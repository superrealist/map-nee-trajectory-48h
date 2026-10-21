"""Source-derived K=3 joint MAP/logNEE Gaussian mixture with AR(1).

No data, fitted parameters, or database access are included.  Inputs are arrays
prepared under the supplement's source and eligibility rules.  This module
retains the published likelihood/optimizer but removes study-specific I/O.
"""
from copy import deepcopy
import os

for _key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[_key] = "1"

import torch
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import linear_sum_assignment, minimize

torch.set_num_threads(1)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DTYPE = torch.float64
WINDOW = 48
K = 3
HOURS = np.arange(WINDOW) + 0.5
NODES = np.array([0.0, 16.0, 32.0, 48.0])
BASIS = CubicSpline(NODES, np.eye(4), bc_type="natural")(HOURS)
TX = torch.as_tensor(BASIS, device=DEVICE, dtype=DTYPE)
MAP_BOUNDS = (46.0, 121.0)
CENTER = np.array([74.34399599811758, -3.435071696071601])
SCALE = np.array([12.534761087258028, 1.3442879939441836])


def transform_hourly(map_values, nee_values):
    """Return N x 48 x 2 standardized values; preserve missingness.

    MAP: mmHg. NEE: microgram/kg/min. NEE=0 remains an observed value.
    No missing-to-zero conversion or interpolation is performed here.
    """
    map_values, nee_values = np.asarray(map_values, float), np.asarray(nee_values, float)
    if map_values.shape != nee_values.shape or map_values.ndim != 2 or map_values.shape[1] != 48:
        raise ValueError("Expected equally sized N x 48 MAP and NEE arrays.")
    if np.any(nee_values[np.isfinite(nee_values)] < 0):
        raise ValueError("Negative NEE is invalid; resolve source records before fitting.")
    y = np.stack([np.clip(map_values, *MAP_BOUNDS), np.log(nee_values + 0.01)], axis=2)
    return (y - CENTER) / SCALE


def tensors(z):
    """Sufficient structures for observed-data likelihood with actual gaps."""
    z = np.asarray(z, float)
    if z.ndim != 3 or z.shape[1:] != (48, 2):
        raise ValueError("z must have shape (N,48,2).")
    answer = []
    n, nt, nm = z.shape
    for marker in range(nm):
        ok = np.isfinite(z[:, :, marker])
        last = np.maximum.accumulate(np.where(ok, np.arange(nt), -1), axis=1)
        prev = np.c_[np.full(n, -1), last[:, :-1]]
        hasprev = (prev >= 0) & ok
        idx = np.maximum(prev, 0)
        gap = np.where(hasprev, np.arange(nt)[None, :] - idx, 1)
        values = dict(y=np.nan_to_num(z[:, :, marker]), ok=ok, hasprev=hasprev, idx=idx, gap=gap)
        answer.append({key: torch.as_tensor(value, device=DEVICE,
                                            dtype=torch.long if key == "idx" else DTYPE)
                       for key, value in values.items()})
    return answer


def vector(model):
    """31 parameters: 24 means, 3 log variances, 2 mixing logits, 2 rhos."""
    return np.r_[np.asarray(model["beta"]).ravel(),
                 np.log(np.asarray(model["var"])[0]),
                 np.log(np.asarray(model["var"])[1, 0]),
                 np.log(np.asarray(model["pi"])[:-1] / model["pi"][-1]),
                 np.asarray(model["rho"])]


def unpack(v):
    beta = v[:24].reshape(2, 3, 4)
    var = torch.stack((v[24:26].exp(), torch.stack((v[26].exp(), v.new_tensor(0.0)))))
    pi = torch.softmax(torch.cat((v[27:29], v.new_zeros(1))), 0)
    return beta, var, pi, v[29:31]


def decode(v):
    values = unpack(torch.as_tensor(v, device=DEVICE, dtype=DTYPE))
    return dict(zip(("beta", "var", "pi", "rho"),
                    (value.detach().cpu().numpy() for value in values)))


def component(v, data):
    """Log pi_k plus marker likelihoods, shape N x 3.

    Cov_m(t,s) = residual_var_m * rho_m**abs(t-s) + random_intercept_var_m.
    The MAP random intercept is integrated analytically. logNEE's random
    intercept variance is fixed to zero. Covariances are common across classes.
    """
    beta, var, pi, rho = unpack(v)
    n = data[0]["y"].shape[0]
    result = torch.log(pi)[None, :].expand(n, 3)
    for marker, a in enumerate(data):
        sig, tau = var[marker]
        mu = TX @ beta[marker].T
        res = (a["y"][:, :, None] - mu[None, :, :]) * a["ok"][:, :, None]
        previous = torch.gather(res, 1, a["idx"][:, :, None].expand(-1, -1, 3))
        r = rho[marker] ** a["gap"] * a["hasprev"]
        den = torch.sqrt(1 - r * r)
        rw = (res - r[:, :, None] * previous) / den[:, :, None] * a["ok"][:, :, None]
        zw = (1 - r) / den * a["ok"]
        count = a["ok"].sum(1)
        h = (zw * zw).sum(1)
        rss = (rw * rw).sum(1)
        rs = (rw * zw[:, :, None]).sum(1)
        ld = torch.log1p(-r * r).sum(1)
        result = result - 0.5 * (count[:, None] * torch.log(2 * torch.pi * sig)
            + ld[:, None] + torch.log1p(h[:, None] * tau / sig) + rss / sig
            - tau * rs * rs / (sig * (sig + h[:, None] * tau)))
    return result


def posterior(z, model, batch_size=600):
    q, ll = [], []
    v = torch.as_tensor(vector(model), device=DEVICE, dtype=DTYPE)
    with torch.no_grad():
        for start in range(0, len(z), batch_size):
            log_component = component(v, tensors(z[start:start + batch_size]))
            q.append(torch.softmax(log_component, 1).cpu().numpy())
            ll.append(torch.logsumexp(log_component, 1).cpu().numpy())
    return np.concatenate(q), np.concatenate(ll)


def align(model, reference):
    """Outcome-blind Hungarian matching of standardized mean curves."""
    a = np.einsum("tp,mkp->ktm", BASIS, reference["beta"])
    b = np.einsum("tp,mkp->ktm", BASIS, model["beta"])
    _, order = linear_sum_assignment(((a[:, None] - b[None]) ** 2).mean((2, 3)))
    answer = deepcopy(model)
    answer["beta"] = answer["beta"][:, order]
    answer["pi"] = answer["pi"][order]
    return answer


def information(data, weights, model):
    w = torch.as_tensor(weights, device=DEVICE, dtype=DTYPE)
    w = w / w.sum()
    loss = lambda v: -(torch.logsumexp(component(v, data), 1) * w).sum()
    v = torch.tensor(vector(model), device=DEVICE, dtype=DTYPE, requires_grad=True)
    return torch.autograd.functional.hessian(loss, v).detach().cpu().numpy()


def estimate(data, weights, initial, preconditioner, precise=False):
    """Published preconditioned BFGS plus Newton refinement.

    preconditioner = chol(inv(information_per_stay)) of the reference used to
    initialize this fit. An explicit initial model is required; labels and
    covariance gates are checked by the caller/accepted_model.
    """
    w = torch.as_tensor(weights, device=DEVICE, dtype=DTYPE)
    if float(w.sum()) <= 0 or torch.any(w < 0):
        raise ValueError("Weights must be nonnegative with positive sum.")
    w = w / w.sum()
    anchor, S = vector(initial), np.asarray(preconditioner)
    rawloss = lambda v: -(torch.logsumexp(component(v, data), 1) * w).sum()

    def fg(u):
        v = anchor + S @ u
        if np.any(v[-2:] <= 0.00001) or np.any(v[-2:] >= 0.995):
            dev = v[-2:] - np.clip(v[-2:], 0.00002, 0.99499)
            gradient = np.zeros(31)
            gradient[-2:] = 2e6 * dev
            return 1e5 + 1e6 * np.sum(dev ** 2), S.T @ gradient
        t = torch.tensor(v, device=DEVICE, dtype=DTYPE, requires_grad=True)
        f = rawloss(t)
        f.backward()
        return float(f.detach().cpu()), S.T @ t.grad.detach().cpu().numpy()

    opt = minimize(fg, np.zeros(31), jac=True, method="BFGS",
                   options={"gtol": 2e-7, "maxiter": 200})
    v = anchor + S @ opt.x
    if precise or abs(opt.jac).max() > 2e-6:
        for _ in range(4):
            t = torch.tensor(v, device=DEVICE, dtype=DTYPE, requires_grad=True)
            f = rawloss(t)
            g = torch.autograd.grad(f, t)[0].detach().cpu().numpy()
            if abs(g).max() < 1e-7:
                break
            h = torch.autograd.functional.hessian(rawloss, t).detach().cpu().numpy()
            step = np.linalg.solve(h + np.eye(31) * max(0, 1e-7 - np.linalg.eigvalsh(h).min()), g)
            alpha = 1.0
            for _ in range(20):
                vn = v - alpha * step
                if (all(vn[-2:] > 0.00001) and all(vn[-2:] < 0.995)
                    and float(rawloss(torch.tensor(vn, device=DEVICE, dtype=DTYPE)).cpu()) < float(f.detach().cpu())):
                    break
                alpha *= 0.5
            if alpha < 1e-5:
                break
            v = vn
    t = torch.tensor(v, device=DEVICE, dtype=DTYPE, requires_grad=True)
    f = rawloss(t)
    g = torch.autograd.grad(f, t)[0].detach().cpu().numpy()
    model = decode(v)
    model.update(mean_loglik=-float(f.detach().cpu()), gradient_max=float(abs(g).max()),
                 scaled_gradient_max=float(abs(S.T @ g).max()),
                 converged=bool(abs(S.T @ g).max() < 2e-5), iterations=int(opt.nit))
    return model


def accepted_model(data, weights, initial, reference, preconditioner, precise=False):
    model = align(estimate(data, weights, initial, preconditioner, precise), reference)
    h = information(data, weights, model)
    minimum = float(np.linalg.eigvalsh(h).min())
    model["information_min_eigenvalue"] = minimum
    model["stage1_gate_passed"] = bool(model["converged"] and np.isfinite(h).all() and minimum > 0)
    return model, h


def observed_class_quality(q):
    """Entropy, class proportions, APP and OCC using hard assigned fractions."""
    q = np.asarray(q, float)
    if q.ndim != 2 or q.shape[1] != 3 or not np.allclose(q.sum(1), 1):
        raise ValueError("q must contain three posterior probabilities per stay.")
    label = q.argmax(1)
    entropy = 1 + np.sum(q * np.log(np.clip(q, 1e-300, 1))) / len(q) / np.log(3)
    rows = []
    for k in range(3):
        selected = label == k
        fraction = selected.mean()
        app = q[selected, k].mean() if selected.any() else np.nan
        occ = (app / (1 - app)) / (fraction / (1 - fraction)) if 0 < fraction < 1 and app < 1 else np.nan
        rows.append(dict(class_id=k + 1, n=int(selected.sum()), fraction=float(fraction),
                         APP=float(app), OCC=float(occ), entropy=float(entropy)))
    return rows


def paired_cluster_weights(mimic_patient, mimic_training, eicu_patient, eicu_hospital,
                           bootstrap_index, seed=2026100917):
    """Integer multiplicities used identically in measurement and association."""
    if bootstrap_index < 1:
        raise ValueError("bootstrap_index starts at one.")
    rng = np.random.default_rng(seed + 100000 + bootstrap_index)
    _, mc = np.unique(np.asarray(mimic_patient).astype(str), return_inverse=True)
    tr = np.asarray(mimic_training, bool)
    tp, vp = np.unique(mc[tr]), np.unique(mc[~tr])
    if np.intersect1d(tp, vp).size:
        raise ValueError("A patient crosses training/holdout partitions.")
    counts = np.zeros(mc.max() + 1, int)
    for group in (tp, vp):
        np.add.at(counts, rng.choice(group, len(group), replace=True), 1)
    _, ec = np.unique(np.asarray(eicu_patient).astype(str), return_inverse=True)
    _, hc = np.unique(np.asarray(eicu_hospital), return_inverse=True)
    ewp = np.bincount(rng.integers(ec.max() + 1, size=ec.max() + 1), minlength=ec.max() + 1)[ec]
    ewh = np.bincount(rng.integers(hc.max() + 1, size=hc.max() + 1), minlength=hc.max() + 1)[hc]
    return dict(mimiciv=counts[mc], eicu=ewp * ewh, eicu_patient=ewp, eicu_hospital=ewh)
