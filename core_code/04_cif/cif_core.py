"""Core Figure 3 estimators. No data or fitted models are included.

Adapted without estimator changes from cohortpub_prepare.py and
local_30_local_curves.py in the canonical publication analysis.
"""
import numpy as np


def cluster_weights(patient_id, B=200, seed=20261009141, hospital_id=None):
    """Patient multinomial weights; crossed patient x hospital for eICU."""
    patient_id = np.asarray(patient_id)
    if patient_id.ndim != 1 or len(patient_id) == 0:
        raise ValueError("One patient identifier is required per stay.")
    rng = np.random.default_rng(seed)
    _, pc = np.unique(patient_id.astype(str), return_inverse=True)
    npat = pc.max() + 1
    if hospital_id is None:
        return rng.multinomial(npat, np.full(npat, 1 / npat), size=B)[:, pc]
    hospital_id = np.asarray(hospital_id)
    if hospital_id.shape != patient_id.shape:
        raise ValueError("Hospital and patient vectors must be aligned.")
    _, hc = np.unique(hospital_id.astype(str), return_inverse=True)
    nh = hc.max() + 1
    wh = rng.multinomial(nh, np.full(nh, 1 / nh), size=B)
    wp = rng.multinomial(npat, np.full(npat, 1 / npat), size=B)
    return wh[:, hc] * wp[:, pc]


def estimators(time, death, weights, grid):
    """Aalen-Johansen hospital-death CIF and discharge-censored KM survival.

    All endpoints must be observed: death=1 is hospital death, death=0 is
    live hospital discharge, NOT right censoring. Follow-up starts at T0+48h.
    For administrative censoring use a general multi-state estimator instead.
    Events at the same time share the risk set immediately before that time.
    """
    time, death, grid = (np.asarray(x, float) for x in (time, death, grid))
    ww = np.atleast_2d(np.asarray(weights, float))
    if (time.ndim != 1 or death.shape != time.shape or not len(time)
            or ww.shape[1] != len(time) or not np.isfinite(time).all()
            or not (time > 0).all() or not np.isin(death, [0, 1]).all()
            or not np.isfinite(ww).all() or (ww < 0).any()
            or not (ww.sum(axis=1) > 0).all()):
        raise ValueError("Invalid endpoints or zero-weight class replicate.")
    ix = np.argsort(time)
    tt, ee, ww = time[ix], death[ix], ww[:, ix]
    unique, first = np.unique(tt, return_index=True)
    total = np.add.reduceat(ww, first, axis=1)
    deaths = np.add.reduceat(ww * ee, first, axis=1)
    risk = ww.sum(axis=1, keepdims=True) - np.c_[
        np.zeros(len(ww)), np.cumsum(total, axis=1)[:, :-1]]
    hd = np.divide(deaths, risk, out=np.zeros_like(deaths), where=risk > 0)
    ha = np.divide(total, risk, out=np.zeros_like(total), where=risk > 0)
    km = np.cumprod(np.clip(1 - hd, 0, 1), axis=1)
    survivor = np.cumprod(np.clip(1 - ha, 0, 1), axis=1)
    prev = np.c_[np.ones(len(ww)), survivor[:, :-1]]
    cif = np.cumsum(prev * hd, axis=1)
    jj = np.searchsorted(unique, grid, side="right") - 1
    return (np.where(jj < 0, 0, cif[:, np.maximum(jj, 0)]),
            np.where(jj < 0, 1, km[:, np.maximum(jj, 0)]), cif[:, -1])


def class_curves(time, death, class_id, patient_id, hospital_id=None,
                 B=200, seed=20261009141, horizon=28):
    """Pointwise percentile intervals conditional on fixed class assignments.

    Use B=200 for MIMIC/eICU and B=500, seed=20261009071 for Local.
    Draw weights once in the full cohort, then subset by class.
    """
    time, death, class_id = map(np.asarray, (time, death, class_id))
    weights = cluster_weights(patient_id, B, seed, hospital_id)
    result = {}
    for k in (1, 2, 3):
        take = class_id == k
        tt, ee = time[take], death[take]
        grid = np.unique(np.r_[0, tt[tt <= horizon], 7, 14, 21, horizon])
        cif, km, terminal = estimators(tt, ee, np.ones(len(tt)), grid)
        boot, bkm, _ = estimators(tt, ee, weights[:, take], grid)
        if not np.isclose(terminal[0], ee.mean()):
            raise AssertionError("Terminal CIF does not match observed death fraction.")
        result[k] = dict(time=grid, CIF=cif[0], KM=km[0],
                         CIF_CI=np.quantile(boot, [.025, .975], axis=0),
                         KM_CI=np.quantile(bkm, [.025, .975], axis=0),
                         at_risk={t: int((tt >= t).sum()) for t in (0, 7, 14, 21, 28)})
    return result
