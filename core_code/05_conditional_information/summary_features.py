"""Observed 0-48 h summaries for Figure 4; no smoothing or imputation."""
import numpy as np


def summary_features(map_values, nee_values):
    """Inputs: stay x 48 matrices; unknown hours remain NaN.

    Uses the production MIMIC/eICU rule: >=3 observed hours and >=4 h span
    for slope; sample SD requires >=3. Local data have all 48 hours, so its
    original >=2-hour slope/SD implementation gives the same estimates.
    """
    maps, nee = np.asarray(map_values, float), np.asarray(nee_values, float)
    if maps.ndim != 2 or maps.shape[1] != 48 or nee.shape != maps.shape:
        raise ValueError("Aligned stay x 48 matrices required.")
    if (nee[np.isfinite(nee)] < 0).any():
        raise ValueError("NEE must be non-negative when recorded.")
    time = np.arange(48)[None, :] + .5
    out = {}
    for name, values in [('MAP', maps), ('logNEE', np.log(nee + .01))]:
        observed = np.isfinite(values)
        count = observed.sum(1)
        mean_time = np.divide((observed * time).sum(1), count,
                              out=np.zeros(len(values)), where=count > 0)
        centered = (time - mean_time[:, None]) * observed
        variance = (centered ** 2).sum(1)
        slope = np.divide((centered * np.nan_to_num(values)).sum(1), variance,
                          out=np.full(len(values), np.nan), where=variance > 0) * 24
        span = np.max(np.where(observed, time, -np.inf), 1) - np.min(np.where(observed, time, np.inf), 1)
        slope[(count < 3) | (span < 4)] = np.nan
        sd = np.nanstd(values, axis=1, ddof=1);sd[count < 3] = np.nan
        out[name + '_mean'] = np.nanmean(values, axis=1)
        out[name + '_slope24'] = slope
        out[name + '_SD'] = sd
    denom = np.isfinite(nee).sum(1)
    out['NEE_positive_fraction'] = np.divide(((nee > 0) & np.isfinite(nee)).sum(1), denom,
                                            out=np.full(len(nee), np.nan), where=denom > 0)
    return out
