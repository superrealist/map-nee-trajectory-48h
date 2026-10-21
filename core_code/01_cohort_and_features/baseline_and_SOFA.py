"""Core bounded baseline selection and explicitly labelled SOFA proxies.

Organ-specific source preparation (arterial gas verification, GCS convention,
urine coverage, drug thresholds) is described in the supplement. This module
does not replace those checks with an unqualified sum of laboratory columns.
"""
import numpy as np
import pandas as pd


def latest_T0_lab(events):
    """Latest valid event in [-24,0]h, result/revision available by T0.

    Caller has standardized units and physiological/source validity.
    Columns: stay_id, variable, value, event_h, available_h, source_row.
    """
    a = events[events['event_h'].between(-24, 0)
               & events['available_h'].le(0)].copy()
    return (a.sort_values(['event_h', 'available_h', 'source_row'])
             .drop_duplicates(['stay_id', 'variable'], keep='last'))


def retrospective_T0_vital_mean(events):
    """Mean pre-T0 physiology documented before the 48h landmark.

    This is retrospective initial-state information; it is NOT a guarantee
    that every record was electronically available at the instant T0.
    Caller resolves duplicate channels/times and physiological validity.
    """
    a = events[events['event_h'].between(-24, 0)
               & events['available_h'].lt(48)].copy()
    return a.groupby(['stay_id', 'variable'])['value'].mean()


def score_from_observed_value(value, organ):
    """Single-value SOFA domains; use worst qualified score per window."""
    x = np.asarray(value, dtype=float)
    good = np.isfinite(x)
    if organ == 'coagulation':
        good &= x > 0
        score = np.select([x < 20, x < 50, x < 100, x < 150], [4, 3, 2, 1], 0)
    elif organ == 'liver':
        good &= x >= 0
        score = np.select([x >= 12, x >= 6, x >= 2, x >= 1.2], [4, 3, 2, 1], 0)
    elif organ == 'renal_creatinine':
        good &= x > 0
        score = np.select([x >= 5, x >= 3.5, x >= 2, x >= 1.2], [4, 3, 2, 1], 0)
    elif organ == 'CNS_GCS':
        good &= (x >= 3) & (x <= 15)
        score = np.select([x < 6, x < 10, x < 13, x < 15], [4, 3, 2, 1], 0)
    else:
        raise ValueError('Unsupported domain; respiratory/CV require source context')
    return np.where(good, score, np.nan)


def proxy_SOFA(observed_domains, known_CV_unscorable=None):
    """Latest publication convention; it explicitly defaults absent domains.

    Columns: respiratory, coagulation, liver, CNS, renal, cardiovascular.
    A non-CV score needs >=1 supported domain; a full score also remains
    missing for known positive cardiovascular support that cannot be graded.
    Domain-complete alternative sums and support counts remain distinguishable.
    """
    noncv = ['respiratory', 'coagulation', 'liver', 'CNS', 'renal']
    full = noncv + ['cardiovascular']
    x = observed_domains[full].apply(pd.to_numeric, errors='coerce')
    if np.any((x.to_numpy() < 0) | (x.to_numpy() > 4)):
        raise ValueError('Each supported SOFA domain must range from 0 to 4')
    n5, n6 = x[noncv].notna().sum(axis=1), x.notna().sum(axis=1)
    block = pd.Series(False, index=x.index) if known_CV_unscorable is None else pd.Series(known_CV_unscorable, index=x.index).fillna(True)
    return pd.DataFrame({
        'nonCV_proxy': x[noncv].fillna(0).sum(axis=1).where(n5 > 0),
        'full_proxy': x.fillna(0).sum(axis=1).where((n6 > 0) & ~block),
        'nonCV_complete': x[noncv].sum(axis=1, min_count=5),
        'full_complete': x.sum(axis=1, min_count=6),
        'nonCV_domains_supported': n5,
        'full_domains_supported': n6,
        'nonCV_default_normal_domains': 5-n5,
        'full_default_normal_domains': 6-n6,
    })
