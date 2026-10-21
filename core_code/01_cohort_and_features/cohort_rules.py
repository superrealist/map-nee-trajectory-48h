"""Core eligibility and observation-quality rules; no database connection.

Inputs are investigator-prepared, source-audited records. These functions do
not assert that supplied local infection/quality flags have been adjudicated.
They do not reproduce upstream extraction or change any published membership.
"""
import numpy as np
import pandas as pd


def official_mimic_by_landmark(frame):
    """Official Sepsis-3 flag and ALL four clinical event times by T0+48h.

    Columns: sepsis3 (boolean), T0 and antibiotic_time, culture_time,
    suspected_infection_time, sofa_time. Times are clinical event times,
    not necessarily the electronic availability times of their records.
    """
    fields = ['antibiotic_time', 'culture_time',
              'suspected_infection_time', 'sofa_time']
    t = frame[fields].apply(pd.to_datetime, errors='coerce')
    landmark = pd.to_datetime(frame['T0']) + pd.Timedelta(hours=48)
    return (frame['sepsis3'].eq(True) & t.notna().all(axis=1)
            & t.le(landmark, axis=0).all(axis=1))


def common_eligibility(frame):
    """Clinical-rule contract for MIMIC and local adjudication.

    The LOCAL release was verified against a supplied final model frame;
    it was not newly selected by this function. Local flags must come from
    the hospital investigators' case-level screening ledger. Unknown flags
    fail the contract; they must not be replaced by True to force agreement.
    """
    flags = ['ICU_eligible', 'T0_supported', 'sepsis_evidence_by48',
             'alive_and_hospitalized48', 'quality48', 'outcome_known',
             'chronology_valid']
    return frame['age'].ge(18) & frame[flags].eq(True).all(axis=1)


def local_final_frame_verification(frame):
    """The exact scope of the current local result-verification gate.

    No new claim of upstream screening is inferred from passing this gate.
    Columns contain supplied-frame status and independently checkable hard
    chronological/age/outcome properties.
    """
    return (frame['supplied_final_frame'].eq(True) & frame['age'].ge(18)
            & frame['T0_available'].eq(True)
            & frame['alive_and_hospitalized48'].eq(True)
            & frame['hospital_mortality'].isin([0, 1]))


def quality48(map_hour, known_minutes, exposure_minutes,
              positive_convertible_first2):
    """Balanced 48h quality rule, normalized to actual ICU observation.

    Three length-48 arrays. known_minutes denotes quantifiable TOTAL NEE.
    exposure_minutes is the observed portion of each hour before ICU exit.
    Trajectory data after ICU exit remain structurally missing.
    """
    m, k, e = map(np.asarray, (map_hour, known_minutes, exposure_minutes))
    if any(v.shape != (48,) for v in (m, k, e)):
        raise ValueError('Exactly 48 hourly values are required')
    if np.any(e < 0) or np.any(e > 60) or np.any(k < 0) or np.any(k > e + 1e-6):
        raise ValueError('Invalid hour-level observation durations')
    valid_nee = (k >= .75 * e - 1e-7) & (e > 0) & (k > 0)
    pair = valid_nee & np.isfinite(m)
    bins = int((pair.reshape(12, 4).sum(axis=1) >= 2).sum())
    coverage = float((pair * e).sum() / e.sum()) if e.sum() else 0.
    longest = run = 0.
    for yes, duration in zip(pair, e):
        if duration <= 0:
            continue
        run = 0. if yes else run + duration / 60
        longest = max(longest, run)
    phase_ok = all(e[s:s+24].sum()/60 < 8 or pair[s:s+24].sum() >= 2
                   for s in (0, 24))
    keep = (bool(positive_convertible_first2) and bins >= 4
            and coverage >= .5 and phase_ok and longest <= 12 + 1e-6)
    return dict(keep=keep, paired_hours=int(pair.sum()), adequate_4h_bins=bins,
                paired_fraction=coverage, longest_unpaired_hours=longest,
                phase_coverage_ok=bool(phase_ok), reliable_nee=valid_nee)


def eicu_sustained_T0(records):
    """First positive record with next same-agent positive record <=120min.

    Caller resolves same-agent/same-time duplicates according to original
    unit/source priorities. Required columns: drug_name, offset_min, rate.
    Offset is ICU-admission minutes; records must already be restricted to
    the index ICU stay. This returns the FIRST record's time, not the later
    confirmation record's time, and does not establish true clinical onset.
    """
    a = records.sort_values(['drug_name', 'offset_min']).copy()
    by = a.groupby('drug_name', sort=False)
    later = by['offset_min'].shift(-1)
    positive_next = by['rate'].shift(-1).gt(0)
    ok = a['rate'].gt(0) & positive_next & (later-a['offset_min']).le(120)
    return a.loc[ok, 'offset_min'].min()
