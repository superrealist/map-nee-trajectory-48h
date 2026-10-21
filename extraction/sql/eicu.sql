-- Core raw-source queries for eICU-CRD v2.0 / PostgreSQL.
-- Offsets are minutes from the index ICU admission, not timestamps.

-- name: cohort_candidates
SELECT p.uniquepid AS patient_id, p.patientunitstayid AS icu_stay_id,
       p.patienthealthsystemstayid, p.hospitalid, p.unittype AS icu_type,
       p.age AS age_raw,
       CASE WHEN p.age IN ('> 89','>89') THEN 90
            WHEN p.age ~ '^[0-9]+$' THEN p.age::integer END AS age,
       p.gender AS sex, p.ethnicity AS race, p.admissionweight,
       p.unitadmitsource, p.hospitaladmitsource, p.hospitaldischargeyear,
       p.hospitaladmitoffset, p.hospitaldischargeoffset,
       p.unitdischargeoffset, p.hospitaldischargestatus, p.unitdischargestatus
FROM eicu_crd.patient p
WHERE CASE WHEN p.age IN ('> 89','>89') THEN 90
           WHEN p.age ~ '^[0-9]+$' THEN p.age::integer END >= 18;

-- name: vasopressors_raw
SELECT d.patientunitstayid AS icu_stay_id, d.infusiondrugid,
       d.infusionoffset, d.drugname, d.drugrate, d.infusionrate,
       d.drugamount, d.volumeoffluid, d.patientweight
FROM eicu_crd.infusiondrug d
JOIN eicu_crd.patient p ON p.patientunitstayid = d.patientunitstayid
WHERE d.patientunitstayid = ANY(%(stay_ids)s::bigint[])
  AND d.infusionoffset >= 0 AND d.infusionoffset < p.unitdischargeoffset
  AND lower(d.drugname) ~
      'norepinephrine|levophed|epinephrine|adrenalin|vasopressin|dopamine|phenylephrine|neosynephrine|neo-synephrine';

-- name: vital_periodic
SELECT v.patientunitstayid AS icu_stay_id, v.vitalperiodicid,
       v.observationoffset, v.systemicmean, v.systemicsystolic,
       v.systemicdiastolic, v.heartrate, v.respiration, v.sao2, v.temperature
FROM eicu_crd.vitalperiodic v
WHERE v.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: laboratory_events
SELECT l.patientunitstayid AS icu_stay_id, l.labid, l.labname,
       l.labresultoffset, l.labresultrevisedoffset, l.labresult,
       l.labmeasurenamesystem, l.labmeasurenameinterface
FROM eicu_crd.lab l
WHERE l.patientunitstayid = ANY(%(stay_ids)s::bigint[])
  AND lower(l.labname) IN ('lactate','creatinine','total bilirubin',
      'platelets x 1000','wbc x 1000','hgb','sodium','ph','pao2','fio2','albumin');

-- name: fluid_balance_raw
SELECT u.patientunitstayid AS icu_stay_id, u.intakeoutputid,
       u.intakeoutputoffset, u.intakeoutputentryoffset, u.cellpath,
       u.celllabel, u.cellvaluenumeric, u.cellvaluetext
FROM eicu_crd.intakeoutput u
WHERE u.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: diagnoses
SELECT d.patientunitstayid AS icu_stay_id, d.diagnosisid,
       d.diagnosisoffset, d.diagnosisstring, d.icd9code, d.diagnosispriority
FROM eicu_crd.diagnosis d
WHERE d.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: history
SELECT h.patientunitstayid AS icu_stay_id, h.pasthistoryid,
       h.pasthistoryoffset, h.pasthistorypath, h.pasthistoryvalue
FROM eicu_crd.pasthistory h
WHERE h.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: medications
SELECT m.patientunitstayid AS icu_stay_id, m.medicationid,
       m.drugorderoffset, m.drugstartoffset, m.drugstopoffset,
       m.drugname, m.dosage, m.routeadmin, m.drugordercancelled
FROM eicu_crd.medication m
WHERE m.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: microbiology
SELECT m.patientunitstayid AS icu_stay_id, m.microlabid,
       m.culturetakenoffset, m.culturesite, m.organism
FROM eicu_crd.microlab m
WHERE m.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: neurologic_chart
SELECT n.patientunitstayid AS icu_stay_id, n.nursingchartoffset,
       n.nursingchartentryoffset, n.nursingchartcelltypevallabel,
       n.nursingchartcelltypevalname, n.nursingchartvalue
FROM eicu_crd.nursecharting n
WHERE n.patientunitstayid = ANY(%(stay_ids)s::bigint[])
  AND (lower(n.nursingchartcelltypevallabel) LIKE '%%glasgow%%'
       OR lower(n.nursingchartcelltypevalname) LIKE '%%gcs%%');

-- name: respiratory_chart
SELECT r.patientunitstayid AS icu_stay_id, r.respchartoffset,
       r.respchartentryoffset, r.respchartvaluelabel, r.respchartvalue
FROM eicu_crd.respiratorycharting r
WHERE r.patientunitstayid = ANY(%(stay_ids)s::bigint[])
  AND r.respchartvaluelabel IN ('FiO2','RT Vent On/Off','Ventilator Mode',
      'Vent Mode','O2 Device','Oxygen Delivery Method','Airway Type');

-- name: respiratory_care
SELECT r.patientunitstayid AS icu_stay_id, r.respcarestatusoffset,
       r.airwaytype, r.ventstartoffset, r.ventendoffset
FROM eicu_crd.respiratorycare r
WHERE r.patientunitstayid = ANY(%(stay_ids)s::bigint[]);

-- name: treatments
SELECT t.patientunitstayid AS icu_stay_id, t.treatmentid,
       t.treatmentoffset, t.treatmentstring
FROM eicu_crd.treatment t
WHERE t.patientunitstayid = ANY(%(stay_ids)s::bigint[]);
