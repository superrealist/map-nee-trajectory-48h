-- Core raw-source queries for MIMIC-IV v3.1 / PostgreSQL.
-- Named blocks are executed separately by extract.py.
-- Event queries preserve source values; final time/quality rules are applied later.

-- name: cohort_candidates
SELECT i.subject_id, i.hadm_id, i.stay_id AS icu_stay_id,
       i.intime, i.outtime, i.first_careunit AS icu_type,
       a.age, p.gender AS sex, h.race, h.admission_type,
       h.admission_location, h.admittime, h.dischtime,
       h.deathtime, h.hospital_expire_flag
FROM mimiciv_icu.icustays i
JOIN mimiciv_derived.age a ON a.hadm_id = i.hadm_id
JOIN mimiciv_hosp.patients p ON p.subject_id = i.subject_id
JOIN mimiciv_hosp.admissions h ON h.hadm_id = i.hadm_id
WHERE a.age >= 18;

-- name: sepsis3
SELECT s.*
FROM mimiciv_derived.sepsis3 s
WHERE s.stay_id = ANY(%(stay_ids)s::bigint[]) AND s.sepsis3 IS TRUE;

-- name: vasopressors_raw
SELECT i.subject_id, i.hadm_id, i.stay_id AS icu_stay_id,
       e.orderid, e.linkorderid, e.itemid,
       CASE e.itemid WHEN 221906 THEN 'norepinephrine'
         WHEN 221289 THEN 'epinephrine' WHEN 222315 THEN 'vasopressin'
         WHEN 221662 THEN 'dopamine' WHEN 221749 THEN 'phenylephrine' END AS drug_name,
       e.starttime, e.endtime, e.storetime, e.rate, e.rateuom,
       e.originalrate, e.amount, e.amountuom, e.patientweight,
       e.statusdescription, e.ordercategoryname, e.ordercategorydescription
FROM mimiciv_icu.icustays i
JOIN mimiciv_icu.inputevents e ON e.stay_id = i.stay_id
WHERE i.stay_id = ANY(%(stay_ids)s::bigint[])
  AND e.itemid IN (221906,221289,222315,221662,221749)
  AND e.starttime < i.outtime AND e.endtime > i.intime;

-- name: chart_events
SELECT c.subject_id, c.hadm_id, c.stay_id AS icu_stay_id,
       c.itemid, d.label, c.charttime, c.storetime,
       c.value, c.valuenum, c.valueuom, c.warning
FROM mimiciv_icu.chartevents c
JOIN mimiciv_icu.d_items d ON d.itemid = c.itemid
WHERE c.stay_id = ANY(%(stay_ids)s::bigint[])
  AND c.itemid IN (220052,220181,220050,220179,220051,220180,
                  220045,220210,224690,220277,223761,223762,
                  224639,226512,226531,223835,226732,223849,
                  220739,223900,223901);

-- name: laboratory_events
SELECT i.stay_id AS icu_stay_id, l.labevent_id, l.specimen_id,
       l.subject_id, l.hadm_id, l.itemid, d.label, d.fluid,
       l.charttime, l.storetime, l.value, l.valuenum, l.valueuom,
       CASE WHEN l.hadm_id = i.hadm_id THEN 'same_hadm'
            ELSE 'null_hadm_subject_window' END AS admission_link
FROM mimiciv_icu.icustays i
JOIN mimiciv_hosp.admissions a ON a.hadm_id = i.hadm_id
JOIN mimiciv_hosp.labevents l ON l.subject_id = i.subject_id
 AND (l.hadm_id = i.hadm_id OR
      (l.hadm_id IS NULL AND l.charttime >= a.admittime - INTERVAL '6 hours'))
JOIN mimiciv_hosp.d_labitems d ON d.itemid = l.itemid
WHERE i.stay_id = ANY(%(stay_ids)s::bigint[])
  AND l.charttime >= a.admittime - INTERVAL '6 hours'
  AND l.charttime < i.outtime
  AND l.itemid IN (50813,52442,53154,50912,52546,50885,53089,
                  51265,53189,51301,51300,51222,50811,50983,50824,
                  50820,52022,50821,52033,50816,50862);

-- name: fluid_input_raw
SELECT e.subject_id, e.hadm_id, e.stay_id AS icu_stay_id,
       e.orderid, e.itemid, d.label, d.category,
       e.starttime, e.endtime, e.storetime, e.amount, e.amountuom,
       e.rate, e.rateuom, e.statusdescription, e.ordercategorydescription
FROM mimiciv_icu.inputevents e
JOIN mimiciv_icu.d_items d ON d.itemid = e.itemid
WHERE e.stay_id = ANY(%(stay_ids)s::bigint[]);

-- name: fluid_output_raw
SELECT e.subject_id, e.hadm_id, e.stay_id AS icu_stay_id,
       e.itemid, d.label, e.charttime, e.storetime, e.value, e.valueuom
FROM mimiciv_icu.outputevents e
JOIN mimiciv_icu.d_items d ON d.itemid = e.itemid
WHERE e.stay_id = ANY(%(stay_ids)s::bigint[]);

-- name: diagnoses
SELECT i.stay_id AS icu_stay_id, d.hadm_id, d.seq_num,
       d.icd_code, d.icd_version, c.long_title
FROM mimiciv_icu.icustays i
JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id = i.hadm_id
JOIN mimiciv_hosp.d_icd_diagnoses c
  ON c.icd_code = d.icd_code AND c.icd_version = d.icd_version
WHERE i.stay_id = ANY(%(stay_ids)s::bigint[]);

-- name: charlson
SELECT i.stay_id AS icu_stay_id, c.*
FROM mimiciv_icu.icustays i
JOIN mimiciv_derived.charlson c ON c.hadm_id = i.hadm_id
WHERE i.stay_id = ANY(%(stay_ids)s::bigint[]);

-- name: sofa_intervals
SELECT s.*
FROM mimiciv_derived.sofa s
WHERE s.stay_id = ANY(%(stay_ids)s::bigint[]);

-- name: ventilation_intervals
SELECT v.*
FROM mimiciv_derived.ventilation v
WHERE v.stay_id = ANY(%(stay_ids)s::bigint[]);
