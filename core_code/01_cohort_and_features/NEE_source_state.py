"""Core source-state reconstruction, extracted from the locked 48h pipeline.

API: records use relative MINUTES; length is min(ICU exit-T0, 2880).
Returns known_minutes, dose_integral, zero_minutes, positive_minutes.
For the 48h analysis slice the returned arrays to [:48].

MIMIC records are (start, end, drug_index, NEE_contribution, end_status).
Pre-filter: valid Continuous Med, positive rate, duration >1min, not
Rewritten/Bolus, deduplicate exact intervals; unknown converted contribution
is NaN. Source-qualified NEE contribution >5 is withheld pending verification.

EICU records are (offset, drug_index, NEE_contribution), unique per agent/time;
conflicting simultaneous doses become NaN. Main: cap=120, stop_forever=True.
Positive dose expires at 2h; explicit zero persists until next record.
Never-documented agents initially equal zero: stream-completeness assumption.
The bounded-zero comparator uses stop_forever=False; it is not the main rule.

Hourly mean = area/known, only if known>=0.75*ICU-exposure minutes and known>0.
No default-zero replacement is applied to UNKNOWN total-dose intervals.
No raw records, model objects, database endpoints or credentials are included.
"""
import re
import math
import numpy as np
DRUGS=['norepinephrine','epinephrine','vasopressin','dopamine','phenylephrine']
FACT=[1,1,2.5,.01,.1]
def drug(s):
 s=str(s).lower()
 for k,p in enumerate(['norepinephrine|levophed','epinephrine|adrenalin','vasopressin','dopamine','phenylephrine|neosynephrine|neo-synephrine']):
  if re.search(p,s):return k
 return -1

def num(x):
 try:return float(x)
 except:return np.nan

def unit(s):
 s=str(s).lower().replace('μg','mcg').replace('ug/','mcg/').replace('hour','hr').replace('units','unit')
 for x in ['mcg/kg/min','mcg/kg/hr','mcg/min','mcg/hr','mg/kg/min','mg/hr','mg/min','unit/kg/min','unit/min','unit/hr','ml/hr']:
  if x in s:return x
 return 'unknown'

def dose(r,u,k,weight):
 if not np.isfinite(r) or r<0:return np.nan
 if r==0:return 0.
 if k==2:
  if u=='unit/min':return r*2.5
  if u=='unit/hr':return r*2.5/60
  if u=='unit/kg/min' and np.isfinite(weight):return r*weight*2.5
  return np.nan
 if u=='mcg/kg/min':return r*FACT[k]
 if u=='mcg/kg/hr':return r*FACT[k]/60
 if not np.isfinite(weight) or not 20<=weight<=400:return np.nan
 z={'mcg/min':1,'mcg/hr':1/60,'mg/min':1000,'mg/hr':1000/60}.get(u,np.nan)
 return r*z/weight*FACT[k]

def bag_conc(s,k):
 # Only explicit dilute infusion bags, not injection-vial strength (e.g. 1 mg/mL).
 s=str(s).lower()
 vals=[]
 for m in re.finditer(r'(\d+(?:\.\d+)?)\s*(mg|mcg|units?|u)\s*(?:in\s*|/\s*)(\d+(?:\.\d+)?)\s*ml',s):
  a,u,v=float(m[1]),m[2],float(m[3])
  if v<50 or a<=0:continue
  if k==2 and u in ['unit','units','u']:vals.append(a/v)
  if k!=2 and u in ['mg','mcg']:vals.append(a/v*(1000 if u=='mg' else 1))
 return vals[0] if len(set(vals))==1 else np.nan

def integrate(start,end,value,length,known,area,zero,positive):
 a=max(0.,start);b=min(length,end)
 if b<=a or not np.isfinite(value):return
 h0=max(0,int(a//60));h1=min(72,int(math.ceil(b/60)))
 hh=np.arange(h0,h1);dur=np.minimum(b,(hh+1)*60)-np.maximum(a,hh*60)
 known[hh]+=dur;area[hh]+=dur*value
 if value==0:zero[hh]+=dur
 else:positive[hh]+=dur

def sweep_mimic(records,length):
 # End status applies to this segment, not all drugs. Track overlapping same-drug records.
 events={};ever=False
 for j,(s,e,k,v,status) in enumerate(records):
  if e<=s or e<=0 or s>=length:continue
  events.setdefault(s,[]).append(('start',k,j,v,status))
  events.setdefault(e,[]).append(('end',k,j,v,status))
 active=[{} for _ in DRUGS];rest=np.zeros(5);known=np.zeros(72);area=known.copy();zero=known.copy();pos=known.copy();last=min([0.]+list(events))
 for t in sorted(set(events)|{0.,float(length)}):
  values=rest.copy()
  for k,d in enumerate(active):
   if d:
    z=np.array(list(d.values()));values[k]=z[0] if np.isfinite(z).all() and np.allclose(z,z[0],rtol=1e-6,atol=1e-10) else np.nan
  if ever:integrate(last,t,values.sum(),length,known,area,zero,pos)
  ev=events.get(t,[])
  for kind,k,j,v,status in ev:
   if kind=='end':active[k].pop(j,None);rest[k]=0. if status in ['Stopped','Paused','explicit_zero'] else np.nan
  for kind,k,j,v,status in ev:
   if kind=='start':active[k][j]=v;ever=True
  last=t
 return known,area,zero,pos

def sweep_eicu(records,length,cap,stop_forever=False):
 events={};versions=np.zeros(5,int);state=np.zeros(5);ever=False
 for j,(s,k,v) in enumerate(records):
  if s>=length:continue
  events.setdefault(s,[]).append(('new',k,j,v))
  expiry=s+(max(length-s,0)+1 if stop_forever and v==0 else cap)
  events.setdefault(expiry,[]).append(('expiry',k,j,v))
 known=np.zeros(72);area=known.copy();zero=known.copy();pos=known.copy();last=min([0.]+list(events))
 for t in sorted(set(events)|{0.,float(length)}):
  if ever:integrate(last,t,state.sum(),length,known,area,zero,pos)
  ev=events.get(t,[])
  for kind,k,j,v in ev:
   if kind=='expiry' and versions[k]==j+1:state[k]=np.nan
  for kind,k,j,v in ev:
   if kind=='new':state[k]=v;versions[k]=j+1;ever=True
  last=t
 return known,area,zero,pos
