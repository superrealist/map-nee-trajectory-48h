"""Audited analytic conditional three-step ML; extracted without old run paths."""
import numpy as np
from scipy.optimize import minimize
from scipy.special import logsumexp, expit
from scipy.stats import t

def error_matrix(q):
 lab=q.argmax(1);D=np.stack([q[lab==w].sum(0) for w in range(3)]);D/=q.sum(0)
 assert np.allclose(D.sum(0),1)
 return D

def parts(v,X,y,w,D):
 n,p=X.shape;eta=X@v[:p,None]+np.r_[0,v[p:p+2]][None,:]
 logits=np.c_[X@v[p+2:].reshape(p,2),np.zeros(n)];lp=logits-logsumexp(logits,axis=1,keepdims=True);pi=np.exp(lp);pr=expit(eta)
 joint=lp+np.log(np.maximum(D[w],1e-300))+y[:,None]*eta-np.logaddexp(0,eta);ll=logsumexp(joint,axis=1);r=np.exp(joint-ll[:,None]);er=y[:,None]-pr
 score=np.c_[X*(r*er).sum(1)[:,None],r[:,1:]*er[:,1:],(X[:,:,None]*(r[:,:2]-pi[:,:2])[:,None,:]).reshape(n,2*p)]
 return ll,r,pi,pr,er,score

def fg(v,X,y,w,D):
 ll,*_,S=parts(v,X,y,w,D);return -ll.mean(),-S.mean(0)

def information(v,X,y,w,D):
 ll,r,pi,pr,er,S=parts(v,X,y,w,D);n,p=X.shape;dim=3*p+2;A=np.zeros((dim,dim));vv=pr*(1-pr)
 A[:p,:p]=X.T@(X*(r*vv).sum(1)[:,None]);A[:p,p:p+2]=X.T@(r[:,1:]*vv[:,1:]);A[p:p+2,:p]=A[:p,p:p+2].T;A[p:p+2,p:p+2]=np.diag((r[:,1:]*vv[:,1:]).sum(0))
 for j in range(2):
  for k in range(2):A[p+2+j::2,p+2+k::2]=X.T@(X*(pi[:,j]*((j==k)-pi[:,k]))[:,None])
 for k in range(3):
  delta=np.zeros((n,2))
  if k>0:delta[:,k-1]=er[:,k]
  gate=-pi[:,:2].copy()
  if k<2:gate[:,k]+=1
  s=np.c_[X*er[:,k,None],delta,(X[:,:,None]*gate[:,None,:]).reshape(n,2*p)]
  A-=s.T@(s*r[:,k,None])
 A+=S.T@S
 return (A+A.T)/2,S

def fit(X,y,w,D,start=None):
 p=X.shape[1]
 if start is None:
  v=np.zeros(3*p+2);v[0]=np.log((y.mean()+1e-4)/(1-y.mean()+1e-4))
  gate_counts=np.maximum(np.bincount(w,minlength=3),1);v[p+2:p+4]=np.log(gate_counts[:2]/gate_counts[2])
 else:v=start.copy()
 opt=minimize(fg,v,args=(X,y,w,D),jac=True,method='L-BFGS-B',options={'gtol':2e-7,'ftol':1e-12,'maxiter':800,'maxls':40})
 v=opt.x;A,S=information(v,X,y,w,D);eig=np.linalg.eigvalsh(A/len(y))
 for _ in range(6):
  loss,g=fg(v,X,y,w,D)
  if max(abs(g))<2e-8:break
  step=np.linalg.solve(A/len(y)+np.eye(len(v))*max(0,1e-7-eig.min()),g);alpha=1.
  for _ in range(18):
   trial=v-alpha*step
   if fg(trial,X,y,w,D)[0]<loss:break
   alpha*=.5
  if alpha<1e-5:break
  v=trial;A,S=information(v,X,y,w,D);eig=np.linalg.eigvalsh(A/len(y))
 grad=max(abs(fg(v,X,y,w,D)[1]));return dict(beta=v,A=A,scores=S,converged=bool(grad<1e-6 and eig.min()>0),gradient=grad,eigenvalue=eig.min(),loglik=-fg(v,X,y,w,D)[0]*len(y))

def cluster_cov(result,cluster1,cluster2=None):
 S=result['scores'];n,p=S.shape
 def meat(ids):
  _,inv=np.unique(ids,return_inverse=True);G=inv.max()+1;u=np.zeros((G,p));np.add.at(u,inv,S);return u.T@u*G/(G-1)
 M=meat(cluster1)
 if cluster2 is not None:M+=meat(cluster2)-meat(np.array([str(a)+'|'+str(b) for a,b in zip(cluster1,cluster2)]))
 Ainv=np.linalg.inv(result['A']);return Ainv@M@Ainv*(n-1)/(n-p)

def pool(Q,U,clusters):
 m=len(Q);q=np.mean(Q);u=np.mean(U);b=np.var(Q,ddof=1);total=u+(1+1/m)*b;r=(1+1/m)*b/u;dfold=(m-1)*(1+1/r)**2 if r>0 else np.inf;lam=(1+1/m)*b/total;dfcom=clusters-1;dfobs=(dfcom+1)/(dfcom+3)*dfcom*(1-lam);df=1/(1/dfold+1/dfobs)
 se=np.sqrt(total);cut=t.ppf(.975,df);return dict(OR=np.exp(q),lower=np.exp(q-cut*se),upper=np.exp(q+cut*se),logOR=q,SE=se,df=df,p=2*t.sf(abs(q/se),df),MCSE_over_SE=np.sqrt(b/m)/se)

def robust_fit(X,y,labels,D,start):
 a=fit(X,y,labels,D,start)
 if a['converged']:return a
 # Escape near-stationary saddles along the negative-curvature direction.
 for attempt in range(3):
  val,vec=np.linalg.eigh(a['A']/len(y));direction=vec[:,0];candidates=[]
  for sign in (-1,1):
   b=a['beta']+sign*direction*(1+attempt)
   opt=minimize(fg,b,args=(X,y,labels,D),jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':600})
   z=fit(X,y,labels,D,opt.x);candidates.append(z)
  candidates.append(a);a=max(candidates,key=lambda z:z['loglik'])
  good=[z for z in candidates if z['converged']]
  if good:
   best=max(good,key=lambda z:z['loglik'])
   if best['loglik']>=a['loglik']-1e-5:return best
 # Floating-point zero eigenvalues can occur for a separated nuisance gate.
 # Do not invert that matrix or pretend it supplies a conditional Wald SE.
 # Accept only numerically PSD information with identifiable, stable contrasts;
 # target-specific profiled covariance must be used rather than inverting the full matrix.
 p=X.shape[1];target=np.array([p,p+1]);nuisance=np.setdiff1d(np.arange(3*p+2),target)
 close=[z for z in candidates if z['loglik']>=a['loglik']-1e-5];stable=max(np.max(abs(z['beta'][target]-a['beta'][target]))for z in close)<1e-3
 for z in sorted(close,key=lambda z:z['gradient']):
  H=z['A']/len(y);profile=H[np.ix_(target,target)]-H[np.ix_(target,nuisance)]@np.linalg.pinv(H[np.ix_(nuisance,nuisance)],rcond=1e-10,hermitian=True)@H[np.ix_(nuisance,target)]
  if z['gradient']<1e-7 and z['eigenvalue']>=-1e-10*np.linalg.eigvalsh(H).max() and np.linalg.eigvalsh(profile).min()>1e-5 and stable:
   a=z;a['converged']=True;a['semidefinite_nuisance_gate']=True;break
 return a


def profile_target_covariance(result, p, patient_id, hospital_id=None):
    """Canonical conditional contrast covariance with nuisance-gate profiling.

    p is the column count of X, including its intercept. Returns the 2 x 2
    class-contrast covariance and diagnostics. This conditions on fixed D.
    """
    H = result['A']; S = result['scores']; n = len(patient_id)
    target = np.array([p, p+1])
    nuisance = np.setdiff1d(np.arange(H.shape[0]), target)
    ev = np.linalg.eigvalsh(H)
    relative = ev.min()/ev.max()
    assert relative >= -1e-9
    covariances, minima = [], []
    for tol in (1e-8, 1e-10, 1e-12):
        inv = np.linalg.pinv(H[np.ix_(nuisance,nuisance)], rcond=tol, hermitian=True)
        J = H[np.ix_(target,target)] - H[np.ix_(target,nuisance)] @ inv @ H[np.ix_(nuisance,target)]
        J = (J+J.T)/2
        minimum = np.linalg.eigvalsh(J/n).min()
        minima.append(minimum)
        assert minimum > 1e-5
        scores = S[:,target] - S[:,nuisance] @ inv @ H[np.ix_(nuisance,target)]
        V = cluster_cov({'A':J,'scores':scores},patient_id,hospital_id)
        rank = int((ev > tol*ev.max()).sum())
        V *= float(n-2)/(n-rank)
        covariances.append(V)
    spread = max(np.max(np.abs(np.sqrt(np.diag(v))/np.sqrt(np.diag(covariances[1]))-1)) for v in covariances)
    assert spread < .01 and np.all(np.diag(covariances[1]) > 0)
    return covariances[1], dict(profile_min_eigenvalue=min(minima),
        profile_SE_tolerance_max_relative_change=spread,
        information_min_relative_eigenvalue=relative)
