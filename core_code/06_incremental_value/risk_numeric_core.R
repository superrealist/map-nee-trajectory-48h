# Pure numerical functions extracted from the validated risk implementation.
# No source data, fitted coefficients, or filesystem paths are included.
ridge<-function(A,y){
 fit<-glmnet(A[,-1,drop=FALSE],y,family='binomial',alpha=0,lambda=1/length(y),standardize=FALSE,intercept=TRUE,thresh=1e-13,maxit=1000000);stopifnot(fit$jerr==0)
 b<-as.numeric(coef(fit));pen<-c(0,rep(1,ncol(A)-1));loss<-function(b){e<-drop(A%*%b);sum(pmax(e,0)+log1p(exp(-abs(e)))-y*e)+.5*sum(pen*b*b)}
 for(it in 1:20){p<-plogis(drop(A%*%b));sc<-drop(crossprod(A,p-y))+pen*b;gr<-max(abs(sc))/length(y);if(gr<1e-8)break;H<-crossprod(A,A*(p*(1-p)))+diag(pen);st<-solve(H,sc);al<-1;ll<-loss(b);for(j in 1:30){bb<-b-al*st;if(loss(bb)<ll)break;al<-al/2};stopifnot(al>1e-9);b<-bb}
 p<-plogis(drop(A%*%b));gr<-max(abs(drop(crossprod(A,p-y))+pen*b))/length(y);stopifnot(gr<1e-7,all(is.finite(b)))
 HH<-crossprod(A,A*(p*(1-p)));edf<-sum(diag(solve(HH+diag(pen),HH)))
 list(coefficients=b,gradient=gr,lambda=1,effective_df_including_intercept=edf,iterations=it)
}
auc<-function(y,p,w){o<-order(p);y<-y[o];p<-p[o];w<-w[o];g<-cumsum(c(TRUE,diff(p)!=0));a<-as.numeric(rowsum(w*y,g,reorder=FALSE));b<-as.numeric(rowsum(w*(1-y),g,reorder=FALSE));if(sum(a)==0||sum(b)==0)return(NA_real_);sum(a*(cumsum(b)-.5*b))/(sum(a)*sum(b))}
cal<-function(y,p,w){lp<-qlogis(pmin(pmax(p,1e-9),1-1e-9));b<-c(0,1);for(j in 1:80){q<-plogis(b[1]+b[2]*lp);v<-w*q*(1-q);g<-c(sum(w*(y-q)),sum(w*(y-q)*lp));H<-matrix(c(sum(v),sum(v*lp),sum(v*lp),sum(v*lp^2)),2);st<-solve(H,g);b<-b+st;if(max(abs(st))<1e-8)break};stopifnot(j<80);a<-0;for(j in 1:80){q<-plogis(a+lp);st<-sum(w*(y-q))/sum(w*q*(1-q));a<-a+st;if(abs(st)<1e-8)break};stopifnot(j<80);c(CITL=a,slope=b[2])}
metric<-function(y,p,w,doc=TRUE){a<-c(AUC=auc(y,p,w),Brier=sum(w*(y-p)^2)/sum(w));if(doc)a<-c(a,cal(y,p,w));a}
