# Figure 5 core: train once in MIMIC, freeze preprocessing and all coefficients,
# then evaluate the same paired models in each evaluation cohort.
# Source risk_numeric_core.R and attach glmnet before using these functions.

fit_numeric_block <- function(v, name, bounds=NULL) {
  v <- as.numeric(v); observed <- v[is.finite(v)]
  stopifnot(length(observed)>0)
  bd <- if(is.null(bounds))range(observed) else bounds
  kn <- as.numeric(quantile(observed,c(1/3,2/3)))
  kn <- unique(kn[kn>bd[1]&kn<bd[2]])
  if(length(kn)!=2L)kn <- numeric()
  list(name=name,median=median(observed),knots=kn,bounds=bd)
}
numeric_block <- function(v, spec) {
  v <- as.numeric(v); missing <- !is.finite(v); v[missing] <- spec$median
  a <- if(length(spec$knots)==2L)splines::ns(v,knots=spec$knots,Boundary.knots=spec$bounds) else matrix(v,ncol=1)
  colnames(a) <- paste0(spec$name,'_s',seq_len(ncol(a)))
  a <- cbind(a,as.numeric(missing)); colnames(a)[ncol(a)] <- paste0(spec$name,'_missing'); a
}
fit_categorical_block <- function(v,name) {
  v <- as.character(v); v[is.na(v)|!nzchar(v)] <- 'Unknown'
  list(name=name,levels=sort(unique(v)),modal=names(sort(table(v),decreasing=TRUE))[1])
}
categorical_block <- function(v,spec) {
  v <- as.character(v); v[is.na(v)|!nzchar(v)] <- 'Unknown'
  v[!v%in%spec$levels] <- if('Unknown'%in%spec$levels)'Unknown' else spec$modal
  if(length(spec$levels)<2L)return(matrix(numeric(),nrow=length(v),ncol=0))
  a <- model.matrix(~factor(v,levels=spec$levels))
  colnames(a) <- c('(Intercept)',paste0(spec$name,'_',spec$levels[-1]));a[,-1,drop=FALSE]
}
risk_inputs <- function(d) {
  required <- c('age','sex','icu','MAP','HR','lactate','NEE','joint_p2','joint_p3',
                paste0(rep(c('S0','S1','S2'),each=2),'_',rep(c('full','nonCV'),3),'_proxy'))
  stopifnot(all(required%in%names(d)),all(is.na(d$NEE)|d$NEE>=0),all(is.na(d$lactate)|d$lactate>=0))
  x <- d; x$log1p_lactate <- log1p(d$lactate);x$logNEE <- log(d$NEE+.01)
  x$clinical_nonCV <- d$S0_nonCV_proxy
  J <- as.matrix(d[,c('joint_p2','joint_p3')]);stopifnot(all(is.finite(J)),all(J>=0),all(rowSums(J)<=1+1e-10))
  list(x=x,J=J)
}
raw_risk_blocks <- function(x,spec) {
  nb <- function(v,nm)numeric_block(v,spec$numeric[[nm]])
  cb <- function(v,nm)categorical_block(v,spec$categorical[[nm]])
  out <- list(T0=cbind(nb(x$age,'age'),cb(x$sex,'sex'),cb(x$icu,'ICU'),
    nb(x$MAP,'MAP'),nb(x$HR,'HR'),nb(x$log1p_lactate,'log1p_lactate'),
    nb(x$logNEE,'logNEE'),nb(x$clinical_nonCV,'clinical_nonCV')))
  for(kind in c('full','nonCV')) {
    z <- lapply(c('S0','S1','S2'),function(wn)nb(x[[paste0(wn,'_',kind,'_proxy')]],paste0(kind,'_',wn)))
    out[[paste0('static_',kind)]] <- z[[1]]
    out[[paste0('dynamic_',kind)]] <- do.call(cbind,z)
  };out
}
fit_pair_scaling <- function(raw,J) {
  a <- cbind('(Intercept)'=1,raw);stopifnot(all(is.finite(a)),!anyDuplicated(colnames(a)))
  mu <- colMeans(raw);ss <- apply(raw,2,sd);ok <- ss>1e-10
  a <- a[,c(TRUE,ok),drop=FALSE];mu <- mu[ok];ss <- ss[ok]
  a[,-1] <- sweep(sweep(a[,-1,drop=FALSE],2,mu,'-'),2,ss,'/')
  retained <- 1L
  for(j in seq.int(2L,ncol(a)))if(qr(a[,c(retained,j),drop=FALSE],tol=1e-10)$rank==length(retained)+1L)retained<-c(retained,j)
  a <- a[,retained,drop=FALSE];js <- apply(J,2,sd);stopifnot(all(js>1e-10))
  list(columns=colnames(a),center=mu[colnames(a)[-1]],scale=ss[colnames(a)[-1]],
       joint_center=colMeans(J),joint_scale=js)
}
apply_pair_scaling <- function(raw,J,spec) {
  a <- cbind('(Intercept)'=1,raw)[,spec$columns,drop=FALSE]
  a[,-1] <- sweep(sweep(a[,-1,drop=FALSE],2,spec$center,'-'),2,spec$scale,'/')
  jj <- sweep(sweep(J,2,spec$joint_center,'-'),2,spec$joint_scale,'/')
  list(a=a,b=cbind(a,jj))
}
fit_risk_models <- function(training) {
  z <- risk_inputs(training);x<-z$x;J<-z$J;y<-training$hospital_mortality
  stopifnot(all(y%in%c(0,1)),length(unique(y))==2L)
  spec <- list(numeric=list(),categorical=list())
  for(nm in c('age','MAP','HR','log1p_lactate','logNEE','clinical_nonCV'))
    spec$numeric[[nm]] <- fit_numeric_block(x[[nm]],nm,if(nm=='clinical_nonCV')c(0,20)else NULL)
  for(kind in c('full','nonCV'))for(wn in c('S0','S1','S2')) {
    nm <- paste0(kind,'_',wn)
    spec$numeric[[nm]] <- fit_numeric_block(x[[paste0(wn,'_',kind,'_proxy')]],nm,if(kind=='full')c(0,24)else c(0,20))
  }
  spec$categorical$sex <- fit_categorical_block(x$sex,'sex')
  spec$categorical$ICU <- fit_categorical_block(x$icu,'ICU')
  raw <- raw_risk_blocks(x,spec);scaling<-fits<-list()
  for(key in names(raw)) {
    scaling[[key]] <- fit_pair_scaling(raw[[key]],J)
    a <- apply_pair_scaling(raw[[key]],J,scaling[[key]])
    for(ver in c('a','b')) {
      stopifnot(qr(a[[ver]])$rank==ncol(a[[ver]]))
      fits[[paste0(key,ver)]] <- ridge(a[[ver]],y)
    }
  }
  list(preprocessing=spec,scaling=scaling,models=fits)
}
predict_risk_models <- function(fit,evaluation) {
  z <- risk_inputs(evaluation);raw<-raw_risk_blocks(z$x,fit$preprocessing);out<-list()
  for(key in names(raw)) {
    a <- apply_pair_scaling(raw[[key]],z$J,fit$scaling[[key]])
    for(ver in c('a','b')) {
      nm <- paste0(key,ver);out[[nm]] <- plogis(drop(a[[ver]]%*%fit$models[[nm]]$coefficients))
    }
  };do.call(cbind,out)
}
paired_cluster_weights <- function(patient_id,hospital_id=NULL,B=1000,seed=2026100911) {
  set.seed(seed);pt<-match(patient_id,unique(patient_id));stopifnot(!anyNA(pt))
  hc<-if(is.null(hospital_id))NULL else match(hospital_id,unique(hospital_id))
  W<-matrix(0L,B,length(pt))
  for(b in seq_len(B)) {
    ww<-tabulate(sample.int(max(pt),max(pt),replace=TRUE),max(pt))[pt]
    if(!is.null(hc))ww<-ww*tabulate(sample.int(max(hc),max(hc),replace=TRUE),max(hc))[hc]
    W[b,]<-ww
  };W
}
evaluate_pairs <- function(y,predictions,patient_id,hospital_id=NULL,B=1000,seed=2026100911) {
  W<-paired_cluster_weights(patient_id,hospital_id,B,seed)
  stopifnot(nrow(predictions)==length(y),all(W%*%y>0),all(W%*%(1-y)>0),all(is.finite(predictions)))
  points<-boots<-list();metrics<-list()
  for(nm in colnames(predictions)) {
    pp<-predictions[,nm];points[[nm]]<-metric(y,pp,rep(1,length(y)),TRUE)
    bb<-t(vapply(seq_len(B),function(b)metric(y,pp,W[b,],TRUE),numeric(4)))
    colnames(bb)<-names(points[[nm]]);stopifnot(all(is.finite(bb)));boots[[nm]]<-bb
    for(me in names(points[[nm]])) {
      ci<-quantile(bb[,me],c(.025,.975))
      metrics[[length(metrics)+1L]]<-data.frame(model=nm,metric=me,estimate=points[[nm]][me],lower=ci[1],upper=ci[2])
    }
  }
  pairs<-c('T0','static_full','static_nonCV','dynamic_full','dynamic_nonCV');inc<-list();deltas<-list()
  for(key in pairs)for(me in c('AUC','Brier')) {
    a<-paste0(key,'a');b<-paste0(key,'b');sg<-if(me=='AUC')1 else -1
    dv<-sg*(boots[[b]][,me]-boots[[a]][,me]);ci<-quantile(dv,c(.025,.975))
    inc[[length(inc)+1L]]<-data.frame(pair=key,metric=paste0('Delta_',me),estimate=sg*(points[[b]][me]-points[[a]][me]),lower=ci[1],upper=ci[2])
    deltas[[length(deltas)+1L]]<-dv
  }
  inc<-do.call(rbind,inc);V<-do.call(cbind,deltas);se<-apply(V,2,sd);stopifnot(all(se>0))
  Z<-sweep(sweep(V,2,inc$estimate,'-'),2,se,'/')
  critical<-as.numeric(quantile(apply(abs(Z),1,max),.95))
  inc$simultaneous_lower<-inc$estimate-critical*se;inc$simultaneous_upper<-inc$estimate+critical*se
  inc$family_size<-10L;inc$critical<-critical
  list(metrics=do.call(rbind,metrics),increments=inc,replicates=boots)
}
