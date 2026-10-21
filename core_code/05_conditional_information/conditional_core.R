# Source robust_MI_core.R first. Requires splines and sandwich.
# Completed datasets must come from the shared current-cohort association MI.

conditional_basis <- function(training_features) {
  out<-lapply(c('MAP_mean','logNEE_mean'),function(nm) {
    v<-training_features[[nm]];stopifnot(all(is.finite(v)))
    list(knots=as.numeric(quantile(v,c(1/3,2/3))),bounds=range(v))
  });names(out)<-c('MAP','NEE');out
}
component_design <- function(features,basis_spec) {
  basis<-function(v,key){s<-basis_spec[[key]];splines::ns(v,knots=s$knots,Boundary.knots=s$bounds)}
  a<-basis(features$MAP_mean,'MAP');n<-basis(features$logNEE_mean,'NEE')
  colnames(a)<-paste0('MAP_mean_ns',1:3);colnames(n)<-paste0('NEE_mean_ns',1:3)
  d<-data.frame(a,MAP_slope10=features$MAP_slope24/10,MAP_sd10=features$MAP_SD/10,
                n,NEE_slope_log=features$logNEE_slope24,NEE_SD_log=features$logNEE_SD)
  stopifnot(all(is.finite(as.matrix(d))));d
}
conditional_analysis <- function(completed,features,basis_spec,patient_id,hospital_id=NULL,
                                 local=FALSE,positive_fraction_sensitivity=FALSE) {
  stopifnot(length(completed)==40L,length(patient_id)==nrow(features))
  MB<-c(paste0('MAP_mean_ns',1:3),'MAP_slope10','MAP_sd10')
  NB<-c(paste0('NEE_mean_ns',1:3),'NEE_slope_log','NEE_SD_log')
  B<-c('age_per10','female','noncv_sofa','baseline_lactate','baseline_MAP','heart_rate10',
       'baseline_logNEE','charlson_retrospective','infection_site')
  O<-character()
  if(!local){B<-c(B,'icu_group');O<-c('ICU_fraction','MAP_coverage_ICU','NEE_coverage_ICU')}
  h<-component_design(features,basis_spec)
  if(positive_fraction_sensitivity) {
    stopifnot(!local,all(is.finite(features$NEE_positive_fraction)))
    NB<-c(NB,'NEE_positive_fraction');h$NEE_positive_fraction<-features$NEE_positive_fraction
  }
  Q<-U<-audit<-list();dfcom<-length(unique(if(is.null(hospital_id))patient_id else hospital_id))-1L
  cluster<-if(is.null(hospital_id))patient_id else list(patient=patient_id,hospital=hospital_id)
  for(i in seq_along(completed)) {
    xx<-as.data.frame(completed[[i]]);stopifnot(nrow(xx)==nrow(h))
    for(nm in names(h))xx[[nm]]<-h[[nm]]
    g<-glm(reformulate(c(B,O,NB,MB),'hospital_mortality'),data=xx,family=binomial(),
           control=glm.control(maxit=100,epsilon=1e-9))
    stopifnot(g$converged,g$rank==length(coef(g)),nobs(g)==length(patient_id),all(is.finite(coef(g))))
    v<-audit_cluster_vcov(g,cluster);Q[[i]]<-coef(g);U[[i]]<-v$V;audit[[i]]<-v$audit
  }
  tests<-rbind(cbind(direction='MAP_given_NEE',robust_D1(Q,U,MB,dfcom)),
               cbind(direction='NEE_given_MAP',robust_D1(Q,U,NB,dfcom)))
  tests$p_Holm<-p.adjust(tests$p,'holm')
  list(tests=tests,coefficients=Q,covariance=U,covariance_audit=do.call(rbind,audit),
       dfcom=dfcom,MAP_block=MB,NEE_block=NB,basis=basis_spec)
}
conditional_contrast <- function(result,contrast) {
  # contrast is a named vector, with one element for each model coefficient.
  Q<-result$coefficients;U<-result$covariance;dfcom<-result$dfcom
  stopifnot(identical(names(contrast),names(Q[[1]])))
  qq<-vapply(Q,function(q)sum(q*contrast),numeric(1))
  uu<-vapply(U,function(v)drop(t(contrast)%*%v%*%contrast),numeric(1))
  qb<-mean(qq);ub<-mean(uu);bt<-var(qq);tt<-ub+(1+1/length(Q))*bt
  r<-(1+1/length(Q))*bt/ub;old_df<-if(r>0)(length(Q)-1)*(1+1/r)^2 else Inf
  lambda<-(1+1/length(Q))*bt/tt;obs_df<-(dfcom+1)/(dfcom+3)*dfcom*(1-lambda)
  df<-1/(1/old_df+1/obs_df);ci<-qb+c(-1,1)*qt(.975,df)*sqrt(tt)
  data.frame(OR=exp(qb),lower=exp(ci[1]),upper=exp(ci[2]),df=df,p=2*pt(-abs(qb/sqrt(tt)),df))
}
reported_component_contrasts <- function(result) {
  basis<-function(v,key){s<-result$basis[[key]];splines::ns(v,knots=s$knots,Boundary.knots=s$bounds)}
  cv<-function(name,amount=1){x<-setNames(rep(0,length(result$coefficients[[1]])),names(result$coefficients[[1]]));x[name]<-amount;x}
  contrasts<-list(
    MAP_mean_80_vs70={x<-cv(result$MAP_block[1],0);x[result$MAP_block[1:3]]<-drop(basis(80,'MAP')-basis(70,'MAP'));x},
    MAP_slope_plus10=cv('MAP_slope10'),MAP_SD_plus5=cv('MAP_sd10',.5),
    NEE_log_location_010_vs005={x<-cv(result$NEE_block[1],0);x[result$NEE_block[1:3]]<-drop(basis(log(.10+.01),'NEE')-basis(log(.05+.01),'NEE'));x},
    NEE_slope_plus_doubling_per_day=cv('NEE_slope_log',log(2)),
    NEE_logSD_plus05=cv('NEE_SD_log',.5))
  do.call(rbind,lapply(names(contrasts),function(nm)cbind(contrast=nm,conditional_contrast(result,contrasts[[nm]]))))
}
