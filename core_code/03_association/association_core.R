# Source-derived core for nested mortality association and boot-then-impute.
# No study data, paths, fitted objects, or orchestration are included.
# Required packages: mice, sandwich. Outcome: hospital_mortality after landmark.

adjustment_spec <- function(local = FALSE) {
  spec <- list(M0=character(), M1=c('age_per10','female'),
    M2=c('age_per10','female','charlson_retrospective','infection_site'),
    M3=c('age_per10','female','charlson_retrospective','infection_site',
         'noncv_sofa','baseline_lactate','baseline_MAP','heart_rate10','baseline_logNEE'))
  if (!local) spec$M3 <- c(spec$M3,'icu_group')
  spec
}

analysis_input <- function(d, local = FALSE) {
  d <- as.data.frame(d)
  spec <- adjustment_spec(local)
  aux <- c('MAP_mean_ns1','MAP_mean_ns2','MAP_mean_ns3',
    'NEE_mean_ns1','NEE_mean_ns2','NEE_mean_ns3','MAP_slope10','MAP_sd10',
    'NEE_slope_log','NEE_SD_log')
  if (!local) aux <- c(aux,'ICU_fraction','MAP_coverage_ICU','NEE_coverage_ICU',
                      'NEE_positive_fraction','p_C1','p_C2')
  variables <- unique(c('hospital_mortality','Jclass',spec$M3,aux))
  stopifnot(all(variables %in% names(d)))
  x <- as.data.frame(d[,variables,drop=FALSE])
  x$Jclass <- factor(x$Jclass,levels=1:3)
  if (local) x$infection_site <- factor(x$infection_site,
    levels=c('undocumented','respiratory','urinary','multiple','other_documented'))
  else {
    x$infection_site <- factor(x$infection_site,levels=sort(unique(x$infection_site)))
    x$icu_group <- factor(x$icu_group,levels=sort(unique(x$icu_group)))
  }
  complete <- c('hospital_mortality','Jclass','infection_site',if (!local) 'icu_group')
  stopifnot(!anyNA(x[,complete,drop=FALSE]))
  x
}

fresh_imputation <- function(x, m=40L, seed=2026109000L, local=FALSE) {
  # Outcome and current assignments are association-imputation predictors;
  # neither the outcome nor a trajectory time series is imputed.
  stopifnot(m>=5L)
  missing <- colSums(is.na(x))
  method <- mice::make.method(x); method[] <- ''
  method[missing>0] <- 'pmm'
  active <- names(method)[method!='']
  stopifnot(all(vapply(x[active],is.numeric,logical(1))))
  predictor <- mice::make.predictorMatrix(x)
  predictor[method=='',] <- 0L; diag(predictor) <- 0L
  init <- NULL
  if (local) {
    init <- x
    for (v in active) init[[v]][is.na(init[[v]])] <- median(init[[v]],na.rm=TRUE)
  }
  imputation <- mice::mice(x,m=m,maxit=20,method=method,
    predictorMatrix=predictor,seed=seed,donors=5,printFlag=FALSE,data.init=init)
  history <- list()
  repeat {
    cv <- do.call(rbind,lapply(c('mean','sd'),function(par) {
      value <- as.data.frame(mice::convergence(imputation,parameter=par))
      value$parameter <- par; value
    }))
    final <- cv[cv$.it==imputation$iteration & cv$vrb %in% active,,drop=FALSE]
    final$missing_values <- as.integer(missing[final$vrb])
    final$applicable <- final$parameter=='mean' | final$missing_values>=2L
    evaluated <- final[final$applicable,,drop=FALSE]
    maximum <- if (nrow(evaluated)) max(evaluated$psrf) else 0
    complete_diagnostics <- nrow(final)==2L*length(active) &&
      !anyDuplicated(final[,c('vrb','parameter'),drop=FALSE]) && setequal(final$vrb,active)
    passed <- complete_diagnostics && all(is.finite(evaluated$psrf)) && is.finite(maximum) && maximum<=1.1
    history[[length(history)+1]] <- data.frame(iterations=imputation$iteration,
      max_psrf=maximum,applicable_statistics=nrow(evaluated),passed=passed)
    if (passed) break
    if (imputation$iteration>=120) stop('MI convergence gate failed: retain diagnostics and do not release estimates.')
    imputation <- mice::mice.mids(imputation,maxit=20,printFlag=FALSE)
  }
  list(mi=imputation,final_diagnostics=final,history=do.call(rbind,history),
    logged_events=imputation$loggedEvents,requested_method=method,
    requested_predictorMatrix=predictor)
}

cluster_covariance <- function(fit, patient_id, hospital_id=NULL) {
  cluster <- if (is.null(hospital_id)) patient_id else list(patient=patient_id,hospital=hospital_id)
  raw <- sandwich::vcovCL(fit,cluster=cluster,type='HC1',fix=FALSE)
  raw <- (raw+t(raw))/2
  ev <- eigen(raw,symmetric=TRUE,only.values=TRUE)$values
  negative <- any(ev<0)
  V <- if (negative) sandwich::vcovCL(fit,cluster=cluster,type='HC1',fix=TRUE) else raw
  list(V=V,audit=data.frame(min_eigen_raw=min(ev),negative_eigenvalues=sum(ev<0),
    PSD_projection_used=negative,
    relative_Frobenius_change=sqrt(sum((V-raw)^2))/max(sqrt(sum(raw^2)),.Machine$double.eps)))
}

rubin_scalar <- function(Q,U,dfcom) {
  m <- length(Q); q <- mean(Q); ub <- mean(U); bt <- var(Q)
  tt <- ub+(1+1/m)*bt; r <- (1+1/m)*bt/ub
  old <- if (r>0) (m-1)*(1+1/r)^2 else Inf
  lambda <- (1+1/m)*bt/tt
  observed <- (dfcom+1)/(dfcom+3)*dfcom*(1-lambda)
  df <- 1/(1/old+1/observed); se <- sqrt(tt)
  data.frame(logOR=q,SE=se,OR=exp(q),lower=exp(q-qt(.975,df)*se),
    upper=exp(q+qt(.975,df)*se),df=df,p=2*pt(-abs(q/se),df),
    within_variance=ub,between_MI_variance=bt,MI_MCSE=sqrt(bt/m))
}

fit_nested <- function(mi,patient_id,hospital_id=NULL,local=FALSE,outer_bootstrap=FALSE) {
  stopifnot(length(patient_id)==nrow(mi$data),!anyNA(patient_id))
  if (!is.null(hospital_id)) stopifnot(length(hospital_id)==length(patient_id),!anyNA(hospital_id))
  spec <- adjustment_spec(local)
  estimates <- diagnostics <- covariance_audit <- list()
  for (i in seq_len(mi$m)) {
    x <- mice::complete(mi,i)
    stopifnot(!anyNA(x),all(x$noncv_sofa>=0 & x$noncv_sofa<=20),
              all(x$noncv_sofa==round(x$noncv_sofa)))
    for (nm in names(spec)) {
      fit <- glm(reformulate(c('Jclass',spec[[nm]]),'hospital_mortality'),data=x,
        family=binomial(),control=glm.control(maxit=100,epsilon=1e-9))
      stopifnot(fit$converged,fit$rank==length(coef(fit)),all(is.finite(coef(fit))))
      if (outer_bootstrap) V <- vcov(fit) else {
        checked <- cluster_covariance(fit,patient_id,hospital_id)
        V <- checked$V
        covariance_audit[[length(covariance_audit)+1]] <- cbind(model=nm,imputation=i,checked$audit)
      }
      for (term in c('Jclass2','Jclass3')) estimates[[length(estimates)+1]] <- data.frame(
        imputation=i,model=nm,contrast=term,logOR=unname(coef(fit)[term]),within_U=V[term,term])
      diagnostics[[length(diagnostics)+1]] <- data.frame(model=nm,imputation=i,
        n=nobs(fit),events=sum(x$hospital_mortality),parameters=length(coef(fit)),
        rank=fit$rank,converged=fit$converged,iterations=fit$iter,
        min_fitted=min(fitted(fit)),max_fitted=max(fitted(fit)))
    }
  }
  dfcom <- length(unique(if (is.null(hospital_id)) patient_id else hospital_id))-1
  rows <- do.call(rbind,estimates)
  groups <- split(rows,interaction(rows$model,rows$contrast,drop=TRUE))
  pooled <- do.call(rbind,lapply(groups,function(a)
    cbind(model=a$model[1],contrast=a$contrast[1],rubin_scalar(a$logOR,a$within_U,dfcom))))
  list(individual_estimates=rows,reference_pooled=pooled,diagnostics=do.call(rbind,diagnostics),
       covariance_audit=do.call(rbind,covariance_audit))
}

boot_then_impute <- function(Q) {
  # Q[b,m] is the log-OR from fresh imputation m inside paired outer bootstrap b.
  # Model-based within-U is NOT added again; do not percentile-pool B*M values.
  stopifnot(is.matrix(Q),all(is.finite(Q)),nrow(Q)>1,ncol(Q)>1)
  B <- nrow(Q); M <- ncol(Q); means <- rowMeans(Q)
  MSW <- sum((Q-means)^2)/(B*(M-1))
  MSB <- M*sum((means-mean(Q))^2)/(B-1)
  between <- (MSB-MSW)/M
  stopifnot(between>0)
  variance <- (1+1/B)*between+MSW/(B*M)
  df <- max(3,variance^2/((((B+1)/(B*M))^2*MSB^2/(B-1)) + MSW^2/(B*M^2*(M-1))))
  q <- mean(Q); se <- sqrt(variance)
  data.frame(logOR=q,OR=exp(q),lower=exp(q-qt(.975,df)*se),upper=exp(q+qt(.975,df)*se),
    SE=se,df=df,p=2*pt(-abs(q/se),df),between_bootstrap_variance=between,
    within_MI_variance=MSW,outer_MCSE=sqrt((between+MSW/M)/B))
}
