# Audited numerical implementation of mitml 0.4-5 .D1, with named extraction
# of user-supplied cluster-robust complete-data covariance matrices.
# Not a call to mice::D1 and not its default model-based covariance estimate.
# Reference: https://amices.org/mice/reference/D1.html
# Installed reference source archived by probe_methods.R -> source_probe.txt.
robust_D1 <- function(q_list, u_list, block, dfcom = NULL) {
  m <- length(q_list); k <- length(block)
  stopifnot(m >= 2L, length(u_list) == m, k >= 1L,
            !anyDuplicated(block), all(vapply(q_list, is.numeric, logical(1))))
  Q <- matrix(NA_real_, k, m, dimnames=list(block,NULL))
  U <- array(NA_real_, c(k,k,m), dimnames=list(block,block,NULL))
  for (i in seq_len(m)) {
    q <- q_list[[i]]; v <- u_list[[i]]
    stopifnot(!is.null(names(q)), !anyDuplicated(names(q)),
              !is.null(rownames(v)), !is.null(colnames(v)),
              !anyDuplicated(rownames(v)), !anyDuplicated(colnames(v)),
              all(block %in% names(q)), all(block %in% rownames(v)),
              all(block %in% colnames(v)))
    Q[,i] <- q[block]
    V <- v[block,block,drop=FALSE]
    stopifnot(all(is.finite(V)), all(is.finite(Q[,i])),
              max(abs(V-t(V))) < 1e-7 * max(1,max(abs(V))))
    U[,,i] <- (V+t(V))/2
  }
  qbar <- rowMeans(Q)
  ubar <- apply(U,c(1,2),mean)
  b <- cov(t(Q))
  ev <- eigen(ubar,symmetric=TRUE,only.values=TRUE)$values
  # A pooled block covariance must be invertible. No silent nearPD/generalized inverse.
  if (min(ev) <= 0 || !is.finite(kappa(ubar)))
    stop('Pooled within-imputation block covariance is not positive definite.')
  r <- (1+1/m)*sum(diag(b %*% solve(ubar)))/k
  if (r < -1e-10) stop('Negative RIV exceeds numerical tolerance.')
  r <- max(0,r)
  F <- as.numeric(crossprod(qbar,solve((1+r)*ubar,qbar))/k)
  t <- k*(m-1)
  if (!is.null(dfcom)) {
    stopifnot(length(dfcom)==1L,is.finite(dfcom),dfcom>4,t>4)
    a <- r*t/(t-2)
    vs <- (dfcom+1)/(dfcom+3)*dfcom
    c0 <- 1/(t-4); c1 <- vs-2*(1+a); c2 <- vs-4*(1+a)
    if (c1 <= 0 || c2 <= 0)
      stop('Finite-df D1 approximation unstable: too few effective clusters/high missing information.')
    z <- 1/c2 + c0*(a^2*c1/((1+a)^2*c2)) +
      c0*(8*a^2*c1/((1+a)*c2^2) + 4*a^2/((1+a)*c2)) +
      c0*(4*a^2/(c2*c1) + 16*a^2*c1/c2^3) + c0*(8*a^2/c2^2)
    df2 <- 4+1/z
  } else if (r==0) {
    df2 <- Inf
  } else if (t>4) {
    df2 <- 4+(t-4)*(1+(1-2/t)/r)^2
  } else {
    df2 <- t*(1+1/k)*(1+1/r)^2/2
  }
  stopifnot(is.finite(F),F>=0,is.finite(df2)||is.infinite(df2),df2>0)
  data.frame(F=F,df1=k,df2=df2,p=pf(F,k,df2,lower.tail=FALSE),
             RIV=r,m=m,dfcom=if(is.null(dfcom))NA_real_ else dfcom,
             min_eigen_ubar=min(ev),condition_ubar=kappa(ubar),
             method='MI D1-type multivariate Wald using explicit cluster sandwich U')
}

# Use after fitting one common full GLM in each imputation, with fixed factor
# levels and spline knots, verifying full rank and convergence in every fit.
# q_list <- lapply(fits,coef)
# u_list <- lapply(fits,function(g) sandwich::vcovCL(g,cluster=...,type='HC1',fix=TRUE))
# MAP|NEE <- robust_D1(q_list,u_list,MAP_block_names,dfcom)
# NEE|MAP <- robust_D1(q_list,u_list,NEE_block_names,dfcom)
# The same B and observation-variable coefficients remain in both hypotheses.

audit_cluster_vcov <- function(g,cluster,type='HC1') {
  raw <- sandwich::vcovCL(g,cluster=cluster,type=type,fix=FALSE)
  raw <- (raw+t(raw))/2
  ev <- eigen(raw,symmetric=TRUE,only.values=TRUE)$values
  neg <- any(ev<0)
  V <- if(neg)sandwich::vcovCL(g,cluster=cluster,type=type,fix=TRUE)else raw
  rel <- sqrt(sum((V-raw)^2))/max(sqrt(sum(raw^2)),.Machine$double.eps)
  list(V=V,audit=data.frame(min_eigen_raw=min(ev),negative_eigenvalues=sum(ev<0),
       materially_negative=sum(ev < -1e-8*max(abs(ev))),PSD_projection_used=neg,
       relative_Frobenius_change=rel,type=type))
}
