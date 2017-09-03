############################
# Various helper functions #
##########################################################
# Partition n into grps with 'best' balance between grps #
##########################################################
nvec <- function(n, grps) {
  ni <- trunc(n/grps)
  nv <- rep.int(ni, times=grps)
  rest <- n-grps*ni
  if(rest!=0){
    nv <- nv + c(rep.int(1,rest), rep.int(0,grps-rest))
  }
  return(nv)
}

####################################################################
# Calculates the CV from the CI in the final analysis taking the   #
# loss of one df for the factor 'stage' in the model into account. #
# Simplified from PowerTOST's CVfromCI()                           #
####################################################################
CVfromFinalCI <- function(lower, upper, alpha=0.0294, n, design="2x2x2") {
  pe <- sqrt(lower*upper)
  if (length(n) == 1) {
    # n given as ntotal
    n <- nvec(n = n, grps = 2)
    if (n[1] != n[length(n)]) {
      message("Unbalanced ", design, " design. n(i)= ",
              paste(n, collapse = "/"), " assumed.")
    }
  } else {
    # n given as vector of No. of subjects in (sequence) groups
    if (length(n) != 2) stop("Length of n vector must be 2!")
  }
  nc     <- sum(1/n)
  n      <- sum(n)
  if (design == "parallel") bkni <- 1 else bkni <- 0.5
  se.fac <- sqrt(bkni * nc)
  df     <- n-3 # One df less than usual ('stage' in the model).
  tval   <- qt(1 - alpha, df)
  s1     <- (log(pe) - log(lower))/se.fac/tval
  s2     <- (log(upper) - log(pe))/se.fac/tval
  sw     <- 0.5 * (s1 + s2)
  if (abs(s1 - s2)/sw > 0.1) {
    warning(paste("sigma based on pe & lower CL more than 10% different than\n",
                  "sigma based on pe & upper CL. Check input."))
  }
  return(se2CV(sw))
}

#########################################################
# Function to estimate the TIE for the given conditions #
# in stage 1.                                           #
#########################################################
TIE1.est <- function(meth, alpha0, alpha1, alpha2, n1, GMR, CV1, target,
                     pmethod, usePE, int.pwr, min.n2, max.n, Nmax, fCrit,
                     fClow, theta2, nsims, setseed, KM, asym, Xover) {
  if (asym) alpha12 <- c(alpha1, alpha2) else alpha12 <- rep(alpha2, 2)
  if (!KM) { # Common methods.
    if (Xover) { # Crossover.
      power.2stage.fC(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                      GMR=GMR, CV=CV1, targetpower=target,
                      pmethod=pmethod, usePE=usePE, powerstep=int.pwr,
                      min.n2=min.n2, max.n=max.n, fCrit=fCrit,
                      fClower=fClow, theta0=theta2, nsims=nsims,
                      setseed=setseed)$pBE_s1
    } else {     # Parallel.
      power.2stage.p(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                     GMR=GMR, CV=CV1, targetpower=target,
                     pmethod=pmethod, usePE=usePE, Nmax=Nmax,
                     test="welch", theta0=theta2, nsims=nsims,
                     setseed=setseed)$pBE_s1
    }
  } else {       # Adaptive (Karalis/Macheras and Karalis).
    power.2stage.KM(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                    CV=CV1, targetpower=target, pmethod=pmethod,
                    Nmax=Nmax, theta0=theta2, nsims=nsims,
                    setseed=setseed)$pBE_s1
  }
}

#########################################################
# Function to estimate the TIE for the given conditions #
# in the final analysis.                                #
#########################################################
TIE.est <- function(meth, alpha0, alpha1, alpha2, n1, GMR, CV1, target,
                    pmethod, usePE, int.pwr, min.n2, max.n, Nmax, fCrit,
                    fClow, theta2, nsims, setseed, KM, asym, Xover) {
  if (asym) alpha12 <- c(alpha1, alpha2) else alpha12 <- rep(alpha2, 2)
  if (!KM) { # Common methods.
    if (Xover) { # Crossover.
      power.2stage.fC(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                      GMR=GMR, CV=CV1, targetpower=target,
                      pmethod=pmethod, usePE=usePE, powerstep=int.pwr,
                      min.n2=min.n2, max.n=max.n, fCrit=fCrit,
                      fClower=fClow, theta0=theta2, nsims=nsims,
                      setseed=setseed)$pBE
    } else {     # Parallel.
      power.2stage.p(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                     GMR=GMR, CV=CV1, targetpower=target,
                     pmethod=pmethod, usePE=usePE, Nmax=Nmax,
                     test="welch", theta0=theta2, nsims=nsims,
                     setseed=setseed)$pBE
    }
  } else {       # Adaptive (Karalis/Macheras and Karalis).
    power.2stage.KM(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                    CV=CV1, targetpower=target, pmethod=pmethod,
                    Nmax=Nmax, theta0=theta2, nsims=nsims,
                    setseed=setseed)$pBE
  }
}

#######################################################
# Function to estimate power for the given conditions #
# in the final analysis.                              #
#######################################################
pwr.est <- function(meth, alpha0, alpha1, alpha2, n1, GMR, CV1, target,
                    pmethod, usePE, int.pwr, min.n2, max.n, Nmax, fCrit,
                    fClow, nsims, setseed, KM, asym, Xover) {
  if (asym) alpha12 <- c(alpha1, alpha2) else alpha12 <- rep(alpha2, 2)
  if (!KM) { # Common methods.
    if (Xover) { # Crossover.
      power.2stage.fC(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                      GMR=GMR, CV=CV1, targetpower=target, pmethod=pmethod,
                      usePE=usePE, powerstep=int.pwr, min.n2=min.n2,
                      max.n=max.n, fCrit=fCrit, fClower=fClow, theta0=GMR,
                      nsims=1e5, setseed=setseed)
    } else {     # Parallel.
      power.2stage.p(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                     GMR=GMR, CV=CV1, targetpower=target, pmethod=pmethod,
                     usePE=usePE, Nmax=Nmax, test="welch", theta0=GMR,
                     nsims=1e5, setseed=setseed)
    }
  } else {       # Adaptive (Karalis/Macheras and Karalis).
    power.2stage.KM(method=meth, alpha0=alpha0, alpha=alpha12, n1=n1,
                    CV=CV1, targetpower=target, pmethod=pmethod,
                    Nmax=Nmax, theta0=GMR, nsims=1e5,
                    setseed=setseed)
  }
}

##########################################
# Function for BE in the final analysis. #
##########################################
BEpooled <- function(alpha2, PE, CV, N, design) {
  N  <- nvec(N, 2) # vectorize
  df <- N[1]+N[2]-3  # one df less than usual (stage in the model!)
  nc <- sum(1/N)
  N  <- sum(N)
  if (design == "2x2x2") {
    bk   <- 2   # 2x2x2 crossover design const
    bkni <- 0.5 # design constant in terms of n(i)
  }
  if (design == "parallel") {
    bk   <- 4   # 2-group parallel design constant
    bkni <- 1   # design constant in terms of n(i)
  }
  mse       <- CV2mse(CV)
  diff      <- log(PE)
  CI        <- exp(diff + c(-1, +1)*qt(1-alpha2, df=df)*sqrt(mse*bkni*nc))
  names(CI) <- c("lower", "upper")
  return(CI)
}
