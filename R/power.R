#------------------------------------------------------------------------------
# functions for power calculation
# Author: dlabes
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# helper function to allow partial match of the method

.powerMethod <- function(method){
  meth <- tolower(method[1])
  if (method=="") method <- "exact"
  methods <- c("exact","owenq","noncentral","nct","shifted","central","mvt")
  #                         ^ = match at start
  meth <- methods[grep(paste("^",meth,sep=""), methods)]
  if (length(meth)==0) meth <- tolower(method)
  if (length(meth)>1)  meth <- "nct" # happens if only "n" is given

  return(meth)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks,
# does not vectorize propperly!
# to be used in sampleN.TOST avoiding overhead of redundant calculations
# in case of multiplicative model:
# diffm=log(null ratio), theta1=log(lower BE limit), theta2=log(upper BE limit)
# in case of additive model:
# diffm=1-null ratio, theta1=lower BE limit-1, theta2=upper BE limit -1
# Jan 2015: Interface changed to sem
# so call it with sem= se*sqrt(bk/n) if balanced or se*sqrt(bkni*sum(1/n))
.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df)
{

  OwensQ <- NULL
  # Allow alpha to be scalar or a matrix. If matrix, the interpretation is
  # - Require 2 columns
  # - Column 1 = alpha value for left hypothesis
  # - Column 2 = alpha value for right hypothesis
  # - Convention: If multiple rows, require diffm & sem to be of the same length
  #   and evaluate power element-wise for each combination (alpha, diffm, sem)
  if (is.atomic(alpha) && !is.matrix(alpha)) {
    # We enforce matrix structure -> recycling will not work, so do it manually
    if (length(alpha) != max(length(diffm), length(sem))) {
      alpha <- rep.int(alpha, max(length(diffm), length(sem)))
    }
    alpha <- matrix(alpha, ncol = 1)
  } else {
    if (nrow(alpha) != length(diffm) || length(diffm) != length(sem))
      stop("Number of rows of alpha must match length of diffm and sem.")
  }
  dl <- ncol(alpha)
  if (dl > 2)
    stop("Number of columns of alpha should be 1 or 2.")
  if (length(df) == 1) {
    df <- rep.int(df, max(length(diffm), length(sem)))
  } else {
    if (length(df) != max(length(diffm), length(sem)))
      stop("df must be of the same length as diffm and sem.")
  }
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)
  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2 and sem=0!
  delta1 <- (diffm-ltheta1)/sem
  delta2 <- (diffm-ltheta2)/sem
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0
  # R is infinite in case of alpha=0.5 where tval is == 0
  R <- (delta1-delta2)*sqrt(df)/(tval[,1]+tval[,dl])
  # in case of se=0 it results: delta1=Inf, delta2=inf if diffm>ltheta2
  # Inf - Inf is NaN
  R[is.nan(R)] <- 0

  # if alpha>0.5 (very unusual!) t(1-alpha,df) is <0 and then R is negative
  # i.e in OwensQ the upper integration limit is lower then the lower limit!
  # SAS OwenQ gives missings if b or a are negative!
  # On the other hand SAS Proc Power gives values which are seemingly calculated
  # with abs(R).
  # Correct acc. to Fig. 1 given in K.Philips
  # "Power for Testing Multiple Instances of the Two One-Sided Tests Procedure"
  # The International Journal of Biostatistics: Vol. 5: Iss. 1, Article 15.
  # DOI: 10.2202/1557-4679.1169
  # should be R=Inf, i.e. unlimited integration with respect to sigma.
  # This gives the same values (within certain precision) as Ben's power.1TOST
  # aka power.TOST(..., method="mvt").
  # Can also be checked via function power.TOST.sim().
  R[R<=0] <- Inf
  # to check SAS Proc power values comment above out and write
  # R <- abs(R)

  # to avoid numerical errors in OwensQ implementation
  if (min(df)>10000){
    # 'shifted' normal approximation Jan 2015
    # former Julious formula (57)/(58) doesn't work
    tval <- qnorm(1-alpha)
    p1   <- pnorm(tval[,1]-delta1)
    p2   <- pnorm(-tval[,dl]-delta2)
    # may give negative values
    # thus set to zero
    pwr <- p2-p1
    pwr[pwr<0] <- 0
    return(pwr)
  }
  if (min(df)>=5000 & min(df<=10000)) {
    # approximation via non-central t-distribution
    return(.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df))
  }

  # attempt to vectorize (it vectorizes properly if diffm is a vector
  # OR se OR n,df are vectors)
  p1 <- vector("numeric", length(delta1))
  p2 <- vector("numeric", length(delta1))
  for (i in seq_along(delta1)) {
    # get correct values for p1 and p2:
    # p1 is for left hypothesis (right-tailed)
    # -> critical value is always first column
    # p2 for right hypothesis (left-tailed)
    # -> critical value is first column if alpha has 1 column,
    #    second column if alpha has 2 columns
    p1[i] <- OwensQ(df[i], tval[i, 1], delta1[i], 0, R[i])
    p2[i] <- OwensQ(df[i], -tval[i, dl], delta2[i], 0, R[i])
  }
  pwr <- p2-p1
  # due to numeric inaccuracies power < 0
  # paranoia
  pwr[pwr<0] <- 0
  return( pwr )
}

#------------------------------------------------------------------------------
# Ben's implementation of power via integration of the bivariate t-distribution
# with correlation == 1, also exact
# does'nt vectorize in any respect!
.power.1TOST <- function(alpha, ltheta1, ltheta2, diffm, sem, df, setseed = TRUE)
{
  pmvt <- GenzBretz <- OwensQ <- NULL
  if (is.atomic(alpha) && !is.matrix(alpha)) {
    # We enforce matrix structure -> recycling will not work, so do it manually
    if (length(alpha) != max(length(diffm), length(sem))) {
      alpha <- rep.int(alpha, max(length(diffm), length(sem)))
    }
    alpha <- matrix(alpha, ncol = 1)
  } else {
    if (nrow(alpha) != length(diffm) || length(diffm) != length(sem))
      stop("Number of rows of alpha must match length of diffm and sem.")
  }
  dl <- ncol(alpha)
  if (dl > 2)
    stop("Number of columns of alpha should be 1 or 2.")
  if (length(df) == 1) {
    df <- rep.int(df, max(length(diffm), length(sem)))
  } else {
    if (length(df) != max(length(diffm), length(sem)))
      stop("df must be of the same length as diffm and sem.")
  }
  if (setseed) set.seed(123456)

  corr  <- matrix(1, ncol = 2, nrow = 2)
  tval  <- qt(1 - alpha, df)
  delta1 <- (diffm - ltheta1) / sem
  delta2 <- (diffm - ltheta2) / sem
  pow <- rep(0, times=length(delta1))
  # attempt to vectorize if ltheta0 OR se is a vector
  for(i in seq_along(delta1)){
    lower <- c(tval[i, 1], -Inf)
    upper <- c(Inf, -tval[i, dl])
    delta <- c(delta1[i], delta2[i])
    prob  <- pmvt(lower = lower, upper = upper, delta = delta, df = df[i],
                  corr = corr,
                  algorithm = GenzBretz(maxpts=100000, abseps = 1e-05))#[1]
    # abseps=1e-6 gives often "Completion with error > abseps"
    # give a warning if attr(prob,"msg") not equal "Normal completion"?
    if(attr(prob, which="msg")!="Normal Completion")
      warning("pmvt returned message ", attr(prob, which="msg"), call.=FALSE)
    pow[i] <- prob[1]
  }
  pow
}

#------------------------------------------------------------------------------
# 'raw' approximate power function without any error checks,
# approximation based on non-central t
# this vectorizes ok
.approx.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df)
{
  if (is.atomic(alpha) && !is.matrix(alpha)) {
    # We enforce matrix structure -> recycling will not work, so do it manually
    if (length(alpha) != max(length(diffm), length(sem))) {
      alpha <- rep.int(alpha, max(length(diffm), length(sem)))
    }
    alpha <- matrix(alpha, ncol = 1)
  } else {
    if (nrow(alpha) != length(diffm) || length(diffm) != length(sem))
      stop("Number of rows of alpha must match length of diffm and sem.")
  }
  dl <- ncol(alpha)
  if (dl > 2)
    stop("Number of columns of alpha should be 1 or 2.")
  tval <- qt(1 - alpha, df, lower.tail = TRUE, log.p = FALSE)

  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
  # and sem=0!
  delta1 <- (diffm-ltheta1)/sem
  delta2 <- (diffm-ltheta2)/sem
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0

  # suppress warnings with regard to insufficient precision of nct
  pow <- suppressWarnings(pt(-tval[, dl], df, ncp=delta2) -
                            pt(tval[, 1], df, ncp=delta1))
  pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)

  return(pow)
}
#------------------------------------------------------------------------------
# 'raw' power function without any error checks,
# approximation based on central 'shifted' central t distribution
# according to Chow, Liu "Design and Analysis of Bioavailability ..."
# Chapter 9.6 and implemented in PASS 2008
# where does this all come from?
.approx2.power.TOST <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df)
{
  if (is.atomic(alpha) && !is.matrix(alpha)) {
    # We enforce matrix structure -> recycling will not work, so do it manually
    if (length(alpha) != max(length(diffm), length(sem))) {
      alpha <- rep.int(alpha, max(length(diffm), length(sem)))
    }
    alpha <- matrix(alpha, ncol = 1)
  } else {
    if (nrow(alpha) != length(diffm) || length(diffm) != length(sem))
      stop("number of rows of alpha must match length of diffm and sem.")
  }
  dl <- ncol(alpha)
  if (dl > 2)
    stop("Number of columns of alpha should be 1 or 2.")
  tval   <- qt(1 - alpha, df, lower.tail = TRUE)

  # 0/0 -> NaN in case diffm=ltheta1 or diffm=ltheta2
  # and se=0!
  delta1 <- (diffm-ltheta1)/sem
  delta2 <- (diffm-ltheta2)/sem
  # is this correct?
  delta1[is.nan(delta1)] <- 0
  delta2[is.nan(delta2)] <- 0

  pow <- pt(-tval[, dl]-delta2, df) - pt(tval[, 1]-delta1, df)
  pow[pow<0] <- 0 # this is to avoid neg. power due to approx. (vector form)

	return(pow)
}

#------------------------------------------------------------------------------
# function for merging the various power calculations
.calc.power <- function(alpha=0.05, ltheta1, ltheta2, diffm, sem, df, method="exact")
{
  pow <- switch(
      method,
      exact=.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      owenq=.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      mvt=  .power.1TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      nct=  .approx.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      noncentral=.approx.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      shifted=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      central=.approx2.power.TOST(alpha, ltheta1, ltheta2, diffm, sem, df),
      stop("Method '", method, "' unknown!\n", call.=TRUE)
  )
  return(pow)
}
