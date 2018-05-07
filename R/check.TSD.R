check.TSD <- function(Var1, PE1, n1, Var, PE, N, type = 1, usePE = FALSE,
                      GMR, alpha0 = 0.05, alpha1 = 0.0294, alpha2 = 0.0294,
                      theta1, theta2, target = 0.80,
                      pmethod = c("shifted", "nct", "exact"),
                      int.pwr = TRUE, min.n2 = 0, max.n = Inf, Nmax = Inf,
                      fCrit = c("PE", "CI"), fClow = 0, nsims = 1e6,
                      setseed = TRUE, tol = 1e-8, pa = FALSE, skip = TRUE,
                      algo = 1, plot.it = FALSE, valid = FALSE, expl = 3,
                      Xover = TRUE, stop1 = FALSE, KM = FALSE,
                      KM.des = c("TSD", "TSD-1", "TSD-2"), CIs = FALSE)
{
  ######################################################################
  ## Check.TSD.R                                                      ##
  ## 1. Assess the empiric Type I Error (TIE) and power of crossover  ##
  ##    or parallel Two-Stage (Sequential) Designs based on the con-  ##
  ##    ditions given in the study protocol and stage 1 data by simu- ##
  ##    lations. If the TIE is not inflated, the applied method and   ##
  ##    alpha(s) should be considered justified and the result of the ##
  ##    study acceptable.                                             ##
  ## 2. Optimize the alpha to obtain a TIE as close as possible to    ##
  ##    the fixed consumer risk of alpha0 (generally 0.05) by the R-  ##
  ##    function uniroot().                                           ##
  ##    Notes:                                                        ##
  ##    a. It is futile to expect a TIE of /exactly/ 0.05 in simula-  ##
  ##       tions. Generally in one million simulations a nonsignifi-  ##
  ##       cant inflation of the TIE (one-sided limit of the binomial ##
  ##       test for one mio sim's ~0.05036) is considered sufficient. ##
  ##    b. Following an internal (!) statement of the BSWP optimiza-  ##
  ##       tion will be performed for /any/ TIE >0.05.                ##
  ##    c. In an asymmetric split (alpha1 != alpha2), alpha1 is kept  ##
  ##       and alpha2 optimized.                                      ##
  ## ---------------------------------------------------------------- ##
  ## Author:                                                          ##
  ##   Helmut Schütz                                                  ##
  ##   BEBAC, Neubaugasse 36/11, 1070 Vienna, Austria                 ##
  ##   helmut.schuetz@bebac.at                                        ##
  ## ---------------------------------------------------------------- ##
  ## License:                                                         ##
  ##   GPL-3 (http://www.gnu.org/licenses/gpl-3.0.en.html)            ##
  ## ---------------------------------------------------------------- ##
  ## Required:                                                        ##
  ##   Power2Stage, PowerTOST                                         ##
  ## ---------------------------------------------------------------- ##
  ## Tested on Windows 7 Pro SP1 64bit:                               ##
  ##   R 3.2.5 64bit (2016-04-14)                                     ##
  ##   R 3.3.0 64bit (2016-05-03)                                     ##
  ##   R 3.3.1 64bit (2016-06-21)                                     ##
  ##   R 3.4.1 64bit (2017-06-30)                                     ##
  ##   R 3.4.4 64bit (2018-03-15)                                     ##
  ##   Power2Stage 0.4-3 (2015-11-24)                                 ##
  ##   Power2Stage 0.4-5 (2017-08-24) beta-release                    ##
  ##   Power2Stage 0.5-1 (2018-04-03)                                 ##
  ##   PowerTOST 1.3-5 (2016-04-12)                                   ##
  ##   PowerTOST 1.3-6 (2016-06-06)                                   ##
  ##   PowerTOST 1.3-7 (2016-06-07) beta release                      ##
  ##   PowerTOST 1.4-6 (2017-08-19) beta-release                      ##
  ##   PowerTOST 1.4-7 (2018-04-12)                                   ##
  ## ---------------------------------------------------------------- ##
  ## History:                                                         ##
  ## 2016-05-02 v0.0: New                                             ##
  ## 2016-05-03 v0.1: - Merged the two optimizations functions.       ##
  ##                  - Added power of specified alphas.              ##
  ##                  - Enhanced the sample size plots.               ##
  ##                  - Windows-conform line-endings (CR/LF)          ##
  ##                  - Documents: System, user, OS, versions of R,   ##
  ##                    Power2Stage, PowerTOST, and this code.        ##
  ##                  - Calculates CI in the final analysis based on  ##
  ##                    specified and - if applicable - adjusted      ##
  ##                    alpha.                                        ##
  ## 2016-05-04 v0.2: - The power analysis is now optional.           ##
  ##                    Follows the boolean argument pa (TRUE|FALSE). ##
  ##                  - Added the interim analysis.                   ##
  ##                  - Tries to reconstruct the applied method from  ##
  ##                    the specified type and alpha(s).              ##
  ##                  - The alpha-optimzation is now optional. Skips  ##
  ##                    this time-consuming step if the TIE with the  ##
  ##                    specified alpha(s) is <=alpha0 and the argu-  ##
  ##                    ment skip is set to TRUE.                     ##
  ## 2016-05-11 v0.3: - Added parallel TSDs (by Welch-test only).     ##
  ##                  - Variabilities can be given as CV or MSE.      ##
  ##                  - PE can be given as a ratio (backtransformed)  ##
  ##                    or as a difference of log-means.              ##
  ##                  - In the interim all specified conditions (N,   ##
  ##                    futility conditions(s), alpha0 or alpha1 in   ##
  ##                    Type 1 (based on power) are assessed.         ##
  ##                  - Corrected the final analysis (one degree of   ##
  ##                    less due to stage as a fixed effect in the    ##
  ##                    model).                                       ##
  ##                  - If adjustment was done, gives a statement     ##
  ##                    about agreement and regulatory suggestions.   ##
  ##                  - Included validation examples.                 ##
  ## 2016-05-16 v0.4: - Added new algorithm: Simulate a set of TIEs   ##
  ##                    around the optimzed alpha. Fit two models     ##
  ##                    (linear and quadratic), select the better one ##
  ##                    by minimum AIC and solve for TIE = alpha0.    ##
  ##                    algo <- 1: Like in previous versions.         ##
  ##                    algo <- 2: Fitted model.                      ##
  ##                  - Started to work on vectorized CV and n1       ##
  ##                    (needed for Welsh-test). Buggy!               ##
  ## 2016-06-14 v0.5: - Changed to a function and added a fair amount ##
  ##                    of input / consistency checking.              ##
  ##                  - Convergence tolerance for uniroot() can be    ##
  ##                    given by the argument tol. Defaults to 1e-8.  ##
  ##                  - Added a plot of PE and CI (for specified and  ##
  ##                    adjusted alpha).                              ##
  ##                  - In case the CVs are not given in the report,  ##
  ##                    a four-element vector can be given in Var1    ##
  ##                    and/or Var. 1st two elements = lower and      ##
  ##                    upper CI, 3rd element = alpha, 4th element =  ##
  ##                    "CI".                                         ##
  ##                    Warning: The CV based on the final CI (from   ##
  ##                    Var) is only approximate since the factor     ##
  ##                    'stage' is not taken into account. Hence, the ##
  ##                    calculated CI will disagree with the report!  ##
  ##                    However, this has no influence on the esti-   ##
  ##                    mated TIE (based on the /interim/ data).      ##
  ## 2016-06-17 v0.6: - Additional information supporting validation. ##
  ##                  - Estimate N in sampleN.TOST(method=pmethod)    ##
  ##                    for consistency.                              ##
  ##                  - Bug corrected in assessing CI-futility.       ##
  ##                  - CV is now correctly calculated from the CI in ##
  ##                    in the final analysis (one df less due to     ##
  ##                    factor 'stage' in the model).                 ##
  ##                  - First package-version.                        ##
  ##                  - Replaced windows() by dev.new() for cross-    ##
  ##                    platform compatibility.                       ##
  ## 2016-06-21 v0.7: - If the study's specifications match any of    ##
  ##                    the published methods check whether n1, CV1,  ##
  ##                    max.n, GMR, target are within the validated   ##
  ##                    range. If not issue a warning.                ##
  ## 2016-09-19 v0.8: - Changed the specification of PEs:             ##
  ##                    "log" => "difflog", "lin" = "ratio".          ##
  ## 2016-09-27 v0.8.1: - Added the disclaimer to the output.         ##
  ## 2016-10-14 v0.8.2: - Added (optional) 95% CI of the TIE(s).      ##
  ## 2017-08-29 v0.8.3: - Changed the lower interval in uniroot()     ##
  ##                      from 0 to tol. Does not need any more to    ##
  ##                      force pmethod to "nct". However,            ##
  ##                      pmethod="exact" is  /extremely/ slow!       ##
  ##                    - Published on GitHub.                        ##
  ##                    - n2 with sampleN2.TOST of Power2Stage.       ##
  ## 2018-04-23 v0.8.4: - Adapted syntax to Power2Stage v0.5-1        ##
  ######################################################################
  ## PROGRAM OFFERED FOR USE WITHOUT ANY GUARANTEES AND ABSOLUTELY NO ##
  ## WARRANTY. NO LIABILITY IS ACCEPTED FOR ANY LOSS AND RISK TO      ##
  ## PUBLIC HEALTH RESULTING FROM USE OF THIS R-CODE.                 ##
  ######################################################################
  ## Arguments:                                                       ##
  ## ================================================================ ##
  ## Stage 1 (interim) data:                                          ##
  ## ---------------------------------------------------------------- ##
  ##   Var1:    Vector of observed variability.                       ##
  ##            - If two elements:                                    ##
  ##              1: value, 2: "CV" or "MSE".                         ##
  ##              Use the one with the highest numeric precision      ##
  ##              given in the report.                                ##
  ##            - If four elements:                                   ##
  ##              1: lower CL, 2: upper CL, 3: alpha, 4: "CI".        ##
  ##   PE1:     Observed PE (two element vector).                     ##
  ##            1: value, 2: "ratio" or "difflog".                    ##
  ##            "ratio": back-transformed (T/R-ratio /not/ in %).     ##
  ##            "difflog": difference of means in log-scale.          ##
  ##            Use the one with the highest numeric precision given  ##
  ##            in the report.                                        ##
  ##            Not needed if Var1 contains four elements.            ##
  ##  n1:       Sample size.                                          ##
  ##            Note: In the case of parallel designs with unequal    ##
  ##            group sizes, n1 can be given as a two-element vector, ##
  ##            where the 1st element is the group under T and the    ##
  ##            2nd one under R.                                      ##
  ##  stop1:    TRUE|FALSE. Set to TRUE if the study stopped in the   ##
  ##            1st stage (BE demonstrated). The TIE will be assessed ##
  ##            for the 1st stage only and all other calculations     ##
  ##            skipped.                                              ##
  ##            Defaults to FALSE (study proceeded to the 2nd stage). ##
  ## ---------------------------------------------------------------- ##
  ## Final (pooled) data: Not needed if stop1 == TRUE.                ##
  ## ---------------------------------------------------------------- ##
  ##   Var:     Observed variability (see also above).                ##
  ##   PE:      Observed PE (see also above).                         ##
  ##   N:       Total sample size.                                    ##
  ##            Currently unequal group sizes of parallel designs are ##
  ##            not implemented.                                      ##
  ## ---------------------------------------------------------------- ##
  ## Conditions (hopefully) stated in the protocol:                   ##
  ## ---------------------------------------------------------------- ##
  ##   type:    'Type' of design (for the definition see Schütz 2015) ##
  ##            1: Potvin et al. "B", Fuglsang 2013 "B",              ##
  ##               Fuglsang 2014 "B" (parallel), Xu et al. "E"        ##
  ##            2: Potvin et al. "C", Montague et al. "D",            ##
  ##               Fuglsang 2013 "C/D", Fuglsang 2014 "B" (parallel), ##
  ##               Xu et al. "F"                                      ##
  ##   usePE:   if FALSE (default) a fixed GMR (specified below) is   ##
  ##            used.                                                 ##
  ##            if TRUE the PE observed in stage 1 is used.           ##
  ##   GMR:     Fixed GMR for sample size estimation; will /not/ be   ##
  ##            used if above usePE = TRUE or in TSDs of Karalis /    ##
  ##            Macheras and Karalis.                                 ##
  ##   alpha0:  Nominal alpha and /primary/ level in stage 1 of all   ##
  ##           'Type 1' designs. Defaults to 0.05.                    ##
  ##   alpha1:  Mandatory adjustment for stage 1 of 'Type 1', 2nd ad- ##
  ##            justment in stage 1 of 'Type 2'. Defaults to the most ##
  ##            commonly used 0.0294.                                 ##
  ##   alpha2:  Adjustment in the final analysis. Defaults to the     ##
  ##            most commonly used 0.0294.                            ##
  ##   theta1:  Lower limit of the BE acceptance range; defaults to   ##
  ##            0.80.                                                 ##
  ##   theta2:  Upper limit of the BE acceptance range; defaults to   ##
  ##            1.25. If not given will be set to 1/theta1. The TIE   ##
  ##            will be simulated here.                               ##
  ##   target:  Target (desired) power for the estimation of the      ##
  ##            total sample size. Defaults to 0.80.                  ##
  ##   Xover:   If TRUE, 2x2x2 crossover design (default).            ##
  ##            if FALSE, two-group parallel (only Welch-test imple-  ##
  ##            mented so far).                                       ##
  ##   KM:      FALSE (default): /not/ Karalis/Macheras or Karalis;   ##
  ##            if set to TRUE additional input is needed below.      ##
  ##   KM.des:  "TSD"  : Karalis/Macheras (modified "C")              ##
  ##            "TSD-1": Karalis (modified "B")                       ##
  ##            "TSD-2": Karalis (modified "C")                       ##
  ##   pmethod: Power calculation method                              ##
  ##            "nct"    : Approximate calculations via non-central   ##
  ##                       t-distribution (default for speed reasons) ##
  ##            "shifted": Approximate calculation via shifted cent-  ##
  ##                       ral t-distribution like in most papers.    ##
  ##            "exact"  : Exact calculations via Owen's Q-functions  ##
  ##                       Warning: ~40times slower than the others.  ##
  ##   int.pwr: If TRUE (the default) uses the interim power step     ##
  ##            like in all [sic] papers.                             ##
  ##            FALSE: Do /not/ use interim power (Dutch MEB?).       ##
  ##            Possible only for 'Type 1' designs.                   ##
  ##   min.n2:  Minimum sample size in stage 2 (defaults to 0).       ##
  ##            Set to 2 according to Q&A-document (Rev. 7).          ##
  ##            Nonsense: Any study failing in the interim will       ##
  ##            /always/ have less than the target power. Therefore,  ##
  ##            the minimum stage 2 sample size will /always/ be 2    ##
  ##            anyway. If you like 'cosmetics' in the output, go     ##
  ##            ahead...                                              ##
  ##   max.n:   Defaults to unconstrainted total sample size (Inf).   ##
  ##            If set to a number, all studies with an estimated     ##
  ##            /higher/ total sample size will be forced to this     ##
  ##            number. Might degrade the power to show BE in the     ##
  ##            final analyis (see Fuglsang 2014) but once the study  ##
  ##            is done it is not of regulatory concern.              ##
  ##            Note: This is /not/ a futility criterion like in the  ##
  ##            TSDs of Karalis/Macheras and Karalis.                 ##
  ##   Nmax:    Defaults to unconstrainted total sample size (Inf).   ##
  ##            Can be set as a futility criterion for parallel       ##
  ##            designs and is mandatory for KM-designs.              ##
  ##   fCrit:   Futility criterion (early stopping) in the interim;   ##
  ##            Can be set to "PE" for the point estimate or "CI" for ##
  ##            the CI.                                               ##
  ##   fClow:   Lower futility limit of the PE or CI given above.     ##
  ##            Defaults to 0 (i.e., does /not/ stop for futility).   ##
  ##            Common values for "PE":                               ##
  ##              0.80: Armitage 1991, Karalis/Macheras, Karalis      ##
  ##              0.85: Charles Bon (AAPS Annual Meeting 2007)        ##
  ##            Common values for "CI":                               ##
  ##              0.925 : Personal communication Diane Potvin         ##
  ##              0.9374: Xu et al. "Method E" for CV 10-30%          ##
  ##              0.9305: Xu et al. "Method E" for CV 30-55%          ##
  ##              0.9492: Xu et al. "Method F" for CV 10-30%          ##
  ##              0.9350: Xu et al. "Method F" for CV 30-55%          ##
  ## ---------------------------------------------------------------- ##
  ## Simulation settings:                                             ##
  ## ---------------------------------------------------------------- ##
  ##   nsims:   Number of simulated studies. Don't use a lower number ##
  ##            than the default of 1e6 in a productive environment;  ##
  ##            only for testing!                                     ##
  ##            Note: Increasing nsims /rarely/ gives an adjusted     ##
  ##            alpha and respective TIE substantially closer to      ##
  ##            alpha0.                                               ##
  ##   setseed: If TRUE (default) uses a fixed seed of 123456.        ##
  ##            If FALSE, a random seed is used to check robustness.  ##
  ##   tol:     Convergence tolerance of uniroot. Defaults to 1e-8.   ##
  ##            Lower values /rarely/ give an adjusted alpha with TIE ##
  ##            closer to alpha0.                                     ##
  ##   pa:      Should the power analysis be printed/plotted?         ##
  ##            Defaults to FALSE.                                    ##
  ##   skip:    If set to TRUE (default) the alpha-optimization is    ##
  ##            /only/ performed if there is inflation of the TIE.    ##
  ##            The optimization will take some time. Be  patient;    ##
  ##            5-10 minutes are absolutely normal!                   ##
  ##   algo:    1 Optimize alpha2 by uniroot (default for speed).     ##
  ##            2 Set of alphas/TIES around the optimized alpha. Fit  ##
  ##              two models (linear, quadratic), select the better   ##
  ##              one (minimum AIC), and solve for TIE=alpha0. Doubles##
  ##              the computation time. IMHO, of doubtful value.      ##
  ##   plot.it: TRUE|FALSE. Should PE and CI of ajusted alpha be      ##
  ##            plotted? Defaults to FALSE.                           ##
  ##            Note: If no inflation of the TIE is observed with the ##
  ##            specified alpha(s) plotting will /not/ be done even   ##
  ##            if requested!                                         ##
  ##   valid:   TRUE|FALSE. Defaults to FALSE (run your own data).    ##
  ##            If set to TRUE, one of the built-in datasets must be  ##
  ##            specified in the next argument.                       ##
  ##   expl:    Specifies the validation dataset:                     ##
  ##            Defaults to 3 (Potvin 2008, Method B, Example 2)      ##
  ##            14 datasets are provided.                             ##
  ##   CIs:     If set to TRUE, the two-sided 95% CI of the TIE(s)    ##
  ##            will be shown. Defaults to FALSE.                     ##
  ## ================================================================ ##
  ## Known bugs:                                                      ##
  ##  - None so far.                                                  ##
  ## ================================================================ ##
  ## TODO (not covered yet):                                          ##
  ##  - Parallel groups:                                              ##
  ##    - Vectorize CVs of groups (currently assumed to be equal).    ##
  ##    - Unequal group sizes in the final analysis / Welch-test.     ##
  ##  - The code does not "know" whether a study was unbalanced       ##
  ##    (Xover) or performed in unequal group sizes (parallel). For   ##
  ##    odd sample sizes a best guess will be done. If f.i. N is 23   ##
  ##    the code will assume 12|11. However the code does not have a  ##
  ##    crystal ball. If N is 24 the true sequences might have been   ##
  ##    14|10. However, the CI (if more balance is incorrectly        ##
  ##    assumed) will be wider and thus conservative.                 ##
  ##    Remedy: Vectorized input of sample sizes.                     ##
  ##  - No alpha-spending in the first stage = (blinded) sample size  ##
  ##    re-estimation acc. to Golkowski et al. 2014. Available in     ##
  ##    function power.tsd.ssr(). Generally needs a /lot/ of ad-   ##
  ##    justment (contrary to what the BWSP 'believes').              ##
  ##    Too many different arguments for my taste - likely better to  ##
  ##    write specific code instead of incorporating it here.         ##
  ######################################################################
  # Warning: Change below only if you know what you are doing! #
  ##############################################################
  exec.start <- strftime(Sys.time(), usetz=TRUE) # Timestamp.
  graphics.off()                                 # Close eventual plots.
  # Check input for completeness and plausibility.
  if (!valid) { # Only if not validating.
    # Variabilities
    if (missing(Var1) | (length(Var1) != 2 & length(Var1) != 4))
      stop("Var1 must be a two- or four-element vector.")
    if (Var1[length(Var1)] != "CV" &
        Var1[length(Var1)] != "MSE" &
        Var1[length(Var1)] != "CI")
      stop("last element of Var1 must be \'CV\', \'MSE\', or \'CI\'.")
    if (stop1) Var <- PE <- N <- NA # Not needed; don't check.
    if (!stop1) { # Not needed if the study stopped in stage 1.
      if (missing(Var) | (length(Var) != 2 & length(Var) != 4))
        stop("Var must be a two- or four-element vector.")
      if (Var[length(Var)] != "CV" &
          Var[length(Var)] != "MSE" &
          Var[length(Var)] != "CI")
        stop("last element of Var must be \'CV\', \'MSE\', or \'CI\'.")
    }
    # PEs.
    if (missing(PE1) & length(Var1) != 4)
      stop("PE1 or stage 1 CI must be given.")
    if (!stop1) { # Not needed if the study stopped in stage 1.
      if (missing(PE) & length(Var) != 4)
        stop("PE or final CI must be given.")
    }
    # Acceptance range defaults.
    if (missing(theta1) & missing(theta2)) theta1 <- 0.8
    if (missing(theta2)) theta2 <- 1/theta1
    # Check GMR.
    if (missing(GMR)) GMR <- 0.95
    if (GMR < theta1 | GMR > theta2)
      stop("GMR must be within [theta1, theta2].")
    if (missing(pmethod)) pmethod <- "nct"
    pmethod <- match.arg(pmethod)
    # Check type of design.
    type <- as.character(type)
    if (!type %in% c("1", "2", "MSDBE")) stop("type must be 1, 2, or \'MSDBE\'.")
    # Futility criterion default and checking.
    if (missing(fCrit)) fCrit <- "PE"
    fCrit <- match.arg(fCrit)
    # Check target power.
    if (target > 1) stop("target must be <=1.")
    if (target < 0.5) warning("target <0.5 does not make sense.")
    # check algo
    if (algo < 1 | algo > 2) stop("algo must be 1 or 2.")
    # Check conditions of KM-designs.
    if (KM) {
      # Check PE futility.
      if (fClow == 0) {
        fCrit <- "PE"
        fClow <- 0.8
        warning("Lower futility criterion for PE set to 0.80.")
      }
      # Check Nmax futility in KM-designs.
      if (Nmax == Inf)
        warning("Should give a finite total sample size for futility in K/M-designs.")
      # Check KM type.
      KM.des <- match.arg(KM.des)
      # Check pmethod.
      if (pmethod == "shifted") {
        pmethod <- "nct"
        warning("pmethod 'shifted' not supported; changed to 'nct'.")
      }
    }
  }

  ########################################################
  # The workhorse for optimization. Called by uniroot(). #
  ########################################################
  opt <- function(x) {
    # Optimizes adjusted alpha(s) based on stage 1 data.
    # Objective function: TIE(x) - alpha0 = 0
    # Arguments         : Study specifications as given above.
    # Returns           : Adjusted alpha(s) giving TIE ~alpha0.
    # asymmetric split of alphas?
    if (asym) alpha <- c(alpha1, x) else alpha <- rep(x, 2)
    if (!KM) {  # The common methods.
      if (Xover) { # Crossover.
        power.tsd.fC(method=meth, alpha0=alpha0, alpha=alpha,
                    n1=n1, GMR=GMR, CV=CV1, targetpower=target,
                    pmethod=pmethod, usePE=usePE,
                    powerstep=int.pwr, min.n2=min.n2,
                    max.n=max.n, fCrit=fCrit, fClower=fClow,
                    theta0=theta2, nsims=nsims,
                    setseed=setseed)$pBE - alpha0
      } else {     # Parallel.
        power.tsd.p(method=meth, alpha0=alpha0, alpha=alpha,
                    n1=n1, GMR=GMR, CV=CV1, targetpower=target,
                    pmethod=pmethod, usePE=usePE, Nmax=Nmax,
                    test="welch", theta0=theta2, nsims=nsims,
                    setseed=setseed)$pBE - alpha0
      }
    } else {    # Adaptive (Karalis/Macheras and Karalis).
      power.tsd.KM(method=meth, alpha0=alpha0, alpha=alpha,
                   n1=n1, CV=CV1, targetpower=target,
                   pmethod=pmethod, Nmax=Nmax, theta0=theta2,
                   nsims=nsims, setseed=setseed)$pBE - alpha0
    }
  }

  ########################
  # Validation examples  #
  ########################
  if (valid) {
    # retrieve data
    # allow overruling some arguments to explore impact
    res <- example(expl)
    alpha0  <- res$alpha0
    if (missing(alpha1))  alpha1  <- res$alpha1
    if (missing(alpha2))  alpha2  <- res$alpha2
    if (missing(GMR))     GMR     <- res$GMR
    if (missing(target))  target  <- res$target
    theta1  <- res$theta1
    theta2  <- res$theta2
    if (missing(pmethod)) pmethod <- res$pmethod
    if (missing(type))     type    <- res$type
    Xover   <- res$Xover
    if (missing(usePE))   usePE   <- res$usePE
    if (missing(int.pwr)) int.pwr <- res$int.pwr
    if (missing(min.n2))  min.n2  <- res$min.n2
    if (missing(max.n))   max.n   <- res$max.n
    if (missing(Nmax))    Nmax    <- res$Nmax
    if (missing(fCrit))   fCrit   <- res$fCrit
    if (missing(fClow))   fClow   <- res$fClow
    if (missing(stop1))   stop1   <- res$stop1
    if (missing(setseed)) setseed <- res$setseed
    if (missing(nsims))   nsims   <- res$nsims
    if (missing(tol))     tol     <- res$tol
    if (missing(algo))    algo    <- res$algo
    if (missing(plot.it)) plot.it <- res$plot.it
    if (missing(skip))    skip    <- res$skip
    KM      <- res$KM
    KM.des  <- res$KM.des
    if (missing(CIs))     CIs     <- res$CIs
    Var1    <- res$Var1
    PE1     <- res$PE1
    n1      <- res$n1
    Var     <- res$Var
    PE      <- res$PE
    N       <- res$N
    info    <- res$info
  } # End of validation examples.

  # Create the "method"-argument for PowerTOST and Power2Stage.
  if (Xover) {
    design <- "2x2x2"
    if (length(n1) > 1) {
      n1 <- 2*round(mean(n1), 0)
      warning("n1 should be given as a scalar.\n  Mean will be used.")
    }
  } else {
    if (KM) stop("KM only for crossover designs.")
    design <- "parallel"
    if (length(n1) < 2) {
      n1 <- nvec(n1, 2)
      warning("n1 should be given for each group (T and R).\n  n1 was vectorized.")
    }
  }

  # CVs from variabilities (given as CV, MSE, or CI).
  # Stage 1 data.
  if (Var1[length(Var1)] == "CV") {
    CV1 <- as.double(Var1[1:length(Var1)-1])
  } else if (Var1[length(Var1)] == "MSE") {
    CV1 <- mse2CV(as.double(Var1[1:length(Var1)-1]))
  } else {
    # Suppress messages regarding unbalanced designs.
    suppressMessages(
      CV1 <- CVfromCI(lower=as.double(Var1[1]), upper=as.double(Var1[2]),
                      alpha=as.double(Var1[3]), n=n1, design=design) )
  }
  # Pooled data.
  if(!stop1) { # Not needed if stopped.
    if (Var[length(Var)] == "CV") {
      CV <- as.double(Var[1:length(Var)-1])
    } else if (Var[length(Var)] == "MSE") {
      CV <- mse2CV(as.double(Var[1:length(Var)-1]))
    } else { # from CI
      # Suppress messages regarding unbalanced designs.
      suppressMessages(
        CV <- CVfromFinalCI(lower=as.double(Var[1]), upper=as.double(Var[2]),
                            alpha=as.double(Var[3]), n=N, design=design) )
    }
  }

  # Transform PEs or calculate from the CIs.
  # Stage 1 data.
  if (!missing(PE1)) {
    if (PE1[2] == "ratio") {
      PE1 <- as.double(PE1[1])
    } else {
      PE1 <- exp(as.double(PE1[1]))
    }
  } else {
    PE1 <- sqrt(as.double(Var1[1]) * as.double(Var1[2]))
  }
  if (!stop1) { # Not needed if stopped in the interim.
    if (!missing(PE)) {
      if (PE[2] == "ratio") {
        PE <- as.double(PE[1])
      } else {
        PE <- exp(as.double(PE[1]))
      }
    } else {
      PE <- sqrt(as.double(Var[1]) * as.double(Var[2]))
    }
  }

  # Significance limit of the binomial test for alpha0 and No. of sim's.
  sig    <- binom.test(alpha0*nsims, nsims, alternative="less")$conf.int[2]

  # Information about the computing environment.
  system <- Sys.info()
  node   <- system["nodename"]
  user   <- system["user"]
  OS     <- system["sysname"]
  OSrel  <- system["release"]
  OSver  <- system["version"]
  rver   <- sessionInfo()$R.version$version.string
  rver   <- substr(rver, which(strsplit(rver, "")[[1]]=="n")+2, nchar(rver))
  env    <- as.character(OS)
  flushable <- FALSE
  if (env == "Windows" | env == "Darwin") flushable <- TRUE
  cit    <- citation("Power2Stage")
  year1  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  cit    <- citation("PowerTOST")
  year2  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  cit    <- citation("AdaptiveBE")
  year3  <- paste0(" (", substr(cit, regexpr("year", cit)+8,
                                regexpr("year", cit)+11), ")")
  hr     <- paste0(rep("\u2500", 61), collapse="")
  if (pmethod == "nct") pverbose <- "(approx. via non-central t)"
  if (pmethod == "shifted") pverbose <- "(approx. via shifted central t)"
  if (pmethod == "exact") pverbose <- "(exact via Owen\u2019s Q)"
  # Splash screens.
  if (stop1) {
    info1 <- paste(
      "\n====================================================",
      "\n One million BE studies in a TSD will be simulated.",
      "\n Should take less than one minute on most systems.",
      "\n====================================================\n")
  } else {
    if (pmethod != "exact") {
      info1 <- paste(
        "\n====================================================",
        "\n 1.1 million BE studies in a TSD will be simulated.",
        "\n Should take less than one minute on most systems.",
        "\n====================================================\n")
    } else {
      info1 <- paste(
        "\n====================================================",
        "\n 1.1 million BE studies in a TSD will be simulated.",
        "\n Should take less than one minute on most systems.",
        "\n Will take a couple of minutes on most systems.   ",
        "\n====================================================\n")
    }
  }
  if (pmethod != "exact") {
    info2 <- paste(
      "\n=============================================================",
      "\n Be patient. Adjusting \u03B1 requires simulating of up to 20",
      "\n million BE studies in a TSD. Will take a couple of minutes\u2026 ",
      "\n=============================================================\n")
  } else {
    info2 <- paste(
      "\n=============================================================",
      "\n Be patient. Adjusting \u03B1 requires simulating of up to 20",
      "\n million BE studies in a TSD. Can take many (!) hours\u2026       ",
      "\n=============================================================\n")
  }

  # Check for asymmetric split of alphas (set a boolean variable).
  if (alpha1 == alpha2) {
    asym <- FALSE # Same alphas in both stages.
  } else {
    asym <- TRUE  # Haybittle/Peto, O'Brien/Fleming, Zheng et al., Xu et al.
  }

  # Convert type 1 or 2 to argument "B" or "C" for the function-calls.
  if (type == "MSDBE") {
    meth <- "B0"
  } else {
    meth <- intToUtf8(65 + as.integer(type))
  }

  # Assessment text.
  pass <- "(BE concluded)"
  fail <- "(failed to demonstrate BE)"

  # Try to reconstruct methods from type and alpha levels.
  des.type <- type
  outside.range <- FALSE
  if (KM) {
    des.type <- paste0(type, " (KM ", KM.des, ")")
  } else {
    if (type == 1) {
      if (asym) {
        if (alpha1 == 0.001 & alpha2 == 0.049)
          des.type <- paste(type, "(Haybittle/Peto)")
        if (alpha1 == 0.005 & alpha2 == 0.048)
          des.type <- paste(type, "(O\u2019Brien/Fleming)")
        if (alpha1 == 0.01 & alpha2 == 0.04 & min.n2 == 0)
          des.type <- paste(type, "(Zheng et al. 2015, MSDBE)")
        if (alpha1 == 0.0249 & alpha2 == 0.0363 & min.n2 == 0 &
            fClow == 0.9374) {
          des.type <- paste(type, "(Xu et al. 2015, Method E, low CVs)")
          if (CV1 > 0.3 | n1 < 12 | GMR != 0.95 | target != 0.8 |
              n1 > 30 | max.n > 42) outside.range <- TRUE
        }
        if (alpha1 == 0.0254 & alpha2 == 0.0357 & fClow == 0.9305) {
          des.type <- paste(type, "(Xu et al. 2015, Method E, high CVs)")
          if (CV1 < 0.3 | CV1 > 0.55 | n1 < 48 | GMR != 0.95 | target != 0.8 |
              n1 > 60 | max.n > 180) outside.range <- TRUE
        }
        if (nchar(des.type) == 1) des.type <- paste(type, "(user specified)")
      } else {
        if (unique(c(c(alpha1, alpha2))) == 0.0294 & min.n2 == 0) {
          if (design == "2x2x2") {
            des.type <- paste(type, "(Potvin et al. 2008, Method B)")
            if (CV1 < 0.1 | CV1 > 1 | n1 < 12 | n1 > 60 | GMR != 0.95 |
                target != 0.8) outside.range <- TRUE
          } else {
            des.type <- paste(type, "(Fuglsang 2014, Method B)")
            if (CV1 < 0.1 | CV1 > 1 | sum(n1) < 48 | sum(n1) > 120 |
              GMR != 0.95 | target != 0.8) outside.range <- TRUE
          }
        }
        if (unique(c(c(alpha1, alpha2))) == 0.0304 & min.n2 == 0)
          des.type <- paste(type, "(Pocock equivalence, Method B)")
        if (unique(c(c(alpha1, alpha2))) == 0.0284 & min.n2 == 0) {
          des.type <- paste(type, "(Fuglsang 2013, Method B)")
          if (CV1 < 0.1 | CV1 > 0.8 | n1 < 12 | n1 > 60 | GMR != 0.95 |
              target != 0.9) outside.range <- TRUE
        }
        if (nchar(des.type) == 1) des.type <- paste(type, "(user specified)")
      }
    } else {
      if (asym) {
        if (alpha1 == 0.0248 & alpha2 == 0.0364 & min.n2 == 0 &
            fClow == 0.9492) {
          des.type <- paste(type, "(Xu et al. 2015, Method F, low CVs)")
          if (CV1 > 0.3 | n1 < 18 | GMR != 0.95 | target != 0.8 |
              n1 > 30 | max.n > 42) outside.range <- TRUE
        }
        if (alpha1 == 0.0259 & alpha2 == 0.0349 & min.n2 == 0 &
            fClow == 0.9350) {
          des.type <- paste(type, "(Xu et al. 2016, Method F, high CVs)")
          if (CV1 < 0.3 | CV1 > 0.55 | n1 < 48 | GMR != 0.95 | target != 0.8 |
              n1 > 60 | max.n > 180) outside.range <- TRUE
        }
        if (nchar(des.type) == 1) des.type <- paste(type, "(user specified)")
      } else {
        if (unique(c(c(alpha1, alpha2))) == 0.0294 & min.n2 == 0) {
          if (design == "2x2x2") {
            des.type <- paste(type, "(Potvin et al. 2008, Method C)")
            if (CV1 < 0.1 | CV1 > 1 | n1 < 12 | n1 > 60 | GMR != 0.95 |
                target != 0.8) outside.range <- TRUE
          } else {
            des.type <- paste(type, "(Fuglsang 2014, Method C)")
            if (CV1 < 0.1 | CV1 > 1 | sum(n1) < 48 | sum(n1) > 120 |
              GMR != 0.95 | target != 0.8) outside.range <- TRUE
          }
        }
        if (unique(c(c(alpha1, alpha2))) == 0.0304 & min.n2 == 0)
          des.type <- paste(type, "(Pocock equivalence, Method C)")
        if (unique(c(c(alpha1, alpha2))) == 0.0280 & min.n2 == 0) {
          des.type <- paste(type, "(Montague et al. 2011, Method D)")
          if (CV1 < 0.1 | CV1 > 1 | n1 < 12 | n1 > 60 | GMR != 0.90 |
              target != 0.8) outside.range <- TRUE
        }
        if ((unique(c(c(alpha1, alpha2))) == 0.0274 |
             unique(c(c(alpha1, alpha2))) == 0.0269) & min.n2 == 0) {
          des.type <- paste(type, "(Fuglsang 2013, Method C/D)")
          if (CV1 < 0.1 | CV1 > 0.8 | n1 < 12 | n1 > 60 | GMR != 0.95 |
              target != 0.9) outside.range <- TRUE
        }
        if (nchar(des.type) == 1) des.type <- paste(type, "(user specified)")
      }
    }
  }
  if (type == "MSDBE") {
    if (alpha1 == 0.01 & alpha2 == 0.04)
      des.type <- paste(type, "(Zheng 2015)")
  }

  # Fixed GMR or fully adaptive?
  if (!usePE | !KM) { # fixed
    GMR.used <- paste(sprintf("%.2f", GMR), "(fixed)")
  } else {             # adaptive
    GMR.used <- "PE observed in stage 1 (fully adaptive)"
  }

  # Interim power check performed?
  if (int.pwr) { # default
    pwr.check <- "yes"
  } else {       # Dutch MEB?
    pwr.check <- "no"
  }

  # Futility criterion used?
  if (fClow == 0) {
    fC.check <- "none"
  } else {
    if (fCrit == "PE" | KM) {
      fC.check <- paste0("PE outside [",
                         sprintf("%.4f, %.4f", fClow, 1/fClow), "]")
    } else {
      fC.check <- paste0("CI outside [",
                         sprintf("%.4f, %.4f", fClow, 1/fClow), "]")
    }
  }

  # Minimum n2 specified?
  if (min.n2 == 0) {   # all (!) papers.
    n2.check <- "not specified"
  } else {
    if (min.n2 == 2) { # Q&A.
      n2.check <- "2 (EMA Q&A document Rev. 7)"
    } else {
      n2.check <- min.n2
    }
  }

  # Maximum sample size specified?
  if (max.n == Inf & Nmax == Inf) { # default
    maxN.check <- "not specified"
  } else {
    if (!KM) {
      if (Xover) {
        maxN.check <- max.n
      } else {
        maxN.check <- Nmax
      }
    } else {
      maxN.check <- Nmax
    }
  }
  txt <- paste(paste0("\n", hr,
                      "\nSystem             : ", node,
                      "\nUser               : ", user,
                      "\nOperating System   : ", OS, " ", OSrel))
  if (OS == "Darwin") { # special treatment (long system[["version"]])
    tmp <- strwrap(OSver, width=61, prefix="\n                     ")
    for (j in 1:length(tmp)) {
      txt <- paste0(txt, tmp[[j]])
    }
  } else {
    txt <- paste0(txt, " ", OSver)
  }
  txt <- paste(paste0(txt, "\nR version          : ", rver,
                      "\nPower2Stage version: ", packageVersion("Power2Stage"), year1,
                      "\nPowerTOST version  : ", packageVersion("PowerTOST"), year2,
                      "\nAdaptiveBE version : ", packageVersion("AdaptiveBE"), year3,
                      "\ncheck.TSD() started: ", exec.start, "\n"))

  if (valid) {
    txt <- paste(txt, "\nValidation;")
    txt <- paste(txt, info)
  }
  cond <- paste(txt,
                "\nStudy conditions and assessment of empiric Type I Error",
                paste0("\n", paste0(rep("\u2500", 61), collapse="")),
                "\nDesign             :")
  if (design == "2x2x2") {
    cond <- paste(cond, "2\u00D72\u00D72 crossover")
  } else {
    cond <- paste(cond, "parallel")
  }
  cond <- paste(cond, "\nTSD Type           :", des.type)
  if (outside.range) {
    cond <- paste(cond,
                  "\nWarning            : At least one of the study conditions",
                  "\n                     is outside the validated range of the",
                  "\n                     reference\u2019s method.")
  }
  cond <- paste(cond,
                "\nTarget power       :", sprintf("%.2f", target),
                "\nGMR used           :", GMR.used,
                "\nInterim power check:", pwr.check,
                "\nFutility criterion :", fC.check,
                "\nMinimum n2         :", n2.check,
                "\nMaximum N          :", maxN.check)
  if (asym) {
    if (type == 1) {
      cond <- paste(cond, "\nSpecified \u03B1 1, 2   :",
                    sprintf("%.4f, %.4f", alpha1, alpha2),
                    "\nSpecified CIs      :",
                    sprintf("%.2f%%, %.2f%%",
                            100*(1-2*alpha1), 100*(1-2*alpha2)))
    } else {
      cond <- paste(cond, "\nSpecified \u03B1 1, 2   :",
                    paste0(sprintf("%.3f", alpha0), " | ",
                           sprintf("%.4f, %.4f", alpha1, alpha2)),
                    "\nSpecified CIs      :",
                    paste0(sprintf("%.2f%%", 100*(1-2*alpha0)), " | ",
                           sprintf("%.2f%%, %.2f%%",
                                   100*(1-2*alpha1), 100*(1-2*alpha2))))
    }
  } else {
    if (type == 1) {
      cond <- paste(cond, "\nSpecified \u03B1 1, 2   :",
                    sprintf("%.4f, %.4f", alpha1, alpha2),
                    "\nSpecified CIs      :",
                    sprintf("%.2f%%, %.2f%%",
                            100*(1-2*alpha1), 100*(1-2*alpha2)))
    } else {
      cond <- paste(cond, "\nSpecified \u03B1 1, 2   :",
                    paste0(sprintf("%.3f", alpha0), " | ",
                           sprintf("%.4f, %.4f", alpha1, alpha2)),
                    "\nSpecified CIs      :",
                    paste0(sprintf("%.2f%%", 100*(1-2*alpha0)), " | ",
                           sprintf("%.2f%%, %.2f%%",
                                   100*(1-2*alpha1), 100*(1-2*alpha2))))
    }
  }

  CI2.txt <- CI1.txt <- ""
  ptm <- proc.time() # Start timer.
  # Estimate the empiric TIE and power for the specified alpha(s)
  # based on the interim data.
  if (flushable) {
    cat(info1)
    flush.console()
  }
  if (!stop1) { # Study proceeded to the second stage.
    TIE.interim.est <- TIE.est(meth, alpha0, alpha1, alpha2, n1, GMR,
                               CV1, target, pmethod, usePE, int.pwr,
                               min.n2, max.n, Nmax, fCrit, fClow,
                               theta2, nsims, setseed, KM, asym, Xover)
    # 95% CI of empiric TIE
    TIE.CI1 <- binom.test(TIE.interim.est*nsims, nsims,
                          alternative="two.sided")$conf.int[1:2]
    if (pa) {
      pwr.interim.est <- pwr.est(meth, alpha0, alpha1, alpha2, n1, GMR,
                                 CV1, target, pmethod, usePE, int.pwr,
                                 min.n2, max.n, Nmax, fCrit, fClow,
                                 nsims, setseed, KM, asym, Xover)
    }
  } else {      # Study stopped in the interim.
    TIE.interim.est <- TIE1.est(meth, alpha0, alpha1, alpha2, n1, GMR,
                                CV1, target, pmethod, usePE, int.pwr,
                                min.n2, max.n, Nmax, fCrit, fClow,
                                theta2, nsims, setseed, KM, asym, Xover)
    # 95% CI of empiric TIE
    TIE.CI1 <- binom.test(TIE.interim.est*nsims, nsims,
                          alternative="two.sided")$conf.int[1:2]
    # Inflated TIE for specified alpha(s)?
    justif <- paste("\n\nTIE for specified \u03B1:", sprintf("%1.5f",
                                                             TIE.interim.est))

    if (TIE.interim.est <= alpha0) { # yes
      justif <- paste(justif, "(\u22640.05)", "\n                    ",
                      "Applied adjustment is justified.")
    } else {                         # no (BSWP rule: >0.05 not acceptable)
      justif <- paste(justif, "(>0.05)", "\n                    ",
                      "Applied adjustment is not justified.")
      if (TIE.interim.est <= sig) {  # significant inflation of the TIE?
        justif <- paste(justif, "\n                    ",
                        "However, the TIE is n.s. >0.05.")
      }
    }
    if (CIs) {
      CI1.txt <- paste("\n                     95% CI:",
                       sprintf("%1.5f%s%1.5f", TIE.CI1[1], "\u2013", TIE.CI1[2]))
    }
    if (type == 1) { # Method B, E, KM TSD-1, MSDBE, Haybittle/Peto,...
      interim <- paste0("\n\nInterim analysis (specified \u03B11 ", alpha1, ")")
      interim <- paste0(interim, paste0("\n",
                                        paste0(rep("\u2500", 52), collapse="")))
      # In Type 1 BE is based on alpha1
      BE.interim.spec <- round(100*CI.BE(alpha=alpha1, pe=PE1, CV=CV1,
                                         n=n1, design=design), 2)
      if (BE.interim.spec[["lower"]] >= 80 &
          BE.interim.spec[["upper"]] <= 125) {
        BE.interim.spec.assess <- pass
      } else {
        BE.interim.spec.assess <- fail
      }
      interim <- paste(interim,
                       sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha1), " CI:"),
                       sprintf("%.2f%s%.2f%%", BE.interim.spec[["lower"]], "\u2013",
                               BE.interim.spec[["upper"]]),
                       BE.interim.spec.assess,
                       "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
    } else {         # Method C, D, C/D, F, KM TSD, KM TSD-2
      if (!KM) { # Conventional (fixed).
        pwr.interim.calc <- power.TOST(alpha=alpha0, theta0=GMR, CV=CV1,
                                       n=n1, method=pmethod, design=design)
      } else {   # Adaptive.
        pwr.interim.calc <- power.TOST(alpha=alpha0, theta0=PE1, CV=CV1,
                                       n=n1, method=pmethod, design=design)
      }
      interim <- "\n\nInterim analysis (specified \u03B1"
      if (pwr.interim.calc >= target) {
        interim <- paste0(interim, "0 ", alpha0, ")")
      } else {
        interim <- paste0(interim, "1 ", alpha1, ")")
      }
      interim <- paste0(interim, paste0("\n",
                                        paste0(rep("\u2500", 52), collapse="")))
      # If at least the target power, assess BE with alpha0.
      if (pwr.interim.calc >= target) {
        BE.interim.spec0 <- round(100*CI.BE(alpha=alpha0, pe=PE1, CV=CV1,
                                            n=n1, design=design), 2)
        if (BE.interim.spec0[["lower"]] >= 80 &
            BE.interim.spec0[["upper"]] <= 125) {
          BE.interim.spec.assess <- pass
        } else {
          BE.interim.spec.assess <- fail
        }
        interim <- paste(interim,
                         sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha0), " CI:"),
                         sprintf("%.2f%s%.2f%%", BE.interim.spec0[["lower"]], "\u2013",
                                 BE.interim.spec0[["upper"]]),
                         BE.interim.spec.assess,
                         "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
        # If less than the target power, assess BE with alpha1.
      } else {
        BE.interim.spec1 <- round(100*CI.BE(alpha=alpha1, pe=PE1, CV=CV1,
                                            n=n1, design=design), 2)
        if (BE.interim.spec1[["lower"]] >= 80 &
            BE.interim.spec1[["upper"]] <= 125) {
          BE.interim.spec.assess <- pass
        } else {
          BE.interim.spec.assess <- fail
        }
        interim <- paste(interim,
                         sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha1), " CI:"),
                         sprintf("%.2f%s%.2f%%", BE.interim.spec1[["lower"]], "\u2013",
                                 BE.interim.spec1[["upper"]]),
                         BE.interim.spec.assess,
                         "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
      }
    }
    txt <- paste(cond, "\n\nData for the interim analysis",
                 paste0("\n", paste0(rep("\u2500", 41), collapse="")),
                 "\nCV (MSE)           :", sprintf("%6.2f%% (%.6g)", 100*CV1,
                                                   CV2mse(CV1)),
                 "\nPE (ln(T)\u2013ln(R))   :", sprintf("%6.2f%% (%.6g)", 100*PE1,
                                                        log(PE1)),
                 "\nSample size        :", sprintf("%3i", sum(n1)))
    txt <- paste(txt, "\n\nStudy stopped in stage 1 (interim \u2192 pivotal).",
                 justif)
    txt <- paste(txt, interim,
                 paste0("\n\n\u250C", paste0(rep("\u2500", 53), collapse=""), "\u2510",
                        "\n\u2502 Since no inflation of the Type I Error is expected, \u2502",
                        "\n\u2502 can accept the reported analysis.                   \u2502",
                        "\n\u2514", paste0(rep("\u2500", 53), collapse=""), "\u2518\n"))
    cat(txt, "\n")
    return(invisible("done"))
  }

  # Only if needed / requested.
  if (TIE.interim.est > alpha0 | (TIE.interim.est <= alpha0 & !skip)) {
    if (flushable) {
      cat(info2)
      flush.console()
    }
    # Optimize alpha in the interval {tol, alpha0}.
    # Convergence tolerance in the function-call: 1e-8 (or smaller).
    x <- uniroot(opt, interval=c(tol, alpha0), tol=tol)
    run.time1 <- proc.time() - ptm
    # Vector of adjusted alphas.
    if (!asym) { # Same adjusted alphas in both stages.
      adj <- rep(x$root, 2)
    } else {     # Alpha1 in stage 1; set adj. alpha2 for stage 2.
      adj <- c(alpha1, x$root)
    }
    # Empiric TIE for adjusted alpha(s).
    TIE <- alpha0 + x$f.root
    if (algo == 2) {
      ptm       <- proc.time() # Start timer.
      delta     <- abs(alpha2-adj[2]) / 4
      alpha.set <- seq(adj[2]-delta, adj[2]+delta, length.out=24)
      TIE.set   <- rep(NA, length(alpha.set))
      pb <- txtProgressBar(0, 1, 0, char="\u2588", width=NA, style=3)
      for (j in seq_along(alpha.set)) {
        if (asym) {
          alpha12 <- c(alpha1, alpha.set[j])
        } else {
          alpha12 <- rep(alpha.set[j], 2)
        }
        TIE.set[j] <- TIE.est(meth, alpha0, alpha12[1], alpha12[2], n1,
                              GMR, CV1, target, pmethod, usePE, int.pwr,
                              min.n2, max.n, Nmax, fCrit, fClow, theta2,
                              nsims, setseed, KM, asym, Xover)
        setTxtProgressBar(pb, j/length(alpha.set))
      }
      close(pb)
      alpha.set <- c(alpha.set, adj[2]) # Add the result
      TIE.set   <- c(TIE.set, TIE)      # we already have.
      dev.new(record=TRUE)
      op <- par(no.readonly=TRUE) # save par() options
      par(mar=c(c(4, 4, 2, 0))+0.1, cex.main=1.1, cex.axis=0.95)
      plot(alpha.set, TIE.set, type="n", xlab="adjusted \u03B1",
           ylab="empiric TIE",
           main="Type I Errors close to optimized \u03B1")
      grid()
      abline(h=c(alpha0, sig), lty=c(1, 2), col=c("black", "red"))
      abline(v=x$root, lty=2, col="blue")
      # Fits (linear: Fuglsang 2011; quadratic sometimes better)
      mod1 <- lm(TIE.set ~ alpha.set)
      mod2 <- lm(TIE.set ~ alpha.set + I(alpha.set^2))
      # Select the better model based on lower AIC.
      if (extractAIC(mod1, k=2)[2] <= extractAIC(mod2, k=2)[2]) {
        fit    <- "linear"
        best   <- (alpha0-coef(mod1)[[1]])/coef(mod1)[[2]]
        points(sort(alpha.set), sort(fitted(mod1)), type="l", col="red")
      } else {
        fit    <- "quadratic"
        det    <- sqrt((coef(mod2)[[2]]/2/coef(mod2)[[3]])^2-
                         (coef(mod2)[[1]]-alpha0)/coef(mod2)[[3]])
        if (coef(mod2)[[3]] < 0) { # convex shape.
          best <- -(coef(mod2)[[2]]/2/coef(mod2)[[3]] + det)
        } else {                   # concave shape.
          best <- -(coef(mod2)[[2]]/2/coef(mod2)[[3]] - det)
        }
        points(sort(alpha.set), sort(fitted(mod2)), type="l", col="red")
      }
      abline(v=best)
      legend("topleft", inset=0.02, bg="white", box.lty=0, y.intersp=1.25,
             legend=c(
               paste0("\u03b1-range (fixed seed; ", fit, " fit)"),
               "optimized \u03b1",
               "fitted \u03b1",
               "significance limit",
               "target"),
             pch=c(rep(21, 3), NA, NA),
             pt.cex=c(1.25, rep(2, 2), 0, 0),
             col=c("#FF0000AA", "blue", "black", "red", "black"),
             pt.bg=c("#FFD700AA", rep(NA, 4)),
             lty=c(1, 1, 1, 2, 1))
      points(alpha.set, TIE.set, cex=1.25, pch=21,
             col="#FF0000AA", bg="#FFD700AA")
      points(x$root, TIE, cex=2, pch=21, col="blue", bg=NA)
      if (!asym) { # Same adjusted alphas in both stages.
        adj <- rep(best, 2)
      } else {     # Alpha1 in stage 1 and set adj. alpha2 for stage 2.
        adj <- c(alpha1, best)
      }
      TIE <- TIE.est(meth, alpha0, adj[1], adj[2], n1, GMR, CV1, target,
                     pmethod, usePE, int.pwr, min.n2, max.n, Nmax, fCrit,
                     fClow, theta2, nsims, setseed, KM, asym, Xover)
      points(adj[2], TIE, cex=2, pch=21, col="black", bg=NA)
      par(op) # reset options
      run.time2 <- proc.time() - ptm
    }
    if (TIE <= sig) { # no
      infl <- "(n.s. >0.05)"
    } else {          # yes (Unlikely! Try lower tol-argument.)
      infl <- "(significantly >0.05)"
    }
    # 95% CI of empiric TIE for adjusted alpha
    TIE.CI2 <- binom.test(TIE*nsims, nsims, alternative="two.sided")$conf.int[1:2]
    if (CIs) {
      CI2.txt <- paste("\n                     95% CI:",
                       sprintf("%1.5f%s%1.5f", TIE.CI2[1], "\u2013", TIE.CI2[2]))
    }
    pwr.interim.adj <- pwr.est(meth, alpha0, adj[1], adj[2], n1, GMR, CV1,
                               target, pmethod, usePE, int.pwr, min.n2,
                               max.n, Nmax, fCrit, fClow, nsims, setseed,
                               KM, asym, Xover)
  }
  # Collect data for console output.
  # Inflated TIE for specified alpha(s)?
  justif <- paste("\n\nTIE for specified \u03B1:", sprintf("%1.5f",
                                                           TIE.interim.est))
  if (TIE.interim.est <= alpha0) { # yes
    justif <- paste(justif, "(\u22640.05)", "\n                    ",
                    "Applied adjustment is justified.")
  } else {                         # no (BSWP rule: >0.05 not acceptable)
    justif <- paste(justif, "(>0.05)", "\n                    ",
                    "Applied adjustment is not justified.")
    if (TIE.interim.est <= sig) {  # significant inflation of the TIE?
      justif <- paste(justif, "\n                    ",
                      "However, the TIE is n.s. >0.05.")
    }
  }
  if (CIs) {
    CI1.txt <- paste("\n                     95% CI:",
                     sprintf("%1.5f%s%1.5f", TIE.CI1[1], "\u2013", TIE.CI1[2]))
  }

  #####################
  # Interim analysis. #
  ################################################
  # First part: based on the specified alpha(s). #
  ################################################
  if (type == 1) {
    pwr.interim.calc <- power.TOST(alpha=alpha2, theta0=GMR, CV=CV1, n=n1,
                                   method=pmethod, design=design)
    interim <- paste0("\n\nInterim analysis (specified \u03B11 ", alpha1, ")")
    interim <- paste0(interim, paste0("\n",
                                      paste0(rep("\u2500", 52), collapse="")))
    # In Type 1 BE is based on alpha1
    BE.interim.spec <- round(100*CI.BE(alpha=alpha1, pe=PE1, CV=CV1,
                                       n=n1, design=design), 2)
    if (BE.interim.spec[["lower"]] >= 80 &
        BE.interim.spec[["upper"]] <= 125) {
      BE.interim.spec.assess <- pass
    } else {
      BE.interim.spec.assess <- fail
    }
    interim <- paste(interim,
                     sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha1), " CI:"),
                     sprintf("%.2f%s%.2f%%", BE.interim.spec[["lower"]], "\u2013",
                             BE.interim.spec[["upper"]]),
                     BE.interim.spec.assess,
                     "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
    if (pwr.interim.calc >= target & BE.interim.spec.assess == fail) {
      interim <- paste(interim,
                       "\nStudy should have been stopped (\u2265 target power) and",
                       "\nconclusions based on the interim.\n")
    } else { # Lower than target power and not BE. Proceed to stage 2.
      if (!usePE) { # Use fixed GMR.
        n2.spec <- sampleN2.TOST(alpha=alpha2, CV=CV1, n1=n1, theta0=GMR,
                                 theta1=theta1, theta2=theta2,
                                 targetpower=target, design=design,
                                 method=pmethod)[["Sample size"]]
      } else {      # Use the PE of stage 1.
        n2.spec <- sampleN2.TOST(alpha=alpha2, CV=CV1, n1=n1, theta0=PE1,
                                 theta1=theta1, theta2=theta2,
                                 targetpower=target, design=design,
                                 method=pmethod)[["Sample size"]]
      }
      interim <- paste(interim, sprintf("%s %i %s%i%s",
                                        "\nSecond stage with", n2.spec,
                                        "subjects (N=", n1 + n2.spec,
                                        ") is justified.\n"))
    }
  } else { # Type 2
    if (!KM) { # Conventional (fixed).
      pwr.interim.calc <- power.TOST(alpha=alpha0, theta0=GMR, CV=CV1,
                                     n=n1, method=pmethod, design=design)
    } else {   # Adaptive.
      pwr.interim.calc <- power.TOST(alpha=alpha0, theta0=PE1, CV=CV1,
                                     n=n1, method=pmethod, design=design)
    }
    interim <- "\n\nInterim analysis (specified \u03B1"
    if (pwr.interim.calc >= target) {
      interim <- paste0(interim, "0 ", alpha0, ")")
    } else {
      interim <- paste0(interim, "1 ", alpha1, ")")
    }
    interim <- paste0(interim, paste0("\n",
                                      paste0(rep("\u2500", 52), collapse="")))
    # If at least the target power, assess BE with alpha0.
    if (pwr.interim.calc >= target) {
      BE.interim.spec0 <- round(100*CI.BE(alpha=alpha0, pe=PE1, CV=CV1,
                                          n=n1, design=design), 2)
      if (BE.interim.spec0[["lower"]] >= 80 &
          BE.interim.spec0[["upper"]] <= 125) {
        BE.interim.spec.assess <- pass
      } else {
        BE.interim.spec.assess <- fail
      }
      interim <- paste(interim,
                       sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha0), " CI:"),
                       sprintf("%.2f%s%.2f%%", BE.interim.spec0[["lower"]], "\u2013",
                               BE.interim.spec0[["upper"]]),
                       BE.interim.spec.assess,
                       "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
      # If less than the target power, assess BE with alpha1.
    } else {
      BE.interim.spec1 <- round(100*CI.BE(alpha=alpha1, pe=PE1, CV=CV1,
                                          n=n1, design=design), 2)
      if (BE.interim.spec1[["lower"]] >= 80 &
          BE.interim.spec1[["upper"]] <= 125) {
        BE.interim.spec.assess <- pass
      } else {
        BE.interim.spec.assess <- fail
      }
      interim <- paste(interim,
                       sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha1), " CI:"),
                       sprintf("%.2f%s%.2f%%", BE.interim.spec1[["lower"]], "\u2013",
                               BE.interim.spec1[["upper"]]),
                       BE.interim.spec.assess,
                       "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
    }
    if (pwr.interim.calc >= target | BE.interim.spec.assess == pass) {
      if (pwr.interim.calc >= target) {
      interim <- paste(interim,
                       "\nStudy should have been stopped (\u2265 target power) and",
                       "\nconclusions based on the interim.\n")
      } else {
      interim <- paste(interim,
                       "\nStudy should have been stopped (BE in the interim).\n")
      }
    } else { # Lower than target power and not BE: Proceed to stage 2.
      if (!usePE) { # Use fixed GMR.
        n2.spec <- sampleN2.TOST(alpha=alpha2, CV=CV1, n1=n1, theta0=GMR,
                                 theta1=theta1, theta2=theta2,
                                 targetpower=target, design=design,
                                 method=pmethod)[["Sample size"]]
      } else {      # Use the PE of stage 1.
        n2.spec <- sampleN2.TOST(alpha=alpha2, CV=CV1, n1=n1, theta0=PE1,
                                 theta1=theta1, theta2=theta2,
                                 targetpower=target, design=design,
                                 method=pmethod)[["Sample size"]]
      }
      interim <- paste(interim, sprintf("%s %i %s%i%s",
                                        "\nSecond stage with", n2.spec,
                                        "subjects (N=", n1 + n2.spec,
                                        ") is justified.\n"))
    }
  }

  # Check whether the study should have been stopped for futility.
  if (fC.check != "none") {
    if ((KM | design == "parallel") & (n1+n2.spec > Nmax)) {
      interim <- paste(interim, paste0(
        "\nTotal sample size (", n1+n2.spec, ") > Nmax (", Nmax, "); the",
        "\nstudy should have been stopped for futility.\n"))
    }
    if ((fCrit == "PE" | KM) & (PE1 < fClow | PE1 > 1/fClow)) {
      interim <- paste(interim, paste0(
        "\nPE outside specified range of [",
        sprintf("%.4f, %.4f", fClow, 1/fClow), "];",
        "\nthe study should have been stopped for futility.\n"))
    }
    if (fCrit == "CI") {
      if (type == 1) {
        if (BE.interim.spec[["upper"]] < 100*fClow |
            BE.interim.spec[["lower"]] > 100/fClow) {
          interim <- paste(interim, paste0(
            "\nCI outside specified range of [",
            sprintf("%.4f, %.4f", fClow, 1/fClow), "];",
            "\nthe study should have been stopped for futility.\n"))
        }
      } else {
        if (BE.interim.spec1[["upper"]] < 100*fClow |
            BE.interim.spec1[["lower"]] > 100/fClow) {
          interim <- paste(interim, paste0(
            "\nCI outside specified range of [",
            sprintf("%.4f, %.4f", fClow, 1/fClow), "];",
            "\nthe study should have been stopped for futility.\n"))
        }
      }
    }
  }

  #########################################################
  # Optional second part: based on the adjusted alpha(s). #
  #########################################################
  if (!skip | TIE.interim.est > alpha0) {
    if (type == 1) {
      pwr.interim.calc <- power.TOST(alpha=adj[1], theta0=GMR, CV=CV1,
                                     n=n1, method=pmethod, design=design)
      interim.adj <- paste0("\n\nInterim analysis (adjusted \u03B11 ",
                            round(adj[1], 5), ")")
      interim.adj <- paste0(interim.adj, paste0("\n",
                                                paste0(rep("\u2500", 52), collapse="")))
      # In Type 1 BE is based on alpha1
      BE.interim.adj <- round(100*CI.BE(alpha=adj[1], pe=PE1, CV=CV1,
                                        n=n1, design=design), 2)
      if (BE.interim.adj[["lower"]] >= 80 &
          BE.interim.adj[["upper"]] <= 125) {
        BE.interim.adj.assess <- pass
      } else {
        BE.interim.adj.assess <- fail
      }
      interim.adj <- paste(interim.adj,
                           sprintf("%s%.2f%%%s", "\n", 100*(1-2*adj[1]), " CI:"),
                           sprintf("%.2f%s%.2f%%", BE.interim.adj[["lower"]], "\u2013",
                                   BE.interim.adj[["upper"]]),
                           BE.interim.adj.assess,
                           "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
      if (pwr.interim.calc >= target &
          BE.interim.adj.assess == fail) {
        interim.adj <- paste(interim.adj,
                             "\nStudy should have been stopped (\u2265 target power) and",
                             "\nconclusions based on the interim.\n")
      } else { # Lower than target power and not BE. Proceed to stage 2.
        ifelse (design != "parallel", n1.tmp <- n1, n1.tmp <- sum(n1))
        if (!usePE) { # Use fixed GMR.
          n2.adj <- sampleN2.TOST(alpha=adj[2], CV=CV1, n1=sum(n1), theta0=GMR,
                                  theta1=theta1, theta2=theta2,
                                  targetpower=target, design=design,
                                  method=pmethod)[["Sample size"]]
        } else {      # Use the PE of stage 1.
          n2.adj <- sampleN2.TOST(alpha=adj[2], CV=CV1, n1=sum(n1), theta0=PE1,
                                  theta1=theta1, theta2=theta2,
                                  targetpower=target, design=design,
                                  method=pmethod)[["Sample size"]]
        }
        interim.adj <- paste(interim.adj, sprintf("%s %i %s%i%s",
                                                  "\nSecond stage with", n2.adj,
                                                  "subjects (N=", n1 + n2.adj,
                                                  ") is justified.\n"))
      }
    } else { # Type 2
      pwr.interim.calc <- power.TOST(alpha=alpha0, theta0=GMR, CV=CV1, n=n1,
                                     method=pmethod, design=design)
      interim.adj <- "\n\nInterim analysis (adjusted \u03B1"
      if (pwr.interim.calc >= target) {
        interim.adj <- paste0(interim.adj, "0 ", alpha0, ")")
      } else {
        interim.adj <- paste0(interim.adj, "1 ", round(adj[1], 5), ")")
      }
      interim.adj <- paste0(interim.adj, paste0("\n",
                                                paste0(rep("\u2500", 52), collapse="")))
      # If at least the target power, assess BE with alpha0.
      if (pwr.interim.calc >= target) {
        BE.interim.adj0 <- round(100*CI.BE(alpha=alpha0, pe=PE1, CV=CV1,
                                           n=n1, design=design), 2)
        if (BE.interim.adj0[["lower"]] >= 80 &
            BE.interim.adj0[["upper"]] <= 125) {
          BE.interim.adj.assess <- pass
        } else {
          BE.interim.adj.assess <- fail
        }
        interim.adj <- paste(interim.adj,
                             sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha0), " CI:"),
                             sprintf("%.2f%s%.2f%%", BE.interim.adj0[["lower"]], "\u2013",
                                     BE.interim.adj0[["upper"]]),
                             BE.interim.adj.assess,
                             "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
        # If less than the target power, assess BE with adjusted alpha1.
      } else {
        BE.interim.adj1 <- round(100*CI.BE(alpha=adj[1], pe=PE1, CV=CV1,
                                           n=n1, design=design), 2)
        if (BE.interim.adj1[["lower"]] >= 80 &
            BE.interim.adj1[["upper"]] <= 125) {
          BE.interim.adj.assess <- pass
        } else {
          BE.interim.adj.assess <- fail
        }
        interim.adj <- paste(interim.adj,
                             sprintf("%s%.2f%%%s", "\n", 100*(1-2*adj[1]), " CI:"),
                             sprintf("%.2f%s%.2f%%", BE.interim.adj1[["lower"]], "\u2013",
                                     BE.interim.adj1[["upper"]]),
                             BE.interim.adj.assess,
                             "\nPower    :", sprintf("%1.4f", pwr.interim.calc), pverbose)
      }
      if (pwr.interim.calc >= target | BE.interim.adj.assess == pass) {
        if (pwr.interim.calc >= target) {
          interim <- paste(interim,
                           "\nStudy should have been stopped (\u2265 target power) and",
                         "\nconclusions based on the interim.\n")
        } else {
          interim <- paste(interim,
                           "\nStudy should have been stopped (BE in the interim).\n")
        }
      } else { # Lower than target power and not BE. Proceed to stage 2.
        if (!usePE) { # Use fixed GMR.
          n2.adj <- sampleN2.TOST(alpha=adj[2], CV=CV1, n1=sum(n1), theta0=GMR,
                                  theta1=theta1, theta2=theta2,
                                  targetpower=target, design=design,
                                  method=pmethod)[["Sample size"]]
        } else {      # Use the PE of stage 1.
          n2.adj <- sampleN2.TOST(alpha=adj[2], CV=CV1, n1=sum(n1), theta0=PE1,
                                  theta1=theta1, theta2=theta2,
                                  targetpower=target, design=design,
                                  method=pmethod)[["Sample size"]]
        }
        interim.adj <- paste(interim.adj, sprintf("%s %i %s%i%s",
                                                  "\nSecond stage with", n2.adj,
                                                  "subjects (N=", n1 + n2.adj,
                                                  ") is justified.\n"))
      }
    }

    # Check whether the study should have been stopped for futility.
    if (fC.check != "none") {
      if ((KM | design == "parallel") & (n1+n2.adj > Nmax)) {
        interim.adj <- paste(interim.adj, paste0(
          "\nTotal sample size (", n1+n2.adj, ") > Nmax (", Nmax, "); the",
          "\nstudy should have been stopped for futility.\n"))
      }
      if ((fCrit == "PE" | KM) & (PE1 < fClow | PE1 > 1/fClow)) {
        interim.adj <- paste(interim.adj, paste0(
          "\nPE outside specified range of [",
          sprintf("%.4f, %.4f", fClow, 1/fClow), "];",
          "\nthe study should have been stopped for futility.\n"))
      }
      if (fCrit == "CI") {
        if (type == 1) {
          if (BE.interim.adj[["upper"]] < 100*fClow |
              BE.interim.adj[["lower"]] > 100/fClow) {
            interim.adj <- paste(interim.adj, paste0(
              "\nCI outside specified range of [",
              sprintf("%.4f, %.4f", fClow, 1/fClow), "];",
              "\nthe study should have been stopped for futility.\n"))
          }
        } else {
          if (BE.interim.adj1[["upper"]] < 100*fClow |
              BE.interim.adj1[["lower"]] > 100/fClow) {
            interim.adj <- paste(interim.adj, paste0(
              "\nCI outside specified range of [",
              sprintf("%.4f, %.4f", fClow, 1/fClow), "];",
              "\nthe study should have been stopped for futility.\n"))
          }
        }
      }
    }
  } # End of interim based on adjusted alpha(s).

  # BE in the final analysis.
  # Caution: Welch-test not implemented yet!
  # Original alpha(s).
  BE.spec <- round(100*BEpooled(alpha2, PE, CV=CV, N, design), 2)
  if (BE.spec[["lower"]] >= 80 & BE.spec[["upper"]] <= 125) {
    BE.spec.assess <- pass
  } else {
    BE.spec.assess <- fail
  }
  BE.adj.assess <- BE.spec.assess # Workaround to allow final
                                  # assessment in the output.
  if (valid) { # Only in validation: irrelevant post-hoc power
    ifelse (design == "parallel", bk <- 4, bk <- 2)
    ph.spec <- .calc.power(alpha=alpha2, ltheta1=log(theta1),
                           ltheta2=log(theta2), diffm=log(c(GMR, PE)),
                           sem=CV2se(CV)*sqrt(bk/N), df=N-3,
                           method=pmethod)
  }
  if (!skip | TIE.interim.est > alpha0) {
    # Adjusted alpha(s).
    BE.adj <- round(100*BEpooled(adj[2], PE, CV, N, design), 2)
    if (BE.adj[["lower"]] >= 80 & BE.adj[["upper"]] <= 125) {
      BE.adj.assess <- pass
    } else {
      BE.adj.assess <- fail
    }
    if (valid) { # Only in validation: irrelevant post-hoc power
      ifelse (design == "parallel", bk <- 4, bk <- 2)
      ph.adj <- .calc.power(alpha=adj[2], ltheta1=log(theta1),
                            ltheta2=log(theta2), diffm=log(c(GMR, PE)),
                            sem=CV2se(CV)*sqrt(bk/N), df=N-3,
                            method=pmethod)
    }
  }

  # Needle plots of total sample size N.
  if (pa & (TIE.interim.est > alpha0)) {
    dev.new(record=TRUE)
    op <- par(no.readonly=TRUE) # save par() options
    par(mar=c(c(3, 3, 2.5, 0))+0.1, cex.main=1.1, cex.axis=0.95,
        mgp=c(2, 0.75, 0), las=1, lend="square", tcl=-0.3)
    par(mfrow=c(1, 2), oma=c(0, 0, 2, 0 )) # split layout into 3 parts
    # top (entire width) for title
    # 50:50 below for N plots
    # 2 lines for the title
    # Axes common to both plots.
    xlim <- range(as.integer(c(names(pwr.interim.est$ntable),
                               names(pwr.interim.adj$ntable))))
    ylim <- c(0, max(c(pwr.interim.est$ntable/sum(pwr.interim.est$ntable),
                       pwr.interim.adj$ntable/sum(pwr.interim.adj$ntable))))
    # First plot: Study data.
    plot(pwr.interim.est$ntable/sum(pwr.interim.est$ntable), col="blue",
         xlim=xlim, ylim=ylim, xlab="N", lwd=3, ylab="Density",
         main="Study\u2019s approach")
    axis(2, labels=FALSE, tcl=-0.15,
         at=seq(ylim[1], ylim[2], signif(median(diff(pretty(ylim))), 3)/5))
    abline(v=c(pwr.interim.est$nmean, pwr.interim.est$nperc[2],
               range(pwr.interim.est$nperc)),
           lty=c("solid", "dashed", rep("dotted", 2)))
    legend("topright", legend=
             c(paste0("average (", round(pwr.interim.est$nmean, 1), ")"),
               paste0("median (", pwr.interim.est$nperc[[2]], ")"),
               paste0("5, 95 percentiles (",
                      pwr.interim.est$nperc[[1]], ", ",
                      pwr.interim.est$nperc[[3]], ")")),
           bg="white", box.lty=0, inset=0.005, cex=0.95, y.intersp=1.15,
           lty=c(1, 2, 3))
    # Second plot: Optimized alpha(s).
    plot(pwr.interim.adj$ntable/sum(pwr.interim.adj$ntable), col="blue",
         xlim=xlim, ylim=ylim, xlab="N", lwd=3, ylab="",
         main="Optimized \u03B1")
    axis(2, labels=FALSE, tcl=-0.15,
         at=seq(ylim[1], ylim[2], signif(median(diff(pretty(ylim))), 3)/5))
    abline(v=c(pwr.interim.adj$nmean, pwr.interim.adj$nperc[2],
               range(pwr.interim.adj$nperc)),
           lty=c("solid", "dashed", rep("dotted", 2)))
    legend("topright", legend=
             c(paste0("average (", round(pwr.interim.adj$nmean, 1), ")"),
               paste0("median (", pwr.interim.adj$nperc[[2]], ")"),
               paste0("5, 95 percentiles (",
                      pwr.interim.adj$nperc[[1]], ", ",
                      pwr.interim.adj$nperc[[3]], ")")),
           bg="white", box.lty=0, inset=0.005,
           cex=0.95, y.intersp=1.15, lty=c(1, 2, 3))
    # Title on top.
    title("Distribution of expected total sample size", cex.main=1.15,
          outer=TRUE)
    par(op) # reset options
  }

  # Collect results.
  txt <- paste(cond, "\n\nData for the interim analysis",
               paste0("\n", paste0(rep("\u2500", 41), collapse="")),
               "\nCV (MSE)           :", sprintf("%6.2f%% (%.6g)", 100*CV1,
                                                 CV2mse(CV1)),
               "\nPE (ln(T)\u2013ln(R))   :", sprintf("%6.2f%% (%.6g)", 100*PE1,
                                                      log(PE1)),
               "\nSample size        :", sprintf("%3i", sum(n1)),
               "\n\nData for the final (pooled) analysis",
               paste0("\n", paste0(rep("\u2500", 41), collapse="")),
               "\nCV (MSE)           :", sprintf("%6.2f%% (%.6g)", 100*CV,
                                                 CV2mse(CV)),
               "\nPE (ln(T)\u2013ln(R))   :", sprintf("%6.2f%% (%.6g)", 100*PE,
                                                      log(PE)),
               "\nTotal sample size  :", sprintf("%3i", sum(N)))
  txt <- paste(txt, justif, CI1.txt, interim)
  if (pa) {
    txt <- paste(txt, "\nPower based on interim data (specified \u03B1)")
    txt <- paste0(txt, paste0("\n",
                              paste0(rep("\u2500", 50), collapse="")))
    txt <- paste(txt,
                 "\nMethod             :", substr(pverbose, 2, nchar(pverbose)-1),
                 "\nStage 1            :", sprintf("%1.4f", pwr.interim.est$pBE_s1),
                 "\nExpected (based on simulations)",
                 "\nBoth stages        :", sprintf("%1.4f", pwr.interim.est$pBE),
                 "\nStudies in stage 2 :", sprintf("%.1f%%", pwr.interim.est$pct_s2),
                 "\nTotal sample size (N)",
                 "\n  Average          :", round(pwr.interim.est$nmean, 1),
                 "\n  Median           :", pwr.interim.est$nperc[[2]],
                 "\n  5, 95 percentiles:",
                 paste0(pwr.interim.est$nperc[[1]], ", ",
                        pwr.interim.est$nperc[[3]]), "\n")
  }
  txt <- paste(txt, "\nFinal analysis of pooled data (specified \u03B12",
               paste0(alpha2, ")"),
               paste0("\n", paste0(rep("\u2550", 52), collapse="")),
               sprintf("%s%.2f%%%s", "\n", 100*(1-2*alpha2), " CI:"),
               sprintf("%.2f%s%.2f%%", BE.spec[["lower"]], "\u2013",
                       BE.spec[["upper"]]), BE.spec.assess)
  if (valid) {
    txt <- paste(txt, "\nPost hoc power (irrelevant; for validation purposes)",
                 "\nBased on GMR       :", sprintf("%1.4f", ph.spec[1]),
                 "\nBased on PE        :", sprintf("%1.4f", ph.spec[2]))
  }
  if (BE.adj.assess == BE.spec.assess & TIE.interim.est <= alpha0) {
    txt <- paste(txt,
                 paste0("\n\n\u250C", paste0(rep("\u2500", 53), collapse=""), "\u2510",
                        "\n\u2502 Since no inflation of the Type I Error is expected, \u2502",
                        "\n\u2502 can accept the reported analysis.                   \u2502",
                        "\n\u2514", paste0(rep("\u2500", 53), collapse=""), "\u2518\n"))
  }
  if (TIE.interim.est > alpha0 | (TIE.interim.est <= alpha0 & !skip)) {
    txt <- paste(txt,
                 "\n\n\u03B1-optimization (objective function: TIE \u2013 0.05 \u2192 0)",
                 paste0("\n", paste0(rep("\u2500", 54), collapse="")),
                 "\nMethod             :", substr(pverbose, 2, nchar(pverbose)-1),
                 "\nConvergence        :", x$iter, "iterations",
                 paste("(run-time", signif(run.time1[3]/60, 3), "min)"),
                 "\nEstimated precision:", sprintf("%.2E", x$estim.prec))
    if (algo == 2) {
      txt <- paste(txt,
                   "\nAlgorithm 2        :", paste(fit, "fit near the estimate"),
                   "\n                    ", paste("(run-time", signif(run.time2[3]/60, 3),
                                                   "min)"))
    }
    if (type == 1 | asym) {
      txt <- paste(txt, "\nAdjusted \u03B1 1, 2    :",
                   sprintf("%1.5f, %1.5f", adj[1], adj[2]),
                   "\nAdjusted CIs       :",
                   sprintf("%.2f%%, %.2f%%",
                           100*(1-2*adj[1]), 100*(1-2*adj[2])))
    } else {
      txt <- paste(txt, "\nAdjusted \u03B1 1, 2    :",
                   paste0(sprintf("%.3f", alpha0), " | ",
                          sprintf("%.5f, %.5f", adj[1], adj[2])),
                   "\nAdjusted CIs       :",
                   paste0(sprintf("%.2f%%", 100*(1-2*alpha0)), " | ",
                          sprintf("%.2f%%, %.2f%%",
                                  100*(1-2*adj[1]), 100*(1-2*adj[2]))))
    }
    txt <- paste(txt,
                 "\n\nTIE for adjusted \u03B1 :", sprintf("%1.5f", TIE), infl, CI2.txt,
                 interim.adj)
    if (pa) {
      txt <- paste(txt, "\nPower based on interim data (adjusted \u03B1)")
      txt <- paste0(txt, paste0("\n", paste0(rep("\u2500", 50), collapse="")))
      txt <- paste(txt,
                   "\nMethod             :", substr(pverbose, 2, nchar(pverbose)-1),
                   "\nStage 1            :", sprintf("%1.4f", pwr.interim.adj$pBE_s1),
                   "\nExpected (based on simulations)",
                   "\nBoth stages        :", sprintf("%1.4f", pwr.interim.adj$pBE),
                   "\nStudies in stage 2 :", sprintf("%.1f%%", pwr.interim.adj$pct_s2),
                   "\nTotal sample size (N)",
                   "\n  Average          :", round(pwr.interim.adj$nmean, 1),
                   "\n  Median           :", pwr.interim.adj$nperc[[2]],
                   "\n  5, 95 percentiles:",
                   paste0(pwr.interim.adj$nperc[[1]], ", ",
                          pwr.interim.adj$nperc[[3]]), "\n")
    }
    txt <- paste(txt, "\nFinal analysis of pooled data (adjusted \u03B12",
                 paste0(round(adj[2], 5), ")"),
                 paste0("\n", paste0(rep("\u2550", 52), collapse="")),
                 sprintf("%s%.2f%%%s", "\n", 100*(1-2*adj[2]), " CI:"),
                 sprintf("%.2f%s%.2f%%", BE.adj[["lower"]], "\u2013",
                         BE.adj[["upper"]]), BE.adj.assess)
    if (valid) {
      txt <- paste(txt, "\nPost hoc power (irrelevant; for validation purposes)",
                   "\nBased on GMR       :", sprintf("%1.4f", ph.adj[1]),
                   "\nBased on PE        :", sprintf("%1.4f", ph.adj[2]))
    }
    if (BE.adj.assess == BE.spec.assess) {
      txt <- paste(txt,
                   paste0("\n\n\u250C", paste0(rep("\u2500", 43), collapse=""), "\u2510",
                          "\n\u2502 Since conclusions of both analyses agree, \u2502",
                          "\n\u2502 can accept the original analysis.         \u2502",
                          "\n\u2514", paste0(rep("\u2500", 43), collapse=""), "\u2518\n"))
    } else {
      if (BE.spec.assess == pass & BE.adj.assess == fail) {
        risk <- sprintf("%.1f%%.", 100*(TIE.interim.est-alpha0)/alpha0)
        txt <- paste(txt,
                     paste0("\n\n\u250C", paste0(rep("\u2500", 40+nchar(risk)), collapse=""), "\u2510",
                            "\n\u2502 Accepting the reported analysis could in-",
                            paste0(rep(" ", nchar(risk)-2), collapse=""), "\u2502",
                            "\n\u2502 crease the relative consumer risk by ~", risk, " \u2502",
                            "\n\u2514", paste0(rep("\u2500", 40+nchar(risk)), collapse=""), "\u2518\n"))
      } else {
        txt <- paste(txt,
                     paste0("\n\n\u250C", paste0(rep("\u2500", 55), collapse=""), "\u2510",
                            "\n\u2502 The study could have demonstrated BE with adjusted \u03B1. \u2502",
                            "\n\u2502 However, this is not of a regulatory concern.         \u2502",
                            "\n\u2514", paste0(rep("\u2500", 55), collapse=""), "\u2518\n"))
      }
    }
  }
  txt <- paste(txt, "\nProgram offered for Use without any Guarantees and Absolutely",
               "\nNo Warranty. No Liability is accepted for any Loss and Risk",
               "\nto Public Health Resulting from Use of this Code.")

  if (plot.it & (TIE.interim.est > alpha0)) {
    ylim <- c(0.95*min(80, BE.spec, BE.adj, 125),
              max(80, BE.spec, BE.adj, 125)/0.95)
    xlim <- c(0, 3)
    plot(1, 1, type="n", log="y", axes=FALSE, xlim=xlim, ylim=ylim,
         xlab="adjusted \u03B1", ylab="% of reference", cex.axis=0.95)
    axis(1, labels=sprintf("%.5f", c(alpha2, adj[2])), at=1:2,
         tick=FALSE)
    axis(2, las=1)
    axis(2, labels=FALSE, tcl=-0.15,
         at=seq(ylim[1], ylim[2], signif(median(diff(pretty(ylim))), 3)/5))
    axis(3, labels=sprintf("%.2f%%", 100*(1-2*c(alpha2, adj[2]))),
         at=1:2, tick=FALSE)
    mtext("adjusted confidence interval", 3, line=3)
    axis(4, labels=FALSE)
    axis(4, labels=FALSE, tcl=-0.15,
         at=seq(ylim[1], ylim[2], signif(median(diff(pretty(ylim))), 3)/5))
    abline(h=100*c(theta1, 1, theta2), col=c("red", "blue", "red"),
           lty=c(2, 1, 2))
    box()
    points(1:2, rep(100*PE, 2), pch=22, cex=1.5)
    text(1.5, 100*PE, sprintf("%.2f%%", 100*PE), cex=0.9)
    text(1, 100*PE, gsub("[()]", "", BE.spec.assess), pos=2)
    text(2, 100*PE, gsub("[()]", "", BE.adj.assess), pos=4)
    arrows(1, BE.spec[["lower"]], 1, BE.spec[["upper"]], angle=90, code=3)
    text(1, BE.spec[["lower"]], pos=1, cex=0.9,
         labels=sprintf("%.2f%%", BE.spec[["lower"]]))
    text(1, BE.spec[["upper"]], pos=3, cex=0.9,
         labels=sprintf("%.2f%%", BE.spec[["upper"]]))
    arrows(2, BE.adj[["lower"]], 2, BE.adj[["upper"]], angle=90, code=3)
    text(2, BE.adj[["lower"]], pos=1, cex=0.9,
         labels=sprintf("%.2f%%", BE.adj[["lower"]]))
    text(2, BE.adj[["upper"]], pos=3, cex=0.9,
         labels=sprintf("%.2f%%", BE.adj[["upper"]]))
  }
  # insert exec.end into txt
  exec.end <- strftime(Sys.time(), usetz=TRUE) # Timestamp.
  if (valid) {
    cut.pos <- unlist(gregexpr(pattern="\nValidation", txt))
  } else {
    cut.pos <- unlist(gregexpr(pattern="\nStudy conditions ", txt))
  }
  left.str   <- substr(txt, 1, cut.pos-1)
  right.str  <- substr(txt, cut.pos, nchar(txt))
  insert.str <- paste0("         completed: ", exec.end, "\n")
  txt <- paste0(left.str, insert.str, hr, "\n", right.str, "\n\n")

  # Results to console.
  cat(txt, "\n")
}
