example <- function(expl = 3) { # returns expl's data
  if(expl < 1 | expl > 14) stop("\'expl\' must within 1 and 14.")
  ########################
  # defaults
  theta1  <- 0.80
  theta2  <- 1.25
  GMR     <- 0.95
  target  <- 0.80
  pmethod <- "shifted"
  type    <- 1
  alpha0  <- 0.05
  alpha1  <- alpha2 <- 0.0294
  usePE   <- pa <- plot.it <- stop1 <- KM <-  CIs <- FALSE
  int.pwr <- setseed <- skip <- Xover <- TRUE
  min.n2  <- 0
  max.n   <- Nmax <- Inf
  fCrit   <- "PE"
  fClow   <- 0
  nsims   <- 1e6
  tol     <- 1e-8
  algo    <- 1
  hr      <- paste0(rep("\u2500", 61), collapse="")
  res     <- list(Var1=NA, PE1=NA, n1=NA, Var=NA, PE=NA, N=NA, type=type,
                  usePE=usePE, GMR=GMR, alpha0=alpha0, alpha1=alpha1,
                  alpha2=alpha2, theta1=theta1, theta2=theta2, target=target,
                  pmethod=pmethod, int.pwr=int.pwr, min.n2=min.n2, max.n=max.n,
                  Nmax=Nmax, fCrit=fCrit, fClow=fClow, nsims=nsims,
                  setseed=setseed, tol=tol, pa=pa, skip=skip, algo=algo,
                  plot.it=plot.it, Xover=Xover, stop1=stop1, KM=KM, KM.des="",
                  CIs=CIs, info="")
  ########################

  if (expl == 1 | expl == 2) { # Potvin's Example 1
    res$Var1 <- c(0.020977, "MSE")
    res$PE1  <- c(0.16785, "difflog")
    res$n1   <- 12
    res$Var  <- c(0.021224, "MSE")
    res$PE   <- c(0.14401, "difflog")
    res$N    <- 14
    if (expl == 1) {
      res$info <- paste(paste0("expl = 1\n", hr,
                    "\nPotvin et al. 2008 (doi 10.1002/pst.294) target power 80%"),
                    "\nExample 1, Method B. Reported results (power for GMR 0.95):",
                    "\nInterim analyis: 94.12% CI: 104.27\u2013134.17% (power 75.6%)",
                    "\nFinal analysis : 94.12% CI: 102.83\u2013129.71% (power NR)\n")
    } else {
      res$type  <- 2
      res$stop1 <- TRUE
      res$info <- paste(paste0("expl = 2\n", hr,
                    "\nPotvin et al. 2008 (doi 10.1002/pst.294) target power 80%"),
                    "\nExample 1, Method C. Reported results (power for GMR 0.95):",
                    "\nInterim analyis: 90.00% CI: 106.26\u2013131.66% (power 84.1%)",
                    "\nFinal analysis : NA (stopped in the interim)\n")
    }
  }
  if (expl == 3 | expl == 4) { # Potvin's Example 2
    res$Var1 <- c(0.032634, "MSE")
    res$PE1  <- c(0.08396, "difflog")
    res$n1   <- 12
    res$Var  <- c(0.045896, "MSE")
    res$PE   <- c(0.014439, "difflog")
    res$N    <- 20
    if (expl == 3) {
      res$info <- paste(paste0("expl = 3\n", hr,
                    "\nPotvin et al. 2008 (doi 10.1002/pst.294) target power 80%"),
                    "\nExample 2, Method B. Reported results (power for GMR 0.95):",
                    "\nInterim analyis: 94.12% CI: 92.93\u2013127.28% (power 50.5%)",
                    "\nFinal analysis : 94.12% CI: 88.45\u2013116.38% (power 66.3%)\n")
    } else {
      res$type  <- 2
      res$info <- paste(paste0("expl = 4\n", hr,
                    "\nPotvin et al. 2008 (doi 10.1002/pst.294) target power 80%"),
                    "\nExample 2, Method C. Reported results (power for GMR 0.95):",
                    "\nInterim analyis: 94.12% CI: 92.93\u2013127.28% (power 64.9%)",
                    "\nFinal analysis : 94.12% CI: 88.45\u2013116.38% (power 66.3%)\n")
    }
  }
  if (expl == 5) { # Montague et al. D at the location of the maximum TIE.
    res$GMR     <- 0.90
    res$type    <- 2
    res$alpha1  <- res$alpha2 <- 0.0280
    res$Var1    <- c(0.20, "CV")
    res$PE1     <- c(0.92, "ratio")
    res$n1      <- 12
    res$Var     <- c(0.23315, "CV")
    res$PE      <- c(0.88, "ratio")
    res$N       <- 45 # one dropout
    res$info    <- paste(paste0("expl = 5\n", hr,
                         "\nMontague et al. 2011 (doi 10.1002/pst.483) target power 80%"),
                         "\nMethod D. GMR 0.90:",
                         "\nTable I. Maximum TIE at n1 12 and CV 20%: 0.0518\n")
  }
  if (expl == 6) { # Haybittle/Peto.
    res$pmethod <- "nct"
    res$pa     <- TRUE
    res$alpha1 <- 0.001
    res$alpha2 <- 0.049
    res$Var1   <- c(0.20, "CV")
    res$PE1    <- c(0.94, "ratio")
    res$n1     <- 12
    res$Var    <- c(0.24, "CV")
    res$PE     <- c(0.92, "ratio")
    res$N      <- 19 # one dropout
    res$info   <- paste(paste0("expl = 6\n", hr,
                        "\nHaybittle/Peto, GMR 0.95, target power 80%\n"))
  }
  if (expl == 7) { # Zheng et al.
    res$type    <- "MSDBE"
    res$pa      <- TRUE
    res$alpha1  <- 0.01
    res$alpha2  <- 0.04
    res$nsims   <- 20000
    res$CIs     <- TRUE
    res$Var1    <- c(0.30, "CV")
    res$PE1     <- c(0.95, "ratio")
    res$n1      <- 12
    res$Var     <- c(0.30, "CV")
    res$PE      <- c(0.95, "ratio")
    res$N       <- 42
    res$info    <- paste(paste0("expl = 7\n", hr,
                         "\nZheng et al. 2015 (doi 10.1002/pst.1672)"),
                         "\nMSDBE (GMR 0.95, target power 80%)",
                         "\nn1 12, CV 30%, 20,000 simulations",
                         "\nTable I  : TIE   0.046",
                         "\nTable III: power 0.775\n")
  }
  if (expl == 8 | expl == 9) { # Xu et al. high CVs
    res$fCrit   <- "CI"
    res$Var1    <- c(0.483, "CV")
    res$PE1     <- c(0.943, "ratio")
    res$n1      <- 48
    res$Var     <- c(0.429, "CV")
    res$PE      <- c(1.01, "ratio")
    res$N       <- 104
    if (expl == 8) { # Method E
      res$alpha1 <- 0.0254
      res$alpha2 <- 0.0357
      res$fClow  <- 0.9305
      res$info   <- paste(paste0("expl =  = 8\n", hr,
                          "\nXu et al. 2015 (doi 10.1002/pst.1721) target power 80%"),
                          "\nMethod E for ISCV 30\u201350%. Reported results (power for GMR 0.95):",
                          "\nInterim analyis: 94.92% CI: 78\u2013114% (power 35.8%)",
                          "\nFinal analysis : 92.86% CI: 91\u2013112% (power NR)\n")
    } else {         # Method F
      res$type   <- 2
      res$alpha1 <- 0.0259
      res$alpha2 <- 0.0349
      res$fClow  <- 0.9350
      res$info   <- paste(paste0("expl =  = 9\n", hr,
                          "\nXu et al. 2015 (doi 10.1002/pst.1721) target power 80%"),
                          "\nMethod F for ISCV 30\u201350%. Reported results (power for GMR 0.95):",
                          "\nInterim analyis: 94.82% CI: 78\u2013114% (power 48.7%)",
                          "\nFinal analysis : 93.02% CI: 91\u2013112% (power NR)\n")
    }
  }
  if (expl == 10) { # Fuglsang 2013, Method B, GMR 0.95, power 0.90
    res$target  <- 0.90
    res$alpha1  <- res$alpha2 <- 0.0284
    res$Var1    <- c(0.5, "CV")
    res$PE1     <- c(0.93, "ratio")
    res$n1      <- 60
    res$Var     <- c(0.55, "CV")
    res$PE      <- c(0.92, "ratio")
    res$N       <- 156
    res$info    <- paste(paste0("expl = 10\n", hr,
                         "\nFuglsang 2013 (doi 10.1208/s12248-013-9475-5) target power 90%"),
                         "\nMethod B. GMR 0.95:",
                         "\nTable I. Maximum TIE at n1 60 and CV 50%: 0.0501\n")
  }
  if (expl == 11) { # Fuglsang 2013, Method C/D, GMR 0.90, power=0.90
    res$GMR     <- 0.90
    res$target  <- 0.90
    res$type    <- 2
    res$alpha1  <- res$alpha2 <- 0.0269
    res$Var1    <- c(0.2, "CV")
    res$PE1     <- c(0.88, "ratio")
    res$n1      <- 12
    res$Var     <- c(0.22, "CV")
    res$PE      <- c(0.89, "ratio")
    res$N       <- 62
    res$info    <- paste(paste0("expl = 11\n", hr,
                                "\nFuglsang 2013 (doi 10.1208/s12248-013-9475-5) target power 90%"),
                                "\nMethod C/D. GMR 0.90:",
                                "\nTable I. Maximum TIE at n1 12 and CV 20%: 0.0501\n")
  }
  if (expl == 12) { # parallel res$type 2, GMR 0.95, power=0.80
    res$type    <- 2
    res$pa      <- TRUE
    res$Xover   <- FALSE
    res$Var1    <- c(0.40, "CV")
    res$PE1     <- c(0.95, "ratio")
    res$n1      <- c(60, 60)
    res$Var     <- c(0.40, "CV")
    res$PE      <- c(0.95, "ratio")
    res$N       <- 156
    res$info    <- paste(paste0("expl = 12\n", hr,
                         "\nFuglsang 2014 (doi 10.1208/s12248-014-9571-1) target power 80%"),
                         "\nParallel groups. Method C. GMR 0.95:",
                         "\nTable II. At n1 120 and CV 40%: TIE 0.0454, power 0.830\n")
  }
  if (expl == 13) { # Karalis/Macheras TSD, GMR 0.95, power=0.80
    res$pmethod <- "nct"
    res$usePE   <- TRUE
    res$Nmax    <- 150
    res$fClow   <- 0.8
    res$type    <- 2
    res$nsims   <- 100000
    res$CIs     <- TRUE
    res$KM      <- TRUE
    res$KM.des  <- "TSD"
    res$Var1    <- c(0.2, "CV")
    res$PE1     <- c(0.95,  "ratio")
    res$n1      <- 12
    res$Var     <- c(0.21672, "CV")
    res$PE      <- c(0.9195, "ratio")
    res$N       <- 24
    res$info    <- paste(paste0("expl = 13\n", hr,
                         "\nKaralis/Macheras 2013 (doi 10.1007/s11095-013-1026-3) target power 80%"),
                         "\nAdaptive TSD (modified C):",
                         "\nTable II. TIE at n1 12 and CV 20%: 0.0446 (100,000 simulations)\n")
  }
  if (expl == 14) { # one BEBAC's studies (stopped in the interim for BE)
    res$pmethod <- "exact"
    res$type    <- 2
    res$stop1   <- TRUE
    res$Var1    <- c(0.0108779, "MSE")
    res$PE1     <- c(0.05132263, "difflog")
    res$n1      <- 18
    res$info    <- paste(paste0("expl = 14\n", hr,
                         "\nBEBAC 0110804, target power 80%, BE in the interim."),
                         "\nPhoenix/WinNonlin (subject=random), PowerTOST (method=\"exact\").",
                         "\nMethod C. Reported results (power for GMR 0.95):",
                         "\nInterim analyis: 90.00% CI: 99.07\u2013111.85% (power 99.90%)",
                         "\nFinal analysis : NA (stopped in the interim)\n")
  }
  return(res)
}
