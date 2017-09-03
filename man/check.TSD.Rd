\encoding{utf-8}
\name{check.TSD}
\alias{check.TSD}
\title{
Checks Acceptability of Two-Stage Sequential Designs in Bioequivalence
}
\description{
  This function assesses the pre-specified alpha(s) for the \acronym{BE} decision
  in Adaptive Sequential Two-Stage Designs (\acronym{TSD}s) with sample size
  re-estimation based on simulations in order to control the consumer risk
  \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} at the nominal level.
}
\usage{
check.TSD(Var1, PE1, n1, Var, PE, N, type = 1, usePE = FALSE, GMR,
          alpha0 = 0.05, alpha1 = 0.0294, alpha2 = 0.0294, theta1, theta2,
          target = 0.80, pmethod = c("shifted", "nct", "exact"),
          int.pwr = TRUE, min.n2 = 0, max.n = Inf, Nmax = Inf,
          fCrit = c("PE", "CI"), fClow = 0, nsims = 1e6, setseed = TRUE,
          tol = 1e-8, pa = FALSE, skip = TRUE, algo = 1, plot.it = FALSE,
          valid = FALSE, expl = 3, Xover = TRUE, stop1 = FALSE,
          KM = FALSE, KM.des = c("TSD", "TSD-1", "TSD-2"), CIs = FALSE)
}
\arguments{
  \item{Var1}{
    Vector of observed variability in the first stage.\cr\cr
    If two elements are given:
    \enumerate{
      \item numeric value
      \item \verb{"CV"} coefficient of variation (use \emph{e.g.}, 0.3 for
    30\%) or \verb{"MSE"} (residual error from \acronym{ANOVA})
    }

    If four elements are given, the CV is calculated from the confidence
    interval by PowerTOST\enc{’}{'}s function \verb{CVfromCI()}:
    \enumerate{
      \item lower confidence limit (use \emph{e.g.}, 0.81 for
    81\%)
      \item upper confidence limit (use \emph{e.g.}, 1.15 for
    115\%)
      \item alpha
      \item \verb{"CI"}
    }
  }
  \item{PE1}{
    Observed point estimate in the first stage.\cr
    Must be given if a two-element vector of \verb{Var1} is given.
    \enumerate{
      \item numeric value (use \emph{e.g.}, 0.96 for
    96\%)
      \item \verb{"ratio"} or \verb{"difflog"} (difference of LSMs in log-scale)
    }
    In the case of a four-element vector of \verb{Var1}, \acronym{PE1} will
    be calculated as \eqn{\sqrt(Var1[1]*Var1[2])}.
  }
  \item{n1}{
    Sample size of the first stage of the study. In a parallel design with
    unequal group sizes may be given as a two-element vector where the first
    element is the group-size under Test and the second one under Reference.
    If \verb{n1} leads to an unbalanced design (\emph{i.e.}, is not a multiple of
    two), the code tries to keep sequences (or groups) as balanced as possible.\cr
  }
  \item{Var}{
    Vector of observed variability in the final (pooled) analysis.\cr\cr
    If two elements are given:
    \enumerate{
      \item numeric value
      \item \verb{"CV"} coefficient of variation (use \emph{e.g.}, 0.3 for
    30\%) or \verb{"MSE"} (residual error from \acronym{ANOVA})
    }

    If four elements are given, the CV is calculated from the confidence
    interval (taking the degrees of freedom of the stage-term into account):
    \enumerate{
      \item lower confidence limit (use \emph{e.g.}, 0.83 for
    83\%)
      \item upper confidence limit (use \emph{e.g.}, 1.13 for
    113\%)
      \item alpha
      \item \verb{"CI"}
    }
  }
  \item{PE}{
   Observed point estimate in the final (pooled) analysis.\cr
    Must be given if a two-element vector of \verb{Var} is given.
    \enumerate{
      \item numeric value (use \emph{e.g.}, 0.96 for
    96\%)
      \item \verb{"ratio"} or \verb{"difflog"} (difference of LSMs in log-scale)
    }
    In the case of a four-element vector of \verb{Var}, \acronym{PE} will
    be calculated as \eqn{\sqrt(Var[1]*Var[2])}.
  }
  \item{N}{
    Total sample size of the the study. If \verb{N} leads to an unbalanced
    design (\emph{i.e.}, is not a multiple of two), the code tries to keep
    sequences (or groups) as balanced as possible.
  }
  \item{type}{
    \sQuote{Type} of the Two-Stage Design (\verb{1}, \verb{2}, or \verb{"MSDBE"}).
    For the definition see \enc{Schütz}{Schuetz}.\cr
    Defaults to \verb{1} if not given explicitly.
  }
  \item{usePE}{
    If \verb{TRUE} the sample size estimation step is done with \verb{PE1}
    and the MSE of stage 1 (\emph{i.e.}, fully adaptive).\cr
    Defaults to \verb{FALSE}, (\emph{i.e.}, the sample size is estimated with
    the fixed \verb{GMR} and MSE of stage 1.\cr
    NB: The power inspection steps in the Potvin methods are always done with
    the \verb{GMR} argument and MSE (CV) of stage 1.
  }
  \item{GMR}{
    \sQuote{True} or assumed bioavailability ratio for the sample size re-estimation.\cr
    \verb{0.95}: Potvin \emph{et al.}, Zheng \emph{et al.}, Fuglsang B and C/D (1),
    Fuglsang (parallel), Karalis & Macheras, Karalis, Xu \emph{et al.}\cr
    \verb{0.90}: Montague D, Fuglsang C/D (2)
  }
  \item{alpha0}{
    Type I error (\acronym{TIE}) probability (\ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}}; nominal level of the test). Per convention
    commonly set to \verb{0.05}.
  }
  \item{alpha1}{
    Specified adjusted \ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}} of the test in the interim. Defaults to \verb{0.0294}.\cr
    Common values:\cr
    \verb{0.0010}: Haybittle/Peto\cr
    \verb{0.0050}: O\enc{’}{'}Brien/Fleming\cr
    \verb{0.0100}: Zheng \emph{et al.}\cr
    \verb{0.0248}: Xu \emph{et al.} Method F (\sQuote{Type 2}) for CV 10--30\%\cr
    \verb{0.0249}: Xu \emph{et al.} Method E (\sQuote{Type 1}) for CV 10--30\%\cr
    \verb{0.0254}: Xu \emph{et al.} Method E (\sQuote{Type 1}) for CV 30--55\%\cr
    \verb{0.0259}: Xu \emph{et al.} Method F (\sQuote{Type 2}) for CV 30--55\%\cr
    \verb{0.0269}: Fulgsang C/D (GMR 0.90, power 0.90)\cr
    \verb{0.0274}: Fuglsang C/D (GMR 0.95, power 0.90)\cr
    \verb{0.0280}: Montague D, Karalis TSD-1\cr
    \verb{0.0284}: Fuglsang B (GMR 0.95, power 0.90)\cr
    \verb{0.0294}: Potvin B, Karalis & Macheras TSD, Karalis TSD-2, Fuglsang (parallel)\cr
    \verb{0.0304}: Kieser & Rauch
  }
  \item{alpha2}{
    Specified adjusted \ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}} of the test in the final (pooled) analysis. Defaults to \verb{0.0294}.\cr
    Common values:\cr
    \verb{0.0400}: Zheng \emph{et al.}\cr
    \verb{0.0480}: O\enc{’}{'}Brien/Fleming\cr
    \verb{0.0490}: Haybittle/Peto\cr
    \verb{0.0269}: Fuglsang C/D (GMR 0.90, power 0.90)\cr
    \verb{0.0274}: Fuglsang C/D (GMR 0.95, power 0.90)\cr
    \verb{0.0280}: Montague D, Karalis TSD-1\cr
    \verb{0.0284}: Fulgsang B (GMR 0.95, power 0.90)\cr
    \verb{0.0294}: Potvin B, Karalis & Macheras TSD, Karalis TSD-2, Fuglsang (parallel)\cr
    \verb{0.0304}: Kieser & Rauch\cr
    \verb{0.0349}: Xu \emph{et al.} Method F (\sQuote{Type 2}) for CV 30--55\%\cr
    \verb{0.0357}: Xu \emph{et al.} Method E (\sQuote{Type 1}) for CV 30--55\%\cr
    \verb{0.0363}: Xu \emph{et al.} Method E (\sQuote{Type 1}) for CV 10--30\%\cr
    \verb{0.0364}: Xu \emph{et al.} Method F (\sQuote{Type 2}) for CV 10--30\%
  }
  \item{theta1}{
    Lower acceptance limit for \acronym{BE}. Defaults to \verb{0.80} if not given explicitly.
  }
  \item{theta2}{
    Upper acceptance limit for \acronym{BE}. Defaults to \verb{1.25} if not given explicitly.
  }
  \item{target}{
    Power threshold in the first step of \sQuote{Type 1} designs and power to achieve
    in the sample size re-estimation step.
  }
  \item{pmethod}{
    Power calculation method; also to be used in the sample size re-estimation
    for stage 2.\cr
    Implemented are \verb{"nct"} (approximate calculations via the non-central
    \emph{t}-distribution) and \verb{"exact"} (exact calculations via
    Owen\enc{’}{'}s Q-functions), as well as \verb{"shifted"} (approximation via
    the shifted central \emph{t}-distribution, as used in most references).\cr
    Defaults to \verb{"nct"} as a reasonable compromise between speed and accuracy
    in the power and sample size estimation step.\cr
    \verb{"shifted"} is not state of the art and should be mainly used in validation.\cr
    It is likely that in the study one of the better methods was used -- which
    are available in current software;\cr
    \verb{"nct"} is available in \emph{e.g.}, \acronym{FARTSSIE}, \acronym{PASS},
    \acronym{SAS}, and package \verb{PowerTOST}, whereas \verb{"exact"} is available
    in \acronym{NQuery}, \acronym{SAS}, and package \verb{PowerTOST}.\cr
    Caution: Adjusting \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} with
    \verb{pmethod="exact"} can take \emph{hours}.
  }
  \item{int.pwr}{
    If \verb{TRUE} (the default) the interim power monitoring step in the
    stage 1 evaluation of \sQuote{Type 1} TSDs will be done as described in all
    [\emph{sic}] references.
    Setting this argument to \verb{FALSE} will omit this step (to satisfy some
    assessors of the Dutch \acronym{MEB}).\cr
    Has no effect if \verb{type=2} is choosen.
  }
  \item{min.n2}{
    Minimum sample size of stage 2.  Defaults to 0 (as in all references).\cr
    If the sample size estimation step gives \verb{N < n1} the sample size for
    stage 2 will be set to \verb{min.n2}, \emph{i.e.}, the total sample size to
    \verb{n1+min.n2}.
  }
  \item{max.n}{
    If \verb{max.n} is set to a finite value the re-estimated total sample
    size (\verb{N}) is set to \verb{min(max.n, N)}. Defaults to \verb{Inf} which
    is equivalent to not constrain the re-estimated sample size.\cr
    Caution: \verb{max.n} here is \bold{not} a futility criterion like \verb{Nmax}.
  }
  \item{Nmax}{
    Futility criterion.  If set to a finite value all simulated studies in which
    a sample size \verb{>Nmax} is obtained will be regarded as failed in the
    interim. Mandatory if \verb{KM=TRUE} (Karalis & Macheras recommend 150).\cr
    Defaults to \verb{Inf} (unconstrained total sample size).
  }
  \item{fCrit}{
    Futility criterion to use for the point estimate or confidence interval in
    the interim. Acceptable values are \verb{"PE"} or \verb{"CI"}.\cr
    A suitable value must be given in the argument \verb{fClow}.
  }
  \item{fClow}{
    Lower futility limit for \verb{fCrit} in the interim. If the point estimate
    or confidence interval (as specified by \verb{fCrit}) is outside \verb{fClow}
    \ldots \verb{1/fClow} the study is stopped in the interim with the result FAIL
    (not BE).\cr
    May be missing. Defaults then to 0, \emph{i.e.}, no futility rule is applied.\cr
    Values if \verb{fCrit == "PE"}:\cr
    \verb{0.80}: Armitage, Karalis & Macheras, Karalis\cr
    \verb{0.85}: Charles Bon (AAPS Annual Meeting 2007)\cr
    Values if \verb{fCrit == "CI"}:\cr
    \verb{0.9250}: Detlew Labes (personal communication with Diane Potvin)\cr
    \verb{0.9374}: Xu \emph{et al.} Method E (\sQuote{Type 1}) for CV 10--30\%\cr
    \verb{0.9305}: Xu \emph{et al.} Method E (\sQuote{Type 1}) for CV 30--55\%\cr
    \verb{0.9492}: Xu \emph{et al.} Method F (\sQuote{Type 2}) for CV 10--30\%\cr
    \verb{0.9350}: Xu \emph{et al.} Method F (\sQuote{Type 2}) for CV 30--55\%
  }
  \item{nsims}{
    Number of simulations to be performed to estimate the (\sQuote{empiric})
    \acronym{TIE} error and in optimizing adjusting \ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}}. The default value 1,000,000 = 1e+6 should not be lowered.
  }
  \item{setseed}{
    Simulations are dependent on the starting point of the (pseudo)random number
    generator.To avoid differences in power for different runs, \verb{set.seed(1234567)}
    is issued if \verb{TRUE} (default).\cr
    Set to \verb{FALSE} in order to assess robustness (\emph{i.e.}, a different
    seed is issued in every call of the function).
  }
  \item{tol}{
    Desired accuracy (convergence tolerance) of \verb{uniroot()}.\cr
    The objective function is \eqn{|TIE-alpha0| \le tol}. Defaults to 1e-8.
  }
  \item{pa}{
    Should results of the power analysis for adjusted \ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}} be shown? Defaults to \verb{FALSE}.
  }
  \item{skip}{
    Should optimization of \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} be
    skipped if the \acronym{TIE} with the specified \ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}} already preserves the consumer risk (\verb{alpha0})?
    Defaults to \verb{TRUE} for speed reasons.
  }
  \item{algo}{
    Defaults to \verb{1}: Optimization by \verb{uniroot()}.\cr
    If set to \verb{2} an attempt will be made to further improve the estimate
    by assessing a set of adjusted alphas close to the first estimated
    \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}}. Easily doubles the
    computation time.
  }
  \item{plot.it}{
    Should a comparative plot of PE and CI int the final analysis be made?
    Defaults to \verb{FALSE}.\cr
    Only applicable if \ifelse{html}{\out{&alpha;}}{\eqn{\alpha}{alpha}} was
    optimized (\emph{i.e.}, if the specified \ifelse{html}{\out{&alpha;}}{
    \eqn{\alpha}{alpha}} inflates the \acronym{TIE}).
  }
  \item{valid}{
    Should one of the validation examples be assessed? Defaults to \verb{FALSE}.\cr
    If set to \verb{TRUE} one of the built-in examples must be specified in
    the argument \verb{expl}. Any of the other arguments can be additionally
    specified in order to assess their impact.
  }
  \item{expl}{
    Number of the validation examples (\verb{1} to \verb{14}).
    \enumerate{
      \item Potvin \emph{et al.} Example 1, Method B (\sQuote{Type 1}), alpha1=alpha2=0.0294, GMR=0.95, target=0.80
      \item Potvin \emph{et al.} Example 1, Method C (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0294, GMR=0.95, target=0.80, stopped in the interim
      \item Potvin \emph{et al.} Example 2, Method B (\sQuote{Type 1}), alpha1=alpha2=0.0294, GMR=0.95, target=0.80
      \item Potvin \emph{et al.} Example 2, Method C (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0294, GMR=0.95, target=0.80
      \item Montague \emph{et al.} Method D (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0280, GMR=0.90, target=0.80
      \item Haybittle/Peto (\sQuote{Type 1}), alpha1=0.001, alpha2=0.049, GMR=0.95, target=0.80
      \item Zheng \emph{et al.} (\acronym{MSDBE}), alpha1=0.01, alpha2=0.04, GMR=0.95, target=0.80
      \item Xu \emph{et al.} Method E, high CV (\sQuote{Type 1}), alpha1=0.0254, alpha2=0.0357, GMR=0.95, target=0.80, n.max=180, futility on CI 0.9305\ldots1/0.9305
      \item Xu \emph{et al.} Method F, high CV (\sQuote{Type 2}), alpha0=0.05, alpha1=0.0259, alpha2=0.0349, GMR=0.95, target=0.80, n.max=180, futility on CI 0.9350\ldots1/0.9350
      \item Fuglsang 2013, Method B (\sQuote{Type 1}), alpha1=alpha2=0.0284, GMR=0.95, target=0.90
      \item Fuglsang 2013, Method C/D (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0269, GMR=0.90, target=0.90
      \item Fuglsang 2014, parallel, Method C (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0294, GMR=0.95, target=0.80 (parallel design)
      \item Karalis/Macheras, TSD (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0294, target=0.80, futility on PE 0.80\ldots1.25, Nmax=150
      \item \acronym{BEBAC} internal, Method C (\sQuote{Type 2}), alpha0=0.05, alpha1=alpha2=0.0294, GMR=0.95, target=0.80, stopped in the interim
    }
  }
  \item{Xover}{
    Defaults to \verb{TRUE}.Set to \verb{FALSE} if the study was performed in a parallel design.
  }
  \item{stop1}{
    Defaults to \verb{FALSE} (\emph{i.e.}, the study proceeded to the second stage).\cr
    Set to \verb{TRUE} if the study stopped in the interim.  Naturally in this case arguments
    \verb{Var}, \verb{PE}, and \verb{N} are not required.
  }
  \item{KM}{
    Should a fully adaptive design according to Karalis & Macheras or Karalis be
    used? Defaults to \verb{FALSE}.\cr
    If set to \verb{TRUE} the upper sample size for futility (\verb{Nmax}) and the
    additional argument \verb{KM.des} must be given.
  }
  \item{KM.des}{
    Design-specification of Karalis & Macheras or Karalis. Acceptable values:\cr
    \verb{TSD} \sQuote{Type 2} (Karalis & Macheras).\cr
    \verb{TSD-1} \sQuote{Type 2} (Karalis).\cr
    \verb{TSD-2} \sQuote{Type 1} (Karalis).
  }
  \item{CIs}{
    Defaults to \verb{FALSE}. Set to \verb{TRUE} to show the two-sided 95\% CI
    of the empiric TIE(s).
  }
}
\details{
  The calculations follow in principle the simulations as described in the references.\cr
  The underlying subject data are assumed to be evaluated after log-transformation.
  Instead of -- time-consuming -- simulating subject data, the statistics PE1, PE, and
  MSE1, SS2 are simulated via their associated distributions (normal and
  \ifelse{html}{\out{&chi;&sup2;}}{\eqn{\chi}{chi^2}}) as suggested by Zheng \emph{et al.}
}

\note{
  The code tries to reconstruct the applied method from the study\enc{’}{'}s
  conditions. In the strict sense applying a specified \ifelse{html}{\out{&alpha;}}{
  \eqn{\alpha}{alpha}} is only justified if the study\enc{’}{'}s conditions were
  within the validated range (\emph{e.g.}, n1, CV1, GMR, target power, etc) given
  in the references.

  Hoping to get in simulations an empiric Type I Error of \emph{exactly} 0.05 (as
  -- \emph{unoffically} -- demanded by the \acronym{EMA}\enc{’}{'}s Biostatistics
  Working Party) is futile. Already rounding the confidence interval in fixed sample
  designs to two digits according to the guideline (\emph{i.e.}, 79.995\% upwards
  to 80.00\% or 125.004\% downwards to 125.00\%) may transform into a \acronym{TIE}
  of up to 0.0508.\cr
  For more details see a description in the doc subdirectory of the package.\cr
  So far the author did not came across a \acronym{TIE} after optimization of
  >0.05001. The question whether this is acceptable to the \acronym{BSWP}
  remains open.

  \bold{Program offered for Use without any Guarantees and Absolutely No Warranty.\cr
  No Liability is accepted for any Loss and Risk to Public Health Resulting from Use of this Code}
}

\author{
Helmut \enc{Schütz}{Schuetz}
}

\references{
O\enc{’}{'}Brien PC, Fleming TR. \emph{A Multiple Testing Procedure for Clinical Trials.}\cr
Biometrics. 1979;35(3):549--56.

Armitage P. \emph{Interim Analysis in Clinical Trials.}\cr
Stat Med. 1991;10(6):925--37. \doi{10.1002/sim.4780100613}

Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ, Smith RA. \emph{Sequential design approaches for bioequivalence studies with crossover designs}.\cr
Pharm Stat. 2008;7(4):245--62. \doi{10.1002/pst.294}\cr

Montague TH, Potvin D, DiLiberti CE, Hauck WW, Parr AF, Schuirmann DJ. \emph{Additional results for \sQuote{Sequential design approaches for bioequivalence studies with crossover designs}.}\cr
Pharm Stat. 2012;11(1):8--13. \doi{10.1002/pst.483}\cr

Zheng Ch, Wang J, Zhao L. \emph{Testing bioequivalence for multiple formulations with power and sample size calculations.}\cr
Pharm Stat. 2012;11(4):334--41. \doi{10.1002/pst.1522}\cr

Fuglsang A. \emph{Sequential Bioequivalence Trial Designs with Increased Power and Controlled Type I Error Rates.}\cr
AAPS J. 2013;15(3):659--61. \doi{10.1208/s12248-013-9475-5}\cr

Karalis V, Macheras P. \emph{An Insight into the Properties of a Two-Stage Design in Bioequivalence Studies.}\cr
Pharm Res. 2013;30(7):1824--35. \doi{10.1007/s11095-013-1026-3}\cr

Karalis V. \emph{The role of the upper sample size limit in two-stage bioequivalence designs.}\cr
Int J Pharm. 2013;456(1):87--94. \doi{10.1016/j.ijpharm.2013.08.013}\cr

Fuglsang A. \emph{Futility Rules in Bioequivalence Trials with Sequential Designs.}\cr
AAPS J. 2014;16(1):79--82. \doi{10.1208/s12248-013-9540-0}\cr

Fuglsang A. \emph{Sequential Bioequivalence Approaches for Parallel Designs}.\cr
AAPS J. 2014;16(3):373--8. \doi{10.1208/s12248-014-9571-1}\cr

Zheng Ch, Zhao L, Wang J. \emph{Modifications of sequential designs in bioequivalence trials.}\cr
Pharm Stat. 2015;14(3):180--8. \doi{10.1002/pst.1672}\cr

\enc{Schütz}{Schuetz} H. \emph{Two-stage designs in bioequivalence trials.}\cr
Eur J Clin Pharmacol. 2015;71(3):271--81. \doi{10.1007/s00228-015-1806-2}\cr

Kieser M, Rauch G. \emph{Two-stage designs for cross-over bioequivalence trials.}\cr
Stat Med. 2015;34(16):2403--16. \doi{10.1002/sim.6487}\cr

Xu J, Audet C, DiLiberti CE, Hauck WW, Montague TH, Parr TH, Potvin D, Schuirmann DJ. \emph{Optimal adaptive sequential designs for crossover bioequivalence studies.}\cr
Pharm Stat. 2015;15(1):15--27. \doi{10.1002/pst.1721}
}

\examples{
# Potvin et al., Example 2, 'Method B'
# Not run: Due to timing policy of CRAN for examples.
\dontrun{
check.TSD(Var1=c(0.032634, "MSE"), PE1=c(0.08396, "difflog"), n1=12,
          Var=c(0.045896, "MSE"), PE=c(0.014439, "difflog"), N=20,
          pmethod="shifted", type=1)
# or simply using a built-in dataset
check.TSD(valid=TRUE, expl=3)
# Should give TIE 0.04307 and acceptance of the reported results.
#
# As above but variabilities and PEs calculated from CIs and alpha(s).
# Varibility as a 4-element vector, and PEs not given.
# Not run: Due to timing policy of CRAN for examples.
check.TSD(Var1=c(0.9293, 1.2728, 0.0294, "CI"), n1=12,
          Var=c(0.8845, 1.1638, 0.0294, "CI"), N=20,
          pmethod="shifted", type=1)
# Slightly different TIE (0.04300) since the precision of the input
# is limited.
#
# Potvin et al., Example 2, 'Method C'
check.TSD(Var1=c(0.032634, "MSE"), PE1=c(0.08396, "difflog"), n1=12,
          Var=c(0.045896, "MSE"), PE=c(0.014439, "difflog"), N=20,
          pmethod="shifted", type=2, plot.it=TRUE)
# or using the built-in dataset
check.TSD(valid=TRUE, expl=4)
# Should give TIE 0.05062 and optimized alpha 0.02858 (TIE 0.04992).
# Original and adjusted results agree; can accept the reported results.
#
# Xu et al., 'Method E' for high CVs
check.TSD(Var1=c(0.483, "CV"), PE1=c(0.943, "ratio"), n1=48,
          Var=c(0.429, "CV"), PE=c(1.01, "ratio"), N=104,
          pmethod="shifted", type=1, fCrit="CI", fClow=0.9305,
          alpha1=0.0254, alpha2=0.0357)
# or using the built-in dataset
check.TSD(valid=TRUE, expl=8)}
# Should give TIE 0.04953 and acceptance of reported results.
}
