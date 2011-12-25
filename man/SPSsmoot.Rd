\name{SPSsmoot}
\alias{SPSsmoot}
\title{Spline Smoothing of Outcome by Treatment in Supervised Propensiy Scoring}
\description{Express Expected Outcome by Treatment as Spline Functions of Fitted
 Propensity Scores.
}
\usage{ 
  SPSsmoot(dframe, trtm, pscr, yvar, faclev=3, df=5, spar=NULL, cv=FALSE, penalty=1)
}
\arguments{
  \item{dframe}{data.frame of the form returned by SPSlogit().}
  \item{trtm}{the two-level factor on the left-hand-side in the formula argument to SPSlogit().}
  \item{pscr}{fitted propensity scores of the form returned by SPSlogit().}
  \item{yvar}{continuous outcome measure or result unknown at the time patient was assigned
   (possibly non-randomly) to treatment; "NA"s are allowed in yvar.}
  \item{faclev}{optional; maximum number of distinct numerical values a variable can assume
   and yet still be converted into a factor variable; faclev=1 causes a binary indicator to be
   treated as a continuous variable determining a proportion.}
  \item{cv}{optional; ordinary cross-validation (T) or generalized cross-validation, GCV (F).}
  \item{df}{optional; degrees-of-freedom of B-spline fit.}
  \item{spar}{spar argument for smooth.spline() function.}
  \item{penalty}{coefficient of penalty for df in the GCV criterion.}
}
\details{Once one has fitted a somewhat smooth curve through scatters of observed outcomes, Y,
 versus the fitted propensity scores, X, for the patients in each of the two treatment groups,
 one can consider the question: "Over the range where both smooth curves are defined (i.e. their
 common support), what is the (weighted) average signed difference between these two curves?"

 If the distribution of patients (either treated or untreated) were UNIFORM over this range, the
 (unweighted) average signed difference (treated minus untreated) would be an appropriate
 estimate of the overall difference in outcome due to choice of treatment.

 Histogram patient counts within 100 cells of width 0.01 provide a naive "non-parametric density
 estimate" for the distribution of total patients (treated or untreated) along the propensity
 score axis.  The weighted average difference (and standard error) displayed by SPSsmoot() are
 based on an R density() smooth of these counts.

 In situations where the propensity scoring distribution for all patients in a therapeutic class
 is known to differ from that of the patients within the current study, that population weighted
 average would also be of interest.  Thus the SPSsmoot() output object contains two data frames,
 ssgrid and ssfit, useful in further computations.
}
\value{ An output list object of class SPSsmoot:
  \item{ssgrid}{spline grid data.frame containing 11 variables and 100 observations. The PS
   variable contains propensity score "cell means" of 0.005 to 0.995 in steps of 0.010.
   Variables F0, S0 and C0 for treatment 0 and variables F1, S1 and C1 for treatment 1 contain
   fitted smooth spline values, standard error estimates and patient counts, respectively.  The
   DIF variable is simply (F1\-F0), the SED variable is sqrt(S1\^2+S0\^2), the HST variable is
   proportional to (C0+C1), and the DEN variable is the estimated probability density of patients
   along the PS axis.}
  \item{spsub0, spsub1}{spline fit data.frames containing 4 variables for each distinct PS value in ssfit.
   These 4 variables are named PS, YAVG, TRT==0 and 1, respectively, and FIT = spline prediction.}
  \item{df}{smooth.spline() degrees-of-freedom}
  \item{sptdif}{outcome treatment difference mean.}
  \item{sptsde}{outcome treatment difference standard deviation.}
}
\references{
  Chambers JM, Hastie T. (1992) \bold{Statistical Models in S} \emph{Wadsworth & Brooks/Cole}.

  Green PJ, Silverman BW. (1994) \bold{Nonparametric Regression and Generalized Linear Models:
  A Roughness Penalty Approach}. \emph{Chapman and Hall}.

  Hastie TJ, Tibshirani RJ. (1990) \bold{Generalized Additive Models}. \emph{Chapman and Hall}.

  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.

  Sheather SJ, Jones MC. (1991) A reliable data-based bandwidth selection method for kernel
  density estimation. \emph{J Roy Statist Soc B} \bold{53}: 683-690.

  R implementation of smooth.spline() by Ripley BD and Maechler M. ('spar/lambda', etc).
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{SPSloess}}, \code{\link{SPSbalan}} and \code{\link{SPSoutco}}.}
\examples{
  data(lindner)
  PStreat <- abcix~stent+height+female+diabetic+acutemi+ejecfrac+ves1proc
  logtSPS <- SPSlogit(lindner, PStreat, PSfit, PSrnk, PSbin, appn="lindSPS")

  SPScbss7 <- SPSsmoot(lindSPS, abcix, PSfit, cardbill, df=7)
  SPScbss7
  plot(SPScbss7)
}
\keyword{smooth}