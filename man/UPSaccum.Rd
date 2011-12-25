\name{UPSaccum}
\alias{UPSaccum}
\title{Prepare for Accumulation of (Outcome,Treatment) Results in Unsupervised Propensity
 Scoring.
}
\description{Specify key result accumulation parameters: Treatment t-Factor, Outcome
  Y-variable, faclev setting, scedasticity assumption, and name of the UPSgraph() data
  accumulation object.
}
\usage{UPSaccum(hiclus, dframe, trtm, yvar, faclev=3, scedas="hete", accobj="UPSframe")}
\arguments{
  \item{hiclus}{Name of UPShclus() output object created using the diana, agnes or hclust
     method.}
  \item{dframe}{Name of data.frame containing the X, t & Y variables.}
  \item{trtm}{Name of treatment factor variable.}
  \item{yvar}{Name of outcome Y variable.}
  \item{faclev}{Maximum number of different numerical values an outcome variable can assume
   without automatically being converted into a "factor" variable; faclev=1 causes a binary
   indicator to be treated as a continuous variable determining an average or proportion.}
  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
  \item{accobj}{Name of the object for accumulation of I-plots to be ultimately displayed
     using UPSgraph().}
}
\details{The second phase in an Unsupervised Propensity Scoring analysis is to prepare to
  accumulate results over a wide range of values for "Number of Clusters."  As the number of
  such clusters increases, individual clusters will tend to become smaller and smaller and,
  thus, more and more compact in covariate X-space. 
}
\value{The output object will automaticlly be named UPSaccum.pars:
  \item{hiclus}{Name of a diana, agnes or hclust object created by UPShclus().}
  \item{dframe}{Name of data.frame containing the X, t & Y variables.}
  \item{trtm}{Name of treatment factor variable.}
  \item{yvar}{Name of outcome Y variable.}
  \item{faclev}{Maximum number of different numerical values an outcome variable can assume
   without automatically being converted into a "factor" variable; faclev=1 causes a binary
   indicator to be treated as a continuous variable determining a proportion.}
  \item{scedas}{Scedasticity assumption: "homo" or "hete"}
  \item{accobj}{Name of the object for accumulation of I-plots to be ultimately displayed
     using UPSgraph().}
  \item{nnymax}{Maximum NN LTD Standard Error observed; Upper NN plot limit;
     initialized to zero.}
  \item{nnxmin}{Minimum NN LTD observed; Left NN plot limit; initialized to zero.}
  \item{nnxmax}{Maximum NN LTD observed; Right NN plot limit; initialized to zero.}
}
\references{
  Obenchain RL. (2004) Unsupervised Propensity Scoring: NN and IV Plots.
  \emph{Proceedings of the American Statistical Association (on CD)} 8 pages.

  Obenchain RL. (2011) \bold{USPSinR.pdf}  USPS R-package vignette, 40 pages.
}
\author{Bob Obenchain <wizbob@att.net>}
\seealso{\code{\link{UPSnnltd}}, \code{\link{UPSivadj}} and \code{\link{UPShclus}}.}
\examples{
  data(lindner)
  UPSxvars <- c("stent", "height", "female", "diabetic", "acutemi", "ejecfrac", "ves1proc")
  UPSharch <- UPShclus(lindner, UPSxvars)

  UPSaccum(UPSharch, lindner, abcix, lifepres, faclev=1, scedas="homo", accobj="ABClife")
  UPSaccum.pars
    
  lif001nn <- UPSnnltd(1)
  lif020nn <- UPSnnltd(20)
  lif070nn <- UPSnnltd(70)
  UPSgraph()

  ABClife
}
\keyword{univar}
\keyword{design}