\name{RcmdrPlugin.EBM-package}
\alias{RcmdrPlugin.EBM-package}
\alias{RcmdrPlugin.EBM}
\docType{package}
\title{
Evidence Based Medicine Rcmdr Plug-In
}
\description{
This package provides an Rcmdr GUI extension plug-in for Evidence Based Medicine calculations. The purpose of this GUI extension was to improve teaching in the use of Evidence Based Medicine calculations. The interface allows for computing diagnostic medical indicators: sensitivity, specificity, diagnostic accuracy, Youden's index, positive predictive value, negative predictive value, positive likelyhood ratio, negative likelyhood ratio, number needed to diagnose; prognosis medical indicators: relative risk, atributable risk, odds ratio; therapeutic medical indicators: absolute risk reduction, number needed to treat, relative risk reduction, relative risk, odds ratio. Also the interface allows computing the post-test probability from a pre-test probability and the likelyhood ratio.

 The medical indicators are computed with epi.2by2.
 For therapy indicators it offers results based AR, OR, RR with their confidence intervals. The following calculations are done to compute the therapeutic indicators:\cr
  ARR <- - AR (ARR CI lower bound =min(-AR$lower, -AR$upper), ARR CI upper bound =max(-AR$lower, - AR$upper))\cr
	NNT <- 1 / ARR (NNT CI lower bound =min(1/ARR.lower, 1/ARR.upper), NNT CI upper bound =max(1/ARR.lower, 1/ARR.upper)) \cr
		if ARR.lower<0 then NNT.lower = NNT.upper and NNT.upper = inf\cr
	RR <- RR\cr
	OR <- OR\cr
	RRR <- 1 - RR (RRR CI lower bound =min(1-RR$lower, 1-RR$upper), RRR CI upper bound =max(1-RR$lower, 1-RR$upper))\cr
}
\details{
\tabular{ll}{
Package: \tab RcmdrPlugin.EBM\cr
Type: \tab Package\cr
Version: \tab 1.0-2\cr
Date: \tab 2012-04-13\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
}

\author{
Daniel-Corneliu Leucuta

Maintainer: Daniel-Corneliu Leucuta <danny.ldc@gmail.com>
}

\keyword{ package }

\seealso{
\code{\link[Rcmdr]{Rcmdr}}
}