\name{fncEBMTherapy}
\alias{fncEBMTherapy}

\title{Computes ARR, NNT, RR, OR, RRR and their confidence intervals}
\description{
 This dialog box permits you to compute several medical indicators ARR, OR, RR, RRR, NNT and their confidence intervals, for comparing two medical interventions. 

 The medical indicators are computed with epi.2by2 that offer calculations for AR, OR, RR with their confidence intervals. Based on these the following calculations are done to compute the therapeutic indicators:\cr
	ARR <- - AR (ARR CI lower bound =min(-AR$lower, -AR$upper), ARR CI upper bound =max(-AR$lower, - AR$upper))\cr
	NNT <- 1 / ARR (NNT CI lower bound =min(1/ARR.lower, 1/ARR.upper), NNT CI upper bound =max(1/ARR.lower, 1/ARR.upper)) \cr
		if ARR.lower<0 then NNT.lower = NNT.upper and NNT.upper = inf\cr
	RR <- RR\cr
	OR <- OR\cr
	RRR <- 1 - RR (RRR CI lower bound =min(1-RR$lower, 1-RR$upper), RRR CI upper bound =max(1-RR$lower, 1-RR$upper))\cr

}


\author{Daniel-Corneliu Leucuta <danny.ldc@gmail.com>}

\keyword{therapeutic medical indicators}
