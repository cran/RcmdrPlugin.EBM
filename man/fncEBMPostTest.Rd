\name{fncEBMPostTest}
\alias{fncEBMPostTest}

\title{Computes post test probability to have a dissease based on the pre-test probability and the likelihood ratio of the diagnostic test}
\description{
 This dialog box permits you to define the probability of a patient to have a disease
 before performing a diagnostic test (pre-test probability), and the likelihood ratio
 of the diagnostic test. Based on these two, the post-test probability to actually 
 have the disease is computed.
 
 This is computed the following way:\cr
 	pre-test odds <- pre-test probability / (1 - pre-test probability)\cr
	post-test odds <- pre-test odds * likelihood ratio\cr
	post-test probability<- post-test odds / (post-test odds + 1)\cr
}

\usage{
fncEBMPostTest(.pretest, .LR)
}

\arguments{
\item{.pretest}{Pre-test probability.}
\item{.LR}{Likekihood ratio of the diagnostic test.}
}

\value{
The post-test probability of having the dissease.
}

\examples{
fncEBMPostTest(.12, 5.7)
}

\author{Daniel-Corneliu Leucuta <danny.ldc@gmail.com>}

\keyword{post-test probability}
