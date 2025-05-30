% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mss.bb.R
\name{mss.bb}
\alias{mss.bb}
\title{Minimum sample size for a single binomial parameter with a beta prior distribution}
\usage{
mss.bb(
  crit,
  c,
  d,
  rho.min = NULL,
  len = NULL,
  rho = NULL,
  len.max = NULL,
  R = 1000,
  n0 = 1
)
}
\arguments{
\item{crit}{A character string specifying the criterion. Criteria: "ACC", "ALC" and "ALCApprox".}

\item{c}{First parameter of the beta prior distribution.}

\item{d}{Second parameter of the beta prior distribution.}

\item{rho.min}{A number in (0, 1) representing the minimum admissible posterior probability for the HPD interval in the ACC.}

\item{len}{A positive real number representing the length of the HPD interval in the ACC.}

\item{rho}{A number in (0, 1) representing the posterior probability of the HPD in the ALC.}

\item{len.max}{A positive real number representing the maximum admissible length for the HPD interval in the ALC.}

\item{R}{Number of replicates used in the simulation. Default is 1000.}

\item{n0}{A positive integer, \code{n0+1} is the initial sample size in which the function will check the criterion. Default is 1.}
}
\value{
An integer representing the sample size.
}
\description{
Minimum sample size for a single binomial parameter with a beta prior distribution
}
\note{
Depending on the fixed values for interval length and probability, the function may take a while to calculate the sample size. ALC tends to be faster than ACC. Since this function uses Monte Carlo simulations, the provided minimum sample sizes may vary from one call to the next. The difference is expected to decrease as the number of replicates (\code{R}) used in the Monte Carlo simulation increases. For the criterion "ALCApprox" it is used the result of Theorem 4.1 of M'Lan et al. (2008).
}
\examples{
mss.bb(crit = "ALC", c = 10, d = 2, rho = 0.9, len.max = 0.25)

mss.bb(crit = "ALCApprox", c = 10, d = 2, rho = 0.9, len.max = 0.25)

\dontrun{mss.bb(crit = "ACC", c = 2, d = 10, rho.min = 0.9, len = 0.25)}
}
\references{
Costa, E. G. (2025). Bayesian Sample Size for Binomial Proportions with Applications in R. In: Awe, O.O., A. Vance, E. (eds) Practical Statistical Learning and Data Science Methods. STEAM-H: Science, Technology, Engineering, Agriculture, Mathematics & Health. Springer, Cham. \doi{10.1007/978-3-031-72215-8_14}.

M’Lan, C.E., Joseph, L., Wolfson, D.B. (2008). Bayesian sample size determination for binomial proportions. Bayesian Analysis, 3, 269–296.
}
