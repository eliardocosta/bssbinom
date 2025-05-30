#' Minimum sample size for a single binomial parameter with a beta prior distribution
#'
#' @param crit A character string specifying the criterion. Criteria: "ACC", "ALC" and "ALCApprox".
#' @param c First parameter of the beta prior distribution.
#' @param d Second parameter of the beta prior distribution.
#' @param rho.min A number in (0, 1) representing the minimum admissible posterior probability for the HPD interval in the ACC.
#' @param len A positive real number representing the length of the HPD interval in the ACC.
#' @param rho A number in (0, 1) representing the posterior probability of the HPD in the ALC.
#' @param len.max A positive real number representing the maximum admissible length for the HPD interval in the ALC.
#' @param R Number of replicates used in the simulation. Default is 1000.
#' @param n0 A positive integer, \code{n0+1} is the initial sample size in which the function will check the criterion. Default is 1.
#'
#' @return An integer representing the sample size.
#' 
#' @note Depending on the fixed values for interval length and probability, the function may take a while to calculate the sample size. ALC tends to be faster than ACC. Since this function uses Monte Carlo simulations, the provided minimum sample sizes may vary from one call to the next. The difference is expected to decrease as the number of replicates (\code{R}) used in the Monte Carlo simulation increases. For the criterion "ALCApprox" it is used the result of Theorem 4.1 of M'Lan et al. (2008).
#' 
#' @references
#' Guimarães da Costa, E. (2025). Bayesian Sample Size for Binomial Proportions with Applications in R. In: Awe, O.O., A. Vance, E. (eds) Practical Statistical Learning and Data Science Methods. STEAM-H: Science, Technology, Engineering, Agriculture, Mathematics & Health. Springer, Cham. <https://doi.org/10.1007/978-3-031-72215-8_14>
#' 
#' M’Lan, C.E., Joseph, L., Wolfson, D.B. (2008). Bayesian sample size determination for binomial proportions. Bayesian Analysis, 3, 269–296.
#' 
#' @export
#'
#' @examples
#' mss.bb(crit = "ALC", c = 10, d = 2, rho = 0.9, len.max = 0.2)
#' 
#' mss.bb(crit = "ALCApprox", c = 10, d = 2, rho = 0.95, len.max = 0.2)
#' 
#' mss.bb(crit = "ACC", c = 2, d = 10, rho.min = 0.9, len = 0.2)
mss.bb <- function(crit, c, d, rho.min = NULL, len = NULL, rho = NULL, len.max = NULL, R = 1E3, n0 = 1) {
  message("The computation may take a while.")
  message("Computing ...", appendLF = FALSE)
  cl <- match.call()
  if (crit == "ACC") {
    cov <- 0 
    n <- n0
    while (mean(cov) < rho.min) {
      n <- n + 1
      cov <- numeric()
      for (i in 1:R) {
        theta <- stats::rbeta(n = 1, c, d)
        xn <- stats::rbinom(n = 1, size = n, prob = theta)
        ab <- hd.beta(c = c + xn, d = d + n - xn, len = len)
        cov[i] <- stats::pbeta(ab[2], c + xn, d + n - xn) - stats::pbeta(ab[1], c + xn, d + n - xn)
      }
    }
  }
  if (crit == "ALCApprox") {
    zrho <- stats::qnorm((1 + rho)/2)
    n <- ceiling(4*(zrho/len.max)^2*(base::beta(c + 0.5, d + 0.5)/base::beta(c, d))^2 - c - d)
  }
  if (crit == "ALC") {
    len <- len.max + 1
    n <- n0
    while (mean(len) > len.max) {
      n <- n + 1
      len <- numeric(R)
      for (i in 1:R) {
        theta <- stats::rbeta(n = 1, c, d)
        xn <- stats::rbinom(n = 1, size = n, prob = theta)
        ab <- hd.beta(c = c + xn, d = d + n - xn, rho = rho)
        len[i] <- ab[2] - ab[1]
      }
    }
  }
  # Output
  message(" done!")
  cat("\nCall:\n")
  print(cl)
  cat("\nMinimum sample size:\n")
  cat("n  = ", n, "\n")
} 
