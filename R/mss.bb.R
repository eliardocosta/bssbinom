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
#' @param n0 A positive integer representing the initial sample size in which the function will check the criterion. Default is 1.
#'
#' @return An integer representing the sample size.
#' @export
#'
#' @examples .
mss.bb <- function(crit, c, d, rho.min = NULL, len = NULL, rho = NULL, len.max = NULL, 
                   R = 1E3, n0 = 1) {
  packageStartupMessage("Computing...")
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
        cov <- append(cov, stats::pbeta(ab[2], c + xn, d + n - xn) - stats::pbeta(ab[1], c + xn, d + n - xn))
      }
    }
  }
  if (crit == "ALCApprox") {
    zrho <- stats::qnorm((1 + rho)/2)
    n <- 4*(zrho/len.max)^2*(base::beta(c + 0.5, d + 0.5)/base::beta(c, d))^2 - c - d
  }
  if (crit == "ALC") {
    len <- len.max + 1
    n <- n0
    while (mean(len) > len.max) {
      n <- n + 1
      len <- numeric()
      for (i in 1:R) {
        theta <- stats::rbeta(n = 1, c, d)
        xn <- stats::rbinom(n = 1, size = n, prob = theta)
        ab <- hd.beta(c = c + xn, d = d + n - xn, rho = rho)
        len <- append(len, ab[2] - ab[1])
      }
    }
  }
  # Output
  cat("\nCall:\n")
  print(cl)
  cat("\nMinimum sample size:\n")
  cat("n  = ", n, "\n")
} 
