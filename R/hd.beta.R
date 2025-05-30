#' HD interval for the beta distribution
#'
#' Computes the highest density interval for the beta distribution by optimization or simulation.
#'
#' @param c First parameter of the beta distribution.
#' @param d Second parameter of the beta distribution.
#' @param rho A number in (0, 1) representing the fixed probability for the HD interval.
#' @param len A positive real number representing the fixed length for the HD interval.
#' @param N Number of replicates used in the simulation when the length is fixed. Default is 1000.
#'
#' @return Lower and upper bounds of the HD interval.
#' 
#' @note For fixed probability (\code{rho}) of the HD interval the function uses \code{betaHPD} 
#' function from \code{pscl} R package. For fixed length of the HD interval the function uses an algorithm proposed in M'Lan et al. (2006).
#'
#' @references 
#' Costa, E. G. (2025). Bayesian Sample Size for Binomial Proportions with Applications in R. In: Awe, O.O., A. Vance, E. (eds) Practical Statistical Learning and Data Science Methods. STEAM-H: Science, Technology, Engineering, Agriculture, Mathematics & Health. Springer, Cham. \doi{10.1007/978-3-031-72215-8_14}.
#' 
#' M’Lan, C.E., Joseph, L., Wolfson, D.B. (2006). Bayesian sample size determination for case-control studies. Journal of the American Statistical Association, 101, 760–772.
#' 
#' @export
#' 
#' @importFrom pscl betaHPD
#'
#' @examples
#' # HD interval with probability 0.95 for Beta(2, 3)
#' hd.beta(c = 2, d = 3, rho = 0.95)
#' 
#' # HD interval with length 0.3 for Beta(2, 3)
#' hd.beta(c = 2, d = 3, len = 0.3)
hd.beta <- function(c, d, rho = NULL, len = NULL, N = 1E3) {
  if (is.null(len)) {
    out <- pscl::betaHPD(c, d, p = rho)
  }
  if (is.null(rho)) {
    theta <- sort(stats::rbeta(N, c, d))
    cov <- numeric()
    for (i in 1:N) {
      cov[i] <- sum(theta > max(theta[i], 0) & theta < min(theta[i] + len, 1))/N
    }
    ind <- which.max(cov)
    out <- c(max(theta[ind], 0), min(theta[ind] + len, 1))
  }
  return(out)  
}
