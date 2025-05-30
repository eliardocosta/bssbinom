#' HD interval for the difference of two beta distributions
#' 
#' Computes the highest density interval for the difference of two beta distributions by simulation.
#'
#' @param c1 First parameter of the beta distribution of the first proportion.
#' @param d1 Second parameter of the beta distribution of the first proportion.
#' @param c2 First parameter of the beta distribution of the second proportion.
#' @param d2 Second parameter of the beta distribution of the second proportion.
#' @param rho A number in (0, 1) representing the fixed probability for the HD interval.
#' @param len A positive real number representing the fixed length for the HD interval.
#' @param N Number of replicates used in the simulation. Default is 1000.
#'
#' @return Lower and upper bounds of the HD interval.
#' 
#' @note For fixed probability (\code{rho}) of the HD interval the function uses \code{emp.hpd} function from \code{TeachingDemos} R package. For fixed length of the HD interval the function uses an algorithm proposed in M'Lan et al. (2006).
#' 
#' @references 
#' Guimarães da Costa, E. (2025). Bayesian Sample Size for Binomial Proportions with Applications in R. In: Awe, O.O., A. Vance, E. (eds) Practical Statistical Learning and Data Science Methods. STEAM-H: Science, Technology, Engineering, Agriculture, Mathematics & Health. Springer, Cham. \doi{10.1007/978-3-031-72215-8_14}.
#' 
#' M’Lan, C.E., Joseph, L., Wolfson, D.B. (2006). Bayesian sample size determination for case-control studies. Journal of the American Statistical Association, 101, 760–772.
#' 
#' @export
#' 
#' @importFrom TeachingDemos emp.hpd
#'
#' @examples
#' hd.diffbetas(c1 = 8, d1 = 2, c2 = 2, d2 = 8, rho = 0.95)
#' 
#' hd.diffbetas(c1 = 8, d1 = 2, c2 = 2, d2 = 8, len = 0.2)
hd.diffbetas <- function(c1, d1, c2, d2, rho = NULL, len = NULL, N = 1E3) {
  theta1 <- stats::rbeta(N, c1, d1)
  theta2 <- stats::rbeta(N, c2, d2)
  thetad <- sort(theta1 - theta2)
  if (is.null(len)) {
    out <- TeachingDemos::emp.hpd(thetad, conf = rho)
  }
  if (is.null(rho)) {
    cov <- numeric()
    for (i in 1:N) {
      cov[i] <- sum(thetad > thetad[i] & thetad < min(thetad[i] + len, 1))/N
    }
    ind <- which.max(cov)
    out <- c(thetad[ind], min(thetad[ind] + len, 1))
  }
  return(out)  
}
