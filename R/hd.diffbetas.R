#' Highest density interval for the difference of two beta distributions
#'
#' @param c1 
#' @param d1 
#' @param c2 
#' @param d2 
#' @param rho 
#' @param len 
#' @param N 
#'
#' @return Lower and upper bounds of the HD interval.
#' @export
#'
#' @examples .
hd.diffbetas <- function(c1, d1, c2, d2, rho = NULL, len = NULL, N = 1E3) {
  theta1 <- rbeta(N, c1, d1)
  theta2 <- rbeta(N, c2, d2)
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
