#' Highest density interval for the beta distribution
#'
#' @param c First parameter of the beta distribution.
#' @param d Second parameter of the beta distribution.
#' @param rho A number in (0, 1) representing the fixed posterior probability of the HPD interval.
#' @param len A positive real number representing the fixed length for the HPD interval.
#'
#' @return Lower and upper bounds of the HD interval.
#' @export
#'
#' @examples .
hd.beta <- function(c, d, rho = NULL, len = NULL, N = 1E3) {
  if (is.null(len)) {
    out <- pscl::betaHPD(c, d, p = rho)
  }
  if (is.null(rho)) {
    theta <- sort(rbeta(N, c, d))
    cov <- numeric()
    for (i in 1:N) {
      cov[i] <- sum(theta > max(theta[i], 0) & theta < min(theta[i] + len, 1))/N
    }
    ind <- which.max(cov)
    out <- c(max(theta[ind], 0), min(theta[ind] + len, 1))
  }
  return(out)  
}
