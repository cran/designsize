#' Sample size determination using power and precision analysis
#'
#' @description It determines the ratio between the sample size of power analysis and
#' precision analysis analysis and also, give the required sample sizes.
#'
#' @usage precsize(pR, pT, sigr, sigt, c, alpha, beta)
#'
#' @details
#' A pre-study power analysis for sample size determination is usually performed to calculate an appropriate
#' sample size for achieving a desired power for detecting a clinically meaningful difference at a prespecified
#' level of significance. In practice, a much larger sample size is expected for detecting a relatively smaller
#' difference, especially for clinical trials with extremely low incidence rate. As a result, sample size
#' determination based on power analysis may not be feasible. So, it is a good suggestion to determine the
#' sample size based on precision analysis.
#'
#'
#' @param pR     Incidence rate for the reference group
#' @param pT     Incidence rate for the test group
#' @param sigr   Variability in the reference group
#' @param sigt   Variability in the test group
#' @param c      Constant value for allowance of maximum error margin
#' @param alpha  Level of significance
#' @param beta   The probability of type-II error
#'
#' @return
#' precsize returns 3 values:
#' \describe{
#'   \item{R}{Ratio between the sample size of power analysis and precision analysis}
#'   \item{n.power}{Sample size required for power analysis}
#'   \item{n.precision}{Sample size required for precision analysis}
#' }
#'
#'
#' @examples
#' # The incidence rate of reference group is pR = 0.8 per thousands and that of test group
#' # is pT = 0.7 per thousands. It is also assumed that  the respective stanadard deviation
#' # of reference and test group are sigr = 2 and sigt = 1 respectively. The constant value
#' # is chosen c = 0.08 to allow the  maximum marginal error. The  level of significance is
#' # alpha = 0.05 and the probability of type-II error is beta = 0.20.
#'
#' precsize(pR = 0.8, pT = 0.7, sigr = 2, sigt = 1, c = 0.08, alpha = 0.05, beta = 0.20)
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.match crt.unmatch phsize prsize crsize
#' @export
#'
precsize <- function(pR, pT, sigr, sigt, c, alpha, beta)
{
  delta = pR - pT
  sigma.sq = sigr^2 + sigt^2
  z.alpha = qnorm(1-alpha/2)
  z.beta = qnorm(1-beta)
  n.power = ((z.alpha + z.beta)^2 * sigma.sq)/delta^2
  w = c * sqrt(sigma.sq) #maximum error margin
  R = ((1 + (z.beta/z.alpha))^2) * (w^2 / delta^2)
  n.precision = n.power/R
  return(cat(" Ratio between the sample size of power analysis and precision analysis:",R,
             "\n","Sample size required for power  analysis:",round(n.power),
             "\n","Sample size required for precision analysis:",round(n.precision)))
}
