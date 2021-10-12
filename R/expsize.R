#' Sample size determination for survival data using exponential assumption
#'
#' @description Sample size determination for control drug and test drug for time to event
#' outcome using exponential assumption
#'
#' @usage expsize(type, k, delta, lambda1, lambda2, sigma1, sigma2,
#'         sigma.lambda, alpha, beta)
#'
#' @details
#' Our aim is to determine the sample size based on the hazard rates for median survival times between
#' control drug and test drug. Since, the hazard function is constant for an exponential distribution, the
#' median survival time is determined by the hazard function. Moreover, comparing the hazard rates between
#' the treatment drugs is our hypothesis of interest.
#'
#' @param type           There are three different types of comparison tests: (1) test for equality,
#'                       (2) test for non-inferiority/superiority, (3) test for equivalence,
#'                       ie. type = c("equal", "noninf.sup", "equiv")
#'
#' @param k              Ratio of sample sizes
#' @param delta          The superiority or non-inferiority margin
#' @param lambda1        Hazard rate of the control drug
#' @param lambda2        Hazard rate of the test drug
#' @param sigma1         Variability in the hazard rate due to using control drug
#' @param sigma2         Variability in the hazard rate due to using test drug
#' @param sigma.lambda   Variability in the hazard rate due to combination of control and test drug
#' @param alpha          Level of significance
#' @param beta           The probability of type-II error
#'
#' @return
#' expsize returns a sample size for control and test drug intervention.
#'
#'
#' @examples
#'
#' # (a) Test for equality:
#'
#' # The exponential assumption is used to determine the  sample size with  null hypothesis
#' # that the hazard rates of a test drug and a reference drug are equal i.e.type ="equal".
#' # The both sample sizes are taken to be equal (k = 1). The  hazard rate of  control drug
#' # is lambda1 = 2 and that of test drug is  lambda2 = 1. The standard deviation (s.d.) in
#' # hazard rate due to using control drug & test drug is 0.97 and 3.94 respectively. Their
#' # combined standard deviation is sigma.lambda = 2.56. The level of significance is alpha
#' # = 0.05 and the probability of type-II error is beta = 0.20.
#'
#'
#' expsize(type = "equal", k = 1, delta = 0, lambda1 = 2, lambda2 = 1, sigma1 = 0.97,
#'         sigma2 = 3.94, sigma.lambda = 2.56, alpha = 0.05, beta = 0.20)
#'
#'
#' # (b) Test for noninferiority/superiority:
#'
#' # The exponential assumption is used to determine sample size by testing null hypothesis
#' # (type = "noninf.sup") that the difference  between the hazard rates of a test drug and
#' # the reference drug is less than or equal to a superiority margin delta = 0.2,where k=1
#' # indicates both the sample  sizes are taken to be equal. The hazard rate of the control
#' # drug is lambda1 = 2 and  that of test drug is  lambda2 = 1. The standard  deviation in
#' # hazard rate due to using control drug & test drug is 0.97 and 3.94 respectively. Their
#' # combined standard deviation is sigma, lambda =2.56. The level of significance is alpha
#' # = 0.05 and the probability of type-II error is beta = 0.20.
#'
#'
#' expsize(type = "noninf.sup", k = 1, delta = 0.2, lambda1 = 2, lambda2 = 1, sigma1 = 0.97,
#'         sigma2 = 3.94, sigma.lambda = 2.56, alpha = 0.05, beta = 0.20)
#'
#'
#'  # (c) Test for equivalence:
#'
#' # The exponential assumption is used to determine sample size by testing null hypothesis
#' # (type = "equiv") that the absolute difference between the hazard rates of a  test drug
#' # and a ref drug is greater than or equal to a superiority margin delta =0.5, where k =1
#' # indicates both the sample  sizes are taken to be equal. The hazard rate of the control
#' # drug is lambda1 = 2 and  that of test drug is  lambda2 = 1.  The standard deviation in
#' # the hazard rate  due to using control drug and test drug is 0.97 and 3.94 respectively.
#' # Their combined standard deviation is sigma.lambda = 2.56. The level of significance is
#' #  alpha = 0.05 and the probability of type-IIqerror is beta = 0.20.
#'
#'
#' expsize(type = "equiv", k = 1, delta = 0.5, lambda1 = 2, lambda2 = 1, sigma1 = 0.97,
#'         sigma2 = 3.94, sigma.lambda = 2.56, alpha = 0.05, beta = 0.20)
#'
#'
#'
#'
#'
#'
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.match crt.unmatch phsize precsize prsize crsize
#' @export
#'
expsize <- function(type, k, delta, lambda1, lambda2, sigma1, sigma2,
                    sigma.lambda, alpha, beta){
  epsilon = lambda1 - lambda2
  if(type == "equal")
  {
    z.alpha = qnorm(1-alpha/2)
    z.beta = qnorm(1-beta)
    n2 = ((z.alpha + z.beta)^2 / epsilon^2)*((sigma1^2/k) + sigma2^2)
    n1 = k*n2
  }else if(type == "noninf.sup"){
    z.alpha = qnorm(1-alpha)
    z.beta = qnorm(1-beta)
    n2 = ((z.alpha + z.beta)^2 / (epsilon - delta)^2)*((sigma1^2/k) + sigma2^2)
    n1 = k*n2
  }else if(type == "equiv"){
    z.alpha = qnorm(1-alpha)
    z.beta = qnorm(1 - beta/2)
    n2 = ((z.alpha + z.beta)^2 / (delta - abs(epsilon))^2)*((sigma1^2/k) + sigma2^2)
    n1 = k*n2
  }else{
    z.alpha = qnorm(1 - alpha/2)
    z.beta = qnorm(1-beta)
    n2 = (((z.alpha * sigma.lambda^2 * (1 + 1/k)) +
             (z.beta * sqrt((sigma1^2/k) + sigma2^2)))^2 / epsilon^2)
    n1 = k*n2
  }
  return(cat("Sample size required for control drug is",ceiling(n1),
             "and that of test drug is",ceiling(n2),"."))
}

