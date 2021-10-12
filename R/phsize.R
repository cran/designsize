#' Sample size determination using proportional hazard assumption
#'
#' @description Determination of sample sizes of the control drug and the test drug intervention for
#' time to event outcome using proportional hazard assumption.
#'
#' @usage phsize(type, lambda1, lambda2, delta, prop, d, alpha, beta)
#'
#' @details
#' The proportional hazards assumption is used for comparing time to event data where, we assume that
#' the hazard function is the product of two components. One component is the non-parametric part which
#' is generally called baseline hazard and another one is the parametric part. The covariates of this
#' regression model are included in the later component. Because of the combination of parametric and
#' non-parametric components, the model is known as semi-parametric model.
#'
#'
#' @param type     There are three different types of comparison tests: (1) test for equality,
#'                 (2) test for non-inferiority/superiority, (3) test for equivalence,
#'                 ie. type = c("equal", "noninf.sup", "equiv")
#' @param lambda1  Hazard rate of the control group
#' @param lambda2  Hazard rate of the test group
#' @param delta    The inferiority or superiority margin
#' @param prop     The proportion of patients in the control group
#' @param d        The probability of observing an event
#' @param alpha    Level of significance
#' @param beta     The probability of type-II error
#'
#' @return
#' phsize returns a sample size for the control and the test drug intervention.
#'
#' @examples
#' # (a) Test for equality:
#'
#' # The phsize function determines the sample size using  proportional hazards assumption.
#' # The type = "equal" denotes the two survival curves are equal under the null hypothesis
#' # and the hazard rate of the control group is lambda1 = 1 and that of the  test group is
#' # lambda2 = 2. The proportion of patients in the control group is prop = 0.5. The proba-
#' # bility of observing an event is 0.8. The level of significance is alpha = 0.05 and the
#' # probability of type-II error is beta = 0.20.
#'
#'
#' phsize(type = "equal", lambda1 = 1, lambda2 = 2, delta = 0, prop = 0.5,
#'        d = 0.8, alpha = 0.05, beta = 0.20)
#'
#'
#' # (b) Test for non-inferiority/superiority:
#'
#' # The phsize function determines the sample size using  proportional hazards assumption.
#' # The type = "noninf.sup" denotes the difference of the two survival curves is less than
#' # or  equal to the marginal value delta = 0.3. The  hazard rate of  the control group is
#' # lambda1 = 1 and  that of  the test group is lambda2 = 2. The proportion of patients in
#' # the control group is prop = 0.5, the probability of observing a event is 0.8 and level
#' # of significance is alpha = 0.05 and the probability of type-II error is beta = 0.20.
#'
#'
#' phsize(type = "noninf.sup", lambda1 = 1, lambda2 = 2, delta = 0.3, prop = 0.5,
#'        d = 0.8, alpha = 0.05, beta = 0.20)
#'
#'
#' # (c) Test for equivalence:
#'
#' # The phsize function determines the sample size using  proportional hazards assumption.
#' # The type = "equiv", denotes whether absolute value of  the differences between the two
#' # survival curves is greater than or equal to the marginal value delta = 0.5. The hazard
#' # rate of the  control group is lambda1 = 1 and that of the test group is lambda2 = 2, &
#' # the proportion of patients in control group is prop = 0.5 and the  probability of obs-
#' # erving a event ais 0.8. The  level of significance is alpha = 0.05 and the probability
#' # of type-II error is beta = 0.20.
#'
#'
#' phsize(type = "equiv", lambda1 = 1, lambda2 = 1, delta = 0.5, prop = 0.5,
#'        d = 0.8, alpha = 0.05, beta = 0.20)
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.match crt.unmatch precsize prsize crsize
#' @export
#'
phsize <- function(type, lambda1, lambda2, delta, prop, d, alpha, beta)
{
  if(lambda1 > 0 && lambda2 > 0 && max(prop,d) <=1 && min(prop,d) > 0)
  {
    b = log(lambda2/lambda1) #logarithm of the ratio of hazard rates
    p2 = 1 - prop #the proportion of patients in the control group
    if(type == "equal")
    {
      if(b != 0)
      {
        z.alpha = qnorm(1-alpha/2)
        z.beta = qnorm(1-beta)
        n = ((z.alpha + z.beta)^2)/((b^2) * prop * p2 * d)
      }else{
        return(print("Logarithm of the ratio of hazard rates cannot be zero for the test of equality."))
      }
    }else if(type == "noninf.sup"){
      z.alpha = qnorm(1-alpha)
      z.beta = qnorm(1-beta)
      n = ((z.alpha + z.beta)^2)/((b-delta)^2 * prop * p2 * d)
    }else{
      z.alpha = qnorm(1-alpha)
      z.beta = qnorm(1-beta/2)
      n = ((z.alpha + z.beta)^2)/((delta-abs(b))^2 * prop * p2 * d)
    }
    return(cat("The required sample size for each treatment group is",round(n),"."))
  }else{
    return(print("Something is wrong! Please check again."))
  }
}
