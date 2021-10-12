#' Cluster number determination for cluster randomized trails (CRT) unmatched case
#'
#' @description Determination of number of clusters per group for unmatched cluster randomized trails
#'
#' @usage crt.unmatch(type, m, u1, u2, sigma1.B, sigma1.W, sigma2.B, sigma2.W,
#'             rho1, rho2, alpha, beta)
#'
#' @details
#' Instead of independent individuals, the unit of randomization is a group of subjects in a cluster
#' randomized trial(s) or group randomized trial(s). CRTs are generally more complex and the investigators
#' consider the selection of the unit of randomization and the unit of inference, matching or stratification
#' to improve treatment balance across clusters. It is also well known that CRTs need more subjects than
#' individually randomized trials to be adequately powered.
#'
#' Under unmatched case, no pair of matching is used to control the balance. The simple randomization
#' technique is generally used.
#'
#' @param type         There are three different types of comparison i.e. type = c("M", "P", "IR"),
#'                     "M" stands for comparison of means,
#'                     "P" stands for comparison of proportions and
#'                     "IR" stands for comparison of incidence rates
#' @param m            A common cluster size
#' @param u1           Mean(for M)/proportion(for P)/incidence rate(for IR) of group-1
#' @param u2           Mean(for M)/proportion(for P)/incidence rate(for IR) of group-2
#' @param sigma1.B     Between-cluster standard deviation of group-1
#' @param sigma1.W     Within-cluster standard deviation of group-1
#' @param sigma2.B     Between-cluster standard deviation of group-2
#' @param sigma2.W     Within-cluster standard deviation of group-2
#' @param rho1         Intra-cluster correlation coefficient(ICC) of group-1
#' @param rho2         Intra-cluster correlation coefficient(ICC) of group-2
#' @param alpha        Level of significance
#' @param beta         The probability of type-II error
#'
#'
#' @return
#' crt.unmatch returns a value indicating the number of clusters needed per group.
#'
#'
#' @examples
#' # (a) Comparison of means:
#'
#' # This is a cluster randomized trials for unmatched cases. The type = "M", indicates the
#' # comparison of means of two groups  taking a common cluster size m = 20. The mean value
#' # of group-1 and group-2 is u1 = 0.5 and u2 = 0.3 respectively. For group-1, the between
#' # -cluster standard deviation is sigma1.B = 0.3 and within-cluster standard deviation is
#' # sigma1.W = 0.3. Similarly, for group-2 those values are sigma2.B=0.3 and sigma2.W=0.3.
#' # The intra-cluster correlation coefficient (ICC) of group - 1 is rho1 = 0.2 and that of
#' # group - 2 is rho2 = 0.2. The level of significance is alpha = 0.05 and the probability
#' # of type-II error is beta = 0.20.
#'
#'
#' crt.unmatch(type = "M", m = 20, u1 = 0.5, u2 = 0.3, sigma1.B = 0.3, sigma1.W = 0.3,
#'             sigma2.B = 0.3, sigma2.W = 0.3, rho1 = 0.2, rho2 = 0.2,
#'             alpha = 0.05, beta = 0.20)
#'
#'
#' # (b) Comparison of proportions:
#'
#' # This is a cluster randomized trials for unmatched cases. The type = "P", indicates the
#' # comparison of proportions of two groups taking a common cluster size m = 20. The prop-
#' # ortion of group-1 andgroup-2 is u1 = 0.5 and u2 = 0.3 respectively. For group - 1, the
#' # between-cluster standard deviation is sigma1.B = 0.3 and within-cluster standard devi-
#' # ation is sigma1.W = 0.3. Similarly, for group-2 the standard deviations are sigma2.B =
#' # 0.3 and sigma2.W = 0.3. The intra-cluster correlation coefficient(ICC) of group - 1 is
#' # rho1 = 0.2 and that of group-2 is rho2 = 0.2. The level of significance is alpha =0.05
#' # and the probability of type-II error is beta = 0.20.
#'
#'
#' crt.unmatch(type = "P", m = 20, u1 = 0.5, u2 = 0.3, sigma1.B = 0.3, sigma1.W = 0.3,
#'             sigma2.B = 0.3, sigma2.W = 0.3, rho1 = 0.2, rho2 = 0.2,
#'             alpha = 0.05, beta = 0.20)
#'
#'
#' # (c) Comparison of incidence rates:
#'
#' # This is a cluster randomized trials for unmatched cases. The type = "IR" indicates the
#' # comparison of incidence rates of two groups taking a total of m = 20 person-years. The
#' # incidence rate of group-1 and group-2 is u1 = 0.5 and u2 = 0.3 respectively. For group
#' # -1, the between-cluster standard deviation is sigma1.B = 0.3 and within-cluster stand-
#' # ard deviation is sigma1.W = 0.3. Similarly, for  group-2  the  standard deviations are
#' # sigma2.B = 0.3 and sigma2.W = 0.3. The  intra-cluster correlation coefficient (ICC) of
#' # group-1 is rho1 = 0.2 and  that of group-2 is rho2 = 0.2. The level of significance is
#' # alpha = 0.05 and the probability of type-II error is beta = 0.20.
#'
#'
#' crt.unmatch(type = "IR", m = 20, u1 = 0.5, u2 = 0.3, sigma1.B = 0.3, sigma1.W = 0.3,
#'             sigma2.B = 0.3, sigma2.W = 0.3, rho1 = 0.2, rho2 = 0.2,
#'             alpha = 0.05, beta = 0.20)
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.match expsize phsize precsize prsize crsize
#' @export
#'
crt.unmatch <- function(type, m, u1, u2, sigma1.B, sigma1.W, sigma2.B,
                        sigma2.W,rho1, rho2, alpha, beta){
  sigma1.sq = sigma1.B^2 + sigma1.W^2
  sigma2.sq = sigma2.B^2 + sigma2.W^2
  z.alpha = qnorm(1-alpha/2)
  z.beta = qnorm(1-beta)
  delta = u1 - u2
  if(type == "M")
  {
    if(sigma1.sq > 0 && sigma2.sq > 0 && rho1 <= 1 && rho1 > -1 && rho2 <= 1 && rho2 > -1)
    {
      sigma.sq = sigma1.sq + sigma2.sq
      sigma.rho = sigma1.sq*rho1 + sigma2.sq*rho2
      c = ((z.alpha + z.beta)^2/delta^2)*((sigma.sq/m) + ((m-1)/m)*sigma.rho)
    }
  }else if(type == "P"){
    if(u1 <= 1 && u1 > 0 && u2 <= 1 && u2 > 0)
    {
      if(rho1 <= 0 && rho1 > -1 && rho2 <= 0 && rho2 > -1)
      {
        r1 = -rho1
        r2 = -rho2
      }else if(rho1 <= 0 && rho1 > -1 && rho2 <= 1 && rho2 > 0){
        r1 = -rho1
        r2 = rho2
      }else if(rho1 <= 1 && rho1 > 0 && rho2 <= 0 && rho2 > -1){
        r1 = rho1
        r2 = -rho2
      }else{
        r1 = rho1
        r2 = rho2
      }
      p1.star = u1 * (1-u1) * (1 + r1*(m-1))
      p2.star = u2 * (1-u2) * (1 + r2*(m-1))
      c = ((z.alpha + z.beta)^2/(m * delta^2)) * (p1.star + p2.star)
    }
  }else if(type == "IR"){
    k1 = sqrt(sigma1.sq)/u1
    k2 = sqrt(sigma2.sq)/u2
    k.star = (k1^2 * u1^2 + k2^2 * u2^2)/(u1 + u2)
    c = (((z.alpha + z.beta)^2 * (u1 + u2))/(m * delta^2))*(1 + m*k.star)
  }else{
    return(cat("Some entries are misleading. Please check again!!"))
  }
  return(cat("Number of clusters needed per group:",round(c)))
}
