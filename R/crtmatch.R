#' Cluster numbers determination for cluster randomized trails (CRT) matched case
#'
#' @description Determine the number of clusters need per group for matched cluster randomized trails
#'
#' @usage crt.match(type, mu1, mu2, alpha, beta, sig.w, sig.bm, m, k)
#'
#' @details
#' In cluster-randomized trials (CRTs), matching is a technique that can be used to improve
#' covariate balance. Matching protects against chance imbalances in baseline covariate distributions and is
#' thought to improve study credibility. Matching is also implemented to increase study power. Pairs of similar
#' clusters are formed and then one cluster from the pair is randomized to group 1 while the other is assigned
#' to group 2. Now we are going to determine the number of cluster in each group.
#'
#' @param type    There are three types of comparison i.e. type = c("M", "P", "IR"),
#'                      "M" stands for comparison of means
#'                      "P" stands for comparison of proportions
#'                      "IR" stands for comparison of incidence rates
#' @param mu1 The mean/proportion/incidence rate value of the 1st group
#' @param mu2 The mean/proportion/incidence rate value of the 2nd group
#' @param alpha Level of significance
#' @param beta The probability of type-II error
#' @param sig.w Standard deviation of within cluster
#' @param sig.bm Standard deviation of between cluster
#' @param m Number of subject in each cluster which is person-years in the case of incidence rates
#' @param k Common value of the coefficient of variation for each group
#'
#' @return
#' crt.match returns a value indicating the number of clusters needed per group
#'
#' @importFrom stats qnorm
#'
#' @examples
#'
#' # (a) Comparison of means:
#'
#' # This is a matched cluster randomized trials. The type ="M" indicates the comparison of
#' # means. The mean responses of a test group is mu1 = 0.06 and a reference group is mur =
#' # 0.4. The standard deviationof within cluster and between  cluster are sig.w = 0.69 and
#' # sig.bm = 0.224 respectively, m = 20 indicates number of subject in each cluster. alpha
#' # =.05 is the level of significance and the probability of type-II error is beta = 0.10.
#'
#'
#' crt.match(type="M",mu1=0.6,mu2=0.4,alpha=0.05,beta=0.20,sig.w=0.69,sig.bm=0.224,m=20)
#'
#'
#' # (b) Comparison of proportions:
#'
#' # This is a matched  cluster randomized trials. Where type = "P" indicates the tests for
#' # comparison of proportions. The proportion of a test group is mu1 =0.01 and a reference
#' # group is mur = 0.0075. The Standard deviation of between cluster is sig.bm=0.224 and m
#' # = 2750 indicates number of subject in each cluster, alpha = 0.05 is the level of signi
#' # -ficance and probability of type-II error is beta = 0.10.
#'
#'
#' crt.match(type="P",mu1=0.01,mu2=0.0075,alpha=0.05,beta=0.10,sig.bm=0.0075,m=2750)
#'
#'
#' #(c) Comparison of incidence rates:
#'
#' # This is a matched cluster randomized trials. Where type = "IR" indicates the tests the
#' # comparison of incidence rates. The incidence rate of a test group is mu1 = 4.5 and for
#' # reference group is mur = 3.6. A total of m = 50 person years is considered, alpha =.05
#' # is the level ofsignificance and the probability of type-II error is beta = 0.10. k=0.3
#' # indicates the common value of the coefficient of variation for each group.
#'
#'
#' crt.match(type="IR",mu1=4.5,mu2=3.6,alpha=0.05,beta=0.20,m=50,k=0.3)
#'
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.unmatch expsize phsize precsize prsize crsize
#' @export
#'
crt.match <- function(type,mu1,mu2,alpha,beta,sig.w,sig.bm,m,k)
{
  z.a=qnorm(1-(alpha/2));
  z.b=qnorm(1-beta);
  if(type=="M")
  {
    d=mu1-mu2
    e=round((2*((sig.w^2/m)+sig.bm^2)*((z.a+z.b)^2))/d^2)

  }
  else if(type=="P")
  {
    e=round((((z.a+z.b)^2)*(mu1*(1-mu1)+mu2*(1-mu2)+2*sig.bm^2*(m-1)))/(m*(mu1-mu2)^2));

  }
  else
  {
    z.a=qnorm(1-(alpha/2));
    z.b=qnorm(1-beta);
    a=((z.a+z.b)^2)*(mu1+mu2);
    b=m*(mu1-mu2)^2;
    d=m*k^2*(mu1^2+mu2^2)/(mu1+mu2)
    e=round(((a/b)*(1+d)))

  }
  return(cat("The number of cluster in each group: ",e,"\n"))

}
crt.match(type="IR",mu1=4.5,mu2=3.6,alpha=0.05,beta=0.20,
          m=50,k=0.3)

crt.match(type="P",mu1=0.01,mu2=0.0075,alpha=0.05,beta=0.10,
          sig.bm=0.0075,m=2750)
crt.match(type="M",mu1=0.6,mu2=0.4,alpha=0.05,beta=0.20,
          sig.w=0.69,sig.bm=0.224,m=20)
