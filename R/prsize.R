#' Sample size determination for parallel study design.
#'
#' @description Determination of sample sizes of two factors of each of the two groups using one of the tests for
#' equality, non-inferiority/superiority or equivalence.
#'
#' @usage prsize(type, mu1, mu2, s, alpha, beta, k, r1, r2, del)
#'
#' @details
#' Parallel arm design is the most commonly used study design where subjects are randomized
#' to one or more study arms. Each study arm will be allocated a different intervention.
#' After randomization each subject will stay in their assigned arm during the whole study.
#' The randomized subjects should not inadvertently contaminate with the other group. A major
#' characteristic of a parallel study is randomization, which ensures accuracy of the results
#' and lower risk of being biased.
#'
#' @param type  	There are three types of test, (1) test of equality, (2) test for non-inferiority/superiority,
#          			  (3) test for equivalence i.e. type=c("equal","noninf.sup","equiv")
#' @param mu1   	The mean value of 1st group
#' @param mu2   	The mean value of 2nd group
#' @param s     	The common standard deviation
#' @param alpha 	The level of significance
#' @param beta  	The probability of the type II error i.e. 1 - power
#' @param k     	The ratio of 1st sample size(n1) and 2nd sample size(n2) i.e k=n1/n2
#' @param r1    	The ratio of n1fac1 (sample size of the 1st factor for 1st group) and n1 i.e r1=n1fac1/n1
#' @param r2    	The ratio of n2fac2 (sample size of the 1st factor for 2nd group) and n2 i.e r2=n2fac1/n2
#' @param del   	The superiority or non-inferiority margin
#'
#' @importFrom stats pt qt
#'
#' @return
#' prsize returns  returns the required sample sizes for each groups and their factors
#' in a 2x2 contingency table.
#'
#' @examples
#'
#' # (a) Test for equality:
#'
#' # This is a parallel study design. The type = "equal" tests the equality of mean respon-
#' # ses of a test drug (mu1 = 12) and a reference drug (mur = 8). The common standard dev-
#' # iation of the drugs is s = 5. k = 2 indicates the ratio of the sample sizes of the two
#' # groups. alpha = 0.05 is the level of significance and the probability of type-II error
#' # is beta = 0.10. The proportion  of factor- 1 and factor-2 are taken to be r1 = 0.6 and
#' # r2 = 0.6 respectively.
#'
#'
#'  prsize(type="equal", mu1=12, mu2=8, s=5, alpha=0.05, beta=0.10, k=2, r1=0.6, r2=0.6)
#'
#'
#'  # (b) Test for superiority/noninferiority:
#'
#' # This is a Parallel design. The type = "noninf.sup" test whether the difference of mean
#' # responses  of a test drug (mu1 = 12) and a reference drug (mu2 = 8) being greater than
#' # or equal to the marginal value delta = 0.8. s = 5 is the  common standard deviation of
#' # the drugs. The value k = 2 indicates the ratio of the  sample sizes of the two groups.
#' # alpha = 0.05 is the level of significance and the probability of type-II error is beta
#' # = 0.10. The proportion  of factor-1 and factor-2 are taken to be r1 = 0.6 and r2 = 0.6
#' # respectively.
#'
#'
#'  prsize(type="noninf.sup", mu1=12, mu2=8, s=5, alpha=0.05, beta=0.10, k=2, r1=0.6,
#'         r2=0.6, del=0.8)
#'
#'
#' # (c) Test for equivalence:
#'
#' # This is a Parallel design. The type = "equiv" tests  whether the absolute value of the
#' #  difference of mean responses of a test drug (mu1 = 12) and a reference drug (mu2 = 8)
#' # being less  than or equal to the  marginal value delta = 0.8. The  number of responses
#' # are m = 4 observed from each subject in each sequence.The s = 5 is the common standard
#' # deviation of the drugs. The value k = 2 indicates the ratio of the sample sizes of the
#' # two groups. The alpha = 0.05 is the level of significance and the  probability of type
#' # -II error is  beta = 0.10. The proportion  of factor-1 (r1) and factor-2 (r2) both are
#' # taken to be equal to 0.6.
#'
#'
#'  prsize(type="equiv", mu1=12, mu2=8, s=5, alpha=0.05, beta=0.10, k=2, r1=0.6,
#'         r2=0.6, del=0.8)
#'
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.match crt.unmatch phsize precsize crsize
#' @export
#'
prsize <- function(type,mu1,mu2,s,alpha,beta,k,r1,r2,del)
{ n2 = 1;
  ch = 0;
  if(type=="equal")
  {
    while(ch!=1)
    {
      n2=n2+1
      theta=(abs(mu1-mu2))/s
      lamda=(sqrt(n2)*theta)/sqrt(1+1/k)
      t.a=qt(1-alpha/2,df=((1+k)*n2-2))
      T.a=pt(t.a,df=((1+k)*n2-2),ncp=lamda)
      if(T.a<=beta){ch=1}
    }
  }
  else if(type=="noninf.sup")
  {
    while(ch!=1)
    {
      n2=n2+1
      ef=mu1-mu2;
      theta=(ef-del)/s
      lamda=(theta*sqrt(n2))/sqrt(1+(1/k))
      t.a=qt(1-alpha,df=((1+k)*n2-2))
      T.a=pt(t.a,df=((1+k)*n2-2),ncp=lamda)
      if(T.a<=beta){ch=1}
    }
  }

  else
  {
    while(ch!=1)
    {
      n2=n2+1
      ef=(abs(mu1-mu2))
      lamda=(sqrt(n2)*ef)/s*sqrt(1+1/k)
      t.a=qt(1-alpha,df=((1+k)*n2-2))
      T.a=pt(t.a,df=((1+k)*n2-2),ncp=lamda)
      if(T.a<=(beta/2)){ch=1}
    }
  }
  n2
  n1=k*n2
  n1f=round(r1*n1);
  n1m=round(n1*(1-r1));
  n2f=round(r2*n2);
  n2m=round((1-r2)*n2);
  print("Sample Size Calculation under parallel Design")
  P= matrix(c(n1f, n2f, n1f+n2f,n1m, n2m, n1m+n2m, n1, n2, n1+n2),ncol=3,byrow =T);
  colnames(P) = c("1st Group", "2nd Group", "Total");
  rownames(P) = c("factor-1", "factor-2", "Total");
  return(P)

}
