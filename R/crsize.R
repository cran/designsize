#' Sample size determination for crossover study design
#'
#' @description Determination of sample sizes for two factors of each group using one of the tests for
#' equality, non-inferiority/superiority or equivalence
#'
#' @usage crsize(type, delta, m, k, mur, mut, sigbr, sigbt, rho, sigwr, sigwt,
#'        alpha, beta, r1, r2)
#'
#' @details
#' Consider a 2x2m replicated crossover design for comparing mean responses of a test drug and a reference
#' drug. Under both treatments the design consists of two sequences with m subjects each.
#'
#'
#' @param type         The three different types of tests are (1) test for equality, (2) test for non-inferiority/
#'                     superiority, (3) test for equivalence i.e. type = c("equal", "noninf.sup", "equiv")
#' @param delta        Non-inferiority/Superiority margin
#' @param m            Number of responses observed from each subject in each sequence
#'                     under a fixed treatment
#' @param k            Ratio of the sample sizes of the two sequences
#' @param mur          Mean value of reference therapy
#' @param mut          Mean value of test therapy
#' @param sigbr        Between standard deviation due to the effect of reference therapy
#' @param sigbt        Between standard deviation due to the effect of test therapy
#' @param rho          Correlation between reference and test therapy
#' @param sigwr        Within standard deviation due to the effect of reference therapy
#' @param sigwt        Within standard deviation due to the effect of test therapy
#' @param alpha        Level of significance
#' @param beta         The probability of type-II error
#' @param r1           Proportion of factor-1
#' @param r2           Proportion of factor-2
#'
#' @importFrom stats qnorm
#'
#' @return
#' crsize returns the required sample sizes for each sequence and their factors
#' in a 2x2 contingency table.
#'
#' @examples
#'
#' # (a) Test for equality:
#'
#' # This is a crossover design. The type = "equal" tests the equality of mean responses of
#' # a test drug (mut = 9) and a reference drug (mur = 8.5) and the number of responses are
#' # m = 4 observed from  each subject in each  sequence. k = 1 indicates the  ratio of the
#' # sample sizes of the two sequences are equal. The between standard deviation due to the
#' # effect of reference therapy is sigbr = 1.5 and that of test therapy is 1.5. The corre-
#' # lation between reference and test  therapy is rho = 0.7. The within standard deviation
#' # due to the effect of reference therapy is sigwr = 1 as well as test therapy is sigwt =
#' # 1. The alpha = 0.05 is level of significance and the probability of type - II error is
#' # beta = 0.10. The proportion  of factor - 1 and factor - 2 are taken to be r1 = 0.5 and
#' # r2 = 0.5 respectively.
#'
#'
#' crsize(type= "equal", delta = 0.4, m = 4, k = 1, mur = 8.5, mut = 9, sigbr = 1.5,
#'        sigbt = 1.5, rho = 0.7, sigwr = 1, sigwt = 1, alpha = 0.05, beta = 0.10,
#'        r1 = 0.5, r2 = 0.5)
#'
#'
#' # (b) Test for non-inferiority/superiority:
#'
#' # This is a crossover design. The type = "noninf.sup", tests  whether the  difference of
#' # mean responses of a test drug (mut = 9) and a reference drug (mur = 8.5) being greater
#' # than or equal to  the marginal value delta = 0.4. The  number of  responses are m = 4,
#' # observed from each subject in each sequence. The value of k = 1 indicates the ratio of
#' # the sample sizes of the two sequences are equal. The between standard deviation due to
#' # the effect of  reference  therapy is sigbr = 1.5 and that of  test therapy is 1.5. The
#' # correlation between reference and test therapy is rho = 0.7. The within standard devi-
#' # ation due to  the effect of reference therapy is sigwr = 1, as well as test therapy is
#' # sigwt = 1. A alpha = 0.05 is the level of significance and  the probability of type-II
#' # error is beta = 0.10. The proportion of factor-1 (r1) and factor-2 (r2) both are taken
#' # to be 0.5.
#'
#'
#' crsize(type = "noninf.sup", delta = 0.4, m = 4, k = 1, mur = 8.5, mut = 9, sigbr = 1.5,
#'        sigbt = 1.5, rho = 0.7, sigwr = 1, sigwt = 1, alpha = 0.05, beta = 0.10,
#'        r1 = 0.5, r2 = 0.5)
#'
#'
#'
#' #(c) Test for equivalence:
#'
#' # This is a crossover design. The type = "equiv" tests whether the absolute value of the
#' # difference of mean responses of a test drug (mut = 9) and a reference drug (mur = 8.5)
#' # being  less than or equal to  the marginal value delta = 0.6. The number  of responses
#' # are m = 4 observed from each subject in each sequence. k = 1, indicates that the ratio
#' # of the sample sizes of the two sequences are equal. The between standard deviation due
#' # to the effect of reference therapy is sigbr = 1.5 and that of test therapy is 1.5. The
#' # correlation between reference and test therapy is rho = 0.7. The within standard devi-
#' # ation due to  the effect of reference therapy is sigwr = 1 as well as  test therapy is
#' # sigwt = 1. alpha = 0.05 is  the level of significance and the probability of type - II
#' # error is beta = 0.10. The proportion of factor-1 (r1) and factor-2 (r2) both are taken
#' # to be 0.5.
#'
#'
#' crsize(type = "equiv", delta = 0.6, m = 4, k = 1, mur = 8.5, mut = 9, sigbr = 1.5,
#'        sigbt = 1.5,rho = 0.7, sigwr = 1, sigwt = 1, alpha = 0.05, beta = 0.10,
#'        r1 = 0.5, r2 = 0.5)
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso ABdesign crt.match crt.unmatch phsize precsize
#' @export
#'
crsize <- function(type, delta, m, k, mur, mut, sigbr, sigbt, rho, sigwr,
                   sigwt,alpha, beta, r1, r2){
  epsilon = mut - mur  #difference between test and reference therapy
  sigmaD.sq = (sigbr)^2 + (sigbt)^2 - 2*rho*sigbr*sigbt
  #variability due to the effect of subject-by-treatment interaction
  sigmaM.sq = sigmaD.sq + (sigwr^2 + sigwt^2)/m  #random error variance
  error = 0.01
  temp = 2
  if(type == "equal")
  {
    while(error >= 0.0001 && r1+r2 == 1)
    {
      N = temp*(1+k)
      theta = sqrt(N)*(abs(epsilon)/sqrt(sigmaM.sq))
      t.star = qt(1-alpha/2, df = N-2)
      T.cdf = pt(t.star, df = N-2, ncp = theta)
      error = T.cdf - beta
      n2 = temp
      temp = n2 + 1
    }
  }else if(type == "noninf.sup"){
    if(epsilon <= delta)
    {
      temp = ((qnorm(1-alpha) + qnorm(1-beta))^2 * sigmaM.sq)/(2*(epsilon-delta)^2)
      n2 = temp
      N = n2*(1+k)
    }else{
      while(error >= 0.0001 && r1+r2 == 1)
      {
        N = temp*(1+k)
        theta = sqrt(N)*((epsilon - delta)/sqrt(sigmaM.sq))
        t.star = qt(1-alpha, df = N-2)
        T.cdf = pt(t.star, df = N-2, ncp = theta)
        error = T.cdf - beta
        n2 = temp
        temp = n2 + 1
      }
    }
  }else{
    if(abs(epsilon) <= delta)
    {
      while(error >= 0.0001 && r1+r2 == 1)
      {
        N = temp*(1+k)
        theta = sqrt(N)*((delta - abs(epsilon))/sqrt(sigmaM.sq))
        t.star = qt(1-alpha, df = N-2)
        T.cdf = pt(t.star, df = N-2, ncp = theta)
        error = T.cdf - (beta/2)
        n2 = temp
        temp = n2 + 1
      }
    }else{
      return(cat("The equivalence test is not possible for this crossover design."))
    }
  }
  n2 = round(n2)
  N = round(N) #Total sample size
  n1 = N - n2
  n1.factor1 = round(n1*r1)
  n1.factor2 = n1 - n1.factor1
  n2.factor1 = round(n2*r1)
  n2.factor2 = n2 - n2.factor1
  print("Sample size determination under Crossover Design")
  B = matrix(c(n1.factor1, n2.factor1, n1.factor1+n2.factor1,
               n1.factor2, n2.factor2, n1.factor2+n2.factor2, n1, n2, N), 3, 3)
  colnames(B) = c("Factor-1", "Factor-2", "Total")
  rownames(B) = c("Sequence-1", "Sequence-2", "Total")
  return(B)
}
