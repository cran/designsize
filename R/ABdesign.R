#' Sample size determination for A + B escalation design without dose de-escalation
#'
#' @description Determination of sample size for each dose level using A + B escalation design without dose de-escalation
#'
#' @usage ABdesign(A, B, C, D, E, prop=c())
#'
#' @details
#' Let there are "A" patients at the dose level "i" and also consider "C" and "D" as
#' predetermined value where (D >= C).
#'
#' If less than C patients have DLTs out of A patints then we escalate the dose at
#' (i+1)th level and if more than D have the DLT's out of A then we consider the
#' previous dose level (i-1)th as MTD(maximum dose level with toxicity rates occurring
#' no more than a predetermined value).
#'
#' If more than C and less than D
#' patients have DLTs then we add B more patients at the ith dose level and then if more than E (E >= D) out
#' of (A+B) patients have the DLTs then we consider the previous dose level as MTD.
#'
#' Now we are going to determine the expected number of sample size at the jth dose level.
#'
#' # prop = Vector of DLT rates at different dose level
#'
#' # n = Total number of doses
#'
#' # N = Vector of expected number of patients at different level of doses
#'
#' @param A    	Number patients at the dose level i
#' @param B    	Number of patients add at the dose level i when more than D number of patients have DLT
#' @param C    	Predetermined number patients out of A.
#' @param D    	Predetermined number patients out of A and D>=C
#' @param E    	Predetermined number patients out of A+B and E>=D
#' @param prop 	Vector of DLT rates at different dose level
#'
#' @return
#' The expected number of patients at dose levels
#'
#' @importFrom stats pbinom
#'
#'
#' @examples
#'
#' # This is A+B escalation design without dose de-escalation. Here A=3 , B=3 indicates the
#' # number of patients at the  dose level i and taking C=D=E=1 the predetermined number of
#' # patients with DLT. Prop indicates the vector of the DLT rates at different dose level.
#'
#' ABdesign(A = 3,B = 3,C = 1,D = 1,E = 1,prop= c(0.01,0.014,0.025,0.056,0.177,0.594,0.963))
#'
#' @author Atanu Bhattacharjee, Rajashree Dey ,Soutik Halder and Akash Pawar
#' @seealso crt.match crt.unmatch phsize precsize
#' @export
#'
#'
ABdesign<-function(A,B,C,D,E,prop=c())
{
  if(C<=D &&  D<=E)
  {
    n=length(prop)
    P0=rep(0,n)
    Q0=rep(0,n)
    h=rep(0,n)
    z=rep(0,C)
    P1=rep(0,n)
    Q1=rep(0,n)
    Nji=matrix(c(rep(0)),byrow=T, ncol=n, nrow=n)
    for(j in 1:n)
    {
      P0[j]=pbinom(C-1,A,prop[j])
      P1[j]=pbinom(D,A,prop[j])-pbinom((C-1),A,prop[j])
      Q0[j]=0;
      Q1[j]=0;
      for(k in C:D)
      {
        h[k]=pbinom(E-k,B,prop[j])
        Q0[j]= Q0[j]+h[k]*choose(A,k)*(prop[j])^k*(1-prop[j])^(A-k)

      }
      for(k in 1:C)
      {
        z[k]=pbinom(E-(k-1),B,prop[j])
        Q1[j]=Q1[j]+z[k]*choose(A,k-1)*(prop[j])^(k-1)*(1-prop[j])^(A-k+1)
      }

    }

    for(j in 1:n)
    {
      for( i in 1:n)
      {
        if(j<= i+1)
        {
          Nji[j,i]=(A*P0[j]+(A+B)*Q0[j])/(P0[j]+Q0[j])
        }
        else if(j==i+1)
        {
          Nji[j,i]=(A*(1-P0[j]-P1[j])+(A+B)*(P1[j]-Q0[j]))/(1-P0[j]-Q0[j])
        }
        else
        {
          Nji[j,i]=0
        }
      }

    }
    g=c()
    PP=c()
    for( i in 1:n-1)
    {
      g[i]=1
      for(j in 1:i)
      {
        g[i]=g[i]*(P0[j]+Q0[j])

      }
      PP[i]=(1-P0[i+1]-Q0[i+1])*(g[i])
    }
    PP[n]=1
    for(j in 1:n)
    {
      PP[n]=PP[n]*(P0[j]+Q0[j])
    }
    N=c()
    for( j in 1:n)
    {
      N[j]=sum(Nji[j,]*PP)
    }

    return(cat("The expected number of patients:",round(N,2)))
  }
  else
  {
    return(print("something is wrong"))
  }
}

