#' Weighted Sum test
#'
#' @description 
#' Performs the weighted sum test described by Madsen and Browning (2009).
#'
#' @param Y a numeric vector of phenotypes. Affected individuals are coded 1 and unaffected individuals are coded 0.
#' @param X a numeric matrix of genotypes (row: individual, column: variant). Genotypes are coded 0,1 or 2 corresponding to the number of minor alleles.
#' @param Z optional numeric matrix of covariates. See Details. 
#' @param correction type of covariates adjustment. See Details.
#' @param perm number of permutations. If not NULL,a "standard permutation procedure" is performed to compute the significance. See Details. 
#' @param alpha error level. If not NULL,an "adaptive permutation procedure" is performed to compute the significance. See Details.
#' @param c precision of the p-value. If not NULL,an "adaptive permutation procedure" is performed to compute the significance. See Details.
#' @param weight.type onderation type. By default the ponderation is "Madsen-Browning". For a beta distribution  of parameters 1, 25, choose "beta".
#' @param weights numeric vector of weights. By default weights=NULL and argument weight.type will be used.
#' @param autosomal boolean. If TRUE, autosomal chromosome; FALSE, X chromosome.
#' @param gender numeric vector. 1=male; 2=female.
#'
#' @details 
#' \subsection{Permutation procedures}{
#'  Two types of permutations procedures can be defined in the function: the standard permutation procedure and the adaptive permutation procedure.
#'  
#'  In the standard permutation procedure, the user specifies, in the argument "perm", the number of permutations to be done. The p-value will be \eqn{(R+1)(B+1)}. With \eqn{R} the number of permutation statistics superior to the observed statistic and \eqn{B} the number of permutations.
#'  
#'  In the adaptive permutation procedure, the user specifies, in the argument "alpha", the significance to achieve after multiple testing correction. In the argument "c", the estimation precision of the p-value. In function of these two paremeters, the maximal number of permutations and the maximal number of success to achieve will be computed. If the maximal number of success is reached, the p-value will be \eqn{R/B}. If not, the maximal number of permutations will be used to compute the p-value \eqn{(R+1)(B+1)}.
#' }
#'
#' \subsection{Adjustment for covariates}{
#'  The adjustment for covariates is based on the permutation procedure described by Epstein et al. (2012).
#'  
#'  This package incorporates code from the package R BiasedUrn for the Multivariate Fisher NonCentral Hypergeometric distribution. The modification resides in a different installation settings that enables the analysis of thousands of individuals,
#' }
#' 
#' @return 
#' \item{p.value}{the p-value obtained by the phenotype permutation procedure. }
#' \item{stat}{the test statistic. }
#' 
#' @references 
#' TODO : REF of CAST test
#'
#' Che R, Jack JR, Motsinger-Reif AA, Brown CC (2014) An adaptive permutation approach for genome-wide association study: evaluation and recommendations for use. BioData Min 7:9 . doi: 10.1186/1756-0381-7-9
#'
#' Epstein MP, Duncan R, Broadaway KA, et al (2012) Stratification Score Matching Improves Correction for Confounding by Population Stratification in Case-Control Association Studies. Genet Epidemiol 36:195-205 . doi: 10.1002/gepi.21611
#'
#' Agner Fog (2015). BiasedUrn: Biased Urn Model Distributions. R package version 1.07. https://CRAN.R-project.org/package=BiasedUrn
#'
#' 
#' @export
#' @useDynLib DoEstRare wSum_stat
#'
#' @examples
#' X=matrix(sample(c(0,1,2), prob=c(0.999,0.001,0),replace=TRUE, 50000), nrow=500)
#' m=apply(X, MARGIN=2, 'sum')
#' X=X[,-which(m==0), drop=FALSE]
#' Y=sample(c(0,1), replace=TRUE, 500)
#' Z=NULL
#'
#' alpha=0.05
#' c=0.2
#' 
#' res=wSum.test(X=X, Y=Y, alpha=alpha, c=c); res
#' 
wSum.test=function(Y, X, Z=NULL, correction="Price", perm=100, alpha=NULL, c=NULL, weight.type="Madsen-Browning", weights=NULL, autosomal=TRUE, gender=rep(0, length(Y))){
  
  N=nrow(X)
  Ncase=sum(Y)
  P=ncol(X)
  
  Madsen=0
  wbeta=0
  Madsentot=0
  
  
  # weights to use
  if(is.null(weights)){
    weights=rep(0,P)   #to change
    if(weight.type=="Madsen-Browning"){
      Madsen=1
    }else if(weight.type=="beta"){
      wbeta=1
    }else if(weight.type=="Madsen-Browning_tot"){
      Madsentot=1
    }
  }
  
  if(is.null(Z)){
    Ymu=rep(mean(Y), N)
    converged=T
  }else{
    res.glm=glm(Y~Z, family="binomial")
    Ymu=res.glm$fitted.values
    converged=res.glm$converged
    
    if(correction=="Epstein"){
      d.odds <- exp (res.glm$linear.predictors)
      m1 <- c(rep(1, length(Y)))
    }
  }
  
  if(converged==T & sum(abs(Y-Ymu))!=0){
    stat=.C("wSum_stat",
            as.numeric(c(X)),
            as.numeric(Y),
            as.numeric(Ymu),
            as.numeric(sum(Y)),
            as.numeric(N-sum(Y)),
            as.numeric(P),
            as.numeric(autosomal),
            as.numeric(gender),
            as.numeric(Madsen),
            as.numeric(wbeta),
            as.numeric(Madsentot),
            as.numeric(weights),
            numeric(1)
    )[[13]]
    
    #adaptive permutations
    if(!is.null(alpha) & !is.null(c)){
      b=choose_b(alpha, c)
      r=choose_r(alpha, c)
      
      j=1
      Rij=0
      while(Rij<r & j<b){
        
        if(is.null(Z)){
          Yperm=sample(Y)
          converged=T
        }else if (correction=="Price"){
          Yperm=sample(Y)
          res.glm=glm(Yperm~Z, family="binomial")
          converged=res.glm$converged
          Ymu=res.glm$fitted.values
        }else if (correction=="Epstein"){
          Yperm=rMFNCHypergeo(1, m1, Ncase, d.odds)
        }
        
        if(converged==T & sum(abs(Yperm-Ymu))!=0){
          
          stat_perm=.C("wSum_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       as.numeric(autosomal),
                       as.numeric(gender),
                       as.numeric(Madsen),
                       as.numeric(wbeta),
                       as.numeric(Madsentot),
                       as.numeric(weights),
                       numeric(1)
          )[[13]]
          
          
          if(stat_perm>=stat){
            Rij=Rij+1
          }
          
          j=j+1
        }
      }
      
      if(j<b){
        p.value=Rij/(j-1)
      }else{
        p.value=(Rij+1)/(b+1)
      }
      #standard permutation procedure
    }else{
      j=1
      Rij=0
      while(j<perm){
        
        if(is.null(Z)){
          Yperm=sample(Y)
          converged=T
        }else if (correction=="Price"){
          Yperm=sample(Y)
          res.glm=glm(Yperm~Z, family="binomial")
          converged=res.glm$converged
          Ymu=res.glm$fitted.values
        }else if (correction=="Epstein"){
          Yperm=rMFNCHypergeo(1, m1, Ncase, d.odds)
        }
        
        
        if(converged==T & sum(abs(Yperm-Ymu))!=0){
          
          stat_perm=.C("wSum_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       as.numeric(autosomal),
                       as.numeric(gender),
                       as.numeric(Madsen),
                       as.numeric(wbeta),
                       as.numeric(Madsentot),
                       as.numeric(weights),
                       numeric(1)
          )[[13]]
          
          if(stat_perm>=stat){
            Rij=Rij+1
          }
          
          j=j+1
        }
      }
      p.value=(Rij+1)/(perm+1)
    }
    
    #if no convergence of ajustment for covariables in the test statistic
  }else{
    print("algorithm did not converge")
    stat=NA
    p.value=NA
  }
  
  return(list(stat=stat, p.value=p.value))
}