#' DoEstRare
#' 
#' @description Rare variant association test comparing position density functions and mutation counts between cases and controls.
#' 
#' @param Y a numeric vector of phenotypes. Affected individuals are coded 1 and unaffected individuals are coded 0.
#' @param X a numeric matrix of genotypes (row: individual, column: variant). Genotypes are coded 0,1 or 2 corresponding to the number of minor alleles.
#' @param position a numeric vector of variant positions. 
#' @param genome.size a numeric value corresponding to the length of the analyzed region. 
#' @param Z optional numeric matrix of covariates. See Details. 
#' @param perm number of permutations. If not NULL,a "standard permutation procedure" is performed to compute the significance. See Details. 
#' @param alpha error level. If not NULL,an "adaptive permutation procedure" is performed to compute the significance. See Details.
#' @param c precision of the p-value. If not NULL,an "adaptive permutation procedure" is performed to compute the significance. See Details.
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
#' 
#' @return 
#' \item{p.value}{the p-value obtained by the phenotype permutation procedure. }
#' \item{stat}{the test statistic. }
#'
#' @references 
#' Persyn E, Karakachoff M, Le Scouarnec S, et al (2017) DoEstRare: A statistical test to identify local enrichments in rare genomic variants associated with disease. PLOS ONE 12:e0179364 . doi: 10.1371/journal.pone.0179364
#'
#' Che R, Jack JR, Motsinger-Reif AA, Brown CC (2014) An adaptive permutation approach for genome-wide association study: evaluation and recommendations for use. BioData Min 7:9 . doi: 10.1186/1756-0381-7-9
#'
#' Epstein MP, Duncan R, Broadaway KA, et al (2012) Stratification Score Matching Improves Correction for Confounding by Population Stratification in Case-Control Association Studies. Genet Epidemiol 36:195-205 . doi: 10.1002/gepi.21611
#'
#' Agner Fog (2015). BiasedUrn: Biased Urn Model Distributions. R package version 1.07. https://CRAN.R-project.org/package=BiasedUrn
#'
#'
#'
#' @export
#' @useDynLib DoEstRare DoEstRare_stat
#' @importFrom stats bw.nrd0 qnbinom binomial glm density pbinom
#' 
#'
#' @examples
#' Y=rep(c(0,1), 500)
#' X=matrix(sample(c(0,1),prob=c(0.7,0.3) ,1000*30, replace=TRUE), ncol=30)
#' position=sample(1:500,30)
#' genome.size=500
#' perm=200
#'
#' #Autosomal gene
#' #standard phenotype permutation procedure
#' DoEstRare(Y, X, position, genome.size, perm=perm)
#' #adaptive phenotype permutation procedure
#' DoEstRare(Y, X, position, genome.size, alpha=0.05, c=0.2)
#'
#' #X gene
#' gender=rep(c(1,2), each=500)
#' #standard phenotype permutation procedure
#' DoEstRare(Y, X, position, genome.size, perm=perm, autosomal=FALSE, gender=gender)
#' #adaptive phenotype permutation procedure
#' DoEstRare(Y, X, position, genome.size, alpha=0.05, c=0.2, autosomal=FALSE, gender=gender)
#' 
DoEstRare=function(Y, X, position, genome.size,  Z=NULL, perm=NULL, alpha=NULL, c=NULL, autosomal=TRUE, gender=NULL){

  #check
  #phenotype
  if(!is.numeric(Y)){
    stop("argument \"Y\" requires a numeric vector")
  }
  tab=as.data.frame(table(Y))
  if(nrow(tab)==2){
  }else{
    stop("The phenotype is not a binary trait")
  }

  if(length(Y)!=nrow(X)){
    stop("the \"Y\" vector length is not equal to the number of rows in the \"X\" matrix")
  }

  #position
  if(!is.numeric(position)){
    stop("argument \"position argument\" requires a numeric vector")
  }

  if(length(position)!=ncol(X)){
    stop("the \"position\" vector length is not equal to the number of columns in the \"X\" matrix")
  }

  #genome.size
  if(!is.numeric(genome.size)){
    stop("argument \"genome.size\" requires a positive numeric value")
  }else{
    if(genome.size<0){
      stop("argument \"genome.size\" requires a positive numeric value")
    }
  }
  #check if genome.size is correct

  #permutations
  if(is.null(perm)){
    if(is.null(alpha)| is.null(c)){
      stop("Please use arguments \"perm\" or \"alpha\" + \"c\"  for permutations")
    }
  }

  if(length(unique(position))<=1){
    stat_obs=NA
    p.value=NA
    warning("Only one rare variant position, NA value returned")
  }else{

    P=ncol(X)   #number of variants
    N=nrow(X)   #number of individuals
    from=1
    to=genome.size
    bw=bw.nrd0(position)
    Ncase=sum(Y)

    stat_obs=.C("DoEstRare_stat",
                as.numeric(0),
                as.numeric(c(X)),
                as.numeric(Y),
                as.numeric(N),
                as.numeric(P),
                as.numeric(autosomal),
                as.numeric(gender),
                as.numeric(position),
                as.numeric(from),
                as.numeric(to) ,
                as.numeric(bw))[[1]]

    if(!is.null(Z)){
      model <- glm (Y~Z, family= binomial())
      d.odds <- exp (model$linear.predictors)
      m1 <- c(rep(1, length(Y)))
    }

    if(!is.null(perm)){
      stat_perm=rep(NA,perm)
      for(i in 1:perm){
        if(is.null(Z)){
          Y_perm=sample(Y)
        }else{
          Y_perm <- rMFNCHypergeo(1, m1, Ncase, d.odds)
        }

        stat_perm[i]=.C("DoEstRare_stat",
                        as.numeric(0),
                        as.numeric(c(X)),
                        as.numeric(Y_perm),
                        as.numeric(N),
                        as.numeric(P),
                        as.numeric(autosomal),
                        as.numeric(gender),
                        as.numeric(position),
                        as.numeric(from),
                        as.numeric(to) ,
                        as.numeric(bw))[[1]]
      }
      p.value=(length(which(stat_perm>=stat_obs))+1)/(perm+1)

    }else{
      b=choose_b(alpha, c)
      r=choose_r(alpha, c)

      j=1
      Rij=0
      while(Rij<r & j<b){
        stat_perm=NA

        if(is.null(Z)){
          Y_perm=sample(Y)
        }else{
          Y_perm <- rMFNCHypergeo(1, m1, Ncase, d.odds)
        }

        stat_perm=.C("DoEstRare_stat",
                     as.numeric(0),
                     as.numeric(c(X)),
                     as.numeric(Y_perm),
                     as.numeric(N),
                     as.numeric(P),
                     as.numeric(autosomal),
                     as.numeric(gender),
                     as.numeric(position),
                     as.numeric(from),
                     as.numeric(to) ,
                     as.numeric(bw))[[1]]

        if(stat_perm>=stat_obs){
          Rij=Rij+1
        }

        j=j+1
      }

      if(j<b){
        p.value=Rij/(j-1)
      }else {
        p.value=(Rij+1)/(b+1)
      }
    }
  }


  res_DoEstRare=list(stat=stat_obs, p.value=p.value)
  return(res_DoEstRare)
}
