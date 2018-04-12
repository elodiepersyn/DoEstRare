#author: Elodie Persyn

#Implementation of burden rare variant association tests


#----------------------------------------------
#           Sum                           ####
#----------------------------------------------

Sum.test=function(X, Y, Z=NULL, correction="Price", alpha=NULL, c=NULL, perm=100){

  N=nrow(X)
  Ncase=sum(Y)
  P=ncol(X)

  if(is.null(Z)){
    Ymu=rep(mean(Y), N)
    converged=T
  }else{
    res.glm=glm(Y~Z, family = "binomial")
    Ymu=res.glm$fitted.values
    converged=res.glm$converged

    if(correction=="Epstein"){
      d.odds <- exp (res.glm$linear.predictors)
      m1 <- c(rep(1, length(Y)))
    }
  }

  if(converged==T & sum(abs(Y-Ymu))!=0){
    stat=.C("Sum_stat",
            as.numeric(c(X)),
            as.numeric(Y),
            as.numeric(Ymu),
            as.numeric(sum(Y)),
            as.numeric(N-sum(Y)),
            as.numeric(P),
            numeric(1)
    )[[7]]

    #adaptive permutations
    if(!is.null(alpha) & !is.null(c)){
      b=choose_b(alpha, c)
      r=choose_r(alpha, c)

      j=1
      Rij=0
      while(Rij<r & j<b){
        #sampling
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

          stat_perm=.C("Sum_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       numeric(1)
          )[[7]]


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

          stat_perm=.C("Sum_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       numeric(1)
          )[[7]]


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


#----------------------------------------------
#           CAST                           ####
#----------------------------------------------

CAST.test=function(X, Y, Z=NULL, correction="Price", alpha=NULL, c=NULL, perm=100){

  N=nrow(X)
  Ncase=sum(Y)
  P=ncol(X)

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
    stat=.C("CAST_stat",
            as.numeric(c(X)),
            as.numeric(Y),
            as.numeric(Ymu),
            as.numeric(sum(Y)),
            as.numeric(N-sum(Y)),
            as.numeric(P),
            numeric(1)
    )[[7]]

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

          stat_perm=.C("CAST_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       numeric(1)
          )[[7]]


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

          stat_perm=.C("CAST_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       numeric(1)
          )[[7]]


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


#----------------------------------------------
#           wSum                           ####
#----------------------------------------------


wSum.test=function(X, Y, Z=NULL, correction="Price", alpha=NULL, c=NULL, perm=100, weight.type="Madsen-Browning", weights=NULL, autosomal=TRUE, gender=rep(0, length(Y))){

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


#----------------------------------------------
#           aSum                           ####
#----------------------------------------------

aSum.test=function(X, Y, Z=NULL, correction="Price", alpha=NULL, c=NULL, perm=100, alpha0=0.10){

  N=nrow(X)
  Ncase=sum(Y)
  P=ncol(X)

  Madsen=0
  wbeta=0

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
    stat=.C("aSum_stat",
           as.numeric(c(X)),
           as.numeric(Y),
           as.numeric(Ymu),
           as.numeric(sum(Y)),
           as.numeric(N-sum(Y)),
           as.numeric(P),
           as.numeric(alpha0),
           as.numeric(0)
    )[[8]]

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

          stat_perm=.C("aSum_stat",
                  as.numeric(c(X)),
                  as.numeric(Yperm),
                  as.numeric(Ymu),
                  as.numeric(sum(Y)),
                  as.numeric(N-sum(Y)),
                  as.numeric(P),
                  as.numeric(alpha0),
                  as.numeric(0)
          )[[8]]


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

          stat_perm=.C("aSum_stat",
                       as.numeric(c(X)),
                       as.numeric(Yperm),
                       as.numeric(Ymu),
                       as.numeric(sum(Y)),
                       as.numeric(N-sum(Y)),
                       as.numeric(P),
                       as.numeric(alpha0),
                       as.numeric(0)
          )[[8]]

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



