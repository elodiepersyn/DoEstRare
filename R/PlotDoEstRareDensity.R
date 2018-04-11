#' PlotDoEstRareDensity
#'
#' @description Rare variant association test comparing position density functions and mutation counts between cases and controls.
#'
#' @param pheno a numeric vector of phenotypes. Affected individuals are coded 1 and unaffected individuals are coded 0.
#' @param geno a numeric matrix of genotypes (row: individual, column: variant). Genotypes are coded 0,1 or 2 corresponding to the number of minor alleles.
#' @param position a numeric vector of variant positions.
#' @param genome.size a numeric value corresponding to the length of the analyzed region. 
#' @param autosomal boolean. If TRUE, autosomal chromosome; FALSE, X chromosome.
#' @param gender numeric vector. 1=male; 2=female.
#' @param plot.density By default, the density of rare variant positions is plotted (plot.density="density"). If you want the density function being multiplied by the weighted allele frequency, use  "burden-density".
#' @param plot.counts If TRUE, rare allele counts are indicated.
#' @param plot.legend If TRUE, the legend is added.
#' @param col.case color to use for cases.
#' @param col.ctrl color to use for controls.
#' @param args.plot a list of additional graphical parameters for the 'plot' function.
#' @param args.legend a list of additional graphical parameters for the 'legend' function.
#' @param args.text a list of additional graphical parameters for the 'text' function.
#'
#' @export
#' 
#' @importFrom graphics legend lines plot text
#'
#' @examples
#' pheno=rep(c(0,1), 500)
#' geno=matrix(sample(c(0,1),prob=c(0.99,0.01) ,1000*30, replace=TRUE), ncol=30)
#' position=sample(1:500,30)
#' genome.size=500
#'
#' # plot.density="burden-density"
#' PlotDoEstRareDensity(pheno, geno, position, genome.size, plot.legend=TRUE,
#'                     plot.density="burden-density",col.case="blue",
#'                     args.plot=list(lwd=2, xaxt="no", yaxt="no", main="Gene"),
#'                     args.legend=list(x="topleft", lwd=2),
#'                     args.text=list(cex=1))
#'
#' # plot.density="density"
#' PlotDoEstRareDensity(pheno, geno, position, genome.size, plot.legend=TRUE,
#'                     plot.density="density",col.case="blue",
#'                     args.plot=list(lwd=2, xaxt="no", yaxt="no", main="Gene"),
#'                     args.legend=list(x="topleft", lwd=2),
#'                     args.text=list(cex=1))
#' 
#' 
#' 
PlotDoEstRareDensity=function(pheno,
                              geno,
                              position,
                              genome.size,
                              autosomal=TRUE,
                              gender=NULL,
                              plot.density="density",
                              plot.counts=TRUE,
                              plot.legend=TRUE,
                              col.case="red",
                              col.ctrl="black",
                              args.plot=list(),
                              args.legend=list(),
                              args.text=list()){

  noyau="gaussian"
  lissage=1

  P=ncol(geno)   #number of variants
  N=nrow(geno)   #number of individuals
  Na=sum(pheno)  #number of cases
  Nu=N-Na        #number of controls

  #counts of rare alleles in controls
  geno_ctrl=geno[pheno==0,,drop=F]
  nombre_mutations_ctrl=apply(geno_ctrl, MARGIN=2, "sum")

  #counts of rare alleles in cases
  geno_cas=geno[pheno==1,,drop=F]
  nombre_mutations_cas=apply(geno_cas, MARGIN=2, "sum")

  #PART 1 : estimation of mutation position density functions in controls and cases
  if(sum(nombre_mutations_ctrl)==0){
    res_ctrl=rep(0,genome.size)
  } else {
    res_ctrl=density(x=position, kernel=noyau, adjust=lissage,weights=nombre_mutations_ctrl/sum(nombre_mutations_ctrl), from=1, to=genome.size, n=genome.size)$y
  }

  if(sum(nombre_mutations_cas)==0){
    res_cas=rep(0,genome.size)
  }else{
    res_cas=density(x=position, kernel=noyau, adjust=lissage,weights=nombre_mutations_cas/sum(nombre_mutations_cas), from=1, to=genome.size, n=genome.size)$y
  }

  # PLOT
  if(!is.null(args.plot$col)){
    args.plot=args.plot[-which(names(args.plot)=="col")]
    warning("Unused argument 'col'. Please use col.case and col.ctrl arguments.")
  }
  if(!is.null(args.plot$xlim)| !is.null(args.plot$ylim)){
    args.plot=args.plot[-which(names(args.plot)%in%c("xlim", "ylim"))]
    warning("Unused argument 'xlim' or 'ylim'. These parameters are fixed.")
  }

  if(plot.density=="density"){
  }else if(plot.density=="burden-density"){
    #PART 2 : computation of burden components
    if(autosomal==TRUE){
      weight=apply(cbind(nombre_mutations_cas, nombre_mutations_ctrl), MARGIN=1, function(x) pbinom(q=x[1], size=2*Na,prob=(x[2]+1)/(2*Nu+2),lower.tail=T))
      nombre_mutations_cas_weight=nombre_mutations_cas*weight
      nombre_mutations_ctrl_weight=nombre_mutations_ctrl*weight
      pa=sum(nombre_mutations_cas_weight)/(2*Na*P)/sum(weight)
      pu=sum(nombre_mutations_ctrl_weight)/(2*Nu*P)/sum(weight)
    }else{
      weight=apply(cbind(nombre_mutations_cas, nombre_mutations_ctrl), MARGIN=1, function(x) pbinom(q=x[1], size=sum(gender[pheno==1]),prob=(x[2]+1)/(sum(gender[pheno==0])+2),lower.tail=T))
      nombre_mutations_cas_weight=nombre_mutations_cas*weight
      nombre_mutations_ctrl_weight=nombre_mutations_ctrl*weight
      pa=sum(nombre_mutations_cas_weight)/(sum(gender[pheno==1])*P)/sum(weight)
      pu=sum(nombre_mutations_ctrl_weight)/(sum(gender[pheno==0])*P)/sum(weight)
    }
    res_cas=res_cas*pa
    res_ctrl=res_ctrl*pu
  }else{
    warning("The 'toplot' argument is not correct. The value by default 'toplot='density'' is used.")
  }

  if(is.null(args.plot$xlab)){
    args.plot=append(list(xlab="position"), args.plot)
  }
  if(is.null(args.plot$ylab)){
    if(plot.density=="density"){
      args.plot=append(list(ylab="position density"), args.plot)
    }else if(plot.density=="burden-density"){
      args.plot=append(list(ylab="position density x weighted average frequency"), args.plot)
    }else{
      args.plot=append(list(ylab="position density"), args.plot)
    }
  }

  do.call(plot, args=append(list(res_cas, col=col.case, type="l", ylim=c(-max(c(res_cas,res_ctrl))/5,max(c(res_cas,res_ctrl)))), args.plot))
  do.call(lines, args=append(list(res_ctrl, col=col.ctrl), args.plot))

  if(plot.counts==TRUE){
    if(!is.null(args.text$col)){
      args.text=args.text[-which(names(args.text)=="col")]
      warning("Unused argument 'col'. Please use col.case and col.ctrl arguments.")
      }
    if(sum(nombre_mutations_cas)!=0){
      text_vector_x=position[nombre_mutations_cas!=0]
      text_vector_y=rep(-max(c(res_cas,res_ctrl))/20,length(text_vector_x))
      do.call(text, args=append(list(x=text_vector_x,y=text_vector_y, labels=nombre_mutations_cas[nombre_mutations_cas!=0],pos=4, col=col.case), args.text))
    }

    if(sum(nombre_mutations_ctrl)!=0){
      text_vector_x=position[nombre_mutations_ctrl!=0]
      text_vector_y=rep(-max(c(res_cas,res_ctrl))/10,length(text_vector_x))
      do.call(text, args=append(list(x=text_vector_x,y=text_vector_y, labels=nombre_mutations_ctrl[nombre_mutations_ctrl!=0],pos=1, col=col.ctrl), args.text))
    }
  }

  if(plot.legend==TRUE){
    if(!is.null(args.legend$col)){
      args.legend=args.legend[-which(names(args.legend)=="col")]
      warning("Unused argument 'col'. Please use col.case and col.ctrl arguments.")
    }
    if(is.null(args.legend$lwd)){args.legend=append(args.legend, list(lwd=1))}
    if(is.null(args.legend$x)){args.legend=append(args.legend, list(x="topright"))}
    do.call(legend, args=append(list(col=c(col.case,col.ctrl),legend=c("cases","controls")), args.legend))
  }

}
