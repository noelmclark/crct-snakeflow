#redirect output and messages/erros to the log
#log <- file(snakemake@log[[1]], open="wt")
#sink(log, type = "output")
#sink(log, type = "message")

library(data.table)
library(stringr)
library(colorspace)
library(pbapply)
library(parallel)
library(scales)

## define the functions from the rCNV source code first

readVCF <- function(vcf.file.path,verbose=FALSE){
  tt <- fread(vcf.file.path,sep="\t",skip="#CHROM",verbose=verbose,stringsAsFactors=FALSE)
  return(list(vcf=tt))
}

h.zygosity<-function(vcf,plot=FALSE,pops=NA,verbose=TRUE,parallel=FALSE){
  if(inherits(vcf,"list")) {
    vcf<-vcf$vcf
    gtt<-hetTgen(vcf,"GT",verbose=verbose)
  } else {
    if(any(colnames(vcf)=="REF")){gtt<-hetTgen(vcf,"GT",verbose=verbose)}
    else {gtt <-vcf}
  }
  if(parallel){
    numCores<-detectCores()-1
    cl<-makeCluster(numCores)
    clusterExport(cl, c("het.sity1","het.sity2"))
    eh<-unlist(t(parApply(cl=cl,gtt[,-c(1:4)],1,het.sity1)))
    hh<-t(parApply(cl=cl,gtt[,-c(1:4)],2,het.sity2,eh=eh))
    stopCluster(cl)
  } else {
    if(verbose){
      message("assessing per sample homozygosity")
      eh<-unlist(apply_pb(gtt[,-c(1:4)],1,het.sity1))
      hh<-t(apply_pb(gtt[,-c(1:4)],2,het.sity2,eh=eh))
    } else {
      eh<-unlist(t(apply(gtt[,-c(1:4)],1,het.sity1)))
      hh<-t(apply(gtt[,-c(1:4)],2,het.sity2,eh=eh))
    }
  }
  hh<-data.frame(rownames(hh),hh)
  colnames(hh)<-c("ind","O(Hom)","E(Hom)","total","Fis")
  rownames(hh)<-NULL
  if(plot){
    if(is.na(pops[1])){
      warning("Please provide a population list")
    } else {
      hh$pop<-pops
      if(length(unique(hh$pop))<=1){
        warning("individuals must come from at least two populations")
      }
      boxplot(hh$Fis~hh$pop,ann=F)
    }
  }
  return(hh)
}

hetTgen <- function(vcf, info.type = c("AD", "AD-tot", "GT", "GT-012", "GT-AB", "DP"), verbose = TRUE, parallel = FALSE) {
  if (inherits(vcf, "list")) { vcf <- vcf$vcf }
  if (inherits(vcf, "data.frame")) { vcf <- data.table::data.table(vcf) }
  if (any(nchar(vcf$ALT) > 1)) {
    warning("vcf file contains multi-allelic variants: \nonly bi-allelic SNPs allowed\nUse maf() to remove non-bi-allilic snps or drop minimum frequency alleles")
  }

  info.type <- match.arg(info.type)
  itype <- substr(info.type, 1, 2)

  adn <- sapply(strsplit(unlist(vcf[,"FORMAT"], use.names = FALSE), ":"), function(x) match(itype, x))
  max_adn <- max(adn) + 1L
  ind <- cbind(seq_along(adn), adn)
  xx <- data.frame(vcf[,-c(1:9)])

  h.table <- matrix(NA, nrow(xx), ncol(xx))

  process_column_0 <- function(i) {
    if (info.type == "AD-tot") {
      tmp <- stringr::str_split_fixed(xx[,i], ":", max_adn)[ind]
      tmp <- stringr::str_split_fixed(tmp, ",", 2L)
      as.numeric(tmp[,1]) + as.numeric(tmp[,2])
    } else {
      tmp <- stringr::str_split_fixed(xx[,i], ":", max_adn)[ind]
      if(info.type!="DP"){tmp[is.na(tmp) | tmp==".,."] <- "./."}
      tmp
    }
  }

  if (verbose & parallel==FALSE) {
    message("generating table")
    pb <- txtProgressBar(min = 0, max = ncol(xx), style = 3, width = 50, char = "=")
  }

  if (parallel) {
    numCores <- detectCores() - 1
    cl <- makeCluster(numCores)
    clusterExport(cl, varlist = c("xx", "max_adn", "ind", "info.type", "process_column"), envir = environment())

    results <- parLapply(cl, seq_len(ncol(xx)), function(i) {
      res <- process_column_0(i)
      res
    })
    stopCluster(cl)

    h.table <- do.call(cbind, results)
  } else {
    for (i in seq_len(ncol(xx))) {
      h.table[, i] <- process_column_0(i)
      if (verbose) {
        setTxtProgressBar(pb, i)
      }
    }
  }

  if (verbose) {
    close(pb)
  }
  if (info.type == "GT-012") {
    h.table[h.table == "0/0"] <- 0
    h.table[h.table == "1/1"] <- 1
    h.table[h.table == "1/0" | h.table == "0/1"] <- 2
    h.table[h.table == "./." | h.table == "."] <- NA
  }
  if (info.type == "GT-AB") {
    h.table[h.table == "0/0"] <- "AA"
    h.table[h.table == "1/1"] <- "BB"
    h.table[h.table == "1/0" | h.table == "0/1"] <- "AB"
    h.table[h.table == "./." | h.table == "."] <- -9
  }
  if (info.type == "AD") {
    h.table[h.table == "./." | h.table == "." | is.na(h.table)] <- "0,0"
  }
  if (info.type == "DP") {
    h.table[is.character(h.table)] <- 0
    h.table[is.na(h.table)] <- 0
  }

  het.table <- as.data.frame(cbind(vcf[, c(1:3, 5)], h.table))
  colnames(het.table) <- c("CHROM", colnames(vcf)[c(2, 3, 5, 10:ncol(vcf))])
  return(het.table)
}

allele.info<-function(X,x.norm=NULL,Fis,method=c("MedR","QN","pca","TMM","TMMex"),logratioTrim = 0.3,sumTrim = 0.05,Weighting = TRUE,Acutoff = -1e+10,plot.allele.cov=TRUE,verbose = TRUE,parallel=FALSE,...){
  method=match.arg(method)
  if(is.null(x.norm)){
    x.norm<-cpm.normal(X,method=method,logratioTrim=logratioTrim,sumTrim = sumTrim,Weighting = Weighting,Acutoff = Acutoff,verbose = verbose)
  }
  if(!inherits(x.norm,"list")){x.norm<-list(AD=x.norm)}
  if(inherits(x.norm,"list")){x.norm<-x.norm$AD}

  if(parallel){
    numCores<-detectCores()-1
    cl<-makeCluster(numCores)

    p.cal<-parApply(cl=cl,x.norm[,-c(1:4)],1,function(snp1){
      if(is.character(unname(unlist(snp1[1])))){
        y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
        y[,1]<-as.integer(y[,1])
        y[,2]<-as.integer(y[,2])} else {y<-snp1}
      rs<-rowSums(y)
      rs[rs==0]<-NA
      cv<-sd(unlist(rs),na.rm = TRUE)/mean(unlist(rs),na.rm = TRUE)
      rr1<-y[,2]/rs
      snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
      homalt<-sum(rr1==1,na.rm=TRUE)
      homref<-sum(rr1==0,na.rm=TRUE)
      covrefhomo<-sum(y[c(rr1 == 0,na.rm = TRUE),],na.rm = TRUE)
      covalthomo<-sum(y[c(rr1 == 1,na.rm = TRUE),],na.rm = TRUE)
      covalt<-sum(y[,2],na.rm = TRUE)
      covref<-sum(y[,1],na.rm = TRUE)
      NHet<-nrow(snp1het)
      if(NHet>3){
        p.all<-(covalt/(NHet+(2*homalt)))/((covalt/(NHet+(2*homalt)))+(covref/(NHet+(2*homref))))
        p.sum<-sum(snp1het[,2])/sum(snp1het)
        ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/NHet,mean.r.het=sum(snp1het[,1])/NHet,cv=cv)
      } else {
        ll<-NA
      }
      return(ll)
    })

  } else {
    if(verbose){
      message("calculating probability values of alleles")
      p.cal<-apply_pb(x.norm[,-c(1:4)],1,function(snp1){
        if(is.character(unname(unlist(snp1[1])))){
          y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
          y[,1]<-as.integer(y[,1])
          y[,2]<-as.integer(y[,2])} else {y<-snp1}
        rs<-rowSums(y)
        rs[rs==0L]<-NA
        cv<-sd(unlist(rs),na.rm = TRUE)/mean(unlist(rs),na.rm = TRUE)
        rr1<-y[,2]/rs
        snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
        homalt<-sum(rr1==1,na.rm=TRUE)
        homref<-sum(rr1==0,na.rm=TRUE)
        covrefhomo<-sum(y[c(rr1 == 0,na.rm = TRUE),],na.rm = TRUE)
        covalthomo<-sum(y[c(rr1 == 1,na.rm = TRUE),],na.rm = TRUE)
        covalt<-sum(y[,2],na.rm = TRUE)
        covref<-sum(y[,1],na.rm = TRUE)
        NHet<-nrow(snp1het)
        if(NHet>3){
          p.all<-(covalt/(NHet+(2*homalt)))/((covalt/(NHet+(2*homalt)))+(covref/(NHet+(2*homref))))
          p.sum<-sum(snp1het[,2])/sum(snp1het)
          ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/NHet,mean.r.het=sum(snp1het[,1])/NHet,cv=cv)
        } else {
          ll<-NA
        }
        return(ll)
      })
    } else {
      p.cal<-apply(x.norm[,-c(1:4)],1,function(snp1){
        if(is.character(unname(unlist(snp1[1])))){
          y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
          y[,1]<-as.integer(y[,1])
          y[,2]<-as.integer(y[,2])} else {y<-snp1}
        rs<-rowSums(y)
        rs[rs==0]<-NA
        cv<-sd(unlist(rs),na.rm = TRUE)/mean(unlist(rs),na.rm = TRUE)
        rr1<-y[,2]/rs
        snp1het<-y[-which(rr1 == 0 | rr1 == 1 | is.na(rr1)==TRUE),]
        homalt<-sum(rr1==1,na.rm=TRUE)
        homref<-sum(rr1==0,na.rm=TRUE)
        covrefhomo<-sum(y[c(rr1 == 0,na.rm = TRUE),],na.rm = TRUE)
        covalthomo<-sum(y[c(rr1 == 1,na.rm = TRUE),],na.rm = TRUE)
        covalt<-sum(y[,2],na.rm = TRUE)
        covref<-sum(y[,1],na.rm = TRUE)
        NHet<-nrow(snp1het)
        if(NHet>3){
          p.all<-(covalt/(NHet+(2*homalt)))/((covalt/(NHet+(2*homalt)))+(covref/(NHet+(2*homref))))
          p.sum<-sum(snp1het[,2])/sum(snp1het)
          ll <-data.frame(p.all,p.sum,mean.a.homo=covalthomo/(2*homalt),mean.r.homo=covrefhomo/(2*homref),mean.a.het=sum(snp1het[,2])/NHet,mean.r.het=sum(snp1het[,1])/NHet,cv=cv)
        } else {
          ll<-NA
        }
        return(ll)
      })
    }
  }

  if(is.list(p.cal)){
    p.cal<-do.call(rbind,p.cal)
  } else {
    p.cal<-t(p.cal)
  }
  p.cal[p.cal=="NaN"]<-0
  if(plot.allele.cov){
    p.list<-list(...)
    if(is.null(p.list$pch)) p.list$pch=19
    if(is.null(p.list$cex)) p.list$cex=0.6
    if(is.null(p.list$col)) p.list$col<-makeTransparent(colorspace::heat_hcl(12,h=c(0,-100),c=c(40,80),l=c(75,40),power=1)[11])
    if(is.null(p.list$lcol)) p.list$lcol="tomato"
    par(mfrow=c(2,2))
    par(mar=c(4,5,2,2))
    plot(p.cal$mean.a.homo,p.cal$mean.r.homo,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nalt. allele in homozygotes",
         ylab="Mean coverage of \n ref. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.a.het,p.cal$mean.r.het,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nalt. allele in heterozygotes",
         ylab="Mean coverage of \n ref. allele in heterozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.a.het,p.cal$mean.a.homo,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nalt. allele in heterozygotes",
         ylab="Mean coverage of \n alt. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    plot(p.cal$mean.r.het,p.cal$mean.r.homo,pch=p.list$pch,cex=p.list$cex,
         col=p.list$col,xlab="Mean coverage of \nref. allele in heterozygotes",
         ylab="Mean coverage of \n ref. allele in homozygotes",cex.lab=0.8)
    abline(0,1,col=p.list$lcol)
    par(mfrow=c(1,1))
  }

  if(parallel){
    pvals<-parLapply(1:nrow(X),get.pvals,df=X,p.cal=p.cal,cl=cl)
    stopCluster(cl)
  } else {
    if(verbose){
      message("calculating chi-square significance")
      pvals<-lapply_pb(1:nrow(X),get.pvals,df=X,p.cal=p.cal)
    } else {
      pvals<-lapply(1:nrow(X),get.pvals,df=X,p.cal=p.cal)
    }
  }
  pvals<-do.call(rbind,pvals)
  pvals<-cbind(X[,1:3],pvals)
  pvals<-na.omit(pvals)
  ht<-sig.hets(pvals,Fis,plot = FALSE, verbose = verbose)
  pvals<-data.frame(pvals,eH.pval=ht[,"eH.pval"],eH.delta=ht[,"eH.delta"],cv=na.omit(p.cal[,"cv"]))

  return(pvals)
}

het.sity1 <- function(ind){  ### calculate expected number of hets
  tab <- table(factor(ind, levels=c("0/0", "1/1", "1/0", "0/1", "./.", ".")))
  N <- tab["0/0"] + tab["1/1"] + tab["0/1"] + tab["1/0"]
  p <- (2 * tab["0/0"] + tab["0/1"] + tab["1/0"])/ (2*N)
  q <- 1 - p
  E <-(p^2+q^2)
  return(E)
}

het.sity2 <- function(ind,eh){
  tab <- table(factor(ind, levels=c("0/0", "1/1", "1/0", "0/1", "./.", ".")))
  O <- tab["0/0"] + tab["1/1"]
  N <- tab["0/0"] + tab["1/1"] + tab["0/1"] + tab["1/0"]
  E <- sum(eh[which(ind == "0/0" | ind == "1/1" | ind == "1/0" | ind == "0/1" )])
  FF <-(O-E)/(N-E)
  return(c(O,E,N,FF))
}

apply_pb <- function(X, MARGIN, FUN, ...)
{
  env <- environment()
  pb_Total <- sum(dim(X)[MARGIN])
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total,width = 50,
                       style = 3)

  wrapper <- function(...)
  {
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir= env)
    setTxtProgressBar(get("pb", envir= env),
                      curVal +1)
    FUN(...)
  }
  res <- apply(X, MARGIN, wrapper, ...)
  close(pb)
  res
}

lapply_pb <- function(X, FUN, ...)
{
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, width = 50,style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}

get.pvals<-function(x,df,p.cal){
  snp1<-unlist(df[x,-c(1:4)])
  y<-data.frame(stringr::str_split_fixed(snp1,",",n=2L))
  y[,1]<-as.integer(y[,1]);y[,2]<-as.integer(y[,2])
  rr1<-y[,2]/rowSums(y,na.rm = TRUE)
  nonhet<-which(rr1 == 0L | rr1 == 1L | is.na(rr1)==TRUE)
  if(length(nonhet)>0){snp1het<-y[-which(rr1 == 0L | rr1 == 1L | is.na(rr1)==TRUE),]} else {snp1het<-y}
  homalt<-sum(rr1==1L,na.rm=TRUE)
  homref<-sum(rr1==0L,na.rm=TRUE)
  NHet<-nrow(na.omit(snp1het))
  Nsamp <- NHet+homalt+homref
  if(NHet>3L){
    propHet<-NHet/length(na.omit(rr1))
    medRatio<-median(proportions(as.matrix(snp1het),margin = 1)[,2],na.rm = TRUE)
    p.sum<-p.cal[x,2]
    p.05<-0.5
    p.all<-p.cal[x,1]
    n<-unname(rowSums(snp1het,na.rm = TRUE))

    chi.het<-sum((((n*p.sum)-snp1het[,2])^2L/n*p.sum)+(((n*(1L-p.sum))-snp1het[,1])^2L/n*(1L-p.sum)),na.rm = TRUE)
    chi.het.sum<-chi.het
    chi.05<-sum((((n*p.05)-snp1het[,2])^2L/n*p.05)+(((n*(1L-p.05))-snp1het[,1])^2L/n*(1L-p.05)),na.rm = TRUE)
    chi.05.sum<-chi.05
    chi.all<-sum((((n*p.all)-snp1het[,2])^2L/n*p.all)+(((n*(1L-p.all))-snp1het[,1])^2L/n*(1L-p.all)),na.rm = TRUE)
    chi.all.sum<-chi.all

    z <- (n*p.sum-snp1het[,2])/sqrt(n*p.sum*(1L-p.sum))
    z.het.sum<-sum(z,na.rm = TRUE)
    z<-pnorm(z.het.sum,0,sqrt(NHet))
    z.05 <- (n*p.05-snp1het[,2])/sqrt(n*p.05*(1L-p.05))
    z.05.sum<-sum(z.05,na.rm = TRUE)
    z.05<-pnorm(z.05.sum,0,sqrt(NHet))
    z.all<- (n*p.all-snp1het[,2])/sqrt(n*p.all*(1L-p.all))
    z.all.sum<-sum(z.all,na.rm = TRUE)
    z.all<-pnorm(z.all.sum,0,sqrt(NHet))
    ll<-data.frame(NHet=NHet,propHet,medRatio,NHomRef=homref,NHomAlt=homalt,propHomAlt=homalt/Nsamp,Nsamp,
                   pAll=p.all,pHet=p.sum,fis=1L-(NHet/(2L*(homref+(NHet/2L))*(homalt+(NHet/2L)))),
                   z.het=ifelse(z>0.5, (1L-z)*2L, z*2L),
                   z.05=ifelse(z.05>0.5, (1L-z.05)*2L, z.05*2L),
                   z.all=ifelse(z.all>0.5, (1L-z.all)*2L, z.all*2L),
                   chi.het=pchisq(chi.het,NHet-1L,lower.tail=F),
                   chi.05=pchisq(chi.05,NHet-1L,lower.tail = F),
                   chi.all=pchisq(chi.all,NHet-1L,lower.tail = F),
                   z.het.sum,z.05.sum,z.all.sum,chi.het.sum,chi.05.sum,chi.all.sum)
  } else {
    ll<-NA
  }
  return(ll)
}

sig.hets<-function(a.info,Fis,method=c("chi.sq","fisher"),plot=TRUE,verbose=TRUE,...){
  if(!any(colnames(a.info)=="NHomRef")){
    if(verbose){message("assessing excess of heterozygotes")
      df<-apply(a.info,1,get.eHpvals,method=method,Fis=Fis)
    } else {df<-apply(a.info,1,get.eHpvals,method=method,Fis=Fis)}
    df<-data.frame(do.call(rbind,df))
    colnames(df)<-c("medRatio","propHomAlt","propHet","p2","het","q2","eH.pval","eH.delta")
    #df<-na.omit(df)
  } else {
    d<-a.info[,c("NHomRef","NHet","NHomAlt","Nsamp")]
    colnames(d)<-c("h1","het","h2","truNsample")
    method<-match.arg(method)
    if(verbose){
      message("assessing excess of heterozygotes")
      df<-data.frame(t(apply(d,1,ex.prop,method=method,Fis=Fis)))
    } else {
      df<-data.frame(t(apply(d,1,ex.prop,method=method,Fis=Fis)))
    }
    colnames(df)<-c("p2","het","q2","eH.pval","eH.delta")
    df$propHomAlt <- a.info$propHomAlt
    df$propHet<-a.info$propHet
    df$medRatio<-a.info$medRatio
  }

  df$dup.stat<-"non-deviant";df$dup.stat[which(df$eH.pval < 0.05/nrow(df) & df$eH.delta > 0 )]<-"deviant"
  d<-na.omit(df[,c("p2","het","q2")])
  df<-na.omit(data.frame(cbind(a.info[,c(1:4)],df[,c("medRatio","propHomAlt","propHet","eH.pval","eH.delta","dup.stat")]),row.names = NULL))

  if(plot){
    l<-list(...)
    if(is.null(l$cex)) l$cex=0.2
    if(is.null(l$pch)) l$pch=19
    if(is.null(l$xlim)) l$xlim=c(0,1)
    if(is.null(l$ylim)) l$ylim=c(0,1)
    if(is.null(l$col)) cols<-makeTransparent(colorspace::rainbow_hcl(2),alpha=0.3) else cols<-makeTransparent(l$col,alpha=0.3)

    d$Color <- cols[1]
    d$Color [which(df$dup.stat=="deviant")]<- cols[2]#& df$delta > 0
    plot(df$propHet~df$propHomAlt, pch=l$pch, cex=l$cex,col=d$Color,xlim=l$xlim,ylim=l$ylim,
         xlab="Proportion of Alternate Homozygotes",ylab="Proportion of Heterozygotes")
    lines((smm<-smooth.spline(d$het~d$q2)),col="blue")
    legend("bottomright", c("non-deviants","deviants","expected"), col = c(cols,"blue"), lty = c(0, 0, 1), lwd = c(0, 0, 1),pch = c(l$pch, l$pch, NA),
           cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")
  }

  return(df)
}

cpm.normal <- function(het.table, method=c("MedR","QN","pca","TMM","TMMex"),
                       logratioTrim=.3, sumTrim=0.05, Weighting=TRUE,
                       Acutoff=-1e10, verbose=TRUE, plot=TRUE){
  method<-match.arg(method)
  if(length(method)>1){method="MedR"}
  dm <- dim(het.table)
  pb_Total <- dm[2] - 4L
  if(verbose){
    message("calculating normalization factor")
    pb <- txtProgressBar(min = 0, max = pb_Total, width = 50, style = 3)
  }
  #  tdep<-parApply(cl=cl,(het.table[,-c(1:4)],2,function(tmp){
  y1 <- y2 <- matrix(NA_integer_, dm[1], pb_Total)
  for(i in seq_len(pb_Total)){
    if (verbose) setTxtProgressBar(pb, i)
    tmp <- stringr::str_split_fixed(het.table[,i+4L], ",", 2L)
    y1[, i] <- as.integer(tmp[,1])
    y2[, i] <- as.integer(tmp[,2])
  }
  if(verbose) close(pb)
  tdep <- y1 + y2
  colnames(tdep) <- colnames(het.table[-c(1:4)])

  #find and warn about outliers
  ot<-boxplot.stats(colSums(tdep,na.rm = T))$out
  cl<-rep("dodgerblue",ncol(tdep))
  ot.ind<-which(colnames(tdep) %in% names(ot))
  cl[ot.ind]<-2
  if(length(ot)>0){
    if(plot){
      barplot(colSums(tdep,na.rm = T), col=cl, border=NA, xlab="sample",
              ylab="total depth")
    }
    message("OUTLIERS DETECTED\nConsider removing the samples:")
    cat(colnames(tdep)[ot.ind])
  }

  if(method=="TMM" | method=="TMMex"){
    if(verbose)   message("\ncalculating normalized depth")
    nf<-norm.fact(tdep, method=method, logratioTrim=logratioTrim,
                  sumTrim=sumTrim, Weighting=Weighting, Acutoff=Acutoff)
    sc <- median(colSums(tdep)) / (nf[,1]*nf[,2])
    y1 <- round(y1 * rep(sc, each=dm[1]), 2)
    y2 <- round(y2 * rep(sc, each=dm[1]), 2)
    out <- paste0(y1, ",", y2)
    attributes(out) <- attributes(tdep)
  } else if(method=="MedR") {
    pseudo <- apply(tdep,1,function(xx){exp(mean(log(as.numeric(xx)[as.numeric(xx)>0])))})
    nf <- apply(tdep,2,function(xx){median(as.numeric(xx)/pseudo,na.rm=T)})
    if(verbose)   message("\ncalculating normalized depth")
    y1 <- round(y1 / rep(nf, each=dm[1]), 0)
    y2 <- round(y2 / rep(nf, each=dm[1]), 0)
    out <- paste0(y1, ",", y2)
    attributes(out) <- attributes(tdep)
  } else if(method=="QN"){
    out <- do.call(cbind,quantile_normalisation(tdep,het.table,verbose=verbose,cl=cl))
  } else if(method=="pca"){
    if(verbose){message("\ncalculating normalized depth")}
    new.mat <- t(tdep) ### check the direction to confirm if this step need to be done
    colmean <- colMeans(new.mat)
    colsd <- apply(new.mat, 2, sd)
    new.mat <- apply(new.mat, 2, scale,scale = TRUE) #### essential before SVD
    test.la.svd <- La.svd(new.mat)
    u <- test.la.svd$u
    d <- test.la.svd$d
    vt <- test.la.svd$vt
    ## optimal PCs to remove with d values
    ddl<-NULL
    for(i in seq_along(1:50)){
      if(i<50) ddl[i]<-d[i+1]-d[i]
    }
    ## modified Kaiser's Rule: Sample variation of eigen values smaller than 0.7 should be kept (i.e., the first eigen value < 0.7)
    rmpc<-min(which(abs(ddl)<0.7))
    #plot(d[1:50],pch=19,cex=0.5)
    #points(rmpc,d[rmpc],cex=1.5,col="red")
    d[1:rmpc] <- 0
    out <- u %*% diag(d) %*% vt
    out <- apply(out, 1, FUN = function(x){round(x*colsd + colmean,0)})#### back transform to depth matrix
    out[out<0]<-0
    out<-t(out)

    ht<-het.table[,-c(1:4)]
    tout<-NULL
    for(i in 1:ncol(ht)){
      if(verbose){
        pb <- txtProgressBar(min = 0, max = ncol(ht), style = 3, width = 50, char = "=")
        setTxtProgressBar(pb, i)
      }
      tmp <- stringr::str_split_fixed(ht[,i], ",", 2L)
      tt<-matrix(NA,nrow = nrow(tmp),ncol = 2)
      tt[,1]<-as.integer(tmp[,1])
      tt[,2]<-as.integer(tmp[,2])
      tt <- proportions(tt,margin = 1)
      tt[is.na(tt)]<-0
      tt<-tt*out[i,]
      tt<-paste0(round(tt[,1],0),",",round(tt[,2],0))
      tout<-cbind(tout,tt)
    }
    out<- out #t(tout)
  }
  #  browser()
  out<-data.frame(het.table[,c(1:4)], out)
  colnames(out)<-colnames(het.table)
  return(list(AD=out,outliers=data.frame(column=(ot.ind+4),sample=colnames(tdep)[ot.ind])))
}


makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}

allele.info.WGS <- function (ad, gt, fis = NULL, vcf = NULL, parallel = FALSE, numCores = NULL, ...)
{
  ll <- list(...)
  # <---per sample total depth
  message("step 1/5: calculating depth values")
  if (parallel) {
    if (is.null(numCores)) {
      numCores <- detectCores() - 1
      cl <- makeCluster(numCores)
    }
    tot <- parApply(cl = cl, ad[, -(1:4)], 2, FUN = function(x) {
      unlist(lapply(x, FUN = function(x) {
        sum(as.numeric(unlist(strsplit(x, split = ","))))
      }))
    })
  }
  else {
    tot <- apply(ad[, -(1:4)], 2, FUN = function(x) {
      unlist(lapply(x, FUN = function(x) {
        sum(as.numeric(unlist(strsplit(x, split = ","))))
      }))
    })
  }
  tot <- apply(tot, 2, as.numeric)
  tot <- as.data.frame(tot)

  # <--- global inbreeding coefficient
  if (is.null(fis)) {
    if (is.null(vcf)) {
      stop("No fis (global inbreeding coefficient) or vcf provided \n\n           Either calculate fis using h.zygosity() function or provide a vcf")
    }
    else {
      fis <- h.zygosity(vcf)
      fis <- mean(fis$Fis)
    }
  }

  # <--- statistics for depth simulation, works only for WGS
  if (parallel) {
    colstat <- nb_stats(tot,cl=cl)
  }
  else{
    colstat <- nb_stats(tot,cl=NULL)
  }

  # <--- likelihood for observing the allele-specific depth under N = 2 and N = 4
  message("step 2/5: calculating allele specific depth-likelihood")
  if (parallel) {
    gene_prob <- t(parApply(ad[, -(1:4)], 1, FUN = cal_geno_lld2,
                            inb = fis, nb_stat = colstat, cl = cl))
  }
  else {
    gene_prob <- t(apply(ad[, -(1:4)], 1, FUN = cal_geno_lld2,
                         inb = fis, nb_stat = colstat))
  }

  # <--- likelihood for observing the total depth under N = 2 and N = 4
  message("step 3/5: calculating total depth-specific likelihood")
  if (parallel) {
    depth_prob <- parApply(tot, 2, FUN = cal_depth_lld_indi2,
                           nb_stat = colstat, cl = cl)
  }
  else {
    depth_prob <- apply(tot, 2, FUN = cal_depth_lld_indi2,
                        nb_stat = colstat)
  }

  # <--- likelihood ratio: N=4 against N=2 (allele-specific depth + total depth)
  lhr <- 2 * (gene_prob + depth_prob)
  lhr[lhr < 0] <- 0
  lhr <- rowSums(lhr)

  # <--- permutation test
  message("step 4/5: performing permutation accross SNPs and samples")
  # likelihood ratio for each SNP each individual
  lld.ratio.geno <- na.omit(gene_prob[gt[, -(1:4)] == "0/1"])
  lld.ratio.geno[lld.ratio.geno < 0] <- 0
  lld.ratio.dep <- depth_prob
  lld.ratio.dep[lld.ratio.dep < 0] <- 0
  # the 0.95 significant threshold considering both allelic ratio and depth by permutation,
  # significant level can be adjusted
  nhet <- rowSums(gt[, -(1:4)] == "0/1")
  # 0.95 threshold
  lld.ratio.dep.0.95 <- quantile(replicate(10000, sum(sample(lld.ratio.dep,
                                                             ncol(tot)))), 0.95)
  if (parallel) {
    lld.ratio.geno.0.95 <- unlist(parLapply(cl = cl, unique(nhet),
                                            sample_and_quantile, samples = sample(lld.ratio.geno,
                                                                                  1e+05), nrep = 10000, quan = 0.95))
  }
  else {
    lld.ratio.geno.0.95 <- unlist(lapply(unique(nhet), sample_and_quantile,
                                         samples = sample(lld.ratio.geno, 1e+05), nrep = 10000,
                                         quan = 0.95))
  }
  thre0.95 <- 2 * (lld.ratio.dep.0.95 + lld.ratio.geno.0.95)
  names(thre0.95) <- unique(nhet)
  # 0.99 threshold
  lld.ratio.dep.0.99 <- quantile(replicate(10000, sum(sample(lld.ratio.dep,
                                                             ncol(tot)))), 0.99)
  if (parallel) {
    lld.ratio.geno.0.99 <- unlist(parLapply(cl = cl, unique(nhet),
                                            sample_and_quantile, samples = sample(lld.ratio.geno,
                                                                                  1e+05), nrep = 10000, quan = 0.99))
  }
  else {
    lld.ratio.geno.0.99 <- unlist(lapply(unique(nhet), sample_and_quantile,
                                         samples = sample(lld.ratio.geno, 1e+05), nrep = 10000,
                                         quan = 0.99))
  }
  thre0.99 <- 2 * (lld.ratio.dep.0.99 + lld.ratio.geno.0.99)
  names(thre0.99) <- unique(nhet)
  if (parallel) {
    stopCluster(cl)
  }

  # <--- create output info table
  message("step 5/5: finishing....")
  info <- data.frame(ad[, 1:4], mdep = rowMeans(tot), Nhet = nhet,
                     lhr = lhr, perm.95 = thre0.95[as.character(nhet)], perm.99 = thre0.99[as.character(nhet)])
  info$dup.stat <- info$lhr > info$perm.95
  info$dup.stat <- ifelse(info$dup.stat, "duplicated", "non-duplicated")
  return(info)
}

##### internal functions of likelihood ratio calculations ######
nb_stats <- function(tot.tab,cl = NULL){
  ##### correlation between mean depth and sd, calculate parameter for nb distribution
  if (is.null(cl)){
    colstat <- apply(tot.tab, 2, FUN = function(x){
      mean <- mean(x[x> quantile(x,0.05) & x < quantile(x,0.95)])
      sd <- sd(x[x> quantile(x,0.05) & x < quantile(x,0.95)])
      size <- MASS::fitdistr(x, "Negative Binomial")$estimate["size"]
      return(data.frame(mean,sd,size))
    })
  }else {
    colstat <- parApply(tot.tab, 2, FUN = function(x){
      mean <- mean(x[x> quantile(x,0.05) & x < quantile(x,0.95)])
      sd <- sd(x[x> quantile(x,0.05) & x < quantile(x,0.95)])
      size <- MASS::fitdistr(x, "Negative Binomial")$estimate["size"]
      return(data.frame(mean,sd,size))
    },cl=cl)
  }

  par(mfrow =c(4,5))
  colstat <- do.call(rbind,colstat)
  plot(colstat$mean,colstat$sd,xlab = "mean",ylab = "sd") ## check correlation


  fit <- lm(sd~mean,data = colstat)
  colstat$mean2 <- 2*colstat$mean              ## expected mean and depth for N=4
  colstat$sd2 <- predict(fit,newdata = data.frame(mean = colstat$mean2))
  colstat$size2 <- colstat$mean2^2/(colstat$sd2^2-colstat$mean2)       ## expected mean and depth for N=4

  ############# check if nb distribution fit well, output first 20, can be disabled
  for (i in 1:19) {
    hist(tot.tab[,i], breaks = max(tot.tab[,i]), freq = FALSE, col = "gray", main = "Negative Binomial Fit",
         xlab = "Counts", xlim = c(0, 200))
    suppressWarnings(lines(0:200, dnbinom(0:200, size = colstat$size[i], mu = colstat$mean[i]),
                           col = "red", lwd = 2))
    suppressWarnings(lines(0:200, dnbinom(0:200, size = colstat$size2[i], mu = colstat$mean2[i]),
                           col = "black", lwd = 2))
  }
  return(colstat)
}








#############################################################
## run the rCNV flow using defined functions

#load the paths from the snakemake rule
vcf_file <- snakemake@input$vcf
outfile <- snakemake@output$tsv
outwgsfile <- snakemake@output$wgs

vcf<-readVCF(vcf_file)

fis<-mean(h.zygosity(vcf)$Fis)
gt<-hetTgen(vcf,"GT")
ad<-hetTgen(vcf,"AD")
ad[ad==""] <- "0,0"
adnor<-cpm.normal(ad, method = "MedR")
out<-allele.info(X = ad, x.norm = adnor, Fis = fis)
out.wgs<-allele.info.WGS(ad, gt, fis = fis)

# X is the corrected non-normalized allele depth table and x.norm is the normalized allele depth table
# out is a data frame with duplication status of each snp, and other stats.

#write the output to the file path
write_tsv(out, file = outfile)
write_tsv(out.wgs, file = outwgsfile)