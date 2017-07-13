#16 Jan 2017
#TNS
#script to look for presence of alternative paths under a different evolutionary scheme, that incorporates the quantitative nature of our meanF estimates rather than the threshold function applied up until now (plus a couple simpler alternative schemes in script 9)

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)
library(rgexf)
library(igraph)
library(Hmisc)

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)

#before I can test for paths in which quantitative measure of SRE-binding meanF goes up, need to determine SEMs

#approach: for variants with similar numbers of cfus in each repliacte, bin by average cfu of two rep (or min cfu?); in each bin, estimate standard deviation of sampling distribution from distribution of meanF.repX - meanF.pooled; then, for each variant, SEM is the s.d. of the sampling distribution for its relevant coverage bin given the number of pooled cfus with which it was determined
#may have to do separately for e.g. nulls and positives, as SEM might vary as a function of meanF
#can also compare to s.d. of null sampling distribution as a function of coverage

breaks.ERE.11P <- cut2(dt.11P.coding[ERE.rep1.cfu > 0.8*ERE.rep2.cfu & ERE.rep1.cfu < 1.2*ERE.rep2.cfu, (ERE.rep1.cfu+ERE.rep2.cfu)/2],m=1000,g=31,onlycuts=T);breaks.ERE.11P[31] <- breaks.ERE.11P[31]+500000
breaks.SRE.11P <- cut2(dt.11P.coding[SRE.rep1.cfu > 0.8*SRE.rep2.cfu & SRE.rep1.cfu < 1.2*SRE.rep2.cfu, (SRE.rep1.cfu+SRE.rep2.cfu)/2],m=1000,g=31,onlycuts=T);breaks.SRE.11P[31] <- breaks.SRE.11P[31]+500000

ERE.sample.dists.11P <- data.frame(bin=1:30)
for(i in 1:nrow(ERE.sample.dists.11P)){
  ERE.sample.dists.11P$range.cfu[i] <- list(c(breaks.ERE.11P[i],breaks.ERE.11P[i+1]))
  data <- dt.11P.coding[ERE.rep1.cfu > 0.8*ERE.rep2.cfu & ERE.rep1.cfu < 1.2*ERE.rep2.cfu & (ERE.rep1.cfu+ERE.rep2.cfu)/2 > ERE.sample.dists.11P$range.cfu[i][[1]][[1]] & (ERE.rep1.cfu+ERE.rep2.cfu)/2 < ERE.sample.dists.11P$range.cfu[i][[1]][[2]] & !is.na(ERE.rep1.meanF) & !is.na(ERE.rep2.meanF) & !is.na(ERE.pooled.meanF),]
  ERE.sample.dists.11P$n.var[i] <- nrow(data)
  ERE.sample.dists.11P$median.cfu[i] <- median((data$ERE.rep1.cfu + data$ERE.rep2.cfu)/2)
  ERE.sample.dists.11P$residuals[i] <- list(c((data$ERE.rep1.meanF - data$ERE.pooled.meanF),(data$ERE.rep2.meanF - data$ERE.pooled.meanF)))
  ERE.sample.dists.11P$sample.sd[i] <- sd(ERE.sample.dists.11P$residuals[i][[1]])
}

SRE.sample.dists.11P <- data.frame(bin=1:30)
for(i in 1:nrow(SRE.sample.dists.11P)){
  SRE.sample.dists.11P$range.cfu[i] <- list(c(breaks.SRE.11P[i],breaks.SRE.11P[i+1]))
  data <- dt.11P.coding[SRE.rep1.cfu > 0.8*SRE.rep2.cfu & SRE.rep1.cfu < 1.2*SRE.rep2.cfu & (SRE.rep1.cfu+SRE.rep2.cfu)/2 > SRE.sample.dists.11P$range.cfu[i][[1]][[1]] & (SRE.rep1.cfu+SRE.rep2.cfu)/2 < SRE.sample.dists.11P$range.cfu[i][[1]][[2]] & !is.na(SRE.rep1.meanF) & !is.na(SRE.rep2.meanF) & !is.na(SRE.pooled.meanF),]
  SRE.sample.dists.11P$n.var[i] <- nrow(data)
  SRE.sample.dists.11P$median.cfu[i] <- median((data$SRE.rep1.cfu + data$SRE.rep2.cfu)/2)
  SRE.sample.dists.11P$residuals[i] <- list(c((data$SRE.rep1.meanF - data$SRE.pooled.meanF),(data$SRE.rep2.meanF - data$SRE.pooled.meanF)))
  SRE.sample.dists.11P$sample.sd[i] <- sd(SRE.sample.dists.11P$residuals[i][[1]])
}

#how do these compare to st.dev of sampling distribution of stop-codon-containing variants?
load("3_classify-variants_define-distributions_out/ERE.pooled.null.11P.Rda")
load("3_classify-variants_define-distributions_out/SRE.pooled.null.11P.Rda")

#for purposes of supplemental figure, want this same plot for the SR1 data
dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

breaks.ERE.SR1 <- cut2(dt.SR1.coding[ERE.rep1.cfu > 0.8*ERE.rep2.cfu & ERE.rep1.cfu < 1.2*ERE.rep2.cfu, (ERE.rep1.cfu+ERE.rep2.cfu)/2],m=1000,g=31,onlycuts=T);breaks.ERE.SR1[31] <- breaks.ERE.SR1[31]+500000
breaks.SRE.SR1 <- cut2(dt.SR1.coding[SRE.rep1.cfu > 0.8*SRE.rep2.cfu & SRE.rep1.cfu < 1.2*SRE.rep2.cfu, (SRE.rep1.cfu+SRE.rep2.cfu)/2],m=1000,g=31,onlycuts=T);breaks.SRE.SR1[31] <- breaks.SRE.SR1[31]+500000

ERE.sample.dists.SR1 <- data.frame(bin=1:30)
for(i in 1:nrow(ERE.sample.dists.SR1)){
  ERE.sample.dists.SR1$range.cfu[i] <- list(c(breaks.ERE.SR1[i],breaks.ERE.SR1[i+1]))
  data <- dt.SR1.coding[ERE.rep1.cfu > 0.8*ERE.rep2.cfu & ERE.rep1.cfu < 1.2*ERE.rep2.cfu & (ERE.rep1.cfu+ERE.rep2.cfu)/2 > ERE.sample.dists.SR1$range.cfu[i][[1]][[1]] & (ERE.rep1.cfu+ERE.rep2.cfu)/2 < ERE.sample.dists.SR1$range.cfu[i][[1]][[2]] & !is.na(ERE.rep1.meanF) & !is.na(ERE.rep2.meanF) & !is.na(ERE.pooled.meanF),]
  ERE.sample.dists.SR1$n.var[i] <- nrow(data)
  ERE.sample.dists.SR1$median.cfu[i] <- median((data$ERE.rep1.cfu + data$ERE.rep2.cfu)/2)
  ERE.sample.dists.SR1$residuals[i] <- list(c((data$ERE.rep1.meanF - data$ERE.pooled.meanF),(data$ERE.rep2.meanF - data$ERE.pooled.meanF)))
  ERE.sample.dists.SR1$sample.sd[i] <- sd(ERE.sample.dists.SR1$residuals[i][[1]])
}

SRE.sample.dists.SR1 <- data.frame(bin=1:30)
for(i in 1:nrow(SRE.sample.dists.SR1)){
  SRE.sample.dists.SR1$range.cfu[i] <- list(c(breaks.SRE.SR1[i],breaks.SRE.SR1[i+1]))
  data <- dt.SR1.coding[SRE.rep1.cfu > 0.8*SRE.rep2.cfu & SRE.rep1.cfu < 1.2*SRE.rep2.cfu & (SRE.rep1.cfu+SRE.rep2.cfu)/2 > SRE.sample.dists.SR1$range.cfu[i][[1]][[1]] & (SRE.rep1.cfu+SRE.rep2.cfu)/2 < SRE.sample.dists.SR1$range.cfu[i][[1]][[2]] & !is.na(SRE.rep1.meanF) & !is.na(SRE.rep2.meanF) & !is.na(SRE.pooled.meanF),]
  SRE.sample.dists.SR1$n.var[i] <- nrow(data)
  SRE.sample.dists.SR1$median.cfu[i] <- median((data$SRE.rep1.cfu + data$SRE.rep2.cfu)/2)
  SRE.sample.dists.SR1$residuals[i] <- list(c((data$SRE.rep1.meanF - data$SRE.pooled.meanF),(data$SRE.rep2.meanF - data$SRE.pooled.meanF)))
  SRE.sample.dists.SR1$sample.sd[i] <- sd(SRE.sample.dists.SR1$residuals[i][[1]])
}

#how do these compare to st.dev of sampling idstribution of stop-codon-containing variants?
load("3_classify-variants_define-distributions_out/ERE.pooled.null.SR1.Rda")
load("3_classify-variants_define-distributions_out/SRE.pooled.null.SR1.Rda")

pdf(file="9_quant-scale-landscape_out/plots_standard-error-v-coverage.pdf",5,4,useDingbats=F)
par(mfrow=c(2,2),mar=c(3,4,0,0),oma=c(2,2,2,2))
#11P
plot(ERE.sample.dists.11P$median.cfu, ERE.sample.dists.11P$sample.sd, log="x",type="l",lty=1,lwd=2,col=rgb(100,69,155,maxColorValue=255),ylim=c(0.05,0.4),xlim=c(1,10000),ylab="Estimated standard error",xlab="")
points(ERE.pooled.null.11P$median.cfu, ERE.pooled.null.11P$sd.meanF, type="l",lty=2,lwd=2,col=rgb(100,69,155,maxColorValue=255))
points(SRE.sample.dists.11P$median.cfu, SRE.sample.dists.11P$sample.sd,type="l",lty=1,lwd=2,col=rgb(5,172,72,maxColorValue=255))
points(SRE.pooled.null.11P$median.cfu, SRE.pooled.null.11P$sd.meanF, type="l",lty=2,lwd=2,col=rgb(5,172,72,maxColorValue=255))

#SR1
plot(ERE.sample.dists.SR1$median.cfu, ERE.sample.dists.SR1$sample.sd, log="x",type="l",lty=1,lwd=2,col=rgb(100,69,155,maxColorValue=255),ylim=c(0.05,0.4),xlim=c(1,10000),ylab="",xlab="")
points(ERE.pooled.null.SR1$median.cfu, ERE.pooled.null.SR1$sd.meanF, type="l",lty=2,lwd=2,col=rgb(100,69,155,maxColorValue=255))
points(SRE.sample.dists.SR1$median.cfu, SRE.sample.dists.SR1$sample.sd,type="l",lty=1,lwd=2,col=rgb(5,172,72,maxColorValue=255))
points(SRE.pooled.null.SR1$median.cfu, SRE.pooled.null.SR1$sd.meanF, type="l",lty=2,lwd=2,col=rgb(5,172,72,maxColorValue=255))

#11P
plot(ecdf(log10(dt.11P.coding$ERE.pooled.cfu[dt.11P.coding$ERE.pooled.cfu>0])),xlim=c(log10(1),log10(10000)),col=rgb(100,69,155,maxColorValue=255),xlab="Number of cells sampled",ylab="Cumulative density",main="")
plot(ecdf(log10(dt.11P.coding$SRE.pooled.cfu[dt.11P.coding$SRE.pooled.cfu>0])),xlim=c(log10(1),log10(10000)),col=rgb(5,172,72,maxColorValue=255),lwd=2,add=T)

#SR1
plot(ecdf(log10(dt.SR1.coding$ERE.pooled.cfu[dt.SR1.coding$ERE.pooled.cfu>0])),xlim=c(log10(1),log10(10000)),col=rgb(100,69,155,maxColorValue=255),lwd=2,main="",xlab="",ylab="")
plot(ecdf(log10(dt.SR1.coding$SRE.pooled.cfu[dt.SR1.coding$SRE.pooled.cfu>0])),xlim=c(log10(1),log10(10000)),col=rgb(5,172,72,maxColorValue=255),lwd=2,add=T)
dev.off()

##############################################################################
#assign SE to each variant via binned SE estimates given above
assign.SE <- function(cfu, RE, bg){
  if(cfu < 0.5){
    return(as.numeric(NA))
  }else if(cfu > 0.5){
    if(RE=="ERE"){
      if(bg=="11P"){bins<-ERE.sample.dists.11P}else if(bg=="SR1"){bins<-ERE.sample.dists.SR1}
    }else if(RE=="SRE"){
      if(bg=="11P"){bins<-SRE.sample.dists.11P}else if(bg=="SR1"){bins<-SRE.sample.dists.SR1}
    }
    x <- 1
    while(bins[x,"range.cfu"][[1]][2] < cfu){x <- x+1}
    return(bins[x,"sample.sd"])
  }
}

dt.11P.coding[,ERE.SE.meanF := assign.SE(ERE.pooled.cfu, "ERE", "11P"),by=AAseq]
dt.11P.coding[,SRE.SE.meanF := assign.SE(SRE.pooled.cfu, "SRE", "11P"),by=AAseq]

dt.SR1.coding[,ERE.SE.meanF := assign.SE(ERE.pooled.cfu, "ERE", "SR1"),by=AAseq]
dt.SR1.coding[,SRE.SE.meanF := assign.SE(SRE.pooled.cfu, "SRE", "SR1"),by=AAseq]

save(dt.11P.coding, file="9_quant-scale-landscape_out/dt.11P.coding.Rda")
save(dt.SR1.coding, file="9_quant-scale-landscape_out/dt.SR1.coding.Rda")

##############################################################################
#first, are there paths from EGKA to SRE-specificity that a) maintain "functionality" (strong on at least one RE) *and* are monotonically increasing in SRE-activity?
#11P space
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,ERE.SE.meanF,SRE.SE.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,ERE.SE.meanF,SRE.SE.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,ERE.SE.meanF,SRE.SE.meanF,specificity)]

#read in list of all 1-letter amino acid codes, and a table that gives all amino acids that can be encoded in adjacent (single-nt mutation) codons
AAs <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
AA.nt.neighbors <- read.csv(file="7_assess-connectivity_in/AA_single-mutant-codon-neighbors.csv",header=T,row.names=1,quote="")

#define function that takes a 4-aa RH sequence and returns all possible single-mutant neighbors (nt space)
get.Hamming1.nt <- function(AAseq){
  seq <- unlist(strsplit(as.character(AAseq),split=""))
  mutant <- vector()
  for(i in AAs[AA.nt.neighbors[seq[1],]==1]){
    mutant <- c(mutant,paste(i,seq[2],seq[3],seq[4],sep=""))
  }
  for(i in AAs[AA.nt.neighbors[seq[2],]==1]){
    mutant <- c(mutant,paste(seq[1],i,seq[3],seq[4],sep=""))
  }
  for(i in AAs[AA.nt.neighbors[seq[3],]==1]){
    mutant <- c(mutant,paste(seq[1],seq[2],i,seq[4],sep=""))
  }
  for(i in AAs[AA.nt.neighbors[seq[4],]==1]){
    mutant <- c(mutant,paste(seq[1],seq[2],seq[3],i,sep=""))
  }
  return(mutant)
}

#for each ERE.specific variant, output 3mut paths to SRE with monotonically increasing SRE-activity, defined as when the 95% confidence intervals for SRE activity between two variants don't overlap
#do two different ways: first, ask about the ones EGKA can access in <=3mut steps, and then do in the context of an entire adjacency matrix, ask how many places EGKA can end up
#old function for delimiting 3mut paths, in which any step among functional variants is allowed (barring functional reversions)
get.3mut.paths.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P$AAseq),prom=as.character(promiscuous.11P$AAseq),SRE.sp=as.character(SRE.specific.11P$AAseq)){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  if(seq %in% ERE.sp){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(ERE.sp,prom,SRE.sp)]
  }else if(seq %in% prom){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(prom,SRE.sp)]
  }else{
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(SRE.sp)]
  }
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    if(i %in% ERE.sp){
      mutant2 <- c(mutant2, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)])
    }else if(i %in% prom){
      mutant2 <- c(mutant2, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)])
    }else{
      mutant2 <- c(mutant2, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)])
    }
  }
  mutant2 <- unique(mutant2)
  mutant3 <- vector(mode="character")
  for(i in mutant2){
    if(i %in% ERE.sp){
      mutant3 <- c(mutant3, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)])
    }else if(i %in% prom){
      mutant3 <- c(mutant3, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)])
    }else{
      mutant3 <- c(mutant3, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)])
    }
  }
  mutant3 <- unique(mutant3)
  return(list(mut0=seq,mut1=mutant1, mut2=mutant2, mut3=mutant3))
}

#new function for delimiting 3mut paths, in which a step is only allowed if it involves an increase in SRE-binding activity (no overlap in SRE meanF 90% confidence inerval)
get.3mut.paths.sig.increasing.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P$AAseq),prom=as.character(promiscuous.11P$AAseq),SRE.sp=as.character(SRE.specific.11P$AAseq)){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  if(seq %in% ERE.sp){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(ERE.sp,prom,SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[seq,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq])
  }else if(seq %in% prom){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(prom,SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[seq,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq])
  }else{
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[seq,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq])
  }
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    if(i %in% ERE.sp){
      muts2 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)]
      mutant2 <- c(mutant2,as.character(dt.11P.coding[AAseq %in% muts2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq]))
    }else if(i %in% prom){
      muts2 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)]
      mutant2 <- c(mutant2,as.character(dt.11P.coding[AAseq %in% muts2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq]))
    }else{
      muts2 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)]
      mutant2 <- c(mutant2,as.character(dt.11P.coding[AAseq %in% muts2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq]))
    }
  }
  mutant2 <- unique(mutant2)
  mutant3 <- vector(mode="character")
  for(i in mutant2){
    if(i %in% ERE.sp){
      muts3 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)]
      mutant3 <- c(mutant3,as.character(dt.11P.coding[AAseq %in% muts3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq]))
    }else if(i %in% prom){
      muts3 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)]
      mutant3 <- c(mutant3,as.character(dt.11P.coding[AAseq %in% muts3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq]))
    }else{
      muts3 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)]
      mutant3 <- c(mutant3,as.character(dt.11P.coding[AAseq %in% muts3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)],AAseq]))
    }
  }
  mutant3 <- unique(mutant3)
  return(list(mut0=seq,mut1=mutant1, mut2=mutant2, mut3=mutant3))
}

muts <- get.3mut.paths.sig.increasing.11P("EGKA")
all.muts <- unique(unlist(muts))
ERE.sp <- unique(unlist(muts)[unlist(muts) %in% ERE.specific.11P$AAseq])
SRE.sp <- unique(unlist(muts)[unlist(muts) %in% SRE.specific.11P$AAseq])
prom <- unique(unlist(muts)[unlist(muts) %in% promiscuous.11P$AAseq])

#on spoke-o-gram, need to highlight paths from variants that are allowed (sig increase in SRE binding)
#change the RH under x, gives all single-nt-mutant neighbors with significan increases in SRE meanF
x <- "GAKA";as.character(dt.11P.coding[get.Hamming1.nt(x)][(SRE.pooled.meanF - 1.645*SRE.SE.meanF) > (dt.11P.coding[x, SRE.pooled.meanF] + 1.645*dt.11P.coding[x,SRE.SE.meanF]),AAseq])[as.character(dt.11P.coding[get.Hamming1.nt(x)][(SRE.pooled.meanF - 1.645*SRE.SE.meanF) > (dt.11P.coding[x, SRE.pooled.meanF] + 1.645*dt.11P.coding[x,SRE.SE.meanF]),AAseq]) %in% unlist(get.3mut.paths.11P("EGKA"))]

#pdf(file="9_quant-scale-landscape_out/EGKA-3-mut-paths_sig-increase-SRE_meanF-plot.pdf",4,4.5,useDingbats=F)
plot(dt.11P.coding$SRE.pooled.meanF,dt.11P.coding$ERE.pooled.meanF,pch="",xlab="SRE meanF (a.u.), pooled reps",ylab="ERE meanF (a.u.), pooled reps",xlim=c(4.5,10.2),ylim=c(4.5,10.2))
#points(dt.11P.coding[ERE.pooled.cfu > 15 & SRE.pooled.cfu > 15, SRE.pooled.meanF],dt.11P.coding[ERE.pooled.cfu > 15 & SRE.pooled.cfu > 15, ERE.pooled.meanF],pch=20,cex=0.6,col="#00000010")
abline(v=dt.11P.coding["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
abline(h=7.104872,lwd=2,lty=2,col=rgb(100,69,155,maxColorValue=255)) #mean(ERE.pooled.egka)
points(mean(SRE.pooled.null.11P$mean.meanF),mean(ERE.pooled.null.11P$mean.meanF),pch="*",cex=1.5,col=rgb(211,211,211,maxColorValue=255))
arrows(dt.11P.coding[all.muts,SRE.pooled.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF]-1.645*dt.11P.coding[all.muts,ERE.SE.meanF], dt.11P.coding[all.muts,SRE.pooled.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF]+1.645*dt.11P.coding[all.muts,ERE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
arrows(dt.11P.coding[all.muts,SRE.pooled.meanF]-1.645*dt.11P.coding[all.muts,SRE.SE.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF], dt.11P.coding[all.muts,SRE.pooled.meanF]+1.645*dt.11P.coding[all.muts,SRE.SE.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF],length=0.01,angle=90,code=3,col="gray60")
points(dt.11P.coding[ERE.sp,SRE.pooled.meanF],dt.11P.coding[ERE.sp,ERE.pooled.meanF],col=rgb(100,69,155,maxColorValue=255),pch=19,cex=0.75)
points(dt.11P.coding[prom,SRE.pooled.meanF],dt.11P.coding[prom,ERE.pooled.meanF],col=rgb(111,204,255,maxColorValue=255),pch=19,cex=0.75)
points(dt.11P.coding[SRE.sp,SRE.pooled.meanF],dt.11P.coding[SRE.sp,ERE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.75)
#dev.off()
x <- "AGRV";points(dt.11P.coding[x,SRE.pooled.meanF],dt.11P.coding[x,ERE.pooled.meanF],pch=19,col="orange",cex=1.2)

#alternate: show just SRE activity on y-axis, number steps on x-axis (with some jitter?)
#pdf(file="9_quant-scale-landscape_out/EGKA-3-mut-paths_sig-increase-SRE_reaction-coordinate.pdf",5,4.5,useDingbats=F)
plot(NA,NA,xlim=c(-0.5,3.5),ylim=c(4.5, 10.2),xlab="Number of mutations",ylab="SRE mean fluorescence (a.u.)")
abline(h=dt.11P.coding["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
arrows(0,dt.11P.coding["EGKA",SRE.pooled.meanF]-1.645*dt.11P.coding["EGKA",SRE.SE.meanF],0,dt.11P.coding["EGKA",SRE.pooled.meanF]+1.645*dt.11P.coding["EGKA",SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(0,dt.11P.coding["EGKA",SRE.pooled.meanF],col=rgb(111,204,255,maxColorValue=255),pch=19,cex=0.75)
set.seed(105);x <- prom[prom %in% muts$mut1];jitter <- 1+runif(length(x),min=-0.3,max=0.3); arrows(jitter,dt.11P.coding[x,SRE.pooled.meanF]-1.645*dt.11P.coding[x,SRE.SE.meanF],jitter,dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[x,SRE.pooled.meanF],col=rgb(111,204,255,maxColorValue=255),pch=19,cex=0.75)
set.seed(105);x <- prom[prom %in% muts$mut2];jitter <- 2+runif(length(x),min=-0.3,max=0.3); arrows(jitter,dt.11P.coding[x,SRE.pooled.meanF]-1.645*dt.11P.coding[x,SRE.SE.meanF],jitter,dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[x,SRE.pooled.meanF],col=rgb(111,204,255,maxColorValue=255),pch=19,cex=0.75)
set.seed(105);x <- prom[prom %in% muts$mut3];jitter <- 3+runif(length(x),min=-0.3,max=0.3); arrows(jitter,dt.11P.coding[x,SRE.pooled.meanF]-1.645*dt.11P.coding[x,SRE.SE.meanF],jitter,dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[x,SRE.pooled.meanF],col=rgb(111,204,255,maxColorValue=255),pch=19,cex=0.75)
set.seed(104);x <- SRE.sp[SRE.sp %in% muts$mut1];jitter <- 1+runif(length(x),min=-0.3,max=0.3); arrows(jitter,dt.11P.coding[x,SRE.pooled.meanF]-1.645*dt.11P.coding[x,SRE.SE.meanF],jitter,dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[x,SRE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.75)
set.seed(102);x <- SRE.sp[SRE.sp %in% muts$mut2];jitter <- 2+runif(length(x),min=-0.3,max=0.3); arrows(jitter,dt.11P.coding[x,SRE.pooled.meanF]-1.645*dt.11P.coding[x,SRE.SE.meanF],jitter,dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[x,SRE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.75)
set.seed(105);x <- SRE.sp[SRE.sp %in% muts$mut3];jitter <- 3+runif(length(x),min=-0.3,max=0.3); arrows(jitter,dt.11P.coding[x,SRE.pooled.meanF]-1.645*dt.11P.coding[x,SRE.SE.meanF],jitter,dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[x,SRE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.75)
#dev.off()
x <- "AAKE";abline(h=dt.11P.coding[x,SRE.pooled.meanF])

#probability of occurrence of different SRE-specific outcomes, two ways: either each SRE-increasing step equally likely, or probability of a step proportional to the relative increase in SRE activity
#from EGKA (6 options)
#a <- b <- c <- d <- e <- f <- 1/6
a <- (dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])/((dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF]))
b <- (dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])/((dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF]))
c <- (dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])/((dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF]))
d <- (dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])/((dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF]))
e <- (dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])/((dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF]))
f <- (dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])/((dt.11P.coding["EAKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["QGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF])+(dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["EGKA",SRE.pooled.meanF]))

#from DGKA (3 options)
#g <- h <- i <- 1/3
g <- (dt.11P.coding["DAKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])/((dt.11P.coding["DAKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF]))
h <- (dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])/((dt.11P.coding["DAKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF]))
i <- (dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])/((dt.11P.coding["DAKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])+(dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["DGKA",SRE.pooled.meanF]))

#from DGKV (2 options)
#j <- k <- 1/2
j <- (dt.11P.coding["DAKV",SRE.pooled.meanF]-dt.11P.coding["DGKV",SRE.pooled.meanF])/((dt.11P.coding["DAKV",SRE.pooled.meanF]-dt.11P.coding["DGKV",SRE.pooled.meanF])+(dt.11P.coding["DGKM",SRE.pooled.meanF]-dt.11P.coding["DGKV",SRE.pooled.meanF]))
k <- (dt.11P.coding["DGKM",SRE.pooled.meanF]-dt.11P.coding["DGKV",SRE.pooled.meanF])/((dt.11P.coding["DAKV",SRE.pooled.meanF]-dt.11P.coding["DGKV",SRE.pooled.meanF])+(dt.11P.coding["DGKM",SRE.pooled.meanF]-dt.11P.coding["DGKV",SRE.pooled.meanF]))

#from AGKA (6 options)
#l <- m <- n <- o <- p <- q <- 1/6
l <- (dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])/((dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF]))
m <- (dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])/((dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF]))
n <- (dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])/((dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF]))
o <- (dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])/((dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF]))
p <- (dt.11P.coding["AAKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])/((dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF]))
q <- (dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])/((dt.11P.coding["VGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKE",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["AGKV",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["SGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF])+(dt.11P.coding["GGKA",SRE.pooled.meanF]-dt.11P.coding["AGKA",SRE.pooled.meanF]))

#from VGKA (1 option)
r <- 1

#from AGKV (6 options)
#s <- t <- u <- v <- w <- x <- 1/6
s <- (dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])/((dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF]))
t <- (dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])/((dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF]))
u <- (dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])/((dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF]))
v <- (dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])/((dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF]))
w <- (dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])/((dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF]))
x <- (dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])/((dt.11P.coding["DGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AAKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGKM",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["ASKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["AGRV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF])+(dt.11P.coding["GGKV",SRE.pooled.meanF]-dt.11P.coding["AGKV",SRE.pooled.meanF]))

#from AGKE (3 options)
#y <- z <- A <- 1/3
y <- (dt.11P.coding["SGKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])/((dt.11P.coding["SGKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])+(dt.11P.coding["AAKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])+(dt.11P.coding["AGKQ",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF]))
z <- (dt.11P.coding["AAKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])/((dt.11P.coding["SGKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])+(dt.11P.coding["AAKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])+(dt.11P.coding["AGKQ",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF]))
A <- (dt.11P.coding["AGKQ",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])/((dt.11P.coding["SGKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])+(dt.11P.coding["AAKE",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF])+(dt.11P.coding["AGKQ",SRE.pooled.meanF]-dt.11P.coding["AGKE",SRE.pooled.meanF]))

VAKA <- f*r + e*l*r
DAKV <- c*h*j + e*n*s*j
DGKM <- c*h*k + e*n*s*k
SGKE <- e*m*y
AAKE <- e*m*z
AGKQ <- e*m*A
AAKV <- e*n*t
AGKM <- e*n*u
ASKV <- e*n*v
AGRV <- e*n*w
GGKV <- e*n*x

#make another related figure, for Extended Data Fig. 7: the path illustrated from EGKA in the AncSR1 space, how does SRE meanF change along the trajectory?
pdf(file="9_quant-scale-landscape_out/SR1_EGKA-paths-to-SRE_reaction-coordinate.pdf",3.75,3.5,useDingbats=F)
plot(NA,NA,xlim=c(-0.5,5.5),ylim=c(5, 9),xlab="Number of mutations",ylab="SRE mean fluorescence (a.u.)")
abline(h=dt.11P.coding["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
x <- "EGKA";y <- 0; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(100, 69, 155,maxColorValue=255),pch=19,cex=0.75)
x <- "EGRA";y <- 1; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(100, 69, 155,maxColorValue=255),pch=19,cex=0.75)
x <- "EARA";y <- 2; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(100, 69, 155,maxColorValue=255),pch=19,cex=0.75)
x <- "DARA";y <- 3; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(111, 204, 255,maxColorValue=255),pch=19,cex=0.75)
x <- "DARV";y <- 4; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(111, 204, 255,maxColorValue=255),pch=19,cex=0.75)
x <- "DARI";y <- 5; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(5, 172, 72,maxColorValue=255),pch=19,cex=0.75)
x <- "HARV";y <- 5; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(5, 172, 72,maxColorValue=255),pch=19,cex=0.75)
x <- "DARL";y <- 5; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(5, 172, 72,maxColorValue=255),pch=19,cex=0.75)
x <- "DARM";y <- 5; arrows(y,dt.SR1.coding[x,SRE.pooled.meanF]-1.645*dt.SR1.coding[x,SRE.SE.meanF],y,dt.SR1.coding[x,SRE.pooled.meanF]+1.645*dt.SR1.coding[x,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(y,dt.SR1.coding[x,SRE.pooled.meanF],col=rgb(5, 172, 72,maxColorValue=255),pch=19,cex=0.75)
dev.off()

#are there single-mut neighbors of GSKV with significantly higher SRE mean fluorescence?
x <- "GSRV";as.character(dt.11P.coding[get.Hamming1.nt(x),][specificity=="SRE-specific" & SRE.pooled.meanF - 1.645*SRE.SE.meanF > dt.11P.coding[x,SRE.pooled.meanF]+1.645*dt.11P.coding[x,SRE.SE.meanF],AAseq])

#illustrate on a plot: y axis, SRE mean fluor, x axis, mut distance from GSKV
step0 <- "GSKV"

step1 <- as.character(dt.11P.coding[get.Hamming1.nt("GSKV")][specificity=="SRE-specific",AAseq])

step2 <- vector()
for(i in step1){
  step2 <- c(step2, as.character(dt.11P.coding[get.Hamming1.nt(i)][specificity=="SRE-specific",AAseq]))
}
step2 <- unique(step2); step2 <- step2[!(step2 %in% c(step1, step0))]

pdf(file="9_quant-scale-landscape_out/11P_GSKV-paths-to-addl-SRE_reaction-coordinate.pdf",2.75,3.5,useDingbats=F)
plot(NA,NA,xlim=c(-0.25,2.5),ylim=c(4.5, 10.2),xlab="Number of steps from GSKV",ylab="SRE mean fluorescence (a.u.)")
#abline(h=dt.11P.coding["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
arrows(0,dt.11P.coding["GSKV",SRE.pooled.meanF]-1.645*dt.11P.coding["GSKV",SRE.SE.meanF],0,dt.11P.coding["GSKV",SRE.pooled.meanF]+1.645*dt.11P.coding["GSKV",SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60")
points(0,dt.11P.coding["GSKV",SRE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.5)
set.seed(104);jitter <- 1+runif(length(step1),min=-0.2,max=0.2); arrows(jitter,dt.11P.coding[step1,SRE.pooled.meanF]-1.645*dt.11P.coding[step1,SRE.SE.meanF],jitter,dt.11P.coding[step1,SRE.pooled.meanF]+1.645*dt.11P.coding[step1,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[step1,SRE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.5)
set.seed(104);jitter <- 2+runif(length(step2),min=-0.2,max=0.2); arrows(jitter,dt.11P.coding[step2,SRE.pooled.meanF]-1.645*dt.11P.coding[step2,SRE.SE.meanF],jitter,dt.11P.coding[step2,SRE.pooled.meanF]+1.645*dt.11P.coding[step2,SRE.SE.meanF],length=0.01,angle=90,code=3,col="gray60");points(jitter,dt.11P.coding[step2,SRE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19,cex=0.5)
dev.off()

##############################################################################
#alternative: allow a step if it increases SRE activity OR decreases ERE activity (and maintains functionality); first, just with values, then add some accounting for SE
#fixget.3mut.paths.increasingSP.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P$AAseq),prom=as.character(promiscuous.11P$AAseq),SRE.sp=as.character(SRE.specific.11P$AAseq)){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  if(seq %in% ERE.sp){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(ERE.sp,prom,SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[seq,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[seq,ERE.pooled.meanF]),AAseq])
  }else if(seq %in% prom){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(prom,SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[seq,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[seq,ERE.pooled.meanF]),AAseq])
  }else{
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[seq,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[seq,ERE.pooled.meanF]),AAseq])
  }
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    if(i %in% ERE.sp){
      mutant2 <- c(mutant2, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)])
      mutant2 <- as.character(dt.11P.coding[AAseq %in% mutant2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[i,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[i,ERE.pooled.meanF]),AAseq])
    }else if(i %in% prom){
      mutant2 <- c(mutant2, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)])
      mutant2 <- as.character(dt.11P.coding[AAseq %in% mutant2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[i,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[i,ERE.pooled.meanF]),AAseq])
    }else{
      mutant2 <- c(mutant2, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)])
      mutant2 <- as.character(dt.11P.coding[AAseq %in% mutant2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[i,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[i,ERE.pooled.meanF]),AAseq])
    }
  }
  mutant2 <- unique(mutant2)
  mutant3 <- vector(mode="character")
  for(i in mutant2){
    if(i %in% ERE.sp){
      mutant3 <- c(mutant3, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)])
      mutant3 <- as.character(dt.11P.coding[AAseq %in% mutant3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[i,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[i,ERE.pooled.meanF]),AAseq])
    }else if(i %in% prom){
      mutant3 <- c(mutant3, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)])
      mutant3 <- as.character(dt.11P.coding[AAseq %in% mutant3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[i,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[i,ERE.pooled.meanF]),AAseq])
    }else{
      mutant3 <- c(mutant3, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)])
      mutant3 <- as.character(dt.11P.coding[AAseq %in% mutant3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & (SRE.pooled.meanF > dt.11P.coding[i,SRE.pooled.meanF] | ERE.pooled.meanF < dt.11P.coding[i,ERE.pooled.meanF]),AAseq])
    }
  }
  mutant3 <- unique(mutant3)
  return(list(mut0=seq,mut1=mutant1, mut2=mutant2, mut3=mutant3))
}

get.3mut.paths.sig.increasingSP.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P$AAseq),prom=as.character(promiscuous.11P$AAseq),SRE.sp=as.character(SRE.specific.11P$AAseq)){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  if(seq %in% ERE.sp){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(ERE.sp,prom,SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[seq,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[seq,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq])
  }else if(seq %in% prom){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(prom,SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[seq,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[seq,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq])
  }else{
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(SRE.sp)]
    mutant1 <- as.character(dt.11P.coding[AAseq %in% mutant1 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[seq,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[seq,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq])
  }
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    if(i %in% ERE.sp){
      muts2 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)]
      mutant2 <- c(mutant2,as.character(dt.11P.coding[AAseq %in% muts2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[i,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq]))
    }else if(i %in% prom){
      muts2 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)]
      mutant2 <- c(mutant2,as.character(dt.11P.coding[AAseq %in% muts2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[i,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq]))
    }else{
      muts2 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)]
      mutant2 <- c(mutant2,as.character(dt.11P.coding[AAseq %in% muts2 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[i,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq]))
    }
  }
  mutant2 <- unique(mutant2)
  mutant3 <- vector(mode="character")
  for(i in mutant2){
    if(i %in% ERE.sp){
      muts3 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)]
      mutant3 <- c(mutant3,as.character(dt.11P.coding[AAseq %in% muts3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[i,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq]))
    }else if(i %in% prom){
      muts3 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)]
      mutant3 <- c(mutant3,as.character(dt.11P.coding[AAseq %in% muts3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[i,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq]))
    }else{
      muts3 <- get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)]
      mutant3 <- c(mutant3,as.character(dt.11P.coding[AAseq %in% muts3 & !is.na(SRE.pooled.meanF) & !is.na(SRE.SE.meanF) & ((SRE.pooled.meanF - 1.645*SRE.SE.meanF) > dt.11P.coding[i,(SRE.pooled.meanF + 1.645*SRE.SE.meanF)] | (ERE.pooled.meanF + 1.645*ERE.SE.meanF) < dt.11P.coding[i,(ERE.pooled.meanF - 1.645*ERE.SE.meanF)]),AAseq]))
    }
  }
  mutant3 <- unique(mutant3)
  return(list(mut0=seq,mut1=mutant1, mut2=mutant2, mut3=mutant3))
}

muts <- get.3mut.paths.sig.increasingSP.11P("EGKA")
all.muts <- unlist(muts)
ERE.sp <- unlist(muts)[unlist(muts) %in% ERE.specific.11P$AAseq]
SRE.sp <- unlist(muts)[unlist(muts) %in% SRE.specific.11P$AAseq]
prom <- unlist(muts)[unlist(muts) %in% promiscuous.11P$AAseq]

pdf(file="9_quant-scale-landscape_out/EGKA-3-mut-paths_sig-increase-sp_meanF-plot.pdf",3,3.5,useDingbats=F)
plot(dt.11P.coding$SRE.pooled.meanF,dt.11P.coding$ERE.pooled.meanF,pch="",xlab="SRE meanF (a.u.), pooled reps",ylab="ERE meanF (a.u.), pooled reps",xlim=c(4.5,10.2),ylim=c(4.5,10.2))
#points(dt.11P.coding[ERE.pooled.cfu > 15 & SRE.pooled.cfu > 15, SRE.pooled.meanF],dt.11P.coding[ERE.pooled.cfu > 15 & SRE.pooled.cfu > 15, ERE.pooled.meanF],pch=20,cex=0.6,col="#00000010")
abline(v=dt.11P.coding["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
abline(h=7.104872,lwd=2,lty=2,col=rgb(100,69,155,maxColorValue=255)) #mean(ERE.pooled.egka)
points(mean(SRE.pooled.null.11P$mean.meanF),mean(ERE.pooled.null.11P$mean.meanF),pch="*",cex=1.5,col=rgb(211,211,211,maxColorValue=255))
arrows(dt.11P.coding[all.muts,SRE.pooled.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF]-1.645*dt.11P.coding[all.muts,ERE.SE.meanF], dt.11P.coding[all.muts,SRE.pooled.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF]+1.645*dt.11P.coding[all.muts,ERE.SE.meanF],length=0.025,angle=90,code=3,col="gray60")
arrows(dt.11P.coding[all.muts,SRE.pooled.meanF]-1.645*dt.11P.coding[all.muts,SRE.SE.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF], dt.11P.coding[all.muts,SRE.pooled.meanF]+1.645*dt.11P.coding[all.muts,SRE.SE.meanF],dt.11P.coding[all.muts,ERE.pooled.meanF],length=0.025,angle=90,code=3,col="gray60")
points(dt.11P.coding[ERE.sp,SRE.pooled.meanF],dt.11P.coding[ERE.sp,ERE.pooled.meanF],col=rgb(100,69,155,maxColorValue=255),pch=19)
points(dt.11P.coding[prom,SRE.pooled.meanF],dt.11P.coding[prom,ERE.pooled.meanF],col=rgb(111,204,255,maxColorValue=255),pch=19)
points(dt.11P.coding[SRE.sp,SRE.pooled.meanF],dt.11P.coding[SRE.sp,ERE.pooled.meanF],col=rgb(5,172,72,maxColorValue=255),pch=19)
dev.off()