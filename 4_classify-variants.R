#16 Jan 2017
#TNS
#script to classify variants as strong, weak, null on each of ERE and SRE in AncSR1 and AncSR1+11P backgrounds

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)

#read in tables with calculated meanF for each variant, output of calc-meanF.R scripts
dt.11P <- read.table(file="1_calc-meanF_11P_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.SR1 <- read.table(file="2_calc-meanF_SR1_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

#define function to get p value for test of superiority versus null
get.p.null <- function(cfu,mean,null.bins){
  if(is.na(mean)){
    return(as.numeric(NA))
  } else {
    for(i in 1:nrow(null.bins)){
      if(cfu >= null.bins[i,"range.cfu"][[1]][1] & cfu < null.bins[i,"range.cfu"][[1]][2]){
        return(length(which(null.bins[i,"list.means"][[1]]>=mean))/length(null.bins[i,"list.means"][[1]]))
      }
    }
  }
}

#load sampling distributions of null variants across depth of # cells, delimited in 3_classify-variants_define-distributions.R
load(file="./3_classify-variants_define-distributions_out/ERE.rep1.null.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.rep2.null.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.pooled.null.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep1.null.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep2.null.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.pooled.null.11P.Rda")

#calculate p-values
dt.11P[,ERE.rep1.p.null := get.p.null(ERE.rep1.cfu,ERE.rep1.meanF,ERE.rep1.null.11P),by=AAseq]
dt.11P[,ERE.rep2.p.null := get.p.null(ERE.rep2.cfu,ERE.rep2.meanF,ERE.rep2.null.11P),by=AAseq]
dt.11P[,ERE.pooled.p.null := get.p.null(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.null.11P),by=AAseq]
dt.11P[,SRE.rep1.p.null := get.p.null(SRE.rep1.cfu,SRE.rep1.meanF,SRE.rep1.null.11P),by=AAseq]
dt.11P[,SRE.rep2.p.null := get.p.null(SRE.rep2.cfu,SRE.rep2.meanF,SRE.rep2.null.11P),by=AAseq]
dt.11P[,SRE.pooled.p.null := get.p.null(SRE.pooled.cfu,SRE.pooled.meanF,SRE.pooled.null.11P),by=AAseq]

#load samplingdistributions of null variants across depth of # cells
load(file="./3_classify-variants_define-distributions_out/ERE.rep1.null.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.rep2.null.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.pooled.null.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep1.null.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep2.null.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.pooled.null.SR1.Rda")

#calculate p-values
dt.SR1[,ERE.rep1.p.null := get.p.null(ERE.rep1.cfu,ERE.rep1.meanF,ERE.rep1.null.SR1),by=AAseq]
dt.SR1[,ERE.rep2.p.null := get.p.null(ERE.rep2.cfu,ERE.rep2.meanF,ERE.rep2.null.SR1),by=AAseq]
dt.SR1[,ERE.pooled.p.null := get.p.null(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.null.SR1),by=AAseq]
dt.SR1[,SRE.rep1.p.null := get.p.null(SRE.rep1.cfu,SRE.rep1.meanF,SRE.rep1.null.SR1),by=AAseq]
dt.SR1[,SRE.rep2.p.null := get.p.null(SRE.rep2.cfu,SRE.rep2.meanF,SRE.rep2.null.SR1),by=AAseq]
dt.SR1[,SRE.pooled.p.null := get.p.null(SRE.pooled.cfu,SRE.pooled.meanF,SRE.pooled.null.SR1),by=AAseq]

####################################################################################################

#define function to get p value for test of non-inferiority versus evolutionary wildtype (with 20% margin of indistinguishability)
get.p.weak <- function(cfu,mean,pos.bins){
  if(is.na(mean)){
    return(as.numeric(NA))
  } else {
    for(i in 1:nrow(pos.bins)){
      if(cfu >= pos.bins[i,"range.cfu"][[1]][1] & cfu < pos.bins[i,"range.cfu"][[1]][2]){
        return(length(which(pos.bins[i,"list.means"][[1]]>=mean))/length(pos.bins[i,"list.means"][[1]]))
      }
    }
  } 
}

#load sampling distributions of positive reference versus depth of reads, delimited in 3_classify-variants_define-distributions.R
load(file="./3_classify-variants_define-distributions_out/ERE.rep1.pos.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.rep2.pos.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.pooled.pos.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep1.pos.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep2.pos.11P.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.pooled.pos.11P.Rda")

#calculate p-values
dt.11P[,ERE.rep1.p.weak := get.p.weak(ERE.rep1.cfu,ERE.rep1.meanF,ERE.rep1.pos.11P),by=AAseq]
dt.11P[,ERE.rep2.p.weak := get.p.weak(ERE.rep2.cfu,ERE.rep2.meanF,ERE.rep2.pos.11P),by=AAseq]
dt.11P[,ERE.pooled.p.weak := get.p.weak(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.pos.11P),by=AAseq]
dt.11P[,SRE.rep1.p.weak := get.p.weak(SRE.rep1.cfu,SRE.rep1.meanF,SRE.rep1.pos.11P),by=AAseq]
dt.11P[,SRE.rep2.p.weak := get.p.weak(SRE.rep2.cfu,SRE.rep2.meanF,SRE.rep2.pos.11P),by=AAseq]
dt.11P[,SRE.pooled.p.weak := get.p.weak(SRE.pooled.cfu,SRE.pooled.meanF,SRE.pooled.pos.11P),by=AAseq]

#load sampling distributions of positive reference versus depth of reads, delimited in 3_classify-variants_define-distributions.R
load(file="./3_classify-variants_define-distributions_out/ERE.rep1.pos.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.rep2.pos.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/ERE.pooled.pos.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep1.pos.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.rep2.pos.SR1.Rda")
load(file="./3_classify-variants_define-distributions_out/SRE.pooled.pos.SR1.Rda")

#calculate p-values
dt.SR1[,ERE.rep1.p.weak := get.p.weak(ERE.rep1.cfu,ERE.rep1.meanF,ERE.rep1.pos.SR1),by=AAseq]
dt.SR1[,ERE.rep2.p.weak := get.p.weak(ERE.rep2.cfu,ERE.rep2.meanF,ERE.rep2.pos.SR1),by=AAseq]
dt.SR1[,ERE.pooled.p.weak := get.p.weak(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.pos.SR1),by=AAseq]
dt.SR1[,SRE.rep1.p.weak := get.p.weak(SRE.rep1.cfu,SRE.rep1.meanF,SRE.rep1.pos.SR1),by=AAseq]
dt.SR1[,SRE.rep2.p.weak := get.p.weak(SRE.rep2.cfu,SRE.rep2.meanF,SRE.rep2.pos.SR1),by=AAseq]
dt.SR1[,SRE.pooled.p.weak := get.p.weak(SRE.pooled.cfu,SRE.pooled.meanF,SRE.pooled.pos.SR1),by=AAseq]

####################################################################################################################

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#determine p-value cutoff for p.null (alt hyp: weak or strong) that controls FDR to be 5% false discoveries, as judged by false discovery of stop codon containing variants
#ERE rep1, 11P
p.values.coding <- c(dt.11P.coding$ERE.rep1.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.11P.stop$ERE.rep1.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #848
p.null.ERE.rep1.cutoff.11P <- p.values.coding[i] # significant p.value < 0.0008176615

#ERE rep2, 11P
p.values.coding <- c(dt.11P.coding$ERE.rep2.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.11P.stop$ERE.rep2.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #1334
p.null.ERE.rep2.cutoff.11P <- p.values.coding[i] # significant p.value < 0.0007874016

#ERE pooled, 11P
p.values.coding <- c(dt.11P.coding$ERE.pooled.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.11P.stop$ERE.pooled.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #1263
p.null.ERE.pooled.cutoff.11P <- p.values.coding[i] # significant p.value < 0.0007513148

#SRE rep1, 11P
p.values.coding <- c(dt.11P.coding$SRE.rep1.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.11P.stop$SRE.rep1.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #4274
p.null.SRE.rep1.cutoff.11P <- p.values.coding[i] # significant p.value < 0.001605136

#SRE rep2, 11P
p.values.coding <- c(dt.11P.coding$SRE.rep2.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.11P.stop$SRE.rep2.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #2253
p.null.SRE.rep2.cutoff.11P <- p.values.coding[i] # significant p.value < 0.000770416

#SRE pooled, 11P
p.values.coding <- c(dt.11P.coding$SRE.pooled.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.11P.stop$SRE.pooled.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #3352
p.null.SRE.pooled.cutoff.11P <- p.values.coding[i] # significant p.value < 0.00148258


dt.11P[,ERE.rep1.class := as.character(NA)];dt.11P[ERE.rep1.p.null<p.null.ERE.rep1.cutoff.11P,ERE.rep1.class:="weak"];dt.11P[ERE.rep1.p.null>=p.null.ERE.rep1.cutoff.11P,ERE.rep1.class:="null"]
dt.11P[,ERE.rep2.class := as.character(NA)];dt.11P[ERE.rep2.p.null<p.null.ERE.rep2.cutoff.11P,ERE.rep2.class:="weak"];dt.11P[ERE.rep2.p.null>=p.null.ERE.rep2.cutoff.11P,ERE.rep2.class:="null"]
dt.11P[,ERE.pooled.class := as.character(NA)];dt.11P[ERE.pooled.p.null<p.null.ERE.pooled.cutoff.11P,ERE.pooled.class:="weak"];dt.11P[ERE.pooled.p.null>=p.null.ERE.pooled.cutoff.11P,ERE.pooled.class:="null"]
dt.11P[,SRE.rep1.class := as.character(NA)];dt.11P[SRE.rep1.p.null<p.null.SRE.rep1.cutoff.11P,SRE.rep1.class:="weak"];dt.11P[SRE.rep1.p.null>=p.null.SRE.rep1.cutoff.11P,SRE.rep1.class:="null"]
dt.11P[,SRE.rep2.class := as.character(NA)];dt.11P[SRE.rep2.p.null<p.null.SRE.rep2.cutoff.11P,SRE.rep2.class:="weak"];dt.11P[SRE.rep2.p.null>=p.null.SRE.rep2.cutoff.11P,SRE.rep2.class:="null"]
dt.11P[,SRE.pooled.class := as.character(NA)];dt.11P[SRE.pooled.p.null<p.null.SRE.pooled.cutoff.11P,SRE.pooled.class:="weak"];dt.11P[SRE.pooled.p.null>=p.null.SRE.pooled.cutoff.11P,SRE.pooled.class:="null"]

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

#ERE rep1, SR1
p.values.coding <- c(dt.SR1.coding$ERE.rep1.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.SR1.stop$ERE.rep1.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #491
p.null.ERE.rep1.cutoff.SR1 <- p.values.coding[i] # significant p.value < 0.0007739938

#ERE rep2, SR1
p.values.coding <- c(dt.SR1.coding$ERE.rep2.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.SR1.stop$ERE.rep2.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #449
p.null.ERE.rep2.cutoff.SR1 <- p.values.coding[i] # significant p.value < 0.000798722

#ERE pooled, SR1
p.values.coding <- c(dt.SR1.coding$ERE.pooled.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.SR1.stop$ERE.pooled.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #597
p.null.ERE.pooled.cutoff.SR1 <- p.values.coding[i] # significant p.value < 0.0007412898

#SRE rep1, SR1
p.values.coding <- c(dt.SR1.coding$SRE.rep1.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.SR1.stop$SRE.rep1.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #294
p.null.SRE.rep1.cutoff.SR1 <- p.values.coding[i] # significant p.value < 0.0007342144

#SRE rep2, SR1
p.values.coding <- c(dt.SR1.coding$SRE.rep2.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.SR1.stop$SRE.rep2.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #294
p.null.SRE.rep2.cutoff.SR1 <- p.values.coding[i] # significant p.value < 0.0007535795

#SRE pooled, SR1
p.values.coding <- c(dt.SR1.coding$SRE.pooled.p.null); p.values.coding <- p.values.coding[!is.na(p.values.coding)]; p.values.coding <- sort(p.values.coding)
p.values.stop <- c(dt.SR1.stop$SRE.pooled.p.null); p.values.stop <- p.values.stop[!is.na(p.values.stop)]; p.values.stop <- sort(p.values.stop)
i <- 1
while(length(which(p.values.stop <= p.values.coding[i]))*(length(p.values.coding)/length(p.values.stop))/length(which(p.values.coding <= p.values.coding[i])) <= 0.05){
  i <- i+1
}
i #370
p.null.SRE.pooled.cutoff.SR1 <- p.values.coding[i] # significant p.value < 0.0007283321

dt.SR1[,ERE.rep1.class := as.character(NA)];dt.SR1[ERE.rep1.p.null<p.null.ERE.rep1.cutoff.SR1,ERE.rep1.class:="weak"];dt.SR1[ERE.rep1.p.null>=p.null.ERE.rep1.cutoff.SR1,ERE.rep1.class:="null"]
dt.SR1[,ERE.rep2.class := as.character(NA)];dt.SR1[ERE.rep2.p.null<p.null.ERE.rep2.cutoff.SR1,ERE.rep2.class:="weak"];dt.SR1[ERE.rep2.p.null>=p.null.ERE.rep2.cutoff.SR1,ERE.rep2.class:="null"]
dt.SR1[,ERE.pooled.class := as.character(NA)];dt.SR1[ERE.pooled.p.null<p.null.ERE.pooled.cutoff.SR1,ERE.pooled.class:="weak"];dt.SR1[ERE.pooled.p.null>=p.null.ERE.pooled.cutoff.SR1,ERE.pooled.class:="null"]
dt.SR1[,SRE.rep1.class := as.character(NA)];dt.SR1[SRE.rep1.p.null<p.null.SRE.rep1.cutoff.SR1,SRE.rep1.class:="weak"];dt.SR1[SRE.rep1.p.null>=p.null.SRE.rep1.cutoff.SR1,SRE.rep1.class:="null"]
dt.SR1[,SRE.rep2.class := as.character(NA)];dt.SR1[SRE.rep2.p.null<p.null.SRE.rep2.cutoff.SR1,SRE.rep2.class:="weak"];dt.SR1[SRE.rep2.p.null>=p.null.SRE.rep2.cutoff.SR1,SRE.rep2.class:="null"]
dt.SR1[,SRE.pooled.class := as.character(NA)];dt.SR1[SRE.pooled.p.null<p.null.SRE.pooled.cutoff.SR1,SRE.pooled.class:="weak"];dt.SR1[SRE.pooled.p.null>=p.null.SRE.pooled.cutoff.SR1,SRE.pooled.class:="null"]

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#################################################################################################

#determine p-value cutoff for hypothesis: weak (alt hypothesis: strong binder) that controls FDR to 5%
#ERE rep1, 11P
p.values.coding <- sort(c(dt.11P.coding[ERE.rep1.class=="weak",ERE.rep1.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 160
p.weak.ERE.rep1.cutoff.11P <- p.values.coding[i] # 0.0096

#ERE rep2, 11P
p.values.coding <- sort(c(dt.11P.coding[ERE.rep2.class=="weak",ERE.rep2.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 217
p.weak.ERE.rep2.cutoff.11P <- p.values.coding[i] # 0.0084

#ERE pooled, 11P
p.values.coding <- sort(c(dt.11P.coding[ERE.pooled.class=="weak",ERE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 299
p.weak.ERE.pooled.cutoff.11P <- p.values.coding[i] # 0.0126

#SRE rep1, 11P
p.values.coding <- sort(c(dt.11P.coding[SRE.rep1.class=="weak",SRE.rep1.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 1067
p.weak.SRE.rep1.cutoff.11P <- p.values.coding[i] # 0.0125

#SRE rep2, 11P
p.values.coding <- sort(c(dt.11P.coding[SRE.rep2.class=="weak",SRE.rep2.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 908
p.weak.SRE.rep2.cutoff.11P <- p.values.coding[i] # 0.0217

#SRE pooled, 11P
p.values.coding <- sort(c(dt.11P.coding[SRE.pooled.class=="weak",SRE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 1085
p.weak.SRE.pooled.cutoff.11P <- p.values.coding[i] # 0.0173

dt.11P[ERE.rep1.p.weak<p.weak.ERE.rep1.cutoff.11P & ERE.rep1.class=="weak",ERE.rep1.class:="strong"]
dt.11P[ERE.rep2.p.weak<p.weak.ERE.rep2.cutoff.11P & ERE.rep2.class=="weak",ERE.rep2.class:="strong"]
dt.11P[ERE.pooled.p.weak<p.weak.ERE.pooled.cutoff.11P & ERE.pooled.class=="weak",ERE.pooled.class:="strong"]

dt.11P[SRE.rep1.p.weak<p.weak.SRE.rep1.cutoff.11P & SRE.rep1.class=="weak",SRE.rep1.class:="strong"]
dt.11P[SRE.rep2.p.weak<p.weak.SRE.rep2.cutoff.11P & SRE.rep2.class=="weak",SRE.rep2.class:="strong"]
dt.11P[SRE.pooled.p.weak<p.weak.SRE.pooled.cutoff.11P & SRE.pooled.class=="weak",SRE.pooled.class:="strong"]

write.table(dt.11P, file="4_classify-variants_out/dt_output_11P.csv",col.names=T,row.names=F,quote=F,sep=",")

#ERE rep1, SR1
p.values.coding <- sort(c(dt.SR1.coding[ERE.rep1.class=="weak",ERE.rep1.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 58
p.weak.ERE.rep1.cutoff.SR1 <- p.values.coding[i] # 0.0067

#ERE rep2, SR1
p.values.coding <- sort(c(dt.SR1.coding[ERE.rep2.class=="weak",ERE.rep2.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 33
p.weak.ERE.rep2.cutoff.SR1 <- p.values.coding[i] # .0037

#ERE pooled, SR1
p.values.coding <- sort(c(dt.SR1.coding[ERE.pooled.class=="weak",ERE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 65
p.weak.ERE.pooled.cutoff.SR1 <- p.values.coding[i] # .0067

#SRE rep1, SR1
p.values.coding <- sort(c(dt.SR1.coding[SRE.rep1.class=="weak",SRE.rep1.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 68
p.weak.SRE.rep1.cutoff.SR1 <- p.values.coding[i] # 0.0156

#SRE rep2, SR1
p.values.coding <- sort(c(dt.SR1.coding[SRE.rep2.class=="weak",SRE.rep2.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 49
p.weak.SRE.rep2.cutoff.SR1 <- p.values.coding[i] # .009

#SRE pooled, SR1
p.values.coding <- sort(c(dt.SR1.coding[SRE.pooled.class=="weak",SRE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i # 70
p.weak.SRE.pooled.cutoff.SR1 <- p.values.coding[i] # .0306

dt.SR1[ERE.rep1.p.weak<p.weak.ERE.rep1.cutoff.SR1 & ERE.rep1.class=="weak",ERE.rep1.class:="strong"]
dt.SR1[ERE.rep2.p.weak<p.weak.ERE.rep2.cutoff.SR1 & ERE.rep2.class=="weak",ERE.rep2.class:="strong"]
dt.SR1[ERE.pooled.p.weak<p.weak.ERE.pooled.cutoff.SR1 & ERE.pooled.class=="weak",ERE.pooled.class:="strong"]

dt.SR1[SRE.rep1.p.weak<p.weak.SRE.rep1.cutoff.SR1 & SRE.rep1.class=="weak",SRE.rep1.class:="strong"]
dt.SR1[SRE.rep2.p.weak<p.weak.SRE.rep2.cutoff.SR1 & SRE.rep2.class=="weak",SRE.rep2.class:="strong"]
dt.SR1[SRE.pooled.p.weak<p.weak.SRE.pooled.cutoff.SR1 & SRE.pooled.class=="weak",SRE.pooled.class:="strong"]

write.table(dt.SR1, file="4_classify-variants_out/dt_output_SR1.csv",col.names=T,row.names=F,quote=F,sep=",")