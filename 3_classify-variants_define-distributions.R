#16 Jan 2017
#TNS
#generate null distributions for determination of p.weak and p.strong

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)
library(Hmisc)
library(fitdistrplus)

dt.11P <- read.table(file="1_calc-meanF_11P_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

dt.SR1 <- read.table(file="2_calc-meanF_SR1_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#to determine probability that a variant is null, compare to distribution of stop-codon-containing variants of equal depth of coverage
#ERE rep1, 11P
ERE.rep1.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(ERE.rep1.meanF),ERE.rep1.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(ERE.rep1.null.11P)){
  ERE.rep1.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[ERE.rep1.cfu>=ERE.rep1.null.11P$range.cfu[i][[1]][[1]] & ERE.rep1.cfu<ERE.rep1.null.11P$range.cfu[i][[1]][[2]],]
  ERE.rep1.null.11P$mean.meanF[i] <- mean(data$ERE.rep1.meanF,na.rm=T)
  ERE.rep1.null.11P$sd.meanF[i] <- sd(data$ERE.rep1.meanF,na.rm=T)
  ERE.rep1.null.11P$median.cfu[i] <- median(data$ERE.rep1.cfu,na.rm=T)
  ERE.rep1.null.11P$max.meanF[i] <- max(data$ERE.rep1.meanF,na.rm=T)
  ERE.rep1.null.11P$list.means[i] <- list(data[!is.na(ERE.rep1.meanF),ERE.rep1.meanF])
}
save(ERE.rep1.null.11P, file="3_classify-variants_define-distributions_out/ERE.rep1.null.11P.Rda")

#ERE rep2, 11P
ERE.rep2.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(ERE.rep2.meanF),ERE.rep2.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(ERE.rep2.null.11P)){
  ERE.rep2.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[ERE.rep2.cfu>=ERE.rep2.null.11P$range.cfu[i][[1]][[1]] & ERE.rep2.cfu<ERE.rep2.null.11P$range.cfu[i][[1]][[2]],]
  ERE.rep2.null.11P$mean.meanF[i] <- mean(data$ERE.rep2.meanF,na.rm=T)
  ERE.rep2.null.11P$sd.meanF[i] <- sd(data$ERE.rep2.meanF,na.rm=T)
  ERE.rep2.null.11P$median.cfu[i] <- median(data$ERE.rep2.cfu,na.rm=T)
  ERE.rep2.null.11P$max.meanF[i] <- max(data$ERE.rep2.meanF,na.rm=T)
  ERE.rep2.null.11P$list.means[i] <- list(data[!is.na(ERE.rep2.meanF),ERE.rep2.meanF])
}
save(ERE.rep2.null.11P, file="3_classify-variants_define-distributions_out/ERE.rep2.null.11P.Rda")

#ERE pooled, 11P
ERE.pooled.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(ERE.pooled.meanF),ERE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(ERE.pooled.null.11P)){
  ERE.pooled.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[ERE.pooled.cfu>=ERE.pooled.null.11P$range.cfu[i][[1]][[1]] & ERE.pooled.cfu<ERE.pooled.null.11P$range.cfu[i][[1]][[2]],]
  ERE.pooled.null.11P$mean.meanF[i] <- mean(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.11P$sd.meanF[i] <- sd(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.11P$median.cfu[i] <- median(data$ERE.pooled.cfu,na.rm=T)
  ERE.pooled.null.11P$max.meanF[i] <- max(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.11P$list.means[i] <- list(data[!is.na(ERE.pooled.meanF),ERE.pooled.meanF])
}
save(ERE.pooled.null.11P, file="3_classify-variants_define-distributions_out/ERE.pooled.null.11P.Rda")

#SRE rep1, 11P
SRE.rep1.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(SRE.rep1.meanF),SRE.rep1.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(SRE.rep1.null.11P)){
  SRE.rep1.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[SRE.rep1.cfu>=SRE.rep1.null.11P$range.cfu[i][[1]][[1]] & SRE.rep1.cfu<SRE.rep1.null.11P$range.cfu[i][[1]][[2]],]
  SRE.rep1.null.11P$mean.meanF[i] <- mean(data$SRE.rep1.meanF,na.rm=T)
  SRE.rep1.null.11P$sd.meanF[i] <- sd(data$SRE.rep1.meanF,na.rm=T)
  SRE.rep1.null.11P$median.cfu[i] <- median(data$SRE.rep1.cfu,na.rm=T)
  SRE.rep1.null.11P$max.meanF[i] <- max(data$SRE.rep1.meanF,na.rm=T)
  SRE.rep1.null.11P$list.means[i] <- list(data[!is.na(SRE.rep1.meanF),SRE.rep1.meanF])
}
save(SRE.rep1.null.11P, file="3_classify-variants_define-distributions_out/SRE.rep1.null.11P.Rda")

#SRE rep2, 11P
SRE.rep2.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(SRE.rep2.meanF),SRE.rep2.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(SRE.rep2.null.11P)){
  SRE.rep2.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[SRE.rep2.cfu>=SRE.rep2.null.11P$range.cfu[i][[1]][[1]] & SRE.rep2.cfu<SRE.rep2.null.11P$range.cfu[i][[1]][[2]],]
  SRE.rep2.null.11P$mean.meanF[i] <- mean(data$SRE.rep2.meanF,na.rm=T)
  SRE.rep2.null.11P$sd.meanF[i] <- sd(data$SRE.rep2.meanF,na.rm=T)
  SRE.rep2.null.11P$median.cfu[i] <- median(data$SRE.rep2.cfu,na.rm=T)
  SRE.rep2.null.11P$max.meanF[i] <- max(data$SRE.rep2.meanF,na.rm=T)
  SRE.rep2.null.11P$list.means[i] <- list(data[!is.na(SRE.rep2.meanF),SRE.rep2.meanF])
}
save(SRE.rep2.null.11P, file="3_classify-variants_define-distributions_out/SRE.rep2.null.11P.Rda")

#SRE pooled, 11P
SRE.pooled.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(SRE.pooled.meanF),SRE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(SRE.pooled.null.11P)){
  SRE.pooled.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[SRE.pooled.cfu>=SRE.pooled.null.11P$range.cfu[i][[1]][[1]] & SRE.pooled.cfu<SRE.pooled.null.11P$range.cfu[i][[1]][[2]],]
  SRE.pooled.null.11P$mean.meanF[i] <- mean(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.11P$sd.meanF[i] <- sd(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.11P$median.cfu[i] <- median(data$SRE.pooled.cfu,na.rm=T)
  SRE.pooled.null.11P$max.meanF[i] <- max(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.11P$list.means[i] <- list(data[!is.na(SRE.pooled.meanF),SRE.pooled.meanF])
}
save(SRE.pooled.null.11P, file="3_classify-variants_define-distributions_out/SRE.pooled.null.11P.Rda")

#ERE rep1, SR1
ERE.rep1.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(ERE.rep1.meanF),ERE.rep1.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(ERE.rep1.null.SR1)){
  ERE.rep1.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[ERE.rep1.cfu>=ERE.rep1.null.SR1$range.cfu[i][[1]][[1]] & ERE.rep1.cfu<ERE.rep1.null.SR1$range.cfu[i][[1]][[2]],]
  ERE.rep1.null.SR1$mean.meanF[i] <- mean(data$ERE.rep1.meanF,na.rm=T)
  ERE.rep1.null.SR1$sd.meanF[i] <- sd(data$ERE.rep1.meanF,na.rm=T)
  ERE.rep1.null.SR1$median.cfu[i] <- median(data$ERE.rep1.cfu,na.rm=T)
  ERE.rep1.null.SR1$max.meanF[i] <- max(data$ERE.rep1.meanF,na.rm=T)
  ERE.rep1.null.SR1$list.means[i] <- list(data[!is.na(ERE.rep1.meanF),ERE.rep1.meanF])
}
save(ERE.rep1.null.SR1, file="3_classify-variants_define-distributions_out/ERE.rep1.null.SR1.Rda")

#ERE rep2, SR1
ERE.rep2.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(ERE.rep2.meanF),ERE.rep2.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(ERE.rep2.null.SR1)){
  ERE.rep2.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[ERE.rep2.cfu>=ERE.rep2.null.SR1$range.cfu[i][[1]][[1]] & ERE.rep2.cfu<ERE.rep2.null.SR1$range.cfu[i][[1]][[2]],]
  ERE.rep2.null.SR1$mean.meanF[i] <- mean(data$ERE.rep2.meanF,na.rm=T)
  ERE.rep2.null.SR1$sd.meanF[i] <- sd(data$ERE.rep2.meanF,na.rm=T)
  ERE.rep2.null.SR1$median.cfu[i] <- median(data$ERE.rep2.cfu,na.rm=T)
  ERE.rep2.null.SR1$max.meanF[i] <- max(data$ERE.rep2.meanF,na.rm=T)
  ERE.rep2.null.SR1$list.means[i] <- list(data[!is.na(ERE.rep2.meanF),ERE.rep2.meanF])
}
save(ERE.rep2.null.SR1, file="3_classify-variants_define-distributions_out/ERE.rep2.null.SR1.Rda")

#ERE pooled, SR1
ERE.pooled.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(ERE.pooled.meanF),ERE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(ERE.pooled.null.SR1)){
  ERE.pooled.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[ERE.pooled.cfu>=ERE.pooled.null.SR1$range.cfu[i][[1]][[1]] & ERE.pooled.cfu<ERE.pooled.null.SR1$range.cfu[i][[1]][[2]],]
  ERE.pooled.null.SR1$mean.meanF[i] <- mean(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$sd.meanF[i] <- sd(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$median.cfu[i] <- median(data$ERE.pooled.cfu,na.rm=T)
  ERE.pooled.null.SR1$max.meanF[i] <- max(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$list.means[i] <- list(data[!is.na(ERE.pooled.meanF),ERE.pooled.meanF])
}
save(ERE.pooled.null.SR1, file="3_classify-variants_define-distributions_out/ERE.pooled.null.SR1.Rda")

#SRE rep1, SR1
SRE.rep1.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(SRE.rep1.meanF),SRE.rep1.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(SRE.rep1.null.SR1)){
  SRE.rep1.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[SRE.rep1.cfu>=SRE.rep1.null.SR1$range.cfu[i][[1]][[1]] & SRE.rep1.cfu<SRE.rep1.null.SR1$range.cfu[i][[1]][[2]],]
  SRE.rep1.null.SR1$mean.meanF[i] <- mean(data$SRE.rep1.meanF,na.rm=T)
  SRE.rep1.null.SR1$sd.meanF[i] <- sd(data$SRE.rep1.meanF,na.rm=T)
  SRE.rep1.null.SR1$median.cfu[i] <- median(data$SRE.rep1.cfu,na.rm=T)
  SRE.rep1.null.SR1$max.meanF[i] <- max(data$SRE.rep1.meanF,na.rm=T)
  SRE.rep1.null.SR1$list.means[i] <- list(data[!is.na(SRE.rep1.meanF),SRE.rep1.meanF])
}
save(SRE.rep1.null.SR1, file="3_classify-variants_define-distributions_out/SRE.rep1.null.SR1.Rda")

#SRE rep2, SR1
SRE.rep2.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(SRE.rep2.meanF),SRE.rep2.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(SRE.rep2.null.SR1)){
  SRE.rep2.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[SRE.rep2.cfu>=SRE.rep2.null.SR1$range.cfu[i][[1]][[1]] & SRE.rep2.cfu<SRE.rep2.null.SR1$range.cfu[i][[1]][[2]],]
  SRE.rep2.null.SR1$mean.meanF[i] <- mean(data$SRE.rep2.meanF,na.rm=T)
  SRE.rep2.null.SR1$sd.meanF[i] <- sd(data$SRE.rep2.meanF,na.rm=T)
  SRE.rep2.null.SR1$median.cfu[i] <- median(data$SRE.rep2.cfu,na.rm=T)
  SRE.rep2.null.SR1$max.meanF[i] <- max(data$SRE.rep2.meanF,na.rm=T)
  SRE.rep2.null.SR1$list.means[i] <- list(data[!is.na(SRE.rep2.meanF),SRE.rep2.meanF])
}
save(SRE.rep2.null.SR1, file="3_classify-variants_define-distributions_out/SRE.rep2.null.SR1.Rda")

#SRE pooled, SR1
SRE.pooled.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(SRE.pooled.meanF),SRE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
#for each bin, make data subset that has # cfus within break range (lower boundary inclusive)
for(i in 1:nrow(SRE.pooled.null.SR1)){
  SRE.pooled.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[SRE.pooled.cfu>=SRE.pooled.null.SR1$range.cfu[i][[1]][[1]] & SRE.pooled.cfu<SRE.pooled.null.SR1$range.cfu[i][[1]][[2]],]
  SRE.pooled.null.SR1$mean.meanF[i] <- mean(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.SR1$sd.meanF[i] <- sd(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.SR1$median.cfu[i] <- median(data$SRE.pooled.cfu,na.rm=T)
  SRE.pooled.null.SR1$max.meanF[i] <- max(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.SR1$list.means[i] <- list(data[!is.na(SRE.pooled.meanF),SRE.pooled.meanF])
}
save(SRE.pooled.null.SR1, file="3_classify-variants_define-distributions_out/SRE.pooled.null.SR1.Rda")

#determine meanF of stop variants
ERE.rep1.stop.meanF.11P <- mean(dt.11P.stop$ERE.rep1.meanF,na.rm=T)
ERE.rep2.stop.meanF.11P <- mean(dt.11P.stop$ERE.rep2.meanF,na.rm=T)
ERE.pooled.stop.meanF.11P <- mean(dt.11P.stop$ERE.pooled.meanF,na.rm=T)
SRE.rep1.stop.meanF.11P <- mean(dt.11P.stop$SRE.rep1.meanF,na.rm=T)
SRE.rep2.stop.meanF.11P <- mean(dt.11P.stop$SRE.rep2.meanF,na.rm=T)
SRE.pooled.stop.meanF.11P <- mean(dt.11P.stop$SRE.pooled.meanF,na.rm=T)
ERE.rep1.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.rep1.meanF,na.rm=T)
ERE.rep2.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.rep2.meanF,na.rm=T)
ERE.pooled.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.pooled.meanF,na.rm=T)
SRE.rep1.stop.meanF.SR1 <- mean(dt.SR1.stop$SRE.rep1.meanF,na.rm=T)
SRE.rep2.stop.meanF.SR1 <- mean(dt.SR1.stop$SRE.rep2.meanF,na.rm=T)
SRE.pooled.stop.meanF.SR1 <- mean(dt.SR1.stop$SRE.pooled.meanF,na.rm=T)

#read in parallel control flow cytometry values for "wt" controls
ERE.rep1.egka.11P <- log(read.csv(file="3_classify-variants_define-distributions_in/l6_egka_ctrl_distribution.csv",header=F)$V1[1:7948])
ERE.rep2.egka.11P <- (log(read.csv(file="3_classify-variants_define-distributions_in/l13_egka_ctrl_distribution.csv",header=F)$V1)-0.01014)/1.00035
SRE.rep1.GSKV.11P <- (log(read.csv(file="3_classify-variants_define-distributions_in/l5_GSKV_ctrl_distribution.csv",header=F)$V1)-0.46994)/0.95213
SRE.rep2.GSKV.11P <- log(read.csv(file="3_classify-variants_define-distributions_in/l8_GSKV_ctrl_distribution.csv",header=F)$V1[1:length(SRE.rep1.GSKV.11P)])
set.seed(103)
ERE.pooled.egka.11P <- sample(c(ERE.rep1.egka.11P,ERE.rep2.egka.11P))
set.seed(201)
SRE.pooled.GSKV.11P <- sample(c(SRE.rep1.GSKV.11P,SRE.rep2.GSKV.11P))

ERE.rep1.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l10_egka_ctrl_distribution.csv",header=F)$V1)+0.59618)/1.01805; ERE.rep1.egka.SR1 <- ERE.rep1.egka.SR1[is.finite(ERE.rep1.egka.SR1)]
ERE.rep2.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l12_egka_ctrl_distribution.csv",header=F)$V1[1:length(ERE.rep1.egka.SR1)])-0.5260)/0.9633
SRE.rep1.GSKV.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l9_GSKV_ctrl_distribution.csv",header=F)$V1)+0.82284)/1.17925
SRE.rep2.GSKV.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l11_GSKV_ctrl_distribution.csv",header=F)$V1[1:length(SRE.rep1.GSKV.SR1)])+0.57030)/1.11911
set.seed(104)
ERE.pooled.egka.SR1 <- sample(c(ERE.rep1.egka.SR1,ERE.rep2.egka.SR1))
set.seed(202)
SRE.pooled.GSKV.SR1 <- sample(c(SRE.rep1.GSKV.SR1,SRE.rep2.GSKV.SR1))

#shift controls over to 80th percentile between null and wildtype, to give most extreme instantiation of the null hypothesis that variant is weak that I want to test
ERE.rep1.egka.p80.11P <- ERE.rep1.egka.11P - (mean(ERE.rep1.egka.11P)-ERE.rep1.stop.meanF.11P)*0.2
ERE.rep2.egka.p80.11P <- ERE.rep2.egka.11P - (mean(ERE.rep2.egka.11P)-ERE.rep2.stop.meanF.11P)*0.2
ERE.pooled.egka.p80.11P <- ERE.pooled.egka.11P - (mean(ERE.pooled.egka.11P)-ERE.pooled.stop.meanF.11P)*0.2
SRE.rep1.GSKV.p80.11P <- SRE.rep1.GSKV.11P + (dt.11P["GSKV",SRE.rep1.meanF]-mean(SRE.rep1.GSKV.11P)) - (dt.11P["GSKV",SRE.rep1.meanF]-SRE.rep1.stop.meanF.11P)*0.2
SRE.rep2.GSKV.p80.11P <- SRE.rep2.GSKV.11P + (dt.11P["GSKV",SRE.rep2.meanF]-mean(SRE.rep2.GSKV.11P)) - (dt.11P["GSKV",SRE.rep2.meanF]-SRE.rep2.stop.meanF.11P)*0.2
SRE.pooled.GSKV.p80.11P <- SRE.pooled.GSKV.11P + (dt.11P["GSKV",SRE.pooled.meanF]-mean(SRE.pooled.GSKV.11P)) - (dt.11P["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.11P)*0.2

ERE.rep1.egka.p80.SR1 <- ERE.rep1.egka.SR1 - (mean(ERE.rep1.egka.SR1)-ERE.rep1.stop.meanF.SR1)*0.2
ERE.rep2.egka.p80.SR1 <- ERE.rep2.egka.SR1 - (mean(ERE.rep2.egka.SR1)-ERE.rep2.stop.meanF.SR1)*0.2
ERE.pooled.egka.p80.SR1 <- ERE.pooled.egka.SR1 - (mean(ERE.pooled.egka.SR1)-ERE.pooled.stop.meanF.SR1)*0.2
SRE.rep1.GSKV.p80.SR1 <- SRE.rep1.GSKV.SR1 + (dt.11P["GSKV",SRE.rep1.meanF]-mean(SRE.rep1.GSKV.SR1)) - (dt.11P["GSKV",SRE.rep1.meanF]-SRE.rep1.stop.meanF.SR1)*0.2
SRE.rep2.GSKV.p80.SR1 <- SRE.rep2.GSKV.SR1 + (dt.11P["GSKV",SRE.rep2.meanF]-mean(SRE.rep2.GSKV.SR1)) - (dt.11P["GSKV",SRE.rep2.meanF]-SRE.rep2.stop.meanF.SR1)*0.2
SRE.pooled.GSKV.p80.SR1 <- SRE.pooled.GSKV.SR1 + (dt.11P["GSKV",SRE.pooled.meanF]-mean(SRE.pooled.GSKV.SR1)) - (dt.11P["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.SR1)*0.2

#give sort bin boundaries
min.l6.b1 <- log(1); min.l6.b2 <- log(129.5); min.l6.b3 <- log(614.5); min.l6.b4 <- log(1284.5); max.l6.b4 <- log(262144)
min.l5.b1 <- log(1); min.l5.b2 <- (log(178.5)-0.46994)/0.95213; min.l5.b3 <- (log(496.5)-0.46994)/0.95213; min.l5.b4 <- (log(2183.5)-0.46994)/0.95213; max.l5.b4 <- log(262144)
min.l13.b1<- log(1); min.l13.b2<- (log(137.5)-0.01014)/1.00035; min.l13.b3<- (log(329.5)-0.01014)/1.00035; min.l13.b4<- (log(938.5)-0.01014)/1.00035; max.l13.b4 <- log(262144)
min.l8.b1 <- log(1); min.l8.b2 <- log(140.5); min.l8.b3 <- log(472.5); min.l8.b4 <- log(1875.5); max.l8.b4 <- log(262144)
min.l10.b1 <- log(1); min.l10.b2 <- (log(73.5)+0.59618)/1.01805; min.l10.b3 <- (log(170.5)+0.59618)/1.01805; min.l10.b4 <- (log(374.5)+0.59618)/1.01805; max.l10.b4 <- log(262144)
min.l9.b1 <- log(1); min.l9.b2 <- (log(149.5)+0.82284)/1.17925; min.l9.b3 <- (log(483.5)+0.82284)/1.17925; min.l9.b4 <- (log(1147.5)+0.82284)/1.17925; max.l9.b4 <- log(262144)
min.l12.b1 <- log(1); min.l12.b2 <- (log(175.5)-0.5260)/0.9633; min.l12.b3 <- (log(374.5)-0.5260)/0.9633; min.l12.b4 <- (log(933.5)-0.5260)/0.9633; max.l12.b4 <- log(262144)
min.l11.b1 <- log(1); min.l11.b2 <- (log(121)+0.57030)/1.11911; min.l11.b3 <- (log(455)+0.57030)/1.11911; min.l11.b4 <- (log(1136.5)+0.57030)/1.11911; max.l11.b4 <- log(262144)

#interval censor wildtype observations
ERE.rep1.egka.p80.11P.cens.left <- sapply(1:length(ERE.rep1.egka.p80.11P), function(x) if(ERE.rep1.egka.p80.11P[x]<min.l6.b2){return(min.l6.b1)}else if(ERE.rep1.egka.p80.11P[x]<min.l6.b3){return(min.l6.b2)}else if(ERE.rep1.egka.p80.11P[x]<min.l6.b4){return(min.l6.b3)}else{return(min.l6.b4)})
ERE.rep1.egka.p80.11P.cens.right <- sapply(1:length(ERE.rep1.egka.p80.11P), function(x) if(ERE.rep1.egka.p80.11P[x]<min.l6.b2){return(min.l6.b2)}else if(ERE.rep1.egka.p80.11P[x]<min.l6.b3){return(min.l6.b3)}else if(ERE.rep1.egka.p80.11P[x]<min.l6.b4){return(min.l6.b4)}else{return(max.l6.b4)})
ERE.rep1.egka.p80.11P.cens <- data.frame(left=ERE.rep1.egka.p80.11P.cens.left,right=ERE.rep1.egka.p80.11P.cens.right)

ERE.rep2.egka.p80.11P.cens.left <- sapply(1:length(ERE.rep2.egka.p80.11P), function(x) if(ERE.rep2.egka.p80.11P[x]<min.l13.b2){return(min.l13.b1)}else if(ERE.rep2.egka.p80.11P[x]<min.l13.b3){return(min.l13.b2)}else if(ERE.rep2.egka.p80.11P[x]<min.l13.b4){return(min.l13.b3)}else{return(min.l13.b4)})
ERE.rep2.egka.p80.11P.cens.right <- sapply(1:length(ERE.rep2.egka.p80.11P), function(x) if(ERE.rep2.egka.p80.11P[x]<min.l13.b2){return(min.l13.b2)}else if(ERE.rep2.egka.p80.11P[x]<min.l13.b3){return(min.l13.b3)}else if(ERE.rep2.egka.p80.11P[x]<min.l13.b4){return(min.l13.b4)}else{return(max.l13.b4)})
ERE.rep2.egka.p80.11P.cens <- data.frame(left=ERE.rep2.egka.p80.11P.cens.left,right=ERE.rep2.egka.p80.11P.cens.right)

#for pooled meanF, want to interval censor some observations by rep1 bins, some by rep2: proportional to total number of reads from each rep
ERE.pooled.egka.p80.11P.cens.left <- c(sapply(1:round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l6.b2){return(min.l6.b1)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b3){return(min.l6.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b4){return(min.l6.b3)}else{return(min.l6.b4)}),sapply(round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)+1):length(ERE.pooled.egka.p80.11P), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l13.b2){return(min.l13.b1)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b3){return(min.l13.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b4){return(min.l13.b3)}else{return(min.l13.b4)}))
ERE.pooled.egka.p80.11P.cens.right <- c(sapply(1:round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l6.b2){return(min.l6.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b3){return(min.l6.b3)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b4){return(min.l6.b4)}else{return(max.l6.b4)}),sapply(round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)+1):length(ERE.pooled.egka.p80.11P), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l13.b2){return(min.l13.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b3){return(min.l13.b3)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b4){return(min.l13.b4)}else{return(max.l13.b4)}))
ERE.pooled.egka.p80.11P.cens <- data.frame(left=ERE.pooled.egka.p80.11P.cens.left,right=ERE.pooled.egka.p80.11P.cens.right)

SRE.rep1.GSKV.p80.11P.cens.left <- sapply(1:length(SRE.rep1.GSKV.p80.11P), function(x) if(SRE.rep1.GSKV.p80.11P[x]<min.l5.b2){return(min.l5.b1)}else if(SRE.rep1.GSKV.p80.11P[x]<min.l5.b3){return(min.l5.b2)}else if(SRE.rep1.GSKV.p80.11P[x]<min.l5.b4){return(min.l5.b3)}else{return(min.l5.b4)})
SRE.rep1.GSKV.p80.11P.cens.right <- sapply(1:length(SRE.rep1.GSKV.p80.11P), function(x) if(SRE.rep1.GSKV.p80.11P[x]<min.l5.b2){return(min.l5.b2)}else if(SRE.rep1.GSKV.p80.11P[x]<min.l5.b3){return(min.l5.b3)}else if(SRE.rep1.GSKV.p80.11P[x]<min.l5.b4){return(min.l5.b4)}else{return(max.l5.b4)})
SRE.rep1.GSKV.p80.11P.cens <- data.frame(left=SRE.rep1.GSKV.p80.11P.cens.left,right=SRE.rep1.GSKV.p80.11P.cens.right)

SRE.rep2.GSKV.p80.11P.cens.left <- sapply(1:length(SRE.rep2.GSKV.p80.11P), function(x) if(SRE.rep2.GSKV.p80.11P[x]<min.l8.b2){return(min.l8.b1)}else if(SRE.rep2.GSKV.p80.11P[x]<min.l8.b3){return(min.l8.b2)}else if(SRE.rep2.GSKV.p80.11P[x]<min.l8.b4){return(min.l8.b3)}else{return(min.l8.b4)})
SRE.rep2.GSKV.p80.11P.cens.right <- sapply(1:length(SRE.rep2.GSKV.p80.11P), function(x) if(SRE.rep2.GSKV.p80.11P[x]<min.l8.b2){return(min.l8.b2)}else if(SRE.rep2.GSKV.p80.11P[x]<min.l8.b3){return(min.l8.b3)}else if(SRE.rep2.GSKV.p80.11P[x]<min.l8.b4){return(min.l8.b4)}else{return(max.l8.b4)})
SRE.rep2.GSKV.p80.11P.cens <- data.frame(left=SRE.rep2.GSKV.p80.11P.cens.left,right=SRE.rep2.GSKV.p80.11P.cens.right)

SRE.pooled.GSKV.p80.11P.cens.left <- c(sapply(1:round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.11P)), function(x) if(SRE.pooled.GSKV.p80.11P[x]<min.l5.b2){return(min.l5.b1)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l5.b3){return(min.l5.b2)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l5.b4){return(min.l5.b3)}else{return(min.l5.b4)}),sapply(round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.11P)+1):length(SRE.pooled.GSKV.p80.11P), function(x) if(SRE.pooled.GSKV.p80.11P[x]<min.l8.b2){return(min.l8.b1)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l8.b3){return(min.l8.b2)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l8.b4){return(min.l8.b3)}else{return(min.l8.b4)}))
SRE.pooled.GSKV.p80.11P.cens.right <- c(sapply(1:round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.11P)), function(x) if(SRE.pooled.GSKV.p80.11P[x]<min.l5.b2){return(min.l5.b2)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l5.b3){return(min.l5.b3)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l5.b4){return(min.l5.b4)}else{return(max.l5.b4)}),sapply(round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.11P)+1):length(SRE.pooled.GSKV.p80.11P), function(x) if(SRE.pooled.GSKV.p80.11P[x]<min.l8.b2){return(min.l8.b2)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l8.b3){return(min.l8.b3)}else if(SRE.pooled.GSKV.p80.11P[x]<min.l8.b4){return(min.l8.b4)}else{return(max.l8.b4)}))
SRE.pooled.GSKV.p80.11P.cens <- data.frame(left=SRE.pooled.GSKV.p80.11P.cens.left,right=SRE.pooled.GSKV.p80.11P.cens.right)

ERE.rep1.egka.p80.SR1.cens.left <- sapply(1:length(ERE.rep1.egka.p80.SR1), function(x) if(ERE.rep1.egka.p80.SR1[x]<min.l10.b2){return(min.l10.b1)}else if(ERE.rep1.egka.p80.SR1[x]<min.l10.b3){return(min.l10.b2)}else if(ERE.rep1.egka.p80.SR1[x]<min.l10.b4){return(min.l10.b3)}else{return(min.l10.b4)})
ERE.rep1.egka.p80.SR1.cens.right <- sapply(1:length(ERE.rep1.egka.p80.SR1), function(x) if(ERE.rep1.egka.p80.SR1[x]<min.l10.b2){return(min.l10.b2)}else if(ERE.rep1.egka.p80.SR1[x]<min.l10.b3){return(min.l10.b3)}else if(ERE.rep1.egka.p80.SR1[x]<min.l10.b4){return(min.l10.b4)}else{return(max.l10.b4)})
ERE.rep1.egka.p80.SR1.cens <- data.frame(left=ERE.rep1.egka.p80.SR1.cens.left,right=ERE.rep1.egka.p80.SR1.cens.right)

ERE.rep2.egka.p80.SR1.cens.left <- sapply(1:length(ERE.rep2.egka.p80.SR1), function(x) if(ERE.rep2.egka.p80.SR1[x]<min.l12.b2){return(min.l12.b1)}else if(ERE.rep2.egka.p80.SR1[x]<min.l12.b3){return(min.l12.b2)}else if(ERE.rep2.egka.p80.SR1[x]<min.l12.b4){return(min.l12.b3)}else{return(min.l12.b4)})
ERE.rep2.egka.p80.SR1.cens.right <- sapply(1:length(ERE.rep2.egka.p80.SR1), function(x) if(ERE.rep2.egka.p80.SR1[x]<min.l12.b2){return(min.l12.b2)}else if(ERE.rep2.egka.p80.SR1[x]<min.l12.b3){return(min.l12.b3)}else if(ERE.rep2.egka.p80.SR1[x]<min.l12.b4){return(min.l12.b4)}else{return(max.l12.b4)})
ERE.rep2.egka.p80.SR1.cens <- data.frame(left=ERE.rep2.egka.p80.SR1.cens.left,right=ERE.rep2.egka.p80.SR1.cens.right)

ERE.pooled.egka.p80.SR1.cens.left <- c(sapply(1:round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l10.b2){return(min.l10.b1)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b3){return(min.l10.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b4){return(min.l10.b3)}else{return(min.l10.b4)}),sapply(round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)+1):length(ERE.pooled.egka.p80.SR1), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l12.b2){return(min.l12.b1)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b3){return(min.l12.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b4){return(min.l12.b3)}else{return(min.l12.b4)}))
ERE.pooled.egka.p80.SR1.cens.right <- c(sapply(1:round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l10.b2){return(min.l10.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b3){return(min.l10.b3)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b4){return(min.l10.b4)}else{return(max.l10.b4)}),sapply(round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)+1):length(ERE.pooled.egka.p80.SR1), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l12.b2){return(min.l12.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b3){return(min.l12.b3)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b4){return(min.l12.b4)}else{return(max.l12.b4)}))
ERE.pooled.egka.p80.SR1.cens <- data.frame(left=ERE.pooled.egka.p80.SR1.cens.left,right=ERE.pooled.egka.p80.SR1.cens.right)

SRE.rep1.GSKV.p80.SR1.cens.left <- sapply(1:length(SRE.rep1.GSKV.p80.SR1), function(x) if(SRE.rep1.GSKV.p80.SR1[x]<min.l9.b2){return(min.l9.b1)}else if(SRE.rep1.GSKV.p80.SR1[x]<min.l9.b3){return(min.l9.b2)}else if(SRE.rep1.GSKV.p80.SR1[x]<min.l9.b4){return(min.l9.b3)}else{return(min.l9.b4)})
SRE.rep1.GSKV.p80.SR1.cens.right <- sapply(1:length(SRE.rep1.GSKV.p80.SR1), function(x) if(SRE.rep1.GSKV.p80.SR1[x]<min.l9.b2){return(min.l9.b2)}else if(SRE.rep1.GSKV.p80.SR1[x]<min.l9.b3){return(min.l9.b3)}else if(SRE.rep1.GSKV.p80.SR1[x]<min.l9.b4){return(min.l9.b4)}else{return(max.l9.b4)})
SRE.rep1.GSKV.p80.SR1.cens <- data.frame(left=SRE.rep1.GSKV.p80.SR1.cens.left,right=SRE.rep1.GSKV.p80.SR1.cens.right)

SRE.rep2.GSKV.p80.SR1.cens.left <- sapply(1:length(SRE.rep2.GSKV.p80.SR1), function(x) if(SRE.rep2.GSKV.p80.SR1[x]<min.l11.b2){return(min.l11.b1)}else if(SRE.rep2.GSKV.p80.SR1[x]<min.l11.b3){return(min.l11.b2)}else if(SRE.rep2.GSKV.p80.SR1[x]<min.l11.b4){return(min.l11.b3)}else{return(min.l11.b4)})
SRE.rep2.GSKV.p80.SR1.cens.right <- sapply(1:length(SRE.rep2.GSKV.p80.SR1), function(x) if(SRE.rep2.GSKV.p80.SR1[x]<min.l11.b2){return(min.l11.b2)}else if(SRE.rep2.GSKV.p80.SR1[x]<min.l11.b3){return(min.l11.b3)}else if(SRE.rep2.GSKV.p80.SR1[x]<min.l11.b4){return(min.l11.b4)}else{return(max.l11.b4)})
SRE.rep2.GSKV.p80.SR1.cens <- data.frame(left=SRE.rep2.GSKV.p80.SR1.cens.left,right=SRE.rep2.GSKV.p80.SR1.cens.right)

SRE.pooled.GSKV.p80.SR1.cens.left <- c(sapply(1:round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.SR1)), function(x) if(SRE.pooled.GSKV.p80.SR1[x]<min.l9.b2){return(min.l9.b1)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l9.b3){return(min.l9.b2)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l9.b4){return(min.l9.b3)}else{return(min.l9.b4)}),sapply(round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.SR1)+1):length(SRE.pooled.GSKV.p80.SR1), function(x) if(SRE.pooled.GSKV.p80.SR1[x]<min.l11.b2){return(min.l11.b1)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l11.b3){return(min.l11.b2)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l11.b4){return(min.l11.b3)}else{return(min.l11.b4)}))
SRE.pooled.GSKV.p80.SR1.cens.right <- c(sapply(1:round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.SR1)), function(x) if(SRE.pooled.GSKV.p80.SR1[x]<min.l9.b2){return(min.l9.b2)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l9.b3){return(min.l9.b3)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l9.b4){return(min.l9.b4)}else{return(max.l9.b4)}),sapply(round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p80.SR1)+1):length(SRE.pooled.GSKV.p80.SR1), function(x) if(SRE.pooled.GSKV.p80.SR1[x]<min.l11.b2){return(min.l11.b2)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l11.b3){return(min.l11.b3)}else if(SRE.pooled.GSKV.p80.SR1[x]<min.l11.b4){return(min.l11.b4)}else{return(max.l11.b4)}))
SRE.pooled.GSKV.p80.SR1.cens <- data.frame(left=SRE.pooled.GSKV.p80.SR1.cens.left,right=SRE.pooled.GSKV.p80.SR1.cens.right)

boot.cens.dist <- function(data, cfu, rep){
  sample.means <- vector()
  while(length(sample.means)<rep){
    boot <- data[sample(1:nrow(data),cfu,replace=T),]
    sample.means <- c(sample.means,tryCatch( summary(fitdistcens(boot,"logis"))$estimate["location"],error=function(e){return(as.numeric(NA))}))
  }
  return(sample.means)
}

ERE.rep1.pos.11P <- data.frame(bin=1:25); ERE.rep1.pos.11P$range.cfu <- ERE.rep1.null.11P$range.cfu; ERE.rep1.pos.11P$median.cfu <- ERE.rep1.null.11P$median.cfu
for(i in 1:nrow(ERE.rep1.pos.11P)){
  ERE.rep1.pos.11P$list.means[i] <- list(boot.cens.dist(ERE.rep1.egka.p80.11P.cens, round(ERE.rep1.pos.11P$median.cfu[i]),10000))
  ERE.rep1.pos.11P$mean.meanF[i] <- mean(ERE.rep1.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.rep1.pos.11P$sd.meanF[i] <- sd(ERE.rep1.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.rep1.pos.11P$max.meanF[i] <- max(ERE.rep1.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(ERE.rep1.pos.11P,file="3_classify-variants_define-distributions_out/ERE.rep1.pos.11P.Rda")

ERE.rep2.pos.11P <- data.frame(bin=1:25); ERE.rep2.pos.11P$range.cfu <- ERE.rep2.null.11P$range.cfu; ERE.rep2.pos.11P$median.cfu <- ERE.rep2.null.11P$median.cfu
for(i in 1:nrow(ERE.rep2.pos.11P)){
  ERE.rep2.pos.11P$list.means[i] <- list(boot.cens.dist(ERE.rep2.egka.p80.11P.cens, round(ERE.rep2.pos.11P$median.cfu[i]),10000))
  ERE.rep2.pos.11P$mean.meanF[i] <- mean(ERE.rep2.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.rep2.pos.11P$sd.meanF[i] <- sd(ERE.rep2.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.rep2.pos.11P$max.meanF[i] <- max(ERE.rep2.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(ERE.rep2.pos.11P,file="3_classify-variants_define-distributions_out/ERE.rep2.pos.11P.Rda")

ERE.pooled.pos.11P <- data.frame(bin=1:25); ERE.pooled.pos.11P$range.cfu <- ERE.pooled.null.11P$range.cfu; ERE.pooled.pos.11P$median.cfu <- ERE.pooled.null.11P$median.cfu
for(i in 1:nrow(ERE.pooled.pos.11P)){
  ERE.pooled.pos.11P$list.means[i] <- list(boot.cens.dist(ERE.pooled.egka.p80.11P.cens, round(ERE.pooled.pos.11P$median.cfu[i]),10000))
  ERE.pooled.pos.11P$mean.meanF[i] <- mean(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.11P$sd.meanF[i] <- sd(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.11P$max.meanF[i] <- max(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(ERE.pooled.pos.11P,file="3_classify-variants_define-distributions_out/ERE.pooled.pos.11P.Rda")

SRE.rep1.pos.11P <- data.frame(bin=1:25); SRE.rep1.pos.11P$range.cfu <- SRE.rep1.null.11P$range.cfu; SRE.rep1.pos.11P$median.cfu <- SRE.rep1.null.11P$median.cfu
for(i in 1:nrow(SRE.rep1.pos.11P)){
  SRE.rep1.pos.11P$list.means[i] <- list(boot.cens.dist(SRE.rep1.GSKV.p80.11P.cens, round(SRE.rep1.pos.11P$median.cfu[i]),10000))
  SRE.rep1.pos.11P$mean.meanF[i] <- mean(SRE.rep1.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.rep1.pos.11P$sd.meanF[i] <- sd(SRE.rep1.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.rep1.pos.11P$max.meanF[i] <- max(SRE.rep1.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(SRE.rep1.pos.11P,file="3_classify-variants_define-distributions_out/SRE.rep1.pos.11P.Rda")

SRE.rep2.pos.11P <- data.frame(bin=1:25); SRE.rep2.pos.11P$range.cfu <- SRE.rep2.null.11P$range.cfu; SRE.rep2.pos.11P$median.cfu <- SRE.rep2.null.11P$median.cfu
for(i in 1:nrow(SRE.rep2.pos.11P)){
  SRE.rep2.pos.11P$list.means[i] <- list(boot.cens.dist(SRE.rep2.GSKV.p80.11P.cens, round(SRE.rep2.pos.11P$median.cfu[i]),10000))
  SRE.rep2.pos.11P$mean.meanF[i] <- mean(SRE.rep2.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.rep2.pos.11P$sd.meanF[i] <- sd(SRE.rep2.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.rep2.pos.11P$max.meanF[i] <- max(SRE.rep2.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(SRE.rep2.pos.11P,file="3_classify-variants_define-distributions_out/SRE.rep2.pos.11P.Rda")

SRE.pooled.pos.11P <- data.frame(bin=1:25); SRE.pooled.pos.11P$range.cfu <- SRE.pooled.null.11P$range.cfu; SRE.pooled.pos.11P$median.cfu <- SRE.pooled.null.11P$median.cfu
for(i in 1:nrow(SRE.pooled.pos.11P)){
  SRE.pooled.pos.11P$list.means[i] <- list(boot.cens.dist(SRE.pooled.GSKV.p80.11P.cens, round(SRE.pooled.pos.11P$median.cfu[i]),10000))
  SRE.pooled.pos.11P$mean.meanF[i] <- mean(SRE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.11P$sd.meanF[i] <- sd(SRE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.11P$max.meanF[i] <- max(SRE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(SRE.pooled.pos.11P,file="3_classify-variants_define-distributions_out/SRE.pooled.pos.11P.Rda")


ERE.rep1.pos.SR1 <- data.frame(bin=1:25); ERE.rep1.pos.SR1$range.cfu <- ERE.rep1.null.SR1$range.cfu; ERE.rep1.pos.SR1$median.cfu <- ERE.rep1.null.SR1$median.cfu
for(i in 1:nrow(ERE.rep1.pos.SR1)){
  ERE.rep1.pos.SR1$list.means[i] <- list(boot.cens.dist(ERE.rep1.egka.p80.SR1.cens, round(ERE.rep1.pos.SR1$median.cfu[i]),10000))
  ERE.rep1.pos.SR1$mean.meanF[i] <- mean(ERE.rep1.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.rep1.pos.SR1$sd.meanF[i] <- sd(ERE.rep1.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.rep1.pos.SR1$max.meanF[i] <- max(ERE.rep1.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(ERE.rep1.pos.SR1,file="3_classify-variants_define-distributions_out/ERE.rep1.pos.SR1.Rda")

ERE.rep2.pos.SR1 <- data.frame(bin=1:25); ERE.rep2.pos.SR1$range.cfu <- ERE.rep2.null.SR1$range.cfu; ERE.rep2.pos.SR1$median.cfu <- ERE.rep2.null.SR1$median.cfu
for(i in 1:nrow(ERE.rep2.pos.SR1)){
  ERE.rep2.pos.SR1$list.means[i] <- list(boot.cens.dist(ERE.rep2.egka.p80.SR1.cens, round(ERE.rep2.pos.SR1$median.cfu[i]),10000))
  ERE.rep2.pos.SR1$mean.meanF[i] <- mean(ERE.rep2.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.rep2.pos.SR1$sd.meanF[i] <- sd(ERE.rep2.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.rep2.pos.SR1$max.meanF[i] <- max(ERE.rep2.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(ERE.rep2.pos.SR1,file="3_classify-variants_define-distributions_out/ERE.rep2.pos.SR1.Rda")

ERE.pooled.pos.SR1 <- data.frame(bin=1:25); ERE.pooled.pos.SR1$range.cfu <- ERE.pooled.null.SR1$range.cfu; ERE.pooled.pos.SR1$median.cfu <- ERE.pooled.null.SR1$median.cfu
for(i in 1:nrow(ERE.pooled.pos.SR1)){
  ERE.pooled.pos.SR1$list.means[i] <- list(boot.cens.dist(ERE.pooled.egka.p80.SR1.cens, round(ERE.pooled.pos.SR1$median.cfu[i]),10000))
  ERE.pooled.pos.SR1$mean.meanF[i] <- mean(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.SR1$sd.meanF[i] <- sd(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.SR1$max.meanF[i] <- max(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(ERE.pooled.pos.SR1,file="3_classify-variants_define-distributions_out/ERE.pooled.pos.SR1.Rda")

SRE.rep1.pos.SR1 <- data.frame(bin=1:25); SRE.rep1.pos.SR1$range.cfu <- SRE.rep1.null.SR1$range.cfu; SRE.rep1.pos.SR1$median.cfu <- SRE.rep1.null.SR1$median.cfu
for(i in 1:nrow(SRE.rep1.pos.SR1)){
  SRE.rep1.pos.SR1$list.means[i] <- list(boot.cens.dist(SRE.rep1.GSKV.p80.SR1.cens, round(SRE.rep1.pos.SR1$median.cfu[i]),10000))
  SRE.rep1.pos.SR1$mean.meanF[i] <- mean(SRE.rep1.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.rep1.pos.SR1$sd.meanF[i] <- sd(SRE.rep1.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.rep1.pos.SR1$max.meanF[i] <- max(SRE.rep1.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(SRE.rep1.pos.SR1,file="3_classify-variants_define-distributions_out/SRE.rep1.pos.SR1.Rda")

SRE.rep2.pos.SR1 <- data.frame(bin=1:25); SRE.rep2.pos.SR1$range.cfu <- SRE.rep2.null.SR1$range.cfu; SRE.rep2.pos.SR1$median.cfu <- SRE.rep2.null.SR1$median.cfu
for(i in 1:nrow(SRE.rep2.pos.SR1)){
  SRE.rep2.pos.SR1$list.means[i] <- list(boot.cens.dist(SRE.rep2.GSKV.p80.SR1.cens, round(SRE.rep2.pos.SR1$median.cfu[i]),10000))
  SRE.rep2.pos.SR1$mean.meanF[i] <- mean(SRE.rep2.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.rep2.pos.SR1$sd.meanF[i] <- sd(SRE.rep2.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.rep2.pos.SR1$max.meanF[i] <- max(SRE.rep2.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(SRE.rep2.pos.SR1,file="3_classify-variants_define-distributions_out/SRE.rep2.pos.SR1.Rda")

SRE.pooled.pos.SR1 <- data.frame(bin=1:25); SRE.pooled.pos.SR1$range.cfu <- SRE.pooled.null.SR1$range.cfu; SRE.pooled.pos.SR1$median.cfu <- SRE.pooled.null.SR1$median.cfu
for(i in 1:nrow(SRE.pooled.pos.SR1)){
  SRE.pooled.pos.SR1$list.means[i] <- list(boot.cens.dist(SRE.pooled.GSKV.p80.SR1.cens, round(SRE.pooled.pos.SR1$median.cfu[i]),10000))
  SRE.pooled.pos.SR1$mean.meanF[i] <- mean(SRE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.SR1$sd.meanF[i] <- sd(SRE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.SR1$max.meanF[i] <- max(SRE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(SRE.pooled.pos.SR1,file="3_classify-variants_define-distributions_out/SRE.pooled.pos.SR1.Rda")
