#16 Jan 2017
#TNS
#script to assess quality of replicates:
# -- assess bulk meanF values compared to variants picked and isolated in monoculture
# -- assess reproducibility in meanF, classification between replicates
# -- figure out minimum # cfu cutoff
# -- censor values below cutoff, output final table
# -- make plots of correlations in meanF, distribution of specificity classes in censored data

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)

#read in tables with calculated meanF for each variant, output of classify-variants.R
dt.11P <- read.table(file="4_classify-variants_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.SR1 <- read.table(file="4_classify-variants_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

#read in meanF of isogenic clones picked from post-sort plates, re-induced, analyzed via flow cytometry and genotyped
controls.ERE.11P <- read.table(file="5_assess-replicates_in/ERE-controls_11P.csv",header=T,sep=",")
controls.SRE.11P <- read.table(file="5_assess-replicates_in/SRE-controls_11P.csv",header=T,sep=",")
controls.ERE.SR1 <- read.table(file="5_assess-replicates_in/ERE-controls_SR1.csv",header=T,sep=",")
controls.SRE.SR1 <- read.table(file="5_assess-replicates_in/SRE-controls_SR1.csv",header=T,sep=",")

#for each set of isogenic clones, fill in bulk meanF estimates
for(i in 1:nrow(controls.ERE.11P)){
  controls.ERE.11P$rep1.meanF[i] <- dt.11P[as.character(controls.ERE.11P[i,"AAseq"]),ERE.rep1.meanF]
  controls.ERE.11P$rep2.meanF[i] <- dt.11P[as.character(controls.ERE.11P[i,"AAseq"]),ERE.rep2.meanF]
  controls.ERE.11P$pooled.meanF[i] <- dt.11P[as.character(controls.ERE.11P[i,"AAseq"]),ERE.pooled.meanF]
  controls.ERE.11P$rep1.cfu[i] <- dt.11P[as.character(controls.ERE.11P[i,"AAseq"]),ERE.rep1.cfu]
  controls.ERE.11P$rep2.cfu[i] <- dt.11P[as.character(controls.ERE.11P[i,"AAseq"]),ERE.rep2.cfu]
}
for(i in 1:nrow(controls.SRE.11P)){
  controls.SRE.11P$rep1.meanF[i] <- dt.11P[as.character(controls.SRE.11P[i,"AAseq"]),SRE.rep1.meanF]
  controls.SRE.11P$rep2.meanF[i] <- dt.11P[as.character(controls.SRE.11P[i,"AAseq"]),SRE.rep2.meanF]
  controls.SRE.11P$pooled.meanF[i] <- dt.11P[as.character(controls.SRE.11P[i,"AAseq"]),SRE.pooled.meanF]
  controls.SRE.11P$rep1.cfu[i] <- dt.11P[as.character(controls.SRE.11P[i,"AAseq"]),SRE.rep1.cfu]
  controls.SRE.11P$rep2.cfu[i] <- dt.11P[as.character(controls.SRE.11P[i,"AAseq"]),SRE.rep2.cfu]
}
for(i in 1:nrow(controls.ERE.SR1)){
  controls.ERE.SR1$rep1.meanF[i] <- dt.SR1[as.character(controls.ERE.SR1[i,"AAseq"]),ERE.rep1.meanF]
  controls.ERE.SR1$rep2.meanF[i] <- dt.SR1[as.character(controls.ERE.SR1[i,"AAseq"]),ERE.rep2.meanF]
  controls.ERE.SR1$pooled.meanF[i] <- dt.SR1[as.character(controls.ERE.SR1[i,"AAseq"]),ERE.pooled.meanF]
  controls.ERE.SR1$rep1.cfu[i] <- dt.SR1[as.character(controls.ERE.SR1[i,"AAseq"]),ERE.rep1.cfu]
  controls.ERE.SR1$rep2.cfu[i] <- dt.SR1[as.character(controls.ERE.SR1[i,"AAseq"]),ERE.rep2.cfu]
}
for(i in 1:nrow(controls.SRE.SR1)){
  controls.SRE.SR1$rep1.meanF[i] <- dt.SR1[as.character(controls.SRE.SR1[i,"AAseq"]),SRE.rep1.meanF]
  controls.SRE.SR1$rep2.meanF[i] <- dt.SR1[as.character(controls.SRE.SR1[i,"AAseq"]),SRE.rep2.meanF]
  controls.SRE.SR1$pooled.meanF[i] <- dt.SR1[as.character(controls.SRE.SR1[i,"AAseq"]),SRE.pooled.meanF]
  controls.SRE.SR1$rep1.cfu[i] <- dt.SR1[as.character(controls.SRE.SR1[i,"AAseq"]),SRE.rep1.cfu]
  controls.SRE.SR1$rep2.cfu[i] <- dt.SR1[as.character(controls.SRE.SR1[i,"AAseq"]),SRE.rep2.cfu]
}

#plot bulk vs independent meanF, give R2 of correlation
pdf(file="5_assess-replicates_out/bulk-quant_v_isolated-clone.pdf",width=6,height=7,useDingbats=F)
par(mar=c(6.1,5.1,4.1,2.1),mfrow=c(2,2))
plot(controls.ERE.SR1$meanF,controls.ERE.SR1$pooled.meanF,xlab="isolated clone,\nmean fluorescence (a.u.)",ylab="bulk quantification,\nmean fluorescence (a.u.)",pch=19,col=rgb(100,69,155,maxColorValue=255),main="ERE, SR1"); abline(lm(controls.ERE.SR1$pooled.meanF~controls.ERE.SR1$meanF),lty=2,lwd=2); summary(lm(controls.ERE.SR1$pooled.meanF~controls.ERE.SR1$meanF))
legend("bottomright",legend="R-squared: 0.88",bty="n",cex=0.8)
plot(controls.SRE.SR1$meanF,controls.SRE.SR1$pooled.meanF,xlab="isolated clone,\nmean fluorescence (a.u.)",ylab="bulk quantification,\nmean fluorescence (a.u.)",pch=19,col=rgb(5,172,72,maxColorValue=255),main="SRE, SR1"); abline(lm(controls.SRE.SR1$pooled.meanF~controls.SRE.SR1$meanF),lwd=2,lty=2); summary(lm(controls.SRE.SR1$pooled.meanF~controls.SRE.SR1$meanF))
legend("bottomright",legend="R-squared: 0.79",bty="n",cex=0.8)
plot(controls.ERE.11P$meanF,controls.ERE.11P$pooled.meanF,xlab="isolated clone,\nmean fluorescence (a.u.)",ylab="bulk quantification,\nmean fluorescence (a.u.)",pch=19,col=rgb(100,69,155,maxColorValue=255),main="ERE, 11P"); abline(lm(controls.ERE.11P$pooled.meanF~controls.ERE.11P$meanF),lty=2,lwd=2);summary(lm(controls.ERE.11P$pooled.meanF~controls.ERE.11P$meanF))
legend("bottomright",legend="R-squared: 0.87",bty="n",cex=0.8)
plot(controls.SRE.11P$meanF,controls.SRE.11P$pooled.meanF,xlab="isolated clone,\nmean fluorescence (a.u.)",ylab="bulk quantification,\nmean fluorescence (a.u.)",pch=19,col=rgb(5,172,72,maxColorValue=255),main="SRE, 11P"); abline(lm(controls.SRE.11P$pooled.meanF~controls.SRE.11P$meanF),lty=2,lwd=2);summary(lm(controls.SRE.11P$pooled.meanF~controls.SRE.11P$meanF))
legend("bottomright",legend="R-squared: 0.92",bty="n",cex=0.8)
dev.off()

################################################################################################
#11P: 
#define cfu.min as the lesser amt of cells a genotype is represented in between the two replicates
dt.11P[,ERE.cfu.min := min(ERE.rep1.cfu,ERE.rep2.cfu),by=AAseq]
dt.11P[,SRE.cfu.min := min(SRE.rep1.cfu,SRE.rep2.cfu),by=AAseq]

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

#ERE
#make ~even bins on the basis of the minimum number of reads between replicates
bins.ERE.11P <- as.numeric(quantile(dt.11P.coding[!is.na(ERE.rep1.class) & !is.na(ERE.rep2.class),ERE.cfu.min],seq(0,1,by=0.05)))
shared.strong.ERE.11P <- shared.positive.ERE.11P <- vector(mode="numeric",length=length(bins.ERE.11P)-1)
for(i in 2:(length(bins.ERE.11P)-1)){
  data <- dt.11P.coding[ERE.cfu.min >= bins.ERE.11P[i] & ERE.cfu.min < bins.ERE.11P[i+1] & !is.na(dt.11P.coding$ERE.rep1.class) & !is.na(dt.11P.coding$ERE.rep2.class),]
  shared.strong.ERE.11P[i] <- (nrow(data[data$ERE.rep1.class %in% c("strong") & data$ERE.rep2.class %in% c("strong"),])/nrow(data[data$ERE.rep1.class %in% c("strong"),]) + nrow(data[data$ERE.rep1.class %in% c("strong") & data$ERE.rep2.class %in% c("strong"),])/nrow(data[data$ERE.rep2.class %in% c("strong"),]))/2
  shared.positive.ERE.11P[i] <- (nrow(data[data$ERE.rep1.class %in% c("weak","strong") & data$ERE.rep2.class %in% c("weak","strong"),])/nrow(data[data$ERE.rep1.class %in% c("weak","strong"),]) + nrow(data[data$ERE.rep1.class %in% c("weak","strong") & data$ERE.rep2.class %in% c("weak","strong"),])/nrow(data[data$ERE.rep2.class %in% c("weak","strong"),]))/2
}

#repeat SRE
bins.SRE.11P <- as.numeric(quantile(dt.11P.coding[!is.na(SRE.rep1.class) & !is.na(SRE.rep2.class),SRE.cfu.min],seq(0,1,by=0.05)))
shared.strong.SRE.11P <- shared.positive.SRE.11P <- vector(mode="numeric",length=length(bins.SRE.11P)-1)
for(i in 2:(length(bins.SRE.11P)-1)){
  data <- dt.11P.coding[SRE.cfu.min >= bins.SRE.11P[i] & SRE.cfu.min < bins.SRE.11P[i+1] & !is.na(dt.11P.coding$SRE.rep1.class) & !is.na(dt.11P.coding$SRE.rep2.class),]
  shared.strong.SRE.11P[i] <- (nrow(data[data$SRE.rep1.class %in% c("strong") & data$SRE.rep2.class %in% c("strong"),])/nrow(data[data$SRE.rep1.class %in% c("strong"),]) + nrow(data[data$SRE.rep1.class %in% c("strong") & data$SRE.rep2.class %in% c("strong"),])/nrow(data[data$SRE.rep2.class %in% c("strong"),]))/2
  shared.positive.SRE.11P[i] <- (nrow(data[data$SRE.rep1.class %in% c("weak","strong") & data$SRE.rep2.class %in% c("weak","strong"),])/nrow(data[data$SRE.rep1.class %in% c("weak","strong"),]) + nrow(data[data$SRE.rep1.class %in% c("weak","strong") & data$SRE.rep2.class %in% c("weak","strong"),])/nrow(data[data$SRE.rep2.class %in% c("weak","strong"),]))/2
}

pdf(file="./5_assess-replicates_out/11P_shared-proportion-positive_v_min-cfu-binned-max-min.pdf",width=8,height=4.5)
par(mfrow=c(1,2))
plot(bins.ERE.11P[1:20],shared.strong.ERE.11P,type="l",lwd=2,log="x",col=rgb(100,69,155,maxColorValue=255),xlab="minimum number of reads, binned",ylab="average shared proportion strong btwn reps",xlim=c(1,1000),main="",ylim=c(0,1))
points(bins.SRE.11P[1:20],shared.strong.SRE.11P,type="l",lty=1,lwd=2,col=rgb(5,172,72,maxColorValue=255))
abline(v=15,lty=2)
legend("bottomright",legend=c("ERE","SRE"),lty=c(1,1),lwd=2,col=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),bty="n")

plot(bins.ERE.11P[1:20],shared.positive.ERE.11P,type="l",lwd=2,log="x",col=rgb(100,69,155,maxColorValue=255),xlab="minimum number of reads, binned",ylab="average shared proportion positive btwn reps",xlim=c(1,1000),main="",ylim=c(0,1))
points(bins.SRE.11P[1:20],shared.positive.SRE.11P,type="l",lty=1,lwd=2,col=rgb(5,172,72,maxColorValue=255))
abline(v=15,lty=2)
legend("bottomright",legend=c("ERE","SRE"),lty=c(1,1),lwd=2,col=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),bty="n")
dev.off()

#censor all classifications based on less than 15 cfu
dt.11P[ERE.rep1.cfu<15, ERE.rep1.class := as.factor(NA)];dt.11P[ERE.rep2.cfu<15, ERE.rep2.class := as.factor(NA)];dt.11P[ERE.pooled.cfu<15, ERE.pooled.class := as.factor(NA)]
dt.11P[SRE.rep1.cfu<15, SRE.rep1.class := as.factor(NA)];dt.11P[SRE.rep2.cfu<15, SRE.rep2.class := as.factor(NA)];dt.11P[SRE.pooled.cfu<15, SRE.pooled.class := as.factor(NA)]

#output table
write.table(dt.11P[,.(AAseq,ERE.rep1.b1,ERE.rep1.b2,ERE.rep1.b3,ERE.rep1.b4,ERE.rep2.b1,ERE.rep2.b2,ERE.rep2.b3,ERE.rep2.b4,SRE.rep1.b1,SRE.rep1.b2,SRE.rep1.b3,SRE.rep1.b4,SRE.rep2.b1,SRE.rep2.b2,SRE.rep2.b3,SRE.rep2.b4,ERE.rep1.cfu,ERE.rep2.cfu,ERE.pooled.cfu,SRE.rep1.cfu,SRE.rep2.cfu,SRE.pooled.cfu,ERE.rep1.meanF,ERE.rep2.meanF,ERE.pooled.meanF,SRE.rep1.meanF,SRE.rep2.meanF,SRE.pooled.meanF,ERE.rep1.class,ERE.rep2.class,ERE.pooled.class,SRE.rep1.class,SRE.rep2.class,SRE.pooled.class)],file="5_assess-replicates_out/dt_output_11P_censor15.csv",col.names=T,row.names=F,quote=F,sep=",")

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

#plot by specificity
pdf(file="5_assess-replicates_out/11P_pooled_ERE-meanF_v_SRE-meanF_points-colored-by-specificity.pdf",5.5,6,useDingbats=F)
plot(dt.11P.coding$SRE.pooled.meanF,dt.11P.coding$ERE.pooled.meanF,pch="",xlab="SRE meanF (a.u.), pooled reps",ylab="ERE meanF (a.u.), pooled reps",xlim=c(3,10),ylim=c(3,10))
points(dt.11P.coding[ERE.pooled.class %in% c("null","weak") & SRE.pooled.class %in% c("null","weak"),SRE.pooled.meanF],dt.11P.coding[ERE.pooled.class %in% c("null","weak") & SRE.pooled.class %in% c("null","weak"),ERE.pooled.meanF],pch=20,cex=0.6)
points(dt.11P.stop[!is.na(as.character(ERE.pooled.class)) & !is.na(as.character(SRE.pooled.class)),SRE.pooled.meanF],dt.11P.stop[!is.na(as.character(ERE.pooled.class)) & !is.na(as.character(SRE.pooled.class)),ERE.pooled.meanF],pch=20,cex=0.6,col="coral")
points(dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="null",SRE.pooled.meanF],dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="null",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(100,69,155,maxColorValue=255))
points(dt.11P.coding[ERE.pooled.class=="null" & SRE.pooled.class=="strong",SRE.pooled.meanF],dt.11P.coding[ERE.pooled.class=="null" & SRE.pooled.class=="strong",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(5,172,72,maxColorValue=255))
points(dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="weak",SRE.pooled.meanF],dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="weak",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(111,204,221,maxColorValue=255))
points(dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="strong",SRE.pooled.meanF],dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="strong",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(111,204,221,maxColorValue=255))
points(dt.11P.coding[ERE.pooled.class=="weak" & SRE.pooled.class=="strong",SRE.pooled.meanF],dt.11P.coding[ERE.pooled.class=="weak" & SRE.pooled.class=="strong",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(111,204,221,maxColorValue=255))
abline(v=dt.11P["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
abline(h=7.104872,lwd=2,lty=2,col=rgb(100,69,155,maxColorValue=255)) #mean(ERE.pooled.egka)
dev.off()

#correlation in meanF for variants with >15 cfu cutoff, ERE
pdf(file="5_assess-replicates_out/11P_ERE_meanF-correlation.pdf",6.25,6.75,useDingbats=F)
plot(dt.11P.coding[ERE.cfu.min>15,ERE.rep1.meanF],dt.11P.coding[ERE.cfu.min>15,ERE.rep2.meanF],pch=19,xlab="ERE mean fluorescence, replicate 1 (a.u.)",ylab="ERE mean fluorescence, replicate 2 (a.u.)",cex.lab=1.4)
points(dt.11P.stop[ERE.cfu.min>15,ERE.rep1.meanF],dt.11P.stop[ERE.cfu.min>15,ERE.rep2.meanF],pch=19,col="coral")
summary(lm(dt.11P.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep2.meanF]~dt.11P.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep1.meanF]))#;abline(lm(dt.11P.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep2.meanF]~dt.11P.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep1.meanF]),lty=2,lwd=3,col="red")
legend("topleft",legend="R-squared_pos: 0.58",bty="n",cex=0.8)
dev.off()

#correlation in meanF for variants with >15 cfu cutoff, SRE
pdf(file="5_assess-replicates_out/11P_SRE_meanF-correlation.pdf",6.25,6.75,useDingbats=F)
plot(dt.11P.coding[SRE.cfu.min>15,SRE.rep1.meanF],dt.11P.coding[SRE.cfu.min>15,SRE.rep2.meanF],pch=19,xlab="SRE mean fluorescence, replicate 1 (a.u.)",ylab="SRE mean fluorescence, replicate 2 (a.u.)",cex.lab=1.4)
points(dt.11P.stop[SRE.cfu.min>15,SRE.rep1.meanF],dt.11P.stop[SRE.cfu.min>15,SRE.rep2.meanF],pch=19,col="coral")
summary(lm(dt.11P.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep2.meanF]~dt.11P.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep1.meanF]))#;abline(lm(dt.11P.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep2.meanF]~dt.11P.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep1.meanF]),lty=2,lwd=3,col="red")
legend("topleft",legend="R-squared_pos: 0.76",bty="n",cex=0.8)
dev.off()

###################################################################################################################
#SR1: 
#define cfu.min as the lesser amt of cells a genotype is represented in between the two replicates
dt.SR1[,ERE.cfu.min := min(ERE.rep1.cfu,ERE.rep2.cfu),by=AAseq]
dt.SR1[,SRE.cfu.min := min(SRE.rep1.cfu,SRE.rep2.cfu),by=AAseq]

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#ERE
#make ~even bins on the basis of the minimum number of reads between replicates
bins.ERE.SR1 <- as.numeric(quantile(dt.SR1.coding[!is.na(ERE.rep1.class) & !is.na(ERE.rep2.class),ERE.cfu.min],seq(0,1,by=0.05)))
shared.strong.ERE.SR1 <- shared.positive.ERE.SR1 <- vector(mode="numeric",length=length(bins.ERE.SR1)-1)
for(i in 2:(length(bins.ERE.SR1)-1)){
  data <- dt.SR1.coding[ERE.cfu.min >= bins.ERE.SR1[i] & ERE.cfu.min < bins.ERE.SR1[i+1] & !is.na(dt.SR1.coding$ERE.rep1.class) & !is.na(dt.SR1.coding$ERE.rep2.class),]
  shared.strong.ERE.SR1[i] <- (nrow(data[data$ERE.rep1.class %in% c("strong") & data$ERE.rep2.class %in% c("strong"),])/nrow(data[data$ERE.rep1.class %in% c("strong"),]) + nrow(data[data$ERE.rep1.class %in% c("strong") & data$ERE.rep2.class %in% c("strong"),])/nrow(data[data$ERE.rep2.class %in% c("strong"),]))/2
  shared.positive.ERE.SR1[i] <- (nrow(data[data$ERE.rep1.class %in% c("weak","strong") & data$ERE.rep2.class %in% c("weak","strong"),])/nrow(data[data$ERE.rep1.class %in% c("weak","strong"),]) + nrow(data[data$ERE.rep1.class %in% c("weak","strong") & data$ERE.rep2.class %in% c("weak","strong"),])/nrow(data[data$ERE.rep2.class %in% c("weak","strong"),]))/2
}

#repeat SRE
bins.SRE.SR1 <- as.numeric(quantile(dt.SR1.coding[!is.na(SRE.rep1.class) & !is.na(SRE.rep2.class),SRE.cfu.min],seq(0,1,by=0.05)))
shared.strong.SRE.SR1 <- shared.positive.SRE.SR1 <- vector(mode="numeric",length=length(bins.SRE.SR1)-1)
for(i in 2:(length(bins.SRE.SR1)-1)){
  data <- dt.SR1.coding[SRE.cfu.min >= bins.SRE.SR1[i] & SRE.cfu.min < bins.SRE.SR1[i+1] & !is.na(dt.SR1.coding$SRE.rep1.class) & !is.na(dt.SR1.coding$SRE.rep2.class),]
  shared.strong.SRE.SR1[i] <- (nrow(data[data$SRE.rep1.class %in% c("strong") & data$SRE.rep2.class %in% c("strong"),])/nrow(data[data$SRE.rep1.class %in% c("strong"),]) + nrow(data[data$SRE.rep1.class %in% c("strong") & data$SRE.rep2.class %in% c("strong"),])/nrow(data[data$SRE.rep2.class %in% c("strong"),]))/2
  shared.positive.SRE.SR1[i] <- (nrow(data[data$SRE.rep1.class %in% c("weak","strong") & data$SRE.rep2.class %in% c("weak","strong"),])/nrow(data[data$SRE.rep1.class %in% c("weak","strong"),]) + nrow(data[data$SRE.rep1.class %in% c("weak","strong") & data$SRE.rep2.class %in% c("weak","strong"),])/nrow(data[data$SRE.rep2.class %in% c("weak","strong"),]))/2
}

pdf(file="./5_assess-replicates_out/SR1_shared-proportion-positive_v_min-cfu-binned-max-min.pdf",width=8,height=4.5)
par(mfrow=c(1,2))
plot(bins.ERE.SR1[1:20],shared.strong.ERE.SR1,type="l",lwd=2,log="x",col=rgb(100,69,155,maxColorValue=255),xlab="minimum number of reads, binned",ylab="average shared proportion strong btwn reps",xlim=c(1,1000),main="",ylim=c(0,1))
points(bins.SRE.SR1[1:20],shared.strong.SRE.SR1,type="l",lty=1,lwd=2,col=rgb(5,172,72,maxColorValue=255))
abline(v=15,lty=2)
legend("bottomright",legend=c("ERE","SRE"),lty=c(1,1),lwd=2,col=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),bty="n")

plot(bins.ERE.SR1[1:20],shared.positive.ERE.SR1,type="l",lwd=2,log="x",col=rgb(100,69,155,maxColorValue=255),xlab="minimum number of reads, binned",ylab="average shared proportion positive btwn reps",xlim=c(1,1000),main="",ylim=c(0,1))
points(bins.SRE.SR1[1:20],shared.positive.SRE.SR1,type="l",lty=1,lwd=2,col=rgb(5,172,72,maxColorValue=255))
abline(v=15,lty=2)
legend("bottomright",legend=c("ERE","SRE"),lty=c(1,1),lwd=2,col=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),bty="n")
dev.off()

#censor all classifications based on less than 15 cfu
dt.SR1[ERE.rep1.cfu<15, ERE.rep1.class := as.factor(NA)];dt.SR1[ERE.rep2.cfu<15, ERE.rep2.class := as.factor(NA)];dt.SR1[ERE.pooled.cfu<15, ERE.pooled.class := as.factor(NA)]
dt.SR1[SRE.rep1.cfu<15, SRE.rep1.class := as.factor(NA)];dt.SR1[SRE.rep2.cfu<15, SRE.rep2.class := as.factor(NA)];dt.SR1[SRE.pooled.cfu<15, SRE.pooled.class := as.factor(NA)]

#output table
write.table(dt.SR1[,.(AAseq,ERE.rep1.b1,ERE.rep1.b2,ERE.rep1.b3,ERE.rep1.b4,ERE.rep2.b1,ERE.rep2.b2,ERE.rep2.b3,ERE.rep2.b4,SRE.rep1.b1,SRE.rep1.b2,SRE.rep1.b3,SRE.rep1.b4,SRE.rep2.b1,SRE.rep2.b2,SRE.rep2.b3,SRE.rep2.b4,ERE.rep1.cfu,ERE.rep2.cfu,ERE.pooled.cfu,SRE.rep1.cfu,SRE.rep2.cfu,SRE.pooled.cfu,ERE.rep1.meanF,ERE.rep2.meanF,ERE.pooled.meanF,SRE.rep1.meanF,SRE.rep2.meanF,SRE.pooled.meanF,ERE.rep1.class,ERE.rep2.class,ERE.pooled.class,SRE.rep1.class,SRE.rep2.class,SRE.pooled.class)],file="5_assess-replicates_out/dt_output_SR1_censor15.csv",col.names=T,row.names=F,quote=F,sep=",")

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#plot by specificity
pdf(file="5_assess-replicates_out/SR1_pooled_ERE-meanF_v_SRE-meanF_points-colored-by-specificity.pdf",5.5,6,useDingbats=F)
plot(dt.SR1.coding$SRE.pooled.meanF,dt.SR1.coding$ERE.pooled.meanF,pch="",xlab="SRE meanF (a.u.), pooled reps",ylab="ERE meanF (a.u.), pooled reps",xlim=c(3,10),ylim=c(3,10))
points(dt.SR1.coding[ERE.pooled.class %in% c("null","weak") & SRE.pooled.class %in% c("null","weak"),SRE.pooled.meanF],dt.SR1.coding[ERE.pooled.class %in% c("null","weak") & SRE.pooled.class %in% c("null","weak"),ERE.pooled.meanF],pch=20,cex=0.6)
points(dt.SR1.stop[!is.na(as.character(ERE.pooled.class)) & !is.na(as.character(SRE.pooled.class)),SRE.pooled.meanF],dt.SR1.stop[!is.na(as.character(ERE.pooled.class)) & !is.na(as.character(SRE.pooled.class)),ERE.pooled.meanF],pch=20,cex=0.6,col="coral")
points(dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="null",SRE.pooled.meanF],dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="null",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(100,69,155,maxColorValue=255))
points(dt.SR1.coding[ERE.pooled.class=="null" & SRE.pooled.class=="strong",SRE.pooled.meanF],dt.SR1.coding[ERE.pooled.class=="null" & SRE.pooled.class=="strong",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(5,172,72,maxColorValue=255))
points(dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="weak",SRE.pooled.meanF],dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="weak",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(111,204,221,maxColorValue=255))
points(dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="strong",SRE.pooled.meanF],dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="strong",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(111,204,221,maxColorValue=255))
points(dt.SR1.coding[ERE.pooled.class=="weak" & SRE.pooled.class=="strong",SRE.pooled.meanF],dt.SR1.coding[ERE.pooled.class=="weak" & SRE.pooled.class=="strong",ERE.pooled.meanF],pch=20,cex=0.6,col=rgb(111,204,221,maxColorValue=255))
abline(v=dt.11P["GSKV",SRE.pooled.meanF],lwd=2,lty=2,col=rgb(5,172,72,maxColorValue=255))
abline(h=7.104872,lwd=2,lty=2,col=rgb(100,69,155,maxColorValue=255)) #mean(ERE.pooled.egka)
dev.off()

#correlation in meanF for variants with >15 cfu cutoff, ERE
pdf(file="5_assess-replicates_out/SR1_ERE_meanF-correlation.pdf",6.25,6.75,useDingbats=F)
plot(dt.SR1.coding[ERE.cfu.min>15,ERE.rep1.meanF],dt.SR1.coding[ERE.cfu.min>15,ERE.rep2.meanF],pch=19,xlab="ERE mean fluorescence, replicate 1 (a.u.)",ylab="ERE mean fluorescence, replicate 2 (a.u.)",cex.lab=1.4)
points(dt.SR1.stop[ERE.cfu.min>15,ERE.rep1.meanF],dt.SR1.stop[ERE.cfu.min>15,ERE.rep2.meanF],pch=19,col="coral")
summary(lm(dt.SR1.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep2.meanF]~dt.SR1.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep1.meanF]))#;abline(lm(dt.SR1.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep2.meanF]~dt.SR1.coding[ERE.cfu.min>15 & (ERE.rep1.class %in% c("weak","strong") | ERE.rep2.class %in% c("weak","strong")),ERE.rep1.meanF]),lty=2,lwd=3,col="red")
legend("topleft",legend="R-squared_pos: 0.50",bty="n",cex=0.8)
dev.off()

#correlation in meanF for variants with >15 cfu cutoff, SRE
pdf(file="5_assess-replicates_out/SR1_SRE_meanF-correlation.pdf",6.25,6.75,useDingbats=F)
plot(dt.SR1.coding[SRE.cfu.min>15,SRE.rep1.meanF],dt.SR1.coding[SRE.cfu.min>15,SRE.rep2.meanF],pch=19,xlab="SRE mean fluorescence, replicate 1 (a.u.)",ylab="SRE mean fluorescence, replicate 2 (a.u.)",cex.lab=1.4)
points(dt.SR1.stop[SRE.cfu.min>15,SRE.rep1.meanF],dt.SR1.stop[SRE.cfu.min>15,SRE.rep2.meanF],pch=19,col="coral")
summary(lm(dt.SR1.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep2.meanF]~dt.SR1.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep1.meanF]))#;abline(lm(dt.SR1.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep2.meanF]~dt.SR1.coding[SRE.cfu.min>15 & (SRE.rep1.class %in% c("weak","strong") | SRE.rep2.class %in% c("weak","strong")),SRE.rep1.meanF]),lty=2,lwd=3,col="red")
legend("topleft",legend="R-squared_pos: 0.69",bty="n",cex=0.8)
dev.off()