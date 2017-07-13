#16 Jan 2017
#TNS
#script to output alignments for weblogo creation; logistic regressions for determinants of binding specificity; analyze foldx output data

library(seqinr)
library(data.table)
library(forestplot)
library(ggplot2)

setwd("path/to/source/directory") #setwd to source file directory

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

#output list of all ERE-specific, SRE-specific variants in each background --> make logo plots with Weblogo
write(as.character(dt.11P.coding$AAseq[which(dt.11P.coding$specificity %in% c("ERE-specific"))]), file="8_bchem-determ_out/fasta/ERE_specific_11P.txt",sep="/n")
write(as.character(dt.11P.coding$AAseq[which(dt.11P.coding$specificity %in% c("SRE-specific"))]), file="8_bchem-determ_out/fasta/SRE_specific_11P.txt",sep="/n")
write(as.character(dt.SR1.coding$AAseq[which(dt.SR1.coding$specificity %in% c("ERE-specific"))]), file="8_bchem-determ_out/fasta/ERE_specific_SR1.txt",sep="/n")
write(as.character(dt.SR1.coding$AAseq[which(dt.SR1.coding$specificity %in% c("SRE-specific"))]), file="8_bchem-determ_out/fasta/SRE_specific_SR1.txt",sep="/n")
# #opened up individual .txt files, added > signifier to make fasta for logo generation with regexp, generated logos with Weblogo given settings described in screenshot in input directory

################################################################################

#bchem determinants regression using a set of biochemical/biophysical parameters that have been used previously in the literature (Ariata, Palzkill, Dal Peraro 2015, PLOS ONE e0118684), shown to have reduced correlation between stats
AA.props <- read.csv(file="8_bchem-determ_in/AA-properties.csv",header=T,stringsAsFactors=F)[1:20,c("AA","Volume","Hydrophobicity","Isoelectric.point","P.helix")]
#get variables ready for a logistic regression: center and normalize continuous variables
AA.props$Volume <- (AA.props$Volume-mean(AA.props$Volume))/sd(AA.props$Volume)
AA.props$Hydrophobicity <- (AA.props$Hydrophobicity-mean(AA.props$Hydrophobicity))/sd(AA.props$Hydrophobicity)
AA.props$Isoelectric.point <- (AA.props$Isoelectric.point-mean(AA.props$Isoelectric.point))/sd(AA.props$Isoelectric.point)
AA.props$P.helix <- (AA.props$P.helix-mean(AA.props$P.helix))/sd(AA.props$P.helix)

#construct tables for logistic regression
dt.ERE.SR1 <- dt.SR1.coding[,.(AAseq,AA1,AA2,AA3,AA4,meanF=ERE.pooled.meanF,class=ERE.full.class,sp.class=specificity)]
dt.ERE.SR1[,c("AA1.volume","AA1.hydrophobicity","AA1.isoelectric.point","AA1.P.helix",
                "AA2.volume","AA2.hydrophobicity","AA2.isoelectric.point","AA2.P.helix",
                "AA3.volume","AA3.hydrophobicity","AA3.isoelectric.point","AA3.P.helix",
                "AA4.volume","AA4.hydrophobicity","AA4.isoelectric.point","AA4.P.helix") :=
              list(AA.props[AA.props$AA==AA1,"Volume"],AA.props[AA.props$AA==AA1,"Hydrophobicity"],AA.props[AA.props$AA==AA1,"Isoelectric.point"],AA.props[AA.props$AA==AA1,"P.helix"],
                   AA.props[AA.props$AA==AA2,"Volume"],AA.props[AA.props$AA==AA2,"Hydrophobicity"],AA.props[AA.props$AA==AA2,"Isoelectric.point"],AA.props[AA.props$AA==AA2,"P.helix"],
                   AA.props[AA.props$AA==AA3,"Volume"],AA.props[AA.props$AA==AA3,"Hydrophobicity"],AA.props[AA.props$AA==AA3,"Isoelectric.point"],AA.props[AA.props$AA==AA3,"P.helix"],
                   AA.props[AA.props$AA==AA4,"Volume"],AA.props[AA.props$AA==AA4,"Hydrophobicity"],AA.props[AA.props$AA==AA4,"Isoelectric.point"],AA.props[AA.props$AA==AA4,"P.helix"]),by=AAseq]

dt.ERE.11P <- dt.ERE.SR1
dt.ERE.11P$meanF <- dt.11P.coding$ERE.pooled.meanF
dt.ERE.11P$class <- dt.11P.coding$ERE.full.class
dt.ERE.11P$sp.class <- dt.11P.coding$specificity

dt.SRE.SR1 <- dt.SR1.coding[,.(AAseq,AA1,AA2,AA3,AA4,meanF=SRE.pooled.meanF,class=SRE.full.class,sp.class=specificity)]
dt.SRE.SR1[,c("AA1.volume","AA1.hydrophobicity","AA1.isoelectric.point","AA1.P.helix",
              "AA2.volume","AA2.hydrophobicity","AA2.isoelectric.point","AA2.P.helix",
              "AA3.volume","AA3.hydrophobicity","AA3.isoelectric.point","AA3.P.helix",
              "AA4.volume","AA4.hydrophobicity","AA4.isoelectric.point","AA4.P.helix") :=
             list(AA.props[AA.props$AA==AA1,"Volume"],AA.props[AA.props$AA==AA1,"Hydrophobicity"],AA.props[AA.props$AA==AA1,"Isoelectric.point"],AA.props[AA.props$AA==AA1,"P.helix"],
                  AA.props[AA.props$AA==AA2,"Volume"],AA.props[AA.props$AA==AA2,"Hydrophobicity"],AA.props[AA.props$AA==AA2,"Isoelectric.point"],AA.props[AA.props$AA==AA2,"P.helix"],
                  AA.props[AA.props$AA==AA3,"Volume"],AA.props[AA.props$AA==AA3,"Hydrophobicity"],AA.props[AA.props$AA==AA3,"Isoelectric.point"],AA.props[AA.props$AA==AA3,"P.helix"],
                  AA.props[AA.props$AA==AA4,"Volume"],AA.props[AA.props$AA==AA4,"Hydrophobicity"],AA.props[AA.props$AA==AA4,"Isoelectric.point"],AA.props[AA.props$AA==AA4,"P.helix"]),by=AAseq]

dt.SRE.11P <- dt.SRE.SR1
dt.SRE.11P$meanF <- dt.11P.coding$SRE.pooled.meanF
dt.SRE.11P$class <- dt.11P.coding$SRE.full.class
dt.SRE.11P$sp.class <- dt.11P.coding$specificity

dt.ERE.SR1[,specific := 0];dt.ERE.SR1[sp.class == "ERE-specific",specific:=1];dt.ERE.SR1$specific <- as.factor(dt.ERE.SR1$specific)

dt.ERE.11P[,specific := 0];dt.ERE.11P[sp.class == "ERE-specific",specific:=1];dt.ERE.11P$specific <- as.factor(dt.ERE.11P$specific)

dt.SRE.SR1[,specific := 0];dt.SRE.SR1[sp.class == "SRE-specific",specific:=1];dt.SRE.SR1$specific <- as.factor(dt.SRE.SR1$specific)

dt.SRE.11P[,specific := 0];dt.SRE.11P[sp.class == "SRE-specific",specific:=1];dt.SRE.11P$specific <- as.factor(dt.SRE.11P$specific)

#define RE factor: 0 if ERE, 1 if SRE
dt.ERE.SR1[,RE := 0]
dt.ERE.11P[,RE := 0]
dt.SRE.SR1[,RE := 1]
dt.SRE.11P[,RE := 1]

#define perm factor: 0 if SR1 background, 1 if 11P
dt.ERE.SR1[,perm := 0]
dt.ERE.11P[,perm := 1]
dt.SRE.SR1[,perm := 0]
dt.SRE.11P[,perm := 1]

#logistic regressions for specific/not specific as a function of underlying biophysical properties
logit.ERE.SR1.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                                AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                                AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                                AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix,
                              data=dt.ERE.SR1,family="binomial")

logit.ERE.11P.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                                AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                                AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                                AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix,
                              data=dt.ERE.11P,family="binomial")

#to test for coefficients that significantly differ in SR1 and 11P background, built joint model with *perm interaction term
dt.ERE <- rbind(dt.ERE.SR1, dt.ERE.11P)
dt.ERE$perm <- as.factor(dt.ERE$perm)

logit.ERE.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                            AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                            AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                            AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix+
                            AA1.volume*perm+AA1.hydrophobicity*perm+AA1.isoelectric.point*perm+AA1.P.helix*perm+
                            AA2.volume*perm+AA2.hydrophobicity*perm+AA2.isoelectric.point*perm+AA2.P.helix*perm+
                            AA3.volume*perm+AA3.hydrophobicity*perm+AA3.isoelectric.point*perm+AA3.P.helix*perm+
                            AA4.volume*perm+AA4.hydrophobicity*perm+AA4.isoelectric.point*perm+AA4.P.helix*perm,
                          data=dt.ERE,family="binomial")

logit.SRE.SR1.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                                  AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                                  AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                                  AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix,
                                data=dt.SRE.SR1,family="binomial")

logit.SRE.11P.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                                  AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                                  AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                                  AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix,
                                data=dt.SRE.11P,family="binomial")

#to test for coefficients that significantly differ in SR1 and 11P background, built joint model with *perm interaction term
dt.SRE <- rbind(dt.SRE.SR1, dt.SRE.11P)
dt.SRE$perm <- as.factor(dt.SRE$perm)

logit.SRE.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                              AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                              AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                              AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix+
                              AA1.volume*perm+AA1.hydrophobicity*perm+AA1.isoelectric.point*perm+AA1.P.helix*perm+
                              AA2.volume*perm+AA2.hydrophobicity*perm+AA2.isoelectric.point*perm+AA2.P.helix*perm+
                              AA3.volume*perm+AA3.hydrophobicity*perm+AA3.isoelectric.point*perm+AA3.P.helix*perm+
                              AA4.volume*perm+AA4.hydrophobicity*perm+AA4.isoelectric.point*perm+AA4.P.helix*perm,
                            data=dt.SRE,family="binomial")

#to test for coefficients that significantly differ for ERE- veruss SRE-specificity within a background, built joint model with *RE interaction term
dt.SR1 <- rbind(dt.ERE.SR1,dt.SRE.SR1)
dt.SR1$RE <- as.factor(dt.SR1$RE)

dt.11P <- rbind(dt.ERE.11P,dt.SRE.11P)
dt.11P$RE <- as.factor(dt.11P$RE)

logit.SR1.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                            AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                            AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                            AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix+
                            AA1.volume*RE+AA1.hydrophobicity*RE+AA1.isoelectric.point*RE+AA1.P.helix*RE+
                            AA2.volume*RE+AA2.hydrophobicity*RE+AA2.isoelectric.point*RE+AA2.P.helix*RE+
                            AA3.volume*RE+AA3.hydrophobicity*RE+AA3.isoelectric.point*RE+AA3.P.helix*RE+
                            AA4.volume*RE+AA4.hydrophobicity*RE+AA4.isoelectric.point*RE+AA4.P.helix*RE,
                          data=dt.SR1,family="binomial")

logit.11P.specific <- glm(specific ~ AA1.volume+AA1.hydrophobicity+AA1.isoelectric.point+AA1.P.helix+
                            AA2.volume+AA2.hydrophobicity+AA2.isoelectric.point+AA2.P.helix+
                            AA3.volume+AA3.hydrophobicity+AA3.isoelectric.point+AA3.P.helix+
                            AA4.volume+AA4.hydrophobicity+AA4.isoelectric.point+AA4.P.helix+
                            AA1.volume*RE+AA1.hydrophobicity*RE+AA1.isoelectric.point*RE+AA1.P.helix*RE+
                            AA2.volume*RE+AA2.hydrophobicity*RE+AA2.isoelectric.point*RE+AA2.P.helix*RE+
                            AA3.volume*RE+AA3.hydrophobicity*RE+AA3.isoelectric.point*RE+AA3.P.helix*RE+
                            AA4.volume*RE+AA4.hydrophobicity*RE+AA4.isoelectric.point*RE+AA4.P.helix*RE,
                          data=dt.11P,family="binomial")

coefs <- names(summary(logit.ERE.SR1.specific)$coefficients[,"Estimate"])

mean.ERE.SR1.specific <- as.vector(summary(logit.ERE.SR1.specific)$coefficients[,"Estimate"])
lower.ERE.SR1.specific <- as.vector(summary(logit.ERE.SR1.specific)$coefficients[,"Estimate"])-1.96*as.vector(summary(logit.ERE.SR1.specific)$coefficients[,"Std. Error"])
upper.ERE.SR1.specific <- as.vector(summary(logit.ERE.SR1.specific)$coefficients[,"Estimate"])+1.96*as.vector(summary(logit.ERE.SR1.specific)$coefficients[,"Std. Error"])

mean.ERE.11P.specific <- as.vector(summary(logit.ERE.11P.specific)$coefficients[,"Estimate"])
lower.ERE.11P.specific <- as.vector(summary(logit.ERE.11P.specific)$coefficients[,"Estimate"])-1.96*as.vector(summary(logit.ERE.11P.specific)$coefficients[,"Std. Error"])
upper.ERE.11P.specific <- as.vector(summary(logit.ERE.11P.specific)$coefficients[,"Estimate"])+1.96*as.vector(summary(logit.ERE.11P.specific)$coefficients[,"Std. Error"])

pdf(file="8_bchem-determ_out/forest-plots_ERE-specific_SR1-v-11P.pdf",height=5,width=5)
forestplot(labeltext=list(coefs[-1]),mean=cbind((mean.ERE.SR1.specific[-1]),(mean.ERE.11P.specific[-1])),lower=cbind((lower.ERE.SR1.specific[-1]),(lower.ERE.11P.specific[-1])),upper=cbind((upper.ERE.SR1.specific[-1]),(upper.ERE.11P.specific[-1])),xlog=F,xlab="log2(odds ratio)",lwd.zero=1,lwd.ci=1,boxsize=0.3,col=fpColors(box=c("yellow","orange"),lines=c("black","black"),zero="black"),clip=c(-8,8),xticks=seq(-8,8,by=2),new_page=F)
dev.off()

#significant differences in coefficeints, ERE specific, SR1 versus 11P given in:
summary(logit.ERE.specific)

mean.SRE.SR1.specific <- as.vector(summary(logit.SRE.SR1.specific)$coefficients[,"Estimate"])
lower.SRE.SR1.specific <- as.vector(summary(logit.SRE.SR1.specific)$coefficients[,"Estimate"])-1.96*as.vector(summary(logit.SRE.SR1.specific)$coefficients[,"Std. Error"])
upper.SRE.SR1.specific <- as.vector(summary(logit.SRE.SR1.specific)$coefficients[,"Estimate"])+1.96*as.vector(summary(logit.SRE.SR1.specific)$coefficients[,"Std. Error"])

mean.SRE.11P.specific <- as.vector(summary(logit.SRE.11P.specific)$coefficients[,"Estimate"])
lower.SRE.11P.specific <- as.vector(summary(logit.SRE.11P.specific)$coefficients[,"Estimate"])-1.96*as.vector(summary(logit.SRE.11P.specific)$coefficients[,"Std. Error"])
upper.SRE.11P.specific <- as.vector(summary(logit.SRE.11P.specific)$coefficients[,"Estimate"])+1.96*as.vector(summary(logit.SRE.11P.specific)$coefficients[,"Std. Error"])

pdf(file="8_bchem-determ_out/forest-plots_SRE-specific_SR1-v-11P.pdf",height=5,width=5)
forestplot(labeltext=list(coefs[-1]),mean=cbind((mean.SRE.SR1.specific[-1]),(mean.SRE.11P.specific[-1])),lower=cbind((lower.SRE.SR1.specific[-1]),(lower.SRE.11P.specific[-1])),upper=cbind((upper.SRE.SR1.specific[-1]),(upper.SRE.11P.specific[-1])),xlog=F,xlab="log2(odds ratio)",lwd.zero=1,lwd.ci=1,boxsize=0.3,col=fpColors(box=c("yellow","orange"),lines=c("black","black"),zero="black"),clip=c(-8,8),xticks=seq(-8,8,by=2),new_page=F)
dev.off()

#significant differences in coefficeints, SRE specific, SR1 versus 11P given in:
summary(logit.SRE.specific)


#forest plots for ERE versus SRE in each of SR1 and 11P backgrounds
pdf(file="8_bchem-determ_out/forest-plots_SR1_ERE-specific-v-SRE-specific.pdf",height=5,width=5)
forestplot(labeltext=list(coefs[-1]),mean=cbind((mean.ERE.SR1.specific[-1]),(mean.SRE.SR1.specific[-1])),lower=cbind((lower.ERE.SR1.specific[-1]),(lower.SRE.SR1.specific[-1])),upper=cbind((upper.ERE.SR1.specific[-1]),(upper.SRE.SR1.specific[-1])),xlog=F,xlab="log2(odds ratio)",lwd.zero=1,lwd.ci=1,boxsize=0.3,col=fpColors(box=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),lines=c("black","black"),zero="black"),clip=c(-8,8),xticks=seq(-8,8,by=2),new_page=F)
dev.off()

#significant differences in coefficeints, SR1 background, ERE specific versus SRE specific given in:
summary(logit.SR1.specific)

pdf(file="8_bchem-determ_out/forest-plots_11P_ERE-specific-v-SRE-specific.pdf",height=5,width=5)
forestplot(labeltext=list(coefs[-1]),mean=cbind((mean.ERE.11P.specific[-1]),(mean.SRE.11P.specific[-1])),lower=cbind((lower.ERE.11P.specific[-1]),(lower.SRE.11P.specific[-1])),upper=cbind((upper.ERE.11P.specific[-1]),(upper.SRE.11P.specific[-1])),xlog=F,xlab="log2(odds ratio)",lwd.zero=1,lwd.ci=1,boxsize=0.3,col=fpColors(box=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),lines=c("black","black"),zero="black"),clip=c(-8,8),xticks=seq(-8,8,by=2),new_page=F)
dev.off()

#significant differences in coefficeints, 11P background, ERE specific versus SRE specific given in:
summary(logit.11P.specific)

#to have all on one scale, make one heatmap with all four regressions
values.all.specific <- data.frame(variable=rep(coefs[-1],2),set=c(rep("ERE, SR1",length(coefs)-1),rep("ERE, 11P",length(coefs)-1),rep("SRE, SR1",length(coefs)-1),rep("SRE, 11P",length(coefs)-1)),values=c(mean.ERE.SR1.specific[-1],mean.ERE.11P.specific[-1],mean.SRE.SR1.specific[-1],mean.SRE.11P.specific[-1]))
pdf(file="8_bchem-determ_out/heatmap_all-specific_raw.pdf",width=6,height=3.5)
p <- ggplot(values.all.specific,aes(variable,set))+geom_tile(aes(fill=values))+
  scale_fill_gradientn(colours=c("#A94E35","#F48365","#FFFFFF","#7378B9","#383C6C"),limits=c(-6,6),values=c(0,0.3387,0.5,0.66129,1))+
  labs(x="variable",y="background")+theme_classic()+
  coord_equal()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  scale_x_discrete(limits=values.all.specific$variable)
p
dev.off()

############################################################################################################################################

#example of epistasis: look at difference in logo for SRE-specific and SRE-speicifc with K3, R3: de-enrichment of K1 in +3 subsets
#read in fasta-formatted SRE.specific.11P file (after output above, must add > signifier with new line in fasta format manually)
SRE.specific.11P <- read.fasta("./8_bchem-determ_out/fasta/SRE_specific_11P.txt", seqtype="AA", as.string=FALSE)

SRE.specific.K3 <- vector()
for(i in 1:length(SRE.specific.11P)){
  if(SRE.specific.11P[i][[1]][[3]]=="K"){
    SRE.specific.K3 <- c(SRE.specific.K3,paste(SRE.specific.11P[i][[1]][[1]],SRE.specific.11P[i][[1]][[2]],SRE.specific.11P[i][[1]][[3]],SRE.specific.11P[i][[1]][[4]],sep=""))
  }
}
write(SRE.specific.K3, file="8_bchem-determ_out/fasta/SRE.specific.K3.txt",sep="/n")

SRE.specific.R3 <- vector()
for(i in 1:length(SRE.specific.11P)){
  if(SRE.specific.11P[i][[1]][[3]]=="R"){
    SRE.specific.R3 <- c(SRE.specific.R3,paste(SRE.specific.11P[i][[1]][[1]],SRE.specific.11P[i][[1]][[2]],SRE.specific.11P[i][[1]][[3]],SRE.specific.11P[i][[1]][[4]],sep=""))
  }
}
write(SRE.specific.R3, file="8_bchem-determ_out/fasta/SRE.specific.R3.txt",sep="/n")

#opened up SRE.speicfic.R3 and K3, added > with regexp to make fasta format, made logos

SRE.specific.K3 <- read.fasta("./8_bchem-determ_out/fasta/SRE.specific.K3.txt", seqtype="AA", as.string=FALSE)
SRE.specific.R3 <- read.fasta("./8_bchem-determ_out/fasta/SRE.specific.R3.txt", seqtype="AA", as.string=FALSE)

#define function to calculate enrichment/depletion of amino acid states in subsets (query) of a ref set
AAs <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
calc.enrichments <- function(query, ref, reps){
  query.prop.site <- data.frame(site=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20)), AA=rep(AAs,4))
  for(i in 1:4){
    observed.AAs <- vector()
    for(k in 1:length(query)){
      observed.AAs <- c(observed.AAs,query[k][[1]][[i]])
    }
    for(j in AAs){
      query.prop.site[query.prop.site$site==i & query.prop.site$AA==j,"proportion"] <- length(which(observed.AAs==j))/length(observed.AAs)
    }
  }
  boot.prop.site <- data.frame(site=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20)), AA=rep(AAs,4))
  for(b in 1:reps){
    boot <- sample(ref,length(query),replace=T)
    for(i in 1:4){
      observed.AAs <- vector()
      for(k in 1:length(boot)){
        observed.AAs <- c(observed.AAs,boot[k][[1]][[i]])
      }
      for(j in AAs){
        boot.prop.site[boot.prop.site$site==i & boot.prop.site$AA==j, b+2] <- length(which(observed.AAs==j))/length(observed.AAs)
      }
    }
  }
  for(i in 1:nrow(query.prop.site)){
    query.prop.site$p.enriched[i] <- length(which(boot.prop.site[i,-(1:2)]>=query.prop.site[i,"proportion"]))/length(boot.prop.site[i,-(1:2)])
    query.prop.site$p.depleted[i] <- length(which(boot.prop.site[i,-(1:2)]<=query.prop.site[i,"proportion"]))/length(boot.prop.site[i,-(1:2)])
  }
  return(query.prop.site)
}

SRE.sp.K3.v.all <- calc.enrichments(query=SRE.specific.K3,ref=SRE.specific.11P,reps=10000)
SRE.sp.R3.v.all <- calc.enrichments(query=SRE.specific.R3,ref=SRE.specific.11P,reps=10000)

SRE.sp.K3.v.all[SRE.sp.K3.v.all$site==3,"p.enriched"] <- NA;SRE.sp.K3.v.all[SRE.sp.K3.v.all$site==3,"p.depleted"] <- NA
SRE.sp.R3.v.all[SRE.sp.R3.v.all$site==3,"p.enriched"] <- NA;SRE.sp.R3.v.all[SRE.sp.R3.v.all$site==3,"p.depleted"] <- NA

SRE.sp.K3.v.all[SRE.sp.K3.v.all$p.enriched<(0.05/120) & !is.na(SRE.sp.K3.v.all$p.enriched),] #E4, Q4
SRE.sp.K3.v.all[SRE.sp.K3.v.all$p.depleted<(0.05/120) & !is.na(SRE.sp.K3.v.all$p.depleted),] #K1, R4

SRE.sp.R3.v.all[SRE.sp.R3.v.all$p.enriched<(0.05/120) & !is.na(SRE.sp.R3.v.all$p.enriched),] #none
SRE.sp.R3.v.all[SRE.sp.R3.v.all$p.depleted<(0.05/120) & !is.na(SRE.sp.R3.v.all$p.depleted),] #K1, A2

rm(list=ls())
#################################################################################
#used FoldX to model all SRE-strong binders on SRE in the AncSR1 and AncSR2 crystal structures
#read in data, ask whether 11P-dependent versus 11P-independent variants differ in their predicted binding energy in AncSR2 (carrying permissives) backgorund
#corroborate with raw meanF measurements

#read in data from modeling all SRE-sp RH combos in AncSR1 background and measuring predicted binding affinity for SRE via FoldX
SR1.SREsp.11P.dependent <- read.csv(file="./8_bchem-determ_in/SR1.SREsp.11P.dependent.csv", header=T)
SR1.SREsp.11P.independent <- read.csv(file="./8_bchem-determ_in/SR1.SREsp.11P.independent.csv", header=T)
#read in data from modeling all SRE-sp RH combos in AncSR2 background and measuring predicted binding affinity for SRE via FoldX
SR2.SREsp.11P.dependent <- read.csv(file="./8_bchem-determ_in/SR2.SREsp.11P.dependent.csv", header=T)
SR2.SREsp.11P.independent <- read.csv(file="./8_bchem-determ_in/SR2.SREsp.11P.independent.csv", header=T)

#plot total interaction energies on one figure output
pdf(file="8_bchem-determ_out/FoldX-summary_11P-dep-v-ind_SR1-and-SR2-backgrounds.pdf",4.5,4.5,useDingbats=F)
boxplot(SR1.SREsp.11P.dependent$Interaction.Energy, SR1.SREsp.11P.independent$Interaction.Energy,SR2.SREsp.11P.dependent$Interaction.Energy,SR2.SREsp.11P.independent$Interaction.Energy,add=F,lwd=2,notch=T,range=0,at=c(1,2,4,5),ylab="predicted binding energy (kcal/mol)",names=c("11P-dependent","11P-independent","11P-dependent","11P-independent"),col=c(rgb(241,101,33,maxColorValue=255),rgb(255,255,0,maxColorValue=255)))
dev.off()

qqnorm(SR1.SREsp.11P.dependent$Interaction.Energy)
qqnorm(SR1.SREsp.11P.independent$Interaction.Energy)
qqnorm(SR2.SREsp.11P.dependent$Interaction.Energy)
qqnorm(SR2.SREsp.11P.independent$Interaction.Energy)

wilcox.test(SR1.SREsp.11P.dependent$Interaction.Energy, SR1.SREsp.11P.independent$Interaction.Energy)
wilcox.test(SR2.SREsp.11P.dependent$Interaction.Energy, SR2.SREsp.11P.independent$Interaction.Energy)

#do 11P-dependent variants get a significantly greater affinity enhancement by 11P than 11P-dependent variants?
qqnorm(SR2.SREsp.11P.dependent$Interaction.Energy-SR1.SREsp.11P.dependent$Interaction.Energy)
qqnorm(SR2.SREsp.11P.independent$Interaction.Energy-SR1.SREsp.11P.independent$Interaction.Energy)

#is the difference in predicted binding energy in AncSR1 versus AncSR2 background different for 11P-dependent versus 11P-independent variants? Paired t.test
t.test(c(SR2.SREsp.11P.dependent$Interaction.Energy-SR1.SREsp.11P.dependent$Interaction.Energy),c(SR2.SREsp.11P.independent$Interaction.Energy-SR1.SREsp.11P.independent$Interaction.Energy))
#no significant difference in the gain in predicted binding energy due to 35 AncSR2 subs for 11P-dependent versus 11P-independent variants; that is, no evidence that 11P have an enhanced permissive effect for the 11P-dependent variantsr

#make equivalent figure, just from raw FACS-seq estimates?

dt.11P.coding <- read.table(file="./6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)

dt.SR1.coding <- read.table(file="./6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

for(i in 1:nrow(SR1.SREsp.11P.dependent)){
  if(dt.SR1.coding[as.character(SR1.SREsp.11P.dependent$AAseq[i]),SRE.pooled.cfu]>15){
    SR1.SREsp.11P.dependent$meanF[i] <- dt.SR1.coding[as.character(SR1.SREsp.11P.dependent$AAseq[i]),SRE.pooled.meanF]
  }else{
    SR1.SREsp.11P.dependent$meanF[i] <- NA
  }
}

for(i in 1:nrow(SR1.SREsp.11P.independent)){
  if(dt.SR1.coding[as.character(SR1.SREsp.11P.independent$AAseq[i]),SRE.pooled.cfu]>15){
    SR1.SREsp.11P.independent$meanF[i] <- dt.SR1.coding[as.character(SR1.SREsp.11P.independent$AAseq[i]),SRE.pooled.meanF]
  }else{
    SR1.SREsp.11P.independent$meanF[i] <- NA
  }
}


for(i in 1:nrow(SR2.SREsp.11P.dependent)){
  if(dt.11P.coding[as.character(SR2.SREsp.11P.dependent$AAseq[i]),SRE.pooled.cfu]>15){
    SR2.SREsp.11P.dependent$meanF[i] <- dt.11P.coding[as.character(SR2.SREsp.11P.dependent$AAseq[i]),SRE.pooled.meanF]
  }else{
    SR2.SREsp.11P.dependent$meanF[i] <- NA
  }
}

for(i in 1:nrow(SR2.SREsp.11P.independent)){
  if(dt.11P.coding[as.character(SR2.SREsp.11P.independent$AAseq[i]),SRE.pooled.cfu]>15){
    SR2.SREsp.11P.independent$meanF[i] <- dt.11P.coding[as.character(SR2.SREsp.11P.independent$AAseq[i]),SRE.pooled.meanF]
  }else{
    SR2.SREsp.11P.independent$meanF[i] <- NA
  }
}

pdf(file="8_bchem-determ_out/meanF-summary_11P-dep-v-ind_SR1-and-SR2-backgrounds_include-restricted.pdf",4.5,4.5,useDingbats=F)
boxplot(SR1.SREsp.11P.dependent$meanF, SR1.SREsp.11P.independent$meanF,SR2.SREsp.11P.dependent$meanF,SR2.SREsp.11P.independent$meanF,add=F,lwd=2,notch=T,range=0,at=c(1,2,4,5),ylab="FACS-seq meanF (a.u.)",names=c("11P-dependent","11P-independent","11P-dependent","11P-independent"),col=c(rgb(241,101,33,maxColorValue=255),rgb(255,255,0,maxColorValue=255)))
stripchart(c(SR2.SREsp.11P.independent[SR2.SREsp.11P.independent$AAseq %in% c("CARV","HARV","HPRM"),"meanF"]),vertical=T,method="jitter",jitter=0.15,pch=4,col="red",add=T,at=5)
stripchart(c(SR1.SREsp.11P.independent[SR1.SREsp.11P.independent$AAseq %in% c("CARV","HARV","HPRM"),"meanF"]),vertical=T,method="jitter",jitter=0.15,pch=4,col="red",add=T,at=2)
stripchart(c(SR2.SREsp.11P.independent[SR2.SREsp.11P.independent$AAseq %in% c("KASM"),"meanF"]),vertical=T,method="jitter",jitter=0.15,pch=15,col="blue",add=T,at=5)
stripchart(c(SR1.SREsp.11P.independent[SR1.SREsp.11P.independent$AAseq %in% c("KASM"),"meanF"]),vertical=T,method="jitter",jitter=0.15,pch=15,col="blue",add=T,at=2)
stripchart(dt.11P.coding["GSKV",SRE.pooled.meanF],vertical=T,method="jitter",jitter=0.15,pch=19,col="brown",add=T,at=4)
stripchart(dt.SR1.coding["GSKV",SRE.pooled.meanF],vertical=T,method="jitter",jitter=0.15,pch=19,col="brown",add=T,at=1)
stripchart(dt.11P.coding[c("YGKQ","SPKM"),SRE.pooled.meanF],vertical=T,method="jitter",jitter=0.15,pch=19,col="brown",add=T,at=4)
stripchart(dt.SR1.coding[c("YGKQ","SPKM"),SRE.pooled.meanF],vertical=T,method="jitter",jitter=0.15,pch=19,col="brown",add=T,at=1)
dev.off()

qqnorm(SR1.SREsp.11P.dependent$meanF)
qqnorm(SR1.SREsp.11P.independent$meanF)
qqnorm(SR2.SREsp.11P.dependent$meanF)
qqnorm(SR2.SREsp.11P.independent$meanF)
#definitely non-normal, as expected given floor and ceiling of quantitative range of the assay

wilcox.test(SR1.SREsp.11P.dependent$meanF, SR1.SREsp.11P.independent$meanF)#$p.value
wilcox.test(SR2.SREsp.11P.dependent$meanF,SR2.SREsp.11P.independent$meanF)

#re-do, but take out the three red x points, the putatively restricted genotypes; of these, all are predicted strong, FACS-seq not strong in AncSR1+11P; however, we re-tested these in isogenic cultures and they were false negatives in the FACS-seq assay, likely due to a noticeable detrimental effect on growth due to yEGFP overexpression
for(i in 1:nrow(SR2.SREsp.11P.independent)){
  if(dt.11P.coding[as.character(SR2.SREsp.11P.independent$AAseq[i]),SRE.pooled.cfu]>15 & !(as.character(SR2.SREsp.11P.independent$AAseq[i]) %in% c("CARV","HARV","HPRM"))){
    SR2.SREsp.11P.independent$meanF[i] <- dt.11P.coding[as.character(SR2.SREsp.11P.independent$AAseq[i]),SRE.pooled.meanF]
  }else{
    SR2.SREsp.11P.independent$meanF[i] <- NA
    #print(SR2.SREsp.11P.independent$AAseq[i])
  }
}

pdf(file="8_bchem-determ_out/meanF-summary_11P-dep-v-ind_SR1-and-SR2-backgrounds_exclude-restricted.pdf",4.5,4.5,useDingbats=F)
boxplot(SR1.SREsp.11P.dependent$meanF, SR1.SREsp.11P.independent$meanF,SR2.SREsp.11P.dependent$meanF,SR2.SREsp.11P.independent$meanF,add=F,lwd=2,notch=T,range=0,at=c(1,2,4,5),ylab="FACS-seq meanF (a.u.)",names=c("11P-dependent","11P-independent","11P-dependent","11P-independent"),col=c(rgb(241,101,33,maxColorValue=255),rgb(255,255,0,maxColorValue=255)))
dev.off()

wilcox.test(SR1.SREsp.11P.dependent$meanF, SR1.SREsp.11P.independent$meanF)
wilcox.test(SR2.SREsp.11P.dependent$meanF,SR2.SREsp.11P.independent$meanF)
