#16 Jan 2017
#TNS
#script to build models to predict missing variants

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)
library(Matrix)
library(glmnetcr)

dt.11P <- read.table(file="5_assess-replicates_out/dt_output_11P_censor15.csv",header=TRUE,sep=",");dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)
dt.11P[,c("AA1","AA2","AA3","AA4") := list(strsplit(as.character(AAseq),split="")[[1]][1],strsplit(as.character(AAseq),split="")[[1]][2],strsplit(as.character(AAseq),split="")[[1]][3],strsplit(as.character(AAseq),split="")[[1]][4]),by=AAseq]
dt.11P[,AA1 := as.factor(AA1)];dt.11P[,AA2 := as.factor(AA2)];dt.11P[,AA3 := as.factor(AA3)];dt.11P[,AA4 := as.factor(AA4)]

#set order of classification factors for ERE, SRE classes
dt.11P[,ERE.rep1.class:=factor(ERE.rep1.class,levels=c("null","weak","strong"))]
dt.11P[,ERE.rep2.class:=factor(ERE.rep2.class,levels=c("null","weak","strong"))]
dt.11P[,ERE.pooled.class:=factor(ERE.pooled.class,levels=c("null","weak","strong"))]
dt.11P[,SRE.rep1.class:=factor(SRE.rep1.class,levels=c("null","weak","strong"))]
dt.11P[,SRE.rep2.class:=factor(SRE.rep2.class,levels=c("null","weak","strong"))]
dt.11P[,SRE.pooled.class:=factor(SRE.pooled.class,levels=c("null","weak","strong"))]

dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]
rm(dt.11P)

##############################################################################################################
#use glmnet.cr to infer ordinal logistic model with penalization, including main effect and pairwise epistasis terms
dt.ERE.11P.pooled <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.pooled.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled,
                                  y=dt.11P.coding$ERE.pooled.class[!is.na(dt.11P.coding$ERE.pooled.class)],
                                  maxit=500)
save(dt.ERE.11P.pooled, file="6_predict-missing_out/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="6_predict-missing_out/glmnet.cr.ERE.11P.pooled.Rda")
rm(dt.ERE.11P.pooled)
rm(logit.ERE.11P.pooled)

dt.SRE.11P.pooled <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$SRE.pooled.class)])
logit.SRE.11P.pooled <- glmnet.cr(x=dt.SRE.11P.pooled,
                                  y=dt.11P.coding$SRE.pooled.class[!is.na(dt.11P.coding$SRE.pooled.class)],
                                  maxit=500)
save(dt.SRE.11P.pooled, file="6_predict-missing_out/dt.SRE.11P.pooled")
save(logit.SRE.11P.pooled, file="6_predict-missing_out/glmnet.cr.SRE.11P.pooled.Rda")
rm(dt.SRE.11P.pooled)
rm(logit.SRE.11P.pooled)

#####################################################################################################
#cross validation
#ERE pooled, 11P:
#get path of lambdas used in full data inference
load(file="./6_predict-missing_out/glmnet.cr.ERE.11P.pooled.Rda")
lambda <- logit.ERE.11P.pooled$lambda
rm(logit.ERE.11P.pooled)

#break data into ten groups of ~equal numbers
#randomly order rows in dt
set.seed(1001)
dt.ERE <- dt.11P.coding[!is.na(ERE.pooled.class)]
dt.ERE <- dt.ERE[sample(nrow(dt.ERE)),]
folds.ERE <- cut(seq(1,nrow(dt.ERE)),breaks=10,labels=F)
save(dt.ERE, file="./6_predict-missing_out/cross-validate_11P/ERE/dt.ERE.Rda")
save(folds.ERE, file="./6_predict-missing_out/cross-validate_11P/ERE/folds.ERE.Rda")

#dropping out test1 set, infer model
test1.i <- which(folds.ERE==1,arr.ind=T)
train1.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test1.i,])
train1.logit.ERE.11P.pooled <- glmnet.cr(x=train1.d,
                                            y=dt.ERE[-test1.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train1.d, file="6_predict-missing_out/cross-validate_11P/ERE/train1.d.Rda")
save(train1.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train1.logit.ERE.11P.pooled.Rda")
rm(train1.d)
rm(train1.logit.ERE.11P.pooled)
rm(test1.i)

#dropping out test2 set, infer model
test2.i <- which(folds.ERE==2,arr.ind=T)
train2.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test2.i,])
train2.logit.ERE.11P.pooled <- glmnet.cr(x=train2.d,
                                            y=dt.ERE[-test2.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train2.d, file="6_predict-missing_out/cross-validate_11P/ERE/train2.d.Rda")
save(train2.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train2.logit.ERE.11P.pooled.Rda")
rm(train2.d)
rm(train2.logit.ERE.11P.pooled)
rm(test2.i)

#dropping out test3 set, infer model
test3.i <- which(folds.ERE==3,arr.ind=T)
train3.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test3.i,])
train3.logit.ERE.11P.pooled <- glmnet.cr(x=train3.d,
                                            y=dt.ERE[-test3.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train3.d, file="6_predict-missing_out/cross-validate_11P/ERE/train3.d.Rda")
save(train3.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train3.logit.ERE.11P.pooled.Rda")
rm(train3.d)
rm(train3.logit.ERE.11P.pooled)
rm(test3.i)

#dropping out test4 set, infer model
test4.i <- which(folds.ERE==4,arr.ind=T)
train4.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test4.i,])
train4.logit.ERE.11P.pooled <- glmnet.cr(x=train4.d,
                                            y=dt.ERE[-test4.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train4.d, file="6_predict-missing_out/cross-validate_11P/ERE/train4.d.Rda")
save(train4.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train4.logit.ERE.11P.pooled.Rda")
rm(train4.d)
rm(train4.logit.ERE.11P.pooled)
rm(test4.i)

#dropping out test5 set, infer model
test5.i <- which(folds.ERE==5,arr.ind=T)
train5.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test5.i,])
train5.logit.ERE.11P.pooled <- glmnet.cr(x=train5.d,
                                            y=dt.ERE[-test5.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train5.d, file="6_predict-missing_out/cross-validate_11P/ERE/train5.d.Rda")
save(train5.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train5.logit.ERE.11P.pooled.Rda")
rm(train5.d)
rm(train5.logit.ERE.11P.pooled)
rm(test5.i)

#dropping out test6 set, infer model
test6.i <- which(folds.ERE==6,arr.ind=T)
train6.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test6.i,])
train6.logit.ERE.11P.pooled <- glmnet.cr(x=train6.d,
                                            y=dt.ERE[-test6.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train6.d, file="6_predict-missing_out/cross-validate_11P/ERE/train6.d.Rda")
save(train6.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train6.logit.ERE.11P.pooled.Rda")
rm(train6.d)
rm(train6.logit.ERE.11P.pooled)
rm(test6.i)

#dropping out test7 set, infer model
test7.i <- which(folds.ERE==7,arr.ind=T)
train7.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test7.i,])
train7.logit.ERE.11P.pooled <- glmnet.cr(x=train7.d,
                                            y=dt.ERE[-test7.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train7.d, file="6_predict-missing_out/cross-validate_11P/ERE/train7.d.Rda")
save(train7.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train7.logit.ERE.11P.pooled.Rda")
rm(train7.d)
rm(train7.logit.ERE.11P.pooled)
rm(test7.i)

#dropping out test8 set, infer model
test8.i <- which(folds.ERE==8,arr.ind=T)
train8.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test8.i,])
train8.logit.ERE.11P.pooled <- glmnet.cr(x=train8.d,
                                            y=dt.ERE[-test8.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train8.d, file="6_predict-missing_out/cross-validate_11P/ERE/train8.d.Rda")
save(train8.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train8.logit.ERE.11P.pooled.Rda")
rm(train8.d)
rm(train8.logit.ERE.11P.pooled)
rm(test8.i)

#dropping out test9 set, infer model
test9.i <- which(folds.ERE==9,arr.ind=T)
train9.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test9.i,])
train9.logit.ERE.11P.pooled <- glmnet.cr(x=train9.d,
                                            y=dt.ERE[-test9.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train9.d, file="6_predict-missing_out/cross-validate_11P/ERE/train9.d.Rda")
save(train9.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train9.logit.ERE.11P.pooled.Rda")
rm(train9.d)
rm(train9.logit.ERE.11P.pooled)
rm(test9.i)

#dropping out test10 set, infer model
test10.i <- which(folds.ERE==10,arr.ind=T)
train10.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test10.i,])
train10.logit.ERE.11P.pooled <- glmnet.cr(x=train10.d,
                                            y=dt.ERE[-test10.i,ERE.pooled.class],
                                            maxit=500,
                                            lambda=lambda)
save(train10.d, file="6_predict-missing_out/cross-validate_11P/ERE/train10.d.Rda")
save(train10.logit.ERE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/ERE/train10.logit.ERE.11P.pooled.Rda")
rm(train10.d)
rm(train10.logit.ERE.11P.pooled)
rm(test10.i)

rm(lambda)
rm(folds.ERE)
rm(dt.ERE)

######################################################################
#SRE pooled, 11P:
#get path of lambdas used in full data infSREnce
load(file="./6_predict-missing_out/glmnet.cr.SRE.11P.pooled.Rda")
lambda <- logit.SRE.11P.pooled$lambda
rm(logit.SRE.11P.pooled)

#break data into ten groups of ~equal numbers
#randomly order rows in dt
set.seed(1001)
dt.SRE <- dt.11P.coding[!is.na(SRE.pooled.class)]
dt.SRE <- dt.SRE[sample(nrow(dt.SRE)),]
folds.SRE <- cut(seq(1,nrow(dt.SRE)),breaks=10,labels=F)
save(dt.SRE, file="./6_predict-missing_out/cross-validate_11P/SRE/dt.SRE.Rda")
save(folds.SRE, file="./6_predict-missing_out/cross-validate_11P/SRE/folds.SRE.Rda")

#dropping out test1 set, infer model
test1.i <- which(folds.SRE==1,arr.ind=T)
train1.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test1.i,])
train1.logit.SRE.11P.pooled <- glmnet.cr(x=train1.d,
                                         y=dt.SRE[-test1.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train1.d, file="6_predict-missing_out/cross-validate_11P/SRE/train1.d.Rda")
save(train1.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train1.logit.SRE.11P.pooled.Rda")
rm(train1.d)
rm(train1.logit.SRE.11P.pooled)
rm(test1.i)

#dropping out test2 set, infer model
test2.i <- which(folds.SRE==2,arr.ind=T)
train2.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test2.i,])
train2.logit.SRE.11P.pooled <- glmnet.cr(x=train2.d,
                                         y=dt.SRE[-test2.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train2.d, file="6_predict-missing_out/cross-validate_11P/SRE/train2.d.Rda")
save(train2.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train2.logit.SRE.11P.pooled.Rda")
rm(train2.d)
rm(train2.logit.SRE.11P.pooled)
rm(test2.i)

#dropping out test3 set, infer model
test3.i <- which(folds.SRE==3,arr.ind=T)
train3.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test3.i,])
train3.logit.SRE.11P.pooled <- glmnet.cr(x=train3.d,
                                         y=dt.SRE[-test3.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train3.d, file="6_predict-missing_out/cross-validate_11P/SRE/train3.d.Rda")
save(train3.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train3.logit.SRE.11P.pooled.Rda")
rm(train3.d)
rm(train3.logit.SRE.11P.pooled)
rm(test3.i)

#dropping out test4 set, infer model
test4.i <- which(folds.SRE==4,arr.ind=T)
train4.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test4.i,])
train4.logit.SRE.11P.pooled <- glmnet.cr(x=train4.d,
                                         y=dt.SRE[-test4.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train4.d, file="6_predict-missing_out/cross-validate_11P/SRE/train4.d.Rda")
save(train4.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train4.logit.SRE.11P.pooled.Rda")
rm(train4.d)
rm(train4.logit.SRE.11P.pooled)
rm(test4.i)

#dropping out test5 set, infer model
test5.i <- which(folds.SRE==5,arr.ind=T)
train5.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test5.i,])
train5.logit.SRE.11P.pooled <- glmnet.cr(x=train5.d,
                                         y=dt.SRE[-test5.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train5.d, file="6_predict-missing_out/cross-validate_11P/SRE/train5.d.Rda")
save(train5.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train5.logit.SRE.11P.pooled.Rda")
rm(train5.d)
rm(train5.logit.SRE.11P.pooled)
rm(test5.i)

#dropping out test6 set, infer model
test6.i <- which(folds.SRE==6,arr.ind=T)
train6.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test6.i,])
train6.logit.SRE.11P.pooled <- glmnet.cr(x=train6.d,
                                         y=dt.SRE[-test6.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train6.d, file="6_predict-missing_out/cross-validate_11P/SRE/train6.d.Rda")
save(train6.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train6.logit.SRE.11P.pooled.Rda")
rm(train6.d)
rm(train6.logit.SRE.11P.pooled)
rm(test6.i)

#dropping out test7 set, infer model
test7.i <- which(folds.SRE==7,arr.ind=T)
train7.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test7.i,])
train7.logit.SRE.11P.pooled <- glmnet.cr(x=train7.d,
                                         y=dt.SRE[-test7.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train7.d, file="6_predict-missing_out/cross-validate_11P/SRE/train7.d.Rda")
save(train7.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train7.logit.SRE.11P.pooled.Rda")
rm(train7.d)
rm(train7.logit.SRE.11P.pooled)
rm(test7.i)

#dropping out test8 set, infer model
test8.i <- which(folds.SRE==8,arr.ind=T)
train8.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test8.i,])
train8.logit.SRE.11P.pooled <- glmnet.cr(x=train8.d,
                                         y=dt.SRE[-test8.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train8.d, file="6_predict-missing_out/cross-validate_11P/SRE/train8.d.Rda")
save(train8.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train8.logit.SRE.11P.pooled.Rda")
rm(train8.d)
rm(train8.logit.SRE.11P.pooled)
rm(test8.i)

#dropping out test9 set, infer model
test9.i <- which(folds.SRE==9,arr.ind=T)
train9.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test9.i,])
train9.logit.SRE.11P.pooled <- glmnet.cr(x=train9.d,
                                         y=dt.SRE[-test9.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train9.d, file="6_predict-missing_out/cross-validate_11P/SRE/train9.d.Rda")
save(train9.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train9.logit.SRE.11P.pooled.Rda")
rm(train9.d)
rm(train9.logit.SRE.11P.pooled)
rm(test9.i)

#dropping out test10 set, infer model
test10.i <- which(folds.SRE==10,arr.ind=T)
train10.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test10.i,])
train10.logit.SRE.11P.pooled <- glmnet.cr(x=train10.d,
                                          y=dt.SRE[-test10.i,SRE.pooled.class],
                                          maxit=500,
                                          lambda=lambda)
save(train10.d, file="6_predict-missing_out/cross-validate_11P/SRE/train10.d.Rda")
save(train10.logit.SRE.11P.pooled, file="6_predict-missing_out/cross-validate_11P/SRE/train10.logit.SRE.11P.pooled.Rda")
rm(train10.d)
rm(train10.logit.SRE.11P.pooled)
rm(test10.i)

rm(lambda)
rm(folds.SRE)
rm(dt.SRE)

##############################################################################################################
#look at CV stats as a function of lambda: 11P data
#ERE
load("./6_predict-missing_out/cross-validate_11P/ERE/dt.ERE.Rda");load("./6_predict-missing_out/cross-validate_11P/ERE/folds.ERE.Rda")

load(file="./6_predict-missing_out/glmnet.cr.ERE.11P.pooled.Rda")
lambda <- logit.ERE.11P.pooled$lambda
rm(logit.ERE.11P.pooled)

test1.i <- which(folds.ERE==1,arr.ind=T)
test1.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test1.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train1.logit.ERE.11P.pooled.Rda")
test1.p <- predict(train1.logit.ERE.11P.pooled,newx=test1.d)
test1.class <- dt.ERE[test1.i,ERE.pooled.class]
rm(test1.d);rm(train1.logit.ERE.11P.pooled)

test2.i <- which(folds.ERE==2,arr.ind=T)
test2.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test2.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train2.logit.ERE.11P.pooled.Rda")
test2.p <- predict(train2.logit.ERE.11P.pooled,newx=test2.d)
test2.class <- dt.ERE[test2.i,ERE.pooled.class]
rm(test2.d);rm(train2.logit.ERE.11P.pooled)

test3.i <- which(folds.ERE==3,arr.ind=T)
test3.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test3.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train3.logit.ERE.11P.pooled.Rda")
test3.p <- predict(train3.logit.ERE.11P.pooled,newx=test3.d)
test3.class <- dt.ERE[test3.i,ERE.pooled.class]
rm(test3.d);rm(train3.logit.ERE.11P.pooled)

test4.i <- which(folds.ERE==4,arr.ind=T)
test4.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test4.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train4.logit.ERE.11P.pooled.Rda")
test4.p <- predict(train4.logit.ERE.11P.pooled,newx=test4.d)
test4.class <- dt.ERE[test4.i,ERE.pooled.class]
rm(test4.d);rm(train4.logit.ERE.11P.pooled)

test5.i <- which(folds.ERE==5,arr.ind=T)
test5.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test5.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train5.logit.ERE.11P.pooled.Rda")
test5.p <- predict(train5.logit.ERE.11P.pooled,newx=test5.d)
test5.class <- dt.ERE[test5.i,ERE.pooled.class]
rm(test5.d);rm(train5.logit.ERE.11P.pooled)

test6.i <- which(folds.ERE==6,arr.ind=T)
test6.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test6.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train6.logit.ERE.11P.pooled.Rda")
test6.p <- predict(train6.logit.ERE.11P.pooled,newx=test6.d)
test6.class <- dt.ERE[test6.i,ERE.pooled.class]
rm(test6.d);rm(train6.logit.ERE.11P.pooled)

test7.i <- which(folds.ERE==7,arr.ind=T)
test7.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test7.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train7.logit.ERE.11P.pooled.Rda")
test7.p <- predict(train7.logit.ERE.11P.pooled,newx=test7.d)
test7.class <- dt.ERE[test7.i,ERE.pooled.class]
rm(test7.d);rm(train7.logit.ERE.11P.pooled)

test8.i <- which(folds.ERE==8,arr.ind=T)
test8.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test8.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train8.logit.ERE.11P.pooled.Rda")
test8.p <- predict(train8.logit.ERE.11P.pooled,newx=test8.d)
test8.class <- dt.ERE[test8.i,ERE.pooled.class]
rm(test8.d);rm(train8.logit.ERE.11P.pooled)

test9.i <- which(folds.ERE==9,arr.ind=T)
test9.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test9.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train9.logit.ERE.11P.pooled.Rda")
test9.p <- predict(train9.logit.ERE.11P.pooled,newx=test9.d)
test9.class <- dt.ERE[test9.i,ERE.pooled.class]
rm(test9.d);rm(train9.logit.ERE.11P.pooled)

test10.i <- which(folds.ERE==10,arr.ind=T)
test10.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test10.i,])
load("./6_predict-missing_out/cross-validate_11P/ERE/train10.logit.ERE.11P.pooled.Rda")
test10.p <- predict(train10.logit.ERE.11P.pooled,newx=test10.d)
test10.class <- dt.ERE[test10.i,ERE.pooled.class]
rm(test10.d);rm(train10.logit.ERE.11P.pooled)

save(test1.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test1.p.Rda");save(test2.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test2.p.Rda");save(test3.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test3.p.Rda");save(test4.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test4.p.Rda");save(test5.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test5.p.Rda");save(test6.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test6.p.Rda");save(test7.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test7.p.Rda");save(test8.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test8.p.Rda");save(test9.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test9.p.Rda");save(test10.p, file="./6_predict-missing_out/cross-validate_11P/ERE/test10.p.Rda")
save(test1.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test1.class.Rda");save(test2.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test2.class.Rda");save(test3.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test3.class.Rda");save(test4.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test4.class.Rda");save(test5.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test5.class.Rda");save(test6.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test6.class.Rda");save(test7.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test7.class.Rda");save(test8.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test8.class.Rda");save(test9.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test9.class.Rda");save(test10.class, file="./6_predict-missing_out/cross-validate_11P/ERE/test10.class.Rda")
save(test1.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test1.i.Rda");save(test2.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test2.i.Rda");save(test3.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test3.i.Rda");save(test4.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test4.i.Rda");save(test5.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test5.i.Rda");save(test6.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test6.i.Rda");save(test7.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test7.i.Rda");save(test8.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test8.i.Rda");save(test9.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test9.i.Rda");save(test10.i, file="./6_predict-missing_out/cross-validate_11P/ERE/test10.i.Rda")
save(lambda, file="./6_predict-missing_out/cross-validate_11P/ERE/lambda.Rda")

#record TPR (proportion of true positives called positive in prediction) and PPV (proportion of predicted positives that are true positives)
test1.TPR.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.TPR.strong <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.strong <- vector(mode="numeric",length=ncol(test1.p$class))
for(i in 1:ncol(test1.p$class)){
  model.pred <- test1.p$class[,i]
  test1.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class%in% c("weak","strong"))/sum(test1.class %in% c("weak","strong"));test1.TPR.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(test1.class=="strong")
  test1.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test1.PPV.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(model.pred=="strong")
}

test2.TPR.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.TPR.strong <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.strong <- vector(mode="numeric",length=ncol(test2.p$class))
for(i in 1:ncol(test2.p$class)){
  model.pred <- test2.p$class[,i]
  test2.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class%in% c("weak","strong"))/sum(test2.class %in% c("weak","strong"));test2.TPR.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(test2.class=="strong")
  test2.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test2.PPV.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(model.pred=="strong")
}

test3.TPR.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.TPR.strong <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.strong <- vector(mode="numeric",length=ncol(test3.p$class))
for(i in 1:ncol(test3.p$class)){
  model.pred <- test3.p$class[,i]
  test3.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class%in% c("weak","strong"))/sum(test3.class %in% c("weak","strong"));test3.TPR.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(test3.class=="strong")
  test3.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test3.PPV.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(model.pred=="strong")
}

test4.TPR.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.TPR.strong <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.strong <- vector(mode="numeric",length=ncol(test4.p$class))
for(i in 1:ncol(test4.p$class)){
  model.pred <- test4.p$class[,i]
  test4.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class%in% c("weak","strong"))/sum(test4.class %in% c("weak","strong"));test4.TPR.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(test4.class=="strong")
  test4.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test4.PPV.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(model.pred=="strong")
}

test5.TPR.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.TPR.strong <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.strong <- vector(mode="numeric",length=ncol(test5.p$class))
for(i in 1:ncol(test5.p$class)){
  model.pred <- test5.p$class[,i]
  test5.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class%in% c("weak","strong"))/sum(test5.class %in% c("weak","strong"));test5.TPR.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(test5.class=="strong")
  test5.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test5.PPV.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(model.pred=="strong")
}

test6.TPR.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.TPR.strong <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.strong <- vector(mode="numeric",length=ncol(test6.p$class))
for(i in 1:ncol(test6.p$class)){
  model.pred <- test6.p$class[,i]
  test6.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class%in% c("weak","strong"))/sum(test6.class %in% c("weak","strong"));test6.TPR.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(test6.class=="strong")
  test6.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test6.PPV.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(model.pred=="strong")
}

test7.TPR.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.TPR.strong <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.strong <- vector(mode="numeric",length=ncol(test7.p$class))
for(i in 1:ncol(test7.p$class)){
  model.pred <- test7.p$class[,i]
  test7.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class%in% c("weak","strong"))/sum(test7.class %in% c("weak","strong"));test7.TPR.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(test7.class=="strong")
  test7.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test7.PPV.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(model.pred=="strong")
}

test8.TPR.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.TPR.strong <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.strong <- vector(mode="numeric",length=ncol(test8.p$class))
for(i in 1:ncol(test8.p$class)){
  model.pred <- test8.p$class[,i]
  test8.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class%in% c("weak","strong"))/sum(test8.class %in% c("weak","strong"));test8.TPR.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(test8.class=="strong")
  test8.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test8.PPV.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(model.pred=="strong")
}

test9.TPR.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.TPR.strong <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.strong <- vector(mode="numeric",length=ncol(test9.p$class))
for(i in 1:ncol(test9.p$class)){
  model.pred <- test9.p$class[,i]
  test9.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class%in% c("weak","strong"))/sum(test9.class %in% c("weak","strong"));test9.TPR.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(test9.class=="strong")
  test9.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test9.PPV.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(model.pred=="strong")
}

test10.TPR.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.TPR.strong <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.strong <- vector(mode="numeric",length=ncol(test10.p$class))
for(i in 1:ncol(test10.p$class)){
  model.pred <- test10.p$class[,i]
  test10.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class%in% c("weak","strong"))/sum(test10.class %in% c("weak","strong"));test10.TPR.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(test10.class=="strong")
  test10.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test10.PPV.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(model.pred=="strong")
}

TPR.strong <- vector(mode="numeric",length=length(test1.TPR.strong));TPR.pos <- vector(mode="numeric",length=length(test1.TPR.pos));PPV.strong <- vector(mode="numeric",length=length(test1.PPV.strong));PPV.pos <- vector(mode="numeric",length=length(test1.PPV.pos))
for(i in 1:length(test1.TPR.strong)){
  TPR.strong[i] <- mean(c(test1.TPR.strong[i],test2.TPR.strong[i],test3.TPR.strong[i],test4.TPR.strong[i],test5.TPR.strong[i],test6.TPR.strong[i],test7.TPR.strong[i],test8.TPR.strong[i],test9.TPR.strong[i],test10.TPR.strong[i]),na.rm=T)
  PPV.strong[i] <- mean(c(test1.PPV.strong[i],test2.PPV.strong[i],test3.PPV.strong[i],test4.PPV.strong[i],test5.PPV.strong[i],test6.PPV.strong[i],test7.PPV.strong[i],test8.PPV.strong[i],test9.PPV.strong[i],test10.PPV.strong[i]),na.rm=T)
  TPR.pos[i] <- mean(c(test1.TPR.pos[i],test2.TPR.pos[i],test3.TPR.pos[i],test4.TPR.pos[i],test5.TPR.pos[i],test6.TPR.pos[i],test7.TPR.pos[i],test8.TPR.pos[i],test9.TPR.pos[i],test10.TPR.pos[i]),na.rm=T)
  PPV.pos[i] <- mean(c(test1.PPV.pos[i],test2.PPV.pos[i],test3.PPV.pos[i],test4.PPV.pos[i],test5.PPV.pos[i],test6.PPV.pos[i],test7.PPV.pos[i],test8.PPV.pos[i],test9.PPV.pos[i],test10.PPV.pos[i]),na.rm=T)
}

pdf(file="./6_predict-missing_out/cross-validate_11P/ERE/TPR-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.TPR.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, strong",main="ERE, 11P")
points(lambda,test2.TPR.strong,type="l",lty=3)
points(lambda,test3.TPR.strong,type="l",lty=3)
points(lambda,test4.TPR.strong,type="l",lty=3)
points(lambda,test5.TPR.strong,type="l",lty=3)
points(lambda,test6.TPR.strong,type="l",lty=3)
points(lambda,test7.TPR.strong,type="l",lty=3)
points(lambda,test8.TPR.strong,type="l",lty=3)
points(lambda,test9.TPR.strong,type="l",lty=3)
points(lambda,test10.TPR.strong,type="l",lty=3)
points(lambda,TPR.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_11P/ERE/TPR-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.TPR.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, positive",main="ERE, 11P")
points(lambda,test2.TPR.pos,type="l",lty=3)
points(lambda,test3.TPR.pos,type="l",lty=3)
points(lambda,test4.TPR.pos,type="l",lty=3)
points(lambda,test5.TPR.pos,type="l",lty=3)
points(lambda,test6.TPR.pos,type="l",lty=3)
points(lambda,test7.TPR.pos,type="l",lty=3)
points(lambda,test8.TPR.pos,type="l",lty=3)
points(lambda,test9.TPR.pos,type="l",lty=3)
points(lambda,test10.TPR.pos,type="l",lty=3)
points(lambda,TPR.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_11P/ERE/PPV-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.PPV.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, strong",main="ERE, 11P")
points(lambda,test2.PPV.strong,type="l",lty=3)
points(lambda,test3.PPV.strong,type="l",lty=3)
points(lambda,test4.PPV.strong,type="l",lty=3)
points(lambda,test5.PPV.strong,type="l",lty=3)
points(lambda,test6.PPV.strong,type="l",lty=3)
points(lambda,test7.PPV.strong,type="l",lty=3)
points(lambda,test8.PPV.strong,type="l",lty=3)
points(lambda,test9.PPV.strong,type="l",lty=3)
points(lambda,test10.PPV.strong,type="l",lty=3)
points(lambda,PPV.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_11P/ERE/PPV-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.PPV.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, positive",main="ERE, 11P")
points(lambda,test2.PPV.pos,type="l",lty=3)
points(lambda,test3.PPV.pos,type="l",lty=3)
points(lambda,test4.PPV.pos,type="l",lty=3)
points(lambda,test5.PPV.pos,type="l",lty=3)
points(lambda,test6.PPV.pos,type="l",lty=3)
points(lambda,test7.PPV.pos,type="l",lty=3)
points(lambda,test8.PPV.pos,type="l",lty=3)
points(lambda,test9.PPV.pos,type="l",lty=3)
points(lambda,test10.PPV.pos,type="l",lty=3)
points(lambda,PPV.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

step.ERE.11P <- which(abs(lambda-(1e-5))==min(abs(lambda-(1e-5))))
save(step.ERE.11P,file="./6_predict-missing_out/step.ERE.11P.Rda")

lambda[step.ERE.11P] 

#bring together all cross-validation predictions: in dt.ERE, add column giving prediction.CV -- use this to compute PPV and TPR, since above these shouldn't be averaged completely, should they? different subsets have different #'s of positives etc.
dt.ERE[test1.i,prediction.CV:=test1.p$class[,step.ERE.11P]]
dt.ERE[test2.i,prediction.CV:=test2.p$class[,step.ERE.11P]]
dt.ERE[test3.i,prediction.CV:=test3.p$class[,step.ERE.11P]]
dt.ERE[test4.i,prediction.CV:=test4.p$class[,step.ERE.11P]]
dt.ERE[test5.i,prediction.CV:=test5.p$class[,step.ERE.11P]]
dt.ERE[test6.i,prediction.CV:=test6.p$class[,step.ERE.11P]]
dt.ERE[test7.i,prediction.CV:=test7.p$class[,step.ERE.11P]]
dt.ERE[test8.i,prediction.CV:=test8.p$class[,step.ERE.11P]]
dt.ERE[test9.i,prediction.CV:=test9.p$class[,step.ERE.11P]]
dt.ERE[test10.i,prediction.CV:=test10.p$class[,step.ERE.11P]]

write.table(data.frame(table(dt.ERE$ERE.pooled.class,dt.ERE$prediction.CV,useNA="always",dnn=c("expt","pred"))),file="6_predict-missing_out/cross-validate_11P/ERE/CV-summary.txt")

#SRE
load("./6_predict-missing_out/cross-validate_11P/SRE/dt.SRE.Rda");load("./6_predict-missing_out/cross-validate_11P/SRE/folds.SRE.Rda")

load(file="./6_predict-missing_out/glmnet.cr.SRE.11P.pooled.Rda")
lambda <- logit.SRE.11P.pooled$lambda
rm(logit.SRE.11P.pooled)

test1.i <- which(folds.SRE==1,arr.ind=T)
test1.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test1.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train1.logit.SRE.11P.pooled.Rda")
test1.p <- predict(train1.logit.SRE.11P.pooled,newx=test1.d)
test1.class <- dt.SRE[test1.i,SRE.pooled.class]
rm(test1.d);rm(train1.logit.SRE.11P.pooled)

test2.i <- which(folds.SRE==2,arr.ind=T)
test2.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test2.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train2.logit.SRE.11P.pooled.Rda")
test2.p <- predict(train2.logit.SRE.11P.pooled,newx=test2.d)
test2.class <- dt.SRE[test2.i,SRE.pooled.class]
rm(test2.d);rm(train2.logit.SRE.11P.pooled)

test3.i <- which(folds.SRE==3,arr.ind=T)
test3.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test3.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train3.logit.SRE.11P.pooled.Rda")
test3.p <- predict(train3.logit.SRE.11P.pooled,newx=test3.d)
test3.class <- dt.SRE[test3.i,SRE.pooled.class]
rm(test3.d);rm(train3.logit.SRE.11P.pooled)

test4.i <- which(folds.SRE==4,arr.ind=T)
test4.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test4.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train4.logit.SRE.11P.pooled.Rda")
test4.p <- predict(train4.logit.SRE.11P.pooled,newx=test4.d)
test4.class <- dt.SRE[test4.i,SRE.pooled.class]
rm(test4.d);rm(train4.logit.SRE.11P.pooled)

test5.i <- which(folds.SRE==5,arr.ind=T)
test5.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test5.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train5.logit.SRE.11P.pooled.Rda")
test5.p <- predict(train5.logit.SRE.11P.pooled,newx=test5.d)
test5.class <- dt.SRE[test5.i,SRE.pooled.class]
rm(test5.d);rm(train5.logit.SRE.11P.pooled)

test6.i <- which(folds.SRE==6,arr.ind=T)
test6.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test6.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train6.logit.SRE.11P.pooled.Rda")
test6.p <- predict(train6.logit.SRE.11P.pooled,newx=test6.d)
test6.class <- dt.SRE[test6.i,SRE.pooled.class]
rm(test6.d);rm(train6.logit.SRE.11P.pooled)

test7.i <- which(folds.SRE==7,arr.ind=T)
test7.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test7.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train7.logit.SRE.11P.pooled.Rda")
test7.p <- predict(train7.logit.SRE.11P.pooled,newx=test7.d)
test7.class <- dt.SRE[test7.i,SRE.pooled.class]
rm(test7.d);rm(train7.logit.SRE.11P.pooled)

test8.i <- which(folds.SRE==8,arr.ind=T)
test8.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test8.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train8.logit.SRE.11P.pooled.Rda")
test8.p <- predict(train8.logit.SRE.11P.pooled,newx=test8.d)
test8.class <- dt.SRE[test8.i,SRE.pooled.class]
rm(test8.d);rm(train8.logit.SRE.11P.pooled)

test9.i <- which(folds.SRE==9,arr.ind=T)
test9.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test9.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train9.logit.SRE.11P.pooled.Rda")
test9.p <- predict(train9.logit.SRE.11P.pooled,newx=test9.d)
test9.class <- dt.SRE[test9.i,SRE.pooled.class]
rm(test9.d);rm(train9.logit.SRE.11P.pooled)

test10.i <- which(folds.SRE==10,arr.ind=T)
test10.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test10.i,])
load("./6_predict-missing_out/cross-validate_11P/SRE/train10.logit.SRE.11P.pooled.Rda")
test10.p <- predict(train10.logit.SRE.11P.pooled,newx=test10.d)
test10.class <- dt.SRE[test10.i,SRE.pooled.class]
rm(test10.d);rm(train10.logit.SRE.11P.pooled)

save(test1.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test1.p.Rda");save(test2.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test2.p.Rda");save(test3.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test3.p.Rda");save(test4.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test4.p.Rda");save(test5.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test5.p.Rda");save(test6.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test6.p.Rda");save(test7.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test7.p.Rda");save(test8.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test8.p.Rda");save(test9.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test9.p.Rda");save(test10.p, file="./6_predict-missing_out/cross-validate_11P/SRE/test10.p.Rda")
save(test1.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test1.class.Rda");save(test2.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test2.class.Rda");save(test3.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test3.class.Rda");save(test4.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test4.class.Rda");save(test5.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test5.class.Rda");save(test6.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test6.class.Rda");save(test7.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test7.class.Rda");save(test8.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test8.class.Rda");save(test9.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test9.class.Rda");save(test10.class, file="./6_predict-missing_out/cross-validate_11P/SRE/test10.class.Rda")
save(test1.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test1.i.Rda");save(test2.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test2.i.Rda");save(test3.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test3.i.Rda");save(test4.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test4.i.Rda");save(test5.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test5.i.Rda");save(test6.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test6.i.Rda");save(test7.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test7.i.Rda");save(test8.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test8.i.Rda");save(test9.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test9.i.Rda");save(test10.i, file="./6_predict-missing_out/cross-validate_11P/SRE/test10.i.Rda")
save(lambda, file="./6_predict-missing_out/cross-validate_11P/SRE/lambda.Rda")

#record TPR (proportion of true positives called positive in prediction) and PPV (proportion of predicted positives that are true positives)
test1.TPR.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.TPR.strong <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.strong <- vector(mode="numeric",length=ncol(test1.p$class))
for(i in 1:ncol(test1.p$class)){
  model.pred <- test1.p$class[,i]
  test1.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class%in% c("weak","strong"))/sum(test1.class %in% c("weak","strong"));test1.TPR.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(test1.class=="strong")
  test1.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test1.PPV.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(model.pred=="strong")
}

test2.TPR.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.TPR.strong <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.strong <- vector(mode="numeric",length=ncol(test2.p$class))
for(i in 1:ncol(test2.p$class)){
  model.pred <- test2.p$class[,i]
  test2.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class%in% c("weak","strong"))/sum(test2.class %in% c("weak","strong"));test2.TPR.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(test2.class=="strong")
  test2.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test2.PPV.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(model.pred=="strong")
}

test3.TPR.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.TPR.strong <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.strong <- vector(mode="numeric",length=ncol(test3.p$class))
for(i in 1:ncol(test3.p$class)){
  model.pred <- test3.p$class[,i]
  test3.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class%in% c("weak","strong"))/sum(test3.class %in% c("weak","strong"));test3.TPR.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(test3.class=="strong")
  test3.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test3.PPV.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(model.pred=="strong")
}

test4.TPR.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.TPR.strong <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.strong <- vector(mode="numeric",length=ncol(test4.p$class))
for(i in 1:ncol(test4.p$class)){
  model.pred <- test4.p$class[,i]
  test4.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class%in% c("weak","strong"))/sum(test4.class %in% c("weak","strong"));test4.TPR.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(test4.class=="strong")
  test4.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test4.PPV.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(model.pred=="strong")
}

test5.TPR.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.TPR.strong <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.strong <- vector(mode="numeric",length=ncol(test5.p$class))
for(i in 1:ncol(test5.p$class)){
  model.pred <- test5.p$class[,i]
  test5.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class%in% c("weak","strong"))/sum(test5.class %in% c("weak","strong"));test5.TPR.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(test5.class=="strong")
  test5.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test5.PPV.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(model.pred=="strong")
}

test6.TPR.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.TPR.strong <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.strong <- vector(mode="numeric",length=ncol(test6.p$class))
for(i in 1:ncol(test6.p$class)){
  model.pred <- test6.p$class[,i]
  test6.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class%in% c("weak","strong"))/sum(test6.class %in% c("weak","strong"));test6.TPR.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(test6.class=="strong")
  test6.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test6.PPV.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(model.pred=="strong")
}

test7.TPR.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.TPR.strong <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.strong <- vector(mode="numeric",length=ncol(test7.p$class))
for(i in 1:ncol(test7.p$class)){
  model.pred <- test7.p$class[,i]
  test7.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class%in% c("weak","strong"))/sum(test7.class %in% c("weak","strong"));test7.TPR.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(test7.class=="strong")
  test7.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test7.PPV.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(model.pred=="strong")
}

test8.TPR.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.TPR.strong <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.strong <- vector(mode="numeric",length=ncol(test8.p$class))
for(i in 1:ncol(test8.p$class)){
  model.pred <- test8.p$class[,i]
  test8.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class%in% c("weak","strong"))/sum(test8.class %in% c("weak","strong"));test8.TPR.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(test8.class=="strong")
  test8.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test8.PPV.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(model.pred=="strong")
}

test9.TPR.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.TPR.strong <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.strong <- vector(mode="numeric",length=ncol(test9.p$class))
for(i in 1:ncol(test9.p$class)){
  model.pred <- test9.p$class[,i]
  test9.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class%in% c("weak","strong"))/sum(test9.class %in% c("weak","strong"));test9.TPR.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(test9.class=="strong")
  test9.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test9.PPV.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(model.pred=="strong")
}

test10.TPR.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.TPR.strong <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.strong <- vector(mode="numeric",length=ncol(test10.p$class))
for(i in 1:ncol(test10.p$class)){
  model.pred <- test10.p$class[,i]
  test10.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class%in% c("weak","strong"))/sum(test10.class %in% c("weak","strong"));test10.TPR.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(test10.class=="strong")
  test10.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test10.PPV.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(model.pred=="strong")
}

TPR.strong <- vector(mode="numeric",length=length(test1.TPR.strong));TPR.pos <- vector(mode="numeric",length=length(test1.TPR.pos));PPV.strong <- vector(mode="numeric",length=length(test1.PPV.strong));PPV.pos <- vector(mode="numeric",length=length(test1.PPV.pos))
for(i in 1:length(test1.TPR.strong)){
  TPR.strong[i] <- mean(c(test1.TPR.strong[i],test2.TPR.strong[i],test3.TPR.strong[i],test4.TPR.strong[i],test5.TPR.strong[i],test6.TPR.strong[i],test7.TPR.strong[i],test8.TPR.strong[i],test9.TPR.strong[i],test10.TPR.strong[i]),na.rm=T)
  PPV.strong[i] <- mean(c(test1.PPV.strong[i],test2.PPV.strong[i],test3.PPV.strong[i],test4.PPV.strong[i],test5.PPV.strong[i],test6.PPV.strong[i],test7.PPV.strong[i],test8.PPV.strong[i],test9.PPV.strong[i],test10.PPV.strong[i]),na.rm=T)
  TPR.pos[i] <- mean(c(test1.TPR.pos[i],test2.TPR.pos[i],test3.TPR.pos[i],test4.TPR.pos[i],test5.TPR.pos[i],test6.TPR.pos[i],test7.TPR.pos[i],test8.TPR.pos[i],test9.TPR.pos[i],test10.TPR.pos[i]),na.rm=T)
  PPV.pos[i] <- mean(c(test1.PPV.pos[i],test2.PPV.pos[i],test3.PPV.pos[i],test4.PPV.pos[i],test5.PPV.pos[i],test6.PPV.pos[i],test7.PPV.pos[i],test8.PPV.pos[i],test9.PPV.pos[i],test10.PPV.pos[i]),na.rm=T)
}

pdf(file="./6_predict-missing_out/cross-validate_11P/SRE/TPR-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.TPR.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, strong",main="SRE, 11P")
points(lambda,test2.TPR.strong,type="l",lty=3)
points(lambda,test3.TPR.strong,type="l",lty=3)
points(lambda,test4.TPR.strong,type="l",lty=3)
points(lambda,test5.TPR.strong,type="l",lty=3)
points(lambda,test6.TPR.strong,type="l",lty=3)
points(lambda,test7.TPR.strong,type="l",lty=3)
points(lambda,test8.TPR.strong,type="l",lty=3)
points(lambda,test9.TPR.strong,type="l",lty=3)
points(lambda,test10.TPR.strong,type="l",lty=3)
points(lambda,TPR.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_11P/SRE/TPR-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.TPR.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, positive",main="SRE, 11P")
points(lambda,test2.TPR.pos,type="l",lty=3)
points(lambda,test3.TPR.pos,type="l",lty=3)
points(lambda,test4.TPR.pos,type="l",lty=3)
points(lambda,test5.TPR.pos,type="l",lty=3)
points(lambda,test6.TPR.pos,type="l",lty=3)
points(lambda,test7.TPR.pos,type="l",lty=3)
points(lambda,test8.TPR.pos,type="l",lty=3)
points(lambda,test9.TPR.pos,type="l",lty=3)
points(lambda,test10.TPR.pos,type="l",lty=3)
points(lambda,TPR.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_11P/SRE/PPV-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.PPV.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, strong",main="SRE, 11P")
points(lambda,test2.PPV.strong,type="l",lty=3)
points(lambda,test3.PPV.strong,type="l",lty=3)
points(lambda,test4.PPV.strong,type="l",lty=3)
points(lambda,test5.PPV.strong,type="l",lty=3)
points(lambda,test6.PPV.strong,type="l",lty=3)
points(lambda,test7.PPV.strong,type="l",lty=3)
points(lambda,test8.PPV.strong,type="l",lty=3)
points(lambda,test9.PPV.strong,type="l",lty=3)
points(lambda,test10.PPV.strong,type="l",lty=3)
points(lambda,PPV.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_11P/SRE/PPV-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.PPV.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, positive",main="SRE, 11P")
points(lambda,test2.PPV.pos,type="l",lty=3)
points(lambda,test3.PPV.pos,type="l",lty=3)
points(lambda,test4.PPV.pos,type="l",lty=3)
points(lambda,test5.PPV.pos,type="l",lty=3)
points(lambda,test6.PPV.pos,type="l",lty=3)
points(lambda,test7.PPV.pos,type="l",lty=3)
points(lambda,test8.PPV.pos,type="l",lty=3)
points(lambda,test9.PPV.pos,type="l",lty=3)
points(lambda,test10.PPV.pos,type="l",lty=3)
points(lambda,PPV.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

step.SRE.11P <- which(abs(lambda-(1e-5))==min(abs(lambda-(1e-5))))
save(step.SRE.11P,file="./6_predict-missing_out/step.SRE.11P.Rda")

lambda[step.SRE.11P] 

#bring together all cross-validation predictions: in dt.SRE, add column giving prediction.CV -- use this to compute PPV and TPR, since above these shouldn't be averaged completely, should they? diffSREnt subsets have diffSREnt #'s of positives etc.
dt.SRE[test1.i,prediction.CV:=test1.p$class[,step.SRE.11P]]
dt.SRE[test2.i,prediction.CV:=test2.p$class[,step.SRE.11P]]
dt.SRE[test3.i,prediction.CV:=test3.p$class[,step.SRE.11P]]
dt.SRE[test4.i,prediction.CV:=test4.p$class[,step.SRE.11P]]
dt.SRE[test5.i,prediction.CV:=test5.p$class[,step.SRE.11P]]
dt.SRE[test6.i,prediction.CV:=test6.p$class[,step.SRE.11P]]
dt.SRE[test7.i,prediction.CV:=test7.p$class[,step.SRE.11P]]
dt.SRE[test8.i,prediction.CV:=test8.p$class[,step.SRE.11P]]
dt.SRE[test9.i,prediction.CV:=test9.p$class[,step.SRE.11P]]
dt.SRE[test10.i,prediction.CV:=test10.p$class[,step.SRE.11P]]

write.table(data.frame(table(dt.SRE$SRE.pooled.class,dt.SRE$prediction.CV,useNA="always",dnn=c("expt","pred"))),file="6_predict-missing_out/cross-validate_11P/SRE/CV-summary.txt")

rm(list=ls())

##############################################################################################################
#repeat for SR1 data
dt.SR1 <- read.table(file="5_assess-replicates_out/dt_output_SR1_censor15.csv",header=TRUE,sep=",");dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)
dt.SR1[,c("AA1","AA2","AA3","AA4") := list(strsplit(as.character(AAseq),split="")[[1]][1],strsplit(as.character(AAseq),split="")[[1]][2],strsplit(as.character(AAseq),split="")[[1]][3],strsplit(as.character(AAseq),split="")[[1]][4]),by=AAseq]
dt.SR1[,AA1 := as.factor(AA1)];dt.SR1[,AA2 := as.factor(AA2)];dt.SR1[,AA3 := as.factor(AA3)];dt.SR1[,AA4 := as.factor(AA4)]

#set order of classification factors for ERE, SRE classes
dt.SR1[,ERE.rep1.class:=factor(ERE.rep1.class,levels=c("null","weak","strong"))]
dt.SR1[,ERE.rep2.class:=factor(ERE.rep2.class,levels=c("null","weak","strong"))]
dt.SR1[,ERE.pooled.class:=factor(ERE.pooled.class,levels=c("null","weak","strong"))]
dt.SR1[,SRE.rep1.class:=factor(SRE.rep1.class,levels=c("null","weak","strong"))]
dt.SR1[,SRE.rep2.class:=factor(SRE.rep2.class,levels=c("null","weak","strong"))]
dt.SR1[,SRE.pooled.class:=factor(SRE.pooled.class,levels=c("null","weak","strong"))]

dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]
rm(dt.SR1)

##############################################################################################################
#use glmnet.cr to infer ordinal logistic model with penalization, including main effect and pairwise epistasis terms
dt.ERE.SR1.pooled <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.pooled.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled,
                                  y=dt.SR1.coding$ERE.pooled.class[!is.na(dt.SR1.coding$ERE.pooled.class)],
                                  maxit=500)
save(dt.ERE.SR1.pooled, file="6_predict-missing_out/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="6_predict-missing_out/glmnet.cr.ERE.SR1.pooled.Rda")
rm(dt.ERE.SR1.pooled)
rm(logit.ERE.SR1.pooled)

dt.SRE.SR1.pooled <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$SRE.pooled.class)])
logit.SRE.SR1.pooled <- glmnet.cr(x=dt.SRE.SR1.pooled,
                                  y=dt.SR1.coding$SRE.pooled.class[!is.na(dt.SR1.coding$SRE.pooled.class)],
                                  maxit=500)
save(dt.SRE.SR1.pooled, file="6_predict-missing_out/dt.SRE.SR1.pooled")
save(logit.SRE.SR1.pooled, file="6_predict-missing_out/glmnet.cr.SRE.SR1.pooled.Rda")
rm(dt.SRE.SR1.pooled)
rm(logit.SRE.SR1.pooled)

#####################################################################################################
#cross validation
#ERE pooled, SR1:
#get path of lambdas used in full data inference
load(file="./6_predict-missing_out/glmnet.cr.ERE.SR1.pooled.Rda")
lambda <- logit.ERE.SR1.pooled$lambda
rm(logit.ERE.SR1.pooled)

#break data into ten groups of ~equal numbers
#randomly order rows in dt
set.seed(1003)
dt.ERE <- dt.SR1.coding[!is.na(ERE.pooled.class)]
dt.ERE <- dt.ERE[sample(nrow(dt.ERE)),]
folds.ERE <- cut(seq(1,nrow(dt.ERE)),breaks=10,labels=F)
save(dt.ERE, file="./6_predict-missing_out/cross-validate_SR1/ERE/dt.ERE.Rda")
save(folds.ERE, file="./6_predict-missing_out/cross-validate_SR1/ERE/folds.ERE.Rda")

#dropping out test1 set, infer model
test1.i <- which(folds.ERE==1,arr.ind=T)
train1.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test1.i,])
train1.logit.ERE.SR1.pooled <- glmnet.cr(x=train1.d,
                                         y=dt.ERE[-test1.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train1.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train1.d.Rda")
save(train1.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train1.logit.ERE.SR1.pooled.Rda")
rm(train1.d)
rm(train1.logit.ERE.SR1.pooled)
rm(test1.i)

#dropping out test2 set, infer model
test2.i <- which(folds.ERE==2,arr.ind=T)
train2.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test2.i,])
train2.logit.ERE.SR1.pooled <- glmnet.cr(x=train2.d,
                                         y=dt.ERE[-test2.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train2.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train2.d.Rda")
save(train2.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train2.logit.ERE.SR1.pooled.Rda")
rm(train2.d)
rm(train2.logit.ERE.SR1.pooled)
rm(test2.i)

#dropping out test3 set, infer model
test3.i <- which(folds.ERE==3,arr.ind=T)
train3.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test3.i,])
train3.logit.ERE.SR1.pooled <- glmnet.cr(x=train3.d,
                                         y=dt.ERE[-test3.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train3.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train3.d.Rda")
save(train3.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train3.logit.ERE.SR1.pooled.Rda")
rm(train3.d)
rm(train3.logit.ERE.SR1.pooled)
rm(test3.i)

#dropping out test4 set, infer model
test4.i <- which(folds.ERE==4,arr.ind=T)
train4.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test4.i,])
train4.logit.ERE.SR1.pooled <- glmnet.cr(x=train4.d,
                                         y=dt.ERE[-test4.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train4.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train4.d.Rda")
save(train4.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train4.logit.ERE.SR1.pooled.Rda")
rm(train4.d)
rm(train4.logit.ERE.SR1.pooled)
rm(test4.i)

#dropping out test5 set, infer model
test5.i <- which(folds.ERE==5,arr.ind=T)
train5.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test5.i,])
train5.logit.ERE.SR1.pooled <- glmnet.cr(x=train5.d,
                                         y=dt.ERE[-test5.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train5.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train5.d.Rda")
save(train5.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train5.logit.ERE.SR1.pooled.Rda")
rm(train5.d)
rm(train5.logit.ERE.SR1.pooled)
rm(test5.i)

#dropping out test6 set, infer model
test6.i <- which(folds.ERE==6,arr.ind=T)
train6.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test6.i,])
train6.logit.ERE.SR1.pooled <- glmnet.cr(x=train6.d,
                                         y=dt.ERE[-test6.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train6.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train6.d.Rda")
save(train6.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train6.logit.ERE.SR1.pooled.Rda")
rm(train6.d)
rm(train6.logit.ERE.SR1.pooled)
rm(test6.i)

#dropping out test7 set, infer model
test7.i <- which(folds.ERE==7,arr.ind=T)
train7.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test7.i,])
train7.logit.ERE.SR1.pooled <- glmnet.cr(x=train7.d,
                                         y=dt.ERE[-test7.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train7.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train7.d.Rda")
save(train7.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train7.logit.ERE.SR1.pooled.Rda")
rm(train7.d)
rm(train7.logit.ERE.SR1.pooled)
rm(test7.i)

#dropping out test8 set, infer model
test8.i <- which(folds.ERE==8,arr.ind=T)
train8.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test8.i,])
train8.logit.ERE.SR1.pooled <- glmnet.cr(x=train8.d,
                                         y=dt.ERE[-test8.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train8.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train8.d.Rda")
save(train8.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train8.logit.ERE.SR1.pooled.Rda")
rm(train8.d)
rm(train8.logit.ERE.SR1.pooled)
rm(test8.i)

#dropping out test9 set, infer model
test9.i <- which(folds.ERE==9,arr.ind=T)
train9.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test9.i,])
train9.logit.ERE.SR1.pooled <- glmnet.cr(x=train9.d,
                                         y=dt.ERE[-test9.i,ERE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train9.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train9.d.Rda")
save(train9.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train9.logit.ERE.SR1.pooled.Rda")
rm(train9.d)
rm(train9.logit.ERE.SR1.pooled)
rm(test9.i)

#dropping out test10 set, infer model
test10.i <- which(folds.ERE==10,arr.ind=T)
train10.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[-test10.i,])
train10.logit.ERE.SR1.pooled <- glmnet.cr(x=train10.d,
                                          y=dt.ERE[-test10.i,ERE.pooled.class],
                                          maxit=500,
                                          lambda=lambda)
save(train10.d, file="6_predict-missing_out/cross-validate_SR1/ERE/train10.d.Rda")
save(train10.logit.ERE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/ERE/train10.logit.ERE.SR1.pooled.Rda")
rm(train10.d)
rm(train10.logit.ERE.SR1.pooled)
rm(test10.i)

rm(lambda)
rm(folds.ERE)
rm(dt.ERE)

######################################################################
#SRE pooled, SR1:
#get path of lambdas used in full data infSREnce
load(file="./6_predict-missing_out/glmnet.cr.SRE.SR1.pooled.Rda")
lambda <- logit.SRE.SR1.pooled$lambda
rm(logit.SRE.SR1.pooled)

#break data into ten groups of ~equal numbers
#randomly order rows in dt
set.seed(1004)
dt.SRE <- dt.SR1.coding[!is.na(SRE.pooled.class)]
dt.SRE <- dt.SRE[sample(nrow(dt.SRE)),]
folds.SRE <- cut(seq(1,nrow(dt.SRE)),breaks=10,labels=F)
save(dt.SRE, file="./6_predict-missing_out/cross-validate_SR1/SRE/dt.SRE.Rda")
save(folds.SRE, file="./6_predict-missing_out/cross-validate_SR1/SRE/folds.SRE.Rda")

#dropping out test1 set, infer model
test1.i <- which(folds.SRE==1,arr.ind=T)
train1.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test1.i,])
train1.logit.SRE.SR1.pooled <- glmnet.cr(x=train1.d,
                                         y=dt.SRE[-test1.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train1.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train1.d.Rda")
save(train1.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train1.logit.SRE.SR1.pooled.Rda")
rm(train1.d)
rm(train1.logit.SRE.SR1.pooled)
rm(test1.i)

#dropping out test2 set, infer model
test2.i <- which(folds.SRE==2,arr.ind=T)
train2.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test2.i,])
train2.logit.SRE.SR1.pooled <- glmnet.cr(x=train2.d,
                                         y=dt.SRE[-test2.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train2.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train2.d.Rda")
save(train2.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train2.logit.SRE.SR1.pooled.Rda")
rm(train2.d)
rm(train2.logit.SRE.SR1.pooled)
rm(test2.i)

#dropping out test3 set, infer model
test3.i <- which(folds.SRE==3,arr.ind=T)
train3.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test3.i,])
train3.logit.SRE.SR1.pooled <- glmnet.cr(x=train3.d,
                                         y=dt.SRE[-test3.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train3.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train3.d.Rda")
save(train3.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train3.logit.SRE.SR1.pooled.Rda")
rm(train3.d)
rm(train3.logit.SRE.SR1.pooled)
rm(test3.i)

#dropping out test4 set, infer model
test4.i <- which(folds.SRE==4,arr.ind=T)
train4.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test4.i,])
train4.logit.SRE.SR1.pooled <- glmnet.cr(x=train4.d,
                                         y=dt.SRE[-test4.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train4.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train4.d.Rda")
save(train4.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train4.logit.SRE.SR1.pooled.Rda")
rm(train4.d)
rm(train4.logit.SRE.SR1.pooled)
rm(test4.i)

#dropping out test5 set, infer model
test5.i <- which(folds.SRE==5,arr.ind=T)
train5.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test5.i,])
train5.logit.SRE.SR1.pooled <- glmnet.cr(x=train5.d,
                                         y=dt.SRE[-test5.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train5.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train5.d.Rda")
save(train5.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train5.logit.SRE.SR1.pooled.Rda")
rm(train5.d)
rm(train5.logit.SRE.SR1.pooled)
rm(test5.i)

#dropping out test6 set, infer model
test6.i <- which(folds.SRE==6,arr.ind=T)
train6.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test6.i,])
train6.logit.SRE.SR1.pooled <- glmnet.cr(x=train6.d,
                                         y=dt.SRE[-test6.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train6.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train6.d.Rda")
save(train6.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train6.logit.SRE.SR1.pooled.Rda")
rm(train6.d)
rm(train6.logit.SRE.SR1.pooled)
rm(test6.i)

#dropping out test7 set, infer model
test7.i <- which(folds.SRE==7,arr.ind=T)
train7.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test7.i,])
train7.logit.SRE.SR1.pooled <- glmnet.cr(x=train7.d,
                                         y=dt.SRE[-test7.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train7.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train7.d.Rda")
save(train7.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train7.logit.SRE.SR1.pooled.Rda")
rm(train7.d)
rm(train7.logit.SRE.SR1.pooled)
rm(test7.i)

#dropping out test8 set, infer model
test8.i <- which(folds.SRE==8,arr.ind=T)
train8.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test8.i,])
train8.logit.SRE.SR1.pooled <- glmnet.cr(x=train8.d,
                                         y=dt.SRE[-test8.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train8.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train8.d.Rda")
save(train8.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train8.logit.SRE.SR1.pooled.Rda")
rm(train8.d)
rm(train8.logit.SRE.SR1.pooled)
rm(test8.i)

#dropping out test9 set, infer model
test9.i <- which(folds.SRE==9,arr.ind=T)
train9.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test9.i,])
train9.logit.SRE.SR1.pooled <- glmnet.cr(x=train9.d,
                                         y=dt.SRE[-test9.i,SRE.pooled.class],
                                         maxit=500,
                                         lambda=lambda)
save(train9.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train9.d.Rda")
save(train9.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train9.logit.SRE.SR1.pooled.Rda")
rm(train9.d)
rm(train9.logit.SRE.SR1.pooled)
rm(test9.i)

#dropping out test10 set, infer model
test10.i <- which(folds.SRE==10,arr.ind=T)
train10.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[-test10.i,])
train10.logit.SRE.SR1.pooled <- glmnet.cr(x=train10.d,
                                          y=dt.SRE[-test10.i,SRE.pooled.class],
                                          maxit=500,
                                          lambda=lambda)
save(train10.d, file="6_predict-missing_out/cross-validate_SR1/SRE/train10.d.Rda")
save(train10.logit.SRE.SR1.pooled, file="6_predict-missing_out/cross-validate_SR1/SRE/train10.logit.SRE.SR1.pooled.Rda")
rm(train10.d)
rm(train10.logit.SRE.SR1.pooled)
rm(test10.i)

rm(lambda)
rm(folds.SRE)
rm(dt.SRE)

##############################################################################################################
#look at CV stats as a function of lambda: SR1 data
#ERE
load("./6_predict-missing_out/cross-validate_SR1/ERE/dt.ERE.Rda");load("./6_predict-missing_out/cross-validate_SR1/ERE/folds.ERE.Rda")

load(file="./6_predict-missing_out/glmnet.cr.ERE.SR1.pooled.Rda")
lambda <- logit.ERE.SR1.pooled$lambda
rm(logit.ERE.SR1.pooled)

test1.i <- which(folds.ERE==1,arr.ind=T)
test1.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test1.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train1.logit.ERE.SR1.pooled.Rda")
test1.p <- predict(train1.logit.ERE.SR1.pooled,newx=test1.d)
test1.class <- dt.ERE[test1.i,ERE.pooled.class]
rm(test1.d);rm(train1.logit.ERE.SR1.pooled)

test2.i <- which(folds.ERE==2,arr.ind=T)
test2.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test2.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train2.logit.ERE.SR1.pooled.Rda")
test2.p <- predict(train2.logit.ERE.SR1.pooled,newx=test2.d)
test2.class <- dt.ERE[test2.i,ERE.pooled.class]
rm(test2.d);rm(train2.logit.ERE.SR1.pooled)

test3.i <- which(folds.ERE==3,arr.ind=T)
test3.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test3.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train3.logit.ERE.SR1.pooled.Rda")
test3.p <- predict(train3.logit.ERE.SR1.pooled,newx=test3.d)
test3.class <- dt.ERE[test3.i,ERE.pooled.class]
rm(test3.d);rm(train3.logit.ERE.SR1.pooled)

test4.i <- which(folds.ERE==4,arr.ind=T)
test4.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test4.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train4.logit.ERE.SR1.pooled.Rda")
test4.p <- predict(train4.logit.ERE.SR1.pooled,newx=test4.d)
test4.class <- dt.ERE[test4.i,ERE.pooled.class]
rm(test4.d);rm(train4.logit.ERE.SR1.pooled)

test5.i <- which(folds.ERE==5,arr.ind=T)
test5.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test5.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train5.logit.ERE.SR1.pooled.Rda")
test5.p <- predict(train5.logit.ERE.SR1.pooled,newx=test5.d)
test5.class <- dt.ERE[test5.i,ERE.pooled.class]
rm(test5.d);rm(train5.logit.ERE.SR1.pooled)

test6.i <- which(folds.ERE==6,arr.ind=T)
test6.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test6.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train6.logit.ERE.SR1.pooled.Rda")
test6.p <- predict(train6.logit.ERE.SR1.pooled,newx=test6.d)
test6.class <- dt.ERE[test6.i,ERE.pooled.class]
rm(test6.d);rm(train6.logit.ERE.SR1.pooled)

test7.i <- which(folds.ERE==7,arr.ind=T)
test7.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test7.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train7.logit.ERE.SR1.pooled.Rda")
test7.p <- predict(train7.logit.ERE.SR1.pooled,newx=test7.d)
test7.class <- dt.ERE[test7.i,ERE.pooled.class]
rm(test7.d);rm(train7.logit.ERE.SR1.pooled)

test8.i <- which(folds.ERE==8,arr.ind=T)
test8.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test8.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train8.logit.ERE.SR1.pooled.Rda")
test8.p <- predict(train8.logit.ERE.SR1.pooled,newx=test8.d)
test8.class <- dt.ERE[test8.i,ERE.pooled.class]
rm(test8.d);rm(train8.logit.ERE.SR1.pooled)

test9.i <- which(folds.ERE==9,arr.ind=T)
test9.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test9.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train9.logit.ERE.SR1.pooled.Rda")
test9.p <- predict(train9.logit.ERE.SR1.pooled,newx=test9.d)
test9.class <- dt.ERE[test9.i,ERE.pooled.class]
rm(test9.d);rm(train9.logit.ERE.SR1.pooled)

test10.i <- which(folds.ERE==10,arr.ind=T)
test10.d <- model.matrix(ERE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.ERE[test10.i,])
load("./6_predict-missing_out/cross-validate_SR1/ERE/train10.logit.ERE.SR1.pooled.Rda")
test10.p <- predict(train10.logit.ERE.SR1.pooled,newx=test10.d)
test10.class <- dt.ERE[test10.i,ERE.pooled.class]
rm(test10.d);rm(train10.logit.ERE.SR1.pooled)

save(test1.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test1.p.Rda");save(test2.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test2.p.Rda");save(test3.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test3.p.Rda");save(test4.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test4.p.Rda");save(test5.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test5.p.Rda");save(test6.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test6.p.Rda");save(test7.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test7.p.Rda");save(test8.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test8.p.Rda");save(test9.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test9.p.Rda");save(test10.p, file="./6_predict-missing_out/cross-validate_SR1/ERE/test10.p.Rda")
save(test1.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test1.class.Rda");save(test2.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test2.class.Rda");save(test3.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test3.class.Rda");save(test4.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test4.class.Rda");save(test5.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test5.class.Rda");save(test6.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test6.class.Rda");save(test7.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test7.class.Rda");save(test8.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test8.class.Rda");save(test9.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test9.class.Rda");save(test10.class, file="./6_predict-missing_out/cross-validate_SR1/ERE/test10.class.Rda")
save(test1.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test1.i.Rda");save(test2.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test2.i.Rda");save(test3.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test3.i.Rda");save(test4.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test4.i.Rda");save(test5.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test5.i.Rda");save(test6.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test6.i.Rda");save(test7.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test7.i.Rda");save(test8.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test8.i.Rda");save(test9.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test9.i.Rda");save(test10.i, file="./6_predict-missing_out/cross-validate_SR1/ERE/test10.i.Rda")
save(lambda, file="./6_predict-missing_out/cross-validate_SR1/ERE/lambda.Rda")

#record TPR (proportion of true positives called positive in prediction) and PPV (proportion of predicted positives that are true positives)
test1.TPR.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.TPR.strong <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.strong <- vector(mode="numeric",length=ncol(test1.p$class))
for(i in 1:ncol(test1.p$class)){
  model.pred <- test1.p$class[,i]
  test1.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class%in% c("weak","strong"))/sum(test1.class %in% c("weak","strong"));test1.TPR.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(test1.class=="strong")
  test1.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test1.PPV.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(model.pred=="strong")
}

test2.TPR.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.TPR.strong <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.strong <- vector(mode="numeric",length=ncol(test2.p$class))
for(i in 1:ncol(test2.p$class)){
  model.pred <- test2.p$class[,i]
  test2.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class%in% c("weak","strong"))/sum(test2.class %in% c("weak","strong"));test2.TPR.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(test2.class=="strong")
  test2.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test2.PPV.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(model.pred=="strong")
}

test3.TPR.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.TPR.strong <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.strong <- vector(mode="numeric",length=ncol(test3.p$class))
for(i in 1:ncol(test3.p$class)){
  model.pred <- test3.p$class[,i]
  test3.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class%in% c("weak","strong"))/sum(test3.class %in% c("weak","strong"));test3.TPR.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(test3.class=="strong")
  test3.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test3.PPV.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(model.pred=="strong")
}

test4.TPR.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.TPR.strong <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.strong <- vector(mode="numeric",length=ncol(test4.p$class))
for(i in 1:ncol(test4.p$class)){
  model.pred <- test4.p$class[,i]
  test4.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class%in% c("weak","strong"))/sum(test4.class %in% c("weak","strong"));test4.TPR.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(test4.class=="strong")
  test4.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test4.PPV.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(model.pred=="strong")
}

test5.TPR.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.TPR.strong <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.strong <- vector(mode="numeric",length=ncol(test5.p$class))
for(i in 1:ncol(test5.p$class)){
  model.pred <- test5.p$class[,i]
  test5.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class%in% c("weak","strong"))/sum(test5.class %in% c("weak","strong"));test5.TPR.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(test5.class=="strong")
  test5.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test5.PPV.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(model.pred=="strong")
}

test6.TPR.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.TPR.strong <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.strong <- vector(mode="numeric",length=ncol(test6.p$class))
for(i in 1:ncol(test6.p$class)){
  model.pred <- test6.p$class[,i]
  test6.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class%in% c("weak","strong"))/sum(test6.class %in% c("weak","strong"));test6.TPR.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(test6.class=="strong")
  test6.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test6.PPV.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(model.pred=="strong")
}

test7.TPR.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.TPR.strong <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.strong <- vector(mode="numeric",length=ncol(test7.p$class))
for(i in 1:ncol(test7.p$class)){
  model.pred <- test7.p$class[,i]
  test7.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class%in% c("weak","strong"))/sum(test7.class %in% c("weak","strong"));test7.TPR.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(test7.class=="strong")
  test7.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test7.PPV.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(model.pred=="strong")
}

test8.TPR.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.TPR.strong <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.strong <- vector(mode="numeric",length=ncol(test8.p$class))
for(i in 1:ncol(test8.p$class)){
  model.pred <- test8.p$class[,i]
  test8.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class%in% c("weak","strong"))/sum(test8.class %in% c("weak","strong"));test8.TPR.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(test8.class=="strong")
  test8.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test8.PPV.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(model.pred=="strong")
}

test9.TPR.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.TPR.strong <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.strong <- vector(mode="numeric",length=ncol(test9.p$class))
for(i in 1:ncol(test9.p$class)){
  model.pred <- test9.p$class[,i]
  test9.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class%in% c("weak","strong"))/sum(test9.class %in% c("weak","strong"));test9.TPR.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(test9.class=="strong")
  test9.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test9.PPV.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(model.pred=="strong")
}

test10.TPR.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.TPR.strong <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.strong <- vector(mode="numeric",length=ncol(test10.p$class))
for(i in 1:ncol(test10.p$class)){
  model.pred <- test10.p$class[,i]
  test10.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class%in% c("weak","strong"))/sum(test10.class %in% c("weak","strong"));test10.TPR.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(test10.class=="strong")
  test10.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test10.PPV.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(model.pred=="strong")
}

TPR.strong <- vector(mode="numeric",length=length(test1.TPR.strong));TPR.pos <- vector(mode="numeric",length=length(test1.TPR.pos));PPV.strong <- vector(mode="numeric",length=length(test1.PPV.strong));PPV.pos <- vector(mode="numeric",length=length(test1.PPV.pos))
for(i in 1:length(test1.TPR.strong)){
  TPR.strong[i] <- mean(c(test1.TPR.strong[i],test2.TPR.strong[i],test3.TPR.strong[i],test4.TPR.strong[i],test5.TPR.strong[i],test6.TPR.strong[i],test7.TPR.strong[i],test8.TPR.strong[i],test9.TPR.strong[i],test10.TPR.strong[i]),na.rm=T)
  PPV.strong[i] <- mean(c(test1.PPV.strong[i],test2.PPV.strong[i],test3.PPV.strong[i],test4.PPV.strong[i],test5.PPV.strong[i],test6.PPV.strong[i],test7.PPV.strong[i],test8.PPV.strong[i],test9.PPV.strong[i],test10.PPV.strong[i]),na.rm=T)
  TPR.pos[i] <- mean(c(test1.TPR.pos[i],test2.TPR.pos[i],test3.TPR.pos[i],test4.TPR.pos[i],test5.TPR.pos[i],test6.TPR.pos[i],test7.TPR.pos[i],test8.TPR.pos[i],test9.TPR.pos[i],test10.TPR.pos[i]),na.rm=T)
  PPV.pos[i] <- mean(c(test1.PPV.pos[i],test2.PPV.pos[i],test3.PPV.pos[i],test4.PPV.pos[i],test5.PPV.pos[i],test6.PPV.pos[i],test7.PPV.pos[i],test8.PPV.pos[i],test9.PPV.pos[i],test10.PPV.pos[i]),na.rm=T)
}

pdf(file="./6_predict-missing_out/cross-validate_SR1/ERE/TPR-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.TPR.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, strong",main="ERE, SR1")
points(lambda,test2.TPR.strong,type="l",lty=3)
points(lambda,test3.TPR.strong,type="l",lty=3)
points(lambda,test4.TPR.strong,type="l",lty=3)
points(lambda,test5.TPR.strong,type="l",lty=3)
points(lambda,test6.TPR.strong,type="l",lty=3)
points(lambda,test7.TPR.strong,type="l",lty=3)
points(lambda,test8.TPR.strong,type="l",lty=3)
points(lambda,test9.TPR.strong,type="l",lty=3)
points(lambda,test10.TPR.strong,type="l",lty=3)
points(lambda,TPR.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_SR1/ERE/TPR-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.TPR.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, positive",main="ERE, SR1")
points(lambda,test2.TPR.pos,type="l",lty=3)
points(lambda,test3.TPR.pos,type="l",lty=3)
points(lambda,test4.TPR.pos,type="l",lty=3)
points(lambda,test5.TPR.pos,type="l",lty=3)
points(lambda,test6.TPR.pos,type="l",lty=3)
points(lambda,test7.TPR.pos,type="l",lty=3)
points(lambda,test8.TPR.pos,type="l",lty=3)
points(lambda,test9.TPR.pos,type="l",lty=3)
points(lambda,test10.TPR.pos,type="l",lty=3)
points(lambda,TPR.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_SR1/ERE/PPV-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.PPV.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, strong",main="ERE, SR1")
points(lambda,test2.PPV.strong,type="l",lty=3)
points(lambda,test3.PPV.strong,type="l",lty=3)
points(lambda,test4.PPV.strong,type="l",lty=3)
points(lambda,test5.PPV.strong,type="l",lty=3)
points(lambda,test6.PPV.strong,type="l",lty=3)
points(lambda,test7.PPV.strong,type="l",lty=3)
points(lambda,test8.PPV.strong,type="l",lty=3)
points(lambda,test9.PPV.strong,type="l",lty=3)
points(lambda,test10.PPV.strong,type="l",lty=3)
points(lambda,PPV.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_SR1/ERE/PPV-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda,test1.PPV.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, positive",main="ERE, SR1")
points(lambda,test2.PPV.pos,type="l",lty=3)
points(lambda,test3.PPV.pos,type="l",lty=3)
points(lambda,test4.PPV.pos,type="l",lty=3)
points(lambda,test5.PPV.pos,type="l",lty=3)
points(lambda,test6.PPV.pos,type="l",lty=3)
points(lambda,test7.PPV.pos,type="l",lty=3)
points(lambda,test8.PPV.pos,type="l",lty=3)
points(lambda,test9.PPV.pos,type="l",lty=3)
points(lambda,test10.PPV.pos,type="l",lty=3)
points(lambda,PPV.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

step.ERE.SR1 <- which(abs(lambda-(1e-5))==min(abs(lambda-(1e-5))))
save(step.ERE.SR1,file="./6_predict-missing_out/step.ERE.SR1.Rda")

lambda[step.ERE.SR1] 

#bring together all cross-validation predictions: in dt.ERE, add column giving prediction.CV -- use this to compute PPV and TPR, since above these shouldn't be averaged completely, should they? different subsets have different #'s of positives etc.
dt.ERE[test1.i,prediction.CV:=test1.p$class[,step.ERE.SR1]]
dt.ERE[test2.i,prediction.CV:=test2.p$class[,step.ERE.SR1]]
dt.ERE[test3.i,prediction.CV:=test3.p$class[,step.ERE.SR1]]
dt.ERE[test4.i,prediction.CV:=test4.p$class[,step.ERE.SR1]]
dt.ERE[test5.i,prediction.CV:=test5.p$class[,step.ERE.SR1]]
dt.ERE[test6.i,prediction.CV:=test6.p$class[,step.ERE.SR1]]
dt.ERE[test7.i,prediction.CV:=test7.p$class[,step.ERE.SR1]]
dt.ERE[test8.i,prediction.CV:=test8.p$class[,step.ERE.SR1]]
dt.ERE[test9.i,prediction.CV:=test9.p$class[,step.ERE.SR1]]
dt.ERE[test10.i,prediction.CV:=test10.p$class[,step.ERE.SR1]]

write.table(data.frame(table(dt.ERE$ERE.pooled.class,dt.ERE$prediction.CV,useNA="always",dnn=c("expt","pred"))),file="6_predict-missing_out/cross-validate_SR1/ERE/CV-summary.txt")

#SRE
load("./6_predict-missing_out/cross-validate_SR1/SRE/dt.SRE.Rda");load("./6_predict-missing_out/cross-validate_SR1/SRE/folds.SRE.Rda")

load(file="./6_predict-missing_out/glmnet.cr.SRE.SR1.pooled.Rda")
lambda <- logit.SRE.SR1.pooled$lambda
rm(logit.SRE.SR1.pooled)

test1.i <- which(folds.SRE==1,arr.ind=T)
test1.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test1.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train1.logit.SRE.SR1.pooled.Rda")
test1.p <- predict(train1.logit.SRE.SR1.pooled,newx=test1.d)
test1.class <- dt.SRE[test1.i,SRE.pooled.class]
rm(test1.d);rm(train1.logit.SRE.SR1.pooled)

test2.i <- which(folds.SRE==2,arr.ind=T)
test2.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test2.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train2.logit.SRE.SR1.pooled.Rda")
test2.p <- predict(train2.logit.SRE.SR1.pooled,newx=test2.d)
test2.class <- dt.SRE[test2.i,SRE.pooled.class]
rm(test2.d);rm(train2.logit.SRE.SR1.pooled)

test3.i <- which(folds.SRE==3,arr.ind=T)
test3.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test3.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train3.logit.SRE.SR1.pooled.Rda")
test3.p <- predict(train3.logit.SRE.SR1.pooled,newx=test3.d)
test3.class <- dt.SRE[test3.i,SRE.pooled.class]
rm(test3.d);rm(train3.logit.SRE.SR1.pooled)

test4.i <- which(folds.SRE==4,arr.ind=T)
test4.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test4.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train4.logit.SRE.SR1.pooled.Rda")
test4.p <- predict(train4.logit.SRE.SR1.pooled,newx=test4.d)
test4.class <- dt.SRE[test4.i,SRE.pooled.class]
rm(test4.d);rm(train4.logit.SRE.SR1.pooled)

test5.i <- which(folds.SRE==5,arr.ind=T)
test5.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test5.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train5.logit.SRE.SR1.pooled.Rda")
test5.p <- predict(train5.logit.SRE.SR1.pooled,newx=test5.d)
test5.class <- dt.SRE[test5.i,SRE.pooled.class]
rm(test5.d);rm(train5.logit.SRE.SR1.pooled)

test6.i <- which(folds.SRE==6,arr.ind=T)
test6.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test6.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train6.logit.SRE.SR1.pooled.Rda")
test6.p <- predict(train6.logit.SRE.SR1.pooled,newx=test6.d)
test6.class <- dt.SRE[test6.i,SRE.pooled.class]
rm(test6.d);rm(train6.logit.SRE.SR1.pooled)

test7.i <- which(folds.SRE==7,arr.ind=T)
test7.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test7.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train7.logit.SRE.SR1.pooled.Rda")
test7.p <- predict(train7.logit.SRE.SR1.pooled,newx=test7.d)
test7.class <- dt.SRE[test7.i,SRE.pooled.class]
rm(test7.d);rm(train7.logit.SRE.SR1.pooled)

test8.i <- which(folds.SRE==8,arr.ind=T)
test8.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test8.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train8.logit.SRE.SR1.pooled.Rda")
test8.p <- predict(train8.logit.SRE.SR1.pooled,newx=test8.d)
test8.class <- dt.SRE[test8.i,SRE.pooled.class]
rm(test8.d);rm(train8.logit.SRE.SR1.pooled)

test9.i <- which(folds.SRE==9,arr.ind=T)
test9.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test9.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train9.logit.SRE.SR1.pooled.Rda")
test9.p <- predict(train9.logit.SRE.SR1.pooled,newx=test9.d)
test9.class <- dt.SRE[test9.i,SRE.pooled.class]
rm(test9.d);rm(train9.logit.SRE.SR1.pooled)

test10.i <- which(folds.SRE==10,arr.ind=T)
test10.d <- model.matrix(SRE.pooled.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SRE[test10.i,])
load("./6_predict-missing_out/cross-validate_SR1/SRE/train10.logit.SRE.SR1.pooled.Rda")
test10.p <- predict(train10.logit.SRE.SR1.pooled,newx=test10.d)
test10.class <- dt.SRE[test10.i,SRE.pooled.class]
rm(test10.d);rm(train10.logit.SRE.SR1.pooled)

save(test1.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test1.p.Rda");save(test2.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test2.p.Rda");save(test3.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test3.p.Rda");save(test4.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test4.p.Rda");save(test5.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test5.p.Rda");save(test6.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test6.p.Rda");save(test7.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test7.p.Rda");save(test8.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test8.p.Rda");save(test9.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test9.p.Rda");save(test10.p, file="./6_predict-missing_out/cross-validate_SR1/SRE/test10.p.Rda")
save(test1.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test1.class.Rda");save(test2.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test2.class.Rda");save(test3.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test3.class.Rda");save(test4.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test4.class.Rda");save(test5.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test5.class.Rda");save(test6.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test6.class.Rda");save(test7.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test7.class.Rda");save(test8.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test8.class.Rda");save(test9.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test9.class.Rda");save(test10.class, file="./6_predict-missing_out/cross-validate_SR1/SRE/test10.class.Rda")
save(test1.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test1.i.Rda");save(test2.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test2.i.Rda");save(test3.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test3.i.Rda");save(test4.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test4.i.Rda");save(test5.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test5.i.Rda");save(test6.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test6.i.Rda");save(test7.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test7.i.Rda");save(test8.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test8.i.Rda");save(test9.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test9.i.Rda");save(test10.i, file="./6_predict-missing_out/cross-validate_SR1/SRE/test10.i.Rda")
save(lambda, file="./6_predict-missing_out/cross-validate_SR1/SRE/lambda.Rda")

#record TPR (proportion of true positives called positive in prediction) and PPV (proportion of predicted positives that are true positives)
test1.TPR.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.TPR.strong <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.pos <- vector(mode="numeric",length=ncol(test1.p$class));test1.PPV.strong <- vector(mode="numeric",length=ncol(test1.p$class))
for(i in 1:ncol(test1.p$class)){
  model.pred <- test1.p$class[,i]
  test1.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class%in% c("weak","strong"))/sum(test1.class %in% c("weak","strong"));test1.TPR.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(test1.class=="strong")
  test1.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test1.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test1.PPV.strong[i] <- sum(model.pred=="strong" & test1.class=="strong")/sum(model.pred=="strong")
}

test2.TPR.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.TPR.strong <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.pos <- vector(mode="numeric",length=ncol(test2.p$class));test2.PPV.strong <- vector(mode="numeric",length=ncol(test2.p$class))
for(i in 1:ncol(test2.p$class)){
  model.pred <- test2.p$class[,i]
  test2.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class%in% c("weak","strong"))/sum(test2.class %in% c("weak","strong"));test2.TPR.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(test2.class=="strong")
  test2.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test2.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test2.PPV.strong[i] <- sum(model.pred=="strong" & test2.class=="strong")/sum(model.pred=="strong")
}

test3.TPR.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.TPR.strong <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.pos <- vector(mode="numeric",length=ncol(test3.p$class));test3.PPV.strong <- vector(mode="numeric",length=ncol(test3.p$class))
for(i in 1:ncol(test3.p$class)){
  model.pred <- test3.p$class[,i]
  test3.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class%in% c("weak","strong"))/sum(test3.class %in% c("weak","strong"));test3.TPR.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(test3.class=="strong")
  test3.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test3.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test3.PPV.strong[i] <- sum(model.pred=="strong" & test3.class=="strong")/sum(model.pred=="strong")
}

test4.TPR.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.TPR.strong <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.pos <- vector(mode="numeric",length=ncol(test4.p$class));test4.PPV.strong <- vector(mode="numeric",length=ncol(test4.p$class))
for(i in 1:ncol(test4.p$class)){
  model.pred <- test4.p$class[,i]
  test4.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class%in% c("weak","strong"))/sum(test4.class %in% c("weak","strong"));test4.TPR.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(test4.class=="strong")
  test4.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test4.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test4.PPV.strong[i] <- sum(model.pred=="strong" & test4.class=="strong")/sum(model.pred=="strong")
}

test5.TPR.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.TPR.strong <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.pos <- vector(mode="numeric",length=ncol(test5.p$class));test5.PPV.strong <- vector(mode="numeric",length=ncol(test5.p$class))
for(i in 1:ncol(test5.p$class)){
  model.pred <- test5.p$class[,i]
  test5.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class%in% c("weak","strong"))/sum(test5.class %in% c("weak","strong"));test5.TPR.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(test5.class=="strong")
  test5.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test5.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test5.PPV.strong[i] <- sum(model.pred=="strong" & test5.class=="strong")/sum(model.pred=="strong")
}

test6.TPR.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.TPR.strong <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.pos <- vector(mode="numeric",length=ncol(test6.p$class));test6.PPV.strong <- vector(mode="numeric",length=ncol(test6.p$class))
for(i in 1:ncol(test6.p$class)){
  model.pred <- test6.p$class[,i]
  test6.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class%in% c("weak","strong"))/sum(test6.class %in% c("weak","strong"));test6.TPR.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(test6.class=="strong")
  test6.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test6.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test6.PPV.strong[i] <- sum(model.pred=="strong" & test6.class=="strong")/sum(model.pred=="strong")
}

test7.TPR.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.TPR.strong <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.pos <- vector(mode="numeric",length=ncol(test7.p$class));test7.PPV.strong <- vector(mode="numeric",length=ncol(test7.p$class))
for(i in 1:ncol(test7.p$class)){
  model.pred <- test7.p$class[,i]
  test7.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class%in% c("weak","strong"))/sum(test7.class %in% c("weak","strong"));test7.TPR.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(test7.class=="strong")
  test7.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test7.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test7.PPV.strong[i] <- sum(model.pred=="strong" & test7.class=="strong")/sum(model.pred=="strong")
}

test8.TPR.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.TPR.strong <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.pos <- vector(mode="numeric",length=ncol(test8.p$class));test8.PPV.strong <- vector(mode="numeric",length=ncol(test8.p$class))
for(i in 1:ncol(test8.p$class)){
  model.pred <- test8.p$class[,i]
  test8.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class%in% c("weak","strong"))/sum(test8.class %in% c("weak","strong"));test8.TPR.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(test8.class=="strong")
  test8.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test8.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test8.PPV.strong[i] <- sum(model.pred=="strong" & test8.class=="strong")/sum(model.pred=="strong")
}

test9.TPR.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.TPR.strong <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.pos <- vector(mode="numeric",length=ncol(test9.p$class));test9.PPV.strong <- vector(mode="numeric",length=ncol(test9.p$class))
for(i in 1:ncol(test9.p$class)){
  model.pred <- test9.p$class[,i]
  test9.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class%in% c("weak","strong"))/sum(test9.class %in% c("weak","strong"));test9.TPR.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(test9.class=="strong")
  test9.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test9.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test9.PPV.strong[i] <- sum(model.pred=="strong" & test9.class=="strong")/sum(model.pred=="strong")
}

test10.TPR.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.TPR.strong <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.pos <- vector(mode="numeric",length=ncol(test10.p$class));test10.PPV.strong <- vector(mode="numeric",length=ncol(test10.p$class))
for(i in 1:ncol(test10.p$class)){
  model.pred <- test10.p$class[,i]
  test10.TPR.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class%in% c("weak","strong"))/sum(test10.class %in% c("weak","strong"));test10.TPR.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(test10.class=="strong")
  test10.PPV.pos[i] <- sum(model.pred %in% c("weak","strong") & test10.class %in% c("weak","strong"))/sum(model.pred %in% c("weak","strong"));test10.PPV.strong[i] <- sum(model.pred=="strong" & test10.class=="strong")/sum(model.pred=="strong")
}

TPR.strong <- vector(mode="numeric",length=length(test1.TPR.strong));TPR.pos <- vector(mode="numeric",length=length(test1.TPR.pos));PPV.strong <- vector(mode="numeric",length=length(test1.PPV.strong));PPV.pos <- vector(mode="numeric",length=length(test1.PPV.pos))
for(i in 1:length(test1.TPR.strong)){
  TPR.strong[i] <- mean(c(test1.TPR.strong[i],test2.TPR.strong[i],test3.TPR.strong[i],test4.TPR.strong[i],test5.TPR.strong[i],test6.TPR.strong[i],test7.TPR.strong[i],test8.TPR.strong[i],test9.TPR.strong[i],test10.TPR.strong[i]),na.rm=T)
  PPV.strong[i] <- mean(c(test1.PPV.strong[i],test2.PPV.strong[i],test3.PPV.strong[i],test4.PPV.strong[i],test5.PPV.strong[i],test6.PPV.strong[i],test7.PPV.strong[i],test8.PPV.strong[i],test9.PPV.strong[i],test10.PPV.strong[i]),na.rm=T)
  TPR.pos[i] <- mean(c(test1.TPR.pos[i],test2.TPR.pos[i],test3.TPR.pos[i],test4.TPR.pos[i],test5.TPR.pos[i],test6.TPR.pos[i],test7.TPR.pos[i],test8.TPR.pos[i],test9.TPR.pos[i],test10.TPR.pos[i]),na.rm=T)
  PPV.pos[i] <- mean(c(test1.PPV.pos[i],test2.PPV.pos[i],test3.PPV.pos[i],test4.PPV.pos[i],test5.PPV.pos[i],test6.PPV.pos[i],test7.PPV.pos[i],test8.PPV.pos[i],test9.PPV.pos[i],test10.PPV.pos[i]),na.rm=T)
}

pdf(file="./6_predict-missing_out/cross-validate_SR1/SRE/TPR-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda[1:length(test1.TPR.strong)],test1.TPR.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, strong",main="SRE, SR1")
points(lambda[1:length(test2.TPR.strong)],test2.TPR.strong,type="l",lty=3)
points(lambda[1:length(test3.TPR.strong)],test3.TPR.strong,type="l",lty=3)
points(lambda[1:length(test4.TPR.strong)],test4.TPR.strong,type="l",lty=3)
points(lambda[1:length(test5.TPR.strong)],test5.TPR.strong,type="l",lty=3)
points(lambda[1:length(test6.TPR.strong)],test6.TPR.strong,type="l",lty=3)
points(lambda[1:length(test7.TPR.strong)],test7.TPR.strong,type="l",lty=3)
points(lambda[1:length(test8.TPR.strong)],test8.TPR.strong,type="l",lty=3)
points(lambda[1:length(test9.TPR.strong)],test9.TPR.strong,type="l",lty=3)
points(lambda[1:length(test10.TPR.strong)],test10.TPR.strong,type="l",lty=3)
points(lambda[1:length(TPR.strong)],TPR.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_SR1/SRE/TPR-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda[1:length(test1.TPR.pos)],test1.TPR.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="TPR, pos",main="SRE, SR1")
points(lambda[1:length(test2.TPR.pos)],test2.TPR.pos,type="l",lty=3)
points(lambda[1:length(test3.TPR.pos)],test3.TPR.pos,type="l",lty=3)
points(lambda[1:length(test4.TPR.pos)],test4.TPR.pos,type="l",lty=3)
points(lambda[1:length(test5.TPR.pos)],test5.TPR.pos,type="l",lty=3)
points(lambda[1:length(test6.TPR.pos)],test6.TPR.pos,type="l",lty=3)
points(lambda[1:length(test7.TPR.pos)],test7.TPR.pos,type="l",lty=3)
points(lambda[1:length(test8.TPR.pos)],test8.TPR.pos,type="l",lty=3)
points(lambda[1:length(test9.TPR.pos)],test9.TPR.pos,type="l",lty=3)
points(lambda[1:length(test10.TPR.pos)],test10.TPR.pos,type="l",lty=3)
points(lambda[1:length(TPR.pos)],TPR.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_SR1/SRE/PPV-strong-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda[1:length(test1.PPV.strong)],test1.PPV.strong,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, strong",main="SRE, SR1")
points(lambda[1:length(test2.PPV.strong)],test2.PPV.strong,type="l",lty=3)
points(lambda[1:length(test3.PPV.strong)],test3.PPV.strong,type="l",lty=3)
points(lambda[1:length(test4.PPV.strong)],test4.PPV.strong,type="l",lty=3)
points(lambda[1:length(test5.PPV.strong)],test5.PPV.strong,type="l",lty=3)
points(lambda[1:length(test6.PPV.strong)],test6.PPV.strong,type="l",lty=3)
points(lambda[1:length(test7.PPV.strong)],test7.PPV.strong,type="l",lty=3)
points(lambda[1:length(test8.PPV.strong)],test8.PPV.strong,type="l",lty=3)
points(lambda[1:length(test9.PPV.strong)],test9.PPV.strong,type="l",lty=3)
points(lambda[1:length(test10.PPV.strong)],test10.PPV.strong,type="l",lty=3)
points(lambda[1:length(PPV.strong)],PPV.strong,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

pdf(file="./6_predict-missing_out/cross-validate_SR1/SRE/PPV-pos-v-lambda.pdf",5,5,useDingbats=F)
plot(lambda[1:length(test1.PPV.pos)],test1.PPV.pos,log="x",type="l",ylim=c(0,1),lty=3,xlab="lambda",ylab="PPV, pos",main="SRE, SR1")
points(lambda[1:length(test2.PPV.pos)],test2.PPV.pos,type="l",lty=3)
points(lambda[1:length(test3.PPV.pos)],test3.PPV.pos,type="l",lty=3)
points(lambda[1:length(test4.PPV.pos)],test4.PPV.pos,type="l",lty=3)
points(lambda[1:length(test5.PPV.pos)],test5.PPV.pos,type="l",lty=3)
points(lambda[1:length(test6.PPV.pos)],test6.PPV.pos,type="l",lty=3)
points(lambda[1:length(test7.PPV.pos)],test7.PPV.pos,type="l",lty=3)
points(lambda[1:length(test8.PPV.pos)],test8.PPV.pos,type="l",lty=3)
points(lambda[1:length(test9.PPV.pos)],test9.PPV.pos,type="l",lty=3)
points(lambda[1:length(test10.PPV.pos)],test10.PPV.pos,type="l",lty=3)
points(lambda[1:length(PPV.pos)],PPV.pos,type="l",lty=1,lwd=2)
abline(v=1e-5)
dev.off()

step.SRE.SR1 <- which(abs(lambda-(1e-5))==min(abs(lambda-(1e-5))))
save(step.SRE.SR1,file="./6_predict-missing_out/step.SRE.SR1.Rda")

lambda[step.SRE.SR1] 

#bring together all cross-validation predictions: in dt.SRE, add column giving prediction.CV -- use this to compute PPV and TPR, since above these shouldn't be averaged completely, should they? diffSREnt subsets have diffSREnt #'s of positives etc.
dt.SRE[test1.i,prediction.CV:=test1.p$class[,step.SRE.SR1]]
dt.SRE[test2.i,prediction.CV:=test2.p$class[,step.SRE.SR1]]
dt.SRE[test3.i,prediction.CV:=test3.p$class[,step.SRE.SR1]]
dt.SRE[test4.i,prediction.CV:=test4.p$class[,step.SRE.SR1]]
dt.SRE[test5.i,prediction.CV:=test5.p$class[,step.SRE.SR1]]
dt.SRE[test6.i,prediction.CV:=test6.p$class[,step.SRE.SR1]]
dt.SRE[test7.i,prediction.CV:=test7.p$class[,step.SRE.SR1]]
dt.SRE[test8.i,prediction.CV:=test8.p$class[,step.SRE.SR1]]
dt.SRE[test9.i,prediction.CV:=test9.p$class[,step.SRE.SR1]]
dt.SRE[test10.i,prediction.CV:=test10.p$class[,step.SRE.SR1]]

write.table(data.frame(table(dt.SRE$SRE.pooled.class,dt.SRE$prediction.CV,useNA="always",dnn=c("expt","pred"))),file="6_predict-missing_out/cross-validate_SR1/SRE/CV-summary.txt")

rm(list=ls())

#################################################################################
#make figure: number of non-zero parameters versus lambda
load(file="./6_predict-missing_out/glmnet.cr.ERE.11P.pooled.Rda")
load(file="./6_predict-missing_out/glmnet.cr.SRE.11P.pooled.Rda")
load(file="./6_predict-missing_out/glmnet.cr.ERE.SR1.pooled.Rda")
load(file="./6_predict-missing_out/glmnet.cr.SRE.SR1.pooled.Rda")

num.coef.ERE.11P <- sapply(1:length(logit.ERE.11P.pooled$lambda), function(x) return(length(nonzero.glmnet.cr(logit.ERE.11P.pooled,s=x)$beta)-2))
lambda.ERE.11P <- logit.ERE.11P.pooled$lambda

num.coef.SRE.11P <- sapply(1:length(logit.SRE.11P.pooled$lambda), function(x) return(length(nonzero.glmnet.cr(logit.SRE.11P.pooled,s=x)$beta)-2))
lambda.SRE.11P <- logit.SRE.11P.pooled$lambda

num.coef.ERE.SR1 <- sapply(1:length(logit.ERE.SR1.pooled$lambda), function(x) return(length(nonzero.glmnet.cr(logit.ERE.SR1.pooled,s=x)$beta)-2))
lambda.ERE.SR1 <- logit.ERE.SR1.pooled$lambda

num.coef.SRE.SR1 <- sapply(1:length(logit.SRE.SR1.pooled$lambda), function(x) return(length(nonzero.glmnet.cr(logit.SRE.SR1.pooled,s=x)$beta)-2))
lambda.SRE.SR1 <- logit.SRE.SR1.pooled$lambda

pdf(file="6_predict-missing_out/number-nonzero-params_v_lambda.pdf",4,4)
plot(lambda.ERE.11P,num.coef.ERE.11P,log="x",type="l",lty=1,col=rgb(100,69,155,maxColorValue=255),lwd=2,xlim=c(min(c(lambda.ERE.11P,lambda.SRE.11P,lambda.ERE.SR1,lambda.SRE.SR1)),max(c(lambda.ERE.11P,lambda.SRE.11P,lambda.ERE.SR1,lambda.SRE.SR1))),ylim=c(min(c(num.coef.ERE.11P,num.coef.SRE.11P,num.coef.ERE.SR1,num.coef.SRE.SR1)),max(c(num.coef.ERE.11P,num.coef.SRE.11P,num.coef.ERE.SR1,num.coef.SRE.SR1))),xlab="lambda",ylab="number of non-zero parameters")
points(lambda.SRE.11P,num.coef.SRE.11P,type="l",lty=1,col=rgb(5,172,72,maxColorValue=255),lwd=2)
points(lambda.ERE.SR1,num.coef.ERE.SR1,type="l",lty=2,col=rgb(100,69,155,maxColorValue=255),lwd=2)
points(lambda.SRE.SR1,num.coef.SRE.SR1,type="l",lty=2,col=rgb(5,172,72,maxColorValue=255),lwd=2)
abline(v=1e-5,lty=2)
legend("topright",legend=c("ERE 11P","SRE 11P","ERE SR1","SRE SR1"),lwd=2,lty=c(1,1,2,2),col=c(rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255),rgb(100,69,155,maxColorValue=255),rgb(5,172,72,maxColorValue=255)),bty="n")
dev.off()

#################################################################################
#predict class of missing genos
#chose lambda = 1e-5 from cross-validation scheme
#11P data:
dt.11P <- read.table(file="5_assess-replicates_out/dt_output_11P_censor15.csv",header=TRUE,sep=",");dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)
dt.11P[,c("AA1","AA2","AA3","AA4") := list(strsplit(as.character(AAseq),split="")[[1]][1],strsplit(as.character(AAseq),split="")[[1]][2],strsplit(as.character(AAseq),split="")[[1]][3],strsplit(as.character(AAseq),split="")[[1]][4]),by=AAseq]
dt.11P[,AA1 := as.factor(AA1)];dt.11P[,AA2 := as.factor(AA2)];dt.11P[,AA3 := as.factor(AA3)];dt.11P[,AA4 := as.factor(AA4)]

#set order of classification factors for ERE, SRE classes
dt.11P[,ERE.rep1.class:=factor(ERE.rep1.class,levels=c("null","weak","strong"))]
dt.11P[,ERE.rep2.class:=factor(ERE.rep2.class,levels=c("null","weak","strong"))]
dt.11P[,ERE.pooled.class:=factor(ERE.pooled.class,levels=c("null","weak","strong"))]
dt.11P[,SRE.rep1.class:=factor(SRE.rep1.class,levels=c("null","weak","strong"))]
dt.11P[,SRE.rep2.class:=factor(SRE.rep2.class,levels=c("null","weak","strong"))]
dt.11P[,SRE.pooled.class:=factor(SRE.pooled.class,levels=c("null","weak","strong"))]

dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]
rm(dt.11P)

#SR1
dt.SR1 <- read.table(file="5_assess-replicates_out/dt_output_SR1_censor15.csv",header=TRUE,sep=",");dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)
dt.SR1[,c("AA1","AA2","AA3","AA4") := list(strsplit(as.character(AAseq),split="")[[1]][1],strsplit(as.character(AAseq),split="")[[1]][2],strsplit(as.character(AAseq),split="")[[1]][3],strsplit(as.character(AAseq),split="")[[1]][4]),by=AAseq]
dt.SR1[,AA1 := as.factor(AA1)];dt.SR1[,AA2 := as.factor(AA2)];dt.SR1[,AA3 := as.factor(AA3)];dt.SR1[,AA4 := as.factor(AA4)]

#set order of classification factors for ERE, SRE classes
dt.SR1[,ERE.rep1.class:=factor(ERE.rep1.class,levels=c("null","weak","strong"))]
dt.SR1[,ERE.rep2.class:=factor(ERE.rep2.class,levels=c("null","weak","strong"))]
dt.SR1[,ERE.pooled.class:=factor(ERE.pooled.class,levels=c("null","weak","strong"))]
dt.SR1[,SRE.rep1.class:=factor(SRE.rep1.class,levels=c("null","weak","strong"))]
dt.SR1[,SRE.rep2.class:=factor(SRE.rep2.class,levels=c("null","weak","strong"))]
dt.SR1[,SRE.pooled.class:=factor(SRE.pooled.class,levels=c("null","weak","strong"))]

dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]
rm(dt.SR1)

dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)

load(file="./6_predict-missing_out/glmnet.cr.ERE.11P.pooled.Rda")
load(file="./6_predict-missing_out/step.ERE.11P.Rda")
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m, s=step.ERE.11P)
save(ERE.11P.predictions, file="./6_predict-missing_out/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.prediction <- ERE.11P.predictions$class

load(file="./6_predict-missing_out/glmnet.cr.SRE.11P.pooled.Rda")
load(file="./6_predict-missing_out/step.SRE.11P.Rda")
SRE.11P.predictions <- fitted(logit.SRE.11P.pooled, newx=dt.m, s=step.SRE.11P)
save(SRE.11P.predictions, file="./6_predict-missing_out/SRE.11P.predictions.Rda")
dt.11P.coding$SRE.prediction <- SRE.11P.predictions$class

load(file="./6_predict-missing_out/glmnet.cr.ERE.SR1.pooled.Rda")
load(file="./6_predict-missing_out/step.ERE.SR1.Rda")
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m, s=step.ERE.SR1)
save(ERE.SR1.predictions, file="./6_predict-missing_out/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.prediction <- ERE.SR1.predictions$class

load(file="./6_predict-missing_out/glmnet.cr.SRE.SR1.pooled.Rda")
load(file="./6_predict-missing_out/step.SRE.SR1.Rda")
SRE.SR1.predictions <- fitted(logit.SRE.SR1.pooled, newx=dt.m, s=step.SRE.SR1)
save(SRE.SR1.predictions, file="./6_predict-missing_out/SRE.SR1.predictions.Rda")
dt.SR1.coding$SRE.prediction <- SRE.SR1.predictions$class

#define new columns in dt.11P.coding and dt.SR1.coding: ERE.full.class and SRE.full.class:
#these columns are the exptl classification if determined (>= 15 cfu), prediction if otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]
dt.SR1.coding[,ERE.full.class := ERE.pooled.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]

#define columns for "specificity", given ERE.full.class and SRE.full.class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]

write.table(dt.11P.coding,file="6_predict-missing_out/dt_output_11P.csv",col.names=T,row.names=F,quote=F,sep=",")

write.table(dt.SR1.coding,file="6_predict-missing_out/dt_output_SR1.csv",col.names=T,row.names=F,quote=F,sep=",")