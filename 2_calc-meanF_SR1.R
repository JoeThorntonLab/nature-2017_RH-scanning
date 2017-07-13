#16 Jan 2017
#TNS
#read in demultiplexed sequence reads, eliminate possible PCR artefacts (?), calculate ML meanF for each library genotype in AncSR1 background for ERE, SRE binding
#library codes: l10 = ERE, rep1 ; l12 = ERE, rep 2; l9 = SRE, rep 1; l11 = SRE, rep 2

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)
library(fitdistrplus)

#read input file, which is a table of all RH variants and the number of sequence reads mapped to that RH for each library/barcode combo
dt <- read.csv(file="2_calc-meanF_SR1_in/data.in.SR1.csv",header=T)

# #check for any variants with abnormal variation in # reads across duplicate bc PCRs == evidence for PCR bias? (only visible in b1/b2 which have multiple bcs/bin)
# #simply visually ID outliers, remove reads for that variant across all bins for that replicate
# #l9 (SRE rep1) bin 1, dt
# plot(dt$l9_b1_bc01+.1,dt$l9_b1_bc06+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# plot(dt$l9_b1_bc01+.1,dt$l9_b1_bc10+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# plot(dt$l9_b1_bc06+.1,dt$l9_b1_bc10+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers

# #l9 bin2, dt
# plot(dt$l9_b2_bc03+.1,dt$l9_b2_bc09+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# plot(dt$l9_b2_bc03+.1,dt$l9_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# plot(dt$l9_b2_bc09+.1,dt$l9_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers

# #l9 bin2, dt
# plot(dt$l9_b2_bc02+.1,dt$l9_b2_bc07+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# plot(dt$l9_b2_bc02+.1,dt$l9_b2_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# plot(dt$l9_b2_bc07+.1,dt$l9_b2_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers

# #l10 (ERE rep1) bin1, dt
# plot(dt$l10_b1_bc05+.1,dt$l10_b1_bc12+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 5 & dt$l10_b1_bc12 > 10),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 5 & dt$l10_b1_bc12 > 10),"l10_b1_bc12"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 8 & dt$l10_b1_bc12 > 4),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 8 & dt$l10_b1_bc12 > 4),"l10_b1_bc12"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 18 & dt$l10_b1_bc12 > 1),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 18 & dt$l10_b1_bc12 > 1),"l10_b1_bc12"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 50 & dt$l10_b1_bc12 == 1),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 50 & dt$l10_b1_bc12 == 1),"l10_b1_bc12"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05 > 50 & dt$l10_b1_bc12 == 0),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05 > 50 & dt$l10_b1_bc12 == 0),"l10_b1_bc12"]+.1,pch=20,col="coral")
AAseqs.l10 <- unique(as.character(dt[c(which(dt$l10_b1_bc05 > 50 & dt$l10_b1_bc12 == 0),which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 5 & dt$l10_b1_bc12 > 10),which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 8 & dt$l10_b1_bc12 > 4),which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 18 & dt$l10_b1_bc12 > 1),which(dt$l10_b1_bc05/dt$l10_b1_bc12 > 50 & dt$l10_b1_bc12 == 1)),1]))

# plot(dt$l10_b1_bc05+.1,dt$l10_b1_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 4 & dt$l10_b1_bc14 > 20),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 4 & dt$l10_b1_bc14 > 20),"l10_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 6 & dt$l10_b1_bc14 > 10),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 6 & dt$l10_b1_bc14 > 10),"l10_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 10 & dt$l10_b1_bc14 > 5),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 10 & dt$l10_b1_bc14 > 5),"l10_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 12 & dt$l10_b1_bc14 > 2),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 12 & dt$l10_b1_bc14 > 2),"l10_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 20 & dt$l10_b1_bc14 == 2),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 20 & dt$l10_b1_bc14 == 2),"l10_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 30 & dt$l10_b1_bc14 == 1),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 30 & dt$l10_b1_bc14 == 1),"l10_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b1_bc05 > 50 & dt$l10_b1_bc14 == 0),"l10_b1_bc05"]+.1,dt[which(dt$l10_b1_bc05 > 50 & dt$l10_b1_bc14 == 0),"l10_b1_bc14"]+.1,pch=20,col="coral")
AAseqs.l10 <- unique(c(AAseqs.l10,as.character(dt[c(which(dt$l10_b1_bc05 > 50 & dt$l10_b1_bc14 == 0),which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 30 & dt$l10_b1_bc14 == 1),which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 20 & dt$l10_b1_bc14 == 2),which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 12 & dt$l10_b1_bc14 > 2),which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 10 & dt$l10_b1_bc14 > 5),which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 6 & dt$l10_b1_bc14 > 10),which(dt$l10_b1_bc05/dt$l10_b1_bc14 > 4 & dt$l10_b1_bc14 > 20)),1])))

# plot(dt$l10_b1_bc12+.1,dt$l10_b1_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# points(dt[dt$AAseq %in% AAseqs.l10, "l10_b1_bc12"]+.1,dt[dt$AAseq %in% AAseqs.l10, "l10_b1_bc14"]+.1,pch=19,col="coral") #outlier points already selected

# #l10 bin2, dt
# plot(dt$l10_b2_bc08+.1,dt$l10_b2_bc11.x+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# points(dt[which(dt$l10_b2_bc08/dt$l10_b2_bc11.x > 2.5 & dt$l10_b2_bc11.x > 6),"l10_b2_bc08"]+.1,dt[which(dt$l10_b2_bc08/dt$l10_b2_bc11.x > 2.5 & dt$l10_b2_bc11.x > 6),"l10_b2_bc11.x"]+.1,pch=20,col="coral")
# points(dt[which(dt$l10_b2_bc08/dt$l10_b2_bc11.x > 5 & dt$l10_b2_bc11.x == 5),"l10_b2_bc08"]+.1,dt[which(dt$l10_b2_bc08/dt$l10_b2_bc11.x > 5 & dt$l10_b2_bc11.x == 5),"l10_b2_bc11.x"]+.1,pch=20,col="coral")
AAseqs.l10 <- unique(c(AAseqs.l10,as.character(dt[c(which(dt$l10_b2_bc08/dt$l10_b2_bc11.x > 5 & dt$l10_b2_bc11.x == 5),which(dt$l10_b2_bc08/dt$l10_b2_bc11.x > 2.5 & dt$l10_b2_bc11.x > 6)),1])))

# plot(dt$l10_b2_bc08+.1,dt$l10_b2_bc15+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# points(dt[dt$AAseq %in% AAseqs.l10, "l10_b2_bc08"]+.1,dt[dt$AAseq %in% AAseqs.l10, "l10_b2_bc15"]+.1,pch=19,col="coral",cex=0.5) #outliers already accounted for

# #l10 bin2, dt
# plot(dt$l10_b2_bc01+.1,dt$l10_b2_bc05+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(a=0,b=1,col="red")
# plot(dt$l10_b2_bc01+.1,dt$l10_b2_bc11.y+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(a=0,b=1,col="red")
# #bc01 has weird behavior... if I don't use bc01 reads entirely, I should still have plenty sufficient reads from this bin to supplement dt.
# plot(dt$l10_b2_bc05+.1,dt$l10_b2_bc11.y+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# points(dt[dt$AAseq %in% AAseqs.l10, "l10_b2_bc05"]+.1,dt[dt$AAseq %in% AAseqs.l10, "l10_b2_bc11.y"]+.1,pch=19,col="coral") #outliers already accounted for

# #l11 (SRE rep1) bin 1
# plot(dt$l11_b1_bc05+.1,dt$l11_b1_bc07+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 5 & dt$l11_b1_bc07 > 10),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 5 & dt$l11_b1_bc07 > 10),"l11_b1_bc07"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 7 & dt$l11_b1_bc07 > 5),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 7 & dt$l11_b1_bc07 > 5),"l11_b1_bc07"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 12 & dt$l11_b1_bc07 > 1),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 12 & dt$l11_b1_bc07 > 1),"l11_b1_bc07"]+.1,pch=20,col="coral")
AAseqs.l11 <- unique(as.character(dt[c(which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 12 & dt$l11_b1_bc07 > 1),which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 7 & dt$l11_b1_bc07 > 5),which(dt$l11_b1_bc05/dt$l11_b1_bc07 > 5 & dt$l11_b1_bc07 > 10)),1]))

# plot(dt$l11_b1_bc05+.1,dt$l11_b1_bc13+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# abline(h=10);abline(v=10)
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 3 & dt$l11_b1_bc13 > 20),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 3 & dt$l11_b1_bc13 > 20),"l11_b1_bc13"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 5 & dt$l11_b1_bc13 > 10),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 5 & dt$l11_b1_bc13 > 10),"l11_b1_bc13"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 8 & dt$l11_b1_bc13 > 3),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 8 & dt$l11_b1_bc13 > 3),"l11_b1_bc13"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 25 & dt$l11_b1_bc13 > 0),"l11_b1_bc05"]+.1,dt[which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 25 & dt$l11_b1_bc13 > 0),"l11_b1_bc13"]+.1,pch=20,col="coral")
AAseqs.l11 <- unique(c(AAseqs.l11,as.character(dt[c(which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 25 & dt$l11_b1_bc13 > 0),which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 8 & dt$l11_b1_bc13 > 3),which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 5 & dt$l11_b1_bc13 > 10),which(dt$l11_b1_bc05/dt$l11_b1_bc13 > 3 & dt$l11_b1_bc13 > 20)),1])))

# plot(dt$l11_b1_bc07+.1,dt$l11_b1_bc13+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no clear outliers
# points(dt[dt$AAseq %in% AAseqs.l11,"l11_b1_bc07"]+.1,dt[dt$AAseq %in% AAseqs.l11,"l11_b1_bc13"]+.1,pch=20,col="coral")

# #l11 (SRE rep1) bin 2
# plot(dt$l11_b2_bc08+.1,dt$l11_b2_bc09+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 3 & dt$l11_b2_bc09 > 30),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 3 & dt$l11_b2_bc09 > 30),"l11_b2_bc09"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 5 & dt$l11_b2_bc09 > 20),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 5 & dt$l11_b2_bc09 > 20),"l11_b2_bc09"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 7 & dt$l11_b2_bc09 > 5),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 7 & dt$l11_b2_bc09 > 5),"l11_b2_bc09"]+.1,pch=20,col="coral")
AAseqs.l11 <- unique(c(AAseqs.l11,as.character(dt[c(which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 7 & dt$l11_b2_bc09 > 5),which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 5 & dt$l11_b2_bc09 > 20),which(dt$l11_b2_bc08/dt$l11_b2_bc09 > 3 & dt$l11_b2_bc09 > 30)),1])))

# plot(dt$l11_b2_bc08+.1,dt$l11_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 3 & dt$l11_b2_bc16 > 30),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 3 & dt$l11_b2_bc16 > 30),"l11_b2_bc16"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 5 & dt$l11_b2_bc16 > 20),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 5 & dt$l11_b2_bc16 > 20),"l11_b2_bc16"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 8 & dt$l11_b2_bc16 > 5),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 8 & dt$l11_b2_bc16 > 5),"l11_b2_bc16"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 15 & dt$l11_b2_bc16 > 1),"l11_b2_bc08"]+.1,dt[which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 15 & dt$l11_b2_bc16 > 1),"l11_b2_bc16"]+.1,pch=20,col="coral")
AAseqs.l11 <- unique(c(AAseqs.l11,as.character(dt[c(which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 15 & dt$l11_b2_bc16 > 1),which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 8 & dt$l11_b2_bc16 > 5),which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 5 & dt$l11_b2_bc16 > 20),which(dt$l11_b2_bc08/dt$l11_b2_bc16 > 3 & dt$l11_b2_bc16 > 30)),1])))

# plot(dt$l11_b2_bc09+.1,dt$l11_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 3.5 & dt$l11_b2_bc16 > 30),"l11_b2_bc09"]+.1,dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 3.5 & dt$l11_b2_bc16 > 30),"l11_b2_bc16"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 5 & dt$l11_b2_bc16 > 20),"l11_b2_bc09"]+.1,dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 5 & dt$l11_b2_bc16 > 20),"l11_b2_bc16"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 8 & dt$l11_b2_bc16 > 5),"l11_b2_bc09"]+.1,dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 8 & dt$l11_b2_bc16 > 5),"l11_b2_bc16"]+.1,pch=20,col="coral")
# points(dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 20 & dt$l11_b2_bc16 > 1),"l11_b2_bc09"]+.1,dt[which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 20 & dt$l11_b2_bc16 > 1),"l11_b2_bc16"]+.1,pch=20,col="coral")
AAseqs.l11 <- unique(c(AAseqs.l11,as.character(dt[c(which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 20 & dt$l11_b2_bc16 > 1),which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 8 & dt$l11_b2_bc16 > 5),which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 5 & dt$l11_b2_bc16 > 20),which(dt$l11_b2_bc09/dt$l11_b2_bc16 > 3.5 & dt$l11_b2_bc16 > 30)),1])))

# #l12 (ERE rep2) bin 1
# plot(dt$l12_b1_bc01+.1,dt$l12_b1_bc03+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no evidence for PCR bias

# plot(dt$l12_b1_bc01+.1,dt$l12_b1_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10)
# points(dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 3 & dt$l12_b1_bc14 > 30),"l12_b1_bc01"]+.1,dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 3 & dt$l12_b1_bc14 > 30),"l12_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 4.5 & dt$l12_b1_bc14 > 20),"l12_b1_bc01"]+.1,dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 4.5 & dt$l12_b1_bc14 > 20),"l12_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 7 & dt$l12_b1_bc14 > 10),"l12_b1_bc01"]+.1,dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 7 & dt$l12_b1_bc14 > 10),"l12_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 10 & dt$l12_b1_bc14 > 5),"l12_b1_bc01"]+.1,dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 10 & dt$l12_b1_bc14 > 5),"l12_b1_bc14"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 20 & dt$l12_b1_bc14 > 0),"l12_b1_bc01"]+.1,dt[which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 20 & dt$l12_b1_bc14 > 0),"l12_b1_bc14"]+.1,pch=20,col="coral")
AAseqs.l12 <- unique(as.character(dt[c(which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 20 & dt$l12_b1_bc14 > 0),which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 10 & dt$l12_b1_bc14 > 5),which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 7 & dt$l12_b1_bc14 > 10),which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 4.5 & dt$l12_b1_bc14 > 20),which(dt$l12_b1_bc01/dt$l12_b1_bc14 > 3 & dt$l12_b1_bc14 > 30)),1]))

# plot(dt$l12_b1_bc03+.1,dt$l12_b1_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# points(dt[dt$AAseq %in% AAseqs.l12,"l12_b1_bc03"]+.1,dt[dt$AAseq %in% AAseqs.l12,"l12_b1_bc14"]+.1,pch=20,col="coral")

# #l12 (ERE rep1) bin 2
# plot(dt$l12_b2_bc02+.1,dt$l12_b2_bc06+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no evidence for PCR bias
# abline(h=10);abline(v=10)
# points(dt[which(dt$l12_b2_bc02/dt$l12_b2_bc06 > 5 & dt$l12_b2_bc06 > 10),"l12_b2_bc02"]+.1,dt[which(dt$l12_b2_bc02/dt$l12_b2_bc06 > 5 & dt$l12_b2_bc06 > 10),"l12_b2_bc06"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b2_bc02/dt$l12_b2_bc06 > 15 & dt$l12_b2_bc06 > 1),"l12_b2_bc02"]+.1,dt[which(dt$l12_b2_bc02/dt$l12_b2_bc06 > 15 & dt$l12_b2_bc06 > 1),"l12_b2_bc06"]+.1,pch=20,col="coral")
AAseqs.l12 <- unique(c(AAseqs.l12,as.character(dt[c(which(dt$l12_b2_bc02/dt$l12_b2_bc06 > 5 & dt$l12_b2_bc06 > 10),which(dt$l12_b2_bc02/dt$l12_b2_bc06 > 15 & dt$l12_b2_bc06 > 1)),1])))

# plot(dt$l12_b2_bc02+.1,dt$l12_b2_bc11+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no evidence for PCR bias
# abline(h=10);abline(v=10)
# points(dt[which(dt$l12_b2_bc02/dt$l12_b2_bc11 > 5 & dt$l12_b2_bc11 > 10),"l12_b2_bc02"]+.1,dt[which(dt$l12_b2_bc02/dt$l12_b2_bc11 > 5 & dt$l12_b2_bc11 > 10),"l12_b2_bc11"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b2_bc02/dt$l12_b2_bc11 > 15 & dt$l12_b2_bc11 > 1),"l12_b2_bc02"]+.1,dt[which(dt$l12_b2_bc02/dt$l12_b2_bc11 > 15 & dt$l12_b2_bc11 > 1),"l12_b2_bc11"]+.1,pch=20,col="coral")
AAseqs.l12 <- unique(c(AAseqs.l12,as.character(dt[c(which(dt$l12_b2_bc02/dt$l12_b2_bc11 > 15 & dt$l12_b2_bc11 > 1),which(dt$l12_b2_bc02/dt$l12_b2_bc11 > 5 & dt$l12_b2_bc11 > 10)),1])))

# plot(dt$l12_b2_bc06+.1,dt$l12_b2_bc11+.1,pch=20,cex=0.5,log="xy",col="#00000066") #no evidence for PCR bias
# abline(h=10);abline(v=10)
# points(dt[which(dt$l12_b2_bc06/dt$l12_b2_bc11 > 4 & dt$l12_b2_bc11 > 10),"l12_b2_bc06"]+.1,dt[which(dt$l12_b2_bc06/dt$l12_b2_bc11 > 4 & dt$l12_b2_bc11 > 10),"l12_b2_bc11"]+.1,pch=20,col="coral")
# points(dt[which(dt$l12_b2_bc06/dt$l12_b2_bc11 > 15 & dt$l12_b2_bc11 > 1),"l12_b2_bc06"]+.1,dt[which(dt$l12_b2_bc06/dt$l12_b2_bc11 > 15 & dt$l12_b2_bc11 > 1),"l12_b2_bc11"]+.1,pch=20,col="coral")
AAseqs.l12 <- unique(c(AAseqs.l12,as.character(dt[c(which(dt$l12_b2_bc06/dt$l12_b2_bc11 > 15 & dt$l12_b2_bc11 > 1),which(dt$l12_b2_bc06/dt$l12_b2_bc11 > 4 & dt$l12_b2_bc11 > 10)),1])))


#remove variants with evidence of PCR bias
dt[dt$AAseq %in% AAseqs.l10,c("l10_b1_bc05","l10_b1_bc12","l10_b1_bc14","l10_b2_bc08","l10_b2_bc11.x","l10_b2_bc15","l10_b3_bc02","l10_b4_bc04","l10_b1_bc06","l10_b2_bc01","l10_b2_bc05","l10_b2_bc11.y","l10_b3_bc12")] <- 0
dt[dt$AAseq %in% AAseqs.l11,c("l11_b1_bc05","l11_b1_bc07","l11_b1_bc13","l11_b2_bc08","l11_b2_bc09","l11_b2_bc16","l11_b3_bc15","l11_b4_bc12")] <- 0
dt[dt$AAseq %in% AAseqs.l12,c("l12_b1_bc01","l12_b1_bc03","l12_b1_bc14","l12_b2_bc02","l12_b2_bc06","l12_b2_bc11","l12_b3_bc04","l12_b4_bc10")] <- 0

rm(AAseqs.l10);rm(AAseqs.l11);rm(AAseqs.l12)

#colony counts, estimated from serial dilution and plating after cell sorting
cfu.l10.b1 <- 20403126; cfu.l10.b2 <- 31752091; cfu.l10.b3 <- 5254508; cfu.l10.b4 <- 218800
cfu.l9.b1 <- 20909000; cfu.l9.b2 <- 35685000; cfu.l9.b3 <- 5404565; cfu.l9.b4 <- 132735
cfu.l12.b1 <- 15797500; cfu.l12.b2 <- 16107000; cfu.l12.b3 <- 2817900; cfu.l12.b4 <- 156000
cfu.l11.b1 <- 20344400; cfu.l11.b2 <- 30681600; cfu.l11.b3 <- 5538600; cfu.l11.b4 <- 585900

#pool together individual barcodes representing identical bins
dt$l10.b1 <- rowSums(dt[,c("l10_b1_bc05","l10_b1_bc12","l10_b1_bc14","l10_b1_bc06")]);dt$l10.b2 <- rowSums(dt[,c("l10_b2_bc08","l10_b2_bc11.x","l10_b2_bc15","l10_b2_bc05","l10_b2_bc11.y")]);dt$l10.b3 <- rowSums(dt[,c("l10_b3_bc02","l10_b3_bc12")]);dt$l10.b4 <- dt[,"l10_b4_bc04"] #no bc01
dt$l9.b1 <- rowSums(dt[,c("l9_b1_bc01","l9_b1_bc06","l9_b1_bc10","l9_b1_bc15")]);dt$l9.b2 <- rowSums(dt[,c("l9_b2_bc03","l9_b2_bc09","l9_b2_bc16","l9_b2_bc02","l9_b2_bc07","l9_b2_bc14")]);dt$l9.b3 <- rowSums(dt[,c("l9_b3_bc07","l9_b3_bc13")]);dt$l9.b4 <- dt[,"l9_b4_bc13"]
dt$l12.b1 <- rowSums(dt[,c("l12_b1_bc01","l12_b1_bc03","l12_b1_bc14")]);dt$l12.b2 <- rowSums(dt[,c("l12_b2_bc02","l12_b2_bc06","l12_b2_bc11")]);dt$l12.b3 <- dt[,"l12_b3_bc04"]; dt$l12.b4 <- dt[,"l12_b4_bc10"]
dt$l11.b1 <- rowSums(dt[,c("l11_b1_bc05","l11_b1_bc07","l11_b1_bc13")]);dt$l11.b2 <- rowSums(dt[,c("l11_b2_bc08","l11_b2_bc09","l11_b2_bc16")]);dt$l11.b3 <- dt[,"l11_b3_bc15"];dt$l11.b4 <- dt[,"l11_b4_bc12"]

#normalize # of reads by number of cfus in bin
dt$ERE.rep1.b1 <- dt$l10.b1/sum(dt$l10.b1)*cfu.l10.b1;dt$ERE.rep1.b2 <- dt$l10.b2/sum(dt$l10.b2)*cfu.l10.b2;dt$ERE.rep1.b3 <- dt$l10.b3/sum(dt$l10.b3)*cfu.l10.b3;dt$ERE.rep1.b4 <- dt$l10.b4/sum(dt$l10.b4)*cfu.l10.b4
dt$SRE.rep1.b1 <- dt$l9.b1/sum(dt$l9.b1)*cfu.l9.b1;dt$SRE.rep1.b2 <- dt$l9.b2/sum(dt$l9.b2)*cfu.l9.b2;dt$SRE.rep1.b3 <- dt$l9.b3/sum(dt$l9.b3)*cfu.l9.b3;dt$SRE.rep1.b4 <- dt$l9.b4/sum(dt$l9.b4)*cfu.l9.b4
dt$ERE.rep2.b1 <- dt$l12.b1/sum(dt$l12.b1)*cfu.l12.b1;dt$ERE.rep2.b2 <- dt$l12.b2/sum(dt$l12.b2)*cfu.l12.b2;dt$ERE.rep2.b3 <- dt$l12.b3/sum(dt$l12.b3)*cfu.l12.b3;dt$ERE.rep2.b4 <- dt$l12.b4/sum(dt$l12.b4)*cfu.l12.b4
dt$SRE.rep2.b1 <- dt$l11.b1/sum(dt$l11.b1)*cfu.l11.b1;dt$SRE.rep2.b2 <- dt$l11.b2/sum(dt$l11.b2)*cfu.l11.b2;dt$SRE.rep2.b3 <- dt$l11.b3/sum(dt$l11.b3)*cfu.l11.b3;dt$SRE.rep2.b4 <- dt$l11.b4/sum(dt$l11.b4)*cfu.l11.b4

#sum total number of weighted reads per variant in each rep
dt$ERE.rep1.cfu <- rowSums(dt[,c("ERE.rep1.b1","ERE.rep1.b2","ERE.rep1.b3","ERE.rep1.b4")])
dt$SRE.rep1.cfu <- rowSums(dt[,c("SRE.rep1.b1","SRE.rep1.b2","SRE.rep1.b3","SRE.rep1.b4")])
dt$ERE.rep2.cfu <- rowSums(dt[,c("ERE.rep2.b1","ERE.rep2.b2","ERE.rep2.b3","ERE.rep2.b4")])
dt$SRE.rep2.cfu <- rowSums(dt[,c("SRE.rep2.b1","SRE.rep2.b2","SRE.rep2.b3","SRE.rep2.b4")])
dt$ERE.pooled.cfu <- rowSums(dt[,c("ERE.rep1.cfu","ERE.rep2.cfu")])
dt$SRE.pooled.cfu <- rowSums(dt[,c("SRE.rep1.cfu","SRE.rep2.cfu")])

#convert to data.table for faster computation and indexing
dt <- data.table(dt[,c("AAseq","ERE.rep1.b1","ERE.rep1.b2","ERE.rep1.b3","ERE.rep1.b4","SRE.rep1.b1","SRE.rep1.b2","SRE.rep1.b3","SRE.rep1.b4","ERE.rep2.b1","ERE.rep2.b2","ERE.rep2.b3","ERE.rep2.b4","SRE.rep2.b1","SRE.rep2.b2","SRE.rep2.b3","SRE.rep2.b4","ERE.rep1.cfu","ERE.rep2.cfu","ERE.pooled.cfu","SRE.rep1.cfu","SRE.rep2.cfu","SRE.pooled.cfu")]);setkey(dt,AAseq)

#rescale fluorescence scales to +11P reference scales (ERE rep1=l6; SRE rep2=l8), so I can combine observations to calculate pooled meanF and compare meanF between experiments
ERE.l6.ctrl <- read.csv(file="./1_calc-meanF_11P_in/FITC_lib6-controls.csv",header=T);ERE.l10.ctrl <- read.csv(file="./2_calc-meanF_SR1_in/FITC_lib10-controls.csv",header=T);ERE.l12.ctrl <- read.csv(file="./2_calc-meanF_SR1_in/FITC_lib12-controls.csv",header=T)
SRE.l8.ctrl <- read.csv(file="./1_calc-meanF_11P_in/FITC_lib8-controls.csv",header=T);SRE.l9.ctrl <- read.csv(file="./2_calc-meanF_SR1_in/FITC_lib9-controls.csv",header=T);SRE.l11.ctrl <- read.csv(file="./2_calc-meanF_SR1_in/FITC_lib11-controls.csv",header=T)

ERE.ctrls <- data.frame(geno=c("lib","null","egka","EGKA","GGKA","GSKV","GGKV","GSKA"))
for(i in 1:nrow(ERE.ctrls)){
  ERE.ctrls$l6.meanF[i] <- mean(log(ERE.l6.ctrl[ERE.l6.ctrl[,as.character(ERE.ctrls$geno[i])] > 0,as.character(ERE.ctrls$geno[i])]),na.rm=T)
  ERE.ctrls$l10.meanF[i] <- mean(log(ERE.l10.ctrl[ERE.l10.ctrl[,as.character(ERE.ctrls$geno[i])] > 0,as.character(ERE.ctrls$geno[i])]),na.rm=T)
  ERE.ctrls$l12.meanF[i] <- mean(log(ERE.l12.ctrl[ERE.l12.ctrl[,as.character(ERE.ctrls$geno[i])] > 0,as.character(ERE.ctrls$geno[i])]),na.rm=T)
}
plot(ERE.ctrls$l6.meanF,ERE.ctrls$l10.meanF)
lm.ERE.ctrl <- lm(ERE.ctrls$l10.meanF~ERE.ctrls$l6.meanF);abline(lm.ERE.ctrl);summary(lm.ERE.ctrl)
#l10 = 1.01805*l6 - 0.59618

plot(ERE.ctrls$l6.meanF,ERE.ctrls$l12.meanF)
lm.ERE.ctrl <- lm(ERE.ctrls$l12.meanF~ERE.ctrls$l6.meanF);abline(lm.ERE.ctrl);summary(lm.ERE.ctrl)
#l12 = 0.9633*l6 + 0.5260


SRE.ctrls <- data.frame(geno=c("lib","null","GSKV","GGKA","EGKA","egka","GGKV","GSKA"))
for(i in 1:nrow(SRE.ctrls)){
  SRE.ctrls$l8.meanF[i] <- mean(log(SRE.l8.ctrl[SRE.l8.ctrl[,as.character(SRE.ctrls$geno[i])] > 0,as.character(SRE.ctrls$geno[i])]),na.rm=T)
  SRE.ctrls$l9.meanF[i] <- mean(log(SRE.l9.ctrl[SRE.l9.ctrl[,as.character(SRE.ctrls$geno[i])] > 0,as.character(SRE.ctrls$geno[i])]),na.rm=T)
  SRE.ctrls$l11.meanF[i] <- mean(log(SRE.l11.ctrl[SRE.l11.ctrl[,as.character(SRE.ctrls$geno[i])] > 0,as.character(SRE.ctrls$geno[i])]),na.rm=T)
}
plot(SRE.ctrls$l8.meanF,SRE.ctrls$l9.meanF)
lm.SRE.ctrl <- lm(SRE.ctrls$l9.meanF~SRE.ctrls$l8.meanF);abline(lm.SRE.ctrl);summary(lm.SRE.ctrl)
#l9 = 1.17925*l8 - 0.82284

plot(SRE.ctrls$l8.meanF,SRE.ctrls$l11.meanF)
lm.SRE.ctrl <- lm(SRE.ctrls$l11.meanF~SRE.ctrls$l8.meanF);abline(lm.SRE.ctrl);summary(lm.SRE.ctrl)
#l11 = 1.11911*l8 - 0.57030

#give fluor boundaries of sort bins, rescaled to 11P reference scales
min.l10.b1 <- log(1); min.l10.b2 <- (log(73.5)+0.59618)/1.01805; min.l10.b3 <- (log(170.5)+0.59618)/1.01805; min.l10.b4 <- (log(374.5)+0.59618)/1.01805; max.l10.b4 <- log(262144)
min.l9.b1 <- log(1); min.l9.b2 <- (log(149.5)+0.82284)/1.17925; min.l9.b3 <- (log(483.5)+0.82284)/1.17925; min.l9.b4 <- (log(1147.5)+0.82284)/1.17925; max.l9.b4 <- log(262144)
min.l12.b1 <- log(1); min.l12.b2 <- (log(175.5)-0.5260)/0.9633; min.l12.b3 <- (log(374.5)-0.5260)/0.9633; min.l12.b4 <- (log(933.5)-0.5260)/0.9633; max.l12.b4 <- log(262144)
min.l11.b1 <- log(1); min.l11.b2 <- (log(121)+0.57030)/1.11911; min.l11.b3 <- (log(455)+0.57030)/1.11911; min.l11.b4 <- (log(1136.5)+0.57030)/1.11911; max.l11.b4 <- log(262144)

#define functions to estimate ML mean and sd from logistic distribution for individual reps, pooled reps
fit.bin.logistic <- function(b1,b2,b3,b4,min.b1,min.b2,min.b3,min.b4,max.b4){
  data <- data.frame(left=c(rep(min.b1,round(b1)),rep(min.b2,round(b2)),rep(min.b3,round(b3)),rep(min.b4,round(b4))),right=c(rep(min.b2,round(b1)),rep(min.b3,round(b2)),rep(min.b4,round(b3)),rep(max.b4,round(b4))))
  if(nrow(unique(data))>1){
    fit <- fitdistcens(data,"logis")
    return(list(as.numeric(summary(fit)$estimate["location"]),sqrt((as.numeric(summary(fit)$estimate["scale"])^2*pi^2)/3)))
  } else {
    return(list(as.numeric(NA),as.numeric(NA)))
  }
}

fit.bin.logistic.pooled <- function(b1.1,b2.1,b3.1,b4.1,b1.2,b2.2,b3.2,b4.2,min.b1.1,min.b2.1,min.b3.1,min.b4.1,max.b4.1,min.b1.2,min.b2.2,min.b3.2,min.b4.2,max.b4.2){
  data <- data.frame(left=c(rep(min.b1.1,round(b1.1)),rep(min.b2.1,round(b2.1)),rep(min.b3.1,round(b3.1)),rep(min.b4.1,round(b4.1)),rep(min.b1.2,round(b1.2)),rep(min.b2.2,round(b2.2)),rep(min.b3.2,round(b3.2)),rep(min.b4.2,round(b4.2))),right=c(rep(min.b2.1,round(b1.1)),rep(min.b3.1,round(b2.1)),rep(min.b4.1,round(b3.1)),rep(max.b4.1,round(b4.1)),rep(min.b2.2,round(b1.2)),rep(min.b3.2,round(b2.2)),rep(min.b4.2,round(b3.2)),rep(max.b4.2,round(b4.2))))
  fit <- fitdistcens(data,"logis")
  return(list(as.numeric(summary(fit)$estimate["location"]),sqrt((as.numeric(summary(fit)$estimate["scale"])^2*pi^2)/3)))
}

#calculate ML meanF, sdF for each replicate, pooled replicates
dt[,c("ERE.rep1.meanF","ERE.rep1.sdF") := fit.bin.logistic(ERE.rep1.b1,ERE.rep1.b2,ERE.rep1.b3,ERE.rep1.b4,min.l10.b1,min.l10.b2,min.l10.b3,min.l10.b4,max.l10.b4),by=AAseq]
dt[,c("ERE.rep2.meanF","ERE.rep2.sdF") := fit.bin.logistic(ERE.rep2.b1,ERE.rep2.b2,ERE.rep2.b3,ERE.rep2.b4,min.l12.b1,min.l12.b2,min.l12.b3,min.l12.b4,max.l12.b4),by=AAseq]
dt[,c("ERE.pooled.meanF","ERE.pooled.sdF") := tryCatch(fit.bin.logistic.pooled(ERE.rep1.b1,ERE.rep1.b2,ERE.rep1.b3,ERE.rep1.b4,ERE.rep2.b1,ERE.rep2.b2,ERE.rep2.b3,ERE.rep2.b4,min.l10.b1,min.l10.b2,min.l10.b3,min.l10.b4,max.l10.b4,min.l12.b1,min.l12.b2,min.l12.b3,min.l12.b4,max.l12.b4),error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=AAseq]
dt[,c("SRE.rep1.meanF","SRE.rep1.sdF") := fit.bin.logistic(SRE.rep1.b1,SRE.rep1.b2,SRE.rep1.b3,SRE.rep1.b4,min.l9.b1,min.l9.b2,min.l9.b3,min.l9.b4,max.l9.b4),by=AAseq]
dt[,c("SRE.rep2.meanF","SRE.rep2.sdF") := fit.bin.logistic(SRE.rep2.b1,SRE.rep2.b2,SRE.rep2.b3,SRE.rep2.b4,min.l11.b1,min.l11.b2,min.l11.b3,min.l11.b4,max.l11.b4),by=AAseq]
dt[,c("SRE.pooled.meanF","SRE.pooled.sdF") := tryCatch(fit.bin.logistic.pooled(SRE.rep1.b1,SRE.rep1.b2,SRE.rep1.b3,SRE.rep1.b4,SRE.rep2.b1,SRE.rep2.b2,SRE.rep2.b3,SRE.rep2.b4,min.l9.b1,min.l9.b2,min.l9.b3,min.l9.b4,max.l9.b4,min.l11.b1,min.l11.b2,min.l11.b3,min.l11.b4,max.l11.b4),error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=AAseq]
write.table(dt, file="2_calc-meanF_SR1_out/dt_output_SR1.csv",col.names=T,row.names=F,quote=F,sep=",")

#what variants are undetermined in pooled rep, but determined in either or both indiivdual reps?
#919 variants determined in ERE rep1, undetermined pooled ERE --> 
#910  have >1cfu in l10.b1, l10.b2, and l12.b1, 0 in all others --> assign pooled.meanF to be same as ERE.rep1.meanF
dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF & ERE.rep1.b1>0.5 & ERE.rep1.b2>0.5 & ERE.rep2.b1>0.5 & ERE.rep1.b3<0.5&ERE.rep1.b4<0.5&ERE.rep2.b2<0.5 & ERE.rep2.b3<0.5 & ERE.rep2.b4<0.5),c(1:5,10:13,18:20,24,26,28),with=FALSE];dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF & ERE.rep1.b1>0.5 & ERE.rep1.b2>0.5 & ERE.rep2.b1>0.5 & ERE.rep1.b3<0.5&ERE.rep1.b4<0.5&ERE.rep2.b2<0.5 & ERE.rep2.b3<0.5 & ERE.rep2.b4<0.5),ERE.pooled.meanF:=ERE.rep1.meanF]
#8 have >1cfu in l10.b2, l10.b3, and l12.b3, 0 in all others --> assign pooled.meanF to be same as ERE.rep1.meanF
dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF & ERE.rep1.b1<0.5 & ERE.rep1.b2>0.5 & ERE.rep1.b3>0.5 & ERE.rep1.b4<0.5&ERE.rep2.b1<0.5&ERE.rep2.b2<0.5 & ERE.rep2.b3>0.5 & ERE.rep2.b4<0.5),c(1:5,10:13,18:20,24,26,28),with=FALSE];dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF & ERE.rep1.b1<0.5 & ERE.rep1.b2>0.5 & ERE.rep1.b3>0.5 & ERE.rep1.b4<0.5&ERE.rep2.b1<0.5&ERE.rep2.b2<0.5 & ERE.rep2.b3>0.5 & ERE.rep2.b4<0.5),ERE.pooled.meanF:=ERE.rep1.meanF]
#remaining one is determined for both ERE.rep1.meanF and ERE.rep2.meanF --> give ERE.pooled.meanF as mean of two determinations, weighted by cfu from each rep
dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF),ERE.pooled.meanF := (ERE.rep1.meanF * ERE.rep1.cfu + ERE.rep2.meanF * ERE.rep2.cfu)/(ERE.rep1.cfu+ERE.rep2.cfu)]

#1319 variants determined in ERE rep2, undetermined pooled ERE -->
#1254 have >1cfu in l10.b2, l12.b1, l12.b2, <1 in all others --> assign pooled.meanF to be same as ERE.rep2.meanF
dt[!is.na(ERE.rep2.meanF) & is.na(ERE.pooled.meanF) & ERE.rep1.b1<0.5 & ERE.rep1.b2>0.5 & ERE.rep1.b3<0.5 & ERE.rep1.b4<0.5&ERE.rep2.b1>0.5&ERE.rep2.b2>0.5 & ERE.rep2.b3<0.5 & ERE.rep2.b4<0.5,c(1:5,10:13,18:20,24,26,28),with=FALSE];dt[!is.na(ERE.rep2.meanF) & is.na(ERE.pooled.meanF)& ERE.rep1.b1<0.5 & ERE.rep1.b2>0.5 & ERE.rep1.b3<0.5 & ERE.rep1.b4<0.5&ERE.rep2.b1>0.5&ERE.rep2.b2>0.5 & ERE.rep2.b3<0.5 & ERE.rep2.b4<0.5,ERE.pooled.meanF := ERE.rep2.meanF]
#remaining 65 have >1cfu in l10.b2, l12.b2, l12.b3, <1 in all others --> assign pooled.meanF to be same as ERE.rep2.meanF
dt[!is.na(ERE.rep2.meanF) & is.na(ERE.pooled.meanF) & ERE.rep1.b1<0.5 & ERE.rep1.b2>0.5 & ERE.rep1.b3<0.5 & ERE.rep1.b4<0.5&ERE.rep2.b1<0.5&ERE.rep2.b2>0.5 & ERE.rep2.b3>0.5 & ERE.rep2.b4<0.5,c(1:5,10:13,18:20,24,26,28),with=FALSE];dt[!is.na(ERE.rep2.meanF) & is.na(ERE.pooled.meanF)& ERE.rep1.b1<0.5 & ERE.rep1.b2>0.5 & ERE.rep1.b3<0.5 & ERE.rep1.b4<0.5&ERE.rep2.b1<0.5&ERE.rep2.b2>0.5 & ERE.rep2.b3>0.5 & ERE.rep2.b4<0.5,ERE.pooled.meanF := ERE.rep2.meanF]

#1 variant determined in SRE rep1, undetermined pooled SRE --> >1 cfu l9.b2, l9.b3, l11.b2 --> assign SRE.pooled.meanF to be same as SRE.rep1.meanF
dt[!is.na(SRE.rep1.meanF) & is.na(SRE.pooled.meanF),c(1,6:9,14:17,21:23,30,32,34),with=FALSE];dt[!is.na(SRE.rep1.meanF) & is.na(SRE.pooled.meanF),SRE.pooled.meanF:=SRE.rep1.meanF]

#2 variants determined in SRE rep2, undetermined pooled SRE --> >1 cfu l9.b3, l1.b2, l11.b3 --> assign SRE.pooled.meanF to be same as SRE.rep2.meanF
dt[!is.na(SRE.rep2.meanF) & is.na(SRE.pooled.meanF),c(1,6:9,14:17,21:23,30,32,34),with=FALSE];dt[!is.na(SRE.rep2.meanF) & is.na(SRE.pooled.meanF),SRE.pooled.meanF:=SRE.rep2.meanF]

write.table(dt, file="2_calc-meanF_SR1_out/dt_output_SR1.csv",col.names=T,row.names=F,quote=F,sep=",")
