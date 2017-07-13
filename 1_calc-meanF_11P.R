#16 Jan 2017
#TNS
#read in demultiplexed sequence reads, eliminate possible PCR artefacts (?), calculate ML meanF for each library genotype in AncSR1+11P background for ERE, SRE binding
#library codes: l6 = ERE, rep1 ; l13 = ERE, rep 2; l5 = SRE, rep 1; l8 = SRE, rep 2

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)
library(fitdistrplus)

#read input file, which is a table of all RH variants and the number of sequence reads mapped to that RH for each library/barcode combo
dt <- read.csv(file="1_calc-meanF_11P_in/data.in.11P.csv",header=T)

# #check for any variants with abnormal variation in # reads across duplicate bc PCRs == evidence for PCR bias? (only visible in b1/b2 which have multiple bcs/bin)
# #visually ID outliers, remove reads for that variant across all bins for that replicate
# #lib6 (ERE rep1) bin 1
# plot(dt$l6_b1_bc01+.1,dt$l6_b1_bc06+.1,pch=20,cex=0.5,log="xy",col="#00000066")

# plot(dt$l6_b1_bc01+.1,dt$l6_b1_bc10+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=30);abline(v=30)
# points(dt[which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 3 & dt$l6_b1_bc10>=30),"l6_b1_bc01"],dt[which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 3 & dt$l6_b1_bc10>=30),"l6_b1_bc10"],pch=20,col="coral")
# points(dt[which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 4 & dt$l6_b1_bc10>=10),"l6_b1_bc01"],dt[which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 4 & dt$l6_b1_bc10>=10),"l6_b1_bc10"],pch=20,col="coral")
# points(dt[which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 8 & dt$l6_b1_bc10>=3),"l6_b1_bc01"],dt[which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 8 & dt$l6_b1_bc10>=3),"l6_b1_bc10"],pch=20,col="coral")
# points(dt[which(dt$l6_b1_bc01 >= 20 & dt$l6_b1_bc10==0),"l6_b1_bc01"]+0.1,dt[which(dt$l6_b1_bc01 >= 20 & dt$l6_b1_bc10==0),"l6_b1_bc10"]+0.1,pch=20,col="coral")
AAseqs.l6 <- unique(as.character(dt[c(which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 3 & dt$l6_b1_bc10>=30),which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 4 & dt$l6_b1_bc10>=10),which(dt$l6_b1_bc01/dt$l6_b1_bc10 > 8 & dt$l6_b1_bc10>=3),which(dt$l6_b1_bc01 >= 20 & dt$l6_b1_bc10==0)),1]))

# plot(dt$l6_b1_bc06+.1,dt$l6_b1_bc10+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=25);abline(v=25)
# points(dt[which(dt$l6_b1_bc06/dt$l6_b1_bc10 > 3 & dt$l6_b1_bc10>=24),"l6_b1_bc06"],dt[which(dt$l6_b1_bc06/dt$l6_b1_bc10 > 3 & dt$l6_b1_bc10>=24),"l6_b1_bc10"],pch=20,col="coral")
# points(dt[which(dt$l6_b1_bc06/dt$l6_b1_bc10 > 8 & dt$l6_b1_bc10>=3),"l6_b1_bc06"],dt[which(dt$l6_b1_bc06/dt$l6_b1_bc10 > 8 & dt$l6_b1_bc10>=3),"l6_b1_bc10"],pch=20,col="coral")
# points(dt[which(dt$l6_b1_bc06 >= 20 & dt$l6_b1_bc10==0),"l6_b1_bc06"]+0.1,dt[which(dt$l6_b1_bc06 >= 20 & dt$l6_b1_bc10==0),"l6_b1_bc10"]+0.1,pch=20,col="coral")
AAseqs.l6 <- unique(c(AAseqs.l6,as.character(dt[c(which(dt$l6_b1_bc06/dt$l6_b1_bc10 > 3 & dt$l6_b1_bc10>=24),which(dt$l6_b1_bc06/dt$l6_b1_bc10 > 8 & dt$l6_b1_bc10>=3),which(dt$l6_b1_bc06 >= 20 & dt$l6_b1_bc10==0)),1])))

# #l6 (ERE rep1) bin2
# plot(dt$l6_b2_bc03+.1,dt$l6_b2_bc09+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=25);abline(v=25)
# points(dt[which(dt$l6_b2_bc03/dt$l6_b2_bc09 > 4 & dt$l6_b2_bc09>=24),"l6_b2_bc03"],dt[which(dt$l6_b2_bc03/dt$l6_b2_bc09 > 4 & dt$l6_b2_bc09>=24),"l6_b2_bc09"],pch=20,col="coral")
# points(dt[which(dt$l6_b2_bc03/dt$l6_b2_bc09 > 6 & dt$l6_b2_bc09>=10),"l6_b2_bc03"],dt[which(dt$l6_b2_bc03/dt$l6_b2_bc09 > 6 & dt$l6_b2_bc09>=10),"l6_b2_bc09"],pch=20,col="coral")
# points(dt[which(dt$l6_b2_bc03 >= 20 & dt$l6_b2_bc09==0),"l6_b2_bc03"]+0.1,dt[which(dt$l6_b2_bc03 >= 20 & dt$l6_b2_bc09==0),"l6_b2_bc09"]+0.1,pch=20,col="coral")
AAseqs.l6 <- unique(c(AAseqs.l6,as.character(dt[c(which(dt$l6_b2_bc03/dt$l6_b2_bc09 > 4 & dt$l6_b2_bc09>=24),which(dt$l6_b2_bc03/dt$l6_b2_bc09 > 6 & dt$l6_b2_bc09>=10),which(dt$l6_b2_bc03 >= 20 & dt$l6_b2_bc09==0)),1])))

# plot(dt$l6_b2_bc03+.1,dt$l6_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=25);abline(v=25)
# points(dt[which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 3 & dt$l6_b2_bc16>=35),"l6_b2_bc03"],dt[which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 3 & dt$l6_b2_bc16>=35),"l6_b2_bc16"],pch=20,col="coral")
# points(dt[which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 4 & dt$l6_b2_bc16>=10),"l6_b2_bc03"],dt[which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 4 & dt$l6_b2_bc16>=10),"l6_b2_bc16"],pch=20,col="coral")
# points(dt[which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 9 & dt$l6_b2_bc16>=4),"l6_b2_bc03"],dt[which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 9 & dt$l6_b2_bc16>=4),"l6_b2_bc16"],pch=20,col="coral")
# points(dt[which(dt$l6_b2_bc03 >= 20 & dt$l6_b2_bc16==0),"l6_b2_bc03"]+0.1,dt[which(dt$l6_b2_bc03 >= 20 & dt$l6_b2_bc16==0),"l6_b2_bc16"]+0.1,pch=20,col="coral")
AAseqs.l6 <- unique(c(AAseqs.l6,as.character(dt[c(which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 3 & dt$l6_b2_bc16>=35),which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 4 & dt$l6_b2_bc16>=10),which(dt$l6_b2_bc03/dt$l6_b2_bc16 > 9 & dt$l6_b2_bc16>=4),which(dt$l6_b2_bc03 >= 20 & dt$l6_b2_bc16==0)),1])))

# plot(dt$l6_b2_bc09+.1,dt$l6_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=25);abline(v=25)
# points(dt[which(dt$l6_b2_bc09/dt$l6_b2_bc16 > 3 & dt$l6_b2_bc16>=20),"l6_b2_bc09"],dt[which(dt$l6_b2_bc09/dt$l6_b2_bc16 > 3 & dt$l6_b2_bc16>=20),"l6_b2_bc16"],pch=20,col="coral")
# points(dt[which(dt$l6_b2_bc09 >= 20 & dt$l6_b2_bc16==0),"l6_b2_bc09"]+0.1,dt[which(dt$l6_b2_bc09 >= 20 & dt$l6_b2_bc16==0),"l6_b2_bc16"]+0.1,pch=20,col="coral")
AAseqs.l6 <- unique(c(AAseqs.l6,as.character(dt[c(which(dt$l6_b2_bc09 >= 20 & dt$l6_b2_bc16==0),which(dt$l6_b2_bc09/dt$l6_b2_bc16 > 3 & dt$l6_b2_bc16>=20)),1])))


# #l5 bin1 (SRE rep1)
# plot(dt$l5_b1_bc01+.1,dt$l5_b1_bc07+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=15);abline(v=15)
# points(dt[which(dt$l5_b1_bc01/dt$l5_b1_bc07 > 5 & dt$l5_b1_bc07>=15),"l5_b1_bc01"]+.1,dt[which(dt$l5_b1_bc01/dt$l5_b1_bc07 > 5 & dt$l5_b1_bc07>=15),"l5_b1_bc07"]+.1,pch=20,col="coral")
AAseqs.l5 <- unique(as.character(dt[which(dt$l5_b1_bc01/dt$l5_b1_bc07 > 5 & dt$l5_b1_bc07>=15),1]))

# plot(dt$l5_b1_bc01+.1,dt$l5_b1_bc13+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=15);abline(v=15)
# points(dt[which(dt$l5_b1_bc01/dt$l5_b1_bc13 > 4 & dt$l5_b1_bc13>=15),"l5_b1_bc01"]+.1,dt[which(dt$l5_b1_bc01/dt$l5_b1_bc13 > 4 & dt$l5_b1_bc13>=15),"l5_b1_bc13"]+.1,pch=20,col="coral")
# points(dt[which(dt$l5_b1_bc01 > 20 & dt$l5_b1_bc13==0),"l5_b1_bc01"]+.1,dt[which(dt$l5_b1_bc01 > 20 & dt$l5_b1_bc13==0),"l5_b1_bc13"]+.1,pch=20,col="coral")
AAseqs.l5 <- unique(c(AAseqs.l5,as.character(dt[which(dt$l5_b1_bc01/dt$l5_b1_bc13 > 3.5 & dt$l5_b1_bc13>=15),1])))

# plot(dt$l5_b1_bc07+.1,dt$l5_b1_bc13+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=15);abline(v=15)
# points(dt[which(dt$l5_b1_bc07/dt$l5_b1_bc13 > 3 & dt$l5_b1_bc13>=15),"l5_b1_bc07"]+.1,dt[which(dt$l5_b1_bc07/dt$l5_b1_bc13 > 3 & dt$l5_b1_bc13>=15),"l5_b1_bc13"]+.1,pch=20,col="coral")
# points(dt[which(dt$l5_b1_bc07/dt$l5_b1_bc13 > 10 & dt$l5_b1_bc13>=5),"l5_b1_bc07"]+.1,dt[which(dt$l5_b1_bc07/dt$l5_b1_bc13 > 10 & dt$l5_b1_bc13>=5),"l5_b1_bc13"]+.1,pch=20,col="coral")
# points(dt[which(dt$l5_b1_bc07 > 20 & dt$l5_b1_bc13==0),"l5_b1_bc07"]+.1,dt[which(dt$l5_b1_bc07 > 20 & dt$l5_b1_bc13==0),"l5_b1_bc13"]+.1,pch=20,col="coral")
AAseqs.l5 <- unique(c(AAseqs.l5,as.character(dt[c(which(dt$l5_b1_bc07/dt$l5_b1_bc13 > 3 & dt$l5_b1_bc13>=15),which(dt$l5_b1_bc07/dt$l5_b1_bc13 > 10 & dt$l5_b1_bc13>=5),which(dt$l5_b1_bc07 > 20 & dt$l5_b1_bc13==0)),1])))

# #lib5 bin2
# plot(dt$l5_b2_bc06+.1,dt$l5_b2_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=5);abline(v=5)
# points(dt[which(dt$l5_b2_bc06/dt$l5_b2_bc14 > 8 & dt$l5_b2_bc14 >= 5),"l5_b2_bc06"]+.1,dt[which(dt$l5_b2_bc06/dt$l5_b2_bc14 > 8 & dt$l5_b2_bc14 >= 5),"l5_b2_bc14"]+.1,pch=20,col="coral")
AAseqs.l5 <- unique(c(AAseqs.l5,as.character(dt[which(dt$l5_b2_bc06/dt$l5_b2_bc14 > 8 & dt$l5_b2_bc14 >= 5),1])))

# #l13 (ERE rep2) bin1
# plot(dt$l13_b1_bc03+.1,dt$l13_b1_bc09+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=25);abline(v=25)
# points(dt[which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 2.5 & dt$l13_b1_bc09>=100),"l13_b1_bc03"],dt[which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 2.5 & dt$l13_b1_bc09>=100),"l13_b1_bc09"],pch=20,col="coral")
# points(dt[which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 4 & dt$l13_b1_bc09>=20),"l13_b1_bc03"],dt[which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 4 & dt$l13_b1_bc09>=20),"l13_b1_bc09"],pch=20,col="coral")
# points(dt[which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 6 & dt$l13_b1_bc09>=5),"l13_b1_bc03"],dt[which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 6 & dt$l13_b1_bc09>=5),"l13_b1_bc09"],pch=20,col="coral")
# points(dt[which(dt$l13_b1_bc03 >= 50 & dt$l13_b1_bc09==2),"l13_b1_bc03"]+0.1,dt[which(dt$l13_b1_bc03 >= 50 & dt$l13_b1_bc09==2),"l13_b1_bc09"]+0.1,pch=20,col="coral")
AAseqs.l13 <- unique(as.character(dt[c(which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 2.5 & dt$l13_b1_bc09>=100),which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 4 & dt$l13_b1_bc09>=20),which(dt$l13_b1_bc03/dt$l13_b1_bc09 > 6 & dt$l13_b1_bc09>=5),which(dt$l13_b1_bc03 >= 50 & dt$l13_b1_bc09==2)),1]))

# #l13 (ERE rep2) bin2
# plot(dt$l13_b2_bc08+.1,dt$l13_b2_bc16+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=10);abline(v=10);abline(h=25);abline(v=25)
# points(dt[which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 2.5 & dt$l13_b2_bc16>=100),"l13_b2_bc08"],dt[which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 2.5 & dt$l13_b2_bc16>=100),"l13_b2_bc16"],pch=20,col="coral")
# points(dt[which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 3.5 & dt$l13_b2_bc16>=20),"l13_b2_bc08"],dt[which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 3.5 & dt$l13_b2_bc16>=20),"l13_b2_bc16"],pch=20,col="coral")
# points(dt[which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 8 & dt$l13_b2_bc16>=5),"l13_b2_bc08"],dt[which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 8 & dt$l13_b2_bc16>=5),"l13_b2_bc16"],pch=20,col="coral")
AAseqs.l13 <- unique(c(AAseqs.l13,as.character(dt[c(which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 2.5 & dt$l13_b2_bc16>=100),which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 3.5 & dt$l13_b2_bc16>=20),which(dt$l13_b2_bc08/dt$l13_b2_bc16 > 8 & dt$l13_b2_bc16>=5)),1])))


# #lib8 bin 1 (SRE rep2)
# plot(dt$l8_b1_bc05+.1,dt$l8_b1_bc12+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=20);abline(v=20)
# points(dt[which(dt$l8_b1_bc05/dt$l8_b1_bc12 > 6 & dt$l8_b1_bc12>=20),"l8_b1_bc05"]+.1,dt[which(dt$l8_b1_bc05/dt$l8_b1_bc12 > 6 & dt$l8_b1_bc12>=20),"l8_b1_bc12"]+.1,pch=20,col="coral")
# points(dt[which(dt$l8_b1_bc05/dt$l8_b1_bc12 > 10 & dt$l8_b1_bc12>=10),"l8_b1_bc05"]+.1,dt[which(dt$l8_b1_bc05/dt$l8_b1_bc12 > 10 & dt$l8_b1_bc12>=10),"l8_b1_bc12"]+.1,pch=20,col="coral")
AAseqs.l8 <- unique(as.character(dt[c(which(dt$l8_b1_bc05/dt$l8_b1_bc12 > 6 & dt$l8_b1_bc12>=20),which(dt$l8_b1_bc05/dt$l8_b1_bc12 > 10 & dt$l8_b1_bc12>=10)),1]))

# plot(dt$l8_b1_bc05+.1,dt$l8_b1_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=12);abline(v=12)
# points(dt[which(dt$l8_b1_bc05/dt$l8_b1_bc14 > 7 & dt$l8_b1_bc14>=12),"l8_b1_bc05"],dt[which(dt$l8_b1_bc05/dt$l8_b1_bc14 > 7 & dt$l8_b1_bc14>=12),"l8_b1_bc14"],pch=20,col="coral")
AAseqs.l8 <- unique(c(AAseqs.l8,as.character(dt[which(dt$l8_b1_bc05/dt$l8_b1_bc14 > 7 & dt$l8_b1_bc14>=12),1])))

# plot(dt$l8_b1_bc12+.1,dt$l8_b1_bc14+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=12);abline(v=12)
# points(dt[which(dt$l8_b1_bc12/dt$l8_b1_bc14 > 6 & dt$l8_b1_bc14>=12),"l8_b1_bc12"]+.1,dt[which(dt$l8_b1_bc12/dt$l8_b1_bc14 > 6 & dt$l8_b1_bc14>=12),"l8_b1_bc14"]+.1,pch=20,col="coral")
AAseqs.l8 <- unique(c(AAseqs.l8,as.character(dt[which(dt$l8_b1_bc12/dt$l8_b1_bc14 > 6 & dt$l8_b1_bc14>=12),1])))

# #l8 b2
# plot(dt$l8_b2_bc08+.1,dt$l8_b2_bc11+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=15);abline(v=15)
# points(dt[which(dt$l8_b2_bc08/dt$l8_b2_bc11 > 4 & dt$l8_b2_bc11>=15),"l8_b2_bc08"]+.1,dt[which(dt$l8_b2_bc08/dt$l8_b2_bc11 > 4 & dt$l8_b2_bc11>=15),"l8_b2_bc11"]+.1,pch=20,col="coral")
# points(dt[which(dt$l8_b2_bc08/dt$l8_b2_bc11 > 3 & dt$l8_b2_bc11>=50),"l8_b2_bc08"]+.1,dt[which(dt$l8_b2_bc08/dt$l8_b2_bc11 > 3 & dt$l8_b2_bc11>=50),"l8_b2_bc11"]+.1,pch=20,col="coral")
AAseqs.l8 <- unique(c(AAseqs.l8,as.character(dt[c(which(dt$l8_b2_bc08/dt$l8_b2_bc11 > 3 & dt$l8_b2_bc11>=50),which(dt$l8_b2_bc08/dt$l8_b2_bc11 > 4 & dt$l8_b2_bc11>=15)),1])))

# plot(dt$l8_b2_bc08+.1,dt$l8_b2_bc15+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=15);abline(v=15)
# points(dt[which(dt$l8_b2_bc08/dt$l8_b2_bc15 > 4 & dt$l8_b2_bc15>=25),"l8_b2_bc08"]+.1,dt[which(dt$l8_b2_bc08/dt$l8_b2_bc15 > 4 & dt$l8_b2_bc15>=25),"l8_b2_bc15"]+.1,pch=20,col="coral")
# points(dt[which(dt$l8_b2_bc08/dt$l8_b2_bc15 > 6 & dt$l8_b2_bc15>=12),"l8_b2_bc08"]+.1,dt[which(dt$l8_b2_bc08/dt$l8_b2_bc15 > 6 & dt$l8_b2_bc15>=12),"l8_b2_bc15"]+.1,pch=20,col="coral")
AAseqs.l8 <- unique(c(AAseqs.l8,as.character(dt[c(which(dt$l8_b2_bc08/dt$l8_b2_bc15 > 4 & dt$l8_b2_bc15>=25),which(dt$l8_b2_bc08/dt$l8_b2_bc15 > 6 & dt$l8_b2_bc15>=12)),1])))

# plot(dt$l8_b2_bc11+.1,dt$l8_b2_bc15+.1,pch=20,cex=0.5,log="xy",col="#00000066")
# abline(h=15);abline(v=15)
# points(dt[which(dt$l8_b2_bc11/dt$l8_b2_bc15 > 4.5 & dt$l8_b2_bc15>15),"l8_b2_bc11"]+.1,dt[which(dt$l8_b2_bc11/dt$l8_b2_bc15 > 4.5 & dt$l8_b2_bc15>15),"l8_b2_bc15"]+.1,pch=20,col="coral")
AAseqs.l8 <- unique(c(AAseqs.l8,as.character(dt[which(dt$l8_b2_bc11/dt$l8_b2_bc15 > 4.5 & dt$l8_b2_bc15>15),1])))

#remove reads for potentially biased variants from each replicate
dt[dt$AAseq %in% AAseqs.l6,c("l6_b1_bc01","l6_b1_bc06","l6_b1_bc10","l6_b2_bc03","l6_b2_bc09","l6_b2_bc16","l6_b3_bc07","l6_b4_bc04")] <- 0
dt[dt$AAseq %in% AAseqs.l5,c("l5_b1_bc01","l5_b1_bc07","l5_b1_bc13","l5_b2_bc06","l5_b2_bc14","l5_b3_bc02","l5_b4_bc04")] <- 0
dt[dt$AAseq %in% AAseqs.l13,c("l13_b1_bc03","l13_b1_bc09","l13_b2_bc08","l13_b2_bc16","l13_b3_bc10","l13_b4_bc04")] <- 0
dt[dt$AAseq %in% AAseqs.l8,c("l8_b1_bc05","l8_b1_bc12","l8_b1_bc14","l8_b2_bc08","l8_b2_bc11","l8_b2_bc15","l8_b3_bc02","l8_b4_bc13")] <- 0
rm(AAseqs.l6);rm(AAseqs.l5);rm(AAseqs.l13);rm(AAseqs.l8)

#colony counts, estimated from serial dilution and plating after cell sorting
cfu.l6.b1 <- 20212500; cfu.l6.b2 <- 33810000; cfu.l6.b3 <- 3220988; cfu.l6.b4 <- 474340
cfu.l5.b1 <- 1.575e7;cfu.l5.b2 <- 1.313e7;cfu.l5.b3 <- 1.785e6;cfu.l5.b4 <- 1.86e5
cfu.l13.b1 <- 16883700; cfu.l13.b2 <- 15939000; cfu.l13.b3 <- 2891125; cfu.l13.b4 <- 108000
cfu.l8.b1 <- 20246000; cfu.l8.b2 <- 31683000; cfu.l8.b3 <- 4909250; cfu.l8.b4 <- 410872

#pool together individual barcodes representing identical bins
dt$l6.b1 <- rowSums(dt[,c("l6_b1_bc01","l6_b1_bc06","l6_b1_bc10")]);dt$l6.b2 <- rowSums(dt[,c("l6_b2_bc03","l6_b2_bc09","l6_b2_bc16")]);dt$l6.b3 <- dt$l6_b3_bc07; dt$l6.b4 <- dt$l6_b4_bc04
dt$l5.b1 <- rowSums(dt[,c("l5_b1_bc01","l5_b1_bc07","l5_b1_bc13")]);dt$l5.b2 <- rowSums(dt[,c("l5_b2_bc06","l5_b2_bc14")]);dt$l5.b3 <- dt$l5_b3_bc02; dt$l5.b4 <- dt$l5_b4_bc04
dt$l13.b1 <- rowSums(dt[,c("l13_b1_bc03","l13_b1_bc09")]);dt$l13.b2 <- rowSums(dt[,c("l13_b2_bc08","l13_b2_bc16")]);dt$l13.b3 <- dt$l13_b3_bc10;dt$l13.b4 <- dt$l13_b4_bc04
dt$l8.b1 <- rowSums(dt[,c("l8_b1_bc05","l8_b1_bc12","l8_b1_bc14")]);dt$l8.b2 <- rowSums(dt[,c("l8_b2_bc08","l8_b2_bc11","l8_b2_bc15")]);dt$l8.b3 <- dt$l8_b3_bc02;dt$l8.b4 <- dt$l8_b4_bc13

#normalize # of reads by number of cfus in bin
dt$ERE.rep1.b1 <- dt$l6.b1/sum(dt$l6.b1)*cfu.l6.b1;dt$ERE.rep1.b2 <- dt$l6.b2/sum(dt$l6.b2)*cfu.l6.b2;dt$ERE.rep1.b3 <- dt$l6.b3/sum(dt$l6.b3)*cfu.l6.b3;dt$ERE.rep1.b4 <- dt$l6.b4/sum(dt$l6.b4)*cfu.l6.b4
dt$SRE.rep1.b1 <- dt$l5.b1/sum(dt$l5.b1)*cfu.l5.b1;dt$SRE.rep1.b2 <- dt$l5.b2/sum(dt$l5.b2)*cfu.l5.b2;dt$SRE.rep1.b3 <- dt$l5.b3/sum(dt$l5.b3)*cfu.l5.b3;dt$SRE.rep1.b4 <- dt$l5.b4/sum(dt$l5.b4)*cfu.l5.b4
dt$ERE.rep2.b1 <- dt$l13.b1/sum(dt$l13.b1)*cfu.l13.b1;dt$ERE.rep2.b2 <- dt$l13.b2/sum(dt$l13.b2)*cfu.l13.b2;dt$ERE.rep2.b3 <- dt$l13.b3/sum(dt$l13.b3)*cfu.l13.b3;dt$ERE.rep2.b4 <- dt$l13.b4/sum(dt$l13.b4)*cfu.l13.b4
dt$SRE.rep2.b1 <- dt$l8.b1/sum(dt$l8.b1)*cfu.l8.b1;dt$SRE.rep2.b2 <- dt$l8.b2/sum(dt$l8.b2)*cfu.l8.b2;dt$SRE.rep2.b3 <- dt$l8.b3/sum(dt$l8.b3)*cfu.l8.b3;dt$SRE.rep2.b4 <- dt$l8.b4/sum(dt$l8.b4)*cfu.l8.b4

#sum total number of weighted reads per variant in each rep
dt$ERE.rep1.cfu <- rowSums(dt[,c("ERE.rep1.b1","ERE.rep1.b2","ERE.rep1.b3","ERE.rep1.b4")])
dt$SRE.rep1.cfu <- rowSums(dt[,c("SRE.rep1.b1","SRE.rep1.b2","SRE.rep1.b3","SRE.rep1.b4")])
dt$ERE.rep2.cfu <- rowSums(dt[,c("ERE.rep2.b1","ERE.rep2.b2","ERE.rep2.b3","ERE.rep2.b4")])
dt$SRE.rep2.cfu <- rowSums(dt[,c("SRE.rep2.b1","SRE.rep2.b2","SRE.rep2.b3","SRE.rep2.b4")])
dt$ERE.pooled.cfu <- rowSums(dt[,c("ERE.rep1.cfu","ERE.rep2.cfu")])
dt$SRE.pooled.cfu <- rowSums(dt[,c("SRE.rep1.cfu","SRE.rep2.cfu")])

#convert to data.table for faster computation and indexing
dt <- data.table(dt[,c("AAseq","ERE.rep1.b1","ERE.rep1.b2","ERE.rep1.b3","ERE.rep1.b4","SRE.rep1.b1","SRE.rep1.b2","SRE.rep1.b3","SRE.rep1.b4","ERE.rep2.b1","ERE.rep2.b2","ERE.rep2.b3","ERE.rep2.b4","SRE.rep2.b1","SRE.rep2.b2","SRE.rep2.b3","SRE.rep2.b4","ERE.rep1.cfu","ERE.rep2.cfu","ERE.pooled.cfu","SRE.rep1.cfu","SRE.rep2.cfu","SRE.pooled.cfu")]);setkey(dt,AAseq)

#rescale fluorescence scales btwn reps, so I can combine observations to calculate pooled meanF, and also compare SR1 bkgrd data
#rescale to ERE rep1, SRE rep2 -- use SRE rep2 b/c I ran more controls with this expt, so I have more points to compare to with the SR1 data in the accompanying script
ERE.rep1.ctrl <- read.csv(file="1_calc-meanF_11P_in/FITC_lib6-controls.csv",header=T);ERE.rep2.ctrl <- read.csv(file="./1_calc-meanF_11P_in/FITC_lib13-controls.csv",header=T)
SRE.rep1.ctrl <- read.csv(file="1_calc-meanF_11P_in/FITC_lib5-controls.csv",header=T);SRE.rep2.ctrl <- read.csv(file="./1_calc-meanF_11P_in/FITC_lib8-controls.csv",header=T)

ERE.ctrls <- data.frame(geno=c("lib","null","egka","EGKA","GGKA","GSKV","GGKV"))
for(i in 1:nrow(ERE.ctrls)){
  ERE.ctrls$rep1.meanF[i] <- mean(log(ERE.rep1.ctrl[ERE.rep1.ctrl[,as.character(ERE.ctrls$geno[i])] > 0,as.character(ERE.ctrls$geno[i])]),na.rm=T)
  ERE.ctrls$rep2.meanF[i] <- mean(log(ERE.rep2.ctrl[ERE.rep2.ctrl[,as.character(ERE.ctrls$geno[i])] > 0,as.character(ERE.ctrls$geno[i])]),na.rm=T)
}
plot(ERE.ctrls$rep1.meanF,ERE.ctrls$rep2.meanF)
lm.ERE.ctrl <- lm(ERE.ctrls$rep2.meanF~ERE.ctrls$rep1.meanF);abline(lm.ERE.ctrl);summary(lm.ERE.ctrl)
#rep2 = 1.00035*rep1 + 0.01014

SRE.ctrls <- data.frame(geno=c("lib","null","GSKV","GGKA","EGKA"))
for(i in 1:nrow(SRE.ctrls)){
  SRE.ctrls$rep1.meanF[i] <- mean(log(SRE.rep1.ctrl[SRE.rep1.ctrl[,as.character(SRE.ctrls$geno[i])] > 0,as.character(SRE.ctrls$geno[i])]),na.rm=T)
  SRE.ctrls$rep2.meanF[i] <- mean(log(SRE.rep2.ctrl[SRE.rep2.ctrl[,as.character(SRE.ctrls$geno[i])] > 0,as.character(SRE.ctrls$geno[i])]),na.rm=T)
}
plot(SRE.ctrls$rep2.meanF,SRE.ctrls$rep1.meanF)
lm.SRE.ctrl <- lm(SRE.ctrls$rep1.meanF~SRE.ctrls$rep2.meanF);abline(lm.SRE.ctrl);summary(lm.SRE.ctrl)
#rep1 = 0.95213*rep2 + 0.46994

#give fluor boundaries of sort bins, with rescaled rep2 values
min.l6.b1 <- log(1); min.l6.b2 <- log(129.5); min.l6.b3 <- log(614.5); min.l6.b4 <- log(1284.5); max.l6.b4 <- log(262144)
min.l5.b1 <- log(1); min.l5.b2 <- (log(178.5)-0.46994)/0.95213; min.l5.b3 <- (log(496.5)-0.46994)/0.95213; min.l5.b4 <- (log(2183.5)-0.46994)/0.95213; max.l5.b4 <- log(262144)
min.l13.b1<- log(1); min.l13.b2<- (log(137.5)-0.01014)/1.00035; min.l13.b3<- (log(329.5)-0.01014)/1.00035; min.l13.b4<- (log(938.5)-0.01014)/1.00035; max.l13.b4 <- log(262144)
min.l8.b1 <- log(1); min.l8.b2 <- log(140.5); min.l8.b3 <- log(472.5); min.l8.b4 <- log(1875.5); max.l8.b4 <- log(262144)


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
dt[,c("ERE.rep1.meanF","ERE.rep1.sdF") := fit.bin.logistic(ERE.rep1.b1,ERE.rep1.b2,ERE.rep1.b3,ERE.rep1.b4,min.l6.b1,min.l6.b2,min.l6.b3,min.l6.b4,max.l6.b4),by=AAseq]
dt[,c("ERE.rep2.meanF","ERE.rep2.sdF") := fit.bin.logistic(ERE.rep2.b1,ERE.rep2.b2,ERE.rep2.b3,ERE.rep2.b4,min.l13.b1,min.l13.b2,min.l13.b3,min.l13.b4,max.l13.b4),by=AAseq]
dt[,c("ERE.pooled.meanF","ERE.pooled.sdF") := tryCatch(fit.bin.logistic.pooled(ERE.rep1.b1,ERE.rep1.b2,ERE.rep1.b3,ERE.rep1.b4,ERE.rep2.b1,ERE.rep2.b2,ERE.rep2.b3,ERE.rep2.b4,min.l6.b1,min.l6.b2,min.l6.b3,min.l6.b4,max.l6.b4,min.l13.b1,min.l13.b2,min.l13.b3,min.l13.b4,max.l13.b4),error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=AAseq]
dt[,c("SRE.rep1.meanF","SRE.rep1.sdF") := fit.bin.logistic(SRE.rep1.b1,SRE.rep1.b2,SRE.rep1.b3,SRE.rep1.b4,min.l5.b1,min.l5.b2,min.l5.b3,min.l5.b4,max.l5.b4),by=AAseq]
dt[,c("SRE.rep2.meanF","SRE.rep2.sdF") := fit.bin.logistic(SRE.rep2.b1,SRE.rep2.b2,SRE.rep2.b3,SRE.rep2.b4,min.l8.b1,min.l8.b2,min.l8.b3,min.l8.b4,max.l8.b4),by=AAseq]
dt[,c("SRE.pooled.meanF","SRE.pooled.sdF") := tryCatch(fit.bin.logistic.pooled(SRE.rep1.b1,SRE.rep1.b2,SRE.rep1.b3,SRE.rep1.b4,SRE.rep2.b1,SRE.rep2.b2,SRE.rep2.b3,SRE.rep2.b4,min.l5.b1,min.l5.b2,min.l5.b3,min.l5.b4,max.l5.b4,min.l8.b1,min.l8.b2,min.l8.b3,min.l8.b4,max.l8.b4),error=function(e){return(list(as.numeric(NA),as.numeric(NA)))}),by=AAseq]

#what variants are undetermined in pooled rep, but determined in either or both indiivdual reps?
#19 variants determined in ERE rep1, undetermined pooled ERE --> all have >1 read in l6.b1 and l6.b2, >1 read in l13.b1, <1 in all other bins --> assign pooled.meanF to be same as ERE.rep1.meanF
dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF),c(1:5,10:13,18:20,24,26,28),with=FALSE]
dt[!is.na(ERE.rep1.meanF) & is.na(ERE.pooled.meanF),ERE.pooled.meanF:=ERE.rep1.meanF]
#32 variants determined in ERE rep2 (lib13), undetermined pooled ERE --> all have >1 read in l13.b1 and l13.b2, and >1 read in l6.b2, with <1 in all other bins --> assign pooled.meanF to be same as ERE.rep2.meanF
dt[!is.na(ERE.rep2.meanF) & is.na(ERE.pooled.meanF),c(1:5,10:13,18:20,24,26,28),with=FALSE]
dt[!is.na(ERE.rep2.meanF) & is.na(ERE.pooled.meanF),ERE.pooled.meanF:=ERE.rep2.meanF]

#1615 variants determined in SRE rep1 (lib5), undetermined in pooled SRE --> all have >1 read in l5.b1 and l5.b2 and l8.b2, <1 read in l5.b3, l5.b4, l8.b3, l8.b4, some of which have >1 read in l8.b1
#--> for 0 reads in l8.b1 (undetermined rep2), assign pooled.meanF to be rep1.meanF. For 130 variants with bin1+2 reads in both replicates, assign pooled meanF to be average meanF of two replicates, both of which are individually determined
dt[!is.na(SRE.rep1.meanF) & is.na(SRE.pooled.meanF) & SRE.rep2.b1<0.5,c(1,6:9,14:17,21:23,30,32,34),with=FALSE]
dt[!is.na(SRE.rep1.meanF) & is.na(SRE.pooled.meanF) & SRE.rep2.b1<0.5,SRE.pooled.meanF:=SRE.rep1.meanF]
dt[!is.na(SRE.rep1.meanF) & is.na(SRE.pooled.meanF) & SRE.rep2.b1>0.5 & SRE.rep2.b2>0.5,c(1,6:9,14:17,21:23,30,32,34),with=FALSE]
dt[!is.na(SRE.rep1.meanF) & is.na(SRE.pooled.meanF) & SRE.rep2.b1>0.5 & SRE.rep2.b2>0.5,SRE.pooled.meanF := (SRE.rep1.meanF+SRE.rep2.meanF)/2,by=AAseq]
#2242 variants determined in SRE rep2 (lib8), undetermined in pooled SRE --> all have >1 read in l8.b1, l8.b2, l5.b1 <1 read in l8.b3, l8.b4, l5.b2 through b4 --> assign pooled meanF to be rep2 meanF
dt[!is.na(SRE.rep2.meanF) & is.na(SRE.pooled.meanF),c(1,6:9,14:17,21:23,30,32,34),with=FALSE];dt[!is.na(SRE.rep2.meanF) & is.na(SRE.pooled.meanF),SRE.pooled.meanF:=SRE.rep2.meanF]

write.table(dt, file="1_calc-meanF_11P_out/dt_output_11P.csv",col.names=T,row.names=F,quote=F,sep=",")
