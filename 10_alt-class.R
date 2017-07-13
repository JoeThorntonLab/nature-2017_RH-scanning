#16 Jan 2017
#TNS
#script to output alignments for weblogo creation; logistic regressions for determinants of binding specificity

library(data.table)
library(rgexf)
library(igraph)
library(Matrix)
library(glmnetcr)
library(fitdistrplus)
library(Hmisc)

setwd("path/to/source/directory") #setwd to source file directory

#define functions, background info...
AAs <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
AA.nt.neighbors <- read.csv(file="7_assess-connectivity_in/AA_single-mutant-codon-neighbors.csv",header=T,row.names=1,quote="")

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

get.3mut.paths.SR1 <- function(AAseq,ERE.sp=c(as.character(ERE.specific.SR1$AAseq),"EGKA"),prom=as.character(promiscuous.SR1$AAseq),SRE.sp=as.character(SRE.specific.SR1$AAseq)){
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

###########################################################################################
#scheme 1: use only model-predicted classifications for all variants, even if >15 cfu

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
#define new specificity class, given entirely by ERE and SRE.prediction classes
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.prediction=="strong" & SRE.prediction=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.prediction=="null" & SRE.prediction=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.prediction=="strong" & SRE.prediction=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.prediction=="weak" & SRE.prediction=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.prediction=="strong" & SRE.prediction=="weak",specificity:="promiscuous"]

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)
#define new specificity class, given entirely by ERE and SRE.prediction classes
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.prediction=="strong" & SRE.prediction=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.prediction=="null" & SRE.prediction=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.prediction=="strong" & SRE.prediction=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.prediction=="weak" & SRE.prediction=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.prediction=="strong" & SRE.prediction=="weak",specificity:="promiscuous"]

#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme1/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme1/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme1/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme1/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme1/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme1/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme1/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme1/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme1/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme1/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-5,5*round(max(ERE.specific.11P$num.SRE.outcomes)/5+5),5),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-1,max(SRE.specific.11P$accessed.by),1),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme1/all.pos.11P.Rda")
load("10_alt-class_out/scheme1/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme1/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme1/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme1/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme1/all.pos.SR1.Rda")
load("10_alt-class_out/scheme1/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme1/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme1/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme1/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme1/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])

###########################################################################################
#scheme 2: I expect my text and scheme 1 classifications are more prone to false negatives than false positives. First attempt
#at a method that decrease false negatives (at expense of false positives): for ERE and SRE categories, take highest category
#of either expt OR prediction

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
#define new specificity class, given by best-of prediction and experimental classes
dt.11P.coding[,ERE.full:="null"]
dt.11P.coding[ERE.pooled.class=="weak" | ERE.prediction=="weak", ERE.full:="weak"]
dt.11P.coding[ERE.pooled.class=="strong" | ERE.prediction=="strong", ERE.full:="strong"]
dt.11P.coding[,SRE.full:="null"]
dt.11P.coding[SRE.pooled.class=="weak" | SRE.prediction=="weak", SRE.full:="weak"]
dt.11P.coding[SRE.pooled.class=="strong" | SRE.prediction=="strong", SRE.full:="strong"]
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full=="strong" & SRE.full=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full=="null" & SRE.full=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full=="strong" & SRE.full=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full=="weak" & SRE.full=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full=="strong" & SRE.full=="weak",specificity:="promiscuous"]

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)
#define new specificity class, given by best-of prediction and experimental classes
names(dt.SR1.coding)
dt.SR1.coding[,ERE.full:="null"]
dt.SR1.coding[ERE.pooled.class=="weak" | ERE.prediction=="weak", ERE.full:="weak"]
dt.SR1.coding[ERE.pooled.class=="strong" | ERE.prediction=="strong", ERE.full:="strong"]
dt.SR1.coding[,SRE.full:="null"]
dt.SR1.coding[SRE.pooled.class=="weak" | SRE.prediction=="weak", SRE.full:="weak"]
dt.SR1.coding[SRE.pooled.class=="strong" | SRE.prediction=="strong", SRE.full:="strong"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full=="strong" & SRE.full=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full=="null" & SRE.full=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full=="strong" & SRE.full=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full=="weak" & SRE.full=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full=="strong" & SRE.full=="weak",specificity:="promiscuous"]

#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme2/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme2/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme2/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme2/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme2/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme2/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme2/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme2/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme2/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme2/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-5,5*round(max(ERE.specific.11P$num.SRE.outcomes)/5+5),5),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-1,max(SRE.specific.11P$accessed.by),1),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme2/all.pos.11P.Rda")
load("10_alt-class_out/scheme2/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme2/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme2/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme2/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme2/all.pos.SR1.Rda")
load("10_alt-class_out/scheme2/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme2/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme2/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme2/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme2/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])

###########################################################################################
#scheme3:
#try most extreme reduction of stringency of calling "strongs": anything with measureable activity above background (weak+strong)

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]

#define new specificity class, where any activity (weak+strong)is functional
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class%in%c("weak","strong") & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class%in%c("weak","strong"),specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class%in%c("weak","strong") & SRE.full.class%in%c("weak","strong"),specificity:="promiscuous"]

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)
#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.SR1.coding[,ERE.full.class := ERE.pooled.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]

#define new specificity class, where any activity (weak+strong)is functional
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class%in%c("weak","strong") & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class%in%c("weak","strong"),specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class%in%c("weak","strong") & SRE.full.class%in%c("weak","strong"),specificity:="promiscuous"]

#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme3/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme3/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme3/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme3/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme3/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme3/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme3/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme3/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme3/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme3/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5),5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme3/all.pos.11P.Rda")
load("10_alt-class_out/scheme3/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme3/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme3/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme3/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme3/all.pos.SR1.Rda")
load("10_alt-class_out/scheme3/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme3/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme3/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme3/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme3/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])



###########################################################################################
#scheme4:
#instead of determining "strong" as "significantly" better than the 80th percentile of evolutionary wildtype, define as simply greater than this threshold
#(that is, reduce problem of things close to this threshold being above it but not significantly so)
#I will have to re-do the glmnet.cr model predictions for this one to fill in missing genos after new classification

#get 80%ile cutoffs for ERE, SRE strong: this is EGKA in the SR1 library, GSKV in the 11P library
dt.11P <- read.table(file="1_calc-meanF_11P_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

dt.SR1 <- read.table(file="2_calc-meanF_SR1_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

SRE.pooled.stop.meanF.11P <- mean(dt.11P.stop$SRE.pooled.meanF,na.rm=T)
ERE.pooled.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.pooled.meanF,na.rm=T)

ERE.rep1.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l10_egka_ctrl_distribution.csv",header=F)$V1)+0.59618)/1.01805; ERE.rep1.egka.SR1 <- ERE.rep1.egka.SR1[is.finite(ERE.rep1.egka.SR1)]
ERE.rep2.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l12_egka_ctrl_distribution.csv",header=F)$V1[1:length(ERE.rep1.egka.SR1)])-0.5260)/0.9633
set.seed(104)
ERE.pooled.egka.SR1 <- sample(c(ERE.rep1.egka.SR1,ERE.rep2.egka.SR1))
ERE.cutoff <- (mean(ERE.pooled.egka.SR1)-ERE.pooled.stop.meanF.SR1)*0.8+ERE.pooled.stop.meanF.SR1
SRE.cutoff <- (dt.11P.coding["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.11P)*0.8+SRE.pooled.stop.meanF.11P
rm(dt.SR1);rm(dt.SR1.coding);rm(dt.SR1.stop);rm(dt.11P);rm(dt.11P.coding);rm(dt.11P.stop);rm(ERE.rep1.egka.SR1);rm(ERE.rep2.egka.SR1);rm(ERE.pooled.egka.SR1);rm(ERE.pooled.stop.meanF.SR1);rm(SRE.pooled.stop.meanF.11P)

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
dt.11P.coding[,ERE.pooled.new.class := as.character(NA)]
dt.11P.coding[,SRE.pooled.new.class := as.character(NA)]
dt.11P.coding[ERE.pooled.class=="null",ERE.pooled.new.class:="null"]
dt.11P.coding[ERE.pooled.class=="strong",ERE.pooled.new.class:="strong"]
dt.11P.coding[ERE.pooled.class=="weak" & ERE.pooled.meanF > ERE.cutoff,ERE.pooled.new.class:="strong"]
dt.11P.coding[ERE.pooled.class=="weak" & ERE.pooled.meanF <= ERE.cutoff,ERE.pooled.new.class:="weak"]
dt.11P.coding[SRE.pooled.class=="null",SRE.pooled.new.class:="null"]
dt.11P.coding[SRE.pooled.class=="strong",SRE.pooled.new.class:="strong"]
dt.11P.coding[SRE.pooled.class=="weak" & SRE.pooled.meanF > SRE.cutoff,SRE.pooled.new.class:="strong"]
dt.11P.coding[SRE.pooled.class=="weak" & SRE.pooled.meanF <= SRE.cutoff,SRE.pooled.new.class:="weak"]

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)
dt.SR1.coding[,ERE.pooled.new.class := as.character(NA)]
dt.SR1.coding[,SRE.pooled.new.class := as.character(NA)]
dt.SR1.coding[ERE.pooled.class=="null",ERE.pooled.new.class:="null"]
dt.SR1.coding[ERE.pooled.class=="strong",ERE.pooled.new.class:="strong"]
dt.SR1.coding[ERE.pooled.class=="weak" & ERE.pooled.meanF > ERE.cutoff,ERE.pooled.new.class:="strong"]
dt.SR1.coding[ERE.pooled.class=="weak" & ERE.pooled.meanF <= ERE.cutoff,ERE.pooled.new.class:="weak"]
dt.SR1.coding[SRE.pooled.class=="null",SRE.pooled.new.class:="null"]
dt.SR1.coding[SRE.pooled.class=="strong",SRE.pooled.new.class:="strong"]
dt.SR1.coding[SRE.pooled.class=="weak" & SRE.pooled.meanF > SRE.cutoff,SRE.pooled.new.class:="strong"]
dt.SR1.coding[SRE.pooled.class=="weak" & SRE.pooled.meanF <= SRE.cutoff,SRE.pooled.new.class:="weak"]

dt.11P.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.11P.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]

#build ordinal logistic regression models for each of the four experiments with the new classifications to predict missing genos
#ERE, 11P
#infer model
dt.ERE.11P.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.pooled.new.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled, 
                                     y=dt.11P.coding$ERE.pooled.new.class[!is.na(dt.11P.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.11P.pooled, file="10_alt-class_out/scheme4/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="10_alt-class_out/scheme4/glmnet.cr.ERE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.11P <- which(abs(logit.ERE.11P.pooled$lambda-(1e-5))==min(abs(logit.ERE.11P.pooled$lambda-(1e-5))))

#predict missing
dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m,s=step.ERE.11P)
save(ERE.11P.predictions, file="10_alt-class_out/scheme4/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.new.prediction <- ERE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme4/dt.11P.coding.Rda")
rm(ERE.11P.predictions);rm(step.ERE.11P);rm(logit.ERE.11P.pooled);rm(dt.ERE.11P.pooled)

#SRE, 11P
#infer model
dt.SRE.11P.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$SRE.pooled.new.class)])
logit.SRE.11P.pooled <- glmnet.cr(x=dt.SRE.11P.pooled, 
                                     y=dt.11P.coding$SRE.pooled.new.class[!is.na(dt.11P.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.11P.pooled, file="10_alt-class_out/scheme4/dt.SRE.11P.pooled")
save(logit.SRE.11P.pooled, file="10_alt-class_out/scheme4/glmnet.cr.SRE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.11P <- which(abs(logit.SRE.11P.pooled$lambda-(1e-5))==min(abs(logit.SRE.11P.pooled$lambda-(1e-5))))

#predict missing
SRE.11P.predictions <- fitted(logit.SRE.11P.pooled, newx=dt.m,s=step.SRE.11P)
save(SRE.11P.predictions, file="10_alt-class_out/scheme4/SRE.11P.predictions.Rda")
dt.11P.coding$SRE.new.prediction <- SRE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme4/dt.11P.coding.Rda")
rm(SRE.11P.predictions);rm(step.SRE.11P);rm(logit.SRE.11P.pooled);rm(dt.SRE.11P.pooled)

#ERE, SR1
#infer model
dt.ERE.SR1.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.pooled.new.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled, 
                                     y=dt.SR1.coding$ERE.pooled.new.class[!is.na(dt.SR1.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.SR1.pooled, file="10_alt-class_out/scheme4/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="10_alt-class_out/scheme4/glmnet.cr.ERE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.SR1 <- which(abs(logit.ERE.SR1.pooled$lambda-(1e-5))==min(abs(logit.ERE.SR1.pooled$lambda-(1e-5))))

#predict missing
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m,s=step.ERE.SR1)
save(ERE.SR1.predictions, file="10_alt-class_out/scheme4/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.new.prediction <- ERE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme4/dt.SR1.coding.Rda")
rm(ERE.SR1.predictions);rm(step.ERE.SR1);rm(logit.ERE.SR1.pooled);rm(dt.ERE.SR1.pooled)

#SRE, SR1
#infer model
dt.SRE.SR1.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$SRE.pooled.new.class)])
logit.SRE.SR1.pooled <- glmnet.cr(x=dt.SRE.SR1.pooled, 
                                     y=dt.SR1.coding$SRE.pooled.new.class[!is.na(dt.SR1.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.SR1.pooled, file="10_alt-class_out/scheme4/dt.SRE.SR1.pooled")
save(logit.SRE.SR1.pooled, file="10_alt-class_out/scheme4/glmnet.cr.SRE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.SR1 <- which(abs(logit.SRE.SR1.pooled$lambda-(1e-5))==min(abs(logit.SRE.SR1.pooled$lambda-(1e-5))))

#predict missing
SRE.SR1.predictions <- fitted(logit.SRE.SR1.pooled, newx=dt.m,s=step.SRE.SR1)
save(SRE.SR1.predictions, file="10_alt-class_out/scheme4/SRE.SR1.predictions.Rda")
dt.SR1.coding$SRE.new.prediction <- SRE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme4/dt.SR1.coding.Rda")
rm(SRE.SR1.predictions);rm(step.SRE.SR1);rm(logit.SRE.SR1.pooled);rm(dt.SRE.SR1.pooled);rm(dt.m)


#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.new.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.new.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.new.prediction]
dt.SR1.coding[,ERE.full.class := ERE.pooled.new.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.new.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.new.prediction]


#define specificity class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]


#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme4/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme4/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme4/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme4/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme4/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme4/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme4/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme4/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme4/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme4/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5),5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme4/all.pos.11P.Rda")
load("10_alt-class_out/scheme4/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme4/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme4/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme4/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme4/all.pos.SR1.Rda")
load("10_alt-class_out/scheme4/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme4/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme4/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme4/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme4/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])


###########################################################################################
#scheme5:
#for p.strong, use larger margin of indistinguishability: make the cutoff 50%ile between wt and null

dt.11P <- read.table(file="1_calc-meanF_11P_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

dt.SR1 <- read.table(file="2_calc-meanF_SR1_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#define null distributions for hypothesis: inactive (pretty much only need for defining the positive bins, since I don't need to re-do these p.values)
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
#SRE pooled, 11P
SRE.pooled.null.11P <- data.frame(bin=1:25)
breaks <- cut2(dt.11P.stop[!is.na(SRE.pooled.meanF),SRE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
for(i in 1:nrow(SRE.pooled.null.11P)){
  SRE.pooled.null.11P$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.11P.stop[SRE.pooled.cfu>=SRE.pooled.null.11P$range.cfu[i][[1]][[1]] & SRE.pooled.cfu<SRE.pooled.null.11P$range.cfu[i][[1]][[2]],]
  SRE.pooled.null.11P$mean.meanF[i] <- mean(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.11P$sd.meanF[i] <- sd(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.11P$median.cfu[i] <- median(data$SRE.pooled.cfu,na.rm=T)
  SRE.pooled.null.11P$max.meanF[i] <- max(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.11P$list.means[i] <- list(data[!is.na(SRE.pooled.meanF),SRE.pooled.meanF])
}
#ERE pooled, SR1
ERE.pooled.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(ERE.pooled.meanF),ERE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
for(i in 1:nrow(ERE.pooled.null.SR1)){
  ERE.pooled.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[ERE.pooled.cfu>=ERE.pooled.null.SR1$range.cfu[i][[1]][[1]] & ERE.pooled.cfu<ERE.pooled.null.SR1$range.cfu[i][[1]][[2]],]
  ERE.pooled.null.SR1$mean.meanF[i] <- mean(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$sd.meanF[i] <- sd(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$median.cfu[i] <- median(data$ERE.pooled.cfu,na.rm=T)
  ERE.pooled.null.SR1$max.meanF[i] <- max(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$list.means[i] <- list(data[!is.na(ERE.pooled.meanF),ERE.pooled.meanF])
}
#SRE pooled, SR1
SRE.pooled.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(SRE.pooled.meanF),SRE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
for(i in 1:nrow(SRE.pooled.null.SR1)){
  SRE.pooled.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[SRE.pooled.cfu>=SRE.pooled.null.SR1$range.cfu[i][[1]][[1]] & SRE.pooled.cfu<SRE.pooled.null.SR1$range.cfu[i][[1]][[2]],]
  SRE.pooled.null.SR1$mean.meanF[i] <- mean(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.SR1$sd.meanF[i] <- sd(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.SR1$median.cfu[i] <- median(data$SRE.pooled.cfu,na.rm=T)
  SRE.pooled.null.SR1$max.meanF[i] <- max(data$SRE.pooled.meanF,na.rm=T)
  SRE.pooled.null.SR1$list.means[i] <- list(data[!is.na(SRE.pooled.meanF),SRE.pooled.meanF])
}

#determine meanF of stop variants
ERE.pooled.stop.meanF.11P <- mean(dt.11P.stop$ERE.pooled.meanF,na.rm=T)
SRE.pooled.stop.meanF.11P <- mean(dt.11P.stop$SRE.pooled.meanF,na.rm=T)
ERE.pooled.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.pooled.meanF,na.rm=T)
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

#shift controls over to 50th percentile between null and wildtype, to give most extreme instantiation of the null hypothesis that variant is weak that I want to test
ERE.pooled.egka.p50.11P <- ERE.pooled.egka.11P - (mean(ERE.pooled.egka.11P)-ERE.pooled.stop.meanF.11P)*0.5
SRE.pooled.GSKV.p50.11P <- SRE.pooled.GSKV.11P + (dt.11P["GSKV",SRE.pooled.meanF]-mean(SRE.pooled.GSKV.11P)) - (dt.11P["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.11P)*0.5
ERE.pooled.egka.p50.SR1 <- ERE.pooled.egka.SR1 - (mean(ERE.pooled.egka.SR1)-ERE.pooled.stop.meanF.SR1)*0.5
SRE.pooled.GSKV.p50.SR1 <- SRE.pooled.GSKV.SR1 + (dt.11P["GSKV",SRE.pooled.meanF]-mean(SRE.pooled.GSKV.SR1)) - (dt.11P["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.SR1)*0.5

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
ERE.pooled.egka.p50.11P.cens.left <- c(sapply(1:round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.11P)), function(x) if(ERE.pooled.egka.p50.11P[x]<min.l6.b2){return(min.l6.b1)}else if(ERE.pooled.egka.p50.11P[x]<min.l6.b3){return(min.l6.b2)}else if(ERE.pooled.egka.p50.11P[x]<min.l6.b4){return(min.l6.b3)}else{return(min.l6.b4)}),sapply(round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.11P)+1):length(ERE.pooled.egka.p50.11P), function(x) if(ERE.pooled.egka.p50.11P[x]<min.l13.b2){return(min.l13.b1)}else if(ERE.pooled.egka.p50.11P[x]<min.l13.b3){return(min.l13.b2)}else if(ERE.pooled.egka.p50.11P[x]<min.l13.b4){return(min.l13.b3)}else{return(min.l13.b4)}))
ERE.pooled.egka.p50.11P.cens.right <- c(sapply(1:round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.11P)), function(x) if(ERE.pooled.egka.p50.11P[x]<min.l6.b2){return(min.l6.b2)}else if(ERE.pooled.egka.p50.11P[x]<min.l6.b3){return(min.l6.b3)}else if(ERE.pooled.egka.p50.11P[x]<min.l6.b4){return(min.l6.b4)}else{return(max.l6.b4)}),sapply(round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.11P)+1):length(ERE.pooled.egka.p50.11P), function(x) if(ERE.pooled.egka.p50.11P[x]<min.l13.b2){return(min.l13.b2)}else if(ERE.pooled.egka.p50.11P[x]<min.l13.b3){return(min.l13.b3)}else if(ERE.pooled.egka.p50.11P[x]<min.l13.b4){return(min.l13.b4)}else{return(max.l13.b4)}))
ERE.pooled.egka.p50.11P.cens <- data.frame(left=ERE.pooled.egka.p50.11P.cens.left,right=ERE.pooled.egka.p50.11P.cens.right)

SRE.pooled.GSKV.p50.11P.cens.left <- c(sapply(1:round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.11P)), function(x) if(SRE.pooled.GSKV.p50.11P[x]<min.l5.b2){return(min.l5.b1)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l5.b3){return(min.l5.b2)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l5.b4){return(min.l5.b3)}else{return(min.l5.b4)}),sapply(round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.11P)+1):length(SRE.pooled.GSKV.p50.11P), function(x) if(SRE.pooled.GSKV.p50.11P[x]<min.l8.b2){return(min.l8.b1)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l8.b3){return(min.l8.b2)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l8.b4){return(min.l8.b3)}else{return(min.l8.b4)}))
SRE.pooled.GSKV.p50.11P.cens.right <- c(sapply(1:round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.11P)), function(x) if(SRE.pooled.GSKV.p50.11P[x]<min.l5.b2){return(min.l5.b2)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l5.b3){return(min.l5.b3)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l5.b4){return(min.l5.b4)}else{return(max.l5.b4)}),sapply(round(sum(dt.11P$SRE.rep1.cfu)/sum(dt.11P$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.11P)+1):length(SRE.pooled.GSKV.p50.11P), function(x) if(SRE.pooled.GSKV.p50.11P[x]<min.l8.b2){return(min.l8.b2)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l8.b3){return(min.l8.b3)}else if(SRE.pooled.GSKV.p50.11P[x]<min.l8.b4){return(min.l8.b4)}else{return(max.l8.b4)}))
SRE.pooled.GSKV.p50.11P.cens <- data.frame(left=SRE.pooled.GSKV.p50.11P.cens.left,right=SRE.pooled.GSKV.p50.11P.cens.right)

ERE.pooled.egka.p50.SR1.cens.left <- c(sapply(1:round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.SR1)), function(x) if(ERE.pooled.egka.p50.SR1[x]<min.l10.b2){return(min.l10.b1)}else if(ERE.pooled.egka.p50.SR1[x]<min.l10.b3){return(min.l10.b2)}else if(ERE.pooled.egka.p50.SR1[x]<min.l10.b4){return(min.l10.b3)}else{return(min.l10.b4)}),sapply(round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.SR1)+1):length(ERE.pooled.egka.p50.SR1), function(x) if(ERE.pooled.egka.p50.SR1[x]<min.l12.b2){return(min.l12.b1)}else if(ERE.pooled.egka.p50.SR1[x]<min.l12.b3){return(min.l12.b2)}else if(ERE.pooled.egka.p50.SR1[x]<min.l12.b4){return(min.l12.b3)}else{return(min.l12.b4)}))
ERE.pooled.egka.p50.SR1.cens.right <- c(sapply(1:round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.SR1)), function(x) if(ERE.pooled.egka.p50.SR1[x]<min.l10.b2){return(min.l10.b2)}else if(ERE.pooled.egka.p50.SR1[x]<min.l10.b3){return(min.l10.b3)}else if(ERE.pooled.egka.p50.SR1[x]<min.l10.b4){return(min.l10.b4)}else{return(max.l10.b4)}),sapply(round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p50.SR1)+1):length(ERE.pooled.egka.p50.SR1), function(x) if(ERE.pooled.egka.p50.SR1[x]<min.l12.b2){return(min.l12.b2)}else if(ERE.pooled.egka.p50.SR1[x]<min.l12.b3){return(min.l12.b3)}else if(ERE.pooled.egka.p50.SR1[x]<min.l12.b4){return(min.l12.b4)}else{return(max.l12.b4)}))
ERE.pooled.egka.p50.SR1.cens <- data.frame(left=ERE.pooled.egka.p50.SR1.cens.left,right=ERE.pooled.egka.p50.SR1.cens.right)

SRE.pooled.GSKV.p50.SR1.cens.left <- c(sapply(1:round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.SR1)), function(x) if(SRE.pooled.GSKV.p50.SR1[x]<min.l9.b2){return(min.l9.b1)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l9.b3){return(min.l9.b2)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l9.b4){return(min.l9.b3)}else{return(min.l9.b4)}),sapply(round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.SR1)+1):length(SRE.pooled.GSKV.p50.SR1), function(x) if(SRE.pooled.GSKV.p50.SR1[x]<min.l11.b2){return(min.l11.b1)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l11.b3){return(min.l11.b2)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l11.b4){return(min.l11.b3)}else{return(min.l11.b4)}))
SRE.pooled.GSKV.p50.SR1.cens.right <- c(sapply(1:round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.SR1)), function(x) if(SRE.pooled.GSKV.p50.SR1[x]<min.l9.b2){return(min.l9.b2)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l9.b3){return(min.l9.b3)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l9.b4){return(min.l9.b4)}else{return(max.l9.b4)}),sapply(round(sum(dt.SR1$SRE.rep1.cfu)/sum(dt.SR1$SRE.pooled.cfu)*length(SRE.pooled.GSKV.p50.SR1)+1):length(SRE.pooled.GSKV.p50.SR1), function(x) if(SRE.pooled.GSKV.p50.SR1[x]<min.l11.b2){return(min.l11.b2)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l11.b3){return(min.l11.b3)}else if(SRE.pooled.GSKV.p50.SR1[x]<min.l11.b4){return(min.l11.b4)}else{return(max.l11.b4)}))
SRE.pooled.GSKV.p50.SR1.cens <- data.frame(left=SRE.pooled.GSKV.p50.SR1.cens.left,right=SRE.pooled.GSKV.p50.SR1.cens.right)

boot.cens.dist <- function(data, cfu, rep){
  sample.means <- vector()
  while(length(sample.means)<rep){
    boot <- data[sample(1:nrow(data),cfu,replace=T),]
    sample.means <- c(sample.means,tryCatch( summary(fitdistcens(boot,"logis"))$estimate["location"],error=function(e){return(as.numeric(NA))}))
  }
  return(sample.means)
}

ERE.pooled.pos.11P <- data.frame(bin=1:25); ERE.pooled.pos.11P$range.cfu <- ERE.pooled.null.11P$range.cfu; ERE.pooled.pos.11P$median.cfu <- ERE.pooled.null.11P$median.cfu
for(i in 1:nrow(ERE.pooled.pos.11P)){
  ERE.pooled.pos.11P$list.means[i] <- list(boot.cens.dist(ERE.pooled.egka.p50.11P.cens, round(ERE.pooled.pos.11P$median.cfu[i]),10000))
  ERE.pooled.pos.11P$mean.meanF[i] <- mean(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.11P$sd.meanF[i] <- sd(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.11P$max.meanF[i] <- max(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(ERE.pooled.pos.11P,file="10_alt-class_out/scheme5/ERE.pooled.pos.11P.Rda")

SRE.pooled.pos.11P <- data.frame(bin=1:25); SRE.pooled.pos.11P$range.cfu <- SRE.pooled.null.11P$range.cfu; SRE.pooled.pos.11P$median.cfu <- SRE.pooled.null.11P$median.cfu
for(i in 1:nrow(SRE.pooled.pos.11P)){
  SRE.pooled.pos.11P$list.means[i] <- list(boot.cens.dist(SRE.pooled.GSKV.p50.11P.cens, round(SRE.pooled.pos.11P$median.cfu[i]),10000))
  SRE.pooled.pos.11P$mean.meanF[i] <- mean(SRE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.11P$sd.meanF[i] <- sd(SRE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.11P$max.meanF[i] <- max(SRE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(SRE.pooled.pos.11P,file="10_alt-class_out/scheme5/SRE.pooled.pos.11P.Rda")

ERE.pooled.pos.SR1 <- data.frame(bin=1:25); ERE.pooled.pos.SR1$range.cfu <- ERE.pooled.null.SR1$range.cfu; ERE.pooled.pos.SR1$median.cfu <- ERE.pooled.null.SR1$median.cfu
for(i in 1:nrow(ERE.pooled.pos.SR1)){
  ERE.pooled.pos.SR1$list.means[i] <- list(boot.cens.dist(ERE.pooled.egka.p50.SR1.cens, round(ERE.pooled.pos.SR1$median.cfu[i]),10000))
  ERE.pooled.pos.SR1$mean.meanF[i] <- mean(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.SR1$sd.meanF[i] <- sd(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.SR1$max.meanF[i] <- max(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(ERE.pooled.pos.SR1,file="10_alt-class_out/scheme5/ERE.pooled.pos.SR1.Rda")

SRE.pooled.pos.SR1 <- data.frame(bin=1:25); SRE.pooled.pos.SR1$range.cfu <- SRE.pooled.null.SR1$range.cfu; SRE.pooled.pos.SR1$median.cfu <- SRE.pooled.null.SR1$median.cfu
for(i in 1:nrow(SRE.pooled.pos.SR1)){
  SRE.pooled.pos.SR1$list.means[i] <- list(boot.cens.dist(SRE.pooled.GSKV.p50.SR1.cens, round(SRE.pooled.pos.SR1$median.cfu[i]),10000))
  SRE.pooled.pos.SR1$mean.meanF[i] <- mean(SRE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.SR1$sd.meanF[i] <- sd(SRE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  SRE.pooled.pos.SR1$max.meanF[i] <- max(SRE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(SRE.pooled.pos.SR1,file="10_alt-class_out/scheme5/SRE.pooled.pos.SR1.Rda")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt","ERE.pooled.pos.11P","SRE.pooled.pos.11P","ERE.pooled.pos.SR1","SRE.pooled.pos.SR1"))])


dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

#define function to get p value for test of non-inferiority versus evolutionary wildtype (with 50% margin of indistinguishability)
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

dt.11P.coding[,ERE.pooled.p.weak := get.p.weak(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.pos.11P),by=AAseq]
dt.11P.coding[,SRE.pooled.p.weak := get.p.weak(SRE.pooled.cfu,SRE.pooled.meanF,SRE.pooled.pos.11P),by=AAseq]
dt.SR1.coding[,ERE.pooled.p.weak := get.p.weak(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.pos.SR1),by=AAseq]
dt.SR1.coding[,SRE.pooled.p.weak := get.p.weak(SRE.pooled.cfu,SRE.pooled.meanF,SRE.pooled.pos.SR1),by=AAseq]

dt.11P.coding[,ERE.pooled.new.class:=ERE.pooled.class]
dt.11P.coding[ERE.pooled.class %in% c("weak","strong"),ERE.pooled.new.class:="weak"] #revert all strongs to weak
dt.11P.coding[,SRE.pooled.new.class:=SRE.pooled.class]
dt.11P.coding[SRE.pooled.class %in% c("weak","strong"),SRE.pooled.new.class:="weak"] #revert all strongs to weak
dt.SR1.coding[,ERE.pooled.new.class:=ERE.pooled.class]
dt.SR1.coding[ERE.pooled.class %in% c("weak","strong"),ERE.pooled.new.class:="weak"] #revert all strongs to weak
dt.SR1.coding[,SRE.pooled.new.class:=SRE.pooled.class]
dt.SR1.coding[SRE.pooled.class %in% c("weak","strong"),SRE.pooled.new.class:="weak"] #revert all strongs to weak


#determine p-value cutoff for hypothesis: weak (alt hypothesis: strong binder) that controls FDR to 5%
#ERE pooled, 11P
p.values.coding <- sort(c(dt.11P.coding[ERE.pooled.new.class=="weak",ERE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i
p.weak.ERE.pooled.cutoff.11P <- p.values.coding[i] 

#SRE pooled, 11P
p.values.coding <- sort(c(dt.11P.coding[SRE.pooled.new.class=="weak",SRE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i
p.weak.SRE.pooled.cutoff.11P <- p.values.coding[i] 

#ERE pooled, SR1
p.values.coding <- sort(c(dt.SR1.coding[ERE.pooled.new.class=="weak",ERE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i
p.weak.ERE.pooled.cutoff.SR1 <- p.values.coding[i] 

#SRE pooled, SR1
p.values.coding <- sort(c(dt.SR1.coding[SRE.pooled.new.class=="weak",SRE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i
p.weak.SRE.pooled.cutoff.SR1 <- p.values.coding[i] 

#assign as strong if p.value below p.value.cutoff
dt.11P.coding[ERE.pooled.p.weak<p.weak.ERE.pooled.cutoff.11P & ERE.pooled.new.class=="weak",ERE.pooled.new.class:="strong"]
dt.11P.coding[SRE.pooled.p.weak<p.weak.SRE.pooled.cutoff.11P & SRE.pooled.new.class=="weak",SRE.pooled.new.class:="strong"]
dt.SR1.coding[ERE.pooled.p.weak<p.weak.ERE.pooled.cutoff.SR1 & ERE.pooled.new.class=="weak",ERE.pooled.new.class:="strong"]
dt.SR1.coding[SRE.pooled.p.weak<p.weak.SRE.pooled.cutoff.SR1 & SRE.pooled.new.class=="weak",SRE.pooled.new.class:="strong"]

dt.11P.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.11P.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]

#build ordinal logistic regression models for each of the four experiments with the new classifications to predict missing genos
#ERE, 11P
#infer model
dt.ERE.11P.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.pooled.new.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled, 
                                     y=dt.11P.coding$ERE.pooled.new.class[!is.na(dt.11P.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.11P.pooled, file="10_alt-class_out/scheme5/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="10_alt-class_out/scheme5/glmnet.cr.ERE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.11P <- which(abs(logit.ERE.11P.pooled$lambda-(1e-5))==min(abs(logit.ERE.11P.pooled$lambda-(1e-5))))

#predict missing
dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m,s=step.ERE.11P)
save(ERE.11P.predictions, file="10_alt-class_out/scheme5/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.new.prediction <- ERE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme5/dt.11P.coding.Rda")
rm(ERE.11P.predictions);rm(step.ERE.11P);rm(logit.ERE.11P.pooled);rm(dt.ERE.11P.pooled)

#SRE, 11P
#infer model
dt.SRE.11P.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$SRE.pooled.new.class)])
logit.SRE.11P.pooled <- glmnet.cr(x=dt.SRE.11P.pooled, 
                                     y=dt.11P.coding$SRE.pooled.new.class[!is.na(dt.11P.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.11P.pooled, file="10_alt-class_out/scheme5/dt.SRE.11P.pooled")
save(logit.SRE.11P.pooled, file="10_alt-class_out/scheme5/glmnet.cr.SRE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.11P <- which(abs(logit.SRE.11P.pooled$lambda-(1e-5))==min(abs(logit.SRE.11P.pooled$lambda-(1e-5))))

#predict missing
SRE.11P.predictions <- fitted(logit.SRE.11P.pooled, newx=dt.m,s=step.SRE.11P)
save(SRE.11P.predictions, file="10_alt-class_out/scheme5/SRE.11P.predictions.Rda")
dt.11P.coding$SRE.new.prediction <- SRE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme5/dt.11P.coding.Rda")
rm(SRE.11P.predictions);rm(step.SRE.11P);rm(logit.SRE.11P.pooled);rm(dt.SRE.11P.pooled)

#ERE, SR1
#infer model
dt.ERE.SR1.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.pooled.new.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled, 
                                     y=dt.SR1.coding$ERE.pooled.new.class[!is.na(dt.SR1.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.SR1.pooled, file="10_alt-class_out/scheme5/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="10_alt-class_out/scheme5/glmnet.cr.ERE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.SR1 <- which(abs(logit.ERE.SR1.pooled$lambda-(1e-5))==min(abs(logit.ERE.SR1.pooled$lambda-(1e-5))))

#predict missing
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m,s=step.ERE.SR1)
save(ERE.SR1.predictions, file="10_alt-class_out/scheme5/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.new.prediction <- ERE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme5/dt.SR1.coding.Rda")
rm(ERE.SR1.predictions);rm(step.ERE.SR1);rm(logit.ERE.SR1.pooled);rm(dt.ERE.SR1.pooled)

#SRE, SR1
#infer model
dt.SRE.SR1.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$SRE.pooled.new.class)])
logit.SRE.SR1.pooled <- glmnet.cr(x=dt.SRE.SR1.pooled, 
                                     y=dt.SR1.coding$SRE.pooled.new.class[!is.na(dt.SR1.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.SR1.pooled, file="10_alt-class_out/scheme5/dt.SRE.SR1.pooled")
save(logit.SRE.SR1.pooled, file="10_alt-class_out/scheme5/glmnet.cr.SRE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.SR1 <- which(abs(logit.SRE.SR1.pooled$lambda-(1e-5))==min(abs(logit.SRE.SR1.pooled$lambda-(1e-5))))

#predict missing
SRE.SR1.predictions <- fitted(logit.SRE.SR1.pooled, newx=dt.m,s=step.SRE.SR1)
save(SRE.SR1.predictions, file="10_alt-class_out/scheme5/SRE.SR1.predictions.Rda")
dt.SR1.coding$SRE.new.prediction <- SRE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme5/dt.SR1.coding.Rda")
rm(SRE.SR1.predictions);rm(step.SRE.SR1);rm(logit.SRE.SR1.pooled);rm(dt.SRE.SR1.pooled);rm(dt.m)


#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.new.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.new.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.new.prediction]
dt.SR1.coding[,ERE.full.class := ERE.pooled.new.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.new.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.new.prediction]


#define specificity class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]


#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme5/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme5/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme5/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme5/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme5/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme5/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme5/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme5/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme5/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme5/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5)+5,5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme5/all.pos.11P.Rda")
load("10_alt-class_out/scheme5/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme5/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme5/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme5/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme5/all.pos.SR1.Rda")
load("10_alt-class_out/scheme5/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme5/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme5/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme5/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme5/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])


###########################################################################################
#scheme6:
#for the AncSR1 wildtype, use the level of activation of SR1:EGKA in the bulk assay rather than the side-by-side ctrl
#SRE classes stay the same

dt.11P <- read.table(file="1_calc-meanF_11P_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

dt.SR1 <- read.table(file="2_calc-meanF_SR1_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

#define null distributions for hypothesis: inactive (pretty much only need for defining the positive bins, since I don't need to re-do these p.values)
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

#ERE pooled, SR1
ERE.pooled.null.SR1 <- data.frame(bin=1:25)
breaks <- cut2(dt.SR1.stop[!is.na(ERE.pooled.meanF),ERE.pooled.cfu],m=1000,g=25,onlycuts=T);breaks[26] <- breaks[26]+500000
for(i in 1:nrow(ERE.pooled.null.SR1)){
  ERE.pooled.null.SR1$range.cfu[i] <- list(c(breaks[i],breaks[i+1]))
  data <- dt.SR1.stop[ERE.pooled.cfu>=ERE.pooled.null.SR1$range.cfu[i][[1]][[1]] & ERE.pooled.cfu<ERE.pooled.null.SR1$range.cfu[i][[1]][[2]],]
  ERE.pooled.null.SR1$mean.meanF[i] <- mean(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$sd.meanF[i] <- sd(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$median.cfu[i] <- median(data$ERE.pooled.cfu,na.rm=T)
  ERE.pooled.null.SR1$max.meanF[i] <- max(data$ERE.pooled.meanF,na.rm=T)
  ERE.pooled.null.SR1$list.means[i] <- list(data[!is.na(ERE.pooled.meanF),ERE.pooled.meanF])
}


#determine meanF of stop variants
ERE.pooled.stop.meanF.11P <- mean(dt.11P.stop$ERE.pooled.meanF,na.rm=T)
ERE.pooled.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.pooled.meanF,na.rm=T)

#read in parallel control flow cytometry values for "wt" controls
ERE.rep1.egka.11P <- log(read.csv(file="3_classify-variants_define-distributions_in/l6_egka_ctrl_distribution.csv",header=F)$V1[1:7948])
ERE.rep2.egka.11P <- (log(read.csv(file="3_classify-variants_define-distributions_in/l13_egka_ctrl_distribution.csv",header=F)$V1)-0.01014)/1.00035
set.seed(103)
ERE.pooled.egka.11P <- sample(c(ERE.rep1.egka.11P,ERE.rep2.egka.11P))

ERE.rep1.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l10_egka_ctrl_distribution.csv",header=F)$V1)+0.59618)/1.01805; ERE.rep1.egka.SR1 <- ERE.rep1.egka.SR1[is.finite(ERE.rep1.egka.SR1)]
ERE.rep2.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l12_egka_ctrl_distribution.csv",header=F)$V1[1:length(ERE.rep1.egka.SR1)])-0.5260)/0.9633
set.seed(104)
ERE.pooled.egka.SR1 <- sample(c(ERE.rep1.egka.SR1,ERE.rep2.egka.SR1))

#shift controls over to 80th percentile between null and wildtype, to give most extreme instantiation of the null hypothesis that variant is weak that I want to test
ERE.pooled.egka.p80.11P <- ERE.pooled.egka.11P + (dt.SR1["EGKA",ERE.pooled.meanF]-mean(ERE.pooled.egka.11P)) - (dt.SR1["EGKA",ERE.pooled.meanF]-ERE.pooled.stop.meanF.11P)*0.2
ERE.pooled.egka.p80.SR1 <- ERE.pooled.egka.SR1 + (dt.SR1["EGKA",ERE.pooled.meanF]-mean(ERE.pooled.egka.SR1)) - (dt.SR1["EGKA",ERE.pooled.meanF]-ERE.pooled.stop.meanF.SR1)*0.2

#give sort bin boundaries
min.l6.b1 <- log(1); min.l6.b2 <- log(129.5); min.l6.b3 <- log(614.5); min.l6.b4 <- log(1284.5); max.l6.b4 <- log(262144)
min.l13.b1<- log(1); min.l13.b2<- (log(137.5)-0.01014)/1.00035; min.l13.b3<- (log(329.5)-0.01014)/1.00035; min.l13.b4<- (log(938.5)-0.01014)/1.00035; max.l13.b4 <- log(262144)
min.l10.b1 <- log(1); min.l10.b2 <- (log(73.5)+0.59618)/1.01805; min.l10.b3 <- (log(170.5)+0.59618)/1.01805; min.l10.b4 <- (log(374.5)+0.59618)/1.01805; max.l10.b4 <- log(262144)
min.l12.b1 <- log(1); min.l12.b2 <- (log(175.5)-0.5260)/0.9633; min.l12.b3 <- (log(374.5)-0.5260)/0.9633; min.l12.b4 <- (log(933.5)-0.5260)/0.9633; max.l12.b4 <- log(262144)

#interval censor wildtype observations
ERE.pooled.egka.p80.11P.cens.left <- c(sapply(1:round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l6.b2){return(min.l6.b1)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b3){return(min.l6.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b4){return(min.l6.b3)}else{return(min.l6.b4)}),sapply(round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)+1):length(ERE.pooled.egka.p80.11P), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l13.b2){return(min.l13.b1)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b3){return(min.l13.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b4){return(min.l13.b3)}else{return(min.l13.b4)}))
ERE.pooled.egka.p80.11P.cens.right <- c(sapply(1:round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l6.b2){return(min.l6.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b3){return(min.l6.b3)}else if(ERE.pooled.egka.p80.11P[x]<min.l6.b4){return(min.l6.b4)}else{return(max.l6.b4)}),sapply(round(sum(dt.11P$ERE.rep1.cfu)/sum(dt.11P$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.11P)+1):length(ERE.pooled.egka.p80.11P), function(x) if(ERE.pooled.egka.p80.11P[x]<min.l13.b2){return(min.l13.b2)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b3){return(min.l13.b3)}else if(ERE.pooled.egka.p80.11P[x]<min.l13.b4){return(min.l13.b4)}else{return(max.l13.b4)}))
ERE.pooled.egka.p80.11P.cens <- data.frame(left=ERE.pooled.egka.p80.11P.cens.left,right=ERE.pooled.egka.p80.11P.cens.right)

ERE.pooled.egka.p80.SR1.cens.left <- c(sapply(1:round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l10.b2){return(min.l10.b1)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b3){return(min.l10.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b4){return(min.l10.b3)}else{return(min.l10.b4)}),sapply(round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)+1):length(ERE.pooled.egka.p80.SR1), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l12.b2){return(min.l12.b1)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b3){return(min.l12.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b4){return(min.l12.b3)}else{return(min.l12.b4)}))
ERE.pooled.egka.p80.SR1.cens.right <- c(sapply(1:round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l10.b2){return(min.l10.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b3){return(min.l10.b3)}else if(ERE.pooled.egka.p80.SR1[x]<min.l10.b4){return(min.l10.b4)}else{return(max.l10.b4)}),sapply(round(sum(dt.SR1$ERE.rep1.cfu)/sum(dt.SR1$ERE.pooled.cfu)*length(ERE.pooled.egka.p80.SR1)+1):length(ERE.pooled.egka.p80.SR1), function(x) if(ERE.pooled.egka.p80.SR1[x]<min.l12.b2){return(min.l12.b2)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b3){return(min.l12.b3)}else if(ERE.pooled.egka.p80.SR1[x]<min.l12.b4){return(min.l12.b4)}else{return(max.l12.b4)}))
ERE.pooled.egka.p80.SR1.cens <- data.frame(left=ERE.pooled.egka.p80.SR1.cens.left,right=ERE.pooled.egka.p80.SR1.cens.right)

boot.cens.dist <- function(data, cfu, rep){
  sample.means <- vector()
  while(length(sample.means)<rep){
    boot <- data[sample(1:nrow(data),cfu,replace=T),]
    sample.means <- c(sample.means,tryCatch( summary(fitdistcens(boot,"logis"))$estimate["location"],error=function(e){return(as.numeric(NA))}))
  }
  return(sample.means)
}

ERE.pooled.pos.11P <- data.frame(bin=1:25); ERE.pooled.pos.11P$range.cfu <- ERE.pooled.null.11P$range.cfu; ERE.pooled.pos.11P$median.cfu <- ERE.pooled.null.11P$median.cfu
for(i in 1:nrow(ERE.pooled.pos.11P)){
  ERE.pooled.pos.11P$list.means[i] <- list(boot.cens.dist(ERE.pooled.egka.p80.11P.cens, round(ERE.pooled.pos.11P$median.cfu[i]),10000))
  ERE.pooled.pos.11P$mean.meanF[i] <- mean(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.11P$sd.meanF[i] <- sd(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.11P$max.meanF[i] <- max(ERE.pooled.pos.11P[i,"list.means"][[1]],na.rm=T)
}
save(ERE.pooled.pos.11P,file="10_alt-class_out/scheme6/ERE.pooled.pos.11P.Rda")

ERE.pooled.pos.SR1 <- data.frame(bin=1:25); ERE.pooled.pos.SR1$range.cfu <- ERE.pooled.null.SR1$range.cfu; ERE.pooled.pos.SR1$median.cfu <- ERE.pooled.null.SR1$median.cfu
for(i in 1:nrow(ERE.pooled.pos.SR1)){
  ERE.pooled.pos.SR1$list.means[i] <- list(boot.cens.dist(ERE.pooled.egka.p80.SR1.cens, round(ERE.pooled.pos.SR1$median.cfu[i]),10000))
  ERE.pooled.pos.SR1$mean.meanF[i] <- mean(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.SR1$sd.meanF[i] <- sd(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
  ERE.pooled.pos.SR1$max.meanF[i] <- max(ERE.pooled.pos.SR1[i,"list.means"][[1]],na.rm=T)
}
save(ERE.pooled.pos.SR1,file="10_alt-class_out/scheme6/ERE.pooled.pos.SR1.Rda")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt","ERE.pooled.pos.11P","ERE.pooled.pos.SR1"))])


dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

#define function to get p value for test of non-inferiority versus evolutionary wildtype (with 50% margin of indistinguishability)
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

dt.11P.coding[,ERE.pooled.p.weak := get.p.weak(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.pos.11P),by=AAseq]
dt.SR1.coding[,ERE.pooled.p.weak := get.p.weak(ERE.pooled.cfu,ERE.pooled.meanF,ERE.pooled.pos.SR1),by=AAseq]

dt.11P.coding[,ERE.pooled.new.class:=ERE.pooled.class]
dt.11P.coding[ERE.pooled.class %in% c("weak","strong"),ERE.pooled.new.class:="weak"] #revert all strongs to weak
dt.11P.coding[,SRE.pooled.new.class:=SRE.pooled.class]
dt.SR1.coding[,ERE.pooled.new.class:=ERE.pooled.class]
dt.SR1.coding[ERE.pooled.class %in% c("weak","strong"),ERE.pooled.new.class:="weak"] #revert all strongs to weak
dt.SR1.coding[,SRE.pooled.new.class:=SRE.pooled.class]


#determine p-value cutoff for hypothesis: weak (alt hypothesis: strong binder) that controls FDR to 5%
#ERE pooled, 11P
p.values.coding <- sort(c(dt.11P.coding[ERE.pooled.new.class=="weak",ERE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i
p.weak.ERE.pooled.cutoff.11P <- p.values.coding[i] 

#ERE pooled, SR1
p.values.coding <- sort(c(dt.SR1.coding[ERE.pooled.new.class=="weak",ERE.pooled.p.weak]))
i <- 1
while( p.values.coding[i]/(length(which(p.values.coding <= p.values.coding[i]))/length(p.values.coding)) <= 0.05 ){
  i <- i + 1
}
i
p.weak.ERE.pooled.cutoff.SR1 <- p.values.coding[i] 

#assign as strong if p.value below p.value.cutoff
dt.11P.coding[ERE.pooled.p.weak<p.weak.ERE.pooled.cutoff.11P & ERE.pooled.new.class=="weak",ERE.pooled.new.class:="strong"]
dt.SR1.coding[ERE.pooled.p.weak<p.weak.ERE.pooled.cutoff.SR1 & ERE.pooled.new.class=="weak",ERE.pooled.new.class:="strong"]

dt.11P.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.11P.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]

#build ordinal logistic regression models for each of the ERE experiments with the new classifications to predict missing genos
#ERE, 11P
#infer model
dt.ERE.11P.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.pooled.new.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled, 
                                     y=dt.11P.coding$ERE.pooled.new.class[!is.na(dt.11P.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.11P.pooled, file="10_alt-class_out/scheme6/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="10_alt-class_out/scheme6/glmnet.cr.ERE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.11P <- which(abs(logit.ERE.11P.pooled$lambda-(1e-5))==min(abs(logit.ERE.11P.pooled$lambda-(1e-5))))

#predict missing
dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m,s=step.ERE.11P)
save(ERE.11P.predictions, file="10_alt-class_out/scheme6/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.new.prediction <- ERE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme6/dt.11P.coding.Rda")
rm(ERE.11P.predictions);rm(step.ERE.11P);rm(logit.ERE.11P.pooled);rm(dt.ERE.11P.pooled)

#ERE, SR1
#infer model
dt.ERE.SR1.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.pooled.new.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled, 
                                     y=dt.SR1.coding$ERE.pooled.new.class[!is.na(dt.SR1.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.SR1.pooled, file="10_alt-class_out/scheme6/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="10_alt-class_out/scheme6/glmnet.cr.ERE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.SR1 <- which(abs(logit.ERE.SR1.pooled$lambda-(1e-5))==min(abs(logit.ERE.SR1.pooled$lambda-(1e-5))))

#predict missing
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m,s=step.ERE.SR1)
save(ERE.SR1.predictions, file="10_alt-class_out/scheme6/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.new.prediction <- ERE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme6/dt.SR1.coding.Rda")
rm(ERE.SR1.predictions);rm(step.ERE.SR1);rm(logit.ERE.SR1.pooled);rm(dt.ERE.SR1.pooled)

#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.new.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.new.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]
dt.SR1.coding[,ERE.full.class := ERE.pooled.new.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.new.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]


#define specificity class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]


#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme6/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme6/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme6/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme6/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme6/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme6/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme6/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme6/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme6/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme6/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5)+5,5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme6/all.pos.11P.Rda")
load("10_alt-class_out/scheme6/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme6/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme6/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme6/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme6/all.pos.SR1.Rda")
load("10_alt-class_out/scheme6/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme6/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme6/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme6/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme6/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])

###########################################################################################
#scheme7:
#use the main-text data classification scheme, but only keep "high-quality" classifications: classifications that were the same between both replicates on an RE.
#if the two replicates differ, classify as "NA" and fill in via a newly derived regression

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

dt.11P.coding[,ERE.pooled.new.class:=as.character(NA)]
dt.11P.coding[ERE.rep1.class=="null" & ERE.rep2.class=="null",ERE.pooled.new.class:="null"]
dt.11P.coding[ERE.rep1.class=="weak" & ERE.rep2.class=="weak",ERE.pooled.new.class:="weak"]
dt.11P.coding[ERE.rep1.class=="strong" & ERE.rep2.class=="strong",ERE.pooled.new.class:="strong"]
dt.11P.coding[,SRE.pooled.new.class:=as.character(NA)]
dt.11P.coding[SRE.rep1.class=="null" & SRE.rep2.class=="null",SRE.pooled.new.class:="null"]
dt.11P.coding[SRE.rep1.class=="weak" & SRE.rep2.class=="weak",SRE.pooled.new.class:="weak"]
dt.11P.coding[SRE.rep1.class=="strong" & SRE.rep2.class=="strong",SRE.pooled.new.class:="strong"]
dt.SR1.coding[,ERE.pooled.new.class:=as.character(NA)]
dt.SR1.coding[ERE.rep1.class=="null" & ERE.rep2.class=="null",ERE.pooled.new.class:="null"]
dt.SR1.coding[ERE.rep1.class=="weak" & ERE.rep2.class=="weak",ERE.pooled.new.class:="weak"]
dt.SR1.coding[ERE.rep1.class=="strong" & ERE.rep2.class=="strong",ERE.pooled.new.class:="strong"]
dt.SR1.coding[,SRE.pooled.new.class:=as.character(NA)]
dt.SR1.coding[SRE.rep1.class=="null" & SRE.rep2.class=="null",SRE.pooled.new.class:="null"]
dt.SR1.coding[SRE.rep1.class=="weak" & SRE.rep2.class=="weak",SRE.pooled.new.class:="weak"]
dt.SR1.coding[SRE.rep1.class=="strong" & SRE.rep2.class=="strong",SRE.pooled.new.class:="strong"]

dt.11P.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.11P.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong"))]

#build ordinal logistic regression models for each of the ERE experiments with the new classifications to predict missing genos
#ERE, 11P
#infer model
dt.ERE.11P.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.pooled.new.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled, 
                                     y=dt.11P.coding$ERE.pooled.new.class[!is.na(dt.11P.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.11P.pooled, file="10_alt-class_out/scheme7/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="10_alt-class_out/scheme7/glmnet.cr.ERE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.11P <- which(abs(logit.ERE.11P.pooled$lambda-(1e-5))==min(abs(logit.ERE.11P.pooled$lambda-(1e-5))))

#predict missing
dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m,s=step.ERE.11P)
save(ERE.11P.predictions, file="10_alt-class_out/scheme7/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.new.prediction <- ERE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme7/dt.11P.coding.Rda")
rm(ERE.11P.predictions);rm(step.ERE.11P);rm(logit.ERE.11P.pooled);rm(dt.ERE.11P.pooled)

#SRE, 11P
#infer model
dt.SRE.11P.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$SRE.pooled.new.class)])
logit.SRE.11P.pooled <- glmnet.cr(x=dt.SRE.11P.pooled, 
                                     y=dt.11P.coding$SRE.pooled.new.class[!is.na(dt.11P.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.11P.pooled, file="10_alt-class_out/scheme7/dt.SRE.11P.pooled")
save(logit.SRE.11P.pooled, file="10_alt-class_out/scheme7/glmnet.cr.SRE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.11P <- which(abs(logit.SRE.11P.pooled$lambda-(1e-5))==min(abs(logit.SRE.11P.pooled$lambda-(1e-5))))

#predict missing
SRE.11P.predictions <- fitted(logit.SRE.11P.pooled, newx=dt.m,s=step.SRE.11P)
save(SRE.11P.predictions, file="10_alt-class_out/scheme7/SRE.11P.predictions.Rda")
dt.11P.coding$SRE.new.prediction <- SRE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme7/dt.11P.coding.Rda")
rm(SRE.11P.predictions);rm(step.SRE.11P);rm(logit.SRE.11P.pooled);rm(dt.SRE.11P.pooled)


#ERE, SR1
#infer model
dt.ERE.SR1.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.pooled.new.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled, 
                                     y=dt.SR1.coding$ERE.pooled.new.class[!is.na(dt.SR1.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.SR1.pooled, file="10_alt-class_out/scheme7/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="10_alt-class_out/scheme7/glmnet.cr.ERE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.SR1 <- which(abs(logit.ERE.SR1.pooled$lambda-(1e-5))==min(abs(logit.ERE.SR1.pooled$lambda-(1e-5))))

#predict missing
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m,s=step.ERE.SR1)
save(ERE.SR1.predictions, file="10_alt-class_out/scheme7/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.new.prediction <- ERE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme7/dt.SR1.coding.Rda")
rm(ERE.SR1.predictions);rm(step.ERE.SR1);rm(logit.ERE.SR1.pooled);rm(dt.ERE.SR1.pooled)

#SRE, SR1
#infer model
dt.SRE.SR1.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$SRE.pooled.new.class)])
logit.SRE.SR1.pooled <- glmnet.cr(x=dt.SRE.SR1.pooled, 
                                     y=dt.SR1.coding$SRE.pooled.new.class[!is.na(dt.SR1.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.SR1.pooled, file="10_alt-class_out/scheme7/dt.SRE.SR1.pooled")
save(logit.SRE.SR1.pooled, file="10_alt-class_out/scheme7/glmnet.cr.SRE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.SR1 <- which(abs(logit.SRE.SR1.pooled$lambda-(1e-5))==min(abs(logit.SRE.SR1.pooled$lambda-(1e-5))))

#predict missing
SRE.SR1.predictions <- fitted(logit.SRE.SR1.pooled, newx=dt.m,s=step.SRE.SR1)
save(SRE.SR1.predictions, file="10_alt-class_out/scheme7/SRE.SR1.predictions.Rda")
dt.SR1.coding$SRE.new.prediction <- SRE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme7/dt.SR1.coding.Rda")
rm(SRE.SR1.predictions);rm(step.SRE.SR1);rm(logit.SRE.SR1.pooled);rm(dt.SRE.SR1.pooled)

#fill in ERE, SRE full, which is exptl class if >15 cfu and reps agree, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.new.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.new.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]
dt.SR1.coding[,ERE.full.class := ERE.pooled.new.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.new.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]


#define specificity class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]


#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme7/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme7/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme7/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme7/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme7/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme7/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme7/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme7/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme7/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme7/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5)+5,5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme7/all.pos.11P.Rda")
load("10_alt-class_out/scheme7/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme7/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme7/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme7/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme7/all.pos.SR1.Rda")
load("10_alt-class_out/scheme7/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme7/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme7/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme7/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme7/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])


###########################################################################################
#scheme8:
#add an upper margin on the level of binding that can be tolerated as "acceptable" -- that is, too good of binding is also not allowed

#get 80%ile cutoffs for ERE, SRE strong: this is EGKA in the SR1 library, GSKV in the 11P library
dt.11P <- read.table(file="1_calc-meanF_11P_out/dt_output_11P.csv",header=TRUE,sep=",")
dt.11P <- data.table(dt.11P);setkey(dt.11P,AAseq)

dt.11P.stop <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[grep("\\*",dt.11P$AAseq)],]
dt.11P.coding <- dt.11P[dt.11P$AAseq %in% dt.11P$AAseq[-grep("\\*",dt.11P$AAseq)],]

dt.SR1 <- read.table(file="2_calc-meanF_SR1_out/dt_output_SR1.csv",header=TRUE,sep=",")
dt.SR1 <- data.table(dt.SR1);setkey(dt.SR1,AAseq)

dt.SR1.stop <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[grep("\\*",dt.SR1$AAseq)],]
dt.SR1.coding <- dt.SR1[dt.SR1$AAseq %in% dt.SR1$AAseq[-grep("\\*",dt.SR1$AAseq)],]

SRE.pooled.stop.meanF.11P <- mean(dt.11P.stop$SRE.pooled.meanF,na.rm=T)
ERE.pooled.stop.meanF.SR1 <- mean(dt.SR1.stop$ERE.pooled.meanF,na.rm=T)

ERE.rep1.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l10_egka_ctrl_distribution.csv",header=F)$V1)+0.59618)/1.01805; ERE.rep1.egka.SR1 <- ERE.rep1.egka.SR1[is.finite(ERE.rep1.egka.SR1)]
ERE.rep2.egka.SR1 <- (log(read.csv(file="3_classify-variants_define-distributions_in/l12_egka_ctrl_distribution.csv",header=F)$V1[1:length(ERE.rep1.egka.SR1)])-0.5260)/0.9633
set.seed(104)
ERE.pooled.egka.SR1 <- sample(c(ERE.rep1.egka.SR1,ERE.rep2.egka.SR1))
ERE.lower.cutoff <- (mean(ERE.pooled.egka.SR1)-ERE.pooled.stop.meanF.SR1)*0.8+ERE.pooled.stop.meanF.SR1
SRE.lower.cutoff <- (dt.11P.coding["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.11P)*0.8+SRE.pooled.stop.meanF.11P
ERE.upper.cutoff <- (mean(ERE.pooled.egka.SR1)-ERE.pooled.stop.meanF.SR1)*0.2+mean(ERE.pooled.egka.SR1)
SRE.upper.cutoff <- (dt.11P.coding["GSKV",SRE.pooled.meanF]-SRE.pooled.stop.meanF.11P)*0.2+dt.11P.coding["GSKV",SRE.pooled.meanF]

rm(dt.SR1);rm(dt.SR1.coding);rm(dt.SR1.stop);rm(dt.11P);rm(dt.11P.coding);rm(dt.11P.stop);rm(ERE.rep1.egka.SR1);rm(ERE.rep2.egka.SR1);rm(ERE.pooled.egka.SR1);rm(ERE.pooled.stop.meanF.SR1);rm(SRE.pooled.stop.meanF.11P)

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
dt.11P.coding[,ERE.pooled.new.class := as.character(NA)]
dt.11P.coding[ERE.pooled.class=="null",ERE.pooled.new.class:="null"]
dt.11P.coding[ERE.pooled.class%in%c("weak","strong"),ERE.pooled.new.class:="weak"]
dt.11P.coding[ERE.pooled.class%in%c("weak","strong") & ERE.pooled.meanF > ERE.lower.cutoff & ERE.pooled.meanF <= ERE.upper.cutoff,ERE.pooled.new.class:="strong"]
dt.11P.coding[ERE.pooled.class%in%c("weak","strong") & ERE.pooled.meanF > ERE.upper.cutoff,ERE.pooled.new.class:="null2"]

dt.11P.coding[,SRE.pooled.new.class := as.character(NA)]
dt.11P.coding[SRE.pooled.class=="null",SRE.pooled.new.class:="null"]
dt.11P.coding[SRE.pooled.class%in%c("weak","strong"),SRE.pooled.new.class:="weak"]
dt.11P.coding[SRE.pooled.class%in%c("weak","strong") & SRE.pooled.meanF > SRE.lower.cutoff & SRE.pooled.meanF <= SRE.upper.cutoff,SRE.pooled.new.class:="strong"]
dt.11P.coding[SRE.pooled.class%in%c("weak","strong") & SRE.pooled.meanF > SRE.upper.cutoff,SRE.pooled.new.class:="null2"]

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)
dt.SR1.coding[,ERE.pooled.new.class := as.character(NA)]
dt.SR1.coding[ERE.pooled.class=="null",ERE.pooled.new.class:="null"]
dt.SR1.coding[ERE.pooled.class%in%c("weak","strong"),ERE.pooled.new.class:="weak"]
dt.SR1.coding[ERE.pooled.class%in%c("weak","strong") & ERE.pooled.meanF > ERE.lower.cutoff & ERE.pooled.meanF <= ERE.upper.cutoff,ERE.pooled.new.class:="strong"]
dt.SR1.coding[ERE.pooled.class%in%c("weak","strong") & ERE.pooled.meanF > ERE.upper.cutoff,ERE.pooled.new.class:="null2"]

dt.SR1.coding[,SRE.pooled.new.class := as.character(NA)]
dt.SR1.coding[SRE.pooled.class=="null",SRE.pooled.new.class:="null"]
dt.SR1.coding[SRE.pooled.class%in%c("weak","strong"),SRE.pooled.new.class:="weak"]
dt.SR1.coding[SRE.pooled.class%in%c("weak","strong") & SRE.pooled.meanF > SRE.lower.cutoff & SRE.pooled.meanF <= SRE.upper.cutoff,SRE.pooled.new.class:="strong"]
dt.SR1.coding[SRE.pooled.class%in%c("weak","strong") & SRE.pooled.meanF > SRE.upper.cutoff,SRE.pooled.new.class:="null2"]

dt.11P.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong","null2"))]
dt.11P.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong","null2"))]
dt.SR1.coding[,ERE.pooled.new.class:=factor(ERE.pooled.new.class,levels=c("null","weak","strong","null2"))]
dt.SR1.coding[,SRE.pooled.new.class:=factor(SRE.pooled.new.class,levels=c("null","weak","strong","null2"))]

#build ordinal logistic regression models for each of the four experiments with the new classifications to predict missing genos
#ERE, 11P
#infer model
dt.ERE.11P.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.pooled.new.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled, 
                                     y=dt.11P.coding$ERE.pooled.new.class[!is.na(dt.11P.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.11P.pooled, file="10_alt-class_out/scheme8/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="10_alt-class_out/scheme8/glmnet.cr.ERE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.11P <- which(abs(logit.ERE.11P.pooled$lambda-(1e-5))==min(abs(logit.ERE.11P.pooled$lambda-(1e-5))))

#predict missing
dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m,s=step.ERE.11P)
save(ERE.11P.predictions, file="10_alt-class_out/scheme8/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.new.prediction <- ERE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme8/dt.11P.coding.Rda")
rm(ERE.11P.predictions);rm(step.ERE.11P);rm(logit.ERE.11P.pooled);rm(dt.ERE.11P.pooled)

#SRE, 11P
#infer model
dt.SRE.11P.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$SRE.pooled.new.class)])
logit.SRE.11P.pooled <- glmnet.cr(x=dt.SRE.11P.pooled, 
                                     y=dt.11P.coding$SRE.pooled.new.class[!is.na(dt.11P.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.11P.pooled, file="10_alt-class_out/scheme8/dt.SRE.11P.pooled")
save(logit.SRE.11P.pooled, file="10_alt-class_out/scheme8/glmnet.cr.SRE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.11P <- which(abs(logit.SRE.11P.pooled$lambda-(1e-5))==min(abs(logit.SRE.11P.pooled$lambda-(1e-5))))

#predict missing
SRE.11P.predictions <- fitted(logit.SRE.11P.pooled, newx=dt.m,s=step.SRE.11P)
save(SRE.11P.predictions, file="10_alt-class_out/scheme8/SRE.11P.predictions.Rda")
dt.11P.coding$SRE.new.prediction <- SRE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme8/dt.11P.coding.Rda")
rm(SRE.11P.predictions);rm(step.SRE.11P);rm(logit.SRE.11P.pooled);rm(dt.SRE.11P.pooled)

#ERE, SR1
#infer model
dt.ERE.SR1.pooled <- model.matrix(ERE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.pooled.new.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled, 
                                     y=dt.SR1.coding$ERE.pooled.new.class[!is.na(dt.SR1.coding$ERE.pooled.new.class)],
                                     maxit=500)
save(dt.ERE.SR1.pooled, file="10_alt-class_out/scheme8/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="10_alt-class_out/scheme8/glmnet.cr.ERE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.SR1 <- which(abs(logit.ERE.SR1.pooled$lambda-(1e-5))==min(abs(logit.ERE.SR1.pooled$lambda-(1e-5))))

#predict missing
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m,s=step.ERE.SR1)
save(ERE.SR1.predictions, file="10_alt-class_out/scheme8/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.new.prediction <- ERE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme8/dt.SR1.coding.Rda")
rm(ERE.SR1.predictions);rm(step.ERE.SR1);rm(logit.ERE.SR1.pooled);rm(dt.ERE.SR1.pooled)

#SRE, SR1
#infer model
dt.SRE.SR1.pooled <- model.matrix(SRE.pooled.new.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$SRE.pooled.new.class)])
logit.SRE.SR1.pooled <- glmnet.cr(x=dt.SRE.SR1.pooled, 
                                     y=dt.SR1.coding$SRE.pooled.new.class[!is.na(dt.SR1.coding$SRE.pooled.new.class)],
                                     maxit=500)
save(dt.SRE.SR1.pooled, file="10_alt-class_out/scheme8/dt.SRE.SR1.pooled")
save(logit.SRE.SR1.pooled, file="10_alt-class_out/scheme8/glmnet.cr.SRE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.SR1 <- which(abs(logit.SRE.SR1.pooled$lambda-(1e-5))==min(abs(logit.SRE.SR1.pooled$lambda-(1e-5))))

#predict missing
SRE.SR1.predictions <- fitted(logit.SRE.SR1.pooled, newx=dt.m,s=step.SRE.SR1)
save(SRE.SR1.predictions, file="10_alt-class_out/scheme8/SRE.SR1.predictions.Rda")
dt.SR1.coding$SRE.new.prediction <- SRE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme8/dt.SR1.coding.Rda")
rm(SRE.SR1.predictions);rm(step.SRE.SR1);rm(logit.SRE.SR1.pooled);rm(dt.SRE.SR1.pooled);rm(dt.m)


#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.pooled.new.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.11P.coding[,SRE.full.class := SRE.pooled.new.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.new.prediction]
dt.SR1.coding[,ERE.full.class := ERE.pooled.new.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.new.prediction]
dt.SR1.coding[,SRE.full.class := SRE.pooled.new.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.new.prediction]


#define specificity class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null2",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="null2" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null2",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="null2" & SRE.full.class=="strong",specificity:="promiscuous"]


#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme8/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.11P$AAseq)==FALSE){
  ERE.specific.11P <- data.table(rbind(as.data.frame(ERE.specific.11P),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.11P,AAseq)
}


#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme8/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme8/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme8/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme8/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme8/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme8/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme8/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme8/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme8/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5),5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme8/all.pos.11P.Rda")
load("10_alt-class_out/scheme8/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme8/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme8/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme8/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme8/all.pos.SR1.Rda")
load("10_alt-class_out/scheme8/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme8/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme8/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme8/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme8/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])


###########################################################################################
#scheme 9: use only experimental classifications for all variants, undetermined variants remain missing

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)
#define new specificity class, given entirely by ERE and SRE experimental classes
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.pooled.class=="null" & SRE.pooled.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.pooled.class=="weak" & SRE.pooled.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="weak",specificity:="promiscuous"]

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)
#define new specificity class, given entirely by ERE and SRE experimental classes
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.pooled.class=="null" & SRE.pooled.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.pooled.class=="weak" & SRE.pooled.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.pooled.class=="strong" & SRE.pooled.class=="weak",specificity:="promiscuous"]

#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme9/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to promiscuous if not already there
if("EGKA"%in%as.character(promiscuous.11P$AAseq)==FALSE){
  promiscuous.11P <- data.table(rbind(as.data.frame(promiscuous.11P),c("EGKA",NA,NA,"promiscuous")));setkey(promiscuous.11P,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme9/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme9/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme9/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme9/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme9/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme9/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme9/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme9/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme9/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-5,5*round(max(ERE.specific.11P$num.SRE.outcomes)/5),5),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-1,max(SRE.specific.11P$accessed.by),1),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme9/all.pos.11P.Rda")
load("10_alt-class_out/scheme9/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme9/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme9/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme9/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme9/all.pos.SR1.Rda")
load("10_alt-class_out/scheme9/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme9/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme9/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme9/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme9/egka-SR1-requires-permissives.txt")

rm(list=ls()[!(ls() %in% c("AA.nt.neighbors","AAs","get.3mut.paths.11P","get.3mut.paths.SR1","get.Hamming1.nt"))])



############################################################################################
#scheme 10: classify variants as null/weak/strong on each RE given the per-variant estimates of SE generated in script 9
load(file="9_quant-scale-landscape_out/dt.11P.coding.Rda")
load(file="9_quant-scale-landscape_out/dt.SR1.coding.Rda")

#clean up table, removing unneeded columns
dt.11P.coding[,c("ERE.rep1.b1","ERE.rep1.b2","ERE.rep1.b3","ERE.rep1.b4","ERE.rep2.b1","ERE.rep2.b2","ERE.rep2.b3","ERE.rep2.b4","SRE.rep1.b1","SRE.rep1.b2","SRE.rep1.b3","SRE.rep1.b4","SRE.rep2.b1","SRE.rep2.b2","SRE.rep2.b3","SRE.rep2.b4","ERE.rep1.class","ERE.rep2.class","ERE.pooled.class","SRE.rep1.class","SRE.rep2.class","SRE.pooled.class","ERE.prediction","SRE.prediction","ERE.full.class","SRE.full.class","specificity") := NULL]
dt.SR1.coding[,c("ERE.rep1.b1","ERE.rep1.b2","ERE.rep1.b3","ERE.rep1.b4","ERE.rep2.b1","ERE.rep2.b2","ERE.rep2.b3","ERE.rep2.b4","SRE.rep1.b1","SRE.rep1.b2","SRE.rep1.b3","SRE.rep1.b4","SRE.rep2.b1","SRE.rep2.b2","SRE.rep2.b3","SRE.rep2.b4","ERE.rep1.class","ERE.rep2.class","ERE.pooled.class","SRE.rep1.class","SRE.rep2.class","SRE.pooled.class","ERE.prediction","SRE.prediction","ERE.full.class","SRE.full.class","specificity") := NULL]
#make NA any meanF/SE from a variant with <15 cfu
dt.11P.coding[ERE.pooled.cfu<15, c("ERE.pooled.meanF","ERE.SE.meanF") := NA]
dt.11P.coding[SRE.pooled.cfu<15, c("SRE.pooled.meanF","SRE.SE.meanF") := NA]
dt.SR1.coding[ERE.pooled.cfu<15, c("ERE.pooled.meanF","ERE.SE.meanF") := NA]
dt.SR1.coding[SRE.pooled.cfu<15, c("SRE.pooled.meanF","SRE.SE.meanF") := NA]

#ask whether, for each variant, its mean is significantly greater than the average null variant (meaning it's positive, weak or strong); and for each positive variants, whether its significantly greater than the relevant ancestral reference (meaning it's strong)

#from script 3, pulled out reference points for comparison
null.11P.ERE <- 5.070092
null.11P.SRE <- 5.043342
null.SR1.ERE <- 4.982878
null.SR1.SRE <- 5.139408

pos.11P.ERE <- 6.765023
pos.11P.SRE <- 6.538804
pos.SR1.ERE <- 6.773614
pos.SR1.SRE <- 6.495246

dt.11P.coding[,ERE.p.null := pnorm(null.11P.ERE, ERE.pooled.meanF, ERE.SE.meanF),by=AAseq]
dt.11P.coding[,SRE.p.null := pnorm(null.11P.SRE, SRE.pooled.meanF, SRE.SE.meanF),by=AAseq]
dt.SR1.coding[,ERE.p.null := pnorm(null.SR1.ERE, ERE.pooled.meanF, ERE.SE.meanF),by=AAseq]
dt.SR1.coding[,SRE.p.null := pnorm(null.SR1.SRE, SRE.pooled.meanF, SRE.SE.meanF),by=AAseq]

dt.11P.coding[,ERE.p.weak := pnorm(pos.11P.ERE, ERE.pooled.meanF, ERE.SE.meanF),by=AAseq]
dt.11P.coding[,SRE.p.weak := pnorm(pos.11P.SRE, SRE.pooled.meanF, SRE.SE.meanF),by=AAseq]
dt.SR1.coding[,ERE.p.weak := pnorm(pos.SR1.ERE, ERE.pooled.meanF, ERE.SE.meanF),by=AAseq]
dt.SR1.coding[,SRE.p.weak := pnorm(pos.SR1.SRE, SRE.pooled.meanF, SRE.SE.meanF),by=AAseq]

#fdr adjust p-values for 5% FDR
#p.null pvalues
dt.11P.coding$ERE.q.null <- p.adjust(dt.11P.coding$ERE.p.null, method="BH")
dt.11P.coding$SRE.q.null <- p.adjust(dt.11P.coding$SRE.p.null, method="BH")
dt.SR1.coding$ERE.q.null <- p.adjust(dt.SR1.coding$ERE.p.null, method="BH")
dt.SR1.coding$SRE.q.null <- p.adjust(dt.SR1.coding$SRE.p.null, method="BH")

#for p.weak, make NA any variant that does not reject null
dt.11P.coding[ERE.q.null > 0.05, ERE.p.weak := NA]
dt.11P.coding[SRE.q.null > 0.05, SRE.p.weak := NA]
dt.SR1.coding[ERE.q.null > 0.05, ERE.p.weak := NA]
dt.SR1.coding[SRE.q.null > 0.05, SRE.p.weak := NA]

dt.11P.coding$ERE.q.weak <- p.adjust(dt.11P.coding$ERE.p.weak, method="BH")
dt.11P.coding$SRE.q.weak <- p.adjust(dt.11P.coding$SRE.p.weak, method="BH")
dt.SR1.coding$ERE.q.weak <- p.adjust(dt.SR1.coding$ERE.p.weak, method="BH")
dt.SR1.coding$SRE.q.weak <- p.adjust(dt.SR1.coding$SRE.p.weak, method="BH")

#classify based on q-values at 5% FDR
dt.11P.coding$ERE.class <- "null"
dt.11P.coding$SRE.class <- "null"
dt.SR1.coding$ERE.class <- "null"
dt.SR1.coding$SRE.class <- "null"

dt.11P.coding[is.na(ERE.pooled.meanF), ERE.class := NA]
dt.11P.coding[ERE.q.null < 0.05 & ERE.q.weak > 0.05, ERE.class := "weak"]
dt.11P.coding[ERE.q.weak < 0.05, ERE.class := "strong"]

dt.11P.coding[is.na(SRE.pooled.meanF), SRE.class := NA]
dt.11P.coding[SRE.q.null < 0.05 & SRE.q.weak > 0.05, SRE.class := "weak"]
dt.11P.coding[SRE.q.weak < 0.05, SRE.class := "strong"]

dt.SR1.coding[is.na(ERE.pooled.meanF), ERE.class := NA]
dt.SR1.coding[ERE.q.null < 0.05 & ERE.q.weak > 0.05, ERE.class := "weak"]
dt.SR1.coding[ERE.q.weak < 0.05, ERE.class := "strong"]

dt.SR1.coding[is.na(SRE.pooled.meanF), SRE.class := NA]
dt.SR1.coding[SRE.q.null < 0.05 & SRE.q.weak > 0.05, SRE.class := "weak"]
dt.SR1.coding[SRE.q.weak < 0.05, SRE.class := "strong"]

dt.11P.coding[,ERE.class:=factor(ERE.class,levels=c("null","weak","strong"))]
dt.11P.coding[,SRE.class:=factor(SRE.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,ERE.class:=factor(ERE.class,levels=c("null","weak","strong"))]
dt.SR1.coding[,SRE.class:=factor(SRE.class,levels=c("null","weak","strong"))]

#build ordinal logistic regression models for each of the four experiments with the new classifications to predict missing genos
#ERE, 11P
#infer model
dt.ERE.11P.pooled <- model.matrix(ERE.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$ERE.class)])
logit.ERE.11P.pooled <- glmnet.cr(x=dt.ERE.11P.pooled, 
                                  y=dt.11P.coding$ERE.class[!is.na(dt.11P.coding$ERE.class)],
                                  maxit=500)
save(dt.ERE.11P.pooled, file="10_alt-class_out/scheme10/dt.ERE.11P.pooled")
save(logit.ERE.11P.pooled, file="10_alt-class_out/scheme10/glmnet.cr.ERE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.11P <- which(abs(logit.ERE.11P.pooled$lambda-(1e-5))==min(abs(logit.ERE.11P.pooled$lambda-(1e-5))))

#predict missing
dt.m <- model.matrix(rep(1,160000) ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding)
ERE.11P.predictions <- fitted(logit.ERE.11P.pooled, newx=dt.m,s=step.ERE.11P)
save(ERE.11P.predictions, file="10_alt-class_out/scheme10/ERE.11P.predictions.Rda")
dt.11P.coding$ERE.prediction <- ERE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme10/dt.11P.coding.Rda")
rm(ERE.11P.predictions);rm(step.ERE.11P);rm(logit.ERE.11P.pooled);rm(dt.ERE.11P.pooled)

#SRE, 11P
#infer model
dt.SRE.11P.pooled <- model.matrix(SRE.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.11P.coding[!is.na(dt.11P.coding$SRE.class)])
logit.SRE.11P.pooled <- glmnet.cr(x=dt.SRE.11P.pooled, 
                                  y=dt.11P.coding$SRE.class[!is.na(dt.11P.coding$SRE.class)],
                                  maxit=500)
save(dt.SRE.11P.pooled, file="10_alt-class_out/scheme10/dt.SRE.11P.pooled")
save(logit.SRE.11P.pooled, file="10_alt-class_out/scheme10/glmnet.cr.SRE.11P.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.11P <- which(abs(logit.SRE.11P.pooled$lambda-(1e-5))==min(abs(logit.SRE.11P.pooled$lambda-(1e-5))))

#predict missing
SRE.11P.predictions <- fitted(logit.SRE.11P.pooled, newx=dt.m,s=step.SRE.11P)
save(SRE.11P.predictions, file="10_alt-class_out/scheme10/SRE.11P.predictions.Rda")
dt.11P.coding$SRE.prediction <- SRE.11P.predictions$class
save(dt.11P.coding, file="10_alt-class_out/scheme10/dt.11P.coding.Rda")
rm(SRE.11P.predictions);rm(step.SRE.11P);rm(logit.SRE.11P.pooled);rm(dt.SRE.11P.pooled)

#ERE, SR1
#infer model
dt.ERE.SR1.pooled <- model.matrix(ERE.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$ERE.class)])
logit.ERE.SR1.pooled <- glmnet.cr(x=dt.ERE.SR1.pooled, 
                                  y=dt.SR1.coding$ERE.class[!is.na(dt.SR1.coding$ERE.class)],
                                  maxit=500)
save(dt.ERE.SR1.pooled, file="10_alt-class_out/scheme10/dt.ERE.SR1.pooled")
save(logit.ERE.SR1.pooled, file="10_alt-class_out/scheme10/glmnet.cr.ERE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.ERE.SR1 <- which(abs(logit.ERE.SR1.pooled$lambda-(1e-5))==min(abs(logit.ERE.SR1.pooled$lambda-(1e-5))))

#predict missing
ERE.SR1.predictions <- fitted(logit.ERE.SR1.pooled, newx=dt.m,s=step.ERE.SR1)
save(ERE.SR1.predictions, file="10_alt-class_out/scheme10/ERE.SR1.predictions.Rda")
dt.SR1.coding$ERE.prediction <- ERE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme10/dt.SR1.coding.Rda")
rm(ERE.SR1.predictions);rm(step.ERE.SR1);rm(logit.ERE.SR1.pooled);rm(dt.ERE.SR1.pooled)

#SRE, SR1
#infer model
dt.SRE.SR1.pooled <- model.matrix(SRE.class ~ AA1+AA2+AA3+AA4+AA1*AA2+AA1*AA3+AA1*AA4+AA2*AA3+AA2*AA4+AA3*AA4,data=dt.SR1.coding[!is.na(dt.SR1.coding$SRE.class)])
logit.SRE.SR1.pooled <- glmnet.cr(x=dt.SRE.SR1.pooled, 
                                  y=dt.SR1.coding$SRE.class[!is.na(dt.SR1.coding$SRE.class)],
                                  maxit=500)
save(dt.SRE.SR1.pooled, file="10_alt-class_out/scheme10/dt.SRE.SR1.pooled")
save(logit.SRE.SR1.pooled, file="10_alt-class_out/scheme10/glmnet.cr.SRE.SR1.pooled.Rda")

#pull out step closest to lambda=1e-5
step.SRE.SR1 <- which(abs(logit.SRE.SR1.pooled$lambda-(1e-5))==min(abs(logit.SRE.SR1.pooled$lambda-(1e-5))))

#predict missing
SRE.SR1.predictions <- fitted(logit.SRE.SR1.pooled, newx=dt.m,s=step.SRE.SR1)
save(SRE.SR1.predictions, file="10_alt-class_out/scheme10/SRE.SR1.predictions.Rda")
dt.SR1.coding$SRE.prediction <- SRE.SR1.predictions$class
save(dt.SR1.coding, file="10_alt-class_out/scheme10/dt.SR1.coding.Rda")
rm(SRE.SR1.predictions);rm(step.SRE.SR1);rm(logit.SRE.SR1.pooled);rm(dt.SRE.SR1.pooled);rm(dt.m)


#fill in ERE, SRE full, which is exptl class if >15 cfu, predicted class otherwise
dt.11P.coding[,ERE.full.class := ERE.class]
dt.11P.coding[is.na(ERE.full.class),ERE.full.class := ERE.prediction]
dt.11P.coding[,SRE.full.class := SRE.class]
dt.11P.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]
dt.SR1.coding[,ERE.full.class := ERE.class]
dt.SR1.coding[is.na(ERE.full.class),ERE.full.class := ERE.prediction]
dt.SR1.coding[,SRE.full.class := SRE.class]
dt.SR1.coding[is.na(SRE.full.class),SRE.full.class := SRE.prediction]


#define specificity class
dt.11P.coding[,specificity:="null"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.11P.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.11P.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[,specificity:="null"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="null",specificity:="ERE-specific"]
dt.SR1.coding[ERE.full.class=="null" & SRE.full.class=="strong",specificity:="SRE-specific"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="strong",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="strong" & SRE.full.class=="weak",specificity:="promiscuous"]
dt.SR1.coding[ERE.full.class=="weak" & SRE.full.class=="strong",specificity:="promiscuous"]


#output table giving number of variants falling into each specificity class
write.csv(data.frame(AncSR1=table(dt.SR1.coding$specificity),AncSR1.11P=table(dt.11P.coding$specificity))[,c(1,2,4)],file="10_alt-class_out/scheme10/class-freqs.csv")

#re-analyses related to Figure 2: stochasticity and contingency in the 11P network
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq)),label=c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.11P))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.11P))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.11P)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.11P)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme10/11P_all-functional-positives.gexf")

#how many SRE-sp outcomes can EGKA access in <=3 nt mutations?
egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.11P),"unique outcomes in 11P space","\n","accesses GSKV?","GSKV"%in%egka.3mut.SRE.11P),file="10_alt-class_out/scheme10/egka-stochasticity.txt")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
  ERE.specific.11P$num.SRE.outcomes[i] <- length(ERE.specific.11P[i,all.SRE.muts][[1]])
}
save(ERE.specific.11P, file="10_alt-class_out/scheme10/ERE-specific.11P.Rda")
load(file="10_alt-class_out/scheme10/ERE-specific.11P.Rda")

#how many ERE-sp starting points access each SRE-specific outcome in 3 or less mutations?
for(i in 1:nrow(SRE.specific.11P)){
  print(i)
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count
}
save(SRE.specific.11P, file="10_alt-class_out/scheme10/SRE-specific.11P.Rda")
load(file="10_alt-class_out/scheme10/SRE-specific.11P.Rda")


# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap)[j]
    ERE.sp.overlap[geno.i,geno.j] <- sum(ERE.specific.11P[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap)){
  ERE.sp.overlap[i,i] <- NA
}
#take out cells involving an EREsp that accesses no SRE-sp outcomes
ERE.sp.overlap[as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap[,as.character(ERE.specific.11P[num.SRE.outcomes==0,AAseq])] <- NA
save(ERE.sp.overlap, file="10_alt-class_out/scheme10/ERE.sp.overlap.Rda")
load(file="10_alt-class_out/scheme10/ERE.sp.overlap.Rda")


pdf(file="./10_alt-class_out/scheme10/Fig2-hists.pdf",8,2.5)
par(mfrow=c(1,3))
par(mar=c(5,5,1,1))
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-10,10*round(max(ERE.specific.11P$num.SRE.outcomes)/10)+10,10),xlab="number of SRE-specific outcomes\naccessed in <= 3-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
abline(v=length(egka.3mut.SRE.11P),lty=2)
hist(SRE.specific.11P$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P$accessed.by)/5)+5,5),xlab="number of ERE-specific starting points\naccessing in <= 3-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(ERE.sp.overlap[which(ERE.specific.11P$num.SRE.outcomes>0),which(ERE.specific.11P$num.SRE.outcomes>0)],breaks=50,col="gray75",xlab="proportion of SRE-specific outcomes\naccessible from i also accessible from j",ylab="number of pairs of ERE-\nspecific starting points",main="")
dev.off()


#re-analyses related to Figure 3: effects of 11P on paths and evolvability
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.11P <- matrix(data=rep(0,(length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq))^2),nrow=length(ERE.specific.11P$AAseq)+length(SRE.specific.11P$AAseq)+length(promiscuous.11P$AAseq));colnames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq)); rownames(adj.m.11P) <- c(as.character(ERE.specific.11P$AAseq),as.character(SRE.specific.11P$AAseq),as.character(promiscuous.11P$AAseq))
for(i in 1:nrow(adj.m.11P)){
  if(row.names(adj.m.11P)[i] %in% as.character(ERE.specific.11P$AAseq)){
    adj.m.11P[i,colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i])] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(promiscuous.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% c(as.character(promiscuous.11P$AAseq),as.character(SRE.specific.11P$AAseq)))] <- 1
  }else if(row.names(adj.m.11P)[i] %in% as.character(SRE.specific.11P$AAseq)){
    adj.m.11P[i,(colnames(adj.m.11P) %in% get.Hamming1.nt(row.names(adj.m.11P)[i]) & colnames(adj.m.11P) %in% as.character(SRE.specific.11P$AAseq))] <- 1
  }else{print(i)}
}

all.pos.11P <- graph.adjacency(adj.m.11P,mode="directed")
save(all.pos.11P,file="10_alt-class_out/scheme10/all.pos.11P.Rda")
load("10_alt-class_out/scheme10/all.pos.11P.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  geno <- as.character(ERE.specific.11P$AAseq[i])
  distances <- distances(all.pos.11P, v=geno,to=as.character(SRE.specific.11P$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.11P$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.11P,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.11P$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.11P$shortest.path.length.to.SRE <- NA
ERE.specific.11P$num.shortest.paths <- NA
ERE.specific.11P$perm.discrete <- 0
ERE.specific.11P$prom <- 0
ERE.specific.11P$perm.prom <- 0
ERE.specific.11P$discrete <- 0
for(i in 1:nrow(ERE.specific.11P)){
  if(sum(!is.na(ERE.specific.11P[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.11P$shortest.path.length.to.SRE[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.11P$num.shortest.paths[i] <- length(ERE.specific.11P[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.11P$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.11P$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.11P[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.11P$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.11P$perm.discrete[i] <- ERE.specific.11P$perm.discrete[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.11P$prom[i] <- ERE.specific.11P$prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.11P$perm.prom[i] <- ERE.specific.11P$perm.prom[i] + 1/ERE.specific.11P$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.11P$discrete[i] <- ERE.specific.11P$discrete[i]+1/ERE.specific.11P$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.11P, file="10_alt-class_out/scheme10/ERE.specific.11P.Rda")
load("10_alt-class_out/scheme10/ERE.specific.11P.Rda")

#SR1 space
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#add EGKA to ERE.specific if not already there
if("EGKA"%in%as.character(ERE.specific.SR1$AAseq)==FALSE){
  ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)
}

#make file for visualization in gephi
nodes <- data.frame(id=1:(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq)),label=c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[1],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[1],nrow(promiscuous.SR1))),g=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[2],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[2],nrow(promiscuous.SR1))),b=c(rep(col2rgb(rgb(100,69,155,maxColorValue=255))[3],nrow(ERE.specific.SR1)),rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.SR1)),rep(col2rgb(rgb(111,204,221,maxColorValue=255))[3],nrow(promiscuous.SR1))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="10_alt-class_out/scheme10/SR1_all-functional-positives.gexf")

#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions
adj.m.SR1 <- matrix(data=rep(0,(length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq))^2),nrow=length(ERE.specific.SR1$AAseq)+length(SRE.specific.SR1$AAseq)+length(promiscuous.SR1$AAseq));colnames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq)); rownames(adj.m.SR1) <- c(as.character(ERE.specific.SR1$AAseq),as.character(SRE.specific.SR1$AAseq),as.character(promiscuous.SR1$AAseq))
for(i in 1:nrow(adj.m.SR1)){
  if(row.names(adj.m.SR1)[i] %in% as.character(ERE.specific.SR1$AAseq)){
    adj.m.SR1[i,colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i])] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(promiscuous.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% c(as.character(promiscuous.SR1$AAseq),as.character(SRE.specific.SR1$AAseq)))] <- 1
  }else if(row.names(adj.m.SR1)[i] %in% as.character(SRE.specific.SR1$AAseq)){
    adj.m.SR1[i,(colnames(adj.m.SR1) %in% get.Hamming1.nt(row.names(adj.m.SR1)[i]) & colnames(adj.m.SR1) %in% as.character(SRE.specific.SR1$AAseq))] <- 1
  }else{print(i)}
}

all.pos.SR1 <- graph.adjacency(adj.m.SR1,mode="directed")
save(all.pos.SR1,file="10_alt-class_out/scheme10/all.pos.SR1.Rda")
load("10_alt-class_out/scheme10/all.pos.SR1.Rda")

#shortest paths from ERE-specific to ANY SRE-specific variant?
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
  print(i)
  geno <- as.character(ERE.specific.SR1$AAseq[i])
  distances <- distances(all.pos.SR1, v=geno,to=as.character(SRE.specific.SR1$AAseq))
  min <- min(distances,na.rm=FALSE)
  if(is.finite(min)){
    outcomes <- colnames(distances)[distances==min]
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- list(all_shortest_paths(all.pos.SR1,from=geno,to=outcomes)$res)
  }else{
    ERE.specific.SR1$shortest.paths.to.SRE[i] <- NA
  }
}

#indicate, for each ERE-sp starting point, whether its closest SRE-sp neighbors can be accessed through shortest paths with prom int and/or discrete shifts, as well as whether it requires permissive RH steps
#types of trajectories: permissive no promiscuous (1), promiscuous no permissive (2), promiscuous and permissive (3), no promiscuous or permissive (4), or no path available (5)
ERE.specific.SR1$shortest.path.length.to.SRE <- NA
ERE.specific.SR1$num.shortest.paths <- NA
ERE.specific.SR1$perm.discrete <- 0
ERE.specific.SR1$prom <- 0
ERE.specific.SR1$perm.prom <- 0
ERE.specific.SR1$discrete <- 0
for(i in 1:nrow(ERE.specific.SR1)){
  if(sum(!is.na(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]]))>0){
    ERE.specific.SR1$shortest.path.length.to.SRE[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[1]])-1
    ERE.specific.SR1$num.shortest.paths[i] <- length(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]])
    for(k in 1:ERE.specific.SR1$num.shortest.paths[i]){
      promiscuous <- F
      permissive <- F
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(promiscuous.SR1$AAseq))>0){
        promiscuous <- T
      }
      if(sum(as_ids(ERE.specific.SR1[i,shortest.paths.to.SRE][[1]][[k]]) %in% as.character(ERE.specific.SR1$AAseq))>1){
        permissive <- T
      }
      if(promiscuous==F & permissive==T){
        ERE.specific.SR1$perm.discrete[i] <- ERE.specific.SR1$perm.discrete[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==F){
        ERE.specific.SR1$prom[i] <- ERE.specific.SR1$prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==T & permissive==T){
        ERE.specific.SR1$perm.prom[i] <- ERE.specific.SR1$perm.prom[i] + 1/ERE.specific.SR1$num.shortest.paths[i]
      }else if(promiscuous==F & permissive==F){
        ERE.specific.SR1$discrete[i] <- ERE.specific.SR1$discrete[i]+1/ERE.specific.SR1$num.shortest.paths[i]
      }
    }
  }
}
save(ERE.specific.SR1, file="10_alt-class_out/scheme10/ERE.specific.SR1.Rda")
load("10_alt-class_out/scheme10/ERE.specific.SR1.Rda")


#Fig3d-f equivalents
pdf(file="10_alt-class_out/scheme10/Fig3-hists.pdf",height=5,width=4,useDingbats=F)
par(mfrow=c(2,2))
#shortest path lengths
ERE.specific.11P.shortest.paths <- factor(ERE.specific.11P$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
ERE.specific.SR1.shortest.paths <- factor(ERE.specific.SR1$shortest.path.length.to.SRE,levels=1:max(c(ERE.specific.11P$shortest.path.length.to.SRE,ERE.specific.SR1$shortest.path.length.to.SRE),na.rm=T))
par(las=1)
par(mar=c(5,6,4,2))
barplot(rev(table(ERE.specific.SR1.shortest.paths,useNA="always")),horiz=T,xlab="number of ERE-specific starting points",ylab="Length of shortest path\nto SRE-specificity",main="AncSR1")
par(mar=c(5,2,4,6))
barplot(rev(table(ERE.specific.11P.shortest.paths,useNA="always")),horiz=T,main="AncSR1+11P",names.arg="")
#requirement for additional permissive steps within RH and promiscuous intermediates
par(mar=c(5,6,2,2))
barplot(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),names.arg=c("no path","neither","permissive+\npromiscuous","promiscuous","permissive"),horiz=T,xlab="number of ERE-specific starting points")
par(mar=c(5,2,2,6))
barplot(c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)),names.arg=c(""),horiz=T)
dev.off()

egka.3mut.paths.SR1 <- get.3mut.paths.SR1("EGKA")
labels.unordered.SR1 <- unique(c(egka.3mut.paths.SR1[[1]],egka.3mut.paths.SR1[[2]],egka.3mut.paths.SR1[[3]],egka.3mut.paths.SR1[[4]]))

egka.3mut.ERE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(ERE.specific.SR1$AAseq)]
egka.3mut.prom.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(promiscuous.SR1$AAseq)]
egka.3mut.SRE.SR1 <- labels.unordered.SR1[labels.unordered.SR1 %in% as.character(SRE.specific.SR1$AAseq)]

cat(c("EGKA accesses",length(egka.3mut.SRE.SR1),"unique outcomes in SR1 space in three steps\n","EGKA requires permissive:",ERE.specific.SR1["EGKA",prom]==F & ERE.specific.SR1["EGKA",discrete]==F,"\nEGKA shortest path:",ERE.specific.SR1["EGKA",shortest.path.length.to.SRE]),file="10_alt-class_out/scheme10/egka-SR1-requires-permissives.txt")


