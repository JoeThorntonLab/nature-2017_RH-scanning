#16 Jan 2017
#TNS
#script to make figures and diagrams illustrating connectivity of functional variants and trajectories through RH sequence space

setwd("path/to/source/directory") #setwd to source file directory

library(data.table)
library(rgexf)
library(igraph)

dt.11P.coding <- read.table(file="6_predict-missing_out/dt_output_11P.csv",header=TRUE,sep=",");dt.11P.coding <- data.table(dt.11P.coding);setkey(dt.11P.coding,AAseq)

dt.SR1.coding <- read.table(file="6_predict-missing_out/dt_output_SR1.csv",header=TRUE,sep=",");dt.SR1.coding <- data.table(dt.SR1.coding);setkey(dt.SR1.coding,AAseq)

#read in list of all 1-letter amino acid codes, and a table that gives all amino acids that can be encoded in adjacent (single-nt mutation) codons
AAs <- c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
AA.nt.neighbors <- read.csv(file="7_assess-connectivity_in/AA_single-mutant-codon-neighbors.csv",header=T,row.names=1,quote="")

#define function that takes a 4-aa RH sequence and returns all possible single-mutant neighbors
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

#11P space
ERE.specific.11P <- dt.11P.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.11P <- dt.11P.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.11P <- dt.11P.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#look at functional annotation of nearest neighbors:
ERE.specific.11P[,n.ERE.nt := sum(get.Hamming1.nt(AAseq) %in% ERE.specific.11P$AAseq),by=AAseq]
ERE.specific.11P[,n.prom.nt := sum(get.Hamming1.nt(AAseq) %in% promiscuous.11P$AAseq),by=AAseq]
ERE.specific.11P[,n.SRE.nt := sum(get.Hamming1.nt(AAseq) %in% SRE.specific.11P$AAseq),by=AAseq]

SRE.specific.11P[,n.ERE.nt := sum(get.Hamming1.nt(AAseq) %in% ERE.specific.11P$AAseq),by=AAseq]
SRE.specific.11P[,n.prom.nt := sum(get.Hamming1.nt(AAseq) %in% promiscuous.11P$AAseq),by=AAseq]
SRE.specific.11P[,n.SRE.nt := sum(get.Hamming1.nt(AAseq) %in% SRE.specific.11P$AAseq),by=AAseq]

promiscuous.11P[,n.ERE.nt := sum(get.Hamming1.nt(AAseq) %in% ERE.specific.11P$AAseq),by=AAseq]
promiscuous.11P[,n.prom.nt := sum(get.Hamming1.nt(AAseq) %in% promiscuous.11P$AAseq),by=AAseq]
promiscuous.11P[,n.SRE.nt := sum(get.Hamming1.nt(AAseq) %in% SRE.specific.11P$AAseq),by=AAseq]

#make node and edge data frames to convert to gexf format for visualization in gephi, for all variants that are strongly positive on at least one RE
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
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="7_assess-connectivity_out/gephi-files/11P_all-functional-positives_filled-genos.gexf")

##############################################################################################################
#SR1
ERE.specific.SR1 <- dt.SR1.coding[specificity=="ERE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
SRE.specific.SR1 <- dt.SR1.coding[specificity=="SRE-specific",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]
promiscuous.SR1 <- dt.SR1.coding[specificity=="promiscuous",.(AAseq,ERE.pooled.meanF,SRE.pooled.meanF,specificity)]

#Add back in EGKA, false negative in SR1:
ERE.specific.SR1 <- data.table(rbind(as.data.frame(ERE.specific.SR1),c("EGKA",NA,NA,"ERE-specific")));setkey(ERE.specific.SR1,AAseq)

#look at functional annotation of nearest neighbors:
ERE.specific.SR1[,n.ERE.nt := sum(get.Hamming1.nt(AAseq) %in% ERE.specific.SR1$AAseq),by=AAseq]
ERE.specific.SR1[,n.prom.nt := sum(get.Hamming1.nt(AAseq) %in% promiscuous.SR1$AAseq),by=AAseq]
ERE.specific.SR1[,n.SRE.nt := sum(get.Hamming1.nt(AAseq) %in% SRE.specific.SR1$AAseq),by=AAseq]

SRE.specific.SR1[,n.ERE.nt := sum(get.Hamming1.nt(AAseq) %in% ERE.specific.SR1$AAseq),by=AAseq]
SRE.specific.SR1[,n.prom.nt := sum(get.Hamming1.nt(AAseq) %in% promiscuous.SR1$AAseq),by=AAseq]
SRE.specific.SR1[,n.SRE.nt := sum(get.Hamming1.nt(AAseq) %in% SRE.specific.SR1$AAseq),by=AAseq]

promiscuous.SR1[,n.ERE.nt := sum(get.Hamming1.nt(AAseq) %in% ERE.specific.SR1$AAseq),by=AAseq]
promiscuous.SR1[,n.prom.nt := sum(get.Hamming1.nt(AAseq) %in% promiscuous.SR1$AAseq),by=AAseq]
promiscuous.SR1[,n.SRE.nt := sum(get.Hamming1.nt(AAseq) %in% SRE.specific.SR1$AAseq),by=AAseq]

#make node and edge data frames to convert to gexf format for visualization in gephi, for all variants that are strongly positive on at least one RE
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
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="7_assess-connectivity_out/gephi-files/SR1_all-functional-positives_filled-genos.gexf")


###########################################################################################
#11P igraph
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions, convert to igraph object
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

#identify shortest paths from each ERE-specific to it's closest SRE-specific outcomes
ERE.specific.11P[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.11P)){
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

###########################################################################################
#SR1 igraph
#for analyzing paths in igraph, create adjacency matrix for trajectories disallowing functional reversions, convert to igraph object
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

#identify shortest paths from each ERE-specific to it's closest SRE-specific outcomes
ERE.specific.SR1[,shortest.paths.to.SRE := NA]
for(i in 1:nrow(ERE.specific.SR1)){
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

######################################################################################################

#test for difference in robustness, evolvability in network
t.test(ERE.specific.11P$n.ERE.nt,ERE.specific.SR1$n.ERE.nt) #1.1e-8, goes from 1.67 ERE.sp neighbors to 3.44 ERE.sp neighbors
t.test(ERE.specific.11P$n.SRE.nt,ERE.specific.SR1$n.SRE.nt) #7.6e-7, goes from 0 to 0.27 SRE.sp neighbors
t.test(ERE.specific.11P$n.SRE.nt+ERE.specific.11P$n.prom.nt,ERE.specific.SR1$n.SRE.nt+ERE.specific.SR1$n.prom.nt) #4.368554e-18 #goes from 0.49 to 2.58 SRE-active neighbors
t.test(ERE.specific.11P$n.ERE.nt+ERE.specific.11P$n.SRE.nt+ERE.specific.11P$n.prom.nt,ERE.specific.SR1$n.ERE.nt+ERE.specific.SR1$n.SRE.nt+ERE.specific.SR1$n.prom.nt) #1.135911e-20 goes from 2.16 to 6.02 functional neighbors

#charts for comparing shortest path length, dependence on permissives, and promiscuous intermediates along shortest trajectories between AncSR1 and AncSR1+11P starting points
pdf(file="7_assess-connectivity_out/comparison-of-shortest-trajectories_SR1-v-11P.pdf",height=5,width=4,useDingbats=F)
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

chisq.test(cbind(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete))),simulate.p.value=T,B=1000) #simulate B=1e8 or more replicates to get more exact p-value than < x

chisq.test(cbind(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete))),simulate.p.value=T,B=1000)$expected
chisq.test(cbind(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete))),simulate.p.value=T,B=1000)$observed

cbind(c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete)),c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete)))
c(sum(is.na(ERE.specific.SR1$shortest.path.length.to.SRE)),sum(ERE.specific.SR1$discrete),sum(ERE.specific.SR1$perm.prom),sum(ERE.specific.SR1$prom),sum(ERE.specific.SR1$perm.discrete))
c(sum(is.na(ERE.specific.11P$shortest.path.length.to.SRE)),sum(ERE.specific.11P$discrete),sum(ERE.specific.11P$perm.prom),sum(ERE.specific.11P$prom),sum(ERE.specific.11P$perm.discrete))

#how many SRE-sp outcomes can EGKA and other starting points access in <=3 steps while passing througn functional intermediates?
#define function for 11P, to delineate all variants accessible in 3-step trajectories (length of historical trajectory)
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

#second function: output all genotypes that in principle could be reached in three steps from a genotype, if all steps were allowed (including nonfunctional)
get.all.3muts <- function(AAseq){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  mutant1 <- get.Hamming1.nt(seq)
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    mutant2 <- c(mutant2, get.Hamming1.nt(i))
  }
  mutant2 <- unique(mutant2)
  mutant3 <- vector(mode="character")
  for(i in mutant2){
    mutant3 <- c(mutant3, get.Hamming1.nt(i))
  }
  mutant3 <- unique(mutant3)
  return(unique(c(mutant1, mutant2, mutant3)))
}

egka.3mut.paths.11P <- get.3mut.paths.11P("EGKA")
labels.unordered.11P <- unique(c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]],egka.3mut.paths.11P[[4]]))

egka.3mut.ERE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(ERE.specific.11P$AAseq)]
egka.3mut.prom.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(promiscuous.11P$AAseq)]
egka.3mut.SRE.11P <- labels.unordered.11P[labels.unordered.11P %in% as.character(SRE.specific.11P$AAseq)]

#output gephi file for making "wheel" of 3mut paths from EGKA
nodes <- data.frame(id=1:(length(egka.3mut.ERE.11P)+length(egka.3mut.prom.11P)+length(egka.3mut.SRE.11P)),label=c(as.character(egka.3mut.ERE.11P),as.character(egka.3mut.prom.11P),as.character(egka.3mut.SRE.11P)),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  if(nodes[i,"label"] %in% c(egka.3mut.paths.11P[[1]],egka.3mut.paths.11P[[2]],egka.3mut.paths.11P[[3]])){
    if(nodes[i,"label"] %in% as.character(ERE.specific.11P$AAseq)){
      for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
        source <- c(source, i)
        target <- c(target, j)
      }
    }else if(nodes[i,"label"] %in% as.character(promiscuous.11P$AAseq)){
      for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]) & (nodes$label %in% as.character(promiscuous.11P$AAseq) | nodes$label %in% as.character(SRE.specific.11P$AAseq)))){
        source <- c(source, i)
        target <- c(target, j)
      }
    }else if(nodes[i,"label"] %in% as.character(SRE.specific.11P$AAseq)){
      for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]) & nodes$label %in% as.character(SRE.specific.11P$AAseq))){
        source <- c(source, i)
        target <- c(target, j)
      }
    }
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=c(rep(100,length(egka.3mut.ERE.11P)),rep(col2rgb("cyan")[1],length(egka.3mut.prom.11P)),rep(5,length(egka.3mut.SRE.11P))),g=c(rep(69,length(egka.3mut.ERE.11P)),rep(col2rgb("cyan")[2],length(egka.3mut.prom.11P)),rep(172,length(egka.3mut.SRE.11P))),b=c(rep(155,length(egka.3mut.ERE.11P)),rep(col2rgb("cyan")[3],length(egka.3mut.prom.11P)),rep(72,length(egka.3mut.SRE.11P))),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
sizes <- sapply(1:nrow(nodes), function(x) if(nodes$label[x] %in% egka.3mut.paths.11P[[1]]){return(50)}else if(nodes$label[x] %in% egka.3mut.paths.11P[[2]]){return(40)}else if(nodes$label[x] %in% egka.3mut.paths.11P[[3]]){return(20)}else if(nodes$label[x] %in% egka.3mut.paths.11P[[4]]){return(10)})
write.gexf(nodes,edges,nodesVizAtt=list(color=colors, size=sizes),output="7_assess-connectivity_out/gephi-files/EGKA-3-step-wheel_11P.gexf")

#to illustrate that GSKV is densely connected to other SREsp variants in the 11P space, give wheel of single and double mutant neighbors
get.Hamming2.nt <- function(AAseq){
  seq <- as.character(AAseq)
  mutant1 <- get.Hamming1.nt(seq)
  mutant2 <- vector(mode="character")
  for(i in 1:length(mutant1)){
    mutant2 <- unique(c(mutant2, get.Hamming1.nt(mutant1[i])))
  }
  return(unique(c(mutant1,mutant2)))
}

Hamming2.nt.gskv <- get.Hamming2.nt("GSKV")

ERE.specific.gskv <- ERE.specific.11P[AAseq %in% Hamming2.nt.gskv,]
SRE.specific.gskv <- SRE.specific.11P[AAseq %in% Hamming2.nt.gskv,]
promiscuous.gskv <- promiscuous.11P[AAseq %in% Hamming2.nt.gskv,]

#make node and edge data frames to convert to gexf format for visualization in gephi, for all variants that are SRE-specific and 1 or 2 hamming nt from GSKV
nodes <- data.frame(id=1:length(SRE.specific.gskv$AAseq),label=as.character(SRE.specific.gskv$AAseq),stringsAsFactors = F)
source <- vector()
target <- vector()
for(i in 1:nrow(nodes)){
  for(j in which(nodes$label %in% get.Hamming1.nt(nodes$label[i]))){
    source <- c(source, i)
    target <- c(target, j)
  }
}
edges <- data.frame(source=source,target=target,stringsAsFactors = F)
colors <- data.frame(r=rep(col2rgb(rgb(5,172,72,maxColorValue=255))[1],nrow(SRE.specific.gskv)),g=rep(col2rgb(rgb(5,172,72,maxColorValue=255))[2],nrow(SRE.specific.gskv)),b=rep(col2rgb(rgb(5,172,72,maxColorValue=255))[3],nrow(SRE.specific.gskv)),alpha=rep(1.0,nrow(nodes)),stringsAsFactors=F)
write.gexf(nodes,edges,nodesVizAtt=list(color=colors),output="7_assess-connectivity_out/gephi-files/gskv-2-neighbors_graph.gexf")

#how many SRE-sp outcomes are accessible in <=3 nt mutational steps from other ERE-sp starting points?
#for each ERE-specific genotype, output the mut1-mut3 genotypes it can reach
ERE.specific.11P[,c("mut0","mut1","mut2","mut3","all.muts","all.SRE.muts"):=NA]
for(i in 1:nrow(ERE.specific.11P)){
  geno <- as.character(ERE.specific.11P$AAseq[i])
  geno.3mut.paths <- get.3mut.paths.11P(geno)
  ERE.specific.11P$mut0[i] <- list(geno.3mut.paths[[1]])
  ERE.specific.11P$mut1[i] <- list(geno.3mut.paths[[2]])
  ERE.specific.11P$mut2[i] <- list(geno.3mut.paths[[3]])
  ERE.specific.11P$mut3[i] <- list(geno.3mut.paths[[4]])
  ERE.specific.11P$all.muts[i] <- list(unique(c(geno.3mut.paths[[1]],geno.3mut.paths[[2]],geno.3mut.paths[[3]],geno.3mut.paths[[4]])))
  ERE.specific.11P$all.SRE.muts[i] <- list(ERE.specific.11P$all.muts[i][[1]][which(ERE.specific.11P$all.muts[i][[1]] %in% SRE.specific.11P$AAseq)])
}

for(i in 1:nrow(ERE.specific.11P)){
  geno.outcomes <- ERE.specific.11P[i,all.muts][[1]]
  geno.SRE.outcomes <- ERE.specific.11P[i,all.SRE.muts][[1]]
  ERE.specific.11P$num.outcomes[i] <- length(geno.outcomes)
  ERE.specific.11P$num.outcomes.shared.egka[i] <- sum(geno.outcomes %in% labels.unordered.11P)
  ERE.specific.11P$num.SRE.outcomes[i] <- length(geno.SRE.outcomes)
  ERE.specific.11P$num.SRE.outcomes.shared.egka[i] <- sum(geno.SRE.outcomes %in% egka.3mut.SRE.11P)
  if("GSKV" %in% ERE.specific.11P$all.muts[i][[1]]){
    ERE.specific.11P$reaches.GSKV[i] <- TRUE
  }else{
    ERE.specific.11P$reaches.GSKV[i] <- FALSE
  }
}

#output genotypes and outcomes in principle accessible from each ERE-sp starting point if no constriant for functional intermediates
for(i in 1:nrow(ERE.specific.11P)){
  print(i)
  all.3muts <- get.all.3muts(ERE.specific.11P$AAseq[i])
  ERE.specific.11P$num.3mut.SREs[i] <- sum(all.3muts %in% as.character(SRE.specific.11P$AAseq))
  ERE.specific.11P$all.3muts[i] <- list(all.3muts)
  ERE.specific.11P$SRE.3muts[i] <- list(all.3muts[all.3muts %in% as.character(SRE.specific.11P$AAseq)])
}

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point.pdf",4,4)
hist(ERE.specific.11P$num.SRE.outcomes,breaks=seq(-5,170,5),xlab="number of SRE-specific outcomes accessed in <= 3-nt steps",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P[num.SRE.outcomes==0,])),breaks=seq(-5,170,5),col="black",add=T) #label black any ERE not accessing any outcomes; after addition of white bar, the remaining black will be ones that could in principle access an SRE but did not because direct paths blocked by nonfunc intemrediate (epistasis)
hist(rep(0,nrow(ERE.specific.11P[num.SRE.outcomes==0 & num.3mut.SREs==0,])),breaks=seq(-5,170,5),col="white",add=T) #white is ERE-variants that didn't access SRE-outcomes because more than 3 muts away from SRE to begin with, even if able to pass through nonfunc intermediates
abline(v=length(egka.3mut.SRE.11P))
dev.off()


#how many different ERE-sp starting points can access each SRE-sp outcome in 3 or fewer steps?
for(i in 1:nrow(SRE.specific.11P)){
  geno <- as.character(SRE.specific.11P[i,AAseq])
  count <- 0
  for(j in 1:nrow(ERE.specific.11P)){
    if(geno %in% ERE.specific.11P[j,all.muts][[1]]){
      count <- count+1
    }
  }
  SRE.specific.11P$accessed.by[i] <- count #how many ERE-sp starting poitns access each outcome through allowed trajectories
  SRE.specific.11P$num.3mut.EREs[i] <- sum(get.all.3muts(SRE.specific.11P$AAseq[i]) %in% as.character(ERE.specific.11P$AAseq)) #how many ERE-sp starting points are within three mutations barring functional constrain, and in principle could access?
}


pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome.pdf",4,4)
hist(SRE.specific.11P$accessed.by,breaks=seq(-1,43,1),xlab="number of ERE-specific starting poitns accessing in <= 3-nt steps",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P[accessed.by==0,])),breaks=seq(-1,43,1),col="black",add=T) #label 0 bar black; after overlaying white, remaining black will represent SRE-sp outcomes that could in principle be reached by at least one ERE-sp starting point, but the direct path(s) were blocked by nonfunctional intermediates (epistasis)
hist(rep(0,nrow(SRE.specific.11P[accessed.by==0 & num.3mut.EREs==0,])),breaks=seq(-1,43,1),col="white",add=T) #color white those variants that are not reached, and could never be reached even if passing through nonfunctional intermediates because they're more than 3 raw steps away from any EREsp starting points
abline(v=SRE.specific.11P["GSKV",accessed.by])
dev.off()


# what fraction of outcomes are shared between each pair of ERE-sp variants in 11P space?
ERE.sp.overlap <- matrix(data=rep(0,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq);colnames(ERE.sp.overlap) <- as.character(ERE.specific.11P$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap)){
  geno.i <- rownames(ERE.sp.overlap)[i]
  for(j in 1:ncol(ERE.sp.overlap)){
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

#see many variants have zero overlap:
sum(ERE.sp.overlap==0,na.rm=T)
#three reasons there could be zero overlap:
#1 because the starting points began >6 mut steps away from each other to begin with, so there's no chance of them converging on a genotype outcome in 3 steps each (n.possible=0: there are no shared genotypes in the 3-mut neighborhood from get.all.3muts)
#2 because the startign poitns began <= 6 mut steps away from each other to begin with (n.possible > 0), but there are no SRE-sp variants in that potential overlap (n.SRE.possible = 0)
#3 starting points had at least one SRE in their joint 3 mut potential neighborhood, but one or both took more than 3 steps to get there (epistasis): n.SRE.possible > 0
n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap)){
  for(j in 1:nrow(ERE.sp.overlap)){
    if(!is.na(ERE.sp.overlap[i,j]) & ERE.sp.overlap[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P[i,all.3muts][[1]] %in% ERE.specific.11P[j,all.3muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P[i,SRE.3muts][[1]] %in% ERE.specific.11P[j,SRE.3muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}
n.possible #4381 pairs of startign points are 6 or fewer nonsyn muts away from one another, so feasibly could converge on an outcome in 3 steps
n.impossible #523 pairs of staritng points are >6 nonsyn muts away from one another, could not converge on an outcome
n.SRE.possible #1218 pairs of starting points share at least one possible SRE outcomes in each of their 3 nonsyn mut radii (that is, compared to n.possible, many starting points are within striking distance from one another, but there are no SRE.specific outcomes in this overlap in sequence space)
n.SRE.impossible #3686 pairs of starting points do not share at least one SRE outcome in their shared 3 nonsyn mut radii


pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap.pdf",width=4,height=4)
hist(ERE.sp.overlap,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants)
dev.off()

#beyond observing that different starting points access different sets of outcomes, I want to address
#whether the state frequencies among outcomes accessed by different starting points differ
#to do this, measure the Jensen-Shannon distance at each variable position between all sets of outcomes accessed
#by each pair of starting points; to compare to a null distribution, sample the numbers of SRE-specific
#variants found from each starting point randomly from all SRE-specific outcomes, compute JS distance between bootstraps

H <- function(v){
  v <- v[v>0]
  return(sum(-v*log(v)))
}
#function that takes the dataframe entry of the list of reached outcomes, calculates proportion per site of each AA
JS.dist <- function(query1,query2){
  if(length(query1)<16 | length(query2)<16){
    return(list(NA,NA,NA,NA))
  }else{
    data1 <- data.frame(row.names=query1)
    data2 <- data.frame(row.names=query2)
    for(i in 1:length(query1)){
      data1$AA1[i] <- strsplit(query1[i],split="")[[1]][1]
      data1$AA2[i] <- strsplit(query1[i],split="")[[1]][2]
      data1$AA3[i] <- strsplit(query1[i],split="")[[1]][3]
      data1$AA4[i] <- strsplit(query1[i],split="")[[1]][4]
    }
    for(i in 1:length(query2)){
      data2$AA1[i] <- strsplit(query2[i],split="")[[1]][1]
      data2$AA2[i] <- strsplit(query2[i],split="")[[1]][2]
      data2$AA3[i] <- strsplit(query2[i],split="")[[1]][3]
      data2$AA4[i] <- strsplit(query2[i],split="")[[1]][4]
    }
    prop.site1 <- data.frame(site=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20)), AA=rep(AAs,4))
    prop.site2 <- data.frame(site=c(rep(1,20),rep(2,20),rep(3,20),rep(4,20)), AA=rep(AAs,4))
    for(i in 1:nrow(prop.site1)){
      prop.site1$prop[i] <- sum(data1[,prop.site1$site[i]]==as.character(prop.site1$AA[i]))/length(data1[,prop.site1$site[i]])
    }
    for(i in 1:nrow(prop.site2)){
      prop.site2$prop[i] <- sum(data2[,prop.site2$site[i]]==as.character(prop.site2$AA[i]))/length(data2[,prop.site2$site[i]])
    }
    JSD.1 <- sqrt(H((prop.site1[1:20,3]+prop.site2[1:20,3])/2)-(H(prop.site1[1:20,3])+H(prop.site2[1:20,3]))/2)
    JSD.2 <- sqrt(H((prop.site1[21:40,3]+prop.site2[21:40,3])/2)-(H(prop.site1[21:40,3])+H(prop.site2[21:40,3]))/2)
    JSD.3 <- sqrt(H((prop.site1[41:60,3]+prop.site2[41:60,3])/2)-(H(prop.site1[41:60,3])+H(prop.site2[41:60,3]))/2)
    JSD.4 <- sqrt(H((prop.site1[61:80,3]+prop.site2[61:80,3])/2)-(H(prop.site1[61:80,3])+H(prop.site2[61:80,3]))/2)
    return(list(JSD.1,JSD.2,JSD.3,JSD.4))
  }
}

#calculate JS distances for all pairs of starting points that each access at least 16 outcomes (top 50%)
JSdists1 <- matrix(data=rep(NA,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(JSdists1) <- as.character(ERE.specific.11P$AAseq);colnames(JSdists1) <- as.character(ERE.specific.11P$AAseq)
JSdists2 <- matrix(data=rep(NA,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(JSdists2) <- as.character(ERE.specific.11P$AAseq);colnames(JSdists2) <- as.character(ERE.specific.11P$AAseq)
JSdists3 <- matrix(data=rep(NA,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(JSdists3) <- as.character(ERE.specific.11P$AAseq);colnames(JSdists3) <- as.character(ERE.specific.11P$AAseq)
JSdists4 <- matrix(data=rep(NA,nrow(ERE.specific.11P)*nrow(ERE.specific.11P)),nrow=nrow(ERE.specific.11P)); rownames(JSdists4) <- as.character(ERE.specific.11P$AAseq);colnames(JSdists4) <- as.character(ERE.specific.11P$AAseq)
for(i in 1:nrow(JSdists1)){
  for(j in 1:ncol(JSdists1)){
    if(i<j){
      dists <- JS.dist( ERE.specific.11P[rownames(JSdists1)[i],all.SRE.muts[[1]]] , ERE.specific.11P[colnames(JSdists1)[j],all.SRE.muts[[1]]] )
      JSdists1[i,j] <- dists[[1]]
      JSdists2[i,j] <- dists[[2]]
      JSdists3[i,j] <- dists[[3]]
      JSdists4[i,j] <- dists[[4]]
    }
    if(i==j){
      JSdists1[i,j]<-NA
      JSdists2[i,j]<-NA
      JSdists3[i,j]<-NA
      JSdists4[i,j]<-NA
    }
    if(j<i){
      dists.random <- JS.dist(sample(as.character(SRE.specific.11P$AAseq),length(ERE.specific.11P[rownames(JSdists1)[i],all.SRE.muts[[1]]]),replace=F),sample(as.character(SRE.specific.11P$AAseq),length(ERE.specific.11P[rownames(JSdists1)[j],all.SRE.muts[[1]]]),replace=F))
      JSdists1[i,j] <- dists.random[[1]]
      JSdists2[i,j] <- dists.random[[2]]
      JSdists3[i,j] <- dists.random[[3]]
      JSdists4[i,j] <- dists.random[[4]]
    }
  }
}


pdf(file="7_assess-connectivity_out/distributions_Jensen-Shannon-distance-of-outcomes-between-pairs-of-ERE.pdf",7,7.5,useDingbats=F)
par(mfrow=c(2,2))
plot(density(JSdists1[upper.tri(JSdists1)][!is.na(JSdists1[upper.tri(JSdists1)])]),xlim=c(0,1),xlab="Jensen-Shannon distance",main="site 25",ylim=c(0,5))
polygon(density(JSdists1[upper.tri(JSdists1)][!is.na(JSdists1[upper.tri(JSdists1)])]),col=rgb(0,0,255,200,maxColorValue=255))
lines(density(JSdists1[lower.tri(JSdists1)][!is.na(JSdists1[lower.tri(JSdists1)])]))
polygon(density(JSdists1[lower.tri(JSdists1)][!is.na(JSdists1[lower.tri(JSdists1)])]),col=rgb(255,0,0,200,maxColorValue=255))
legend(legend=c("observed outcomes","randomly sampled outcomes"),"topright", pch=15,col=c("blue","red"),bty="n",cex=0.7)

plot(density(JSdists2[upper.tri(JSdists2)][!is.na(JSdists2[upper.tri(JSdists2)])]),xlim=c(0,1),xlab="Jensen-Shannon distance",main="site 26",ylim=c(0,6))
polygon(density(JSdists2[upper.tri(JSdists2)][!is.na(JSdists2[upper.tri(JSdists2)])]),col=rgb(0,0,255,200,maxColorValue=255))
lines(density(JSdists2[lower.tri(JSdists2)][!is.na(JSdists2[lower.tri(JSdists2)])]))
polygon(density(JSdists2[lower.tri(JSdists2)][!is.na(JSdists2[lower.tri(JSdists2)])]),col=rgb(255,0,0,200,maxColorValue=255))
legend(legend=c("observed outcomes","randomly sampled outcomes"),"topright", pch=15,col=c("blue","red"),bty="n",cex=0.7)

plot(density(JSdists3[upper.tri(JSdists3)][!is.na(JSdists3[upper.tri(JSdists3)])]),xlim=c(0,1),xlab="Jensen-Shannon distance",main="site 28",ylim=c(0,6))
polygon(density(JSdists3[upper.tri(JSdists3)][!is.na(JSdists3[upper.tri(JSdists3)])]),col=rgb(0,0,255,200,maxColorValue=255))
lines(density(JSdists3[lower.tri(JSdists3)][!is.na(JSdists3[lower.tri(JSdists3)])]))
polygon(density(JSdists3[lower.tri(JSdists3)][!is.na(JSdists3[lower.tri(JSdists3)])]),col=rgb(255,0,0,200,maxColorValue=255))
legend(legend=c("observed outcomes","randomly sampled outcomes"),"topright", pch=15,col=c("blue","red"),bty="n",cex=0.7)

plot(density(JSdists4[upper.tri(JSdists4)][!is.na(JSdists4[upper.tri(JSdists4)])]),xlim=c(0,1),xlab="Jensen-Shannon distance",main="site 29",ylim=c(0,5))
polygon(density(JSdists4[upper.tri(JSdists4)][!is.na(JSdists4[upper.tri(JSdists4)])]),col=rgb(0,0,255,200,maxColorValue=255))
lines(density(JSdists4[lower.tri(JSdists4)][!is.na(JSdists4[lower.tri(JSdists4)])]))
polygon(density(JSdists4[lower.tri(JSdists4)][!is.na(JSdists4[lower.tri(JSdists4)])]),col=rgb(255,0,0,200,maxColorValue=255))
legend(legend=c("observed outcomes","randomly sampled outcomes"),"topright", pch=15,col=c("blue","red"),bty="n",cex=0.7)
dev.off()

JSdistsavg <- (JSdists1+JSdists2+JSdists3+JSdists4)/4
pdf(file="7_assess-connectivity_out/distributions_Jensen-Shannon-distance-of-outcomes-between-pairs-of-ERE_average-across-sites.pdf",3,3.5,useDingbats=F)
plot(density(JSdistsavg[upper.tri(JSdistsavg)][!is.na(JSdistsavg[upper.tri(JSdistsavg)])]),xlim=c(0,1),ylim=c(0,6.5),xlab="Jensen-Shannon distance",main="average of 4 sites")
polygon(density(JSdistsavg[upper.tri(JSdistsavg)][!is.na(JSdistsavg[upper.tri(JSdistsavg)])]),col=rgb(0,0,255,200,maxColorValue=255))
lines(density(JSdistsavg[lower.tri(JSdistsavg)][!is.na(JSdistsavg[lower.tri(JSdistsavg)])]))
polygon(density(JSdistsavg[lower.tri(JSdistsavg)][!is.na(JSdistsavg[lower.tri(JSdistsavg)])]),col=rgb(255,0,0,200,maxColorValue=255))
dev.off()

#as easily seen, all differences statistically significant
t.test(JSdists1[upper.tri(JSdists1)][!is.na(JSdists1[upper.tri(JSdists1)])],JSdists1[lower.tri(JSdists1)][!is.na(JSdists1[lower.tri(JSdists1)])])
t.test(JSdists2[upper.tri(JSdists2)][!is.na(JSdists2[upper.tri(JSdists2)])],JSdists2[lower.tri(JSdists2)][!is.na(JSdists2[lower.tri(JSdists2)])])
t.test(JSdists3[upper.tri(JSdists3)][!is.na(JSdists3[upper.tri(JSdists3)])],JSdists3[lower.tri(JSdists3)][!is.na(JSdists3[lower.tri(JSdists3)])])
t.test(JSdists4[upper.tri(JSdists4)][!is.na(JSdists4[upper.tri(JSdists4)])],JSdists4[lower.tri(JSdists4)][!is.na(JSdists4[lower.tri(JSdists4)])])

##################################################################################################################################
#how often are randomly drawn networks of 1351 nodes from all possible RH seqs as or more connected than the observed 11P space?
sum(components(all.pos.11P)$csize[-1]) #nodes isolated from primary subnetwork in actual 11P space
components(all.pos.11P)$no #3 components
set.seed(10101010)
bootstraps.num.components <- vector()
for(k in 1:1000){
  print(k)
  AAseqs <- dt.11P.coding$AAseq[sample(1:160000,nrow(ERE.specific.11P)+nrow(SRE.specific.11P)+nrow(promiscuous.11P))]
  adj.m.bootstrap <- matrix(data=rep(0,length(AAseqs)^2),nrow=length(AAseqs));colnames(adj.m.bootstrap) <- as.character(AAseqs); rownames(adj.m.bootstrap) <- as.character(AAseqs)
  for(i in 1:nrow(adj.m.bootstrap)){
    adj.m.bootstrap[i,colnames(adj.m.bootstrap) %in% get.Hamming1.nt(row.names(adj.m.bootstrap)[i])] <- 1
  }
  graph.bootstrap <- graph.adjacency(adj.m.bootstrap,mode="undirected")
  bootstraps.num.components <- c(bootstraps.num.components,components(graph.bootstrap)$no)
  print(components(graph.bootstrap)$no)
}
save(bootstraps.num.components, file="7_assess-connectivity_out/bootstraps.numcomponents.11P.Rda")
pdf(file="7_assess-connectivity_out/boostrap-networks-num-components_11P.pdf",4,4,useDingbats=F)
hist(bootstraps.num.components,xlab="Number of components",ylab="Bootstrap frequency",xlim=c(0,max(bootstraps.num.components)),main="")
abline(v=components(all.pos.11P)$no,lty=2)
dev.off()
median(bootstraps.num.components) #1181, range 1136-1216, mean 1180.7

#how often are randomly drawn networks of 129 nodes from all possible RH seqs as or more connected than the observed SR1 space?
sum(components(all.pos.SR1)$csize[-3]) #nodes isolated from primary subnetwork in actual SR1 space = 22
components(all.pos.SR1)$no #14 components
set.seed(20202020)
bootstraps.num.components <- vector()
for(k in 1:1000){
  print(k)
  AAseqs <- dt.SR1.coding$AAseq[sample(1:160000,nrow(ERE.specific.SR1)+nrow(SRE.specific.SR1)+nrow(promiscuous.SR1))]
  adj.m.bootstrap <- matrix(data=rep(0,length(AAseqs)^2),nrow=length(AAseqs));colnames(adj.m.bootstrap) <- as.character(AAseqs); rownames(adj.m.bootstrap) <- as.character(AAseqs)
  for(i in 1:nrow(adj.m.bootstrap)){
    adj.m.bootstrap[i,colnames(adj.m.bootstrap) %in% get.Hamming1.nt(row.names(adj.m.bootstrap)[i])] <- 1
  }
  graph.bootstrap <- graph.adjacency(adj.m.bootstrap,mode="undirected")
  bootstraps.num.components <- c(bootstraps.num.components,components(graph.bootstrap)$no)
  print(components(graph.bootstrap)$no)
}
save(bootstraps.num.components, file="7_assess-connectivity_out/bootstraps.numcomponents.SR1.Rda")
pdf(file="7_assess-connectivity_out/boostrap-networks-num-components_SR1.pdf",4,4,useDingbats=F)
hist(bootstraps.num.components,xlab="Number of components",ylab="Bootstrap frequency",xlim=c(0,140),main="")
abline(v=components(all.pos.SR1)$no,lty=2)
dev.off()
median(bootstraps.num.components) #128, range 122-129, mean 127.5, 240 of the 1000 have all 129 isolated with absolutely no connection


##################################################################################################################################
#test how contingency stuff in 11P network breaks down with increasing path lengths
#note: this takes a couple days to run on a laptop. Could probably be sped up with sapply's instead of for loops...
#length 1
ERE.specific.11P.1 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.1 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.1 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

get.1mut.paths.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P.1$AAseq),prom=as.character(promiscuous.11P.1$AAseq),SRE.sp=as.character(SRE.specific.11P.1$AAseq)){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  if(seq %in% ERE.sp){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(ERE.sp,prom,SRE.sp)]
  }else if(seq %in% prom){
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(prom,SRE.sp)]
  }else{
    mutant1 <- get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(SRE.sp)]
  }
  all.muts <- unique(mutant1)
  return(all.muts)
}

get.all.1muts <- function(AAseq){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  mutant1 <- get.Hamming1.nt(seq)
  return(unique(mutant1))
}

egka.1mut.paths.11P <- get.1mut.paths.11P("EGKA")

egka.1mut.SRE.11P <- egka.1mut.paths.11P[egka.1mut.paths.11P %in% as.character(SRE.specific.11P.1$AAseq)]

#how many SRE-sp outcomes are accessible in <=1 nonsyn mutational steps from other ERE-sp starting points, and in 1 steps retaining functional intermediates?
#for each ERE-specific genotype, output the mut1-mut4 genotypes it can reach
ERE.specific.11P.1[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.1muts","SRE.1muts","num.1mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.1)){
  print(i)
  geno <- get.1mut.paths.11P(as.character(ERE.specific.11P.1$AAseq[i]))
  all.1muts <- get.all.1muts(as.character(ERE.specific.11P.1$AAseq[i]))
  ERE.specific.11P.1$all.muts[i] <- list(geno)
  ERE.specific.11P.1$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.1$AAseq])
  ERE.specific.11P.1$num.SRE.outcomes[i] <- length(ERE.specific.11P.1[i,all.SRE.muts][[1]])
  ERE.specific.11P.1$all.1muts[i] <- list(all.1muts)
  ERE.specific.11P.1$SRE.1muts[i] <- list(all.1muts[all.1muts %in% as.character(SRE.specific.11P.1$AAseq)])
  ERE.specific.11P.1$num.1mut.SREs[i] <- sum(all.1muts %in% as.character(SRE.specific.11P.1$AAseq))
}
save(ERE.specific.11P.1, file="7_assess-connectivity_out/ERE-specific.11P.1.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.1.Rda")

for(i in 1:nrow(SRE.specific.11P.1)){
  print(i)
  geno <- as.character(SRE.specific.11P.1[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.1)){
    if(geno %in% ERE.specific.11P.1[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.1[j,all.1muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.1$accessed.by[i] <- count
  SRE.specific.11P.1$num.1mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.1, file="7_assess-connectivity_out/SRE-specific.11P.1.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.1.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_1-step-paths.pdf",4,4)
hist(ERE.specific.11P.1$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 1-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.1[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.1[num.SRE.outcomes==0 & num.1mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.1mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_1-step-paths.pdf",4,4)
hist(SRE.specific.11P.1$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 1-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.1[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.1[accessed.by==0 & num.1mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.1["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.1 <- matrix(data=rep(0,nrow(ERE.specific.11P.1)*nrow(ERE.specific.11P.1)),nrow=nrow(ERE.specific.11P.1)); rownames(ERE.sp.overlap.1) <- as.character(ERE.specific.11P.1$AAseq);colnames(ERE.sp.overlap.1) <- as.character(ERE.specific.11P.1$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.1)){
  geno.i <- rownames(ERE.sp.overlap.1)[i]
  for(j in 1:ncol(ERE.sp.overlap.1)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.1)[j]
    ERE.sp.overlap.1[geno.i,geno.j] <- sum(ERE.specific.11P.1[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.1[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.1[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.1)){
  ERE.sp.overlap.1[i,i] <- NA
}
#upper right diagonal puts overlap as zero when j has no outcomes (whereas lower left gives NA like desired, b/c divide by zero)
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.1[as.character(ERE.specific.11P.1[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.1[,as.character(ERE.specific.11P.1[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.1))
save(ERE.sp.overlap.1, file="7_assess-connectivity_out/ERE.sp.overlap.1.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.1.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.1)){
  for(j in 1:nrow(ERE.sp.overlap.1)){
    if(!is.na(ERE.sp.overlap.1[i,j]) & ERE.sp.overlap.1[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.1[i,all.1muts][[1]] %in% ERE.specific.11P.1[j,all.1muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.1[i,SRE.1muts][[1]] %in% ERE.specific.11P.1[j,SRE.1muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}
n.possible 
n.impossible
n.SRE.possible
n.SRE.impossible

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_1-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.1,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 2
ERE.specific.11P.2 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.2 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.2 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

get.2mut.paths.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P.2$AAseq),prom=as.character(promiscuous.11P.2$AAseq),SRE.sp=as.character(SRE.specific.11P.2$AAseq)){
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
  all.muts <- unique(c(mutant1,mutant2))
  return(all.muts)
}

get.all.2muts <- function(AAseq){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  mutant1 <- get.Hamming1.nt(seq)
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    mutant2 <- c(mutant2, get.Hamming1.nt(i))
  }
  mutant2 <- unique(mutant2)
  return(unique(c(mutant1, mutant2)))
}

egka.2mut.paths.11P <- get.2mut.paths.11P("EGKA")

egka.2mut.SRE.11P <- egka.2mut.paths.11P[egka.2mut.paths.11P %in% as.character(SRE.specific.11P.2$AAseq)]

#how many SRE-sp outcomes are accessible in <=4 nonsyn mutational steps from other ERE-sp starting points, and in 4 steps retaining functional intermediates?
#for each ERE-specific genotype, output the mut1-mut4 genotypes it can reach
ERE.specific.11P.2[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.2muts","SRE.2muts","num.2mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.2)){
  print(i)
  geno <- get.2mut.paths.11P(as.character(ERE.specific.11P.2$AAseq[i]))
  all.2muts <- get.all.2muts(as.character(ERE.specific.11P.2$AAseq[i]))
  ERE.specific.11P.2$all.muts[i] <- list(geno)
  ERE.specific.11P.2$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.2$AAseq])
  ERE.specific.11P.2$num.SRE.outcomes[i] <- length(ERE.specific.11P.2[i,all.SRE.muts][[1]])
  ERE.specific.11P.2$all.2muts[i] <- list(all.2muts)
  ERE.specific.11P.2$SRE.2muts[i] <- list(all.2muts[all.2muts %in% as.character(SRE.specific.11P.2$AAseq)])
  ERE.specific.11P.2$num.2mut.SREs[i] <- sum(all.2muts %in% as.character(SRE.specific.11P.2$AAseq))
}
save(ERE.specific.11P.2, file="7_assess-connectivity_out/ERE-specific.11P.2.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.2.Rda")

for(i in 1:nrow(SRE.specific.11P.2)){
  print(i)
  geno <- as.character(SRE.specific.11P.2[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.2)){
    if(geno %in% ERE.specific.11P.2[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.2[j,all.2muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.2$accessed.by[i] <- count
  SRE.specific.11P.2$num.2mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.2, file="7_assess-connectivity_out/SRE-specific.11P.2.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.2.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_2-step-paths.pdf",4,4)
hist(ERE.specific.11P.2$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 2-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.2[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.2[num.SRE.outcomes==0 & num.2mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.2mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_2-step-paths.pdf",4,4)
hist(SRE.specific.11P.2$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 2-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.2[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.2[accessed.by==0 & num.2mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.2["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.11P.2 <- matrix(data=rep(0,nrow(ERE.specific.11P.2)*nrow(ERE.specific.11P.2)),nrow=nrow(ERE.specific.11P.2)); rownames(ERE.sp.11P.2) <- as.character(ERE.specific.11P.2$AAseq);colnames(ERE.sp.11P.2) <- as.character(ERE.specific.11P.2$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.11P.2)){
  geno.i <- rownames(ERE.sp.11P.2)[i]
  for(j in 1:ncol(ERE.sp.11P.2)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.11P.2)[j]
    ERE.sp.11P.2[geno.i,geno.j] <- sum(ERE.specific.11P.2[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.2[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.2[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.11P.2)){
  ERE.sp.11P.2[i,i] <- NA
}
#upper right diagonal puts overlap as zero when j has no outcomes (whereas lower left gives NA like desired, b/c divide by zero)
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.11P.2[as.character(ERE.specific.11P.2[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.11P.2[,as.character(ERE.specific.11P.2[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.11P.2))
save(ERE.sp.11P.2, file="7_assess-connectivity_out/ERE.sp.11P.2.Rda")
load(file="7_assess-connectivity_out/ERE.sp.11P.2.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.11P.2)){
  for(j in 1:nrow(ERE.sp.11P.2)){
    if(!is.na(ERE.sp.11P.2[i,j]) & ERE.sp.11P.2[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.2[i,all.2muts][[1]] %in% ERE.specific.11P.2[j,all.2muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.2[i,SRE.2muts][[1]] %in% ERE.specific.11P.2[j,SRE.2muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_2-step-paths.pdf",width=4,height=4)
hist(ERE.sp.11P.2,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 4
ERE.specific.11P.4 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.4 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.4 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

get.4mut.paths.11P <- function(AAseq,ERE.sp=as.character(ERE.specific.11P.4$AAseq),prom=as.character(promiscuous.11P.4$AAseq),SRE.sp=as.character(SRE.specific.11P.4$AAseq)){
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
  mutant4 <- vector(mode="character")
  for(i in mutant3){
    if(i %in% ERE.sp){
      mutant4 <- c(mutant4, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(ERE.sp,prom,SRE.sp)])
    }else if(i %in% prom){
      mutant4 <- c(mutant4, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(prom,SRE.sp)])
    }else{
      mutant4 <- c(mutant4, get.Hamming1.nt(i)[get.Hamming1.nt(i) %in% c(SRE.sp)])
    }
  }
  all.muts <- unique(c(mutant1,mutant2,mutant3,mutant4))
  return(all.muts)
}

get.all.4muts <- function(AAseq){
  seq <- as.character(AAseq)
  mutant1 <- vector(mode="character")
  mutant1 <- get.Hamming1.nt(seq)
  mutant2 <- vector(mode="character")
  for(i in mutant1){
    mutant2 <- c(mutant2, get.Hamming1.nt(i))
  }
  mutant2 <- unique(mutant2)
  mutant3 <- vector(mode="character")
  for(i in mutant2){
    mutant3 <- c(mutant3, get.Hamming1.nt(i))
  }
  mutant3 <- unique(mutant3)
  mutant4 <- vector(mode="character")
  for(i in mutant3){
    mutant4 <- c(mutant4, get.Hamming1.nt(i))
  }
  mutant4 <- unique(mutant4)
  return(unique(c(mutant1, mutant2, mutant3, mutant4)))
}

egka.4mut.paths.11P <- get.4mut.paths.11P("EGKA")

egka.4mut.SRE.11P <- egka.4mut.paths.11P[egka.4mut.paths.11P %in% as.character(SRE.specific.11P.4$AAseq)]

#how many SRE-sp outcomes are accessible in <=4 nonsyn mutational steps from other ERE-sp starting points, and in 4 steps retaining functional intermediates?
#for each ERE-specific genotype, output the mut1-mut4 genotypes it can reach
ERE.specific.11P.4[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.4muts","SRE.4muts","num.4mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.4)){
  print(i)
  geno <- get.4mut.paths.11P(as.character(ERE.specific.11P.4$AAseq[i]))
  all.4muts <- get.all.4muts(as.character(ERE.specific.11P.4$AAseq[i]))
  ERE.specific.11P.4$all.muts[i] <- list(geno)
  ERE.specific.11P.4$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.4$AAseq])
  ERE.specific.11P.4$num.SRE.outcomes[i] <- length(ERE.specific.11P.4[i,all.SRE.muts][[1]])
  ERE.specific.11P.4$all.4muts[i] <- list(all.4muts)
  ERE.specific.11P.4$SRE.4muts[i] <- list(all.4muts[all.4muts %in% as.character(SRE.specific.11P.4$AAseq)])
  ERE.specific.11P.4$num.4mut.SREs[i] <- sum(all.4muts %in% as.character(SRE.specific.11P.4$AAseq))
}
save(ERE.specific.11P.4, file="7_assess-connectivity_out/ERE-specific.11P.4.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.4.Rda")

for(i in 1:nrow(SRE.specific.11P.4)){
  print(i)
  geno <- as.character(SRE.specific.11P.4[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.4)){
    if(geno %in% ERE.specific.11P.4[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.4[j,all.4muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.4$accessed.by[i] <- count
  SRE.specific.11P.4$num.4mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.4, file="7_assess-connectivity_out/SRE-specific.11P.4.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.4.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_4-step-paths.pdf",4,4)
hist(ERE.specific.11P.4$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 4-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.4[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.4[num.SRE.outcomes==0 & num.4mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.4mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_4-step-paths.pdf",4,4)
hist(SRE.specific.11P.4$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 4-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.4[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.4[accessed.by==0 & num.4mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.4["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.4 <- matrix(data=rep(0,nrow(ERE.specific.11P.4)*nrow(ERE.specific.11P.4)),nrow=nrow(ERE.specific.11P.4)); rownames(ERE.sp.overlap.4) <- as.character(ERE.specific.11P.4$AAseq);colnames(ERE.sp.overlap.4) <- as.character(ERE.specific.11P.4$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.4)){
  geno.i <- rownames(ERE.sp.overlap.4)[i]
  for(j in 1:ncol(ERE.sp.overlap.4)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.4)[j]
    ERE.sp.overlap.4[geno.i,geno.j] <- sum(ERE.specific.11P.4[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.4[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.4[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.4)){
  ERE.sp.overlap.4[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.4[as.character(ERE.specific.11P.4[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.4[,as.character(ERE.specific.11P.4[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.4))
save(ERE.sp.overlap.4, file="7_assess-connectivity_out/ERE.sp.overlap.4.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.4.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.4)){
  for(j in 1:nrow(ERE.sp.overlap.4)){
    if(!is.na(ERE.sp.overlap.4[i,j]) & ERE.sp.overlap.4[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.4[i,all.4muts][[1]] %in% ERE.specific.11P.4[j,all.4muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.4[i,SRE.4muts][[1]] %in% ERE.specific.11P.4[j,SRE.4muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_4-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.4,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()


########################################################################################################################
#length 5
ERE.specific.11P.5 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.5 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.5 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

#define function that takes a vector of AAseqs, returns all additional next functional steps
get.next.func.step <- function(AAseqs,ERE.sp=as.character(ERE.specific.11P.5$AAseq),prom=as.character(promiscuous.11P.5$AAseq),SRE.sp=as.character(SRE.specific.11P.5$AAseq)){
  outcomes <- vector(mode="character")
  for(seq in AAseqs){
    if(seq %in% ERE.sp){
      outcomes <- c(outcomes,get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(ERE.sp,prom,SRE.sp)])
    }else if(seq %in% prom){
      outcomes <- c(outcomes,get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(prom,SRE.sp)])
    }else{
      outcomes <- c(outcomes,get.Hamming1.nt(seq)[get.Hamming1.nt(seq) %in% c(SRE.sp)])
    }
  }
  outcomes <- unique(c(AAseqs, outcomes))
  return(outcomes)
}

#define function that takes a vector of AAseqs, returns all aditional neighbors one more nt mutation away
get.next.step <- function(AAseq){
  outcomes <- vector(mode="character")
  for(seq in AAseq){
    outcomes <- unique(c(outcomes, get.Hamming1.nt(seq)))
  }
  outcomes <- unique(c(AAseq, outcomes))
  return(outcomes)
}

#for each ERE-specific genotype, output the next set of genotypes it can reach given its 4mut trajectories
ERE.specific.11P.5[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.5muts","SRE.5muts","num.5mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.5)){
  print(i)
  geno <- get.next.func.step(ERE.specific.11P.4[i,all.muts][[1]])
  all.5muts <- get.next.step(ERE.specific.11P.4[i,all.4muts][[1]])
  ERE.specific.11P.5$all.muts[i] <- list(geno)
  ERE.specific.11P.5$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P$AAseq])
  ERE.specific.11P.5$num.SRE.outcomes[i] <- length(ERE.specific.11P.5[i,all.SRE.muts][[1]])
  ERE.specific.11P.5$all.5muts[i] <- list(all.5muts)
  ERE.specific.11P.5$SRE.5muts[i] <- list(all.5muts[all.5muts %in% as.character(SRE.specific.11P.5$AAseq)])
  ERE.specific.11P.5$num.5mut.SREs[i] <- sum(all.5muts %in% as.character(SRE.specific.11P.5$AAseq))
}
save(ERE.specific.11P.5, file="7_assess-connectivity_out/ERE-specific.11P.5.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.5.Rda")

for(i in 1:nrow(SRE.specific.11P.5)){
  print(i)
  geno <- as.character(SRE.specific.11P.5[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.5)){
    if(geno %in% ERE.specific.11P.5[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.5[j,all.5muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.5$accessed.by[i] <- count
  SRE.specific.11P.5$num.5mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.5, file="7_assess-connectivity_out/SRE-specific.11P.5.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.5.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_5-step-paths.pdf",4,4)
hist(ERE.specific.11P.5$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 5-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.5[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.5[num.SRE.outcomes==0 & num.5mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.5mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_5-step-paths.pdf",4,4)
hist(SRE.specific.11P.5$accessed.by,breaks=seq(-5,5*round(max(SRE.specific.11P.5$accessed.by)/5)+5,5),xlab="number of ERE-specific starting poitns accessing in <= 5-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.5[accessed.by==0,])),breaks=seq(-5,5*round(max(SRE.specific.11P.5$accessed.by)/5)+5,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.5[accessed.by==0 & num.5mut.EREs==0,])),breaks=seq(-5,5*round(max(SRE.specific.11P.5$accessed.by)/5)+5,5),col="white",add=T)
#abline(v=SRE.specific.11P.5["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.5 <- matrix(data=rep(0,nrow(ERE.specific.11P.5)*nrow(ERE.specific.11P.5)),nrow=nrow(ERE.specific.11P.5)); rownames(ERE.sp.overlap.5) <- as.character(ERE.specific.11P.5$AAseq);colnames(ERE.sp.overlap.5) <- as.character(ERE.specific.11P.5$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.5)){
  geno.i <- rownames(ERE.sp.overlap.5)[i]
  for(j in 1:ncol(ERE.sp.overlap.5)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.5)[j]
    ERE.sp.overlap.5[geno.i,geno.j] <- sum(ERE.specific.11P.5[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.5[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.5[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.5)){
  ERE.sp.overlap.5[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.5[as.character(ERE.specific.11P.5[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.5[,as.character(ERE.specific.11P.5[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.5))
save(ERE.sp.overlap.5, file="7_assess-connectivity_out/ERE.sp.overlap.5.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.5.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.5)){
  for(j in 1:nrow(ERE.sp.overlap.5)){
    if(!is.na(ERE.sp.overlap.5[i,j]) & ERE.sp.overlap.5[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.5[i,all.5muts][[1]] %in% ERE.specific.11P.5[j,all.5muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.5[i,SRE.5muts][[1]] %in% ERE.specific.11P.5[j,SRE.5muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_5-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.5,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 6
ERE.specific.11P.6 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.6 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.6 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]


#for each ERE-specific genotype, output the mut1-mut6 genotypes it can reach
ERE.specific.11P.6[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.6muts","SRE.6muts","num.6mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.6)){
  print(i)
  geno <- get.next.func.step(ERE.specific.11P.5[i,all.muts][[1]])
  all.6muts <- get.next.step(ERE.specific.11P.5[i,all.5muts][[1]])
  ERE.specific.11P.6$all.muts[i] <- list(geno)
  ERE.specific.11P.6$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.6$AAseq])
  ERE.specific.11P.6$num.SRE.outcomes[i] <- length(ERE.specific.11P.6[i,all.SRE.muts][[1]])
  ERE.specific.11P.6$all.6muts[i] <- list(all.6muts)
  ERE.specific.11P.6$SRE.6muts[i] <- list(all.6muts[all.6muts %in% as.character(SRE.specific.11P.6$AAseq)])
  ERE.specific.11P.6$num.6mut.SREs[i] <- sum(all.6muts %in% as.character(SRE.specific.11P.6$AAseq))
}
save(ERE.specific.11P.6, file="7_assess-connectivity_out/ERE-specific.11P.6.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.6.Rda")

for(i in 1:nrow(SRE.specific.11P.6)){
  print(i)
  geno <- as.character(SRE.specific.11P.6[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.6)){
    if(geno %in% ERE.specific.11P.6[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.6[j,all.6muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.6$accessed.by[i] <- count
  SRE.specific.11P.6$num.6mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.6, file="7_assess-connectivity_out/SRE-specific.11P.6.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.6.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_6-step-paths.pdf",4,4)
hist(ERE.specific.11P.6$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 6-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.6[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.6[num.SRE.outcomes==0 & num.6mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.6mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_6-step-paths.pdf",4,4)
hist(SRE.specific.11P.6$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 6-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.6[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.6[accessed.by==0 & num.6mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.6["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.6 <- matrix(data=rep(0,nrow(ERE.specific.11P.6)*nrow(ERE.specific.11P.6)),nrow=nrow(ERE.specific.11P.6)); rownames(ERE.sp.overlap.6) <- as.character(ERE.specific.11P.6$AAseq);colnames(ERE.sp.overlap.6) <- as.character(ERE.specific.11P.6$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.6)){
  geno.i <- rownames(ERE.sp.overlap.6)[i]
  for(j in 1:ncol(ERE.sp.overlap.6)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.6)[j]
    ERE.sp.overlap.6[geno.i,geno.j] <- sum(ERE.specific.11P.6[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.6[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.6[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.6)){
  ERE.sp.overlap.6[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.6[as.character(ERE.specific.11P.6[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.6[,as.character(ERE.specific.11P.6[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.6))
save(ERE.sp.overlap.6, file="7_assess-connectivity_out/ERE.sp.overlap.6.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.6.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.6)){
  for(j in 1:nrow(ERE.sp.overlap.6)){
    if(!is.na(ERE.sp.overlap.6[i,j]) & ERE.sp.overlap.6[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.6[i,all.6muts][[1]] %in% ERE.specific.11P.6[j,all.6muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.6[i,SRE.6muts][[1]] %in% ERE.specific.11P.6[j,SRE.6muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_6-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.6,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 7
ERE.specific.11P.7 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.7 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.7 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

#for each ERE-specific genotype, output the mut1-mut7 genotypes it can reach
ERE.specific.11P.7[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.7muts","SRE.7muts","num.7mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.7)){
  print(i)
  last.geno <- ERE.specific.11P.6[i,all.muts][[1]][!(ERE.specific.11P.6[i,all.muts][[1]] %in% ERE.specific.11P.5[i,all.muts][[1]])] #output only those variants that were just reached in 6 steps -- no need to get hamming 1 of variants that are already accessed in 1st round, accounted for in round 2
  last.mut <- ERE.specific.11P.6[i,all.6muts][[1]][!(ERE.specific.11P.6[i,all.6muts][[1]] %in% ERE.specific.11P.5[i,all.5muts][[1]])] #output only those genotypes that were just reached in 6 steps -- no need to get hamming 1 of genotypes that are already accessed in 1st round, accounted for in round 2
  geno <- unique(c(get.next.func.step(last.geno),ERE.specific.11P.6[i,all.muts][[1]])) # get all new next steps, add to list of all previously reached steps
  all.7muts <- unique(c(get.next.step(last.mut),ERE.specific.11P.6[i,all.6muts][[1]]))
  ERE.specific.11P.7$all.muts[i] <- list(geno)
  ERE.specific.11P.7$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.7$AAseq])
  ERE.specific.11P.7$num.SRE.outcomes[i] <- length(ERE.specific.11P.7[i,all.SRE.muts][[1]])
  ERE.specific.11P.7$all.7muts[i] <- list(all.7muts)
  ERE.specific.11P.7$SRE.7muts[i] <- list(all.7muts[all.7muts %in% as.character(SRE.specific.11P.7$AAseq)])
  ERE.specific.11P.7$num.7mut.SREs[i] <- sum(all.7muts %in% as.character(SRE.specific.11P.7$AAseq))
}
save(ERE.specific.11P.7, file="7_assess-connectivity_out/ERE-specific.11P.7.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.7.Rda")

for(i in 1:nrow(SRE.specific.11P.7)){
  print(i)
  geno <- as.character(SRE.specific.11P.7[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.7)){
    if(geno %in% ERE.specific.11P.7[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.7[j,all.7muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.7$accessed.by[i] <- count
  SRE.specific.11P.7$num.7mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.7, file="7_assess-connectivity_out/SRE-specific.11P.7.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.7.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_7-step-paths.pdf",4,4)
hist(ERE.specific.11P.7$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 7-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.7[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.7[num.SRE.outcomes==0 & num.7mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.7mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_7-step-paths.pdf",4,4)
hist(SRE.specific.11P.7$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 7-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.7[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.7[accessed.by==0 & num.7mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.7["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.7 <- matrix(data=rep(0,nrow(ERE.specific.11P.7)*nrow(ERE.specific.11P.7)),nrow=nrow(ERE.specific.11P.7)); rownames(ERE.sp.overlap.7) <- as.character(ERE.specific.11P.7$AAseq);colnames(ERE.sp.overlap.7) <- as.character(ERE.specific.11P.7$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.7)){
  geno.i <- rownames(ERE.sp.overlap.7)[i]
  for(j in 1:ncol(ERE.sp.overlap.7)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.7)[j]
    ERE.sp.overlap.7[geno.i,geno.j] <- sum(ERE.specific.11P.7[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.7[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.7[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.7)){
  ERE.sp.overlap.7[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.7[as.character(ERE.specific.11P.7[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.7[,as.character(ERE.specific.11P.7[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.7))
save(ERE.sp.overlap.7, file="7_assess-connectivity_out/ERE.sp.overlap.7.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.7.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.7)){
  for(j in 1:nrow(ERE.sp.overlap.7)){
    if(!is.na(ERE.sp.overlap.7[i,j]) & ERE.sp.overlap.7[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.7[i,all.7muts][[1]] %in% ERE.specific.11P.7[j,all.7muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.7[i,SRE.7muts][[1]] %in% ERE.specific.11P.7[j,SRE.7muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_7-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.7,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 7 or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 8
ERE.specific.11P.8 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.8 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.8 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

#for each ERE-specific genotype, output the 1mut-8mut genotypes it can reach
ERE.specific.11P.8[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.8muts","SRE.8muts","num.8mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.8)){
  print(i)
  last.geno <- ERE.specific.11P.7[i,all.muts][[1]][!(ERE.specific.11P.7[i,all.muts][[1]] %in% ERE.specific.11P.6[i,all.muts][[1]])] #output only those variants that were just reached in 7 steps -- no need to get hamming 1 of variants that are already accessed in 1st round, accounted for in round 2
  last.mut <- ERE.specific.11P.7[i,all.7muts][[1]][!(ERE.specific.11P.7[i,all.7muts][[1]] %in% ERE.specific.11P.6[i,all.6muts][[1]])] #output only those genotypes that were just reached in 6 steps -- no need to get hamming 1 of genotypes that are already accessed in 1st round, accounted for in round 2
  geno <- unique(c(get.next.func.step(last.geno),ERE.specific.11P.7[i,all.muts][[1]])) # get all new next steps, add to list of all previously reached steps
  all.8muts <- unique(c(get.next.step(last.mut),ERE.specific.11P.7[i,all.7muts][[1]]))
  ERE.specific.11P.8$all.muts[i] <- list(geno)
  ERE.specific.11P.8$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.8$AAseq])
  ERE.specific.11P.8$num.SRE.outcomes[i] <- length(ERE.specific.11P.8[i,all.SRE.muts][[1]])
  ERE.specific.11P.8$all.8muts[i] <- list(all.8muts)
  ERE.specific.11P.8$SRE.8muts[i] <- list(all.8muts[all.8muts %in% as.character(SRE.specific.11P.8$AAseq)])
  ERE.specific.11P.8$num.8mut.SREs[i] <- sum(all.8muts %in% as.character(SRE.specific.11P.8$AAseq))
}
save(ERE.specific.11P.8, file="7_assess-connectivity_out/ERE-specific.11P.8.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.8.Rda")

for(i in 1:nrow(SRE.specific.11P.8)){
  print(i)
  geno <- as.character(SRE.specific.11P.8[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.8)){
    if(geno %in% ERE.specific.11P.8[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.8[j,all.8muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.8$accessed.by[i] <- count
  SRE.specific.11P.8$num.8mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.8, file="7_assess-connectivity_out/SRE-specific.11P.8.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.8.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_8-step-paths.pdf",4,4)
hist(ERE.specific.11P.8$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 8-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.8[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.8[num.SRE.outcomes==0 & num.8mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.8mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_8-step-paths.pdf",4,4)
hist(SRE.specific.11P.8$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 8-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.8[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.8[accessed.by==0 & num.8mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.8["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.8 <- matrix(data=rep(0,nrow(ERE.specific.11P.8)*nrow(ERE.specific.11P.8)),nrow=nrow(ERE.specific.11P.8)); rownames(ERE.sp.overlap.8) <- as.character(ERE.specific.11P.8$AAseq);colnames(ERE.sp.overlap.8) <- as.character(ERE.specific.11P.8$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.8)){
  geno.i <- rownames(ERE.sp.overlap.8)[i]
  for(j in 1:ncol(ERE.sp.overlap.8)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.8)[j]
    ERE.sp.overlap.8[geno.i,geno.j] <- sum(ERE.specific.11P.8[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.8[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.8[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.8)){
  ERE.sp.overlap.8[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.8[as.character(ERE.specific.11P.8[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.8[,as.character(ERE.specific.11P.8[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.8))
save(ERE.sp.overlap.8, file="7_assess-connectivity_out/ERE.sp.overlap.8.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.8.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.8)){
  for(j in 1:nrow(ERE.sp.overlap.8)){
    if(!is.na(ERE.sp.overlap.8[i,j]) & ERE.sp.overlap.8[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.8[i,all.8muts][[1]] %in% ERE.specific.11P.8[j,all.8muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.8[i,SRE.8muts][[1]] %in% ERE.specific.11P.8[j,SRE.8muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_8-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.8,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 8-nt or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 9
ERE.specific.11P.9 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.9 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.9 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

#for each ERE-specific genotype, output the 1mut-9mut genotypes it can reach
ERE.specific.11P.9[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.9muts","SRE.9muts","num.9mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.9)){
  print(i)
  last.geno <- ERE.specific.11P.8[i,all.muts][[1]][!(ERE.specific.11P.8[i,all.muts][[1]] %in% ERE.specific.11P.7[i,all.muts][[1]])] #output only those variants that were just reached in 8 steps -- no need to get hamming 1 of variants that are already accessed in 1st round, accounted for in round 2
  last.mut <- ERE.specific.11P.8[i,all.8muts][[1]][!(ERE.specific.11P.8[i,all.8muts][[1]] %in% ERE.specific.11P.7[i,all.7muts][[1]])] #output only those genotypes that were just reached in 6 steps -- no need to get hamming 1 of genotypes that are already accessed in 1st round, accounted for in round 2
  geno <- unique(c(get.next.func.step(last.geno),ERE.specific.11P.8[i,all.muts][[1]])) # get all new next steps, add to list of all previously reached steps
  if(length(ERE.specific.11P.8[i,all.8muts][[1]]) < 160000){
    all.9muts <- unique(c(get.next.step(last.mut),ERE.specific.11P.8[i,all.8muts][[1]]))
  }else{
    all.9muts <- ERE.specific.11P.8[i,all.8muts][[1]]
  }
  ERE.specific.11P.9$all.muts[i] <- list(geno)
  ERE.specific.11P.9$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.9$AAseq])
  ERE.specific.11P.9$num.SRE.outcomes[i] <- length(ERE.specific.11P.9[i,all.SRE.muts][[1]])
  ERE.specific.11P.9$all.9muts[i] <- list(all.9muts)
  ERE.specific.11P.9$SRE.9muts[i] <- list(all.9muts[all.9muts %in% as.character(SRE.specific.11P.9$AAseq)])
  ERE.specific.11P.9$num.9mut.SREs[i] <- sum(all.9muts %in% as.character(SRE.specific.11P.9$AAseq))
}
save(ERE.specific.11P.9, file="7_assess-connectivity_out/ERE-specific.11P.9.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.9.Rda")

for(i in 1:nrow(SRE.specific.11P.9)){
  print(i)
  geno <- as.character(SRE.specific.11P.9[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.9)){
    if(geno %in% ERE.specific.11P.9[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.9[j,all.9muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.9$accessed.by[i] <- count
  SRE.specific.11P.9$num.9mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.9, file="7_assess-connectivity_out/SRE-specific.11P.9.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.9.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_9-step-paths.pdf",4,4)
hist(ERE.specific.11P.9$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 9-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.9[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.9[num.SRE.outcomes==0 & num.9mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.9mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_9-step-paths.pdf",4,4)
hist(SRE.specific.11P.9$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 9-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.9[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.9[accessed.by==0 & num.9mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.9["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.9 <- matrix(data=rep(0,nrow(ERE.specific.11P.9)*nrow(ERE.specific.11P.9)),nrow=nrow(ERE.specific.11P.9)); rownames(ERE.sp.overlap.9) <- as.character(ERE.specific.11P.9$AAseq);colnames(ERE.sp.overlap.9) <- as.character(ERE.specific.11P.9$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.9)){
  geno.i <- rownames(ERE.sp.overlap.9)[i]
  for(j in 1:ncol(ERE.sp.overlap.9)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.9)[j]
    ERE.sp.overlap.9[geno.i,geno.j] <- sum(ERE.specific.11P.9[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.9[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.9[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.9)){
  ERE.sp.overlap.9[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.9[as.character(ERE.specific.11P.9[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.9[,as.character(ERE.specific.11P.9[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.9))
save(ERE.sp.overlap.9, file="7_assess-connectivity_out/ERE.sp.overlap.9.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.9.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.9)){
  for(j in 1:nrow(ERE.sp.overlap.9)){
    if(!is.na(ERE.sp.overlap.9[i,j]) & ERE.sp.overlap.9[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.9[i,all.9muts][[1]] %in% ERE.specific.11P.9[j,all.9muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.9[i,SRE.9muts][[1]] %in% ERE.specific.11P.9[j,SRE.9muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_9-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.9,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 9-nt or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()

########################################################################################################################
#length 10
ERE.specific.11P.10 <- ERE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt,shortest.path.length.to.SRE)]
SRE.specific.11P.10 <- SRE.specific.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]
promiscuous.11P.10 <- promiscuous.11P[,.(AAseq,n.ERE.nt,n.prom.nt,n.SRE.nt)]

#for each ERE-specific genotype, output the 1mut-10mut genotypes it can reach
ERE.specific.11P.10[,c("all.muts","all.SRE.muts","num.SRE.outcomes","all.10muts","SRE.10muts","num.10mut.SREs"):=NA]
for(i in 1:nrow(ERE.specific.11P.10)){
  print(i)
  last.geno <- ERE.specific.11P.9[i,all.muts][[1]][!(ERE.specific.11P.9[i,all.muts][[1]] %in% ERE.specific.11P.8[i,all.muts][[1]])] #output only those variants that were just reached in 9 steps -- no need to get hamming 1 of variants that are already accessed in 1st round, accounted for in round 2
  last.mut <- ERE.specific.11P.9[i,all.9muts][[1]][!(ERE.specific.11P.9[i,all.9muts][[1]] %in% ERE.specific.11P.8[i,all.8muts][[1]])] #output only those genotypes that were just reached in 6 steps -- no need to get hamming 1 of genotypes that are already accessed in 1st round, accounted for in round 2
  geno <- unique(c(get.next.func.step(last.geno),ERE.specific.11P.9[i,all.muts][[1]])) # get all new next steps, add to list of all previously reached steps
  if(length(ERE.specific.11P.9[i,all.9muts][[1]]) < 160000){
    all.10muts <- unique(c(get.next.step(last.mut),ERE.specific.11P.9[i,all.9muts][[1]]))
  }else{
    all.10muts <- ERE.specific.11P.9[i,all.9muts][[1]]
  }
  ERE.specific.11P.10$all.muts[i] <- list(geno)
  ERE.specific.11P.10$all.SRE.muts[i] <- list(geno[geno %in% SRE.specific.11P.10$AAseq])
  ERE.specific.11P.10$num.SRE.outcomes[i] <- length(ERE.specific.11P.10[i,all.SRE.muts][[1]])
  ERE.specific.11P.10$all.10muts[i] <- list(all.10muts)
  ERE.specific.11P.10$SRE.10muts[i] <- list(all.10muts[all.10muts %in% as.character(SRE.specific.11P.10$AAseq)])
  ERE.specific.11P.10$num.10mut.SREs[i] <- sum(all.10muts %in% as.character(SRE.specific.11P.10$AAseq))
}
save(ERE.specific.11P.10, file="7_assess-connectivity_out/ERE-specific.11P.10.Rda")
load(file="7_assess-connectivity_out/ERE-specific.11P.10.Rda")

for(i in 1:nrow(SRE.specific.11P.10)){
  print(i)
  geno <- as.character(SRE.specific.11P.10[i,AAseq])
  count <- 0
  count.possible <- 0
  for(j in 1:nrow(ERE.specific.11P.10)){
    if(geno %in% ERE.specific.11P.10[j,all.muts][[1]]){
      count <- count+1
    }
    if(geno %in% ERE.specific.11P.10[j,all.10muts][[1]]){
      count.possible <- count.possible+1
    }
  }
  SRE.specific.11P.10$accessed.by[i] <- count
  SRE.specific.11P.10$num.10mut.EREs[i] <- count.possible
}
save(SRE.specific.11P.10, file="7_assess-connectivity_out/SRE-specific.11P.10.Rda")
load(file="7_assess-connectivity_out/SRE-specific.11P.10.Rda")

pdf(file="./7_assess-connectivity_out/hist_11P_number-SRE-outcomes-accessible-per-ERE-starting-point_10-step-paths.pdf",4,4)
hist(ERE.specific.11P.10$num.SRE.outcomes,breaks=seq(-20,840,20),xlab="number of SRE-specific outcomes accessed in <= 10-nt mutations",ylab="number ERE-specific starting points",main="",col="gray75")
hist(rep(0,nrow(ERE.specific.11P.10[num.SRE.outcomes==0,])),breaks=seq(-20,840,20),col="black",add=T)
hist(rep(0,nrow(ERE.specific.11P.10[num.SRE.outcomes==0 & num.10mut.SREs==0,])),breaks=seq(-20,840,20),col="white",add=T)
#abline(v=length(egka.10mut.SRE.11P))
dev.off()

pdf(file="./7_assess-connectivity_out/hist_11P_number-ERE-starting-points-accessing-per-SRE-outcome_10-step-paths.pdf",4,4)
hist(SRE.specific.11P.10$accessed.by,breaks=seq(-5,145,5),xlab="number of ERE-specific starting poitns accessing in <= 10-nt mutations",ylab="number SRE-specific outcomes",main="",col="gray75")
hist(rep(0,nrow(SRE.specific.11P.10[accessed.by==0,])),breaks=seq(-5,145,5),col="black",add=T)
hist(rep(0,nrow(SRE.specific.11P.10[accessed.by==0 & num.10mut.EREs==0,])),breaks=seq(-5,145,5),col="white",add=T)
#abline(v=SRE.specific.11P.10["GSKV",accessed.by])
dev.off()

# what fraction of outcomes are shared between each pair of ERE-sp variants?
ERE.sp.overlap.10 <- matrix(data=rep(0,nrow(ERE.specific.11P.10)*nrow(ERE.specific.11P.10)),nrow=nrow(ERE.specific.11P.10)); rownames(ERE.sp.overlap.10) <- as.character(ERE.specific.11P.10$AAseq);colnames(ERE.sp.overlap.10) <- as.character(ERE.specific.11P.10$AAseq)
#for each cell [i,j], give the fraction of SRE-sp outcomes found from i also found by j
for(i in 1:nrow(ERE.sp.overlap.10)){
  geno.i <- rownames(ERE.sp.overlap.10)[i]
  for(j in 1:ncol(ERE.sp.overlap.10)){
    print(c(i,j))
    geno.j <- colnames(ERE.sp.overlap.10)[j]
    ERE.sp.overlap.10[geno.i,geno.j] <- sum(ERE.specific.11P.10[geno.i,all.SRE.muts][[1]] %in% ERE.specific.11P.10[geno.j,all.SRE.muts][[1]])/length(ERE.specific.11P.10[geno.i,all.SRE.muts][[1]])
  }
}
#take out diagonals
for(i in 1:nrow(ERE.sp.overlap.10)){
  ERE.sp.overlap.10[i,i] <- NA
}
#make NA any combo i,j that includes an EREsp variant that accesses no SRE
ERE.sp.overlap.10[as.character(ERE.specific.11P.10[num.SRE.outcomes==0,AAseq]),] <- NA
ERE.sp.overlap.10[,as.character(ERE.specific.11P.10[num.SRE.outcomes==0,AAseq])] <- NA
sum(!is.na(ERE.sp.overlap.10))
save(ERE.sp.overlap.10, file="7_assess-connectivity_out/ERE.sp.overlap.10.Rda")
load(file="7_assess-connectivity_out/ERE.sp.overlap.10.Rda")

n.possible <- 0
n.impossible <- 0
n.SRE.possible <- 0
n.SRE.impossible <- 0
for(i in 1:nrow(ERE.sp.overlap.10)){
  for(j in 1:nrow(ERE.sp.overlap.10)){
    if(!is.na(ERE.sp.overlap.10[i,j]) & ERE.sp.overlap.10[i,j] == 0){
      if(i < j){
        print(c(i,j))
        if(sum(ERE.specific.11P.10[i,all.10muts][[1]] %in% ERE.specific.11P.10[j,all.10muts][[1]])>1){
          n.possible <- n.possible + 1
        }else{
          n.impossible <- n.impossible + 1
        }
        if(sum(ERE.specific.11P.10[i,SRE.10muts][[1]] %in% ERE.specific.11P.10[j,SRE.10muts][[1]])>1){
          n.SRE.possible <- n.SRE.possible+1
        }else{
          n.SRE.impossible <- n.SRE.impossible+1
        }
      }
    }
  }
}

pdf(file="7_assess-connectivity_out/hist_11P_ERE-outcomes-overlap_10-step-paths.pdf",width=4,height=4)
hist(ERE.sp.overlap.10,breaks=seq(-0.02,1,0.02),col="gray75",xlab="proportion of SRE-specific outcomes accessible\nfrom i also accessible from j",main="")
hist(rep(0,(n.possible+n.impossible)*2),breaks=seq(-0.02,1,0.02),col="black",add=T) #black bar represents pairs that had at least one SRE outcome in their joint 3mut neighborhoods, but one or both did not access each of these possible outcomes
hist(rep(0,(n.SRE.impossible)*2),breaks=seq(-0.02,1,0.02),col="gray40",add=T) #dark gray bar represents pairs that have shared genotypes in their joint 3mut neighborhoods, but no SRE outcomes are in this shared space so no chance to converge on an outcome
hist(rep(0,(n.impossible)*2),breaks=seq(-0.02,1,0.02),col="white",add=T) #white bar represents pairs that do not share any genotypes in their joint 3mut neighborhoods (the genos are 10-nt or more steps away in raw sequence space without considering funcitonal variants); these are a subset of the dark gray bar, because they are the subset of those without a joint accessible SRE that also don't have ANYTHING jointly accessible
dev.off()


