##################################################################################################################################
# Bean plots of metabolites
# Partho Sen
# Created: 03.11.2020
# Only for the paired samples 
#################################################################################################################################
rm(list=ls())
library("qgraph")
library('qpgraph')
library('ppcor')
library('Hmisc')
library('corrr')
library('readxl')
library('igraph')
library('RColorBrewer')
library('RcmdrMisc')
library('beanplot')

setwd() # Your Wd


edge.weights <- function(community, network, weight.within = 100, wei7ght.between = 1) {
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}


## ** function ** ##
diffExp <- function(A,B)
{
  pval=list()
  for(i in 1:dim(A)[1]){
    pval[[i]] = t.test(A[i,],B[i,],paired = FALSE, var.equal = FALSE)$p.value; ## unpaired
  }
  pval = unlist(pval);   pval.adj = p.adjust(pval, method = "fdr", n = length(pval))
  return(list(pval = pval, pval.adj = pval.adj));
}

######## ! ****** Your datasets ***** ! ############
print('Loading samples.....')
#Samp.nom.multiblock 
# load(file="../Data/metab.multiblock.rda")
load("../Data/Metab.rda") #for all samples
# load(file="../Data/prot.multiblock.rda")

# load(file="../Data/Samp.nom.multiblock.rda")
load(file="../Data/Samp.nom.rda") #for all samples
load(file="../Data/VIPs.metab.multiblock.rda"); ## First component only
# load(file="../Data/VIPs.prot.multiblock.rda"); ## First component only

print('Selecting particular/significant metabolites & proteins ...')

stop()

################# Bean Plots of the Metabolites ##################################
# Mets.DF = metab.multiblock;
Mets.DF = t(Metab) #for all samples

HC =  Mets.DF[grep('HC', rownames(Mets.DF)),];
Bipolar =  Mets.DF[grep('Bipolar', rownames(Mets.DF)),]




DF1  = diffExp(t(HC),t(Bipolar));
pvals = round(DF1$pval,3); names(pvals) <- colnames(HC);
nom.sign.pval.ttest = print(pvals[which(pvals<0.05)]);
save(nom.sign.pval.ttest,file="../Data/Snom.sign.pval.ttest_all_samples.rda")#for all samples

####### Make the bean plots ###################################################

A = Bipolar;
B = HC;

pdf(paste0("../Results/Integratomics_Plots/",'Bean_Plots_Significant_Mets_all_samples2',".pdf"), width = 25, height = 15, paper = 'special')
par(mfrow = c (4,4))

for (uu in 1:dim(A)[2]){
  bh = beanplot (A[,uu], B[,uu],
                 main=colnames(A)[uu],col=list('darkgoldenrod2','chartreuse4'),ll=0.1, method="jitter",overallline="mean",
                 beanlines="median",innerborder="black",log="auto",ylab="log2(intensity)",cex=3, what=c(1,1,1,1,0),boxwex =1, xaxt="n",las=2,
                 at=c(1,2), xlim = c(0,3), jitter = T,cex.main=2,cex.lab=1.2,cex.axis=1.5)
  
  legs = paste('p =', pvals[uu]);
  
  legend('topright',legend = legs, cex = 1.5, border = NA, fill = NULL, bty = "n")
  axis(side = 1, at=c(1,2),labels = c('Bipolar','HC'),cex.axis=2,las=1,cex.lab=1.2)
}

dev.off()

################# Bean Plots of the Proteins ##################################
Mets.DF = prot.multiblock;

HC =  Mets.DF[grep('HC', rownames(Mets.DF)),];
Bipolar =  Mets.DF[grep('Bipolar', rownames(Mets.DF)),]

DF1  = diffExp(t(HC),t(Bipolar)); 
pvals = round(DF1$pval,5); names(pvals) <- colnames(HC);
nom.sign.pval.ttest = print(pvals[which(pvals<0.05)]);


####### Make the bean plots ###################################################

A = Bipolar;
B = HC;

pdf(paste0("../Results/Integratomics_Plots/",'Bean_Plots_Significant_Prots',".pdf"), width = 25, height = 15, paper = 'special')
par(mfrow = c (4,4))

for (uu in 1:dim(A)[2]){
  bh = beanplot (A[,uu], B[,uu],
                 main=colnames(A)[uu],col=list('darkgoldenrod2','chartreuse4'),ll=0.1, method="jitter",overallline="mean",
                 beanlines="mean",innerborder="black",log="auto",ylab="log2(intensity)",cex=3, what=c(1,1,1,1,0),boxwex =1, xaxt="n",las=2,
                 at=c(1,2), xlim = c(0,3), jitter = T,cex.main=2,cex.lab=1.2,cex.axis=1.5)
  
  legs = paste('p =', pvals[uu]);
  
  legend('topright',legend = legs, cex = 1.5, border = NA, fill = NULL, bty = "n")
  axis(side = 1, at=c(1,2),labels = c('Bipolar','HC'),cex.axis=2,las=1,cex.lab=1.2)
}

dev.off()

