########################################################################################
#  For Henique
#  Partho: November 4th, 2020
#  Update: 28.09.2021 (Medication types and dosage indicated!)
######################################################################################

rm(list=ls())

# Load the package

library(gplots)
library(ggplot2)
library(nlme) ## gives lme & nlme
library(lme4) ## gives lmer
library(gplots)
library(ggplot2)
library(ade4)
library(corrplot)
library(reshape)
library(Hmisc) ## Also using for imputed fumction
library(plyr)
library(car)
library(maptools)
library(rgeos)
library(plotrix)
library(factoextra)
library(magrittr)
library(robustbase)
library(MASS)
library(scater)
library(SingleCellExperiment)
library(knitr)
options(stringsAsFactors = FALSE)


setwd()

print('Loading samples.....')
load(file="../Data/metab.multiblock.rda")
load(file="../Data/prot.multiblock.rda")
load(file="../Data/Samp.nom.multiblock.rda")
load(file="../Data/VIPs.metab.multiblock.rda"); ## First component only
load(file="../Data/VIPs.prot.multiblock.rda"); ## First component only

print('Reading clinical data .....')
CD.DF <- read.csv("../Data/Demographic_data.csv",header = T, sep = ','); 

Antipsy = matrix(0,1,dim(CD.DF)[1]); Antipsy[grep('Antipsychotic',CD.DF$Medication.type.dummies..remove.)] = 1;
MoodS = matrix(0,1,dim(CD.DF)[1]); MoodS[grep('Mood stabilizer',CD.DF$Medication.type.dummies..remove.)] = 1;

CD.DF = cbind(CD.DF, t(Antipsy), t(MoodS));


rownames(CD.DF) <- CD.DF[,1]; CD.DF = CD.DF[,-1];

print('Removing height, weight and medication names ... as BMI is estimated f(height & weight)....')
CD.DF <- CD.DF[,-c(3,4,6,7)]; 


print('Selecting particular/significant metabolites & proteins ...')
Mets.DF = metab.multiblock;
Prots.DF = prot.multiblock;

mets.sel = VIPs.metab.multiblock
mets.sel <- mets.sel[abs(mets.sel)>0.15]; Mets.DF = Mets.DF[ ,which(colnames(Mets.DF) %in% names(mets.sel))]

prots.sel <- VIPs.prot.multiblock
prots.sel <- prots.sel[abs(prots.sel)>0.15]; Prots.DF = Prots.DF[ ,which(colnames(Prots.DF) %in% names(prots.sel))]

rownames(Mets.DF) <- Samp.nom.multiblock;
rownames(Prots.DF) <- Samp.nom.multiblock;

########################### Combine Matrix for partial correlation analysis ###################################
Map <- read.table("../Data/Sample_mapping.csv", header =F, sep = '\t', stringsAsFactors = FALSE)

idp = match(rownames(CD.DF),Map[,2])
rownames(CD.DF) = Map[idp,1];

Comm = Reduce(intersect, list( rownames(CD.DF),rownames(Mets.DF),rownames(Prots.DF) ) )

idx.Cli = match(Comm, rownames(CD.DF)); print(cbind(Comm, rownames(CD.DF)[idx.Cli])); CD.DF = CD.DF[idx.Cli, ];
idx.Mets = match(Comm, rownames(Mets.DF)); print(cbind(Comm, rownames(Mets.DF)[idx.Mets])); Mets.DF = Mets.DF[idx.Mets, ];
idx.Prots = match(Comm, rownames(Prots.DF)); print(cbind(Comm, rownames(Prots.DF)[idx.Prots])); Prots.DF = Prots.DF[idx.Prots, ];

######### Factor Analysis ##############

### ****** METABOLOMICS ****** #####

umi <- SingleCellExperiment(assays = list(counts = 2^(t(Mets.DF)), logcounts=t(Mets.DF) ), colData = CD.DF )
exprs(umi)

drop_metabs <- apply(exprs(umi), 1, function(x) {var(x) == 0})
umi <- umi[!drop_metabs, ];

vars <- getVarianceExplained(umi)

vars.nom = colnames(vars)
print(vars.nom)

umi <- perCellQCMetrics(umi); ## Calculates QC for each feature



pdf("../Results/Integratomics_Plots/Factor_Analysis_Metabolomics.pdf", width = 12, height = 8, paper = 'special')
oripar=par()$mar
par(xpd=T, mar=oripar+c(2,2,2,2))
par(mfrow=c(1,1))

g=plotExplanatoryVariables(vars, theme_size = 12, nvars_to_plot = length(colnames(vars)), min_marginal_r2 = 0.5)
print(g)

dev.off()


### ****** Proteomics ****** #####

umi <- SingleCellExperiment(assays = list(counts = 2^(t(Prots.DF)), logcounts=t(Prots.DF) ), colData = CD.DF )
exprs(umi)

drop_metabs <- apply(exprs(umi), 1, function(x) {var(x) == 0})
umi <- umi[!drop_metabs, ];
vars <- getVarianceExplained(umi)
vars.nom = colnames(vars)
print(vars.nom)

umi <- perCellQCMetrics(umi); ## Calculates QC for each feature

pdf("../Results/Integratomics_Plots/Factor_Analysis_Proteomics_ver2.pdf", width = 12, height = 8, paper = 'special')
oripar=par()$mar
par(xpd=T, mar=oripar+c(2,2,2,2))
par(mfrow=c(1,1))

g=plotExplanatoryVariables(vars, theme_size = 12, nvars_to_plot = length(colnames(vars)), min_marginal_r2 = 0.5)
print(g)


dev.off()
