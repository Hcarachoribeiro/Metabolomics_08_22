##################################################################################################################################
# Multivariate analysis of the whole metabolomics dataset
# For Henrique Caracho Ribeiro
# Partho Sen, created March 03,2020
# Edited on 19.03.2020: PS
# Script for Proteome Analysis #
#################################################################################################################################
rm(list=ls())

library(mixOmics)
library(readr)

print('Change accordingly, the $PATH to the script folder in your computer :-')
setwd("") # Your working directory

# Loading variables
DF <- read.table("../Data/ProteomicsBDxHC_multiblock_VIP.csv", header =T, sep = ',', stringsAsFactors = FALSE)
Protein <- DF[,1]; DF = data.matrix(DF[,-1]); rownames(DF) <- Protein;

Samp.nom = colnames(DF);

DF = DF[,-which(colnames(DF) %in% 'CN3')];
Samp.nom = Samp.nom[-which(Samp.nom %in% 'CN3')];

Grps = c( rep('Bipolar',1,length(grep('B',colnames(DF)))),
          rep('HC',1,length(grep('C',colnames(DF)))) );

print('Assigning groups to the column of the data.frame (DF)')

colnames(DF) <- Grps;  

## log2 transform and autoscale the data ##
print('log2 transform')
DF = log2(DF);
DF[is.infinite(DF)] = 0;
DF[is.na(DF)] = 0;


## Auto-scale the data if required! ## Not recommended unless needed!

############################### Multivariate dimensional reduction tecniques ########################################
print('Loading scripts/functions and necessary dependencies for PCA analysis...')
source('ellipseplot.R')
source('PCAvarAxis.R')
source('range.R')

PCs <- prcomp(t(DF), scale. = TRUE, center = TRUE); ## Already centered before ##

pdf(paste0("../Results/Proteomics_Plots/PCA_proteomics_no_QCs_knownOnly.pdf"), width = 14, height = 12, paper = 'special')
oripar=par()$mar
par(xpd=F, mar=oripar+c(2,2,0,6))

explainPCAvar <- PCAvarAxis(PCs)

ellipseplot(PCs$x[,1],                              # data for x-axis
            PCs$x[,2],                              # data for y-axis
            as.factor(Grps),                        # factor with classes
            pcol= c('#3795E3','#E3378A'), # colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
            pbgcol=TRUE,                           # point borders black?
            cexsize=6,                              # size of points
            ppch=rep(21,2),                         # shape of points (must match # of factors)
            legpos="topright",                      # position of legend          
            legcexsize=2,                          # legend text size
            legptsize=3,                           # legend point size
            axissize=2,                            # Set axis text size
            linewidth=2,
            Samp.nom = Samp.nom,
            elev = 0.95,
            legs= as.factor(c('Bipolar','HC')))


title(xlab=explainPCAvar$PC1,    # % variance explained on PC1
      ylab=explainPCAvar$PC2,    # % variance explained on PC2
      main='PCA',                  # Title
      cex.lab=2,                    # size of label text
      cex.main=2                    # size of title text
)
dev.off()

####### Next: Do PLS-DA or sparse-PLSDA ################################

if(TRUE){
  B = t(DF); rownames(B) <- NULL;

  Grps[Grps=='Bipolar'] =1;
  Grps[Grps=='HC'] = 0;
  Grps = as.numeric(Grps);
  
  plsda <- plsda(B,Grps, ncomp = 3);

  Scores = plsda$variates$X;  rownames(Scores) = Grps;
  loads =  plsda$loadings$X
}

####### Model diagnostics & performance ################################

perf.plsda <- perf(plsda, validation = c("Mfold"), folds = 5, progressBar = TRUE, auc = TRUE, nrepeat=100)

print(perf.plsda$error.rate)
print(perf.plsda$auc)


###### Comment: Q2 and R2 are poor :: better go for PLSDA with AUC estimates 
if(0){
  Y.mat <- unmap(Grps) # creates a dummy matrix
  res <- plsda(B,Y.mat, ncomp = 3)
  val <- perf(res, criterion = c("R2", "Q2"), validation = c("Mfold"), folds = 5, progressBar = TRUE, auc = TRUE, nrepeat=100)
  print(sum(val$R2));
  print(sum(val$Q2.total));
}
############################################################################
## % of expained variance ##
if(1){
  print('Estimate the % of explained variance for each PLS component')
  Rd.YvsU = cor(as.numeric(as.factor(Grps)), plsda$variates$X[, 1:3])
  Rd.YvsU = apply(Rd.YvsU^2, 2, sum)
  Rd.Y = cbind(Rd.YvsU, cumsum(Rd.YvsU))
  colnames(Rd.Y) = c("Proportion", "Cumulative")
  
}
############################################################################

pdf(paste0("../Results/Proteomics_Plots/PLSDA_proteomics_no_QCs_knownOnly.pdf"), width = 14, height = 12, paper = 'special')
oripar=par()$mar
par(xpd=F, mar=oripar+c(2,2,0,6))

ellipseplot(Scores[,1],                              # data for x-axis
            Scores[,2],                              # data for y-axis
            as.factor(Grps),                        # factor with classes
            pcol= c('#E3378A','#3795E3'), # colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
            pbgcol=TRUE,                           # point borders black?
            cexsize=6,                              # size of points
            ppch=rep(21,2),                         # shape of points (must match # of factors)
            legpos="topright",                      # position of legend          
            legcexsize=2,                          # legend text size
            legptsize=3,                           # legend point size
            axissize=2,                            # Set axis text size
            linewidth=2,
            Samp.nom = Samp.nom,
            elev = 0.95,
            legs= as.factor(c('HC','Bipolar')))


title(xlab=paste0('PLS1 (', round(Rd.Y[1,1]*100,1)[1], '%)'),    # % variance explained on PLS1
      ylab=paste0('PLS2 (', round(Rd.Y[2,1]*100,1)[1], '%)'),    # % variance explained on PLS2
      main= paste0('PLSDA ','[AUC = ',round(mean(c(perf.plsda$auc$comp1[1],perf.plsda$auc$comp2[1])),2), ']'),                  # Title
      cex.lab=2,                    # size of label text
      cex.main=2                    # size of title text
      
)

dev.off()


################### Next: Variable Selection ################################

print('Getting the VIPs and Regression coefficients from model objects')
RCV = plsda$mat.c
VIP = vip(plsda);

### getting all variables with more than 1 VIPs ####
vip1 <- sort(VIP[,1]); vip1  <- vip1[vip1>=1];
rcv <- sort(RCV[,1]); rcv  <- rcv[vip1>=1];

library("viridis")

cols   = (colorRampPalette(c("steelblue3","steelblue2", "steelblue1","lightskyblue", "royalblue3"))(n = length(vip1)));
cols = rev(viridis(30:50));

pdf(paste0("../Results/Proteomics_Plots/PLSDA_proteomics_VIP_variable_selection.pdf"), width = 10, height = 8, paper = 'special')
oripar=par()$mar
par(xpd=F, mar=oripar+c(2,15,0,2));
bh = barplot(vip1, horiz = T, col = cols , las=2, xlim=c(0,2.5), width = 0.1, space = 0.3, xlab = 'VIPs', main = 'PLSDA:VIPs', las=1, cex.names = 1, cex.axis = 1, border = 'black' );
legend('bottomright', legend = paste0('PLSDA ','[AUC = ', round(perf.plsda$auc$comp1[1],2),']')) ## note this AUC is from 1st component ony ##
box(lty = 1, lwd=2, col = 'black')
dev.off()


################### save the object for multiblock analysis ###################################

Prots = DF;
print("Saving R object ... ")
save(Prots,file="../Data/Prots.rda")
Samp.nom.prots = Samp.nom;
save(Samp.nom.prots,file="../Data/Samp.nom.prots.rda")



########################### The END ############################################################
