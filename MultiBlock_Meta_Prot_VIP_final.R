# #########################################################################################################################
# Multiomics integration of Proteomics and Metabolomics 
# Partho: March 5th, 2020
# Edited on : 30.05.2020 --PS
##########################################################################################################################
rm(list=ls())
setwd() #Your working directory
library(mixOmics)
library(readr)

print('Loading scripts/functions and necessary dependencies for PCA analysis...')
source('ellipseplot.R')
source('PCAvarAxis.R')
source('range.R')

######### ****** Perform multiblock PLS-DAs ******** ################################
print('Loading previously saved data ..... ')
load(file="../Data/Metab.rda")
load(file="../Data/Prots.rda") 
load(file="../Data/Samp.nom.rda")

Map <- read.table("../Data/Sample_mapping.csv", header =F, sep = '\t', stringsAsFactors = FALSE)

# Selecting metab samples
metabsamples = c(3,6,8,11,12,15,16,17,19,20,23,26,27)

metabmultiblock = Metab[ ,metabsamples]
Samp.nom.multiblock = Samp.nom[metabsamples];


metab = t(metabmultiblock); 
prot = t(Prots)


Grps <- rownames(metab);

rownames(metab) <- paste0(rownames(metab),1:length(rownames(metab)) )
rownames(prot) <- paste0(rownames(prot), 1:length(rownames(prot)) )

# extract training data and name each data frame
X <- list(Metabolites = metab, 
          Proteins = prot);

Y =  as.factor(Grps);

print('Fit a multiblock splsda model...')
res <- block.splsda(X,Y,ncomp = 2, scale = T);


####### Model diagnostics & performance ################################
perf.plsda <- perf(res, validation = c("Mfold"), folds = 5, progressBar = TRUE, auc = TRUE, nrepeat=10)
print(perf.plsda$error.rate)
print(perf.plsda$auc)
##############################################################################

###### Score Plots for different Blcks #########################

print('Getting variables...')

variates.Mets = res$variates$Metabolites; ## Score PCA plot ##
variates.Prots = res$variates$Proteins; ## Score PCA plot ##
variates.Response = res$variates$Y; ## Score PCA plot ##

EV = print(res$explained_variance);

AUC = auroc(res);


print('Make score plots...')

pdf(paste0("../Results/Integratomics_Plots/Integ_meta_prot_MB.pdf"), width = 20, height = 14, paper = 'special')
oripar=par()$mar
par(xpd=F, mar=oripar+c(2,10,0,2))
par(mfrow=c(2,2))

ellipseplot(variates.Mets[,1],                              # data for x-axis
            variates.Mets[,2],                              # data for y-axis
            as.factor(Grps),                        # factor with classes
            pcol= c('#E3378A','#3795E3'), # colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
            pbgcol=TRUE,                           # point borders black?
            cexsize=3,                              # size of points
            ppch=rep(21,2),                         # shape of points (must match # of factors)
            legpos="topright",                      # position of legend          
            legcexsize=2,                          # legend text size
            legptsize=2,                           # legend point size
            axissize=2,                            # Set axis text size
            linewidth=2,
            Samp.nom = NULL,
            elev = 0.95,
            legs= as.factor(c('Bipolar','HC')))



title(xlab=paste0('Component1 (', round(EV$Metabolites[1]*100,1)[1], '%)'),    # % variance explained on PLS1
      ylab=paste0('Component2 (', round(EV$Metabolites[2]*100,1)[1], '%)'),    # % variance explained on PLS2
      main= paste0('Metabolome'),                  # Title
      cex.lab=2,                    # size of label text
      cex.main=2                    # size of title text
      
)

legend('bottomright',
       legend= paste0('AUC','=', AUC$Metabolites$comp1[1]),
       col = 'black',
       cex = 1.5,
       border = T)



ellipseplot(variates.Prots[,1],                     # data for x-axis
            variates.Prots[,2],                    # data for y-axis
            as.factor(Grps),                        # factor with classes
            pcol= c('#E3378A','#3795E3'), # colors for plotting (must match # of factors) ; '#F01F02','#04F1F9''#999999','#FFE4C4)
            pbgcol=TRUE,                           # point borders black?
            cexsize=3,                              # size of points
            ppch=rep(21,2),                         # shape of points (must match # of factors)
            legpos="topright",                      # position of legend          
            legcexsize=2,                          # legend text size
            legptsize=2,                           # legend point size
            axissize=2,                            # Set axis text size
            linewidth=2,
            Samp.nom = NULL,
            elev = 0.95,
            legs= as.factor(c('Bipolar','HC')))


title(xlab=paste0('Component1 (', round(EV$Proteins[1]*100,1)[1], '%)'),    # % variance explained on PLS1
      ylab=paste0('Component2 (', round(EV$Proteins[2]*100,1)[1], '%)'),    # % variance explained on PLS2
      main= paste0('Proteome'), # Title
      cex.lab=2,                    # size of label text
      cex.main=2                    # size of title text
      
)


legend('bottomright',
       legend= paste0('AUC','=', '0.97'),
       col = 'black',
       cex = 1.5,
       border = T)


#### ** Main Contributorsn (loadings to maximize covariance between Y & X): just like VIPs ** ####

loading.Mets = res$loadings$Metabolites;
loading.Prots = res$loadings$Proteins;
loading.Response = res$loadings$Y;

p1 <- sort( abs(loading.Mets[,1]) ); p1 = p1[p1>0.15]; p1 = p1[-grep('Unknown',names(p1))];
p2 <- sort( abs(loading.Prots[,1]) ); p2 = p2[p2>0.15];

library(RColorBrewer)


col = colorRampPalette(brewer.pal(n = 5, name = "GnBu"))(n = length(p1))
bh = barplot(p1, horiz = T, col = col , las=2, xlim=c(0,max(p1)+0.03), width = 0.1, space = 0.3, xlab = 'Loadings', main = 'Metabolites', las=1, cex.names = 1.5, cex.axis = 2, cex.main=2, cex.lab=2, border = 'black' );
box(lty = 1, lwd=2, col = 'black')

col = colorRampPalette(brewer.pal(n = 5, name = "OrRd"))(n = length(p2))
bh = barplot(p2, horiz = T, col = col , las=2, xlim=c(0,max(p2)+0.1), width = 0.1, space = 0.3, xlab = 'Loadings', main = 'Proteins', las=1, cex.names = 1.5, cex.axis = 2, cex.main=2, cex.lab=2, border = 'black' );
box(lty = 1, lwd=2, col = 'black')

dev.off()


AUC <- print(auroc(res, roc.block = "protein", roc.comp = 2))

print("Saving R object ... ")

metab.multiblock = metab;
prot.multiblock = prot;
save(metab.multiblock,file="../Data/metab.multiblock.rda")
save(prot.multiblock,file="../Data/prot.multiblock.rda")
save(Samp.nom.multiblock,file="../Data/Samp.nom.multiblock.rda")

VIPs.metab.multiblock = p1;
VIPs.prot.multiblock = p2;

save(VIPs.metab.multiblock,file="../Data/VIPs.metab.multiblock.rda"); ## First component only
save(VIPs.prot.multiblock,file="../Data/VIPs.prot.multiblock.rda"); ## First component only
