##################################################################################################################################
# Partial correlation analysis: first order
# Partho Sen
# Created: 17.08.2020
# Edited on 03.11.2020
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



edge.weights <- function(community, network, weight.within = 100, weight.between = 1) {
  bridges <- crossing(communities = community, graph = network)
  weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
  return(weights) 
}


######## ! ****** Your dataset ***** ! ############

setwd() #Your working directory

print('Loading samples.....')
load(file="../Data/metab.multiblock.rda")
load(file="../Data/prot.multiblock.rda")
load(file="../Data/Samp.nom.multiblock.rda")
load(file="../Data/VIPs.metab.multiblock.rda"); ## First component only
load(file="../Data/VIPs.prot.multiblock.rda"); ## First component only

print('Reading clinical data .....')
CD.DF <- read.csv("../Data/Demographic_data_ABNT.csv",header = T, sep = ','); 
CD.DF <- CD.DF[,-ncol(CD.DF)]; rownames(CD.DF) <- CD.DF[,1]; CD.DF = CD.DF[,-1]
CD.DF <- CD.DF[,-c(3:4)];

print('Selecting particular/significant metabolites & proteins ...')

Mets.DF <- Mets.DF [,-1]; Mets.DF <- t(Mets.DF);
Mets.DF = metab.multiblock;



Prots.DF = prot.multiblock;

### select significant metabolites ###

mets.sel = VIPs.metab.multiblock
mets.sel <- mets.sel[abs(mets.sel)>0.15]; Mets.DF = Mets.DF[ ,which(colnames(Mets.DF) %in% names(mets.sel))]

### select significant proteins ###

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

CN.CD = colnames(CD.DF) ; 
CN.Mets = colnames(Mets.DF); colnames(Mets.DF) = paste0('M',1:length(CN.Mets));
CN.Prots = colnames(Prots.DF); colnames(Prots.DF) = paste0('P',1:length(CN.Prots));
CNs = c(CN.CD,CN.Mets,CN.Prots);

stop()

Groups <- c(rep("Clinical_Info",length(CN.CD)),rep("Metabolomics",length(CN.Mets)),rep("Proteomics",length(CN.Prots)))

print('Make sure they are from same subjects')
print(cbind(rownames(CD.DF), rownames(Mets.DF), rownames(Prots.DF) ))

print('If required we can also do the variable (significant Mets & Prots) selection in the Column of each matrix before merging ....')
Mat = cbind(CD.DF,Mets.DF, Prots.DF);

print('Omit rows with NAs')
Mat <- na.omit(Mat)


print(dim(Mat))
print(head(Mat))


############## Do correlation ######################################################
R1 <- rcorr.adjust(data.matrix(Mat), type = 'spearman',use=c("complete.obs"))
corMat <- R1$R$r 
p.Mat <- R1$R$P; ## Adjusted p-values
corMat[!is.finite(corMat)] <- 0

########## Perform partial correlation analysis ####################

print("Performing partial correlation analysis ...")

Map = cbind(colnames(Mat),Groups,CNs); ## Species and groups 

if(1){
  print('Estimating non-rejection rates')
  nRR = data.matrix(qpNrr(corMat, q=1, I=NULL, restrict.Q=NULL, fix.Q=NULL, nTests=1000,
                                alpha=0.01, pairup.i=NULL, pairup.j=NULL,
                                long.dim.are.variables=TRUE, verbose=TRUE, identicalQs=FALSE,
                                exact.test=TRUE, use=c("complete.obs", "em"), tol=0.1,
                                R.code.only=FALSE, clusterSize=1, estimateTime=FALSE,
                                nAdj2estimateTime=10))
  
  save(nRR,file='../Results/Integratomics_Plots/nRR.rda')
  ## NOTE: nRR should be high to remove spurious correlations ##
  
}



load(file='../Results/Integratomics_Plots/nRR.rda')

print('Filter based on non-rejection rates: Spurious correlation: 1st order')

col <- brewer.pal(n = 11, name = "RdYlBu")
col <- colorRampPalette(col);


corMat[abs(corMat)<0.50] <- 0; ### remove this filter to see positive correlations

diag(corMat) <- 0;


## Change your path here :-

pdf("../Results/Integratomics_Plots/PartialCorrNetwork_.pdf", width = 30, height = 30, paper = 'special')
par(xpd=T,mar=c(2,1,2,4))

### Graphical Visualization and filter ###

g1 <- graph_from_adjacency_matrix(corMat, weighted=T, mode="undirected", diag = F)

E(g1)$color <- ifelse(E(g1)$weight > 0,'steelblue1','maroon') #You had this as "V3 > 0.6" which I guess works but it is more readable as 0. that way if you decide to lower the correlation threshold you do not have to change this line too.

coul <- brewer.pal(nlevels(as.factor(Groups)), "Set2")
coul <- adjustcolor(coul,alpha.f = 0.9)

# Map the color to cylinders
my_color <- coul[as.numeric(as.factor(Groups))]

# plot
par(bg="white", mar=c(4,0,0,2))
set.seed(123)

### Grouping several groups together ####

V(g1)$Groups = Groups;
for(i in unique(V(g1)$Groups)) {
  GroupV = which(V(g1)$Groups == i)
  g1 = add_edges(g1, combn(GroupV, 2), attr=list(weight=50))
} 

## Now create a layout based on G_Grouped
set.seed(56789)
      l = layout_with_lgl(g1)

E(g1)$weight[which((E(g1)$weight!=50))]  = E(g1)$weight[which((E(g1)$weight!=50))]*2;

plot(g1, 
     vertex.size=8,
     vertex.color=my_color, 
     vertex.label.cex=2,
     vertex.label.color="black",
     vertex.frame.color="white",
     vertex.shape="circle",
     vertex.label.family="Helvetica",                  # Font family of the label (e.g.?Times?, ?Helvetica?)
     vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
     vertex.label.dist=0,                          # Distance between the label and the vertex
     vertex.label.degree=0,                       # The position of the label in relation to the vertex (use pi)
     #edge.width=2,                                # Edge width, defaults to 1
     edge.arrow.size=1,                            # Arrow size, defaults to 1
     edge.arrow.width=1,                           # Arrow width, defaults to 1
     edge.lty=c("solid"),                          # Line type, could be 0 or ?blank?, 1 or ?solid?, 2 or ?dashed?, 3 or ?dotted?, 4 or ?dotdash?, 5 or ?longdash?, 6 or ?twodash?
     edge.curved=0,
     edge.width=abs(E(g1)$weight)*2, ## adjust thickness
     #rescale=T, 
     asp = 0,
     layout=l
)

# title and legend
legend('topright', 
       legend=paste(levels(as.factor(Groups)), sep=""), 
       col = coul, 
       title = 'Groups',
       pch=15, pt.cex = 3, cex = 1.5,
       text.col="black", horiz = F,border = T)

legend('bottomright',
       legend= paste0(V(g1)$name,': ',Map[,3]),
       col = my_color,
       bty = "n", pch=19, pt.cex = 1.1, cex = 1.5,
       text.col="black", horiz = F,border = T)

legend('topleft', 
       legend=c('positive','negative'), 
       col = c('steelblue1','maroon'), 
       title = 'Partial correlation',
       lty=c(1,1),lwd=c(6,6),pch=c(NA,NA), pt.cex = 3, cex = 1.5,
       text.col="black", horiz = F,border = T)

dev.off()
