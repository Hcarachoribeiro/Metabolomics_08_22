##################################################################################################################################
# Logistic regression
# For Henrique Caracho Ribeiro
# Partho Sen, created March 03,2020
# Edited on 04.11.2020: PS
# Script for Metabolome Analysis #
# After the RSD correction and identification by ALEX using Golm database
#################################################################################################################################
rm(list=ls())


library('nlme') ## gives lme & nlme
library('lme4') ## gives lmer
library('gplots')
library('ggplot2')
library('ade4')
library('corrplot')
library('reshape')
library('Hmisc') ## Also using for imputed fumction
library('plyr')
library('car')
library('maptools')
library('plotrix')
library('ade4')
library('factoextra')
library('magrittr')
library('ropls')
library('beanplot')
library('cluster')
library('ggfortify')
library('MASS')
library('glmnet')
require('glmnet')
library('ROCR')
library('glmnet')
library('caret')
library('plotmo') # for plot_glmnet
library('broom')
library('unbalanced')
library('ROCR')
library('nnet')
library('caret')
library('pROC')
library('klaR')


r1 = adjustcolor( "lightgreen", alpha.f = 0.2)
r2 = adjustcolor( "blue", alpha.f = 0.1)

print('Change accordingly, the $PATH to the script folder in your computer :-')
setwd() # Your WD

load("../Data/Metab.rda")
DF = Metab;
load("../Data/Snom.sign.pval.ttest.rda")
load(file="../Data/Samp.nom.rda")
noms = names(nom.sign.pval.ttest);
noms = noms[-grep("Unknown",noms)];

### Pre-select metabolites based on the univariate t-test and VIPs #####

DF.sig = DF[which(rownames(DF) %in% noms), ]; 
colnames(DF.sig) <- Samp.nom;

################### *** Prepare data for Logistic Regression *** #######################
x = t(DF.sig); 
y = colnames(Metab); y[which(y=="HC")] = 0; y[which(y=="Bipolar")] = 1; y = as.numeric(y)


print('Reading clinical data .....')
CD.DF <- read.csv("../Data/Demographic data_full.csv",header = T, sep = ','); 
CD.DF <- CD.DF[,-ncol(CD.DF)]; rownames(CD.DF) <- CD.DF[,1]; CD.DF = CD.DF[,-1]
CD.DF <- CD.DF[,-c(3:4)];

Map <- read.table("../Data/Sample_mapping.csv", header =F, sep = '\t', stringsAsFactors = FALSE)
idp = match(rownames(CD.DF),Map[,2]);
CD.DF = CD.DF[-which(is.na(idp)),]; idp = idp[-which(is.na(idp))];

print(cbind(rownames(CD.DF),Map[idp,]))
rownames(CD.DF) = Map[idp,1];

idx = match(rownames(CD.DF),rownames(x));
CD.DF = CD.DF[-which(is.na(idx)),]; idx = idx[-which(is.na(idx))];
print(cbind(rownames(CD.DF),rownames(x)[idx]))

x = x[idx,]; CD.DF[,1] = as.factor(CD.DF[,1]);

if(0){ ## Turn this on or off to compute
  
  # Xx = x[,-1]; ## Remove gender ###
  # meta.Xx = x[,1];  ## keep the metaData vector ###
  # 
  Xx = x
  
  Stor.fit.all <- list();
  Lips.combo <- list();
  auc.dist = list();
  TSTl = list();
  TST1l = list();
  ROCX.Lipids = list();
  ROCY.Lipids = list();
  
  AUC.matrix.max = matrix(0,dim(Xx)[2]-1,dim(Xx)[2]+1); ## Max ###
  AUC.matrix.mean = matrix(0,dim(Xx)[2]-1,dim(Xx)[2]+1); ## Mean of the samples or bootstrap ###
  
  ## gg1 = c(ncol(Xx):2); ## indices of columns of remove
  
  gg1 = c(ncol(Xx):1)
  
  #ll = (dim(Xx)[2]-1); ## initialize counters ##
  #vec = c(ll:1);
  
  for(gg in 1:(dim(Xx)[2])){ ## For e.g 20 loops for 21 selected lipids
    
    if(!is.null(dim(Xx))){
      
    Xx1 = Xx; ## Fit all the lipids at once: whole model ##
    print(dim(Xx1));
    
    mm=1;
    #pp1= ncol(Xx1)+1;
    
    temp=list();
    temp1=list();
    
    TST = list();
    TST1 = list();
    
    
    for(pp in 1:((dim(Xx1)[2]+1) )){
      
      print(paste0('case:',gg,'::',pp)); ## counters
      
      auc.dist = list();
      ROCX.mean = list();
      ROCY.mean = list();
      
      
      Samples = 10; ## Change here ##
      
      
      ROCX.Samps = list();
      ROCY.Samps= list();
      Stor.fit = list();
      tst1 = list();
      tst2 = list();
      
      for(tt in 1:Samples){
        
        Train_Rows <- createDataPartition(y, p=0.70, list=FALSE); ## Caret package ##
        
        mydata.train = cbind(Xx1,CD.DF)[Train_Rows, ]; ## make training sets
        mydata.test  = cbind(Xx1,CD.DF)[-Train_Rows, ]; ## make tests sets
        
        y.train = y[Train_Rows];
        y.test  = y[-Train_Rows];
        
        mydata.train = apply(mydata.train,2,as.numeric);
        mydata.test  = apply(mydata.test,2,as.numeric);
        
        training.set.x = cbind(mydata.train);
        test.set.x = cbind(mydata.test);
        
        test.set.y = y.test;
        training.set.y = y.train;
        
        if(length(which(training.set.y==1))>1 & length(which(training.set.y==0))>1 & length(which(test.set.y==1))>1 & length(which(test.set.y==0))>1 ) {
          
          ### ** ridge is alpha = 0 ; ** lasso is alpha=1 ** ### Note: random folds assigned by cv.glmnet ##
          y.fit0 = cv.glmnet(training.set.x,training.set.y,family="binomial", alpha=0, nfolds=5, type.measure = "auc", standardize=TRUE); ## fit and 10 fold cross-validate one for all #
          
          if(is.numeric(y.fit0$lambda.min)){
            
            yhat0 = predict(y.fit0, s=c("lambda.min"), newx=test.set.x, type = "response"); ## predict by response or probalities
            
            pred0 = ROCR::prediction(yhat0,test.set.y)
            nbperf0 = ROCR::performance(pred0, "tpr", "fpr"); 
            AUC0 = ROCR::performance(pred0, "auc")
            AUC0 = unlist(slot(AUC0, "y.values"))
            
            auc.dist[tt] = AUC0;
            Stor.fit[[tt]] =  y.fit0;
            tst1[[tt]] = test.set.x;
            tst2[[tt]] = test.set.y;
            ROCX.Samps[[tt]] = unlist(nbperf0@x.values);
            ROCY.Samps[[tt]] = unlist(nbperf0@y.values);
          }
        }
        
      }
      
      
      
      ROCX.Samps.temps = Filter(Negate(is.null), ROCX.Samps); ## The code removes the NULL list before doing the mean ##
      ROCY.Samps.temps = Filter(Negate(is.null), ROCY.Samps); 
      
      ## Only select profiles with AUC over 0.5 for ROCX and ROCY curves ## 
      
      if(length(ROCX.Samps.temps) > 1 & length(ROCY.Samps.temps) > 1 & length(unique(lapply(ROCX.Samps.temps,length))) ==1 & length(unique(lapply(ROCY.Samps.temps,length))) ==1 &  mean(unlist(auc.dist)) >= 0.5 ){
        ROCX.mean[[pp]] =  apply(simplify2array(ROCX.Samps.temps),1,mean) ### mean of the ROC curves
        ROCY.mean[[pp]] =  apply(simplify2array(ROCY.Samps.temps),1,mean) ### mean of the ROC curves
      }else{
         
        auc.dist[[pp]] = 0; ## Incase you donot have a proper ROC profile put the AUC to zero
        ROCX.mean[[pp]] = lapply(ROCX.mean,mean);
        ROCY.mean[[pp]] = lapply(ROCY.mean,mean);
      }
      
      AVG = mean(unlist(auc.dist));
      #OR
      #AVG = 0.60; # choose ROCs curve close to 0.7;
      AUC.matrix.mean[gg,pp] = AVG; ## Matrix of AUCs and so on but account for each rows;
      

      temp[[pp]] <- Stor.fit[[ which.min( abs(unlist(auc.dist) - AVG) ) ]]; ## NOTE list are sorted; temp[[1]][[22]]; ## Has all the lipids but temp[[1]][[1]]; just elimates the first one;
      temp1[[pp]] <- colnames(Xx1);
      
      TST[[pp]]   = tst1[[ which.min( abs(unlist(auc.dist) - AVG) ) ]]; ## Store the  model closest of mean AUC ##
      TST1[[pp]]  = tst2[[ which.min( abs(unlist(auc.dist) - AVG) ) ]]; ## Store the  model closest of mean AUC ##
      
      if(!is.null(dim(Xx1))){
       Xx1 = Xx1[,-pp]; #print(mm) same as pp
       #Xx1 = Xx[,-mm]; 
       #mm=mm+1; ## re-initialize Xx1
      }
      
    
  }
    
    Stor.fit.all[[gg]] <- temp; ## List of structure Stor.fit.all[[1]][[22]] and so on;
    Lips.combo[[gg]] <- temp1;  ## List of Lipids combo and so on;
    TSTl[[gg]] <- TST;
    TST1l[[gg]] <- TST1;
    
    ROCX.Lipids[[gg]] = ROCX.mean;
    ROCY.Lipids[[gg]] = ROCY.mean;
    
    Xx = Xx[,-gg1[gg]]; ##

    print(paste0('Remove Colmn',ncol(Xx)));
    print(paste0(gg,'::', gg1[gg]));
    print(dim(Xx));
    
  } ##  if(!is.null(dim(Xx1))){
    
  } ### TURN calculation off ###
  
  
  RN = dim(x)[2]:2
  colnames(AUC.matrix.mean) = c('All.PM.Combined.rows',paste0(rev(colnames(x)),'.Eliminate' ));
  rownames(AUC.matrix.mean) <- paste0('PolarMetabolites:',RN);
  
  write.table(AUC.matrix.mean, "../Results/Metabolite_Plots/AUC.matrix.mean.csv",col.names = NA, quote = F,sep=',')
  
  #################### Variable selection ##################################################
  
  x.ind = which(AUC.matrix.mean == AUC.matrix.mean[AUC.matrix.mean!=max(AUC.matrix.mean)], arr.ind = TRUE)[1]
  y.ind = which(AUC.matrix.mean == AUC.matrix.mean[AUC.matrix.mean!=max(AUC.matrix.mean)], arr.ind = TRUE)[2]
  
  model.sel    = Stor.fit.all[[ x.ind ]][[ y.ind ]];
  Lips.sel     = Lips.combo[[ x.ind ]][[ y.ind ]];
  AUC.binomial = AUC.matrix.mean[x.ind, y.ind];
  Test.set.x   = TSTl[[x.ind]][[ y.ind ]];
  Test.set.y   = TST1l[[x.ind]][[ y.ind ]];
  
  ROCXX = ROCX.Lipids[[x.ind]][[ y.ind ]]; ## Mean of the ROCs computed within the program ##
  ROCYY = ROCY.Lipids[[x.ind]][[ y.ind ]]; ## Mean of the ROCs computed within the program ##
  
  ### Store best models and parameters LISTS  ##
  if(1){
    save(model.sel, file = "../Data/model.sel.mets.rda"); 
    save(Lips.sel, file = "../Data/Lips.sel.mets.rda"); 
    save(AUC.binomial, file = "../Data/AUC.binomial.mets.rda"); 
    save(Test.set.x, file = "../Data/Test.set.x.mets.rda"); 
    save(Test.set.y, file = "../Data/Test.set.y.mets.rda"); 
    save(ROCXX, file = "../Data/ROCXX.mets.rda"); 
    save(ROCYY, file = "../Data/ROCYY.mets.rda"); 
  }
  
  
} ## END of if true loop ##

#################### Plot the model ##########################################

### load best models and parameters Lists  ##
if(1){
  load(file = "../Data/model.sel.mets.rda"); 
  load(file = "../Data/Lips.sel.mets.rda"); 
  load(file = "../Data/AUC.binomial.mets.rda"); 
  load(file = "../Data/Test.set.x.mets.rda"); 
  load(file = "../Data/Test.set.y.mets.rda"); 
  load(file = "../Data/ROCXX.mets.rda"); 
  load(file = "../Data/ROCYY.mets.rda");  
}

###################### Generate ROCs ########################################


print('Model fitting and performances')

yhat3 = predict(model.sel, s=c("lambda.min"), newx=Test.set.x, type = "response"); ## predict by response or probalities
yhat3.link = predict(model.sel, s=c("lambda.min"), newx=Test.set.x, type = "link"); ## predict by response or probalities
yhat3.auc = predict(model.sel, s=c("lambda.min"), newx=Test.set.x, type = "coefficients"); ## predict by response or probalities
yhat3.class = predict(model.sel, s=c("lambda.min"), newx=Test.set.x, type = "class"); ## predict by response or probalities


pdf(paste0("../Results/Metabolite_Plots/Logistic_Predictions_Bionomial_meanAUC_meanROCs2",".pdf"), width = 8, height = 8, paper = 'special')
oripar=par()$mar

if(0){
  obj <- roc(Test.set.y,as.vector(yhat3), of="se",ci=F, auc=T, plot=FALSE, smooth=FALSE, percent=F, boot.n=100, boot.stratified=T);
  obj.smooth <- smooth(obj,method="fitdistr", density = "normal", reuse.ci=F)
  
  ciobj <- ci.se(obj.smooth, specificities=1-obj.smooth$specificities, conf.level=0.50, boot.n=100, boot.stratified=T)
  ciobj[,2] <- rowMeans(ciobj[,c(1,3)]);
  
  save(obj, file = "../Data/obj.rda"); ## save roc smoothended object previous fit ###
  save(obj.smooth, file = "../Data/obj.smooth.rda"); ## save roc smoothended object previous fit ###
  save(ciobj, file = "../Data/ciobj.rda"); ## save roc smoothended object previous fit ###
  
}

load (file = "../Data/ciobj.rda")
load (file = "../Data/obj.smooth.rda")
load (file = "../Data/obj.rda")


print("Confidence interval")
print(colMeans(ciobj))

ci.interval <- colMeans(ciobj);


plot(1-obj.smooth$specificities,rev(ciobj[,2]),type="l", xlab = "1 - Specificity [FPR]", ylab = "Sensitivity [TPR]")

par(new=TRUE)

x1 = c(1-obj.smooth$specificities,rev(1-obj.smooth$specificities));
y1 = c(rev(ciobj[,1]), rev(rev(ciobj[,3])) )

polygon(x=x1,
        y=y1, col=adjustcolor("blue", alpha.f = 0.2), border = NA)
par(new=TRUE)
lines(x = c(0,1), y = c(0,1),lty=2, lwd=1.2,col='black',type="l", pch=21)
legend('bottom',legend = paste0('AUC = ',round(obj$auc[1],3),' [95% CI: ',round(ci.interval[1],3),' - ',round(ci.interval[3],3),']'), cex = 0.8);
legend('bottomright',legend = Lips.sel, cex = 0.8);

dev.off()





########################### The END ############################################################