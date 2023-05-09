## Axis legends for PCA output using prcomp() function
PCAvarAxis <- function(PCA, decimal=1) {
  pcavar <- round((PCA$sdev^2)/sum((PCA$sdev^2)),3)*100   #Calculate % variance explained
  PC1var <- paste("Principal Component 1 (", pcavar[1], "%)", sep="")
  PC2var <- paste("Principal Component 2 (", pcavar[2], "%)", sep="")
  PC3var <- paste("Principal Component 3 (", pcavar[3], "%)", sep="")
  PC4var <- paste("Principal Component 4 (", pcavar[4], "%)", sep="") 
  PC5var <- paste("Principal Component 5 (", pcavar[5], "%)", sep="")     
  return(list(PC1=PC1var, PC2=PC2var, PC3=PC3var, PC4=PC4var, PC5=PC5var))
}   