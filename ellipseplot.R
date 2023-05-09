## ellipseplot function
ellipseplot <- function(x, y, factr, 
                        elev=elev, # Ellipse probability level
                        legpos=c("topright","topleft","bottomleft","bottomleft"), # Legend position
                        pcol=pcol, # manual addition of colors, must meet length of factors
                        cexsize=cexsize, # point size
                        ppch=ppch, # Point type, must meet length of factors
                        legcexsize=legcexsize, # legend font size
                        legptsize=legptsize, # legend point size
                        pbgcol=0,
                        axissize=axissize, 
                        linewidth=linewidth, 
                        font=1,
                        Samp.nom=Samp.nom,
                        legs=legs) {
  require(plyr)
  require(car)
  ## Set factor levels
  if(is.factor(factr)) {
    f <- factr
  } else {
    f <- factor(factr, levels=unique(as.character(factr)))
  }
  intfactr <- as.integer(f) # Set integer vector that matches factor levels
  # Checking to make sure length of ppch equals number of factor levels
  if((length(ppch) > 1 & length(unique(intfactr)) != length(ppch))) stop("Can only increase point shape if equal to factor levels")
  
  ## Get data for ellipses
  edf <- data.frame(LV1 = x, LV2=y, factr = f) # create data frame with data and factor
  print(edf)
  print(f)
  print(factr)
  
  ellipses <- dlply(edf, .(factr), function(x) {
    LV1 <- x[,1]
    LV2 <- x[,2]
    dataEllipse(LV1, LV2, levels=elev, robust=TRUE, draw=FALSE) # Get confidence ellipse points from dataEllipse() function by factor level
  })
  
  ## Get range of x and y data
  xrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,1])), max(abs(x))*(-1), max(abs(x)) )))
  yrange <- plotat(range(c(as.vector(sapply(ellipses, function(x) x[,2])), max(abs(y))*(-1), max(abs(y)) )))
  
  ## Set colors for plots
  if(is.null(pcol) != TRUE) { # If colors are supplied by user
    ptcol <- pcol
    pgcol <- paste(pcol, "7e", sep="") # adds opaqueness
  } else { # Default
    pgcol <- c("#e41a1c7e","#377eb87e","#4daf4a7e","#984ea37e","#807f7d7e") # Defaults at 5 colors
    ptcol <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#807f7d") # For opaqueness
  }
  # Plotting graphic
  plot(x,y, type="n", xlab="", ylab="", main="", xlim=range(xrange), ylim=range(yrange), axes=FALSE)
  axis(1, at=xrange, labels=xrange, cex.axis=axissize,lwd=linewidth, font=font)
  axis(2, las=2, cex.axis=axissize,lwd=linewidth, font=font)
  box(lwd=linewidth, font=font)
  abline(h=0, v=0, col="maroon", lty=1, lwd=2) # Adds lines at 0
  legpch <- c() # vector to collect legend pch data
  legcol <- c() # vector to collect legend col data
  ## Not sure why I split this up, might have been an artifact of an older version.
  ## Adds points, ellipse, and determines color specifications for legend 
  
  if(pbgcol==TRUE)  {
    for(i in 1:length(unique(intfactr))){
      #polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
      polygon(ellipses[[i]], col = adjustcolor(pgcol[i], alpha.f = 0.2), border='black')
      legpch[i] <- ppch[i]
      legcol[i] <- ptcol[i]
      points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col='black', bg=ptcol[i],cex=cexsize)
      
    }
    
  } else {
    for(i in 1:length(unique(intfactr))){
      #polygon(ellipses[[i]], col=pgcol[i], border=ptcol[i])
      polygon(ellipses[[i]], col=adjustcolor(pgcol[i], alpha.f = 0.2), border='black')
      legpch[i] <- ppch[i]
      legcol[i] <- ptcol[i]   
      points(x[intfactr==i], y[intfactr==i], pch=ppch[i], col='black', bg=ptcol[i],cex=cexsize)
      
    }
  }
  

  ### Changed here ## PS
  text(x,y,labels = Samp.nom, cex = 1)
  
  ## Legend
  ## legend(x=legpos, legend=levels(f), pch=legpch,
         legend(x=legpos, legend=legs, pch=legpch,
         pt.bg=legcol, col=legcol, bty=1, border=T, pt.cex=legptsize, cex=legcexsize)
}   
