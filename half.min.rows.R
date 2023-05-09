half.minimum.rows <- function(DF)
{
  
  DF = t(DF);
  DF[DF== 0] <- NA
  mins = apply(DF, 1, FUN = function(x) {min(x[x > 0], na.rm = T)}); mins = mins/2;
  
  for (tt in 1:dim(DF)[1]){
    DF[tt,][is.na(DF[tt,])] <- mins[tt];
  }
  #mins = apply(DFF, 1, FUN=min); mins = mins/2;
  DF = t(DF);
  
  return(DF);

}