plotat <- function(RANGE) {
  if(length(RANGE) != 2) stop("RANGE argument must have a length of 2")
  if(RANGE[1] > RANGE[2]) stop("First element in RANGE must be smaller then second element")
  prettyres <- pretty(sprintf("%.2f",RANGE[1]):sprintf("%.2f",RANGE[2]), 7)
  while((min(prettyres) < RANGE[1]) == FALSE) {
    prdiff <- prettyres[2] - prettyres[1]
    prettyres[length(prettyres) + 1] <- prettyres[1] - prdiff
    prettyres <- sort(prettyres)
  } 
  while((max(prettyres) > RANGE[2]) == FALSE) {
    prdiff <- prettyres[2] - prettyres[1]
    prettyres[length(prettyres) + 1] <- prettyres[length(prettyres)] + prdiff
    prettyres <- sort(prettyres)    
  }   
  plotticks <- as.numeric(sprintf("%.2f",prettyres))
  plotticks
}