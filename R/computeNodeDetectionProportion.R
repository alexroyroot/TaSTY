# A function to calculate the proportion of nodes in a path that are detected

computeNodeDetectionProportion <- function(pdb, sty) {

  styDetected <- unique(sty[!is.na(sty$norm), ]$pSite)

  ## Add indicator variables if the nodes are detected
  pdb$i2 <- 0
  pdb$i3 <- 0
  pdb$i4 <- 0
  pdb$i5 <- 0
  # pdb$i6 <- 0 not needed - only score intermediate nodes

  pdb[pdb$Node2 %in% styDetected, ]$i2 <- 1
  pdb[pdb$Node3 %in% styDetected & pdb$pathLength %in% c(4,5,6), ]$i3 <- 1
  pdb[pdb$Node4 %in% styDetected & pdb$pathLength %in% c(5,6), ]$i4 <- 1
  pdb[pdb$Node5 %in% styDetected & pdb$pathLength %in% c(6), ]$i5 <- 1
  #pdb[pdb$Node6 %in% styDetected, ]$i6 <- 1

  ## Determine how many of the path nodes are detected
  pdb$iSum <- 0
  pdb$iSum <- rowSums(pdb[ ,c('i2', 'i3', 'i4', 'i5')], na.rm=TRUE)
  pdb$detectedPercent <- pdb$iSum / (pdb$pathLength-2)
  return(pdb)
}
