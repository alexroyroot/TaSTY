# scorePaths, the main function in the TaSTY package

scorePaths <- function(pdb, tdb, # input databases
                       pro,sty, #input data
                       c_pro, c_sty, c_tdb,# input constants
                       detectionCutoff, # filter for the proportion of a pathway detected
                       nPerms) {  # number of permutations to perform

  # Filter the pdb based on cut-off for PathwayDetectionProportion
  # Filter out the paths where the receptor wasn't detection
  proteinsDetected <- unique(pro[!is.na(pro$norm), ]$gene)
  pdb <- pdb[pdb$detectedPercent > detectionCutoff &
               pdb$Node1 %in% proteinsDetected &
               pdb$TF %in% proteinsDetected, ]

  # For node1 only, it's sty value is taken as the max or min
  # Add a column to pdb for max phosphorylation event
  sty2 <- separate(sty, col="pSite", sep="_", into=c('gene','residue'), remove=FALSE)
  sty3 <- ddply(sty2, .(gene), summarize, max.norm=max(norm), min.norm=min(norm))
  sty3$norm <- sty3$max.norm

  # if the abs value of the minimum is greater than the maximum, replace it
  sty3[abs(sty3$min.norm) > sty3$max.norm, ]$norm <- sty3[abs(sty3$min.norm) > sty3$max.norm, ]$min.norm
  sty3 <- sty3[ ,c('gene', 'norm')]
  sty3 <- as.data.table(sty3)
  setkey(sty3, gene)


  ## Add columns to pdb to store the results of path scores
  # DT[i, (colvector) := val] basic syntax
  colvector = c('proN1', 'proG2', 'proG3','proG4','proG5','proG6',
                'styN1', 'styN2', 'styN3','styN4','styN5','styN6',
                'TFscore') # create the vector of columns to add
  val = numeric(nrow(pdb)) # initialize with zeros
  pdb[ ,(colvector) := val]

  # Add weights to each datatype
  pro[ ,2] <- c_pro * pro[ ,2]
  sty[ ,2] <- c_sty * sty[ ,2]
  tdb[ ,2] <- c_tdb * tdb[ ,2]

  pdb[ ,'proN1'] <- pro[pdb[ ,c('Node1')], ]$norm
  pdb[ ,'proG2'] <- pro[pdb[ ,c('Gene2')], ]$norm
  pdb[ ,'proG3'] <- pro[pdb[ ,c('Gene3')], ]$norm
  pdb[ ,'proG4'] <- pro[pdb[ ,c('Gene4')], ]$norm
  pdb[ ,'proG5'] <- pro[pdb[ ,c('Gene5')], ]$norm
  pdb[ ,'proG6'] <- pro[pdb[ ,c('Gene6')], ]$norm

  #pdb[ ,'styN1'] <- sty[pdb[ ,c('Node1')], ]$norm  ## THIS IS THE PROBLEMATIC ONE...GET NO MATCHES
  pdb[ ,'styN1'] <- sty3[pdb[ ,c('Node1')], ]$norm
  pdb[ ,'styN2'] <- sty[pdb[ ,c('Node2')], ]$norm
  pdb[ ,'styN3'] <- sty[pdb[ ,c('Node3')], ]$norm
  pdb[ ,'styN4'] <- sty[pdb[ ,c('Node4')], ]$norm
  pdb[ ,'styN5'] <- sty[pdb[ ,c('Node5')], ]$norm
  pdb[ ,'styN6'] <- sty[pdb[ ,c('Node6')], ]$norm

  pdb[ ,'TFscore'] <- tdb[pdb[ ,c('TF')], ]$TFscore
  pdb[ ,'ReceptorProteinExpression'] <- pro[pdb[ ,c('Node1')], ]$norm

  # Compute total scores by path
  #tmp <- rowSums(pdb[ ,c('cnaN1', 'cnaG2', 'cnaG3','cnaG4','cnaG5','cnaG6',
  #              'rnaN1', 'rnaG2', 'rnaG3','rnaG4','rnaG5','rnaG6',
  #             'mutN1', 'mutG2', 'mutG3','mutG4','mutG5','mutG6',
  #            'proN1', 'proG2', 'proG3','proG4','proG5','proG6',
  #           'styN1', 'styN2', 'styN3','styN4','styN5','styN6',
  #          'TFscore')], na.rm=TRUE)

  #pdb[,('sumScore'):=tmp]

  pdb$sumScore <- rowSums(pdb[ ,c('proN1', 'proG2', 'proG3','proG4','proG5','proG6',
                                  'styN1', 'styN2', 'styN3','styN4','styN5','styN6',
                                  'TFscore')], na.rm=TRUE) / pdb$pathLength

  # Perform a large number of permutation tests
  # randomly swap rows in the data and calculate path scores
  # Create a matrix of null scores that is initialized with 0
  # if the null matrix isn't initialized with 0 then NAs will dominate
  # The presence of NAs will throw off the p-value calculation
  nullScores <- matrix(0, nrow = nrow(pdb), ncol=nPerms)
  pdb_tmp <- pdb
  for (i in 1 : nPerms){
    print(c("On permutation ",i, " of ", nPerms))
    pro[ ,2] <- pro[sample(1:nrow(pro),nrow(pro), replace=FALSE), 2]
    sty[ ,2] <- sty[sample(1:nrow(sty),nrow(sty), replace=FALSE), 2]
    tdb[ ,2] <- tdb[sample(1:nrow(tdb), nrow(tdb), replace=FALSE), 2]

    pdb_tmp[ ,'proN1'] <- pro[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'proG2'] <- pro[pdb[ ,c('Gene2')], ]$norm
    pdb_tmp[ ,'proG3'] <- pro[pdb[ ,c('Gene3')], ]$norm
    pdb_tmp[ ,'proG4'] <- pro[pdb[ ,c('Gene4')], ]$norm
    pdb_tmp[ ,'proG5'] <- pro[pdb[ ,c('Gene5')], ]$norm
    pdb_tmp[ ,'proG6'] <- pro[pdb[ ,c('Gene6')], ]$norm

    #pdb_tmp[ ,'styN1'] <- sty[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'styN1'] <- sty3[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'styN2'] <- sty[pdb[ ,c('Node2')], ]$norm
    pdb_tmp[ ,'styN3'] <- sty[pdb[ ,c('Node3')], ]$norm
    pdb_tmp[ ,'styN4'] <- sty[pdb[ ,c('Node4')], ]$norm
    pdb_tmp[ ,'styN5'] <- sty[pdb[ ,c('Node5')], ]$norm
    pdb_tmp[ ,'styN6'] <- sty[pdb[ ,c('Node6')], ]$norm

    pdb_tmp[ ,'TFscore'] <- tdb[pdb[ ,c('TF')], ]$TFscore

    # Compute total scores by path
    nullScores[ ,i] <- rowSums(pdb_tmp[ ,c('proN1', 'proG2', 'proG3','proG4','proG5','proG6',
                                           'styN1', 'styN2', 'styN3','styN4','styN5','styN6',
                                           'TFscore')], na.rm=TRUE) / pdb_tmp$pathLength
  }
  # APPLY VARIOUS CUTOFFS TO GET THE SCORE SIGNIFICANCE LEVELS
  print("Now calculating statistical significance")
  cut95 <- apply(nullScores, 1, quantile, 0.95)
  cut99 <- apply(nullScores, 1, quantile, 0.99)
  cut999 <- apply(nullScores, 1, quantile, 0.999)
  cut9999 <- apply(nullScores, 1, quantile, 0.9999)

  cut05 <- apply(nullScores, 1, quantile, 0.05)
  cut01 <- apply(nullScores, 1, quantile, 0.01)
  cut001 <- apply(nullScores, 1, quantile, 0.001)
  cut0001 <- apply(nullScores, 1, quantile, 0.0001)
  #testfun <- apply(nullScores, 1, function(x) ecdf(x))


  pdb[ ,('pVal') := numeric(nrow(pdb))]
  pdb[ ,pVal := pVal + 1]
  pdb[pdb$sumScore > cut95 | pdb$sumScore < cut05, ]$pVal <- 0.05
  pdb[pdb$sumScore > cut99 | pdb$sumScore < cut01, ]$pVal <- 0.01
  pdb[pdb$sumScore > cut999 | pdb$sumScore < cut001, ]$pVal <- 0.001
  pdb[pdb$sumScore > cut9999 | pdb$sumScore < cut0001, ]$pVal <- 0.0001
  #for (i in 1 : length(testfun)) {
  # tmpecdf <- testfun[[i]]
  #  pdb[i, ]$pVal <- tmpecdf(pdb[i, ]$sumScore)
  #}
  pdb[ ,('padj') := p.adjust(pdb$pVal, method='BH')]

  # Useful to return the nullScores for debugging
  #fncout <- list(nullScores, pdb)
  #return(fncout)
  #return(pdb[pdb$pVal <= 0.05, ])
  return(pdb)
}






