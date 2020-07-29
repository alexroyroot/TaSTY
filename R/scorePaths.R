# function to determine path scores

scorePaths <- function(pdb, tdb, # input databases
                       cna,
                       rna,
                       mut,
                       pro,
                       sty, #input data
                       c_cna, c_rna, c_mut, c_pro, c_sty, c_tdb,# input constants
                       nPerms) {  # number of permutations to perform

  ## Add columns to pdb to store the results of path scores
  # DT[i, (colvector) := val] basic syntax
  colvector = c('cnaN1', 'cnaG2', 'cnaG3','cnaG4','cnaG5','cnaG6',
                'rnaN1', 'rnaG2', 'rnaG3','rnaG4','rnaG5','rnaG6',
                'mutN1', 'mutG2', 'mutG3','mutG4','mutG5','mutG6',
                'proN1', 'proG2', 'proG3','proG4','proG5','proG6',
                'styN1', 'styN2', 'styN3','styN4','styN5','styN6',
                'TFscore') # create the vector of columns to add
  val = numeric(nrow(pdb)) # initialize with zeros
  pdb[ ,(colvector) := val]

  # Add weights to each datatype
  cna[ ,2] <- c_cna * cna[ ,2]
  rna[ ,2] <- c_rna * rna[ ,2]
  mut[ ,2] <- c_mut * mut[ ,2]
  pro[ ,2] <- c_pro * pro[ ,2]
  sty[ ,2] <- c_sty * sty[ ,2]
  tdb[ ,2] <- c_tdb * tdb[ ,2]

  pdb[ ,'cnaN1'] <- cna[pdb[ ,c('Node1')], ]$norm
  pdb[ ,'cnaG2'] <- cna[pdb[ ,c('Gene2')], ]$norm
  pdb[ ,'cnaG3'] <- cna[pdb[ ,c('Gene3')], ]$norm
  pdb[ ,'cnaG4'] <- cna[pdb[ ,c('Gene4')], ]$norm
  pdb[ ,'cnaG5'] <- cna[pdb[ ,c('Gene5')], ]$norm
  pdb[ ,'cnaG6'] <- cna[pdb[ ,c('Gene6')], ]$norm

  pdb[ ,'rnaN1'] <- rna[pdb[ ,c('Node1')], ]$norm
  pdb[ ,'rnaG2'] <- rna[pdb[ ,c('Gene2')], ]$norm
  pdb[ ,'rnaG3'] <- rna[pdb[ ,c('Gene3')], ]$norm
  pdb[ ,'rnaG4'] <- rna[pdb[ ,c('Gene4')], ]$norm
  pdb[ ,'rnaG5'] <- rna[pdb[ ,c('Gene5')], ]$norm
  pdb[ ,'rnaG6'] <- rna[pdb[ ,c('Gene6')], ]$norm

  pdb[ ,'mutN1'] <- mut[pdb[ ,c('Node1')], ]$norm
  pdb[ ,'mutG2'] <- mut[pdb[ ,c('Gene2')], ]$norm
  pdb[ ,'mutG3'] <- mut[pdb[ ,c('Gene3')], ]$norm
  pdb[ ,'mutG4'] <- mut[pdb[ ,c('Gene4')], ]$norm
  pdb[ ,'mutG5'] <- mut[pdb[ ,c('Gene5')], ]$norm
  pdb[ ,'mutG6'] <- mut[pdb[ ,c('Gene6')], ]$norm

  pdb[ ,'proN1'] <- pro[pdb[ ,c('Node1')], ]$norm
  pdb[ ,'proG2'] <- pro[pdb[ ,c('Gene2')], ]$norm
  pdb[ ,'proG3'] <- pro[pdb[ ,c('Gene3')], ]$norm
  pdb[ ,'proG4'] <- pro[pdb[ ,c('Gene4')], ]$norm
  pdb[ ,'proG5'] <- pro[pdb[ ,c('Gene5')], ]$norm
  pdb[ ,'proG6'] <- pro[pdb[ ,c('Gene6')], ]$norm

  pdb[ ,'styN1'] <- sty[pdb[ ,c('Node1')], ]$norm
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

  pdb$sumScore <- rowSums(pdb[ ,c('cnaN1', 'cnaG2', 'cnaG3','cnaG4','cnaG5','cnaG6',
                                  'rnaN1', 'rnaG2', 'rnaG3','rnaG4','rnaG5','rnaG6',
                                  'mutN1', 'mutG2', 'mutG3','mutG4','mutG5','mutG6',
                                  'proN1', 'proG2', 'proG3','proG4','proG5','proG6',
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
    print(i)

    cna[ ,2] <- cna[sample(1:nrow(cna),nrow(cna), replace=FALSE), 2]
    rna[ ,2] <- rna[sample(1:nrow(rna),nrow(rna), replace=FALSE), 2]
    mut[ ,2] <- mut[sample(1:nrow(mut),nrow(mut), replace=FALSE), 2]
    pro[ ,2] <- pro[sample(1:nrow(pro),nrow(pro), replace=FALSE), 2]
    sty[ ,2] <- sty[sample(1:nrow(sty),nrow(sty), replace=FALSE), 2]
    tdb[ ,2] <- tdb[sample(1:nrow(tdb), nrow(tdb), replace=FALSE), 2]

    pdb_tmp[ ,'cnaN1'] <- cna[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'cnaG2'] <- cna[pdb[ ,c('Gene2')], ]$norm
    pdb_tmp[ ,'cnaG3'] <- cna[pdb[ ,c('Gene3')], ]$norm
    pdb_tmp[ ,'cnaG4'] <- cna[pdb[ ,c('Gene4')], ]$norm
    pdb_tmp[ ,'cnaG5'] <- cna[pdb[ ,c('Gene5')], ]$norm
    pdb_tmp[ ,'cnaG6'] <- cna[pdb[ ,c('Gene6')], ]$norm

    pdb_tmp[ ,'rnaN1'] <- rna[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'rnaG2'] <- rna[pdb[ ,c('Gene2')], ]$norm
    pdb_tmp[ ,'rnaG3'] <- rna[pdb[ ,c('Gene3')], ]$norm
    pdb_tmp[ ,'rnaG4'] <- rna[pdb[ ,c('Gene4')], ]$norm
    pdb_tmp[ ,'rnaG5'] <- rna[pdb[ ,c('Gene5')], ]$norm
    pdb_tmp[ ,'rnaG6'] <- rna[pdb[ ,c('Gene6')], ]$norm

    pdb_tmp[ ,'mutN1'] <- mut[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'mutG2'] <- mut[pdb[ ,c('Gene2')], ]$norm
    pdb_tmp[ ,'mutG3'] <- mut[pdb[ ,c('Gene3')], ]$norm
    pdb_tmp[ ,'mutG4'] <- mut[pdb[ ,c('Gene4')], ]$norm
    pdb_tmp[ ,'mutG5'] <- mut[pdb[ ,c('Gene5')], ]$norm
    pdb_tmp[ ,'mutG6'] <- mut[pdb[ ,c('Gene6')], ]$norm

    pdb_tmp[ ,'proN1'] <- pro[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'proG2'] <- pro[pdb[ ,c('Gene2')], ]$norm
    pdb_tmp[ ,'proG3'] <- pro[pdb[ ,c('Gene3')], ]$norm
    pdb_tmp[ ,'proG4'] <- pro[pdb[ ,c('Gene4')], ]$norm
    pdb_tmp[ ,'proG5'] <- pro[pdb[ ,c('Gene5')], ]$norm
    pdb_tmp[ ,'proG6'] <- pro[pdb[ ,c('Gene6')], ]$norm

    pdb_tmp[ ,'styN1'] <- sty[pdb[ ,c('Node1')], ]$norm
    pdb_tmp[ ,'styN2'] <- sty[pdb[ ,c('Node2')], ]$norm
    pdb_tmp[ ,'styN3'] <- sty[pdb[ ,c('Node3')], ]$norm
    pdb_tmp[ ,'styN4'] <- sty[pdb[ ,c('Node4')], ]$norm
    pdb_tmp[ ,'styN5'] <- sty[pdb[ ,c('Node5')], ]$norm
    pdb_tmp[ ,'styN6'] <- sty[pdb[ ,c('Node6')], ]$norm

    pdb_tmp[ ,'TFscore'] <- tdb[pdb[ ,c('TF')], ]$TFscore

    # Compute total scores by path
    nullScores[ ,i] <- rowSums(pdb_tmp[ ,c('cnaN1', 'cnaG2', 'cnaG3','cnaG4','cnaG5','cnaG6',
                                           'rnaN1', 'rnaG2', 'rnaG3','rnaG4','rnaG5','rnaG6',
                                           'mutN1', 'mutG2', 'mutG3','mutG4','mutG5','mutG6',
                                           'proN1', 'proG2', 'proG3','proG4','proG5','proG6',
                                           'styN1', 'styN2', 'styN3','styN4','styN5','styN6',
                                           'TFscore')], na.rm=TRUE)/pdb_tmp$pathLength
  }
  # APPLY VARIOUS CUTOFFS TO GET THE SCORE SIGNIFICANCE LEVELS
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
