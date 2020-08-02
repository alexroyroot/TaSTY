# Compute a score for transcription factor targets


computeTFscore <- function(tdb, rna) {
  ### Prepare the TF scores

  tdb2 <- merge(tdb, rna, by.y='gene', by.x = "Target", all.x=TRUE, all.y=FALSE)

  # Switch the sign when a TF represses a target
  tdb2[tdb2$ControlType %in% 'Repression', ]$norm <- -1 * tdb2[tdb2$ControlType %in% 'Repression', ]$norm


  ## TRANSCRIPTION FACTOR SCORE - MEAN OR SUM?
  t1 <- ddply(tdb2, .(TF), summarize, TFscore = mean(norm, na.rm=TRUE))

  # Take the TFs with non-zero target expression
  t2 <- as.data.table(t1[!is.na(t1$TFscore) & !t1$TFscore == 0, ])

  setkey(t2, TF)
  return(t2)

}
