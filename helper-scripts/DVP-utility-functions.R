#### DVP Utility Functions ##########

## DVP function 
run_DVP <- function(mvals, design){
  
  # fit the model
  fitvar <- varFit(mvals, design = design, coef = c(1,2))
  
  # extract standard errors
  std_err <- fitvar$stdev.unscaled[, 2] * fitvar$sigma
  std_err_df <- data.frame(std_err)
  std_err_df$cpg <- rownames(std_err_df)
  
  # identify top var CpG sites
  topDV <- topVar(fitvar, coef = 2, number = nrow(mvals))
  topDV <- data.frame(Probe_ID = rownames(topDV), topDV)
  
  # adjust p-values with Bonferroni corrections
  topDV$Bonferroni <- p.adjust(topDV$P.Value, method = "bonferroni")
  
  # join results table with SE table
  topDV$cpg <- rownames(topDV)
  topDV <- merge(topDV, std_err_df, by = 'cpg')
  rownames(topDV) <- topDV$cpg
  
  return(topDV)
}
