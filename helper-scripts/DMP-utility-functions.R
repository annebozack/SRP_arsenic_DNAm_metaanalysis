#### DMP Utility Functions ##########

run_DMP <- function(mvals, design){
  
  # fit model
  l_fit <- limma::lmFit(object = mvals, design = design)
  
  # extract standard errors
  std_err <- l_fit$stdev.unscaled[,2]*l_fit$sigma
  std_err_df <- data.frame(std_err)
  std_err_df$cpg <- rownames(std_err_df)
  
  e_fit <- limma::eBayes(l_fit, robust = TRUE)
  
  # extract results and add Bonferroni correction
  p_top <- limma::topTable(e_fit, adjust = "BH", coef = 2,
                           num = Inf, confint = TRUE)
  p_top <- p_top[order(p_top$P.Value), , drop = FALSE]
  p_top$adj.P.Val.bonf <- topTable(e_fit, adjust="bonferroni", coef=2,
                                   number = Inf)$adj.P.Val
  
  # merge results and standard errors
  p_top$cpg <- rownames(p_top)
  p_top <- merge(p_top, std_err_df, by = 'cpg')
  rownames(p_top) <- p_top$cpg
  
  return(p_top)
}