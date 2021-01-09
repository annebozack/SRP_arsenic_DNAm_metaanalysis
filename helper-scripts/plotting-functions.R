# Plotting functions for EWAS and DVP
manhattan <- function(probe, locations, colors = c("#2D728F", "#2B4162")) {
  
  probe <- cbind(
    probe, locations[, 'chr'][match(rownames(probe), rownames(locations))]
  )
  probe <- cbind(
    probe, locations[, 'pos'][match(rownames(probe), rownames(locations))]
  )
  colnames(probe)[c(ncol(probe) - 1, ncol(probe))] = c('chr', 'pos')
  probe$chr = as.numeric(gsub("chr", "", probe$chr))
  
  don <- probe %>% 
    # Compute chromosome size
    group_by(chr) %>% summarise(chr_len=as.numeric(max(pos))) %>% 
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(chr_len)-chr_len) %>% dplyr::select(-chr_len) %>%
    # Add this info to the initial dataset
    left_join(probe, ., by=c("chr"="chr")) %>%
    # Add a cumulative position of each site
    arrange(chr, pos) %>% mutate(poscum=pos+tot) # %>%
  # Add highlight and annotation information
  # mutate(is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  # mutate(is_annotate_fdr=ifelse(adj.P.Val<0.05 & adj.P.Val.bonf>0.05, "yes", "no")) %>%
  # mutate(is_annotate_bonf=ifelse(adj.P.Val.bonf <0.05, "yes", "no"))
  # don$alpha[don$is_annotate_fdr == 'no'] = 0.6
  # don$alpha[don$is_annotate_fdr == 'yes'] = 0
  # Prepare X axis
  axisdf <- don %>%
    group_by(chr) %>%
    summarize(center = (max(poscum) + min(poscum))/2)
  don <- merge(don, probe[, c(1,11)], on = 'cpg', all.x = T)
  
  manhattan <- ggplot(don, aes(x = poscum, y = -log10(P.Value))) +
    geom_point(aes(color = as.factor(chr)), size=0.8, alpha = 0.55) +
    scale_color_manual(values = rep(colors, dim(table(don$chr)))) +
    # p-value cutoffs
    geom_hline(yintercept = -log10(0.05/nrow(don)),
               colour = '#AB3428', size=.5) +
    #   geom_hline(yintercept=-log10(max(don$P.Value[don$adj.P.Val < 0.05])), colour='#AB3428', size=.4) +
    # custom X axis:
    scale_x_continuous(expand = c(0, 0),
                       limits = c(min(don$poscum),max(don$poscum)),
                       label = axisdf$chr, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0)) + 
    # Custom the theme:
    theme_minimal() + theme( 
      legend.position="none", panel.border = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.2, color = 'gray65')) +
    theme(axis.title.y = element_blank()) +
    theme(text = element_text(size = 7.5)) + 
    labs(y=expression(-log[10]*"(P-value)"), x='Chromosome') 
  return(manhattan)
}

gg_qqplot = function(pvector) {
  
  l = round(lambda(pvector), 3)
  o = -log10(sort(pvector, decreasing = FALSE))
  e = -log10(ppoints(length(pvector)))
  df = data.frame(o = o, e = e)
  ggplot(df, aes(e, o)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(intercept = 0, slope = 1, color = '#AB3428') +
    labs(y = expression(Observed ~ ~-log[10](italic(p))),
         x = expression(Expected ~ ~-log[10](italic(p)))) +
    theme_classic() +
    annotate("text", x = 1, y = 5,
             label = paste0('lambda = ', l))
  
}


volcano = function(probe) {
  volcano = ggplot(probe, aes(x = logFC, y  = -log10(P.Value))) + 
    geom_point(size = 0.8, alpha = 0.4) + 
      geom_hline(aes(yintercept = -log10(0.05/nrow(probe)),
                     linetype = 'Bonferroni threshold'),
                 color = "#AB3428", size = 0.5) + 
      # geom_hline(aes(yintercept = -log10(max(P.Value[adj.P.Val < 0.05])),
      #                linetype = 'FDR threshold'), color = "#AB3428", size = 0.3) + 
      scale_linetype_manual(
        name = '', values = c(1,2),
        guide = guide_legend(
          override.aes = list(
          color = c("#AB3428", "#AB3428")
          )
        )
      ) +
      theme_minimal() + 
      labs(y=expression(-log[10]*"(P-value)"), x='Effect estimate') +
      theme(panel.grid.minor.y = element_blank()) +
      theme(text = element_text(size = 8)) + 
      theme(legend.position="none") +
    scale_y_continuous(expand = c(0, 0)) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(size = 0.2, color = 'gray65'),
          panel.grid.major.x = element_line(size = 0.2, color = 'gray65'))
  return(volcano)
}