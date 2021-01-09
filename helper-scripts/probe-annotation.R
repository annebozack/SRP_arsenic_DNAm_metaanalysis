# create arsenic annotation data file

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
data(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# 
# data(Locations)
# data(Other)
# annot <- data.frame(Other)
# loc <- data.frame(Locations)
# probes <- merge(loc, annot, by = "row.names")
# rownames(probes) <- probes$Row.names
# probes <- probes[, c(2:4, 9:18)]
# 
# data(SNPs.Illumina)
# data(Islands.UCSC)
# snps <- data.frame(SNPs.Illumina)
# islands <- data.frame(Islands.UCSC)
# probes2 <- merge(snps, islands, by = "row.names")
# rownames(probes2) <- probes2$Row.names
# probes2 <- probes2[, c(2:5)]
# 
# annotation <- merge(probes, probes2, by = "row.names")
# colnames(annotation)[1] <- "Probe_ID"
# 
# save(
#   annotation,
#   file = "../arsenic-epigenetics-meta/Chile/arsenic_annot.RData"
# )

load("../../arsenic-epigenetics-meta/Chile/arsenic_annot.Rdata")

annotate_probes <- function(probeID_column = 1, sig_dat) {
  
  sig_probes <- as.character(sig_dat[, 1])
  sig_dat$Probe_ID <- sig_dat[, 1]
  annot <- annotation[(annotation$Probe_ID %in% sig_probes), ]
  dat <- merge(sig_dat, annot, by = "Probe_ID")
  
  return(dat)
}

library(TxDb.Hsapiens.UCSC.hg19.knownGene)

annotate_regions <- function(regions_data) {
  
  GR <- makeGRangesFromDataFrame(regions_data)
  GR1 <- subsetByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR)
  sym <- as.data.frame(org.Hs.egSYMBOL)
  merged <- merge(GR1, sym, by = "gene_id")
  GR2 <- mergeByOverlaps(genes(TxDb.Hsapiens.UCSC.hg19.knownGene), GR)
  m <- DataFrame(geneid = GR2[, 2], GR_region = GR2[, 3])
  m$gene_id <- make.names(m$geneid, unique = TRUE)
  m$gene_id <- substring(m$gene_id, 2)
  #m1 <- names(m$GR_region)
  #names(m$GR_region) <- make.names(m1, unique=TRUE)
  merged2 <- merge(m, merged,by = "gene_id")
  annot <- data.frame(merged2)
  #annot <- data.frame(merged)
  annot1 <- annot[,c(2:5,13,6,9:11)]
  colnames(annot1) <- c("gene_id", "chr", "start", "end", "gene_symbol",
                        "width", "gene_start", "gene_end", "gene_width")
  join <- left_join(regions_data, annot1, by = c("chr", "start", "end"))
  join$width <- join$end - join$start
  join1 <- filter(join, width != 0)
  
  return(join1)
}

annotate_probes_DMR <- function(metaresults, annotation) {
  
  # extract CpG sites and nominal p-values
  metaresults <- metaresults %>%
    dplyr::select(Name, P.value)
  
  # format annotation data
  annotation <- annotation %>%
    as.data.frame %>%
    dplyr::select(Probe_ID, chr, pos, strand)
  
  # merge with annotation data
  metaresults <- metaresults %>%
    left_join(annotation, by = c("Name" = "Probe_ID"))
  
  # cleanup to bedfile column ordering
  metaresults %>%
    mutate(
      chrom = chr,
      start = pos,
      end = pos,
      name = Name,
      strand = strand
    ) %>%
    arrange(chrom, start) %>%
    dplyr::select(chrom, start, end, name, strand, P.value)
}
