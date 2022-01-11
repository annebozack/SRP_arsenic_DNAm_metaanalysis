---
title: "FACT 450k processing"
output:
  html_document:
    toc: true
    toc_float: true
---

## Required Packages
```{r library,eval=FALSE}
library(here)
library(tidyverse)
library(minfi)
library(ggplot2)
library(ENmix)
library(stringr)
library(geneplotter)
library(ChAMP)
library(wateRmelon)
library(FlowSorted.Blood.450k)
```

## Load Data
```{r data, warning=FALSE,message=FALSE,eval=FALSE}

# setting project and data directories

change_here = function(new_path){
  new_root = here:::.root_env
  new_root$f = function(...){file.path(new_path, ...)}
  assignInNamespace(".root_env", new_root, ns = "here")
}

change_here('/S2_121611_run 2/')

sample_info = read.csv(here("Gamble_meth450_SampleSheet_121611_2.csv"))

# paths to idats

files_idat_raw = paste0(here(), 'idats/', sample_info$Sentrix_ID, '/', sample_info$Sentrix_ID, '_', sample_info$Sentrix_Position)

# combining with pheno file

pheno = read.csv(here('phenodata.csv'))

pheno = cbind(pheno, files_idat_raw)
colnames(pheno)[ncol(pheno)] = 'Basename'
```

Reading idats

```{r,warning=FALSE,eval=FALSE}
RGset = read.metharray.exp(targets=pheno, extended = TRUE)
```

Cleaning pheno data
```{r,warning=FALSE,eval=FALSE}
pheno = pData(RGset)
pheno = data.frame(pheno)
pheno$Sample_Name = rownames(pheno)

pheno$expGroup[pheno$Group == 'High'] = 1
pheno$expGroup[pheno$Group == 'Low'] = 0
pheno$expGroup = as.factor(pheno$expGroup)

pheno$LandOwn = factor(pheno$LandOwn, levels = c('0', '1'))
pheno$FolDef = factor(pheno$FolDef, levels = c('0', '1'))
pheno$Betel = factor(pheno$Betel, levels = c('0', '1'))

# Adding smoking data
fact = read.csv('factDataset_with_bAsMetabs_corrected2.csv')
colnames(fact)[3] = 'factid'
pheno = merge(pheno, fact[,c(3,192)], by = 'factid')
pheno$Cig.Ever = factor(pheno$Cig.Ever, levels = c(0,1))
ids = colnames(RGset)
pheno = pheno[match(ids, pheno$Sample_Name),]
table(pheno$Sample_Name == colnames(RGset))
# TRUE 
# 48 
```

## Internal Quality Control
Produces a PDF QC report for Illumina Infinium Human Methylation 450k arrays, which is useful for identifying failed samples.

```{r,warning=FALSE,cache=TRUE,eval = FALSE}
qcReport(RGset, sampGroups = pheno$expGroup, pdf="fact_ControlProbes.pdf")
```

The ControlProbesPlot.pdf file is a 19-page report of the internal quality control information.

## Outlier information

A number of functions are run sequentially on the RGset object.First outlier values are thresholded using fixMethOutliers. Then qc is performed using getQC and then sample specific sex is estimated using getSex.

```{r,warning=FALSE,cache=TRUE,eval = FALSE}
# Creation Of A MethylSet Without Normalization (minifi)

Mset = preprocessRaw(RGset)

QC_minfi = minfiQC(Mset, fixOutliers = TRUE)

pdf(file = here("QC", "PotentialOutliers.pdf"))

# Produce QC diagnostic plots

plotQC(QC_minfi$qc)
dev.off()
```

### Potential Outliers

```{r, out.width = '100%'}
knitr::include_graphics("QC/PotentialOutliers.png")
```

## Pre Filtering

```{r, warning=FALSE,eval = FALSE}
#Plot quality control plot, mdsplot, densityPlot, dendrogram.

champ.QC(beta = getBeta(Mset), pheno=pheno$expGroup, mdsPlot=TRUE, densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE,Feature.sel="None", resultsDir=here::here("QC", "QC_preFiltering"))
```

### Pre Filtering Density 
```{r, out.width = '100%'}
knitr::include_graphics("QC/QC_preFiltering/raw_densityPlot.png")
```

### Pre Filtering Beta Values
```{r, out.width = '100%'}
knitr::include_graphics("QC/QC_preFiltering/raw_mdsPlot.png")
```

### Pre Filtering Dendogram
```{r, out.width = '100%'}
knitr::include_graphics("QC/QC_preFiltering/raw_SampleCluster.png")
```

### Pre-Filtering SVD
Imputing missing beta values before performing SVD
```{r, warning=FALSE,eval = FALSE}
table(is.na(getBeta(Mset)))
# FALSE     TRUE 
# 23303867  709 

beta = getBeta(Mset)
betaImpute = champ.impute(beta=beta, pd=pheno, k=5, ProbeCutoff=0.2, SampleCutoff=0.1)

champ.SVD(beta = betaImpute$beta, rgSet=RGset, pd=pheno[-c(1,2,3,9,18,19,20)], RGEffect=TRUE, Rplot=FALSE, resultsDir=here::here("QC_preFiltering"))
```

```{r, out.width = '100%'}
knitr::include_graphics("QC_preFilteringSVDsummary.png")
```

## Post Filtering
Extract informations for data quanlity controls: detection P values, number of beads and averaged bisulfite conversion intensity. The function can also identify low quality samples and probes, as well as outlier samples based on total intensity or beta value distribution.


```{r, warning=FALSE,eval = FALSE}
qc = ENmix::QCinfo(RGset, distplot=FALSE)

# 0  samples with percentage of low quality CpG value greater than  0.05  or bisulfite intensity less than  12774.17 
# 9286  CpGs with percentage of low quality value greater than  0.05 
# Ploting qc_sample.jpg ...Done
# Ploting qc_CpG.jpg ...Done
# Identifying ourlier samples based on beta or total intensity values...
# After excluding low quality samples and CpGs
# 0  samples are outliers based on averaged total intensity value 
# 0  samples are outliers in beta value distribution 
# 0  outlier samples were added into badsample list


filtered = ChAMP::champ.filter(beta=getBeta(Mset), pd = pheno, beadcount = qc$nbead, detP = qc$detP, arraytype = "450K", SampleCutoff=0.05)
              
  # You have inputed beta for Analysis.

  # pd file provided, checking if it's in accord with Data Matrix...
    # pd file check success.

  # Parameter filterDetP is TRUE, checking if detP in accord with Data Matrix...
    # detP check success.

  # Parameter filterBeads is TRUE, checking if beadcount in accord with Data Matrix...
    # beadcount check success.

  # parameter autoimpute is TRUE. Checking if the conditions are fulfilled...
    # !!! ProbeCutoff is 0, which means you have no needs to do imputation. autoimpute has been reset FALSE.

  # Checking Finished :filterDetP,filterBeads,filterMultiHit,filterSNPs,filterNoCG,filterXY would be done on beta.
  # You also provided :detP,beadcount .
# [ Section 1: Check Input Done ]


# [ Section 2: Filtering Start >>

  # Filtering Detect P value Start
    # The fraction of failed positions per sample
    # You may need to delete samples with high proportion of failed probes:

                  # Failed CpG Fraction.
# 6229025012_R01C01         1.606551e-04
# 6229025012_R02C01         2.224456e-04
# 6229025012_R03C01         1.565358e-04
# 6229025012_R04C01         1.421180e-04
# 6229025012_R05C01         2.245053e-04
# 6229025012_R06C01         5.005026e-04
# 6229025012_R01C02         1.750729e-04
# 6229025012_R02C02         1.791923e-04
# 6229025012_R03C02         1.297599e-04
# 6229025012_R04C02         2.595198e-04
# 6229025012_R05C02         2.451021e-04
# 6229025012_R06C02         3.748620e-04
# 6229025035_R01C01         1.833116e-04
# 6229025035_R02C01         1.482971e-04
# 6229025035_R03C01         1.997891e-04
# 6229025035_R04C01         2.018488e-04
# 6229025035_R05C01         2.656989e-04
# 6229025035_R06C01         3.027732e-04
# 6229025035_R01C02         2.018488e-04
# 6229025035_R02C02         1.833116e-04
# 6229025035_R03C02         2.492214e-04
# 6229025035_R04C02         1.318196e-04
# 6229025035_R05C02         1.544761e-04
# 6229025035_R06C02         2.759973e-04
# 6229050114_R01C01         2.348037e-04
# 6229050114_R02C01         2.801167e-04
# 6229050114_R03C01         2.286246e-04
# 6229050114_R04C01         2.636392e-04
# 6229050114_R05C01         2.677586e-04
# 6229050114_R06C01         6.549787e-04
# 6229050114_R01C02         2.080278e-04
# 6229050114_R02C02         1.565358e-04
# 6229050114_R03C02         7.620821e-05
# 6229050114_R04C02         1.688939e-04
# 6229050114_R05C02         2.492214e-04
# 6229050114_R06C02         2.986538e-04
# 6229050144_R01C01         3.460265e-04
# 6229050144_R02C01         3.048328e-04
# 6229050144_R03C01         2.698183e-04
# 6229050144_R04C01         1.771326e-04
# 6229050144_R05C01         2.306843e-04
# 6229050144_R06C01         4.304734e-04
# 6229050144_R01C02         1.565358e-04
# 6229050144_R02C02         1.297599e-04
# 6229050144_R03C02         1.730132e-04
# 6229050144_R04C02         1.668342e-04
# 6229050144_R05C02         1.565358e-04
# 6229050144_R06C02         2.904151e-04

    # Filtering probes with a detection p-value above 0.01.
    # Removing 2152 probes.
    # If a large number of probes have been removed, ChAMP suggests you to identify potentially bad samples

  # Filtering BeadCount Start
    # Filtering probes with a beadcount <3 in at least 5% of samples.
    # Removing 0 probes

  # Filtering NoCG Start
    # Only Keep CpGs, removing 3081 probes from the analysis.

  # Filtering SNPs Start
    # Using general 450K SNP list for filtering.
    # Filtering probes with SNPs as identified in Zhou's Nucleic Acids Research Paper 2016.
    # Removing 58923 probes from the analysis.

  # Filtering MultiHit Start
    # Filtering probes that align to multiple locations as identified in Nordlund et al
    # Removing 11 probes from the analysis.

  # Filtering XY Start
    # Filtering probes located on X,Y chromosome, removing 9945 probes from the analysis.

  # Updating PD file

  # Fixing Outliers Start
    # Replacing all value smaller/equal to 0 with smallest positive value.
    # Replacing all value greater/equal to 1 with largest value below 1..
# [ Section 2: Filtering Done ]

 # All filterings are Done, now you have 411400 probes and 48 samples.


 # All filterings are Done, now you have 414676 probes and 48 samples.
 
pheno = filtered$pd
betas = filtered$beta
  
mvals = filtered$M
filtered_CpG = rownames(betas)
unfiltered_CpG = rownames(getBeta(Mset))
outCpG = unfiltered_CpG[!(unfiltered_CpG %in% filtered_CpG)]
  
champ.QC(beta = betas, pheno=pheno$expGroup, mdsPlot=TRUE, densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE, Feature.sel="None", resultsDir=here::here("QC", "QC_postFiltering"))
  
save(qc, filtered, RGset, outCpG, pheno, mvals, betas, file = here("FACT_450k_QC.RData"))

```

### Post-filtering Density 
```{r, out.width = '100%'}
knitr::include_graphics("QC/QC_postFiltering/raw_densityPlot.png")
```

### Post-filtering Beta Values
```{r, out.width = '100%'}
knitr::include_graphics("QC/QC_postFiltering/raw_mdsPlot.png")
```

### Post-filtering Dendogram
```{r, out.width = '100%'}
knitr::include_graphics("QC/QC_postFiltering/raw_SampleCluster.png")
```

### Post-Filtering SVD
```{r, warning=FALSE,eval = FALSE}
champ.SVD(beta = betas, rgSet=RGset, pd=pheno[-c(1,2,3,9,18,19,20)], RGEffect=TRUE, Rplot=FALSE, resultsDir=here::here("QC_postFiltering"))

```

```{r, out.width = '100%'}
knitr::include_graphics("QC_postFilteringSVDsummary.png")
```

### Probe-level filtering
```{r, warning=FALSE,eval = FALSE}
RGset.noob <- preprocessNoob(RGset)
save(RGset.noob,file="RGsetnoob.RData")

# plot probes density 
densityPlot(RGset, main = "density plots before and after preprocessing", pal="#440154FF", ylim=c(0,4.5))
densityPlot(RGset.noob, add = F, pal = "#FDE725FF")

# Add legend
legend("topleft", c("Noob","Raw"), lty=c(1,1), title="Normalization", bty='n', cex=1.3, col=c("#FDE725FF","#440154FF"))
 
quartz.save('noobDensity.png', type = "png", dpi = 300)
```

```{r, out.width = '100%'}
knitr::include_graphics('noobDensity.png')
```

### Estimate cell counts

```{r, warning=FALSE,eval = FALSE}
cell_counts = estimateCellCounts(RGset, compositeCellType = "Blood", referencePlatform = "IlluminaHumanMethylation450k", returnAll = TRUE, meanPlot = FALSE, verbose = TRUE)

write.csv(cell_counts$counts, file = "fact_cellcounts.csv", row.names = F)

save(cell_counts,file = "fact_cell_counts.RData")

# add cell estimates to pheno file
cellprop = cell_counts$counts
cellprop = data.frame(cellprop)
cellprop$Sample_Name = rownames(cellprop)

pheno = merge(pheno, cellprop, by = 'Sample_Name')

ids = colnames(RGset)

pheno = pheno[match(ids, pheno$Sample_Name),]
table(pheno$Sample_Name == colnames(RGset))
# TRUE 
# 48 

boxplot(cellprop[,c(1:6)]*100, col=1:ncol(cellprop),xlab="Cell type",ylab="Estimated %",main="Cell type distribution")
quartz.save('cellprop', type = "png", dpi = 300)
```

```{r, out.width = '75%'}
knitr::include_graphics('cellprop.png')
```


## Normalization

### FunNorm

```{r, warning=FALSE,eval = FALSE}

fun = preprocessFunnorm(RGset)
betas_fun = getBeta(fun)
betas_fun = betas_fun[!rownames(betas_fun) %in% outCpG,]

dim(betas_fun)
# 411400     48
betas_fun = Harman::shiftBetas(betas_fun, shiftBy=1e-4)
# No beta values were found to be 0 or 1. No shifts made.

mvals_fun <- B2M(betas_fun)

champ.QC(beta = betas_fun, pheno= pheno$expGroup, mdsPlot=TRUE, densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE, Feature.sel="None", resultsDir=here::here("postFunnorm"))

champ.SVD(beta = betas_fun, rgSet=RGset, pd=pheno[-c(1,2,3,4,10,19,20)], RGEffect=TRUE, Rplot=FALSE, resultsDir=here::here("postFunnorm"))

save(betas_fun, mvals_fun, pheno, file = here("fact_funnorm_data.RData"))
```

### Post-normalization Density 
```{r, out.width = '100%'}
knitr::include_graphics("postFunnorm/raw_densityPlot.png")
```

### Post-normalization Beta Values
```{r, out.width = '100%'}
knitr::include_graphics("postFunnorm/raw_mdsPlot.png")
```

### Post-normalization Dendogram
```{r, out.width = '100%'}
knitr::include_graphics("postFunnorm/raw_SampleCluster.png")
```

```{r, out.width = '100%'}
knitr::include_graphics("postFunnormSVDsummary.png")
```


### Quantile normalizaiton

```{r, warning=FALSE,eval = FALSE}
quan = preprocessQuantile(RGset, fixOutliers = TRUE, quantileNormalize=TRUE, stratified=TRUE)
quan = quan[!rownames(quan) %in% outCpG]
quan = mapToGenome(quan)
mvals = assays(quan)$M
betas = M2B(mvals)
newBetas = Harman::shiftBetas(betas, shiftBy=1e-4)
# No beta values were found to be 0 or 1. No shifts made.
mvals_Quantile = B2M(newBetas)
betas_Quantile = M2B(mvals_Quantile)
dim(betas_Quantile)
# 411400     48

champ.QC(beta = betas_Quantile, pheno= pheno$expGroup, mdsPlot=TRUE,
           densityPlot=TRUE, dendrogram=TRUE, PDFplot=TRUE, Rplot= FALSE,
           Feature.sel="None", resultsDir=here::here("postNorm"))

champ.SVD(beta = betas_Quantile, rgSet=RGset, pd=pheno[-c(1,2,3,4,10,19,20)], RGEffect=TRUE,
            Rplot=FALSE, resultsDir=here::here("postNorm"))

save(betas_Quantile, mvals_Quantile, pheno, file = here("fact_quantNorm_data.RData"))
```


### Post-normalization Density 
```{r, out.width = '100%'}
knitr::include_graphics("postNorm/raw_densityPlot.png")
```

### Post-normalization Beta Values
```{r, out.width = '100%'}
knitr::include_graphics("postNorm/raw_mdsPlot.png")
```

### Post-normalization Dendogram
```{r, out.width = '100%'}
knitr::include_graphics("postNorm/raw_SampleCluster.png")
```

```{r, out.width = '100%'}
knitr::include_graphics("postNormSVDsummary.png")
```


## Batch correction
```{r, warning=FALSE,eval = FALSE}

table(pheno$Sample_Name == colnames(mvals_fun))
# TRUE 
#  48 

# PCA 
fact.pca = prcomp(mvals_fun, scale = TRUE)

pairs(fact.pca$rotation[,1:5], pch = 19, cex = 0.5, col = factor(pheno$Batch))

quartz.save('pca_batch.png', type = "png", dpi = 300)
```

```{r, out.width = '100%'}
knitr::include_graphics("pca_batch.png")
```

```{r, warning=FALSE,eval = FALSE}
batch = factor(pheno$Batch)
modcombat = model.matrix(~1, data=pheno)

Mcombat = ComBat(dat = as.matrix(mvals_fun), batch = batch, mod = modcombat)

# PCA 
fact.pca.combat = prcomp(Mcombat, scale = TRUE)

pairs(fact.pca.combat$rotation[,1:5], pch = 19, cex = 0.5, col = factor(pheno$Batch))

quartz.save('pca_batch_combat.png', type = "png", dpi = 300)

save(Mcombat, pheno, file = here("fact_funnorm_combat_data.RData"))
```

```{r, out.width = '100%'}
knitr::include_graphics("pca_batch_combat.png")
```
