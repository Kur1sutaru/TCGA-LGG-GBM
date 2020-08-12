# TCGA-LGG-GBM
################################################################################################
# Single Nucleotide Variation Analysis - TCGA - Lower Grade Glioma and Glioblastoma Multiforme #
################################################################################################

setwd("~/path_to_the_directory")

#Load dependencies

library(TCGAbiolinks)
library(pvclust)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(DT)
library(readr)
library(SummarizedExperiment)
library(maftools)
library(downloader)
library(readr)
library(gaia)
library(GenVisR)


query.maf.hg19 <- GDCquery(project = "TCGA-GBM", 
                           
                           data.category = "Simple nucleotide variation", 
                           
                           data.type = "Simple somatic mutation",
                           
                           access = "open", 
                           
                           legacy = TRUE)
                           
                           
 # Check maf availables
datatable(dplyr::select(getResults(query.maf.hg19),-contains("cases")),
                           filter = 'top',
                           options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
                           rownames = FALSE)
                           
                           
query.maf.hg19 <- GDCquery(project = "TCGA-GBM", 
                           data.category = "Simple nucleotide variation",
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "gsc_GBM_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf",
                           legacy = TRUE)

GDCdownload(query.maf.hg19)


maf <- GDCprepare(query.maf.hg19)


# Only first 50 to make render faster

datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


library(maftools)

library(dplyr)
maf <- GDCquery_Maf("GBM", pipelines = "muse") %>% read.maf

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

pdf('oncoplot_variant_GBM_topgenes.pdf')
oncoplot(maf = maf, top = 30, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)
dev.off()


#######################################################
#                                                     #
#   Lysosomal storage related genes                   #
#                                                     #
#######################################################

pdf('oncoplot_variant_GBM_topgenes_genes.pdf')
oncoplot(maf = maf, 
         genes = c('AGA', 'ARSA', 'ARSB', 'ASAH1', 'CLN3', 'CTNS', 'CTSA', 'CTSK', 'FUCA1', 'GAA', 'GALC', 'GALNS', 'GBA', 'GLA', 'GLB1', 'GM2A', 'GNPTAB', 'GNPTG', 'GNS', 'GUSB', 
                   'HEXA', 'HEXB', 'HGSNAT', 'HYAL1', 'IDS', 'IDUA', 'LAMP2', 'LIPA', 'MAN2B1', 'MANBA', 'MCOLN1', 'NAGA', 'NAGLU', 'NEU1', 'NPC1', 'NPC2', 'PPT1', 'PSAP', 'SGSH', 'SMPD1', 'SUMF1', 'TPP1'),
         removeNonMutated= TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)
dev.off()




#######################################################
#                                                     #
#   Other genes to compare                            #
#                                                     #
#######################################################

pdf('oncoplot_variant_GBM_topgenes_genes_oxphosp.pdf')
oncoplot(maf = maf, 
         genes = c('Chchd10', 'MTND1', 'MTCO2', 'COX14', 'COX10', 'FXN', 'MSH', 'Nipsnap2', 'SURF1', 'TEFM', 'MSH2'),
                   removeNonMutated= TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)
dev.off()

query.maf.hg19 <- GDCquery(project = "TCGA-LGG", 
                           
                           data.category = "Simple nucleotide variation", 
                           
                           data.type = "Simple somatic mutation",
                           
                           access = "open", 
                           
                           legacy = TRUE)


# Check maf availables
datatable(dplyr::select(getResults(query.maf.hg19),-contains("cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
          rownames = FALSE)

# URL with all the curated .maf files in the TCGA cohorts
# https://docs.cancergenomicscloud.org/discuss/572d0663b2c3d1170088fabd



query.maf.hg19 <- GDCquery(project = "TCGA-LGG", 
                           data.category = "Simple nucleotide variation",
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "LGG_FINAL_ANALYSIS.aggregated.capture.tcga.uuid.curated.somatic.maf",
                           legacy = TRUE)

GDCdownload(query.maf.hg19)


maf <- GDCprepare(query.maf.hg19)


# Only first 50 to make render faster

datatable(maf[1:20,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)


library(maftools)
library(dplyr)

# Other example - https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/DNA_Seq_Variant_Calling_Pipeline/
maf <- GDCquery_Maf("LGG", pipelines = "muse") %>% read.maf

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, showBarcodes = FALSE, top = 20)

pdf('oncoplot_variant_LGG_FINAL_ANALYSIS_topgenes.pdf')
oncoplot(maf = maf, top = 30, removeNonMutated = TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)
dev.off()


#######################################################
#                                                     #
#   LSD related genes                                 #
#                                                     #
#######################################################

pdf('oncoplot_variant_LGG_FINAL_ANALYSIS_lysogenes.pdf')
oncoplot(maf = maf, 
         genes = c('AGA', 'ARSA', 'ARSB', 'ASAH1', 'CLN3', 'CTNS', 'CTSA', 'CTSK', 'FUCA1', 'GAA', 'GALC', 'GALNS', 'GBA', 'GLA', 'GLB1', 'GM2A', 'GNPTAB', 'GNPTG', 'GNS', 'GUSB', 
                   'HEXA', 'HEXB', 'HGSNAT', 'HYAL1', 'IDS', 'IDUA', 'LAMP2', 'LIPA', 'MAN2B1', 'MANBA', 'MCOLN1', 'NAGA', 'NAGLU', 'NEU1', 'NPC1', 'NPC2', 'PPT1', 'PSAP', 'SGSH', 'SMPD1', 'SUMF1', 'TPP1'),
         removeNonMutated= TRUE)
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)
dev.off()




############################################################
#                                                          #
#                                                          #
# Draw two oncoplots side by side for cohort comparision   #
#                                                          #
############################################################


genes = c('AGA', 'ARSA', 'ARSB', 'ASAH1', 'CLN3', 'CTNS', 'CTSA', 'CTSK', 'FUCA1', 'GAA', 'GALC', 'GALNS', 'GBA', 'GLA', 'GLB1', 'GM2A', 'GNPTAB', 'GNPTG', 'GNS', 'GUSB', 
          'HEXA', 'HEXB', 'HGSNAT', 'HYAL1', 'IDS', 'IDUA', 'LAMP2', 'LIPA', 'MAN2B1', 'MANBA', 'MCOLN1', 'NAGA', 'NAGLU', 'NEU1', 'NPC1', 'NPC2', 'PPT1', 'PSAP', 'SGSH', 'SMPD1', 'SUMF1', 'TPP1')
coOncoplot(m1 = maf2gbm, m2 = maf3, m1Name = 'gbm', m2Name = 'lgg', genes = genes, removeNonMutated = TRUE)




coOncoplot(maf2gbm, maf3, m1Name = 'GBM', m2Name = 'LGG',
           showSampleNames = TRUE)

#########################################################################
#                                                                       #
# MULTIVARIATE SURVIVAL ANALYSIS                                        #
#                                                                       #
# I WILL SURVIVE!                                                       #
########################################################################

library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)

# More info https://cran.r-project.org/web/packages/survivalAnalysis/vignettes/multivariate.html

library(monocle)
biocLite(c("DDRTree", "pheatmap"))

# helper function to set the root state correctly 
get_correct_root_state <- function(cds){
  T0_counts <- table(pData(cds)$State, pData(cds)$Type)[,"Lsk"]
  as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])
}


# To perform Basic Differential Analysis in Monocle, please visit the url: http://cole-trapnell-lab.github.io/monocle-release/docs/#differential-expression-analysis


suppressMessages(library(monocle))
suppressMessages(library(stringr))
suppressMessages(library(plyr))
suppressMessages(library(netbiov))
