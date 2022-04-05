###############################################################
#                        DEPENDENCIAS                         #

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install()
}
if (!require("devtools")) {
  install.packages("devtools")
  devtools::install_github(repo = "BioinformaticsFMRP/TCGAbiolinks")
}
if (!require("SummarizedExperiment")) {
  install.packages("SummarizedExperiment", dependencies = TRUE)
  library(SummarizedExperiment)
}
if (!require("TCGAbiolinks")) {
  install.packages("TCGAbiolinks", dependencies = TRUE)
  library(TCGAbiolinks)
}
if (!require("pvclust")) {
  install.packages("pvclust", dependencies = TRUE)
  library(pvclust)
}
if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
if (!require("DT")) {
  install.packages("DT", dependencies = TRUE)
  library(DT)
}
if (!require("dplyr")) {
  install.packages("dplyr", dependencies = TRUE)
  library(dplyr)
}
if (!require("maftools")) {
  install.packages("maftools", dependencies = TRUE)
  library(maftools)
}
if (!require("downloader")) {
  install.packages("downloader", dependencies = TRUE)
  library(downloader)
}
if (!require("readr")) {
  install.packages("readr", dependencies = TRUE)
  library(readr)
}
if (!require("gaia")) {
  install.packages("gaia", dependencies = TRUE)
  library(gaia)
}

setwd("/Users/martielafreitas/Documents/Colaboracoes/Cristal/HG19")

###############################################################
#             ANALISE DE EXPRESSAO DIFERENCIAL                #

# Buscando dados LGG
query.lgg <- GDCquery(project = "TCGA-LGG", 
                          legacy = TRUE,
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq", 
                          file.type = "results")
GDCdownload(query.lgg)
data.lgg <- GDCprepare(query = query.lgg, save = TRUE, save.filename = "lggExp.rda")

# Buscando dados GBM
query.gbm <- GDCquery(project = "TCGA-GBM", 
                          legacy = TRUE,
                          data.category = "Gene expression",
                          data.type = "Gene expression quantification",
                          platform = "Illumina HiSeq", 
                          file.type = "results")
GDCdownload(query.gbm)
data.gbm <- GDCprepare(query = query.gbm, save = TRUE, save.filename = "gbmExp.rda")

# Unindo dados brutos LGG e GBM
data.gbm.lgg <- SummarizedExperiment::cbind(data.lgg, data.gbm)

rowRanges(data.gbm)
rowRanges(data.lgg)
rowRanges(data.gbm.lgg)

# ------------- CRIANDO MATRIZ - RAW DATA --------------------
# cria matriz a partir de: data
dataMatrix.gbm <- assay(data.gbm,"raw_count")
dataMatrix.lgg <- assay(data.lgg,"raw_count")
dataMatrix.gbm.lgg <- assay(data.gbm.lgg,"raw_count")

# ------------------ SALVANDO RAW DATA -----------------------
#save raw_data
write.table(dataMatrix.gbm, 'gbm_raw_dataMatrix.txt',sep='\t')
write.table(dataMatrix.lgg, 'lgg_raw_dataMatrix.txt',sep='\t')
write.table(dataMatrix.gbm.lgg, 'gbm_lgg_raw_dataMatrix.txt',sep='\t')

# ------------------------ OUTLIERS --------------------------
# remove outliers
data_CorOutliers_gbm <- TCGAanalyze_Preprocessing(data.gbm)
data_CorOutliers_lgg <- TCGAanalyze_Preprocessing(data.lgg)
data_CorOutliers_gbm_lgg <- TCGAanalyze_Preprocessing(data.gbm.lgg)

# ------------------- SALVANDO TABELAS -----------------------
# Downstream analysis using gene expression data  
# TCGA samples from IlluminaHiSeq_RNASeqV2 with type rsem.genes.results
save(data.gbm, geneInfo , file = "dataGeneExpressionGBM.rda")
save(data.lgg, geneInfo , file = "dataGeneExpressionLGG.rda")
save(data.gbm.lgg, geneInfo , file = "dataGeneExpressionGBM_LGG.rda")


# ----------------------- NORMALIZACAO -----------------------
# normalization of genes using dataMatrix
# apenas essa, para garantir que nao vamos perder
# nenhum gene. Por GC content excluiriamos varios genes.
dataNorm.gbm <- TCGAanalyze_Normalization(tabDF = dataMatrix.gbm,
                                          geneInfo,
                                          method = "geneLength")
dataNorm.lgg <- TCGAanalyze_Normalization(tabDF = dataMatrix.lgg,
                                          geneInfo,
                                          method = "geneLength")
dataNorm.gbm.lgg <- TCGAanalyze_Normalization(tabDF = dataMatrix.gbm.lgg,
                                              geneInfo,
                                              method = "geneLength")

# ---------------- FILTRAGEM POR QUARTIL ----------------------
# quantile filter of genes using normalized data
# general expression: differential and non differential expressed
dataFilt.gbm <- TCGAanalyze_Filtering(tabDF = dataNorm.gbm,
                                      method = "quantile", 
                                      qnt.cut =  0.25)
dataFilt.lgg <- TCGAanalyze_Filtering(tabDF = dataNorm.lgg,
                                      method = "quantile", 
                                      qnt.cut =  0.25)
dataFilt.lgg.blood <- TCGAanalyze_Filtering(tabDF = dataNorm.lgg,
                                      method = "quantile", 
                                      qnt.cut =  0.75)
dataFilt.gbm.lgg <- TCGAanalyze_Filtering(tabDF = dataNorm.gbm.lgg,
                                      method = "quantile", 
                                      qnt.cut =  0.25)

# ------------------------- NORMAL --------------------------
# normal tissue sample
samplesNT.gbm <- TCGAquery_SampleTypes(barcode = colnames(dataFilt.gbm),
                                   typesample = c("NT"))
samplesNT.lgg <- TCGAquery_SampleTypes(barcode = colnames(dataFilt.lgg),
                                       typesample = c("NT")) 
samplesNT.gbm.lgg <- TCGAquery_SampleTypes(barcode = colnames(dataFilt.gbm.lgg),
                                       typesample = c("NT"))

# -------------------------- TUMOR --------------------------
# tumoral tissue sample
samplesTP.gbm <- TCGAquery_SampleTypes(barcode = colnames(dataFilt.gbm), 
                                       typesample = c("TP"))
samplesTP.lgg <- TCGAquery_SampleTypes(barcode = colnames(dataFilt.lgg), 
                                       typesample = c("TP"))
samplesTP.gbm.lgg <- TCGAquery_SampleTypes(barcode = colnames(dataFilt.gbm.lgg), 
                                       typesample = c("TP"))


# ------------------------- DEAll ----------------------------
# Esse nao deve encontrar diferencas, pois os parametros
# FDR e logFC nao irao filtrar nada. Apenas teremos um panorama 
# geral da expressao genica.
# usa pra o volcano plot
# FDR = 2 pra fazer o que nao fiz solteira: PEGAR GERAL. \o/
#
# ----------------------- DEAll GBM --------------------------
#  there are Cond1 type Normal in  5 samples
#  there are Cond2 type Tumor in  156 samples
#  there are  14893 features as miRNA or genes
dataDEGs_GBM_allgenes <- TCGAanalyze_DEA(mat1 = dataFilt.gbm[,samplesNT.gbm],
                                         mat2 = dataFilt.gbm[,samplesTP.gbm],
                                         Cond1type = "Normal",
                                         Cond2type = "Tumor",
                                         fdr.cut = 2 ,
                                         logFC.cut = 0,
                                         method = "glmLRT")

# --------------------- DEAll LGG ---------------------------
#  there are Cond1 type Normal in  0 samples
#  there are Cond2 type Tumor in  516 samples
#  there are  14893 features as miRNA or genes
dataDEGs_LGG_allgenes <- TCGAanalyze_DEA(mat1 = dataFilt.lgg[,samplesNT.lgg],
                                         mat2 = dataFilt.lgg[,samplesTP.lgg],
                                         Cond1type = "Normal",
                                         Cond2type = "Tumor",
                                         fdr.cut = 2 ,
                                         logFC.cut = 0,
                                         method = "glmLRT")
# OBS: Dá erro por não ter amostra de tecido normal.
# Rodar de novo com amostra de sangue(?)

# ------------------- DEAll GBM e LGG ------------------------
#  there are Cond1 type Normal in  5 samples
#  there are Cond2 type Tumor in  672 samples
#  there are  14893 features as miRNA or genes 
dataDEGs_GBM_LGG_allgenes <- TCGAanalyze_DEA(mat1 = dataFilt.gbm.lgg[,samplesNT.gbm.lgg],
                                         mat2 = dataFilt.gbm.lgg[,samplesTP.gbm.lgg],
                                         Cond1type = "Normal",
                                         Cond2type = "Tumor",
                                         fdr.cut = 2 ,
                                         logFC.cut = 0,
                                         method = "glmLRT")
# OBS: NÃO PRECISA continuar a análise com LGG + GBM.
# A partir daqui vamos usar apenas GBM, pois LGG não
# tem amostras normais no HG19. Nem para bloodinho.

# -------------------------- DEA -----------------------------
#
# Esse sim deve achar DE. E MUITA!
# Usa para enriquecimento.
#
# ------------------------ DEA GBM ---------------------------
#  there are Cond1 type Normal in  5 samples
#  there are Cond2 type Tumor in  156 samples
#  there are  14893 features as miRNA or genes
dataDEGs_GBM <- TCGAanalyze_DEA(mat1 = dataFilt.gbm[,samplesNT.gbm],
                                mat2 = dataFilt.gbm[,samplesTP.gbm],
                                Cond1type = "Normal",
                                Cond2type = "Tumor",
                                fdr.cut = 0.01 ,
                                logFC.cut = 2,
                                method = "glmLRT")
# RESERVA!
dataDEGs_GBM_allgenes_um <- TCGAanalyze_DEA(mat1 = dataFilt.gbm[,samplesNT.gbm],
                                      mat2 = dataFilt.gbm[,samplesTP.gbm],
                                      Cond1type = "Normal",
                                      Cond2type = "Tumor",
                                      fdr.cut = 0.999 ,
                                      logFC.cut = 0,
                                      method = "glmLRT")

# ------------------------ GRAFICOS GBM ---------------------------
#
# DEAll GBM
pdf("volcanoplot_GBM_allgenes.pdf")
with(dataDEGs_GBM_allgenes,
     plot(logFC, -log10(FDR), pch=20, cex=0.5, main="Azul LogFC+TP -- Vermelho LogFC-NT", cex.main=0.5))
with(subset(dataDEGs_GBM_allgenes,((logFC>2) & (FDR <0.01))),
     points(logFC, -log10(FDR), pch=20, cex=0.3, col="blue"))
with(subset(dataDEGs_GBM_allgenes,((logFC<(-2)) & (FDR <0.01))),
     points(logFC, -log10(FDR), pch=20, cex=0.3, col="red"))
abline(v=2,col="gray60")
abline(v=-2,col="gray60")
abline(h=(-1*(log10(0.01))), col="purple")
dev.off()

# DEA GBM
#volcano
pdf("volcanoplot_GBM.pdf")
with(dataDEGs_GBM,
     plot(logFC, -log10(FDR), pch=20, cex=0.5, main="Azul LogFC+TP -- Vermelho LogFC-NT", cex.main=0.5))
with(subset(dataDEGs_GBM, ((logFC>2) & (FDR <0.05))),
     points(logFC, -log10(FDR), pch=20, cex=0.3, col="blue"))
with(subset(dataDEGs_GBM, ((logFC<(-2)) & (FDR <0.05))),
     points(logFC, -log10(FDR), pch=20, cex=0.3, col="red"))
abline(v=2,col="gray60")
abline(v=-2,col="gray60")
abline(h=(-1*(log10(0.05))),col="purple")
dev.off()

# ------------------- SALVANDO TABELAS -----------------------
#
write.table(dataDEGs_GBM_allgenes, 'dataDEGs_GBM_allgenes.txt', sep='\t')
write.table(dataDEGs_GBM, 'dataDEGs_GBM.txt', sep='\t')
write.table(dataDEGs_GBM_allgenes_um, 'dataDEGs_GBM_allgenes_fdr999.txt', sep='\t')

###############################################################
#                 ANALISE DE ENRIQUECIMENTO                   #

# Usa a tabela dos diferencialmente expressos
# DEGs table with expression values in normal and tumor samples
dataDEGsFiltLevel_GBM <- TCGAanalyze_LevelTab(dataDEGs_GBM,"Tumor","Normal",
                                          dataFilt.gbm[,samplesTP.gbm],
                                          dataFilt.gbm[,samplesNT.gbm])

write.table(dataDEGsFiltLevel_GBM, 'dataDEGsFiltLevel_GBM.txt', sep='\t')

# Enrichment Analysis EA
# Gene Ontology (GO) and Pathway enrichment by DEGs list
Genelist.GBM <- rownames(dataDEGsFiltLevel_GBM)

system.time(ansEA.gbm <- TCGAanalyze_EAcomplete(TFname="DEA GBM genes Normal Vs Tumor",
                                                Genelist.GBM))

# Enrichment Analysis EA (TCGAVisualize)
# Gene Ontology (GO) and Pathway enrichment barPlot
TCGAvisualize_EAbarplot(tf = rownames(ansEA.gbm$ResBP), 
                        GOBPTab = ansEA.gbm$ResBP,
                        GOCCTab = ansEA.gbm$ResCC,
                        GOMFTab = ansEA.gbm$ResMF,
                        PathTab = ansEA.gbm$ResPat,
                        nRGTab = Genelist.GBM, 
                        nBar = 10)

###PCA plot
pdf('pca_top_50_miRNA.pdf')
TCGAvisualize_PCA(dataFilt,
                  dataDEGsFiltLevel,
                  ntopgenes = 50,
                  colnames(dataFilt[,samplesTP]),
                  colnames(dataFilt[,samplesNT]))
###azul ficou TP e vermelho NT como no volcano
dev.off()

# Se der erro precisa excluir os TM
# descomentar a linha abaixo e rodar.
# dataFilt2 <- dataFilt[,!(colnames(dataFilt)%in%samplesTM)]

### PCA plot -> apos remocao do TM
# descomentar as linhas abaixo e rodar.
#pdf('pca_top_50.pdf')
#TCGAvisualize_PCA(dataFilt2,
#                  dataDEGsFiltLevel,
#                  ntopgenes = 50,
#                  colnames(dataFilt2[,samplesTP]),
#                  colnames(dataFilt2[,samplesNT])) ### azul ficou TP e vermelho NT como no volcano
#dev.off()

### PCA plot -> apos remocao do TM
# descomentar as linhas abaixo e rodar.
#pdf('pca_top_100.pdf')
#TCGAvisualize_PCA(dataFilt2,
#                  dataDEGsFiltLevel,
#                  ntopgenes = 100,
#                  colnames(dataFilt2[,samplesTP]),
#                  colnames(dataFilt2[,samplesNT])) ###azul ficou TP e vermelho NT como no volcano
#dev.off()

### Salvando tabelas -> apos remocao do TM
# descomentar as linhas abaixo e rodar.
#write.table(dataFilt2,"samplesNT_TP.txt", sep="\t")
#write.table(dataDEGs_001,"difexp_allgenes.txt", sep="\t")

###############################################################
#                       ANALISE DE SNVs                       #

query.maf.hg19 <- GDCquery(project = "TCGA-GBM", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           legacy = TRUE)

datatable(select(getResults(query.maf.hg19),-contains("cases")),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 10), 
          rownames = FALSE)

query.maf.hg19 <- GDCquery(project = "TCGA-GBM", 
                           data.category = "Simple nucleotide variation", 
                           data.type = "Simple somatic mutation",
                           access = "open", 
                           file.type = "gsc_GBM_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf",
                           legacy = TRUE)

GDCdownload(query.maf.hg19,
            directory = "/Users/martielafreitas/Documents/Colaboracoes/Cristal/HG19")

maf <- GDCprepare(query.maf.hg19,
                directory = "/Users/martielafreitas/Documents/Colaboracoes/Cristal/HG19")

datatable(maf[1:50,],
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

library(maftools)
library(tibble)

maf2 <- read.maf(maf)

datatable(getSampleSummary(maf2),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)

pdf('variant_summary_data.pdf')
plotmafSummary(maf = maf2, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)
dev.off()

pdf('oncoplot_variant_data.pdf')
oncoplot(maf = maf2, top = 10, removeNonMutated = TRUE)
dev.off()

titv = titv(maf = maf2, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)