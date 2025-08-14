#open seurat
library(dplyr)
library(Seurat)
library(patchwork)

#load mtx, import, edit and add metadata (see Seurat_20241226.rtf)

# check files
readLines("lung_metaData.txt", n=10)

# open scp1219.RData
load("~/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/20240610/Projects/Pdlim2 in lung infection/scRNAseq/BroadInstitute_SingleCellPortal/SCP1219_GSE171524/scp1219.RData")

# subset object with PDLIM2>0
scp1219_PDLIM2 <- subset (scp1219, subset = PDLIM2 >0)

#random check subset results
scp1219_PDLIM2@assays$RNA@data ["PDLIM2", 1:6]

#colnames of metadata
colnames(scp1219_PDLIM2@meta.data)

# FetchData with vars of interest
scp1219_PDLIM2_metadata <- FetchData (scp1219_PDLIM2, vars = c("donor_id", "disease", "group", "cell_type_main", "cell_type_intermediate", "cell_type_fine", "intubation_days", "interval_death_symptoms_onset_days", "pmi_h", "PDLIM2"))

#write.csv
write.csv (scp1219_PDLIM2_metadata, file= "scp1219_PDLIM2_metadata.csv")

# Install Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install MAST
BiocManager::install("MAST")

# load MAST
library(MAST)

# DEG All with MAST
DEG <- FindMarkers(scp1219_PDLIM2,
ident.1 = "Control",
ident.2 = "COVID-19",
group.by = "group",
test.use = "MAST")  # Default is Wilcoxon test

# write.csv
write.csv (DEG, file= "DEG_All_MAST.csv")

#DEG All with wilcox
DEG <- FindMarkers(scp1219_PDLIM2,
                   ident.1 = "Control",
                   ident.2 = "COVID-19",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test

# write.csv
write.csv (DEG, file= "DEG_All_wilcox.csv")


# subset AT1 for DEG
scp1219_PDLIM2_AT1 <- subset(scp1219_PDLIM2, subset = cell_type_intermediate == "AT1")

# DEG AT1 with MAST
DEG <- FindMarkers(scp1219_PDLIM2_AT1,
                       ident.1 = "Control",
                       ident.2 = "COVID-19",
                       group.by = "group",
                       test.use = "MAST")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT1_MAST.csv")

# DEG AT1 with wilcox
DEG <- FindMarkers(scp1219_PDLIM2_Fibroblasts,
                       ident.1 = "COVID-19",
                       ident.2 = "Control",
                       group.by = "group",
                       test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Fibroblast_Covid-Ctrl.csv")

# subset AT2 for DEG
scp1219_PDLIM2_AT2 <- subset(scp1219_PDLIM2, subset = cell_type_intermediate == "AT2")

# DEG AT2
DEG <- FindMarkers(scp1219_PDLIM2_AT2,
                       ident.1 = "Control",
                       ident.2 = "COVID-19",
                       group.by = "group",
                       test.use = "MAST")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT2_MAST.csv")

# DEG AT2
DEG <- FindMarkers(scp1219_PDLIM2_AT2,
                       ident.1 = "Control",
                       ident.2 = "COVID-19",
                       group.by = "group",
                       test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT2_wilcox.csv")

# subset Fibroblasts for DEG
scp1219_PDLIM2_Fibroblasts <- subset(scp1219_PDLIM2, subset = cell_type_intermediate == "Fibroblasts")

# DEG Fibroblasts
DEG <- FindMarkers(scp1219_PDLIM2_Fibroblasts,
                   ident.1 = "Control",
                   ident.2 = "COVID-19",
                   group.by = "group",
                   test.use = "MAST")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Fibroblasts_MAST.csv")

# DEG Fibroblasts
DEG <- FindMarkers(scp1219_PDLIM2_Fibroblasts,
                   ident.1 = "Control",
                   ident.2 = "COVID-19",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Fibroblasts_wilcox.csv")

# subset "Airway epithelial cells" for DEG
scp1219_PDLIM2_AirwayEpi <- subset(scp1219_PDLIM2, subset = cell_type_intermediate == "Airway epithelial cells")

# DEG AirwayEpi
DEG <- FindMarkers(scp1219_PDLIM2_AirwayEpi,
                   ident.1 = "Control",
                   ident.2 = "COVID-19",
                   group.by = "group",
                   test.use = "MAST")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AirwayEpi_MAST.csv")

# DEG AirwayEpi
DEG <- FindMarkers(scp1219_PDLIM2_AirwayEpi,
                   ident.1 = "Control",
                   ident.2 = "COVID-19",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AirwayEpi_wilcox.csv")

# subset "Macrophages" for DEG
scp1219_PDLIM2_Macrophages <- subset(scp1219_PDLIM2, subset = cell_type_intermediate == "Macrophages")

# DEG Macrophages
DEG <- FindMarkers(scp1219_PDLIM2_Macrophages,
                   ident.1 = "Control",
                   ident.2 = "COVID-19",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Macrophages_wilcox.csv")

# test subset - worked
test1 <- subset(scp1219, subset = cell_type_intermediate %in% c("AT1", "AT2"))

test1 <- subset(scp1219, subset = cell_type_intermediate %in% c("AT1", "AT2") & PDLIM2>0)

# save RData
save (list=ls(), file="GSE171524.RData")

# open GSE171524.RData
load("GSE171524.RData")

 # classify PDLIM2+ AT1 to High/Low level Pdlim2 groups
# Fetch PDLIM2 expression data and calculate its level at certain percentile
# Choose the gene of interest
gene_of_interest <- "PDLIM2"  # Replace with your gene

# Extract expression levels
gene_expr <- FetchData(scp1219_PDLIM2_AT1, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
scp1219_PDLIM2_AT1$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(scp1219_PDLIM2_AT1$expression_group)

colnames(scp1219_PDLIM2_AT1[[]])

table(scp1219_PDLIM2_AT1[[]]$group)

# Fetch above High/Low PDLIM2 AT1 groups
scp1219_PDLIM2_AT1_PDLIM2LowHigh <- FetchData(object = scp1219_PDLIM2_AT1, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AT1_PDLIM2LowHigh, file="scp1219_PDLIM2_AT1_PDLIM2LowHigh.csv")

#check active.ident
head(scp1219_PDLIM2_AT1@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AT1 <- SetIdent(scp1219_PDLIM2_AT1, value = scp1219_PDLIM2_AT1@meta.data$expression_group)

# DEG AT1 with High vs Low on PDLIM2 level
DEG <- FindMarkers(scp1219_PDLIM2_AT1,
                   ident.1 = "Low",
                   ident.2 = "High",
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT1_PDLIM2_Low-High.csv")

## AT2
# Extract expression levels
gene_expr <- FetchData(scp1219_PDLIM2_AT2, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
scp1219_PDLIM2_AT2$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(scp1219_PDLIM2_AT2$expression_group)

colnames(scp1219_PDLIM2_AT2[[]])

table(scp1219_PDLIM2_AT2[[]]$group)

# Fetch above High/Low PDLIM2 AT2 groups
scp1219_PDLIM2_AT2_PDLIM2LowHigh <- FetchData(object = scp1219_PDLIM2_AT2, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AT2_PDLIM2LowHigh, file="scp1219_PDLIM2_AT2_PDLIM2LowHigh.csv")

#check active.ident
head(scp1219_PDLIM2_AT2@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AT2 <- SetIdent(scp1219_PDLIM2_AT2, value = scp1219_PDLIM2_AT2@meta.data$expression_group)

# DEG AT2 with High vs Low on PDLIM2 level
DEG <- FindMarkers(scp1219_PDLIM2_AT2,
                   ident.1 = "High",
                   ident.2 = "Low",
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT2_PDLIM2_High-Low.csv")

## Airway Epi
# Extract expression levels
gene_expr <- FetchData(scp1219_PDLIM2_AirwayEpi, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
scp1219_PDLIM2_AirwayEpi$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(scp1219_PDLIM2_AirwayEpi$expression_group)

colnames(scp1219_PDLIM2_AirwayEpi[[]])

table(scp1219_PDLIM2_AirwayEpi[[]]$group)

# Fetch above High/Low PDLIM2 AirwayEpi groups
scp1219_PDLIM2_AirwayEpi_PDLIM2LowHigh <- FetchData(object = scp1219_PDLIM2_AirwayEpi, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AirwayEpi_PDLIM2LowHigh, file="scp1219_PDLIM2_AirwayEpi_PDLIM2LowHigh.csv")

#check active.ident
head(scp1219_PDLIM2_AirwayEpi@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AirwayEpi <- SetIdent(scp1219_PDLIM2_AirwayEpi, value = scp1219_PDLIM2_AirwayEpi@meta.data$expression_group)

# DEG AirwayEpi with High vs Low on PDLIM2 level
DEG <- FindMarkers(scp1219_PDLIM2_AirwayEpi,
                   ident.1 = "High",
                   ident.2 = "Low",
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AirwayEpi_PDLIM2_High-Low.csv")

## Fibroblast
# Extract expression levels
gene_expr <- FetchData(scp1219_PDLIM2_Fibroblasts, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
scp1219_PDLIM2_Fibroblasts$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(scp1219_PDLIM2_Fibroblasts$expression_group)

colnames(scp1219_PDLIM2_Fibroblasts[[]])

table(scp1219_PDLIM2_Fibroblasts[[]]$group)

# Fetch above High/Low PDLIM2 AirwayEpi groups
scp1219_PDLIM2_Fibroblasts_PDLIM2LowHigh <- FetchData(object = scp1219_PDLIM2_Fibroblasts, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_Fibroblasts_PDLIM2LowHigh, file="scp1219_PDLIM2_Fibroblasts_PDLIM2LowHigh.csv")

#check active.ident
head(scp1219_PDLIM2_Fibroblasts@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_Fibroblasts <- SetIdent(scp1219_PDLIM2_Fibroblasts, value = scp1219_PDLIM2_Fibroblasts@meta.data$expression_group)

# DEG Fibroblasts with High vs Low on PDLIM2 level
DEG <- FindMarkers(scp1219_PDLIM2_Fibroblasts,
                   ident.1 = "Low",
                   ident.2 = "High",
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Fibroblasts_PDLIM2_Low-High.csv")

# load data
load("GSE171524.RData")

colnames(scp1219_PDLIM2_AT1[[]])
table(scp1219_PDLIM2_AT1[[]]$group)

##scp1219_PDLIM2_AT1_Covid
# subset scp1219_PDLIM2_AT1_Covid
scp1219_PDLIM2_AT1_Covid <- subset(scp1219_PDLIM2_AT1, subset = group %in% c("COVID-19"))

table(scp1219_PDLIM2_AT1_Covid[[]]$group)

# Extract expression levels
gene_of_interest <- "PDLIM2"
gene_expr <- FetchData(scp1219_PDLIM2_AT1_Covid, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
scp1219_PDLIM2_AT1_Covid$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(scp1219_PDLIM2_AT1_Covid$expression_group)

colnames(scp1219_PDLIM2_AT1_Covid[[]])

table(scp1219_PDLIM2_AT1_Covid[[]]$group)

# Fetch above High/Low PDLIM2 AT1_Covid groups
scp1219_PDLIM2_AT1_Covid_PDLIM2LowHigh <- FetchData(object = scp1219_PDLIM2_AT1_Covid, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AT1_Covid_PDLIM2LowHigh, file="scp1219_PDLIM2_AT1_Covid_PDLIM2LowHigh.csv")

#check active.ident
head(scp1219_PDLIM2_AT1_Covid@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AT1_Covid <- SetIdent(scp1219_PDLIM2_AT1_Covid, value = scp1219_PDLIM2_AT1_Covid@meta.data$expression_group)

# DEG AT1_Covid with High vs Low on PDLIM2 level
DEG <- FindMarkers(scp1219_PDLIM2_AT1_Covid,
                   ident.1 = "Low",
                   ident.2 = "High",
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT1_Covid_PDLIM2_Low-High.csv")

## scp1219_PDLIM2_AT1_Ctrl
# subset scp1219_PDLIM2_AT1_Ctrl
table(scp1219_PDLIM2_AT1[[]]$group)
scp1219_PDLIM2_AT1_Ctrl <- subset(scp1219_PDLIM2_AT1, subset = group %in% c("Control"))
table(scp1219_PDLIM2_AT1_Ctrl[[]]$group)

# Extract expression levels
gene_of_interest <- "PDLIM2"
gene_expr <- FetchData(scp1219_PDLIM2_AT1_Ctrl, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
scp1219_PDLIM2_AT1_Ctrl$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(scp1219_PDLIM2_AT1_Ctrl$expression_group)

colnames(scp1219_PDLIM2_AT1_Ctrl[[]])

table(scp1219_PDLIM2_AT1_Ctrl[[]]$group)

# Fetch above High/Low PDLIM2 AT1_Ctrl groups
scp1219_PDLIM2_AT1_Ctrl_PDLIM2LowHigh <- FetchData(object = scp1219_PDLIM2_AT1_Ctrl, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AT1_Ctrl_PDLIM2LowHigh, file="scp1219_PDLIM2_AT1_Ctrl_PDLIM2LowHigh.csv")

#check active.ident
head(scp1219_PDLIM2_AT1_Ctrl@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AT1_Ctrl <- SetIdent(scp1219_PDLIM2_AT1_Ctrl, value = scp1219_PDLIM2_AT1_Ctrl@meta.data$expression_group)

# DEG AT1_Covid with High vs Low on PDLIM2 level
DEG <- FindMarkers(scp1219_PDLIM2_AT1_Ctrl,
                   ident.1 = "Low",
                   ident.2 = "High",
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT1_Ctrl_PDLIM2_Low-High.csv")

## scp1219_PDLIM2_AT1_Pdlim2Low
# subset scp1219_PDLIM2_AT1_Pdlim2Low
colnames(scp1219_PDLIM2_AT1[[]])
table(scp1219_PDLIM2_AT1[[]]$expression_group)
scp1219_PDLIM2_AT1_Pdlim2Low <- subset(scp1219_PDLIM2_AT1, subset = expression_group %in% c("Low"))
table(scp1219_PDLIM2_AT1_Pdlim2Low[[]]$group)

# Fetch from above seurat object Pdlim2 expression from Covid and Ctrl groups
scp1219_PDLIM2_AT1_Pdlim2Low_PDLIM2 <- FetchData(object = scp1219_PDLIM2_AT1_Pdlim2Low, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AT1_Pdlim2Low_PDLIM2, file="scp1219_PDLIM2_AT1_Pdlim2Low_PDLIM2.csv")

#check active.ident
head(scp1219_PDLIM2_AT1_Pdlim2Low@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AT1_Pdlim2Low <- SetIdent(scp1219_PDLIM2_AT1_Pdlim2Low, value = scp1219_PDLIM2_AT1_Pdlim2Low@meta.data$group)

# DEG AT1_Pdlim2Low_Covid-Ctrl
DEG <- FindMarkers(scp1219_PDLIM2_AT1_Pdlim2Low,
                   ident.1 = "COVID-19",
                   ident.2 = "Control",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT1_Pdlim2Low_Covid-Ctrl.csv")

## scp1219_PDLIM2_AT1_Pdlim2High
# subset scp1219_PDLIM2_AT1_Pdlim2High
colnames(scp1219_PDLIM2_AT1[[]])
table(scp1219_PDLIM2_AT1[[]]$expression_group)
scp1219_PDLIM2_AT1_Pdlim2High <- subset(scp1219_PDLIM2_AT1, subset = expression_group %in% c("High"))
table(scp1219_PDLIM2_AT1_Pdlim2High[[]]$group)

# Fetch from above seurat object the Pdlim2 expression data from Covid and Ctrl groups
scp1219_PDLIM2_AT1_Pdlim2High_PDLIM2 <- FetchData(object = scp1219_PDLIM2_AT1_Pdlim2High, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (scp1219_PDLIM2_AT1_Pdlim2High_PDLIM2, file="scp1219_PDLIM2_AT1_Pdlim2High_PDLIM2.csv")

#check active.ident
head(scp1219_PDLIM2_AT1_Pdlim2High@active.ident)

#Changing active.ident in Seurat
scp1219_PDLIM2_AT1_Pdlim2High <- SetIdent(scp1219_PDLIM2_AT1_Pdlim2High, value = scp1219_PDLIM2_AT1_Pdlim2High@meta.data$group)

# DEG AT1_Pdlim2High_Covid-Ctrl
DEG <- FindMarkers(scp1219_PDLIM2_AT1_Pdlim2High,
                   ident.1 = "COVID-19",
                   ident.2 = "Control",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT1_Pdlim2High_Covid-Ctrl02.csv")

## sum expression of NF-kB targetd genes in AT1/AT2/Airway/Fibroblast
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_Fibroblasts, layer = "data")

# Sum the expression for your gene set per cell
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add gene-set-sum to scp1219_PDLIM2_AT1 metadata
scp1219_PDLIM2_Fibroblasts$NFkBSum <- gene_set_sum
colnames(scp1219_PDLIM2_Fibroblasts[[]])
head(scp1219_PDLIM2_Fibroblasts[[]]$NFkBSum)

# fetch data
scp1219_PDLIM2_Fibroblasts_NFkBSum <- FetchData(object = scp1219_PDLIM2_Fibroblasts, vars = c("group", "expression_group", "PDLIM2", "NFkBSum")) 
write.csv (scp1219_PDLIM2_Fibroblasts_NFkBSum, file="scp1219_PDLIM2_Fibroblasts_NFkBSum.csv")

#### average expression of NF-kB targetd genes in AT1/AT2/Airway/Fibroblast
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_Fibroblasts, layer = "data")

# average expression for your gene set per cell
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to scp1219_PDLIM2_AT1 metadata
scp1219_PDLIM2_Fibroblasts$NFkBSum <- gene_set_sum
scp1219_PDLIM2_Fibroblasts$NFkBAverage <- avg_expr_per_cell


colnames(scp1219_PDLIM2_Fibroblasts[[]])
head(scp1219_PDLIM2_Fibroblasts[[]]$NFkBSum)
head(scp1219_PDLIM2_Fibroblasts[[]]$NFkBAverage)

# fetch data
scp1219_PDLIM2_Fibroblasts_NFkB <- FetchData(object = scp1219_PDLIM2_Fibroblasts, vars = c("group", "expression_group", "PDLIM2", "NFkBSum","NFkBAverage")) 
write.csv (scp1219_PDLIM2_Fibroblasts_NFkB, file="scp1219_PDLIM2_Fibroblasts_NFkB.csv")

#divide disease cells into Pdlim2 low and High groups, and calculate their DEGs to Ctrl, and sum expression of STAT signaling target genes
# Subset disease cells
disease_cells <- WhichCells(scp1219_PDLIM2_AT2, expression = group == "COVID-19")

# Get PDLIM2 expression in disease cells
pdl_expr <- FetchData(scp1219_PDLIM2_AT2, vars = "PDLIM2")[disease_cells, , drop = FALSE]

# Determine median expression
pdl_median <- median(pdl_expr$PDLIM2)

# Classify disease cells into 'low' and 'high'
low_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 <= pdl_median]
high_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 > pdl_median]

# Assign new identities
scp1219_PDLIM2_AT2$DiseasePdlim2 <- "other"  # default
scp1219_PDLIM2_AT2$DiseasePdlim2[low_cells] <- "disease_low"
scp1219_PDLIM2_AT2$DiseasePdlim2[high_cells] <- "disease_high"
scp1219_PDLIM2_AT2$DiseasePdlim2[scp1219_PDLIM2_AT2$group == "Control"] <- "Control"

# fetch data
scp1219_PDLIM2_AT2_DiseasePdlim2 <- FetchData(object = scp1219_PDLIM2_AT2, vars = c("group", "expression_group", "PDLIM2", "DiseasePdlim2")) 
write.csv (scp1219_PDLIM2_AT2$DiseasePdlim2, file="scp1219_PDLIM2_AT2$DiseasePdlim2.csv")

# Set identities
Idents(scp1219_PDLIM2_AT2) <- scp1219_PDLIM2_AT2$DiseasePdlim2

# Run differential expression
DEG <- FindMarkers(scp1219_PDLIM2_AT2,
                   ident.1 = "disease_low",
                   ident.2 = "Control",
                   group.by = "DiseasePdlim2",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_AT2_Pdlim2Low-Ctrl.csv")

#sum expression of STAT target genes
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AT2, layer = "data")

# sum/average expression for your gene set per cell
gene_set <- c("CDKN1A", "CISH", "EGF", "ERAS", "HRAS", "IL6", "IL6R", "JAK2", "KRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK14", "MAPK3", "MRAS", "MYC", "NRAS", "PDGFB", "PIAS3", "PTPN6", "RAF1", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RRAS", "RRAS2", "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "STAT3", "TYK2", "VEGFA", "BCL2", "BMP6", "BMPR1A", "BMPR1B", "BMPR2", "CDC25A", "CSF2RB", "CXCR1", "CXCR2", "EGFR", "FGF2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT1", "FLT4", "GHR", "HGF", "IFNAR1", "IFNLR1", "IGF1", "IGF1R", "IGF2R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RL1", "IL1RL2", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2", "IL27RA", "IL2RA", "IL2RB", "IL2RG", "IL31RA", "IL3RA", "IL4R", "IL5RA", "IL6ST", "IL7R", "IL9", "IL9R", "INSR", "KDR", "MAP2K4", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K20", "MAP3K21", "MAP3K9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK8", "MAPK9", "NDUFA13", "NGFR", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PDGFRB", "PIM1", "PTPN2", "RAC1", "SRC", "TGFA", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TGFBR3", "TNFRSF11A", "IFNG", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IL12B", "TNF", "NOS2", "IRF1", "GBP1-5", "IDO1", "PSMB8", "TAP1", "TAP2", "HLA-DRA", "ICAM1", "FCGR1A", "BST2", "IL10", "LIF", "CSF2", "PDGFA", "PDGFC", "PDGFD", "CCL2", "SETDB2", "MMP9", "CCL17", "CCL22", "ARG1", "MRC1", "CHIT1", "IL4", "PPARG", "GATA3", "AKT1", "AKT2", "AKT3", "BCL2L1", "CCKBR", "CEBPB", "FOS", "GAST", "GNAQ", "GRB2", "JAK1", "JAK3", "JUN", "MTOR", "NFKB1", "NFKB2", "PIAS1", "PIAS2", "PIAS4", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PTPN1", "PTPN11", "REL", "RELA", "RELB", "SHC1", "SOS1", "SOS2", "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AT2$STATSum <- gene_set_sum
scp1219_PDLIM2_AT2$STATAverage <- avg_expr_per_cell


colnames(scp1219_PDLIM2_AT2[[]])
head(scp1219_PDLIM2_AT2[[]]$STATSum)
head(scp1219_PDLIM2_AT2[[]]$STATAverage)

# fetch data
scp1219_PDLIM2_AT2_STAT <- FetchData(object = scp1219_PDLIM2_AT2, vars = c("group", "expression_group", "PDLIM2", "STATSum","STATAverage", "DiseasePdlim2")) 
write.csv (scp1219_PDLIM2_AT2_STAT, file="scp1219_PDLIM2_AT2_STAT-PDLIM2.csv")

# sum expression of MHCI antigen presentation genes
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AT1, layer = "data")

#reactome MHCI antigen presentation molecules
gene_set <- c("ADRM1", "ANAPC1", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "AREL1", "ARIH2", "ASB1", "ASB10", "ASB11", "ASB12", "ASB13", "ASB14", "ASB15", "ASB16", "ASB17", "ASB18", "ASB2", "ASB3", "ASB4", "ASB5", "ASB6", "ASB7", "ASB8", "ASB9", "ATG14", "ATG7", "B2M", "BCAP31", "BECN1", "BLMH", "BTBD1", "BTBD6", "BTK", "BTRC", "CALR", "CANX", "CBLB", "CBLL2", "CCNF", "CD14", "CD207", "CD36", "CDC16", "CDC20", "CDC23", "CDC26", "CDC27", "CDC34", "CHUK", "CTSL", "CTSS", "CTSV", "CUL1", "CUL2", "CUL3", "CUL5", "CUL7", "CYBA", "CYBB", "DCAF1", "DET1", "DTX3L", "DZIP3", "ELOB", "ELOC", "ERAP1", "ERAP2", "FBXL12", "FBXL13", "FBXL14", "FBXL15", "FBXL16", "FBXL18", "FBXL19", "FBXL20", "FBXL22", "FBXL3", "FBXL4", "FBXL5", "FBXL7", "FBXL8", "FBXO10", "FBXO11", "FBXO15", "FBXO17", "FBXO2", "FBXO21", "FBXO22", "FBXO27", "FBXO30", "FBXO31", "FBXO32", "FBXO4", "FBXO40", "FBXO41", "FBXO44", "FBXO6", "FBXO7", "FBXO9", "FBXW10", "FBXW11", "FBXW12", "FBXW2", "FBXW4", "FBXW5", "FBXW7", "FBXW8", "FBXW9", "FCGR1A", "FGA", "FGB", "FGG", "FZR1", "GAN", "GLMN", "HACE1", "HECTD1", "HECTD2", "HECTD3", "HECW2", "HERC1", "HERC2", "HERC3", "HERC4", "HERC5", "HERC6", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HMGB1", "HSPA5", "HUWE1", "IKBKB", "IKBKG", "ITCH", "ITGAV", "ITGB5", "KBTBD13", "KBTBD6", "KBTBD7", "KBTBD8", "KCTD6", "KCTD7", "KEAP1", "KLHL11", "KLHL13", "KLHL2", "KLHL20", "KLHL21", "KLHL22", "KLHL25", "KLHL3", "KLHL41", "KLHL42", "KLHL5", "KLHL9", "LMO7", "LNPEP", "LNX1", "LONRF1", "LRR1", "LRRC41", "LRSAM1", "LTN1", "LY96", "MEX3C", "MGRN1", "MIB2", "MKRN1", "MRC1", "MRC2", "MYD88", "MYLIP", "NCF1", "NCF2", "NCF4", "NEDD4", "NEDD4L", "NPEPPS", "PDIA3", "PIK3C3", "PIK3R4", "PJA1", "PJA2", "PRKN", "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB10", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD11", "PSMD12", "PSMD13", "PSMD14", "PSMD2", "PSMD3", "PSMD6", "PSMD7", "PSMD8", "PSME1", "PSME2", "RBBP6", "RBCK1", "RBX1", "RCHY1", "RLIM", "RNF111", "RNF114", "RNF115", "RNF123", "RNF126", "RNF130", "RNF138", "RNF14", "RNF144B", "RNF182", "RNF19A", "RNF19B", "RNF213", "RNF217", "RNF220", "RNF25", "RNF34", "RNF4", "RNF41", "RNF6", "RNF7", "RPS27A", "S100A1", "S100A8", "S100A9", "SAR1B", "SEC13", "SEC22B", "SEC23A", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC31A", "SEC61A1", "SEC61A2", "SEC61B", "SEC61G", "SEM1", "SH3RF1", "SIAH1", "SIAH2", "SKP1", "SKP2", "SMURF1", "SMURF2", "SNAP23", "SOCS1", "SOCS3", "SPSB1", "SPSB2", "SPSB4", "STUB1", "STX4", "TAP1", "TAP2", "TAPBP", "THOP1", "TIRAP", "TLR1", "TLR2", "TLR4", "TLR6", "TPP2", "TRAF7", "TRAIP", "TRIM11", "TRIM21", "TRIM32", "TRIM36", "TRIM37", "TRIM39", "TRIM4", "TRIM41", "TRIM50", "TRIM63", "TRIM69", "TRIM71", "TRIM9", "TRIP12", "UBA1", "UBA3", "UBA5", "UBA52", "UBA6", "UBA7", "UBAC1", "UBB", "UBC", "UBE2A", "UBE2B", "UBE2C", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2D4", "UBE2E1", "UBE2E2", "UBE2E3", "UBE2F", "UBE2G1", "UBE2G2", "UBE2H", "UBE2J1", "UBE2J2", "UBE2K", "UBE2L3", "UBE2L6", "UBE2M", "UBE2N", "UBE2O", "UBE2Q1", "UBE2Q2", "UBE2R2", "UBE2S", "UBE2U", "UBE2V1", "UBE2V2", "UBE2W", "UBE2Z", "UBE3A", "UBE3B", "UBE3C", "UBE3D", "UBE4A", "UBOX5", "UBR1", "UBR2", "UBR4", "UFL1", "UNKL", "VAMP3", "VAMP8", "VHL", "WSB1", "WWP1", "ZBTB16", "ZNRF1", "ZNRF2")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AT1$MHCI <- gene_set_sum

colnames(scp1219_PDLIM2_AT1[[]])
head(scp1219_PDLIM2_AT1[[]]$MHCI)

# fetch data
scp1219_PDLIM2_AT1_MHCI <- FetchData(object = scp1219_PDLIM2_AT1, vars = c("group", "expression_group", "PDLIM2", "MHCI")) 
write.csv (scp1219_PDLIM2_AT1_MHCI, file="scp1219_PDLIM2_AT1_MHCI-PDLIM2.csv")

# sum expression of Qiagen HIF1a signaling genelist in AT1
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AT1, layer = "data")

# Qiagen HIF1a signaling genelist
gene_set<- c("ADM", "ADRA1B", "AKT1", "AKT2", "AKT3", "APEX1", "ARAF", "ARNT", "BMP6", "BRAF", "CAMK1", "CAMK1D", "CAMK1G", "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G", "CAMK4", "CCNG2", "CDKN1A", "CHP1", "COPS5", "CREBBP", "CUL2", "CYBB", "EDN1", "EGF", "EGLN1", "EGLN2", "EGLN3", "EIF4E", "EIF4EBP1", "ELOB", "ELOC", "EP300", "EPO", "ERAS", "FGF2", "FLT1", "FLT4", "FOXP3", "GCK", "GPI", "HGF", "HIF1A", "HIF1AN", "HK1", "HK2", "HK3", "HMOX1", "HRAS", "HSP90AA1", "HSPA14", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", "HSPA5", "HSPA6", "HSPA8", "HSPA9", "HSPH1", "IGF1", "IGF2", "IL17A", "IL6", "IL6R", "JUN", "KDM1A", "KDR", "KRAS", "LDHA", "LDHB", "LDHC", "MAP2K1", "MAP2K2", "MAP2K3", "MAP2K4", "MAP2K5", "MAP2K6", "MAP2K7", "MAPK1", "MAPK14", "MAPK3", "MDM2", "MET", "MIR135A1", "MIR3662", "MKNK1", "MKNK2", "MMP1", "MMP10", "MMP11", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16", "MMP17", "MMP19", "MMP2", "MMP20", "MMP21", "MMP23B", "MMP24", "MMP25", "MMP26", "MMP27", "MMP28", "MMP3", "MMP7", "MMP8", "MMP9", "MRAS", "MTOR", "NAA10", "NCF1", "NCF2", "NCOA1", "NOS2", "NOX1", "NOX3", "NOX4", "NRAS", "P4HTM", "PDGFB", "PDGFC", "PGF", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PKM", "PLCG1", "PLCG2", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1", "PPP3R2", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCI", "PRKCQ", "PRKCZ", "PRKD1", "PRKD3", "PROK1", "RAC1", "RAC2", "RAC3", "RACK1", "RAF1", "RALA", "RALB", "RAN", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RBX1", "RORC", "RPS6", "RPS6KB1", "RPS6KB2", "RRAS", "RRAS2", "SAT1", "SAT2", "SERPINE1", "SLC2A1", "SLC2A14", "SLC2A2", "SLC2A3", "SLC2A4", "SLC2A5", "SLC2A8", "STAT3", "STUB1", "TF", "TGFA", "TGFB1", "TGFB2", "TP53", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VHL", "VIM")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AT1$HIF1a <- gene_set_sum

colnames(scp1219_PDLIM2_AT1[[]])
head(scp1219_PDLIM2_AT1[[]]$HIF1a)

# fetch data
scp1219_PDLIM2_AT1_HIF1a <- FetchData(object = scp1219_PDLIM2_AT1, vars = c("group", "expression_group", "PDLIM2", "HIF1a")) 
write.csv (scp1219_PDLIM2_AT1_HIF1a, file="scp1219_PDLIM2_AT1_HIF1a-PDLIM2.csv")

# sum expression of Qiagen Virus entry via endocytic genelist in AT2
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AT2, layer = "data")

# Qiagen Virus entry via endocytic genelist
gene_set <- c("ABL1", "ACTA1", "ACTA2", "ACTB", "ACTC1", "ACTG1", "ACTG2", "AP1B1", "AP1G1", "AP1G2", "AP1M1", "AP1M2", "AP1S1", "AP1S2", "AP1S3", "AP2A1", "AP2A2", "AP2B1", "AP2M1", "AP2S1", "AP3B1", "AP3B2", "AP3D1", "AP3M1", "AP3M2", "AP3S1", "AP3S2", "AP4B1", "AP4E1", "AP4M1", "AP4S1", "AP5B1", "AP5M1", "AP5S1", "AP5Z1", "CAV1", "CD55", "CDC42", "CLTA", "CLTB", "CLTC", "CLTCL1", "CXADR", "DNM1", "DNM2", "ERAS", "FLNA", "FLNB", "FLNC", "FOLR1", "FYN", "HRAS", "ITGA2", "ITGA5", "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8", "ITSN1", "KRAS", "MRAS", "NRAS", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PLCG1", "PLCG2", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCI", "PRKCQ", "PRKCZ", "PRKD1", "PRKD3", "RAC1", "RAC2", "RAC3", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RRAS", "RRAS2", "SRC", "TFRC")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AT2$Endocytic <- gene_set_sum

colnames(scp1219_PDLIM2_AT2[[]])
head(scp1219_PDLIM2_AT2[[]]$Endocytic)

# fetch data
scp1219_PDLIM2_AT2_Endocytic <- FetchData(object = scp1219_PDLIM2_AT2, vars = c("group", "expression_group", "PDLIM2", "Endocytic")) 
write.csv (scp1219_PDLIM2_AT2_Endocytic, file="scp1219_PDLIM2_AT2_Endocytic-PDLIM2.csv")

#sum expression of Immunogenic cell death genes in AT2
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AT2, layer = "data")

# Immunogenic cell death genes
gene_set <- c("ARHGAP10", "TLR4", "TNFRSF10A", "PSMD11", "PSMD12", "DFFA", "DNM1L", "PSMD14", "SDCBP", "APAF1", "TNFRSF10B", "PSMA7", "OGT", "TP73", "CFLAR", "SEPTIN4", "PSMD3", "DAPK3", "BCL2L11", "CHMP2A", "GAS2", "PRKN", "OPA1", "GSDME", "FLOT1", "DFFB", "IL1A", "IL1B", "LMNA", "", "TP53", "GSN", "H1-0", "HSP90AA1", "ELANE", "", "CD14", "VIM", "HMGB1", "UBB", "UBC", "GZMB", "H1-4", "BCL2", "MAPT", "IRF1", "CDH1", "IRF2", "DSP", "H1-5", "H1-3", "H1-2", "PSMC3", "PSMB1", "LMNB1", "APC", "FAS", "PSMA1", "PSMA2", "PSMA3", "PSMA4", "HMGB2", "YWHAQ", "MAPK3", "PSMA5", "PSMB4", "PSMB6", "PSMB5", "MAPK1", "CASP1", "NMT1", "AKT1", "AKT2", "YWHAB", "SFN", "DSG3", "CTNNB1", "ADD1", "PSMC2", "STAT3", "CASP3", "DCC", "PSMC4", "MAPK8", "FASLG", "PPP3CC", "PSMD8", "FNTA", "CASP4", "PSMB3", "PSMB2", "TNFSF10", "BCAP31", "PSMD7", "BMX", "CASP5", "KPNA1", "DAPK1", "CASP7", "CASP9", "CASP6", "BID", "GSDMD", "SEM1", "PSMA6", "YWHAG", "PSMC1", "PSMC5", "YWHAE", "PSMC6", "RPS27A", "UBA52", "PPP3R1", "YWHAZ", "DYNLL1", "UBE2L3", "XIAP", "CYCS", "E2F1", "SATB1", "DSG1", "H1-1", "PRKCQ", "YWHAH", "PTK2", "PRKCD", "C1QBP", "TJP1", "BAX", "BCL2L1", "TRAF2", "FADD", "PAK2", "PSMD2", "ROCK1", "BIRC3", "BIRC2", "RIPK1", "TP53BP2", "PMAIP1", "SPTAN1", "PKP1", "IL18", "DSG2", "TFDP1", "TFDP2", "FLOT2", "CASP8", "KPNB1", "PSMD6", "PLEC", "TRADD", "ADRM1", "CDC37", "BAK1", "OCLN", "UNC5A", "", "TMED7-TICAM2", "TICAM1", "UNC5B", "CDKN2A", "MLKL", "PDCD6IP", "CHMP7", "BAD", "CHMP4C", "OMA1", "PELI1", "DYNLL2", "CHMP6", "APIP", "ITCH", "PPP1R13B", "BMF", "", "PSMB7", "PSMD1", "BBC3", "CHMP4A", "UACA", "TP63", "CHMP4B", "CLSPN", "AVEN", "DIABLO", "STK26", "TJP2", "DAPK2", "DBNL", "APPL1", "ACIN1", "STUB1", "PSMD13", "CHMP2B", "AKT3", "CARD8", "CHMP3", "RIPK3", "MAGED1", "STK24", "LY96")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AT2$Immunogenic <- gene_set_sum

colnames(scp1219_PDLIM2_AT2[[]])
head(scp1219_PDLIM2_AT2[[]]$Immunogenic)

# fetch data
scp1219_PDLIM2_AT2_Immunogenic <- FetchData(object = scp1219_PDLIM2_AT2, vars = c("group", "expression_group", "PDLIM2", "Immunogenic")) 
write.csv (scp1219_PDLIM2_AT2_Immunogenic, file="scp1219_PDLIM2_AT2_Immunogenic-PDLIM2.csv")

# sum expression of MHCI antigen presentation genes in Airway
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AirwayEpi, layer = "data")

#reactome MHCI antigen presentation molecules
gene_set <- c("ADRM1", "ANAPC1", "ANAPC10", "ANAPC11", "ANAPC13", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "AREL1", "ARIH2", "ASB1", "ASB10", "ASB11", "ASB12", "ASB13", "ASB14", "ASB15", "ASB16", "ASB17", "ASB18", "ASB2", "ASB3", "ASB4", "ASB5", "ASB6", "ASB7", "ASB8", "ASB9", "ATG14", "ATG7", "B2M", "BCAP31", "BECN1", "BLMH", "BTBD1", "BTBD6", "BTK", "BTRC", "CALR", "CANX", "CBLB", "CBLL2", "CCNF", "CD14", "CD207", "CD36", "CDC16", "CDC20", "CDC23", "CDC26", "CDC27", "CDC34", "CHUK", "CTSL", "CTSS", "CTSV", "CUL1", "CUL2", "CUL3", "CUL5", "CUL7", "CYBA", "CYBB", "DCAF1", "DET1", "DTX3L", "DZIP3", "ELOB", "ELOC", "ERAP1", "ERAP2", "FBXL12", "FBXL13", "FBXL14", "FBXL15", "FBXL16", "FBXL18", "FBXL19", "FBXL20", "FBXL22", "FBXL3", "FBXL4", "FBXL5", "FBXL7", "FBXL8", "FBXO10", "FBXO11", "FBXO15", "FBXO17", "FBXO2", "FBXO21", "FBXO22", "FBXO27", "FBXO30", "FBXO31", "FBXO32", "FBXO4", "FBXO40", "FBXO41", "FBXO44", "FBXO6", "FBXO7", "FBXO9", "FBXW10", "FBXW11", "FBXW12", "FBXW2", "FBXW4", "FBXW5", "FBXW7", "FBXW8", "FBXW9", "FCGR1A", "FGA", "FGB", "FGG", "FZR1", "GAN", "GLMN", "HACE1", "HECTD1", "HECTD2", "HECTD3", "HECW2", "HERC1", "HERC2", "HERC3", "HERC4", "HERC5", "HERC6", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "HMGB1", "HSPA5", "HUWE1", "IKBKB", "IKBKG", "ITCH", "ITGAV", "ITGB5", "KBTBD13", "KBTBD6", "KBTBD7", "KBTBD8", "KCTD6", "KCTD7", "KEAP1", "KLHL11", "KLHL13", "KLHL2", "KLHL20", "KLHL21", "KLHL22", "KLHL25", "KLHL3", "KLHL41", "KLHL42", "KLHL5", "KLHL9", "LMO7", "LNPEP", "LNX1", "LONRF1", "LRR1", "LRRC41", "LRSAM1", "LTN1", "LY96", "MEX3C", "MGRN1", "MIB2", "MKRN1", "MRC1", "MRC2", "MYD88", "MYLIP", "NCF1", "NCF2", "NCF4", "NEDD4", "NEDD4L", "NPEPPS", "PDIA3", "PIK3C3", "PIK3R4", "PJA1", "PJA2", "PRKN", "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMB1", "PSMB10", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6", "PSMD1", "PSMD11", "PSMD12", "PSMD13", "PSMD14", "PSMD2", "PSMD3", "PSMD6", "PSMD7", "PSMD8", "PSME1", "PSME2", "RBBP6", "RBCK1", "RBX1", "RCHY1", "RLIM", "RNF111", "RNF114", "RNF115", "RNF123", "RNF126", "RNF130", "RNF138", "RNF14", "RNF144B", "RNF182", "RNF19A", "RNF19B", "RNF213", "RNF217", "RNF220", "RNF25", "RNF34", "RNF4", "RNF41", "RNF6", "RNF7", "RPS27A", "S100A1", "S100A8", "S100A9", "SAR1B", "SEC13", "SEC22B", "SEC23A", "SEC24A", "SEC24B", "SEC24C", "SEC24D", "SEC31A", "SEC61A1", "SEC61A2", "SEC61B", "SEC61G", "SEM1", "SH3RF1", "SIAH1", "SIAH2", "SKP1", "SKP2", "SMURF1", "SMURF2", "SNAP23", "SOCS1", "SOCS3", "SPSB1", "SPSB2", "SPSB4", "STUB1", "STX4", "TAP1", "TAP2", "TAPBP", "THOP1", "TIRAP", "TLR1", "TLR2", "TLR4", "TLR6", "TPP2", "TRAF7", "TRAIP", "TRIM11", "TRIM21", "TRIM32", "TRIM36", "TRIM37", "TRIM39", "TRIM4", "TRIM41", "TRIM50", "TRIM63", "TRIM69", "TRIM71", "TRIM9", "TRIP12", "UBA1", "UBA3", "UBA5", "UBA52", "UBA6", "UBA7", "UBAC1", "UBB", "UBC", "UBE2A", "UBE2B", "UBE2C", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2D4", "UBE2E1", "UBE2E2", "UBE2E3", "UBE2F", "UBE2G1", "UBE2G2", "UBE2H", "UBE2J1", "UBE2J2", "UBE2K", "UBE2L3", "UBE2L6", "UBE2M", "UBE2N", "UBE2O", "UBE2Q1", "UBE2Q2", "UBE2R2", "UBE2S", "UBE2U", "UBE2V1", "UBE2V2", "UBE2W", "UBE2Z", "UBE3A", "UBE3B", "UBE3C", "UBE3D", "UBE4A", "UBOX5", "UBR1", "UBR2", "UBR4", "UFL1", "UNKL", "VAMP3", "VAMP8", "VHL", "WSB1", "WWP1", "ZBTB16", "ZNRF1", "ZNRF2")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AirwayEpi$MHCI <- gene_set_sum

colnames(scp1219_PDLIM2_AirwayEpi[[]])
head(scp1219_PDLIM2_AirwayEpi[[]]$MHCI)

# fetch data
scp1219_PDLIM2_AirwayEpi_MHCI <- FetchData(object = scp1219_PDLIM2_AirwayEpi, vars = c("group", "expression_group", "PDLIM2", "MHCI")) 
write.csv (scp1219_PDLIM2_AirwayEpi_MHCI, file="scp1219_PDLIM2_AirwayEpi_MHCI-PDLIM2.csv")

# sum expression of IL-1 signaling genelist
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_AirwayEpi, layer = "data")

# Qiagen IL-1 signaling genelist
gene_set <- c("ADCY1", "ADCY10", "ADCY2", "ADCY3", "ADCY4", "ADCY5", "ADCY6", "ADCY7", "ADCY8", "ADCY9", "CAMP", "CHUK", "ECSIT", "ELP1", "FOS", "GNA11", "GNA12", "GNA13", "GNA14", "GNA15", "GNAI1", "GNAI2", "GNAI3", "GNAL", "GNAO1", "GNAQ", "GNAS", "GNAT1", "GNAT2", "GNAZ", "GNB1", "GNB1L", "GNB2", "GNB3", "GNB4", "GNB5", "GNG10", "GNG11", "GNG12", "GNG13", "GNG2", "GNG3", "GNG4", "GNG5", "GNG7", "GUCY1A1", "GUCY1B1", "IKBKB", "IKBKE", "IKBKG", "IL1A", "IL1R1", "IL1RAP", "IRAK1", "IRAK2", "IRAK3", "IRAK4", "JUN", "MAP2K3", "MAP2K4", "MAP2K6", "MAP2K7", "MAP3K1", "MAP3K14", "MAP3K7", "MAPK1", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "MAPK8", "MAPK9", "MRAS", "MYD88", "NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "NFKBID", "NFKBIE", "PRKACA", "PRKACB", "PRKACG", "PRKAG1", "PRKAG2", "PRKAR1A", "PRKAR1B", "PRKAR2A", "PRKAR2B", "REL", "RELA", "RELB", "TAB1", "TAB2", "TOLLIP", "TRAF6")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_AirwayEpi$IL1 <- gene_set_sum

colnames(scp1219_PDLIM2_AirwayEpi[[]])
head(scp1219_PDLIM2_AirwayEpi[[]]$IL1)

# fetch data
scp1219_PDLIM2_AirwayEpi_IL1 <- FetchData(object = scp1219_PDLIM2_AirwayEpi, vars = c("group", "expression_group", "PDLIM2", "IL1")) 
write.csv (scp1219_PDLIM2_AirwayEpi_IL1, file="scp1219_PDLIM2_AirwayEpi_IL1-PDLIM2.csv")

# sum expression of Interferon a/b signaling genes
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_Fibroblasts, layer = "data")

# reactome IFNa&b signaling and antiviral mechanism by ISGs molecules
gene_set <- c("ABCE1", "ADAR", "BST2", "EGR1", "EIF2AK2", "GBP2", "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "IFI27", "IFI35", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3", "IFNA1", "IFNA10", "IFNA14", "IFNA16", "IFNA17", "IFNA2", "IFNA21", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", "IFNAR1", "IFNAR2", "IFNB1", "IP6K2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9", "ISG15", "ISG20", "JAK1", "KPNA1", "KPNB1", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "PSMB8", "PTPN1", "PTPN11", "PTPN6", "RNASEL", "RPS27A", "RSAD2", "SAMHD1", "SOCS1", "SOCS3", "STAT1", "STAT2", "TYK2", "UBA52", "UBB", "UBC", "USP18", "XAF1", "AAAS", "ABCE1", "ADAR", "ARIH1", "BECN1", "CDK1", "CENPS-CORT", "CENPX", "CHUK", "DHX9", "DNAJC3", "DUS2", "EIF2AK2", "EIF2S1", "EIF2S2", "EIF2S3", "EIF4A1", "EIF4A2", "EIF4A3", "EIF4E", "EIF4E2", "EIF4E3", "EIF4G1", "EIF4G2", "EIF4G3", "FAAP100", "FAAP20", "FAAP24", "FANCA", "FANCB", "FANCC", "FANCE", "FANCF", "FANCG", "FANCL", "FANCM", "FLNA", "FLNB", "HERC5", "HSPA1A", "HSPA1L", "HSPA2", "HSPA8", "IFIT1", "IKBKB", "IKBKG", "ILF2", "ILF3", "IRF3", "ISG15", "JAK1", "KPNA1", "KPNA2", "KPNA3", "KPNA4", "KPNA5", "KPNA7", "KPNB1", "MAP2K6", "MAPK3", "MAPT", "MAVS", "MX1", "MX2", "NCK1", "NDC1", "NEDD4", "NPM1", "NUP107", "NUP133", "NUP153", "NUP155", "NUP160", "NUP188", "NUP205", "NUP210", "NUP214", "NUP35", "NUP37", "NUP42", "NUP43", "NUP50", "NUP54", "NUP58", "NUP62", "NUP85", "NUP88", "NUP93", "NUP98", "OAS1", "OAS2", "OAS3", "OASL", "PDE12", "PIN1", "PLCG1", "POM121", "PPM1B", "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B", "PPP2R5A", "PRKRA", "PTPN2", "RAE1", "RANBP2", "RIGI", "RNASEL", "RPS27A", "SEC13", "SEH1L", "SNCA", "SPHK1", "STAT1", "STAT3", "SUMO1", "TARBP2", "TP53", "TPR", "TRIM25", "TUBA1A", "TUBA1B", "TUBA1C", "TUBA3C", "TUBA3D", "TUBA3E", "TUBA4A", "TUBA4B", "TUBA8", "TUBAL3", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3", "TUBB4A", "TUBB4B", "TUBB6", "TUBB8", "TUBB8B", "UBA52", "UBA7", "UBB", "UBC", "UBE2E1", "UBE2I", "UBE2L6", "UBE2N", "USP18")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_Fibroblasts$IFN <- gene_set_sum

colnames(scp1219_PDLIM2_Fibroblasts[[]])
head(scp1219_PDLIM2_Fibroblasts[[]]$IFN)

# fetch data
scp1219_PDLIM2_Fibroblasts_IFN <- FetchData(object = scp1219_PDLIM2_Fibroblasts, vars = c("group", "expression_group", "PDLIM2", "IFN")) 
write.csv (scp1219_PDLIM2_Fibroblasts_IFN, file="scp1219_PDLIM2_Fibroblasts_IFN-PDLIM2.csv")

# sum expression of pulmonary Fibrosis genelist
# Extract the normalized expression values
norm_expr <- GetAssayData(scp1219_PDLIM2_Fibroblasts, layer = "data")

# Qiagen pathways&Lists Pulmonary fibrosis genelist
gene_set <- c("ACTA2","ADAM10","JAG1","BIRC5","AREG","RHOA","AXL","BAX","CCND1","BCL2","CAV1","CDH1","CDH2","COL1A1","COL3A1","ATF2","CREBBP","CCN2","CTNNB1","LPAR1","EDN1","EDNRA","EFNB2soluble","EGFR","EGR1","EPHB3","EPHB4","EZH1","EZH2","F2","F2R","FGF2","FGF9","FN1","FOS","MTOR","GLI1","GLI2","GSK3B","HES1","RBPJ","IL1B","IL4","IL6","IL6ST","IL11","IL13","IL13RA1","IL17A","ILK","ITGB1","ITGB6","JAK2","JUN","SMAD2","SMAD3","SMAD4","MMP2","MMP9","MUC1","NOTCHIC","SERPINE1","PRKN","PEX13","PLAU","PLG","PMAIP1","POLR2A","MAP2K1","PTEN","RELA","CXCL12","MAP2K4","SFTPC","SNAI1","SRF","STAT3","STAT6","MAP3K7","TERC","TERT","TGFA","THBS1","TP53","TYK2","VIM","WNT1","FGF18","CCN4","ADAMTS1","IL17RA","BBC3","NOX4","RTEL1","MRTFA","PINK1","AJUBA","WNT3A","SFTPA2","MIR410")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
scp1219_PDLIM2_Fibroblasts$Fibrosis <- gene_set_sum

colnames(scp1219_PDLIM2_Fibroblasts[[]])
head(scp1219_PDLIM2_Fibroblasts[[]]$Fibrosis)

# fetch data
scp1219_PDLIM2_Fibroblasts_Fibrosis <- FetchData(object = scp1219_PDLIM2_Fibroblasts, vars = c("group", "expression_group", "PDLIM2", "Fibrosis")) 
write.csv (scp1219_PDLIM2_Fibroblasts_Fibrosis, file="scp1219_PDLIM2_Fibroblasts_Fibrosis-PDLIM2.csv")

load("GSE171524.RData")
save (list=ls(), file="GSE171524.RData")
 