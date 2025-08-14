#Python error when openning RStudioo on Macbook
#check the running version of RStudio
#on Thinkpad: [1] ‘2023.6.1.524’ _ subset() works
#on Macbook: [1] ‘2024.12.0.467’  _ subset() not working
RStudio.Version()$version

#problem solved with steps below
#First, verify that Python is installed correctly. Open Terminal and run:
/usr/bin/python3 --version
# error found and Python installation occured automatically
# open RStudio and no more Python error
# however subset() still not working with error: Error in slot(object = object, name = s) : 
no slot of name "images" for this object of class "Seurat"


library(dplyr)
library(Seurat)
library(patchwork)

load ("/Users/fremontgao/Library/CloudStorage/OneDrive-UniversityofSouthernCalifornia/20240610/Projects/Pdlim2 in lung infection/scRNAseq/GSE145926_RAW/GSE145926.RData")

load("C:/Users/Fremont Gao/OneDrive - University of Southern California/20240610/Projects/Pdlim2 in lung infection/scRNAseq/GSE145926_RAW/GSE145926.RData")

load("GSE145926.RData")

# seurat subset with PDLIM2>0
nCoV_PDLIM2 <- subset(nCoV, subset = PDLIM2 >0)
nCoV_PDLIM2@assays$RNA@data["PDLIM2", 1:10]

unique(nCoV_PDLIM2[[]]$cluster)

unique(nCoV_PDLIM2[[]]$hasnCoV)

unique(nCoV_PDLIM2[[]]$group)

# DEG Total Covid vs Healthy
DEG <- FindMarkers(nCoV_PDLIM2,
                   ident.1 = "Y",
                   ident.2 = "N",
                   group.by = "hasnCoV",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Total_wilcox.csv")

# DEG Total HC vs O vs S/C
# to use FindAllMarkers, active.ident has to be preset in advance
head(nCoV_PDLIM2@active.ident)
nCoV_PDLIM2 <- SetIdent(nCoV_PDLIM2, value = nCoV_PDLIM2@meta.data$group)
head(nCoV_PDLIM2@active.ident)
DEG <- FindAllMarkers(nCoV_PDLIM2, group.by = "group")
write.csv (DEG, file = "DEG_Total_ANOVA.csv")


# subset Macrophage clusters for DEG
nCoV_PDLIM2_Macrophage <- subset(nCoV_PDLIM2, subset = cluster %in% c("0", "1", "2", "3", "4", "5", "7", "8", "10", "11", "12", "18", "21", "22", "23", "26"))

# DEG Macrophage Covid vs Healthy
DEG <- FindMarkers(nCoV_PDLIM2_Macrophage,
                   ident.1 = "N",
                   ident.2 = "Y",
                   group.by = "disease",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Macrophage_Ctrl-Covid_wilcox.csv")

# DEG Macrophage CovidMod/Svre vs Healthy
DEG <- FindMarkers(nCoV_PDLIM2_Macrophage,
                   ident.1 = "O",
                   ident.2 = "HC",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Macrophage_CovidMod-Ctrl_wilcox.csv")

# check to see if DEG (Covid-Ctrl) is equal to DEG((O+S/C)-HC), confirmed to be true
DEG <- FindMarkers(nCoV_PDLIM2_Macrophage,
                   ident.1 = c("O", "S/C"),
                   ident.2 = "HC",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Macrophage_CovidModSvre-Ctrl_wilcox.csv")

# DEG Macrophage HC vs O vs S/C
DEG <- FindAllMarkers(nCoV_PDLIM2_Macrophage, group.by = "group")

write.csv (DEG, file = "DEG_Macrophage_ANOVA.csv")

# subset T cell clusters for DEG
nCoV_PDLIM2_T <- subset(nCoV_PDLIM2, subset = cluster %in% c("6", "9", "14"))

# DEG T cell between Covid/Ctrl
DEG <- FindMarkers(nCoV_PDLIM2_T,
                   ident.1 = "N",
                   ident.2 = "Y",
                   group.by = "disease",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_T_Ctrl-Covid.csv")

# DEG T cell between CovidMod/Ctrl, CovidSvre/Ctrl
DEG <- FindMarkers(nCoV_PDLIM2_T,
                   ident.1 = "S/C",
                   ident.2 = "HC",
                   group.by = "group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_T_CovidSvre-Ctrl.csv")

# DEG T cell between PDLIM2H/LM, LM/H, M/H, L/H
DEG <- FindMarkers(nCoV_PDLIM2_T,
                   ident.1 = "High",
                   ident.2 = c("Low","Middle"),
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_T_Pdlim2H-LM.csv")


# DEG T with HC vs O vs S/C
head(nCoV_PDLIM2_T@active.ident)
nCoV_PDLIM2_T <- SetIdent(nCoV_PDLIM2_T, value = nCoV_PDLIM2_T@meta.data$group)
head(nCoV_PDLIM2_T@active.ident)
DEG <- FindAllMarkers(nCoV_PDLIM2_T, group.by = "group")

write.csv (DEG, file = "DEG_T_ANOVA.csv")

# subset PDLIM2+ Macrophage cells to High Middle Low level Pdlim2 groups
# Fetch PDLIM2 expression data and calculate its level at certain percentile
# Choose the gene of interest
gene_of_interest <- "PDLIM2"  # Replace with your gene

# Extract expression levels
gene_expr <- FetchData(nCoV_PDLIM2_Macrophage, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.67)
low_threshold <- quantile(gene_expr[,1], 0.33)

# Assign group labels
nCoV_PDLIM2_Macrophage$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High",
  ifelse(gene_expr[,1] <= low_threshold, "Low", "Middle")
)

# Check new metadata
table(nCoV_PDLIM2_Macrophage$expression_group)

# Fetch above High Middle Low PDLIM2 macrophage groups
nCoV_PDLIM2_Macrophage_PDLIM2Percentile <- FetchData(object = nCoV_PDLIM2_Macrophage, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (nCoV_PDLIM2_Macrophage_PDLIM2Percentile, file="nCoV_PDLIM2_Macrophage_PDLIM2Percentile.csv")

# DEG Macrophage with High vs Middle vs Low on PDLIM2 level
head(nCoV_PDLIM2_Macrophage@active.ident)
nCoV_PDLIM2_Macrophage <- SetIdent(nCoV_PDLIM2_Macrophage, value = nCoV_PDLIM2_Macrophage@meta.data$expression_group)
head(nCoV_PDLIM2_Macrophage@active.ident)
DEG <- FindAllMarkers(nCoV_PDLIM2_Macrophage, group.by = "expression_group")

write.csv (DEG, file = "DEG_Macrophage_PDLIM2_H-M-L_ANOVA.csv")

# DEG Macrophage with High/(Middle+Low), (middle+low)/high, middle/high, low/high on PDLIM2 level

DEG <- FindMarkers(nCoV_PDLIM2_Macrophage,
                   ident.1 = "High",
                   ident.2 = c("Low","Middle"),
                   group.by = "expression_group",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Macrophage_Pdlim2High-LowMid_wilcox.csv")





# subset PDLIM2+ T cells to High Middle Low level Pdlim2 groups
# Fetch PDLIM2 expression data and calculate its level at certain percentile
# Choose the gene of interest
gene_of_interest <- "PDLIM2"  # Replace with your gene

# Extract expression levels
gene_expr <- FetchData(nCoV_PDLIM2_T, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.67)
low_threshold <- quantile(gene_expr[,1], 0.33)

# Assign group labels
nCoV_PDLIM2_T$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High",
  ifelse(gene_expr[,1] <= low_threshold, "Low", "Middle")
)

# Check new metadata
table(nCoV_PDLIM2_T$expression_group)

# Fetch above High Middle Low PDLIM2 T groups
nCoV_PDLIM2_T_PDLIM2Percentile <- FetchData(object = nCoV_PDLIM2_T, vars = c("group", "expression_group", "PDLIM2")) 
write.csv (nCoV_PDLIM2_T_PDLIM2Percentile, file="nCoV_PDLIM2_T_PDLIM2Percentile.csv")

# DEG T with High vs Middle vs Low on PDLIM2 level
head(nCoV_PDLIM2_T@active.ident)
nCoV_PDLIM2_T <- SetIdent(nCoV_PDLIM2_T, value = nCoV_PDLIM2_T@meta.data$expression_group)
head(nCoV_PDLIM2_T@active.ident)
DEG <- FindAllMarkers(nCoV_PDLIM2_T, group.by = "expression_group")

write.csv (DEG, file = "DEG_T_PDLIM2_H-M-L_ANOVA.csv")

save (list=ls(), file="GSE145926.RData")

## check correlation b/w PDLIM2 and nCoV_mean
nCoV_PDLIM2_Macrophage_PDLIM2Percentile_nCoV <- FetchData(object = nCoV_PDLIM2_Macrophage, vars = c("group", "nCoV_mean", "expression_group", "PDLIM2")) 
write.csv (nCoV_PDLIM2_Macrophage_PDLIM2Percentile_nCoV, file="nCoV_PDLIM2_Macrophage_PDLIM2Percentile_nCoV.csv")

###### sum/average expression of all NF-kB targeted genes in macrophage/T
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T, layer = "data")

# average expression for your gene set per cell
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T$NFkBSum <- gene_set_sum
nCoV_PDLIM2_T$NFkBAverage <- avg_expr_per_cell


colnames(nCoV_PDLIM2_T[[]])
head(nCoV_PDLIM2_T[[]]$NFkBSum)
head(nCoV_PDLIM2_T[[]]$NFkBAverage)

# fetch data
nCoV_PDLIM2_T_NFkB <- FetchData(object = nCoV_PDLIM2_T, vars = c("group", "expression_group", "PDLIM2", "NFkBSum","NFkBAverage")) 
write.csv (nCoV_PDLIM2_T_NFkB, file="nCoV_PDLIM2_T_NFkB.csv")

######## sum/average expression of top upregulated NF-kB targetd genes in macrophage/T
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T, layer = "data")

# average expression for your gene set per cell
gene_set <- c("CXCL10", "CCL2", "CD38", "MX1", "IRF7", "IRF4", "IFI44L", "TNFSF10", "CCL4", "GBP1", "BCL3", "PRF1", "HIF1A", "RELB", "NFKBIZ", "CCL3", "FAS", "FOS", "NFKBIA", "FASLG", "SAT1", "TNFAIP3", "NFKB1", "CEBPD", "TAP1", "GZMB")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T$NFkBSum2 <- gene_set_sum
nCoV_PDLIM2_T$NFkBAverage2 <- avg_expr_per_cell


colnames(nCoV_PDLIM2_T[[]])
head(nCoV_PDLIM2_T[[]]$NFkBSum2)
head(nCoV_PDLIM2_T[[]]$NFkBAverage2)

# fetch data
nCoV_PDLIM2_T_NFkB_Sum <- FetchData(object = nCoV_PDLIM2_T, vars = c("group", "expression_group", "PDLIM2", "NFkBSum2","NFkBAverage2")) 
write.csv (nCoV_PDLIM2_T_NFkB_Sum, file="nCoV_PDLIM2_T_NFkB_Sum.csv")

# subset HC+O to a new seurat object, divide Covid to Pdlim2 low and high group, calculate DE genes comparing to HC
# subset Macrophage to MacrophageHCO
nCoV_PDLIM2_MacrophageHCO <- subset(nCoV_PDLIM2_Macrophage, subset = group %in% c("HC", "O"))

#Error in validObject(object = x) :  invalid class "Assay" object: slots in class definition but not in object: "assay.orig"

# update the nCoV_PDLIM2_T object
nCoV_PDLIM2_T_update <- UpdateSeuratObject(nCoV_PDLIM2_T)

# subset T_update to T_update_HCO
nCoV_PDLIM2_T_update_HCO <- subset(nCoV_PDLIM2_T_update, subset = group %in% c("HC", "O"))

# Subset O cells
O_cells <- WhichCells(nCoV_PDLIM2_T_update_HCO, expression = group == "O")

# Get PDLIM2 expression in O cells
pdl_expr <- FetchData(nCoV_PDLIM2_T_update_HCO, vars = "PDLIM2")[O_cells, , drop = FALSE]

# Determine median expression
pdl_median <- median(pdl_expr$PDLIM2)

# Classify disease cells into 'low' and 'high'
low_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 <= pdl_median]
high_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 > pdl_median]

# Assign new identities
nCoV_PDLIM2_T_update_HCO$HCO <- "other"  # default
nCoV_PDLIM2_T_update_HCO$HCO[low_cells] <- "disease_low"
nCoV_PDLIM2_T_update_HCO$HCO[high_cells] <- "disease_high"
nCoV_PDLIM2_T_update_HCO$HCO[nCoV_PDLIM2_T_update_HCO$group == "HC"] <- "HC"

# fetch data
nCoV_PDLIM2_T_update_HCO_PDLIM2 <- FetchData(object = nCoV_PDLIM2_T_update_HCO, vars = c("group", "expression_group", "PDLIM2", "HCO")) 
write.csv (nCoV_PDLIM2_T_update_HCO_PDLIM2, file="nCoV_PDLIM2_T_update_HCO_PDLIM2.csv")

# Set identities
Idents(nCoV_PDLIM2_T_update_HCO) <- nCoV_PDLIM2_T_update_HCO$HCO

# Run differential expression
DEG <- FindMarkers(nCoV_PDLIM2_T_update_HCO,
                   ident.1 = "disease_low",
                   ident.2 = "HC",
                   group.by = "HCO",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_HCO_Pdlim2Low-HC.csv")

#sum expression of STAT target genes in nCoV_PDLIM2_T_update_HCO
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T_update_HCO, layer = "data")

# sum/average expression for your gene set per cell
gene_set <- c("CDKN1A", "CISH", "EGF", "ERAS", "HRAS", "IL6", "IL6R", "JAK2", "KRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK14", "MAPK3", "MRAS", "MYC", "NRAS", "PDGFB", "PIAS3", "PTPN6", "RAF1", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RRAS", "RRAS2", "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "STAT3", "TYK2", "VEGFA", "BCL2", "BMP6", "BMPR1A", "BMPR1B", "BMPR2", "CDC25A", "CSF2RB", "CXCR1", "CXCR2", "EGFR", "FGF2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT1", "FLT4", "GHR", "HGF", "IFNAR1", "IFNLR1", "IGF1", "IGF1R", "IGF2R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RL1", "IL1RL2", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2", "IL27RA", "IL2RA", "IL2RB", "IL2RG", "IL31RA", "IL3RA", "IL4R", "IL5RA", "IL6ST", "IL7R", "IL9", "IL9R", "INSR", "KDR", "MAP2K4", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K20", "MAP3K21", "MAP3K9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK8", "MAPK9", "NDUFA13", "NGFR", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PDGFRB", "PIM1", "PTPN2", "RAC1", "SRC", "TGFA", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TGFBR3", "TNFRSF11A", "IFNG", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IL12B", "TNF", "NOS2", "IRF1", "GBP1-5", "IDO1", "PSMB8", "TAP1", "TAP2", "HLA-DRA", "ICAM1", "FCGR1A", "BST2", "IL10", "LIF", "CSF2", "PDGFA", "PDGFC", "PDGFD", "CCL2", "SETDB2", "MMP9", "CCL17", "CCL22", "ARG1", "MRC1", "CHIT1", "IL4", "PPARG", "GATA3", "AKT1", "AKT2", "AKT3", "BCL2L1", "CCKBR", "CEBPB", "FOS", "GAST", "GNAQ", "GRB2", "JAK1", "JAK3", "JUN", "MTOR", "NFKB1", "NFKB2", "PIAS1", "PIAS2", "PIAS4", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PTPN1", "PTPN11", "REL", "RELA", "RELB", "SHC1", "SOS1", "SOS2", "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T_update_HCO$STATSum <- gene_set_sum
nCoV_PDLIM2_T_update_HCO$STATAverage <- avg_expr_per_cell


colnames(nCoV_PDLIM2_T_update_HCO[[]])
head(nCoV_PDLIM2_T_update_HCO[[]]$STATSum)
head(nCoV_PDLIM2_T_update_HCO[[]]$STATAverage)

# fetch data
nCoV_PDLIM2_T_update_HCO_STAT <- FetchData(object = nCoV_PDLIM2_T_update_HCO, vars = c("group", "expression_group", "PDLIM2", "STATSum","STATAverage", "HCO")) 
write.csv (nCoV_PDLIM2_T_update_HCO_STAT, file="nCoV_PDLIM2_T_update_HCO_STAT-PDLIM2.csv")

# Phagosome formation 
###### sum expression of phagosome formation genes in macrophage
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# sum expression for your gene set per cell
gene_set1 <- c("ABHD3", "ACKR1", "ACKR3", "ACKR4", "ACTR2", "ACTR3", "ADCYAP1R1", "ADGRA1", "ADGRA2", "ADGRA3", "ADGRB1", "ADGRB2", "ADGRB3", "ADGRD1", "ADGRD2", "ADGRE1", "ADGRE2", "ADGRE3", "ADGRE5", "ADGRF1", "ADGRF2", "ADGRF3", "ADGRF4", "ADGRF5", "ADGRG1", "ADGRG2", "ADGRG3", "ADGRG4", "ADGRG5", "ADGRG6", "ADGRG7", "ADGRL1", "ADGRL3", "ADGRL4", "ADGRV1", "ADORA1", "ADORA2A", "ADORA2B", "ADORA3", "ADRA1A", "ADRA1B", "ADRA1D", "ADRA2A", "ADRA2B", "ADRA2C", "ADRB1", "ADRB2", "ADRB3", "AGTR1", "AGTR2", "AKT1", "AKT2", "AKT3", "AP1B1", "AP1G1", "AP1G2", "AP1M1", "AP1M2", "AP1S1", "AP1S2", "AP1S3", "APBB1IP", "APLNR", "ARF6", "ARFIP2", "ARPC1A", "ARPC1B", "ARPC2", "ARPC3", "ARPC4", "ARPC5", "ARPC5L", "AVPR1A", "AVPR1B", "AVPR2", "BCAR1", "BCL10", "BDKRB1", "BDKRB2", "BRS3", "C3AR1", "C5AR1", "C5AR2", "CALCR", "CALCRL", "CAMP", "CARD9", "CASR", "CCKAR", "CCKBR", "CCR1", "CCR10", "CCR2", "CCR3", "CCR4", "CCR5", "CCR6", "CCR7", "CCR8", "CCR9", "CCRL2", "CD14", "CD209", "CD36", "CDC42", "CELSR1", "CELSR2", "CELSR3", "CFL1", "CFL2", "CHRM1", "CHRM2", "CHRM3", "CHRM4", "CHRM5", "CLEC4D", "CLEC4E", "CLEC7A", "CLIP1", "CMKLR1", "CNR1", "CNR2", "COLEC12", "CR1", "CR2", "CRHR1", "CRHR2", "CRKL", "CX3CR1", "CXCR1", "CXCR6", "CYSLTR1", "CYSLTR2", "DIAPH1", "DOCK1", "DRD1", "DRD2", "DRD3", "DRD4", "DRD5", "EDNRA", "EDNRB", "ELMO1", "ELMO2", "ELMO3", "ERAS", "F2R", "F2RL1", "F2RL2", "F2RL3", "FCAMR", "FCAR", "FCER1A", "FCER1G", "FCER2", "FCGR1A", "FCGR2A", "FCGR2B", "FCGR2C", "FFAR1", "FFAR2", "FFAR3", "FFAR4", "FGR", "FN1", "FPR1", "FPR2", "FPR3", "FSHR", "FYN", "FZD1", "FZD10", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6", "FZD7", "FZD8", "FZD9", "GAB2", "GABBR1", "GABBR2", "GALR1", "GALR2", "GALR3", "GCGR", "GHRHR", "GHSR", "GIPR", "GLP1R", "GLP2R", "GNRHR", "GPER1", "GPLD1", "GPR101", "GPR107", "GPR108", "GPR119", "GPR12", "GPR132", "GPR135", "GPR137", "GPR137B", "GPR137C", "GPR139", "GPR141", "GPR142", "GPR143", "GPR146", "GPR148", "GPR149", "GPR15", "GPR150", "GPR151", "GPR152", "GPR153", "GPR155", "GPR156", "GPR157", "GPR158", "GPR160", "GPR161", "GPR162", "GPR17", "GPR171", "GPR173", "GPR174", "GPR176", "GPR179", "GPR18", "GPR180", "GPR182", "GPR183", "GPR19", "GPR20", "GPR21", "GPR22", "GPR25", "GPR26", "GPR27", "GPR3", "GPR31", "GPR32", "GPR33", "GPR34", "GPR35", "GPR37", "GPR37L1", "GPR39", "GPR4", "GPR45", "GPR50", "GPR52", "GPR55", "GPR6", "GPR61", "GPR62", "GPR63", "GPR65", "GPR68", "GPR75", "GPR78", "GPR82", "GPR83", "GPR84", "GPR85", "GPR87", "GPR88", "GPRC5A", "GPRC5B", "GPRC5C", "GPRC5D", "GPRC6A", "GRB2", "GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8", "GRPR")

gene_set2 <- c("HCAR1", "HCAR2", "HCAR3", "HCK", "HCRTR1", "HCRTR2", "HMOX1", "HRAS", "HRH1", "HRH2", "HRH3", "HRH4", "HTR1A", "HTR1B", "HTR1D", "HTR1E", "HTR1F", "HTR2A", "HTR2B", "HTR2C", "HTR4", "HTR5A", "HTR6", "HTR7", "IGHA1", "IGHE", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM", "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC7", "ITGA1", "ITGA10", "ITGA11", "ITGA2", "ITGA2B", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8", "ITGA9", "ITGAD", "ITGAE", "ITGAL", "ITGAM", "ITGAV", "ITGAX", "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8", "ITPR1", "ITPR2", "ITPR3", "JCHAIN", "KISS1R", "KRAS", "LAT", "LBP", "LCAT", "LCK", "LGR4", "LGR5", "LHCGR", "LIMK1", "LIMK2", "LPAR1", "LPAR2", "LPAR3", "LPAR4", "LPAR5", "LPAR6", "LTB4R", "LTB4R2", "LYN", "MALT1", "MAP2K1", "MAP2K2", "MAP2K3", "MAP2K4", "MAP2K5", "MAP2K6", "MAP2K7", "MAPK1", "MAPK12", "MAPK15", "MAPK3", "MAPK4", "MAPK6", "MAPK7", "MARCKS", "MARCO", "MAS1", "MAS1L", "MC1R", "MC2R", "MC3R", "MC4R", "MC5R", "MCHR1", "MLNR", "MRAS", "MRC1", "MRC2", "MSR1", "MTNR1A", "MTNR1B", "MTOR", "MYD88", "MYH1", "MYH10", "MYH11", "MYH13", "MYH14", "MYH2", "MYH3", "MYH4", "MYH6", "MYH7", "MYH7B", "MYH8", "MYH9", "MYL1", "MYL12A", "MYL12B", "MYL2", "MYL3", "MYL4", "MYL5", "MYL6", "MYL6B", "MYL7", "MYL9", "MYLK", "MYLK2", "MYLK3", "MYO10", "MYO18A", "MYO18B", "MYO1A", "MYO1H", "MYO5C", "NAPEPLD", "NMBR", "NMUR1", "NMUR2", "NPBWR1", "NPBWR2", "NPFFR1", "NPFFR2", "NPY1R", "NPY2R", "NPY4R", "NPY5R", "NRAS", "NTSR1", "NTSR2", "OCRL", "OPN1LW", "OPN1SW", "OPRD1", "OPRK1", "OPRL1", "OPRM1", "OXGR1", "OXTR", "P2RY1", "P2RY10", "P2RY11", "P2RY12", "P2RY13", "P2RY14", "P2RY2", "P2RY4", "P2RY6", "PAK1", "PAK2", "PAK3", "PAK4", "PAK5", "PAK6", "PI4KA", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PIKFYVE", "PIP4K2A", "PIP4K2B", "PIP4K2C", "PIP5K1A", "PIP5K1B", "PIP5K1C", "PIP5KL1", "PLA2G10", "PLA2G12A", "PLA2G12B", "PLA2G1B", "PLA2G2A", "PLA2G2C", "PLA2G2D", "PLA2G2E", "PLA2G2F", "PLA2G3", "PLA2G4A", "PLA2G4B", "PLA2G4C", "PLA2G4D", "PLA2G4E", "PLA2G4F", "PLA2G5", "PLA2G6", "PLA2G7", "PLA2R1", "PLAAT1", "PLAAT2", "PLAAT3", "PLAAT4", "PLAAT5", "PLB1", "PLCG1", "PLCG2", "PLD1", "PLD2", "PLD3", "PLD4", "PLD5", "PLD6", "PNPLA2", "PNPLA3", "PNPLA4", "PNPLA8", "PRDX6", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCI", "PRKCQ", "PRKCZ", "PRKD1", "PRKD3", "PRLHR", "PROKR1", "PTAFR", "PTGDR", "PTGDR2", "PTGER1", "PTGER2", "PTGER3", "PTGER4", "PTGFR", "PTGIR", "PTH1R", "PTK2", "PTK2B", "QRFPR", "RAC1", "RAC2", "RAC3", "RAF1", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RAPGEF1", "RAPGEF3", "RAPGEF4", "RASD1", "RASD2", "RASGRP1", "RGR", "RHO", "RHOA", "ROCK1", "ROCK2", "RPS6KB1", "RPS6KB2", "RRAS", "RRAS2", "RRH", "RXFP1", "RXFP2", "RXFP3", "RXFP4", "S1PR1", "S1PR2", "S1PR3", "S1PR4", "S1PR5", "SAP130", "SCARA3", "SCARA5", "SCTR", "SLC52A1", "SLC52A2", "SMO", "SOS1", "SOS2", "SPHK1", "SPHK2", "SRC", "SSTR1", "SSTR2", "SSTR3", "SSTR4", "SSTR5", "SUCNR1", "SYK", "TAAR1", "TAAR2", "TAAR5", "TAAR8", "TACR1", "TACR2", "TACR3", "TAS1R1", "TAS1R2", "TAS1R3", "TBXA2R", "TIMD4", "TLN1", "TLN2", "TLR1", "TLR10", "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TLR9", "TRHR", "TSHR", "TTN", "UTS2R", "VAV1", "VAV2", "VAV3", "VIPR1", "VIPR2", "VN1R1", "VN1R2", "VN1R4", "VN1R5", "VTN", "WAS", "WASF1", "WASF2", "WASL", "WIPF1", "XCR1", "YES1")  # your list

common_genes1 <- intersect(rownames(norm_expr), gene_set1)

common_genes2 <- intersect(rownames(norm_expr), gene_set2)

gene_set_sum1 <- colSums(norm_expr[common_genes1, ])
gene_set_sum2 <- colSums(norm_expr[common_genes2, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$PhagosomeSum1 <- gene_set_sum1
nCoV_PDLIM2_Macrophage_update_HCO$PhagosomeSum2 <- gene_set_sum2

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$PhagosomeSum1)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_Phagosome <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "GAPDH", "NFkBGillmorSum","NFkBGillmorAverage", "PhagosomeSum1", "PhagosomeSum2")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_Phagosome, file="nCoV_PDLIM2_Macrophage_update_HCO_Phagosome.csv")

save (list=ls(), file="GSE145926.RData")

# phagosome formation gene list from IPA_comparison_pathway_PhagosomeFormation_COPdlim2High-HC
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# sum expression for your gene set per cell
gene_set <- c("ADGRL1", "ADRB2", "ARPC1A", "ARPC4", "C3", "C3AR1", "CCR1", "CD36", "CLEC4D", "CMKLR1", "CYSLTR2", "ELMO2", "F2RL2", "FCAR", "FCGR3B", "FFAR2", "FPR2", "FPR3", "FZD1", "FZD5", "GPR132", "GPR141", "GPR25", "GPR35", "GPR84", "GPRC5C", "HRH1", "HRH2", "ITGB7", "LPAR3", "LTB4R2", "LYN", "MAP2K6", "MARCKS", "MYD88", "MYO10", "MYO18A", "P2RY12", "P2RY14", "P2RY6", "PIK3R3", "PIK3R6", "PLA2G7", "PLCG1", "PTAFR", "PTGER2", "PTK2", "SLC52A1", "SPHK1", "SSTR2", "TIMD4", "TLR1", "TLR2", "TLR3", "TTN", "VAV3")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$PhagosomeSum3 <- gene_set_sum

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$PhagosomeSum3)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_Phagosome <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "GAPDH", "NFkBGillmorSum","NFkBGillmorAverage", "PhagosomeSum1", "PhagosomeSum2", "PhagosomeSum3")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_Phagosome, file="nCoV_PDLIM2_Macrophage_update_HCO_Phagosome.csv")

# Qiagen "Production of Nitric Oxide and Reactive Oxygen Species in Macrophages" gene list
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# sum expression for your gene set per cell
gene_set<- c("AKT1", "AKT2", "AKT3", "ALB", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3", "APOC4", "APOD", "APOE", "APOF", "APOL1", "APOM", "ARG2", "CAT", "CDC42", "CHUK", "CLU", "CREBBP", "CYBA", "CYBB", "DIRAS3", "ELP1", "FNBP1", "FOS", "HOXA10", "IFNG", "IFNGR1", "IFNGR2", "IKBKB", "IKBKE", "IKBKG", "IL4", "IRF1", "IRF8", "JAK1", "JAK2", "JAK3", "JUN", "LPA", "LYZ", "MAP2K1", "MAP2K4", "MAP2K7", "MAP3K1", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K13", "MAP3K14", "MAP3K15", "MAP3K2", "MAP3K3", "MAP3K4", "MAP3K5", "MAP3K6", "MAP3K7", "MAP3K8", "MAP3K9", "MAPK1", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "MAPK3", "MAPK8", "MAPK9", "MPO", "NCF1", "NCF2", "NCF4", "NFKB1", "NFKB2", "NFKBIA", "NFKBIB", "NFKBID", "NFKBIE", "NGFR", "NOS2", "ORM1", "ORM2", "PCYOX1", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PLCG1", "PLCG2", "PON1", "PPARA", "PPM1J", "PPM1L", "PPP1CA", "PPP1CB", "PPP1CC", "PPP1R10", "PPP1R11", "PPP1R12A", "PPP1R14A", "PPP1R14B", "PPP1R14C", "PPP1R14D", "PPP1R3A", "PPP1R3C", "PPP1R3D", "PPP1R7", "PPP2CA", "PPP2CB", "PPP2R1A", "PPP2R1B", "PPP2R2A", "PPP2R2B", "PPP2R2C", "PPP2R3A", "PPP2R3B", "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D", "PPP2R5E", "PRKCA", "PRKCB", "PRKCD", "PRKCE", "PRKCG", "PRKCH", "PRKCI", "PRKCQ", "PRKCZ", "PRKD1", "PRKD3", "PTPA", "PTPN6", "RAC1", "RAC2", "RAC3", "RAP1A", "RAP1B", "RBP4", "REL", "RELA", "RELB", "RHOA", "RHOB", "RHOBTB1", "RHOBTB2", "RHOC", "RHOD", "RHOF", "RHOG", "RHOH", "RHOJ", "RHOQ", "RHOT1", "RHOT2", "RHOU", "RHOV", "RND1", "RND2", "RND3", "S100A8", "SAA4", "SERPINA1", "SIRPA", "SPI1", "STAT1", "TLR2", "TLR4", "TNF", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TYK2")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$NOROSinMac <- gene_set_sum

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$NOROSinMac)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_NOROSinMac <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "NOROSinMac")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_NOROSinMac, file="nCoV_PDLIM2_Macrophage_update_HCO_NOROSinMac.csv")

# Sum expression of apoptosis in macrophage
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# Reactome apoptosis molecules
gene_set <- c("ARHGAP10", "TLR4", "TNFRSF10A", "PSMD11", "PSMD12", "DFFA", "DNM1L", "PSMD14", "SDCBP", "APAF1", "TNFRSF10B", "PSMA7", "OGT", "TP73", "CFLAR", "SEPTIN4", "PSMD3", "DAPK3", "BCL2L11", "CHMP2A", "GAS2", "PRKN", "OPA1", "GSDME", "FLOT1", "DFFB", "IL1A", "IL1B", "LMNA", "", "TP53", "GSN", "H1-0", "HSP90AA1", "ELANE", "", "CD14", "VIM", "HMGB1", "UBB", "UBC", "GZMB", "H1-4", "BCL2", "MAPT", "IRF1", "CDH1", "IRF2", "DSP", "H1-5", "H1-3", "H1-2", "PSMC3", "PSMB1", "LMNB1", "APC", "FAS", "PSMA1", "PSMA2", "PSMA3", "PSMA4", "HMGB2", "YWHAQ", "MAPK3", "PSMA5", "PSMB4", "PSMB6", "PSMB5", "MAPK1", "CASP1", "NMT1", "AKT1", "AKT2", "YWHAB", "SFN", "DSG3", "CTNNB1", "ADD1", "PSMC2", "STAT3", "CASP3", "DCC", "PSMC4", "MAPK8", "FASLG", "PPP3CC", "PSMD8", "FNTA", "CASP4", "PSMB3", "PSMB2", "TNFSF10", "BCAP31", "PSMD7", "BMX", "CASP5", "KPNA1", "DAPK1", "CASP7", "CASP9", "CASP6", "BID", "GSDMD", "SEM1", "PSMA6", "YWHAG", "PSMC1", "PSMC5", "YWHAE", "PSMC6", "RPS27A", "UBA52", "PPP3R1", "YWHAZ", "DYNLL1", "UBE2L3", "XIAP", "CYCS", "E2F1", "SATB1", "DSG1", "H1-1", "PRKCQ", "YWHAH", "PTK2", "PRKCD", "C1QBP", "TJP1", "BAX", "BCL2L1", "TRAF2", "FADD", "PAK2", "PSMD2", "ROCK1", "BIRC3", "BIRC2", "RIPK1", "TP53BP2", "PMAIP1", "SPTAN1", "PKP1", "IL18", "DSG2", "TFDP1", "TFDP2", "FLOT2", "CASP8", "KPNB1", "PSMD6", "PLEC", "TRADD", "ADRM1", "CDC37", "BAK1", "OCLN", "UNC5A", "", "TMED7-TICAM2", "TICAM1", "UNC5B", "CDKN2A", "MLKL", "PDCD6IP", "CHMP7", "BAD", "CHMP4C", "OMA1", "PELI1", "DYNLL2", "CHMP6", "APIP", "ITCH", "PPP1R13B", "BMF", "", "PSMB7", "PSMD1", "BBC3", "CHMP4A", "UACA", "TP63", "CHMP4B", "CLSPN", "AVEN", "DIABLO", "STK26", "TJP2", "DAPK2", "DBNL", "APPL1", "ACIN1", "STUB1", "PSMD13", "CHMP2B", "AKT3", "CARD8", "CHMP3", "RIPK3", "MAGED1", "STK24", "LY96")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$apoptosis <- gene_set_sum

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$apoptosis)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_apoptosis <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "apoptosis")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_apoptosis, file="nCoV_PDLIM2_Macrophage_update_HCO_CellDeathMac.csv")


# IPA diseases and Functions "cell death of macrophages molecules" gene list
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# sum expression for your gene set per cell
gene_set<- c("ADAM8", "AIM2", "BCL2L1", "BIRC3", "C3AR1", "CDKN1A", "CGAS", "DDIT3", "EIF2AK2", "FAS", "GBP1", "GPR132", "IL10", "IL6R", "IRF1", "IRF2", "ISG15", "LDLR", "MEFV", "MUC5AC", "MYD88", "NLRP3", "PFKFB3", "PLAT", "PMAIP1", "PML", "RIPK1", "SOD2", "ST6GAL1", "STAT1", "TICAM1", "TICAM2", "TLR2", "TLR3", "TNF", "TNFSF10", "TRIM21", "ZBP1")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$CellDeathMac <- gene_set_sum

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$CellDeathMac)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_CellDeathMac <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "CellDeathMac")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_CellDeathMac, file="nCoV_PDLIM2_Macrophage_update_HCO_CellDeathMac.csv")

# sum expression of reactome autophagy genes in Macrophage
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

#reactome autophagy molecules
gene_set <- c("AMBRA1","ARL13B","ATG10","ATG101","ATG12","ATG13","ATG14","ATG16L1","ATG16L2","ATG3","ATG4A","ATG4B","ATG4C","ATG4D","ATG5","ATG7","ATG9A","ATG9B","ATM","BECN1","CETN1","CFTR","CHMP2A","CHMP2B","CHMP3","CHMP4A","CHMP4B","CHMP4C","CHMP6","CHMP7","CSNK2A1","CSNK2A2","CSNK2B","DYNC1H1","DYNC1I1","DYNC1I2","DYNC1LI1","DYNC1LI2","DYNLL1","DYNLL2","EEF1A1","EPAS1","FUNDC1","GABARAP","GABARAPL1","GABARAPL2","GFAP","HBB","HDAC6","HSF1","HSP90AA1","HSP90AB1","HSPA8","IFT88","LAMP2","LAMTOR1","LAMTOR2","LAMTOR3","LAMTOR4","LAMTOR5","MAP1LC3A","MAP1LC3B","MAP1LC3C","MFN1","MFN2","MLST8","MTERF3","MTMR14","MTMR3","MTOR","MVB12A","MVB12B","NBR1","OPTN","PARK7","PCNT","PEX5","PGAM5","PIK3C3","PIK3R4","PINK1","PLIN2","PLIN3","PRKAA1","PRKAA2","PRKAB1","PRKAB2","PRKAG1","PRKAG2","PRKAG3","PRKN","RB1CC1","RHEB","RNASE1","RPS27A","RPTOR","RRAGA","RRAGB","RRAGC","RRAGD","SLC38A9","SQSTM1","SRC","TBK1","TOMM20","TOMM22","TOMM40","TOMM5","TOMM6","TOMM7","TOMM70","TSC1","TSC2","TSG101","TUBA1A","TUBA1B","TUBA1C","TUBA3C","TUBA3D","TUBA3E","TUBA4A","TUBA4B","TUBA8","TUBAL3","TUBB1","TUBB2A","TUBB2B","TUBB3","TUBB4A","TUBB4B","TUBB6","TUBB8","TUBB8B","UBA52","UBAP1","UBB","UBC","UBE2D2","UBE2D3","UBE2L3","UBE2N","UBE2V1","ULK1","USP30","UVRAG","VCP","VDAC1","VDAC2","VDAC3","VIM","VPS28","VPS37A","VPS37B","VPS37C","VPS37D","WDR45","WDR45B","WIPI1","WIPI2")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$autophagy <- gene_set_sum

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$autophagy)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_autophagy <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "autophagy")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_autophagy, file="nCoV_HCO_autophagy.csv")

# DEG between Pdlim2LowO and Pdlim2HighO in macrophage
colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
table(nCoV_PDLIM2_Macrophage_update_HCO[[]]$HCO)

DEG <- FindMarkers(nCoV_PDLIM2_Macrophage_update_HCO,
                   ident.1 = "disease_low",
                   ident.2 = "disease_high",
                   group.by = "HCO",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_Macrophage_DiseasePdlim2L-H.csv")

#sum expression of "T cell exhaustion signaling pathway"
# IPA pathway "T cell exhaustion signaling pathway" DEGs
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T_update_HCO, layer = "data")

# sum expression for your gene set per cell
gene_set<- c("CD247", "CD3E", "FCER1G", "FOS", "HLA-DMA", "HLA-DQB1", "HLA-DRA", "HLA-E", "IL12RB2", "IRF4", "JAK1", "JAK2", "LAG3", "LGALS9", "NFATC2", "NFATC3", "PIK3CD", "PLCG1", "PPP2R1A", "PPP2R2A", "PTPN11", "STAT1", "STAT2", "STAT4", "TBX21", "TGFBR2", "TGFBR3", "TNFRSF14", "TRAV8-2", "TRAV8-3", "TRBV20-1", "TRGV2", "TRGV3", "TRGV5", "ZAP70", "CD86", "HLA-DQA1", "TRBV9")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T_update_HCO$exhaustion <- gene_set_sum

colnames(nCoV_PDLIM2_T_update_HCO[[]])
head(nCoV_PDLIM2_T_update_HCO[[]]$exhaustion)

# fetch data
nCoV_PDLIM2_T_update_HCO_exhaustion <- FetchData(object = nCoV_PDLIM2_T_update_HCO, vars = c("group", "expression_group", "PDLIM2", "exhaustion")) 
write.csv (nCoV_PDLIM2_T_update_HCO_exhaustion, file="nCoV_PDLIM2_T_update_HCO_exhaustion.csv")

# sum expression of senescence genes in T
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T_update_HCO, layer = "data")

# Reactome senescence molecules
gene_set <- c("CBX4", "E2F3", "UBE2C", "MAP2K7", "KDM6B", "MDM4", "H2BC12", "NBN", "EED", "VENTX", "CBX6", "MAP4K4", "CCNE2", "FOS", "IFNB1", "IL1A", "TP53", "H2AC4", "IL6", "JUN", "RB1", "H2BC11", "H1-0", "SP1", "H2AB1", "UBB", "UBC", "CXCL8", "H1-4", "TXN", "CDK4", "ETS1", "ETS2", "H2AX", "H1-5", "H1-3", "H1-2", "HMGA1", "CEBPB", "NFKB1", "CCNA2", "H2AC7", "LMNB1", "H2BC17", "CCNE1", "CDK2", "MAPK3", "MAPK1", "CDC27", "H2BC3", "BMI1", "CDKN1A", "STAT3", "ID1", "CDKN2A", "CDKN2B", "CDKN2C", "MAPK8", "MAPK9", "MAP2K4", "CDKN1B", "MAP2K3", "MAPKAPK2", "MRE11", "ERF", "UBE2D1", "RPS6KA3", "UBE2E1", "MAP2K6", "HMGA2", "MAPK10", "HIRA", "TERF1", "CDKN2D", "H2BC12L", "H2BC5", "ANAPC15", "H4C4", "H2BC7", "RPS27A", "UBA52", "H3C3", "PHC1", "CCNA1", "H3-3A", "CDK6", "MDM2", "E2F1", "H1-1", "RELA", "RING1", "RBBP4", "CDC16", "MAPK7", "ATM", "TFDP1", "TFDP2", "E2F2", "CBX2", "SUZ12", "RPS6KA2", "RPS6KA1", "TERF2", "MAPK11", "EZH2", "IGFBP7", "MAPK14", "RBBP7", "MAPKAPK3", "H3-4", "UBE2S", "H2AC20", "H2BC21", "H2AC19", "H3C14", "H2AZ2", "MAPKAPK5", "PHC2", "H2BC26", "MINK1", "CDKN2A", "TNRC6A", "PHC3", "CDC26", "RAD50", "KAT5", "H2AC6", "H2BC9", "H2BC1", "ACD", "ANAPC16", "SCMH1", "EHMT2", "EP400", "RNF2", "MAP3K5", "H2BC15", "H2AC14", "H2BC14", "H2BC13", "TINF2", "H2AJ", "ANAPC1", "EHMT1", "AGO3", "CBX8", "MOV10", "TNRC6C", "AGO4", "UBN1", "POT1", "TERF2IP", "ANAPC11", "CDC23", "ANAPC7", "ANAPC5", "ANAPC4", "ANAPC2", "TNIK", "AGO1", "FZR1", "ANAPC10", "TNRC6B", "ASF1A", "CABIN1")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T_update_HCO$senescence <- gene_set_sum

colnames(nCoV_PDLIM2_T_update_HCO[[]])
head(nCoV_PDLIM2_T_update_HCO[[]]$senescence)

# fetch data
nCoV_PDLIM2_T_update_HCO_senescence <- FetchData(object = nCoV_PDLIM2_T_update_HCO, vars = c("group", "expression_group", "PDLIM2", "senescence")) 
write.csv (nCoV_PDLIM2_T_update_HCO_senescence, file="nCoV_PDLIM2_T_update_HCO_senescence.csv")

# sum expression of qiagen FXR/RXR activation genelist in T
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T_update_HCO, layer = "data")

# qiagen FXR/RXR activation genelist
gene_set <- c("A1BG","ABCB11","ABCB4","ABCC2","ABCG5","ABCG8","AGT","AHSG","AKT1","AKT2","AKT3","ALB","AMBP","APOA1","APOA2","APOA4","APOB","APOC1","APOC2","APOC3","APOC4","APOD","APOE","APOF","APOH","APOL1","APOM","BAAT","C3","C9","CAMP","CETP","CLU","CREBBP","CYP27A1","CYP7A1","CYP8B1","FABP6","FASN","FBP1","FETUB","FGA","FGF19","FGFR4","FOXA1","FOXA2","FOXA3","FOXO1","G6PC1","G6PC2","G6PC3","GC","HNF1A","HNF4A","HPR","HPX","IL18","IL1A","IL1B","IL1F10","IL1RN","IL33","IL36A","IL36B","IL36G","IL36RN","IL37","INS","ITIH4","KNG1","LCAT","LIPC","LPL","MAP2K4","MAPK10","MAPK12","MAPK8","MAPK9","MLXIPL","MTTP","NR0B2","NR1H3","NR1H4","NR1I2","NR5A2","ORM1","ORM2","PCK2","PCYOX1","PKLR","PLTP","PON1","PON3","PPARA","PPARG","PPARGC1A","RARA","RBP4","RXRA","SAA1","SAA2","SAA4","SCARB1","SDC1","SERPINA1","SERPINF1","SERPINF2","SLC10A1","SLC10A2","SLC22A7","SLC27A5","SLC4A2","SLC51A","SLC51B","SLCO1B1","SLCO1B3","SULT2A1","TF","TNF","TTR","UGT2B4","VLDLR","VTN")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T_update_HCO$FXR <- gene_set_sum

colnames(nCoV_PDLIM2_T_update_HCO[[]])
head(nCoV_PDLIM2_T_update_HCO[[]]$FXR)

# fetch data
nCoV_PDLIM2_T_update_HCO_FXR <- FetchData(object = nCoV_PDLIM2_T_update_HCO, vars = c("group", "expression_group", "PDLIM2", "FXR")) 
write.csv (nCoV_PDLIM2_T_update_HCO_FXR, file="nCoV_PDLIM2_T_update_HCO_FXR.csv")

save (list=ls(), file="GSE145926.RData")
 