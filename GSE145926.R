#Python error when openning RStudio on Macbook
#check the running version of RStudio
RStudio.Version()$version

#problem solved with steps below
#First, verify that Python is installed correctly. Open Terminal and run:/usr/bin/python3 --version
# open RStudio and no more Python error

library(dplyr)
library(Seurat)
library(patchwork)

# download the annotated Seurat rds file at http://cells.ucsc.edu/covid19-balf/nCoV.rds.
# convert RDS file to seurat project
nCoV <- readRDS("nCoV.rds")

# check seurat object
nCov[[]]
nCoV@meta.data
nCoV[[]]
head(nCoV, n=2)
head(nCoV@assays$RNA@counts)
head(nCoV@assays$RNA@data)
head(nCoV@assays$RNA@scale.data)
head(nCoV, n=2)
tail(nCoV, n=2)
colnames(nCoV@meta.data)
dim(nCoV)
dim(nCoV@assays$RNA@counts)
tail(nCoV@assays$RNA@counts)
nCoV@assays$RNA@counts["PDLIM2",]
nCoV@assays$RNA@counts["PDLIM2",400:405]
dim(nCoV@assays$RNA@counts)
dim(cell_annotation)

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



# subset HC+O to a new seurat object, divide Covid to Pdlim2 low and high group, calculate DE genes comparing to HC
# subset Macrophage to MacrophageHCO
nCoV_PDLIM2_MacrophageHCO <- subset(nCoV_PDLIM2_Macrophage, subset = group %in% c("HC", "O"))

#Error in validObject(object = x) :  invalid class "Assay" object: slots in class definition but not in object: "assay.orig"

# update the nCoV_PDLIM2_Macrophage object
nCoV_PDLIM2_Macrophage_update <- UpdateSeuratObject(nCoV_PDLIM2_Macrophage)

# subset T_update to T_update_HCO
nCoV_PDLIM2_Macrophage_update_HCO <- subset(nCoV_PDLIM2_Macrophage_update, subset = group %in% c("HC", "O"))

# Subset O cells
O_cells <- WhichCells(nCoV_PDLIM2_Macrophage_update_HCO, expression = group == "O")

# Get PDLIM2 expression in O cells
pdl_expr <- FetchData(nCoV_PDLIM2_Macrophage_update_HCO, vars = "PDLIM2")[O_cells, , drop = FALSE]

# Determine median expression
pdl_median <- median(pdl_expr$PDLIM2)

# Classify disease cells into 'low' and 'high'
low_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 <= pdl_median]
high_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 > pdl_median]

# Assign new identities
nCoV_PDLIM2_Macrophage_update_HCO$HCO <- "other"  # default
nCoV_PDLIM2_Macrophage_update_HCO$HCO[low_cells] <- "disease_low"
nCoV_PDLIM2_Macrophage_update_HCO$HCO[high_cells] <- "disease_high"
nCoV_PDLIM2_Macrophage_update_HCO$HCO[nCoV_PDLIM2_Macrophage_update_HCO$group == "HC"] <- "HC"

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_PDLIM2 <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "HCO")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_PDLIM2, file="nCoV_PDLIM2_Macrophage_update_HCO_PDLIM2.csv")

# Set identities
Idents(nCoV_PDLIM2_Macrophage_update_HCO) <- nCoV_PDLIM2_Macrophage_update_HCO$HCO

# Run differential expression
DEG <- FindMarkers(nCoV_PDLIM2_Macrophage_update_HCO,
                   ident.1 = "disease_low",
                   ident.2 = "HC",
                   group.by = "HCO",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_HCO_Pdlim2Low-HC.csv")

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

# sum/average expression of all NF-kB targeted genes in macrophage
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# NF-kB target genes From the Gillmore Lab of Boston University
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$NFkBSum <- gene_set_sum
nCoV_PDLIM2_Macrophage_update_HCO$NFkBAverage <- avg_expr_per_cell


colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$NFkBSum)
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$NFkBAverage)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_NFkB <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "NFkBSum","NFkBAverage")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_NFkB, file="nCoV_PDLIM2_Macrophage_update_HCO_NFkB.csv")

# sum/average expression of all NF-kB targeted genes in T
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T_update_HCO, layer = "data")

# NF-kB target genes From the Gillmore Lab of Boston University
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_T_update_HCO$NFkBSum <- gene_set_sum
nCoV_PDLIM2_T_update_HCO$NFkBAverage <- avg_expr_per_cell


colnames(nCoV_PDLIM2_T_update_HCO[[]])
head(nCoV_PDLIM2_T_update_HCO[[]]$NFkBSum)
head(nCoV_PDLIM2_T_update_HCO[[]]$NFkBAverage)

# fetch data
nCoV_PDLIM2_T_update_HCO_NFkB <- FetchData(object = nCoV_PDLIM2_T_update_HCO, vars = c("group", "expression_group", "PDLIM2", "NFkBSum","NFkBAverage")) 
write.csv (nCoV_PDLIM2_T_update_HCO_NFkB, file="nCoV_PDLIM2_T_update_HCO_NFkB.csv")

#sum expression of STAT3 target genes in nCoV_PDLIM2_Macrophage_update_HCO
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# STAT3 signaling gene list from Qiagen
gene_set <- c("CDKN1A", "CISH", "EGF", "ERAS", "HRAS", "IL6", "IL6R", "JAK2", "KRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK14", "MAPK3", "MRAS", "MYC", "NRAS", "PDGFB", "PIAS3", "PTPN6", "RAF1", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RRAS", "RRAS2", "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "STAT3", "TYK2", "VEGFA", "BCL2", "BMP6", "BMPR1A", "BMPR1B", "BMPR2", "CDC25A", "CSF2RB", "CXCR1", "CXCR2", "EGFR", "FGF2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT1", "FLT4", "GHR", "HGF", "IFNAR1", "IFNLR1", "IGF1", "IGF1R", "IGF2R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RL1", "IL1RL2", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2", "IL27RA", "IL2RA", "IL2RB", "IL2RG", "IL31RA", "IL3RA", "IL4R", "IL5RA", "IL6ST", "IL7R", "IL9", "IL9R", "INSR", "KDR", "MAP2K4", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K20", "MAP3K21", "MAP3K9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK8", "MAPK9", "NDUFA13", "NGFR", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PDGFRB", "PIM1", "PTPN2", "RAC1", "SRC", "TGFA", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TGFBR3", "TNFRSF11A", "IFNG", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IL12B", "TNF", "NOS2", "IRF1", "GBP1-5", "IDO1", "PSMB8", "TAP1", "TAP2", "HLA-DRA", "ICAM1", "FCGR1A", "BST2", "IL10", "LIF", "CSF2", "PDGFA", "PDGFC", "PDGFD", "CCL2", "SETDB2", "MMP9", "CCL17", "CCL22", "ARG1", "MRC1", "CHIT1", "IL4", "PPARG", "GATA3", "AKT1", "AKT2", "AKT3", "BCL2L1", "CCKBR", "CEBPB", "FOS", "GAST", "GNAQ", "GRB2", "JAK1", "JAK3", "JUN", "MTOR", "NFKB1", "NFKB2", "PIAS1", "PIAS2", "PIAS4", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PTPN1", "PTPN11", "REL", "RELA", "RELB", "SHC1", "SOS1", "SOS2", "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$STATSum <- gene_set_sum
nCoV_PDLIM2_Macrophage_update_HCO$STATAverage <- avg_expr_per_cell


colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$STATSum)
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$STATAverage)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_STAT <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "STATSum","STATAverage", "HCO")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_STAT, file="nCoV_PDLIM2_Macrophage_update_HCO_STAT-PDLIM2.csv")


#sum expression of STAT3 target genes in nCoV_PDLIM2_T_update_HCO
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_T_update_HCO, layer = "data")

# STAT3 signaling gene list from Qiagen
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


# sum expression of Production of NO & ROS in Macrophages gene list
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# Qiagen Production of NO & ROS in Macrophages gene list
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


# sum expression of Coronavirus Pathogenesis Pathway-Infection molecules
# Extract the normalized expression values
norm_expr <- GetAssayData(nCoV_PDLIM2_Macrophage_update_HCO, layer = "data")

# Reactome Coronavirus Pathogenesis Pathway-Infection molecules
gene_set <- c("ACE2", "AGRN", "ANO1", "ANO10", "ANO2", "ANO3", "ANO4", "ANO5", "ANO6", "ANO7", "ANO8", "ANO9", "BECN1", "CANX", "CHMP2A", "CHMP2B", "CHMP3", "CHMP4A", "CHMP4B", "CHMP4C", "CHMP6", "CHMP7", "CSNK1A1", "CTSL", "DAD1", "DDOST", "DDX5", "EDEM2", "FURIN", "FUT8", "GALNT1", "GANAB", "GOLGA7", "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6", "GSK3A", "GSK3B", "HAVCR1", "HSPG2", "ISCU", "MAGT1", "MAN1B1", "MAN2A1", "MAP1LC3B", "MGAT1", "MGAT2", "MGAT4A", "MGAT4B", "MGAT4C", "MGAT5", "MOGS", "NRP1", "OST4", "OSTC", "PARP10", "PARP14", "PARP16", "PARP4", "PARP6", "PARP8", "PARP9", "PIK3C3", "PIK3R4", "PRKCSH", "PRMT1", "RB1", "RPN1", "RPN2", "RPS27A", "SDC1", "SDC2", "SDC3", "SDC4", "SRPK1", "SRPK2", "ST3GAL1", "ST3GAL2", "ST3GAL3", "ST3GAL4", "ST6GAL1", "ST6GALNAC2", "ST6GALNAC3", "ST6GALNAC4", "STT3A", "STT3B", "SUMO1", "TMEM258", "TMPRSS2", "TUSC3", "UBA52", "UBB", "UBC", "UBE2I", "UVRAG", "VCP", "VHL", "ZCRB1", "ZDHHC11", "ZDHHC2", "ZDHHC20", "ZDHHC3", "ZDHHC5", "ZDHHC8", "ZDHHC9")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
nCoV_PDLIM2_Macrophage_update_HCO$SARSCoV2_Infection <- gene_set_sum

colnames(nCoV_PDLIM2_Macrophage_update_HCO[[]])
head(nCoV_PDLIM2_Macrophage_update_HCO[[]]$SARSCoV2_Infection)

# fetch data
nCoV_PDLIM2_Macrophage_update_HCO_SARSCoV2_Infection <- FetchData(object = nCoV_PDLIM2_Macrophage_update_HCO, vars = c("group", "expression_group", "PDLIM2", "SARSCoV2_Infection")) 
write.csv (nCoV_PDLIM2_Macrophage_update_HCO_SARSCoV2_Infection, file="nCoV_PDLIM2_Macrophage_update_HCO_SARSCoV2_Infection.csv")

save (list=ls(), file="GSE145926.RData")
 