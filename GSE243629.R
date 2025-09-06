# download the following sample files at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243629: "IC_1", "IC_2", "IA_1", "IA_2", "PI_1", "PI_2", "PI_3", "PHC_1", "PHC_2", "HC_1", "HC_2"

#Load Required Libraries
library(dplyr)
library(Matrix)
library(Seurat)
library(patchwork)

#Read and Create Seurat Objects for Each Sample
samples <- list("IC_1", "IC_2", "IA_1", "IA_2", "PI_1", "PI_2", "PI_3", "PHC_1", "PHC_2", "HC_1", "HC_2")  # Update with your sample names
seurat_list <- list() # blank Seurat list

for (sample in samples) {
  counts <- ReadMtx(mtx = paste0(sample, "/matrix.mtx"),
                    features = paste0(sample, "/genes.tsv"),
                    cells = paste0(sample, "/barcodes.tsv"))
  seurat_obj <- CreateSeuratObject(counts = counts, project = sample)
  seurat_obj$sample <- sample  # Add sample identity as metadata
  seurat_list[[sample]] <- seurat_obj # store in seurat list
}

#Merge the Seurat Objects
merged_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)

# QC
merged_seurat[["percent.mt"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")
head(merged_seurat[[]])

# Visualize QC metrics as a violin plot
VlnPlot(merged_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#cell filtering
merged_seurat_FilteredCells <- subset(merged_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)

#Normalization
merged_seurat_FilteredCells <- NormalizeData(merged_seurat_FilteredCells)

#FindVariableFeatures
merged_seurat_FilteredCells <- FindVariableFeatures(merged_seurat_FilteredCells, selection.method = "vst", nfeatures = 2000)

#Scaling
merged_seurat_FilteredCells <- ScaleData(merged_seurat_FilteredCells)

#linear dimensional reduction
merged_seurat_FilteredCells <- RunPCA(merged_seurat_FilteredCells, features = VariableFeatures(object = merged_seurat_FilteredCells))
ElbowPlot(merged_seurat_FilteredCells)

#clustering
merged_seurat_FilteredCells <- FindNeighbors(merged_seurat_FilteredCells, dims = 1:13)
merged_seurat_FilteredCells <- FindClusters(merged_seurat_FilteredCells, resolution = 0.5)

#UMAP
merged_seurat_FilteredCells <- RunUMAP(merged_seurat_FilteredCells, dims = 1:13)
DimPlot(merged_seurat_FilteredCells, reduction = "umap", label=TRUE)
head(Idents(merged_seurat_FilteredCells))

#FindAllMarkers
merged_seurat_FilteredCells_AllMarkers <- FindAllMarkers(merged_seurat_FilteredCells) # error: data layers are not joined. Please run JoinLayers

# JoinLayers
Layers(merged_seurat_FilteredCells)
merged_seurat_FilteredCells <- JoinLayers(
  object = merged_seurat_FilteredCells, 
  layers = c("data.IC_1", "data.IC_2", "data.IA_1", "data.IA_2", "data.PI_1", "data.PI_2", "data.PI_3", "data.PHC_1", "data.PHC_2", "data.HC_1", "data.HC_2"),  # Specify layers to merge
  new.layer = "data"     # Name of the new layer
)
#after running script above, no error but layers not joined

# Manually Combine the Layers
# Extract and combine the layers
merged_matrix <- cbind(
  ic_1 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.IC_1"),
  ic_2 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.IC_2"),
  ia_1 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.IA_1"),
  ia_2 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.IA_2"),
  pi_1 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.PI_1"),
  pi_2 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.PI_2"),
  pi_3 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.PI_3"),
  phc_1 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.PHC_1"),
  phc_2 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.PHC_2"),
  hc_1 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.HC_1"),
  hc_2 <- GetAssayData(merged_seurat, assay = "RNA", slot = "counts.HC_2")
)

# Add the combined data back to the Seurat object
merged_seurat[["RNA"]]$counts <- merged_matrix
head(merged_seurat)
merged_seurat@assays$RNA@layers$counts[1:5,1:5]


# check layers
Layers(merged_seurat)

dim(merged_seurat@assays$RNA@layers$counts.IC_1)
dim(merged_seurat@assays$RNA@layers$counts)

# Remove multiple layers
merged_seurat[["RNA"]]$counts.IC_1 <- NULL
merged_seurat[["RNA"]]$counts.IC_2 <- NULL
merged_seurat[["RNA"]]$counts.IA_1 <- NULL
merged_seurat[["RNA"]]$counts.IA_2 <- NULL
merged_seurat[["RNA"]]$counts.PI_1 <- NULL
merged_seurat[["RNA"]]$counts.PI_2 <- NULL
merged_seurat[["RNA"]]$counts.PI_3 <- NULL
merged_seurat[["RNA"]]$counts.PHC_1 <- NULL
merged_seurat[["RNA"]]$counts.PHC_2 <- NULL
merged_seurat[["RNA"]]$counts.HC_1 <- NULL
merged_seurat[["RNA"]]$counts.HC_2 <- NULL
# Repeat for other layers as necessary

# redo some steps above to create merged_seurat_FilteredCells till the script below
# FindAllMarkers
merged_seurat_FilteredCells_AllMarkers <- FindAllMarkers(merged_seurat_FilteredCells) # error: error in evaluating the argument 'x' in selecting a method for function 'rowSums': vector memory limit of 16.0 Gb reached, see mem.maxVSize()

# mem.maxVSize()
mem.maxVSize()
mem.maxVSize(vsize = 50000)
gc()

# redo the FindAllMarkers
merged_seurat_FilteredCells_AllMarkers <- FindAllMarkers(merged_seurat_FilteredCells) # no error
head(merged_seurat_FilteredCells_AllMarkers)
tail(merged_seurat_FilteredCells_AllMarkers)
merged_seurat_FilteredCells_AllMarkers$cluster
unique(merged_seurat_FilteredCells_AllMarkers$cluster)

merged_seurat_FilteredCells_AllMarkers_cluster0 <- 
merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="0",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster0, "merged_seurat_FilteredCells_AllMarkers_cluster0.csv")

merged_seurat_FilteredCells_AllMarkers_cluster1 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="1",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster1, "merged_seurat_FilteredCells_AllMarkers_cluster1.csv")

merged_seurat_FilteredCells_AllMarkers_cluster2 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="2",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster2, "merged_seurat_FilteredCells_AllMarkers_cluster2.csv")

merged_seurat_FilteredCells_AllMarkers_cluster3 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="3",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster3, "merged_seurat_FilteredCells_AllMarkers_cluster3.csv")

merged_seurat_FilteredCells_AllMarkers_cluster4 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="4",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster4, "merged_seurat_FilteredCells_AllMarkers_cluster4.csv")


merged_seurat_FilteredCells_AllMarkers_cluster5 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="5",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster5, "merged_seurat_FilteredCells_AllMarkers_cluster5.csv")

merged_seurat_FilteredCells_AllMarkers_cluster6 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="6",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster6, "merged_seurat_FilteredCells_AllMarkers_cluster6.csv")

merged_seurat_FilteredCells_AllMarkers_cluster7 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="7",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster7, "merged_seurat_FilteredCells_AllMarkers_cluster7.csv")


merged_seurat_FilteredCells_AllMarkers_cluster8 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="8",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster8, "merged_seurat_FilteredCells_AllMarkers_cluster8.csv")

merged_seurat_FilteredCells_AllMarkers_cluster9 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="9",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster9, "merged_seurat_FilteredCells_AllMarkers_cluster9.csv")

merged_seurat_FilteredCells_AllMarkers_cluster10 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="10",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster10, "merged_seurat_FilteredCells_AllMarkers_cluster10.csv")

merged_seurat_FilteredCells_AllMarkers_cluster11 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="11",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster11, "merged_seurat_FilteredCells_AllMarkers_cluster11.csv")

merged_seurat_FilteredCells_AllMarkers_cluster12 <- 
  merged_seurat_FilteredCells_AllMarkers[merged_seurat_FilteredCells_AllMarkers$cluster=="12",]

write.csv (merged_seurat_FilteredCells_AllMarkers_cluster12, "merged_seurat_FilteredCells_AllMarkers_cluster12.csv")




# Save
save (list=ls(), file="GSE243629.RData")

# Load
load("GSE243629.RData")

# subset PDLIM2+ cells from GSE_seurat_FilteredCells
merged_seurat_FilteredCells_Pdlim2Active <- subset(merged_seurat_FilteredCells, subset = PDLIM2 >0)

# FetechData
colnames(merged_seurat_FilteredCells_Pdlim2Active@meta.data)
unique(merged_seurat_FilteredCells_Pdlim2Active@meta.data$sample)
FetchData <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active, vars = c("sample", "seurat_clusters", "PDLIM2")) 
write.csv (FetchData, "merged_seurat_FilteredCells_Pdlim2Active.csv")

# save
save (list=ls(), file="GSE243629.RData")

# subset monocytes from GSE_seurat_FilteredCells_Pdlim2Active
merged_seurat_FilteredCells_Pdlim2Active_Monocyte <- subset(merged_seurat_FilteredCells_Pdlim2Active, subset = seurat_clusters %in% c("6","7","9","10"))
colnames(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]])
unique(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$seurat_clusters)
unique(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$sample)

# DEG Monocyte Influenza vs Healthy
DEG <- FindMarkers(merged_seurat_FilteredCells_Pdlim2Active_Monocyte,
                   ident.1 = c("IC_1","IC_2","IA_1","IA_2","PI_1","PI_2","PI_3"),
                   ident.2 = c("PHC_1", "PHC_2", "HC_1","HC_2"),
                   group.by = "sample",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_monocyte_wilcox.csv")

# subset T cells from GSE_seurat_FilteredCells_Pdlim2Active
merged_seurat_FilteredCells_Pdlim2Active_T <- subset(merged_seurat_FilteredCells_Pdlim2Active, subset = seurat_clusters %in% c("4","5"))
colnames(merged_seurat_FilteredCells_Pdlim2Active_T[[]])
unique(merged_seurat_FilteredCells_Pdlim2Active_T[[]]$seurat_clusters)
unique(merged_seurat_FilteredCells_Pdlim2Active_T[[]]$sample)

# DEG T Influenza vs Healthy
DEG <- FindMarkers(merged_seurat_FilteredCells_Pdlim2Active_T,
                   ident.1 = c("IC_1","IA_1","IA_2","PI_1"),
                   ident.2 = c("PHC_1", "PHC_2", "HC_1","HC_2"),
                   group.by = "sample",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_T_wilcox.csv")

# subset NK cells from GSE_seurat_FilteredCells_Pdlim2Active
merged_seurat_FilteredCells_Pdlim2Active_NK <- subset(merged_seurat_FilteredCells_Pdlim2Active, subset = seurat_clusters %in% c("3"))
colnames(merged_seurat_FilteredCells_Pdlim2Active_NK[[]])
unique(merged_seurat_FilteredCells_Pdlim2Active_NK[[]]$seurat_clusters)
unique(merged_seurat_FilteredCells_Pdlim2Active_NK[[]]$sample)

# DEG T Influenza vs Healthy
DEG <- FindMarkers(merged_seurat_FilteredCells_Pdlim2Active_NK,
                   ident.1 = c("IC_1","IC_2", "IA_1","IA_2","PI_1"),
                   ident.2 = c("PHC_2", "HC_1","HC_2"),
                   group.by = "sample",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_NK_wilcox.csv")

# load
load ("GSE243629.RData")

# subset neutrophil from GSE_seurat_FilteredCells_Pdlim2Active
merged_seurat_FilteredCells_Pdlim2Active_neutrophil <- subset(merged_seurat_FilteredCells_Pdlim2Active, subset = seurat_clusters %in% c("0","1"))
colnames(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]])
unique(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$seurat_clusters)
unique(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$sample)

# DEG Monocyte Influenza vs Healthy
DEG <- FindMarkers(merged_seurat_FilteredCells_Pdlim2Active_Monocyte,
                   ident.1 = c("IC_2","IA_1","IA_2","PI_1","PI_2","PI_3"),
                   ident.2 = c("PHC_1", "PHC_2", "HC_1","HC_2"),
                   group.by = "sample",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_neutrophil_wilcox.csv")

# load GSE234629.RData
load("GSE243629.RData")

# classify PDLIM2+ Neutrophils to High/Low level Pdlim2 groups
# Fetch PDLIM2 expression data and calculate its level at certain percentile
# Choose the gene of interest
gene_of_interest <- "PDLIM2"  # Replace with your gene

# Extract expression levels
gene_expr <- FetchData(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(merged_seurat_FilteredCells_Pdlim2Active_neutrophil$expression_group)

colnames(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]])

unique(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$orig.ident)

# Fetch above High/Low PDLIM2 neutrophil groups
merged_seurat_FilteredCells_Pdlim2Active_neutrophil_PDLIM2Percentile <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = c("sample", "expression_group", "PDLIM2")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_neutrophil_PDLIM2Percentile, file="merged_seurat_FilteredCells_Pdlim2Active_neutrophil_PDLIM2Percentile.csv")

# DEG neutrophils with High vs Low on PDLIM2 level
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil@active.ident)
merged_seurat_FilteredCells_Pdlim2Active_neutrophil <- SetIdent(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, value = merged_seurat_FilteredCells_Pdlim2Active_neutrophil@meta.data$expression_group)
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil@active.ident)
DEG <- FindAllMarkers(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, group.by = "expression_group")

write.csv (DEG, file = "DEG_Neutrophil_PDLIM2_H-L_ANOVA.csv")

# subset PDLIM2+ monocytes to High/Low level Pdlim2 groups
# Fetch PDLIM2 expression data and calculate its level at certain percentile
# Choose the gene of interest
gene_of_interest <- "PDLIM2"  # Replace with your gene

# Extract expression levels
gene_expr <- FetchData(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, vars = gene_of_interest)

# Define thresholds (e.g., using quantiles)
high_threshold <- quantile(gene_expr[,1], 0.50)

# Assign group labels
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$expression_group <- ifelse(
  gene_expr[,1] >= high_threshold, "High", "Low")

# Check new metadata
table(merged_seurat_FilteredCells_Pdlim2Active_Monocyte$expression_group)

colnames(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]])

unique(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$orig.ident)

# Fetch above High/Low PDLIM2 neutrophil groups
merged_seurat_FilteredCells_Pdlim2Active_Monocyte_PDLIM2Percentile <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_Monocyte, vars = c("sample", "expression_group", "PDLIM2")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_Monocyte_PDLIM2Percentile, file="merged_seurat_FilteredCells_Pdlim2Active_Monocyte_PDLIM2Percentile.csv")

# DEG Monocytes with High vs Low on PDLIM2 level
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte@active.ident)
merged_seurat_FilteredCells_Pdlim2Active_Monocyte <- SetIdent(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, value = merged_seurat_FilteredCells_Pdlim2Active_Monocyte@meta.data$expression_group)
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte@active.ident)
DEG <- FindAllMarkers(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, group.by = "expression_group")

write.csv (DEG, file = "Seurat_DEG_Monocyte_PDLIM2_H-L.csv")


#divide disease cells into Pdlim2 low and High groups, and calculate their DEGs to Ctrl, and sum expression of STAT signaling target genes

# Assign new identities
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2 <- "other"  # default
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "PI_1"] <- "Flu"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "PI_2"] <- "Flu"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "PI_3"] <- "Flu"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "IA_1"] <- "Flu"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "IA_2"] <- "Flu"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "IC_2"] <- "Flu"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "HC_1"] <- "Control"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "HC_1"] <- "Control"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "HC_2"] <- "Control"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "PHC_1"] <- "Control"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[merged_seurat_FilteredCells_Pdlim2Active_neutrophil$sample == "PHC_2"] <- "Control"

# Subset disease cells
disease_cells <- WhichCells(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, expression = DiseasePdlim2 == "Flu")

# Get PDLIM2 expression in disease cells
pdl_expr <- FetchData(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = "PDLIM2")[disease_cells, , drop = FALSE]

# Determine median expression
pdl_median <- median(pdl_expr$PDLIM2)

# Classify disease cells into 'low' and 'high'
low_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 <= pdl_median]
high_cells <- rownames(pdl_expr)[pdl_expr$PDLIM2 > pdl_median]

# Assign new identities
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[low_cells] <- "disease_low"
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2[high_cells] <- "disease_high"

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_neutrophil_DiseasePdlim2 <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = c("sample", "expression_group", "PDLIM2", "DiseasePdlim2")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2, file="merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2.csv")

# Set identities
Idents(merged_seurat_FilteredCells_Pdlim2Active_neutrophil) <- merged_seurat_FilteredCells_Pdlim2Active_neutrophil$DiseasePdlim2

# Run differential expression
DEG <- FindMarkers(merged_seurat_FilteredCells_Pdlim2Active_neutrophil,
                   ident.1 = "disease_high",
                   ident.2 = "Control",
                   group.by = "DiseasePdlim2",
                   test.use = "wilcox")  # Default is Wilcoxon test
write.csv (DEG, file = "DEG_neutrophil_Pdlim2High-Ctrl.csv")

###### sum/average expression of NF-kB targetd genes in neutrophil
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, layer = "data")

# NF-kB target genes From the Gillmore Lab of Boston University
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$NFkBSum <- gene_set_sum
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$NFkBAverage <- avg_expr_per_cell


colnames(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$NFkBSum)
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$NFkBAverage)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_neutrophil_NFkB <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = c("sample", "expression_group", "PDLIM2", "NFkBSum","NFkBAverage")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_neutrophil_NFkB, file="merged_seurat_FilteredCells_Pdlim2Active_neutrophil_NFkB.csv")

###### sum/average expression of NF-kB targetd genes in monocyte
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, layer = "data")

# NF-kB target genes From the Gillmore Lab of Boston University
gene_set <- c("ABCA1", "B2M", "BACE1", "BCL2A1", "BCL3", "CCL2", "CCL3", "CCL5", "CCND3", "CD274", "CD38", "CD40", "CD44", "CD83", "CDKN1A", "CREB3", "CXCL10", "DPYD", "EGR1", "FCGRT", "FOS", "GBP1", "GSTP1", "HIF1A", "HLA-B", "HPSE", "IDO1", "IER2", "IFI44L", "IL1RN", "IL27", "IRF1", "IRF7", "JUNB", "LAMB2", "LGALS3", "LYZ", "MBP", "MDK", "MX1", "NFKBIZ", "NRG1", "NUAK2", "PIM1", "PLK3", "PSMB9", "PSME2", "PTEN", "PTPN1", "RIPK2", "S100A4", "S100A6", "SAT1", "SERPINA1", "SERPINB1", "SH3BGRL", "SLC3A2", "SPI1", "TAP1", "TGM1", "TLR2", "TNFAIP2", "TNFAIP3", "TNFSF10", "TNFSF13B", "UPP1", "VIM", "A4", "ABCB1", "ABCB4", "ABCB9", "ABCC6", "ABCG5", "ABCG8", "ADAM19", "ADH1A", "ADORA1", "ADORA2A", "ADRA2B", "AFP", "AGER", "AGT", "AHCTF1", "AICDA", "AKR1C1", "ALOX12", "ALOX5", "AMACR", "AMH", "ANGPT1", "APOBEC2", "APOC3", "APOD", "APOE", "AQP4", "AR", "ARFRP1", "ART1", "ASPH", "ASS1", "ATP1A2", "BAX", "BCL2", "BCL2L1", "BCL2L11", "BDKRB1", "BDNF", "BGN", "BLIMP1 /PRDM1", "BLNK", "BLR1", "BMI1", "BMP2", "BMP4", "BNIP3", "BRCA2", "BTK", "C3", "C4A", "C4BPA", "C69", "CALCB", "CASP4", "CAV1", "CCL1", "CCL15", "CCL17", "CCL19", "CCL20", "CCL22", "CCL23", "CCL28", "CCL4", "CCND1", "CCND2", "CCR5", "CCR7", "CD209", "CD3G", "CD40LG", "CD48", "CD54", "CD80", "CD86", "CDK6", "CDX1", "CEBPD", "CFB", "CFLAR", "CGM3", "CHI3L1", "CIDEA", "CLDN2", "COL1A2", "CR2", "CRP", "CSF1", "CSF2", "CSF3", "CTSB", "CTSL1", "CXCL1", "CXCL11", "CXCL3", "CXCL5", "CXCL6", "CXCL9", "CYP19A1", "CYP27B1", "CYP2C11", "CYP2E1", "CYP7B1", "DCTN4", "DEFB2", "DIO2", "DMP1", "DNASE1L2", "DUSP1", "E2F3", "EBI3", "EDN1", "EGFR", "ELF3", "ENG", "ENO2", "EPHA1", "EPO", "ERBB2", "ERVWE1", "F11R", "F3", "F8", "FABP6", "FAM148A", "FAS", "FASLG", "FCER2", "FCER2/CD23", "FGF8", "FN1", "FSTL3", "FTH1", "G6PC", "G6PD", "GAD1", "GADD45B", "GATA3", "GCLC", "GCLM", "GCNT1", "GJB1", "GNAI2", "GNB2L1", "GNRH2", "GRIN1", "GRIN2A", "GRM2", "GRO-beta", "GRO-gamma", "GUCY1A2", "GZMB", "HAMP", "HAS1", "HBE1", "HBZ", "HGF", "HLA-G", "HMGN1", "HMOX1", "HOXA9", "HSD11B2", "HSD17B8", "HSP90AA1", "ICOS", "IER3", "IFNB1", "IFNG", "IGFBP1", "IGFBP2", "IGHE", "IGHG1", "IGHG2", "IGHG4", "IGKC", "Iigp1", "IL10", "IL11", "IL12A", "IL12B", "IL13", "IL15", "IL17", "IL1A", "IL1B", "IL2", "IL23A", "IL2RA", "IL6", "IL8", "IL8RA", "IL8RB", "IL9", "INHBA", "IRF2", "IRF4", "JMJD3", "KC", "KCNK5", "KCNN2", "KISS1", "KITLG", "KLF10", "KLK3", "KLRA1", "KRT15", "KRT3", "KRT5", "KRT6B", "LBP", "LCN2", "LEF1", "LIPG", "LTA", "LTB", "LTF", "MADCAM1", "MAP4K1", "MMP1", "MMP3", "MMP9", "MT3", "MTHFR", "MUC2", "MYB", "MYC", "MYLK", "MYOZ1", "NCAM", "NFKB1", "NFKB2", "NFKBIA", "NFKBIE", "NGFB", "NK4", "NLRP2", "NOD2", "NOS1", "NOS2A", "NOX1", "NPY1R", "NQO1", "NR3C1", "NR4A2", "OLR1", "OPN1SW", "OPRD1", "OPRM1", "ORM1", "Osterix", "OXTR", "PAFAH2", "PAX8", "PDE7A", "PDGFB", "PDYN", "PENK", "PGK1", "PGLYRP1", "PGR", "PI3KAP1", "PIGF", "pIgR", "PIK3CA", "PLA2", "PLAU", "PLCD1", "POMC", "PPARGC1B", "PPP5C", "PRF1", "PRKACA", "PRKCD", "PRL", "PSME1", "PTAFR", "PTGDS", "PTGES", "PTGIS", "PTGS2", "PTHLH", "PTPN13", "PTS", "PTX3", "PYCARD", "RAG1", "RAG2", "RBBP4", "Rdh1", "Rdh7", "REL", "RELB", "REV3L", "S100A10", "SAA1", "SAA2", "SAA3", "SCNN1A", "SDC4", "SELE", "SELP", "SELS", "SENP2", "SERPINA2", "SERPINA3", "SERPINE1, PAI-1", "SERPINE2", "SKALP, PI3", "SKP2", "SLC11A2", "SLC16A1", "SLC6A6", "Slfn2", "SNAI1", "SOD1", "SOD2", "SOX9", "SPATA19", "SPP1", "ST6GAL1", "ST8SIA1", "STAT5A", "SUPV3L1", "TACR1", "TAPBP", "TCRB", "TERT", "TF", "TFEC", "TFF3", "TFPI2", "TGM2", "THBS1", "THBS2", "TICAM1", "TIFA", "TLR9", "TNC", "TNF", "TNFRSF1B", "TNFRSF4", "TNFRSF9", "TNFSF15", "TNIP1", "TNIP3", "TP53", "TRAF1", "TRAF2", "TREM1", "TRPC1", "TWIST1", "UBE2M", "UCP2", "UGCGL1", "UPK1B", "VCAM1", "VEGFC", "VPS53", "WNT10B", "WT1", "XDH", "XIAP", "YY1", "ZNF366")  # your list
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$NFkBSum <- gene_set_sum
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$NFkBAverage <- avg_expr_per_cell


colnames(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$NFkBSum)
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$NFkBAverage)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_Monocyte_NFkB <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_Monocyte, vars = c("sample", "expression_group", "PDLIM2", "NFkBSum","NFkBAverage")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_Monocyte_NFkB, file="merged_seurat_FilteredCells_Pdlim2Active_Monocyte_NFkB.csv")


#sum expression of STAT3 target genes in neutrophil
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, layer = "data")

# STAT3 signaling gene list from Qiagen
gene_set <- c("CDKN1A", "CISH", "EGF", "ERAS", "HRAS", "IL6", "IL6R", "JAK2", "KRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK14", "MAPK3", "MRAS", "MYC", "NRAS", "PDGFB", "PIAS3", "PTPN6", "RAF1", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RRAS", "RRAS2", "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "STAT3", "TYK2", "VEGFA", "BCL2", "BMP6", "BMPR1A", "BMPR1B", "BMPR2", "CDC25A", "CSF2RB", "CXCR1", "CXCR2", "EGFR", "FGF2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT1", "FLT4", "GHR", "HGF", "IFNAR1", "IFNLR1", "IGF1", "IGF1R", "IGF2R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RL1", "IL1RL2", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2", "IL27RA", "IL2RA", "IL2RB", "IL2RG", "IL31RA", "IL3RA", "IL4R", "IL5RA", "IL6ST", "IL7R", "IL9", "IL9R", "INSR", "KDR", "MAP2K4", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K20", "MAP3K21", "MAP3K9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK8", "MAPK9", "NDUFA13", "NGFR", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PDGFRB", "PIM1", "PTPN2", "RAC1", "SRC", "TGFA", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TGFBR3", "TNFRSF11A", "IFNG", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IL12B", "TNF", "NOS2", "IRF1", "GBP1-5", "IDO1", "PSMB8", "TAP1", "TAP2", "HLA-DRA", "ICAM1", "FCGR1A", "BST2", "IL10", "LIF", "CSF2", "PDGFA", "PDGFC", "PDGFD", "CCL2", "SETDB2", "MMP9", "CCL17", "CCL22", "ARG1", "MRC1", "CHIT1", "IL4", "PPARG", "GATA3", "AKT1", "AKT2", "AKT3", "BCL2L1", "CCKBR", "CEBPB", "FOS", "GAST", "GNAQ", "GRB2", "JAK1", "JAK3", "JUN", "MTOR", "NFKB1", "NFKB2", "PIAS1", "PIAS2", "PIAS4", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PTPN1", "PTPN11", "REL", "RELA", "RELB", "SHC1", "SOS1", "SOS2", "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$STATSum <- gene_set_sum
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$STATAverage <- avg_expr_per_cell


colnames(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$STATSum)
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$STATAverage)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_neutrophil_STAT <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = c("sample", "expression_group", "PDLIM2", "STATSum","STATAverage", "DiseasePdlim2")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_neutrophil_STAT, file="merged_seurat_FilteredCells_Pdlim2Active_neutrophil_STAT-PDLIM2.csv")

#sum expression of STAT3 target genes in Monocyte
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, layer = "data")

# STAT3 signaling gene list from Qiagen
gene_set <- c("CDKN1A", "CISH", "EGF", "ERAS", "HRAS", "IL6", "IL6R", "JAK2", "KRAS", "MAP2K1", "MAP2K2", "MAPK1", "MAPK14", "MAPK3", "MRAS", "MYC", "NRAS", "PDGFB", "PIAS3", "PTPN6", "RAF1", "RALA", "RALB", "RAP1A", "RAP1B", "RAP2A", "RAP2B", "RASD1", "RASD2", "RRAS", "RRAS2", "SOCS1", "SOCS2", "SOCS3", "SOCS4", "SOCS5", "SOCS6", "SOCS7", "STAT3", "TYK2", "VEGFA", "BCL2", "BMP6", "BMPR1A", "BMPR1B", "BMPR2", "CDC25A", "CSF2RB", "CXCR1", "CXCR2", "EGFR", "FGF2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "FLT1", "FLT4", "GHR", "HGF", "IFNAR1", "IFNLR1", "IGF1", "IGF1R", "IGF2R", "IL10RA", "IL10RB", "IL11RA", "IL12RB1", "IL12RB2", "IL13RA1", "IL13RA2", "IL15RA", "IL17RA", "IL17RB", "IL17RC", "IL17RD", "IL17RE", "IL18R1", "IL18RAP", "IL1A", "IL1B", "IL1R1", "IL1R2", "IL1RL1", "IL1RL2", "IL20RA", "IL20RB", "IL21R", "IL22RA1", "IL22RA2", "IL27RA", "IL2RA", "IL2RB", "IL2RG", "IL31RA", "IL3RA", "IL4R", "IL5RA", "IL6ST", "IL7R", "IL9", "IL9R", "INSR", "KDR", "MAP2K4", "MAP3K10", "MAP3K11", "MAP3K12", "MAP3K20", "MAP3K21", "MAP3K9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK8", "MAPK9", "NDUFA13", "NGFR", "NTRK1", "NTRK2", "NTRK3", "PDGFRA", "PDGFRB", "PIM1", "PTPN2", "RAC1", "SRC", "TGFA", "TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TGFBR3", "TNFRSF11A", "IFNG", "CXCL9", "CXCL10", "CXCL11", "CCL5", "IL12B", "TNF", "NOS2", "IRF1", "GBP1-5", "IDO1", "PSMB8", "TAP1", "TAP2", "HLA-DRA", "ICAM1", "FCGR1A", "BST2", "IL10", "LIF", "CSF2", "PDGFA", "PDGFC", "PDGFD", "CCL2", "SETDB2", "MMP9", "CCL17", "CCL22", "ARG1", "MRC1", "CHIT1", "IL4", "PPARG", "GATA3", "AKT1", "AKT2", "AKT3", "BCL2L1", "CCKBR", "CEBPB", "FOS", "GAST", "GNAQ", "GRB2", "JAK1", "JAK3", "JUN", "MTOR", "NFKB1", "NFKB2", "PIAS1", "PIAS2", "PIAS4", "PIK3C2A", "PIK3C2B", "PIK3C2G", "PIK3C3", "PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "PIK3R3", "PIK3R4", "PIK3R5", "PIK3R6", "PTPN1", "PTPN11", "REL", "RELA", "RELB", "SHC1", "SOS1", "SOS2", "STAT1", "STAT2", "STAT4", "STAT5A", "STAT5B", "STAT6")
common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

avg_expr_per_cell <- Matrix::colMeans(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$STATSum <- gene_set_sum
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$STATAverage <- avg_expr_per_cell


colnames(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$STATSum)
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$STATAverage)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_Monocyte_STAT <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_Monocyte, vars = c("sample", "expression_group", "PDLIM2", "STATSum","STATAverage", "DiseasePdlim2")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_Monocyte_STAT, file="merged_seurat_FilteredCells_Pdlim2Active_Monocyte_STAT-PDLIM2.csv")


# sum expression of pathway "neutrophil degranulation" genes
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, layer = "data")

# sum expression for your gene set per cell
# Reactome neutrophil degranulation molecules
gene_set1 <- c("TARM1", "SNAP23", "PSMD11", "PSMD12", "SIRPB1", "PGRMC1", "QSOX1", "MANBA", "PSMD14", "SDCBP", "DDX3X", "RNASET2", "FCN1", "MAN2B1", "PDXK", "ADAM10", "APAF1", "SLC27A2", "CEP290", "AGPAT2", "DEGS1", "SCAMP1", "SURF4", "SIGLEC5", "ARPC5", "PSMD3", "MGAM", "SYNGR1", "TYROBP", "GMFG", "TLR2", "DIAPH1", "TOM1", "FCGR3B", "LILRB3", "CPNE3", "DNAJC13", "PGLYRP1", "CREG1", "ATP6AP2", "IDH1", "DGAT1", "TRPM2", "STK10", "NFASC", "TMEM63A", "AP2A2", "STBD1", "NDUFC2", "ATG7", "BRI3", "VNN1", "RAB3D", "SNAP29", "CYB5R3", "PNP", "HP", "CFD", "PLAU", "SERPINA1", "SERPINA3", "C3", "CST3", "NRAS", "PIGR", "HLA-B", "", "LRG1", "ORM1", "AHSG", "TTR", "PPBP", "LTF", "FTL", "FTH1", "SLPI", "CAT", "FUCA1", "ALDOA", "CSTB", "A1BG", "KRT1", "CYBB", "ARG1", "ITGB2", "S100A8", "MPO", "GLA", "GSN", "S100A9", "PYGL", "GPI", "ITGAV", "CTSD", "ANXA2", "CAPN1", "TUBB", "PRSS2", "PSAP", "HEXB", "APRT", "CTSB", "HSP90AA1", "CD55", "GUSB", "HSP90AB1", "ELANE", "CTSG", "MME", "CD14", "PTPRC", "CD63", "ACAA1", "GSTP1", "CXCL1", "HMGB1", "CTSH", "FGR", "ALOX5", "LTA4H", "ALDOC", "HSPA1A", "HSPA1A", "RNASE2", "GAA", "HLA-C", "CTSA", "MGST1", "HSPA8", "SLC2A3", "ITGAM", "PYGB", "LAMP1", "EPX", "IGF2R", "IMPDH2", "FCGR2A", "RNASE3", "DEFA4", "XRCC6", "XRCC5", "LAMP2", "CYBA", "EEF2", "CEACAM1", "ALAD", "PRG2", "APEH", "CD59", "SELL", "MIF", "PKM", "MMP9", "JUP", "ANPEP", "ARSA", "B4GALT1", "ACP3", "GNS", "ARSB", "DSP", "TIMP2", "CD44", "GLB1", "PECAM1", "CD36", "HSPA6", "BPI", "PFKL", "GM2A", "CR1", "LGALS3", "PSMC3", "VCL", "PGAM1", "CD58", "CD53", "ORM2", "AOC1", "NFKB1", "TCN1", "CD33", "AZU1", "TNFRSF1B", "RAB3A", "RAB6A", "PSMB1", "ITGAL", "ITGAX", "IMPDH1", "AGA", "FPR1", "C5AR1", "NME2", "SLC2A5", "MMP8", "PTPRB", "FCAR", "PRTN3", "CXCR1", "CXCR2", "FPR2", "PTAFR", "CTSS", "PSMA2", "S100P", "PTX3", "STOM", "ATP6V0C", "CFP", "PSMA5", "MAPK1", "GCA", "GRN", "PTPN6", "SERPINB3", "PRDX6", "FCER1G", "HMOX2", "SERPINB1", "CHRNB4", "S100A7", "S100A11", "CEACAM8", "CDA", "GALNS", "CD68", "PRSS3", "SERPINB6", "AGL", "PSMC2")

gene_set2 <- c("CHI3L1", "PGM1", "SRP14", "DDOST", "CEACAM3", "CEACAM6", "MNDA", "FOLR3", "ACTR1B", "FRK", "PRCP", "ALDH3B1", "IQGAP1", "GYG1", "GLIPR1", "SERPINB10", "ADGRE5", "DNASE1L1", "SLC11A1", "PSEN1", "CAMP", "GDI2", "CCT8", "RAB5C", "RAB7A", "RAB27A", "P2RX1", "PSMD7", "HK3", "ACLY", "COPB1", "CTSC", "IST1", "CRISP3", "VCP", "NCKAP1L", "GSDMD", "DEFA1B", "SNAP25", "RAB4B", "RAB5B", "RAB10", "RAB14", "ACTR2", "RAP1B", "RAP2B", "RHOA", "LYZ", "B2M", "NPC2", "YPEL5", "RAP1A", "PPIA", "RAC1", "DYNLL1", "DYNLT1", "CSNK2B", "EEF1A1", "TUBB4B", "PAFAH1B2", "HBB", "SIRPA", "ADAM8", "CCT2", "OLR1", "LCN2", "S100A12", "RHOG", "TNFAIP6", "ATP11A", "AMPD3", "FABP5", "CAP1", "DSG1", "PLAUR", "PRKCD", "CKAP4", "DSC1", "CD47", "SIGLEC14", "BST1", "BST2", "ILF2", "IRAG2", "PTPRJ", "PRDX4", "PSMD2", "DNAJC3", "CHIT1", "KCNAB2", "PLD1", "PDAP1", "ROCK1", "TCIRG1", "ASAH1", "IQGAP2", "RAB31", "SPTAN1", "PKP1", "CDK13", "COTL1", "MLEC", "DYNC1H1", "FGL2", "MVP", "KPNB1", "PSMD6", "MAPK14", "C3AR1", "QPCT", "ANO6", "FLG2", "NHLRC3", "FRMPD3", "CLEC12A", "UBR4", "ATAD3B", "HGSNAT", "LAIR1", "LAMTOR1", "OLFM4", "NAPRT", "NBEAL2", "SLCO4C1", "UNC13D", "CYFIP1", "DOK3", "VPS35L", "TMC6", "C6orf120", "GOLGA7", "RAB44", "HUWE1", "TMEM179B", "ABCA13", "CAND1", "STING1", "TMED7-TICAM2", "ADGRG3", "HRNR", "ARMC8", "TBC1D10C", "SLC44A2", "MCEMP1", "OSCAR", "ORMDL3", "STK11IP", "GHDC", "LILRB2", "SLC15A4", "LILRA3", "CD177", "TXNDC5", "NFAM1", "LPCAT1", "TSPAN14", "SVIP", "MOSPD2", "PLEKHO2", "ATP8B4", "CLEC4C", "CANT1", "CLEC4D", "NCSTN", "DOCK2", "ARHGAP45", "GGH", "OSTF1", "PTPRN2", "ATP6V0A1", "RAB24", "TMBIM1", "RAB37", "ARL8A", "FAF2", "HVCN1", "PGM2", "MS4A3", "", "SERPINB12", "PSMB7", "CNN2", "PSMD1", "NEU1", "VAT1", "CPNE1", "CPPED1", "ARHGAP9", "ERP44", "CRACR2A", "FUCA2", "C1orf35", "VAMP8", "AP1M1", "ADGRE3", "CRISPLD2", "TOLLIP", "MAGT1", "CYSTM1", "DNAJC5", "DSN1", "PTGES2", "TMT1A", "RHOF", "RETN", "RAB18", "RAB9B", "MMP25", "CD93", "NIT2", "GPR84", "TMEM30A", "CMTM6", "CLEC5A", "ACTR10", "PLAC8", "ADA2", "CALML5", "COMMD9", "KCMF1", "VAPA", "COMMD3", "CTSZ", "BIN2", "CD300A", "LAMTOR3", "DPP7", "DBNL", "PYCARD", "PSMD13", "PPIE", "PA2G4", "HPSE", "ATP11B", "PADI2", "ATP8A1", "LAMTOR2", "PRG3", "DERA", "SIGLEC9", "CAB39", "RAP2C", "ATP6V1D", "TRAPPC1", "HEBP2", "DYNC1LI1", "ENPP4")

common_genes1 <- intersect(rownames(norm_expr), gene_set1)
common_genes2 <- intersect(rownames(norm_expr), gene_set2)

gene_set_sum1 <- colSums(norm_expr[common_genes1, ])
gene_set_sum2 <- colSums(norm_expr[common_genes2, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$Degranulation1 <- gene_set_sum1
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$Degranulation2 <- gene_set_sum2

colnames(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$Degranulation1)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_neutrophil_Degranulation <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = c("sample", "expression_group", "PDLIM2", "Degranulation1", "Degranulation2")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_neutrophil_Degranulation, file="merged_seurat_FilteredCells_Pdlim2Active_neutrophil_Degranulation-PDLIM2.csv")


# sum expression of Necroptosis genes in neutrophil
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_neutrophil, layer = "data")

# Reactome pathway necroptosis cell death molecules
gene_set <- c("ARHGAP10", "TLR4", "TNFRSF10A", "PSMD11", "PSMD12", "DFFA", "DNM1L", "PSMD14", "SDCBP", "APAF1", "TNFRSF10B", "PSMA7", "OGT", "TP73", "CFLAR", "SEPTIN4", "PSMD3", "DAPK3", "BCL2L11", "CHMP2A", "GAS2", "PRKN", "OPA1", "GSDME", "FLOT1", "DFFB", "IL1A", "IL1B", "LMNA", "", "TP53", "GSN", "H1-0", "HSP90AA1", "ELANE", "", "CD14", "VIM", "HMGB1", "UBB", "UBC", "GZMB", "H1-4", "BCL2", "MAPT", "IRF1", "CDH1", "IRF2", "DSP", "H1-5", "H1-3", "H1-2", "PSMC3", "PSMB1", "LMNB1", "APC", "FAS", "PSMA1", "PSMA2", "PSMA3", "PSMA4", "HMGB2", "YWHAQ", "MAPK3", "PSMA5", "PSMB4", "PSMB6", "PSMB5", "MAPK1", "CASP1", "NMT1", "AKT1", "AKT2", "YWHAB", "SFN", "DSG3", "CTNNB1", "ADD1", "PSMC2", "STAT3", "CASP3", "DCC", "PSMC4", "MAPK8", "FASLG", "PPP3CC", "PSMD8", "FNTA", "CASP4", "PSMB3", "PSMB2", "TNFSF10", "BCAP31", "PSMD7", "BMX", "CASP5", "KPNA1", "DAPK1", "CASP7", "CASP9", "CASP6", "BID", "GSDMD", "SEM1", "PSMA6", "YWHAG", "PSMC1", "PSMC5", "YWHAE", "PSMC6", "RPS27A", "UBA52", "PPP3R1", "YWHAZ", "DYNLL1", "UBE2L3", "XIAP", "CYCS", "E2F1", "SATB1", "DSG1", "H1-1", "PRKCQ", "YWHAH", "PTK2", "PRKCD", "C1QBP", "TJP1", "BAX", "BCL2L1", "TRAF2", "FADD", "PAK2", "PSMD2", "ROCK1", "BIRC3", "BIRC2", "RIPK1", "TP53BP2", "PMAIP1", "SPTAN1", "PKP1", "IL18", "DSG2", "TFDP1", "TFDP2", "FLOT2", "CASP8", "KPNB1", "PSMD6", "PLEC", "TRADD", "ADRM1", "CDC37", "BAK1", "OCLN", "UNC5A", "", "TMED7-TICAM2", "TICAM1", "UNC5B", "CDKN2A", "MLKL", "PDCD6IP", "CHMP7", "BAD", "CHMP4C", "OMA1", "PELI1", "DYNLL2", "CHMP6", "APIP", "ITCH", "PPP1R13B", "BMF", "", "PSMB7", "PSMD1", "BBC3", "CHMP4A", "UACA", "TP63", "CHMP4B", "CLSPN", "AVEN", "DIABLO", "STK26", "TJP2", "DAPK2", "DBNL", "APPL1", "ACIN1", "STUB1", "PSMD13", "CHMP2B", "AKT3", "CARD8", "CHMP3", "RIPK3", "MAGED1", "STK24", "LY96")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_neutrophil$Necroptosis <- gene_set_sum

colnames(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_neutrophil[[]]$Necroptosis)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_neutrophil_Necroptosis <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_neutrophil, vars = c("sample", "expression_group", "PDLIM2", "Necroptosis")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_neutrophil_Necroptosis, file="GSE243629_Pdlim2Active_Monocyte_Necroptosis-PDLIM2.csv")


# sum expression of reactome oxidative stress induced senescence molecules in monocyte
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, layer = "data")

# Reactome oxidative stress induced senescence molecules
gene_set <- c("AGO1", "AGO3", "AGO4", "BMI1", "CBX2", "CBX4", "CBX6", "CBX8", "CDK4", "CDK6", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "E2F1", "E2F2", "E2F3", "EED", "EZH2", "FOS", "H2AB1", "H2AC14", "H2AC19", "H2AC20", "H2AC4", "H2AC6", "H2AC7", "H2AJ", "H2AX", "H2AZ2", "H2BC1", "H2BC11", "H2BC12", "H2BC12L", "H2BC13", "H2BC14", "H2BC15", "H2BC17", "H2BC21", "H2BC26", "H2BC3", "H2BC5", "H2BC7", "H2BC9", "H3-3A", "H3C14", "H3C3", "H4C4", "IFNB1", "JUN", "KDM6B", "MAP2K3", "MAP2K4", "MAP2K6", "MAP2K7", "MAP3K5", "MAP4K4", "MAPK1", "MAPK10", "MAPK11", "MAPK14", "MAPK3", "MAPK8", "MAPK9", "MAPKAPK2", "MAPKAPK3", "MAPKAPK5", "MDM2", "MDM4", "MINK1", "MOV10", "PHC1", "PHC2", "PHC3", "RBBP4", "RBBP7", "RING1", "RNF2", "RPS27A", "SCMH1", "SUZ12", "TFDP1", "TFDP2", "TNIK", "TNRC6A", "TNRC6B", "TNRC6C", "TP53", "TXN", "UBA52", "UBB", "UBC")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$OxidativeSenescence <- gene_set_sum

colnames(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$OxidativeSenescence)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_Monocyte_OxidativeSenescence <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_Monocyte, vars = c("sample", "expression_group", "PDLIM2", "OxidativeSenescence")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_Monocyte_OxidativeSenescence, file="GSE243629_Pdlim2Active_Monocyte_OxidativeSenescence-PDLIM2.csv")


# sum expression of reactome SASP signaling molecules
# Extract the normalized expression values
norm_expr <- GetAssayData(merged_seurat_FilteredCells_Pdlim2Active_Monocyte, layer = "data")

# reactome SASP molecules
gene_set <- c("ANAPC1", "ANAPC10", "ANAPC11", "ANAPC15", "ANAPC16", "ANAPC2", "ANAPC4", "ANAPC5", "ANAPC7", "CCNA1", "CCNA2", "CDC16", "CDC23", "CDC26", "CDC27", "CDK2", "CDK4", "CDK6", "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "CEBPB", "CXCL8", "EHMT1", "EHMT2", "FOS", "FZR1", "H2AB1", "H2AC14", "H2AC19", "H2AC20", "H2AC4", "H2AC6", "H2AC7", "H2AJ", "H2AX", "H2AZ2", "H2BC1", "H2BC11", "H2BC12", "H2BC12L", "H2BC13", "H2BC14", "H2BC15", "H2BC17", "H2BC21", "H2BC26", "H2BC3", "H2BC5", "H2BC8", "H2BC9", "H3-3A", "H3C1", "H3C14", "H4C11", "IGFBP7", "IL1A", "IL6", "JUN", "MAPK1", "MAPK3", "MAPK7", "NFKB1", "RELA", "RPS27A", "RPS6KA1", "RPS6KA2", "RPS6KA3", "STAT3", "UBA52", "UBB", "UBC", "UBE2C", "UBE2D1", "UBE2E1", "UBE2S", "VENTX")

common_genes <- intersect(rownames(norm_expr), gene_set)

gene_set_sum <- colSums(norm_expr[common_genes, ])

# add them to metadata
merged_seurat_FilteredCells_Pdlim2Active_Monocyte$SASP <- gene_set_sum

colnames(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]])
head(merged_seurat_FilteredCells_Pdlim2Active_Monocyte[[]]$SASP)

# fetch data
merged_seurat_FilteredCells_Pdlim2Active_Monocyte_SASP <- FetchData(object = merged_seurat_FilteredCells_Pdlim2Active_Monocyte, vars = c("sample", "expression_group", "PDLIM2", "SASP")) 
write.csv (merged_seurat_FilteredCells_Pdlim2Active_Monocyte_SASP, file="GSE243629_Pdlim2Active_Monocyte_SASP-PDLIM2.csv")

load("GSE243629.RData")
save (list=ls(), file="GSE243629.RData")
