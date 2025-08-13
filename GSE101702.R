# Version info: R 4.2.2,Biobase 2.58.0,GEOquery 2.66.0,limma 3.54.0
################################################################
# install GEOquery and limma
if (!require("BiocManager",quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")
BiocManager::install("limma")

# install umap
install.packages("umap")

#   Differential expression analysis with limma
library(GEOquery)
library(limma)
library(umap)

# load saved RData
load("GSE101702.RData")

# PDLIM2 Agilent probe names: A_33_P3311371 (used for analysis),A_24_P160380 

# load series and platform data from GEO

gset <- getGEO("GSE101702",GSEMatrix =TRUE,AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL21185",attr(gset,"names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group membership for all samples
gsms <- paste0("11111111111111111111111111111000000000000000000011",
               "11111111111111111111111111111111111111000000000000",
               "00000000000000000000011111111111111111111111111111",
               "111111111")
sml <- strsplit(gsms,split="")[[1]]

# log2 transformation
ex <- exprs(gset)
qx <- as.numeric(quantile(ex,c(0.,0.25,0.5,0.75,0.99,1.0),na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

ex[1:5,1:5]
exprs(gset)[1:5,1:5]

# assign samples to groups and set up design matrix
gs <- factor(sml)
groups <- make.names(c("healthy","Influenza"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0,gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)),] # skip missing values

fit <- lmFit(gset,design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts,levels=design)
fit2 <- contrasts.fit(fit,cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2,0.01)
tT <- topTable(fit2,adjust="fdr",sort.by="B",number=250)

tT <- subset(tT,select=c("ID","adj.P.Val","P.Value","t","B","logFC","GB_ACC","SEQUENCE","SPOT_ID"))
write.table(tT,file=stdout(),row.names=F,sep="\t")

# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.
tT2 <- topTable(fit2,adjust="fdr",sort.by="B",number=Inf)
hist(tT2$adj.P.Val,col = "grey",border = "white",xlab = "P-adj",
     ylab = "Number of genes",main = "P-adj value distribution")

# summarize test results as "up","down" or "not expressed"
dT <- decideTests(fit2,adjust.method="fdr",p.value=0.05,lfc=0)

# Venn diagram of results
vennDiagram(dT,circle.col=palette())

# create Q-Q plot for t-statistic
t.good <- which(!is.na(fit2$F)) # filter out bad probes
qqt(fit2$t[t.good],fit2$df.total[t.good],main="Moderated t statistic")

# volcano plot (log P-value vs log fold change)
colnames(fit2) # list contrast names
ct <- 1        # choose contrast of interest
# Please note that the code provided to generate graphs serves as a guidance to
# the users. It does not replicate the exact GEO2R web display due to multitude
# of graphical options.
# 
# The following will produce basic volcano plot using limma function:
volcanoplot(fit2,coef=ct,main=colnames(fit2)[ct],pch=20,
            highlight=length(which(dT[,ct]!=0)),names=rep('+',nrow(fit2)))

# MD plot (log fold change vs mean log expression)
# highlight statistically significant (p-adj < 0.05) probes
plotMD(fit2,column=ct,status=dT[,ct],legend=F,pch=20,cex=1)
abline(h=0)

################################################################
# General expression data analysis
ex <- exprs(gset)

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6,height=5)
ord <- order(gs)  # order samples by group
palette(c("#1B9E77","#7570B3","#E7298A","#E6AB02","#D95F02",
          "#66A61E","#A6761D","#B32424","#B324B3","#666666"))
par(mar=c(7,4,2,1))
title <- paste ("GSE101702","/",annotation(gset),sep ="")
boxplot(ex[,ord],boxwex=0.6,notch=T,main=title,outline=FALSE,las=2,col=gs[ord])
legend("topleft",groups,fill=palette(),bty="n")
dev.off()

# expression value distribution
par(mar=c(4,4,2,1))
title <- paste ("GSE101702","/",annotation(gset)," value distribution",sep ="")
plotDensities(ex,group=gs,main=title,legend ="topright")

# UMAP plot (dimensionality reduction)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex),]  # remove duplicates
ump <- umap(t(ex),n_neighbors = 15,random_state = 123)
par(mar=c(3,3,2,6),xpd=TRUE)
plot(ump$layout,main="UMAP plot,nbrs=15",xlab="",ylab="",col=gs,pch=20,cex=1.5)
legend("topright",inset=c(-0.15,0),legend=levels(gs),pch=20,
       col=1:nlevels(gs),title="Group",pt.cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout,labels = rownames(ump$layout),method="SANN",cex=0.6)

# mean-variance trend,helps to see if precision weights are needed
plotSA(fit2,main="Mean variance trend,GSE101702")

# save table and data
write.table(tT,file="GSE101702_DE_Healthy-Influ_top250.txt",row.names=F,sep="\t")
write.table(tT,file="GSE101702_DE_Healthy-Influ_top250.csv",row.names=F,sep=",")
save(list=ls(),file="GSE101702.RData")

# Specify the target gene
PDLIM2 <- "A_33_P3311371"

# Extract expression levels of the target gene
PDLIM2_expr <- ex[PDLIM2,]

# Calculate quantiles (33% and 66%)
low_cutoff <- quantile(PDLIM2_expr,0.33)
high_cutoff <- quantile(PDLIM2_expr,0.66)

# Classify samples
PDLIM2_groups <- ifelse(
  PDLIM2_expr <= low_cutoff,"Low",
  ifelse(PDLIM2_expr >= high_cutoff,"High","Middle")
)

# Output the grouping
print(PDLIM2_groups)

# add PDLIM2_groups to gset
gset@phenoData@data$PDLIM2_groups <- PDLIM2_groups
colnames(gset@phenoData@data)

# fetch PDLIM2_groups data
PDLIM2_groups_fetch <- data.frame(Expression = exprs(gset)["A_33_P3311371",],
           Group = pData(gset)$PDLIM2_groups)
write.table(PDLIM2_groups_fetch,file="GSE101702_PDLIM2_groups_fetch.csv",row.names=F,sep=",")

# define fetch_data function to Fetch gene expression for specific genes and metadata columns
fetch_data <- function(gset,genes,meta_columns) {
  # Extract expression matrix
  expression_data <- exprs(gset)[genes,,drop = FALSE]
  
  # Extract metadata
  metadata <- pData(gset)[,meta_columns,drop = FALSE]
  
  # Combine expression and metadata
  data.frame(t(expression_data),metadata)
}

# fetch PDLIM2 gene expression and pdata 
PDLIM2_groups_fetch02<-fetch_data(gset,genes = c("A_33_P3311371"),meta_columns = c("title","PDLIM2_groups","group","diagnosis:ch1","severity:ch1"))
write.table(PDLIM2_groups_fetch02,file="GSE101702_PDLIM2_groups_fetch02.csv",row.names=F,sep=",")


## DEG from  PDLIM2_Low vs PDLIM2_High
library(limma)

# Extract metadata
meta_data <- pData(gset)

# Subset for Low and High groups
subset_samples <- meta_data$PDLIM2_groups %in% c("Low","High")
exprSet_subset <- exprs(gset)[,subset_samples]

# Design matrix
group_factor <- factor(meta_data$PDLIM2_groups[subset_samples],levels = c("Low","High"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("Low","High")

# Fit the linear model
fit <- lmFit(exprSet_subset,design)

# Contrast for Low vs High
contrast.matrix <- makeContrasts(Low_vs_High = Low - High,levels = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEG results
DEG_PDLIM2_LowHigh <- topTable(fit2,number= Inf,coef = "Low_vs_High",adjust = "fdr")
write.table(DEG_PDLIM2_LowHigh,file="GSE101702_DEG_PDLIM2_LowHigh.csv",sep=",")


# Subset for svre and hlty groups
subset_samples <- meta_data$'severity:ch1' %in% c("flu_svre","hlty_ctrl")
exprSet_subset <- exprs(gset)[,subset_samples]

# Design matrix
group_factor <- factor(meta_data$'severity:ch1'[subset_samples],levels = c("flu_svre","hlty_ctrl"))
design <- model.matrix(~ 0 + group_factor)
colnames(design) <- c("flu_svre","hlty_ctrl")

# Fit the linear model
fit <- lmFit(exprSet_subset,design)

# Contrast for Low vs High
contrast.matrix <- makeContrasts(flu_svre_vs_hlty_ctrl = flu_svre - hlty_ctrl,levels = design)
fit2 <- contrasts.fit(fit,contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DEG results
DEG_Disease_SvreHlty <- topTable(fit2,number= Inf,coef = "flu_svre_vs_hlty_ctrl",adjust = "fdr")
write.table(DEG_Disease_SvreHlty,file="GSE101702_DEG_Disease_SvreHlty.csv",sep=",")

# gset@phenoData@data$title(INFG_ID_002)/characteristics_ch1(diagnosis:influenza)/characteristics_ch1.4(severity:flu_mod)/diagnosis:ch1(influenza)/severity:ch1(flu_mod)/group(factor w/ 2 levels "healthy" "Influenza") - sample and meta data

# gset@assayData$exprs - gene expression matrix data with Agilent probes as row and samples as column

# gset@featureData@data$ID (agilent probe names)/GENE_SYMBOL/GENE_NMAE/ENSEMBL_ID - gene ID

###calculate sum expression of all NF-kB targets
NFkB_list1 <- c("A_23_P407096","A_22_P00018083","A_33_P3315600","A_33_P3246990","A_33_P3231187","A_33_P3270863","A_23_P116280","A_23_P162322","A_33_P3394352","A_24_P147398","A_24_P279489","A_33_P3297030","A_23_P161190","A_23_P167096","A_23_P34345","A_23_P147805","A_23_P351275","A_24_P200219","A_23_P47704","A_23_P141847","A_33_P3414591","A_23_P71073","A_24_P28977","A_33_P3319905","A_23_P19333","A_21_P0000178","A_33_P3803639","A_33_P3424295","A_24_P89891","A_33_P3315764","A_23_P26810","A_23_P386478","A_23_P30435","A_23_P94754","A_23_P14174","A_21_P0000113","A_23_P121253","A_33_P3226810","A_23_P51936","A_33_P3286157","A_24_P54174","A_33_P3322274","A_24_P157926","A_23_P421423","A_23_P376488","A_23_P157865","A_33_P3380807","A_23_P92499","A_23_P352619","A_23_P376096","A_23_P90311","A_33_P3227375","A_33_P3365735","A_33_P3415820","A_24_P142118","A_24_P923251","A_32_P86763","A_23_P65618","A_23_P393620","A_33_P3334305","A_23_P393099","NM_003226","A_32_P184394","A_33_P3224710","A_23_P212508","A_23_P110851","A_23_P365719","A_23_P360291","A_23_P259580","A_33_P3325914","A_23_P59005","A_24_P148590","A_23_P63847","A_23_P207367","A_33_P3380751","A_23_P354705","A_33_P3382276","A_24_P388528","A_23_P7313","A_33_P3489646","A_23_P430718","A_23_P26847","A_23_P134176","A_23_P154840","A_23_P131846","A_24_P46093","A_24_P359191","A_23_P75811","A_33_P3282434","A_23_P139820","A_24_P381494","A_23_P156310","A_23_P148297","A_23_P50919","A_23_P214330","A_23_P2920","A_33_P3399788","A_23_P162918","A_23_P44663","A_33_P3289659","A_23_P218111","A_32_P74712","A_33_P3339100","A_23_P137697","A_33_P3421053","A_23_P97112","A_23_P109034","A_23_P128323","A_33_P3371727","A_23_P137016","A_33_P3371718","A_33_P3408913","A_33_P3408918","A_24_P335092","A_23_P201711","A_23_P94800","A_23_P137984","A_23_P252106","A_23_P214139","A_23_P55706","A_23_P56938","A_33_P3232562","A_33_P3338559","A_23_P326691","A_24_P468810","A_21_P0000601","A_23_P64525","A_33_P3260782","A_23_P360744","A_23_P26629","A_23_P121064","A_23_P127579","A_23_P18493","A_23_P338890","A_23_P2271","A_24_P250922","A_24_P48723","A_24_P403417","A_33_P3219811","A_33_P3298159","A_23_P98085","A_24_P913115","A_33_P3393537","A_24_P102821","A_23_P65427","A_23_P151614","A_23_P111000","A_32_P65616","A_23_P144054","A_24_P399630","A_23_P1473","A_23_P153583","A_23_P213959","A_33_P3816688","A_33_P3345036","A_23_P51646","A_23_P80739","A_33_P3306146","A_23_P24104","A_23_P345118","A_23_P92057","A_21_P0012238","A_24_P844984","A_23_P131653","A_23_P138938","A_32_P49199","A_23_P208747","A_23_P125829","A_23_P417918","A_24_P279870","A_24_P339944","A_33_P3301955","A_24_P360529","A_33_P3234809","A_24_P71153","A_23_P137765","A_33_P3413741","A_23_P169494","A_24_P397613","A_33_P3272075","A_23_P419156","A_24_P124624","A_33_P3338166","A_33_P3284345","A_23_P315815","A_33_P3299066","A_23_P131208","A_33_P3216297","A_33_P3767927","A_21_P0000153","A_23_P206661","A_23_P69699","A_23_P217280","A_23_P316410","A_33_P3286552","A_23_P204791","A_33_P3365167","A_23_P420863","A_24_P213161","A_23_P212089","A_23_P30655","A_23_P106002","A_23_P202156","A_23_P30024","A_23_P1320","A_21_P0012447","A_24_P319923","A_23_P143817","A_33_P3245163","A_23_P215956","A_23_P31073","A_33_P3311795","A_23_P17663","A_33_P3412384","A_24_P84657","A_33_P3401571","A_23_P256784","A_23_P404045","A_23_P400078","A_23_P129629","A_23_P40174","A_23_P161698","A_23_P1691","A_23_P116235","A_33_P3227990")
                
NFkB_list2 <- c("A_23_P208085","A_24_P402080","A_23_P5002","A_23_P153616","A_23_P166848","A_33_P3248265","A_23_P156683","A_24_P125469","A_33_P3306163","A_23_P128919","A_24_P20630","A_23_P169437","A_23_P143178","A_33_P3223780","A_33_P3338116","A_23_P76249","A_23_P218047","A_33_P3287611","A_23_P27133","A_33_P3215838","A_33_P3215843","A_23_P168828","A_24_P133253","A_23_P124892","A_23_P500353","A_23_P319423","A_24_P241815","A_24_P378019","A_23_P214360","A_23_P125082","A_23_P41765","A_23_P122924","A_24_P320410","A_33_P3407013","A_23_P71037","A_24_P230563","A_23_P127288","A_33_P3284933","A_23_P76078","A_23_P30122","A_33_P3246829","A_33_P3246833","A_23_P79518","A_23_P72096","A_23_P29953","A_23_P251031","A_23_P7560","A_23_P91943","A_33_P3243887","A_23_P126735","A_23_P119943","A_23_P42868","A_23_P151294","A_23_P71774","A_23_P45871","A_23_P42257","A_33_P3423551","A_32_P36235","A_23_P112026","A_23_P371215","A_23_P162874","A_32_P199252","A_33_P3665777","A_23_P81973","A_23_P14986","A_21_P0000064","A_23_P256107","A_33_P3416231","A_22_P00023586","A_22_P00007805","A_23_P500998","A_23_P120883","A_32_P22257","A_24_P409857","A_24_P311926","A_33_P3399208","A_33_P3424800","A_33_P3379947","A_24_P56388","A_33_P3231277","A_33_P3276718","A_33_P3276713","A_33_P3252534","A_33_P3421827","A_23_P47665","A_23_P27400","A_33_P3295203","A_33_P3295200","A_33_P3336720","A_23_P117602","A_22_P00007487","A_23_P350001","A_23_P202658","A_33_P3358943","A_32_P169114","A_24_P381844","A_33_P3318627","A_33_P3334895","A_23_P257962","A_33_P3348127","A_33_P3418790","A_24_P323072","A_33_P3585268","A_24_P404840","A_33_P3282489","A_33_P3265606","A_23_P145114","A_23_P62890")

NFkB_list3 <- c("A_33_P3360341","A_24_P239606","A_23_P374689","A_33_P3365142","A_23_P34093","A_33_P3253249","A_32_P342064","A_24_P58337","A_33_P3318796","A_23_P106194","A_24_P334130","A_23_P46829","A_23_P55936","A_23_P164773","A_23_P369815","A_23_P63896","A_33_P3381127","A_33_P3332112","A_23_P43846","A_24_P314451","A_33_P3226832","A_24_P319364","A_33_P3238250","A_33_P3292596","A_23_P89249","A_23_P145669","A_23_P157333","A_24_P236091","A_23_P83328","A_23_P104188","A_23_P214080","A_33_P3351944","A_23_P215790","A_33_P3351955","A_33_P3258392","A_23_P214821","A_23_P119478","A_23_P385034","A_33_P3418516","A_23_P110712","A_23_P135548","A_24_P410610","A_23_P3643","A_23_P133153","A_24_P341674","A_23_P48740","A_33_P3325700","A_23_P251945","A_23_P169092","A_24_P394940","A_23_P36397","A_32_P86289","A_24_P920646","A_33_P3351371","A_23_P37410","A_23_P18452","A_23_P155755","A_24_P277367","A_24_P183150","A_23_P125278","A_33_P3343175","A_24_P303091","A_23_P7144","A_33_P3330264","A_33_P3287631","A_23_P501754","A_33_P3222424","A_23_P133408","A_33_P3354945","A_33_P3354935","A_33_P3354940","A_33_P3317670","A_24_P342484","A_23_P423389","A_33_P3303697","A_33_P3281985","A_24_P277934","A_33_P3249046","A_23_P376704","A_23_P137665","A_24_P120115","A_23_P156807","A_23_P209394","A_22_P00003932","A_23_P156687","A_33_P3253804","A_22_P00014880","A_33_P3289848","A_23_P58788","A_24_P89457","A_24_P166663","A_24_P131589","A_23_P70670","A_24_P320033","A_32_P175934","A_21_P0000152","A_33_P3294509","A_33_P3250680","A_23_P57036","A_33_P3245439","A_23_P98410","A_23_P167328","A_23_P338479","A_33_P3381513","A_33_P3272493","A_24_P186539","A_24_P305345","A_23_P343398","A_23_P412321","A_23_P361773","A_24_P278747","A_23_P202837","A_23_P152838","A_23_P207564","A_33_P3354607","A_33_P3316273","A_23_P503072","A_24_P133905","A_24_P313418","A_23_P17065","A_23_P89431","A_33_P3377151","A_23_P123853","A_23_P26325","A_33_P3324004","A_23_P49759","A_23_P134454","A_24_P12626","A_23_P35912","A_23_P97541","A_21_P0014047","A_23_P101407","A_33_P3347869","A_33_P3277514","A_23_P137139","A_23_P99452","A_33_P3419785","A_23_P54144","A_33_P3237150","A_24_P303989","A_23_P314115","A_24_P64344","A_33_P3363637","A_23_P34126","A_23_P127891","A_32_P7316","A_23_P128744","A_23_P4662","A_24_P122921","A_33_P3251932","A_33_P3398526","A_23_P210886","A_23_P321703","A_23_P152002","A_33_P3355185","A_23_P352266","A_23_P208706","A_33_P3341144","A_33_P3240941","A_33_P3240936","A_23_P52806","A_23_P37441","A_33_P3359771","A_23_P148879","A_33_P3234580","A_23_P31921","A_24_P18105","A_23_P216094","A_24_P295245","A_23_P87363","A_33_P3279720","A_24_P159948","A_33_P3242748","A_33_P3242743","A_23_P113111","A_24_P202522","A_33_P3420757","A_33_P3223592","A_23_P136777","A_23_P203183","A_23_P70387","A_23_P216023","A_23_P78944","A_33_P3313245","A_33_P3259736","A_33_P3356007","A_23_P104464","A_23_P152906","A_33_P3365117","A_24_P152968","A_23_P257971","A_23_P36641","A_33_P3666346","A_24_P28657","A_23_P115261","A_33_P3396008","A_23_P93360","A_33_P3396010","A_33_P3421243","A_23_P378926","A_33_P3214035","A_24_P237270","A_23_P74299","A_24_P291658","A_23_P374082","A_23_P81369","A_33_P3361112","A_23_P119763","A_23_P100539","A_33_P3385266","A_33_P3298043","A_23_P15272","A_23_P128094","A_33_P3376947","A_33_P3336760","A_33_P3332414","A_23_P82523","A_33_P3422897","A_24_P235429")

# Make sure gene_list is present in the matrix
genes_present <- intersect(NFkB_list3,rownames(ex))

# Sum expression of the selected genes for each sample (column)
gene_sum <- colSums(ex[genes_present,,drop = FALSE])

# add NFkB_sum to gset
gset@phenoData@data$NFkB_sum3 <- gene_sum
colnames(gset@phenoData@data)

# fetch NFkB_sum and other metadata together with PDLIM2 expression 
NFkB_sum<-fetch_data(gset,genes = c("A_33_P3311371"),meta_columns = c("title","PDLIM2_groups","group","diagnosis:ch1","severity:ch1","NFkB_sum1","NFkB_sum2","NFkB_sum3"))
write.table(NFkB_sum,file="NFkB_sum.csv",row.names=F,sep=",")

#####calculate sum expression of top upregulated NF-kB targets from flu vs healthy
NFkB_list4 <- c("A_23_P45871", "A_23_P169437", "A_23_P166848", "A_23_P17663", "A_23_P40174", "A_24_P378019", "A_23_P167328", "A_23_P351275", "A_23_P152002", "A_33_P3284933", "A_33_P3246833", "A_23_P338479", "A_24_P295245", "A_33_P3226810", "A_24_P235429", "A_23_P14174", "A_23_P147805", "A_23_P208747", "A_23_P62890", "A_33_P3246829", "A_33_P3381513")

#Make sure gene_list is present in the matrix
genes_present <- intersect(NFkB_list4,rownames(ex))

# Sum expression of the selected genes for each sample (column)
gene_sum <- colSums(ex[genes_present,,drop = FALSE])

# add NFkB_sum to gset
gset@phenoData@data$NFkB_sum4 <- gene_sum
colnames(gset@phenoData@data)

# fetch NFkB_sum and other metadata together with PDLIM2 expression 
NFkB_sum<-fetch_data(gset,genes = c("A_33_P3311371"),meta_columns = c("title","PDLIM2_groups","group","diagnosis:ch1","severity:ch1","NFkB_sum1","NFkB_sum2","NFkB_sum3","NFkB_sum4"))
write.table(NFkB_sum,file="NFkB_sum4.csv",row.names=F,sep=",")

#####calculate sum expression of STAT all genes
STAT_list1 <- c("A_24_P12401", "A_23_P70398", "A_23_P141917", "A_23_P390518", "A_23_P376488", "A_23_P200780", "A_23_P211957", "A_33_P3331451", "A_24_P402438", "A_24_P79054", "A_33_P3408203", "A_23_P377291", "A_33_P3214650", "A_23_P368067", "A_23_P59005", "A_23_P47879", "A_33_P3380867", "A_23_P207367", "A_23_P305198", "A_23_P68031", "A_23_P100795", "A_24_P116805", "A_33_P3213064", "A_24_P274270", "A_23_P308603", "A_23_P117546", "A_23_P343808", "A_33_P3356406", "A_24_P367211", "A_23_P207981", "A_24_P328492", "A_33_P3349334", "A_23_P207058", "A_23_P128215", "A_23_P420196", "A_24_P68585", "A_23_P200443", "A_33_P3327245", "A_23_P341532", "A_33_P3297245", "A_33_P3234317", "A_23_P39076", "A_23_P55706", "A_33_P3209433", "A_23_P104689", "A_23_P56938", "A_24_P357100", "A_23_P118392", "A_33_P3303121", "A_24_P81965", "A_23_P2661", "A_33_P3211604", "A_33_P3342285", "A_33_P3356502", "A_24_P192262", "A_23_P40952", "A_23_P215406", "A_23_P162486", "A_23_P309701", "A_32_P1445", "A_23_P99027", "A_23_P105436", "A_23_P338890", "A_23_P250629", "A_33_P3350726", "A_23_P252062", "A_23_P345118", "A_24_P366509", "A_23_P66543", "A_23_P132526", "A_33_P3593774", "A_33_P3367855", "A_24_P29401", "A_33_P3304170", "A_24_P71244", "A_33_P3353941", "A_23_P346969", "A_23_P92057", "A_21_P0012238", "A_23_P164536", "A_23_P76322", "A_33_P3419419", "A_23_P200710", "A_33_P3231878", "A_33_P3260062", "A_33_P3211174", "A_23_P149678", "A_33_P3349269")

STAT_list2 <- c("A_23_P163353", "A_23_P421401", "A_22_P00011688", "A_23_P300033", "A_24_P124349", "A_23_P58396", "A_24_P339944", "A_23_P113701", "A_33_P3372666", "A_23_P205900", "A_33_P3322814", "A_32_P141969", "A_33_P3260134", "A_33_P3322784", "A_33_P3322804", "A_24_P343559", "A_24_P265506", "A_23_P63190", "A_33_P3317825", "A_23_P502464", "A_23_P389897", "A_23_P202156", "A_23_P30024", "A_23_P119627", "A_33_P3245163", "A_23_P215956", "A_23_P34606", "A_23_P12746", "A_24_P88850", "A_23_P40174", "A_33_P3400192", "A_23_P167692", "A_23_P356152", "A_33_P3250289", "A_24_P286898", "A_23_P37910", "A_24_P283288", "A_24_P406132", "A_23_P145376", "A_33_P3250055", "A_33_P3763846", "A_33_P3409513", "A_33_P3409508", "A_23_P45025", "A_33_P3254216", "A_23_P257895", "A_33_P3379571", "A_33_P3323435", "A_24_P265265", "A_33_P3267482", "A_23_P318300", "A_23_P366394", "A_33_P3284557", "A_23_P105409", "A_33_P3301851", "A_24_P284523", "A_33_P3653330", "A_32_P6344", "A_23_P208835", "A_23_P20248", "A_24_P233488", "A_24_P122137", "A_23_P306507", "A_33_P7465707", "A_33_P3317815", "A_24_P71973", "A_33_P3323298", "A_24_P59667", "A_23_P329112", "A_33_P3878772", "A_33_P3784283", "A_23_P41765", "A_33_P3335725", "A_23_P34066", "A_24_P320410", "A_23_P404494", "A_32_P223777", "A_33_P3233834", "A_33_P3287338", "A_33_P3233843", "A_33_P3233841", "A_23_P502470", "A_32_P45168", "A_21_P0000171", "A_33_P3288844", "A_33_P3407013", "A_23_P71037", "A_33_P3328254", "A_23_P129556", "A_33_P3395947", "A_33_P3349045", "A_23_P213706", "A_32_P217750", "A_33_P3691168", "A_33_P3290394", "A_23_P148473", "A_24_P203000", "A_24_P230563", "A_23_P127288", "A_23_P27606", "A_23_P168288", "A_33_P3422124", "A_23_P62607", "A_24_P227927", "A_23_P91850", "A_33_P3421118", "A_23_P145514", "A_23_P56604", "A_33_P3211608", "A_23_P51126", "A_24_P63019", "A_24_P200023", "A_33_P3396389", "A_23_P79518", "A_23_P72096", "A_33_P3221960", "A_23_P28334", "A_24_P208567", "A_33_P3251876", "A_33_P3211666", "A_23_P500206", "A_32_P188860", "A_23_P166775", "A_24_P157370", "A_33_P3540143", "A_33_P3224809", "A_23_P17706", "A_33_P3399268", "A_33_P3399267", "A_33_P3399263", "A_23_P138680", "A_23_P85209", "A_24_P280113", "A_24_P288685", "A_23_P72077", "A_23_P399156", "A_23_P7560", "A_23_P71867", "A_24_P322741", "A_23_P203173", "A_23_P126735", "A_23_P334021", "A_23_P417282", "A_24_P304423", "A_23_P13907", "A_33_P3268612", "A_24_P19677", "A_23_P74112", "A_23_P151294", "A_33_P3277988", "A_23_P112026", "A_23_P153320", "A_23_P98183", "A_32_P87697", "A_33_P3276718", "A_33_P3276713", "A_24_P407717", "A_33_P3343155", "A_24_P72064", "A_33_P3260733", "A_33_P3360341", "A_23_P159191", "A_23_P106194", "A_33_P3409269", "A_23_P356070", "A_33_P3215883", "A_21_P0000024", "A_23_P92754", "A_23_P212830", "A_23_P500501", "A_24_P206624", "A_23_P303145", "A_23_P202334", "A_23_P316447", "A_23_P372923", "A_23_P301304", "A_33_P3379886", "A_33_P3419696", "A_21_P0010727", "A_23_P63390", "A_24_P314585", "A_33_P3351944", "A_23_P215790", "A_33_P3351955", "A_23_P155979", "A_33_P3214550", "A_23_P67932", "A_23_P18452", "A_23_P125278", "A_33_P3343175", "A_24_P303091", "A_23_P120899", "A_33_P3341442", "A_23_P133408", "A_23_P144096", "A_33_P3411333", "A_33_P3388501", "A_23_P126278", "A_23_P411296", "A_24_P89457", "A_24_P397107", "A_23_P152838", "A_24_P313418", "A_23_P89431", "A_23_P26325", "A_23_P162010", "A_23_P39465", "A_33_P3220911", "A_24_P753161", "A_24_P63380", "A_23_P1431", "A_24_P527404", "A_19_P00805548", "A_33_P3219256", "A_23_P19624", "A_23_P210886", "A_33_P3355185", "A_23_P352266", "A_33_P3352382", "A_33_P3319967", "A_23_P160354", "A_24_P110983", "A_23_P208870", "A_33_P3400171", "A_33_P3381647", "A_33_P3275235")

#Make sure gene_list is present in the matrix
genes_present <- intersect(STAT_list2,rownames(ex))

# Sum expression of the selected genes for each sample (column)
gene_sum <- colSums(ex[genes_present,,drop = FALSE])

# add NFkB_sum to gset
gset@phenoData@data$STAT_sum2 <- gene_sum
colnames(gset@phenoData@data)

# fetch NFkB_sum and other metadata together with PDLIM2 expression 
STAT_sum<-fetch_data(gset,genes = c("A_33_P3311371"),meta_columns = c("title","PDLIM2_groups","group","diagnosis:ch1","severity:ch1","NFkB_sum1","NFkB_sum2","NFkB_sum3","NFkB_sum4", "STAT_sum1", "STAT_sum2"))
write.table(STAT_sum,file="STAT_sum.csv",row.names=F,sep=",")

#####calculate sum expression of STAT all genes with DEG (FC>=1,padj<=0.05)
STAT_list3 <- c("A_33_P3319967", "A_33_P3352382", "A_33_P3211666", "A_23_P63390", "A_23_P40174", "A_33_P3251876", "A_24_P63019", "A_24_P208567", "A_21_P0010727", "A_23_P207058", "A_24_P283288", "A_23_P28334", "A_33_P3221960", "A_23_P126278", "A_23_P68031", "A_32_P87697", "A_24_P203000")

# Make sure gene_list is present in the matrix
genes_present <- intersect(STAT_list3,rownames(ex))

# Sum expression of the selected genes for each sample (column)
gene_sum <- colSums(ex[genes_present,,drop = FALSE])

# add NFkB_sum to gset
gset@phenoData@data$STAT_sum3 <- gene_sum
colnames(gset@phenoData@data)

# fetch NFkB_sum and other metadata together with PDLIM2 expression 
STAT_sum<-fetch_data(gset,genes = c("A_33_P3311371"),meta_columns = c("title","PDLIM2_groups","group","diagnosis:ch1","severity:ch1","NFkB_sum1","NFkB_sum2","NFkB_sum3","NFkB_sum4", "STAT_sum1", "STAT_sum2", "STAT_sum3"))
write.table(STAT_sum,file="STAT_sum.csv",row.names=F,sep=",")
