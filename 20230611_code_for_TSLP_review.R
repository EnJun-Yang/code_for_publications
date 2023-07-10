library(Seurat)
library(ggplot2)
library(tidyverse)
library(scRepertoire)
setwd("/Users/enjunyang/Documents/Personal/Taku/20220312_Analysis/")
Cont_Ear.vdj <- read.csv("../data/Cont_Ear_outs/vdj_t/filtered_contig_annotations.csv")
TSLP_Ear.vdj <- read.csv("../data/TSLP_Ear_outs/vdj_t/filtered_contig_annotations.csv")
contig_list <- list(TSLP_Ear.vdj, Cont_Ear.vdj)
combined <- combineTCR(contig_list, samples = c("TSLP_Ear", "Cont_Ear"), ID=c("TSLP", "Control"), cells="T-AB")
# combined1 <- combineTCR(contig_list, samples = c("TSLP_Ear", "Cont_Ear"), ID=c("TSLP", "Control"), cells="T-AB", removeNA=TRUE, removeMulti=TRUE)
### To reduce cell loss, just go with all cells labelled as alpha-betas
dpi=300
png(filename="TSLP_vs_Cont_Clonotype_abundance.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
abundanceContig(combined, cloneCall = "gene", scale = F)
dev.off()
png(filename="TSLP_vs_Cont_Clonotype_fraction.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
clonalHomeostasis(combined, cloneCall="gene", cloneTypes=c(Rare = 1e-04, Small = 0.001,  Medium = 0.01, Large = 0.1,  Hyperexpanded = 1))
dev.off()
### Clonotype expansion by gene usage indicates the Control has more clonality than the TSLP ear
png(filename="TSLP_vs_Cont_TRB_V_usage.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
vizGenes(combined, gene = "V", chain="TRB", plot="bar", order="variance", scale=TRUE)
dev.off()
clonalOverlap(combined, cloneCall = "gene+nt", method="jaccard")
### Clonotype overlap by Jaccard Index suggests no major overlap in gene+nt between Control and TSLP samples
TSLP_Ear.data <- Read10X("../Data/TSLP_Ear_outs/count/sample_feature_bc_matrix/")
Cont_Ear.data <- Read10X("../Data/Cont_Ear_outs/count/sample_feature_bc_matrix/")
### Utilizing recommended Seurat Workflow
TSLP_Ear <- CreateSeuratObject(counts=TSLP_Ear.data, project="TSLP_Ear", min.cells=3, min.features=200)
Cont_Ear <- CreateSeuratObject(counts=Cont_Ear.data, project="Cont_Ear", min.cells=3, min.features=200)
Cont_Ear[["percent.mt"]] <- PercentageFeatureSet(Cont_Ear, pattern="^mt-")
TSLP_Ear[["percent.mt"]] <- PercentageFeatureSet(TSLP_Ear, pattern="^mt-")
### Saving QC Plots
png(filename="TSLP_Ear_QCplots.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
VlnPlot(TSLP_Ear, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()
plot1 <- FeatureScatter(TSLP_Ear, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TSLP_Ear, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="TSLP_Ear_QCplots_2.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
plot1 + plot2
dev.off()
plot1 <- FeatureScatter(Cont_Ear, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Cont_Ear, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(filename="Cont_Ear_QCplots.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
VlnPlot(Cont_Ear, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()
png(filename="Cont_Ear_QCplots_2.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
plot1 + plot2
dev.off()
### Removing Doublets and Dead cells by checking against nFeature_RNA and percent.mt
quantile(Cont_Ear$nFeature_RNA, probs=c(0,0.9,0.95,1))
quantile(Cont_Ear$percent.mt, probs=c(0,0.9,0.95,1))
Cont_Ear1 <- subset (Cont_Ear, subset = nFeature_RNA < 3200 & percent.mt < 3.953)
quantile(TSLP_Ear$nFeature_RNA, probs=c(0,0.9,0.95,1))
quantile(TSLP_Ear$percent.mt, probs=c(0,0.9,0.95,1))
TSLP_Ear1 <- subset (TSLP_Ear, subset = nFeature_RNA < 3616 & percent.mt < 2.954)
### Removed top 5% from each measure to account for doublets and dead cells
Ear <- merge(TSLP_Ear1, y=Cont_Ear1, add.cell.ids=c("TSLP_Ear_TSLP","Cont_Ear_Control"))
Ear <- combineExpression(combined, Ear, cloneCall="gene", group.by = "sample", proportion = FALSE, cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
colnames(Ear@meta.data)
### Adding in Jaccard Plot for lab use
png(filename="TSLP_vs_Cont_jaccard.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
clonalOverlap(combined, cloneCall = "gene+nt", method = "jaccard")
dev.off()
vizGenes(combined, gene = "V", chain="TRB", plot="bar", order="variance", scale=TRUE)
ggsave(filename="TSLP_vs_Cont_TRB_V_usage.png")
### Back to Seurat
### TCR has been merged to seurat object
Ear <- NormalizeData(Ear)
Ear <- FindVariableFeatures(Ear, selection.method="vst", nfeatures=2000)
plot1 <- VariableFeaturePlot(Ear)
top20 <- head(VariableFeatures(Ear), 20)
plot2 <- LabelPoints(plot=plot1, points=top20, repel=TRUE,xnudge=0,ynudge=0)
png(filename="Ear_only_varfeature.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
plot1+plot2
dev.off()
all.genes <- rownames(Ear)
Ear <- ScaleData(Ear, features=all.genes)
Ear <- RunPCA(Ear, features=VariableFeatures(object=Ear))
Ear <- JackStraw(Ear, num.replicate=100)
Ear <- ScoreJackStraw (Ear, dims=1:20)
JackStrawPlot(Ear, dims=1:20)
### Jackstraw ranking isn't coming up with a good discriminant for dimension cutoff
DimHeatmap(Ear, dims=1:20, cells=500, balanced=TRUE)
ElbowPlot(Ear, ndims=50)
### DimHeatmap suggests >20, ElbowPlot looks like mid 20s is best. Select 30 for dimension cutoff
Ear <- FindNeighbors(Ear, dims=1:30)
Ear <- FindClusters(Ear, resolution=0.3)
Ear <- RunUMAP(Ear, dims=1:30)
png(filename="Ear_UMAP_res0_3_split_sample.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
DimPlot(Ear, reduction="umap", split.by="orig.ident", pt.size=1, label=TRUE, label.size=12)
dev.off()
Ear <- FindClusters(Ear, resolution=0.5)
png(filename="Ear_UMAP_res0_5_split_sample.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
DimPlot(Ear, reduction="umap", split.by="orig.ident", pt.size=1, label=TRUE, label.size=12)
dev.off()
### Resolution of 0.5 brings about more distinction between the uppermost clusters
### Subset based on TCR
EarT <- subset(Ear, subset=Frequency >0)
NROW(EarT@meta.data)
NROW(Ear@meta.data)
NROW(EarT@meta.data[which(EarT@meta.data$orig.ident=="Cont_Ear"),])
NROW(EarT@meta.data[which(EarT@meta.data$orig.ident=="TSLP_Ear"),])
NROW(Ear@meta.data[which(Ear@meta.data$orig.ident=="Cont_Ear"),])
NROW(Ear@meta.data[which(Ear@meta.data$orig.ident=="TSLP_Ear"),])
### TSLP ear went from 5476 cellIDs to 4926 single, live cells; 4290 of these cells have at least 1 TCR chain
### Control ear went from 6744 cellIDs to 6067 single, live cells; 3553 of these cells have at least 1 TCR chain
### Final EarT object has 7843 cells
names(EarT@meta.data)
EarT@meta.data <- EarT@meta.data[,1:(ncol(EarT@meta.data)-3)]
EarT <- NormalizeData(EarT)
EarT <- FindVariableFeatures(EarT, selection.method="vst", nFeatures=2000)
plot1 <- VariableFeaturePlot(EarT)
top20 <- head(VariableFeatures(EarT), 20)
plot2 <- LabelPoints(plot=plot1, points=top20, repel=TRUE, xnudge=0, ynudge=0)
png(filename="Ear_Ts_only_varfeature.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
plot1+plot2
dev.off()
all.genes <- rownames(EarT)
EarT <- ScaleData(EarT, features=all.genes)
EarT <- RunPCA(EarT, features=VariableFeatures(object=EarT))
DimHeatmap(EarT, dims=1:20, cells=500, balanced=TRUE)
ElbowPlot(EarT, ndims=50)
### DimHeatmap suggests >20, ElbowPlot looks like early 20s is best. Select 30 for dimension cutoff and consistency
EarT <- FindNeighbors(EarT, dims=1:30)
EarT <- FindClusters(EarT, resolution=0.5)
### Using resolution 0.5 as there were a majority T cells in the entire dataset previously
EarT <- RunUMAP(EarT, dims=1:30)
png(filename="Ear_Ts_UMAP_res0_5_split_sample.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
DimPlot(EarT, reduction="umap", split.by="orig.ident", pt.size=1, label=TRUE, label.size=12)
dev.off()
table(EarT@meta.data[,c("seurat_clusters","orig.ident")])
### Appears that clusters 0 and 3 are highly enriched in TSLP_Ear, and maybe cluster 8/9
### Control has larger cluster 4, 5, and 7, maybe 10/11
TSLP.T.markers <- FindAllMarkers(EarT, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, test.use="MAST")
TSLP.T.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) -> top10
### Grouping clusters that have more than 5x cell numbers in TSLP condition than Control
TSLP.clusters <- c(0,3,8,9)
TSLP.cluster.markers <- FindMarkers(EarT, ident.1=TSLP.clusters, min.pct=0.25, test.use="MAST")
write.csv(TSLP.cluster.markers, file="20220318_TSLP_0_3_8_9_cluster_markers.csv", quote=FALSE)
TSLP.cluster.markers %>% top_n(n=10, wt=avg_log2FC) -> TSLP.top10
DotPlot(EarT, features=rownames(TSLP.top10))
TSLP.T.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC) -> top5
png(filename="Ear_Ts_all_cluster_top5_heatmap.png", width=24*dpi, height = 16*dpi, units="px", type="Cairo")
DoHeatmap(subset(EarT, downsample=100), features=top5$gene, size=12)
dev.off()
TSLP.T.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) -> top10
png(filename="Ear_Ts_all_cluster_top5_heatmap.png", width=24*dpi, height = 16*dpi, units="px", type="cairo")
DoHeatmap(subset(EarT, downsample=100), features=top5$gene, size=48) + theme(axis.text.y=element_text(size=48), text = element_text(size = 24))
dev.off()
png(filename="Ear_Ts_all_cluster_top10_heatmap.png", width=24*dpi, height = 16*dpi, units="px", type="cairo")
DoHeatmap(subset(EarT, downsample=100), features=top10$gene, size=48) + theme(axis.text.y=element_text(size=24), text = element_text(size = 24))
dev.off()

### Making Paper ready Figures
png(filename="Ear_Ts_UMAP_res0_5_split_sample.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
DimPlot(EarT, reduction="umap", split.by="orig.ident", pt.size=1, label=TRUE, label.size=24)
dev.off()
png(filename="TSLP_vs_Cont_Clonotype_fraction.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
clonalHomeostasis(combined, cloneCall="gene", cloneTypes=c(Rare = 1e-04, Small = 0.001,  Medium = 0.01, Large = 0.1,  Hyperexpanded = 1)) + theme(text = element_text(size = 48))
dev.off()
png(filename="Ear_Ts_UMAP_res0_5_Cd4_Cd8a.png", width=8*dpi, height = 4*dpi, units="px", type="cairo")
FeaturePlot(EarT, features=c("Cd4", "Cd8a"), cols=c("lightgrey", "darkblue"), pt.size=1) & theme(plot.title=element_text(size=96), text=element_text(size=24))
dev.off()
write.csv(top10$gene, file="top10gene.csv", quote=FALSE)
# save.image(file="20220318_scRepertoire.RData")
### Reading in Ruth's generalist
# genelist <- read.csv("../genelist.csv", header=FALSE)
# DoHeatmap(subset(EarT, downsample=100), features=genelist$V1)
### Looking specifically at Th2 pathway + Crlf2
Th2genes <- c("Crlf2", "Gata3", "Il4", "Il5", "Il13")
DotPlot(EarT, features=Th2genes, cols=c("lightgrey","darkblue"))
ggsave(file="Ear_Ts_Th2_dotplot.png")

### trying fast gsea
library(fgsea)
library(msigdbr)
### MsigDB doesn't play well with mouse genes, using msigdbr for help
TSLP.clusters <- c(0,3,8,9)
Control.Cd4.clusters <- c(1,2,13)
clustering <- EarT@active.ident
cellsIn <- names(clustering[clustering %in% TSLP.clusters])
cellsIn <- cellsIn[grep("Control",cellsIn,invert=TRUE)]
cellsOut <- names(clustering[clustering %in% Control.Cd4.clusters])
cellsOut <- cellsOut[grep("Control",cellsOut,invert=TRUE)]
fgsea.markers <- FindMarkers(EarT, cellsIn, cellsOut, logfc.threshold=0, min.pct=0, test.use="MAST")
ranks <- fgsea.markers$avg_log2FC
names(ranks) <- rownames(fgsea.markers)
immune_mouse_set = msigdbr(species="mouse", category="C7", subcategory="IMMUNESIGDB")
msigdbr_list = split (x=immune_mouse_set$gene_symbol, f=immune_mouse_set$gs_name)
fgseaRes <- fgsea(pathways=msigdbr_list, stats=ranks, minSize=10, maxSize=500)
str(fgseaRes)
### code for fgsea from scdata taken from https://github.com/ctlab/fgsea/issues/50

head(fgseaRes[order(pval),], n=10)
length(which(str_detect(fgseaRes$pathway, "TH2",)))
Th2.fgseaRes <- fgseaRes[(which(str_detect(fgseaRes$pathway, "TH2"))),]
T.fgseaRes <- fgseaRes[(which(str_detect(fgseaRes$pathway, "TCELL|TCONV|TH1|TH2|TH17"))),]
head(Th2.fgseaRes[order(pval),], n=10)
fwrite(fgseaRes, file="20220322_TSLP_vs_Control_CD4_fgsea_analysis.txt", sep="\t", sep2=c("", " ", ""))
fwrite(T.fgseaRes, file="20220322_TSLP_vs_Control_CD4_fgsea_analysis_Tcell.txt", sep="\t", sep2=c("", " ", ""))
fwrite(Th2.fgseaRes, file="20220322_TSLP_vs_Control_CD4_fgsea_analysis_Th2.txt", sep="\t", sep2=c("", " ", ""))

### wrote tables for analyzing pathways that might be of interest, making gsea plots now
head(Th2.fgseaRes[order(pval),], n=10)
plotEnrichment(msigdbr_list[["GSE14308_TH2_VS_NAIVE_CD4_TCELL_UP"]], ranks) + labs(title="Th2 vs Naive CD4 T")
ggsave(filename="TSLP_vs_Control_Cd4s_Th2_vs_Naive_GSEA.png")
plotEnrichment(msigdbr_list[["GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_DN"]], ranks) + labs(title="Naive vs Memory CD4 T")
ggsave(filename="20220322_TSLP_vs_Control_Cd4s_Naive_vs_Memory_GSEA.png")

### 20230611_edits_for_reviewers
library(Seurat)
library(ggplot2)
library(tidyverse)
library(scRepertoire)
setwd("/Users/enjunyang/Documents/Personal/Taku/20230611_review_analysis/")
load("../20220312_Analysis/20220318_scRepertoire.RData")
### UMAP for all cells
plot1 <- DimPlot(Ear, reduction="umap", split.by="orig.ident", pt.size=1, label=TRUE, label.size=24, cols="alphabet")
plot1 <- plot1 + theme(strip.text.x = element_text(size=96), axis.text=element_text(size=24))
png(filename="20230617_Ear_UMAP_res0_5_split_sample_cols.png", width=8*dpi, height = 8*dpi, units="px", type="cairo")
plot1
dev.off()

### Redrawing figure E1 (S1)
png(filename="20230617_Ear_Ts_UMAP_res0_5_Cd4_Cd8a.png", width=8*dpi, height = 8*dpi, units="px", type="cairo")
FeaturePlot(EarT, features=c("Cd4", "Cd8a"), cols=c("lightgrey", "darkblue"), pt.size=2, split.by="orig.ident") & theme(plot.title=element_text(size=96), text=element_text(size=24))
dev.off()

### Annotating EarT clusters manually
### 1st Genelist from https://www.nature.com/articles/s41467-020-15543-y (CD4s)
cd4markers <- read.table("CD4_markergenes.tsv", sep="\t")
cd4markers$V2 <- stringr::str_to_title(cd4markers$V2)
png(filename="20230617_Ear_Ts_CD4_subtypes_heatmap.png", width=24*dpi, height = 16*dpi, units="px", type="cairo")
DoHeatmap(subset(EarT, downsample=100), features=cd4markers$V2, size=48) + theme(axis.text.y=element_text(size=24), text = element_text(size = 24))
dev.off()
### 2nd gene list from doi: 10.4049/jimmunol.2100581
cd4markers <- c("Cd3d", "Cd4","Cd8b1", "Trdv4","Trdc", "Ccr7","Sell", "Lef1","Timp1","Usp10","Tnfrsf4", "Il2ra","Foxp3","Ccr10","Irf7","Oas1a","Cxcr3","Ctsw","Ccl5","Nkg7","Gzma", "Gzmk","Klrg1","Gzmb","Klrb1","Zbtb16","Rorc","Ifng","Csf2","Stat4","Stat1","Il4","Il5","Il13","Gata3","Il17a","Il22","Stat3","Il10")
plot <- DotPlot(EarT, assay="RNA", features=cd4markers, cols=c("lightgrey","darkblue"))
plot + RotatedAxis()
ggsave(filename="20230617_Ear_Ts_CD4_subtypes_dotplot.png")
# cytokines <- c("Ifng","Csf2","Stat4","Stat1","Il17","Il22","Stat3","Rorc","Il10","Foxp3")
### Manual annotations from literature applied
new.cluster.ids <- c("Th2", "Th2", "Tregs", "Th2", "Naïve CD8", "Naïve CD4", "CD8 TEMRA", "MAIT", "Th2", "Th2", "Gamma-delta T", "Effector T", "Unclassified_1", "Th1", "Unclassified_2")
names(new.cluster.ids) <- levels(EarT)
EarT <- RenameIdents(EarT, new.cluster.ids)
png(filename="20230617_Ear_Ts_UMAP_res0_5_man_annot_cols.png", width=8*dpi, height = 8*dpi, units="px", type="cairo")
DimPlot(EarT, reduction="umap", pt.size=2, label=TRUE, label.size=24, repel=TRUE, cols="polychrome") + theme(axis.text=element_text(size=36), axis.title=element_text(size=48))
dev.off()
Idents(EarT) <- "seurat_clusters"
png(filename="20230617_Ear_Ts_UMAP_res0_5_old_annot.png", width=8*dpi, height = 8*dpi, units="px", type="cairo")
DimPlot(EarT, reduction="umap", pt.size=2, label=TRUE, label.size=24, repel=TRUE) + theme(axis.text=element_text(size=36), axis.title=element_text(size=48))
dev.off()

### checking for gammadelta T cells
TRGD <- c(grep("Trd", EarT@assays$RNA@counts@Dimnames[[1]]), grep("Trg", EarT@assays$RNA@counts@Dimnames[[1]]))
TRGD <- EarT@assays$RNA@counts@Dimnames[[1]][TRGD]
TRGD <- TRGD[c(5,6,7,9)]
FeaturePlot(EarT, features=TRGD)
ggsave(filename="20230617_EarT_gammadeltas.png")

### Rerunning GSEA analysis without preselecting for condition cluster comparisons
library(fgsea)
library(msigdbr)
TSLP.clusters <- c(0,3,8,9)
Control.Cd4.clusters <- c(1,2,13)
clustering <- EarT@active.ident
cellsIn <- names(clustering[clustering %in% TSLP.clusters])
cellsOut <- names(clustering[clustering %in% Control.Cd4.clusters])
fgsea.markers <- FindMarkers(EarT, cellsIn, cellsOut, logfc.threshold=0, min.pct=0, test.use="MAST")
ranks <- fgsea.markers$avg_log2FC
names(ranks) <- row.names(fgsea.markers)
immune_mouse_set = msigdbr(species="mouse", category="C7", subcategory="IMMUNESIGDB")
msigdbr_list = split (x=immune_mouse_set$gene_symbol, f=immune_mouse_set$gs_name)
fgseaRes <- fgsea(pathways=msigdbr_list, stats=ranks, minSize=10, maxSize=500)
str(fgseaRes)
head(fgseaRes[order(pval),], n=10)
length(which(str_detect(fgseaRes$pathway, "TH2",)))
Th2.fgseaRes <- fgseaRes[(which(str_detect(fgseaRes$pathway, "TH2"))),]
T.fgseaRes <- fgseaRes[(which(str_detect(fgseaRes$pathway, "TCELL|TCONV|TH1|TH2|TH17"))),]
head(Th2.fgseaRes[order(pval),], n=10)
data.table::fwrite(fgseaRes, file="20230617_TSLP_vs_Control_CD4_fgsea_analysis.txt", sep="\t", sep2=c("", " ", ""))
data.table::fwrite(T.fgseaRes, file="20230617_TSLP_vs_Control_CD4_fgsea_analysis_Tcell.txt", sep="\t", sep2=c("", " ", ""))
data.table::fwrite(Th2.fgseaRes, file="20230617_TSLP_vs_Control_CD4_fgsea_analysis_Th2.txt", sep="\t", sep2=c("", " ", ""))
head(Th2.fgseaRes[order(pval),], n=20)
plotEnrichment(msigdbr_list[["GSE14308_TH2_VS_NAIVE_CD4_TCELL_UP"]], ranks) + labs(title="Th2 vs Naive CD4 T UP")
ggsave(filename="20230617_TSLP_vs_Control_Cd4s_Th2_vs_Naive_Up_GSEA.png")
plotEnrichment(msigdbr_list[["GSE11924_TFH_VS_TH2_CD4_TCELL_DN"]], ranks) + labs(title="TFH vs TH2 CD4 T DN")
ggsave(filename="20230617_TSLP_vs_Control_Cd4s_TFH_vs_TH2_DN_GSEA.png")

### redrawing CD4,CD8 plots
png(filename="20230617_Ear_Ts_UMAP_res0_5_Cd4_Cd8a.png", width=8*dpi, height = 8*dpi, units="px", type="cairo")
FeaturePlot(EarT, features=c("Cd4", "Cd8a"), cols=c("lightgrey", "darkblue"), split.by="orig.ident", pt.size=1, label=TRUE, label.size=24) & theme(plot.title=element_text(size=96), axis.title= element_text(size=48), axis.text=element_text(size=36))
dev.off()

### Trying scType for annotating automatically for all cells and Ear T cells
### https://github.com/IanevskiAleksandr/sc-type/
library(HGNChelper)
library(openxlsx)
library(dplyr)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system")
es.max <- sctype_score(scRNAseqData = Ear[["RNA"]]@scale.data, scaled=TRUE, gs=gs_list$gs_positive, gs2=gs_list$gs_negative)
cL_results = do.call("rbind", lapply(unique(Ear@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(Ear@meta.data[Ear@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(Ear@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
Ear@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  Ear@meta.data$customclassif[Ear@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
png(filename="20230617_scType_Ear_all_cells_labels_updated.png", width=4*dpi, height = 4*dpi, units="px", type="cairo")
DimPlot(Ear, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif', pt.size=1, label.size=12)  + ggtitle("Custom_Classifier") + theme(plot.title = element_text(face="bold", hjust = 0.5, size=36))
dev.off()

ob.list <- SplitObject(Ear, split.by = "orig.ident")
plot.list <- lapply(X = ob.list, FUN = function(x) {
  DimPlot(x, reduction = "umap", label = TRUE, label.size = 12, pt.size=1) + ggtitle(x@meta.data$orig.ident[1]) + theme(plot.title = element_text(face="bold", hjust = 0.5, size=36))
})
head(ob.list[[1]]@meta.data$orig.ident)
png(filename="20230617_Ear_UMAP_res0_5_Ctrl_sample.png", width=4*dpi, height = 4*dpi, units="px", type="cairo")
plot.list[[2]]
dev.off()
png(filename="20230617_Ear_UMAP_res0_5_TSLP_sample.png", width=4*dpi, height = 4*dpi, units="px", type="cairo")
plot.list[[1]]
dev.off()
saveRDS(Ear, file="20230617_Ear_bak.RDS")

es.max <- sctype_score(scRNAseqData = EarT[["RNA"]]@scale.data, scaled=TRUE, gs=gs_list$gs_positive, gs2=gs_list$gs_negative)
cL_results = do.call("rbind", lapply(unique(EarT@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(EarT@meta.data[EarT@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(EarT@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
EarT@meta.data$sctype_class = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  EarT@meta.data$sctype_class[EarT@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}
DimPlot(EarT, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'sctype_class')
ggsave(filename="20230617_scType_EarT_fail_labels.png")
save.image(file="20230617_scRepertoire_reviews.RData")

