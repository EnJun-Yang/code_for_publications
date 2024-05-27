setwd("~/Desktop/Weiyee/20220527_Rhapsody/20240122_TAC/")
library(Seurat)
library(ggplot2)
library(dplyr)
library(reticulate)
library(knitr)

### Processing for 0h timepoint ###

meta <- read.table(file="../Rawdata/20220507_Rhapsody_metadata.csv", sep=",", header=TRUE)
df_0h <- read.table(file="../Rawdata/Combined_20220527-WYO-0h-_RSEC_MolsPerCell.csv", sep = ",", header=TRUE, row.names=1)
df_0h <- data.frame(t(df_0h), check.names=FALSE)
dim(df_0h)
DC_0h <- CreateSeuratObject(counts=df_0h, project="DC_0h")
DC_0h[["percent.mt"]] <- PercentageFeatureSet(DC_0h, pattern="^MT.")
VlnPlot(DC_0h, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(filename="DC_0h_QC_violin_all.png")
### First three rows are abseq [1] "CD141.Thbd.AMM2163.pAbO" "CD1c.CD1C.AHS0088.pAbO" [3] "Neuropilin.1.NRP1.AHS0315.pAbO"
### Appending sample tag calls to dataset
DC_0h_meta <- read.table("../Rawdata/20220527-WYO-0h-_Sample_Tag_Calls.csv", skip = 0, sep = ",", header = TRUE, row.names = 1)
DC_0h <- AddMetaData(object = DC_0h, metadata = DC_0h_meta)
### Dropping SMK multiplets
DC_0h <- subset(x=DC_0h, subset = Sample_Tag != "Multiplet")
DC_0h <- subset(x=DC_0h, subset = Sample_Tag != "Undetermined")
NCOL(DC_0h)
# left with 6140 cells out of 7120
VlnPlot(DC_0h, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(filename="DC_0h_QC_violin_smk_singlets.png")
### Trim off cells with <200 or >3000 bioproducts, > 10000 RNAs or > 20% mito genes
DC_0h <- subset(DC_0h, subset = percent.mt < 25 & nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 10000)
NCOL(DC_0h)
# left with 5766 cells out of 7120
DC_0h <- NormalizeData(DC_0h)
DC_0h <- FindVariableFeatures(DC_0h, selection.method = "vst", nfeatures = 2000)
top10variable <- head(VariableFeatures(DC_0h), 10)
vfp <- VariableFeaturePlot(DC_0h)
vfp_labels <- LabelPoints(plot=vfp, points=top10variable, repel = TRUE, xnudge=0, ynudge=0)
vfp_labels
ggsave(filename="DC_0h_variable_features.png")
# all.rnagenes <- rownames(DC_0h)[4:NROW(DC_0h)]
all.genes <- rownames(DC_0h)
DC_0h <- ScaleData(DC_0h, features = all.genes, model.use="linear")
DC_0h <- RunPCA(DC_0h, features=VariableFeatures(object=DC_0h), verbose=F)
VizDimLoadings(DC_0h, dims=1:4, reduction ="pca")
ggsave(filename="DC_0h_PCA_top4.png", width=8.26, height=10)
ElbowPlot(DC_0h, ndims=50)
ggsave(filename="DC_0h_elbow.png")
### Use about 30 Dimensions to preserve variability
DC_0h <- FindNeighbors(DC_0h, dims = 1:30)
DC_0h <- FindClusters(DC_0h, resolution=0.5, random.seed = 0)
DC_0h <- RunUMAP(DC_0h, dims=1:30, seed.use = 42)
DimPlot(DC_0h, reduction="umap", pt.size=0.1)
ggsave(filename="DC_0h_umap.png")
DimPlot(DC_0h, group.by="Sample_Tag", pt.size=0.1)
ggsave(filename="DC_0h_umap_by_ST.png")
FeaturePlot(DC_0h, features=rownames(DC_0h)[c(1:3)], pt.size=0.1)
ggsave(filename="DC_0h_abseq.png")
### Adding donor and condition information
str(DC_0h@meta.data)
DC_0h$donor <- plyr::mapvalues(
  x = DC_0h$Sample_Name, 
  from = meta$Sample_Name,
  to = meta$donor
)
DC_0h$condition <- plyr::mapvalues(
  x = DC_0h$Sample_Name, 
  from = meta$Sample_Name,
  to = meta$condition
)
RidgePlot(DC_0h, features=rownames(DC_0h)[c(1:3)])
ggsave(filename="DC_0h_abseq_ridge.png")
saveRDS(DC_0h, file="DC_0h.RDS")
rm(DC_0h)
rm(df_0h)

### Repeat for 4h ###

df_4h <- read.table(file="../Rawdata/Combined_20220527-WYO-4h-_RSEC_MolsPerCell.csv", sep = ",", header=TRUE, row.names=1)
df_4h <- data.frame(t(df_4h), check.names=FALSE)
dim(df_4h)
DC_4h <- CreateSeuratObject(counts=df_4h, project="DC_4h")
DC_4h[["percent.mt"]] <- PercentageFeatureSet(DC_4h, pattern="^MT.")
VlnPlot(DC_4h, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(filename="DC_4h_QC_violin_all.png")
### Appending sample tag calls to dataset
DC_4h_meta <- read.table("../Rawdata/20220527-WYO-4h-_Sample_Tag_Calls.csv", skip = 0, sep = ",", header = TRUE, row.names = 1)
DC_4h <- AddMetaData(object = DC_4h, metadata = DC_4h_meta)
DC_4h <- subset(x=DC_4h, subset = Sample_Tag != "Multiplet")
DC_4h <- subset(x=DC_4h, subset = Sample_Tag != "Undetermined")
NCOL(DC_4h)
# left with 13082 cells out of 17164
VlnPlot(DC_4h, features=c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
ggsave(filename="DC_4h_QC_violin_smk_singlets.png")
### Trim off cells with <200 or >3000 bioproducts, > 10000 RNAs or > 20% mito genes
DC_4h <- subset(DC_4h, subset = percent.mt < 25 & nFeature_RNA > 200 & nFeature_RNA < 3000 & nCount_RNA < 20000)
NCOL(DC_4h)
# left with 12795 cells out of 17164
DC_4h <- NormalizeData(DC_4h)
DC_4h <- FindVariableFeatures(DC_4h, selection.method = "vst", nfeatures = 2000)
top10variable <- head(VariableFeatures(DC_4h), 10)
vfp <- VariableFeaturePlot(DC_4h)
vfp_labels <- LabelPoints(plot=vfp, points=top10variable, repel = TRUE, xnudge=0, ynudge=0)
vfp_labels
ggsave(filename="DC_4h_variable_features.png")
all.genes <- rownames(DC_4h)
DC_4h <- ScaleData(DC_4h, features = all.genes, model.use="linear")
DC_4h <- RunPCA(DC_4h, features=VariableFeatures(object=DC_4h), verbose=F)
VizDimLoadings(DC_4h, dims=1:4, reduction ="pca")
ggsave(filename="DC_4h_PCA_top4.png", width=8.26, height=10)
ElbowPlot(DC_4h, ndims=50)
ggsave(filename="DC_4h_elbow.png")
### Use about 30 Dimensions to preserve variability
DC_4h <- FindNeighbors(DC_4h, dims = 1:30)
DC_4h <- FindClusters(DC_4h, resolution=0.5, random.seed = 0)
DC_4h <- RunUMAP(DC_4h, dims=1:30, seed.use = 42)
DimPlot(DC_4h, reduction="umap", pt.size=0.1, label=TRUE)
ggsave(filename="DC_4h_umap.png")
DimPlot(DC_4h, group.by="Sample_Tag", pt.size=0.1)
ggsave(filename="DC_4h_umap_by_ST.png")
FeaturePlot(DC_4h, features=rownames(DC_4h)[c(1:3)], pt.size=0.1)
ggsave(filename="DC_4h_abseq.png")
### Adding condition information to metadata
str(DC_4h@meta.data)
DC_4h$donor <- plyr::mapvalues(
  x = DC_4h$Sample_Name, 
  from = meta$Sample_Name,
  to = meta$donor
)
DC_4h$condition <- plyr::mapvalues(
  x = DC_4h$Sample_Name, 
  from = meta$Sample_Name,
  to = meta$condition
)
DimPlot(DC_4h, group.by="donor", pt.size=0.1)
ggsave(filename="DC_4h_umap_by_donor.png")
DimPlot(DC_4h, group.by="condition", pt.size=0.1)
ggsave(filename="DC_4h_umap_by_condition.png")
RidgePlot(DC_4h, features=rownames(DC_4h)[c(1:3)])
ggsave(filename="DC_4h_abseq_ridge.png")
saveRDS(DC_4h, file="DC_4h.RDS")

###Integrating data ###

DC_0h <- readRDS("DC_0h.RDS")
DC_4h <- readRDS("DC_4h.RDS")
DC_list <- list(DC_0h, DC_4h)
features <- SelectIntegrationFeatures(object.list=DC_list)
anchors <- FindIntegrationAnchors(object.list= DC_list, anchor.features = features)
DC.combined <- IntegrateData(anchorset=anchors)
DefaultAssay(DC.combined) <- "integrated"
### Standard integrated workflow
DC.combined <- ScaleData(DC.combined, verbose = FALSE)
DC.combined <- RunPCA(DC.combined, verbose=FALSE)
ElbowPlot(DC.combined, ndims=50)
ggsave(filename="DC_0h_4h_elbow.png")
DC.combined <- RunUMAP(DC.combined, reduction="pca", dims=1:30, seed.use = 42)
DC.combined <- FindNeighbors(DC.combined, reduction="pca", dims=1:30)
DC.combined <- FindClusters(DC.combined, resolution=1.2, random.seed = 0)
p1 <- DimPlot(DC.combined, reduction = "umap", group.by = "condition")
p2 <- DimPlot(DC.combined, reduction = "umap", label=TRUE, repel=FALSE)
p3 <- DimPlot(DC.combined, reduction = "umap", group.by = "orig.ident")
p4 <- DimPlot(DC.combined, reduction = "umap", group.by = "donor")
p1+p2+p3+p4
dpi=300
ggsave(file="DC_0h_4h_umaps.png", height=16, width=16)
FeaturePlot(DC.combined, features=c(rownames(DC_0h)[c(1:3)], "CD1C"), pt.size=0.1)
ggsave(file="DC_0h_4h_abseq_umaps.png", height=8, width=8)
RidgePlot(DC.combined, features=rownames(DC_0h)[c(1:3)])
ggsave(filename="DC_0h_4h_abseq_ridge.png")
DimPlot(DC.combined, reduction="umap", split.by="condition", pt.size = 0.1)
ggsave(filename="DC_0h_4h_split_condition.png")
DimPlot(DC.combined, reduction="umap", split.by="orig.ident", pt.size = 0.1)
ggsave(filename="DC_0h_4h_split_timepoint.png")
saveRDS(DC.combined, file = "DC_0h_4h.RDS")
#DC.combined <- readRDS("../20230308_clean/DC_0h_4h.RDS")

### Check for B cells and other Monos/DCs/Baso ###
FeaturePlot(DC.combined, features=c("MS4A1","CD14","FCGR3A","FCER1A"), label=TRUE, pt.size=0.1)
ggsave(filename="DC_0h_4h_cell_markers.png")
all.markers <- FindAllMarkers(DC.combined, test.use="MAST", assay="RNA", logfc.threshold = 0.585)
write.table(all.markers, file="0h_4h_all_DEgenes_RNA.csv", sep=",", quote = FALSE)
DC.combined.small <- subset(DC.combined, downsample=50)
DefaultAssay(DC.combined.small) <- "RNA"
DC.combined.small <- ScaleData(DC.combined.small, verbose = FALSE)
all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(DC.combined.small, features = top5$gene) + NoLegend() + theme(axis.text.y=element_text(size = 8))
ggsave(file="0h_4h_heatmap.png", height=16, width=16)
### high (but not highest) CD1C staining is in clusters 6,5,3. Maybe 0,1
### Comparing highest level of CD1c expression by average
abseq <- rownames(DC_0h)[c(2,3)]
AverageExpression(DC.combined, features=abseq, assays="RNA")
AverageExpression(DC.combined, features=abseq, assays="RNA", group.by="donor")
avg.dat <- as.data.frame(AverageExpression(DC.combined, features=abseq, assays="RNA")[[1]])
boxplot(t(avg.dat[1,]))
avg.dat[,which(avg.dat[1,] > 150)]
FeaturePlot(DC.combined, features=c("CD1C", "MS4A1"), label=TRUE, slot="counts")
### For supplementary figure 20240112
DefaultAssay(DC.combined) <- "RNA"
DimPlot(DC.combined, label=TRUE, pt.size=0.2, label.size=12)
ggsave(file="20240112_0h_4h_UMAP.png", height=16, width=16, dpi=400)
feats <- c("MS4A1","CD14","FCGR3A","FCER1A")
VlnPlot(DC.combined, features=feats)
ggsave(file="20240112_0h_4h_other_marker_vlns.png", height=8, width=16, dpi=400)
feats <- rownames(DC.combined)[grep("pAb", rownames(DC.combined))]
VlnPlot(DC.combined, features=feats)
ggsave(file="20240112_0h_4h_Abseq_vlns.png", height=4, width=16, dpi=400)
feats <- c("THBD", "CD1C", "NRP1")
VlnPlot(DC.combined, features=feats)
ggsave(file="20240112_0h_4h_RNA_vlns.png", height=4, width=16, dpi=400)
feats <- rownames(DC.combined)[grep("pAb", rownames(DC.combined))]
AverageExpression(DC.combined, features=feats, assays="RNA")
write.table(AverageExpression(DC.combined, features=feats, assays="RNA"), file="20240112_abseq_av_expression.txt", sep="\t", quote=FALSE)
### Based on RNA and Abseq, cluster of interest is likely bottom right. Apply sanity cutoff of CD1C abseq 150 < x < 500
### Also, not MS4A1+ (cluster 9), FCGR3A/FCER1A high (cluster 0,1) or NRP1 high (12, 14, 16, 19, 23, 24) --> cluster 4, 5, 7, 11,
DimPlot(DC.combined, label=TRUE)
DC.interest <- subset(DC.combined, idents=c(4,5,7,11))
### dropped from 18561 to 3656 cells
table(DC.interest@meta.data$orig.ident)
table(DC.interest@meta.data$donor)
table(DC.interest@meta.data$condition)
### fewer cells from 0h. Donor 2 has more cells, and even split of treatments (taking into account 0h)
DC.interest <- RunUMAP(DC.interest, reduction="pca", dims=1:30, seed.use = 42)
DC.interest <- FindNeighbors(DC.interest, reduction="pca", dims=1:30)
DC.interest <- FindClusters(DC.interest, resolution=0.3, random.seed = 0)
p1 <- DimPlot(DC.interest, reduction = "umap", group.by = "condition", cols="Set1")
p2 <- DimPlot(DC.interest, reduction = "umap", label=TRUE, repel=FALSE, cols="Accent")
p3 <- DimPlot(DC.interest, reduction = "umap", group.by = "orig.ident")
p4 <- DimPlot(DC.interest, reduction = "umap", group.by = "donor")
p1+p2+p3+p4
ggsave(file="DC_0h_4h_clust4_5_7_11_umaps.png", dpi=dpi, height=8, width=8)
FeaturePlot(DC.interest, features = c("TANK", "TBK1", "HLA.DMB", "HLA.DMA", "IFIH1", "DDX58", "MAVS", "IKBKE", "TLR7", "TLR3", "IRF3", "IRF7"))
ggsave(file="DC_0h_4h_clust4_5_7_11_featureplots.png", dpi=300, height=12, width=16)
CD1c.markers <- FindAllMarkers(DC.interest, assay="RNA", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, test.use="MAST")
CD1c.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) -> top10
#write.table(top10$gene, file="20240122_Heatmap_redrawn_legend.csv", sep=",", quote=FALSE, col.names=FALSE, row.names = FALSE)
write.table(top10, file="20240122_Heatmap_redrawn_legend.csv", sep=",", quote=FALSE, row.names=FALSE)
DC.interest1 <- ScaleData(object=DC.interest, assay="RNA")
DoHeatmap(DC.interest1, features=top10$gene, assay="RNA", slot="scale.data")
DefaultAssay(DC.interest) <- "RNA"
FeaturePlot(DC.interest, features="HLA.DMB", reduction="umap", split.by="condition", pt.size=0.5)
ggsave(file="DC_0h_4h_clust4_5_7_11_HLADMB_split.png", height=8, width=16)
DotPlot(DC.interest, features=c("HLA.DMB", "HLA.DMA", "HLA.DOA", "HLA.DOB"), group.by="condition")
ggsave(file="DC_0h_4h_clust4_5_7_11_HLADMB_dotplot.png", height=8, width=8)
DotPlot(DC.interest, features=c("HLA.DMB", "HLA.DMA", "HLA.DOA", "HLA.DOB"))
ggsave(file="DC_0h_4h_clust4_5_7_11_HLADMB_dotplot_bycluster.png", height=8, width=8)

### DEenrichment plots
library(enrichR)
DEenrichRPlot(DC.interest, ident.1="2", ident.2="3", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA")
ggsave(file="DC_0h_4h_clust4_5_7_11_RSV_vs_Flu_GOplot.png", height=8, width=16)
DEenrichRPlot(DC.interest, ident.1="2", ident.2="0", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA")
ggsave(file="DC_0h_4h_clust4_5_7_11_RSV_vs_UT_GOplot.png", height=8, width=16)
DEenrichRPlot(DC.interest, ident.1="3", ident.2="0", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA")
ggsave(file="DC_0h_4h_clust4_5_7_11_Flu_vs_UT_GOplot.png", height=8, width=16)

### Digging into condition enriched clusters. Cluster numbers for flu is 3; untreated is 0; RSV is 2
interest.markers <- FindMarkers(DC.interest, ident.1="0", ident.2="3", test.use="MAST", assay="RNA")
write.table(interest.markers, file="0h_4h_CD1C_UT_vs_Flu_cluster_DEgenes_RNA.csv", sep=",", quote = FALSE)
UT.Flu.markers.4h.RNA <- interest.markers
interest.markers <- FindMarkers(DC.interest, ident.1="2", ident.2="3", test.use="MAST", assay="RNA")
write.table(interest.markers, file="0h_4h_CD1C_RSV_vs_Flu_cluster_DEgenes_RNA.csv", sep=",", quote = FALSE)
RSV.Flu.markers.4h.RNA <- interest.markers
interest.markers <- FindMarkers(DC.interest, ident.1="0", ident.2="2", test.use="MAST", assay="RNA")
write.table(interest.markers, file="0h_4h_CD1C_UT_vs_RSV_cluster_DEgenes_RNA.csv", sep=",", quote = FALSE)
UT.RSV.markers.4h.RNA <- interest.markers
FeaturePlot(DC.interest, features="TANK", pt.size=1.5)
# ggsave("0h_4h_CD1C_TANK_RNA.png")
FeaturePlot(DC.interest, features=c("TANK","TBK1","IKBKE"), pt.size=1, label=TRUE, label.size= 8)
# ggsave("0h_4h_CD1C_TANK_TBK_IKK_RNA.png")


### Dropping RSV HK from the interest list and rerunning for TAC
### Best to restart integration
DC_0h <- readRDS("DC_0h.RDS")
DC_4h <- readRDS("DC_4h.RDS")
keep <- c("Flu", "RSV", "Untreated")
DC_4h <- subset(x=DC_4h, subset = condition %in% keep)
DC_list <- list(DC_0h, DC_4h)
features <- SelectIntegrationFeatures(object.list=DC_list)
anchors <- FindIntegrationAnchors(object.list= DC_list, anchor.features = features)
DC.combined <- IntegrateData(anchorset=anchors)
DefaultAssay(DC.combined) <- "integrated"
### Standard integrated workflow
DC.combined <- ScaleData(DC.combined, verbose = FALSE)
DC.combined <- RunPCA(DC.combined, verbose=FALSE)
ElbowPlot(DC.combined, ndims=50)
outdir="minus_RSVHK/"
ggsave(filename=paste0(outdir,"DC_0h_4h_elbow.png"))
DC.combined <- RunUMAP(DC.combined, reduction="pca", dims=1:30, seed.use = 42)
DC.combined <- FindNeighbors(DC.combined, reduction="pca", dims=1:30)
DC.combined <- FindClusters(DC.combined, resolution=1.2, random.seed = 0)
p1 <- DimPlot(DC.combined, reduction = "umap", group.by = "condition")
p2 <- DimPlot(DC.combined, reduction = "umap", label=TRUE, repel=FALSE)
p3 <- DimPlot(DC.combined, reduction = "umap", group.by = "orig.ident")
p4 <- DimPlot(DC.combined, reduction = "umap", group.by = "donor")
p1+p2+p3+p4
dpi=300
ggsave(file=paste0(outdir,"DC_0h_4h_umaps.png"), height=16, width=16, dpi=dpi)
FeaturePlot(DC.combined, features=c(rownames(DC_0h)[c(1:3)], "CD1C"), pt.size=0.1)
ggsave(file=paste0(outdir,"DC_0h_4h_abseq_umaps.png"), height=8, width=8, dpi=dpi)
RidgePlot(DC.combined, features=rownames(DC_0h)[c(1:3)])
ggsave(filename=paste0(outdir,"DC_0h_4h_abseq_ridge.png"), dpi=dpi)
DimPlot(DC.combined, reduction="umap", split.by="condition", pt.size = 0.1)
ggsave(filename=paste0(outdir,"DC_0h_4h_split_condition.png"), dpi=dpi)
DimPlot(DC.combined, reduction="umap", split.by="orig.ident", pt.size = 0.1)
ggsave(filename=paste0(outdir,"DC_0h_4h_split_timepoint.png"), dpi=dpi)
saveRDS(DC.combined, file = paste0(outdir,"DC_0h_4h_minus_RSVHK.RDS"))

### Check for B cells and other Monos/DCs/Baso ###
FeaturePlot(DC.combined, features=c("MS4A1","CD14","FCGR3A","FCER1A"), label=TRUE, pt.size=0.1)
ggsave(filename=paste0(outdir,"DC_0h_4h_cell_markers.png"))
all.markers <- FindAllMarkers(DC.combined, test.use="MAST", assay="RNA", logfc.threshold = 0.585)
write.table(all.markers, file=paste0(outdir,"0h_4h_all_DEgenes_RNA.csv"), sep=",", quote = FALSE)
DC.combined.small <- subset(DC.combined, downsample=50)
DefaultAssay(DC.combined.small) <- "RNA"
DC.combined.small <- ScaleData(DC.combined.small, verbose = FALSE)
all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(DC.combined.small, features = top5$gene) + NoLegend() + theme(axis.text.y=element_text(size = 8))
ggsave(file=paste0(outdir,"0h_4h_heatmap.png"), height=16, width=16, dpi=dpi)
### Comparing highest level of CD1c expression by average
abseq <- rownames(DC_0h)[c(2,3)]
AverageExpression(DC.combined, features=abseq, assays="RNA")
AverageExpression(DC.combined, features=abseq, assays="RNA", group.by="donor")
avg.dat <- as.data.frame(AverageExpression(DC.combined, features=abseq, assays="RNA")[[1]])
boxplot(t(avg.dat[1,]))
avg.dat[,which(avg.dat[1,] > 130)]
FeaturePlot(DC.combined, features=c("CD1C", "MS4A1"), label=TRUE, slot="counts")
### high CD1C staining is in clusters 0,1,3,5,6,7,12,21,23 (>130).21 is also NRP1 high. 6 is B (MS4A1) cells, 13,17 are FCGR3A high.  0,1,2,8,21 are FCGR1A high. 3,5,7,12 left. 23 may be background staining?
### Probs 3,5,7,12
### For supplementary figure 20240112
DefaultAssay(DC.combined) <- "RNA"
DimPlot(DC.combined, label=TRUE, pt.size=0.2, label.size=12)
ggsave(file=paste0(outdir, "20240112_0h_4h_UMAP.png"), height=16, width=16, dpi=400)
feats <- c("MS4A1","CD14","FCGR3A","FCER1A")
VlnPlot(DC.combined, features=feats)
ggsave(file=paste0(outdir, "20240112_0h_4h_other_marker_vlns.png"), height=8, width=16, dpi=400)
feats <- rownames(DC.combined)[grep("pAb", rownames(DC.combined))]
VlnPlot(DC.combined, features=feats)
ggsave(file=paste0(outdir, "20240112_0h_4h_Abseq_vlns.png"), height=4, width=16, dpi=400)
feats <- c("THBD", "CD1C", "NRP1")
VlnPlot(DC.combined, features=feats)
ggsave(file=paste0(outdir, "20240112_0h_4h_RNA_vlns.png"), height=4, width=16, dpi=400)
feats <- rownames(DC.combined)[grep("pAb", rownames(DC.combined))]
AverageExpression(DC.combined, features=feats, assays="RNA")
write.table(AverageExpression(DC.combined, features=feats, assays="RNA"), file=paste0(outdir, "20240112_abseq_av_expression.txt"), sep="\t", quote=FALSE)
### Based on RNA and Abseq, apply sanity cutoff of CD1C abseq 130 < x < 350
### Also, not MS4A1+ (cluster 6), CD14 (cluster 7), FCGR3A/FCER1A high (cluster 5,12,13,17,21) or NRP1 high (11,14,21) --> cluster 0,1,3
DimPlot(DC.combined, label=TRUE)
DC.interest <- subset(DC.combined, idents=c(3,5,7,12))
length(Cells(DC.interest))
### dropped from 14997 to 2958 cells
DC.combined.old <- readRDS("../20230308_clean/DC_0h_4h.RDS")
DC.interest.old <- subset(DC.combined.old, idents=c(4,5,7,11))
length(which(Cells(DC.interest.old) %in% Cells(DC.interest)))
### Of which 2591 were in the original subset analysis (out of 3656)
table(DC.interest@meta.data$orig.ident)
table(DC.interest@meta.data$donor)
table(DC.interest@meta.data$condition)
### fewer cells from 0h. Donor 2 has more cells, and even split of treatments (taking into account 0h)
# DC.interest <- subset(DC.combined, idents=c(3,5,7,12))
DefaultAssay(DC.interest) <- "integrated"
DC.interest <- RunUMAP(DC.interest, reduction="pca", dims=1:30, seed.use = 42)
DC.interest <- FindNeighbors(DC.interest, reduction="pca", dims=1:30)
DC.interest <- FindClusters(DC.interest, resolution=0.1, random.seed = 0)
p1 <- DimPlot(DC.interest, reduction = "umap", group.by = "condition", cols="Set1")
p2 <- DimPlot(DC.interest, reduction = "umap", label=TRUE, repel=FALSE, cols="Accent")
p3 <- DimPlot(DC.interest, reduction = "umap", group.by = "orig.ident")
p4 <- DimPlot(DC.interest, reduction = "umap", group.by = "donor")
p1+p2+p3+p4
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_umaps.png"), dpi=dpi, height=8, width=8)
DefaultAssay(DC.interest) <- "RNA"
FeaturePlot(DC.interest, features = c("TANK", "TBK1", "HLA.DMB", "HLA.DMA", "IFIH1", "DDX58", "MAVS", "IKBKE", "TLR7", "TLR3", "IRF3", "IRF7"))
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_featureplots.png"), dpi=300, height=12, width=16)
FeaturePlot(DC.interest, features = c("TRAF3","TRAF6", "TICAM1", "MYD88", "TMEM173"))
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_featureplots_2.png"), dpi=300, height=12, width=8)
FeaturePlot(DC.interest, features = c("SEMA3A", "FLT1", "FGFR1", "FLT3", "FLT4", "EXOC1", "EXOC3", "EXOC4", "EXOC6", "EXOC8", "SLC7A5", "LAT2", "TEAD1", "TEAD2", "TEAD3", "TEAD4"))
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_featureplots_3.png"), dpi=300, height=16, width=16)
cols <- RColorBrewer::brewer.pal(n = 3, name = "Set1") # get colors
RColorBrewer::display.brewer.pal(n = 3, name = "Set1") # see colors
levels(Idents(DC.interest))
index <- c(2,3,1)
cols <- cols[order(index)]
VlnPlot(DC.interest, features=c("TBK1", "TANK", "IKBKE", "IRF3", "IRF7"), cols=cols, pt.size=2)
ggsave(file=paste0(outdir, "20240212_DC_0h_4h_noRSVHK_clust3_5_7_12_IFN_Vlnplots.png"), dpi=300, height=12, width=16)
DotPlot(DC.interest, features=c("TBK1", "TANK", "IKBKE", "IRF3", "IRF7"), group.by="condition", dot.scale=24)
ggsave(file=paste0(outdir, "20240212_DC_0h_4h_noRSVHK_clust3_5_7_12_IFN_Dotplots.png"), dpi=300, height=6, width=8)


### Comparing by condition only
DC.interest[["old.ident"]] <- Idents(DC.interest)
Idents(DC.interest) <- "condition"
CD1c.markers <- FindAllMarkers(DC.interest, assay="RNA", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, test.use="MAST")
CD1c.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_log2FC) -> top10
#write.table(top10$gene, file="20240122_Heatmap_redrawn_legend.csv", sep=",", quote=FALSE, col.names=FALSE, row.names = FALSE)
write.table(top10, file=paste0(outdir, "20240122_Heatmap_redrawn_legend.csv"), sep=",", quote=FALSE, row.names=FALSE)
DC.interest1 <- ScaleData(object=DC.interest, assay="RNA")
DoHeatmap(DC.interest1, features=top10$gene, assay="RNA", slot="scale.data")
ggsave(file=paste0(outdir,"DC_0h_4h_noRSVHK_clust3_5_7_12_heatmap.png"), height=16, width=16, dpi=dpi)
DefaultAssay(DC.interest) <- "RNA"
FeaturePlot(DC.interest, features="HLA.DMB", reduction="umap", split.by="condition", pt.size=0.5)
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_HLADMB_split.png"), height=8, width=16)
DotPlot(DC.interest, features=c("HLA.DMB", "HLA.DMA", "HLA.DOA", "HLA.DOB"), group.by="condition")
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_HLADMB_dotplot.png"), height=8, width=8)
DotPlot(DC.interest, features=c("HLA.DMB", "HLA.DMA", "HLA.DOA", "HLA.DOB"))
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_HLADMB_dotplot_bycluster.png"), height=8, width=8)



### DEenrichment plots
library(enrichR)
DEenrichRPlot(DC.interest, ident.1="RSV", ident.2="Flu", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA")
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_Flu_GOplot.png"), height=8, width=16)
DEenrichRPlot(DC.interest, ident.1="RSV", ident.2="Untreated", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA")
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_UT_GOplot.png"), height=8, width=16)
DEenrichRPlot(DC.interest, ident.1="Flu", ident.2="Untreated", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA")
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_Flu_vs_UT_GOplot.png"), height=8, width=16)

# Customized enrichment plots
enriched <- DEenrichRPlot(DC.interest, ident.1="RSV", ident.2="Flu", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA", return.gene.list = TRUE)
names(enriched[[1]]) <- gsub("GO_Biological_Process_2021.","",names(enriched[[1]]))
p1 <- plotEnrich(enriched[[1]], showTerms=10, numChar=40, y="Ratio", orderBy="P.value", title="Upregulated in RSV vs Flu")
names(enriched[[2]]) <- gsub("GO_Biological_Process_2021.","",names(enriched[[2]]))
p2 <- plotEnrich(enriched[[2]], showTerms=10, numChar=40, y="Ratio", orderBy="P.value", title="Downregulated in RSV vs Flu")
p1 + p2
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_Flu_GOplot_new.png"), height=8, width=16)
write.table(enriched, file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_Flu_GOplot_new.csv"), sep=",", quote=FALSE, row.names=FALSE)

enriched <- DEenrichRPlot(DC.interest, ident.1="RSV", ident.2="Untreated", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA", return.gene.list = TRUE)
names(enriched[[1]]) <- gsub("GO_Biological_Process_2021.","",names(enriched[[1]]))
p1 <- plotEnrich(enriched[[1]], showTerms=10, numChar=40, y="Ratio", orderBy="P.value", title="Upregulated in RSV vs Untreated")
names(enriched[[2]]) <- gsub("GO_Biological_Process_2021.","",names(enriched[[2]]))
p2 <- plotEnrich(enriched[[2]], showTerms=10, numChar=40, y="Ratio", orderBy="P.value", title="Downregulated in RSV vs Untreated")
p1 + p2
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_UT_GOplot_new.png"), height=8, width=16)
write.table(enriched, file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_UT_GOplot_new.csv"), sep=",", quote=FALSE, row.names=FALSE)

enriched <- DEenrichRPlot(DC.interest, ident.1="Flu", ident.2="Untreated", test.use="MAST", max.genes=3000, enrich.database="GO_Biological_Process_2021", assay="RNA", return.gene.list = TRUE)
names(enriched[[1]]) <- gsub("GO_Biological_Process_2021.","",names(enriched[[1]]))
p1 <- plotEnrich(enriched[[1]], showTerms=10, numChar=40, y="Ratio", orderBy="P.value", title="Upregulated in Flu vs Untreated")
names(enriched[[2]]) <- gsub("GO_Biological_Process_2021.","",names(enriched[[2]]))
p2 <- plotEnrich(enriched[[2]], showTerms=10, numChar=40, y="Ratio", orderBy="P.value", title="Downregulated in Flu vs Untreated")
p1 + p2
ggsave(file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_Flu_vs_UT_GOplot_new.png"), height=8, width=16)
write.table(enriched, file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_Flu_vs_UT_GOplot_new.csv"), sep=",", quote=FALSE, row.names=FALSE)


### Making Gene lists
interest.markers <- FindMarkers(DC.interest, ident.1="RSV", ident.2="Flu", test.use="MAST", assay="RNA")
write.table(interest.markers, file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_Flu_cluster_DEgenes_RNA.csv"), sep=",", quote = FALSE)
RSV.Flu.markers.4h.RNA <- interest.markers
interest.markers <- FindMarkers(DC.interest, ident.1="RSV", ident.2="Untreated", test.use="MAST", assay="RNA")
write.table(interest.markers, file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_RSV_vs_UT_cluster_DEgenes_RNA.csv"), sep=",", quote = FALSE)
RSV.UT.markers.4h.RNA <- interest.markers
interest.markers <- FindMarkers(DC.interest, ident.1="Flu", ident.2="Untreated", test.use="MAST", assay="RNA")
write.table(interest.markers, file=paste0(outdir, "DC_0h_4h_noRSVHK_clust3_5_7_12_Flu_vs_UT_cluster_DEgenes_RNA.csv"), sep=",", quote = FALSE)
Flu.UT.markers.4h.RNA <- interest.markers
FeaturePlot(DC.interest, features="TANK", pt.size=1.5)
# ggsave(paste0(outdir,"0h_4h_CD1C_TANK_RNA.png"))
FeaturePlot(DC.interest, features=c("TANK","TBK1","IKBKE"), pt.size=1, label=TRUE, label.size= 6, cols=c("lightgrey","red"))
ggsave(paste0(outdir,"0h_4h_CD1C_TANK_TBK_IKK_RNA.png"))

