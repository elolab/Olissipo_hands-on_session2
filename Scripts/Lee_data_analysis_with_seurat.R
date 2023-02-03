## Hands-on session for Olissipo Modelling and Analysis of Single Cell Multiple Biological Omics
## Basic analysis workflow for Lee data

## Lee data is publicly available in GEO with accession GSE149689
## Contains PBMC data from 11 COVID patients, 5 flu patients and 4 healthy controls
## Download the data
## Link: https://bioinfoshare.utu.fi/DataTransfer/Olissipo_hands-on_session2/data
## Size: 468M
## Username: btk
## Password: WzfwWGUMo7

## First, read in data

library(Seurat)
library(dplyr)

data_lee <- Read10X("data")
## Read10X expects the files to be named barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz

## Create Seurat object
## You can name the project the way you like
seu_lee <- CreateSeuratObject(data_lee, project = "Lee COVID data")

## Have a look at the Seurat object structure
## Check how metadata looks like
head(seu_lee@meta.data)

## Calculate percentage of mitochondrial genes
seu_lee[["percent.mt"]] <- PercentageFeatureSet(seu_lee, pattern = "^MT-")

## QC plots
qc_vln <- VlnPlot(seu_lee, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
png("Figures/qc_violin_plots.png")
print(qc_vln)
dev.off()

## Filtering cells
## Filtering bad quality cells based on their properties, e.g., `nCount_RNA`, `nFeature_RNA` or/and `percent.mt`, 
## and low abundant genes will not be performed because the authors removed them already. 
## The code below is an example that you could use to filter the bad quality cells based on their properties. 

## seu_lee <- subset(seu_lee, subset = nFeature_RNA > 1000 & nFeature_RNA < 4000 & percent.mt < 10) 
## This would select cells expressing more than 1000 and less than 4000 distinct genes as well as cells expressing less than 10% mitochondrial genes. 

metadata <- readRDS("data/metadata_lee.rds")
seu_lee$cellID <- colnames(seu_lee)
## Include only high-quality cells, based on authors' filtering
seu_lee <- subset(seu_lee, subset = cellID %in% rownames(metadata))

## QC plots for filtered data
qc_vln <- VlnPlot(seu_lee, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) 
png("Figures/qc_violin_plots_filtered.png")
print(qc_vln)
dev.off()

## What is your opinion on the data quality?
## How many cells were filtered out?

## Add metadata
seu_lee <- AddMetaData(seu_lee, metadata)

## In the interest of time, we'll keep only the healthy controls and four severe COVID-19 patients
keep <- c("nCoV_1", "nCoV_3", "nCoV_4", "nCoV_7", "Normal_1", "Normal_2", "Normal_3", "Normal_4")
seu_lee <- subset(seu_lee, subset = individual %in% keep)

## Clustering without integration
seu_lee <- NormalizeData(seu_lee, normalization.method="LogNormalize", 
                     scale.factor=10000)
seu_lee <- FindVariableFeatures(seu_lee, selection.method="vst", 
                            nfeatures=2000)

## Identify and visualise the most highly variable features
top10 <- head(VariableFeatures(seu_lee), 10)
plot1 <- VariableFeaturePlot(seu_lee)
plot2 <- LabelPoints(plot1, points = top10, repel = T)
png("Figures/HVGs_top10.png")
plot2
dev.off()

## Scale data and run PCA
seu_lee <- ScaleData(seu_lee)
seu_lee <- RunPCA(seu_lee) 

## Visualise and select the number of PCs to use (the "dimensionality" of the data)
png("Figures/PCA_plot.png")
DimPlot(seu_lee, reduction = "pca")
dev.off()

png("Figures/PCA_plot_byGroup.png")
DimPlot(seu_lee, reduction = "pca", group.by = "group")
dev.off()

png("Figures/elbowPlot.png")
ElbowPlot(seu_lee, ndims = 30)
dev.off()

## How many PCs would you use? The authors used 15.

## Find neighbours and cluster the data
seu_lee <- FindNeighbors(seu_lee, dims=1:20)
seu_lee <- FindClusters(seu_lee)
seu_lee <- RunUMAP(seu_lee, dims=1:20)
seu_lee <- RunTSNE(seu_lee, dims = 1:20)

## Visualise clustering
png("Figures/unintegrated_UMAP.png")
DimPlot(seu_lee, reduction = "umap")
dev.off()

png("Figures/unintegrated_TSNE.png")
DimPlot(seu_lee, reduction = "tsne")
dev.off()

png("Figures/unintegrated_UMAP_byGroup.png")
DimPlot(seu_lee, reduction = "umap", group.by = "group")
dev.off()

png("Figures/unintegrated_UMAP_byIndividual.png")
DimPlot(seu_lee, reduction = "umap", group.by = "individual")
dev.off()

png("Figures/unintegrated_UMAP_byCelltype.png")
DimPlot(seu_lee, reduction = "umap", group.by = "cell_type")
dev.off()

## Are there batch effects in these data? Do you think integration is necessary?

## Integrating the data
seu_list <- SplitObject(seu_lee, split.by = "individual")

## Finding cell anchors and integrating the data
seu_list <- lapply(X=seu_list, FUN=function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures=2000)
})

## Using RPCA for integration
features <- SelectIntegrationFeatures(object.list = seu_list)
seu_list <- lapply(X=seu_list, FUN=function(x){
  x <- ScaleData(x, features = features, verbose = F)
  x <- RunPCA(x, features = features, verbose = F)
})
seu_anchors <- FindIntegrationAnchors(object.list = seu_list, anchor.features = features, reduction = "rpca")
seu_combined <- IntegrateData(anchorset = seu_anchors)

## Clustering the integrated data
DefaultAssay(seu_combined) <- "integrated"
seu_combined <- ScaleData(seu_combined)
seu_combined <- RunPCA(seu_combined)
seu_combined <- RunTSNE(seu_combined, reduction = "pca", dims = 1:20)
seu_combined <- RunUMAP(seu_combined, reduction = "pca", dims = 1:20)
seu_combined <- FindNeighbors(seu_combined, reduction = "pca", dims = 1:20)
seu_combined <- FindClusters(seu_combined, resolution = 0.5)

## Visualise integrated clustering
png("Figures/integrated_UMAP.png")
DimPlot(seu_combined, reduction = "umap")
dev.off()

png("Figures/integrated_TSNE.png")
DimPlot(seu_combined, reduction = "tsne")
dev.off()

png("Figures/integrated_UMAP_byGroup.png")
DimPlot(seu_combined, reduction = "umap", group.by = "group")
dev.off()

png("Figures/integrated_UMAP_byIndividual.png")
DimPlot(seu_combined, reduction = "umap", group.by = "individual")
dev.off()

png("Figures/integrated_UMAP_byCelltype.png")
DimPlot(seu_combined, reduction = "umap", group.by = "cell_type", label = T)
dev.off()

## Visualising marker genes from the paper
DefaultAssay(seu_combined) <- "RNA"
png("Figures/markerGenes_from_pub.png", width=1200, height=1000)
FeaturePlot(seu_combined, features = c("CD3E", "CD4", "CCR7", "CD8A", "NCAM1", "CD14", "FCGR3A", "NR4A1", "CD19", "FCER1A", "PPBP", "HBB"), order = T)
dev.off()

## Do you think the annotation corresponds well with the marker genes?

## Automatic annotation with SingleR

library(SingleR)
library(celldex)

## We'll use the human primary cell atlas data as the reference,
## other references are available through the celldex package
hpca <- HumanPrimaryCellAtlasData()

## SingleR needs the normalised expression matrix
DefaultAssay(seu_combined) <- "RNA"
seu_expm <- GetAssayData(seu_combined, slot="data")

pred_lee <- SingleR(test = seu_expm, ref = hpca, labels = hpca$label.main) ## takes about 4-5 min

## Each row of the output DataFrame contains prediction results for a single cell. Labels are shown before fine-tuning (first.labels), 
## after fine-tuning (labels) and after pruning (pruned.labels), along with the associated scores.

head(pred_lee)

## Visualise the annotation using a heatmap
## Ideally, each cell (i.e., column of the heatmap) should have one score that is obviously larger than the rest, 
## indicating that it is unambiguously assigned to a single label. A spread of similar scores for a given cell indicates 
## that the assignment is uncertain, though this may be acceptable if the uncertainty is distributed across similar cell types that cannot be easily resolved.

png("Figures/annotation_score_heatmap.png")
plotScoreHeatmap(pred_lee)
dev.off()

## Add the annotation to the Seurat object
seu_combined <- AddMetaData(seu_combined, pred_lee[,4], col.name = "SingleR_label")

## Visualise the SingleR annotation
png("Figures/integrated_UMAP_bySingleR_label.png")
DimPlot(seu_combined, reduction = "umap", group.by = "SingleR_label", label = T)
dev.off()

## How did the automatic annotation work? What are the main differences to the manual annotation?


## Differential expression analysis

## We will use the default method in Seurat (Wilcoxon) and pseudobulk method ROTS for DE analysis
## In the interest of time, we'll only perform the analysis on two cell types: Dendritic cells and Monocytes classical

## DE analysis using Wilcoxon

## Add the correct labels
seu_combined$Wilcoxon_group <- paste(seu_combined$cell_type, seu_combined$group, sep = "_")
Idents(seu_combined) <- "Wilcoxon_group"
DefaultAssay(seu_combined) <- "RNA"

## Run the DE analysis, takes about 5-6 mins
celltypes <- c("Monocyte_classical", "Dendritic_cells")
dge_wilcox <- list()

for (cell in celltypes){
  cat("\nPerforming single-cell DGE for", cell, "cell with Wilcox \n")
  dge_res <- FindMarkers(seu_combined, ident.1=paste0(cell, "_COVID-19"), ident.2=paste0(cell, "_Normal"),
                         logfc.threshold=0, min.pct=0.1, test.use="wilcox")
  dge_wilcox[[cell]] <- dge_res %>% mutate("gene"=row.names(dge_res)) %>%
  dplyr::select(all_of(c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")))
  cat("The no. of significant genes (FDR<0.05) found was:", sum(dge_wilcox[[cell]]$p_val_adj<0.05), "\n")
}


## DE analysis using ROTS

library(scuttle)
library(ROTS)
library(edgeR)

## Pseudobulks
## For pseudobulk analysis, we first need to aggregate the counts
seu_combined$ROTS_group <- paste(seu_combined$cell_type, seu_combined$individual, sep="_")
sce <- as.SingleCellExperiment(seu_combined)
ps <- list()
ps$counts <- aggregateAcrossCells(sce, sce$ROTS_group)
ps$counts <- SummarizedExperiment::assay(ps$counts)

## Normalization
ps$DGEList <- edgeR::DGEList(ps$counts, remove.zeros=T)
ps$DGEList <- edgeR::calcNormFactors(ps$DGEList)
ps$logcpm <- edgeR::cpm(ps$DGEList, normalized.lib.sizes=T, prior.count=1, log=T)

## DE analysis
dge_ps <- list()
for (cell in celltypes) {
  ctrl_samps <- grep(paste0(cell, "_Normal"), colnames(ps$logcpm), value=TRUE)
  trt_samps <- grep(paste0(cell, "_nCoV"), colnames(ps$logcpm), value=TRUE)
  group_comp <- c(rep(0, length(trt_samps)), rep(1, length(ctrl_samps)))
  cat("\nPerforming single-cell DGE for", cell, "cell with ROTS \n")
  dge_res <- ROTS(data=ps$logcpm[,c(trt_samps, ctrl_samps)], 
                    groups=group_comp, B=100, seed=1024)
  dge_ps[[cell]] <- data.frame("gene"=names(dge_res$logfc), 
                                         "logfc"=dge_res$logfc, 
                                         "pvalue"=dge_res$pvalue, 
                                         "FDR"=dge_res$FDR)
    cat("The no. of significant genes (FDR<0.05) found was:", sum(dge_res$FDR<0.05), "\n")
}


## How does the number of DE genes differ?

## Visualisation of DE genes
library(patchwork)

## Select top 10 DE genes from each comparison and draw violin plots
Idents(seu_combined) <- "cell_type"
top_rots <- list()
top_wilcox <- list()

for (cell in celltypes){
  print(cell)
  top_rots[[cell]] <- dge_ps[[cell]][order(dge_ps[[cell]]$FDR, abs(dge_ps[[cell]]$logfc), decreasing = c(F,T), method = "radix"),] %>% slice_head(n=10)
  top_wilcox[[cell]] <- dge_wilcox[[cell]][order(dge_wilcox[[cell]]$p_val_adj, abs(dge_wilcox[[cell]]$avg_log2FC), decreasing = c(F,T), method = "radix"),] %>% slice_head(n=10)
  plots_rots <- VlnPlot(seu_combined, features = top_rots[[cell]]$gene, split.by = "group", combine = F)
  png(paste0("Figures/top_rots_", cell, "_Vln.png"), width = 1000, height = 3000)
  print(wrap_plots(plots = plots_rots, ncol=1))
  dev.off()
  plots_wilcox <- VlnPlot(seu_combined, features = top_wilcox[[cell]]$gene, split.by = "group", combine = F)
  png(paste0("Figures/top_wilcox_", cell, "_Vln.png"), width = 1000, height = 3000)
  print(wrap_plots(plots = plots_wilcox, ncol=1))
  dev.off()
}

## How do the DE results differ between the two methods?

## Save the data
save(seu_combined, file="RData/seu_combined.RData")
save(dge_wilcox, file="RData/dge_wilcox.RData")
save(dge_ps, file="RData/dge_ps.RData")
save(top_wilcox, file="RData/top_wilcox.RData")
save(top_rots, file="RData/top_rots.RData")


