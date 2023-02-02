## Trajectory inference with Totem

library(Seurat)
library(Totem)
library(SingleCellExperiment)
library(dyndimred)

## Read in data
engel <- readRDS("Engel_data/Engel_count_data.rds")

## Create Seurat object and run preprocessing and clustering
seu_engel <- CreateSeuratObject(counts = engel)
seu_engel <- NormalizeData(seu_engel)
seu_engel <- FindVariableFeatures(seu_engel, selection.method = "vst", nfeatures = 2000)
#varfeat <- VariableFeatures(seu_engel)
seu_engel <- ScaleData(seu_engel)
seu_engel <- RunPCA(seu_engel)
ElbowPlot(seu_engel)
seu_engel <- FindNeighbors(seu_engel, dims=1:15)
seu_engel <- FindClusters(seu_engel, resolution = 0.5)
seu_engel <- RunUMAP(seu_engel, dims = 1:15)
png("Totem_figures/seurat_clusters.png")
DimPlot(seu_engel, reduction = "umap")
dev.off()
png("Totem_figures/seurat_clusters_ident.png")
DimPlot(seu_engel, reduction = "umap", group.by = "orig.ident")
dev.off()

## How does the clustering look like?

## Create SingleCellExperiment object and prepare for Totem
sce_engel <- SingleCellExperiment(assays = list(counts = seu_engel@assays$RNA@counts, logcounts = seu_engel@assays$RNA@data))
sce_engel <- PrepareTotem(sce_engel)

## Run dimensionality reduction
## Options for dimensionality reduction are LMDS, MDS, PCA, tSNE and UMAP
#sce_engel <- RunDimRed(sce_engel, dim.red.method = "lmds", dim.red.features = varfeat, dim.reduction.par.list = list(ndim=5))
sce_engel <- RunDimRed(sce_engel, dim.red.method = "lmds", dim.red.features = NULL, dim.reduction.par.list = list(ndim=5))

dim_red_engel <- dimred_mds(t(seu_engel@assays$RNA@data), ndim = 2)

## Run Totem trajectory inference
sce_engel <- RunClustering(sce_engel, k.range = 3:20, min.cluster.size = 5, N.clusterings = 10000)

## Visualise cell connectivity
png("Totem_figures/Engel_cellConnectivity.png")
VizCellConnectivity(sce_engel, viz.dim.red = dim_red_engel)
dev.off()

## How does the cell connectivity look like?

## Select top clusterings and visualise them
sce_engel <- SelectClusterings(sce_engel, selection.method = 1, selection.N.models = 10)
png("Totem_figures/Engel_top_trajectories.png", width=1000, height=1000)
VizMST(sce_engel, clustering.names = ReturnTrajNames(sce_engel), viz.dim.red = dim_red_engel)
dev.off()

## Are there many differences between the clusterings?

## Run smoothing
sce_engel <- RunSmoothing(sce_engel)
png("Totem_figures/Engel_smoothed_trajectories.png", width = 900, height = 1000)
VizSmoothedTraj(sce_engel, traj.names = ReturnTrajNames(sce_engel), viz.dim.red = dim_red_engel, plot.pseudotime = F)
dev.off()

## How do the directions of the clusterings vary?

## Visualise the cell identities
png("Totem_figures/Engel_traj_with_cell_idents.png")
VizClustering(sce_engel, clustering = as.character(seu_engel$orig.ident), viz.dim.red = dim_red_engel)
dev.off()

## Based on the cell identities and cluster directions, which trajectory will you choose?

## Select a trajectory and print its milestone network
selected_traj <- names(sce_engel@metadata[["totem"]][["selected.clustering"]])[6]
ReturnSmoothedTrajNetwork(sce_engel, clustering.name = selected_traj)

## Visualise the selected trajectory with and without pseudotime
png("Totem_figures/Engel_selected_trajectory.png")
VizSmoothedTraj(sce_engel, traj.names = selected_traj, viz.dim.red = dim_red_engel, plot.pseudotime = F)
dev.off()
png("Totem_figures/Engel_selected_trajectory_with_pseudotime.png")
VizSmoothedTraj(sce_engel, traj.names = selected_traj, viz.dim.red = dim_red_engel, plot.pseudotime = T)
dev.off()

## Change the root cluster, if needed
# sce_engel <- ChangeTrajRoot(sce_engel, traj.name = selected_traj, root.cluster = 3)

## What do you think, is changing the root cluster necessary?

## Visualising gene expression
## Rorc is a marker for NKT17 cells 
png("Totem_figures/Engel_Rorc_expression.png")
VizFeatureExpression(sce_engel, traj.name = selected_traj, feature.names = "Rorc", viz.dim.red = dim_red_engel, plot.traj = T)
dev.off()

## Is Rorc expression localised to NKT17 cells?

## save the data
save(sce_engel, file="RData/sce_engel.RData")
