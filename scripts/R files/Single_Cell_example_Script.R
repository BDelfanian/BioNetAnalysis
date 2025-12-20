#########
#title: "Single_Cell"
#author: "Dimitrios Kyriakis and Alexander Skupin"
#date: "10/20/2022" and updated 10/10/2025
#############

# SingleCellExperiment
#Defines a S4 class for storing data from single-cell experiments. This includes specialized methods to store and retrieve spike-in information, dimensionality reduction coordinates and size factors for each cell, along with the usual metadata for genes and libraries.

# Seurat
#Seurat is an R package designed for QC, analysis, and exploration of single-cell RNA-seq data. Seurat aims to enable users to identify and interpret sources of heterogeneity from single-cell transcriptomic measurements, and to integrate diverse types of single-cell data.

# Scater
#scater provides tools for visualization of single-cell transcriptomic data. It is based on the SingleCellExperiment class (from the SingleCellExperiment package). 


# reset environment by

rm(list=ls())


library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(patchwork)


#To demonstrate the use of the various scater functions, we will load in the classic pbmc dataset:
#https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz



# Load the PBMC dataset - unzip data and put it into the folder of the script!
#getwd()
#datadir<- getwd()#"../data" #put the right paths here!
#pbmc.data <- Read10X(data.dir = datadir)
# remember to get function description by e.g.
# ?Read10X

# set here the right path to the folder from the zip folder
pbmc.data <- Read10X(data.dir = '/Users/alexander.skupin/projects/teaching/WS2024/MADS/SC_analysis_material1/filtered_gene_bc_matrices/pbmc_data_set/')

# Lets examine a few genes in the first thirty cells
dim(pbmc.data) # dimension of data frame
head(pbmc.data)
head(pbmc.data@Dimnames[[1]]) # gene names
head(pbmc.data@Dimnames[[2]]) # bar code names
pbmc.data[c("CD3D","TCL1A","MS4A1"), 1:30]
hist(pbmc.data["CD3D",])
hist(pbmc.data["CD3D",], breaks = 50)
hist(pbmc.data["CD3D",], breaks = 500)
hist(pbmc.data[122,1:2230])



# Diagnostic plots for quality control
#Quality control to remove damaged cells and poorly sequenced libraries is a common step in single-cell analysis pipelines. We will use some utilities from the scuttle package (conveniently loaded for us when we load scater) to compute the usual quality control metrics for this dataset.

#Metadata variables can be plotted against each other using the plotColData() function, as shown below. We expect to see an increasing number of detected genes with increasing total count. Each point represents a cell that is coloured according to its tissue of origin.




# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc.sce <- as.SingleCellExperiment(pbmc)
pbmc

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)
# Plot Histogram of mitopercentage
hist(pbmc@meta.data$percent.mt, breaks= 50)

#In the example below, we visualize QC metrics, and use these to filter cells.

#We filter cells that have unique feature counts over 2,500 or less than 200
#We filter cells that have >5% mitochondrial counts


# Visualize QC metrics as a violin plot
VlnPlot(pbmc , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# to save the plot for later comparison
vlnpl1 <- VlnPlot(pbmc , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
 
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot3 <- plot1 + plot2

# here quality filtering:
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# compare disrtibutions before and after QC:
VlnPlot(pbmc , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
vlnpl2 <- VlnPlot(pbmc , features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# for comparison
vlnpl1 # plot before QC
vlnpl2 # plot after QC

plot1 <- FeatureScatter(pbmc , feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc , feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

plot4 <- plot1 + plot2

# compare with befroe QC
plot3

# Normalizing the data
#After removing unwanted cells from the dataset, the next step is to normalize the data. 

## LogNormalize
#By default, we employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in pbmc[["RNA"]]@data.

## SCTransform
#sctransform to perform normalization and variance stabilization of scRNA-seq datasets. Finding also the most variable genes 
 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# the following is only needed if we normalize using the default seurat normalization
# Identification of highly variable features (feature selection)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
hist(pbmc.data["FTL",], breaks = 50)
# scale only the most variable genes


# plot variable features with and without labels
plot5 <- VariableFeaturePlot(pbmc)
plot5 #+ plot2

#scale the data
pbmc <- ScaleData(pbmc)
# or  scale all the detected genes
# all.genes <- rownames(pbmc)
# pbmc <- ScaleData(pbmc, features = all.genes)

# gene names are in the row of pbmc as  Features of pbmc
head(Features(pbmc))
head(rownames(pbmc))
# Cells are in columns
head(Cells(pbmc))
head(colnames(pbmc))

# get variable genes
head(VariableFeatures(pbmc))

# to get metadata from object:
head(pbmc[[]])

# to see count matrix
head(pbmc[["RNA"]]$counts)
GetAssayData(object = pbmc, assay = "RNA", slot = "counts")

# to get data out of the Seurat object:
head(FetchData(object = pbmc, vars = c("nFeature_RNA", "ISG15", "TNFRSF4"), layer = "counts"))

plot(as.data.frame(FetchData(object = pbmc, vars = c("FTL"), layer = "counts"))[,])
plot(as.data.frame(
  FetchData(object = pbmc, vars = c("FTL"), layer = "counts"))[,],as.data.frame(
    FetchData(object = pbmc, vars = c("GAPDH"), layer = "counts"))[,])

summary(as.data.frame(FetchData(object = pbmc, vars = c("FTL"), layer = "counts")))

#Perform linear dimensional reduction
#Next we perform PCA on the scaled data. By default, only the previously determined variable features are used as input, but can be defined using features argument if you wish to choose a different subset.


# Run dimensionality reduction on control dataset
pbmc <- RunPCA(pbmc,npcs = 30, verbose = FALSE)
ElbowPlot(pbmc,ndims = 30)
# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
DimPlot(pbmc, reduction = "pca")
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 2600, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)


# Dimensionality of dataset
pbmc <- JackStraw(pbmc, num.replicate = 50)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

JackStrawPlot(pbmc, dims = 1:10)

#have a look into pbmc
head(pbmc[[]])

# Cell-Cycle Scoring and Regression

#A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can segregate this list into markers of G2/M phase and markers of S phase


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(pbmc[[]])
DimPlot(pbmc,group.by = "Phase")


# Regress out cell cycle scores
#Some times we need to mitigate the effects of cell cycle heterogeneity in scRNA-seq data by calculating cell cycle phase scores based on canonical markers, and regressing these out of the data during pre-processing. 

# Regress out cell cycle scores can be done during data scaling -TAKES TIME !!!

pbmc <- ScaleData(pbmc, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(pbmc))
# Now, a PCA on the variable genes no longer returns components associated with cell cycle
pbmc <- RunPCA(pbmc, features = VariableFeatures(pbmc), nfeatures.print = 10)
DimPlot(pbmc,group.by = "Phase")

#Or as a better alternative regressing out the difference between the G2M and S phase scores.
#The previous procedure removes all signal associated with cell cycle. In some cases, we’ve found that this can negatively impact downstream analysis, particularly in differentiating processes (like murine hematopoiesis), where stem cells are quiescent and differentiated cells are proliferating (or vice versa).


pbmc$CC.Difference <- pbmc$S.Score - pbmc$G2M.Score
pbmc <- ScaleData(pbmc, vars.to.regress = "CC.Difference", features = rownames(pbmc)) # TAKES TIME!
pbmc <- RunPCA(pbmc, features =  VariableFeatures(pbmc),reduction.name = "pca_reg_ccs")
DimPlot(pbmc,reduction="pca_reg_ccs") |DimPlot(pbmc,reduction="pca")



#To cluster the cells, we next apply modularity optimization techniques such as the Louvain algorithm (default), to iteratively group cells together, with the goal of optimizing the standard modularity function.



library(dplyr)
pbmc <- FindNeighbors(pbmc, reduction = "pca", dims = 1:20, verbose = FALSE) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)

# check again pbmc
head(pbmc[[]])

# find clusters
pbmc <- FindClusters(pbmc, resolution = 0.5)

# check again pbmc
head(pbmc[[]])

DimPlot(pbmc,group.by = "seurat_clusters")
# Run non-linear dimensional reduction (UMAP/tSNE)
#UMAP is a stochastic algorithm – it makes use of randomness both to speed up approximation steps, and to aid in solving hard optimization problems. This means that different runs of UMAP can produce different results. So, to get the reproducible results, RunUmap function of seurat has as default parameter 
seed.use =42.
  

#library(ggplot2)
pbmc <- RunUMAP(pbmc,reduction = "pca", dims = 1:20, verbose = FALSE)
p1 <- DimPlot(pbmc, label = T, repel = T,reduction = "pca") #+ ggtitle("PCA")
p2 <- DimPlot(pbmc, label = T, repel = T,reduction = "umap") #+ ggtitle("UMAP")

# Change seed
pbmc <- RunUMAP(pbmc,reduction = "pca", dims = 1:20, verbose = FALSE,seed.use = 1234)
p3 <- DimPlot(pbmc, label = T, repel = T,reduction = "umap")# + ggtitle("Seed 1234")

p1 | p2 |p3 #remember your seed!!!


# Finding differentially expressed features (cluster biomarkers)


# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25,min.diff.pct =0.4)
cluster2.markers$Symbol <- rownames(cluster2.markers)
head(cluster2.markers, n = 5)
VlnPlot(pbmc, features = cluster2.markers$Symbol[1:2])
VlnPlot(pbmc, features = cluster2.markers$Symbol[1:6])
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
    group_by(cluster) %>%
    slice_max(n = 2, order_by = avg_log2FC)


#alternative
#cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)


# feature plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"))

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))


#DoHeatmap() generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

pbmc.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

# based on previous knowledge we can assigne cell types
new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()


# Gene-Gene correlation - towards  network description
#Here we will calculate the gene-gene correlation of the top 10 genes of each cluster.

M_gg = cor(t(as.data.frame(pbmc.data[top10$gene,])))
ComplexHeatmap::Heatmap(M_gg,column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6))

# Cell-Cell correlation - towards  network description
#Running cell cell correlation using all the cell will take too much time. This is why we used a subset of pbmc (pbmc_small). Another solution could be to create aggregations of similar cells.


M_cc = cor(as.data.frame(pbmc.data[top10$gene,c(1:400)]))
ComplexHeatmap::Heatmap(M_cc,column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 6))




# Useful webpages

#ANALYSIS OF SINGLE CELL RNA-SEQ DATA
#https://broadinstitute.github.io/2020_scWorkshop/index.html

#Seurat R toolkit for single cell genomics
#https://satijalab.org/seurat/

#Scanpy – Single-Cell Analysis in Python
#https://scanpy.readthedocs.io/en/stable/

