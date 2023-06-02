# Run notes -------------------------------------------------------------

# Run analysis on single cell RNA data
# Sample data is from [PancDB](https://hpap.pmacs.upenn.edu/)

##NB: 
# Note on samples: based on cellrager QC, going to exclude the following (see T2DPancreas-V1.docx)
# Note on sample HPAP-108_FGC2390: files corrupted in PancDB. Excluding until file sizes exceed 512 bytes
#excluded.samples <- c("HPAP-027_78239", "HPAP087_FGC2276", "HPAP-090_FGC2390", "HPAP-092_FGC2390", "HPAP-093_FGC2332", "HPAP-093_FGC2390", "HPAP-093_2430", "HPAP-099_FGC2390", "HPAP-100_FGC2390", "HPAP-101_FGC2390", "HPAP-108_FGC2390")
# 2023.05.26 - due to limited sample number, comparing AA and EU in T2D cases ONLY

# Global parameters -------------------------------------------------------

rnaProject <- "PancT2D_AAvsEUonly-doPar"
regression.vars <- c("sequencerID", "SampleSex", "SampleAge")
cum.var.thresh <- 90
resolution <- 0.5
comp.type <- "biowulf" # one of macbookPro, biowulf
do.sctransform <- "each" # one of FALSE, each, pooled

## infrequently modified
do.doubletFinder <- TRUE
run.jackstraw <- FALSE
min.cells <- 3
min.features <- 200
doublet.var.thresh <- 90
predicted.doubletRate <- 0.05
excluded.samples <- c("HPAP-027_78239", "HPAP087_FGC2276", "HPAP-090_FGC2390", "HPAP-092_FGC2390", "HPAP-093_FGC2332", "HPAP-093_FGC2390", "HPAP-093_2430", "HPAP-099_FGC2390", "HPAP-100_FGC2390", "HPAP-101_FGC2390", "HPAP-108_FGC2390")


# Directories -------------------------------------------------------------

if(comp.type == "macbookPro"){
  rna.dir <- "/Users/heustonef/Desktop/Obesity/scRNA/"
  path_to_data <- "/Users/heustonef/Desktop/PancDB_data/scRNA_noBams"
  sourceable.functions <- list.files(path = "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/RFunctions", pattern = "*.R$", full.names = TRUE)
  metadata.location <- "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/"
} else if(comp.type == "biowulf"){
  rna.dir <- "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/"
  path_to_data <- "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/scRNA_transfer"
  sourceable.functions <- list.files(path = "/data/CRGGH/heustonef/hpapdata/RFunctions", pattern = "*.R", full.names = TRUE)
  metadata.location <- "/data/CRGGH/heustonef/hpapdata/"
}

# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(foreach)
library(doParallel)
cl <- makeCluster(future::availableCores(), outfile = "")

##load local functions
invisible(sapply(sourceable.functions, source))

# Load data ---------------------------------------------------------------

try(setwd(rna.dir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(rnaProject, "_sessionInfo.txt"))

## load data list
sc.data <- sapply(list.dirs(path = path_to_data, recursive = FALSE, full.names = TRUE), 
                  basename, 
                  USE.NAMES = TRUE)
sc.data <- sc.data[grepl(pattern = "^HPAP", sc.data)]

# Select data based on inclusion criteria
metadata <- read.table(file = paste0(metadata.location, "HPAPMetaData.txt"), header = TRUE, sep = "\t", row.names = 1)
metadata <- metadata %>%
  filter(grepl("Af|Cauc|Black", SampleEthnicity) &
           # filter(grepl("Af", SampleEthnicity) & 
           SimpDisease != "T1DM" & 
           !grepl("Fluidigm", scRNA_Platform) & 
           scRNA > 0)


# Exclude samples that failed CellRanger QC
for(i in sc.data){
  x <- strsplit(i, "_")[[1]][1]
  if(i %in% excluded.samples){
    sc.data <- sc.data[sc.data!=i]}
  if(!(x %in% rownames(metadata))){
    sc.data <- sc.data[sc.data!=i]
  }
}

# clean up metadata table to exclude data not relevant to patient
metadata <- metadata[,1:12]


registerDoParallel(cl)
object.list <- foreach(i=1:length(sc.data), .combine="c", .packages = 'Seurat') %dopar% {
  # for(i in 1:length(sc.data)){
  invisible(sapply(sourceable.functions, source))
  object.item <- Read10X_h5(paste0(names(sc.data)[i], "/outs/filtered_feature_bc_matrix.h5"))
  object.item <- CreateSeuratObject(object.item, 
                                    project = rnaProject, 
                                    min.cells = min.cells, 
                                    min.features = min.features)
  object.item$orig.ident <- sc.data[[i]]
  object.item <- AssignMetadata(metadata.df = metadata, seurat.object = object.item)
  object.item <- PercentageFeatureSet(object.item, pattern = "MT-", col.name = "percent.mt")
  
  object.item <- subset(object.item, 
                        subset = nFeature_RNA >= 200 &
                          nFeature_RNA <= 2500 &
                          percent.mt <= 5)
  if(length(Cells(object.item)) < 100){
    print(paste0("Excluding ", sc.data[[i]], " because cell count = ", length(Cells(object.item))))
    return(NULL)
  } else{
    print(paste("adding", sc.data[[i]], "to list"))
    return(object.item)
  }
}

seurat.object <- merge(object.list[[1]], y = object.list[2:length(object.list)], add.cell.ids = names(object.list))
seurat.object$DonorID <- sapply(seurat.object$orig.ident, sub, pattern = "_.*", replacement = "")
seurat.object$sequencerID <- seurat.object$orig.ident
seurat.object$sequencerID <- sapply(seurat.object$orig.ident, sub, pattern = ".*_", replacement = "")
saveRDS(seurat.object, file = paste0(rnaProject, "-rawMergedSeurat.Object.RDS"))

# QC ----------------------------------------------------------------------

##plot qc stats
VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
plot1<- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


# Normalize and scale data ----------------------------------------------------------


if(do.sctransform == FALSE){ # standard method
  print("Performing standard normalization and scaling")
  if(length(regression.vars) >1){
    print("HEY YOU! You're performing standard scaling on more than 1 regression variable. You should probably be doing SCTransform. Set `do.sctransform` to TRUE")
  }
  
  ##normalize
  seurat.object <- NormalizeData(seurat.object, normalization.method = "LogNormalize", scale.factor = 10000)
  
  ##find HVG
  seurat.object <- FindVariableFeatures(seurat.object, selection.method = "vst", nfeatures = 2000)
  top10hvg <- head(VariableFeatures(seurat.object), 10)
  plot1 <- VariableFeaturePlot(seurat.object)
  plot2	<- LabelPoints(plot = plot1, points = top10hvg, repel = TRUE)
  plot1 + plot2
  top10hvg
  
  ##scale (a linear transformation)
  all.genes <- rownames(seurat.object)
  
  seurat.object <- ScaleData(seurat.object, features = all.genes, vars.to.regress = regression.vars)
  
} else if(do.sctransform == "each"){
  print("Performing SCTransform")
  split.object <- SplitObject(seurat.object, split.by = "orig.ident")
  
  object.list <- foreach(i=1:length(split.object), .combine="c", .packages = c('Seurat', 'DoubletFinder', 'ggplot2', 'dplyr')) %dopar% {
    object.item <- SCTransform(split.object[[i]], assay = "RNA", return.only.var.genes = FALSE, vst.flavor = "v2")
    
    # RUN DOUBLETFINDER AFTER UMAPhttps://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_doublets.R
    object.item <- runDoubletFinder(object.item, sctransformed = TRUE, tot.var = doublet.var.thresh, predicted.doubletRate = predicted.doubletRate)
    object.item <- subset(object.item, subset = DF.classifications == "Singlet")
    if(length(Cells(object.item)) < 100){
      print(paste0("Excluding ", unique(object.item$orig.ident), " because corrected cell count = ", length(Cells(object.item))))
      return(NULL)
    }else{
      return(object.item)
    }
  }
  
  integration.features <- SelectIntegrationFeatures(object.list = object.list, verbose = TRUE, nfeatures = 3000)
  object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = integration.features, verbose = TRUE)
  integration.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = integration.features, normalization.method = "SCT", verbose = TRUE)
  seurat.object <- IntegrateData(anchorset = integration.anchors, verbose = TRUE, preserve.order = FALSE, normalization.method = "SCT")
  print(paste("made progress on", as.character(i)))
  
} else if(do.sctransform == "pooled") {
  seurat.object <- SCTransform(seurat.object, method = "glsGamPoi", vars.to.regress = regression.vars, verbose = TRUE, return.only.var.genes = FALSE, vst.flavor = "v2")
  seurat.object <- runDoubletFinder(seurat.object = seurat.object, sctransformed = TRUE, tot.var = doublet.var.thresh, predicted.doubletRate = predicted.doubletRate)
  seurat.object <- CellCycleScoring(seurat.object, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)
  seurat.object <- subset(seurat.object, subset = DF.classifications == "Singlet")
  
}else {
  print("Must set do.sctransform to one of: FALSE, each, pooled")
}
saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))



# Linear dimensional reduction --------------------------------------------


seurat.object <- RunPCA(seurat.object, features = VariableFeatures(object = seurat.object))

print(seurat.object[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(seurat.object, dims = 1:2, reduction = "pca")
DimPlot(seurat.object, reduction = "pca")
DimHeatmap(seurat.object, dims = 1:2, cells = 500, balanced = TRUE)

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))

# Determine dimensionality ------------------------------------------------

if(run.jackstraw == TRUE){
  seurat.object <- JackStraw(seurat.object, num.replicate = 100)
  seurat.object <- ScoreJackStraw(seurat.object, dims = 1:40)
  JackStrawPlot(seurat.object, dims = 1:40)
}

ElbowPlot(seurat.object)

# account for variance
tot.var <- percent.variance(seurat.object@reductions$pca@stdev, plot.var = FALSE, return.val = TRUE)
paste0("Num pcs for 80% variance: ", length(which(cumsum(tot.var) <= 80)))
paste0("Num pcs for 85% variance: ", length(which(cumsum(tot.var) <= 85)))
paste0("Num pcs for 90% variance: ", length(which(cumsum(tot.var) <= 90)))
paste0("Num pcs for 95% variance: ", length(which(cumsum(tot.var) <= 95)))

cluster.dims <- 0
if(cum.var.thresh > 0){
  cluster.dims <- length(which(cumsum(tot.var) <= cum.var.thresh))
}

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, ".RDS"))

# Louvain cluster ---------------------------------------------------------

##cluster cells
seurat.object <- FindNeighbors(seurat.object, dims = 1:cluster.dims)
seurat.object <- FindClusters(seurat.object, resolution = resolution)
seurat.object <- RunUMAP(seurat.object, dims = 1:cluster.dims)

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, "-", as.character(tot.var), "pctvar.RDS"))


##umap
DimPlot(seurat.object, reduction = "umap", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "Obesity", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "DonorID", cols = color.palette, label = T, label.size = 7, repel = T)


FeaturePlot(seurat.object, reduction = "umap", features = "BMI")



# Find cluster biomarkers -------------------------------------------------

seurat.object <- PrepSCTFindMarkers(seurat.object, assay = "SCT")
##find positively expressed markers for all clusters compared to all remaining clusters

markers.seurat.pos <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.pos <- markers.seurat.pos %>% 
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.pos, file = paste0(rnaProject, "-posmarkers-", as.character(tot.var), "pctvar.rds"))

markers.seurat.all <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.all %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.all, file = paste0(rnaProject, "-allmarkers-", as.character(tot.var), "pctvar.rds"))

##create workbook
markers.table <- openxlsx::createWorkbook()


##write positive markers to table
openxlsx::addWorksheet(markers.table, sheetName = "PosMarkers")
openxlsx::writeData(markers.table, sheet = "PosMarkers", x = markers.seurat.pos,startCol = 1, startRow = 1, colNames = TRUE)

##write all markers to table
openxlsx::addWorksheet(markers.table, sheetName = "AllMarkers")
openxlsx::writeData(markers.table, sheet = "AllMarkers", x = markers.seurat.all, startCol = 1, startRow = 1, colNames = TRUE)

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(rnaProject, "_seuratMarkers-", as.character(tot.var), "pctvar.xlsx"), overwrite = TRUE, returnValue = TRUE)


# Cell Cycle Scoring ------------------------------------------------------

if(!do.sctransform == FALSE){
  DefaultAssay(seurat.object) <- "SCT"
}
seurat.object <- CellCycleScoring(seurat.object, s.features = Seurat::cc.genes$s.genes, g2m.features = Seurat::cc.genes$g2m.genes)

# Local visualization -----------------------------------------------------
# seurat.object <- readRDS("Obesity_scRNA.RDS")

# did doublet removal happen?
colnames(seurat.object@meta.data)

DimPlot(seurat.object, cols = color.palette)
DimPlot(seurat.object, cols = color.palette, group.by = "seurat_clusters", split.by = "Obesity", pt.size = 0.4, ncol = 4)
FeaturePlot(seurat.object, features = "BMI", pt.size = 0.4, cols = c("blue", "red"))
FeaturePlot(seurat.object, features = "BMI", pt.size = 0.4, cols = c("blue", "red"), split.by = "SCT_snn_res.0.5", ncol = 4)

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, "-", as.character(tot.var), "pctvar.RDS"))

# Visualize after biowulf run ---------------------------------------------

# seurat.object <- readRDS("Obesity_scRNA-SCTRegression-NW-OB.RDS")
colnames(seurat.object@meta.data)

# cds <- readRDS("Obesity_scRNA-SCTRegression-NW-OB-monocle3CDS.RDS")

# Differential abundance testing ------------------------------------------
library(speckle)
library(limma)
library(ggplot2)


# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
# Run propeller testing for cell type proportion differences between the two 
# groups
propeller(clusters = Idents(seurat.object), sample = seurat.object$DonorID, group = seurat.object$Obesity, transform = "asin")

# Plot cell type proportions
plotCellTypeProps(clusters=seurat.object$Obesity, sample=seurat.object$SCT_snn_res.0.5)




# Trajectory analysis -----------------------------------------------------

library(SeuratWrappers)
library(monocle3)


cds <- as.cell_data_set(seurat.object)
fData(cds)$gene_short_name <- rownames(fData(cds))


recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions

list.cluster <- seurat.object@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat.object@reductions$umap@cell.embeddings
cds <- learn_graph(cds, use_partition = F)










# https://rdrr.io/github/satijalab/seurat-wrappers/f/docs/monocle3.Rmd
# seurat.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)

# From monocle3 tutorial
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds){
  cell_ids <- which(colData(cds)[, "SCT_snn_res.0.5"] %in% c(6, 7))
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
saveRDS(cds, file = paste0(rnaProject, "-monocle3CDS-90pctvar.RDS"))

