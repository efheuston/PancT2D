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
regression.vars <- c("sequencerID", "SampleSex", "SampleAge", "orig.ident")
cum.var.thresh <- 95
resolution <- 0.5
comp.type <- "macbookPro" # one of macbookPro, biowulf
do.sctransform <- "pooled" # one of FALSE, each, pooled

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
  rna.dir <- "/Users/heustonef/Desktop/PancDB_Data/ASHG2023/"
  # path_to_data <- "/Users/heustonef/Desktop/PancDB_data/scRNA_noBams"
  sourceable.functions <- list.files(path = "/Users/heustonef/OneDrive/SingleCellMetaAnalysis/GitRepositories/RFunctions/", pattern = "*.R$", full.names = TRUE)
  metadata.location <- "/Users/heustonef/OneDrive/SingleCellMetaAnalysis/"
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
# library(foreach)
# library(doParallel)
# cl <- makeCluster(future::availableCores(), outfile = "")

##load local functions
invisible(sapply(sourceable.functions, source))
try(setwd(rna.dir), silent = TRUE)

# Load data ---------------------------------------------------------------

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
           SimpDisease != "NoDM" &
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
stopCluster(cl)

seurat.object <- merge(object.list[[1]], y = object.list[2:length(object.list)], add.cell.ids = names(object.list))
seurat.object$DonorID <- sapply(seurat.object$orig.ident, sub, pattern = "_.*", replacement = "")
seurat.object$sequencerID <- seurat.object$orig.ident
seurat.object$sequencerID <- sapply(seurat.object$orig.ident, sub, pattern = ".*_", replacement = "")
saveRDS(seurat.object, file = paste0(rnaProject, "-rawMergedSeurat.Object.RDS"))

remove(object.list)


# QC ----------------------------------------------------------------------

##plot qc stats
#VlnPlot(seurat.object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
#plot1<- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(seurat.object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2


# Normalize and scale data ----------------------------------------------------------

cl <- makeCluster(future::availableCores(), outfile = "")
registerDoParallel(cl)

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
  save.image(file = paste0(rna.dir, "/", rnaProject, ".RData"))
  integration.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = integration.features, normalization.method = "SCT", verbose = TRUE)
  print("Starting IntegrateData function")
  save.image(file = paste0(rna.dir, "/", rnaProject, ".RData"))
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
stopCluster(cl)
remove(object.list)
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

levels(seurat.object$seurat_clusters)


saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, "-", as.character(cum.var.thresh), "pctvar.RDS"))


##umap
DimPlot(seurat.object, reduction = "umap", cols = color.palette, label = T, label.size = 7, repel = T)

DimPlot(seurat.object, reduction = "umap", group.by = "SampleEthnicity", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "SampleEthnicity", cols = color.palette, label = F, label.size = 7, repel = T, split.by = "integrated_snn_res.0.5", ncol = 6)
DimPlot(seurat.object, reduction = "umap", group.by = "orig.ident", cols = color.palette, label = T, label.size = 7, repel = T)
DimPlot(seurat.object, reduction = "umap", group.by = "DonorID", cols = color.palette, label = T, label.size = 7, repel = T)

ethn.per.cluster <- as.data.frame.matrix(table(seurat.object$integrated_snn_res.0.5, seurat.object$SampleEthnicity))
ethn.per.cluster <- round((ethn.per.cluster/rowSums(ethn.per.cluster)),2)
ethn.per.cluster
head(seurat.object@meta.data)
seurat.object@meta.data$cluster.pct.AF <- seurat.object@meta.data$integrated_snn_res.0.5
seurat.object@meta.data$cluster.pct.EU <- seurat.object@meta.data$integrated_snn_res.0.5
# seurat.object@meta.data <- select(seurat.object@meta.data,-c("cluster.pct.AF", "cluster.pct.EU"))
# 

temp <- seurat.object@meta.data
head(temp)
conditions_df <- data.frame(
  Condition = rownames(ethn.per.cluster),  # Conditions to check
  Replacement.AF = ethn.per.cluster[,"African american/Black"],  # Corresponding replacements
  Replacement.EU = ethn.per.cluster[,"Caucasian"]  # Corresponding replacements
)
# conditions_df$Condition <- as.numeric(conditions_df$Condition)
conditions_df$Replacement.AF <- as.character(conditions_df$Replacement.AF)
conditions_df$Replacement.EU <- as.character(conditions_df$Replacement.EU)

temp <- temp %>%
  mutate_at(vars(cluster.pct.AF),
            funs(case_when(
              . %in% conditions_df$Condition ~ conditions_df$Replacement.AF[match(., conditions_df$Condition)],
              TRUE ~ .
            ))
  )

temp <- temp %>%
  mutate_at(vars(cluster.pct.EU),
            funs(case_when(
              . %in% conditions_df$Condition ~ conditions_df$Replacement.EU[match(., conditions_df$Condition)],
              TRUE ~ .
            ))
  )
temp$cluster.pct.AF<- as.numeric(temp$cluster.pct.AF)
temp$cluster.pct.EU<- as.numeric(temp$cluster.pct.EU)
seurat.object@meta.data <- temp
head(seurat.object@meta.data)
# FeaturePlot(seurat.object, features = c("cluster.pct.AF", "cluster.pct.EU"), blend = TRUE)



FeaturePlot(seurat.object, features = "cluster.pct.AF", cols = c("red", "blue"))





FeaturePlot(seurat.object, reduction = "umap", features = "BMI")



# Find cluster biomarkers -------------------------------------------------

seurat.object <- PrepSCTFindMarkers(seurat.object, assay = "SCT")
##find positively expressed markers for all clusters compared to all remaining clusters

markers.seurat.pos <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.pos <- markers.seurat.pos %>% 
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.pos, file = paste0(rnaProject, "-posmarkers-", as.character(cum.var.thresh), "pctvar.rds"))

markers.seurat.all <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.all %>%
  group_by(cluster) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.all, file = paste0(rnaProject, "-allmarkers-", as.character(cum.var.thresh), "pctvar.rds"))

##create workbook
markers.table <- openxlsx::createWorkbook()


##write positive markers to table
openxlsx::addWorksheet(markers.table, sheetName = "PosMarkers")
openxlsx::writeData(markers.table, sheet = "PosMarkers", x = markers.seurat.pos,startCol = 1, startRow = 1, colNames = TRUE)

##write all markers to table
openxlsx::addWorksheet(markers.table, sheetName = "AllMarkers")
openxlsx::writeData(markers.table, sheet = "AllMarkers", x = markers.seurat.all, startCol = 1, startRow = 1, colNames = TRUE)

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(rnaProject, "_seuratMarkers-", as.character(cum.var.thresh), "pctvar.xlsx"), overwrite = TRUE, returnValue = TRUE)


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

saveRDS(seurat.object, file = paste0(rna.dir, "/", rnaProject, "-", as.character(cum.var.thresh), "pctvar.RDS"))

# Visualize after biowulf run ---------------------------------------------

colnames(seurat.object@meta.data)


# Calc percent_ribo
seurat.object <- PercentageFeatureSet(seurat.object, "^RP[SL]", col.name = "percent_ribo")
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
seurat.object <- PercentageFeatureSet(seurat.object, "^HB[^(P)]", col.name = "percent_hb")

seurat.object <- PercentageFeatureSet(seurat.object, "PECAM1|PF4", col.name = "percent_plat")


feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo", "percent_hb", "percent_plat")
VlnPlot(seurat.object, group.by = "integrated_snn_res.0.5", features = feats, pt.size = 0, ncol = 3, cols = color.palette) +
  NoLegend()


FeaturePlot(seurat.object, features = c("PDX1", "NKX6-1", "SOX2", "POU5F1", "CD34", "CD3G", "CD19", "CD14"))


# subset for presentation -------------------------------------------------
seurat.object <- readRDS("/Users/heustonef/Desktop/PancDB_Data/PancT2D_scRNA/PancT2D_AAvsEUonly-doPar-SCTeach-95pctvar.RDS")
seurat.subset <- subset(seurat.object, idents = c("1", "15"), invert = TRUE)
levels(seurat.subset)
cluster.ids.subset <- c(
  "Exocrine1", #SCT0
  "Alpha2", #SCT2
  "Exocrine2", #SCT3--REG1A
  "Beta1", #SCT4
  "Epithelial", #SCT5--KRT18
  "Endothelial", #SCT6
  "Alpha3", #SCT7
  "Beta2", #SCT8
  "Immune1", #SCT9--IGFBP7
  "Alpha4", #SCT10
  "Immune3", #SCT11--IGFBP7/NEAT1
  "Alpha1", #SCT12
  "Immune4", #SCT13--NEAT1
  "Mast Cells", #SCT14--TPSB2
  "Macrophages", #SCT16
  "Immune2" #SCT17
)
names(cluster.ids.subset) <- levels(seurat.subset)
seurat.subset <- RenameIdents(seurat.subset, cluster.ids.subset)
seurat.subset$cell.ids <- seurat.subset@active.ident

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent_ribo", "percent_hb", "percent_plat")
VlnPlot(subset(seurat.subset, idents = c("Beta1", "Beta2")), group.by = "integrated_snn_res.0.5", features = feats, pt.size = 0, ncol = 3, cols = color.palette) +
  NoLegend()

# features.dotplot <- c("PNLIP", "PPY", "GCG", "REG1A", "INS", "KRT18", "PECAM1", "TTR", "IGFBP7", "C11orf96", "NEAT1", "TPSB2", "SST", "APOE", "CDH19")
features.dotplot <- c(
  "TTR", 
  "GCG", #alpha
  "INS", #beta
  "KRT18", #acinar/epithilial
  "TPSB2",
  "MALAT1",
  "NEAT1",
  "IGFBP7",
  "C11orf96", 
  "CDH19",
  "APOE",
  "PECAM1", 
  "REG1A", 
  "PNLIP")
DotPlot(seurat.subset, features = features.dotplot, cluster.idents = TRUE, col.min = 0.01, dot.min = .51, dot.scale = 12)

png(filename = "features-dotMarkers.png", width=1000, height=800, bg = "transparent")
DotPlot(seurat.subset, features = features.dotplot, cluster.idents = TRUE, col.min = 0.01, dot.min = .51, dot.scale = 12)
dev.off()

classic.dotplot <- c(
  "GCG", #alpha
  "INS", #beta
  "PPY", #gamma
  "SST", #delta
  "AMY2A",
  "KRT19",
  "MAFB",
  # "CD4", 
  # "CD8A", 
  # "CD44",
  "VWF",
  "PECAM1", 
  # "ITGAM", 
  # "IL2RA",
  "APOE") #myeloid

clustermap <- DotPlot(seurat.subset, features = classic.dotplot, cluster.idents = TRUE, col.min = 0.01) & 
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))
plot(clustermap)
ggsave("clusters-dotMarkers.png", plot= clustermap, width=10, height=10, dpi=300, bg = "transparent")


DotPlot(seurat.subset, features = features.dotplot, cluster.idents = TRUE, col.min = 0.01, cols = c("gray94", "blue2"))
png(filename = "dotplot.png", height = 800, width = 800, bg = "transparent")
DotPlot(seurat.subset, features = features.dotplot, cluster.idents = TRUE, col.min = 0.01, cols = c("gray94", "blue2"))
dev.off()

DimPlot(seurat.subset, group.by = "integrated_snn_res.0.5", cols = color.palette, shuffle = TRUE, label = FALSE, label.size = 5) + NoLegend()

temp.palette <- color.palette
temp.palette[c(1, 3, 12)] <- c("maroon", "darkturquoise", "red")
DimPlot(seurat.subset, group.by = "integrated_snn_res.0.5", shuffle = T, cols = color.palette, pt.size = 3, label = T)

png(filename = "seurat-clusteredByident.png", height = 1200, width = 1200, bg = "transparent")
DimPlot(seurat.subset, group.by = "integrated_snn_res.0.5", cols = color.palette, shuffle = TRUE, label = FALSE, label.size = 5) + NoLegend() &
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))
dev.off()
png(filename = "seurat-clusteredByethnicity.png", height = 1200, width = 1200, bg = "transparent")
DimPlot(seurat.subset, group.by = "SampleEthnicity", cols = color.palette, shuffle = TRUE, label = FALSE, label.size = 5) + NoLegend() &
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))
dev.off()



DimPlot(subset(seurat.subset, idents = c("Exocrine1", "Exocrine2", "Beta1")), group.by = "SampleEthnicity", split.by = "cell.ids", ncol = 2, shuffle = T, cols = color.palette, alpha = 0.4)
# png(filename = "clusters-ethnicity.png", height = 800, width = 800, bg = "transparent")
clustermap <- DimPlot(subset(seurat.subset, idents = c("Exocrine1", "Exocrine2", "Beta1")), group.by = "SampleEthnicity", split.by = "cell.ids", ncol = 2, shuffle = T, cols = color.palette, alpha = 0.4) & 
  theme(legend.background = element_rect(fill = "transparent"),
                   legend.box.background = element_rect(fill = "transparent"),
                   panel.background = element_rect(fill = "transparent"),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.background = element_rect(fill = "transparent",
                                                  color = NA))
ggsave("clusters-ethnicity.png", plot= clustermap, width=10, height=10, dpi=300, bg = "transparent")


FeaturePlot(seurat.subset, features = features.dotplot, ncol = 4)


abstract.dotplot <- c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1")
abstract.plot <- FeaturePlot(seurat.subset, features = abstract.dotplot, ncol = 2) & 
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))
ggsave("Abstract-featureMarkers.png", plot= abstract.plot, width=10, height=10, dpi=300, bg = "transparent")


VlnPlot(seurat.subset, features = c("IAPP", "PPP1R1A", "PLCG2"), cols = color.palette, idents = c("Beta1", "Beta2"))
VlnPlot(seurat.subset, features = "percent.mt", cols = color.palette, idents = c("Beta1", "Beta2"))
VlnPlot(seurat.subset, features = "G2M.Score", cols = color.palette, idents = c("Beta1", "Beta2"))

VlnPlot(subset(seurat.subset, idents = c("Exocrine1", "Exocrine2")), features = c("HNF4A", "HNF4G", "TCF7L2"), split.by = "SampleEthnicity", split.plot = TRUE)
png(filename = "Beta1VBeta2-dotplot.png", height = 800, width = 1200)
VlnPlot(subset(seurat.subset, idents = c("Beta1", "Beta2")), features = c("IAPP", "PPP1R1A", "percent.mt"), split.plot = FALSE, cols = c("green", "darkviolet"))
dev.off()

# Define Subpopulation differences -----------------------------------------------
levels(seurat.subset@active.ident)

seurat.subset <- PrepSCTFindMarkers(seurat.subset, assay = "SCT")


markers.exocrine <- FindMarkers(seurat.subset, ident.1 = "Exocrine1", ident.2 = "Exocrine2")
markers.exocrine <- markers.exocrine %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.exocrine, file = paste0(rnaProject, "-Exocrine-markers.rds"))


markers.beta <- FindMarkers(seurat.subset, ident.1 = "Beta1", ident.2 = "Beta2")
markers.beta <- markers.beta %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.beta, file = paste0(rnaProject, "-Beta-markers.rds"))


markers.immune1 <- FindMarkers(seurat.subset, ident.1 = "Immune1", ident.2 = c("Immune2", "Immune3", "Immune4"))
markers.immune1 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.immune1, file = paste0(rnaProject, "-Immune1-markers.rds"))

markers.immune2 <- FindMarkers(seurat.subset, ident.1 = "Immune2", ident.2 = c("Immune1", "Immune3", "Immune4"))
markers.immune2 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.immune2, file = paste0(rnaProject, "-Immune2-markers.rds"))

markers.immune3 <- FindMarkers(seurat.subset, ident.1 = "Immune3", ident.2 = c("Immune1", "Immune2", "Immune4"))
markers.immune3 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.immune3, file = paste0(rnaProject, "-Immune3-markers.rds"))

markers.immune4 <- FindMarkers(seurat.subset, ident.1 = "Immune4", ident.2 = c("Immune1", "Immune3", "Immune2"))
markers.immune4 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.immune4, file = paste0(rnaProject, "-Immune4-markers.rds"))


markers.Alpha1 <- FindMarkers(seurat.subset, ident.1 = "Alpha1", ident.2 = c("Alpha2", "Alpha3", "Alpha4"))
markers.Alpha1 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.Alpha1, file = paste0(rnaProject, "-Alpha1-markers.rds"))

markers.Alpha2 <- FindMarkers(seurat.subset, ident.1 = "Alpha2", ident.2 = c("Alpha1", "Alpha3", "Alpha4"))
markers.Alpha2 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.Alpha2, file = paste0(rnaProject, "-Alpha2-markers.rds"))

markers.Alpha3 <- FindMarkers(seurat.subset, ident.1 = "Alpha3", ident.2 = c("Alpha1", "Alpha2", "Alpha4"))
markers.Alpha3 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.Alpha3, file = paste0(rnaProject, "-Alpha3-markers.rds"))

markers.Alpha4 <- FindMarkers(seurat.subset, ident.1 = "Alpha4", ident.2 = c("Alpha1", "Alpha3", "Alpha2"))
markers.Alpha4 %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.Alpha4, file = paste0(rnaProject, "-Alpha4-markers.rds"))




##create workbook
markers.table <- openxlsx::createWorkbook()


##write positive markers to table
openxlsx::addWorksheet(markers.table, sheetName = "Exocrine")
openxlsx::writeData(markers.table, sheet = "Exocrine", x = markers.exocrine,startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Beta")
openxlsx::writeData(markers.table, sheet = "Beta", x = markers.beta, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)

openxlsx::addWorksheet(markers.table, sheetName = "Immune1")
openxlsx::writeData(markers.table, sheet = "Immune1", x = markers.immune1, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Immune2")
openxlsx::writeData(markers.table, sheet = "Immune2", x = markers.immune2, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Immune3")
openxlsx::writeData(markers.table, sheet = "Immune3", x = markers.immune3, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Immune4")
openxlsx::writeData(markers.table, sheet = "Immune4", x = markers.immune4, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)


openxlsx::addWorksheet(markers.table, sheetName = "Alpha1")
openxlsx::writeData(markers.table, sheet = "Alpha1", x = markers.Alpha1, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Alpha2")
openxlsx::writeData(markers.table, sheet = "Alpha2", x = markers.Alpha2, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Alpha3")
openxlsx::writeData(markers.table, sheet = "Alpha3", x = markers.Alpha3, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::addWorksheet(markers.table, sheetName = "Alpha4")
openxlsx::writeData(markers.table, sheet = "Alpha4", x = markers.Alpha4, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)


##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(rnaProject, "_Pop_diff-Markers.xlsx"), overwrite = TRUE, returnValue = TRUE)



# Explore Beta1 vs Beta2 --------------------------------------------------

markers.beta <- readRDS(paste0(rnaProject, "-Beta-markers.rds"))
head(markers.beta)

beta.subset <- subset(seurat.subset, idents = c("Beta1", "Beta2"))
DimPlot(beta.subset)
VlnPlot(beta.subset, features = c("GCK", "HK1", "SLC2A2", "CD9", "ST8SIA1"))

beta.palette <- c("green", "darkviolet")

# Explore Exocrine --------------------------------------------------------

exocrine.subset <- subset(seurat.subset, idents = c("Exocrine1", "Exocrine2"))

Idents(exocrine.subset) <- "SampleEthnicity"
levels(exocrine.subset)
png(filename = "seurat-exocrineAAvsEA.png", height = 1200, width = 1200, bg = "transparent")
DimPlot(exocrine.subset, group.by = "SampleEthnicity", cols = c("red", "blue"), shuffle = TRUE, label = FALSE, label.size = 5) + NoLegend() &
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))
dev.off()

png(filename = "seurat-Exocrine1vsExocrine2.png", height = 1200, width = 1200, bg = "transparent")
DimPlot(exocrine.subset, group.by = "cell.ids", cols = c("maroon", "darkturquoise"), shuffle = TRUE, label = FALSE, label.size = 5) + NoLegend() &
  theme(legend.background = element_rect(fill = "transparent"),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent",
                                       color = NA))
dev.off()


exocrine.subset <- PrepSCTFindMarkers(exocrine.subset)
exocrine.markers <- FindMarkers(exocrine.subset, ident.1 = "African american/Black", ident.2 = "Caucasian")
exocrine.markers <- exocrine.markers %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(exocrine.markers, file = paste0(rnaProject, "-Exocrine-markers-AAvsEA.rds"))

exocrine.features <- c("MTRNR2L12", "AMY2B", "REG3A", "INS", "TTR", "GCG", "PPY", "SST", "SCG5")
features.heatmap <- exocrine.markers %>%
  filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 2) %>%
  rownames()

DoHeatmap(subset(exocrine.subset, downsample = 1000), features = features.heatmap, size = 3, group.colors = c("red", "blue"))

classic.exocrine.markers <- c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1", "AMY2A")
RidgePlot(exocrine.subset, features = c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1", "AMY2A"), ncol = 2, group.by = "cell.ids")
RidgePlot(exocrine.subset, features = features.heatmap[1:8], ncol = 2, group.by = "cell.ids", cols = c("red", "tan1"), sort = TRUE)


png(filename = "Exocrines-AAvsEA_abstractFeatures_AAvsEA-ridge.png", height = 600, width = 1200, bg = "transparent")
RidgePlot(exocrine.subset, features = exocrine.features, ncol = 3, group.by = "SampleEthnicity", cols = c("red", "blue"))
dev.off()

png(filename = "Exocrines-AAvsEA_abstractFeatures_AAvsEA-VLN.png", height = 800, width = 1000, bg = "transparent")
VlnPlot(exocrine.subset, features = c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1","AMY2A"), split.by = "SampleEthnicity", split.plot = FALSE, pt.size = 0.1, combine = TRUE, cols = c("red", "blue"))
dev.off()

png(filename = "Exocrines_abstractFeatures_1vs2-VLN.png", height = 800, width = 1000, bg = "transparent")
VlnPlot(exocrine.subset, features = c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1","AMY2A"), split.by = "cell.ids", split.plot = FALSE, pt.size = 0.1, combine = TRUE, cols = c("maroon", "darkturquoise")) &
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))
dev.off()



png(filename = "Exocrines_features_1vs2-VLN.png", height = 800, width = 1000, bg = "transparent")
VlnPlot(exocrine.subset, features = exocrine.features, split.by = "cell.ids", split.plot = FALSE, pt.size = 0.1, combine = TRUE, cols = c("maroon", "darkturquoise")) &
  theme(legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.box.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))
dev.off()


exocrine.features <- c("MTRNR2L12", "AMY2B", "REG3A", "INS", "TTR", "GCG", "PPY", "SST", "SCG5")
for(i in exocrine.features){
  gene.name <- i
  print(i)
  p <- VlnPlot(exocrine.subset, features = gene.name, split.by = "cell.ids", split.plot = FALSE, pt.size = 0.1, combine = TRUE, cols = c("maroon", "darkturquoise")) + NoLegend() &
    theme(legend.background = element_rect(fill = "transparent"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          legend.box.background = element_rect(fill = "transparent"),
          panel.background = element_rect(fill = "transparent"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA))
  png(filename = paste0("Exocrines_features_1vs2-VLN-", gene.name, ".png"), height = 400, width = 800, bg = "transparent")
  plot(p)
  dev.off()
}

for(i in exocrine.features){
  gene.name <- i
  print(i)
  p <- RidgePlot(exocrine.subset, features = gene.name, group.by = "SampleEthnicity", cols = c("red", "blue"))+ NoLegend() &
    theme(legend.background = element_rect(fill = "transparent"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA))
  png(filename = paste0("Exocrines_classics_1vs2-ridge-", gene.name, ".png"), height = 400, width = 800, bg = "transparent")
  plot(p)
  dev.off()
}



RidgePlot(exocrine.subset, features = "GCG", ncol = 3, group.by = "SampleEthnicity", cols = c("red", "blue"))+ NoLegend() &
  theme(
    legend.background = element_rect(fill = "transparent"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        # panel.background = element_rect(fill = "transparent"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA))







exocrine.heatmap <- exocrine.markers %>%
  mutate_at("Exocrine", ~(scale(.) %>% as.vector))


exocrine.hm <- DoHeatmap(exocrine.subset, features = exocrine.features, slot = "data")
exocrine.hm <- exocrine.hm$data
exocrine.hm <- na.omit(exocrine.hm)
exocrine.hm <- data.frame(exocrine.hm)
exocrine.hm$cell.id <- exocrine.hm$Cell
for(i in nrow(exocrine.hm)){
  cell.name <- exocrine.hm[i, "Cell"]
  exocrine.hm[cell.name, "cell.id"] <- exocrine.subset@meta.data[cell.name, "cell.ids"]
}
head(exocrine.hm)
exocrine.hm[33057, ]

markers.table <- openxlsx::createWorkbook()
openxlsx::addWorksheet(markers.table, sheetName = "ExocrineAAvsEA")
openxlsx::writeData(markers.table, sheet = "ExocrineAAvsEA", x = exocrine.markers, startCol = 1, startRow = 1, colNames = TRUE, rowNames = TRUE)
openxlsx::saveWorkbook(wb = markers.table, file = paste0(rnaProject, "_ExocrineMarkersAAvsEA.xlsx"), overwrite = FALSE, returnValue = TRUE)



# Downsampling ------------------------------------------------------------

seurat.object <- readRDS("/Users/heustonef/Desktop/PancDB_Data/PancT2D_scRNA/PancT2D_AAvsEUonly-doPar-SCTeach-95pctvar.RDS")
seurat.subset <- subset(seurat.object, idents = c("1", "15"), invert = TRUE)
levels(seurat.subset)
cluster.ids.subset <- c(
  "Exocrine1", #SCT0
  "Alpha2", #SCT2
  "Exocrine2", #SCT3--REG1A
  "Beta1", #SCT4
  "Epithelial", #SCT5--KRT18
  "Endothelial", #SCT6
  "Alpha3", #SCT7
  "Beta2", #SCT8
  "Immune1", #SCT9--IGFBP7
  "Alpha4", #SCT10
  "Immune3", #SCT11--IGFBP7/NEAT1
  "Alpha1", #SCT12
  "Immune4", #SCT13--NEAT1
  "Mast Cells", #SCT14--TPSB2
  "Macrophages", #SCT16
  "Immune2" #SCT17
)
names(cluster.ids.subset) <- levels(seurat.subset)
seurat.subset <- RenameIdents(seurat.subset, cluster.ids.subset)
seurat.subset$cell.ids <- seurat.subset@active.ident
seurat.subset$cell.ids <- factor(seurat.subset$cell.ids)

DimPlot(seurat.subset)
ea.subset <- subset(seurat.subset, subset = SampleEthnicity == "Caucasian")
Idents(ea.subset) <- "SampleEthnicity" # subset downsampling works per ident group
aa.subset <- subset(seurat.subset, subset = SampleEthnicity == "African american/Black")
Idents(aa.subset) <- "SampleEthnicity"

ea.downsampled <- subset(x = ea.subset, downsample = ncol(aa.subset))
# downsample.df <- data.frame(matrix(ncol=length(unique(ea.subset$cell.ids)),nrow=0, dimnames=list(NULL, unique(ea.subset$cell.ids))))


# iterate.downsampling <- function(object.seurat, downsample.size = 100, n.iterations=5, metadata.col = "cell.ids", class.ident = NULL){
#   if(is.null(metadata.col)){
#     print("Please specify metadata.column for dataframe")
#   }
#   if(!is.null(class.ident)){
#     print("This feature hasn't been finished yet. Using active.ident for downsampling.")
#   }
#   metadata.levels <- levels(object.seurat@meta.data[,metadata.col])
#   downsample.df <- data.frame(matrix(ncol=length(metadata.levels),nrow=0))
#   colnames(downsample.df) <- metadata.levels
#   for(i in 1:n.iterations){
#     # print(i)
#     downsampled.object <- subset(object.seurat, downsample = downsample.size)
#     for(cell.id in metadata.levels){
#       # print(cell.id)
#       downsample.df[i, cell.id] <- nrow(subset(x = downsampled.object, subset = cell.ids == cell.id))
#     }
#     print(head(downsample.df))
#   } 
#   return(downsample.df)
# }
# 
# iterate.downsampling(object.seurat = ea.subset, n.iterations = 5, metadata.col = "cell.ids")
# 
# levels (downsampled.object@meta.data[,metadata.col])


metadata.levels <- levels(ea.subset$cell.ids)
downsample.df <- data.frame(matrix(ncol=length(metadata.levels),nrow=0))
colnames(downsample.df) <- metadata.levels
downsample.df
n.iterations <- 100
sample_size <- ncol(aa.subset)

for(i in 1:n.iterations){
  ea.downsampled <- ea.subset[, sample(Cells(ea.subset), ncol(aa.subset)), seed = NULL]
  metadata.df <- ea.downsampled@meta.data
  for(cell.type in metadata.levels){
    n.cells <- nrow(metadata.df[metadata.df$cell.ids == cell.type,])
    downsample.df[i, cell.type] <- n.cells
  }
}
downsample.df

library(tidyr)
downsample.stats <- downsample.df %>%
  summarise(across(where(is.numeric), .fns = list(avg = mean, sd = sd), .names = "{.fn}_{.col}")) |> 
  pivot_longer(everything(), names_to = c(".value", "name"), names_sep = "_") |> 
  mutate(name = gsub("\\.", '_', tolower(name)))

ggplot(downsample.stats, aes(name, avg)) +
  geom_col() +
  geom_pointrange(aes(ymin = avg - sd, ymax = avg + sd))


# Evaluate cluster stability - scclusteval --------------------------------

library(scclusteval)



# WebGestaltR -------------------------------------------------------------

library(WebGestaltR)
library(ggplot2)
library(ggrepel)

#load  markers list from seurat `FindMarkers`
seurat.markers <- readRDS("/Users/heustonef/Desktop/PancDB_Data/PancT2D_scRNA/PancT2D_AAvsEUonly-doPar-SCTeach-allmarkers-95pctvar.rds") 

colnames(seurat.markers)
seurat.markers<- seurat.markers %>%
  select(-pct.1, -pct.2, -p_val) %>%
  filter(!cluster %in% c(1, 15))
seurat.markers$cluster <- factor(seurat.markers$cluster)
seurat.markers.cluster0 <- seurat.markers %>%
  filter(cluster == 0 & p_val_adj < 0.05) %>%
  ungroup() %>%
  select(-p_val_adj, -cluster) %>%
  relocate(gene, avg_log2FC)
seurat.markers.cluster0 <- data.frame(seurat.markers.cluster0)

# Web Gestalt Beta Markers --------------------------------------------

library(WebGestaltR)
library(ggplot2)
library(ggrepel)

web.genesets <- listGeneSet()
View(web.genesets)

#load  markers list from seurat `FindMarkers`
beta.markers <- readRDS("~/Desktop/PancDB_Data/ASHG2023/PancT2D_AAvsEUonly-doPar-Beta-markers.rds")

# reformat to include only significant markers
beta.markers$gene <- rownames(beta.markers)
beta.markers<- beta.markers %>%
  select(-pct.1, -pct.2, -p_val)
beta.markers.trim <- beta.markers %>%
  filter(p_val_adj < 0.05) %>%
  select(-p_val_adj) %>%
  relocate(gene, avg_log2FC)
beta.markers.trim <- data.frame(beta.markers.trim)
rownames(beta.markers.trim) <- 1:nrow(beta.markers.trim)

#Select databases to test
enrich.databases <- c(
                      "pathway_KEGG",
                      "pathway_Reactome",
                      # "pathway_Panther",
                      "pathway_Wikipathway",
                      "geneontology_Biological_Process",
                      # "geneontology_Cellular_Component_noRedundant",
                      "geneontology_Molecular_Function_noRedundant"
                      # "disease_OMIM",
                      # "disease_Disgenet",
                      # "drug_DrugBank",
                      # "drug_GLAD4U"
                      # "disease_GLAD4U",
                      # "phenotype_Human_Phenotype_Ontology"
                      # "network_Transcription_Factor_target",
                      # "network_TCGA_RNASeq_PAAD"
                      )
# Perform WebGestaltR analysis using`enrichMethod="GSEA"` & `organism="hsapiens"`, with FDR  as significant
enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                            enrichDatabase=enrich.databases, interestGene = beta.markers.trim,
                            interestGeneType="genesymbol", sigMethod = "fdr", minNum=10, fdrThr = 1, gseaPlotFormat = "png",
                            isOutput = TRUE, saveRawGseaResult = FALSE, projectName = paste0("Beta1vsBeta2-pwayandgo-", as.character(as.integer(Sys.time()))))

#Manually curate list in excel to limit plottable data to `FDR < 0.1  #< 0.05`. Alternatively, plot volcano using function below.

#GSEA per assay

enrichResult.list<- c()
for(val in 1:length(enrich.databases)){
  i <- enrich.databases[val]
  enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                              enrichDatabase=i, interestGene = beta.markers.trim,
                              interestGeneType="genesymbol", sigMethod = "fdr", minNum=10, fdrThr = 1, gseaPlotFormat = "png",
                              isOutput = TRUE, saveRawGseaResult = FALSE, projectName = paste0("Beta1vsBeta2-", i, "-", as.character(as.integer(Sys.time()))))
  enrichResult.list <- rbind(enrichResult.list, enrichResult)
}
View(enrichResult.list)
#ORA
for(i in enrich.databases){
  enrichResult <- WebGestaltR(enrichMethod="ORA", organism="hsapiens",
                              enrichDatabase=i, interestGene = beta.markers.trim$gene,
                              interestGeneType="genesymbol", 
                              referenceGeneType="genesymbol", 
                              sigMethod = "fdr", 
                              minNum=10, fdrThr = 1, gseaPlotFormat = "png",
                              isOutput = TRUE, saveRawGseaResult = FALSE, projectName = paste0("Beta1vsBeta2-", i, "-", as.character(as.integer(Sys.time()))))
  enrichResult.list <- c(enrichResult.list, enrichResult)
}

enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                            enrichDatabase="pathway_Wikipathway", interestGene = beta.markers.trim,
                            interestGeneType="genesymbol", sigMethod = "fdr", minNum=10, fdrThr = 1, gseaPlotFormat = "png",
                            isOutput = FALSE, saveRawGseaResult = FALSE, projectName = paste0("Beta1vsBeta2-", as.character(as.integer(Sys.time()))))


enrichResult <- enrichResult.list
colnames(enrichResult)
enrich.plot <- enrichResult[,c("normalizedEnrichmentScore", "FDR", "size", "description")]
enrich.plot$diffexpressed <- "NO"
enrich.plot$diffexpressed[enrich.plot$normalizedEnrichmentScore > 1.5 & enrich.plot$FDR < 0.05] <- "UP"
enrich.plot$diffexpressed[enrich.plot$normalizedEnrichmentScore < -1.5 & enrich.plot$FDR < 0.05] <- "DOWN"
enrich.plot$delabel <- NA
enrich.plot$delabel[enrich.plot$diffexpressed != "NO"] <- enrich.plot$description[enrich.plot$diffexpressed != "NO"]
enrich.plot$FDR[enrich.plot$FDR == 0] <- min(enrich.plot$FDR[enrich.plot$FDR > 0]/2) #Change positions where enrich.plot$FDR == 0 to something that's plot-able


mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")


p <- ggplot(data=enrich.plot, aes(x=normalizedEnrichmentScore, y=-log10(FDR), col=diffexpressed, label=delabel)) +
  theme_minimal(axis.text.x = element_text(csize = 20)) +
  # geom_text_repel(max.overlaps = 30) +
  scale_color_manual(values=mycolors) +
  ylim(c(0, max(-log10(enrich.plot$FDR)+1))) +
  xlim(c(min(enrich.plot$normalizedEnrichmentScore)-1, max(enrich.plot$normalizedEnrichmentScore)+1)) +
  geom_vline(xintercept=c(-0.6, 0.6), col="gray", linetype = "dashed", size = 1.1) +
  geom_vline(xintercept=c(0), col="black", linetype = 1, size = 1.1) +
  geom_hline(yintercept=c(0), col="black", linetype = 1, size = 1.1) +
  geom_hline(yintercept=-log10(0.05), col="gray", linetype = "dashed", size = 1.1)+
  geom_point(aes(size = size, alpha = 0.5)) + 
  scale_size(range = c(1, 10))
plot(p)

png(filename = "beta1vbeta2_volcano_all.png", height = 600, width = 1200, bg = "transparent")
plot(p)
dev.off()


# Web Gestalt Exocrine Markers --------------------------------------------

seurat.subset$SampleEthnicity <- factor(seurat.subset$SampleEthnicity)

exocrine.subset <- subset(seurat.subset, idents = c("Exocrine1", "Exocrine2"))
Idents(exocrine.subset) <- "SampleEthnicity"
levels(exocrine.subset@active.ident)



exocrine.subset <- PrepSCTFindMarkers(exocrine.subset, assay = "SCT")
exocrine.ancestry <- FindMarkers(exocrine.subset, ident.1 = "African american/Black", ident.2 = "Caucasian")
exocrine.ancestry$genes <- rownames(exocrine.ancestry)
exocrine.ancestry <- exocrine.ancestry %>%
  filter(p_val_adj < 0.05) %>%
  arrange(-abs(avg_log2FC))
rownames(exocrine.ancestry) <- exocrine.ancestry$genes
head(exocrine.ancestry)

saveRDS(exocrine.ancestry, file = paste0(rnaProject, "-ExocrineAncestry-markers.rds"))
View(exocrine.ancestry)

DotPlot(exocrine.subset, features = rownames(exocrine.ancestry[abs(exocrine.ancestry$avg_log2FC) >= 1.5,]), cols = "RdBu") + coord_flip()
png(filename = "ExocrinePancreas_AAvsEA.png", height = 1200, width = 800, bg = "transparent")
DotPlot(exocrine.subset, features = rownames(exocrine.ancestry[abs(exocrine.ancestry$avg_log2FC) >= 1.5,]), cols = "RdBu") + coord_flip()
dev.off()

png(filename = "Exocrines-AAvsEA_abstractFeatures_AAvsEA-VLN.png", height = 800, width = 1000, bg = "transparent")
VlnPlot(exocrine.subset, features = c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1"), split.by = "SampleEthnicity", split.plot = FALSE, pt.size = 0.1, combine = TRUE, cols = c("red", "blue"))
dev.off()

png(filename = "Exocrines-AAvsEA_abstractFeatures_AAvsEA-DOT.png", height = 1000, width = 800, bg = "transparent")
DotPlot(exocrine.subset, features = features.list, cols = "RdBu", dot.scale = 25) + coord_flip()
dev.off()


seurat.markers <- exocrine.ancestry
seurat.markers$gene <- rownames(seurat.markers)
seurat.markers<- seurat.markers %>%
  select(-pct.1, -pct.2, -p_val)
seurat.markers.cluster0 <- seurat.markers %>%
  filter(p_val_adj < 0.05) %>%
  select(-p_val_adj) %>%
  relocate(gene, avg_log2FC)
seurat.markers.cluster0 <- data.frame(seurat.markers.cluster0)
rownames(seurat.markers.cluster0) <- 1:nrow(seurat.markers.cluster0)

features.list <- c("PNLIP", "PRSS1", "REG1A", "REG1B", "SPINK1")
exocrine.ancestry[exocrine.ancestry$genes %in% features.list,]



enrich.databases <- c(
  "pathway_KEGG",
  # "pathway_Reactome",
  # "pathway_Panther",
  # "pathway_Wikipathway",
  # "geneontology_Biological_Process",
  # "geneontology_Cellular_Component_noRedundant",
  "geneontology_Molecular_Function_noRedundant"
  # "disease_OMIM",
  # "disease_Disgenet",
  # "drug_DrugBank",
  # "drug_GLAD4U"
  # "disease_GLAD4U",
  # "phenotype_Human_Phenotype_Ontology"
)

enrichResult <- WebGestaltR(enrichMethod="GSEA", organism="hsapiens",
                            enrichDatabase=enrich.databases, interestGene = seurat.markers.cluster0,
                            interestGeneType="genesymbol", sigMethod="top", topThr=2000, minNum=10, fdrThr = 1,
                            isOutput = TRUE, saveRawGseaResult = FALSE, projectName = paste0("Exocrines-AAvsEA", as.character(as.integer(Sys.time()))))


enrich.plot <- enrichResult[,c("normalizedEnrichmentScore", "FDR", "size", "description")]
enrich.plot$diffexpressed <- "NO"
enrich.plot$diffexpressed[enrich.plot$normalizedEnrichmentScore > 2 & enrich.plot$FDR < 0.05] <- "UP"
enrich.plot$diffexpressed[enrich.plot$normalizedEnrichmentScore < -2 & enrich.plot$FDR < 0.05] <- "DOWN"
enrich.plot$delabel <- NA
enrich.plot$delabel[enrich.plot$diffexpressed != "NO"] <- enrich.plot$description[enrich.plot$diffexpressed != "NO"]
enrich.plot$FDR[enrich.plot$FDR == 0] <- min(enrich.plot$FDR[enrich.plot$FDR > 0]/100) #Change positions where enrich.plot$FDR == 0 to something that's plot-able


mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

ggplot(data=enrich.plot, aes(x=normalizedEnrichmentScore, y=-log10(FDR), col=diffexpressed, label=delabel)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel(max.overlaps = 15, ) +
  scale_color_manual(values=mycolors) +
  ylim(c(0, max(-log10(enrich.plot$FDR)+1))) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")








# Integrate ATAC ----------------------------------------------------------

atacProject.coldata <- getCellColData(arch.proj)

n_occur <- data.frame(table(atacProject.coldata$predictedCell_Un))
n_occur[n_occur$Freq > 1,]
atacProject.coldata[atacProject.coldata$id %in% n_occur$Var1[n_occur$Freq > 1],]









# cds <- readRDS("Obesity_scRNA-SCTRegression-NW-OB-monocle3CDS.RDS")

# Differential abundance testing - speckle ------------------------------------------
library(speckle)
library(limma)
library(ggplot2)


# Get some example data which has two groups, three cell types and two 
# biological replicates in each group
# Run propeller testing for cell type proportion differences between the two 
# groups

seurat.object <- readRDS("/Users/heustonef/Desktop/PancDB_Data/PancT2D_scRNA/PancT2D_AAvsEUonly-doPar-SCTeach-95pctvar.RDS")
seurat.subset <- subset(seurat.object, idents = c("1", "15"), invert = TRUE)
levels(seurat.subset)
cluster.ids.subset <- c(
  "Exocrine1", #SCT0
  "Alpha2", #SCT2
  "Exocrine2", #SCT3--REG1A
  "Beta1", #SCT4
  "Epithelial", #SCT5--KRT18
  "Endothelial", #SCT6
  "Alpha3", #SCT7
  "Beta2", #SCT8
  "Immune1", #SCT9--IGFBP7
  "Alpha4", #SCT10
  "Immune3", #SCT11--IGFBP7/NEAT1
  "Alpha1", #SCT12
  "Immune4", #SCT13--NEAT1
  "Mast Cells", #SCT14--TPSB2
  "Macrophages", #SCT16
  "Immune2" #SCT17
)
names(cluster.ids.subset) <- levels(seurat.subset)
seurat.subset <- RenameIdents(seurat.subset, cluster.ids.subset)
seurat.subset$cell.ids <- seurat.subset@active.ident
seurat.subset$cell.ids <- factor(seurat.subset$cell.ids)
DimPlot(seurat.subset)

propeller(clusters = seurat.subset$cell.ids, sample = seurat.subset$DonorID, group = seurat.subset$SampleEthnicity, transfor = "asin", robust = F)


# Plot cell type proportions
plotCellTypeProps(clusters=seurat.object$Obesity, sample=seurat.object$SCT_snn_res.0.5)



# Differential abundance testing - scProportionTest ------------------------------------------

library(scProportionTest)


seurat.object <- readRDS("/Users/heustonef/Desktop/PancDB_Data/PancT2D_scRNA/PancT2D_AAvsEUonly-doPar-SCTeach-95pctvar.RDS")
seurat.subset <- subset(seurat.object, idents = c("1", "15"), invert = TRUE)
levels(seurat.subset)
cluster.ids.subset <- c(
  "Exocrine1", #SCT0
  "Alpha2", #SCT2
  "Exocrine2", #SCT3--REG1A
  "Beta1", #SCT4
  "Epithelial", #SCT5--KRT18
  "Endothelial", #SCT6
  "Alpha3", #SCT7
  "Beta2", #SCT8
  "Immune1", #SCT9--IGFBP7
  "Alpha4", #SCT10
  "Immune3", #SCT11--IGFBP7/NEAT1
  "Alpha1", #SCT12
  "Immune4", #SCT13--NEAT1
  "Mast Cells", #SCT14--TPSB2
  "Macrophages", #SCT16
  "Immune2" #SCT17
)
names(cluster.ids.subset) <- levels(seurat.subset)
seurat.subset <- RenameIdents(seurat.subset, cluster.ids.subset)
seurat.subset$cell.ids <- seurat.subset@active.ident

prop_test <- sc_utils(seurat.subset) #sc_utils object required for permutation_test
prop_test <- permutation_test(prop_test, cluster_identity = "cell.ids",
                              sample_1 = "Caucasian", sample_2 = "African american/Black", 
                              sample_identity = "SampleEthnicity")

prop_test

permutation_plot(prop_test, order_clusters = TRUE)

png(filename = "CellIDs-permutation_test.png", height = 800, width = 1200)
permutation_plot(prop_test, order_clusters = TRUE)
dev.off()


general.ids.subset <- c(
  "Exocrine", # "Exocrine1", #SCT0
  "Alpha", # "Alpha2", #SCT2
  "Exocrine", # "Exocrine2", #SCT3--REG1A
  "Beta", # "Beta1", #SCT4
  "Epithelial", # "Epithelial", #SCT5--KRT18
  "Endothelial", # "Endothelial", #SCT6
  "Alpha", # "Alpha3", #SCT7
  "Beta", # "Beta2", #SCT8
  "Immune", # "Immune1", #SCT9--IGFBP7
  "Alpha", # "Alpha4", #SCT10
  "Immune", # "Immune3", #SCT11--IGFBP7/NEAT1
  "Alpha", # "Alpha1", #SCT12
  "Immune", # "Immune4", #SCT13--NEAT1
  "Immune", # "Mast Cells", #SCT14--TPSB2
  "Immune", # "Macrophages", #SCT16
  "Immune" # "Immune2" #SCT17
)

prop_test <- seurat.subset
names(general.ids.subset) <- levels(prop_test)
prop_test <- RenameIdents(prop_test, general.ids.subset)
prop_test$general.ids <- prop_test@active.ident


prop_test <- sc_utils(prop_test) #sc_utils object required for permutation_test
prop_test <- permutation_test(prop_test, cluster_identity = "general.ids",
                              sample_1 = "Caucasian", sample_2 = "African american/Black", 
                              sample_identity = "SampleEthnicity")

permutation_plot(prop_test, order_clusters = TRUE)

png(filename = "GeneralIDs-permutation_test.png", height = 800, width = 1200)
permutation_plot(prop_test, order_clusters = TRUE)
dev.off()

# Differential abundance testing - edgeR ------------------------------------------

#https://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html
library(edgeR)

seurat.object <- readRDS("PancT2D_AAvsEUonly-doPar-SCTeach-95pctvar.RDS")
# For me:

# 
# tomato = Ethnicity
# celltype.mapped = integrated_snn_res.0.5
# sample = DonorID?
# pool = orig.ident?




# Step1: quantify the number of cells in each cluster:

sce <- as.SingleCellExperiment(seurat.subset)

abundances <- table(sce$cell.ids, sce$DonorID)
abundances <- unclass(abundances)
head(abundances)

extra.info <- colData(sce)[match(colnames(abundances), sce$DonorID),]
y.ab <- DGEList(abundances, samples=extra.info)
y.ab  
  
keep <- filterByExpr(y.ab, group=y.ab$samples$SampleEthnicity)
y.ab <- y.ab[keep,]
summary(keep)

design <- model.matrix(~factor(SampleEthnicity), y.ab$samples)
design

y.ab <- estimateDisp(y.ab, design, trend="none")
summary(y.ab$common.dispersion)
plotBCV(y.ab, cex=1)
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
summary(fit.ab$var.prior)
summary(fit.ab$df.prior)
plotQLDisp(fit.ab, cex=1)
res <- glmQLFTest(fit.ab, coef=ncol(design))
summary(decideTests(res))
topTags(res)

# Want to know if ethnicity alters abundances between populations
# In this case, we are aiming to identify clusters that change in abundance among the compartment of injected cells compared to the background.






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

