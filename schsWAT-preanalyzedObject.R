
# Global parameters -------------------------------------------------------

rnaProject <- "schsWAT_DEbyCluster"
regression.vars <- c("mitochondrial_percent", "ribosomal_protein_percent", "cell_cycle__phase")
cum.var.thresh <- 95
resolution <- 0.5
comp.type <- "biowulf" # one of macbookPro, biowulf
do.sctransform <- "each" # one of FALSE, each, pooled

## infrequently modified
do.doubletFinder <- TRUE
doubletFinder.splitObject <- "biosample_id"
run.jackstraw <- FALSE
min.cells <- 3
min.features <- 200
doublet.var.thresh <- 90
predicted.doubletRate <- 0.05
excluded.samples <- c()


# Directories -------------------------------------------------------------

if(grepl("mac", comp.type, ignore.case = TRUE)){
  rna.dir <- "/Users/heustonef/Desktop/PancDB_Data/PancT2D/"
  path_to_data <- "/Users/heustonef/Desktop/PancDB_Data/PancT2D/SingleCell-AdiposeTissue-Broad/"
  sourceable.functions <- list.files(path = "/Users/heustonef/Library/CloudStorage/OneDrive-NationalInstitutesofHealth/SingleCellMetaAnalysis/GitRepositories/RFunctions/",
                                     pattern = "*.R$", full.names = TRUE)
  metadata.location <- "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/"
} else if(grepl("biowulf", comp.type, ignore.case = TRUE)){
  rna.dir <- "/data/CRGGH/heustonef/broadSingleCellPortalDownloads/schsWAT/PancT2D/"
  path_to_data <- "/data/CRGGH/heustonef/broadSingleCellPortalDownloads/schsWAT/"
  sourceable.functions <- list.files(path = "/data/CRGGH/heustonef/hpapdata/RFunctions", pattern = "*.R", full.names = TRUE)
  metadata.location <- "/data/CRGGH/heustonef/broadSingleCellPortalDownloads/schsWAT/"
}

# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
# if(grepl("biowulf", comp.type, ignore.case = TRUE)){
#   library(foreach)
#   library(doParallel)
#   cl <- makeCluster(future::availableCores(), outfile = "")
# }

##load local functions
invisible(sapply(sourceable.functions, source))

# Load data ---------------------------------------------------------------

try(setwd(rna.dir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(rnaProject, "_sessionInfo.txt"))


# explore object -----------------------------------------------------------

seurat.object <- readRDS("human_all.rds")
colnames(seurat.object@meta.data)
unique(seurat.object$orig.ident)
head(seurat.object@meta.data)
unique(seurat.object$individual)
DimPlot(seurat.object, group.by = "seurat_clusters", raster=FALSE)

# test change ident
Idents(seurat.object) <- "seurat_clusters"

seurat.object <- PrepSCTFindMarkers(seurat.object, assay = "SCT")
markers.seurat.pos <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.pos <- markers.seurat.pos %>% 
  group_by(seurat_clusters) %>%
  arrange(desc(abs(avg_log2FC)), .by_group = TRUE)
saveRDS(markers.seurat.pos, file = paste0(rnaProject, "-posmarkers-", as.character(cum.var.thresh), "pctvar.rds"))

markers.seurat.all <- FindAllMarkers(seurat.object, assay = "SCT", only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
markers.seurat.all %>%
  group_by(seurat_clusters) %>%
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








