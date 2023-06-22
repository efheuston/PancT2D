# Run analysis on single nucleus ATAC data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Works through Signac

# Start by testing on HPAP-099_FGC2381

# RUN NOTES
# as of 2023.06.16: must use remotes::install_github("stuart-lab/signac", ref = "seurat5") branch for TSS plots to work

# Global parameters -------------------------------------------------------

## frequently modified
atacProject <- "Test_snATACSeq-Signac"
working.dir <- "./"
records.dir <- "~/OneDrive/SingleCellMetaAnalysis/GitRepositories/PancT2D/"
res <- 0.5
comp.type <- "mac" # one of macbookPro, biowulf, or workPC
# testable.factors <- c("BMI", "obesity") # factors to query during Harmony regression
# regression.param <- 0
# cum.var.thresh <- 80
# resolution <- 0.5

## infrequently modified
min.cells <- 3
min.features <- 200
genome <- "hg38" # Use `registered_UCSC_genomes()`
excluded.samples <- NULL
fig.export <- FALSE
#run.jackstraw <- TRUE


# Directories -------------------------------------------------------------

if(grepl("mac", comp.type, ignore.case = TRUE)){
  atac.dir <- "/Users/heustonef/Desktop/PancDB_Data/PancT2D/"
  path_to_data <- "/Users/heustonef/Desktop/PancDB_data/TestFiles"
  sourceable.functions <- list.files(path = "/Users/heustonef/OneDrive/SingleCellMetaAnalysis/GitRepositories/RFunctions", pattern = "*.R$", full.names = TRUE)
  metadata.location <- "/Users/heustonef/OneDrive/SingleCellMetaAnalysis/"
} else if(grepl("biowulf", comp.type, ignore.case = TRUE)){
  atac.dir <- "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/"
  path_to_data <- "/data/CRGGH/heustonef/hpapdata/cellranger_scRNA/scRNA_transfer"
  sourceable.functions <- list.files(path = "/data/CRGGH/heustonef/hpapdata/RFunctions", pattern = "*.R", full.names = TRUE)
  metadata.location <- "/data/CRGGH/heustonef/hpapdata/"
} else {
  print("Could not match comp.type")
}

# Load libraries ----------------------------------------------------------

library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
# library(foreach)
# library(doParallel)
# cl <- makeCluster(future::availableCores(), outfile = "")


##load local functions
invisible(sapply(sourceable.functions, source))


## load data list
sc.data <- sapply(list.dirs(path = path_to_data, recursive = FALSE, full.names = TRUE), 
                  basename, 
                  USE.NAMES = TRUE)
sc.data <- sc.data[grepl(pattern = "^HPAP", sc.data)]
# sc.data <- sc.data[1]


# Select data based on inclusion criteria
metadata <- read.table(file = paste0(metadata.location, "HPAPMetaData.txt"), header = TRUE, sep = "\t", row.names = 1)

#Exclude for testing
# metadata <- metadata %>%
#   filter(grepl("Af|Cauc|Black", SampleEthnicity) &
#            # filter(grepl("Af", SampleEthnicity) & 
#            SimpDisease != "T1DM" & 
#            SimpDisease != "NoDM" &
#            !grepl("Fluidigm", scRNA_Platform) & 
#            scATAC > 0)


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



# Preprocess data ---------------------------------------------------------
atac.object <- c()
for(i in 1:length(sc.data)){
    atac.object <- c(atac.object, SignacObjectFromCellranger(samples.list = sc.data, list.index = i))
}

# Examine Seurat object ---------------------------------------------------

atac.object[['peaks']]
granges(atac.object)

##extract annotations

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(atac.object) <- annotations


# QC Metrics --------------------------------------------------------------

sample.name <- unique(atac.object$orig.ident) ## change once code is parallelized

atac.object <- NucleosomeSignal(object = atac.object)
atac.object <- TSSEnrichment(object = atac.object, fast = FALSE)
atac.object$pct_reads_in_peaks <- atac.object$peak_region_fragments / atac.object$passed_filters * 100


# Calculate blacklist depending on if included in singlecell.csv file

if(!any(atac.object$blacklist_region_fragments>0)){
  print("blacklist not calculated in singlecell.csv; using Signac::FractionCountsInRegion")
  atac.object$blacklist_ratio <- FractionCountsInRegion(
    object = atac.object,
    assay = 'peaks',
    regions = blacklist_hg38
  )
} else {
  atac.object$blacklist_ratio <- atac.object$blacklist_region_fragments / atac.object$peak_region_fragments
}


DensityScatter(atac.object, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

p1 <- DensityScatter(atac.object, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-DenScat.png"), plot.fig = p1)}else{plot(p1)} 
atac.object$high.tss <- ifelse(atac.object$TSS.enrichment > 3, 'High', 'Low')
p1 <- TSSPlot(atac.object, group.by = 'high.tss') + NoLegend()
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-TSSPlot.png"), plot.fig = p1)}else{plot(p1)} 

atac.object$nucleosome_group <- ifelse(atac.object$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p1 <- FragmentHistogram(object = atac.object, group.by = 'nucleosome_group')
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-FragHist.png"), plot.fig = p1)}else{plot(p1)} 

p1 <- VlnPlot(object = atac.object,
        features = c("nCount_peaks", "TSS.enrichment", "blacklist_ratio", "nucleosome_signal", "pct_reads_in_peaks"),
        pt.size = 0.1,
        ncol = 5
)
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-Vln.png"), plot.fig = p1)}else{plot(p1)} 



# subset object
atac.object <- subset(
  x = atac.object,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

atac.object

atac.object <- RunTFIDF(atac.object)
atac.object <- FindTopFeatures(atac.object, min.cutoff = 'q0')
atac.object <- RunSVD(atac.object)
p1 <- DepthCor(atac.object)
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_DimRed-DepthCor.png"), plot.fig = p1)}else{plot(p1)}

atac.object <- RunUMAP(atac.object, reduction = "lsi", dims = 2:30)
atac.object <- FindNeighbors(atac.object, reduction = "lsi", dims = 2:30)
atac.object <- FindClusters(atac.object, verbose = TRUE, algorithm = 3)
p1 <- DimPlot(atac.object, label = TRUE) + NoLegend()

if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_DimPlot-Clst.png"), plot.fig = p1)}else{plot(p1)}

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_04_clustering.html














