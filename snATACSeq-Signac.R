# Run analysis on single nucleus ATAC data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Works through Signac

# Start by testing on HPAP-099_FGC2381

# RUN NOTES
# as of 2023.06.16: must use remotes::install_github("stuart-lab/signac", ref = "seurat5") branch for TSS plots to work

# Global parameters -------------------------------------------------------

## frequently modified
atacProject <- "HPAP-099_FGC2381_test"
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

for(i in 1:length(sc.data)){
  counts <- Read10X_h5(filename = paste0(names(sc.data)[i], "/outs/filtered_peak_bc_matrix.h5"))
  chrom_assay <- CreateChromatinAssay(counts = counts, 
                                      sep = c(":", "-"),
                                      genome = genome,
                                      fragments = paste0(names(sc.data)[i], "/outs/fragments.tsv.gz"), 
                                      min.cells = 10,
                                      min.features = 200
  )
  
  seurat.atac <- CreateSeuratObject(counts = chrom_assay,
                                    assay = "peaks", 
                                    meta.data = read.csv(file = paste0(names(sc.data)[i], "/outs/singlecell.csv"), 
                                                         header = TRUE,
                                                         row.names = 1),
                                    project = sc.data[i]
  )
  seurat.atac <- AssignMetadata(metadata.df = metadata, seurat.object = seurat.atac)
  # remove(counts, chrom_assay)
}

# Examine Seurat object ---------------------------------------------------

seurat.atac[['peaks']]
granges(seurat.atac)

##extract annotations

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(seurat.atac) <- annotations


# QC Metrics --------------------------------------------------------------

sample.name <- unique(seurat.atac$orig.ident) ## change once code is parallelized

seurat.atac <- NucleosomeSignal(object = seurat.atac)
seurat.atac <- TSSEnrichment(object = seurat.atac, fast = FALSE)
seurat.atac$pct_reads_in_peaks <- seurat.atac$peak_region_fragments / seurat.atac$passed_filters * 100
seurat.atac$blacklist_ratio <- seurat.atac$blacklist_region_fragments / seurat.atac$peak_region_fragments

DensityScatter(seurat.atac, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

p1 <- DensityScatter(seurat.atac, x = 'nCount_peaks', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-DenScat.png"), plot.fig = p1)}else{plot(p1)} 
seurat.atac$high.tss <- ifelse(seurat.atac$TSS.enrichment > 3, 'High', 'Low')
p1 <- TSSPlot(seurat.atac, group.by = 'high.tss') + NoLegend()
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-TSSPlot.png"), plot.fig = p1)}else{plot(p1)} 

seurat.atac$nucleosome_group <- ifelse(seurat.atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
p1 <- FragmentHistogram(object = seurat.atac, group.by = 'nucleosome_group')
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-FragHist.png"), plot.fig = p1)}else{plot(p1)} 

p1 <- VlnPlot(object = seurat.atac,
        features = c("nCount_peaks", "TSS.enrichment", "blacklist_ratio", "nucleosome_signal", "pct_reads_in_peaks"),
        pt.size = 0.1,
        ncol = 5
)

if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_QC-Vln.png"), plot.fig = p1)}else{plot(p1)} 
























