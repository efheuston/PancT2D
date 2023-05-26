# Run analysis on single nucleus ATAC data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Works through Signac

# Start by testing on HPAP-049

# Load libraries ----------------------------------------------------------

library(Signac)
library(Seurat)
library(patchwork)
library(ggplot2)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)

# Global parameters -------------------------------------------------------

## frequently modified
projectName <- "hpap049"
workingdir <- "./"
# regression.param <- 0
# cum.var.thresh <- 80
# resolution <- 0.5

## infrequently modified
path_to_data <- "PancDB_data/cellranger_scATAC/HPAP-049_FGC2095/"
#run.jackstraw <- TRUE
min.cells <- 3
min.features <- 200
genome <- "hg38" # Use `registered_UCSC_genomes()`

##load local functions
sourceable.functions <- list.files(path = "RFunctions", pattern = "*.R", full.names = TRUE)
invisible(sapply(sourceable.functions, source))


# Preprocess data ---------------------------------------------------------

counts <- Read10X_h5(filename = paste0(path_to_data, "outs/filtered_peak_bc_matrix.h5"))
# metadata <- 

chrom_assay <- CreateChromatinAssay(counts = counts,
																		sep = c(":", "-"),
																		genome = genome,
																		fragments = paste0(path_to_data, "outs/fragments.tsv.gz"), 
																		min.cells = 10,
																		min.features = 200
																		)

seurat.atac <- CreateSeuratObject(counts = chrom_assay,
																	assay = "peaks",
																	# meta.data = metadata
																		)

# Examine Seurat object ---------------------------------------------------

seurat.atac[['peaks']]
granges(seurat.atac)

##extract annotations

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(seurat.atac) <- annotations


# QC Metrics --------------------------------------------------------------



















