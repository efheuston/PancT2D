# Run analysis on single nucleus ATAC data
# Sample data is from [PandDB](https://hpap.pmacs.upenn.edu/)
# Works through Signac

# Start by testing on HPAP-099_FGC2381 & HPAP-100_FGC2414

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
fc.cutoff <- 3
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
library(cicero)
library(SeuratWrappers)
library(patchwork)
library(ggplot2)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(GenomicRanges)
# library(foreach)
# library(doParallel)
# cl <- makeCluster(future::availableCores(), outfile = "")


##load local functions
invisible(sapply(sourceable.functions, source))


# Preprocess data ---------------------------------------------------------

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

# Create seurat objects ---------------------------------------------------

#Combined peak set
sample.peaks <- lapply(1:length(sc.data), function(x){
  read.table(
    file = paste0(names(sc.data)[x], "/outs/peaks.bed"),
    col.names = c("chr", "start", "end")
  )
  }
)
gr.data <- sapply(sample.peaks, makeGRangesFromDataFrame)
gr.data <- GRangesList(gr.data)
combined.peaks <-reduce(x = c(unlist(gr.data)))

names(gr.data) <- sc.data

peak.widths <- width(combined.peaks)
combined.peaks <- combined.peaks[peak.widths < 10000 & peak.widths > 20]

seurat.atac <- c()
seurat.atac <- lapply(names(sc.data), function(x) {SignacObjectFromCommonPeakset(x, combined.peaks)})



# Examine Seurat object ---------------------------------------------------

##extract annotations

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
Annotation(seurat.atac) <- annotations


# QC Metrics --------------------------------------------------------------

sample.name <- unique(seurat.atac$orig.ident) ## change once code is parallelized

seurat.atac <- NucleosomeSignal(object = seurat.atac)
seurat.atac <- TSSEnrichment(object = seurat.atac, fast = FALSE)
seurat.atac$pct_reads_in_peaks <- seurat.atac$peak_region_fragments / seurat.atac$passed_filters * 100


# Calculate blacklist depending on if included in singlecell.csv file

if(!any(seurat.atac$blacklist_region_fragments>0)){
  print("blacklist not calculated in singlecell.csv; using Signac::FractionCountsInRegion")
  seurat.atac$blacklist_ratio <- FractionCountsInRegion(
    object = seurat.atac,
    assay = 'peaks',
    regions = blacklist_hg38
  )
} else {
  seurat.atac$blacklist_ratio <- seurat.atac$blacklist_region_fragments / seurat.atac$peak_region_fragments
}


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



# subset object
seurat.atac <- subset(
  x = seurat.atac,
  subset = nCount_peaks > 3000 &
    nCount_peaks < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)

temp <- seurat.atac[[1]]


# Doublet detection  --------------------------------------------
# Can try scrublet, as tested in https://www.nature.com/articles/s41467-021-21583-9#Sec8, using the Seurat's Label transfer methods
# Can import doublet enrichment scores from ArchR
# Can us Amulet method and/or scDblFinder (Doublet identifiation in single-cell ATAC-seq)

#import doublet scores from ArchR
doublet.df <- read.csv(file = paste0(atac.dir, "ArchR_allSampleProj-cellColData.txt"),
                       sep = "\t", header = TRUE,
                       row.names = 1)
doublet.df <- doublet.df[,grepl("Doublet|Sample", colnames(doublet.df), ignore.case = TRUE)]
doublet.df$cellID <- sapply(strsplit(rownames(doublet.df),"#"), `[`, 2)



head(seurat.atac[[1]]@meta.data)




# Linear dimensional reduction --------------------------------------------

seurat.atac <- RunTFIDF(seurat.atac)
seurat.atac <- FindTopFeatures(seurat.atac, min.cutoff = 'q0')
seurat.atac <- RunSVD(seurat.atac)
p1 <- DepthCor(seurat.atac)
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_DimRed-DepthCor.png"), plot.fig = p1)}else{plot(p1)}

# Non-linear dimensional reduction  -------------------------

seurat.atac <- RunUMAP(seurat.atac, reduction = "lsi", dims = 2:30)
seurat.atac <- FindNeighbors(seurat.atac, reduction = "lsi", dims = 2:30)
seurat.atac <- FindClusters(seurat.atac, verbose = TRUE, algorithm = 3)

saveRDS(seurat.atac, file = paste0(atacProject, "_SignacObject.RDS"))

p1 <- DimPlot(seurat.atac, label = TRUE) + NoLegend()

if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_DimPlot-Clst.png"), plot.fig = p1)}else{plot(p1)}

#https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_04_clustering.html


# Gene activity matrix by accessibility frequency -------------------------

#fragment count in promoters and TSS
gene.activities <- GeneActivity(seurat.atac)
seurat.atac[['GeneActivity']] <- CreateAssayObject(counts = gene.activities)
seurat.atac <- NormalizeData(seurat.atac, assay = "GeneActivity", normalization.method = "LogNormalize", scale.factor = median(seurat.atac$nCount_GeneActivity))
saveRDS(seurat.atac, file = paste0(atacProject, "_SignacObject.RDS"))

#Cicero

cicero.atac <- as.cell_data_set(x = seurat.atac)
cicero.atac <- make_cicero_cds(cicero.atac)

genome <- seqlengths(seurat.atac)[1]
genome.df <- data.frame("chr" = names(genome), "length" <- genome)
conns <- run_cicero(cicero.atac, genomic_coords = genome.df, sample_num = 100)
saveRDS(conns, file = paste0(atacProject, "_ciceroConnections.RDS"))

ccans <- generate_ccans(conns)
saveRDS(conns, file = paste0(atacProject, "_ciceroNetworks.RDS"))

seurat.links <- ConnectionsToLinks(conns = conns, ccans = ccans)
Links(seurat.atac) <- seurat.links
saveRDS(seurat.atac, file = paste0(atacProject, "_SignacObject.RDS"))



# Integrate with scRNA data -----------------------------------------------

#Is on the to-do list



# Differentially accessible peaks -----------------------------------------

DefaultAssay(seurat.atac) <- 'peaks'

da.peaks <- FindAllMarkers(seurat.atac, assay = "peaks", test.use = "LR", latent.vars = "nCount_peaks")
saveRDS(da.peaks, file = paste0(atacProject, "_daPeaks-all.RDS"))

p1 <- VlnPlot(seurat.atac, features = rownames(da.peaks)[1], pt.size = 0.1)
p2 <- FeaturePlot(seurat.atac, features = rownmaes(da.peaks)[1], pt.size = 0.1, max.cutoff = 'q95')

if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_daPeaks-all.png"), plot.fig = c(p1, p2))}else{plot(p1, p2)}

gene.fcup <- ClosestFeature(seurat.atac, rownames(da.peaks[da.peaks$avg_log2FC > fc.cutoff, ]))
gene.fcdn <- ClosestFeature(seurat.atac, rownames(da.peaks[da.peaks$avg_log2FC < -1*fc.cutoff, ]))


##create workbook
markers.table <- openxlsx::createWorkbook()


##write positive markers to table
openxlsx::addWorksheet(gene.fcup, sheetName = paste0("FC", as.character(fc.cutoff), "up"))
openxlsx::writeData(gene.fcup, sheet = paste0("FC", as.character(fc.cutoff), "up"), x = gene.fcup, startCol = 1, startRow = 1, colNames = TRUE)

##write all markers to table
openxlsx::addWorksheet(markers.table, sheetName = paste0("FC", as.character(fc.cutoff), "dn"))
openxlsx::writeData(markers.table, sheet = paste0("FC", as.character(fc.cutoff), "dn"), x = gene.fcdn, startCol = 1, startRow = 1, colNames = TRUE)

##save workbook
openxlsx::saveWorkbook(wb = markers.table, file = paste0(atacProject, "_closestGene-", as.character(cum.var.thresh), "pctvar.xlsx"), overwrite = TRUE, returnValue = TRUE)




# Plotting genomic regions ------------------------------------------------

#Can set plotting order of regions
#levels(seurat.atac) <- c()
p1 <- CoveragePlot(seurat.atac, region = c("INS", "GCG"), extend.upsteram = 1000, extend.downstream = 1000, ncol = 1)
if(fig.export == TRUE){export.figs(plot.name = paste0(sample.name, "_CoveragePlot-all.png"), plot.fig = p1)}else{plot(p1)}




