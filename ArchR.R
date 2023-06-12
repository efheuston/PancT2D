# Notes
#

# Set up ------------------------------------------------------------------

atacProject <- "PancT2D_AAvsEUonly-doPar"
res <- 0.5
testable.factors <- c("BMI", "obesity") # factors to query during Harmony regression
comp.type <- "biowulf" # one of macbookPro, biowulf, or workPC
make.arrow.files <- FALSE

excluded.samples <- NULL

# Directories -------------------------------------------------------------

if(comp.type == "macbookPro"){
	working.dir <- "/Users/heustonef/Desktop/Obesity/snATAC/"
	path_to_data <- c(list.dirs("/Users/heustonef/Desktop/PancDB_Data/scATAC_noBams/", full.names = TRUE, recursive = FALSE))
	records.dir <- "~/OneDrive-NIH/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/"
	path_to_arrow_files <- "/Users/heustonef/Desktop/Obesity/snATAC/ArrowFiles/"
	metadata.location <- "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/"
	functions.path <- "/Users/heustonef/OneDrive-NIH/SingleCellMetaAnalysis/GitRepository/scMultiomics_MetaAnalysis/RFunctions/"
} else if(comp.type == "biowulf"){
	working.dir <- "/data/CRGGH/heustonef/hpapdata/cellranger_snATAC"
	records.dir <- working.dir
	path_to_data <- c(list.dirs("/data/CRGGH/heustonef/hpapdata/cellranger_snATAC/cellrangerOuts/", full.names = TRUE, recursive = FALSE))
	path_to_arrow_files <- "/data/CRGGH/heustonef/hpapdata/cellranger_snATAC/"
	metadata.location <- "/data/CRGGH/heustonef/hpapdata/"
	funtions.path <- "/data/CRGGH/heustonef/hpapdata/RFunctions/"
	# library(vctrs, lib.loc = "/data/heustonef/Rlib_local/")
	# library(purrr, lib.loc = "/data/heustonef/Rlib_local/")
}


# Libraries ---------------------------------------------------------------

library(dplyr)
library(ArchR)
library(harmony)
# library(foreach)
# library(doParallel)
# cl <- makeCluster(future::availableCores(), outfile = "")
nThreads <- future::availableCores()
addArchRGenome("hg38")
addArchRThreads(threads = nThreads)


# Load data ---------------------------------------------------------------

try(setwd(working.dir), silent = TRUE)
writeLines(capture.output(sessionInfo()), paste0(atacProject, "_sessionInfo.txt"))


# Generate arrowFiles -----------------------------------------------------


if(make.arrow.files == TRUE){
	for(i in path_to_data){
		ifelse(file.exists(paste0(i, "/outs")),"", path_to_data <-path_to_data[!path_to_data %in% i])
	}
	names(path_to_data) <- sapply(path_to_data, basename)
	arrowfiles <- createArrowFiles(
		inputFiles = paste0(path_to_data, "/outs/fragments.tsv.gz"),
		sampleNames = names(path_to_data),
		minTSS = 4, 
		minFrags = 1000, 
		addTileMat = TRUE,
		addGeneScoreMat = TRUE, force = TRUE
	)
	dbltScores <- addDoubletScores(input = (arrowfiles), k = 10, knnMethod = "UMAP", LSIMethod = 1) #25 samples = 40min
	saveRDS(arrowfiles, paste0(atacProject, "-ArrowFiles.RDS"))
}



# Subset Arrow Files ------------------------------------------------------

arrowfiles <- list.files(path = path_to_arrow_files, pattern = "*.arrow")
arrowfiles <- sort(arrowfiles)

# Select data based on inclusion criteria
metadata <- read.table(file = paste0(metadata.location, "HPAPMetaData.txt"), header = TRUE, sep = "\t", row.names = 1)
metadata <- metadata %>%
  filter(grepl("Af|Cauc|Black", SampleEthnicity) &
           SimpDisease != "T1DM" & 
  			 	 SimpDisease != "NoDM" &
           !grepl("Fluidigm", scRNA_Platform) & 
           scATAC > 0)


# Exclude samples that failed CellRanger QC
for(i in arrowfiles){
  x <- strsplit(i, "_")[[1]][1]
  if(i %in% excluded.samples){
    arrowfiles <- arrowfiles[arrowfiles!=i]}
  if(!(x %in% rownames(metadata))){
    arrowfiles <- arrowfiles[arrowfiles!=i]
  }
}

# clean up metadata table to exclude data not relevant to patient
metadata <- metadata[,1:12]



# Create arch.proj -----------------------------------------------------

arch.proj <- ArchRProject(ArrowFiles = sapply(arrowfiles, function(x){paste0(path_to_arrow_files, x)}), outputDirectory = working.dir, copyArrows = FALSE)

# remove samples with unreasonable numbers of cells
for(i in getSampleNames(arch.proj)){
  # print(paste(i, "-", length(which(arch.proj$Sample %in% i))))
  if(length(which(arch.proj$Sample %in% i)) < 500){
    excluded.samples <- c(excluded.samples, i)
  }
  if(length(which(arch.proj$Sample %in% i)) > 20000){
    excluded.samples <- c(excluded.samples, i)
  }
}
arch.proj <- arch.proj[which(arch.proj$Sample %ni% excluded.samples),]
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


cellcoldata <- getCellColData(arch.proj, select = c("log10(nFrags)", "TSSEnrichment"))

# plot unique nuclear fragments vs TSS enrichment score
ggPoint(
	x = cellcoldata[,1],
	y = cellcoldata[,2],
	colorDensity = TRUE,
	continuousSet = "sambaNight", 
	xlabel = "Log10 Unique Fragments",
	ylabel = "TSS Enrichment",
	xlim = c(log10(500), quantile(cellcoldata[,1], probs = 0.99)),
	ylim = c(0, quantile(cellcoldata[,2], probs = 0.99))) +
	geom_hline(yintercept = 4, lty = "dashed") +
	geom_vline(xintercept = 3, lty = "dashed")


# plot TSSEnrichment ridge plot per sample
plotGroups(
	ArchRProj = arch.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "TSSEnrichment", 
	plotAs = "ridges"
) +
	geom_vline(xintercept = 8, lty = "dashed") 

# -- TSS plots show some samples have 2 TSS enrichment at 2 different scores. Will set cutoff @ 8 to maintain a single enriched peak.

# plot nFrags ridge plot per sample
plotGroups(
	ArchRProj = arch.proj,
	groupBy = "Sample", 
	colorBy = "cellColData",
	name = "log10(nFrags)", 
	plotAs = "ridges") +
	geom_vline(xintercept = log10(1000), lty = "dashed") +
	geom_vline(xintercept = log10(30000), lty = "dashed")
metadata$DonorID <- rownames(metadata)

arch.proj@cellColData[,names(metadata)] <- lapply(names(metadata), function(x){
	arch.proj@cellColData[[x]] <- metadata[match(vapply(strsplit(as.character(arch.proj$Sample), "_"), `[`, 1, FUN.VALUE = character(1)), metadata$DonorID), x]
}
)

# Filter  ---------------------------------------------------------

arch.proj <- filterDoublets(ArchRProj = arch.proj, cutEnrich = 1, filterRatio = 1.5) # see notes
arch.proj <- arch.proj[which(arch.proj$TSSEnrichment > 8 & 
														 	arch.proj$nFrags > 1000 & 
														 	arch.proj$nFrags<40000 &
														 	arch.proj$BlacklistRatio < 0.03)]
# >1000 : https://www.nature.com/articles/s41586-021-03604-1#Sec9

arch.proj <- addIterativeLSI(
	ArchRProj = arch.proj,
	useMatrix = "TileMatrix", 
	name = "IterativeLSI", 
	iterations = 10, 
	clusterParams = list( #See Seurat::FindClusters
		resolution = c(0.3), 
		sampleCells = 10000, 
		n.start = 10
	), 
	varFeatures = 25000, 
	dimsToUse = 1:30,
	force = TRUE)



# Run Harmony -------------------------------------------------------------


# factorize regression columns
arch.proj$SampleAge <- as.factor(arch.proj$SampleAge)
arch.proj <- addHarmony(ArchRProj = arch.proj, 
												reducedDims = "IterativeLSI", 
												name = "Harmony", 
												groupBy = c("Sample", "SampleSex", "SampleAge"), 
												max.iter.harmony = 20, #did not converge after 10
												force = TRUE) # addHarmony "groupby" defines variables to correct for

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
saveRDS(arch.proj, paste0(atacProject, "_Harmony.RDS"))



arch.proj <- addUMAP(ArchRProj = arch.proj,
										 reducedDims = "Harmony",
										 name = "UMAP_harmony",
										 nNeighbors = 30,
										 minDist = 0.5,
										 metric = "cosine",
										 force = TRUE)
arch.proj <- addClusters(input = arch.proj, reducedDims = "Harmony", method = "Seurat", name = paste0("Harmony_res", as.character(res)), resolution = res, force = TRUE)

plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "Obesity", embedding = "UMAP_harmony")
plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "BMI", embedding = "UMAP_harmony", plotAs = "points")
p1 <- plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = "Obesity", embedding = "UMAP_harmony", randomize = TRUE)
p2 <- plotEmbedding(ArchRProj = arch.proj, colorBy = "cellColData", name = paste0("Harmony_res", as.character(res)), embedding = "UMAP_harmony")
ggAlignPlots(p1, p2, type = "h")

table(getCellColData(ArchRProj = arch.proj, select = paste0("Harmony_res", as.character(res))))
cM <- confusionMatrix(paste0(arch.proj$Harmony_res0.5), paste0(arch.proj$SampleEthnicity)) # Could not automate this line
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)

cM <- confusionMatrix(paste0(arch.proj@cellColData[,paste0("Harmony_res", as.character(res))]), paste0(arch.proj$SampleEthnicity))
cM <- cM / Matrix::rowSums(cM)
pheatmap::pheatmap(
	mat = as.matrix(cM),
	color = paletteContinuous("whiteBlue"),
	border_color = "black"
)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)





# Identify marker genes ---------------------------------------------------

#Without MAGIC
markergenes <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "GeneScoreMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
markerList <- getMarkers(markergenes, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markergenes, file = paste0(atacProject, "-markergenes.RDS"))
markerList$C1

#With MAGIC
# arch.proj.magic <- addImputeWeights(ArchRProj = arch.proj, reducedDims = "Harmony")
# saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


# Calling peaks -----------------------------------------------------------
arch.proj <- loadArchRProject(working.dir)
library(BSgenome.Hsapiens.UCSC.hg38)
pathToMacs2 <- findMacs2()
BSgenome.Hsapiens.UCSC.hg18 <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

arch.proj <- addGroupCoverages(arch.proj, groupBy = paste0("Harmony_res", as.character(res)))
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)
arch.proj <- addReproduciblePeakSet(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), pathToMacs2 = pathToMacs2)

arch.proj <- addPeakMatrix(arch.proj)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


markerPeaks <- getMarkerFeatures(arch.proj, groupBy = paste0("Harmony_res", as.character(res)), useMatrix = "PeakMatrix", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon")
saveRDS(markerPeaks, file = paste0(atacProject, "-MarkerPeaks.RDS"))


markerPeaks <- readRDS(paste0(atacProject, "-MarkerPeaks.RDS"))

markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
heatmapPeaks <- plotMarkerHeatmap(seMarker = markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", transpose = TRUE)


arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", annoName = "cisbp", force = TRUE)
arch.proj <- addMotifAnnotations(arch.proj, motifSet = "encode", annoName = "encode", force = TRUE)
arch.proj <- addMotifAnnotations(arch.proj, motifSet = "homer", annoName = "homer", force = TRUE)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)




arch.proj <- addArchRAnnotations(ArchRProj = arch.proj, collection = "EncodeTFBS")
enrichEncode <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "EncodeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapEncode <- plotEnrichHeatmap(enrichEncode, n = 7, transpose = TRUE)
ComplexHeatmap::draw(heatmapEncode, heatmap_legend_side = "bot", annotation_legend_side = "bot")

arch.proj <- addArchRAnnotations(ArchRProj = arch.proj, collection = "CistromeTFBS")
enrichcistrome <- peakAnnoEnrichment(seMarker = markerPeaks, ArchRProj = arch.proj, peakAnnotation = "CistromeTFBS", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
heatmapcistrome <- plotEnrichHeatmap(enrichcistrome, n = 7, transpose = FALSE)
ComplexHeatmap::draw(heatmapcistrome, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# Motif Deviations --------------------------------------------------------

arch.proj <- addBgdPeaks(arch.proj)
if("cisbMotif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "cisbp", annoName = "cisbpMotif")
}
if("encodeMotif" %ni% names(arch.proj@peakAnnotation)){
	arch.proj <- addMotifAnnotations(arch.proj, motifSet = "encode", annoName = "encodeMotif")
}
arch.proj <- addDeviationsMatrix(arch.proj, peakAnnotation = "encodeMotif", force = TRUE)
plotVarDev <- getVarDeviations(arch.proj, name = "encodeMotifMatrix", plot = TRUE)
saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)



# Footprinting ------------------------------------------------------------

motifPositions <- getPositions(arch.proj)


# Deviant Motifs ----------------------------------------------------------


seGroupMotif <- getGroupSE(arch.proj, useMatrix = "MotifMatrix", groupBy = paste0("Harmony_res", as.character(res)))
saveRDS(seGroupMotif, file = paste0(atacProject, "_seGroupMotif.RDS"))
corGSM_MM <- correlateMatrices(arch.proj, useMatrix1 = "GeneScoreMatrix", useMatrix2 = "MotifMatrix", reducedDims = "Harmony")
saveRDS(corGSM_MM, file = paste0(atacProject, "_corGSM_MM.RDS"))

saveArchRProject(ArchRProj = arch.proj, outputDirectory = working.dir, load = TRUE)


corGSM_MM <- readRDS(paste0(atacProject, "_corGSM_MM.RDS"))




# Integrate scRNA object (Seurat) -----------------------------------------

#seurat.object <- readRDS("/Users/heustonef/Desktop/Obesity/scRNA/Obesity_scRNA-Anchored-NW-OB_seuratMarkers-90pctvar.RDS")
seurat.object <- readRDS("/data/CRGGH/heustonef/hpapdata/cellranger_snATAC/Obesity_scRNA-Anchored-NW-OB-95pctvar.RDS")
#check import
colnames(seurat.object@meta.data)
seurat.object$integrated_snn_res.0.5 <- paste0("SCT", seurat.object$integrated_snn_res.0.5)

arch.proj <- addGeneIntegrationMatrix(  # step takes ~95min
	ArchRProj = arch.proj,
	useMatrix = "GeneScoreMatrix",
	matrixName = "GeneIntegrationMatrix", 
	reducedDims = "Harmony",
	seRNA = seurat.object,
	addToArrow = FALSE,
	groupRNA = "integrated_snn_res.0.5",
	nameCell = "predictedCell_Un",
	nameGroup = "predictedGroup_Un",
	nameScore = "predictedScore_Un"
	
)

saveArchRProject(arch.proj, outputDirectory = "archRvsRNApctvar95", load = TRUE)
pal <- paletteDiscrete(values = seurat.object$integrated_snn_res.0.5)
plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "predictedGroup_Un", pal = pal)

cM <- as.matrix(confusionMatrix(arch.proj$Harmony_res0.5, arch.proj$predictedGroup_Un))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM)) #Assignments
saveArchRProject(arch.proj, outputDirectory = "archRvsRNApctvar95", load = TRUE)


# Trajectory --------------------------------------------------------------
arch.proj <- loadArchRProject(working.dir)

pal <- paletteDiscrete(values = seurat.object$integrated_snn_res.0.5)

plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "predictedGroup_Un") +
	theme_ArchR(legendTextSize = 10)
plotEmbedding(arch.proj, embedding = "UMAP_harmony", colorBy = "cellColData", name = "Harmony_res0.5",) +
	theme_ArchR(legendTextSize = 10)



# Heatmaps ----------------------------------------------------------------

arch.markers <- readRDS(paste0(working.dir, "Obesity_scHPAP-markergenes.RDS"))
heatmap.islets <- plotMarkerHeatmap(seMarker = arch.markers, 
																		cutOff = "FDR <= 0.01 & Log2FC >=1.25",
																		labelMarkers = unlist(panc.markers),
																		transpose = TRUE)


heatmap.plot <- ComplexHeatmap::draw(heatmap.islets, heatmap_legend_side = "bot", annotation_legend_side = "bot")

png(filename = paste0(atacProject, "-UMAP_harmony-res", as.character(res), "-AllMarkersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
plot(heatmap.plot)
dev.off()

for(i in 1:length(panc.markers)){
	chart.name <- names(panc.markers[i])
	gene.set <- panc.markers[i]
	
	heatmap.islets <- plotMarkerHeatmap(seMarker = arch.markers, 
																			cutOff = "FDR <= 0.01 & Log2FC >=1.25",
																			labelMarkers = unlist(gene.set),
																			transpose = TRUE)
	
	
	heatmap.plot <- ComplexHeatmap::draw(heatmap.islets, heatmap_legend_side = "bot", annotation_legend_side = "bot")
	
	png(filename = paste0(atacProject, "-UMAP_harmony-res", as.character(res), "-", as.character(chart.name), "Markersheatmap.png"), height= 800, width = 1600, bg = "transparent", res = 100)
	plot(heatmap.plot)
	dev.off()
	
	
}
panc.markers

obnw.markertest <- getMarkerFeatures(
  ArchRProj = arch.proj, 
  useMatrix = "PeakMatrix",
  groupBy = "Obesity",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "OB",
  bgdGroups = "NW"
)
