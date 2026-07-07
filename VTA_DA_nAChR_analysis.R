# ==============================================================================
# VTA Dopaminergic Neuron Ligand-Gated Ion Channel (LGIC) Expression Analysis
# Dataset: GSE235149 (10x Genomics scRNA-seq, mouse VTA, 10 samples:
#          control M1-3, control F1-3, hungry M/F, sated M/F)
#
# Pipeline: download -> per-sample Seurat objects -> merge -> normalize ->
#           filter dopaminergic (DA) neurons by marker expression ->
#           quantify nicotinic ACh receptor / LGIC subunit expression in DA cells
#
# Companion to: "Understanding the role of Neuronal Acetylcholine receptor in
#               neuropsychiatric disorder (Schizophrenia)" - minor project thesis
# ==============================================================================

# ------------------------------------------------------------------------
# 0) Packages
# ------------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GEOquery")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
})

# ------------------------------------------------------------------------
# 1) Download and unpack raw data from GEO
# ------------------------------------------------------------------------
getGEOSuppFiles("GSE235149")
untar("GSE235149/GSE235149_RAW.tar", exdir = "GSE235149/data")

# ------------------------------------------------------------------------
# 2) Organize the 10 samples into per-sample folders (Read10X expects one
#    matrix/features/barcodes triplet per directory)
# ------------------------------------------------------------------------
samples <- c("controlM1", "controlM2", "controlM3",
             "controlF1", "controlF2", "controlF3",
             "hungryM", "hungryF", "satedM", "satedF")

for (s in samples) {
  new_dir <- file.path("GSE235149/data", s)
  if (!dir.exists(new_dir)) dir.create(new_dir)

  files <- list.files("GSE235149/data", pattern = s, full.names = TRUE)
  for (f in files) {
    if (grepl("barcodes", f)) file.copy(f, file.path(new_dir, "barcodes.tsv.gz"))
    else if (grepl("features", f)) file.copy(f, file.path(new_dir, "features.tsv.gz"))
    else if (grepl("matrix", f)) file.copy(f, file.path(new_dir, "matrix.mtx.gz"))
  }
}

# ------------------------------------------------------------------------
# 3) Load each sample and merge into one Seurat object
# ------------------------------------------------------------------------
objs <- lapply(samples, function(s) {
  counts <- Read10X(data.dir = file.path("GSE235149/data", s))
  CreateSeuratObject(counts = counts, project = s)
})

combined <- merge(objs[[1]], y = objs[-1], add.cell.ids = samples, project = "VTA_DA")
combined <- JoinLayers(combined)                          # Seurat v5+ layer merge

# ------------------------------------------------------------------------
# 4) Normalize and annotate sample metadata (condition / sex)
# ------------------------------------------------------------------------
combined <- NormalizeData(combined)                       # LogNormalize
combined <- FindVariableFeatures(combined, nfeatures = 3000)
combined <- ScaleData(combined, features = VariableFeatures(combined))

combined$sample <- sub("^(.*?)_.*$", "\\1", colnames(combined))

map_cond <- c(controlM1 = "control", controlM2 = "control", controlM3 = "control",
              controlF1 = "control", controlF2 = "control", controlF3 = "control",
              hungryM = "hungry", hungryF = "hungry", satedM = "sated", satedF = "sated")
map_sex  <- c(controlM1 = "M", controlM2 = "M", controlM3 = "M",
              controlF1 = "F", controlF2 = "F", controlF3 = "F",
              hungryM = "M", hungryF = "F", satedM = "M", satedF = "F")

combined$condition <- unname(map_cond[combined$sample])
combined$sex       <- unname(map_sex[combined$sample])

# ------------------------------------------------------------------------
# 5) Filter dopaminergic (DA) neurons by canonical marker expression
#    (Th, Slc6a3/DAT, Ddc, Slc18a2/VMAT2 - detected at raw-count > 0)
# ------------------------------------------------------------------------
da_markers <- c("Th", "Slc6a3", "Ddc", "Slc18a2")
present <- intersect(da_markers, rownames(combined))
if (length(present) == 0) stop("DA markers not found; check gene symbols/case.")

expr_cnt <- GetAssayData(combined, slot = "counts")[present, , drop = FALSE]
keep <- Matrix::colSums(expr_cnt > 0) >= 1
da <- subset(combined, cells = colnames(combined)[keep])

da <- NormalizeData(da)
da <- FindVariableFeatures(da, nfeatures = 3000)
da <- ScaleData(da, features = VariableFeatures(da), verbose = FALSE)

# ------------------------------------------------------------------------
# 6) Ligand-gated ion channel (LGIC) gene panel, incl. nicotinic ACh
#    receptor subunits (Chrna/Chrnb family - focus of the thesis)
# ------------------------------------------------------------------------
lgic_genes <- unique(c(
  # nAChR
  "Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7","Chrna9","Chrna10",
  "Chrnb2","Chrnb3","Chrnb4",
  # GABA-A (+ rho)
  "Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6",
  "Gabrb1","Gabrb2","Gabrb3",
  "Gabrg1","Gabrg2","Gabrg3",
  "Gabrd","Gabre","Gabrq","Gabrp",
  "Gabrr1","Gabrr2","Gabrr3",
  # Glycine receptor
  "Glra1","Glra2","Glra3","Glra4","Glrb",
  # Ionotropic glutamate (AMPA / Kainate / NMDA)
  "Gria1","Gria2","Gria3","Gria4",
  "Grik1","Grik2","Grik3","Grik4","Grik5",
  "Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grin3a","Grin3b",
  # 5-HT3
  "Htr3a","Htr3b","Htr3c","Htr3d","Htr3e",
  # P2X purinergic
  "P2rx1","P2rx2","P2rx3","P2rx4","P2rx5","P2rx6","P2rx7",
  # ASIC
  "Asic1","Asic2","Asic3","Asic4","Asic5",
  # Zinc-activated
  "Zacn"
))
lgic_genes <- intersect(lgic_genes, rownames(da))

# ------------------------------------------------------------------------
# 7) Quantify average log-normalized expression and % of DA cells
#    expressing each LGIC subunit
# ------------------------------------------------------------------------
mat_log  <- GetAssayData(da, layer = "data")[lgic_genes, , drop = FALSE]
avg_expr <- Matrix::rowMeans(mat_log)
pct_expr <- Matrix::rowMeans(mat_log > 0) * 100

lgic_tbl <- data.frame(
  Gene = names(avg_expr),
  AvgLogNormExpr = as.numeric(avg_expr),
  PctExpressing = as.numeric(pct_expr),
  row.names = NULL
)

top20_by_avg <- lgic_tbl[order(-lgic_tbl$AvgLogNormExpr, -lgic_tbl$PctExpressing), ][1:min(20, nrow(lgic_tbl)), ]
top20_by_pct <- lgic_tbl[order(-lgic_tbl$PctExpressing, -lgic_tbl$AvgLogNormExpr), ][1:min(20, nrow(lgic_tbl)), ]

write.csv(top20_by_avg, "LGIC_top20_DA_by_avgexpr.csv", row.names = FALSE)
write.csv(top20_by_pct, "LGIC_top20_DA_by_pctexpr.csv", row.names = FALSE)

# ------------------------------------------------------------------------
# 8) Visualization
# ------------------------------------------------------------------------
top20_genes <- top20_by_avg$Gene

DotPlot(da, features = top20_genes, group.by = "condition") +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "darkred") +
  ggtitle("Top 20 Ligand-Gated Ion Channels in VTA Dopamine Neurons") +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "italic"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

DoHeatmap(da, features = top20_genes, group.by = "condition") +
  ggtitle("Expression Heatmap of Top 20 LGICs in VTA Dopamine Neurons") +
  theme(
    axis.text.y = element_text(size = 8),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
  )

# ------------------------------------------------------------------------
# 9) Save the filtered DA-neuron Seurat object
# ------------------------------------------------------------------------
saveRDS(da, file = "VTA_DA_neurons_seurat.rds")
