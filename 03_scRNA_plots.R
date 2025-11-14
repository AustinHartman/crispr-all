# This script builds the initial Seurat object from Cellranger outputs
# performs initial QC filtering of cells
# assigns perturbations to cells using deMULTIplex2 for classification
# and performs donor integration using Harmony

library(Seurat)
library(ggplot2)
library(BuenColors)
library(ggrepel)
library(deMULTIplex2)

options(future.globals.maxSize = 8000 * 1024^2)

colors <- list(
  "Cholla" = "#FF5A32",
  "Torch" = "#D88D21",
  "Saguaro" = "#4A6F28",
  "PricklyPear" = "#4A6F28",
  "Kingcup" = "#6FB5D9",
  "Control" = "#CDD5E0"
)

### build seurat object ###
data.dir <- "/home/hartmana_stanford_edu/perturb_crispr_all/output/cellranger_outs/cellranger_sbatch_scripts"
full.path <- list.dirs(data.dir, recursive = FALSE)
objects.list <- list()
for (p in full.path) {
  tmp <- paste0(p, "/outs/filtered_feature_bc_matrix/")
  expression_data <- Read10X(data.dir = tmp)
  obj <- CreateSeuratObject(counts = expression_data)
  objects.list[[p]] <- obj
}
obj <- merge(
  objects.list[[1]], y = unname(objects.list[2:length(objects.list)]),
  add.cell.ids = c(
    "dRL320_A", "dRL320_B", "dRL320_C", "dRL321_A",
    "dRL321_B", "dRL322_A", "dRL322_B"
  ), project = "modpoki")

### filter low quality cells ###
DefaultAssay(obj) <- "RNA"
obj[["SCT"]] <- NULL
obj <- JoinLayers(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")
VlnPlot(object = obj, group.by="replicate", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)
VlnPlot(object = obj, group.by = "replicate", features = "percent.mt", pt.size=0) + geom_hline(yintercept = 10)
obj <- subset(obj, subset = nFeature_RNA > 1000 & nCount_RNA < 40000 & percent.mt < 10)

# run demultiplex2 to assign cell barcodes to perturbation barcodes
mtx <- read.csv("/home/hartmana_stanford_edu/perturb_crispr_all/output/matrices/full/03_full_matrix.csv")
rownames(mtx) <- mtx$bc
mtx$bc <- NULL
res <- demultiplexTags(mtx, # Required, the tag count matrix from your experiment, can be either dense or sparse
                       plot.path = "/home/hartmana_stanford_edu/perturb_crispr_all/output/matrices/demultiplex2_assignment", # Where to output a summary plot
                       plot.name = "07_demultiplex_assign", # text append to the name of the summary plot file
                       plot.diagnostics = TRUE) # Whether to output diagnostics plots for each tag
table(res$final_assign)
mtx$assignment <- res$final_assign
write.csv(x = mtx, file = "/home/hartmana_stanford_edu/perturb_crispr_all/output/seurat/03_full_matrix_assignments.csv")
df <- read.csv("/home/hartmana_stanford_edu/perturb_crispr_all/output/seurat/03_full_matrix_assignments.csv")
rownames(df) <- df$X
values <- df$assignment
names(values) <- rownames(df)
obj$construct <- values
# add perturbation class
result <- sapply(strsplit(as.character(obj$construct), ".", fixed = TRUE), function(x) x[[1]])
obj$class <- result

# dimentional reduction
obj <- JoinLayers(obj)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:30)

# remove negative and multiplet cells
Idents(obj) <- "class"
remove.cells <- WhichCells(obj, idents=c("negative", "multiplet"))
obj <- subset(obj, cells=remove.cells, invert=T)
Idents(obj) <- "construct"

# donor integration
obj <- JoinLayers(obj)
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$donor)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features=VariableFeatures(obj))    
obj <- RunPCA(obj, features = VariableFeatures(object = obj), npcs = 30)
obj <- IntegrateLayers(
  object = obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "integrated.harmony",
  verbose = FALSE)
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj <- RunUMAP(obj, reduction = "integrated.harmony", dims = 1:30, reduction.name = "umap.harmony")
DimPlot(obj, reduction="umap.harmony", group.by="class") + scale_color_manual(values=colors)
DimPlot(obj, reduction="umap.harmony", group.by="construct", split.by="construct", ncol=10) + NoLegend()

# load Seurat object on Zenodo
obj <- readRDS("/home/hartmana_stanford_edu/projects/perturb_crispr_all/output/seurat/24_seurat_object_assigned_singlets_only.rds")

# feature plots
fp1 <- FeaturePlot(obj, reduction="umap.harmony", feature = "CD4", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp2 <- FeaturePlot(obj, reduction="umap.harmony", feature = "CD8A", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp3 <- FeaturePlot(obj, reduction="umap.harmony", feature = "MKI67", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp4 <- FeaturePlot(obj, reduction="umap.harmony", feature = "GNLY", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp5 <- FeaturePlot(obj, reduction="umap.harmony", feature = "GZMB", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp6 <- FeaturePlot(obj, reduction="umap.harmony", feature = "FOXP3", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp7 <- FeaturePlot(obj, reduction="umap.harmony", feature = "CCR7", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
fp8 <- FeaturePlot(obj, reduction="umap.harmony", feature = "CXCR6", min.cutoff = "q05", max.cutoff="q95", cols=c("gray92", "black"), order = T)
(fp1 | fp2 | fp3 | fp4) / (fp5 | fp6 | fp7 | fp8)

# construct highlights
Idents(obj) <- "construct"
plots <- list()
for (const in c("Kingcup.v1.41BBz", "Saguaro.v1.LTBR", "Torch.v1.miR.155", "Saguaro.v1.FOXO1")) {
  if (startsWith(const,"Cholla")) {
    p <- DimPlot(obj, reduction = "umap.harmony", cells.highlight = WhichCells(obj, idents=const), size = 0.3, alpha=0.5) + ggtitle(const) + scale_color_manual(values=c(colors$Control, colors[["Cholla"]]))
  } else if (startsWith(const,"Torch")) {
    p <- DimPlot(obj, reduction = "umap.harmony", cells.highlight = WhichCells(obj, idents=const), size = 0.5, alpha=0.5) + ggtitle(const) + scale_color_manual(values=c(colors$Control, colors[["Torch"]]))
  } else if (startsWith(const,"Saguaro") || startsWith(const, "PricklyPear")) {
    p <- DimPlot(obj, reduction = "umap.harmony", cells.highlight = WhichCells(obj, idents=const), size = 0.5, alpha=0.5) + ggtitle(const) + scale_color_manual(values=c(colors$Control, colors[["Saguaro"]]))
  } else if (startsWith(const,"Kingcup")) {
    p <- DimPlot(obj, reduction = "umap.harmony", cells.highlight = WhichCells(obj, idents=const), size=1.5) + ggtitle(const) + scale_color_manual(values=c(colors$Control, colors[["Kingcup"]]))
  }
  plots[[const]] <- p
}
wrap_plots(plots)

# unsupervised clustering
obj <- FindNeighbors(obj, reduction = "integrated.harmony", dims = 1:20)
obj <- FindClusters(obj, resolution = 1)
Idents(obj) <- "seurat_clusters"
obj$cluster <- Idents(obj)
cluster_construct_df <- obj@meta.data %>%
  group_by(cluster, construct) %>%
  summarise(count = n()) %>%
  ungroup()
total_cells_per_construct <- obj@meta.data %>%
  group_by(construct) %>%
  summarise(total_count = n())
cluster_construct_df <- merge(cluster_construct_df, total_cells_per_construct, by = "construct")
cluster_construct_df <- cluster_construct_df %>%
  mutate(percentage = (count / total_count) * 100)
p1 <- ggplot(cluster_construct_df, aes(x = construct, y = percentage, fill = as.factor(cluster))) +
  geom_bar(stat = "identity") +
  L_border() + scale_fill_manual(values = jdb_palette("corona")) +
  labs(title = "Cluster Abundance", x = "Construct", y = "Percentage", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), text=element_text(size=15))
p2 <- DimPlot(obj, group.by="cluster", reduction="umap.harmony") + scale_color_manual(values = jdb_palette("corona"))
p1 | p2

cluster.1.markers <- FindMarkers(obj, ident.1 = 1, lfc.threshold=0.5, only.pos=T)
cluster.6.markers <- FindMarkers(obj, ident.1 = 6, lfc.threshold=0.5, only.pos=T)

################################################################################
############# MED12-KO and LTBR-OE comparison w/ existing data #################
################################################################################

# Load CRISPR-All-seq Seurat object
# available from Zenodo (crispr_all_seq_seurat_object.rds):
# https://zenodo.org/records/17594332
obj <- readRDS("/home/hartmana_stanford_edu/projects/perturb_crispr_all/output/seurat/24_seurat_object_assigned_singlets_only.rds")

# OverCITE-seq input matrices from GEO
# Data available at: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA797432
gex.mtx <- read.csv("/home/hartmana_stanford_edu/projects/perturb_crispr_all/data/ltbr_paper_sc/GSM5819660_GEX_counts.csv.gz")
orf.mtx <- read.csv("/home/hartmana_stanford_edu/projects/perturb_crispr_all/data/ltbr_paper_sc/GSM5819659_ORF_counts.csv.gz")
adt.mtx <- read.csv("/home/hartmana_stanford_edu/projects/perturb_crispr_all/data/ltbr_paper_sc/GSM5819658_ADT_counts.csv.gz")
hto.mtx <- read.csv("/home/hartmana_stanford_edu/projects/perturb_crispr_all/data/ltbr_paper_sc/GSM5819657_HTO_counts.csv.gz")

# Use the ORF matrix to re-assign perturbation labels (publication labels not found on GEO)
rownames(orf.mtx) <- orf.mtx$X
orf.mtx$X <- NULL
orf.mtx.t <- t(orf.mtx)
res <- demultiplexTags(orf.mtx.t, # Required, the tag count matrix from your experiment, can be either dense or sparse
                       plot.path = "/home/hartmana_stanford_edu/projects/perturb_crispr_all/output/matrices/ltbr_demultiplex2_assignment", # Where to output a summary plot
                       plot.name = "24_demultiplex_assign", # text append to the name of the summary plot file
                       plot.diagnostics = TRUE) # Whether to output diagnostics plots for each tag
table(res$final_assign)
orf.mtx$assignment <- res$final_assign
write.csv(x = orf.mtx, file = "/home/hartmana_stanford_edu/projects/perturb_crispr_all/output/matrices/ltbr_demultiplex2_assignment/24_full_matrix_assignments.csv")

# Create a Seurat object with the cell x gene matrix and add the perturbation labels as metadata
rownames(gex.mtx) <- gex.mtx$X
gex.mtx$X <- NULL
ltbr.obj <- CreateSeuratObject(counts=gex.mtx)
ltbr.obj <- AddMetaData(ltbr.obj, metadata=res$final_assign, col.name="construct")

# Spot-check visualization of assignments 
ltbr.obj <- NormalizeData(ltbr.obj)
ltbr.obj <- FindVariableFeatures(ltbr.obj, selection.method = "vst", nfeatures = 2000)
ltbr.obj <- ScaleData(ltbr.obj)
ltbr.obj <- RunPCA(ltbr.obj)
ltbr.obj <- RunUMAP(ltbr.obj, dims = 1:30)
DimPlot(ltbr.obj, group.by="construct")

# Identify LTBR DEGs from the OverCITE-seq data
Idents(ltbr.obj) <- "construct"
ltbr.paper.de <- FindMarkers(ltbr.obj, ident.1="LTBR", ident.2 = "NGFR", only.pos=T)
gene.sets <- list("A"=rownames(ltbr.paper.de)[2:51]) # the top DEG is LTBR itself, so we skip it

# Add module scores from OverCITE-seq DEGs to both objects
obj <- AddModuleScore(obj, features=gene.sets)
ltbr.obj <- AddModuleScore(ltbr.obj, features=gene.sets)

# Adjust the idents for plotting
obj$construct_tmp <- obj$construct
obj$construct_tmp <- gsub("Saguaro.v1.", "", obj$construct_tmp)
obj$construct_tmp <- gsub("Control.", "", obj$construct_tmp)
Idents(obj) <- "construct_tmp"
ltbr.obj <- SetIdent(
  object = ltbr.obj,
  value   = factor(Idents(ltbr.obj), levels = c("LTBR", "NGFR"))
)

# Compare an LTBR signature in our data vs. the Legut et al data
p1 <- VlnPlot(obj, features="Cluster1", idents=c('LTBR', "GFP"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Our data")
p2 <- VlnPlot(ltbr.obj, features="Cluster1", idents=c('LTBR', "NGFR"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Legut et al")

# Select genes for an LTBR signature derived from CRISPR-All-seq data
ltbr.de <- FindMarkers(
  obj, ident.1 = "Saguaro.v1.LTBR", ident.2 = "Control.GFP", latent.vars = "donor",
  test.use="LR", min.pct = 0.02, logfc.threshold=0)
ltbr.de <- ltbr.de[ltbr.de$avg_log2FC > 0,]
gene.sets <- list("A"=rownames(ltbr.de)[1:50])

# Add module scores from CRISPR-All-seq DEGs to both objects
obj <- AddModuleScore(obj, features=gene.sets)
ltbr.obj <- AddModuleScore(ltbr.obj, features=gene.sets)

# Compare an LTBR signature in our data vs. the Legut et al data
p3 <- VlnPlot(obj, features="Cluster1", idents=c('LTBR', "GFP"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Our data")
p4 <- VlnPlot(ltbr.obj, features="Cluster1", idents=c('NGFR', "LTBR"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Legut et al")
(p1 | p2) / (p3 | p4)

### MED12 KO comparison

# Identify MED12-KO markers using the CRISPR-All-seq data
Idents(obj) <- "construct"
med12.de <- FindMarkers(
  obj, ident.1 = "Cholla.v1.MED12", ident.2 = "Cholla.v1.Scrambled.gRNA.sequence",
  latent.vars = "donor", test.use="LR", min.pct = 0.1, logfc.threshold=0.25, only.pos=T)

# Load the MED12 KO data from Freitas et al (Perturb-seq)
med12.file <- "/home/hartmana_stanford_edu/projects/perturb_crispr_all/data/med12_paper_sc/GSM6568648_MED12_filtered_feature_bc_matrix.h5"
aavs1.file <- "/home/hartmana_stanford_edu/projects/perturb_crispr_all/data/med12_paper_sc/GSM6568647_AAVS1_filtered_feature_bc_matrix.h5"
med12.h5 <- Read10X_h5(med12.file)
aavs1.h5 <- Read10X_h5(aavs1.file) # control gRNA
med12.ko.obj <- CreateSeuratObject(counts=med12.h5)
med12.ko.obj$construct <- "MED12"
tmp.obj <- CreateSeuratObject(counts = aavs1.h5)
tmp.obj$construct <- "AAVS1"
med12.ko.obj <- merge(med12.ko.obj, tmp.obj)
med12.ko.obj <- NormalizeData(med12.ko.obj)
med12.ko.obj <- JoinLayers(med12.ko.obj)
Idents(med12.ko.obj) <- "construct"

# Identify MED12-KO markers using the Perturb-seq data
med12.freitas.de <- FindMarkers(
  med12.ko.obj, ident.1 = "MED12", ident.2 = "AAVS1", test.use="LR",
  min.pct = 0.1, logfc.threshold=0.25, only.pos=T)

# First pass module scoring of MED12-KO gene sets
gene.sets <- list("A"=rownames(med12.freitas.de)[1:25])
obj <- AddModuleScore(obj, features=gene.sets)
med12.ko.obj <- AddModuleScore(med12.ko.obj, features=gene.sets)
obj$construct_tmp <- obj$construct
obj$construct_tmp <- gsub("Cholla.v1.", "", obj$construct_tmp)
obj$construct_tmp <- gsub("Scrambled.gRNA.sequence", "Ctrl gRNA", obj$construct_tmp)
obj$construct_tmp <- as.character(obj$construct_tmp)
Idents(obj) <- "construct_tmp"
p1 <- VlnPlot(obj, features="Cluster1", idents=c('MED12', "Ctrl gRNA", "Control.GFP"), pt.size=0) + geom_boxplot(fill=NA, outliers = F) + ggtitle("Our data")
p2 <- VlnPlot(med12.ko.obj, features="Cluster1", idents=c('MED12', "AAVS1"), pt.size=0) + geom_boxplot(fill=NA, outliers = F) + ggtitle("Freitas et al")
gene.sets <- list("A"=rownames(med12.de)[1:25])
obj <- AddModuleScore(obj, features=gene.sets)
med12.ko.obj <- AddModuleScore(med12.ko.obj, features=gene.sets)
Idents(obj) <- "construct_tmp"
p3 <- VlnPlot(obj, features="Cluster1", idents=c('MED12', "Ctrl gRNA", "Control.GFP"), pt.size=0) + geom_boxplot(fill=NA, outliers = F) + ggtitle("Our data")
p4 <- VlnPlot(med12.ko.obj, features="Cluster1", idents=c('MED12', "AAVS1"), pt.size=0) + geom_boxplot(fill=NA, outliers = F) + ggtitle("Freitas et al")

# First pass visualization revealed non-T cells to filter out
med12.ko.obj <- NormalizeData(med12.ko.obj)
med12.ko.obj <- FindVariableFeatures(med12.ko.obj, selection.method = "vst", nfeatures = 2000)
med12.ko.obj <- ScaleData(med12.ko.obj)
med12.ko.obj <- RunPCA(med12.ko.obj)
med12.ko.obj <- FindNeighbors(med12.ko.obj)
med12.ko.obj <- FindClusters(med12.ko.obj, resolution=0.5)
med12.ko.obj <- RunUMAP(med12.ko.obj, dims = 1:30)
DimPlot(med12.ko.obj, group.by=c("construct", "seurat_clusters"))
Idents(med12.ko.obj) <- "seurat_clusters"
med12.paper.de <- FindAllMarkers(med12.ko.obj, only.pos=T)
med12.ko.obj <- RunAzimuth(med12.ko.obj, reference = "pbmcref")
DimPlot(med12.ko.obj, group.by = "predicted.celltype.l1")

# Subset to the T and NK cell annotations and rerun the analysis
# Sometimes T cells are annotated as NK cells, so we include both classifications
Idents(med12.ko.obj) <- "predicted.celltype.l1"
cells <- WhichCells(med12.ko.obj, idents=c("CD4 T", "CD8 T", "NK"))
med12.ko.obj.ss <- subset(med12.ko.obj, cells=cells)

# Identify markers once more on the filtered dataset
Idents(med12.ko.obj.ss) <- "construct"
med12.freitas.de <- FindMarkers(
  med12.ko.obj.ss, ident.1 = "MED12", ident.2 = "AAVS1",
  test.use="LR", min.pct = 0.1, logfc.threshold=0.25, only.pos=T)

# Add MED12-KO module score from Perturb-seq gene set
gene.sets <- list("A"=rownames(med12.freitas.de)[1:50])
obj <- AddModuleScore(obj, features=gene.sets)
med12.ko.obj.ss <- AddModuleScore(med12.ko.obj.ss, features=gene.sets)

# Clean up idents for plotting
obj$construct_tmp <- obj$construct
obj$construct_tmp <- gsub("Cholla.v1.", "", obj$construct_tmp)
obj$construct_tmp <- gsub("Scrambled.gRNA.sequence", "Ctrl gRNA", obj$construct_tmp)
obj$construct_tmp <- as.character(obj$construct_tmp)
Idents(obj) <- "construct_tmp"

# Compare a MED12 signature in our data vs. the Freitas et al data
p1 <- VlnPlot(obj, features="Cluster1", idents=c('MED12', "Control.GFP"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Our data")
p2 <- VlnPlot(med12.ko.obj.ss, features="Cluster1", idents=c('MED12', "AAVS1"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Freitas et al")

# Add MED12-KO module score from CRISPR-All-seq gene set
gene.sets <- list("A"=rownames(med12.de)[1:50])
obj <- AddModuleScore(obj, features=gene.sets)
med12.ko.obj.ss <- AddModuleScore(med12.ko.obj.ss, features=gene.sets)

p3 <- VlnPlot(obj, features="Cluster1", idents=c('MED12', "Control.GFP"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Our data")
p4 <- VlnPlot(med12.ko.obj.ss, features="Cluster1", idents=c('MED12', "AAVS1"), pt.size=0) +
  geom_boxplot(fill=NA, outliers = F) + ggtitle("Freitas et al")


################################################################################
############################### DEGs plot ######################################
################################################################################
plot.gene.directions.stacked <- function(res) {
  colors <- list(
    "Cholla" = "#FF5A32",
    "Torch" = "#D88D21",
    "Saguaro" = "#4A6F28",
    "Gene" = "#4A6F28",
    "Kingcup" = "#6FB5D9",
    "Other" = "#CDD5E0"
  )
  df <- res$top[res$top$avg_log2FC > 0,]
  family_counts <- df %>% 
    dplyr::count(cluster, name = "n_degs")  
  family_counts$category <- "construct"
  family_counts$cluster <- gsub("Saguaro", "Gene", family_counts$cluster)
  family_counts$cluster <- gsub("PricklyPear", "Gene", family_counts$cluster)
  print(head(family_counts))
  family_counts$cluster <- factor(family_counts$cluster, levels = c(
    family_counts$cluster[str_starts(family_counts$cluster, "Cholla")][order(family_counts$n_degs[str_starts(family_counts$cluster, "Cholla")], decreasing = TRUE)],
    family_counts$cluster[str_starts(family_counts$cluster, "Gene")][order(family_counts$n_degs[str_starts(family_counts$cluster, "Gene")], decreasing = TRUE)],
    family_counts$cluster[str_starts(family_counts$cluster, "Kingcup")][order(family_counts$n_degs[str_starts(family_counts$cluster, "Kingcup")], decreasing = TRUE)],
    family_counts$cluster[str_starts(family_counts$cluster, "Torch")][order(family_counts$n_degs[str_starts(family_counts$cluster, "Torch")], decreasing = TRUE)]
  ))
  
  total_genes <- sum(family_counts$n_degs)
  family_counts <- family_counts %>%
    mutate(share = n_degs / total_genes * 100) %>%
    mutate(label = ifelse(share > 1.5, as.character(cluster), NA)) %>%
    mutate(family = case_when(
      str_starts(cluster, "Cholla") ~ "Cholla",
      str_starts(cluster, "Gene") ~ "Gene",
      str_starts(cluster, "Kingcup") ~ "Kingcup",
      str_starts(cluster, "Torch") ~ "Torch",
      TRUE ~ "Other"
    )) %>%
    arrange(family, n_degs)  # Order by family and number of DE genes in ascending order
  
  p <- ggplot(family_counts,
              aes(x = category,
                  y = n_degs,
                  fill = family)) +
    geom_col(width = 0.8, color = "black") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
    labs(title = "Top DEGs per construct (Up vs. Down)",
         x = "Construct",
         y = "# DEGs") +
    L_border() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_fill_manual(values = colors)
  return(p)
}

normalized.deg.plot<- function(obj, ncell, min.pct, logfc.threshold, donors) {
  # subset cells by donors
  Idents(obj) <- "donor"
  obj <- subset(obj, cells=WhichCells(obj, idents=donors))
  ss <- obj
  
  # find markers on subsetted cells
  Idents(ss) <- "construct"
  ss <- NormalizeData(ss)
  constructs <- unique(ss$construct)
  options(mc.cores = 15) 
  markers <- FindAllMarkers(
    ss,
    only.pos = FALSE,
    logfc.threshold = logfc.threshold,
    max.cells.per.ident = ncell,
    min.cells.group = 50,
    min.pct = min.pct, parallel = TRUE,)
  
  markers <- markers[markers$p_val_adj < 0.05,]
  
  best_hit <- markers %>%
    group_by(gene) %>%
    slice_min(p_val_adj, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  return(list(top=best_hit, all=markers))
}
res <- normalized.deg.plot(obj, 500, 0.01, 0.5, c("dRL320", "dRL321", "dRL322"))
p <- plot.gene.directions.stacked(res)
p
