# This script re-analyzes the OverCITE-seq data from Legut et al (2022) and the Perturb-seq data from Freitas et al (2022)
# in order to compare the LTBR overexpression signature and the MED12 KO signatures in the CRISPR-All-seq dataset

library(Seurat)
library(ggplot2)
library(Azimuth)
library(BuenColors)
library(patchwork)
library(dplyr)
library(deMULTIplex2)

# Load CRISPR-All-seq Seurat object
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
