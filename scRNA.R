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

### run demultiplex2 ###
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

### dimensional reduction ###
# standard analysis
obj <- JoinLayers(obj)
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)
obj <- RunUMAP(obj, dims = 1:30)

# donor-integrated
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

saveRDS(obj, file = "/home/hartmana_stanford_edu/perturb_crispr_all/output/seurat/04_seurat_object.rds")
