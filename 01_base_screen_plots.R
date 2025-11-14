# Script to regenerate plots of the base CACTUS screen
# Download the 'base_screen_amplicon_counts_matrix.csv' file from Zenodo:
# https://zenodo.org/records/17594332 

library(ggplot2)
library(BuenColors)
library(svglite)
library(patchwork)
library(DESeq2)
library(ggrepel)
library(data.table)
library(dplyr)
library(patchwork)
library(forcats)
library(tidyr)
library(stringr)
library(tibble)

# Load in input matrix available on Zenodo
cactus.base.merged.df <- read.csv("/home/hartmana_stanford_edu/projects/perturb_crispr_all/github/crispr-all/base_screen_amplicon_counts_matrix.csv", row.names=1)

type.colors <- list(
  "Torch-KD" = "firebrick", "Cholla-KO" = "#E7B800",
  "Saguaro-Over" = "royalblue", "PricklyPear-Syn" = "darkgreen",
  "Cardon-CAR" = "purple", "Fishhook-Binder" = "#FC4E07",
  "Kingcup-Signaling" = "#00AFBB")
donor.colors <- list(
  "d294" = "#8B0000", "d295" = "#B38600", "d296" = "#00008B",
  "d300" = "#006400", "d301" = "#4B0082", "d312" = "#B22222",
  "d314" = "#007373", "d318" = "black")

# Input coverage plots by sub-library
df <- cactus.base.merged.df
df <- df[grepl("input", colnames(df))]
classes <- c("Torch", "Cholla", "Saguaro", "PricklyPear", "Cardon", "Fishhook", "Kingcup")
plots <- list()
for (c in classes) {
  input <- colnames(df)[grepl("input", colnames(df))]
  print(head(df))
  plot.ss <- sweep(df[,input], 2, colSums(df[,input]),`/`) * 100
  m <- reshape2::melt(plot.ss[grepl(c, rownames(plot.ss)),])
  m$donor <- substr(m$variable, 1, 4)
  # ss <- m[m$value < 5000,]
  p <- ggplot(data=m, aes(value)) + geom_density(position="identity", alpha=0.5, color="black")+ scale_x_log10() +
    pretty_plot() +
    xlab("% of library") + ylab("density") + theme(axis.text=element_text(color="black", size=15), text=element_text(size=20)) +
    ggtitle(paste0("Readcount distribution at input for ", c))
  plots[[c]] <- p
}
wrap_plots(plots)

# Run DESeq2
# Note: in the library multiple 'backbone barcodes' are included which
# aids in statistical confidence within donors. So the original counts matrix
# contains a row of counts for each construct+backbone barcode combination
# rather than just constructs as is the case here. Counts across backbone
# barcodes are summed here.
# The volcano plot displayed in the paper uses the backbone barcode matrix
# rather than the matrix on Zenodo. Please email hartmana@stanford.edu
# if interested in the other matrix.
alias <- colnames(cactus.base.merged.df)
backbone <- substr(colnames(cactus.base.merged.df), nchar(colnames(cactus.base.merged.df)) - 2, nchar(colnames(cactus.base.merged.df)))
donor <- substr(colnames(cactus.base.merged.df), 1, 4)
conditions <- sub(".*_(input|output|acute)(_.*|$)", "\\1", colnames(cactus.base.merged.df))
desseq.df <- data.frame(Alias=alias, Condition=conditions, Donor=donor, Backbone=backbone)
deseq.tmp <- DESeqDataSetFromMatrix(countData = cactus.base.merged.df, colData = desseq.df, design = ~ Donor + Condition)
deseq.tmp <- DESeq(deseq.tmp)
counts <- counts(deseq.tmp, normalized=TRUE)
deseq.tmp <- results(deseq.tmp, contrast = c("Condition", "output", "input"))
deseq.tmp <- as.data.frame(deseq.tmp)
cactus.base.deseq.df <- data.frame(
  Construct <- row.names(deseq.tmp),
  L2FC <- deseq.tmp$log2FoldChange,
  adjPval <- deseq.tmp$padj)
colnames(cactus.base.deseq.df) <- c("Construct", "deseq_end_log2fc", "deseq_end_padj")
cactus.base.deseq.df$negLog10adjPval <- -log10(cactus.base.deseq.df$deseq_end_padj)
cactus.base.deseq.df$Type <- as.character(sapply(strsplit(cactus.base.deseq.df$Construct, "_"), `[`, 1))
cactus.base.deseq.df$Gene <- as.character(sapply(strsplit(cactus.base.deseq.df$Construct, "_"), `[`, 2))
cactus.base.deseq.df$Type <- factor(cactus.base.deseq.df$Type, levels=c("Fishhook-Binder", "Kingcup-Signaling", "Torch-KD", "Cholla-KO", "Cardon-CAR", "Saguaro-Over", "PricklyPear-Syn"))

# Plot log2fc by perturbation type
ggplot(cactus.base.deseq.df, aes(deseq_end_log2fc, Type, fill=Type)) + geom_boxplot() + geom_jitter() + pretty_plot() + scale_fill_manual(values=type.colors) + theme(text=element_text(size=20))

# Plot volcano
base.volcano <- ggplot(data=cactus.base.deseq.df, aes(x = deseq_end_log2fc, y = negLog10adjPval, fill=Type)) + 
  geom_point(alpha=1, size=2, shape=21, color="black") + 
  geom_text_repel(
    aes(label=ifelse((negLog10adjPval > 2 & deseq_end_log2fc > 0), Gene, "")),
    color="black", max.overlaps=1000) +
  scale_fill_manual(values=type.colors) + 
  geom_vline(xintercept=c(0), col="grey", linetype="dashed") +
  L_border() +
  xlab("log2FC") + ylab("-log10(p-value)") + theme(
    axis.text=element_text(color="black"),
    text=element_text(size=20))
base.volcano

# Highlight FDA CARs on volcano
p <- ggplot(data=cactus.base.deseq.df, aes(x = deseq_end_log2fc, y = negLog10adjPval)) + 
  geom_point(alpha=1, size=2, color="darkgrey") + 
  geom_point(data=cactus.base.deseq.df[cactus.base.deseq.df$Type == "Cardon-CAR",], aes(x = deseq_end_log2fc, y = negLog10adjPval), fill=c("#2a9d8f", "#e9c46a", "#264653", "#f4a261", "#e76f51"), alpha=1, size=3.5, shape=21, color="black") +
  geom_text_repel(
    data=cactus.base.deseq.df[cactus.base.deseq.df$Type == "Cardon-CAR",],
    aes(label=Gene), size=6,
    color="black", max.overlaps=1000) +
  geom_vline(xintercept=c(0), col="grey", linetype="dashed") +
  L_border() +
  xlab("log2FC") + ylab("-log10(p-value)") + theme(
    axis.text=element_text(color="black"),
    text=element_text(size=20)) +
  ylim(0, 15)
p

# TCF7s highlighted on volcano
tcf7s.df <- cactus.base.deseq.df[grepl("TCF7", cactus.base.deseq.df$Gene),]
p <- ggplot(data=cactus.base.deseq.df, aes(x = deseq_end_log2fc, y = negLog10adjPval)) + 
  geom_point(alpha=0.2, size=2, color="darkgrey") + 
  geom_point(data=tcf7s.df, aes(x = deseq_end_log2fc, y = negLog10adjPval), fill=c("#2a9d8f", "#e9c46a", "#e76f51"), alpha=1, size=3.5, shape=21, color="black") +
  geom_text_repel(
    data=tcf7s.df,
    aes(label=Gene), size=6,
    color="black", max.overlaps=1000) +
  L_border() +
  xlab("log2FC") + ylab("-log10(p-value)") + theme(
    axis.text=element_text(color="black"),
    text=element_text(size=20))
p

# LTBR constructs highlighted on volcano
ltbr.df <- cactus.base.deseq.df[grepl("LTBR", cactus.base.deseq.df$Gene),]
p <- ggplot(data=cactus.base.deseq.df, aes(x = deseq_end_log2fc, y = negLog10adjPval)) + 
  geom_point(alpha=0.2, size=2, color="darkgrey") + 
  geom_point(data=ltbr.df, aes(x = deseq_end_log2fc, y = negLog10adjPval), fill=c("#2a9d8f", "#e9c46a"), alpha=1, size=3.5, shape=21, color="black") +
  geom_text_repel(
    data=tcf7s.df,
    aes(label=Gene), size=6,
    color="black", max.overlaps=1000) +
  L_border() +
  xlab("log2FC") + ylab("-log10(p-value)") + theme(
    axis.text=element_text(color="black"),
    text=element_text(size=20))
p

# Perturbation-type highlights on volcano
plots <- list()
for (pt in names(type.colors)) {
  p<-ggplot() +
    geom_point(data=cactus.base.deseq.df, aes(x = deseq_end_log2fc, y = negLog10adjPval), color="darkgrey", size=1, alpha=1) +
    geom_point(
      data=cactus.base.deseq.df[cactus.base.deseq.df$Type == pt,],
      aes(x = deseq_end_log2fc, y = negLog10adjPval),
      fill=type.colors[[pt]], alpha=1, size=2, shape=21, color="black") +
    xlab("log2FC") + ylab("-log10(p-value)") + L_border()
  plots[[pt]] <- p
}
wrap_plots(plots)

# Donor enrichment plots
construct.donor.plot <- function(df, genes, type) {
  log2fc.df <- as.data.frame(df)
  donor <- unique(substr(colnames(df), 1, 4))
  conditions <- unique(sub(".*_(input|output|acute)(_.*|$)", "\\1", colnames(df)))
  for (d in donor) {
    input <- paste0(d, "_input")
    output <- paste0(d, "_output")
    new.column.name <- paste0(d, "_log2fc")
    new.column.vals <- log2( (df[,output] + 1) / (df[,input] + 1) )
    log2fc.df[new.column.name] <- unlist(new.column.vals)
  }
  log2fc.df <- log2fc.df[, grepl("log2fc$", colnames(log2fc.df))]
  m <- reshape2::melt(log2fc.df)
  overall_density <- density(m$value)
  density.df <- data.frame(x = overall_density$x, y = overall_density$y)
  constructs <- c()
  log2fcs <- c()
  for (g in genes) {
    constructs <- c(constructs, rep(g, ncol(log2fc.df)))
    log2fcs <- c(log2fcs, unlist(log2fc.df[g,]))
  }
  donor <- substr(colnames(log2fc.df), 1, 4)
  constructs <- data.frame(
    "Donor"=donor,
    "Construct"=constructs,
    "Log2FC"=log2fcs)
  constructs$Construct <- as.character(sapply(strsplit(constructs$Construct, "_"), `[`, 2))
  p<-ggplot() +
    geom_tile(data = density.df, aes(x = x, y = 0, fill = y), height = 1, alpha = 1) +
    geom_vline(data = constructs, aes(xintercept = Log2FC), color=type.colors[[type]], alpha = 1, linewidth=1) +
    facet_wrap(~ Construct, ncol = 1, strip.position = "left") +
    scale_fill_gradient(low = "white", high = "black") +
    scale_color_manual(values=donor.colors) +
    L_border() +
    theme(strip.text.y = element_text(angle = 0),
          text=element_text(size=20, color="black"),
          axis.text=element_text(color="black"),
          axis.text.y = element_blank(), legend.position = "none",
          strip.text = element_text(size = 8)) +
    labs(x = "Log2 Fold Change", y = NULL, fill = "Density")
  return(p)  
}

donor.plots <- list()
for (t in unique(cactus.base.deseq.df$Type)) {
  ss <- cactus.base.deseq.df[cactus.base.deseq.df$Type == t, ]
  # get top 6 genes by pvalue
  ss <- ss[ss$deseq_end_log2fc > 0,]
  genes <- head(ss[order(-ss$negLog10adjPval),"Construct"], 6)
  # genes <- head(ss[order(-ss$deseq_end_log2fc),"Construct"], 6)
  
  print(genes)
  p <- construct.donor.plot(counts, genes, t)
  donor.plots[[t]] <- p
}
wrap_plots(donor.plots)
