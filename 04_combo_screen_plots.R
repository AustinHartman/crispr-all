# Script to regenerate plots of the combo CACTUS screen
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
library(ggpointdensity)

# Load in input matrix available on Zenodo
df <- read.csv("/home/hartmana_stanford_edu/projects/perturb_crispr_all/github/crispr-all/combo_screen_amplicon_counts_matrix.csv", row.names=1)
df <- as.data.frame(df)

### INPUT TIMEPOINT COVERAGE PLOT ###
input <- c("D323_Input", "D324_Input", "D326_Input", "D327_Input", "D328_Input")
df.input <- as.data.frame(df[,input])

m <- reshape2::melt(df[,input])
m$donor <- substr(m$variable, 1, 4)
ss <- m[m$value < 5000,]
p <- ggplot(data=ss, aes(value, fill=donor)) + geom_density(position="identity", alpha=0.5, color="black") +
  scale_fill_manual(values=donor.colors) + scale_x_log10() +
  geom_vline(xintercept = 50, linetype="dashed", color = "red", size=1) + pretty_plot() +
  xlab("Readcount") + ylab("# Constructs") + theme(axis.text=element_text(color="black", size=15), text=element_text(size=20)) +
  ggtitle("Readcount distribution at input")
p

plot.ss <- sweep(df[,input], 2, colSums(df[,input]),`/`) * 100
m <- reshape2::melt(plot.ss)
m$donor <- substr(m$variable, 1, 4)
p <- ggplot(data=m, aes(value, fill=donor)) + geom_density(position="identity", alpha=0.5, color="black") +
  scale_fill_manual(values=donor.colors) + scale_x_log10() +
  pretty_plot() +
  xlab("% Reacount") + ylab("# Constructs") + theme(axis.text=element_text(color="black", size=15), text=element_text(size=20)) +
  ggtitle("Readcount distribution at input")
p

# run deseq2
columns <- c(
  "D323_Input", "D323_Acute", "D323_Rep1Stim2", "D323_Rep2Stim2", "D323_Rep1Stim4", "D323_Rep2Stim4", "D323_Rep1Stim6", "D323_Rep2Stim6",
  "D324_Input", "D324_Acute", "D324_Stim2", "D324_Stim4", "D324_Stim6",
  "D326_Input", "D326_Acute", "D326_Rep1Stim2", "D326_Rep2Stim2", "D326_Rep1Stim4", "D326_Rep2Stim4", "D326_Rep1Stim6", "D326_Rep2Stim6",
  "D327_Input", "D327_Acute", "D327_Stim2", "D327_Stim4", "D327_Stim6",
  "D328_Input", "D328_Acute", "D328_Stim2", "D328_Stim4", "D328_Stim6")
# df.use <- df[df[,"D323_Input"] > 50 & df[,"D324_Input"] > 50 & df[,"D326_Input"] > 50 & df[,"D327_Input"] > 50 & df[,"D328_Input"] > 50, ]
df.use <- df
print(nrow(df))
print(nrow(df.use))
col.data <- data.frame(
  Alias = colnames(df.use),
  Condition = c(
    "Input", "Acute", "Stim2", "Stim2", "Stim4", "Stim4", "Stim6", "Stim6",
    "Input", "Acute", "Stim2", "Stim4", "Stim6",
    "Input", "Acute", "Stim2", "Stim2", "Stim4", "Stim4", "Stim6", "Stim6",
    "Input", "Acute", "Stim2", "Stim4", "Stim6",
    "Input", "Acute", "Stim2", "Stim4", "Stim6"
  ),
  Donor = c(
    "323", "323", "323", "323", "323", "323", "323", "323",
    "324", "324", "324", "324", "324",
    "326", "326", "326", "326", "326", "326", "326", "326",
    "327", "327", "327", "327", "327",
    "328", "328", "328", "328", "328"
  )
)
deseq2Data <- DESeqDataSetFromMatrix(countData = df.use, colData = col.data, design = ~ Donor + Condition)
deseq2Data <- DESeq(deseq2Data)
deSeq <- results(deseq2Data, contrast = c("Condition", "Stim6", "Input"))
deSeq <- as.data.frame(deSeq)

# plot volcano
deseq_end <- data.frame(Construct <- row.names(deSeq),L2FC <- deSeq$log2FoldChange, adjPval <- deSeq$padj, p <- deSeq$pvalue)
colnames(deseq_end) <- c("Construct", "deseq_end_log2fc", "deseq_end_padj", "p")
deseq_end$negLog10adjPval <- -log10(deseq_end$deseq_end_padj)
p <- ggplot(deseq_end, aes(x=deseq_end_log2fc, y=negLog10adjPval)) + 
  geom_point(aes(
    color=ifelse(deseq_end_log2fc>0 & negLog10adjPval>10, "1", "2"), 
    alpha=ifelse(deseq_end_log2fc>0 & negLog10adjPval>10, 1, 0.25)), size=3) + 
  L_border() +
  geom_text_repel(aes(label=ifelse(deseq_end_log2fc>0 & negLog10adjPval>10, Construct, "")), size=3) +
  xlab("log2(FC)") + ylab("-log10(p)") + 
  theme(text=element_text(size=20), legend.position = "none", axis.text=element_text(color="black")) + 
  scale_color_manual(values=c("firebrick", "black"))
p

# donor enrichment plots
df.use <- as.data.frame(counts(deseq2Data, normalize=TRUE))
df.dens <- data.frame(
  "D323" = log2((df.use$D323_Rep1Stim6+1) / (df.use$D323_Input+1)),
  "D324" = log2((df.use$D324_Stim6+1) / (df.use$D324_Input+1)),
  "D326" = log2((df.use$D326_Rep1Stim6+1) / (df.use$D326_Input+1)),
  "D327" = log2((df.use$D327_Stim6+1) / (df.use$D327_Input+1)),
  "D328" = log2((df.use$D328_Stim6+1) / (df.use$D328_Input+1)))
rownames(df.dens) <- rownames(df.use)
df.dens$avg <- rowMeans(df.dens[, 1:5])
# build density dfs
m <- mean(rowMeans(df.dens))
m <- 0
overall_density <- density(rowMeans(df.dens) - m)
# overall_density$y <- overall_density$y - m
density.df <- data.frame(x = overall_density$x, y = overall_density$y)
constructs <- c("Fas_MED12_OX40_CD28z", "Fas_MED12_OX40_delCD28z", "Fas_TOX2_GFP_delCD28z", "Fas_TOX2_BATF_TACI")
log2fcs <- c()
constsructs.list <- c()
for (c in constructs) {
  log2fcs <- c(log2fcs, unlist(df.dens[c,] - m))
  constsructs.list <- c(constsructs.list, c)
}
plot.data <- data.frame(
  Construct=constsructs.list,
  Log2FC = log2fcs)
p<-ggplot() +
  geom_tile(data = density.df, aes(x = x, y = 0, fill = y), height = 1, alpha = 1) +
  geom_vline(data = plot.data, aes(xintercept = Log2FC), color = "firebrick3", alpha = 0.75, size=1) +  # red lines for construct+donor pairs
  facet_wrap(~ Construct, ncol = 1, strip.position = "left") +  # One row per construct
  scale_fill_gradient(low = "white", high = "black") +
  L_border() +
  theme(strip.text.y = element_text(angle = 0),
        text=element_text(size=20, color="black"),
        axis.text=element_text(color="black"),
        axis.text.y = element_blank(), legend.position = "none",
        strip.text = element_text(size = 4)) +
  labs(x = "Log2 Fold Change", y = NULL, fill = "Density")
p
