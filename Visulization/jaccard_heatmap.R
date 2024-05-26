library(tidyverse)
library(ggplot2)
library(pheatmap)
rm(list=ls())
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"

df_jaccard <- read.delim("G:/Micro_Maxquant_v2/results/All.Peptides.jaccard_dist.tsv", header=FALSE, sep="\t")
df_jaccard_NonHuman <- read.delim("G:/Micro_Maxquant_v2/results/All.NonHUman.Peptides.jaccard_dist.tsv", header=FALSE, sep="\t")
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"
jaccard <- df_jaccard[-c(1:3),-c(1:2)]
jaccard_NonHuman <- df_jaccard_NonHuman[-c(1:3),-c(1:2)]
# rownames(jaccard) <- paste(df_jaccard[-c(1:3),1], df_jaccard[-c(1:3),2], sep='_')
# colnames(jaccard) <- paste(df_jaccard[1,-c(1:2)], df_jaccard[2,-c(1:2)], sep='_')
# rownames(jaccard_NonHuman) <- paste(df_jaccard_NonHuman[-c(1:3),1], df_jaccard_NonHuman[-c(1:3),2], sep='_')
# colnames(jaccard_NonHuman) <- paste(df_jaccard_NonHuman[1,-c(1:2)], df_jaccard_NonHuman[2,-c(1:2)], sep='_')
jaccard_final <- jaccard
# jaccard_final[upper.tri(jaccard_final)] <- jaccard_NonHuman[upper.tri(jaccard_NonHuman)]
jaccard_final <- apply(as.matrix(jaccard_final), 2, as.numeric)
diag(jaccard_final) <- NA
colnames(jaccard_final) <- 1:20
rownames(jaccard_final) <- 1:20
p_totalLigand_heatmap <- pheatmap(jaccard_final,
                                  scale = "none", 
                                  cluster_rows = T,
                                  cluster_cols = T,
                                  cellwidth = 8,
                                  cellheight = 8,
                                  border_color = "white",
                                  display_numbers = F,
                                  show_rownames = F,
                                  show_colnames = F,
                                  annotation_row = data.frame(Sample = factor(df_jaccard[-c(1:3),1]),
                                                              type = factor(df_jaccard[-c(1:3),2])),
                                  annotation_col = data.frame(Sample = factor(df_jaccard[-c(1:3),1]),
                                                              type = factor(df_jaccard[-c(1:3),2])),
                                  annotation_colors = list(
                                      Sample = c(CRC01 = "#9e0142", CRC02 = "#d53e4f", CRC03 = "#f46d43",
                                                 CRC04 = "#fdae61", CRC05 = "#fee08b", CRC06 = "#e6f598",
                                                 CRC07 = "#abdda4", CRC08 = "#66c2a5", CRC09 = "#3288bd",
                                                 CRC10 = "#5e4fa2"),
                                      type = c(normal = "#33cccc", tumor  = "#F3756D")
                                      ),
                                  # output file
                                  # filename = paste(my_outdir, "Total_Ligand.heatmap.pdf", sep = "/"),
                                  # width = 110/25.4,
                                  # height = 110/25.4,
                                  )
pdf(paste(my_outdir, "Total_Ligand.heatmap.pdf", sep = "/"), height = 110/25.4, width=110/25.4)
print(p_totalLigand_heatmap)
dev.off()


jaccard_final <- jaccard_NonHuman
# jaccard_final[upper.tri(jaccard_final)] <- jaccard_NonHuman[upper.tri(jaccard_NonHuman)]
jaccard_final <- apply(as.matrix(jaccard_final), 2, as.numeric)
diag(jaccard_final) <- NA
colnames(jaccard_final) <- 1:20
rownames(jaccard_final) <- 1:20
p_nonHumanLigand_heatmap <- pheatmap(jaccard_final,
                                  scale = "none", 
                                  cluster_rows = T,
                                  cluster_cols = T,
                                  cellwidth = 8,
                                  cellheight = 8,
                                  border_color = "white",
                                  display_numbers = F,
                                  show_rownames = F,
                                  show_colnames = F,
                                  annotation_row = data.frame(Sample = factor(df_jaccard[-c(1:3),1]),
                                                              type = factor(df_jaccard[-c(1:3),2])),
                                  annotation_col = data.frame(Sample = factor(df_jaccard[-c(1:3),1]),
                                                              type = factor(df_jaccard[-c(1:3),2])),
                                  annotation_colors = list(
                                      Sample = c(CRC01 = "#9e0142", CRC02 = "#d53e4f", CRC03 = "#f46d43",
                                                 CRC04 = "#fdae61", CRC05 = "#fee08b", CRC06 = "#e6f598",
                                                 CRC07 = "#abdda4", CRC08 = "#66c2a5", CRC09 = "#3288bd",
                                                 CRC10 = "#5e4fa2"),
                                      type = c(normal = "#33cccc", tumor  = "#F3756D")
                                      ),
                                  # output file
                                  # filename = paste(my_outdir, "NonHuman_Ligand.heatmap.pdf", sep = "/"),
                                  # width = 110/25.4,
                                  # height = 110/25.4,
                                  )

pdf(paste(my_outdir, "NonHuman_Ligand.heatmap.pdf", sep = "/"), height = 110/25.4, width=110/25.4)
print(p_nonHumanLigand_heatmap)
dev.off()

df_HLA <- read.delim("G:/Micro_Maxquant_v2/results/CRC.all.HLAType.tsv", header = T, sep = "\t")
df_HLA %>% pivot_longer(cols = everything()) %>% mutate(n = 1) %>%
    ggplot(aes(x = name, y=value, fill = n)) +
    geom_tile()

df_HLA_jaccard <- read.delim("G:/Micro_Maxquant_v2/results/CRC.HLA.jaccard.tsv", header=T, row.names = 1, sep="\t")
rownames(df_HLA_jaccard)
HLA_jaccard <- as.matrix(apply(df_HLA_jaccard, 2, as.numeric))
row.names(HLA_jaccard) <- row.names(df_HLA_jaccard)
p_HLA <- pheatmap(HLA_jaccard,
         scale = "none", 
         cluster_rows = T,
         cluster_cols = T,
         # cellwidth = 8,
         # cellheight = 8,
         border_color = "white",
         display_numbers = F,
         show_rownames = T,
         show_colnames = T,
         )
pdf(paste(my_outdir, "HLA_jaccard.heatmap.pdf", sep = "/"))
print(p_HLA)
dev.off()
