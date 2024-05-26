getwd()
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(scales)
library(ggsignif)
library(cowplot)
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"

# df_MicroOnly <- read.delim("../../results/Micro_Only_Pep.tsv", header = TRUE, sep = "\t")
df_Peptides <- read.delim("../../results/All.Peptides.filtered.tsv.gz", header=TRUE, sep="\t")
read.delim("../../results/Micro_Only_Pep.tsv", header = TRUE, sep = "\t") %>%
    mutate(Sample = gsub("CRC", "#", Sample)) -> df_MicroOnly 

df_MicroOnly %>% select(Sequence, Sample, Seq_Sample_Type) %>%
    tidyr::separate_longer_delim(cols=Sample, delim = '||') -> df_MicroOnly2Sample

df_MicroOnly %>% select(Sequence, Species, Seq_Sample_Type) %>%
    tidyr::separate_longer_delim(cols=Species, delim = '||') -> df_MicroOnly2Species

# Fig4A

df_Peptides %>% 
    select(Sequence, Sample, Species) %>%
    mutate(Sample = gsub("CRC","#", Sample)) %>%
    filter(Sequence %in% df_MicroOnly$Sequence) %>%
    tidyr::separate_longer_delim(cols=Species, delim=',') %>%
    distinct() %>%
    group_by(Species, Sample) %>%
    summarise(Seq_Count = n()) %>%
    ungroup() %>% 
    tidyr::pivot_wider(names_from = Sample, values_from = Seq_Count) %>%
    tidyr::pivot_longer(cols = starts_with('#'), names_to = "Sample", values_to = "Seq_Count") %>%
    mutate(Seq_Count = factor(Seq_Count, levels=unique(sort(Seq_Count, decreasing=TRUE)))) -> df_MicroOnly_Specices2Sample

ggplot(df_MicroOnly_Specices2Sample, aes(x=Sample,y=Species)) +
    geom_tile(aes(fill=Seq_Count), colour="#dbdcdc") +
    scale_fill_manual(values = c('#a50026','#d73027','#f46d43','#fdae61',
                                 '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
                                 '#66bd63','#1a9850'),
                      na.value = '#7e7e7f',
    ) +
    labs(x="",y="",fill="") +
    scale_x_discrete(position="top") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
          axis.text.x = element_text(angle = 0, vjust=0.5, hjust=0.5),
          axis.ticks.x.top = element_line(color="black",linewidth = 1),
          legend.position = "right",
    ) -> p_Pep2Species_heatmap

df_MicroOnly2Sample %>%
    mutate(Reported = 1) %>%
    tidyr::pivot_wider(names_from = Sample, values_from = Reported) %>%
    tidyr::pivot_longer(cols = starts_with('#'), names_to = "Sample", values_to = "Reported") %>%
    mutate(Reported = factor(Reported, levels=unique(sort(Reported, decreasing=TRUE)))) %>%
    # filter(Sequence %in% unlist(select(filter(df_MicroOnly, Seq_Sample_Type=='tumor', Seq_Sample=='Multi_Sample'), Sequence))) %>%
    filter(Sequence %in% unlist(select(filter(df_MicroOnly, Seq_Sample=='Multi_Sample'), Sequence))) %>%
    left_join(distinct(select(df_MicroOnly, Sequence, Species)), by = "Sequence") %>%
    mutate(Label = paste(Seq_Sample_Type, Sequence, sep = " | ")) %>%
    ggplot(aes(x=Sample,y=Label)) +
    geom_tile(aes(fill=Reported), colour="#dbdcdc") +
    scale_fill_manual(values = c('#66bd63'), na.value = '#7e7e7f') +
    labs(x="",y="",fill="") +
    scale_x_discrete(position="top") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
          axis.text.x = element_blank(),
          axis.ticks = element_line(color="black", linewidth = 1),
          axis.ticks.x.top = element_blank(),
          legend.position = "right",
    ) -> p_Sequence_Seq2Sample

ggplot(df_MicroOnly2Species, aes(x=Species, fill=Seq_Sample_Type)) +
    geom_bar() +
    labs(x="", y="Number of Peptides") +
    coord_flip() +
    scale_fill_manual(values=c('#f2f0ba','#33cccc', '#f3756d')) +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.x = element_line(color="black", linewidth=0.5),
        legend.position = "",
    ) -> p_PepNum2Speceis_barplot

pdf(paste(my_outdir, "TumorOnlyMicroSeq.heatmap.pdf", sep = "/"), height = 130/25.4, width=180/25.4)
# plot_grid(p_Pep2Species_heatmap, p_Sequence_Seq2Sample, align = 'v', 
#           axis = 'l', nrow = 2, byrow=F, rel_heights=c(2,1))
plot_grid(p_Pep2Species_heatmap, p_Sequence_Seq2Sample, p_PepNum2Speceis_barplot, NULL,
          align = 'h', axis = 'l', nrow = 2,
          byrow=F,
          rel_heights=c(2,1),
          rel_widths = c(7,3)
)
dev.off()

# pdf(paste(my_outdir, "TumorOnlyMicroSeq.number.bar.pdf", sep = "/"), height = 80/25.4, width=130/25.4)
# print(p_PepNum2Speceis_barplot)
# dev.off()


# Fig4B
df_MicroOnly %>% 
    mutate(Seq_Length = str_length(Sequence)) %>%
    group_by(Seq_Sample_Type, Seq_Length) %>%
    summarise(Count = n()) %>%
    mutate(Percent=Count/sum(Count)*100) %>%
    ungroup() %>%
    ggplot(aes(x=Seq_Length)) +
    geom_bar(aes(y=Percent, fill=Seq_Sample_Type), stat="identity") +
    labs(x="Length of Peptides", y="Percentage of Peptides (%)", fill="Sample Type") +
    scale_x_continuous(breaks=c(7:16)) +
    scale_fill_manual(values = c('#f2f0ba','#33cccc', '#f3756d')) +
    facet_grid(~Seq_Sample_Type) +
    coord_flip() +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color="black", linewidth = 0.5),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
        strip.text.x = element_blank(),
        legend.position = ""
    ) -> p_MicroOnly_LengthDist_bar

pdf("../../results/visulization/MicroOnly_Pep_length.bar.pdf", height = 112.2729/2/25.4, width=63.4677 /25.4)
print(p_MicroOnly_LengthDist_bar)
dev.off()

# df_MicroOnly %>% 
#     mutate(Seq_Length = str_length(Sequence)) %>%
#     ggplot(aes(x=Seq_Sample_Type, y=Seq_Length)) +
#     geom_violin(aes(fill=Seq_Sample_Type), color=NA) +
#     geom_boxplot(fill="white", alpha=.8, color="black",
#                  outlier.shape = NA, 
#                  show.legend = FALSE, width=.05, size=.1) +
#     geom_jitter(width =.05, size=1) +
#     scale_y_continuous(breaks=c(7:16)) +
#     scale_fill_manual(values = c('#f2f0ba','#33cccc', '#f3756d')) +
#     labs(x="Sample Type", y="Length of Peptides", fill="") +
#     theme(
#         panel.background = element_blank(),
#         panel.grid = element_blank(),
#         axis.line = element_line(color="black", linewidth=0.5),
#         axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = "top"
#     ) -> p_MicroOnly_LengthDist_violin

# Fig 4C
df_Affinity <- read.delim("../../results/All.affinity_features.txt", header = TRUE, sep="\t")
colnames(df_Affinity)
df_Affinity %>% 
    # select(-c(Length, BinderNum, Sample, Sample_Type)) %>% 
    rename(Sequence = X.Peptide, IC50 = Aff.IC50nM,
           Rank_EL = Rank_EL..., Rank_BA = Rank_BA...) %>%
    tidyr::separate_longer_delim(cols=c(MHC, IC50, Rank_EL, Rank_BA, BinderType), delim=',') %>%
    distinct() %>%
    mutate(BinderType = factor(BinderType,
                               levels=c("SB","WB","NB"),
                               # labels=c("Strong Binder","Weak Binder","No Binder"),
                               ),
           MHCType = factor(substring(MHC, 1, 5)),
           MHC = substring(MHC, 7)
           ) %>% as_tibble() -> df_Affinity_seperate

df_Affinity_seperate %>%
    group_by(Sample_Type) %>%
    distinct(Sequence) %>%
    mutate(Seq_Length = str_length(Sequence)) %>%
    ggplot(aes(x=Seq_Length, fill=Sample_Type)) +
    geom_bar()

df_Affinity_seperate %>%
    right_join(select(df_MicroOnly, Sequence, Seq_Sample_Type),
               by = join_by(Sequence),
               ) -> df_Affinity_tumorOnly_MicroPep
    
# df_Affinity_seperate %>%
#     group_by(MHCType, MHC, BinderType)%>%
#     summarise(Count=n()) %>%
#     ungroup() -> df_Affinity_Count
# 
# df_Affinity_seperate %>% group_by(Sequence) %>%
#     filter(Rank_EL == min(Rank_EL)) %>%
#     ungroup() -> df_Affinity_min
# 
# df_Affinity_min %>%
#     group_by(MHCType, MHC, BinderType)%>%
#     summarise(Count=n()) %>%
#     ungroup() -> df_Affinity_min.Count
# 
# df_Affinity_Count %>% 
#     group_by(MHCType, MHC) %>%
#     summarise(Total=sum(Count)) %>%
#     ungroup() %>%
#     arrange(MHCType, desc(Total)) %>% 
#     select(MHC) %>%
#     unlist() -> MHC.sort
# 
# df_Affinity_min.Count %>% 
#     group_by(MHCType, MHC) %>%
#     summarise(Total=sum(Count)) %>%
#     ungroup() %>%
#     arrange(MHCType, desc(Total)) %>% 
#     select(MHC) %>%
#     unlist() -> MHC.min.sort

df_Affinity_tumorOnly_MicroPep %>%
    group_by(Seq_Sample_Type, MHCType, MHC, BinderType)%>%
    summarise(Count=n()) %>%
    ungroup() -> df_Affinity_tumorOnly_MicroPep_Count

df_Affinity_tumorOnly_MicroPep %>% group_by(Sequence) %>%
    filter(Rank_EL == min(Rank_EL)) %>%
    ungroup() -> df_Affinity_tumorOnly_MicroPep_min

df_Affinity_tumorOnly_MicroPep_min %>%
    group_by(MHCType, MHC, BinderType)%>%
    summarise(Count=n()) %>%
    ungroup() -> df_Affinity_tumorOnly_MicroPep_min.Count

distinct(df_Affinity_tumorOnly_MicroPep_min, Sequence, .keep_all = TRUE) %>%
    group_by(Seq_Sample_Type, BinderType) %>%
    summarise(Count=n()) %>%
    group_by(Seq_Sample_Type) %>%
    mutate(Per = Count/sum(Count)*100) %>% 
    ungroup() %>%
    ggplot(aes(x=Seq_Sample_Type, fill=BinderType)) +
    geom_bar(aes(y=Per), stat = "identity") +
    labs(x="Sample Type", y="Percentage of Peptides (%)", fill="Binder Type") +
    scale_fill_manual(values=c('#c83e4c','#63b89d', '#3380b1')) +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color="black", linewidth=0.5),
        axis.ticks = element_line(color="black", linewidth=0.5),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
        legend.background = element_blank(),
        legend.position = "",
    ) -> p_MicroOnly_affinity_bar

pdf(paste(my_outdir, "TumorOnlyMicroSeq.BindingType.bar.pdf", sep = "/"), height = 112.2729/2/25.4, width=63.4677 /25.4)
print(p_MicroOnly_affinity_bar)
dev.off()
# IC50 violine plot
distinct(df_Affinity_tumorOnly_MicroPep_min, Sequence, .keep_all = TRUE) %>%
    mutate(IC50 = as.numeric(IC50)) %>%
    ggplot(aes(x=Seq_Sample_Type, y=IC50, fill=Seq_Sample_Type)) +
    geom_violin() +
    geom_boxplot(fill="white", alpha=.8, color="black",
                 outlier.shape = NA,
                 show.legend = FALSE, width=.05, size=.1) +
    geom_jitter(width =.05, size=1, show.legend = FALSE) +
    geom_signif(
        comparisons = list(c("BOTH","normal"),
                           c("BOTH","tumor"),
                           c("normal","tumor")
                           ),
        y_position = c(50000,55000,60000),
        map_signif_level = TRUE,
        textsize = 6,
        vjust = -0.5,
    ) +
    scale_fill_manual(values = c('#f2f0ba','#33cccc', '#f3756d')) +
    labs(x="",y="Affinity with HLA-I in IC50 (nM)") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.5, color="black"),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
        axis.ticks = element_line(linewidth = 0.5, color="black"),
        legend.position = "",
    ) -> p_MicroOnly_affinityIC50_min_violin
pdf(paste(my_outdir, "TumorOnlyMicroSeq.affinityIC50.violin.pdf", sep = "/"), height = 106.9377/25.4, width=106.9377/25.4)
print(p_MicroOnly_affinityIC50_min_violin)
dev.off()
# # Rank_EL violine plot
# distinct(df_Affinity_tumorOnly_MicroPep_min, Sequence, .keep_all = TRUE) %>%
#     mutate(Rank_EL = as.numeric(Rank_EL)) %>%
#     ggplot(aes(x=Seq_Sample_Type, y=Rank_EL, fill=Seq_Sample_Type)) +
#     geom_violin() +
#     geom_boxplot(fill="white", alpha=.8, color="black",
#                  outlier.shape = NA,
#                  show.legend = FALSE, width=.05, size=.1) +
#     geom_jitter(width =.05, size=1, show.legend = FALSE) +
#     geom_signif(
#         comparisons = list(c("BOTH","normal"),
#                            c("BOTH","tumor"),
#                            c("normal","tumor")
#         ),
#         y_position = c(105,110,115),
#         map_signif_level = TRUE,
#         textsize = 6,
#         vjust = -0.5,
#     ) +
#     scale_fill_manual(values = c('#f2f0ba','#33cccc', '#f3756d'))
# Rank_BA violine plot
# distinct(df_Affinity_tumorOnly_MicroPep_min, Sequence, .keep_all = TRUE) %>%
#     mutate(Rank_BA = as.numeric(Rank_BA)) %>%
#     ggplot(aes(x=Seq_Sample_Type, y=Rank_BA, fill=Seq_Sample_Type)) +
#     geom_violin() +
#     geom_boxplot(fill="white", alpha=.8, color="black",
#                  outlier.shape = NA,
#                  show.legend = FALSE, width=.05, size=.1) +
#     geom_jitter(width =.05, size=1, show.legend = FALSE) +
#     geom_signif(
#         comparisons = list(c("BOTH","normal"),
#                            c("BOTH","tumor"),
#                            c("normal","tumor")
#         ),
#         y_position = c(105,110,115),
#         map_signif_level = TRUE,
#         textsize = 6,
#         vjust = -0.5,
#     ) +
#     scale_fill_manual(values = c('#f2f0ba','#33cccc', '#f3756d'))
# df_Affinity_tumorOnly_MicroPep_Count %>%
#     group_by(MHCType, MHC) %>%
#     summarise(Total=sum(Count)) %>%
#     ungroup() %>%
#     arrange(MHCType, desc(Total)) %>%
#     select(MHC) %>%
#     unlist() -> MHC.sort
# 
# df_Affinity_tumorOnly_MicroPep_min.Count %>%
#     group_by(MHCType, MHC) %>%
#     summarise(Total=sum(Count)) %>%
#     ungroup() %>%
#     arrange(MHCType, desc(Total)) %>%
#     select(MHC) %>%
#     unlist() -> MHC.min.sort
# 
# ggplot(df_Affinity_tumorOnly_MicroPep_Count %>% mutate(MHC = factor(MHC, levels = MHC.sort)), aes(x=MHC, group=MHCType)) +
#     geom_bar(aes(y= Count, fill=BinderType), stat="identity") +
#     labs(x="", y="Number of Peptides", fill="Binder Type") +
#     scale_fill_manual(values = rev(c('#c83e4c','#63b89d', '#3380b1'))) +
#     theme(panel.background = element_blank(),
#           panel.grid = element_blank(),
#           axis.line = element_line(color="black", linewidth=0.5),
#           axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
#           axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
#           axis.ticks = element_line(color="black", linewidth = 0.5),
#           legend.position = "right",
#     )
# ggplot(df_Affinity_tumorOnly_MicroPep_min.Count %>% mutate(MHC = factor(MHC, levels = MHC.min.sort)), aes(x=MHC, group=MHCType)) +
#     geom_bar(aes(y= Count, fill=BinderType), stat="identity") +
#     labs(x="", y="Number of Peptides", fill="Binder Type") +
#     scale_fill_manual(values = rev(c('#c83e4c','#63b89d', '#3380b1'))) +
#     theme(panel.background = element_blank(),
#           panel.grid = element_blank(),
#           axis.line = element_line(color="black", linewidth=0.5),
#           axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
#           axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
#           axis.ticks = element_line(color="black", linewidth = 0.5),
#           legend.position = "right",
#     )
# 
# ggplot(df_Affinity_tumorOnly_MicroPep_min.Count %>% mutate(MHC = factor(MHC, levels = MHC.min.sort)), aes(x=MHCType)) +
#     geom_bar(aes(y= Count, fill=BinderType), stat="identity") +
#     labs(x="", y="Number of Peptides", fill="Binder Type") +
#     scale_fill_manual(values = rev(c('#ca3f4c','#64ba9f', '#3482b2'))) +
#     theme(panel.background = element_blank(),
#           panel.grid = element_blank(),
#           axis.line = element_line(color="black", linewidth=0.5),
#           axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
#           axis.text.x = element_text(angle=90, hjust=1, vjust = 0.5),
#           axis.ticks = element_line(color="black", linewidth = 0.5),
#           legend.position = "right",
#     )
# Fig 4D
library(circlize)

df_Peptides %>%
    filter(Sequence %in% df_F_nucleatu$Sequence) %>%
    tidyr::separate_longer_delim(cols = c(Proteins, Species), delim = ',') %>%
    select(Sequence, Proteins, Species) %>%
    filter(grepl("Fusobacterium nucleatum", Species)) %>%
    distinct(Sequence, Proteins) %>%
    inner_join(df_F_nucleatu, by=join_by(Sequence, Proteins)) -> df_F_nucleatu_filtered

# df_F_nucleatu %>%
df_F_nucleatu_filtered %>%
    mutate(Protein_Location = if_else(grepl("^$", Protein_Location), "Unknow", Protein_Location),
           Protein_Location = if_else(grepl("Cell membrane", Protein_Location), "Cell membrane", Protein_Location)) %>%
    distinct(Sequence, Proteins, .keep_all = TRUE) %>%
    select(Sequence, Protein_Location, Seq_Sample_Type) -> df_test

my.group.All <- c(unique(df_test$Protein_Location), unique(df_test$Seq_Sample_Type), unique(df_test$Sequence))
my.group <- structure(c(rep("Location", 4), rep("Type", 3), rep("Sequence", 20)), names=my.group.All)
# grid.col <- structure(c(rep(2,4), rep(3,3), rep(4,20)), names=my.group.All)
grid.col <- structure(c('#cb3228','#f4aa63','#d2e38b','#a1cd6a',
                        '#eeebb8','#e7726a','#44bbbd',
                        rep('#367ba8',20)), names=my.group.All)

df_test %>% group_by(Sequence, Protein_Location) %>% 
    summarise(Count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Protein_Location, values_from = Count, values_fill = 0) -> df_test.1
# mat1 <- as.matrix(select(df_test.1, -Sequence))
mat1 <- t(as.matrix(apply(select(df_test.1, -Sequence), MARGIN = 1, function(x){x/sum(x)})))
rownames(mat1) <- df_test.1$Sequence

df_test %>% group_by(Sequence, Seq_Sample_Type) %>% 
    summarise(Count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Seq_Sample_Type, values_from = Count, values_fill = 0) -> df_test.2
mat2 <- t(as.matrix(apply(select(df_test.2, -Sequence), MARGIN = 1, function(x){x/sum(x)})))
rownames(mat2) <- df_test.2$Sequence

df_test %>% group_by(Protein_Location, Seq_Sample_Type) %>% 
    summarise(Count = n()) %>%
    ungroup() %>%
    pivot_wider(names_from = Seq_Sample_Type, values_from = Count, values_fill = 0) -> df_test.3
# mat3 <- as.matrix(select(df_test.3, -Protein_Location))
mat3 <- matrix(0, nrow = nrow(df_test.3), ncol = ncol(df_test.3)-1)
rownames(mat3) <- df_test.3$Protein_Location
colnames(mat3) <- colnames(df_test.3)[-1]

mat = matrix(0, nrow = 24, ncol = 7)
rownames(mat) = c(rownames(mat2), rownames(mat3))
colnames(mat) = c(colnames(mat1), colnames(mat2))
mat[rownames(mat1), colnames(mat1)] = mat1
mat[rownames(mat2), colnames(mat2)] = mat2
mat[rownames(mat3), colnames(mat3)] = mat3
# mat
pdf(paste(my_outdir, "F_nucleatum.circular.pdf", sep = "/"), height = 82.7925  /25.4, width=67.7666  /25.4)
# chordDiagram(mat, group = my.group, grid.col = grid.col, 
#              annotationTrack = c("grid"),
#              preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(mat))))))
chordDiagram(mat, group = my.group, grid.col = grid.col, 
             annotationTrack = c("grid"),
             preAllocateTracks = list(track.height = 0.2))
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, 
                cex = 0.3,
                adj = c(0, 0),
                )
}, bg.border = NA)
# highlight.sector(rownames(mat3), track.index = 1, col = "grey50", 
#                  text = "Protein Location", cex = 0.8, text.col = "black", niceFacing = TRUE)
# highlight.sector(colnames(mat3), track.index = 1, col = "grey50", 
#                      text = "Sample Type", cex = 0.8, text.col = "black", niceFacing = TRUE)
circos.clear()
dev.off()

df_F_nucleatu_filtered %>% 
    select(Sequence, Proteins, Proteins_names) %>% 
    distinct(Sequence, Proteins, .keep_all = TRUE) %>%
    group_by(Proteins) %>%
    mutate(Count = n()) %>%
    ungroup() %>%
    arrange(desc(Count), Proteins_names) %>% 
    view()

df_MicroOnly %>% 
    filter(grepl("Elongation factor Tu",Proteins_names)) %>%
    view()

