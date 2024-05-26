getwd()
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(scales)
library(ggsignif)
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"

Micro_DB <- list.files("../../ProteinDB/")
Micro_DB <- Micro_DB[grep("selected_microbes.tsv", Micro_DB)]
df_Micro <- data.frame()
for (each in Micro_DB){
    sample <- str_split_1(each, '\\.')[1]
    df_tmp <- read.delim(paste("../../ProteinDB", each, sep="/"), header=T, sep="\t")
    df_tmp[['Sample']] <- sample
    if(nrow(df_Micro) == 0){
        df_Micro <- df_tmp
    }else{
        df_Micro <- bind_rows(df_Micro, df_tmp)
    }
}
df_Micro %>% dplyr::select(-Pro_DB_Dir) %>% dplyr::distinct() %>%
    dplyr::rename(Label = Species,
                  Species = Species_Name) %>%
    dplyr::mutate(Species = gsub("[^A-Za-z0-9]+", "_", Species))-> df_Micro
my.Lineage <- strsplit(df_Micro$Lineage, ';', fixed = TRUE)
df_Micro$genus <- gsub("unclassified\\s+", "", sapply(my.Lineage, function(x) x[[length(x)]]))
# Aquabacterium olei excluded from final result
df_Micro <- df_Micro[-grep("Aquabacterium",df_Micro$Species),]
df_BIC_genus <- read.delim("../../results/BIC_DB/BIC_genus_expression_matrix.txt", header=TRUE, sep="\t")
df_BIC_meta <- read.delim("../../results/BIC_DB/BIC_metadata_matrix.txt", header=TRUE, sep="\t")

df_BIC_genus %>% 
    pivot_longer(cols = where(is.numeric), names_to = "Aliquot_barcode_miRNAseq", values_to = "Rela_Abun") %>%
    mutate(Aliquot_barcode_miRNAseq = gsub("\\.","-",Aliquot_barcode_miRNAseq)) %>%
    left_join(df_BIC_meta[c('Project_name','Aliquot_barcode_miRNAseq')]) -> df_BIC_genus_pivort
colnames(df_BIC_genus_pivort)
P_BIC <- list()
for (i in unique(df_Micro$genus)){
    df_BIC_genus_pivort %>% filter(Genus == i) -> df_tmp
    if (nrow(df_tmp) == 0){ next }
    ggplot(df_tmp, aes(x=Project_name, y=Rela_Abun, fill=Project_name)) +
        geom_boxplot(outlier.size = 0.05) +
        labs(x="", y="", title=i) +
        theme(
            panel.background = element_blank(),
            panel.grid = element_blank(),
            plot.title = element_text(family = 'sans', face= "plain", color='black', size = 12, hjust=1, vjust=-5),
            plot.title.position = "plot",
            axis.line = element_line(linewidth = 0.5, colour = "black"),
            axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
            legend.position = "none"
        ) +
        coord_flip() -> P_BIC[[i]]
}

for (i in sort(names(P_BIC))){
    pdf(paste(my_outdir, "../BIC_DB", paste(i,"BIC_boxplot.pdf", sep="."), sep = "/"), height = 75/25.4, width=64/25.4)
    print(P_BIC[[i]])
    dev.off()
}
library(cowplot)
themes.noYaxis <- theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank())
P_BIC[[i]] + themes.noYaxis
pdf(paste(my_outdir, "Supplementary Fig 4.pdf", sep = "/"), height = 130/25.4, width=203/25.4)
cowplot::plot_grid(P_BIC[["Bacteroides"]],
                   P_BIC[["Blautia"]] + themes.noYaxis,
                   P_BIC[["Collinsella"]] + themes.noYaxis,
                   P_BIC[["Coprococcus"]] + themes.noYaxis,
                   P_BIC[["Eubacterium"]] + themes.noYaxis,
                   P_BIC[["Faecalibacterium"]],
                   P_BIC[["Parvimonas"]] + themes.noYaxis,
                   P_BIC[["Porphyromonas"]] + themes.noYaxis,
                   P_BIC[["Pseudomonas"]] + themes.noYaxis,
                   P_BIC[["Streptococcus"]] + themes.noYaxis,
                   align = "hv",
                   axis = 'lb',
                   nrow=2)

dev.off()
read.csv("../../results/All.Peptides.tsv.gz", header=T, sep="\t") %>% 
    mutate(Sample = gsub("CRC", "#", Sample)) -> df_Peptides
colnames(df_Peptides)
# Count number of species and proteins for each peptides
df_Peptides %>% tidyr::separate_longer_delim(cols = c(Proteins, Species), delim = ';') %>%
    dplyr::filter(!grepl("Aquabacterium", Species)) %>% 
    dplyr::group_by_at(colnames(df_Peptides)[!(grepl("^Proteins|Species$", colnames(df_Peptides)))]) %>%
    dplyr::summarise(
        Protein.number = length(unique(Proteins)),
        Species.number = length(unique(Species)),
        Proteins = paste(Proteins, collapse = ';'),
        Species = paste(Species, collapse = ';')
    ) %>%
    ungroup() -> df_Peptides_filtered
# df_Peptides %>%
df_Peptides_filtered %>%
    dplyr::select(Sequence, Proteins, Species) %>%
    tidyr::separate_longer_delim(cols=c(Proteins, Species), delim = ';') %>% 
    dplyr::filter(!grepl("^CON_", Proteins)) %>% 
    dplyr::distinct() %>%
    dplyr::group_by(Sequence) %>% 
    dplyr::summarise(
        Protein.number = length(unique(Proteins)),
        Species.number = length(unique(Species)),
        Proteins = paste(Proteins, collapse = ';'),
        Species = paste(Species, collapse = ';')) %>%
    dplyr::ungroup() -> df_Total_Pep
# df_Total_Pep %>% filter(grepl("Aquabacterium", Species), Species.number > 1) %>% view()
# classifying peptide according to the source protein
df_Total_Pep %>%
    dplyr::mutate(Origin_Type = if_else((grepl("Homo sapiens", Species) & Species.number == 1), "Human", "Micro")) -> df_Total_Pep
df_Total_Pep[(grepl("Homo sapiens", df_Total_Pep$Species) & (grepl("Micro", df_Total_Pep$Origin_Type))), 'Origin_Type'] <- "Common"
# Number of different classified peptides
df_Total_Pep %>%
    dplyr::group_by(Origin_Type) %>% dplyr::summarise(n=n()) %>%
    ungroup() %>% mutate(Per = round(n/sum(n)*100,2))
dplyr::filter(df_Total_Pep, Origin_Type == "Human") %>% dplyr::select(Sequence) -> seq_tmp
filter(df_Peptides_filtered, Sequence %in% unique(seq_tmp$Sequence)) -> df_tmp
my_group <- unique(df_tmp$Sample_Type)
venn_list = list()
venn_list[[my_group[1]]] <- unlist(distinct(df_tmp[df_tmp$Sample_Type == my_group[1],'Sequence']))
venn_list[[my_group[2]]] <- unlist(distinct(df_tmp[df_tmp$Sample_Type == my_group[2],'Sequence']))
library(VennDiagram)
venn.plot <- venn.diagram(
    venn_list[my_group],
    filename=NULL,
    # Output features
    # imagetype="svg",
    # height = 10,
    # width = 10,
    # resolution = 300,
    # compression = "lzw",
    # Circles
    # col=c('#cc3333','#0000ff'),
    # col=c("#33cccc", "#F3756D"),
    # fill=c('#3288bd','#fdae61'),
    fill=c("#33cccc", "#F3756D"),
    lwd = 0,
    # lty = 'blank',
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    verbose = FALSE,
)
pdf(paste(my_outdir, "Pep_num.Human.Venn.pdf", sep = "/"), height = 71/25.4, width=71/25.4)
grid.draw(venn.plot)
dev.off()
df_Total_Pep %>%
    dplyr::group_by(Origin_Type) %>% dplyr::summarise(n=n()) %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(x="Total",
                  Origin_Type = factor(Origin_Type, 
                             levels=c("Human","Common","Micro"),
                             labels=c("Human Only", "Common","Microbe Only"))) %>% 
    ggplot(aes(x=x, y=n, fill=Origin_Type)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = c('#d53e4f','#66c2a5','#3288bd')) +
    coord_polar(theta = "y") +
    guides(fill=guide_legend(title="Origin Type")) +
    theme_void() +
    theme(legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
    ) -> p_PepType_count
pdf(paste(my_outdir, "Pep_Type_count.Pie.pdf", sep = "/"), height = 53/25.4, width=68/25.4)
print(p_PepType_count)
dev.off()
# Number of proteins for different classified peptides
df_Total_Pep %>%
    tidyr::separate_longer_delim(cols=c(Proteins, Species), delim=';') %>%
    dplyr::mutate(Proteins = gsub("\\-\\d+$", "", Proteins)) %>% distinct() %>%
    group_by(Origin_Type) %>%
    summarise(Protein.number = length(unique(Proteins)))
# Select and save potential-Microbe origing peptides
df_Total_Pep %>% dplyr::filter(Origin_Type != "Human") %>% distinct() %>%
    write.table(file=paste(my_outdir, "../Micro_Origin_Pep.tsv", sep = "/"),
                quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

df_Total_Pep %>% filter(
    grepl("Fusobacterium nucleatum", Species),
    Origin_Type == "Micro",
    # Species.number == 1,
    ) -> df_Fusobacterium_nucleatum

nrow(df_Fusobacterium_nucleatum)
df_Fusobacterium_nucleatum %>% group_by(Origin_Type, Species.number) %>% summarise(n=n())
df_Peptides_filtered %>% dplyr::left_join(dplyr::select(df_Total_Pep, Sequence, Origin_Type), by="Sequence") -> df_Peptides_filtered
df_exogenous <- df_Peptides_filtered %>% dplyr::filter(Origin_Type != "Human")
df_exogenous %>%
    dplyr::select(Sequence, Sample, Sample_Type) %>%
    dplyr::distinct() %>%
    group_by(Sample, Sample_Type) %>%
    summarise(Count=n()) %>%
    ungroup() %>%
    ggplot(aes(x=Sample, y=Count, fill=Sample_Type)) +
    geom_bar(stat="identity", position = "dodge") +
    geom_text(aes(label=Count)) + 
    scale_fill_manual(values = c("#33cccc", "#F3756D")) +
    labs(x="", y="Number of exogenous peptides", fill="") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.line = element_line(color="black",linewidth = 0.6),
        axis.title = element_text(family = 'sans', face= "plain", color='black', size = 8),
        legend.position = c(0.8, 0.8),
        legend.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
    )-> p_exogenous

pdf(paste(my_outdir, "exogenous_Number.pdf", sep = "/"), height = 80/25.4, width=130/25.4)
print(p_exogenous)
dev.off()

df_exogenous %>%
    dplyr::select(Sequence, Sample, Sample_Type) %>%
    dplyr::distinct() %>%
    dplyr::mutate(x = paste(Sample, Sample_Type, sep="_")) %>%
    group_by(x) %>%
    summarise(Count=n()) %>%
    ungroup() %>%
    ggplot(aes(x=x, y=Count)) +
    geom_bar(stat="identity", fill = "#3288bd", width = .8) +
    labs(x="", y="") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.title = element_blank(),
    )-> p_exogenous_2

pdf(paste(my_outdir, "exogenous_Number.bar.v2.pdf", sep = "/"), height = 20/25.4, width=48/25.4)
print(p_exogenous_2)
dev.off()

# Exogenous peptide Venn
library(RColorBrewer)
my_group <- unique(df_exogenous %>%
                       dplyr::select(Sample_Type) %>% 
                       dplyr::pull()
                   )
venn_list = list()
venn_list[[my_group[1]]] <- unlist(distinct(df_exogenous[df_exogenous$Sample_Type == my_group[1],'Sequence']))
venn_list[[my_group[2]]] <- unlist(distinct(df_exogenous[df_exogenous$Sample_Type == my_group[2],'Sequence']))
venn.plot <- venn.diagram(
    venn_list[my_group],
    filename=NULL,
    # Output features
    # imagetype="svg",
    # height = 10,
    # width = 10,
    # resolution = 300,
    # compression = "lzw",
    # Circles
    # col=c('#cc3333','#0000ff'),
    # col=c("#33cccc", "#F3756D"),
    # fill=c('#3288bd','#fdae61'),
    fill=c("#F3756D", "#33cccc"),
    lwd = 0,
    # lty = 'blank',
    # Numbers
    cex = 1,
    fontface = "bold",
    fontfamily = "sans",
    verbose = FALSE,
)
pdf(paste(my_outdir, "Exogenous_Pep_num.Venn.pdf", sep = "/"), height = 80/25.4, width=75/25.4)
grid.draw(venn.plot)
dev.off()


df_exogenous %>%
    dplyr::select(Sequence, Sample, Sample_Type, Proteins, Species, Origin_Type) %>%
    tidyr::separate_longer_delim(cols=c(Proteins, Species), delim = ';') %>% 
    dplyr::filter(!grepl("^CON_", Proteins)) %>%
    dplyr::filter(Species != "Homo sapiens") %>% 
    dplyr::distinct() %>%
    dplyr::mutate(Species = gsub("[^A-Za-z0-9]+", "_", Species)) %>% 
    dplyr::left_join(dplyr::select(df_Micro, Species, Label) %>% distinct(), by="Species") -> df_Pep_Micro
my.Match <- function(x){
    return(unique(df_Micro[grep(substr(x, 1, 40),df_Micro$Species), 'Label']))
}
df_Pep_Micro[is.na(df_Pep_Micro$Label), "Label"] <- sapply(unlist(df_Pep_Micro[is.na(df_Pep_Micro$Label), "Species"]), my.Match)
# df_Pep_Micro <- df_Pep_Micro %>% dplyr::select(-Species) %>% distinct()
df_Pep_Micro %>% 
    dplyr::filter(Origin_Type != "Human") %>%
    dplyr::select(-Sample, -Sample_Type, -Species) %>%
    dplyr::rename(Species = Label) %>%
    dplyr::distinct() %>%
    group_by(Sequence) %>% 
    summarise(
        Protein.number = length(unique(Proteins)),
        Speceis.number = length(unique(Species)),
        Proteins       = paste(Proteins, collapse = ';'),
        Species        = paste(Species, collapse = ';'),
        Origin_Type    = paste(unique(Origin_Type), collapse = ';'),
    ) %>%
    write.table(file=paste(my_outdir, "../Micro_Origin_Pep.tsv", sep = "/"),
                quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

df_Pep_Micro_ID2Name <- read.delim("../../results/Micro_Origin_Pep.Proteins", header = FALSE, sep="\t")
colnames(df_Pep_Micro_ID2Name) <- c("Proteins", "Protein_Names")
df_Pep_Micro %>% 
    # dplyr::distinct(Sequence, Proteins, Label) %>%
    dplyr::left_join(df_Pep_Micro_ID2Name) -> df_Pep_Micro
df_Pep_Micro %>% 
    dplyr::select(Sequence, Proteins, Label, Protein_Names, Origin_Type) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Sequence) %>%
    write.table(file=paste(my_outdir, "../Micro_Origin_Pep.v2.tsv", sep = "/"),
                             quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
df_Pep_Micro %>%
    dplyr::distinct(Sample,Sample_Type,Label,Sequence, .keep_all = TRUE) %>%
    dplyr::group_by(Sample,Sample_Type,Label) %>% 
    dplyr::summarise(
        Pep_num = length(unique(Sequence))
        ) %>%
    ungroup() %>%
    pivot_wider(names_from = Label, values_from = Pep_num) %>%
    pivot_longer(cols = where(is.numeric), names_to = "Label", values_to = "Pep_num") -> df_Sample_PepNum_BySpecies
level <- c("NA",sort(as.numeric(unique(df_Sample_PepNum_BySpecies$Pep_num[!is.na(df_Sample_PepNum_BySpecies$Pep_num)])),decreasing=T))
df_Sample_PepNum_BySpecies$Pep_num <- factor(df_Sample_PepNum_BySpecies$Pep_num, levels = level)
df_Sample_PepNum_BySpecies %>% mutate(x = paste(Sample, Sample_Type, sep="_")) -> df_Sample_PepNum_BySpecies

df_Pep_Micro %>%
    dplyr::distinct(Label, Sequence) %>%
    dplyr::group_by(Label) %>%
    dplyr::summarise(n=length(Sequence)) %>%
    ggplot(aes(x=Label, y= n)) +
    geom_bar(stat="identity", fill="#3288bd", width = .8) +
    # geom_text(aes(y = n+1, label = n)) +
    labs(x="", y="Numbver of Peptides") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_text(family = 'sans', face= "plain", color='black', size = 8),
        axis.text.y = element_blank(),
        axis.line.x = element_line(linewidth = 0.5, color = "black"),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(family = 'sans', face= "plain", color='black', size = 8),
    ) +
    coord_flip() -> p_PepNum_BySpecies
pdf(paste(my_outdir, "PepNum_BySpecies.bar.pdf", sep = "/"), height = 80/25.4, width = 18/25.4)
print(p_PepNum_BySpecies)
dev.off()
ggplot(df_Sample_PepNum_BySpecies, aes(x=x,y=Label)) +
    geom_tile(aes(fill=Pep_num), colour="grey70") +
    # scale_fill_manual(values =c("#543005","#8c510a","#bf812d","#dfc27d",
    #                             "#f6e8c3","#f5f5f5","#c7eae5","#80cdc1",
    #                             "#35978f","#66c2a4","#01665e","#003c30")) +
    # scale_fill_manual(values = c('#a50026','#d73027','#f46d43','#fdae61',
    #                          '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
    #                          "#66c2a4",'#66bd63','#1a9850','#006837')) +
    scale_fill_manual(values = c('#a50026','#d73027','#f46d43','#fdae61',
                                 '#fee08b','#ffffbf','#d9ef8b','#a6d96a',
                                 '#66bd63','#1a9850')) +
    labs(x="",y="",fill="") +
    scale_x_discrete(position="top") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
          axis.text.x = element_text(angle = 45, vjust=0, hjust=0),
          axis.ticks.x.top = element_line(color="black",linewidth = 1),
          legend.position = "right") -> p_Sample_PepNum_BySpecies
pdf(paste(my_outdir, "Sample_PepNum_BySpecies.count.pdf", sep = "/"), height = 95/25.4, width = 125/25.4)
print(p_Sample_PepNum_BySpecies)
dev.off()

length(unique(df_Sample_PepNum_BySpecies$Label))

read.delim("../../results/DB_Search/ALL_Exogenous.DB_search.stat", header=T, sep="\t") %>%
    select(-Sequnece_I2L) -> df_DBsearch
df_DBsearch['Reported'] <- apply(df_DBsearch[,-1], MARGIN = 1, function(x){if(sum(x)>0){return(1)}else{return(0)}})
sum(df_DBsearch$Reported)
df_Pep_Micro %>% 
    group_by(Sequence) %>% 
    summarise(Species.Num = length(unique(Label)),
              Sample.Num = length(unique(Sample)),
              Species = paste(sort(unique(Label)), collapse = ';'),
              Origin_Type = unique(Origin_Type),
              Sample_Type = paste(sort(unique(Sample_Type)), collapse = ';'),
              Protein_Names = paste(sort(unique(Protein_Names)), collapse = ';'),
              Sample = paste(sort(unique(Sample)), collapse = ';')) %>% 
    left_join(df_DBsearch) %>%
    mutate(Sample_Type = factor(Sample_Type,
                                levels=c("normal","tumor","normal;tumor"),
                                labels = c("Normal", "Tumor","BOTH")),
           Origin_Type = as.factor(Origin_Type),
           Reported = factor(Reported,
                             levels = c(0,1),
                             labels = c("Not Reported", "Reported"))
           )-> df_MicroPep_summ

df_MicroPep_summ %>%
    dplyr::arrange(Origin_Type, Sample_Type,Sequence) %>% dplyr::pull(Sequence) -> sequence_arranged

df_MicroPep_summ %>% group_by(Origin_Type, Sample_Type, Sample.Num) %>%
    summarise(n=n())
df_MicroPep_summ %>% group_by(Origin_Type, Reported) %>% summarise(n=n())


df_MicroPep_summ %>%
    dplyr::filter(Sequence %in% df_Fusobacterium_nucleatum$Sequence) %>% 
    count(Sample_Type, Species.Num, Sample.Num)

df_Peptides %>% filter(Sequence %in% df_Fusobacterium_nucleatum$Sequence) %>% view()


df_MicroPep_summ %>%
    dplyr::mutate(Sample.Num = as.factor(Sample.Num),
                  Species.Num = as.factor(Species.Num)) %>%
    dplyr::select(Sequence, Origin_Type, Sample.Num, Sample_Type, Species.Num, Reported) %>%
    pivot_longer(cols = c(Origin_Type, Sample.Num, Sample_Type, Species.Num, Reported), names_to = "type", values_to = "class") %>%
    dplyr::mutate(class = factor(class,
                                 levels = c("Common", "Micro",
                                            "Normal", "Tumor", "BOTH",
                                            "1", "2", "3", "4", "5", "6", "7", "9", "10",
                                            "Not Reported", "Reported")),
                  type = factor(type,
                                levels = rev(c("Origin_Type", "Sample_Type", "Reported", "Sample.Num", "Species.Num")),
                                labels = rev(c("Origin", "Type of Sample", "Reported", "Number of Sample", "Number of Species"))),
                  Sequence = factor(Sequence, levels=sequence_arranged)) %>%
    dplyr::arrange(Sequence, class, type) %>%
    pivot_wider(names_from = type, values_from = class) %>%
    write.table(file=paste(my_outdir, "../Micro_Origin_pep.classify.tsv", sep='/'),
                quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
df_MicroPep_summ %>%
    dplyr::mutate(Sample.Num = as.factor(Sample.Num),
                  Species.Num = as.factor(Species.Num)) %>%
    dplyr::select(Sequence, Origin_Type, Sample.Num, Sample_Type, Species.Num, Reported) %>%
    pivot_longer(cols = c(Origin_Type, Sample.Num, Sample_Type, Species.Num, Reported), names_to = "type", values_to = "class") %>%
    dplyr::mutate(class = factor(class,
                                 levels = c("Common", "Micro",
                                            "Normal", "Tumor", "BOTH",
                                            "1", "2", "3", "4", "5", "6", "7", "9", "10",
                                            "Not Reported", "Reported")),
                  type = factor(type,
                                levels = rev(c("Origin_Type", "Sample_Type", "Reported", "Sample.Num", "Species.Num")),
                                labels = rev(c("Origin", "Type of Sample", "Reported", "Number of Sample", "Number of Species"))),
                  Sequence = factor(Sequence, levels=sequence_arranged)) %>%
    ggplot(aes(x=Sequence, y=type, fill=class)) +
    geom_tile() +
    scale_fill_manual(values = c("#d53e4f", "#3288bd",
                                 "#f46d43", "#66c2a5", "#ffffbf",
                                 "#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#15315d",
                                 "#ffffff", "#9e0142")) +
    labs(x="",y="", fill="") -> p_MicroPep_Sum
pdf(paste(my_outdir, "Micro_Pep_Summ.pdf", sep = "/"), height = 50/25.4, width = 215/25.4)
p_MicroPep_Sum + 
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_text(family = 'sans', face= "plain", color='black', size = 6, angle = 45, hjust = 0, vjust = 1),
    )
dev.off()

colnames(df_MicroPep_summ)

df_MicroPep_summ %>%
    dplyr::filter(Reported == "Reported") %>%
    dplyr::arrange(Origin_Type, Sample_Type,Sequence) %>% dplyr::pull(Sequence) -> selected_sequence_arranged


df_MicroPep_summ %>% 
    dplyr::filter(Reported == "Reported") %>% 
    dplyr::select(Sequence, IEDB, IEDB_T_Cell, HLA_Ligand_Atlas, HLAtlas_Cryptic, IEAtlas, dbPepNeo, MHC_Motif_Atlas) %>%
    tidyr::pivot_longer(cols = where(is.numeric), names_to = "DB_ID", values_to = "Reported") %>% 
    dplyr::mutate(Reported = factor(Reported, levels = c(0,1), labels = c("Not Reported", "Reported"))) %>%
    dplyr::bind_rows(dplyr::filter(df_MicroPep_summ, Reported == "Reported") %>%
                         dplyr::select(Sequence, Origin_Type, Sample_Type) %>%
                         tidyr::pivot_longer(cols = c(Origin_Type, Sample_Type), names_to = "DB_ID", values_to = "Reported")) %>%
    dplyr::mutate(DB_ID = factor(DB_ID, levels=c("dbPepNeo", "HLA_Ligand_Atlas", "HLAtlas_Cryptic", "IEAtlas",
                                                 "IEDB", "IEDB_T_Cell", "MHC_Motif_Atlas", "Sample_Type", "Origin_Type")),
                  Sequence = factor(Sequence, levels=selected_sequence_arranged),
                  ) %>%
    ggplot(aes(x=Sequence, y=DB_ID, fill=Reported)) +
    geom_tile(color="#CCCCCC", linewidth=0.1) +
    scale_fill_manual(values = c("#ffffff", "#9e0142",
                                 "#d53e4f", "#3288bd",
                                 "#f46d43", "#66c2a5", "#ffffbf")) +
    labs(x="",y="", fill="") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(family = 'sans', face= "plain", color='black', size = 6, angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none"
    ) -> p_DB_report
pdf(paste(my_outdir, "Micro_Pep_DBsearch.pdf", sep = "/"), height = 27/25.4, width = 90/25.4)
print(p_DB_report)
dev.off()

df_MicroPep_summ %>% 
    dplyr::filter(Reported == "Reported") %>%
    dplyr::filter(Origin_Type == "Micro") %>% dplyr::pull(Sequence) %>% length()

df_MicroPep_summ %>% 
    dplyr::filter(Sample.Num>1) %>% 
    dplyr::select(Sequence, Species.Num, Sample.Num, Species, Origin_Type) %>%
    dplyr::left_join(select(df_Pep_Micro, Sequence, Sample, Sample_Type)) %>%
    dplyr::distinct() %>%
    dplyr::mutate(x=paste(Sample,Sample_Type, sep="_"),
                  y=paste(Species, Sequence, sep="\\n")) -> df_MultiSample_MicroPep

df_MicroPep_summ %>% 
    filter(
        Sample.Num > 1,
        ) %>%
    distinct(Sequence, Origin_Type, .keep_all = TRUE) %>% group_by(Origin_Type) %>% summarise(n=n())

df_MicroPep_summ %>% 
    filter(
        Sample.Num > 1,
        Origin_Type == "Micro",
    ) %>%
    distinct(Sequence, Origin_Type, .keep_all = TRUE) %>%
    filter(grepl("Elongation factor", Protein_Names)) %>%
    # filter(grepl("Chaperone protein", Protein_Names)) %>%
    view()
df_MicroPep_summ %>% 
    # dplyr::filter(Sample.Num > 1) %>% group_by(Origin_Type) %>% summarise(n=n())
    dplyr::filter(
        Origin_Type == "Micro",
        Sample_Type == "Tumor",
        Sample.Num > 1,
        ) %>%
    distinct(Sequence, Species) %>%
    tidyr::separate_longer_delim(cols = Species, delim = ';') %>%
    distinct(Sequence, Species) %>% 
    group_by(Species) %>% summarise(n=n())
    
df_MultiSample_MicroPep %>%
    dplyr::filter(
        Origin_Type == "Micro",
        Sample_Type == "tumor",
                  ) %>%
    dplyr::select(x, Sequence, Species.Num) %>%
    tidyr::pivot_wider(names_from = Sequence, values_from = Species.Num) %>%
    tidyr::pivot_longer(cols=where(is.numeric),names_to = "Sequence", values_to = "Species.Num") %>%
    dplyr::mutate(Species.Num = factor(Species.Num,
                                       levels = c(1:6,8:10))) %>% 
    ggplot(aes(x=x,y=Sequence,fill=Species.Num)) +
    geom_tile(color="#CCCCCC", linewidth=0.2) +
    scale_fill_manual(values = c("#f7fbff", "#deebf7", "#c6dbef",
                                 "#9ecae1", "#6baed6", "#4292c6",
                                 "#2171b5", "#08519c", "#08306b")) +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
    ) -> P_MultiSample_MicroPep
# pdf(paste(my_outdir, "MultiSample_MicroPep.pdf", sep = "/"), height = 85/25.4, width = 85/25.4)
pdf(paste(my_outdir, "MultiSample_MicroPep.v2.pdf", sep = "/"), height = 85/25.4, width = 85/25.4)
print(P_MultiSample_MicroPep)
dev.off()

length(unique(df_MultiSample_MicroPep$Sequence))

df_affinity <- read.delim("../../results/All.affinity_features.txt", header = TRUE, sep="\t")
colnames(df_affinity) <- c('Peptide', 'Length', 'BinderNum', 'MHC',
                           'IC50nM', 'Rank_EL', 'Rank_BA', 'BinderType',
                           'Sample', 'Sample_Type')
df_affinity %>%
    dplyr::select(-Length, -BinderNum) %>%
    dplyr::filter(Peptide %in% df_MicroPep_summ$Sequence) %>%
    tidyr::separate_longer_delim(cols = c(MHC, IC50nM, Rank_EL, Rank_BA, BinderType), delim = ',') %>%
    # dplyr::filter(Peptide %in% c("AEIKTGAQAI", "AIPHSPILSL", "ATTAALHLR",
    #                              "FLDEPTNHL", "FLRELISNA", "GVTVGILVGV",
    #                              "KEIFLRELI", "TGPAGASASI", "VMKKLLETV",
    #                              "FLTLLREQTERL", "ITGNSILTV")) %>%
    dplyr::group_by(Sample, Peptide) %>%
    dplyr::arrange(Rank_BA) %>%
    dplyr::distinct(Sample, Peptide, .keep_all = TRUE) %>%
    dplyr::mutate(BinderType = factor(BinderType,
                                      levels = c("SB","WB","NB"))) %>%
    # dplyr::group_by(Sample) %>%
    # summarise(
    #     n = length(unique(Peptide))
    # )
    ggplot(aes(x=Sample)) +
    geom_bar(aes(fill=BinderType))

df_MicroPep_summ %>% 
    dplyr::filter(Sequence %in% c("VMKKLLETV", "ITGNSILTV", "TGPAGASASI", "GVTVGILVGV",
                                  "GKSTLLNTL", "FLTLLREQTERL", "ATTAALHLR", "AIPHSPILSL",
                                  "AGKALEELQL", "AEIKTGAQAI", "FLRELISNA", "KEIFLRELI",
                                  "SLISGFTTA", "FLDEPTNHL")) %>%
    # dplyr::distinct(Sequence) %>% pull()
    view()

df_MultiSample_MicroPep %>% 
    dplyr::filter(Sequence %in% c("AEIKTGAQAI", "AIPHSPILSL", "ATTAALHLR", 
                                  "FLDEPTNHL", "FLRELISNA", "GVTVGILVGV", 
                                  "KEIFLRELI", "TGPAGASASI", "VMKKLLETV", 
                                  "FLTLLREQTERL", "ITGNSILTV")) %>%
    dplyr::distinct(Sequence) %>% nrow()

df_MicroPep_summ %>% write.table(file=paste(my_outdir, "../Micro_Origin_Pep.tsv", sep = "/"),
                                 quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
