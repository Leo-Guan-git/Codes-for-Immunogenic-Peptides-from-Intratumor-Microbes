library(tidyverse)
library(ggplot2)
library(scales)
library(ggsignif)
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"

read.csv("../../results/All.Peptides.tsv.gz", header=T, sep="\t") %>% 
    mutate(Sample = gsub("CRC", "#", Sample)) -> df_Peptides

df_Peptides %>% select(Sample,Sequence) %>%
    distinct() %>% group_by(Sample) %>% summarise(Total=n())-> df_LigandNum_bysample

df_Peptides %>% select(Sample, Sequence, Species.number, Leading.razor.species) %>%
    filter(Species.number == 1, Leading.razor.species == "Homo sapiens") %>%
    distinct(Sample, Sequence) %>% group_by(Sample) %>% summarise(HumanOrigin=n()) %>% 
    right_join(df_LigandNum_bysample) %>%
    mutate(Other = Total - HumanOrigin) %>%
    select(-Total) %>%
    pivot_longer(cols=where(is.numeric)) %>% 
    group_by(Sample) %>% 
    arrange(Sample, desc(name)) %>%
    mutate(Per = signif(value/sum(value)*100,3),
           POS = cumsum(value)-0.5*value,
           ) -> df_LigandNum_byorigin_pivort

ggplot(df_LigandNum_byorigin_pivort, aes(x=Sample, y=value, fill=name)) +
    geom_col(position = "stack") +
    geom_text(data=filter(df_LigandNum_byorigin_pivort, name=='HumanOrigin'), aes(y=POS, label = paste0(Per, '%'))) +
    scale_fill_manual(values = c("#3288bd", "#d53e4f")) +
    labs(x="", y="Number of Peptides", fill="") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.line = element_line(colour = "black", linewidth = .5),
        legend.background = element_blank(),
        legend.position = c(0.8,0.1),
        legend.key.size = unit(0.3, "cm")
    ) +
    coord_flip()-> p_peptides_bySample

pdf(paste(my_outdir, "Ligand_Origin_type.pdf", sep = "/"), height = 70/25.4, width=45/25.4)
print(p_peptides_bySample)
dev.off()

df_Peptides %>% select(Sample, Sample_Type, Sequence) %>%
    distinct() %>% group_by(Sample, Sample_Type) %>% summarise(Num=n())-> df_LigandNum_bySampleType

ggplot(df_LigandNum_bySampleType, aes(x=Sample, y=Num, fill=Sample_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#F3756D", "#33cccc")) +
    labs(x="", y="Number of Peptides", fill="") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
        axis.line = element_line(colour = "black", linewidth = .5),
        axis.title = element_text(family = 'sans', face= "plain", color='black', size = 10),
        legend.position = c(0.8,0.8)
    ) -> p_peptides_bySampleType
pdf(paste(my_outdir, "Peptide_Number_BySample.pdf", sep = "/"), height = 55/25.4, width=75/25.4)
print(p_peptides_bySampleType)
dev.off()
