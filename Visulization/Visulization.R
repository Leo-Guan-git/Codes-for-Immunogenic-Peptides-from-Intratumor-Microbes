getwd()
library(tidyverse)
library(ggplot2)
library(scales)
library(ggsignif)
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"

## Bar plot of number of microbes and proteins in each sample
df_DB <- read.csv('../../results/All.selected_microbe.tsv', header=F, sep="\t")
colnames(df_DB) <- c('Sample','Class','Count')
head(df_DB)

# filtering of Aquabacterium olei
Aquabacterium_olei.Pro_num <- 3665
df_DB[grep("Microbes_Number",df_DB$Class), 'Count'] <- df_DB[grep("Microbes_Number",df_DB$Class), 'Count'] - 1
df_DB[grep("Proteins_Number",df_DB$Class), 'Count'] <- df_DB[grep("Proteins_Number",df_DB$Class), 'Count'] - Aquabacterium_olei.Pro_num
df_DB %>% group_by(Class) %>% reframe(n = quantile(Count))
Arial_regular_9 <- element_text(family = 'sans', face= "plain", color='black', size = 8)
Arial_regular_6 <- element_text(family = 'sans', face= "plain", color='black', size = 6)
df_DB %>% mutate(Sample=paste0('#',substring(Sample,4))) -> df_DB
f <- function(x){
    if_else(x>0,x,-1*x/1000)
}
df_DB %>% group_by(Class) %>%
    mutate(max_count=max(Count),
           min_count=min(Count)) %>%
    select(-c(Sample,Count)) %>% distinct()
p <- ggplot(df_DB %>% mutate(Count=if_else(Class == "Proteins_Number", Count, -1*Count*1000),
                        Class=gsub("_"," ", Class)),
       aes(x=Sample, y=Count, fill=Class)) +
    geom_bar(stat = "identity", position = "identity", width = 0.6) +
    geom_text(aes(y=Count*1.05,label=if_else(Count > 0, Count, -1*Count/1000)),
              size=6*5/14, family = 'sans', fontface= "plain", color='black') +
    scale_fill_manual(values = c("#d53e4f","#3288bd"), 
                      labels=c("Number of species", "Number of Proteins")) +
    scale_y_continuous(labels = function(x){if_else(x>0,x,-1*x/1000)}) +
    labs(x="", y="Numbers", fill="") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_blank(),
          axis.line.y = element_line(color="black",linewidth = 0.5),
          axis.ticks.x = element_blank(),
          legend.position = "top",
          axis.text = Arial_regular_6,
          legend.title = Arial_regular_6,
          plot.margin = margin(0,0,0,0),
          )
pdf(paste(my_outdir, "Micro_Pro_Number.plot.pdf", sep = "/"), height = 56/25.4, width=75/25.4)
print(p)
dev.off()

df_DB %>% mutate(Class=gsub("_"," ",Class)) %>%
    pivot_wider(id_cols = Sample, names_from = 'Class', values_from = 'Count') %>%
    write.csv(paste(my_outdir, "Micro_Pro_Number.csv", sep = "/"), row.names = FALSE)

## Species tile plot
df_Species <- read.csv('../../results/visulization/Species_Selected.csv')
head(df_Species)
df_Species <- df_Species[!grepl("Aquabacterium",df_Species$Species),]
df_Species %>% count(SampleID,Domain)

Domains <- unique(df_Species$Domain)

df_Species %>% select(SampleID, Species) %>% distinct() %>%
    mutate(count=1) %>%
    pivot_wider(id_cols=SampleID,
                names_from = Species, 
                values_from = count) %>%
    pivot_longer(cols = where(is.numeric),
                 names_to = 'Species',
                 values_to = 'Included') -> df_species_tile
p_list <- list()
for (i in 1:length(Domains)){
    domain_tmp <- Domains[i]
    species_tmp <- df_Species %>% filter(Domain == domain_tmp) %>%
        select(Species) %>% distinct()
    ggplot(filter(df_species_tile, Species %in% species_tmp$Species)) +
        geom_tile(aes(x=SampleID, y=Species, fill=as.factor(Included)), colour="grey80") +
        scale_fill_manual(values = c("#F3756D")) +
        theme(panel.background = element_blank(),
              panel.grid = element_blank(),
              axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
              axis.ticks.x.top = element_line(color="black",linewidth = 1),
              legend.position = "None") -> p_list[[Domains[i]]]
}

library(cowplot)
df_type <- df_Species %>% select(Domain, Species) %>% 
    distinct() %>% count(Domain) %>% 
    arrange(n) %>% mutate(Per = n/sum(n))
p_draw = list()
for (n in 1:nrow(df_type)){
    domain_tmp <- df_type[["Domain"]][n]
    sample_num <- df_type[["n"]][n]
    p_tmp <- p_list[[domain_tmp]]
    if (n < nrow(df_type)){
        p_tmp <- p_tmp + theme(axis.line.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.x = element_blank(),
                               )
    }
    p_tmp <- p_tmp + theme(axis.title.y = element_blank(),
                           plot.margin = margin(0,0,0,0))
    p_draw[[n]] <- p_tmp
}

pdf(paste(my_outdir, "Micro_Species.tile.pdf", sep = "/"), height = 80/25.4, width=80/25.4)
plot_grid(p_draw[[1]], p_draw[[2]], align = 'v', 
          axis = 'l', nrow = 2, rel_heights = c(0.05,0.9), byrow=F)
dev.off()

## MaxQuant MS/MS PSM Visulization
read.csv("../../results/All.Peptides.tsv.gz", header=T, sep="\t") %>% 
    mutate(Sample = gsub("CRC", "#", Sample)) -> df_Peptides
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

df_Peptides_filtered %>%
    separate_rows(MS.MS.IDs, sep=';') %>%
    select(Sample,Sample_Type, MS.MS.IDs) %>%
    count(Sample,Sample_Type) %>%
    mutate(Sample_Type=factor(Sample_Type, levels=rev(levels(factor(Sample_Type))))) -> df_PSM_num

p_PSM_Num <- ggplot(df_PSM_num, aes(x=Sample, y=n)) +
    geom_bar(aes(fill=Sample_Type), stat = "identity", position = position_dodge(width = 0.9)) +
    scale_fill_manual(values = c("#F3756D", "#33cccc")) +
    geom_text(aes(label=n), position = position_dodge(width = 0.9), vjust=-.5) +
    labs(x="", y="Number of PSMs", fill="") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.ticks.x = element_blank(),
          axis.line = element_line(color="black",linewidth = 0.6),
          axis.title = element_text(family = 'sans', face= "plain", color='black', size = 8),
          legend.position = c(0.8, 0.8),
          legend.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
          )
pdf(paste(my_outdir, "PSM_Number.pdf", sep = "/"), height = 56/25.4, width=72/25.4)
print(p_PSM_Num)
dev.off()

df_Peptides_filtered %>%
    separate_rows(MS.MS.IDs, sep=';') %>%
    select(Sequence, MS.MS.IDs) %>%
    group_by(Sequence) %>%
    summarise(n=n()) %>%
    ungroup() %>% mutate(n=if_else(n<=10,as.character(n),">10"),) %>%
    group_by(n) %>%
    summarise(SeqNum=n()) %>%
    mutate(n=factor(n, levels = c(as.character(1:10), ">10"))) %>%
    arrange(rev(levels(n))) %>%
    mutate(Per=signif(SeqNum/sum(SeqNum),3),
           Pos=cumsum(SeqNum) -0.5*SeqNum)-> df_PSMnumPerSeq

p_PSM_count <- ggplot(df_PSMnumPerSeq,
       aes(x=2,y=SeqNum,fill=n),
) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(x=2.5, y=Pos, label=paste0(Per*100, '%')), hjust=1, vjust=.5) +
    scale_fill_manual(values=c('#9e0142',
                               '#d53e4f',
                               '#f46d43',
                               '#fdae61',
                               '#fee08b',
                               '#ffffbf',
                               '#e6f598',
                               '#abdda4',
                               '#66c2a5',
                               '#3288bd',
                               '#5e4fa2')) +
    coord_polar(theta = "y") +
    guides(fill=guide_legend(title="PSM Counts")) +
    theme_void() +
    theme(legend.key.size = unit(0.2, "cm"),
          legend.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
    )
pdf(paste(my_outdir, "PSM_Count.Pie.pdf", sep = "/"), height = 56/25.4, width=56/25.4)
print(p_PSM_count)
dev.off()
#  MaxQuant Charges Visulization
df_Peptides_filtered %>% separate_rows(Charges, sep=';') %>%
    select(Sample,Sample_Type, Charges) %>%
    count(Charges) %>%
    mutate(Charges=factor(Charges, levels=unique(Charges), labels=paste0(unique(Charges), "+")),
           Desc=paste0(signif(n/sum(n)*100, 4), '%'),
           ) -> df_Charges
p_Charge <- ggplot(df_Charges, aes(x=Charges,y=n)) +
    # geom_bar(stat = "identity", width = 0.65, fill='grey') +
    geom_bar(stat = "identity", width = 0.65, fill='#3288bd') +
    geom_text(aes(label=Desc), vjust=-.5) +
    labs(x="Precusor Charge", y="Number of Spectrums") +
    coord_flip() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.line = element_line(color="black",linewidth = 0.5),
          axis.title = element_text(family = 'sans', face= "plain", color='black', size = 8),
          legend.position = c(0.8, 0.8),
          legend.text = element_text(family = 'sans', face= "plain", color='black', size = 8),
    )
pdf(paste(my_outdir, "Charges.bar.pdf", sep = "/"), height = 60/25.4, width=40/25.4)
print(p_Charge)
dev.off()
#  MaxQuant Sequence Length Visulization
df_Peptides_filtered %>% select(Sequence,Length) %>% distinct() -> df_SeqLen
nrow(df_SeqLen)
p_PepLen <- ggplot(df_SeqLen, aes(x=Length)) +
    geom_bar(fill="#3288bd") +
    labs(x="Length of Peptides (aa)", y="Number of Peptides") +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
          axis.ticks.x = element_line(color="black", linewidth = 0.5),
          axis.line = element_line(color="black",linewidth = 0.5),
          axis.title = element_text(family = 'sans', face= "plain", color='black', size = 8),
          )
pdf(paste(my_outdir, "Length.bar.pdf", sep = "/"), height = 40/25.4, width=55/25.4)
print(p_PepLen)
dev.off()

df_affinity <- read.delim("G:/Micro_Maxquant_v2/results/All.mostBindingHLA.affinity.txt", header = T, sep="\t")
head(df_affinity)
colnames(df_affinity)
df_affinity %>% 
    filter(Peptide %in% df_Peptides_filtered$Sequence) %>%
    group_by(BinderType) %>% summarise(num=n()) %>% 
    mutate(BinderType=factor(BinderType, levels = c('SB', 'WB', 'NB'))) %>% 
    arrange(desc(BinderType)) %>%
    mutate(Per = signif(num/sum(num)*100, 3),
           Pos = cumsum(num)-0.5*num,
           ) %>%
    ggplot(aes(x="", fill=BinderType, y=num)) +
    geom_bar(width = 1, stat = "identity") +
    geom_text(aes(y=Pos, label=paste0(Per, "%"))) +
    labs(x="", y="Number of Peptides", fill="Bind Type") +
    scale_fill_manual(values = c('#d53e4f','#66c2a5','#3288bd')) +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line(colour = "black", linewidth = .5),
        legend.position = "top",
        legend.key.size = unit(0.3, "cm")
    ) -> p_affinity
pdf(paste(my_outdir, "affinity.bar.pdf", sep = "/"), height = 55/25.4, width=25/25.4)
print(p_affinity)
dev.off()

df_weight <- read.delim("G:/Micro_Maxquant_v2/results/Sample_weight.csv", sep=",")
df_Peptides_filtered %>% select(Sample, Sample_Type, Sequence) %>% distinct() %>%
    group_by(Sample, Sample_Type) %>% summarise(n=n()) %>%
    left_join(df_weight %>% pivot_longer(cols=where(is.numeric)), 
              by=join_by(Sample == ID, Sample_Type == name)) -> df_ligand
df_ligand %>% mutate(n_per_mg = n/value) -> df_ligand
ggplot(df_ligand, aes(x=Sample_Type, y=n)) +
    geom_boxplot(
        width = .5,
        outlier.fill = NA,
        outlier.shape = NA,
        ) +
    geom_point(
        size = 1,
        position = position_jitter(width = 0.2),
        ) +
    geom_signif(
        comparisons = list(c("normal", "tumor")),
        y_position = 2500,
        map_signif_level = TRUE,
        textsize = 6,
        vjust = -0.5,
        ) +
    labs(x="Sample Type", y="Total number of ligands") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .5),
    ) -> p_total_Box

ggplot(df_ligand, aes(x=Sample_Type, y=n_per_mg)) +
    geom_boxplot(
        width = .5,
        outlier.fill = NA,
        outlier.shape = NA,
    ) +
    geom_point(
        size = 1,
        position = position_jitter(width = 0.2),
    ) +
    geom_signif(
        comparisons = list(c("normal", "tumor")),
        y_position = 12,
        map_signif_level = TRUE,
        textsize = 6,
        vjust = -0.5,
    ) +
    labs(x="Sample Type", y="number of ligands (/mg)") +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black", linewidth = .5),
    ) -> p_Permg_Box
pdf(paste(my_outdir, "Ligand_Num_compare.pdf", sep = "/"))
plot_grid(p_total_Box, p_Permg_Box, align = 'h', 
          axis = 'l', nrow = 1, byrow=T)
dev.off()

# fitted curve
m <- lm(n ~ value, df_ligand);
eq <- substitute(italic("y") == a + b %.% italic("x")*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))
ggplot(df_ligand, aes(x=value, y=n)) +
    geom_point() +
    geom_smooth(aes(x=value,y=n), method="lm", formula = 'y~x') +
    geom_text(x= 200, y= 3000, label=as.character(as.expression(eq)),hjust=0, parse = TRUE) +
    theme(
        panel.background = element_blank(),
        axis.text = Arial_regular_9,
        axis.line = element_line(color="black", linewidth = 1)) -> p_pep2weight_lm
pdf(paste(my_outdir, "SampleWeight2PepNum.plot.pdf", sep = "/"))
print(p_pep2weight_lm)
dev.off()

# Peptide Venn
library(VennDiagram)
library(RColorBrewer)

my_group <- unique(df_Peptides_filtered$Sample_Type)
venn_list = list()
venn_list[[my_group[1]]] <- unlist(distinct(df_Peptides_filtered[df_Peptides_filtered$Sample_Type == my_group[1],'Sequence']))
venn_list[[my_group[2]]] <- unlist(distinct(df_Peptides_filtered[df_Peptides_filtered$Sample_Type == my_group[2],'Sequence']))
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
pdf(paste(my_outdir, "Pep_num.Venn.pdf", sep = "/"), height = 71/25.4, width=71/25.4)
grid.draw(venn.plot)
dev.off()
