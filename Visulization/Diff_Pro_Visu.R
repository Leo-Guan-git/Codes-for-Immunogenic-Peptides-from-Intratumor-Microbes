library(tidyverse)
library(ggplot2)
library(scales)
library(ggrepel)
library(ggsignif)


library(ggplot2)#柱状图和点状图
library(stringr)#基因ID转换
library(biomaRt)
library(enrichplot)#GO,KEGG,GSEA
library(clusterProfiler)#GO,KEGG,GSEA
library(GOplot)#弦图，弦表图，系统聚类图
library(DOSE)
library(ggnewscale)
library(topGO)#绘制通路网络图
library(circlize)#绘制富集分析圈图
library(ComplexHeatmap)#绘制图例
my_outdir <- "G:/Micro_Maxquant_v2/results/visulization"

read.csv("../../results/All.Peptides.tsv.gz", header=T, sep="\t") %>% 
    mutate(Sample = gsub("CRC", "#", Sample)) -> df_Peptides
colnames(df_Peptides)
df_Peptides %>%
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
df_Total_Pep %>%
  dplyr::mutate(Origin_Type = if_else((grepl("Homo sapiens", Species) & Species.number == 1), "Human", "Micro")) -> df_Total_Pep
df_Total_Pep[(grepl("Homo sapiens", df_Total_Pep$Species) & (grepl("Micro", df_Total_Pep$Origin_Type))), 'Origin_Type'] <- "Common"
df_Total_Pep %>%
  dplyr::group_by(Origin_Type) %>% dplyr::summarise(n=n())
# 
df_Peptides %>% dplyr::left_join(dplyr::select(df_Total_Pep, Sequence, Origin_Type), by="Sequence") -> df_Peptides
df_Peptides %>%
  tidyr::separate_longer_delim(cols = c(Proteins, Species), delim=';') %>%
  dplyr::filter(!grepl("^CON_", Proteins)) %>%
  # dplyr::mutate(Proteins = gsub("CON_[^;]+;?", "", Proteins)) %>%
  dplyr::select(Sequence, Sample, Sample_Type, Proteins, Species, Origin_Type) %>%
  dplyr::filter(Origin_Type == "Human") %>%
  dplyr::mutate(Proteins = gsub("\\-\\d+$", "", Proteins)) %>%
  distinct() %>%
  group_by(Sample, Sample_Type, Proteins) %>%
  summarise(n = n()) %>%
  arrange(Sample, Sample_Type, desc(n)) %>%
  ungroup() -> df_protein_human
# 
df_protein_human %>% group_by(Proteins) %>%
    summarise(n = sum(n)) %>%
    arrange(desc(n))%>%
    dplyr::select(Proteins) %>%
    write.table(paste(my_outdir, "../CRC.All.Protein.txt", sep="/"), sep="\t", row.names = F, col.names = F, quote  = F)
# 
df_protein_human %>% group_by(Proteins) %>%
    summarise(n = sum(n)) %>%
    filter(n > 10) %>%
    dplyr::select(Proteins) %>%
    write.table(paste(my_outdir, "../CRC.All.Protein.Count_10.txt", sep="/"), sep="\t", row.names = F, col.names = F, quote  = F)

df_Proteins_reshape <- data.frame(Proteins=character(),
                                  tumor_mean=numeric(),
                                  normal_mean=numeric(),
                                  tumor_n=numeric(),
                                  normal_n=numeric(),
                                  Pvalue=numeric())
n_sample <- length(unique(df_protein_human$Sample))
my.t.test <- function(...) {
    obj<-try(t.test(...), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else return(obj)
}
for (i in unique(df_protein_human$Proteins)){
    normal <- df_protein_human %>% filter(Sample_Type == 'normal', Proteins == i) %>% dplyr::pull(n)
    tumor <- df_protein_human %>% filter(Sample_Type == 'tumor', Proteins == i) %>% dplyr::pull(n)
    normal_n <- length(normal)
    tumor_n <- length(tumor)
    if ((normal_n < 2) | (tumor_n < 2)){
        next
    }
    t <- my.t.test(x=tumor, y=normal)
    if (!is(t, "htest")){
        next
    }
    df_Proteins_reshape[nrow(df_Proteins_reshape) + 1, ] <- c(Proteins=i,
                                                             tumor_mean=t$estimate[[1]],
                                                             normal_mean=t$estimate[[2]],
                                                             tumor_n=tumor_n,
                                                             normal_n=normal_n,
                                                             Pvalue=t$p.value)
}
df_Proteins_reshape %>% dplyr::mutate(tumor_mean = as.numeric(tumor_mean),
                                      normal_mean = as.numeric(normal_mean),
                                      tumor_n = as.numeric(tumor_n),
                                      normal_n = as.numeric(normal_n),
                                      Pvalue = as.numeric(Pvalue),
                                      FC=tumor_mean/normal_mean) -> df_Proteins_reshape
#gene ID转换
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# listAttributes(ensembl) %>% view()
result <- getBM(attributes = c("ensembl_gene_id", "uniprot_gn_symbol","uniprotswissprot","entrezgene_id"), filters = "uniprotswissprot",
      values = df_Proteins_reshape$Proteins, mart = ensembl)
result %>% 
    dplyr::distinct(uniprotswissprot, .keep_all = T) %>%
    dplyr::right_join(df_Proteins_reshape,
                      by=join_by(uniprotswissprot==Proteins)) -> df_Proteins_reshape

threshold.p <- 0.05
threshold.log2FC <- 1
DrawValcano <- function(df_t_test, threshold.p, threshold.log2FC){
    rownames(df_t_test) = c(1:nrow(df_t_test))
    # df_t_test %>% mutate(FC = y_mean/x_mean) -> df_t_test
    df_t_test$diff <- "No"
    df_t_test$diff[((log2(df_t_test$FC) > threshold.log2FC) & (df_t_test$Pvalue < threshold.p))] <- "UP"
    df_t_test$diff[((log2(df_t_test$FC) < -threshold.log2FC) & (df_t_test$Pvalue < threshold.p))] <- "DOWN"
    df_t_test$delabel <- NA
    df_t_test$delabel[((abs(log2(df_t_test$FC)) > threshold.log2FC) & (df_t_test$Pvalue < threshold.p))] <- df_t_test$uniprot_gn_symbol[((abs(log2(df_t_test$FC)) > threshold.log2FC) & (df_t_test$Pvalue < threshold.p))]
    mycolors <- c("blue", "red", "black")
    names(mycolors) <- c("DOWN", "UP", "No")
    ggplot(df_t_test, aes(x=log2(FC), y=-log10(Pvalue),col=diff, label=delabel)) +
        geom_vline(xintercept=c(-threshold.log2FC, threshold.log2FC), col="grey") +
        geom_hline(yintercept=-log10(threshold.p), col="grey") + 
        geom_point() +
        scale_colour_manual(values = mycolors) +
        geom_text_repel() -> p
    return(p)
}
p <- DrawValcano(df_Proteins_reshape, 0.05, 1)
write.table(p$data[p$data$diff != "No",], file=paste(my_outdir, "../DE_Protein.tsv", sep = "/"),
            quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
p + 
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.line = element_line(color="black",linewidth = 0.6),
        legend.position = "none"
    ) -> p
pdf(paste(my_outdir, "DE_Protein.volcano.pdf", sep = "/"), height = 76/25.4, width=80/25.4)
print(p)
dev.off()
#指定富集分析的物种库
GO_database <- 'org.Hs.eg.db' #GO分析指定物种，物种缩写索引表详见http://bioconductor.org/packages/release/BiocViews.html#OrgDb
KEGG_database <- 'hsa' #KEGG分析指定物种，物种缩写索引表详见http://www.genome.jp/kegg/catalog/org_list.html

#gene ID转换
# gene <- bitr(df_Proteins_reshape$ensembl_gene_id , fromType = 'ENSEMBL',toType = 'ENTREZID',OrgDb = GO_database)
# bitr(result$ensembl_gene_id , fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = GO_database)
#GO富集分析
gene <- p$data[p$data$diff == "DOWN", "entrezgene_id"]
gene <- gene[!is.na(gene)]
GO <- enrichGO(
    # gene$ENTREZID,
    gene = gene,
    OrgDb = GO_database,
    keyType = "ENTREZID",
    ont = "ALL",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = T,
    )
#KEGG富集分析
KEGG<-enrichKEGG(gene,
                 organism = KEGG_database,
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
# KEGG@result %>% view()
p_GO_bar <- barplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#柱状图
p_GO_bar + 
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(family = 'sans', face= "plain", color='black', size = 2),
    axis.line = element_line(color="black",linewidth = 0.6),
    # legend.position = "none"
  ) -> p_GO_bar
pdf(paste(my_outdir, "DE_Protein.GO.bar.pdf", sep = "/"), height = 76/25.4, width=100/25.4)
print(p_GO_bar)
dev.off()

# barplot(KEGG, title = 'KEGG Pathway')
# dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")#点状图
# dotplot(KEGG)
#基因-通路关联网络图
enrichplot::cnetplot(GO,circular=TRUE,colorEdge = TRUE)
#circluar为指定是否环化，基因过多时建议设置为FALSE
# enrichplot::cnetplot(KEGG,circular=TRUE,colorEdge = TRUE)
