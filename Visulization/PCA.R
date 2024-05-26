library(ggplot2)
library(ggrepel)
library(reshape2)
library(RColorBrewer)

df <- read.csv(normalizePath("G:/Micro_Maxquant_v2/results/CRC.all.abundance.tsv"),sep = '\t', header = T)
df$Sample_type <- paste(df$Sample, df$dtype, sep = '_')
df$Abundance <- as.numeric(df$Abundance)
# df <- df[,c("Sample","Abundance","Name")]

df_dcast <- dcast(select(df, Sample_type, Abundance, Name), Sample_type ~ Name, value.var = "Abundance")
df_dcast <- t(na.omit(t(df_dcast))) %>% as.data.frame()

pca1 <- prcomp(sapply(df_dcast[,-1], as.numeric), center = TRUE, scale. = TRUE)
df1 <- as.data.frame(pca1$x) # 提取PC score # 注意：如果不转成数据框形式后续绘图时会报错
summ1 <- summary(pca1)
# color_list <- setNames(c("#E41A1C", "#E41A1C", "#4A72A6", "#4A72A6", "#48A462" , "#48A462" , "#7E6E85", "#7E6E85", "#D16948", "#D16948", "#FFB716", "#FFB716", "#E1C62F", "#E1C62F", "#B75F49", "#B75F49", "#EC83BA", "#EC83BA", "#999999", "#999999"), rownames(df))
# df1$Shape <- ifelse(grepl("Tumor", rownames(df1)), "Tumor", "Normal")

df1$sample <- sapply(strsplit(df_dcast$Sample_type, "_"), "[", 1)
df1$dtype <- sapply(strsplit(df_dcast$Sample_type, "_"), "[", 2)

ggplot(df1, aes(x = PC1, y = PC2, color = dtype)) +
    geom_point() +
    stat_ellipse(
        aes(fill=dtype),
        color = NA,
        type = "t",
        level = 0.95,
        alpha = 0.2,
        geom = "polygon",
        ) +
    labs(
        x = "PC1",
        y = "PC2",
        title = "PCA Scores with Confidence Intervals",
        )
outdir <- "G:/Micro_Maxquant_v2/results/visulization/PCA"
prefix <- "test"
# 绘制PCA得分图
for (i in 1:4) {
    for (j in (i+1):5) {
        x_axis <- paste0("PC",i)
        y_axis <- paste0("PC",j)
        
        xlab1 <- paste0(paste0("PC",i,"("),round(summ1$importance[2,i]*100,2),"%)")
        ylab1 <- paste0(paste0("PC",j,"("),round(summ1$importance[2,j]*100,2),"%)")
        
        # outfile <- paste(outdir, paste(prefix,x_axis,y_axis,'png', sep='.'), sep = '/')
        # png(file=outfile)
        
        p.pca1 <- ggplot(data = df1,aes(x = .data[[x_axis]], y = .data[[y_axis]], color = dtype))+
            stat_ellipse(aes(fill = dtype),
                         type = "t",geom = "polygon",alpha = 0.1,color = NA,)+ 
            # geom_point(aes(shape = Shape),size = 3.5)+
            geom_point(size = 2) +
            labs(x = xlab1,y = ylab1,color = "Type", title = "PCA Scores Plot")+
            # guides(fill = "")+
            theme_bw()+
            # scale_fill_manual(values = color_list)+
            # scale_colour_manual(values = color_list)+
            theme(plot.title = element_text(hjust = 0.5,size = 15),
                  axis.text = element_text(size = 11),axis.title = element_text(size = 13),
                  legend.text = element_text(size = 11),legend.title = element_text(size = 13),
                  plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
        print(p.pca1)
        # dev.off()
    }
}

i <- 1
j <- 2
x_axis <- paste0("PC",i)
y_axis <- paste0("PC",j)

xlab1 <- paste0(paste0("PC",i,"("),round(summ1$importance[2,i]*100,2),"%)")
ylab1 <- paste0(paste0("PC",j,"("),round(summ1$importance[2,j]*100,2),"%)")

outdir <- "G:/Micro_Maxquant_v2/results/visulization/PCA"
prefix <- "CRC"
outfile <- paste(outdir, paste(prefix,x_axis,y_axis,'pdf', sep='.'), sep = '/')

p.pca1 <- ggplot(data = df1,aes(x = .data[[x_axis]], y = .data[[y_axis]], color = dtype))+
    stat_ellipse(aes(fill = dtype),
                 type = "t",geom = "polygon",alpha = 0.1,color = NA,)+ 
    geom_point() +
    labs(x = xlab1,y = ylab1, color = "Type")+
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0),'cm'),
        # legend.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        # legend.position = "bottom",
        # legend.key.size = unit(0.2, "cm"),
        # plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'),
        )
pdf(file=outfile, height = 57/25.4, width=70/25.4)
print(p.pca1)
dev.off()

df_BetaDiver <- read.delim(normalizePath("G:/Micro_Maxquant_v2/results/CRC.microorganism.s.beta_diversity.txt"), header=T, sep="\t")
df_BetaDiver_dcast <- dcast(df_BetaDiver, sample1 ~ sample2, value.var = "bray_curtis")
df_BetaDiver_dcast <- t(na.omit(t(df_BetaDiver_dcast))) %>% as.data.frame()

pca1_beta <- prcomp(sapply(df_BetaDiver_dcast[,-1], as.numeric), center = TRUE, scale. = TRUE)
df1 <- as.data.frame(pca1_beta$x) # 提取PC score # 注意：如果不转成数据框形式后续绘图时会报错
summ1 <- summary(pca1_beta)
df1$sample <- sapply(strsplit(df_BetaDiver_dcast$sample1,"_"), "[", 1)
df1$dtype <- sapply(strsplit(df_BetaDiver_dcast$sample1,"_"), "[", 2)
ggplot(df1, aes(x = PC1, y = PC2, color = dtype)) +
    geom_point() +
    stat_ellipse(
        aes(fill=dtype),
        color = NA,
        type = "t",
        level = 0.95,
        alpha = 0.2,
        geom = "polygon",
    ) +
    labs(
        x = "PC1",
        y = "PC2",
        title = "PCA Scores with Confidence Intervals",
    )

i <- 1
j <- 2
x_axis <- paste0("PC",i)
y_axis <- paste0("PC",j)

xlab1 <- paste0(paste0("PC",i,"("),round(summ1$importance[2,i]*100,2),"%)")
ylab1 <- paste0(paste0("PC",j,"("),round(summ1$importance[2,j]*100,2),"%)")

outdir <- "G:/Micro_Maxquant_v2/results/visulization/PCA"
prefix <- "CRC.Beta_Diversity.bray_curtis"
outfile <- paste(outdir, paste(prefix,x_axis,y_axis,'pdf', sep='.'), sep = '/')

p.pca1 <- ggplot(data = df1,aes(x = .data[[x_axis]], y = .data[[y_axis]], color = dtype))+
    stat_ellipse(aes(fill = dtype),
                 type = "t",geom = "polygon",alpha = 0.1,color = NA,)+ 
    geom_point() +
    labs(x = xlab1,y = ylab1, color = "Type")+
    scale_fill_manual(values = c("#F3756D", "#33cccc")) +
    scale_color_manual(values = c("#F3756D", "#33cccc")) +
    theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        axis.line = element_line(colour = "black", linewidth = 0.4),
        legend.position = "bottom",
        plot.margin = unit(c(0,0,0,0),'cm'),
        # legend.text = element_text(family = 'sans', face= "plain", color='black', size = 6),
        # legend.position = "bottom",
        # legend.key.size = unit(0.2, "cm"),
        # plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'),
    )
pdf(file=outfile)
print(p.pca1)
dev.off()

library(rgl)
open3d()
plot3d(df1$PC1, df1$PC2, df1$PC3, 
       col = c("red", "blue"), 
       xlab = "PC1", ylab = "PC2", zlab = "PC3", 
       size = 5,
       main = "PCA 3D Plot")

library(scatterplot3d)
color = rep(c("red", "blue"), 10)
scatterplot3d(
    df1[,1:3],
    color = rep(c("red", "blue"), 10),
    pch = 16,
    angle = 40,
    box = T,
    type = "p",
    lty.hide = 2,
    lty.grid = 2,
    )
