library(tidyverse)
# load input list file
classes <- c('R','K','D','P','C','O','F','G','S')
t <- data.table::fread("All.list", header=TRUE)
table(t$class)
filter(t, class %in% classes) %>% select(class) %>% table()
# get microbe species tax matrix
row.table <- t[,1:3]
head(row.table)
# split to Normal and Tumor group
normal.mat <- t %>% select(-ends_with("Tumor"))
tumor.mat <- t %>% select(-ends_with("Normal"))
# exclude microbes exist in any Normal sample
tumor.only.id <- normal.mat$id[apply(normal.mat[,-c(1:3)], 1, sum) == 0]
# exclude microbes whose read counts are lower than lower quantile in each sample
juge <- function(x){
  tmp <- x
  tmp[which(x <= quantile(x[which(x>0)])[2])] <- 0
  return(tmp)
}
tumor.mat <- apply(tumor.mat[,-c(1:3)], 2, juge)
tumor.mat <- cbind(row.table, tumor.mat) %>% filter(id %in% tumor.only.id)
# exclude microbes doesn't appeared in all tumor samples
tumor.mat <- tumor.mat[apply(tumor.mat[,-c(1:3)],1,sum) > 0,]
tumor.only.id <- tumor.mat$id
# species level heat map
library(pheatmap)
library(RColorBrewer)
# tumor.mat.plot <- tumor.mat %>% filter(class=="S")
for (each in classes[-c(1:6)]){
  tumor.mat.plot <- tumor.mat[tumor.mat[[3]] == each,]
  name.row <- tumor.mat.plot$name
  tumor.mat.plot <- select(tumor.mat.plot, ends_with("Tumor"))
  tumor.mat.plot[tumor.mat.plot==0] <- NA
  file_name = paste0("selected_species_",each,".pdf")
  # pdf(file="selected_species.pdf", width=15, height=7)
  pdf(file=file_name, width=15, height=7)
  p <- pheatmap(t(tumor.mat.plot), border_color = "white",
           color = colorRampPalette(brewer.pal(n=5, name = "Reds"))(9),
           kmeans_k = NA, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
           cellwidth = 4, cellheight = 4, labels_row = colnames(tumor.mat.plot), labels_col = name.row,
           fontsize_col = 4, fontsize_row = 5, angle_col = 90, na_col= "white", main=each)
  print(p)
  dev.off()
}
# selected microbes in each sample
select.microbe <- function(x){
  tmp <- cbind(row.table[row.table$id %in% tumor.only.id[which(x>0)],], x[which(x >0)])
  colnames(tmp)[4] <-  "read_counts"
  return(tmp)
}
final_result <- lapply(tumor.mat[,-c(1:3)], select.microbe)
saveRDS(final_result, file = "selected_microbes.rds")
