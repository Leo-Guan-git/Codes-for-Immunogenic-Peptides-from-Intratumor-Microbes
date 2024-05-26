library(tidyverse)
## get all sample results
taxonomy <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
get_data <- function(filename, id){
  t <- data.table::fread(filename, header = FALSE)
  colnames(t) <- c('clade_tree','clade_name', 'NCBI_tax_id', id)
  # t <- separate(t, clade_name, sep="\\|?[kpcofgs]__",
  #               into= c(NA, taxonomy))
  return(t)
}

# metaphlan_file <- list.files(getwd())
metaphlan_file <- commandArgs(T)
# metaphlan_file <- metaphlan_file[grep("mpa.v2.txt", metaphlan_file)]
ids <- sapply(basename(metaphlan_file), strsplit, "[_.]")
id.subsample <- function(x){x[1]}
samples <- levels(as.factor(as.vector(unlist(lapply(ids, id.subsample)))))
for(i in 1:length(metaphlan_file)){
  data <- get_data(metaphlan_file[i], paste(ids[[i]][1],ids[[i]][2],sep='.'))
  if(i == 1){
    all_data <- data
  } else {
    # all_data <- dplyr::full_join(all_data, data, by=all_of(c(taxonomy,'NCBI_tax_id'))) %>%
    #   dplyr::arrange(across(all_of(taxonomy),desc))
    all_data <- dplyr::full_join(all_data, data)
  }
}
# all_data <- all_data[apply(all_data[,-c(1:8)],1,sum) >0,]
all_data <- all_data[apply(all_data[,-c(1:3)],1,sum) >0,]
all_data <- all_data[grep("[dk]__(Bacteria)|(Fungi)|(Archaea)|(Viruses)", all_data$clade_tree),]

tmp <- all_data[,-c(1:3)]
tmp[tmp==0.0] <- NA
all_data <- cbind(all_data[,c(1:3)], tmp)
rm(tmp)
all_species_data <- all_data[grep("^s__", all_data$clade_name),]
write.csv(all_data, "metagenome.csv")
write.csv(all_species_data, "metagenome_species.csv")

## generate heatmap of all taxnomys
library(pheatmap)
library(RColorBrewer)
taxonomy.class <- sapply(taxonomy, function(x){tolower(substr(x,1,1))})
# library(viridis)
pheatmap_draw <- function(input_data, output.name.prefix="all", class=taxonomy.class){
  for(i in 1:length(class)){
    if (i == 1){
      # data <- filter(input_data, !is.na(.data[[class[[i]]]]) & is.na(.data[[class[[i+1]]]])) %>%
      #   # as_tibble() %>%
      #   select(-class[-i], -NCBI_tax_id)
      data <- input_data[grep("[dk]__(Bacteria)|(Fungi)|(Archaea)|(Viruses)", input_data$clade_name),-c(1,3)]
    }else{
      # data <- filter(input_data, !is.na(.data[[class[[i]]]])) %>%
      #   # as_tibble() %>%
      #   select(-class[-i], -NCBI_tax_id)
      data <- input_data[grep(paste0("^",class[[i]], "__"),input_data$clade_name),-c(1,3)]
    }
    my_label <- sapply(data[[1]], substring, 4)
    data <- data[,-1]
    if (nrow(data) > 0){
      # data <- as.data.frame(apply(data, 2,jitter))
      p <- pheatmap(t(data), border_color = "white",
                    color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(100),
                    # color = rev(viridis(7,begin=0.3,end=0.8,direction=-1,option="D" )),
                    kmeans_k = NA, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
                    labels_col = gsub("_", " ", my_label), cellwidth = 8, cellheight = 10,
                    fontsize_col = 8, fontsize_row = 8, angle_col = 90, main=taxonomy[i],na_col= "#FFFFFF")
      pdf(paste(output.name.prefix,taxonomy[i],"relative_abundance.heatmap.pdf", sep='.'), 70,14)
      print(p)
      dev.off()
    }
  }
}
pheatmap_draw(all_data)

## select microbiomes for each sample
selected_microbiome <- list()
for(i in samples){
  tmp <- all_data %>% 
    # select(all_of(taxonomy), NCBI_tax_id,starts_with(i)) %>%
    select(clade_tree,clade_name,NCBI_tax_id,starts_with(i)) %>%
    filter(!is.na(.data[[paste(i,'Tumor',sep='.')]])) %>%
    filter(is.na(.data[[paste(i,'Normal',sep='.')]]) | ((.data[[paste(i,'Tumor',sep='.')]]/.data[[paste(i,'Normal',sep='.')]]) > 2 & .data[[paste(i,'Tumor',sep='.')]] > quantile(.data[[paste(i,'Tumor',sep='.')]], na.rm=TRUE)[2]))
  selected_microbiome[[i]] <- tmp
# generate heatmap of all taxnomys for each sample
  # pheatmap_draw(input_data = t, output.name.prefix = i)
  t.plot <- tmp[,-c(1:3)]
  my_label <- tmp[[2]]
  if (nrow(t.plot) > 0){
    # data <- as.data.frame(apply(data, 2,jitter))
    # p <- pheatmap(t(t.plot), border_color = "white",
    p <- pheatmap(t.plot, border_color = "white",
                  color = colorRampPalette(rev(brewer.pal(n=7, name = "RdYlBu")))(100),
                  # color = rev(viridis(7,begin=0.3,end=0.8,direction=-1,option="D" )),
                  kmeans_k = NA, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
                  # labels_col = my_label, cellwidth = 8, cellheight = 10,
                  cellwidth = 8, cellheight = 10, labels_row = my_label, labels_col = c('Normal', 'Tumor'),
                  fontsize_col = 8, fontsize_row = 8, angle_col = 90, main=i,na_col= "#FFFFFF")
    pdf(paste(i,"relative_abundance.heatmap.pdf", sep='.'), 28,56)
    print(p)
    dev.off()
  }
  write.csv(tmp, paste(i,"selected.microbime.csv", sep='.'))
}
save(all_data, all_species_data, samples, taxonomy, selected_microbiome, file = "kra2mpa_results.Rdata")

selected_species <- lapply(selected_microbiome,function(x){x[grep("^s__",x[[2]]),]})
for (i in names(selected_species)){
  if(i=="CRC01"){
    draw_data <- selected_species[[i]]
  }else{
    draw_data <- dplyr::full_join(draw_data, selected_species[[i]])
  }
}

p <- pheatmap(t(draw_data[,-c(1,2,3)]), border_color = "white",
              color = colorRampPalette(rev(brewer.pal(n=10, name = "RdYlBu")))(100),
              # color = rev(viridis(7,begin=0.3,end=0.8,direction=-1,option="D" )),
              kmeans_k = NA, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none",
              labels_col = gsub("_", " ", substring(draw_data[[2]],4)), cellwidth = 8, cellheight = 10,
              fontsize_col = 8, fontsize_row = 8, angle_col = 90, main="Species",na_col= "#FFFFFF")
pdf("selected_species.heatmap.pdf", 56,14)
print(p)
dev.off()
table(selected_species$detected_sample_num)
# focus_species <- filter(selected_species, detected_sample_num >3) %>% select(Species,detected_sample_num)
# focus_species <- focus_species %>% 
#   left_join(select(all_species_data, colnames.selected_microbiome[1:8])) %>%
#   select(colnames.selected_microbiome[1:8],detected_sample_num)
  