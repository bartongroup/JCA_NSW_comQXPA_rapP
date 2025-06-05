#!/usr/bin/env Rscript

# Creates a heatmap of skani outputs, with associated dendrogram.

# Note that due to scaling issues the height/alignment of the dendrogram
# will require post-processing with i.e. inkscape

renv::restore()

suppressPackageStartupMessages({
    library(tidyr)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(viridis)
    library(ggdendro)
    library(grid)
})

genomes = "data/full/fasta/genomes"
skani_output = "results/skani/skani.out"
output_dir = "results/skani"

skani_df <- read.csv(skani_output, sep='\t', skip = 1, header = FALSE) %>%
    mutate(V1 = str_replace(V1, genomes, '') ) %>%
    mutate(V1 = str_replace(V1, '.fasta', ''))

accessions <- skani_df$V1
colnames(skani_df) <- c('Genome_1', unlist(accessions))
rownames(skani_df) <- accessions
skani_long <- pivot_longer(skani_df, cols = -c('Genome_1'), names_to = 'Genome_2', values_to = 'id')
skani_wide <- skani_df %>%
    select(-c('Genome_1'))

ani_matrix <- as.matrix(skani_wide)
ani_dendro <- as.dendrogram(hclust(d = dist(x = ani_matrix)))
dendro_plot <- ggdendrogram(data = ani_dendro, rotate = TRUE, leaf_labels = FALSE, labels = FALSE) +
    theme(axis.text.y = element_text(size = 2))

# Determine sequence order for plotting heatmap based upon dendrogram
ref_order = order.dendrogram(ani_dendro)
ref_levels = skani_long$Genome_2[ref_order]
skani_long$Genome_1 <- factor(skani_long$Genome_1, levels = ref_levels, ordered = TRUE)
skani_long$Genome_2 <- factor(skani_long$Genome_2, levels = ref_levels, ordered = TRUE)

heatmap<-ggplot(data = skani_long,
         aes(x = Genome_1, y = Genome_2)) + 
         geom_tile(aes(fill = id)) +
         scale_fill_viridis(limits = c(80, 100)) +
         theme_bw() +
         theme(axis.text.x = element_text(size = 4, angle = 45, vjust=1, hjust=1),
             axis.text.y = element_text(size = 4),
             legend.position = 'left')

# Output combined plot of heatmap and dendrogram
pdf(file.path(output_dir, 'skani.pdf'), width = 25, height = 25)
grid.newpage()
print(heatmap, vp = viewport(x=0.4, y=0.5, width = 0.8, height = 1))
print(dendro_plot, vp = viewport(x=0.9, y=0.53, width = 0.2, height = 1))
dev.off()
