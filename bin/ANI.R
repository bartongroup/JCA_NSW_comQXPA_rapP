#!/usr/usr/bin/env Rscript

# Creates a heatmap of skani outputs, with associated dendrogram.

# Note that due to scaling issues the height/alignment of the dendrogram
# will require post-processing with i.e. inkscape

# Additionally creates 'classifications.txt' and 'outlier_classifications.txt' 
# files in gtdbtk results directory


suppressPackageStartupMessages({
    library(argparse)
    library(tidyr)
    library(dplyr)
    library(stringr)
    library(ggplot2)
    library(viridis)
    library(ggdendro)
    library(grid)
})

parser <- ArgumentParser()

parser$add_argument("-s", "--skani", required=TRUE, help="Path to skani .out file")
parser$add_argument("-f", "--fasta", required=TRUE, help="Path to genome fasta files")
parser$add_argument("-g", "--gtdbtk", required=TRUE, help="Path to GTDBTK results")
parser$add_argument("-o", "--out_dir", required=TRUE, help="Path to output directory")

args <- parser$parse_args()

skani_df <- read.csv(args$skani, sep='\t', skip = 1, header = FALSE) %>%
    mutate(V1 = str_replace(V1, args$fasta, '') ) %>%
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
pdf(file.path(args$out_dir, 'skani.pdf'), width = 25, height = 25)
grid.newpage()
print(heatmap, vp = viewport(x=0.4, y=0.5, width = 0.8, height = 1))
print(dendro_plot, vp = viewport(x=0.9, y=0.53, width = 0.2, height = 1))
dev.off()

outliers <- skani_long %>%
    filter(id < 97)
high_outliers <- as.data.frame(table(outliers$Genome_1)) %>%
    filter(Freq > 20)
colnames(high_outliers) <- c('accession', 'frequency')

gtdbtk_results <- read.csv(file.path(args$gtdbtk, 'gtdbtk.bac120.summary.tsv'), sep="\t") %>%
    select(c('user_genome','classification')) 

#outlier_classifications <- left_join(x = high_outliers, 
#                                     y = gtdbtk_results, 
#                                     by = join_by('accession' == 'user_genome')) %>%
#                                        select(-c(frequency))

# # finding the average of ANIs for each of these to add to the table...
#outlier_ANIs <- skani_df %>%
#    filter(Genome_1 %in% outlier_classifications$accession) %>%
#    select(-c(Genome_1)) 

#outlier_classifications$AANI <- rowMeans(outlier_ANIs)
#outlier_classifications$AANI <- formatC(signif(outlier_classifications$AANI,digits=3), digits=3,format="fg", flag="#")

write.table(gtdbtk_results, file = file.path(args$gtdbtk, 'classifications.txt'), sep="\t", row.names = FALSE, quote = FALSE )
#write.table(outlier_classifications, file = file.path(args$gtdbtk, 'outlier_classifications.txt'), sep="\t", row.names = FALSE, quote = FALSE )