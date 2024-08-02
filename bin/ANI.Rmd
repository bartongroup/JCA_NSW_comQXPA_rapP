---
title: "fastANI summary"
author: "James Abbott"
date: "18/02/2021"
output: html_document
---

```{r setup, include=FALSE}
require(knitr)
knitr::opts_chunk$set(echo = TRUE)
opts_knit$set(root.dir = "~/Projects/JCA_NSW_comP")
library(tidyverse)
library(viridis)
library(ggdendro)
library(grid)
library(readxl)
```

The fastANI output is tab delimited...just needs a bit of reformatting of the query/reference names and getting rid of columns we don't need...

```{r prep_fastANI_data}
cols = c('Query', 'Reference', 'id', 'frags', 'aligned')
ani_df <- read.csv('refined/fastANI/output.txt', sep = "\t", header = FALSE, col.names = cols)

ani_df <- ani_df %>%
  mutate(Query = str_replace(Query, 'refined/fasta/genomes/', '')) %>%
  mutate(Query = str_replace(Query, '.fasta', '')) %>%
  mutate(Reference = str_replace(Reference, 'refined/fasta/genomes/', '')) %>%
  mutate(Reference = str_replace(Reference, '.fasta', '')) %>%
  select(Query,Reference,id)

# We need a wide format for creating a matrix for the dendrogram, and long format for the heatmap..

# values_fill = 0 is due to fastANI not reporting any results < 80 %ID which could result in missing values
ani_wide = pivot_wider(ani_df, names_from = Reference, values_from = id, values_fill = 0)  %>% 
  as.data.frame()
ani_long <- pivot_longer(ani_wide, cols = -c(Query), names_to = 'Reference', values_to = 'id')
row.names(ani_wide) = ani_wide$query
ani_wide <- ani_wide %>%
  select(-Query)
```

```{r fastANI_dendrogram}
ani_matrix <- as.matrix(ani_wide)
ani_dendro <- as.dendrogram(hclust(d = dist(x = ani_matrix)))
dendro_plot <- ggdendrogram(data = ani_dendro, rotate = TRUE) +
  theme(axis.text.y = element_text(size = 2))
dendro_plot
```

```{r fastANA_heatmap}
ref_order = order.dendrogram(ani_dendro)
ref_levels = ani_long$Reference[ref_order]
ani_long$Reference <- factor(ani_long$Reference, levels = ref_levels, ordered = TRUE)
ani_long$Query <- factor(ani_long$Query, levels = ref_levels, ordered = TRUE)

heatmap<-ggplot(data = ani_long,
         aes(x = Query, y = Reference)) + 
    geom_tile(aes(fill = id)) +
    scale_fill_viridis(limits = c(80, 100)) +
    theme_bw() +
    theme(axis.text.x = element_text(size = 4, angle = 45, vjust=1, hjust=1),
          axis.text.y = element_text(size = 4),
          legend.position = 'left') 
```

```{r combine}
pdf('refined/fastANI/ANI.pdf', width = 25, height = 25)
grid.newpage()
print(heatmap, vp = viewport(x=0.4, y=0.5, width = 0.8, height = 1))
print(dendro_plot, vp = viewport(x=0.9, y=0.53, width = 0.2, height = 1))
dev.off()
```
We have some outliers - let's work out how they are...
```{r outliers}
outliers <- ani_df %>%
  filter(id < 90)
high_outliers <- as.data.frame(table(outliers$Query)) %>%
  filter(Freq > 10)
high_outliers
```
Just for the sake of completion and to ensure these are real results, skani was also run...

```{r skani_data_prep}
skani_df <- read.csv('refined/skani.out', sep='\t', skip = 1, header = FALSE) %>%
    mutate(V1 = str_replace(V1, 'refined/fasta/genomes/', '') ) %>%
    mutate(V1 = str_replace(V1, '.fasta', ''))
  
accessions <- skani_df$V1
colnames(skani_df) <- c('Genome_1', unlist(accessions))
rownames(skani_df) <- accessions
skani_df <- skani_df
skani_long <- pivot_longer(skani_df, cols = -c('Genome_1'), names_to = 'Genome_2', values_to = 'id')
skani_wide <- skani_df %>%
  select(-c('Genome_1'))
```

```{r skani_dendrogram}
ani_matrix <- as.matrix(skani_wide)
ani_dendro <- as.dendrogram(hclust(d = dist(x = ani_matrix)))
dendro_plot <- ggdendrogram(data = ani_dendro, rotate = TRUE, leaf_labels = FALSE, labels = FALSE) +
  theme(axis.text.y = element_text(size = 2))
dendro_plot
```
```{r skani_heatmap}
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
```

```{r skani_combine}
pdf('refined/skani.pdf', width = 25, height = 25)
grid.newpage()
print(heatmap, vp = viewport(x=0.4, y=0.5, width = 0.8, height = 1))
print(dendro_plot, vp = viewport(x=0.9, y=0.53, width = 0.2, height = 1))
dev.off()
```

We clearly have some outliers here - some are clearly still closely related, but others are <80% ANI. GTDBTK has therefore been used to confirm their identities. Merging these results with our existing isolate details hopefully will show what is going on.

```{r skani_outliers}
outliers <- skani_long %>%
  filter(id < 97)
high_outliers <- as.data.frame(table(outliers$Genome_1)) %>%
  filter(Freq > 20)
colnames(high_outliers) <- c('accession', 'frequency')

bs_isolates <- read_xlsx("Bacillus_subtilis_complete_genomes_25-06-2024.xlsx") %>%
  select(c('accession','isolate'))

gtdbtk_results <- read.csv('gtdbtk/results/classify/gtdbtk.bac120.summary.tsv', sep="\t") %>%
  select(c('user_genome','classification')) %>%
  left_join(y= bs_isolates, by = join_by('user_genome' == 'accession')) 

outlier_classifications <- left_join(x = high_outliers, y = gtdbtk_results, by = join_by('accession' == 'user_genome')) %>%
  select(-c(frequency))

# finding the average of ANIs for each of these to add to the table...
outlier_ANIs <- skani_df %>%
  filter(Genome_1 %in% outlier_classifications$accession) %>%
  select(-c(Genome_1)) 
  
outlier_classifications$AANI <- rowMeans(outlier_ANIs)
outlier_classifications$AANI <- formatC(signif(outlier_classifications$AANI,digits=3), digits=3,format="fg", flag="#")

write.table(gtdbtk_results, file = 'gtdbtk/classifications.txt', sep="\t", row.names = FALSE, quote = FALSE )
write.table(outlier_classifications, file = 'gtdbtk/outlier_classifications.txt', sep="\t", row.names = FALSE, quote = FALSE )
```