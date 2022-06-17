# A script to extract marker genes with Seurat.
library(Matrix)
library(Seurat)
library(tidyverse)
# Load ref data set, should have been normalized by sequencing depth.
hypo.ref <- readRDS('data/ref.cam.hypo.rds') 
table(hypo.ref$taxonomy_lvl2)

lvl1_clusters = names(table(hypo.ref$taxonomy_lvl2))
registerDoParallel(8)
hypo_lvl1_markers = foreach(x = lvl1_clusters) %dopar% {
  FindMarkers(hypo.ref, ident.1 = x, group.by = 'taxonomy_lvl2')
}
names(hypo_lvl1_markers) = lvl1_clusters
head(hypo_lvl1_markers$Astrocyte[order(hypo_lvl1_markers$Astrocyte$avg_log2FC,decreasing = T),])
save(hypo_lvl1_markers, file = '~/fkc259/cell_type/testing/ref_data_marker_genes/marker_db_cam_hypo_seurat.Rdata')

  