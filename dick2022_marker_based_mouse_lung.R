# ================================================================
# Marker-based macrophage analysis pipeline (replicating Dick et al.)
# ================================================================
# Author: Jingshuo Li
# Version: 1.0
# ================================================================

# -------------------- Parameters (edit here) --------------------
project_name <- "MouseLung"
data_dir     <- "data/mouse_lung"                  # raw 10X matrix input folder
dag_path     <- "data/DAG_136_MouseGenes.csv"      # DAG list (van den Brink et al. 2017, first column = gene symbols)
results_dir  <- "results/mouse_lung"               # output root

# Create main results directory
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Libraries --------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(openxlsx)
  library(cluster)
  library(plyr)
})

set.seed(111)

# -------------------- Step 0: Load GEO metadata (optional) --------------------
# Note: This step retrieves metadata from GEO for reference purposes only.
# Mouse lung dataset – replace GSM ID if needed
gsm_mouse_lung <- tryCatch({
  getGEO("GSM5687568", GSEMatrix = FALSE)
}, error = function(e) { message("GEO fetch skipped (no internet in some environments)."); NULL })
if (!is.null(gsm_mouse_lung)) utils::str(gsm_mouse_lung)

# -------------------- Step 1: Load raw 10X gene expression data --------------------
mouse_lung <- Read10X(data.dir = data_dir) %>%
  CreateSeuratObject(project = project_name, min.cells = 3, min.features = 200)

# -------------------- Step 2: Load mouse DAG gene list --------------------
dag_df    <- read.csv(dag_path, stringsAsFactors = FALSE)
dag_genes <- unique(dag_df[[1]])
overlap_dag <- intersect(rownames(mouse_lung), dag_genes)

# -------------------- Step 3: Add QC metrics --------------------
# Note: mouse mitochondrial genes use "mt-"
mouse_lung[["percent.mito"]] <- PercentageFeatureSet(mouse_lung, pattern = "^mt-")
mouse_lung[["percent.dag"]]  <- PercentageFeatureSet(mouse_lung, features = overlap_dag)

# Summaries
print(summary(mouse_lung$percent.mito))
print(summary(mouse_lung$percent.dag))

# QC histograms
qc_hist_dir <- file.path(results_dir, "QC_Histograms")
dir.create(qc_hist_dir, showWarnings = FALSE, recursive = TRUE)

png(file.path(qc_hist_dir, "Histogram_percent.mito.png"), width = 800, height = 600)
hist(mouse_lung$percent.mito, breaks = 50, col = "lightblue", main = "Mito %", xlab = "percent.mito")
dev.off()

png(file.path(qc_hist_dir, "Histogram_percent.dag.png"), width = 800, height = 600)
hist(mouse_lung$percent.dag, breaks = 50, col = "pink", main = "DAG %", xlab = "percent.dag")
dev.off()

# -------------------- Step 4: Save QC scatter plots --------------------
qc_plot_dir <- file.path(results_dir, "QC_Plots")
dir.create(qc_plot_dir, showWarnings = FALSE, recursive = TRUE)

plot1 <- FeatureScatter(mouse_lung, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(mouse_lung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(mouse_lung, feature1 = "nCount_RNA", feature2 = "percent.dag")
plot4 <- FeatureScatter(mouse_lung, feature1 = "percent.dag", feature2 = "percent.mito")

ggsave(file.path(qc_plot_dir, "Scatter_nCount_vs_percent.mito.png"), plot1, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(qc_plot_dir, "Scatter_nCount_vs_nFeature.png"),  plot2, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(qc_plot_dir, "Scatter_nCount_vs_percent.dag.png"),plot3, width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(qc_plot_dir, "Scatter_percent.dag_vs_percent.mito.png"), plot4, width = 8, height = 8, dpi = 300, bg = "white")

rm(plot1, plot2, plot3, plot4); gc()

# -------------------- Step 5: Final QC filtering --------------------
mouse_lung_filtered <- subset(
  mouse_lung,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mito   < 15 &
    percent.dag    < 5
)

saveRDS(mouse_lung_filtered, file = file.path(results_dir, "mouse_lung_filtered.rds"))

# -------------------- Step 6: Normalize data --------------------
mouse_lung_filtered <- NormalizeData(mouse_lung_filtered)

# Seurat v5: join layers
mouse_lung_filtered <- JoinLayers(mouse_lung_filtered)

# -------------------- Step 7: FindVariableFeatures --------------------
mouse_lung_meanvar <- FindVariableFeatures(
  mouse_lung_filtered,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.005, 4),
  dispersion.cutoff = c(0.25, Inf)
)
message("Variable features (mean.var.plot): ", length(VariableFeatures(mouse_lung_meanvar)))

output_dir_meanvar <- file.path(results_dir, "Method_meanvar")
dir.create(output_dir_meanvar, recursive = TRUE, showWarnings = FALSE)

# -------------------- Step 8: Scale and PCA --------------------
mouse_lung_meanvar <- ScaleData(mouse_lung_meanvar, vars.to.regress = c("nCount_RNA", "percent.mito"))
mouse_lung_meanvar <- RunPCA(mouse_lung_meanvar, features = VariableFeatures(mouse_lung_meanvar))

# -------------------- Step 9: Elbow plot --------------------
p_elbow_meanvar <- ElbowPlot(mouse_lung_meanvar) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "ElbowPlot_meanvar.png"),
       plot = p_elbow_meanvar, width = 6, height = 5, dpi = 300)
rm(p_elbow_meanvar)

# -------------------- Step 10: Natural clustering (no forcing) --------------------
mouse_lung_meanvar <- FindNeighbors(mouse_lung_meanvar, dims = 1:20)
mouse_lung_meanvar <- FindClusters(mouse_lung_meanvar, resolution = 1.0)

# set identities to seurat clusters
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$seurat_clusters

# -------------------- Rename clusters to start from 1 --------------------
cluster_meanvar_ids <- Idents(mouse_lung_meanvar)
levels(cluster_meanvar_ids) <- as.character(as.numeric(levels(cluster_meanvar_ids)) + 1)
mouse_lung_meanvar$Cluster1based <- cluster_meanvar_ids
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$Cluster1based

# -------------------- Step 11: Run UMAP --------------------
mouse_lung_meanvar <- RunUMAP(mouse_lung_meanvar, reduction = "pca", dims = 1:20)

# -------------------- Step 12: Run tSNE --------------------
mouse_lung_meanvar <- RunTSNE(mouse_lung_meanvar, reduction = "pca", dims = 1:20)

# -------------------- Step 13: Save renamed UMAP --------------------
p_umap_meanvar <- DimPlot(mouse_lung_meanvar, reduction = "umap", label = TRUE) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "UMAP_FirstClusters_meanvar.png"),
       plot = p_umap_meanvar, width = 8, height = 6, dpi = 300)
rm(p_umap_meanvar)

# -------------------- Step 14: Save renamed tSNE --------------------
p_tsne_meanvar <- DimPlot(mouse_lung_meanvar, reduction = "tsne", label = TRUE) + theme_classic()
ggsave(filename = file.path(output_dir_meanvar, "TSNE_FirstClusters_meanvar.png"),
       plot = p_tsne_meanvar, width = 8, height = 6, dpi = 300)
rm(p_tsne_meanvar)

# -------------------- Step 15: Silhouette Score Calculation --------------------
# 1. Extract PCA embeddings
emb_pca <- Embeddings(mouse_lung_meanvar, reduction = "pca")[, 1:20]
# 2. Cluster labels
cluster_labels <- as.factor(mouse_lung_meanvar$Cluster1based)
# 3. Distances
dists <- dist(emb_pca)
# 4. Silhouette
sil <- silhouette(as.integer(cluster_labels), dists)
# 5. Data frame
sil_df <- as.data.frame(sil)
sil_df$Cluster <- cluster_labels[as.numeric(rownames(sil_df))]
# 6. Per-cluster mean
sil_summary <- sil_df %>%
  group_by(Cluster) %>%
  summarise(mean_sil = mean(sil_width)) %>%
  arrange(desc(mean_sil))
print(sil_summary)

# 7. Plot
output_dir_sil <- file.path(results_dir, "SilhouettePlot")
dir.create(output_dir_sil, showWarnings = FALSE, recursive = TRUE)

p_sil <- ggplot(sil_df, aes(x = Cluster, y = sil_width, fill = Cluster)) +
  geom_violin(trim = FALSE, alpha = 0.4) +
  geom_jitter(width = 0.2, size = 0.3, alpha = 0.5) +
  theme_classic() +
  labs(y = "Silhouette width", title = "Silhouette Score per Cluster")

ggsave(filename = file.path(output_dir_sil, "Silhouette_PerCluster.png"),
       plot = p_sil, width = 8, height = 6, dpi = 300)

# -------------------- Step 16-A.1: Detect cluster markers --------------------
Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$Cluster1based

markers_all <- data.frame()
for (clust in levels(Idents(mouse_lung_meanvar))) {
  cat("Finding markers for cluster", clust, "...\n")
  markers <- FindMarkers(
    mouse_lung_meanvar,
    ident.1 = clust,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.2
  )
  if (nrow(markers) > 0) {
    markers$cluster <- clust
    markers$gene    <- rownames(markers)
    markers_all     <- rbind(markers_all, markers)
  }
}

# -------------------- Step 16-A.2: Filter significant markers --------------------
markers_sig <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^mt-", gene, ignore.case = TRUE)) %>%
  filter(!grepl("^Rp[sl]", gene, ignore.case = TRUE))

top30_markers <- markers_sig %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(markers_all,  file.path(output_dir_meanvar, "DEGs_ClusterMarkers_meanvar.xlsx"), rowNames = FALSE)
write.xlsx(top30_markers, file.path(output_dir_meanvar, "Top30Markers_meanvar.xlsx"), rowNames = FALSE)

# -------------------- Step 16-B: Assign Cell Type Labels Based on Marker Expression --------------------
step16b_dir <- file.path(output_dir_meanvar, "ManualAnnotation_From_TopMarkers")
dir.create(step16b_dir, showWarnings = FALSE, recursive = TRUE)

# Define lung-specific marker list
celltype_markers <- list(
  # --- Myeloid lineage ---
  "Alveolar_Macrophages"   = c("Marco", "Siglec5", "Car4", "Chil3", "Pparg", "Ear1", "Ear2"),
  "Interstitial_Macrophages" = c("C1qa", "Trem2", "Cd68", "Mrc1", "Fcgr1"),
  "Monocytes"              = c("Ly6c2", "Sell", "Plac8", "Il1b", "Cd14"),
  "DCs"                    = c("Flt3", "Xcr1", "Cd24a", "Irf8", "Clec9a"),
  # --- Lymphoid ---
  "T_cells" = c("Cd3d", "Cd3e", "Il7r", "Cd8a", "Trac"),
  "B_cells" = c("Cd79a", "Ms4a1", "Cd19", "Cd22"),
  "NK_cells"= c("Nkg7", "Klrb1c", "Gzmb", "Klrd1", "Prf1"),
  # --- Granulocytes ---
  "Neutrophils" = c("S100a8", "S100a9", "Ngp", "Elane"),
  # --- Epithelial (Lung specific) ---
  "AT1_cells"   = c("Ager", "Pdpn", "Hopx", "Aqp5"),
  "AT2_cells"   = c("Sftpc", "Sftpb", "Napsa", "Lamp3"),
  "Ciliated_cells" = c("Foxj1", "Dynlrb2", "Dnah5"),
  "Club_cells"  = c("Scgb1a1", "Cyp2f2", "Scgb3a2"),
  # --- Structural (non-immune) ---
  "Endothelial" = c("Cdh5", "Pecam1", "Stab2", "Emcn"),
  "Fibroblasts" = c("Col1a1", "Pdgfra", "Dcn", "Lum")
)

# AverageExpression in Seurat v5 returns a list; pick RNA matrix
all_marker_genes <- unique(unlist(celltype_markers))
genes_use <- intersect(all_marker_genes, rownames(mouse_lung_meanvar))
avg_list <- AverageExpression(mouse_lung_meanvar, assays = "RNA", slot = "data", features = genes_use)
cluster_expr <- as.data.frame(avg_list$RNA)  # genes x clusters

cluster_ids <- colnames(cluster_expr)
cluster_labels <- rep("Unassigned", length(cluster_ids))

for (i in seq_along(cluster_ids)) {
  clust <- cluster_ids[i]
  avg_scores <- sapply(celltype_markers, function(glist) {
    matched <- intersect(glist, rownames(cluster_expr))
    if (length(matched) == 0) return(NA_real_)
    mean(cluster_expr[matched, clust], na.rm = TRUE)
  })
  cluster_labels[i] <- if (all(is.na(avg_scores))) "Unassigned" else names(which.max(avg_scores))
}

annotation_df <- data.frame(
  Cluster  = cluster_ids,
  CellType = cluster_labels,
  stringsAsFactors = FALSE
)
write.xlsx(annotation_df, file.path(step16b_dir, "Cluster_Celltype_Annotation.xlsx"), row.names = FALSE)

# Map annotation back to Seurat object
mouse_lung_meanvar$celltype <- plyr::mapvalues(
  x = mouse_lung_meanvar$Cluster1based,
  from = annotation_df$Cluster,
  to   = annotation_df$CellType,
  warn_missing = FALSE
)

Idents(mouse_lung_meanvar) <- mouse_lung_meanvar$celltype

p_anno <- DimPlot(mouse_lung_meanvar, group.by = "celltype", label = TRUE, repel = TRUE) + theme_classic()
ggsave(file.path(step16b_dir, "UMAP_Celltype_Annotated.png"), plot = p_anno, width = 8, height = 6, dpi = 300)

p_tsne <- DimPlot(mouse_lung_meanvar, reduction = "tsne", group.by = "celltype", label = TRUE, repel = TRUE) + theme_classic()
ggsave(file.path(step16b_dir, "tSNE_Celltype_Annotated.png"), plot = p_tsne, width = 8, height = 6, dpi = 300)

# -------------------- Step 17: Macrophage Reclustering (Dick Manual) --------------------
output_dir_manual <- file.path(results_dir, "Manualrecluster")
dir.create(output_dir_manual, recursive = TRUE, showWarnings = FALSE)

Manual_macro_subset <- subset(mouse_lung_meanvar, idents = c("Alveolar_Macrophages", "Interstitial_Macrophages"))
saveRDS(Manual_macro_subset, file = file.path(output_dir_manual, "Mouselung_Manual_Macrophages.rds"))

# Recalc variable features (mean.var.plot)
Manual_macro_subset <- FindVariableFeatures(
  Manual_macro_subset,
  selection.method = "mean.var.plot",
  mean.cutoff = c(0.02, 4),
  dispersion.cutoff = c(0.3, Inf)
)
message("Manual macrophage variable features: ", length(VariableFeatures(Manual_macro_subset)))

p_elbow_Manual_macro <- ElbowPlot(Manual_macro_subset) + theme_classic()
ggsave(filename = file.path(output_dir_manual, "ElbowPlot_Manual_macro_subset.png"),
       plot = p_elbow_Manual_macro, width = 6, height = 5, dpi = 300)
rm(p_elbow_Manual_macro); gc()

Manual_macro_subset <- ScaleData(Manual_macro_subset)
Manual_macro_subset <- RunPCA(Manual_macro_subset, features = VariableFeatures(Manual_macro_subset))

Manual_macro_subset <- FindNeighbors(Manual_macro_subset, dims = 1:20)
Manual_macro_subset <- FindClusters(Manual_macro_subset, resolution = 0.15)

Idents(Manual_macro_subset) <- Manual_macro_subset$seurat_clusters
levels(Idents(Manual_macro_subset)) <- as.character(as.numeric(levels(Idents(Manual_macro_subset))) + 1)
Manual_macro_subset$MF_subtype <- Idents(Manual_macro_subset)
Idents(Manual_macro_subset) <- Manual_macro_subset$MF_subtype

Manual_macro_subset <- RunUMAP(Manual_macro_subset, dims = 1:20)
Manual_macro_subset <- RunTSNE(Manual_macro_subset, dims = 1:20)

p_mf_umap <- DimPlot(Manual_macro_subset, reduction = "umap", group.by = "MF_subtype", label = TRUE, repel = TRUE) +
  ggtitle("UMAP of MF Subtypes") + theme_classic()
ggsave(file.path(output_dir_manual, "UMAP_MF_Subtype_ManualAnnotated.png"),
       plot = p_mf_umap, width = 8, height = 6, dpi = 300)

p_mf_tsne <- DimPlot(Manual_macro_subset, reduction = "tsne", group.by = "MF_subtype", label = TRUE, repel = TRUE) +
  ggtitle("tSNE of MF Subtypes") + theme_classic()
ggsave(file.path(output_dir_manual, "tSNE_MF_Subtype_ManualAnnotated.png"),
       plot = p_mf_tsne, width = 8, height = 6, dpi = 300)

# -------------------- Marker Detection for MF Subtypes --------------------
Manual_mf_markers <- data.frame()
for (clust in levels(Idents(Manual_macro_subset))) {
  cat("Finding MF markers for cluster", clust, "...\n")
  markers <- FindMarkers(
    Manual_macro_subset,
    ident.1 = clust,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.2
  )
  if (nrow(markers) > 0) {
    markers$cluster <- clust
    markers$gene    <- rownames(markers)
    Manual_mf_markers <- rbind(Manual_mf_markers, markers)
  }
}

Manual_mf_sig <- Manual_mf_markers %>%
  filter(p_val_adj < 0.05) %>%
  filter(!grepl("^mt-", gene, ignore.case = TRUE)) %>%
  filter(!grepl("^Rp[sl]", gene, ignore.case = TRUE))

Manual_top30_mf <- Manual_mf_sig %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(Manual_mf_markers, file.path(output_dir_manual, "DEGs_Manual_MFSubtype_Markers.xlsx"), rowNames = FALSE)
write.xlsx(Manual_top30_mf,  file.path(output_dir_manual, "Top30_Manual_MFSubtype_Markers.xlsx"), rowNames = FALSE)

# -------------------- Violin plots of MF subtype markers --------------------
Manual_vln_dir <- file.path(output_dir_manual, "Violin_Manual_MFSubtype_Markers")
dir.create(Manual_vln_dir, showWarnings = FALSE, recursive = TRUE)

Manual_genes_to_plot <- c(
  # General macrophage markers
  "Csf1r", "Cd68", "Adgre1", "Fcgr1", "Mrc1",
  # MF1 (TLF⁺)
  "Timd4", "Folr2", "Lyve1",
  # MF2 (CCR2⁺)
  "Ccr2", "Plac8", "Ly6c2", "S100a8",
  # MF3 (MHC-IIhi / proliferative)
  "H2-Aa", "H2-Ab1", "Mki67", "Top2a", "Nusap1",
  # Alveolar macrophage-specific (lung-focused)
  "Ear2", "Pparg", "Ear1",
  # Questioned macrophage
  "Cd3d", "Cd3e"
)

for (gene in Manual_genes_to_plot) {
  p_Manual_vln <- VlnPlot(
    Manual_macro_subset,
    features = gene,
    group.by = "MF_subtype",
    pt.size = 0
  ) + ggtitle(gene) + theme_classic()
  ggsave(filename = file.path(Manual_vln_dir, paste0("Violin_", gene, ".png")),
         plot = p_Manual_vln, width = 5, height = 4, dpi = 300)
}

# -------------------- Step 18: Heatmap of MF Subtypes Top30 Marker Genes --------------------
Manual_top30_path <- file.path(output_dir_manual, "Top30_Manual_MFSubtype_Markers.xlsx")
Manual_df_top30   <- read.xlsx(Manual_top30_path)
Manual_df_top30$cluster <- as.character(Manual_df_top30$cluster)

Manual_top30_genes <- unique(Manual_df_top30$gene)
Manual_MF_subtypes <- unique(Manual_df_top30$cluster)

set.seed(111)
Manual_cells_use <- c()
for (stype in Manual_MF_subtypes) {
  cells <- WhichCells(Manual_macro_subset, expression = MF_subtype == stype)
  if (length(cells) > 0) {
    sampled <- sample(cells, size = min(50, length(cells)))
    Manual_cells_use <- c(Manual_cells_use, sampled)
  }
}

Manual_subset_heat <- subset(Manual_macro_subset, cells = Manual_cells_use)
DefaultAssay(Manual_subset_heat) <- "RNA"
Manual_subset_heat <- ScaleData(Manual_subset_heat, features = Manual_top30_genes, verbose = FALSE)

Manual_subset_heat$MF_subtype <- factor(Manual_subset_heat$MF_subtype, levels = Manual_MF_subtypes)
Idents(Manual_subset_heat) <- Manual_subset_heat$MF_subtype

Manual_p_heat <- DoHeatmap(
  object = Manual_subset_heat,
  features = Manual_top30_genes
) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 5),
    strip.text.y = element_text(size = 10, face = "bold"),
    plot.title   = element_text(hjust = 0.5, face = "bold")
  ) +
  ggtitle("Top30 Marker Genes per Manually Defined MF Subtype")

ggsave(
  filename = file.path(output_dir_manual, "Heatmap_Top30_Manual_MFSubtypes_50cells.png"),
  plot = Manual_p_heat,
  width = 13, height = 14, dpi = 600
)

# -------------------- Step 19: Assign MF Subtype Labels (Dick) --------------------
Manual_mf_labels <- c(
  "1" = "AM",
  "2" = "MF1 (TLF+)",
  "3" = "MF2 (CCR2+)",
  "4" = "Proliferative MF",
  "5" = "Questioned MF",
  "6" = "MF2 (CCR2+)"
)

Manual_macro_subset$Manual_MF_subtype_label <- plyr::mapvalues(
  x = Idents(Manual_macro_subset),
  from = names(Manual_mf_labels),
  to   = Manual_mf_labels
)

Idents(Manual_macro_subset) <- Manual_macro_subset$Manual_MF_subtype_label

# -------------------- UMAP with final manual labels --------------------
p_umap_dick_manual <- DimPlot(
  Manual_macro_subset,
  group.by = "Manual_MF_subtype_label",
  label = TRUE,
  repel = TRUE
) + theme_classic()

ggsave(
  filename = file.path(output_dir_manual, "UMAP_MF_Subtype_ManualLabel_dick.png"),
  plot = p_umap_dick_manual,
  width = 8, height = 6, dpi = 300
)

# -------------------- tSNE with final manual labels --------------------
p_tsne_dick_manual <- DimPlot(
  Manual_macro_subset,
  reduction = "tsne",
  group.by = "Manual_MF_subtype_label",
  label = TRUE,
  repel = TRUE
) + theme_classic()

ggsave(
  filename = file.path(output_dir_manual, "tSNE_MF_Subtype_ManualLabel_dick.png"),
  plot = p_tsne_dick_manual,
  width = 8, height = 6, dpi = 300
)

# -------------------- FeaturePlots for classification markers --------------------
key_markers_dick <- c("Ear1", "Pparg", "Folr2", "Ccr2", "Plac8", "Ly6c2", "Mki67", "Top2a", "Nusap1", "Cd3d")

p_feat_umap_dick <- FeaturePlot(
  Manual_macro_subset,
  features = key_markers_dick,
  reduction = "umap",
  order = TRUE,
  cols = c("lightgrey", "red")
) & theme_classic()

ggsave(
  filename = file.path(output_dir_manual, "FeaturePlot_MFSubtypeMarkers_UMAP_dick.png"),
  plot = p_feat_umap_dick,
  width = 12, height = 10, dpi = 300
)

p_feat_tsne_dick <- FeaturePlot(
  Manual_macro_subset,
  features = key_markers_dick,
  reduction = "tsne",
  order = TRUE,
  cols = c("lightgrey", "red")
) & theme_classic()

ggsave(
  filename = file.path(output_dir_manual, "FeaturePlot_MFSubtypeMarkers_tSNE_dick.png"),
  plot = p_feat_tsne_dick,
  width = 12, height = 10, dpi = 300
)

# -------------------- Step 20: Identify Top30 Marker Genes Based on Final MF Subtype Labels --------------------
Idents(Manual_macro_subset) <- Manual_macro_subset$Manual_MF_subtype_label

deg_mf_dick <- FindAllMarkers(
  object = Manual_macro_subset,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top30_mf_dick <- deg_mf_dick %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 30)

write.xlsx(deg_mf_dick,   file.path(output_dir_manual, "DEGs_Reassigned_MFSubtype_All_dick.xlsx"), rowNames = FALSE)
write.xlsx(top30_mf_dick, file.path(output_dir_manual, "Top30_Reassigned_MFSubtype_dick.xlsx"), rowNames = FALSE)

# -------------------- Step 21: Heatmap of Top30 Markers per Final MF Subtype --------------------
top30_path_dick <- file.path(output_dir_manual, "Top30_Reassigned_MFSubtype_dick.xlsx")
df_top30_dick   <- read.xlsx(top30_path_dick)
df_top30_dick$cluster <- as.character(df_top30_dick$cluster)
top30_genes_dick <- unique(df_top30_dick$gene)

set.seed(111)
subtypes <- unique(df_top30_dick$cluster)
cells_use <- c()
for (stype in subtypes) {
  cells <- WhichCells(Manual_macro_subset, expression = Manual_MF_subtype_label == stype)
  if (length(cells) > 0) {
    sampled <- sample(cells, size = min(50, length(cells)))
    cells_use <- c(cells_use, sampled)
  }
}

subset_obj_dick <- subset(Manual_macro_subset, cells = cells_use)
DefaultAssay(subset_obj_dick) <- "RNA"
subset_obj_dick <- ScaleData(subset_obj_dick, features = top30_genes_dick, verbose = FALSE)

subset_obj_dick$Manual_MF_subtype_label <- factor(subset_obj_dick$Manual_MF_subtype_label, levels = subtypes)
Idents(subset_obj_dick) <- subset_obj_dick$Manual_MF_subtype_label

p_heat_dick <- DoHeatmap(
  object = subset_obj_dick,
  features = top30_genes_dick
) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 5),
    strip.text.y = element_text(size = 10, face = "bold")
  )

ggsave(
  filename = file.path(output_dir_manual, "Heatmap_Top30_Reassigned_MFSubtype_50cells_dick.png"),
  plot = p_heat_dick,
  width = 12, height = 12, dpi = 600
)

# ============================ END ================================
