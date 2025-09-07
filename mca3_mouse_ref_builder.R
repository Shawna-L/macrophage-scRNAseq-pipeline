# ================================================================
# MCA 3.0 → SingleR reference builder (mouse, all organs)
# ================================================================
# Author: Jingshuo Li
# Version: 1.0
# Description: Build SCE reference for SingleR from MCA 3.0 mouse data.
# Data source: Mouse Cell Atlas 3.0 – https://bis.zju.edu.cn/MCA/atlas3.html
# ================================================================

# -------------------- Step 0: Parameters (EDIT HERE) --------------------
organ      <- "lung"                                  # "lung" | "heart" | "liver" | "kidney" | "brain"
organ_tag  <- "Adult_Lung"                            # change to "Adult_Heart"/"Adult_Liver"/...
organ_key  <- "AdultLung"                             # metadata filter key (no underscore)
mca_dir    <- "data/MCA3/mouse"                       # base directory for MCA organ folders
meta_path  <- "data/MCA3/MCA3.0_data_annotation.xlsx" # MCA 3.0 global annotation
refs_dir   <- "refs/mouse"                            # output directory for SingleR references

expr_path  <- file.path(mca_dir, organ_tag, paste0(organ_tag, "_dge.csv"))
gene_path  <- file.path(mca_dir, organ_tag, paste0(organ_tag, "_gene.csv"))
meta_out   <- file.path(mca_dir, organ_tag, paste0(organ_tag, "_Meta.csv"))
ref_out    <- file.path(refs_dir, paste0("mca_", organ, "_ref.rds"))

dir.create(refs_dir, showWarnings = FALSE, recursive = TRUE)

# -------------------- Step 1: Libraries --------------------
suppressPackageStartupMessages({
  library(readxl)
  library(SingleCellExperiment)
})

# -------------------- Step 2: Load expression matrix --------------------
# Expression matrix: genes x cells
expr <- read.csv(expr_path, header = TRUE, check.names = FALSE)

# Load gene names and assign to rownames
genes <- read.csv(gene_path, header = FALSE, stringsAsFactors = FALSE)
rownames(expr) <- genes[[1]]
expr <- as.matrix(expr)

# -------------------- Step 3: Load annotation metadata --------------------
meta <- read_excel(meta_path)

# Filter to target organ and keep barcodes present in expr
organ_meta <- meta[grepl(organ_key, meta$cell_name), ]
organ_meta <- organ_meta[organ_meta$cell_name %in% colnames(expr), ]
organ_meta <- organ_meta[match(colnames(expr), organ_meta$cell_name), ]

# Save filtered metadata snapshot
write.csv(organ_meta, file = meta_out, row.names = FALSE)

# -------------------- Step 4: Create SingleCellExperiment object --------------------
ref_sce <- SingleCellExperiment(
  assays  = list(counts = expr),
  colData = DataFrame(label.main = organ_meta$cell_type)
)

# -------------------- Step 5: Save as RDS for SingleR --------------------
saveRDS(ref_sce, file = ref_out)
