# ================================================================
# Make SingleR reference (human) from 10X-style matrix + labels
# ================================================================
# Dataset sources (CellxGene, for citation):
# - Heart (Tabula Sapiens – heart subset):
#   https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
# - Lung (An Integrated Cell Atlas of the Human Lung in Health and Disease):
#   https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293
# - Liver (A Single-Cell Atlas of Human Pediatric Liver ...):
#   https://cellxgene.cziscience.com/collections/ff69f0ee-fef6-4895-9f48-6c64a68c8289
# ================================================================

# -------------------- Step 0: Parameters (EDIT HERE) --------------------
organ        <- "lung"  # "lung" | "heart" | "liver"
input_dir    <- "refs/human/lung/HLCA_core"     # folder with matrix.mtx(.gz), features.tsv(.gz) or genes.tsv(.gz), barcodes.tsv(.gz), labels.csv
mtx_file     <- "matrix.mtx.gz"
feat_file    <- "features.tsv.gz"               # or "genes.tsv" / "genes.tsv.gz"
barcode_file <- "barcodes.tsv.gz"
labels_file  <- "labels.csv"
label_main   <- "label.main"                    # column name in labels.csv
label_fine   <- "label.fine"                    # column name in labels.csv
map_to_symbol <- TRUE                           # map Ensembl -> SYMBOL if rownames look like ENSG*
add_logcounts <- TRUE                           # add log1p(counts) as logcounts
out_rds      <- file.path(input_dir, "singler_ref.rds")

# -------------------- Step 1: Libraries --------------------
suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(S4Vectors)
})

# helper: read lines from gz or plain file
.read_lines_any <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  readLines(con)
}
# helper: read table from gz or plain file
.read_table_any <- function(path, header = FALSE, sep = "\t") {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  read.table(con, header = header, sep = sep, quote = "", comment.char = "", stringsAsFactors = FALSE)
}

# -------------------- Step 2: Read counts (Matrix Market) --------------------
mtx_path <- file.path(input_dir, mtx_file)
feat_path <- file.path(input_dir, feat_file)
bar_path <- file.path(input_dir, barcode_file)
lab_path <- file.path(input_dir, labels_file)

stopifnot(file.exists(mtx_path), file.exists(bar_path), file.exists(lab_path))
if (!file.exists(feat_path)) {
  # fallback to genes.tsv(.gz)
  alt <- file.path(input_dir, "genes.tsv")
  alt_gz <- file.path(input_dir, "genes.tsv.gz")
  feat_path <- if (file.exists(alt)) alt else if (file.exists(alt_gz)) alt_gz else stop("features/genes file not found.")
}

expr <- readMM(if (grepl("\\.gz$", mtx_path, TRUE)) gzfile(mtx_path, "rb") else mtx_path)

# -------------------- Step 3: Read features/genes & barcodes --------------------
feat <- .read_table_any(feat_path, header = FALSE, sep = "\t")
# 10X convention: col1 = ID (often Ensembl), col2 = name (symbol), col3 = type
if (ncol(feat) >= 2) {
  gene_id   <- as.character(feat[[1]])
  gene_name <- as.character(feat[[2]])
} else {
  gene_id   <- as.character(feat[[1]])
  gene_name <- gene_id
}
genes <- gene_id  # use IDs first; mapping to SYMBOL happens later if requested
cells <- .read_lines_any(bar_path)

stopifnot(nrow(expr) == length(genes), ncol(expr) == length(cells))

# -------------------- Step 4: Read labels & align with cells --------------------
labels <- read.csv(lab_path, row.names = 1, check.names = FALSE)
stopifnot(all(c(label_main, label_fine) %in% colnames(labels)))

# Align order (use intersection to be safe)
common <- intersect(cells, rownames(labels))
if (length(common) == 0) stop("No overlapping cell barcodes between matrix and labels.")
expr   <- expr[, match(common, cells), drop = FALSE]
labels <- labels[common, , drop = FALSE]
cells  <- common

# -------------------- Step 5: Build SCE --------------------
sce_ref <- SingleCellExperiment(
  assays  = list(counts = expr),
  colData = DataFrame(label.main = labels[[label_main]],
                      label.fine = labels[[label_fine]])
)
rownames(sce_ref) <- genes
colnames(sce_ref) <- cells

# -------------------- Step 6: Optional logcounts --------------------
if (add_logcounts) {
  # simple log1p of counts; for full normalization consider scuttle::logNormCounts
  logcounts(sce_ref) <- log1p(as.matrix(counts(sce_ref)))
}

# -------------------- Step 7: Optional Ensembl -> SYMBOL mapping --------------------
if (map_to_symbol && any(grepl("^ENSG", rownames(sce_ref)))) {
  suppressPackageStartupMessages({
    library(org.Hs.eg.db)
    library(AnnotationDbi)
  })
  ens <- rownames(sce_ref)
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys     = ens,
    column   = "SYMBOL",
    keytype  = "ENSEMBL",
    multiVals = "first"
  )
  keep <- !is.na(symbols) & !duplicated(symbols)
  sce_ref <- sce_ref[keep, ]
  rownames(sce_ref) <- symbols[keep]
}

# -------------------- Step 8: Save --------------------
saveRDS(sce_ref, file = out_rds)
message("Saved SingleR reference: ", normalizePath(out_rds, winslash = "/"))
message("  Cells: ", ncol(sce_ref), " | Genes: ", nrow(sce_ref))
message("  label.main head: ", paste(head(sce_ref$label.main), collapse = ", "))
