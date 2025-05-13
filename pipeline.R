# -----------------------------------------------
# Obese vs. Lean Adipose RNA-Seq  (GSE180016)
# complete DESeq2 pipeline + volcano + enrichment
# -----------------------------------------------

library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)
library(gprofiler2)

# ----- 1. download ready-made counts matrix -----
url_main <- "https://geo-downloads.s3.amazonaws.com/GSE180016_raw_counts.txt.gz"
dest     <- "GSE180016_counts.txt.gz"

if (!file.exists(dest)) {
  message("Downloading counts matrix ...")
  download.file(url_main, destfile = dest, mode = "wb", quiet = FALSE)
}

# ----- 2. load counts -----
counts <- read_delim(dest, delim = "\t")
gene_ids <- counts$Gene; counts$Gene <- NULL
rownames(counts) <- gene_ids

# ----- 3. build metadata (first 6 = lean, last 6 = obese) -----
sample_names <- colnames(counts)
meta <- tibble(
  sample    = sample_names,
  condition = c(rep("lean", 6), rep("obese", 6))
) %>% column_to_rownames("sample")

# ----- 4. DESeq2 -----
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = meta,
                              design    = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]   # keep expressed genes
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "obese", "lean"))
resOrdered <- res[order(res$padj), ]
write_csv(as.data.frame(resOrdered), "DESeq2_results.csv")

# ----- 5. volcano plot -----
png("volcano.png", width = 1800, height = 1200, res = 200)
EnhancedVolcano(resOrdered,
                lab      = rownames(resOrdered),
                x        = "log2FoldChange",
                y        = "padj",
                FCcutoff = 1,
                pCutoff  = 0.05,
                title    = "Obese vs Lean adipose (GSE180016)")
dev.off()

# ----- 6. pathway enrichment (up-regulated genes) -----
sig_up <- rownames(resOrdered[resOrdered$padj < 0.05 &
                                resOrdered$log2FoldChange > 1 , ])
gp <- gost(sig_up, organism = "hsapiens", significant = TRUE)
write_csv(gp$result, "enrichment_upregulated.csv")

message("✓ Pipeline complete — results saved:")
message("  • DESeq2_results.csv")
message("  • volcano.png")
message("  • enrichment_upregulated.csv")
