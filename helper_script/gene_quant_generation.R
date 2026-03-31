
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(dplyr))

tx2gene <- read.delim("/mnt/storage1/data/Alejandro/Stranded_check/Scripts/tx2gene.gencode.v38.csv")

args <- commandArgs(trailingOnly=T)

#path_salmon <- "/mnt/storage1/data/Alejandro/Salmon_test/quants/E-2025-12513_ISR"
path_salmon <- args[1]

samples <- path_salmon
files <- file.path(samples, "quant.sf")
file.exists(files)

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

#Keep all genes from gene code version 36:
all_genes <- unique(tx2gene$GENEID)
missing_genes <- setdiff(all_genes, rownames(txi$counts))

# Create empty rows
zero_mat <- matrix(0, nrow = length(missing_genes), ncol = ncol(txi$counts))
rownames(zero_mat) <- missing_genes
colnames(zero_mat) <- colnames(txi$counts)

# Add to counts matrix
txi$counts <- rbind(txi$counts, zero_mat)
txi$counts <- txi$counts[sort(rownames(txi$counts)), ]

gene_counts <- as.data.frame(txi$counts)
gene_counts$gene_id <- rownames(gene_counts)
  

#Add gene information gene code version 36 (hg38)
genes_v36 <- read.delim("/mnt/storage1/data/Alejandro/Salmon_test/Scripts/gencode.gene.info.v36.tsv")

gene_counts <- left_join(gene_counts,genes_v36,by="gene_id")
colnames(gene_counts)[1] <- "raw_counts"

tpm_calc = function(x,genes){
  rpk = x/(genes$exon_length/1000)
  scl = sum(rpk)/1e6
  normalized = rpk/scl
  sum(normalized)
  return(normalized)
}

gene_counts <- gene_counts[,c("gene_id","gene_name","seqname","start","end","strand","gene_type","full_length","exon_length","exon_num","raw_counts")]
gene_counts$TPM <- round(tpm_calc(x = gene_counts$raw_counts,genes = gene_counts),3)

write.table(gene_counts,file = paste0(path_salmon,"/quant.gene.sf"),row.names = F,sep = "\t",quote = F)

message("Gene quantification by salmon completed")


