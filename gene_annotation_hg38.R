library(AnnotationHub)
ah <- AnnotationHub()

info <- query(ah, c("Homo.sapiens","Ensembl","GRCh38","gtf"))
gene_anno <- info[['AH51014']]

chr_info <- c(1:22)
gene_info_hg38 <- gene_anno[gene_anno$type == 'gene']
gene_info_hg38_chr <- gene_info_hg38[seqnames(gene_info_hg38) %in% chr_info]

gene_info_hg38_chr_gr <- GRanges(seqnames=paste0("chr",seqnames(gene_info_hg38_chr)),ranges=ranges(gene_info_hg38_chr),strand=strand(gene_info_hg38_chr),gene_name=gene_info_hg38_chr$gene_name)

save(gene_info_hg38_chr_gr,file="./rda/gene_info_hg38_chr_gr.rda")
