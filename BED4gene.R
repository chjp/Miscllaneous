library(biomaRt)
library(GenomicRanges)

exonbed <- function(ensemID){
  ensembl <- useMart("ensembl", dataset="hsapiens_gene_ensembl", host="www.ensembl.org") #GRCh38.p13
  exons <- getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','exon_chrom_start','exon_chrom_end'),
                 filters = 'ensembl_gene_id',
                 values = ensemID,
                 mart = ensembl)
  exon_range = reduce(IRanges(start=exons$exon_chrom_start, end = exons$exon_chrom_end))
  return(data.frame(chr=unique(exons$chromosome_name), start=exon_range@start, end=(exon_range@start+exon_range@width)))
}


BRCA1=exonbed("ENSG00000012048")
BRCA2=exonbed("ENSG00000139618")
GRCh38.p13.BRCA12 = rbind(BRCA1,BRCA2)
write.table(GRCh38.p13.BRCA12, file = "./GRCh38.p13.gene.bed", sep = "\t", row.names = F)
