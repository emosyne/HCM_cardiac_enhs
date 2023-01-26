#!/usr/bin/env Rscript

## PART ONE: IMPORT SNV AND GWAS FILES, MAKE BED FILES
library(tidyverse)
library(data.table)


args = commandArgs()

print(args)
EP_list = args[8]
EP_file_hg19 = args[9]

#outputs
(list_out = paste0(EP_list,"_EP_file_hg19.tsv"))

  
SNPlist = fread(file = EP_file_hg19,
                  #seqnames,start,end,ID,REF,A1,MAF,enh,GTEx_variant_id,ensembl_gene_id,OBS_CT,hgnc_symbol,GWAS_SNP_position,block_start,block_end,GWAS_SNP,GWAS_block_top_OR,GWAS_block_top_p,OR_ADD,OR_REC,L95_ADD,L95_REC,U95_ADD,U95_REC,P_ADD,P_REC
                  select=c("seqnames","ID","start","A1","REF","P_REC","OR_REC"),
                  col.names=c("seqnames","SNP","start", "A1", "A2","P","OR"))
print(SNPlist)

SNPlist %>%
      drop_na() %>%
      fwrite(file=list_out, sep="\t")

