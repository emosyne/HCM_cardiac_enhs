#!/usr/bin/env Rscript

library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

#set max CPU processes
nthreads = as.numeric(args[8])
setDTthreads( round(nthreads))

(ENH_list = args[9])
(ENH_bed_file = args[10])
(HCM_GWAS = args[11])
(cohort = args[12])


#OUTPUT
TS_EPs_outfilename = paste0(cohort, "_", ENH_list, "_noclump_TS_ENH_GWAS_compartment.tsv.gz")
residual_GWAS_compartment_outfilename = paste0(cohort, "_", ENH_list, "_noclump_residual_GWAS_compartment.tsv.gz")



pDivBy = 100000




(ENH_bed = fread(file = ENH_bed_file, select=c(1:4), col.names = c("seqnames","start","end","name")) %>% 
      makeGRangesFromDataFrame(keep.extra.columns = T))
seqlevelsStyle(x = ENH_bed) <- "UCSC"

#extract GWAS SNPs within ranges
print("GWAS head")
# Read in the GWAS data
# SNP     A1      A2      Z       N       FRQ     P       POS     CHR     BETA    SE
(HCM_GWAS <- fread(file=HCM_GWAS, select = c("CHR", "SNP", "POS", "A1", "A2",  "P", "BETA")) %>%
    dplyr::select(CHR,SNP,POS,A1,A2,P,BETA) %>%
    mutate(OR=exp(BETA), BETA=NULL)
    ) 
(HCM_GWAS_hg19 = makeGRangesFromDataFrame(HCM_GWAS, keep.extra.columns = T,
                                               seqnames.field = "CHR", start.field = "POS", 
                                               end.field = "POS"))
seqlevelsStyle(HCM_GWAS_hg19) <- "UCSC"

#subset GWAS by bed
################################################################################################################ CAN DIVIDE P BY VALUE TO PRIORITISE ENH SNPS ########################################################
(full_GWAS_overlap_beds = subsetByOverlaps(x = HCM_GWAS_hg19, ranges = ENH_bed, type="any"))
(full_GWAS_overlap_beds = as_tibble(full_GWAS_overlap_beds) %>%
  dplyr::select(CHR=seqnames, SNP, POS=start, A1, A2, P, OR) %>%
  dplyr::mutate(P=P/pDivBy)
  ) 


fwrite(x= full_GWAS_overlap_beds, file = TS_EPs_outfilename, sep="\t")
#R.utils::gzip(TS_EPs_outfilename,destname=paste0(TS_EPs_outfilename, ".gz"))


#merge annotated_OR_E_Ps and original_base, map to LD blocks, clump, separate 2 lists
#merge annotated_OR_E_Ps and original_base
(full_GWAS_NOoverlap_beds = subsetByOverlaps(x = HCM_GWAS_hg19, ranges = ENH_bed, invert = TRUE))
(full_GWAS_NOoverlap_beds = as_tibble(full_GWAS_NOoverlap_beds) %>%
  dplyr::select(CHR=seqnames, SNP, POS=start, A1, A2, P, OR))

fwrite(x = full_GWAS_NOoverlap_beds, file = residual_GWAS_compartment_outfilename, sep="\t")