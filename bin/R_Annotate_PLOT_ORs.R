#!/usr/bin/env Rscript

library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
# library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

#INPUT
args = commandArgs()

print(args)
print("inputs:")
#meta_analysis_results_filelist.txt ${EP_list} ${enhancer_annotations_hg19} ${hg38ToHg19_chain} ${PGC_GWAS_GW_LD_blocks_merge_hg19} PLINK_assoc_filelist.txt 
(meta_analytic_results_filelist = readLines(args[8]))
(EP_list=args[9])
(annotated_EP_list_hg19 = data.table::fread(file=args[10]))
(hg38ToHg19_chain = import.chain(args[11]))
(PGC_GWAS_GW_LD_blocks_merge_hg19 = args[12])
(PLINK_assoc_filelist = readLines(args[13]))


#OUTPUT
annotated_ORs_outfilename = paste0(EP_list, "_", Sys.Date(), "_metaPGC_annotated_ORs.csv")
dir.create(file.path(".", "figs"))
#ORfigure SEE BELOW NAMING
# all_associations_out = paste0(EP_list, "_",Sys.Date(), "_all_EPs_associations.tsv")



##set pval threshold for plotting:
pval=0.001


print("1")
(results_per_snp <- rbind(
  cbind(data.table::fread(meta_analytic_results_filelist[1]), recessive=rep(meta_analytic_results_filelist[1])),
  cbind(data.table::fread(meta_analytic_results_filelist[2]), recessive=rep(meta_analytic_results_filelist[2]))
))

print("2")
#annotate SNPs
SNPcoords<-data.frame(matrix(ncol = 3, nrow = 0))
# all_associations<-data.frame()
colnames(SNPcoords) <- c("seqnames","start","rs_number")
for(filo in PLINK_assoc_filelist){
  print(filo)

  print(coords_file <- data.table::fread(file=filo, select=c("#CHROM","POS","ID"),
                                    col.names = c("seqnames","start","rs_number")))
  SNPcoords=rbind(SNPcoords,coords_file)
  # all_associations=rbind(all_associations, data.table::fread(file=filo, header=T))
}

# all_associations = all_associations %>% group_by(ID) %>% slice_head(n=1) %>% ungroup()
# data.table::fwrite(x=all_associations, file=all_associations_out, sep="\t")

print("3")
(SNPcoords<-SNPcoords %>% group_by(rs_number) %>% slice_head(n=1) %>% ungroup() %>%
    mutate(end=start))

# # filter unique SNVs and then unique E-P_eQTLs by joining BB ORs back to E-P and eQTL SNV data 

names(as_tibble(results_per_snp, .name_repair = "universal"))
# [1] "rs_number"        "reference_allele" "other_allele"     "eaf"             
#  [5] "OR"               "OR_se"            "OR_95L"           "OR_95U"          
#  [9] "z"                "p.value"          "._.log10_p.value" "q_statistic"     
# [13] "q_p.value"        "i2"               "n_studies"        "n_samples"       
# [17] "effects"          "recessive"       
print("4")
(results_per_snp3 <- as_tibble(results_per_snp, .name_repair = "universal") %>% 
  left_join(SNPcoords) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T, na.rm=TRUE))
seqlevelsStyle(results_per_snp3) = "UCSC"
results_per_snp3

print("5")
# #NEED TO MERGE ORS WITH SNV DATA
annotated_EP_list_hg19
table(annotated_EP_list_hg19$tidytype)

print("6")
#if there are overlaps with eQTLs:
if(nrow(annotated_EP_list_hg19[annotated_EP_list_hg19$tidytype=="ESpos_Contact_overlap_eQTLs",])>0){
  print("7a FIRST MERGE THOSE OVERLAPPING GTEX EQTLS")
  (significant_SNVs_withinGTEx_hg38 <- annotated_EP_list_hg19 %>%
    #make granges with GTEx SNV position
    dplyr::select(-start ,-end, start=GTEx_SNV_pos_hg38) %>%
    dplyr::mutate(end=start) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T, na.rm=TRUE))

  seqlevelsStyle(significant_SNVs_withinGTEx_hg38) = "UCSC"
  print("7b")
  (significant_SNVs_withinGTEx_hg19 <- liftOver(significant_SNVs_withinGTEx_hg38, hg38ToHg19_chain) %>%
    unlist() %>%  as_tibble())
  print("7c")
  (significant_SNVs_withinGTEx_hg19 <-
    inner_join(
      as_tibble(results_per_snp3) %>% dplyr::select(-width, -strand) ,
      as_tibble(significant_SNVs_withinGTEx_hg19)  ,#%>% dplyr::select(-width, -strand)
      by = c("seqnames", "start", "end")
    ))
} else {
  print("else 7")
  significant_SNVs_withinGTEx_hg19<-data.frame()
}



print("8")
#THEN annotate THOSE not overlapping eqtls
significant_E_Ps_notGTEx_hg19 <- annotated_EP_list_hg19 %>%
  dplyr::filter() %>%
  #make granges with enh position
  dplyr::select(-GTEx_SNV_pos_hg38) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T, na.rm=TRUE)
seqlevelsStyle(significant_E_Ps_notGTEx_hg19) = "UCSC"
#find overlaps between ENH and ORs
overlaps<-mergeByOverlaps(query = results_per_snp3, subject = significant_E_Ps_notGTEx_hg19)

print("9")
significant_unique_E_Ps_notGTEx_hg19 <-
  cbind(
    as_tibble(overlaps$results_per_snp3) ,
    as_tibble(overlaps$significant_E_Ps_notGTEx_hg19) %>% dplyr::select(-width, -strand)
  ) %>%
  #select cols in other significant dataset to merge
  select(names(significant_SNVs_withinGTEx_hg19))  %>%
  #keep only one per SNP
  group_by(rs_number,enh,ensembl_gene_id,recessive) %>% slice_max(effect_size, n=1, with_ties = F) %>%  ungroup()

print("10")
dim(significant_SNVs_withinGTEx_hg19)
dim(significant_unique_E_Ps_notGTEx_hg19)

print("11")
if(nrow(annotated_EP_list_hg19[annotated_EP_list_hg19$tidytype=="ESpos_Contact_overlap_eQTLs",])>0){
  #### MERGE the 2 lists
  # table(names(significant_SNVs_withinGTEx_hg19)==names(significant_unique_E_Ps_notGTEx_hg19))
  significant_SNVs <-
    rbind (
      significant_SNVs_withinGTEx_hg19,
      #add the SNVs in enhancers, minus the ones already in E_P_eQTLs
      significant_unique_E_Ps_notGTEx_hg19 %>% dplyr::filter(!rs_number %in% significant_SNVs_withinGTEx_hg19$rs_number)
    )
  table(duplicated(significant_SNVs$rs_number))
} else {
  significant_SNVs<-significant_unique_E_Ps_notGTEx_hg19
}
head(significant_SNVs)

print("12")
names(significant_SNVs)
#now deal with duplicated rows due to multiple genes associated with one SNV:
significant_SNVs2<-significant_SNVs %>%
  #select(seqnames,start,snp,enh,ensembl_gene_id, GTEx_variant_id, overlapGTEx_eQTLs, effect_size, OR_2_p, hgnc_symbol) %>%
  #group by SNV vars (all but gene vars)
  group_by(rs_number,recessive) %>%
  #merge all gene ids for each SNV
  mutate(ensembl_gene_id = paste0(ensembl_gene_id, collapse = ","),
         hgnc_symbol = paste0(hgnc_symbol, collapse = ","))  %>%
  #now keep one row per SNV
  slice_max(effect_size, with_ties = F) %>%
  ungroup() %>% arrange(p.value)
head(significant_SNVs2)


# # #what about same E-P, but separate SNPs?
# # results_per_snp2 %>%  group_by(enh, ensembl_gene_id) %>% 
# #   select(seqnames,start,ID,enh,ensembl_gene_id, GTEx_variant_id, overlapGTEx_eQTLs, effect_size, P) %>% 
# #   filter(n()>1) %>% arrange(enh) 
# # (results_per_snp2<- results_per_snp2 %>%  group_by(enh, ensembl_gene_id) %>% 
# #   slice_min(P, with_ties = F) %>% arrange(P) )

print("13")
head(significant_SNVs2[c("rs_number","enh","hgnc_symbol")])
table(significant_SNVs2$tidytype)

# # # #export UCSC track in bed format
# # # names(results_per_snp3)
# # # data.table::fwrite(x =
# # #   results_per_snp3 %>% select ("seqnames", "end", "snp", "OR_2", "OR_2_p" ,"enh","hgnc_symbol") %>% 
# # #     mutate(start=end-1, OR_scaled_p=round(scales::rescale(-log10(OR_2_p), to = c(0, 1000))), CHR=paste0("chr",seqnames), seqnames=NULL, strand=".") %>%
# # #     relocate(CHR, start, end, snp, OR_scaled_p, strand, OR_2, OR_2_p, enh, hgnc_symbol)
# # #   , file = "~/2HH/GWAS_summary_results/brainSNVs_22May23_hg19_UCSC.bed", sep = "\t", row.names = F, col.names = F, scipen = 999
# # # )



# # # import GWAS_LD block overlap data in hg19, by condition
# # ##LD block work from 1000 genomes - self calculated on EUR pop
print("14")
(GWAS_block_merge_hg19 <- data.table::fread(file=PGC_GWAS_GW_LD_blocks_merge_hg19))
GWAS_block_merge_hg19 <- GWAS_block_merge_hg19 %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
seqlevelsStyle(GWAS_block_merge_hg19) = "UCSC"


#annotate GWAS SCZ overlap including significance P VALUE and EXPORT E-P_eQTL list
significant_SNVs2<-makeGRangesFromDataFrame(significant_SNVs2,  keep.extra.columns = T)
seqlevelsStyle(significant_SNVs2) = "UCSC"

significant_SNVs2
GWAS_block_merge_hg19

#for those overlapping a block, annotate block start and end, called start and end, and save SNP location as SNP_position
overlaps<-mergeByOverlaps(query = significant_SNVs2, subject = GWAS_block_merge_hg19)

results_per_snp4 <-
  cbind(
    as_tibble(overlaps$significant_SNVs2),
    as_tibble(overlaps$GWAS_block_merge_hg19) %>%
      mutate( block_start=start, start=NULL, GWAS_SNP_position = SNP_position, SNP_position = NULL,
              block_end=end, end=NULL, seqnames=NULL, width=NULL, strand=NULL,
              GWAS_SNP=SNP, SNP=NULL, GWAS_block_top_OR=OR, OR=NULL, GWAS_block_top_p=P,
              P=NULL)
  )

table(duplicated(results_per_snp4$rs_number))
#add SNPs not matching LD blocks
no_overlaps<-subsetByOverlaps(x = significant_SNVs2, ranges  = GWAS_block_merge_hg19, invert = T) %>% as_tibble() %>%
  mutate(block_start=NA, GWAS_SNP_position=NA, block_end=NA,GWAS_SNP=NA,GWAS_block_top_OR=NA,
         GWAS_block_top_p=NA)
names(results_per_snp4)[names(results_per_snp4) != names(no_overlaps)]
results_per_snp4 <- rbind(results_per_snp4,no_overlaps)
results_per_snp4<-results_per_snp4 %>% arrange(p.value)
table(duplicated(results_per_snp4$rs_number))



# #export SNPs of interest
#remove chr as in base file
results_per_snp4$seqnames<-gsub(pattern = 'chr', replacement = "", x = results_per_snp4$seqnames)
names(results_per_snp4)
table(results_per_snp4$recessive)
results_per_snp4 <- results_per_snp4 %>%
  select(seqnames,start,end,ID=rs_number, 
  recessive,
         reference_allele, other_allele, eaf,
         enh,GTEx_variant_id,ensembl_gene_id,
         OR, OR_95L, OR_95U, P=p.value, n_samples,n_studies,
         #  OR_1,OR_1_lowCI,OR_1_high_CI,OR_1_p,
         #  OR_2, OR_2_lowCI, OR_2_high_CI,OR_2_p,FDR,
         #  OR_A_AA, OR_A_AA_lowCI, OR_A_AA_high_CI, OR_A_AA_p,
         #  OR_AA, OR_AA_lowCI, OR_AA_high_CI,OR_AA_p,
         hgnc_symbol,
         GWAS_SNP_position, block_start, block_end,     GWAS_SNP, GWAS_block_top_OR,GWAS_block_top_p
  ) %>%
  #simplify the recessive name
  mutate(recessive=
    ifelse(test= grepl("ADD",recessive), "ADD", "REC")
  )



table(results_per_snp4$recessive)
head(results_per_snp4)


# data.table::fwrite(x = results_per_snp4, file = "results_per_snp4.csv")


#pivot wider (one row per SNP)
names(results_per_snp4)
(results_per_snp5 <- as_tibble(results_per_snp4) %>% 
  pivot_wider(names_from=recessive, 
              values_from = c(OR:n_studies)))
(results_per_snp5)

data.table::fwrite(x = results_per_snp5, file = annotated_ORs_outfilename)

# seqnames     start       end ID     reference_allele other_allele   eaf enh  
# GTEx_variant_id <chr>,
#   ensembl_gene_id <chr>, hgnc_symbol <chr>, GWAS_SNP_position <int>,
#   block_start <int>, block_end <int>, GWAS_SNP <chr>,
#   GWAS_block_top_OR <dbl>, GWAS_block_top_p <dbl>, OR_ADD <dbl>,
#   OR_REC <dbl>, OR_95L_ADD <dbl>, OR_95L_REC <dbl>, OR_95U_ADD <dbl>,
#   OR_95U_REC <dbl>, P_ADD <dbl>, P_REC <dbl>, n_samples_ADD <int>,
#   n_samples_REC <int>, n_studies_ADD <int>, n_studies_REC <int>


#plot OR figure
# names(as_tibble(meta_analytic_results,    .name_repair = "universal"))

(plot <- as_tibble(results_per_snp5) %>%
   dplyr::filter((P_REC< pval | P_ADD< pval)) %>%
   mutate(SNP=factor(ID), SNP = fct_reorder(SNP, OR_REC, .desc = T))
   )

#reemove inf and NAs
plot=plot[is.finite(rowSums(dplyr::select_if(plot, is.numeric))),]

str(plot)
(max_plot = max(plot$OR_95U_REC, na.rm = TRUE) + 10)

ORfigure= paste0("figs/",EP_list, "_", Sys.Date(),"_p_less_",pval,"_metaPGC_associations.png")
png(filename = ORfigure, width = 4000, height=2800, res = 250)
ggplot(data = plot , aes(x=SNP ) ) + 
  scale_color_gradient(low = "navy blue", high = "white")+
  geom_pointrange(
    aes(y=OR_REC,
        ymin = OR_95L_REC, 
        ymax = OR_95U_REC,
        col=P_REC)) +
    geom_pointrange( 
      color="orange",  alpha=0.25,#fill="white",
      aes(
        y=OR_ADD,
        ymin = OR_95L_ADD, 
        ymax = OR_95U_ADD,         col=P_ADD
      )
    ) +
  ylab("log10(OR)")+ 
  xlab('')+
  scale_y_log10(breaks=seq(0,40,1))+  #limits = c(0.2, max_plot)+
  geom_hline(yintercept =1, linetype=2)+ 
  coord_flip() +
  ggnewscale::new_scale_color() +
  geom_text( 
    aes(label=paste(enh,"gene",hgnc_symbol, "MAF",round(eaf,2), "A1", reference_allele, "n_studies", n_studies_REC), #, "chr", seqnames
        y=0), 
    size=3, hjust = 0) +# scale_colour_manual(values=c("#000000", "#FF5733")) +
  geom_text( aes(label=paste("GWAS SNP: ",GWAS_SNP," GWAS_block_top_OR: ",GWAS_block_top_OR,"| p ",GWAS_block_top_p), 
                  y=max_plot, hjust = "inward"), size=3) + 
  # theme(legend.position = "none") + 
  labs(title =  paste(EP_list,"OR with p<",pval,"in PGC cohort for significant E-Ps in PGC."))
dev.off()
