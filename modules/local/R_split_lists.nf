process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.94'
    label 'process_high'
    tag "${ENH_list}"
    cache "lenient"
    // errorStrategy 'ignore'
    

    input:
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, 6k_CARDIAC_NoFibro_significant_noGRB, 
    // UKBB_6k_CARDIAC_NoFibro_significant_noGRB_noclump_TS_ENH_GWAS_compartment.tsv.gz, UKBB_6k_CARDIAC_NoFibro_significant_noGRB_noclump_residual_GWAS_compartment.tsv.gz, clumped_SNPs.clumped]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs), val(multiplier)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"), val(multiplier),       emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${ENH_list}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_TS_ENH_GWAS_compartment} ${EP_ES_gene_brain_exp} ${multiplier}
    
    
   
    """
}
    