process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high'
    tag "${ENH_list}"
    cache "lenient"
    // errorStrategy 'ignore'
    

    input:
    // // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 18k_PsychENCODE_PFCortex,
        //  SCZ_18k_PsychENCODE_PFCortex_noclump_TS_ENH_GWAS_compartment.tsv.gz, SCZ_18k_PsychENCODE_PFCortex_noclump_residual_GWAS_compartment.tsv.gz, SCZ_18k_PsychENCODE_PFCortex_clumped_SNPs.clumped, SCZ]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs), val(condition)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"), val("1"), val(condition),       emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${ENH_list}_${condition}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_TS_ENH_GWAS_compartment} ${EP_ES_gene_brain_exp} "1"
    
    
   
    """
}
    