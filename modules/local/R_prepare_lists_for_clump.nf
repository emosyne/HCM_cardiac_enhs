process R_prepare_lists_for_clump {
    debug true
    errorStrategy 'terminate'
    container 'emosyne/r_docker:1.94'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${ENH_list}"
    cache "lenient"
    

    input:
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 34k_neg, ./input/enh_bedfiles/34k_neg.bed]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), path (HCM_GWAS_QC), val(ENH_list), path(ENH_bed)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), path (HCM_GWAS_QC), val(ENH_list), path(ENH_bed), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${HCM_GWAS_QC}  "UKBB"
    
   
    """
}
    
    