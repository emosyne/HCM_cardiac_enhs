process R_prepare_lists_for_clump {
    // debug true
    errorStrategy 'terminate'
    container 'emosyne/r_docker:1.96'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${ENH_list}"
    cache "lenient"
    

    input:
    // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, SCZ_GWAS_QC_nodups.tsv.gz, SCZ, NEURAL_8k_GRB_significant_EPs, ./input/enh_bedfiles/NEURAL_8k_GRB_significant_EPs.bed]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), path (HCM_GWAS_QC),  val(condition), val(ENH_list), path(ENH_bed)

    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path("*_noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_noclump_residual_GWAS_compartment.tsv.gz"),  val(condition), emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${HCM_GWAS_QC}  ${condition}
    
   
    """
}
    
    