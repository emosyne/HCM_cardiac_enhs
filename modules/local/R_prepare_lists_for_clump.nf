process R_prepare_lists_for_clump {
    // debug true
    // errorStrategy 'terminate'
    container 'emosyne/r_docker:1.97'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${ENH_list}"
    cache "lenient"
    

    input:
    // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, \
    // SCZ_GWAS_QC_nodups.tsv.gz, SCZ, Neural_significant_enh, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/enh_bedfiles/Neural_significant_enh.bed, 
    // /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/EPWAS/UKBB_ENH_associations_DOM.tsv.gz]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), path (QCed_GWAS),  val(condition), val(ENH_list), path(ENH_bed), path(ENH_EPwas)

    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path("*_noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_noclump_residual_GWAS_compartment.tsv.gz"),  val(condition), emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${QCed_GWAS}  ${condition} ${ENH_EPwas}
    
   
    """
}
    
    