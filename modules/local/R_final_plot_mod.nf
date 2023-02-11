process R_final_plot {
    // debug true
    container 'emosyne/r_docker:1.94'
    label 'process_high'
    tag "$ENH_list"
    cache "lenient"
    // errorStrategy 'ignore'

    input: 
    // [905_HEART_EP_eQTL, 
    // 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_originalOR.summary, 
    // 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_originalOR.prsice,
    // 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, 905_HEART_EP_eQTL_clumped_TS_ENH_GWAS_compartment_originalOR.best, 
    // 905_HEART_EP_eQTL_clumped_residual_GWAS_compartment.summary, 905_HEART_EP_eQTL_clumped_residual_GWAS_compartment.prsice, 905_HEART_EP_eQTL_clumped_residual_GWAS_compartment.best, 
    // 905_HEART_EP_eQTL_clumped_merged_GWAS.summary, 905_HEART_EP_eQTL_clumped_merged_GWAS.prsice, 905_HEART_EP_eQTL_clumped_merged_GWAS.best, 
    // GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, 
    // 905_HEART_EP_eQTL_original_HCM_GWAS.summary, 905_HEART_EP_eQTL_original_HCM_GWAS.prsice, 905_HEART_EP_eQTL_original_HCM_GWAS.best,
    // e_log_OR_X__log_max_ES_perSigEnh__X_10, e_log_OR_X__log_cardiac_FANTOM_enh_tpm__X_10]
    tuple val(ENH_list), \
        path(TS_ENH_GWAS_compartment_OR_by_measure1_summary), path(TS_ENH_GWAS_compartment_OR_by_measure2_summary), path(TS_ENH_GWAS_compartment_originalOR_summary), \
        path(TS_ENH_GWAS_compartment_OR_by_measure1_prsice), path(TS_ENH_GWAS_compartment_OR_by_measure2_prsice), path(TS_ENH_GWAS_compartment_originalOR_prsice), \
        path(TS_ENH_GWAS_compartment_OR_by_measure1_best), path(TS_ENH_GWAS_compartment_OR_by_measure2_best), path (TS_ENH_GWAS_compartment_originalOR_best), \
        path(residual_GWAS_compartment_summary), path(residual_GWAS_compartment_prsice), path (residual_GWAS_compartment_best), \
        path(merged_GWAS_summary), path(merged_GWAS_prsice), path (merged_GWAS_best),\
        path(cohort_fam),\
        path(original_HCM_GWAS_summary), path(original_HCM_GWAS_prsice), path (original_HCM_GWAS_best),\
        val(CTthreshold),\
        val(modif_name_1),val(modif_name_2)

    output:
    // path("*.txt")
    path("*.pdf")

    script:
    """
    
    R_final_plot.R $task.cpus ${ENH_list} ${cohort_fam} \
        ${TS_ENH_GWAS_compartment_originalOR_summary} ${TS_ENH_GWAS_compartment_originalOR_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure1_summary} ${TS_ENH_GWAS_compartment_OR_by_measure1_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure2_summary} ${TS_ENH_GWAS_compartment_OR_by_measure2_best}\
        ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        ${merged_GWAS_summary} ${merged_GWAS_best}\
        ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
        ${original_HCM_GWAS_summary} ${original_HCM_GWAS_prsice} ${original_HCM_GWAS_best}\
        ${modif_name_1} ${modif_name_2} ${CTthreshold}
    """
}
    