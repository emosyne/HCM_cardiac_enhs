process R_R2_and_logistic_and_quantile_compare {
    // debug true
    container 'emosyne/r_docker:1.94'
    label 'process_mid'
    tag "$cohort_ENHpart"
    cache "lenient"
    errorStrategy 'ignore'

    input: 
    // [clz2a_notNeural_20k_100flank_noInternalOverlap, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_originalOR.summary, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment_originalOR.best, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.summary, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.prsice, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.best, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_merged_GWAS.summary, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_merged_GWAS.prsice, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_merged_GWAS.best, 
    //     clz2a_QC.fam, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_original_HCM_GWAS.summary, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_original_HCM_GWAS.prsice, 
    //     clz2a_notNeural_20k_100flank_noInternalOverlap_original_HCM_GWAS.best, 
    //     exp_log_OR_times__log_max_ES_perEnh_contact_1_3_times_5, 
    //     maxESperEnh_contact_X__brainEnhFantomExp_1_7]
    tuple val(cohort_ENHpart), \
        path(TS_ENH_GWAS_compartment_OR_by_measure1_summary), path(TS_ENH_GWAS_compartment_OR_by_measure2_summary), path(TS_ENH_GWAS_compartment_originalOR_summary), \
        path(TS_ENH_GWAS_compartment_OR_by_measure1_prsice), path(TS_ENH_GWAS_compartment_OR_by_measure2_prsice), path(TS_ENH_GWAS_compartment_originalOR_prsice), \
        path(TS_ENH_GWAS_compartment_OR_by_measure1_best), path(TS_ENH_GWAS_compartment_OR_by_measure2_best), path (TS_ENH_GWAS_compartment_originalOR_best), \
        path(residual_GWAS_compartment_summary), path(residual_GWAS_compartment_prsice), path (residual_GWAS_compartment_best), \
        path(merged_GWAS_summary), path(merged_GWAS_prsice), path (merged_GWAS_best),\
        path(cohort_fam),\
        path(original_HCM_GWAS_summary), path(original_HCM_GWAS_prsice), path (original_HCM_GWAS_best),\
        val(modif_name_1),val(modif_name_2)

    output:
    tuple path("*.txt"), path("*.pdf")

    script:
    """
    
    logistic_and_quantile_compare.R $task.cpus ${cohort_ENHpart} ${cohort_fam} \
        ${TS_ENH_GWAS_compartment_originalOR_summary} ${TS_ENH_GWAS_compartment_originalOR_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure1_summary} ${TS_ENH_GWAS_compartment_OR_by_measure1_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure2_summary} ${TS_ENH_GWAS_compartment_OR_by_measure2_best}\
        ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        ${merged_GWAS_summary} ${merged_GWAS_best}\
        ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
        ${original_HCM_GWAS_summary} ${original_HCM_GWAS_prsice} ${original_HCM_GWAS_best}\
        ${modif_name_1} ${modif_name_2}
    """
}
    