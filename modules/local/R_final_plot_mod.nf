process R_final_plot {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high_short'
    tag "$ENH_list"
    cache "lenient"
    // errorStrategy 'ignore'

    input: 
    // [NEURAL_21k_significant_EPs, 
        // SCZ_NEURAL_21k_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, 
        // SCZ_NEURAL_21k_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, 
        // SCZ_NEURAL_21k_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, 
        // SCZ_NEURAL_21k_significant_EPs_0.5_clumped_residual_GWAS_compartment.summary, SCZ_NEURAL_21k_significant_EPs_0.5_clumped_residual_GWAS_compartment.prsice, SCZ_NEURAL_21k_significant_EPs_0.5_clumped_residual_GWAS_compartment.best, 
        // SCZ_NEURAL_21k_significant_EPs_0.5_clumped_merged_GWAS.summary, SCZ_NEURAL_21k_significant_EPs_0.5_clumped_merged_GWAS.prsice, SCZ_NEURAL_21k_significant_EPs_0.5_clumped_merged_GWAS.best, 
        // GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 
        // SCZ_NEURAL_21k_significant_EPs_original_GWAS.summary, SCZ_NEURAL_21k_significant_EPs_original_GWAS.prsice, SCZ_NEURAL_21k_significant_EPs_original_GWAS.best, 
        // 0.5, SCZ, 
        // enh_ES, enh_TS_tpm]

    // // [34k_neg, 
        // SCZ_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, 
        // SCZ_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, 
        // SCZ_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, SCZ_34k_neg_0.5_clumped_residual_GWAS_compartment.summary, SCZ_34k_neg_0.5_clumped_residual_GWAS_compartment.prsice, SCZ_34k_neg_0.5_clumped_residual_GWAS_compartment.best, SCZ_34k_neg_0.5_clumped_merged_GWAS.summary, SCZ_34k_neg_0.5_clumped_merged_GWAS.prsice, SCZ_34k_neg_0.5_clumped_merged_GWAS.best, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, SCZ_34k_neg_original_GWAS.summary, SCZ_34k_neg_original_GWAS.prsice, SCZ_34k_neg_original_GWAS.best, 0.5, SCZ, enh_ES, enh_TS_tpm][18k_PsychENCODE_PFCortex, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_residual_GWAS_compartment.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_residual_GWAS_compartment.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_residual_GWAS_compartment.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_merged_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_merged_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_merged_GWAS.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_original_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_original_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_original_GWAS.best, 0.5, SCZ, enh_ES, enh_TS_tpm]
    tuple val(ENH_list), \
        path(TS_ENH_GWAS_compartment_originalOR_summary), path(TS_ENH_GWAS_compartment_OR_by_measure1_summary), path(TS_ENH_GWAS_compartment_OR_by_measure2_summary),  \
        path(TS_ENH_GWAS_compartment_originalOR_prsice), path(TS_ENH_GWAS_compartment_OR_by_measure1_prsice), path(TS_ENH_GWAS_compartment_OR_by_measure2_prsice),  \
        path (TS_ENH_GWAS_compartment_originalOR_best), path(TS_ENH_GWAS_compartment_OR_by_measure1_best), path(TS_ENH_GWAS_compartment_OR_by_measure2_best),  \
        path(residual_GWAS_compartment_summary), path(residual_GWAS_compartment_prsice), path (residual_GWAS_compartment_best), \
        path(merged_GWAS_summary), path(merged_GWAS_prsice), path (merged_GWAS_best),\
        path(cohort_fam),\
        path(original_GWAS_summary), path(original_GWAS_prsice), path (original_GWAS_best),\
        val(CTthreshold), val(condition),\
        val(modif_name_1),val(modif_name_2)

    output:
    // path("*.txt")
    path("*/*/*.pdf")
    

    script:
    """
    
    R_final_plot.R $task.cpus "${ENH_list}" ${cohort_fam} \
        ${TS_ENH_GWAS_compartment_originalOR_summary} ${TS_ENH_GWAS_compartment_originalOR_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure1_summary} ${TS_ENH_GWAS_compartment_OR_by_measure1_best}\
        ${TS_ENH_GWAS_compartment_OR_by_measure2_summary} ${TS_ENH_GWAS_compartment_OR_by_measure2_best}\
        ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        ${merged_GWAS_summary} ${merged_GWAS_best}\
        ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
        ${original_GWAS_summary} ${original_GWAS_prsice} ${original_GWAS_best}\
        ${modif_name_1} ${modif_name_2} ${CTthreshold} ${condition}
    """
}
    