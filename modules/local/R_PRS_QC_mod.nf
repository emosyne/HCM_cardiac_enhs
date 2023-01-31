process R_PRS_QC {
    // debug true
    // container 'emosyne/r_docker:1.92'
    container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "$cohort"
    cache "lenient"
    


    input: 
    ////// [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0e/ce0ddf79c4e1e43923154fa08368cc/celso.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0e/ce0ddf79c4e1e43923154fa08368cc/celso.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c4/f0be9b4b466c6063bbb943a07381f8/celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c4/f0be9b4b466c6063bbb943a07381f8/celso_GWAS_QC.gz]
    tuple val(cohort), path (prune), path (het), path(cohort_dir), path(LOO_GWAS_QC)
    
    

    output:
    tuple val(cohort), path(cohort_dir), path ("*_het_valid_out_vs_LOO_GWAS*.sample"), path("*_a1_cohort_bim_vs_LOO_GWAS*"), path("*_mismatching_SNPs_vs_LOO_GWAS*"),  emit: QC_het_a1_mismatch
    // path("*fullGWAS_merged_with_EP_lists_vs_UKBB_annotated_mismatching.csv")    

    
    script:
    """
    bimfile="${cohort_dir}/imputed/hardcall_genotypes/*.bgn.bim"

    R_PRS_QC2.R ${het} \$bimfile ${LOO_GWAS_QC}  ${cohort}
    
    """
}
