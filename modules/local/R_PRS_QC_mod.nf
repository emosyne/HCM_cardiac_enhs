process R_PRS_QC {
    // debug true
    // container 'emosyne/r_docker:1.92'
    container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "$cohort"
    cache "lenient"
    


    input: 
    ////[celso, celso.prune.in, celso.het, daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, /home/osimoe/PGC_w3_data/celso]
    tuple val(cohort), path (prune), path (het), path(LOO_GWAS), path(cohort_dir)
    
    

    output:
    tuple val(cohort), path ("*_het_valid_out_vs_LOO_GWAS*.sample"), path("*_a1_cohort_bim_vs_LOO_GWAS*"), path("*_mismatching_SNPs_vs_LOO_GWAS*"),  emit: QC_het_a1_mismatch
    // path("*fullGWAS_merged_with_EP_lists_vs_UKBB_annotated_mismatching.csv")    

    
    script:
    """
    bimfile="${cohort_dir}/imputed/hardcall_genotypes/*.bgn.bim"

    R_PRS_QC2.R ${het} \$bimfile ${LOO_GWAS}  ${cohort}
    
    """
}
