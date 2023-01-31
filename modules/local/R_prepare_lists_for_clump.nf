process R_prepare_lists_for_clump {
    // debug true
    container 'emosyne/r_docker:1.92'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${cohort}_${ENH_list}"
    cache "lenient"
    

    input:
    // [xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234_GWAS_QC.gz, NEURAL_8k_GRB_significant_EPs, ./input/enh_bedfiles/NEURAL_8k_GRB_significant_EPs.bed]
    tuple val(cohort),  path(cohort_dir), path (LOO_GWAS_QC), val(ENH_list), path(ENH_bed)


    output:
    tuple val(cohort), path (LOO_GWAS_QC),  val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${LOO_GWAS_QC} ${cohort}
    
   
    """
}
    
    