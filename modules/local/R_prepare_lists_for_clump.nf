process R_prepare_lists_for_clump {
    // debug true
    container 'emosyne/r_docker:1.92'
    // container 'emosyne/simpler:1.1'
    label 'process_high'
    tag "${cohort}_${ENH_list}"
    cache "lenient"
    

    input:
    // [celso, daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, /home/osimoe/PGC_w3_data/celso, PsychENCODE_DER_03b_PFC_enhancers_18k, ./input/enh_bedfiles/PsychENCODE_DER_03b_PFC_enhancers_18k.bed]
    tuple val(cohort), path (LOO_GWAS),  path(cohort_dir), val(ENH_list), path(ENH_bed)


    output:
    tuple val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump

    
    script:
    """
    R_prepare_lists_for_clump.R $task.cpus ${ENH_list} ${ENH_bed}  ${LOO_GWAS} ${cohort}
    
   
    """
}
    
    