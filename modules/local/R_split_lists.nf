process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.92'
    label 'process_high'
    tag "${cohort}_${ENH_list}"
    cache "lenient"
    errorStrategy 'ignore'
    

    input:
    // tuple val(cohort), path (LOO_GWAS),  path(cohort_dir), val(ENH_list), path(PGC_noclump_TS_ENH_GWAS_compartment),  path(PGC_noclump_residual_GWAS_compartment), path("*_PGC_clumped_SNPs.clumped"), emit: clumped_SNPs_and_noclump_lists
    tuple val(cohort), path (LOO_GWAS),  val(ENH_list), path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"),       emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${cohort}_${ENH_list}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_TS_ENH_GWAS_compartment} ${EP_ES_gene_brain_exp}
    
    
   
    """
}
    