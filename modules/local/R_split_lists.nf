process R_split_lists {
    // debug true
    container 'emosyne/r_docker:1.97'
    label 'process_high'
    tag "${ENH_list}"
    cache "lenient"
    // errorStrategy 'ignore'
    

    input:
    //     tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path(noclump_EPWAS),  path(noclump_residual_GWAS_compartment), path("*clumped_SNPs.clumped"), \
            // val(condition), val(EPWAS_model),   emit: clumped_SNPs_and_noclump_lists
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path(noclump_EPWAS),  path(noclump_residual_GWAS_compartment),  path(clumped_SNPs), val(condition)
    each path(EP_ES_gene_brain_exp)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), \
        path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"), val("1"), val(condition),       emit: partitioned
    
    
    
    script:
    """
    R_split_lists.R "${ENH_list}_${condition}" ${clumped_SNPs} ${noclump_residual_GWAS_compartment} ${noclump_EPWAS} ${EP_ES_gene_brain_exp} "1"
    
    
   
    """
}
    