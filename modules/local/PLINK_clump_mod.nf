process PLINK_clump {
    // debug true
    tag "${ENH_list}"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input:
    // [GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 
        // NEURAL_8k_GRB_significant_EPs, SCZ_NEURAL_8k_GRB_significant_EPs_noclump_TS_ENH_GWAS_compartment.tsv.gz, SCZ_NEURAL_8k_GRB_significant_EPs_noclump_residual_GWAS_compartment.tsv.gz, 
        // SCZ, EUR_phase3_autosomes_hg19.bed, EUR_phase3_autosomes_hg19.bim, EUR_phase3_autosomes_hg19.fam]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), \
        val(ENH_list), path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment), \
        val(condition), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)


    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment), path("clumped_SNPs.clumped"), val(condition),    emit: clumped_SNPs_and_noclump_lists
    path("*.log")


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    plink  \\
       --clump ${noclump_TS_ENH_GWAS_compartment},${noclump_residual_GWAS_compartment} \\
       --clump-p1 1 --clump-p2 1 \\
       --clump-kb 500 --clump-r2 0.1 \\
       --bfile ${LD_ref_bed.baseName} \\
       --out clumped_SNPs  \\
       --threads $task.cpus \\
       --memory $mem_mb

    
    """
}
