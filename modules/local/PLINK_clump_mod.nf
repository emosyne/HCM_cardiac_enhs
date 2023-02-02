process PLINK_clump {
    // debug true
    tag "${ENH_list}"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input:
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, 34k_neg, UKBB_34k_neg_noclump_TS_ENH_GWAS_compartment.tsv.gz, UKBB_34k_neg_noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)



    output:
    tuple path(bed_QC),  path(bim_QC), path(fam_QC), val(ENH_list), path(noclump_TS_ENH_GWAS_compartment),  path(noclump_residual_GWAS_compartment), path("clumped_SNPs.clumped"), emit: clumped_SNPs_and_noclump_lists
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
