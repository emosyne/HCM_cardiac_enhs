process PLINK_PRODUCE_QC_DATASET {
    // debug true
    // tag "$cohort"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"


    input:
    //[GWAS_ENH_SNPs_hg19_ALLCHR.bed, GWAS_ENH_SNPs_hg19_ALLCHR.bim, GWAS_ENH_SNPs_hg19_ALLCHR.fam, GWAS_QC.gz, UKBB_het_valid_out_vs_HCM_GWAS.sample, UKBB_a1_cohort_bim_vs_HCM_GWAS, UKBB_mismatching_SNPs_vs_HCM_GWAS]
    tuple path(bed), path(bim), path(fam), path(HCM_GWAS_QC), path (het_valid), path (a1_bim), path(mismatch), val(condition)

    output:
    tuple path ("*_QC.bed"), path ("*_QC.bim"), path ("*_QC.fam"), path(HCM_GWAS_QC), val(condition),           emit: target_QC
    path ("*.log")


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    plink \\
        --bfile ${bed.simpleName} \\
        --make-bed \\
        --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \\
        --a1-allele ${a1_bim} \\
        --keep ${het_valid} \\
        --exclude ${mismatch} \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out ${bed.simpleName}_${condition}_QC
    

    """
}
//        --mind 0.1
