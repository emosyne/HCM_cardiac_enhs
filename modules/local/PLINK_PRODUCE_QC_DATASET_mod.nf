process PLINK_PRODUCE_QC_DATASET {
    // debug true
    // tag "$cohort"
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"


    input:
    //[/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_het_valid_out_vs_HCM_GWAS.sample, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_a1_cohort_bim_vs_HCM_GWAS, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_mismatching_SNPs_vs_HCM_GWAS, SCZ]
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
