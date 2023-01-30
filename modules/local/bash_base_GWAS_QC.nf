process bash_base_GWAS_QC {
    tag "$cohort"
    // debug true
    label 'process_low'
    // container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    tuple val(cohort), path(LOO_GWAS), path(cohort_dir)
    

    output:
    tuple val(cohort), path(cohort_dir), path ("*_GWAS_QC.gz"),         emit: GWAS_QC
    // path("*.log")
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    #remove SNPs with INFO < 0.8 and MAF < 0.01
    zcat ${LOO_GWAS} | awk 'NR==1 || (\$6 < 0.99) && (\$8 > 0.8) {print}' | gzip > ${cohort}_GWAS_QC.gz
    """
}

