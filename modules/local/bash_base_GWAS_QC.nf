process bash_base_GWAS_QC {
    // tag "$cohort"
    // debug true
    label 'process_low'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    path(HCM_GWAS)
    

    output:
    path ("GWAS_QC_nodups.tsv.gz"),         emit: GWAS_QC
    // path("*.log")
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    #remove SNPs with MAF < 0.01 (in this case main allele freq <0.99) and remove duplicated SNPs
    zcat ${HCM_GWAS} | awk 'NR==1 || (\$6 < 0.99) {print}' | awk '!seen[\$1]++' | sed -E 's/,/\t/g' | gzip > GWAS_QC_nodups.tsv.gz

    """
}

