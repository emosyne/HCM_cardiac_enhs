process PLINK_extract_SNPs {
    // debug true
    label 'process_low'
    // errorStrategy 'ignore' 
    container 'emosyne/plink2:1.23'
    cache "lenient"
    tag "$SNP"

    input: 
    //[fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR_QC.bed, fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR_QC.bim, fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR_QC.fam, EUR_phase3_autosomes_hg19.bed, EUR_phase3_autosomes_hg19.bim, EUR_phase3_autosomes_hg19.fam, rs34380086]
    tuple path(UKBBbed), path(UKBBbim), path(UKBBfam), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam), val(SNP)
    
    output:
    tuple val(SNP), path("*.pheno"),             emit: SNP_pheno
    // path("*.log")

    script:
    """ 
    echo ${SNP} > snp.txt

    plink --bfile ${UKBBbed.simpleName} \\
        --recode AD --extract snp.txt --out ${SNP}

    cat ${SNP}.raw | awk '{ print \$1, \$2, \$7 }' > ${SNP}.pheno

    """
}
