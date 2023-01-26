process PRSice_quant_extract {
    // debug true
    tag "$SNPs"
    label 'vlarge2'
    container 'emosyne/prsice_gwama_exec:1.0'
    cache "lenient"
    errorStrategy 'ignore'

    input:
    // [fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR_QC.bed, fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR_QC.bim, fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR_QC.fam,
        //  EUR_phase3_autosomes_hg19.bed, EUR_phase3_autosomes_hg19.bim, EUR_phase3_autosomes_hg19.fam, 
        // rs4609913, rs4609913.pheno]
    // EUR_phase3_autosomes_hg19.bed, EUR_phase3_autosomes_hg19.bim, EUR_phase3_autosomes_hg19.fam, rs34380086]
    
    tuple path(UKBBbed), path(UKBBbim), path(UKBBfam), path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam), val(SNPs), path(SNP_pheno_file)
    each path (clumped_GWAS_hg19)
    each path (UKBB_covariates)
    each path (SCZ_UKBB_pheno)

    output:
    // tuple  path("*_original_PRS.summary"), path("*_original_PRS.prsice"), path("*_original_PRS.best")
    // tuple path("*.png"), path("*.txt"), path("*.log")//, path("*.all_score") //figures, quantiles text and log
    path ("*")

    
    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    gunzip < ${UKBB_covariates} > covariates.pheno
    head covariates.pheno

    # remove NAs from pheno file, and only keep cases (geno 2 in column 3), then keep only ID and FID cols
    cat ${SNP_pheno_file} | grep -v 'NA' | awk '\$3 ~ /2/' | awk -v OFS='\t' '{print \$1, \$2}' > SNP_pheno.txt
    ## replace zeroes with 1s in caco pheno file
    # cat SNP_pheno.txt | awk 'BEGIN{FS=OFS=" "} {gsub(/0/, "1", \$3)} 1' > SNP_pheno2.txt
    
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --target ${UKBBbed.baseName} \\
        --base ${clumped_GWAS_hg19} \\
        --no-clump \\
        --pheno ${SCZ_UKBB_pheno}\\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P\\
        --stat OR --or \\
        --keep-ambig \\
        --binary-target T --prevalence 0.01 \\
        --quant-extract SNP_pheno.txt --quantile 10 --quant-ref 1\\
        --out ${SNPs}_quantile \\
        --thread $task.cpus \\
        --cov covariates.pheno --cov-factor sex\\
        --memory $mem_mb 
    #   
    #    --ld ${LD_ref_bed.baseName} \\
    #    --clump-kb 3M   \\
    #    --clump-p 1     \\
    #    --clump-r2 0.1  \\
    
    
    """
}


