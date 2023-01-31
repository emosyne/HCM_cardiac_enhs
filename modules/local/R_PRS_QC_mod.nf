process R_PRS_QC {
    // debug true
    // container 'emosyne/r_docker:1.92'
    container 'emosyne/simpler:1.1'
    label 'process_high'
    // tag "$cohort"
    cache "lenient"
    


    input: 
    // [GWAS_ENH_SNPs_hg19_ALLCHR.bed, GWAS_ENH_SNPs_hg19_ALLCHR.bim, GWAS_ENH_SNPs_hg19_ALLCHR.fam, GWAS_ENH_SNPs_hg19_ALLCHR.prune.in, GWAS_ENH_SNPs_hg19_ALLCHR.het, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/0b/036c27bada52d3b859916ac5896889/GWAS_QC.gz]
    tuple path(bed), path(bim), path(fam), path (prune), path (het), path(HCM_GWAS_QC)
    
    

    output:
    tuple path(bed), path(bim), path(fam), path(HCM_GWAS_QC), path ("*_het_valid_out_vs_HCM_GWAS*.sample"), path("*_a1_cohort_bim_vs_HCM_GWAS*"), path("*_mismatching_SNPs_vs_HCM_GWAS*"),  emit: QC_het_a1_mismatch

    
    script:
    """
    R_PRS_QC2.R ${het} ${bim} ${HCM_GWAS_QC} "UKBB"
    
    """
}
