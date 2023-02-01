process PRSice_calculate_PRS_split_partitions {
    debug true
    tag "${ENH_list}"
    label 'process_high_memory'
    // clusterOptions "--partition=shared_52c_384g" //only for LISA
    container 'emosyne/prsice_gwama_exec:1.0'
    cache "lenient"
    // maxForks 5
    // errorStrategy 'ignore'


    input:
        // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 6k_CARDIAC_NoFibro_significant_noGRB, 
        // 6k_CARDIAC_NoFibro_significant_noGRB_clumped_TS_ENH_GWAS_compartment.tsv.gz, 6k_CARDIAC_NoFibro_significant_noGRB_clumped_residual_GWAS_compartment.tsv.gz, 6k_CARDIAC_NoFibro_significant_noGRB_clumped_merged_GWAS.tsv.gz, 
        // /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/input/biobank/non_missing_10PCs_Jun22.covariate.gz, 
        // /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
        
    tuple path(cohort_bed_QC),  path(cohort_bim_QC), path(cohort_fam_QC), path (HCM_GWAS_QC), val(ENH_list), \
        path(clumped_TS_ENH_GWAS_compartment), path(clumped_residual_GWAS_compartment), path(clumped_merged_GWAS), \
        path(UKBB_covariates), \
        path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)
    
    output:
    tuple val("${ENH_list}"), path("*_clumped_TS_ENH_GWAS_compartment_*.summary"), path("*_clumped_TS_ENH_GWAS_compartment_*.prsice"), path("*_clumped_TS_ENH_GWAS_compartment_*.best"), \
             emit: clumped_TS_ENH_GWAS_compartment_PRS
    tuple val("${ENH_list}"), path("*_clumped_residual_GWAS_compartment.summary"), path("*_clumped_residual_GWAS_compartment.prsice"), path("*_clumped_residual_GWAS_compartment.best"), \
             emit: clumped_residual_GWAS_compartment_PRS
    tuple val("${ENH_list}"), path("*_clumped_merged_GWAS.summary"), path("*_clumped_merged_GWAS.prsice"), path("*_clumped_merged_GWAS.best"),  path(cohort_fam_QC),                     \
             emit: clumped_merged_GWAS_PRS
    tuple val("${ENH_list}"), path("*_original_HCM_GWAS.summary"), path("*_original_HCM_GWAS.prsice"), path("*_original_HCM_GWAS.best"),                                                 \
             emit: clumped_original_HCM_GWAS_PRS
    tuple  path("*.png"), path("*.txt"), path("*.log") //figures, quantiles text and log

    script:
    def mem_Gb = (task.memory * 0.95).toGiga()
    def max_cpus = Math.round(task.cpus * 4/5)
    """
    gunzip < ${UKBB_covariates} > covariates.pheno
    head covariates.pheno
    echo memory: ${mem_Gb}Gb
    echo cpus: $max_cpus

    echo clumped_TS_ENH_GWAS_compartment - ORIGINAL OR
    #CHR    POS     SNP     A1      A2      P       OR      OR_by_measure1  OR_by_measure2  measure1       measure2
    # ORIGINAL OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${ENH_list}_clumped_TS_ENH_GWAS_compartment_originalOR

    echo clumped_TS_ENH_GWAS_compartment  - OR by measure 1
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure1 --or \\
        --out ${ENH_list}_clumped_TS_ENH_GWAS_compartment_OR_by_measure1
    
    echo clumped_TS_ENH_GWAS_compartment  - OR by measure 2
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure2 --or \\
        --out ${ENH_list}_clumped_TS_ENH_GWAS_compartment_OR_by_measure2

    echo clumped_residual_GWAS_compartment - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_residual_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${ENH_list}_clumped_residual_GWAS_compartment
        
    
    echo clumped_merged_GWAS - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_merged_GWAS} \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --out ${ENH_list}_clumped_merged_GWAS
    
    echo ORIGINAL GWAS LOO
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${HCM_GWAS_QC} \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P --stat BETA --beta \\
        --target ${cohort_bed_QC.simpleName} \\
        --ld ${LD_ref_bed.baseName} \\
        --clump-p 1 \\
        --clump-kb 500 --clump-r2 0.1 \\
        --keep-ambig \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.002 \\
        --cov covariates.pheno --cov-factor sex \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --out ${ENH_list}_original_HCM_GWAS
    """

}
        