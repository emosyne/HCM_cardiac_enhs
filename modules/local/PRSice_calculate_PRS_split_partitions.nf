process PRSice_calculate_PRS_split_partitions {
    // debug true
    tag "${cohort}_${ENH_list}"
    label 'process_high_memory'
    clusterOptions "--partition=shared_52c_384g"
    container 'emosyne/prsice_gwama_exec:1.0'
    cache "lenient"
    // maxForks 5
    errorStrategy 'ignore'


    input:
        // [clz2a, daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap,
        //  clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, 
        // clz2a_QC.bed, clz2a_QC.bim, clz2a_QC.fam, 
        // /home/osimoe/PGC_w3_data/clz2a, 
        // /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
        
    tuple val(cohort), path (LOO_GWAS),  val(ENH_list), \
        path(clumped_TS_ENH_GWAS_compartment), path(clumped_residual_GWAS_compartment), path(clumped_merged_GWAS), \
        path(cohort_bed_QC), path(cohort_bim_QC), path(cohort_fam_QC), \
        path(cohort_dir), \
        path(LD_ref_bed), path(LD_ref_bim), path(LD_ref_fam)
    
    output:
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_TS_ENH_GWAS_compartment_*.summary"), path("*_clumped_TS_ENH_GWAS_compartment_*.prsice"), path("*_clumped_TS_ENH_GWAS_compartment_*.best"),          emit: clumped_TS_ENH_GWAS_compartment_PRS
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_residual_GWAS_compartment.summary"), path("*_clumped_residual_GWAS_compartment.prsice"), path("*_clumped_residual_GWAS_compartment.best"),          emit: clumped_residual_GWAS_compartment_PRS
    tuple val("${cohort}_${ENH_list}"), path("*_clumped_merged_GWAS.summary"), path("*_clumped_merged_GWAS.prsice"), path("*_clumped_merged_GWAS.best"),  path(cohort_fam_QC),                              emit: clumped_merged_GWAS_PRS
    tuple val("${cohort}_${ENH_list}"), path("*_original_LOO_GWAS.summary"), path("*_original_LOO_GWAS.prsice"), path("*_original_LOO_GWAS.best"),                                                          emit: clumped_original_LOO_GWAS_PRS
    tuple  path("*.png"), path("*.txt"), path("*.log") //figures, quantiles text and log

    script:
    def mem_Gb = (task.memory * 0.95).toGiga()
    def max_cpus = Math.round(task.cpus * 4/5)
    """
    covariates="${cohort_dir}/prin_comp/*.mds"
    #remove text before star from IID
    cat \$covariates | sed 's/*/1/g' | awk '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13 }' > covariates.pheno
    #< \$covariates cut -d'*' -f2 | awk '{ print \$1,\$2,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$11,\$12,\$13 }' > covariates.pheno
    head covariates.pheno
    echo memory: ${mem_Gb}Gb
    echo cpus: $max_cpus
    
    
    echo clumped_TS_ENH_GWAS_compartment - ORIGINAL OR
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    # ORIGINAL OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${cohort}_${ENH_list}_clumped_TS_ENH_GWAS_compartment_originalOR

    echo clumped_TS_ENH_GWAS_compartment  - OR by measure 1
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure1 --or \\
        --out ${cohort}_${ENH_list}_clumped_TS_ENH_GWAS_compartment_OR_by_measure1
    
    echo clumped_TS_ENH_GWAS_compartment  - OR by measure 2
    #CHR	POS	SNP	A1	A2	P	OR	measure1	measure2	OR_by_measure1	OR_by_measure2
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_TS_ENH_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR_by_measure2 --or \\
        --out ${cohort}_${ENH_list}_clumped_TS_ENH_GWAS_compartment_OR_by_measure2

    echo clumped_residual_GWAS_compartment - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_residual_GWAS_compartment} \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --out ${cohort}_${ENH_list}_clumped_residual_GWAS_compartment
        
    
    echo clumped_merged_GWAS - original OR
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${clumped_merged_GWAS} \\
        --snp SNP --chr CHR --bp POS --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --target ${cohort_bed_QC.simpleName} \\
        --no-clump  \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --out ${cohort}_${ENH_list}_clumped_merged_GWAS
    
    echo ORIGINAL GWAS LOO
    PRSice.R \\
        --prsice /usr/local/bin/PRSice_linux \\
        --base ${LOO_GWAS} \\
        --snp SNP --chr CHR --bp BP --A1 A1 --A2 A2 --pvalue P --stat OR --or \\
        --target ${cohort_bed_QC.simpleName} \\
        --ld ${LD_ref_bed.baseName} \\
        --clump-p 1 \\
        --clump-kb 500 --clump-r2 0.1 \\
        --quantile 10 --quant-ref 1 \\
        --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \\
        --lower 5e-08 \\
        --binary-target T --prevalence 0.01 \\
        --cov covariates.pheno \\
        --thread $max_cpus \\
        --memory ${mem_Gb}Gb \\
        --out ${cohort}_${ENH_list}_original_LOO_GWAS
    """

}
        