include { bash_base_GWAS_QC }           from '../modules/local/bash_base_GWAS_QC.nf'
include { extract_GWAS_SNPs_fromBedFiles }      from '../modules/local/R_extract_GWAS_SNPs_fromBedFiles.nf'
include { PLINK2_EXTRACT }              from '../modules/local/PLINK2_extract_mod.nf'
include { PLINK_MERGE }                 from '../modules/local/PLINK_MERGE_mod.nf'
include { PLINK2_QC_PRUNE_HET }         from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                    from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }    from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { R_prepare_lists_for_clump }   from '../modules/local/R_prepare_lists_for_clump.nf'
include {PLINK_clump}                   from '../modules/local/PLINK_clump_mod.nf'
include {R_split_lists}                 from '../modules/local/R_split_lists.nf'
include {PRSice_calculate_PRS_split_partitions} from '../modules/local/PRSice_calculate_PRS_split_partitions.nf'
include {R_R2_and_logistic_and_quantile_compare}from '../modules/local/R_R2_and_logistic_and_quantile_compare.nf'

// // chain file
// hg38ToHg19_chain = Channel
//     .fromPath( "./input/chainfiles/hg38ToHg19.over.chain", checkIfExists: true)




// #### UKBB input files ####

genotype_chr_files = Channel
    .fromFilePairs("$geno_input_dir/*c*_b0*.{bgen,sample}", flat: true, checkIfExists: true)
    .map{ it-> [it[1],it[2]] }
// genotype_chr_files.view()

UKBBethinicityRelatedness = Channel.fromPath( './input/biobank/EIDs_nonBritIrish_includingsecondary_or_related_over_king125.tsv' , checkIfExists: true)
HCM_UKBB_pheno =    Channel.fromPath("./input/biobank/HCM.pheno", checkIfExists: true)
UKBB_covariates =   Channel.fromPath('./input/biobank/non_missing_10PCs_Jun22.covariate.gz', checkIfExists: true)


//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("$ld_ref_dir/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()


full_GWAS_hg19 = Channel
    .fromPath("$GWAS_dir/hcm.gwama.sumstats_hg19_24Feb21.gz", checkIfExists: true) 
    // .fromPath("$GWAS_dir/hcm_GWAS_sample.tsv.gz", checkIfExists: true) 



enhancer_lists_bed_files = 
    Channel.from(
        "34k_neg", 
        "notCardiac_40k",
        "9k_CARDIAC_NoFibro_significant",
        "6k_CARDIAC_NoFibro_significant_noGRB",
        "3k_CARDIAC_NoFibro_significant_GRB",
        "905_HEART_EP_eQTL"
        )
        .map { ENH_list -> ["${ENH_list}", 
            file("./input/enh_bedfiles/${ENH_list}*.bed", checkIfExists: true)]
            } 
    







workflow HCM {
    // BASE =   SEAN HCM GWAS
    // TARGET = UKBB


    // BASE (GWAS) QC: REMOVE LOW MAF AND INFO SCORES
    //produce GWAS_QC
    bash_base_GWAS_QC (
        full_GWAS_hg19
    )

    extract_GWAS_SNPs_fromBedFiles ( 
        // THIS MODULE IMPORTS 
        // GWAS (hg19), and selects all SNPs in input bed files
        
        enhancer_lists_bed_files.map{it -> it[1]}.collect(),
        bash_base_GWAS_QC.out.GWAS_QC
        )
    
    chromosomes_by_condition_plus_SNPs = 
        // PGC_GWAS_plus_allEPlists_SNPs_hg19
        extract_GWAS_SNPs_fromBedFiles.out.GWAS_SNPs_in_initial_bed_files
            .combine(genotype_chr_files) //The combine operator combines (cartesian product) the items emitted by two channels
            
        
    // chromosomes_by_condition_plus_SNPs.view()

    // GENERATE UKBB UNIQUE FILE
    PLINK2_EXTRACT ( 
        // extract genotypes at bed file locations
        chromosomes_by_condition_plus_SNPs

        //out tuple val(meta), path("*.bim"), path("*.bed"), path ("*.fam"),  emit: SNPextracted_by_chromosome
        )
    
    

    PLINK_MERGE( // SETTING TO BE RESTORED FOR RUNNING IN IMPERIAL
        // merge all bed files into one:
        PLINK2_EXTRACT.out.SNPextracted_by_chromosome.collect(),
        UKBBethinicityRelatedness,
        HCM_UKBB_pheno
        //out tuplepath ("*.bed"), path ("*.bim"), path ("*.fam"),  emit: all_chromosomes_extracted
        )
    // PLINK_MERGE.out.all_chromosomes_extracted.view()

    // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
    // produce prune.in and het files
    PLINK2_QC_PRUNE_HET (
        PLINK_MERGE.out.all_chromosomes_extracted
    )
    
    // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    //         .combine(bash_base_GWAS_QC.out.GWAS_QC)
    //         .view()
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.prune.in, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.het, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/0b/036c27bada52d3b859916ac5896889/GWAS_QC.gz]

    // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
    // produce QC_het_a1_mismatch, 
    R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
        PLINK2_QC_PRUNE_HET.out.pruned_variants_het
            .combine(bash_base_GWAS_QC.out.GWAS_QC)
    )
    // R_PRS_QC.out.QC_het_a1_mismatch.view()
    //[/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/GWAS_QC.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/UKBB_het_valid_out_vs_HCM_GWAS.sample, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/UKBB_a1_cohort_bim_vs_HCM_GWAS, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/7d/9374f235cb5d2bb5e641c8a63511ba/UKBB_mismatching_SNPs_vs_HCM_GWAS]

    // TARGET QC 3:  
    // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
    // also remove mismatching SNPs
    // also standard QC --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 
    PLINK_PRODUCE_QC_DATASET ( //   SETTING TO BE RESTORED FOR RUNNING IN IMPERIAL     --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 \\

        R_PRS_QC.out.QC_het_a1_mismatch
    )

    // PLINK_PRODUCE_QC_DATASET.out.target_QC.view()
    //[GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz]
    
    
    PLINK_PRODUCE_QC_DATASET.out.target_QC
        .combine(enhancer_lists_bed_files)
        .map { it.flatten() }
        .set{cohort_GWAS_enh_list}
    
    // cohort_GWAS_enh_list.view()
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 34k_neg, ./input/enh_bedfiles/34k_neg.bed]
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, notCardiac_40k, ./input/enh_bedfiles/notCardiac_40k.bed]
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 9k_CARDIAC_NoFibro_significant, ./input/enh_bedfiles/9k_CARDIAC_NoFibro_significant.bed]
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 6k_CARDIAC_NoFibro_significant_noGRB, ./input/enh_bedfiles/6k_CARDIAC_NoFibro_significant_noGRB.bed]
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 3k_CARDIAC_NoFibro_significant_GRB, ./input/enh_bedfiles/3k_CARDIAC_NoFibro_significant_GRB.bed]
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 905_HEART_EP_eQTL, ./input/enh_bedfiles/905_HEART_EP_eQTL.bed]
    
    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
    // R_prepare_lists_for_clump.out.lists_before_clump
    //     .combine(LD_reference)
    //     .view()
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 6k_CARDIAC_NoFibro_significant_noGRB, UKBB_6k_CARDIAC_NoFibro_significant_noGRB_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, UKBB_6k_CARDIAC_NoFibro_significant_noGRB_PGC__noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/73/bcc4d22977187ac6f9c0f66d4b15bf/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/73/bcc4d22977187ac6f9c0f66d4b15bf/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/73/bcc4d22977187ac6f9c0f66d4b15bf/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/73/bcc4d22977187ac6f9c0f66d4b15bf/GWAS_QC.gz, 9k_CARDIAC_NoFibro_significant, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/73/bcc4d22977187ac6f9c0f66d4b15bf/UKBB_9k_CARDIAC_NoFibro_significant_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/73/bcc4d22977187ac6f9c0f66d4b15bf/UKBB_9k_CARDIAC_NoFibro_significant_PGC__noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/97/d4decfcb7462d5dcb60006ddcf5979/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/97/d4decfcb7462d5dcb60006ddcf5979/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/97/d4decfcb7462d5dcb60006ddcf5979/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/97/d4decfcb7462d5dcb60006ddcf5979/GWAS_QC.gz, 905_HEART_EP_eQTL, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/97/d4decfcb7462d5dcb60006ddcf5979/UKBB_905_HEART_EP_eQTL_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/97/d4decfcb7462d5dcb60006ddcf5979/UKBB_905_HEART_EP_eQTL_PGC__noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/24/1ccce3653b38bcfbe62a487e0a78e3/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/24/1ccce3653b38bcfbe62a487e0a78e3/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/24/1ccce3653b38bcfbe62a487e0a78e3/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/24/1ccce3653b38bcfbe62a487e0a78e3/GWAS_QC.gz, 3k_CARDIAC_NoFibro_significant_GRB, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/24/1ccce3653b38bcfbe62a487e0a78e3/UKBB_3k_CARDIAC_NoFibro_significant_GRB_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/24/1ccce3653b38bcfbe62a487e0a78e3/UKBB_3k_CARDIAC_NoFibro_significant_GRB_PGC__noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/ac/851a645b8a1e41c651f6c3375a6386/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/ac/851a645b8a1e41c651f6c3375a6386/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/ac/851a645b8a1e41c651f6c3375a6386/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/ac/851a645b8a1e41c651f6c3375a6386/GWAS_QC.gz, 34k_neg, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/ac/851a645b8a1e41c651f6c3375a6386/UKBB_34k_neg_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/ac/851a645b8a1e41c651f6c3375a6386/UKBB_34k_neg_PGC__noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/8a/5a024c9f8e0c6de9ed830356a9d18a/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/8a/5a024c9f8e0c6de9ed830356a9d18a/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/8a/5a024c9f8e0c6de9ed830356a9d18a/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/8a/5a024c9f8e0c6de9ed830356a9d18a/GWAS_QC.gz, notCardiac_40k, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/8a/5a024c9f8e0c6de9ed830356a9d18a/UKBB_notCardiac_40k_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/8a/5a024c9f8e0c6de9ed830356a9d18a/UKBB_notCardiac_40k_PGC__noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
    
    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference)
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    //     .view()
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 6k_CARDIAC_NoFibro_significant_noGRB, UKBB_6k_CARDIAC_NoFibro_significant_noGRB_noclump_TS_ENH_GWAS_compartment.tsv.gz, UKBB_6k_CARDIAC_NoFibro_significant_noGRB_noclump_residual_GWAS_compartment.tsv.gz, clumped_SNPs.clumped]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/GWAS_QC.gz, 9k_CARDIAC_NoFibro_significant, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/UKBB_9k_CARDIAC_NoFibro_significant_noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/UKBB_9k_CARDIAC_NoFibro_significant_noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/42/b17481ea7d933ec5c9cf53ca605004/clumped_SNPs.clumped]
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/GWAS_QC.gz, 3k_CARDIAC_NoFibro_significant_GRB, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/UKBB_3k_CARDIAC_NoFibro_significant_GRB_noclump_TS_ENH_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/UKBB_3k_CARDIAC_NoFibro_significant_GRB_noclump_residual_GWAS_compartment.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/27/86924a3ca9eed1f3dc273d77076c55/clumped_SNPs.clumped]
    
    

    R_split_lists (
        // first annotate SNPs with ES of relevant E-P - for ENH SNP list
        // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
        // ##################################################### CAN MULTIPLY P BY VALUE TO RESTORE ENH SNPS P ###########################################################
        // output separate lists to calculate split PRSs and also merged one
        PLINK_clump.out.clumped_SNPs_and_noclump_lists,
        Channel.fromPath( "./input/ES_multipliers/2023-02-01_CARDIAC_NoFibro_significant_ES_significant_contact_EPs_ANNOT_plus_100_noOverlap.csv.gz", checkIfExists: true)
    )

    
    R_split_lists.out.partitioned //path(bed_QC),  path(bim_QC), path(fam_QC), path (HCM_GWAS_QC), val(ENH_list), path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz")
        .combine(UKBB_covariates)
        .combine(LD_reference)
        .set{combined_splitlists_bedfile_QCeddata_LDdata}
    
    // combined_splitlists_bedfile_QCeddata_LDdata.view()
    // [GWAS_ENH_SNPs_hg19_ALLCHR_QC.bed, GWAS_ENH_SNPs_hg19_ALLCHR_QC.bim, GWAS_ENH_SNPs_hg19_ALLCHR_QC.fam, GWAS_QC.gz, 6k_CARDIAC_NoFibro_significant_noGRB, 6k_CARDIAC_NoFibro_significant_noGRB_clumped_TS_ENH_GWAS_compartment.tsv.gz, 6k_CARDIAC_NoFibro_significant_noGRB_clumped_residual_GWAS_compartment.tsv.gz, 6k_CARDIAC_NoFibro_significant_noGRB_clumped_merged_GWAS.tsv.gz, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/input/biobank/non_missing_10PCs_Jun22.covariate.gz, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bed, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.bim, /Users/eosimo/large_files_not_to_back_up/LD_ref/EUR_phase3_autosomes_hg19.fam]
                
    
    PRSice_calculate_PRS_split_partitions(
        combined_splitlists_bedfile_QCeddata_LDdata
    )
    
    // ########################################### CHANGE NAMES OF MULTIPLIERS ###########################################
    PRS_results = 
        PRSice_calculate_PRS_split_partitions.out.clumped_TS_ENH_GWAS_compartment_PRS
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_merged_GWAS_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_HCM_GWAS_PRS)
            .map { [it, "e_log_OR_X__log_max_ES_perSigEnh__X_10",
                        "e_log_OR_X__log_cardiac_FANTOM_enh_tpm__X_10"].flatten() }


    PRS_results.view()
    
    
    R_R2_and_logistic_and_quantile_compare (
        PRS_results
    )

}


// 
// srun --pty -t 1-00:00 -p shared -c1 nextflow run https://github.com/emosyne/HCM -latest -r master -profile lisa -resume 
// --slurmd-debug=error