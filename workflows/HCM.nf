include { PLINK_base_GWAS_QC_and_clump }        from '../modules/local/PLINK_base_GWAS_QC_and_clump.nf'
include { R_extract_GWAS_SNPs_into_bed }        from '../modules/local/R_extract_GWAS_SNPs_into_bed.nf'
include { PLINK2_EXTRACT }                      from '../modules/local/PLINK2_extract_mod.nf'
include { PLINK_MERGE }                         from '../modules/local/PLINK_MERGE_mod.nf'
include { PLINK2_QC_PRUNE_HET }                 from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                            from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }            from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { R_prepare_lists_for_clump }           from '../modules/local/R_prepare_lists_for_clump.nf'
include {PLINK_clump}                           from '../modules/local/PLINK_clump_mod.nf'
include {R_split_lists}                         from '../modules/local/R_split_lists.nf'
include {PRSice_calculate_PRS_split_partitions} from '../modules/local/PRSice_calculate_PRS_split_partitions.nf'
include {R_final_plot}                          from '../modules/local/R_final_plot_mod.nf'

// // chain file
// hg38ToHg19_chain = Channel
//     .fromPath( "./input/chainfiles/hg38ToHg19.over.chain", checkIfExists: true)




// #### UKBB input files ####

genotype_chr_files = Channel
    .fromFilePairs("$geno_input_dir/*c*_b0*.{bgen,sample}", flat: true, checkIfExists: true)
    .map{ it-> [it[1],it[2]] }
// genotype_chr_files.view()

UKBBethinicityRelatedness = Channel.fromPath( './input/biobank/EIDs_nonBritIrish_includingsecondary_or_related_over_king125.tsv' , checkIfExists: true)
UKBB_covariates =   Channel.fromPath('./input/biobank/non_missing_10PCs_Jun22.covariate.gz', checkIfExists: true)


//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("$ld_ref_dir/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()


// full_GWAS_hg19 = Channel
//     .fromPath("$GWAS_dir/hcm.gwama.sumstats_hg19_24Feb21.gz", checkIfExists: true) 
//     // .fromPath("$GWAS_dir/hcm_GWAS_sample.tsv.gz", checkIfExists: true) 
// dx_UKBB_pheno =    Channel.fromPath("./input/biobank/HCM.pheno", checkIfExists: true)
// enhancer_lists_bed_files = 
//     Channel.from(
//         "34k_neg", 
//         "notCardiac_40k",
//         "9k_CARDIAC_NoFibro_significant"//,
//         // "6k_CARDIAC_NoFibro_significant_noGRB",
//         // "3k_CARDIAC_NoFibro_significant_GRB",
//         // "905_HEART_EP_eQTL"
//         )
//         .map { ENH_list -> ["${ENH_list}", 
//             file("./input/enh_bedfiles/${ENH_list}*.bed", checkIfExists: true)]
//             } 
// annotations = Channel.fromPath( "./input/ES_multipliers/2023-02-01_CARDIAC_NoFibro_significant_ES_significant_contact_EPs_ANNOT_plus_100_noOverlap.csv.gz", checkIfExists: true)
// condition = "HCM" // SCZ or HCM

//  SCHIZO and neural lists ##############
full_GWAS_hg19 = Channel
    .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz", checkIfExists: true) 
dx_UKBB_pheno =    Channel.fromPath("./input/biobank/SCZ.pheno", checkIfExists: true)
enhancer_lists_bed_files = 
    Channel.from(
        "18k_PsychENCODE_PFCortex", 
        "NEURAL_21k_significant_EPs",
        "NEURAL_8k_GRB_significant_EPs",
        "20k_notNeural",
        "34k_neg")
            .map { ENH_list -> ["${ENH_list}", 
                file("./input/enh_bedfiles/${ENH_list}*.bed", checkIfExists: true)]
            } 
annotations = Channel.fromPath( "./input/ES_multipliers/2023-01-18_2023-02-17_NEURAL_ENH_EXP_significant_plus_100_noOverlap_HCMformat.csv.gz", checkIfExists: true)
condition = "SCZ" // SCZ or HCM


workflow HCM {
    // BASE =   GWAS
    // TARGET = UKBB


    // BASE (GWAS) QC: REMOVE LOW MAF AND INFO SCORES
    //produce GWAS_QC
    PLINK_base_GWAS_QC_and_clump (
        full_GWAS_hg19
            .combine(LD_reference)
            .map { [it, condition].flatten() }
    )
    

    R_extract_GWAS_SNPs_into_bed ( 
        // THIS MODULE IMPORTS 
        // GWAS (hg19), and selects all SNPs in input bed files and all GWAS clumped SNPs and outputs a bed file
        enhancer_lists_bed_files.map{it -> it[1]}.collect(),
        PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump
            .combine(PLINK_base_GWAS_QC_and_clump.out.clumped_SNPs)
            .map { [it, condition].flatten() }
        
        )
    // R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS_SNPs_plus_those_in_bed_files
    //     .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
    //     .view()
    chromosomes_by_condition_plus_SNPs = 
        // PGC_GWAS_plus_allEPlists_SNPs_hg19
        R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS_SNPs_plus_those_in_bed_files
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
        dx_UKBB_pheno
        //out tuplepath ("*.bed"), path ("*.bim"), path ("*.fam"),  emit: all_chromosomes_extracted
        )
    // PLINK_MERGE.out.all_chromosomes_extracted.view()

    // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
    // produce prune.in and het files
    PLINK2_QC_PRUNE_HET (
        PLINK_MERGE.out.all_chromosomes_extracted
    )
    
    // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    //         .combine(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC)
    //         .view()
    // [/Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.prune.in, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/a1/7a67040dfa3a47d6677a6e9f003f56/GWAS_ENH_SNPs_hg19_ALLCHR.het, /Users/eosimo/GoogleDrive/WORK/CF_PhD/NF_2HH/HCM_cardiac_enhs/work/0b/036c27bada52d3b859916ac5896889/GWAS_QC.gz]

    // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
    // produce QC_het_a1_mismatch, 
    R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
        PLINK2_QC_PRUNE_HET.out.pruned_variants_het
            .combine(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump)
            .map { [it, condition].flatten() }
    )
    // R_PRS_QC.out.QC_het_a1_mismatch.view()
    //[/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/GWAS_ENH_SNPs_hg19_ALLCHR.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_het_valid_out_vs_HCM_GWAS.sample, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_a1_cohort_bim_vs_HCM_GWAS, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/16/7e26f02ebc884d13aa58e154e8c3a7/SCZ_mismatching_SNPs_vs_HCM_GWAS, SCZ]

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
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/SCZ_GWAS_QC_nodups.tsv.gz, SCZ, NEURAL_8k_GRB_significant_EPs, ./input/enh_bedfiles/NEURAL_8k_GRB_significant_EPs.bed]
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/SCZ_GWAS_QC_nodups.tsv.gz, SCZ, 20k_notNeural, ./input/enh_bedfiles/20k_notNeural.bed]
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/3a/ae3a303ab84df51da0f936e223c1d7/SCZ_GWAS_QC_nodups.tsv.gz, SCZ, 34k_neg, ./input/enh_bedfiles/34k_neg.bed]    

    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
//    R_prepare_lists_for_clump.out.lists_before_clump
//         .combine(LD_reference)
//         .view()
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/36/bf081900a69b45cca3a65ece2e46af/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/36/bf081900a69b45cca3a65ece2e46af/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/36/bf081900a69b45cca3a65ece2e46af/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, NEURAL_8k_GRB_significant_EPs, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/36/bf081900a69b45cca3a65ece2e46af/SCZ_NEURAL_8k_GRB_significant_EPs_noclump_TS_ENH_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/36/bf081900a69b45cca3a65ece2e46af/SCZ_NEURAL_8k_GRB_significant_EPs_noclump_residual_GWAS_compartment.tsv.gz, SCZ, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/32/8b1a7419f3343722f9ca4094fdbc39/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/32/8b1a7419f3343722f9ca4094fdbc39/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/32/8b1a7419f3343722f9ca4094fdbc39/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 20k_notNeural, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/32/8b1a7419f3343722f9ca4094fdbc39/SCZ_20k_notNeural_noclump_TS_ENH_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/32/8b1a7419f3343722f9ca4094fdbc39/SCZ_20k_notNeural_noclump_residual_GWAS_compartment.tsv.gz, SCZ, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b2/4a149c4188f7278835fc78b74fe44b/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b2/4a149c4188f7278835fc78b74fe44b/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b2/4a149c4188f7278835fc78b74fe44b/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 34k_neg, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b2/4a149c4188f7278835fc78b74fe44b/SCZ_34k_neg_noclump_TS_ENH_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/b2/4a149c4188f7278835fc78b74fe44b/SCZ_34k_neg_noclump_residual_GWAS_compartment.tsv.gz, SCZ, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam]
    
    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference)
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    //     .view()
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/55/66300e53ff3640d71dbb7b878f322f/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/55/66300e53ff3640d71dbb7b878f322f/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/55/66300e53ff3640d71dbb7b878f322f/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 18k_PsychENCODE_PFCortex, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/55/66300e53ff3640d71dbb7b878f322f/SCZ_18k_PsychENCODE_PFCortex_noclump_TS_ENH_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/55/66300e53ff3640d71dbb7b878f322f/SCZ_18k_PsychENCODE_PFCortex_noclump_residual_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/55/66300e53ff3640d71dbb7b878f322f/SCZ_18k_PsychENCODE_PFCortex_clumped_SNPs.clumped, SCZ]
    

    R_split_lists (
        // first annotate SNPs with ES of relevant E-P - for ENH SNP list
        // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
        // ##################################################### CAN MULTIPLY P BY VALUE TO RESTORE ENH SNPS P ###########################################################
        // output separate lists to calculate split PRSs and also merged one
        PLINK_clump.out.clumped_SNPs_and_noclump_lists,//.map { [it, "1"].flatten() }, //######################## multiplier can be set here ########################
        annotations
    )

    
    R_split_lists.out.partitioned 
        .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
        .combine(UKBB_covariates)
        .combine(LD_reference)
        .map { [it, "0.5"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
        .set{combined_splitlists_bedfile_QCeddata_LDdata_05}
    R_split_lists.out.partitioned 
        .combine(R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS)
        .combine(UKBB_covariates)
        .combine(LD_reference)
        .map { [it, "0.05"].flatten() }         // ######################## SET CT THRESHOLD FOR PRSICE ##################
        .set{combined_splitlists_bedfile_QCeddata_LDdata_005}
    
    combined_splitlists_bedfile_QCeddata_LDdata = combined_splitlists_bedfile_QCeddata_LDdata_05.mix(combined_splitlists_bedfile_QCeddata_LDdata_005)
    // combined_splitlists_bedfile_QCeddata_LDdata.view()
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/52/86efea6cc0d5914f82f009c7775fdd/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/52/86efea6cc0d5914f82f009c7775fdd/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/52/86efea6cc0d5914f82f009c7775fdd/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, 20k_notNeural, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/52/86efea6cc0d5914f82f009c7775fdd/20k_notNeural_SCZ_X_1_clumped_TS_ENH_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/52/86efea6cc0d5914f82f009c7775fdd/20k_notNeural_SCZ_clumped_residual_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/52/86efea6cc0d5914f82f009c7775fdd/20k_notNeural_SCZ_clumped_merged_GWAS.tsv.gz, 1, SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8e/451eaa3242048fc07fb0496d3717cd/SCZ_clumped_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/biobank/non_missing_10PCs_Jun22.covariate.gz, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam, 0.5]
    // [/rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/9c/32cf754c26669da6278362f0e09542/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bed, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/9c/32cf754c26669da6278362f0e09542/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.bim, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/9c/32cf754c26669da6278362f0e09542/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, NEURAL_8k_GRB_significant_EPs, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/9c/32cf754c26669da6278362f0e09542/NEURAL_8k_GRB_significant_EPs_SCZ_X_1_clumped_TS_ENH_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/9c/32cf754c26669da6278362f0e09542/NEURAL_8k_GRB_significant_EPs_SCZ_clumped_residual_GWAS_compartment.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/9c/32cf754c26669da6278362f0e09542/NEURAL_8k_GRB_significant_EPs_SCZ_clumped_merged_GWAS.tsv.gz, 1, SCZ, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/8e/451eaa3242048fc07fb0496d3717cd/SCZ_clumped_GWAS_QC_nodups.tsv.gz, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/input/biobank/non_missing_10PCs_Jun22.covariate.gz, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bed, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.bim, /rds/general/user/eosimo/home/lenhard_prs/LD_ref/EUR_phase3_autosomes_hg19.fam, 0.5]
                
    
    PRSice_calculate_PRS_split_partitions(
        combined_splitlists_bedfile_QCeddata_LDdata
    )
    
    // ########################################### SET NAMES OF MULTIPLIERS ###########################################
    PRS_results = 
        PRSice_calculate_PRS_split_partitions.out.clumped_TS_ENH_GWAS_compartment_PRS
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_merged_GWAS_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_GWAS_PRS)
            .map { [it, "enh_ES", "enh_TS_tpm"].flatten() }


    // PRS_results.view()
        // [NEURAL_21k_significant_EPs, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_residual_GWAS_compartment.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_residual_GWAS_compartment.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_residual_GWAS_compartment.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_merged_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_merged_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_0.5_clumped_merged_GWAS.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_original_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_original_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/60/34063b61da3429a69a573b05acb9c1/SCZ_NEURAL_21k_significant_EPs_original_GWAS.best, 0.5, SCZ, enh_ES, enh_TS_tpm]
        // [34k_neg, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_residual_GWAS_compartment.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_residual_GWAS_compartment.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_residual_GWAS_compartment.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_merged_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_merged_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_0.5_clumped_merged_GWAS.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_original_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_original_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/e4/2c51edd868cb408dd99f2ba96ea03a/SCZ_34k_neg_original_GWAS.best, 0.5, SCZ, enh_ES, enh_TS_tpm][18k_PsychENCODE_PFCortex, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_residual_GWAS_compartment.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_residual_GWAS_compartment.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_residual_GWAS_compartment.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_merged_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_merged_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_0.5_clumped_merged_GWAS.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_original_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_original_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/17/27cee8b916da1b4cb4c56cecad973b/SCZ_18k_PsychENCODE_PFCortex_original_GWAS.best, 0.5, SCZ, enh_ES, enh_TS_tpm]
        // [20k_notNeural, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_TS_ENH_GWAS_compartment_originalOR.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_mult_1_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_residual_GWAS_compartment.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_residual_GWAS_compartment.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_residual_GWAS_compartment.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_merged_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_merged_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_0.5_clumped_merged_GWAS.best, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/GWAS_ENH_SNPs_hg19_ALLCHR_SCZ_QC.fam, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_original_GWAS.summary, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_original_GWAS.prsice, /rds/general/ephemeral/user/eosimo/ephemeral/HCM_cardiac_enhs/work/90/dd3c30856f9c78a17fd5dc0fe552d7/SCZ_20k_notNeural_original_GWAS.best, 0.5, SCZ, enh_ES, enh_TS_tpm]  

    R_final_plot (
        PRS_results
    )

}


// 
// srun --pty -t 1-00:00 -p shared -c1 nextflow run https://github.com/emosyne/HCM -latest -r master -profile lisa -resume 
// --slurmd-debug=error