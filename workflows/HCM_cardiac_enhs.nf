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



// #### UKBB input files ####

genotype_chr_files = Channel
    .fromFilePairs("$geno_input_dir/*c*_b0*.{bgen,sample}", flat: true, checkIfExists: true)
    .map{ it-> [it[1],it[2]] }
// genotype_chr_files.view()

UKBBethinicityRelatedness = Channel.fromPath( '../private_input_files/biobank/EIDs_nonBritIrish_includingsecondary_or_related_over_king125.tsv' , checkIfExists: true)
UKBB_covariates =   Channel.fromPath('../private_input_files/biobank/non_missing_10PCs_Jun22.covariate.gz', checkIfExists: true)


//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("$ld_ref_dir/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()


// ## HCM
full_GWAS_hg19 = Channel
    .fromPath("$GWAS_dir/hcm.gwama.sumstats_hg19_24Feb21.gz", checkIfExists: true) 
    // .fromPath("$GWAS_dir/hcm_GWAS_sample.tsv.gz", checkIfExists: true) 
dx_UKBB_pheno =    Channel.fromPath("../private_input_files/biobank/HCM.pheno", checkIfExists: true)
enhancer_lists_bed_files = 
    Channel.from(
        "Non-associated_enh", 
        "Non-cardiac_enh",
        "Cardiac_significant_enh"
        )
        .map { ENH_list -> ["${ENH_list}", 
            file("./input/enh_bedfiles/${ENH_list}.bed", checkIfExists: true)]
            } 
annotations = Channel.fromPath( "../private_input_files/ES_multipliers/2023-02-01_CARDIAC_NoFibro_significant_ES_significant_contact_EPs_ANNOT_plus_100_noOverlap.csv.gz", checkIfExists: true)
condition = "HCM" // SCZ or HCM

// //  SCHIZO and neural lists ##############
// full_GWAS_hg19 = Channel
//     .fromPath("$GWAS_dir/PGC3_SCZ_wave3.european.autosome.public.v3_HCM_format.tsv.gz", checkIfExists: true) 
// dx_UKBB_pheno =    Channel.fromPath("../private_input_files/biobank/SCZ.pheno", checkIfExists: true)
// enhancer_lists_bed_files = 
//     Channel.from(
//         // "18k_PsychENCODE_PFCortex", 
//         "Neural_significant_enh",
//         "Neural_significant_enh_GRB",
//         "Non-neural_enh",
//         "Non-associated_enh")
//             .map { ENH_list -> ["${ENH_list}", 
//                 file("./input/enh_bedfiles/${ENH_list}.bed", checkIfExists: true)]
//             } 
// annotations = Channel.fromPath( "./ ../private_input_files/ES_multipliers/2023-01-18_2023-02-17_NEURAL_ENH_EXP_significant_plus_100_noOverlap_HCMformat.csv.gz", checkIfExists: true)
// condition = "SCZ" // SCZ or HCM


workflow HCM_cardiac_enhs {

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

    chromosomes_by_condition_plus_SNPs = 
        // PGC_GWAS_plus_allEPlists_SNPs_hg19
        R_extract_GWAS_SNPs_into_bed.out.clumped_GWAS_SNPs_plus_those_in_bed_files
            .combine(genotype_chr_files) //The combine operator combines (cartesian product) the items emitted by two channels
            
        
    // chromosomes_by_condition_plus_SNPs.view()

    // GENERATE UKBB UNIQUE FILE
    PLINK2_EXTRACT ( 
        // extract genotypes at bed file locations
        chromosomes_by_condition_plus_SNPs
        )
    
    

    PLINK_MERGE( // SETTING TO BE RESTORED FOR RUNNING IN IMPERIAL
        // merge all bed files into one:
        PLINK2_EXTRACT.out.SNPextracted_by_chromosome.collect(),
        UKBBethinicityRelatedness,
        dx_UKBB_pheno
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

    // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
    // produce QC_het_a1_mismatch, 
    R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
        PLINK2_QC_PRUNE_HET.out.pruned_variants_het
            .combine(PLINK_base_GWAS_QC_and_clump.out.GWAS_QC_noClump)
            .map { [it, condition].flatten() }
    )
    // R_PRS_QC.out.QC_het_a1_mismatch.view()

    // TARGET QC 3:  
    // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
    // also remove mismatching SNPs
    // also standard QC --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 
    PLINK_PRODUCE_QC_DATASET ( //   SETTING TO BE RESTORED FOR RUNNING IN IMPERIAL     --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 \\

        R_PRS_QC.out.QC_het_a1_mismatch
    )

    // PLINK_PRODUCE_QC_DATASET.out.target_QC.view()
    
    
    PLINK_PRODUCE_QC_DATASET.out.target_QC
        .combine(enhancer_lists_bed_files)
        .map { it.flatten() }
        .set{cohort_GWAS_enh_list}
    
    // cohort_GWAS_enh_list.view()

    // BASE subsetting
    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
//    R_prepare_lists_for_clump.out.lists_before_clump
//         .combine(LD_reference)
//         .view()
    
    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference)
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    //     .view()
    

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

    R_final_plot (
        PRS_results
    )

}

