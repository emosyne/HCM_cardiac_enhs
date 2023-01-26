include { PLINK2_ASSOC_GLM }            from '../modules/local/PLINK2_ASSOC_GLM_mod.nf'
include { PLINK2_QC_PRUNE_HET }         from '../modules/local/PLINK2_QC_PRUNE_HET_mod.nf'
include { R_PRS_QC }                    from '../modules/local/R_PRS_QC_mod.nf'
include { PLINK_PRODUCE_QC_DATASET }    from '../modules/local/PLINK_PRODUCE_QC_DATASET_mod.nf'
include { R_forest_by_SNP_by_sample }   from '../modules/local/R_forest_by_SNP_by_sample.nf'
include { R_prepare_lists_for_clump }   from '../modules/local/R_prepare_lists_for_clump.nf'
include {PLINK_clump}                   from '../modules/local/PLINK_clump_mod.nf'
include {R_split_lists}                 from '../modules/local/R_split_lists.nf'
include {PRSice_calculate_PRS_split_partitions} from '../modules/local/PRSice_calculate_PRS_split_partitions.nf'
include {R_R2_and_logistic_and_quantile_compare}from '../modules/local/R_R2_and_logistic_and_quantile_compare.nf'

// // chain file
// hg38ToHg19_chain = Channel
//     .fromPath( "./input/chainfiles/hg38ToHg19.over.chain", checkIfExists: true)
// // SCZ wave 3 full GWAS merged to LD blocks 1000 genomes
// PGC_GWAS_GW_LD_blocks_merge_hg19 = Channel
//     .fromPath( "./input/LD/PGC_Jul22_GWAS_LD_block.tsv.gz", checkIfExists: true)

// // EP annotations
// enhancer_annotations_hg19 = Channel
//     .fromPath("./input/initial_SNP_lists/annotated_EPs.csv.gz", checkIfExists: true)
// // SCZ wave 3 full GWAS
// full_GWAS_hg19 = Channel
//     .fromPath("http://data.genereg.net/emanuele/GWAS_PGC_summary/fullGWAS_SCZ_PGC3_SCZ_wave3.european.autosome.public.v3_tidy_hg19.tsv.gz", checkIfExists: true)





// INPUTS
test_sample_list = ["celso","clz2a"]//"xirwt","xgras","xjrsa","gawli","xmgs2","mcqul","xclm2","xclo3","gpc2a","xboco","xs234","xswe5","xswe6",
LOO_PGC_GWASES = 
    Channel.fromList( test_sample_list )
    .map{[it, file("/home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.no${it}.gz")]}
   

// LOO_PGC_GWASES.view()


no_PCA_sample_list = ["butr","grtr","lemu","uktr"]
Channel.fromPath( '/home/osimoe/PGC_w3_data/*', type: 'dir' )
    .map{ it -> [it.simpleName , it] }
    .branch { sample_name, sample_path ->
        no_PCA: sample_name in no_PCA_sample_list
            return tuple( tuple( sample_name, sample_path ) )
        with_PCA: true//sample_name in test_sample_list
            return tuple( tuple( sample_name, sample_path ) )
    } \
    .set { inputs }

// sample_list_over = ["celso","xirwt","xgras","xjrsa","gawli","xmgs2","mcqul","xclm2","xclo3","gpc2a","xboco","xs234","xswe5","xswe6","clz2a"]
// inputs.with_PCA
//     .branch { sample_name, sample_path ->
//         overthousand: sample_name in sample_list_over
//             return tuple( tuple( sample_name, sample_path ) )
//         small: true
//             return tuple( tuple( sample_name, sample_path ) )
//     } \
//     .set { branchedinputs }

validation_samples = inputs.with_PCA



enhancer_lists_bed_files = 
    Channel.from(
        "PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap", 
        "NEURAL_21k_significant_EPs",
        "NEURAL_8k_GRB_significant_EPs",
        "NEURAL_14k_noGRB_significant_EPs",
        "notNeural_20k_100flank_noInternalOverlap",
        "sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap"
        )
        .map { ENH_list -> ["${ENH_list}", 
            file("./input/enh_bedfiles/${ENH_list}*.bed", checkIfExists: true)]
            } 
    





//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("/home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()

// UKBB_associations = Channel
//     .fromPath("input/UKBB_res/REC_ALL_BRAIN_EPs_UKBB_eur_associations_GWAMA_format")
//     // .fromPath( ["./input/UKBB_res/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.PHENO1.glm.logistic.hybrid",
//                 // "./input/UKBB_res/SCZ_ORs_PLINK2_logistic_firth_fallback_covar_recessive.frq"])
//     // .collect()
//     // .map{ it -> ["UKBB_eur" , it[0], it[1]] }




workflow lisa_percohort_devel {

    input = LOO_PGC_GWASES.join(validation_samples)
    // input.first().view()
    
    PLINK2_QC_PRUNE_HET (
        input
    )

    // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    //     .join(input)
    //     .view()
    //[celso, celso.prune.in, celso.het, daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, /home/osimoe/PGC_w3_data/celso]

    R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
        PLINK2_QC_PRUNE_HET.out.pruned_variants_het
            .join(input)
        //tuple val(cohort), path ("*het_valid_out*.sample"), path("*a1_UKBBbim*"), path("*mismatching_SNPs*"),  emit: QC_het_a1_mismatch
    )
    
    // R_PRS_QC.out.QC_het_a1_mismatch
    //     .join(input)
    //     .view()
    //     [celso, /home/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/46/8e9146b16c03804507ad3d7dbfa11b/celso_het_valid_out_vs_LOO_GWAS.sample, /home/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/46/8e9146b16c03804507ad3d7dbfa11b/celso_a1_cohort_bim_vs_LOO_GWAS, /home/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/46/8e9146b16c03804507ad3d7dbfa11b/celso_mismatching_SNPs_vs_LOO_GWAS, /home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, /home/osimoe/PGC_w3_data/celso]


    // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
    // also remove mismatching SNPs
    PLINK_PRODUCE_QC_DATASET (
        R_PRS_QC.out.QC_het_a1_mismatch
            .join(input)
        //tuple path ("*.bed"), path ("*.bim"), path ("*.fam"), emit: all_chromosomes_QC
    )

    // PLINK_PRODUCE_QC_DATASET.out.all_chromosomes_QC.view()
    // [celso, celso_QC.bed, celso_QC.bim, celso_QC.fam]
    // [xirwt, xirwt_QC.bed, xirwt_QC.bim, xirwt_QC.fam]
    
    
    
    input //= LOO_PGC_GWASES.join(validation_samples)
        .combine(enhancer_lists_bed_files)
        .map { it.flatten() }
        .set{cohort_GWAS_enh_list}
    // cohort_GWAS_enh_list.view()
    // [celso, /home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, /home/osimoe/PGC_w3_data/celso, PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap, ./input/enh_bedfiles/PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap.bed]
    // [celso, /home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, /home/osimoe/PGC_w3_data/celso, NEURAL_21k_significant_EPs, ./input/enh_bedfiles/NEURAL_21k_significant_EPs.bed]
    // [clz2a, /home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, /home/osimoe/PGC_w3_data/clz2a, PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap, ./input/enh_bedfiles/PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap.bed]
    // [clz2a, /home/osimoe/sumstats_w3/autosomes/PRS/ancestry-specific/study_scores_deduped_eur/training/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, /home/osimoe/PGC_w3_data/clz2a, NEURAL_21k_significant_EPs, ./input/enh_bedfiles/NEURAL_21k_significant_EPs.bed]



    R_prepare_lists_for_clump (
        // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
        // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
        cohort_GWAS_enh_list
    )
    
    
    // R_prepare_lists_for_clump.out.lists_before_clump
    //     .combine(LD_reference)
    //     .view()
    // val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump
    // [xirwt, daner_PGC_SCZ_w3_76_0518d_eur.noxirwt.gz, PsychENCODE_DER_03b_PFC_enhancers_18k, xirwt_PsychENCODE_DER_03b_PFC_enhancers_18k_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, xirwt_PsychENCODE_DER_03b_PFC_enhancers_18k_PGC__noclump_residual_GWAS_compartment.tsv.gz, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.fam]


    PLINK_clump (
        //CLUMPING of enhancer-based SNP compartments together 
        R_prepare_lists_for_clump.out.lists_before_clump
            .combine(LD_reference)
    )
    // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    //     .view()
    // [xirwt, daner_PGC_SCZ_w3_76_0518d_eur.noxirwt.gz, ALL_BRAIN_EPs_1569, xirwt_ALL_BRAIN_EPs_1569_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, xirwt_ALL_BRAIN_EPs_1569_PGC__noclump_residual_GWAS_compartment.tsv.gz, xirwt_PGC_clumped_SNPs.clumped]
    

    R_split_lists (
        // first annotate SNPs with ES of relevant E-P - for ENH SNP list
        // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
        // ##################################################### CAN MULTIPLY P BY VALUE TO RESTORE ENH SNPS P ###########################################################
        // output separate lists to calculate split PRSs and also merged one
        PLINK_clump.out.clumped_SNPs_and_noclump_lists,
        Channel.fromPath( "./input/ES_multipliers/2023-01-18_NEURAL_ENH_EXP_significant_ES_significant_contact_EPs_gene_brain_exp_plus_100_noOverlap.csv.gz", checkIfExists: true)
    )

    
    R_split_lists.out.partitioned//tuple val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"),       emit: partitioned
        .combine(PLINK_PRODUCE_QC_DATASET.out.all_chromosomes_QC, by: [0,0])//[celso, celso_QC.bed, celso_QC.bim, celso_QC.fam]
        .combine(validation_samples, by: [0,0])
        .combine(LD_reference)
        .set{combined_splitlists_bedfile_QCeddata_LDdata}
    
    // combined_splitlists_bedfile_QCeddata_LDdata.view()
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, NEURAL_8k_GRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/clz2a_NEURAL_8k_GRB_significant_EPs_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/clz2a_NEURAL_8k_GRB_significant_EPs_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/clz2a_NEURAL_8k_GRB_significant_EPs_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ed/031a06dc4fc3d0799b8321422c7a01/celso_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ed/031a06dc4fc3d0799b8321422c7a01/celso_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ed/031a06dc4fc3d0799b8321422c7a01/celso_QC.fam, /home/osimoe/PGC_w3_data/celso, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, NEURAL_21k_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/clz2a_NEURAL_21k_significant_EPs_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/clz2a_NEURAL_21k_significant_EPs_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/clz2a_NEURAL_21k_significant_EPs_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, NEURAL_14k_noGRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/clz2a_NEURAL_14k_noGRB_significant_EPs_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/clz2a_NEURAL_14k_noGRB_significant_EPs_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/clz2a_NEURAL_14k_noGRB_significant_EPs_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, notNeural_20k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/clz2a_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/clz2a_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/clz2a_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
                
    PRSice_calculate_PRS_split_partitions(
        combined_splitlists_bedfile_QCeddata_LDdata
    )
    
    // ########################################### CHANGE NAMES OF MULTIPLIERS ###########################################
    PRS_results = 
        PRSice_calculate_PRS_split_partitions.out.clumped_TS_ENH_GWAS_compartment_PRS
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_merged_GWAS_PRS)
            .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_LOO_GWAS_PRS)
            .map { [it, "e_log_OR_X__log_max_ES_perEnh_contact_X_10",
                        "e_log_OR_X__log_neuron_FANTOM_enh_tpm_1_4_X10"].flatten() }


    // PRS_results.view()
    // [celso_ALL_BRAIN_EPs_1569, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_originalOR.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_residual_GWAS_compartment.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_residual_GWAS_compartment.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_residual_GWAS_compartment.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_merged_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_merged_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_merged_GWAS.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_QC.fam, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_original_LOO_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_original_LOO_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_original_LOO_GWAS.best, max_ES_1plus, ES_1plus_X_log_GTEx_geneBrainExp_1plus]
    
    R_R2_and_logistic_and_quantile_compare (
        PRS_results
    )


}


// 
// srun --pty -t 1-00:00 -p shared -c1 nextflow run https://github.com/emosyne/lisa_percohort_devel -latest -r master -profile lisa -resume 
// --slurmd-debug=error