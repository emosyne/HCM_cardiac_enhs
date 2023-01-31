include { bash_base_GWAS_QC }           from '../modules/local/bash_base_GWAS_QC.nf'
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




// #### UKBB input files ####

genotype_chr_files = Channel
    .fromFilePairs("$geno_input_dir/c*_b0*.{bgen,sample}", flat: true, checkIfExists: true)
    .map{ it-> [it[1],it[2]] }

UKBBethinicityRelatedness = Channel.fromPath( './input/biobank/EIDs_nonBritIrish_includingsecondary_or_related_over_king125.tsv' , checkIfExists: true)
HCM_UKBB_pheno =    Channel.fromPath("./input/biobank/HCM.pheno", checkIfExists: true)
UKBB_covariates =   Channel.fromPath('./input/biobank/non_missing_10PCs_Jun22.covariate.gz', checkIfExists: true)


//LD ref
LD_reference = Channel.from("bed","bim","fam") 
    .map { ext -> [file("$ld_ref_dir/EUR_phase3_autosomes_hg19.${ext}")] }
            .collect()


full_GWAS_hg19 = Channel
    // .fromPath("./input/textfiles/fullGWAS_SCZ_PGC3_SCZ_wave3.european.autosome.public.v3_tidy_hg19.tsv.gz", checkIfExists: true) 
    .fromPath("$GWAS_dir/hcm.gwama.sumstats_hg19_24Feb21.gz", checkIfExists: true) 



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
    // // BASE =   SEAN HCM GWAS
    // // TARGET = UKBB


    // // BASE (GWAS) QC: REMOVE LOW MAF AND INFO SCORES
    // //produce GWAS_QC
    // bash_base_GWAS_QC (
    //     input
    // )

    // // TARGET QC 1: PRUNE AND HETEROZIGOSITY CALCULATIONS
    // // produce prune.in and het files
    // PLINK2_QC_PRUNE_HET (
    //     input
    // )

    // // PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    // //         .join(bash_base_GWAS_QC.out.GWAS_QC)
    // //         .view()
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0e/ce0ddf79c4e1e43923154fa08368cc/celso.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/0e/ce0ddf79c4e1e43923154fa08368cc/celso.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c4/f0be9b4b466c6063bbb943a07381f8/celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/c4/f0be9b4b466c6063bbb943a07381f8/celso_GWAS_QC.gz]
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ff/060e35ff0fa0dd4990a129e665205c/clz2a.prune.in, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ff/060e35ff0fa0dd4990a129e665205c/clz2a.het, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/de/b15718f7c923436610bb24cdb9f6e0/clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/de/b15718f7c923436610bb24cdb9f6e0/clz2a_GWAS_QC.gz]

    // // TARGET QC 2:  remove heterogeneity outliers, produced A1 alleles, and mismatching SNPs list to be removed
    // // produce QC_het_a1_mismatch, 
    // R_PRS_QC ( // calculates mismatching SNPs and recodes all alleles to GWAS base
    //     PLINK2_QC_PRUNE_HET.out.pruned_variants_het
    //         .join(bash_base_GWAS_QC.out.GWAS_QC)
        
    // )
    
    // // TARGET QC 3:  
    // // Remove individuals with heterozigosity F coefficients that are more than 3 standard deviation (SD) units from the mean
    // // also remove mismatching SNPs
    // // also standard QC --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 
    // PLINK_PRODUCE_QC_DATASET (
    //     R_PRS_QC.out.QC_het_a1_mismatch

    //     //tuple path ("*.bed"), path ("*.bim"), path ("*.fam"), emit: target_QC
    // )

    // // PLINK_PRODUCE_QC_DATASET.out.target_QC.view()
    // // // [celso, celso_QC.bed, celso_QC.bim, celso_QC.fam]
    // // // [xirwt, xirwt_QC.bed, xirwt_QC.bim, xirwt_QC.fam]
    
    
    
    // bash_base_GWAS_QC.out.GWAS_QC
    //     .combine(enhancer_lists_bed_files)
    //     .map { it.flatten() }
    //     .set{cohort_GWAS_enh_list}
    
    // // cohort_GWAS_enh_list.view()
    // // [xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234_GWAS_QC.gz, NEURAL_8k_GRB_significant_EPs, ./input/enh_bedfiles/NEURAL_8k_GRB_significant_EPs.bed]
    // // [xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234_GWAS_QC.gz, NEURAL_14k_noGRB_significant_EPs, ./input/enh_bedfiles/NEURAL_14k_noGRB_significant_EPs.bed]
    // // [xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234_GWAS_QC.gz, notNeural_20k_100flank_noInternalOverlap, ./input/enh_bedfiles/notNeural_20k_100flank_noInternalOverlap.bed]
    // // [xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/bf/d7f39c678de5caf80ae6980c7984b8/xs234_GWAS_QC.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, ./input/enh_bedfiles/sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap.bed]
    // // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3_GWAS_QC.gz, PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap, ./input/enh_bedfiles/PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap.bed]
    // // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3_GWAS_QC.gz, NEURAL_21k_significant_EPs, ./input/enh_bedfiles/NEURAL_21k_significant_EPs.bed]
    // // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3_GWAS_QC.gz, NEURAL_8k_GRB_significant_EPs, ./input/enh_bedfiles/NEURAL_8k_GRB_significant_EPs.bed]
    // // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3_GWAS_QC.gz, NEURAL_14k_noGRB_significant_EPs, ./input/enh_bedfiles/NEURAL_14k_noGRB_significant_EPs.bed]
    // // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3_GWAS_QC.gz, notNeural_20k_100flank_noInternalOverlap, ./input/enh_bedfiles/notNeural_20k_100flank_noInternalOverlap.bed]
    // // [xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/dead58096a916b16e8f35fac26851e/xclo3_GWAS_QC.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, ./input/enh_bedfiles/sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap.bed]

    // // BASE subsetting
    // R_prepare_lists_for_clump (
    //     // SUBSETS GWAS SNPS INTO ENH COMPARTMENT AND RESIDUAL COMPARTMENT.
    //     // ########################### IN PREPARATION FOR CLUMPING, DIVIDE P VALUES FOR ENH SNPS BY X TO PRESERVE ENH SNPS ###########################
    //     cohort_GWAS_enh_list
    // )
    
    
    // // R_prepare_lists_for_clump.out.lists_before_clump
    // //     .combine(LD_reference)
    // //     .view()
    // // val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz"), path("*_PGC__noclump_residual_GWAS_compartment.tsv.gz"), emit: lists_before_clump
    // // [xirwt, daner_PGC_SCZ_w3_76_0518d_eur.noxirwt.gz, PsychENCODE_DER_03b_PFC_enhancers_18k, xirwt_PsychENCODE_DER_03b_PFC_enhancers_18k_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, xirwt_PsychENCODE_DER_03b_PFC_enhancers_18k_PGC__noclump_residual_GWAS_compartment.tsv.gz, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/.nextflow/assets/emosyne/LD_ref/EUR_phase3_autosomes_hg19.fam]


    // PLINK_clump (
    //     //CLUMPING of enhancer-based SNP compartments together 
    //     R_prepare_lists_for_clump.out.lists_before_clump
    //         .combine(LD_reference)
    // )
    // // PLINK_clump.out.clumped_SNPs_and_noclump_lists
    // //     .view()
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_GWAS_QC.gz, NEURAL_8k_GRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_NEURAL_8k_GRB_significant_EPs_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_NEURAL_8k_GRB_significant_EPs_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/a1/90c446ff1933f735dc31ec989df901/celso_PGC_clumped_SNPs.clumped]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_GWAS_QC.gz, notNeural_20k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_notNeural_20k_100flank_noInternalOverlap_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_notNeural_20k_100flank_noInternalOverlap_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/95/93aeae17e1e7b3b5df6b320b5c8064/celso_PGC_clumped_SNPs.clumped]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_GWAS_QC.gz, NEURAL_14k_noGRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_NEURAL_14k_noGRB_significant_EPs_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_NEURAL_14k_noGRB_significant_EPs_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/54/65433ffc2328c2dec9fe2a5be9431c/celso_PGC_clumped_SNPs.clumped]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_GWAS_QC.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_PGC__noclump_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_PGC__noclump_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/80/8446a04f70e95bc51c4e693a56856b/celso_PGC_clumped_SNPs.clumped]

    

    // R_split_lists (
    //     // first annotate SNPs with ES of relevant E-P - for ENH SNP list
    //     // ##################################################### GENERATE MODIFIED ORS MULT BY ES OR EXP       ###########################################################
    //     // ##################################################### CAN MULTIPLY P BY VALUE TO RESTORE ENH SNPS P ###########################################################
    //     // output separate lists to calculate split PRSs and also merged one
    //     PLINK_clump.out.clumped_SNPs_and_noclump_lists,
    //     Channel.fromPath( "./input/ES_multipliers/2023-01-18_NEURAL_ENH_EXP_significant_ES_significant_contact_EPs_gene_brain_exp_plus_100_noOverlap.csv.gz", checkIfExists: true)
    // )

    
    // R_split_lists.out.partitioned//tuple val(cohort), path (LOO_GWAS),  val(ENH_list), path("*_clumped_TS_ENH_GWAS_compartment.tsv.gz"), path("*_clumped_residual_GWAS_compartment.tsv.gz"), path("*_clumped_merged_GWAS.tsv.gz"),       emit: partitioned
    //     .combine(PLINK_PRODUCE_QC_DATASET.out.target_QC, by: [0,0])//[celso, celso_QC.bed, celso_QC.bim, celso_QC.fam]
    //     .combine(validation_samples, by: [0,0])
    //     .combine(LD_reference)
    //     .set{combined_splitlists_bedfile_QCeddata_LDdata}
    
    // // combined_splitlists_bedfile_QCeddata_LDdata.view()
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/09/107e03de16b78c366c461b5679c249/clz2a_PsychENCODE_DER_03b_PFC_enhancers_18k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, NEURAL_8k_GRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/clz2a_NEURAL_8k_GRB_significant_EPs_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/clz2a_NEURAL_8k_GRB_significant_EPs_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/45/dd0a89067854053b4195019fe91310/clz2a_NEURAL_8k_GRB_significant_EPs_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // // [celso, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/daner_PGC_SCZ_w3_76_0518d_eur.nocelso.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/70/dd6a60633015511974424e179ce389/celso_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ed/031a06dc4fc3d0799b8321422c7a01/celso_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ed/031a06dc4fc3d0799b8321422c7a01/celso_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ed/031a06dc4fc3d0799b8321422c7a01/celso_QC.fam, /home/osimoe/PGC_w3_data/celso, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, NEURAL_21k_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/clz2a_NEURAL_21k_significant_EPs_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/clz2a_NEURAL_21k_significant_EPs_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/e3/18263c3d098e1659610cf97d639259/clz2a_NEURAL_21k_significant_EPs_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, NEURAL_14k_noGRB_significant_EPs, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/clz2a_NEURAL_14k_noGRB_significant_EPs_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/clz2a_NEURAL_14k_noGRB_significant_EPs_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9a/0fbbf370d6d51a4ec67f1becdcc560/clz2a_NEURAL_14k_noGRB_significant_EPs_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, notNeural_20k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f7/81dd7f2a86ba496a7f939f819fee7b/clz2a_notNeural_20k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
    // // [clz2a, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/daner_PGC_SCZ_w3_76_0518d_eur.noclz2a.gz, sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/clz2a_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_TS_ENH_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/clz2a_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_residual_GWAS_compartment.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f3/623a209f5434bb745ac2c6c77c65ec/clz2a_sig_ES_sig_contact_EPs_in_Neural_22k_100flank_noInternalOverlap_clumped_merged_GWAS.tsv.gz, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bed, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.bim, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/f0/c7b8deb7b44c89970443f95440e801/clz2a_QC.fam, /home/osimoe/PGC_w3_data/clz2a, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bed, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.bim, /home/osimoe/project/LD_ref/EUR_phase3_autosomes_hg19.fam]
                
    // PRSice_calculate_PRS_split_partitions(
    //     combined_splitlists_bedfile_QCeddata_LDdata
    // )
    
    // // ########################################### CHANGE NAMES OF MULTIPLIERS ###########################################
    // PRS_results = 
    //     PRSice_calculate_PRS_split_partitions.out.clumped_TS_ENH_GWAS_compartment_PRS
    //         .join(PRSice_calculate_PRS_split_partitions.out.clumped_residual_GWAS_compartment_PRS)
    //         .join(PRSice_calculate_PRS_split_partitions.out.clumped_merged_GWAS_PRS)
    //         .join(PRSice_calculate_PRS_split_partitions.out.clumped_original_LOO_GWAS_PRS)
    //         .map { [it, "e_log_OR_X__log_max_ES_perEnh_contact_X_10",
    //                     "e_log_OR_X__log_neuron_FANTOM_enh_tpm_1_4_X10"].flatten() }


    // // PRS_results.view()
    // // [celso_ALL_BRAIN_EPs_1569, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_originalOR.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_originalOR.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure1.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_OR_by_measure2.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_TS_ENH_GWAS_compartment_originalOR.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_residual_GWAS_compartment.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_residual_GWAS_compartment.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_residual_GWAS_compartment.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_merged_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_merged_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_clumped_merged_GWAS.best, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_QC.fam, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_original_LOO_GWAS.summary, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_original_LOO_GWAS.prsice, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/28/c71a646089a6bd625df05cc916f4ce/celso_ALL_BRAIN_EPs_1569_original_LOO_GWAS.best, max_ES_1plus, ES_1plus_X_log_GTEx_geneBrainExp_1plus]
    
    // R_R2_and_logistic_and_quantile_compare (
    //     PRS_results
    // )

}


// 
// srun --pty -t 1-00:00 -p shared -c1 nextflow run https://github.com/emosyne/HCM -latest -r master -profile lisa -resume 
// --slurmd-debug=error