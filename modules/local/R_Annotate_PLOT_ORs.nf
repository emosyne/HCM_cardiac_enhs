process R_Annotate_PLOT_ORs {
    // debug true
    container 'emosyne/r_docker:1.94'
    // stageInMode 'copy'
    label 'process_high'
    tag "$EP_list"
    cache 'lenient'
    // maxForks 1

    input:
    //[ALL_BRAIN_EPs, [/project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/03/7dcfa63bb1b2d1d68551035960c5b9/ALL_BRAIN_EPs_REC_GWAMA_associations_fixedES.out, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/01/ea2eb007917eb84feb1fd19d85ebc0/ALL_BRAIN_EPs_ADD_GWAMA_associations_fixedES.out], [/project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/08/bb4ca9cfcd5c118de82094e64f66a4/recessive_ALL_BRAIN_EPs_clz2a_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/3b/d8739f5ad8d2a14b49fa7b45bd2b6a/recessive_ALL_BRAIN_EPs_celso_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ad/aadd491121aae0c23506ff2549bf8f/recessive_ALL_BRAIN_EPs_gpc2a_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/1f/be6a28389d9341016402eeb1fedf9d/recessive_ALL_BRAIN_EPs_mcqul_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/85/95f9c5b6b9123c369a28429444307c/recessive_ALL_BRAIN_EPs_xboco_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/94/1e483701ec9126f33c85c10a73ccd0/recessive_ALL_BRAIN_EPs_xclm2_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/3e/f903a90ae2b69883f01f9fef1a77cf/recessive_ALL_BRAIN_EPs_xclo3_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/4b/85bc23cf9d856a4feecfd935f00d7b/recessive_ALL_BRAIN_EPs_xgras_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9d/64e63f2d56f646d2ceb638b9f273b6/recessive_ALL_BRAIN_EPs_xirwt_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/24/6607212bd6c734934887f7ba2553e0/recessive_ALL_BRAIN_EPs_xs234_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/8b/e5e97583dd3c323223c8fef83782dc/recessive_ALL_BRAIN_EPs_xswe5_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9b/633ccb064e3250ceaf74e67103e37b/recessive_ALL_BRAIN_EPs_xswe6_eur.PHENO1.glm.logistic.hybrid], [/project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/08/bb4ca9cfcd5c118de82094e64f66a4/recessive_ALL_BRAIN_EPs_clz2a_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/3b/d8739f5ad8d2a14b49fa7b45bd2b6a/recessive_ALL_BRAIN_EPs_celso_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/ad/aadd491121aae0c23506ff2549bf8f/recessive_ALL_BRAIN_EPs_gpc2a_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/1f/be6a28389d9341016402eeb1fedf9d/recessive_ALL_BRAIN_EPs_mcqul_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/85/95f9c5b6b9123c369a28429444307c/recessive_ALL_BRAIN_EPs_xboco_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/94/1e483701ec9126f33c85c10a73ccd0/recessive_ALL_BRAIN_EPs_xclm2_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/3e/f903a90ae2b69883f01f9fef1a77cf/recessive_ALL_BRAIN_EPs_xclo3_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/4b/85bc23cf9d856a4feecfd935f00d7b/recessive_ALL_BRAIN_EPs_xgras_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9d/64e63f2d56f646d2ceb638b9f273b6/recessive_ALL_BRAIN_EPs_xirwt_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/24/6607212bd6c734934887f7ba2553e0/recessive_ALL_BRAIN_EPs_xs234_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/8b/e5e97583dd3c323223c8fef83782dc/recessive_ALL_BRAIN_EPs_xswe5_eur.frq, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/9b/633ccb064e3250ceaf74e67103e37b/recessive_ALL_BRAIN_EPs_xswe6_eur.frq]]
    tuple val(EP_list),  path(meta_analysis_results), path(PLINK_assoc), path(PLINK_freq)
    each path(enhancer_annotations_hg19)
    each path(hg38ToHg19_chain)
    each path(PGC_GWAS_GW_LD_blocks_merge_hg19)
    each path(full_GWAS_hg19)
    
    output:
    tuple val(EP_list),  path ("*_annotated_ORs.csv")        , emit: annotated_ORs
    // path ("*")
    path ("figs/*.png")

    script:
    """
    echo ${PLINK_assoc} | tr ' ' '\\n' > PLINK_assoc_filelist.txt
    echo ${meta_analysis_results} | tr ' ' '\\n' > meta_analysis_results_filelist.txt

    R_Annotate_PLOT_ORs.R meta_analysis_results_filelist.txt ${EP_list} ${enhancer_annotations_hg19} ${hg38ToHg19_chain} ${PGC_GWAS_GW_LD_blocks_merge_hg19} PLINK_assoc_filelist.txt 
    
    
   
    """
}
    