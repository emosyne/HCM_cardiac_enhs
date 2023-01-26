process CLEAN_PLINK_OP_for_GWAMA {
    // debug true
    container 'perl:5.34'
    stageInMode 'copy'
    label 'process_low'
    // maxForks 1
    tag "${cohort}_${EPlist}"
    cache 'lenient'

    input:
    // [celso_eur, ALL_HEART_EPs, REC, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/eb/67e31edcb504abef8643e35027819b/recessive_ALL_HEART_EPs_celso_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/eb/67e31edcb504abef8643e35027819b/recessive_ALL_HEART_EPs_celso_eur.frq]
    // [celso_eur, ALL_HEART_EPs, ADD, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/eb/67e31edcb504abef8643e35027819b/ADD_ALL_HEART_EPs_celso_eur.PHENO1.glm.logistic.hybrid, /project/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel/work/eb/67e31edcb504abef8643e35027819b/ADD_ALL_HEART_EPs_celso_eur.frq]
    tuple val(cohort), val(EPlist), val(assoc_method), path(associations), path(frequencies)
    // each path(hg38ToHg19_chain)
    // each path(GW_LD_blocks)

    output:
    tuple val(EPlist),  val(assoc_method), path("*_associations_GWAMA_format"),       emit: PLINK_assoc_GWAMA_format
    

    script:
    """
    PLINK2GWAMA.pl ${associations} ${frequencies} ${assoc_method}_${EPlist}_${cohort}_associations_GWAMA_format
    
    """
}
// 
    
    // 
  
    // echo \$PATH
    // IFS=':' read -ra ADDR <<< "\$PATH"
    // ls -lh "\${ADDR[0]}"