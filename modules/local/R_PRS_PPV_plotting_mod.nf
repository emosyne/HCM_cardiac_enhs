process R_PRS_PPV_plotting {
    // debug true
    container 'emosyne/r_docker:1.92'
    // container 'emosyne/simpler:1.1' does not work here
    label 'process_low'
    tag "$EP_list"
    errorStrategy 'ignore'

    input: 
    //[EP_NEGATIVE_metaPGC_ADD, 
    //EP_NEGATIVE_metaPGC_ADD_ADDPRS_clumpedOriginalGWAS_NoEPs_overlap.summary, EP_NEGATIVE_metaPGC_ADD_ADDPRS_clumpedOriginalGWAS_NoEPs_overlap.prsice, EP_NEGATIVE_metaPGC_ADD_ADDPRS_clumpedOriginalGWAS_NoEPs_overlap.best, 
    //EP_NEGATIVE_metaPGC_ADD_ADDPRS_EPs_clumped_NonOverlapOriginal.summary, EP_NEGATIVE_metaPGC_ADD_ADDPRS_EPs_clumped_NonOverlapOriginal.prsice, EP_NEGATIVE_metaPGC_ADD_ADDPRS_EPs_clumped_NonOverlapOriginal.best, 
    //clumped_original_PRS.summary, clumped_original_PRS.prsice, clumped_original_PRS.best]
    tuple val (EP_list), \
        path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_summary), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_prsice), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_best), \
        // path(REC_all_TS_EPs_summary), path(REC_all_TS_EPs_prsice), path (REC_all_TS_EPs_best), \
        path(EPs_clumped_NonOverlapOriginal_summary), path(EPs_clumped_NonOverlapOriginal_prsice), path(EPs_clumped_NonOverlapOriginal_best), \
        path(ADD_clumped_original_PRS_summary), path(ADD_clumped_original_PRS_prsice), path (ADD_clumped_original_PRS_best)
    
    output:
    // tuple val(meta), path("*_modified_GWAS.csv"),       emit: modified_GWAS
    path("*.png")


    script:
    """
    R2_compare_by_PRS.R ${EP_list} ${ADD_clumped_original_PRS_prsice} ${ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_prsice} ${EPs_clumped_NonOverlapOriginal_prsice} 
    
    
    """
}
    