process bash_extract_numSNP {
    // debug true
    // container 'emosyne/r_docker:1.94'
    // container 'emosyne/simpler:1.1'
    label 'process_low'
    tag "$EP_list"
    cache "lenient"
    

    input:
    tuple val (EP_list), \
        path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_summary), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_prsice), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_best), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_005),  path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_005numsnp), \
        path(EPs_clumped_NonOverlapOriginal_summary), path(EPs_clumped_NonOverlapOriginal_prsice), path(EPs_clumped_NonOverlapOriginal_best), path(EPs_clumped_NonOverlapOriginal_005), path(EPs_clumped_NonOverlapOriginal_005numsnp), \
        path(ADD_clumped_original_PRS_summary), path(ADD_clumped_original_PRS_prsice), path (ADD_clumped_original_PRS_best), path (ADD_clumped_original_PRS_005), path (ADD_clumped_original_PRS_005numsnp)

    output:
    tuple  val(EP_list), \
        path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_summary), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_prsice), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_best), path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_005),  path(ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_005numsnp), \
        path(EPs_clumped_NonOverlapOriginal_summary), path(EPs_clumped_NonOverlapOriginal_prsice), path(EPs_clumped_NonOverlapOriginal_best), path(EPs_clumped_NonOverlapOriginal_005), path(EPs_clumped_NonOverlapOriginal_005numsnp), \
        path(ADD_clumped_original_PRS_summary), path(ADD_clumped_original_PRS_prsice), path (ADD_clumped_original_PRS_best), path (ADD_clumped_original_PRS_005), path (ADD_clumped_original_PRS_005numsnp),\
        emit: PRS_results_correct_n_SNPs

    script:
    """
    cat ${ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_prsice} | extractSNPn_fix.sh > ${ADDPRS_clumpedOriginalGWAS_NoEPs_overlap_005numsnp}
    cat ${EPs_clumped_NonOverlapOriginal_prsice} | extractSNPn_fix.sh > ${EPs_clumped_NonOverlapOriginal_005numsnp}
    cat ${ADD_clumped_original_PRS_prsice} | extractSNPn_fix.sh > ${ADD_clumped_original_PRS_005numsnp}
    
    

   
    """
}
