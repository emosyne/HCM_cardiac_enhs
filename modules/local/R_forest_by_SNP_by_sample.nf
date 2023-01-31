process R_forest_by_SNP_by_sample {
    // debug true
    container 'emosyne/r_docker:1.94'
    stageInMode 'copy'
    label 'process_low'
    // tag "$meta"
    cache 'lenient'

    input:
    path(associations_GWAMA_format)
    // val(SNPlist)

    output:
    path("*")
    

    script:
    """
    echo ${associations_GWAMA_format} | tr ' ' '\\n' > associations_GWAMA_format_filelist.txt


    R_forest_by_SNP_by_sample.R associations_GWAMA_format_filelist.txt
    
    
    
    """
}
// 