process GWAMA_meta_analysis {
    // debug true
    container 'emosyne/prsice_gwama_exec:1.0'
    stageInMode 'copy'
    label 'process_low'
    tag "${EP_list}_${assoc_method}"
    cache 'lenient'

    input:
    tuple val(EP_list), val(assoc_method), path(associations_GWAMA_format)

    output:
    tuple val(EP_list), val(assoc_method), path ("*_GWAMA_associations_fixedES.out"), emit: meta_analysis_results
    path("*")
    

    script:
    """
    echo ${associations_GWAMA_format} | tr ' ' '\\n' > gwama.in
    
    GWAMA --filelist gwama.in -o ${EP_list}_${assoc_method}_GWAMA_associations_fixedES
    
    
    """
}
