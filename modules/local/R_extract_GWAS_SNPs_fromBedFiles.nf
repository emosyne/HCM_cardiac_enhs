process extract_GWAS_SNPs_fromBedFiles {
    // debug true
    maxForks 2
    container 'emosyne/r_docker:1.92'
    // container 'emosyne/simpler:1.1'
    label 'process_high_resource_short'
    // tag "$meta"
    cache "lenient"

    
    input:
    path(collected_bed_files_for_enhancers)
    each path(full_GWAS_hg19)

    output:
    path("GWAS_SNPs_in_initial_bed_files.bed"), emit: GWAS_SNPs_in_initial_bed_files
    

    script:
    """
    echo ${collected_bed_files_for_enhancers} | tr ' ' '\\n' > collected_bed_files_for_enhancers.txt

    extract_GWAS_SNPs_fromBedFiles.R ${full_GWAS_hg19} collected_bed_files_for_enhancers.txt
    
    """
}
