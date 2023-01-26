process R_plot_GO {
    // debug true
    container 'emosyne/r_docker:1.92'
    label 'process_low'
    tag "$EP_list"
    cache "lenient"
    

    input:
    tuple val(EP_list), path (annotated_ORs)
    each path(full_GWAS_hg19)
    each path(EP_annotations)
    // tuple val(meta), path(annotated_ORs), path(full_GWAS_hg19), path(ENH_PROM_hg38)

    output:
    path("*")
    path("versions.yml")                        ,       emit: versions

    script:
    """
    R_plot_GO.R ${EP_list} ${annotated_ORs} ${full_GWAS_hg19} ${EP_annotations}
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Bash: \$(echo "\$BASH_VERSION")
        R: \$(R --version | head -1 | awk '{print \$3}')
        R_dplyr: \$(Rscript -e 'packageVersion("dplyr")' | awk '{print \$2}' | tr -d "‘’")
        R_ggplot2: \$(Rscript -e 'packageVersion("ggplot2")' | awk '{print \$2}' | tr -d "‘’")
    END_VERSIONS
   
    """
}
    