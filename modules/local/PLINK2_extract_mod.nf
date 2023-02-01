process PLINK2_EXTRACT {
    // debug true
    label 'process_high_memory'
    // maxForks 1
    tag "${chr_bgenfile.simpleName}"
    // errorStrategy 'ignore' 
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    tuple path(extracted_GWAS_SNPs_bed), path(chr_bgenfile), path(chr_samplefile)
    

    output:
    tuple path("*.bim"), path("*.bed"), path ("*.fam"),  emit: SNPextracted_by_chromosome
    path("*.log")

    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    echo $mem_mb

    plink2 \\
        --bgen ${chr_bgenfile} ref-first \\
        --sample ${chr_samplefile} \\
        --chr 1-22\\
        --rm-dup force-first --make-bed \\
        --extract bed1 ${extracted_GWAS_SNPs_bed} \\
        --out extracted_GWAS_SNPs_${chr_bgenfile.simpleName} \\
        --threads $task.cpus \\
        --memory $mem_mb

    """
}
