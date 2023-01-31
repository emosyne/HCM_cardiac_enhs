process PLINK2_QC_PRUNE_HET {
    // tag "$cohort"
    // debug true
    label 'process_high'
    container 'emosyne/plink2:1.23'
    cache "lenient"

    input: 
    tuple path(bed), path(bim), path(fam)
    

    output:
    tuple path(bed), path(bim), path(fam), path ("*.prune.in"), path ("*.het"), emit: pruned_variants_het
    // path("*.log")
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    plink \\
      --bfile ${bed.simpleName} \\
      --indep-pairwise 200 50 0.25 \\
      --out ${bed.simpleName} \\
      --threads $task.cpus \\
      --memory $mem_mb

    plink \\
      --bfile ${bed.simpleName} \\
      --extract ${bed.simpleName}.prune.in \\
      --het \\
      --out ${bed.simpleName} \\
      --threads $task.cpus \\
      --memory $mem_mb
    """
}

