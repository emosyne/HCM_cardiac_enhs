process PLINK2_ASSOC_GLM {
    tag "$meta"
    // debug true
    label 'process_high'
    maxForks 2
    container 'emosyne/plink2:1.23'
    // errorStrategy 'ignore'
    cache 'lenient'

    input: 
    tuple val(meta), path(cohort_dir)
    each path(EP_list_hg19)
    

    output:
    tuple env(meta2), env(EPlist), val("REC"), path ("recessive_*.glm.logistic.hybrid"), path ("recessive_*.frq"),  emit: associations_REC
    tuple env(meta2), env(EPlist), val("ADD"), path ("ADD_*.glm.logistic.hybrid"), path ("ADD_*.frq"),              emit: associations_ADD
    path("*")

    
    script:
    def mem_mb = (task.memory * 0.95).toMega()
    def EPlist = EP_list_hg19.simpleName
    """ 
    EPlist=`echo ${EP_list_hg19} | sed 's/.bed//'`

    bedfile="${cohort_dir}/imputed/hardcall_genotypes/*.bed"
    bedfile2=`echo \$bedfile | sed 's/.bed//'`
    covariates="${cohort_dir}/prin_comp/*.mds"

    echo \$bedfile2

    IFS='_' read -r -a array <<< \$bedfile2
    popu=`echo "\${array[3]}"`
    IFS='-' read -r -a array <<< \$popu
    pop=`echo "\${array[0]}"`

    
    echo \$pop
    
    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile \$bedfile2 \\
        --extract bed1 ${EP_list_hg19} \\
        --write-snplist \\
        --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-6 \\
        --glm recessive firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar \$covariates --covar-name C1-C10 \\
        --out recessive_${EPlist}_${meta}_\$pop
    
    plink \\
        --bfile \$bedfile2 \\
        --extract recessive_${EPlist}_${meta}_\$pop.snplist \\
        --freq \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out recessive_${EPlist}_${meta}_\$pop

    meta2=${meta}_\$pop
    echo \$meta2

    plink2 \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --bfile \$bedfile2 \\
        --extract bed1 ${EP_list_hg19} \\
        --write-snplist \\
        --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-6 \\
        --glm firth-fallback omit-ref hide-covar \\
        --ci 0.95 \\
        --covar \$covariates --covar-name C1-C10 \\
        --out ADD_${EPlist}_${meta}_\$pop
    plink \\
        --bfile \$bedfile2 \\
        --extract ADD_${EPlist}_${meta}_\$pop.snplist \\
        --freq \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out ADD_${EPlist}_${meta}_\$pop

    """

      

}
