process PLINK2_ASSOC_GLM_nocov {
    tag "$meta"
    // debug true
    label 'process_low'
    container 'emosyne/plink2:1.23'
    // errorStrategy 'ignore'
    cache 'lenient'

    input: 
    tuple val(meta), path(cohort_dir)
    each path(EP_list_hg19)
    

    output:
    tuple env(meta2), env(EPlist), path ("*.glm.logistic.hybrid"), path ("*.frq"), emit: associations
    path("*")
    // path "versions2.yml"           , emit: versions
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args
    //def prefix = task.ext.prefix 
    
    def mem_mb = task.memory.toMega()
    // def meta = cohort_dir.simpleName
    script:
    def EPlist = EP_list_hg19.simpleName
    // publishDir "${params.publishDir}/EPlist", mode: 'copy'
    if(workflow.profile == "lisa") {
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
            --glm recessive firth-fallback omit-ref hide-covar allow-no-covars\\
            --ci 0.95 \\
            --out ${EPlist}_${meta}_\$pop
        
        plink \\
            --bfile \$bedfile2 \\
            --extract ${EPlist}_${meta}_\$pop.snplist \\
            --freq \\
            --out ${EPlist}_${meta}_\$pop

        meta2=${meta}_\$pop
        echo \$meta2

        """
    } else {
        """ 
        EPlist=`echo ${EP_list_hg19} | sed 's/.bed//'`

        bedfile="${cohort_dir}/imputed/hardcall_genotypes/*.bed"
        bedfile2=`echo \$bedfile | sed 's/.bed//'`
        covariates="${cohort_dir}/prin_comp/*.covariate.gz"
        
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
            --write-snplist \\
            --maf 0.01 --mac 20 --geno 0.1 --hwe 1e-6 \\
            --glm recessive firth-fallback omit-ref hide-covar  allow-no-covars\\
            --ci 0.95 \\
            --out ${EPlist}_${meta}_\$pop
        
       plink \\
            --bfile \$bedfile2 \\
            --extract ${EPlist}_${meta}_\$pop.snplist \\
            --freq \\
            --out ${EPlist}_${meta}_\$pop
    
        meta2=${meta}_\$pop
        echo \$meta2
        
        """
        
    }


      

//     cat <<-END_VERSIONS > versions2.yml
//     "${task.process}":
//         plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
//     END_VERSIONS
}
