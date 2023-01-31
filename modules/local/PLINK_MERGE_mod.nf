process PLINK_MERGE {
        
    
    // maxForks 1
    // errorStrategy 'ignore' 
    // scratch '/rds/general/user/eosimo/ephemeral/' TO BE RESTORED FOR RUNNING IN IMPERIAL
    label 'vlarge'
    container 'emosyne/plink2:1.23'
    cache "lenient"


    input:
    path(genomes)
    each path (PLINKethinicityRelatedness)
    each path (HCM_UKBB_pheno)
    

    output:
    tuple path ("*.bed"), path ("*.bim"), path ("*.fam"), emit: all_chromosomes_extracted
    path "chr_file_list.txt"
    path ("*.log")
    // path "versions1.yml"           , emit: versions

    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """
    # IMPORT LIST OF FILES INTO chr_file_list
    echo ${genomes} | \\
        tr -d ',[]' | \\
        tr ' ' '\\n' | grep 'bed' | \\
        sed 's/.bed//g' > \\
        chr_file_list.txt
    # take first chrom file
    first_chr=\$(head -n 1 chr_file_list.txt)
    echo \$first_chr
    # REMOVE FIRST CHR FROM FILE
    sed -i '1d' chr_file_list.txt 
    
    plink --bfile \$first_chr \\
        --merge-list chr_file_list.txt \\
        --make-bed \\
        --remove ${PLINKethinicityRelatedness} \\
        --pheno ${HCM_UKBB_pheno} \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out GWAS_ENH_SNPs_hg19_ALLCHR


    
    """
}
