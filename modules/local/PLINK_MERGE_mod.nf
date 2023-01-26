process PLINK_MERGE {
        
    label 'vlarge'
    // maxForks 1
    // tag "${chr_bgenfile.simpleName}"
    // errorStrategy 'ignore' 
    scratch '/rds/general/user/eosimo/ephemeral/'
    container 'emosyne/plink2:1.23'
    cache "lenient"


    input:
    path(genomes)
    each path (PLINKethinicityRelatedness)
    each path (SCZ_UKBB_pheno)
    

    output:
    tuple path ("*.bed"), path ("*.bim"), path ("*.fam"), emit: all_chromosomes_extracted
    path "chr_file_list.txt"
    path ("*.log")
    // path "versions1.yml"           , emit: versions

    script:
    def mem_mb = (task.memory * 0.9).toMega()

    // if( "$bed" == "${prefix}.bed" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
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
        --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 --mind 0.1 \\
        --remove ${PLINKethinicityRelatedness} \\
        --pheno ${SCZ_UKBB_pheno} \\
        --threads $task.cpus \\
        --memory $mem_mb \\
        --out fullPGC_GWAS_plus_allEPlists_SNPs_hg19_ALLCHR


    
    """
}
