process R_GWAS_QC_SBayes {
    // debug true
    container 'emosyne/r_docker:2.01'
    label 'process_high_long'
    tag "$cohort"
    cache "lenient"
    


    input: 
    // [xs234, /gpfs/home2/osimoe/emanuele_project/lisa_percohort_devel_pub_work_dir/5e/c387425b77165778812bfd909a8217/xs234_GWAS_QC_noclump.gz, /home/osimoe/emanuele_project/private_input_files]
    tuple val(cohort), path (GWAS_QC_noclump), path(private_input_files_path)
    
    

    output:
    // $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
    tuple val(cohort), path ("*"),  emit: SBayes
    
    
    script:
    """
    
    R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump} ${private_input_files_path} 
    
    """
}
