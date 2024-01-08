process R_GWAS_QC_SBayes {
    // debug true
    container 'emosyne/r_docker:2.01'
    label 'process_high'
    tag "$cohort"
    cache "lenient"
    


    input: 
    // [xs234, /path/xs234_GWAS_QC_noclump.gz, annot_baseline2_2_with_continuous_enhancers, /path/annot_baseline2_2_with_continuous_enhancers.txt.gz, /home/osimoe/private_input_files]
    val(cohort), path (GWAS_QC_noclump), path(private_input_files_path)
    
    

    output:
    // $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
    tuple val(cohort), path ("*"),  emit: SBayes
    
    
    script:
    """
    
    R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump} ${private_input_files_path} 
    
    """
}
