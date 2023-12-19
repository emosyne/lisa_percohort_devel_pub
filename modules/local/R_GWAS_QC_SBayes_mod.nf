process R_GWAS_QC_SBayes {
    // debug true
    container 'emosyne/r_docker:2.0'
    label 'process_high'
    tag "$cohort"
    cache "lenient"
    


    input: 
    tuple val(cohort), path (GWAS_QC_noclump), path (private_input_files_path)
    
    

    output:
    tuple val(cohort), path ("*"),  emit: SBayes
    
    
    script:
    """
    R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump}  ${private_input_files_path}
    
    """
}
