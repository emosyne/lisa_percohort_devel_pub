process R_GWAS_QC_SBayes {
    // debug true
    container 'emosyne/r_docker:2.0'
    label 'process_high'
    tag "$cohort"
    cache "lenient"
    


    input: 
    // [xs234, /path/to/xs234_GWAS_QC_noclump.gz, annot_baseline2_2_with_continuous_enhancers, /home/osimoe/private_input_files/SBayes_annots/annot_baseline2_2_with_continuous_enhancers.txt.gz]
    tuple val(cohort), path (GWAS_QC_noclump), \
        val (SBayesRC_annot), path (private_input_files_path)
    
    

    output:
    tuple val(cohort), path ("*"),  emit: SBayes
    
    
    script:
    """
    R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump}  ${private_input_files_path}
    
    """
}
