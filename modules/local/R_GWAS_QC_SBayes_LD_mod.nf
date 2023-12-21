process R_GWAS_QC_SBayes_LD {
    // debug true
    container 'emosyne/r_docker:2.0'
    label 'process_high_memory'
    tag "$cohort"
    cache "lenient"
    


    input: 
    // [xs234, /path/xs234_GWAS_QC_noclump.gz, annot_baseline2_2_with_continuous_enhancers, /path/annot_baseline2_2_with_continuous_enhancers.txt.gz, /home/osimoe/private_input_files]
    tuple val(cohort), path (GWAS_QC_noclump), \
        // val (SBayesRC_annot), path (SBayesRC_annot_path), \
        path(private_input_files_path)
    
    

    output:
    // $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
    tuple val(cohort), path ("${cohort}_LD/ldm.info"), path ("${cohort}_LD/snp.info"), path ("${cohort}_LD/block*.eigen.bin"),  emit: SBayes_LD
    
    
    script:
    """
    
    R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump}  ${private_input_files_path} 
    
    """
}
