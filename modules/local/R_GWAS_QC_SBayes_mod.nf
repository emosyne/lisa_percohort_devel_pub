process R_GWAS_QC_SBayes {
    // debug true
    container 'emosyne/r_docker:2.01'
    label 'process_high_long'
    tag "${cohort}_${SBayesRC_annot}"
    cache "lenient"
    


    input: 
    // [xs234, xs234_GWAS_QC_noclump.gz, /home/osimoe/emanuele_project/private_input_files, annot_baseline2_2_with_continuous_enhancers, /home/osimoe/emanuele_project/private_input_files/SBayes_annots/annot_baseline2_2_with_continuous_enhancers.txt.gz]
    tuple val(cohort), path (GWAS_QC_noclump), path(private_input_files_path), val(SBayesRC_annot), path(SBayesRC_annot_file)
    
    

    output:
    // $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
    tuple val(cohort), path ("*"),  emit: SBayes
    
    
    script:
    """
    
    R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump} ${private_input_files_path} ${SBayesRC_annot} ${SBayesRC_annot_file}
    
    """
}
