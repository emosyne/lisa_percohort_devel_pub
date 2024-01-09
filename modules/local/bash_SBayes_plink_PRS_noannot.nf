process bash_SBayes_plink_PRS_noannot {
    tag "${cohort}_${SBayesRC_annot}"
    // debug true
    label 'process_high'
    container 'emosyne/plink2:7.2_Alpha5.9'
    cache "lenient"

    input: 
    // [xs234, annot_baseline2_2_with_continuous_enhancers, xs234_UKBB-LD_annot_baseline2_2_with_continuous_enhancers_sbrc.txt, /gpfs/work5/0/pgcdac/DWFV2CJb8Piv_0116_pgc_data/scz/wave3/v1/xs234]
    tuple val(cohort), val(SBayesRC_annot), path(sbrc_PRS), path(cohort_dir)
    

    output:
    tuple val(cohort), val(SBayesRC_annot), path ("*.profile"),                 emit: SBayesRC_PRS

    // path("*.log")
    


    script:
    def mem_mb = (task.memory * 0.95).toMega()
    """ 
    bedfile="${cohort_dir}/imputed/hardcall_genotypes/*.bed"
    bimfile="${cohort_dir}/imputed/hardcall_genotypes/*.bim"
    bedfile2=`echo \$bedfile | sed 's/.bed//'`
    covariates="${cohort_dir}/prin_comp/*.mds"

    echo \$bedfile2

    plink  \\
       --bfile \$bedfile2 \\
       --score ${sbrc_PRS} 1 2 3 header sum center \\
       --out ${cohort}_${SBayesRC_annot}_SBayes_PRS \\
       --threads $task.cpus \\
       --memory $mem_mb
    """
}

