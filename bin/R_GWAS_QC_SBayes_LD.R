#!/usr/bin/env Rscript
library(data.table)
library(Rcpp)
library(stringi)
library(BH)
library(RcppEigen)
library(tidyverse)
# library(GenomicRanges)
# library(rtracklayer)

#INPUT
args = commandArgs()

print(args)

#R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump}  ${private_input_files_path} ${SBayesRC_annot} ${SBayesRC_annot_path}
(cohort = args[8])
LOO_GWAS_QC_noclump = fread(args[9], select=c("CHR", "SNP", "BP", "A1", "A2", "FRQ_A_51419", "FRQ_U_74993", "OR", "SE", "P", "Nca", "Nco"))
private_input_files_path = args[10]
# SBayesRC_annot = args[11]
# SBayesRC_annot_path = args[12]

LDdir       <- paste0(private_input_files_path, "/LD_ref/ukbEUR_HM3/")
(LDfile_path <- paste0(private_input_files_path, "/LD_ref/1000g_phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_EUR_annot_GRCh37"))
refblocks_GRCh37 <- paste0(private_input_files_path, "/SBayes_annots/refblocks_GRCh37.txt")

#current GWAS format:
# (/gpfs/home2/osimoe/nf) [osimoe@int6 lisa_percohort_devel_pub]$ zcat  /gpfs/home2/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel_pub/work/a4/b9989df138ae1893b5033a3b5a8e66/xs234_GWAS_QC_noclump.gz | head | column -t
# CHR  SNP         BP         A1  A2  FRQ_A_51419  FRQ_U_74993  INFO   OR       SE      P       ngt  Direction                                                                    HetISqt  HetDf  HetPVa   Nca    Nco    Neff
# 8    rs62513865  101592213  C   T   0.931        0.928        0.964  1.01005  0.0175  0.5677  0    --+-++-+--+-+-++-++-++++-+--+-++-----++--------++++++-+-+-+-++---++-+--++++  2.9      74     0.4075   51419  74993  56643.62
# 8    rs79643588  106973048  G   A   0.908        0.906        0.998  0.99134  0.0151  0.5656  0    +--+-+----+-+++++++++++-+-+-+++--+-----+----+++----++-+---++++++++++---+--+  23.0     74     0.04289  51419  74993  56643.62


### Format GWAS in cojo
# colnames(formatted.gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
LOO_GWAS_QC_noclump <- LOO_GWAS_QC_noclump %>%
    mutate(freq = (( FRQ_A_51419 * Nca) + (FRQ_U_74993 * Nco) ) /(Nca + Nco), FRQ_A_51419 = NULL, FRQ_U_74993 = NULL,
           N    = Nca + Nco, Nca = NULL, Nco = NULL,
           b    = log(OR), OR = NULL) %>%
    select(SNP,A1, A2, freq, b, se=SE, p=P, N)

head(LOO_GWAS_QC_noclump)
  


data.table::fwrite(x = LOO_GWAS_QC_noclump, file = paste0(cohort,"_LOO_GWAS_QC_noclump.cojo"))



##############################################
# Code
print("Step1: generate the LD block information and script")
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
SBayesRC::LDstep1(mafile=paste0(cohort,"_LOO_GWAS_QC_noclump.cojo"), 
                  genoPrefix=LDfile_path,
                  outDir=paste0(cohort,'_LD'), genoCHR='', 
                  blockRef=refblocks_GRCh37, log2file=F)

print("Step2: generate each LD matrix for blocks")
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ld.sh, $outDir/snplist/$idx.snplist
#  Ouput $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
for(idx in (1:591)) {
    print(paste("Step2 loop", idx))
    SBayesRC::LDstep2(outDir=paste0(cohort,'_LD'), blockIndex=idx, log2file=F)
}


print("Step3: eigen decomposition for each LD block")
#  Loop idx from 1 to NUM_BLOCK (591)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ldm.info, $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
#  Output $outDir/block$block.eigen.bin, $outDir/block$block.eigen.bin.log
# export OMP_NUM_THREADS=$threads  # parallel computing supported in this step
for(idx in (1:591)) {
    print(paste("Step3 loop", idx))
    SBayesRC::LDstep3(outDir=paste0(cohort,'_LD'), blockIndex=idx, log2file=F)
    gc() #free memory
}

print("Step4: merge LD information")
SBayesRC::LDstep4(outDir=paste0(cohort,'_LD'), log2file=F)

# Step5: clean if necessary
# Essential for analysis: $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
# Other files could be removed
# Note: before removing, check all blocks were finished.



print("Tidy: optional step, tidy summary data")
## "log2file=TRUE" means the messages will be redirected to a log file 
SBayesRC::tidy(mafile=GWAS_QC_noclump_cojo, LDdir=cohort_LD_path,
               output=paste0(cohort,'_LOO_GWAS_QC_noclump_tidy.ma'), log2file=F)



## Best practice: read the log to check issues in your GWAS summary data.  
print("Impute: optional step if your summary data doesn't cover the SNP panel")
SBayesRC::impute(mafile=paste0(cohort,'_LOO_GWAS_QC_noclump_tidy.ma'), LDdir=cohort_LD_path,
                 output=paste0(cohort,'_LOO_GWAS_QC_noclump_imp.ma'), log2file=F)


print("SBayesRC: main function for SBayesRC  without annotation (for comparison)")
SBayesRC::sbayesrc(mafile=paste0(cohort,'_LOO_GWAS_QC_noclump_imp.ma'), LDdir=cohort_LD_path,
                 outPrefix=paste0(cohort,'_',SBayesRC_annot,'_sbrc_noAnnot'), log2file=F)