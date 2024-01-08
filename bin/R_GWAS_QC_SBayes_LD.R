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
refblocks_GRCh37 = args[11]


LDdir       <- paste0(private_input_files_path, "/LD_ref/ukbEUR_HM3/")


(LDfile_path <- paste0(private_input_files_path, "/LD_ref/1000g_phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_EUR_annot_GRCh37_{CHR}"))

# Code
print("Step1: generate the LD block information and script")
# Output $outDir/ldm.info, $outDir/ld.sh, $outDir/snplist/*.snplist
SBayesRC::LDstep1(mafile=paste0(cohort,"_LOO_GWAS_QC_noclump.cojo"), 
                  genoPrefix=LDfile_path,
                  outDir=paste0(cohort,'_LD'), genoCHR='1-23', 
                  blockRef= refblocks_GRCh37, log2file=T)
gc() #free memory

print("Step2: generate each LD matrix for blocks")
#  Loop idx from 1 to NUM_BLOCK (589)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ld.sh, $outDir/snplist/$idx.snplist
#  Ouput $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
for(idx in (1:589)) {
    print(paste("Step2 loop", idx))
    SBayesRC::LDstep2(outDir=paste0(cohort,'_LD'), blockIndex=idx, log2file=T)
    gc() #free memory
}
gc() #free memory

print("obtain block ids from actual files:")
list <- list.files(paste0(cohort,'_LD'))
list <- grep('.ldm.full.info', list, value = TRUE)
(ids_list <- unique(readr::parse_number(list)))


print("Step3: eigen decomposition for each LD block")
#  Loop idx from 1 to NUM_BLOCK (589)
#  Submit multiple jobs on your cluster / clouds instead of for loop
#  Input depends on $outDir/ldm.info, $outDir/b$idx.ldm.full.info, $outDir/b$idx.ldm.full.bin
#  Output $outDir/block$block.eigen.bin, $outDir/block$block.eigen.bin.log
# export OMP_NUM_THREADS=$threads  # parallel computing supported in this step
for(idx in ids_list) {
    print(paste("Step3 loop", idx))
    SBayesRC::LDstep3(outDir=paste0(cohort,'_LD'), blockIndex=idx, log2file=T)
    gc() #free memory
}
gc() #free memory

print("Step4: merge LD information")
SBayesRC::LDstep4(outDir=paste0(cohort,'_LD'), log2file=T)

# Step5: clean if necessary
# Essential for analysis: $outDir/ldm.info, $outDir/snp.info, $outDir/block*.eigen.bin 
# Other files could be removed
# Note: before removing, check all blocks were finished.



