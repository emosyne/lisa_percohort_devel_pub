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

#R_GWAS_QC_SBayes.R ${cohort} ${cohort_LD} ${GWAS_QC_noclump_cojo} ${SBayesRC_annot} ${SBayesRC_annot_path} ${private_input_files_path} 
(cohort = args[8])
cohort_LD_path = args[9]
GWAS_QC_noclump_cojo = fread(args[10])
SBayesRC_annot = args[11]
SBayesRC_annot_path = args[12]
private_input_files_path = args[13]




head(GWAS_QC_noclump_cojo)
  


print("SBayesRC: main function for SBayesRC")
SBayesRC::sbayesrc(mafile=paste0(cohort,'_LOO_GWAS_QC_noclump_imp.ma'), LDdir=cohort_LD_path,
                  outPrefix=paste0(cohort,'_',SBayesRC_annot,'_sbrc'),
                  annot=SBayesRC_annot_path, 
                  log2file=F)

