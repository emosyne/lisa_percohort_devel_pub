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

#R_GWAS_QC_SBayes.R ${cohort} ${GWAS_QC_noclump} ${private_input_files_path} 
(cohort = args[8])
GWAS_QC_noclump = fread(args[9])
private_input_files_path = args[10]
SBayesRC_annot = args[11]
SBayesRC_annot_path = args[12]


#current GWAS format:
# (/gpfs/home2/osimoe/nf) [osimoe@int6 lisa_percohort_devel_pub]$ zcat  /gpfs/home2/osimoe/.nextflow/assets/emosyne/lisa_percohort_devel_pub/work/a4/b9989df138ae1893b5033a3b5a8e66/xs234_GWAS_QC_noclump.gz | head | column -t
# CHR  SNP         BP         A1  A2  FRQ_A_51419  FRQ_U_74993  INFO   OR       SE      P       ngt  Direction                                                                    HetISqt  HetDf  HetPVa   Nca    Nco    Neff
# 8    rs62513865  101592213  C   T   0.931        0.928        0.964  1.01005  0.0175  0.5677  0    --+-++-+--+-+-++-++-++++-+--+-++-----++--------++++++-+-+-+-++---++-+--++++  2.9      74     0.4075   51419  74993  56643.62
# 8    rs79643588  106973048  G   A   0.908        0.906        0.998  0.99134  0.0151  0.5656  0    +--+-+----+-+++++++++++-+-+-+++--+-----+----+++----++-+---++++++++++---+--+  23.0     74     0.04289  51419  74993  56643.62


### Format GWAS in cojo
# colnames(formatted.gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
LOO_GWAS_QC_noclump <- GWAS_QC_noclump %>%
    mutate(freq = (( FRQ_A_51419 * Nca) + (FRQ_U_74993 * Nco) ) /(Nca + Nco), FRQ_A_51419 = NULL, FRQ_U_74993 = NULL,
           N    = Nca + Nco, Nca = NULL, Nco = NULL,
           b    = log(OR), OR = NULL) %>%
    select(SNP,A1, A2, freq, b, se=SE, p=P, N)

head(LOO_GWAS_QC_noclump)
  


data.table::fwrite(x = LOO_GWAS_QC_noclump, file = paste0(cohort,"_LOO_GWAS_QC_noclump.cojo"))



## pipeline using UKBB-LD ##############################################
# 1.135.332 SNPs in common between GWAS summary and LD

(cohort_LD_path     <- paste0(private_input_files_path, "/LD_ref/ukbEUR_HM3/"))
(cohort_LD          <- paste0(cohort,'_UKBB-LD'))
  
print("Tidy: optional step, tidy summary data")
## "log2file=TRUE" means the messages will be redirected to a log file 
SBayesRC::tidy(mafile=paste0(cohort,"_LOO_GWAS_QC_noclump.cojo"), LDdir=cohort_LD_path,
               output=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_tidy.ma'), log2file=T)



## Best practice: read the log to check issues in your GWAS summary data.  
print("Impute: optional step if your summary data doesn't cover the SNP panel")
SBayesRC::impute(mafile=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_tidy.ma'), LDdir=cohort_LD_path,
                 output=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_imp.ma'), log2file=T)



print("SBayesRC: main function for SBayesRC  without annotation (for comparison)")
SBayesRC::sbayesrc(mafile=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_imp.ma'), LDdir=cohort_LD_path,
                 outPrefix=paste0(cohort_LD,'_sbrc_noAnnot'), log2file=T)

print("SBayesRC: main function for SBayesRC")
SBayesRC::sbayesrc(mafile=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_imp.ma'), LDdir=cohort_LD_path,
                  outPrefix=paste0(cohort_LD,'_',SBayesRC_annot,'_sbrc'),
                  annot=SBayesRC_annot_path, 
                  log2file=F)


# ## pipeline using PGC_xs234_LD ##############################################
# # 1.394.135 SNPs in common between GWAS summary and LD

# (cohort_LD_path     <- paste0(private_input_files_path, "/LD_ref/PGC_xs234_LD/"))
# (cohort_LD = paste0(cohort,'_PGC_xs234_LD'))

# print("Tidy: optional step, tidy summary data")
# ## "log2file=TRUE" means the messages will be redirected to a log file 
# SBayesRC::tidy(mafile=paste0(cohort,"_LOO_GWAS_QC_noclump.cojo"), LDdir=cohort_LD_path,
#                output=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_tidy.ma'), log2file=T)



# ## Best practice: read the log to check issues in your GWAS summary data.  
# print("Impute: optional step if your summary data doesn't cover the SNP panel")
# SBayesRC::impute(mafile=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_tidy.ma'), LDdir=cohort_LD_path,
#                  output=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_imp.ma'), log2file=T)




# print("SBayesRC: main function for SBayesRC  without annotation (for comparison)")
# SBayesRC::sbayesrc(mafile=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_imp.ma'), LDdir=cohort_LD_path,
#                  outPrefix=paste0(cohort_LD,'_sbrc_noAnnot'), log2file=T)

# print("SBayesRC: main function for SBayesRC")
# SBayesRC::sbayesrc(mafile=paste0(cohort_LD,'_LOO_GWAS_QC_noclump_imp.ma'), LDdir=cohort_LD_path,
#                   outPrefix=paste0(cohort_LD,'_',SBayesRC_annot,'_sbrc'),
#                   annot=SBayesRC_annot_path, 
#                   log2file=F)