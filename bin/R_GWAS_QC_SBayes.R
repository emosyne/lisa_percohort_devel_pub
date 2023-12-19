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

orig_GWAS = args[8]
cohort_bimfile = args[9]
LOO_GWAS = args[10]
(cohort = args[11])





### Format GWAS in cojo
# colnames(formatted.gwas) = c("SNP", "A1", "A2", "freq", "b", "se", "p", "N")

head(orig_GWAS)
  
orig_GWAS$freq = ( (orig_GWAS[,FCAS] * orig_GWAS[,NCAS]) + (orig_GWAS[,FCON] * orig_GWAS[,NCON]) ) /(orig_GWAS[,NCAS] + orig_GWAS[,NCON])
orig_GWAS$N =  (orig_GWAS[,NCAS] + orig_GWAS[,NCON])

cojo_GWAS <- orig_GWAS %>% select(SNP, A1, A2, freq, b=BETA, se=SE, p=P, N)

data.table::fwrite(x = cojo_GWAS, file = "gwas.cojo")



# Tidy: optional step, tidy summary data
## "log2file=TRUE" means the messages will be redirected to a log file 

SBayesRC::tidy(mafile="gwas.cojo", LDdir='ukbEUR_HM3/', 
                  output='cojo_GWAS_tidy.ma', log2file=TRUE)



## Best practice: read the log to check issues in your GWAS summary data.  

# Impute: optional step if your summary data doesn't cover the SNP panel

SBayesRC::impute(mafile='cojo_GWAS_tidy.ma', LDdir='ukbEUR_HM3/', 
                  output='cojo_GWAS_imp.ma', log2file=TRUE)




# SBayesRC: main function for SBayesRC
```{bash }
#!/bin/bash
#PBS -lselect=1:ncpus=12:mem=250gb
#PBS -lwalltime=24:0:0
#PBS -N SBayesRC
source /rds/general/user/eosimo/home/.bashrc
module load anaconda3/personal
source activate r422

cd /rds/general/user/eosimo/home/lenhard_prs/SBayesRC

##############################################
# Variables: need to be fixed
ma_file="/rds/general/user/eosimo/home/lenhard_prs/SBayesRC/gwas.cojo"               # GWAS summary in COJO format (the only input)
ld_folder="/rds/general/user/eosimo/home/lenhard_prs/SBayesRC/ukbEUR_Imputed"        # LD reference 
annot="/rds/general/user/eosimo/home/lenhard_prs/SBayesRC/annot_baseline2_2_with_enhancers.txt"         # Functional annotation 
out_prefix="annot_baseline2_2_with_enhancers"   # Output prefix, e.g. "./test"
threads=12                       # Number of CPU cores

##############################################
# Code: usually don't need a change in this section
## Note: Flags were documented in the package, use ?function in R to lookup.
## We suggest to run those in multiple jobs (tasks)
export OMP_NUM_THREADS=$threads # Revise the threads

## Tidy: optional step, tidy summary data
### "log2file=TRUE" means the messages will be redirected to a log file 
#Rscript -e "SBayesRC::tidy(mafile='$ma_file', LDdir='$ld_folder', \
#                  output='${out_prefix}_tidy.ma', log2file=TRUE)"
### Best practice: read the log to check issues in your GWAS summary data.  

## Impute: optional step if your summary data doesn't cover the SNP panel
#Rscript -e "SBayesRC::impute(mafile='${out_prefix}_tidy.ma', LDdir='$ld_folder', \
#                  output='${out_prefix}_imp.ma', log2file=TRUE)"

# SBayesRC: main function for SBayesRC
Rscript -e "SBayesRC::sbayesrc(mafile='output_imp.ma', LDdir='$ld_folder', \
                  outPrefix='output/${out_prefix}_sbrc', annot='$annot', log2file=TRUE)"
# Alternative run, SBayesRC without annotation (similar to SBayesR, not recommended)
# Rscript -e "SBayesRC::sbayesrc(mafile='${out_prefix}_imp.ma', LDdir='$ld_folder', \
#                  outPrefix='${out_prefix}_sbrc_noAnnot', log2file=TRUE)"



```
