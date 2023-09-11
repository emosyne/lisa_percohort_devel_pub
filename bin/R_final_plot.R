#!/usr/bin/env Rscript

#disable all warnings
# options(warn=-1)

library(tidyverse)
library(data.table)
# require(rms)
# library(fst)
library(gridExtra)
library(grid)

#INPUT
args = commandArgs()


# R_final_plot.R $task.cpus ${cohort_ENHpart} ${cohort_fam} \
        # ${TS_ENH_GWAS_compartment_originalOR_summary} ${TS_ENH_GWAS_compartment_originalOR_best}\
        # ${TS_ENH_GWAS_compartment_OR_by_measure1_summary} ${TS_ENH_GWAS_compartment_OR_by_measure1_best}\
        # ${TS_ENH_GWAS_compartment_OR_by_measure2_summary} ${TS_ENH_GWAS_compartment_OR_by_measure2_best}\
        # ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        # ${merged_GWAS_summary} ${merged_GWAS_best}\
        # ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
        # ${original_LOO_GWAS_summary} ${original_LOO_GWAS_prsice} ${original_LOO_GWAS_best}\
        # ${modif_name_1} ${modif_name_2} ${CTthreshold}

print(args)

nthreads = as.numeric(args[8])
#set max CPU processes
setDTthreads(nthreads)
# threads_fst(nr_of_threads = round(nthreads/3*2))

(ENH_list = args[9])
(diagnosis = fread(args[10], header=F, col.names = c("FID", "IID", "IIDf", "IIDm", "sex", "dx" )) %>%
    dplyr::select("FID", "IID", "dx"))

TS_ENH_GWAS_compartment_originalOR_summary = args[11]
TS_ENH_GWAS_compartment_originalOR_best = 
  fread(args[12], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(TS_ENH_GWAS_compartment_originalOR_best_PRS = PRS)

TS_ENH_GWAS_compartment_OR_by_measure1_summary = args[13]
TS_ENH_GWAS_compartment_OR_by_measure1_best = 
  fread(args[14], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(TS_ENH_GWAS_compartment_OR_by_measure1_best_PRS = PRS)

TS_ENH_GWAS_compartment_OR_by_measure2_summary = args[15]
TS_ENH_GWAS_compartment_OR_by_measure2_best = 
  fread(args[16], select=c("FID", "IID", "PRS")) %>% 
  dplyr::rename(TS_ENH_GWAS_compartment_OR_by_measure2_best_PRS = PRS)
# ${residual_GWAS_compartment_summary} ${residual_GWAS_compartment_best}\
        # ${merged_GWAS_summary} ${merged_GWAS_best}\
        # ${TS_ENH_GWAS_compartment_originalOR_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure1_prsice} ${TS_ENH_GWAS_compartment_OR_by_measure2_prsice} 
        # ${residual_GWAS_compartment_prsice} ${merged_GWAS_prsice}  \
        # ${original_LOO_GWAS_summary} ${original_LOO_GWAS_prsice} ${original_LOO_GWAS_best}\
        # ${modif_name_1} ${modif_name_2} ${CTthreshold}
residual_GWAS_compartment_summary = args[17]
residual_GWAS_compartment_best = 
  fread(args[18], select=c("FID", "IID", "PRS"))  %>% 
  dplyr::rename(residual_GWAS_compartment_best_PRS = PRS)


TS_ENH_GWAS_compartment_originalOR_prsice = args[21]
TS_ENH_GWAS_compartment_OR_by_measure1_prsice = args[22]
TS_ENH_GWAS_compartment_OR_by_measure2_prsice = args[23]
residual_GWAS_compartment_prsice = args[24]
merged_GWAS_prsice = args[25]

original_GWAS_summary = args[26]
original_GWAS_prsice = args[27]
(original_GWAS_best = fread(args[28], select=c("FID", "IID", "PRS"))  %>% 
  dplyr::rename(original_GWAS_best_PRS = PRS))

modif_name_1 = args[29]
modif_name_2 = args[30]
threshold = args[31]


#set input variables
number_quantiles = 3
condition_name = "schizophrenia"
# pop_prev = population prevalence
pop_prev = 0.01


#OUTPUT_prefix
OUTPUT_prefix = paste0(threshold,"/")
if (!dir.exists(file.path(paste0(threshold,"")))) {dir.create(file.path(paste0(threshold,"")))}



############### start main ###############



#SUMMARY TABLE
#import thresholds and SNP N for each summary
(summary_table = 
   rbind(
     "EPWAS_originalOR_summary" = 
       data.frame(fread(EPWAS_originalOR_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
     "EPWAS_OR_by_measure1_summary" = 
       data.frame(fread(EPWAS_OR_by_measure1_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
     "EPWAS_OR_by_measure2_summary" = 
       data.frame(fread(EPWAS_OR_by_measure2_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
     "residual_GWAS_compartment_summary"=
       data.frame(fread(residual_GWAS_compartment_summary, select=c("Threshold", "PRS.R2","Num_SNP"))),#"PRS.R2.adj",
     "original_GWAS_summary"=
       data.frame(fread(original_GWAS_summary, select=c("Threshold", "PRS.R2","Num_SNP")))#"PRS.R2.adj",
   ) %>%  rownames_to_column(var = "compartment"))



#BEST TABLE
#create total PRS score
(BEST_PRS_score_per_UKBB_participant <- EPWAS_originalOR_best %>%
    left_join(EPWAS_OR_by_measure1_best) %>%
    left_join(EPWAS_OR_by_measure2_best) %>%
    left_join(residual_GWAS_compartment_best) %>%
    # left_join(merged_GWAS_best) %>%
    left_join(original_GWAS_best) %>%
    left_join(diagnosis) %>%
    mutate(dx=factor(dx), IID=factor(IID)) %>%
    select(-FID) %>%
    remove_missing() #%>% head(n=50000)
)


(scaled_BEST_PRS_score_per_UKBB_participant <- BEST_PRS_score_per_UKBB_participant)
# fwrite(BEST_PRS_score_per_UKBB_participant,"BEST_PRS_score_per_UKBB_participant.txt")
scaled_BEST_PRS_score_per_UKBB_participant[,c(2:6)] <-  data.frame(scale(BEST_PRS_score_per_UKBB_participant[,c(2:6)], center = T, scale = T))+10
# head(scaled_BEST_PRS_score_per_UKBB_participant)
# fwrite(scaled_BEST_PRS_score_per_UKBB_participant,"scaled_BEST_PRS_score_per_UKBB_participant.txt")

(scaled_BEST_PRS_score_per_UKBB_participant<-
    scaled_BEST_PRS_score_per_UKBB_participant %>% 
    # mutate(weight_total_PRS_best = 
    #          (EPWAS_originalOR_best_PRS * summary_table[summary_table$compartment=="EPWAS_originalOR_summary",]$Num_SNP  /summary_table[summary_table$compartment=="merged_GWAS_summary",]$Num_SNP) + 
    #          (residual_GWAS_compartment_best_PRS     * summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP      /summary_table[summary_table$compartment=="merged_GWAS_summary",]$Num_SNP)) %>%
    remove_missing() %>% 
    #generate quantiles
    mutate(original_GWAS_q = 
             factor(ntile(original_GWAS_best_PRS, n = number_quantiles))) %>% 
    # mutate(merged_GWAS_q = 
    #          factor(ntile(merged_GWAS_best_PRS, n = number_quantiles))) %>% 
    mutate(residual_GWAS_compartment_q = 
             factor(ntile(residual_GWAS_compartment_best_PRS, n = number_quantiles))) %>% 
    mutate(TS_ENH_compartment_originalOR_q = 
             factor(ntile(EPWAS_originalOR_best_PRS, number_quantiles))) %>% 
    # mutate(weight_total_q = 
    #          factor(ntile(weight_total_PRS_best, number_quantiles)))%>% 
    mutate(TS_ENH_compartment_OR_by_measure1_q = 
             factor(ntile(EPWAS_OR_by_measure1_best_PRS, number_quantiles)))%>% 
    mutate(TS_ENH_compartment_OR_by_measure2_q = 
             factor(ntile(EPWAS_OR_by_measure2_best_PRS, number_quantiles)))
)





# str(scaled_BEST_PRS_score_per_UKBB_participant)
# scaled_BEST_PRS_score_per_UKBB_participant[rowSums(is.na(scaled_BEST_PRS_score_per_UKBB_participant)) > 0,]


#Write output to file
sink(paste0(OUTPUT_prefix, cohort_ENHlist_thresh, "_", Sys.Date(),"_logfile.log"), split = T)
print("measures of model fit:")
#### measures of model fit 
## https://onlinelibrary.wiley.com/doi/full/10.1002/gepi.21614
# nt = total number of the sample
(nt = NROW(diagnosis))
# ncase = number of cases
(ncase = NROW(diagnosis[diagnosis$dx==2,]))
# ncont = number of controls
(ncont = NROW(diagnosis[diagnosis$dx==1,]))
# case_prev_in_sample = proportion of cases in the case-control samples
(case_prev_in_sample = ncase/nt)
# thd = the threshold on the normal distribution which truncates the proportion of disease prevalence
(thd = -qnorm(pop_prev,0,1))
(zv = dnorm(thd)) #z (normal density)
(mv = zv/pop_prev) #mean liability for case
(mv2 = -mv*pop_prev/(1-pop_prev)) #mean liability for controls

# #R2 on the observed scale
# (theta = mv*(case_prev_in_sample-pop_prev)/(1-pop_prev)*(mv*(case_prev_in_sample-pop_prev)/(1-pop_prev)-thd)) #θ in equation 15
# # theta = mv * (P - K)/(1 - K) * (mv * (P - K)/(1 - K) - thd)
# (cv = pop_prev*(1-pop_prev)/zv^2*pop_prev*(1-pop_prev)/(case_prev_in_sample*(1-case_prev_in_sample))) #C inequation 15
# # cv = K * (1 - K)/zv^2 * K * (1 - K)/(P * (1 - P))
# #Choi
# (e = 1 - (case_prev_in_sample^(2 * case_prev_in_sample)) * ((1 - case_prev_in_sample)^(2 * (1 - case_prev_in_sample))))
# (top = cv * e)
# (bottom = cv * e * theta)
# ChoiMe <- function(x,...){
#   top * x / (1 + bottom * x)
# }



# Start writing to an output file
(CoD_per_SNP = data.frame())

scaled_BEST_PRS_score_per_UKBB_participant$dx <- ifelse(scaled_BEST_PRS_score_per_UKBB_participant$dx=="1", 0, 1)

# ## original_GWAS
##bootstrap to obtain R2 variance
(logit = rms::lrm(factor(dx) ~ original_GWAS_best_PRS, scaled_BEST_PRS_score_per_UKBB_participant)$stats[10])

linear = lm(dx ~ original_GWAS_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant)
summary(linear)[[8]]

(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ original_GWAS_best_PRS, data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2       = boot_lm$t0)
(linear_R2_SE    = sd(boot_lm$t))
(linear_R2_95CI  = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))
(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="0",    
               NSNP=summary_table[summary_table$compartment=="original_GWAS_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))



# ## original_GWAS + squared original GWAS
# A more streamlined approach can be adopted by consolidating the various predictors (any_cov1, any_cov2, …, any_covN) into a single predictor, following this procedure in R: 
#   
#   R 
# mod <- lm(y ~PRS + any_cov1 + any_cov2 + … + any_covN) 
# merged_predictor <- cbind(any_cov1, any_cov2, … , any_covN) %*% mod$coefficients[3:(2+N)] 
#logistic model
linear = lm(dx ~ original_GWAS_best_PRS +
              I(original_GWAS_best_PRS^2), data = scaled_BEST_PRS_score_per_UKBB_participant)
summary(linear)[[8]]

merged_predictor_original_plus_squared_original <- 
  cbind(
    scaled_BEST_PRS_score_per_UKBB_participant$original_GWAS_best_PRS, 
    I(scaled_BEST_PRS_score_per_UKBB_participant$original_GWAS_best_PRS^2)) %*% 
  linear$coefficients[2:(3)] 

summary(
  lm(scaled_BEST_PRS_score_per_UKBB_participant$dx ~ merged_predictor_original_plus_squared_original)
)[[8]]

(boot_lm <- 
    boot(
      scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
        summary(lm(dx ~ original_GWAS_best_PRS +
                     I(original_GWAS_best_PRS^2), data[indices,]))[[8]], 
      R=500, ncpus = nthreads
    )
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="0b",    
               NSNP=summary_table[summary_table$compartment=="original_GWAS_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))





##residual
##bootstrap to obtain R2 variance
(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ residual_GWAS_compartment_best_PRS , data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="1",    
               NSNP=summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))



# TS ENH original OR
##bootstrap to obtain R2 variance
(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ EPWAS_originalOR_best_PRS, data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="2",    
               NSNP=summary_table[summary_table$compartment=="EPWAS_originalOR_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))



# TS ENH  OR by measure 1
##bootstrap to obtain R2 variance
(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ EPWAS_OR_by_measure1_best_PRS, data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="2b",    
               NSNP=summary_table[summary_table$compartment=="EPWAS_OR_by_measure1_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))

# TS ENH  OR by measure 2
##bootstrap to obtain R2 variance
(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ EPWAS_OR_by_measure2_best_PRS, data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="2c",    
               NSNP=summary_table[summary_table$compartment=="EPWAS_OR_by_measure2_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))



##simple additive model
linear = lm(dx ~ residual_GWAS_compartment_best_PRS + EPWAS_originalOR_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant)
summary(linear)[[8]]

merged_predictor_residual_plus_EPWAS_originalOR <- 
  cbind(
    scaled_BEST_PRS_score_per_UKBB_participant$residual_GWAS_compartment_best_PRS, 
    scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS) %*% 
  linear$coefficients[2:(3)] 
summary(
  lm(scaled_BEST_PRS_score_per_UKBB_participant$dx ~ merged_predictor_residual_plus_EPWAS_originalOR)
)[[8]]

##bootstrap to obtain R2 variance
(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ residual_GWAS_compartment_best_PRS + EPWAS_originalOR_best_PRS, data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="3",    
               NSNP=summary_table[summary_table$compartment=="EPWAS_OR_by_measure2_summary",]$Num_SNP +
                 summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))



## full factorial design 
linear = lm(dx ~ residual_GWAS_compartment_best_PRS * EPWAS_originalOR_best_PRS, data = scaled_BEST_PRS_score_per_UKBB_participant)
summary(linear)[[8]]

merged_predictor_residual_plus_EPWAS_originalOR_interactions <- 
  cbind(
    scaled_BEST_PRS_score_per_UKBB_participant$residual_GWAS_compartment_best_PRS, 
    scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS,
    scaled_BEST_PRS_score_per_UKBB_participant$residual_GWAS_compartment_best_PRS*scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS) %*% 
  linear$coefficients[2:4] 
summary(
  lm(scaled_BEST_PRS_score_per_UKBB_participant$dx ~ merged_predictor_residual_plus_EPWAS_originalOR_interactions)
)[[8]]

##bootstrap to obtain R2 variance
(boot_lm <- boot(
  scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
    summary(lm(dx ~ residual_GWAS_compartment_best_PRS * EPWAS_originalOR_best_PRS, data[indices,]))[[8]], 
  R=500, ncpus = nthreads
)
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))

(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="3b",    
               NSNP=summary_table[summary_table$compartment=="EPWAS_OR_by_measure2_summary",]$Num_SNP +
                 summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))



# full factorial design including interactions and non-linear interactions
linear = lm(dx ~ residual_GWAS_compartment_best_PRS*EPWAS_originalOR_best_PRS +
              I(residual_GWAS_compartment_best_PRS^2) + I(EPWAS_originalOR_best_PRS^2), 
            data = scaled_BEST_PRS_score_per_UKBB_participant)
summary(linear)[[8]]

merged_predictor_residual_plus_EPWAS_originalOR_interactions_quadraticterms <- 
  cbind(
    scaled_BEST_PRS_score_per_UKBB_participant$residual_GWAS_compartment_best_PRS, 
    scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS,
    I(scaled_BEST_PRS_score_per_UKBB_participant$residual_GWAS_compartment_best_PRS^2),
    I(scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS^2),
    scaled_BEST_PRS_score_per_UKBB_participant$residual_GWAS_compartment_best_PRS*scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS) %*% 
  linear$coefficients[2:6] 
summary(
  lm(scaled_BEST_PRS_score_per_UKBB_participant$dx ~ merged_predictor_residual_plus_EPWAS_originalOR_interactions_quadraticterms)
)[[8]]

##bootstrap to obtain R2 variance
(boot_lm <- 
    boot(
      scaled_BEST_PRS_score_per_UKBB_participant,  function(data,indices)
        summary(lm(dx ~ residual_GWAS_compartment_best_PRS*EPWAS_originalOR_best_PRS +
                     I(residual_GWAS_compartment_best_PRS^2) + I(EPWAS_originalOR_best_PRS^2), 
                   data[indices,]))[[8]],
      R=500, ncpus = nthreads
    )
)
(linear_R2     = boot_lm$t0)
(linear_R2_SE  = sd(boot_lm$t))
(CI_95         = data.frame(t(matrix(quantile(boot_lm$t,c(0.025,0.975)))) ))

(R2_liability    = data.frame(r2redux::cc_trf(R2 = linear_R2, se = linear_R2_SE, K = pop_prev, P = case_prev_in_sample)))
#A 95% confidence interval (CI) is twice the standard error (also called margin of error) plus or minus the mean.
(R2_liability_CI = cbind(R2_liability$R2l - 1.96*R2_liability$sel,
                         R2_liability$R2l + 1.96*R2_liability$sel))
(results_table = rbind(setNames(cbind(linear_R2,linear_R2_95CI, t="r2"), nm = c("R2","R2_LCI","R2_UCI","t")),
                       setNames(data.frame(cbind(R2_liability$R2l, R2_liability_CI,t="Lee")), nm = c("R2","R2_LCI","R2_UCI","t")))
)
(info = tibble(comp="3c",    
               NSNP=summary_table[summary_table$compartment=="EPWAS_OR_by_measure2_summary",]$Num_SNP +
                 summary_table[summary_table$compartment=="residual_GWAS_compartment_summary",]$Num_SNP, 
               CoD_per_SNP=NA) %>% 
    #double row
    slice(rep(1:n(), each = 2)))
## add to df
(CoD_per_SNP = rbind(
  CoD_per_SNP,
  cbind(info,results_table)
))

# CoD PER SNP df ########

colnames(CoD_per_SNP)=c("partition","Num_SNP","CoD_per_SNP","R2","LCL","UCL","R2type")
CoD_per_SNP[c(2:6)]<-sapply(CoD_per_SNP[c(2:6)],as.numeric)
CoD_per_SNP$CoD_per_SNP = (CoD_per_SNP$R2 / CoD_per_SNP$Num_SNP)*10^7
(CoD_per_SNP = as_tibble(CoD_per_SNP))





######### Paiwise comparisons with r2redux: ##########
comparisons <- data.frame(matrix(ncol = 6, nrow = 0))
#original vs original_plus_squared_original
compar_1=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  scaled_BEST_PRS_score_per_UKBB_participant$original_GWAS_best_PRS,
  merged_predictor_original_plus_squared_original
)
(res_2 = r2redux::r2_diff(dat = compar_1, v1=c(2), v2=c(1), nv=nrow(compar_1)))
comparisons = rbind(
  comparisons,
  cbind(
    "original vs original_plus_squared_original",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
comparisons
# original vs residual+EP
compar_2=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  scaled_BEST_PRS_score_per_UKBB_participant$original_GWAS_best_PRS,
  merged_predictor_residual_plus_EPWAS_originalOR
)
(res_2 = r2redux::r2_diff(dat = compar_2, v1=c(2), v2=c(1), nv=nrow(compar_2)))
comparisons = rbind(
  comparisons,
  cbind(
    "original vs residual + EP",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
# original vs residual+EP + residualxEP
compar_3=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  scaled_BEST_PRS_score_per_UKBB_participant$original_GWAS_best_PRS,
  merged_predictor_residual_plus_EPWAS_originalOR_interactions
)
(res_2 = r2redux::r2_diff(dat = compar_3, v1=c(2), v2=c(1), nv=nrow(compar_3)))
comparisons = rbind(
  comparisons,
  cbind(
    "original vs residual + EP+ residualxEP",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
# original vs residual+EP + interactions+quadraticterms
compar_4=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  scaled_BEST_PRS_score_per_UKBB_participant$original_GWAS_best_PRS,
  merged_predictor_residual_plus_EPWAS_originalOR_interactions_quadraticterms
)
(res_2 = r2redux::r2_diff(dat = compar_4, v1=c(2), v2=c(1), nv=nrow(compar_4)))
comparisons = rbind(
  comparisons,
  cbind(
    "original vs residual+EP + interactions+quadraticterms",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
# original_plus_squared_original vs residual+EP + interactions+quadraticterms
compar_5=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  merged_predictor_original_plus_squared_original,
  merged_predictor_residual_plus_EPWAS_originalOR_interactions_quadraticterms
)
(res_2 = r2redux::r2_diff(dat = compar_5, v1=c(2), v2=c(1), nv=nrow(compar_5)))
comparisons = rbind(
  comparisons,
  cbind(
    "original_plus_squared_original vs residual+EP + interactions+quadraticterms",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
comparisons
# original enhancer partition vs ES x OR enhancer partition
compar_6=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS,
  scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_OR_by_measure1_best_PRS
)
(res_2 = r2redux::r2_diff(dat = compar_6, v1=c(2), v2=c(1), nv=nrow(compar_6)))
comparisons = rbind(
  comparisons,
  cbind(
    "original enhancer partition vs ES x OR enhancer partition",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
comparisons
# original enhancer partition vs TS tpm x OR enhancer partition
compar_7=cbind(
  scaled_BEST_PRS_score_per_UKBB_participant$dx,
  scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_originalOR_best_PRS,
  scaled_BEST_PRS_score_per_UKBB_participant$EPWAS_OR_by_measure2_best_PRS
)
(res_2 = r2redux::r2_diff(dat = compar_7, v1=c(2), v2=c(1), nv=nrow(compar_7)))
comparisons = rbind(
  comparisons,
  cbind(
    "original enhancer partition vs TS tpm x OR enhancer partition",
    r2redux::cc_trf(R2 = res_2$mean_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$lower_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    r2redux::cc_trf(R2 = res_2$upper_diff, se = res_2$var_diff, K = pop_prev, P = case_prev_in_sample)$R2l,
    res_2$r2_based_p,
    res_2$var_diff
  )
)
comparisons


colnames(comparisons) <- c("comparison","mean_diff","lower_diff","upper_diff","r2_based_p","var_diff")
comparisons[c(2:6)] <-  sapply(comparisons[c(2:6)], as.numeric)
comparisons
knitr::kable(comparisons, format = "latex", 
             caption = paste0("r2redux r2diff comparisons between variables for ",cohort_ENHlist_thresh, " ", Sys.Date()), 
             digits = 6)
comparisons[,sapply(comparisons, is.numeric)] <-round(comparisons[,sapply(comparisons, is.numeric)],3)

# r2redux::cc_trf(R2 = res_2$upper_diff, se =  res_2$var_diff, K = pop_prev, P = case_prev_in_sample)



addline_format <- function(x,...){
  gsub('_',' ',x)
}

#create DF for plotting
(df_plot<- data.frame(
  partition=c(factor(        c("0","0b","1","2","2b","2c","3","3b","3c"#,"3d"
  ), ordered = T)),
  partition_name= factor(x = c("0","0b","1","2","2b","2c","3","3b","3c"#,"3d"
  ), labels = c("Original GWAS PRS",
                "Original GWAS + quadratic Original PRS",
                "Residual partition PRS", 
                paste0(ENH_list," partition PRS Original OR"),
                paste0(ENH_list," partition PRS OR \u00D7 ",modif_name_1),
                paste0(ENH_list," partition PRS OR \u00D7 ",modif_name_2),
                paste0("Residual + ",ENH_list," partition PRS"),
                paste0("Residual \u00D7 ",ENH_list," partition PRS"),
                paste0("Residual \u00D7 ",ENH_list," partition PRS + quadratic terms")
                #paste0("Residual \u00D7 ",EnhSNP_to_residual_SNP_ratio[1],ENH_list," partition PRS + quadratic terms")
  ), 
  ordered = T)
) %>% 
    left_join(CoD_per_SNP, by="partition", multiple = "all") %>% 
    mutate(xlabel=factor(stringr::str_wrap(paste0(addline_format(partition_name), " (SNP N=", Num_SNP,")"),
                                           width=30))) %>% 
    #remove Nagel
    dplyr::filter(R2type=="Lee")
  
)

knitr::kable(df_plot, format = "latex", caption = "df_plot")

sink()

## FIGURE 1, COD AND COD PER SNP FOR ORIGINAL, ENHANCER AND RESIDUAL PARTITIONS ###
# pos <- position_jitter(width = 0, height = 0.1, seed = 2345)
# pos = position_dodge(width = 0.5)
pos = position_identity()
(p1 <- ggplot(data = df_plot[df_plot$partition == "0" |
                               df_plot$partition == "1" |
                               df_plot$partition == "2",], 
              aes(y=reorder(xlabel, desc(partition)), 
                  x=R2*100, xmin = LCL*100, xmax=UCL*100, 
                  colour=R2type,
                  label=round(R2*100,2))
) +  
    geom_pointrange(position=pos, size=1, lwd=1) + 
    # geom_linerange(position=pos, 
    #                aes(xmin = 0, xmax=R2*100 )) +
    scale_color_manual(values=MetBrewer::met.brewer("Johnson", 2),
                       name = expression(paste(R^2)))+ #labels=c(a=expression(paste(Delta^2))
    ggrepel::geom_label_repel(nudge_y = 0.1, nudge_x = 0.1, size = rel(4),  show.legend = F,min.segment.length = 0,
                              point.padding = NA, box.padding = 0.5)+ #nudge_y = -0.2, 
    scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .15))) + 
    scale_y_discrete(expand = expansion(add = .7)) + #.6 is default - add a little spacing to the sides
    xlab("") +  ylab("HCM diagnosis ~")+
    theme_bw() +
    theme(
      legend.position = "none",  
      text=element_text(lineheight = 0.8, angle = 0, size = 21, color = "gray8"),
      # axis.text.y = element_text(lineheight = 0.8, angle = 0, size = rel(2), color = "gray8"),
      # axis.title.y = element_text(angle = 90, size = rel(2),
      #                             margin = margin(t = 0, r = 20, b = 0, l = 0), color = "gray8"),
      # panel.grid.major.y = element_blank()
    ))


(f1<-arrangeGrob(textGrob("A)", just = "left",
                          gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                 textGrob(paste("Total CoD %, and 95% CI"), 
                          gp = gpar(fontsize = 18, fontface = "bold", col="darkgreen")), 
                 #textGrob("diagnosis ~ PRS, probit link function \nProportion of the total variance explained by the genetic factor on the liability scale, \ncorrected for ascertainment, as per Lee et al 2012", gp = gpar(fontsize = 10)), 
                 p1,
                 layout_matrix=rbind(c(1,2),
                                     c(3,3)),
                 widths = c(0.1, 1), heights = c(0.1, 1)))


(p2 <-ggplot(data = df_plot[df_plot$partition == "0" |
                              df_plot$partition == "1" |
                              df_plot$partition == "2",], 
             aes(
               y=reorder(xlabel, desc(partition)),
               x=CoD_per_SNP, xmax=CoD_per_SNP, xmin=0,
               colour=R2type, group=R2type,
               label=round(CoD_per_SNP,2))) +  
    geom_point(position=pos, size=3) +     
    #geom_linerange(position=pos) +
    scale_color_manual(values=MetBrewer::met.brewer("Johnson", 2),
                       name = expression(paste(R^2)))+
    ggrepel::geom_label_repel(position = pos, size = rel(4),  show.legend = F,min.segment.length = 0,
                              point.padding = NA, box.padding = 0.5)+
    scale_x_continuous(limits = c(0, max_cod_per_snp), expand = expansion(mult = c(0, .1))) + 
    scale_y_discrete(expand = expansion(add = .7)) + #.6 is default - add a little spacing to the sides
    xlab("") +  ylab("")+
    theme_bw() +
    theme(
      legend.position = "none",  
      axis.text.y = element_blank(),#element_text(lineheight = 0.8, angle = 0, size = rel(1.3)),
      text=element_text(lineheight = 0.8, angle = 0, size = 21, color = "gray8"),
      #plot.margin = margin(t = 0, r = 1, b = 1, l = 0.3, "cm"),
      # panel.grid.major.y = element_blank()
    ))


(f2<-arrangeGrob(textGrob("B)", just = "left",
                          gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
                 textGrob(bquote("Coefficients of determination per SNP \u00D7"~10^7), 
                          gp = gpar(fontsize = 18, fontface = "bold", col="darkgreen")), 
                 #textGrob("diagnosis ~ PRS, probit link function \nProportion of the total variance explained by the genetic factor on the liability scale, \ncorrected for ascertainment, as per Lee et al 2012", gp = gpar(fontsize = 10)), 
                 p2,
                 layout_matrix=rbind(c(1,2),
                                     c(3,3)),
                 widths = c(0.1, 1), heights = c(0.1, 1)))


fig1grob<- arrangeGrob(
  textGrob(paste("Coefficients of determination for the main three partitions: original, enhancer and residual"), 
           gp = gpar(fontsize = 22, fontface = "bold", col="darkgreen")), 
  f1, f2, 
  layout_matrix=rbind(c(1,1),
                      c(2,3)),
  widths = c(1, 0.5), heights = c(0.1,1)
)

ggsave(
  filename = paste0(OUTPUT_prefix, "fig1_", cohort_ENHlist_thresh, "_", Sys.Date(),"_CoD_main_partitions.pdf"), 
  fig1grob,  
  width = 17, height = 5)





### FIGURE 2, COD FOR ENH PARTITIONS ####
(plot2= df_plot[grepl(pattern = "^2", perl = T, x=df_plot$partition) ,] %>% droplevels())
(comparisons )

fig2 <-
  ggplot(data = plot2,
         aes(
           x=reorder(xlabel, desc(partition)), 
           y=R2*100, #xmin = plot2$LCL*100, xmax=plot2$UCL*100, 
           colour=R2type,group=R2type
           # label=round(R2*100,2)
         )) +
  geom_point(position = position_dodge(width = 0.5), size=3.5) +
  geom_errorbar(aes(ymin = LCL*100,
                    ymax = UCL*100),
                position = position_dodge(width = 0.5),
                width = 0, linewidth=1) +
  ## first comparison: original OR ENH partition vs OR * ES
  geom_bracket(data = plot2,
               xmin = reorder(plot2$xlabel, desc(plot2$partition))[[1]],
               xmax = reorder(plot2$xlabel, desc(plot2$partition))[[2]],
               y.position = max(plot2$UCL*100)*1.1,
               label =paste0("MD=",comparisons[6,2]*100,", p=",comparisons[6,5]),
               # tip.length = c(0.02, 0.02),
               coord.flip = TRUE,
               # vjust = -1
  ) +
  ## second comparison: original OR ENH partition vs OR * tpm
  geom_bracket(data = plot2,
               xmin = reorder(plot2$xlabel, desc(plot2$partition))[[1]],
               xmax = reorder(plot2$xlabel, desc(plot2$partition))[[3]],
               y.position = max(plot2$UCL*100)*1.2,
               label = paste0("MD=",comparisons[7,2]*100,", p=",comparisons[7,5]),
               # tip.length = c(0.02, 0.02),
               coord.flip = TRUE,
               # vjust = -1
  ) +
  scale_color_manual(values=MetBrewer::met.brewer("Johnson", 2),
                     expression(paste(R^2)))+
  ggrepel::geom_label_repel(
    data = plot2,
    aes(
      x=reorder(xlabel, desc(partition)), 
      y=R2*100,
      colour=R2type,group=R2type,
      label=round(R2*100,2)
    ), nudge_y = 0.1, nudge_x = 0.1, size = 4,  show.legend = F,#min.segment.length = 0.5,
    point.padding = NA, box.padding = 0.5) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .15)))  +
  scale_x_discrete(expand = expansion(add = .7)) + #.6 is default - add a little spacing to the sides
  ylab("") +  xlab(paste0("HCM diagnosis ~"))+theme_bw() +
  theme(
    legend.position = "none",  
    text=element_text(lineheight = 0.8, angle = 0, size = 21, color = "gray8")
  )+
  coord_flip()


fig2_grob = arrangeGrob(
  textGrob(paste("CoDs % and 95% CIs for the three enhancer partitions:\nOriginal OR, enhanced by ES, enhanced by expression"), gp = gpar(fontsize = 22, fontface = "bold", col="darkblue")), 
  fig2, 
  ncol=1, heights = c(0.2,1)
)

ggsave(
  filename = paste0(OUTPUT_prefix, "fig2_", cohort_ENHlist_thresh, "_", Sys.Date(),"_CoD_enh_partitions.pdf"), 
  fig2_grob,  
  width = 10, height = 4, device = "pdf", scale = 1.5)



## FIGURE 3, COD FOR original,  PARTITIONS ####

(fig3 <- ggplot(data = df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,], 
                aes(
                  x=reorder(xlabel, desc(partition)), 
                  y=R2*100, #ymin = LCL*100, xmax=UCL*100, 
                  colour=R2type, group=R2type,
                  label=round(R2*100,2))) +  
   geom_point(position = position_dodge(width = 0.5), size=3.5) +
   geom_errorbar(aes(ymin = LCL*100,
                     ymax = UCL*100),
                 position = position_dodge(width = 0.5),
                 width = 0, linewidth=1) +
   scale_color_manual(values=MetBrewer::met.brewer("Johnson", 2),
                      expression(paste(R^2)))+
   ggrepel::geom_label_repel(position = pos, size = rel(4),  show.legend = F,min.segment.length = 0,
                             point.padding = NA, box.padding = 0.5)+ 
   scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, .1))) + 
   scale_x_discrete(expand = expansion(add = .7)) + #.6 is default - add a little spacing to the sides
   ylab("") +  xlab(paste0("HCM diagnosis ~"))+
   theme_bw() +
   theme(
     legend.position = "none",  
     text=element_text(lineheight = 0.8, angle = 0, size = 21, color = "gray8")
   )+
   ## first comparison: original vs original_plus_squared_original
   geom_bracket(data = df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,],
                xmin = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[1]],
                xmax = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[2]],
                y.position = max(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$UCL*100)*1.08,
                label =paste0("MD=",comparisons[comparisons$comparison=="original vs original_plus_squared_original",2]*100,
                              ", p=",comparisons[comparisons$comparison=="original vs original_plus_squared_original",5]),
                # tip.length = c(0.02, 0.02),
                coord.flip = TRUE,
                # vjust = -1
   ) +
   ## second comparison: original vs residual + EP
   geom_bracket(data = df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,],
                xmin = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[1]],
                xmax = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[3]],
                y.position = max(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$UCL*100)*1.16,
                label = paste0("MD=",comparisons[comparisons$comparison=="original vs residual + EP",2]*100,
                               ", p=",comparisons[comparisons$comparison=="original vs residual + EP",5]),
                # tip.length = c(0.02, 0.02),
                coord.flip = TRUE,
                # vjust = -1
   )+
   ## third comparison: original vs residual + EP+ residualxEP
   geom_bracket(data = df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,],
                xmin = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[1]],
                xmax = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[4]],
                y.position = max(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$UCL*100)*1.24,
                label = paste0("MD=",comparisons[comparisons$comparison=="original vs residual + EP+ residualxEP",2]*100,
                               ", p=",comparisons[comparisons$comparison=="original vs residual + EP+ residualxEP",5]),
                # tip.length = c(0.02, 0.02),
                coord.flip = TRUE,
                # vjust = -1
   )+
   ## third comparison: original_plus_squared_original vs residual+EP + interactions+quadraticterms
   geom_bracket(data = df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,],
                xmin = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[2]],
                xmax = reorder(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$xlabel, desc(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$partition))[[5]],
                y.position = max(df_plot[grepl(pattern = "^0|^3", perl = T, x=df_plot$partition) ,]$UCL*100)*1.32,
                label = paste0("MD=",comparisons[comparisons$comparison=="original_plus_squared_original vs residual+EP + interactions+quadraticterms",2]*100,
                               ", p=",comparisons[comparisons$comparison=="original_plus_squared_original vs residual+EP + interactions+quadraticterms",5]),
                # tip.length = c(0.02, 0.02),
                coord.flip = TRUE,
                # vjust = -1
   )+
   coord_flip())
 
 

fig3_grob= arrangeGrob(
  textGrob(paste("CoDs % and 95% CIs for the Original GWAS PRS vs\nAdditive Models Including the Residual and Enhancer Partitions"), gp = gpar(fontsize = 22, fontface = "bold", col="darkred")), 
  fig3, 
  ncol=1, heights = c(0.2,1)
)


ggsave(
  filename = paste0(OUTPUT_prefix, "fig3_", cohort_ENHlist_thresh, "_", Sys.Date(),"_CoD_original_vs_partitioned_models.pdf"), 
  fig3_grob,  
  width = 10, height = 4, device = "pdf", scale = 1.5)



## OVERALL PLOT ###
fig2_grob_modif = arrangeGrob(
  textGrob("C)", just = "left",
           gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
  textGrob(paste("CoDs % and 95% CIs for the three enhancer partitions:\nOriginal OR, enhanced by ES, enhanced by expression"), gp = gpar(fontsize = 22, fontface = "bold", col="darkblue")), 
  fig2 + theme(legend.position = "none"), 
  layout_matrix=rbind(c(1,2),
                      c(3,3)),
  widths = c(0.1, 1), heights = c(0.15, 1)
)
fig3_modif <- fig3 + theme(axis.title.y =element_blank(), legend.position = "none")
fig3_grob_modif= arrangeGrob(
  textGrob("D)", just = "left",
           gp = gpar(fontsize = 18, fontface = "bold", col="black")), 
  textGrob(paste("CoDs % and 95% CIs for the Original GWAS PRS vs\nAdditive Models Including the Residual and Enhancer Partitions"), gp = gpar(fontsize = 20, fontface = "bold", col="darkred")), 
  fig3_modif, 
  layout_matrix=rbind(c(1,2),
                      c(3,3)),
  widths = c(0.05, 1), heights = c(0.15, 1)
)
ggsave(filename = paste0(OUTPUT_prefix, "overall_", cohort_ENHlist_thresh, "_", Sys.Date(),"_all_plots.pdf"), 
       arrangeGrob(
         textGrob(addline_format(ENH_list),
                  gp = gpar(fontsize = 24, fontface = "bold",col="navyblue")),
         fig1grob, fig2_grob_modif, fig3_grob_modif, 
         layout_matrix=rbind(c(1,1),
                             c(2,2),
                             c(3,4)),
         heights = c(0.08, 0.5, 0.5, 0.04), widths = c(1,1.2)       ),  
       width = 20, height = 14, device = "pdf", scale = 1)








# OR based plots
# double quantile plot for interactions

#https://cran.r-project.org/web/packages/samplesizeCMH/vignettes/samplesizeCMH-introduction.latex
#https://stats.stackexchange.com/questions/593123/can-i-add-up-ors-for-specific-predictors/593130#593130
#https://sphweb.bumc.bu.edu/otlt/mph-modules/bs/bs704-ep713_confounding-em/BS704-EP713_Confounding-EM7.latex
scaled_BEST_PRS_score_per_UKBB_participant

# CALCULATE ORS
#merged_GWAS_q_OR
summary(logistic<-glm(formula = dx ~ original_GWAS_q, 
                      data = scaled_BEST_PRS_score_per_UKBB_participant, family = binomial, na.action = "na.omit"))
(original_GWAS_q_OR<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
colnames(original_GWAS_q_OR) <- c("ENH_compartment_quantile", "OR", "LCI", "UCI")
original_GWAS_q_OR

(ORs <- original_GWAS_q_OR)
ORs[1,]<-list("1",1,1,1)
ORs[,1]<-list(1:nrow(ORs))
ORs$original_OR_quant <- "All"
ORs

sink(paste0(OUTPUT_prefix, cohort_ENHlist_thresh, "_", Sys.Date(),"_ORs.log"), split = T)
ORs_original_OR <- ORs
for  (i in 1:number_quantiles) {
  print(i)
  logistic<-glm(formula = dx ~ original_GWAS_q, 
                data = scaled_BEST_PRS_score_per_UKBB_participant[scaled_BEST_PRS_score_per_UKBB_participant$TS_ENH_compartment_originalOR_q==i,],
                family = binomial, na.action = "na.omit")
  (OR<-exp(cbind(coef(logistic), confint(logistic))) %>% as_tibble(rownames = "quant"))
  OR[1,]<-list("1",1,1,1)
  OR[,1]<-list(1:nrow(OR))
  colnames(OR) <- c("ENH_compartment_quantile", "OR", "LCI", "UCI")
  OR$original_OR_quant <- paste0("Enh q",i)
  OR
  
  ORs_original_OR<-rbind(ORs_original_OR,OR)
}
ORs_original_OR



(all_ORs<-
    ORs_original_OR)
all_ORs$original_OR_quant = factor(all_ORs$original_OR_quant)
all_ORs$original_OR_quant = relevel(all_ORs$original_OR_quant, ref = "All")
sink()

# pdf(file = PRS_double_QUANTILE_PLOT, width = 11, height = 7)
(p4 = ggplot(data = all_ORs , aes(y= OR, ymin = LCI, ymax=UCI, 
                                  x=factor(ENH_compartment_quantile), colour=original_OR_quant, group=original_OR_quant)) + 
    # facet_wrap(facets = vars((comp)))+
    scale_colour_manual(name="ENH compartment quantile", values = c("tomato",MetBrewer::met.brewer("Hokusai2",number_quantiles)))+
    geom_pointrange(position = position_dodge(width = 0.3))  + 
    ylab(paste0("OR for HCM"))+   xlab('Original PRS quantile')+
    # labs(title =  paste("Participant distribution by HCM OR by original PGC GWAS quantile\nand further by", ENH_list, "quantile"))+ 
    theme_bw() +
    theme(
      strip.text.x = element_text(size = rel(1.3)),
      axis.text = element_text(size = rel(1)),
      axis.title = element_text(size = rel(1.5)),
      plot.margin = margin(t = 0, r = 1, b = 0, l = 1, "cm"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    ))

# dev.off()
f4<-arrangeGrob(
  textGrob(paste0("Participant distribution by OR for HCM, first by original GWAS quantile (in red)\nand further by ", addline_format(ENH_list), " quantile (shades of blue)"), gp = gpar(fontsize = 16, fontface = "bold",col="maroon")), 
  p4,
  ncol=1,
  heights = c(0.1, 1))


ggsave(filename = paste0(OUTPUT_prefix, cohort_ENHlist_thresh, "_", Sys.Date(),"_Quant_by_quant_plot.pdf"),
       f4,  width = 9, height = 7)

