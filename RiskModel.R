#!/usr/bin/env Rscript

# # List of packages to install and load
# packages <- c("survival", "dplyr", "survminer", "survivalROC")
# 
# for (pkg in packages) {
#   if (!requireNamespace(pkg, quietly = TRUE)) {
#     install.packages(pkg)
#   }
#   library(pkg, character.only = TRUE)
# }

library(survival)
library(dplyr)
library(survminer)
library(survivalROC)


# Split the patients according to the specified ratio into training group and testing group
data_sampling <- function(clinical_filter_patient, ratio) {
  # Function to filter training samples and preprocess data
  clinical_filter_patient_1 <- clinical_filter_patient[clinical_filter_patient$OS_STATUS == "1", "PATIENT_ID"]
  clinical_filter_patient_0 <- clinical_filter_patient[clinical_filter_patient$OS_STATUS == "0", "PATIENT_ID"]
  train_patient1_size = quantile(range(0, length(clinical_filter_patient_1)), ratio) %>% round
  train_patient0_size = quantile(range(0, length(clinical_filter_patient_0)), (1 - ratio)) %>% round
  train_patients_1 <- sample(clinical_filter_patient_1, train_patient1_size)
  test_patients_1 <- clinical_filter_patient_1[!(clinical_filter_patient_1 %in% train_patients_1)]
  train_patients_0 <- sample(clinical_filter_patient_0, train_patient0_size)
  test_patients_0 <- clinical_filter_patient_0[!(clinical_filter_patient_0 %in% train_patients_0)]
  train_patients <- c(train_patients_0, train_patients_1)
  test_patients <- c(test_patients_0, test_patients_1)
  
  classifed_patients <- list("train_patients" = train_patients,
                             "test_patients" = test_patients)
  
  return(classifed_patients)
}

# Substrating the survival info according to the training patients and testing patients
sampled_data_preprocess <- function(chosen_patients, clinical_filter_patient, rna_filter_nor_info, chosen_genes){
  
  rownames(clinical_filter_patient) <- clinical_filter_patient$PATIENT_ID
  clinical_info <- clinical_filter_patient[chosen_patients, c("OS_STATUS", "OS_MONTHS")]
  rna_info <- t(rna_filter_nor_info[, chosen_patients])
  survival_info <- cbind(clinical_info, rna_info)
  
  return(survival_info)
}

# Module to select genes
uni_select_genes <- function(survival_info){
  # Function to select genes using univariate Cox regression
  covariates <- colnames(survival_info)[-c(1:2)]
  univ_formulas <- sapply(covariates, 
                          function(x) as.formula(paste('Surv(OS_MONTHS, OS_STATUS)~', x)))
  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = survival_info)})
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=3)
                           wald.test<-signif(x$wald["test"], digits=3)
                           beta<-signif(x$coef[1], digits=3);#coeficient beta
                           HR <-signif(x$coef[2], digits=3);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"], 3)
                           HR <- paste0(HR, " (", 
                                        HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, wald.test, p.value)
                           names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                         "p.value")
                           return(res)
                         })
  res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
  uni_selected_genes <- rownames(res[res$p.value < 0.05,])
  
  uniCox_result <- list("uni_res" = res,
                        "uni_genes" = uni_selected_genes)
  
  return(uniCox_result)
}

# Construct the Model
construct_model <- function(candidate_genes, survival_info){
  genes_formula <- paste(candidate_genes, collapse = "+")
  full_formula <- as.formula(paste("Surv(OS_MONTHS, OS_STATUS) ~", genes_formula))
  res.cox <- coxph(full_formula, data = survival_info)
  multi_cox_res <- summary(res.cox)$coefficients %>% as.data.frame
  survival_ORS_genes <- as.matrix(survival_info[, candidate_genes])
  genes_coef <- as.matrix(multi_cox_res[candidate_genes, "coef"])
  formula_str <- mapply(function(n, s) paste(n, "*", s), round(as.vector(genes_coef), 3), candidate_genes)
  formula_exp <- paste(formula_str, collapse = " + ")
  
  
  model_res <- list("cox_res" = res.cox,
                    "coef" = genes_coef,
                    "formula" = formula_exp)
  
  return(model_res)
  
}

# Validate the Model
validate_model <- function(candidate_genes, survival_info, coef, plot_prefix){
  survival_ORS_genes <- as.matrix(survival_info[, candidate_genes])
  survival_info$ORS <- (as.matrix(survival_ORS_genes) %*% coef)[,1]
  threshold.cut <- survminer::surv_cutpoint(survival_info,
                                 time = "OS_MONTHS",
                                 event = "OS_STATUS",
                                 variables = "ORS")
  ORS_threshold <- summary(threshold.cut)$cutpoint
  survival_info$ORS_group <- ifelse(survival_info$ORS > ORS_threshold,
                                    "High",
                                    "Low")
  
  survplot_name = paste(plot_prefix, "survival_plot.pdf", sep = "_")
  print(survplot_name)
  pdf(file = survplot_name, width = 16, height = 12)
  sfit <- surv_fit(Surv(OS_MONTHS, OS_STATUS) ~ ORS_group,
                  data = survival_info)
  
  plot <- survminer::ggsurvplot(data = survival_info, sfit, 
               conf.int=TRUE, pval=TRUE, risk.table=TRUE,
               legend.labs=c("High", "Low"), legend.title="ORS",
               palette=c("red", "darkblue"),
               risk.table.height=.15)
  print(plot)
  
  dev.off()
  
  return(ORS_threshold)
}

# # Module for generating plots
# generate_plot <- function(data) {
#   # Function to generate plots
#   # ...
# }
# 
# # ROC calculation module
# calculate_ROC <- function(data) {
#   # Function to calculate ROC
#   # ...
#   return(roc_results)
# }

args <- commandArgs(trailingOnly = T)

# Check if arguments are passed
if(length(args)<5) {
  stop("At least five argument must be supplied (input file).n", call.=FALSE)
} else {
  print(paste("Loading input file:", args[1], args[2], args[3]))
  clinical_filter_patient <- read.table(args[1],
                                        header = T, sep = "\t")
  colnames(clinical_filter_patient) <- c("PATIENT_ID", "OS_STATUS", "OS_MONTHS")
  
  rna_filter_nor_info <- read.table(args[2],
                                    header = T, row.names = 1, sep = "\t", check.names = F)
  chosen_genes <- readLines(args[3])
}

seed_num = as.numeric(args[4])
ratio = as.numeric(args[5])

set.seed(seed_num)


# Judge whether the input format is correct.
if (any(!(rownames(rna_filter_nor_info) == chosen_genes))){
  stop("Rows in expression file is inconsistent with the order of genes.
         The format of expression file is rows for genes, columns for patients.
         The order of genes in expression file should be consistent with the order of genes")
}

if (any(!(colnames(rna_filter_nor_info) == clinical_filter_patient$PATIENT_ID))) {
  stop("columns in expression file is inconsistent with the order of patients in clinical file.
         The format of expression file is rows for genes, columns for patients.
         The order of patients in expression file should be consiste with the order of patients in clinical file")
}

#############################################


classified_patients <- data_sampling(clinical_filter_patient = clinical_filter_patient, 
                                     ratio = ratio)

train_patients <- classified_patients$train_patients
test_patients <- classified_patients$test_patients

train_survival_info <- sampled_data_preprocess(chosen_patients = train_patients, 
                                               clinical_filter_patient = clinical_filter_patient,
                                               rna_filter_nor_info = rna_filter_nor_info)  

test_survival_info <- sampled_data_preprocess(chosen_patients = test_patients,
                                              clinical_filter_patient = clinical_filter_patient,
                                              rna_filter_nor_info = rna_filter_nor_info)

uni_res <- uni_select_genes(train_survival_info)
uni_genes <- uni_res$uni_genes
model_details <- construct_model(uni_genes, train_survival_info)
coef <- model_details$coef
train_threshold <- validate_model(uni_genes, train_survival_info, coef, "train")
test_threshold <- validate_model(uni_genes, test_survival_info, coef, "test")

total_survival_info <- sampled_data_preprocess(chosen_patients = clinical_filter_patient$PATIENT_ID,
                                               clinical_filter_patient = clinical_filter_patient,
                                               rna_filter_nor_info = rna_filter_nor_info)
total_threshold <- validate_model(uni_genes, total_survival_info, coef, "total")
