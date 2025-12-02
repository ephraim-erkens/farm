#### script for logistic regression training

# packages for general data handling 
library(tidyverse)
library(dplyr)
library(stringi)
# mlr3 packages
library(mlr3)
library(mlr3learners)


#-------------------------------------------------------------------------------


#### load custom performance measures for evaluation
# load function to calculate precision, recall and f1 score per class for both mlr3 prediction objects and prepared dataframes
perf_rec_prec_f1 = readRDS(file = "../custom_performance_measures/perf_rec_prec_f1.RDS")


#-------------------------------------------------------------------------------


# set working directory (adapt to own directory!)
### 
setwd(".../training/LR/")
###

# create output folder if not existent
path_4_saving0 = paste("LR_output_data/")
dir.create(path_4_saving0)

# read compound names for which models are build individually
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name


for(k in seq_along(compound_names)) {
  
  # select and read compound file (has been shuffled before to make sure data frame has no underlying structure)
  compound_name = compound_names[k] 
  df_compound_k <- read.csv(paste0("../input_data/input_data_", compound_name, ".csv"))
  cat("\n", k, "\n", compound_name, "\n", paste("number of sites:", nrow(df_compound_k)), "\n") # to keep track
  
  # create output folder for compound if not existent
  path_4_saving = paste0(path_4_saving0, compound_name, "/")
  dir.create(path_4_saving)

  
  #---------------------------------------------------------------------------
  ### start of logistic regression training
  
  
  # define and assign task for data set
  tsk = as_task_classif(conc_group ~ ., data = df_compound_k)
  tsk$feature_names
  # save tasks
  saveRDS(tsk, file = paste(path_4_saving, "task_", compound_name, "_conc_group.RDS", sep = ""))
  
  ### loop over different seeds to create several reproducible data set splits 
  # for TFA seeds are slightly adapted to prevent folds where not all feature attribute classes are included
  if (k == 8) { # TFA
    seeds = c(123, 234, 345, 456, 567, 987, 765, 543, 432, 321)
  } else {
    seeds = c(123, 234, 345, 456, 567, 678, 789, 987, 876, 765)
  }
  
  for (w in 1:length(seeds)){
  
    seed = seeds[w]
    set.seed(seed)

    ### set up single cross-validation splits
    outer_resampling = rsmp("cv", folds = 5)
    outer_resampling$instantiate(tsk); outer_resampling
    
    #---------------------------------------------------------------------------
    ### 1. perform single cross-validation
    
    ### loop over different folds of single cross-validation
    HP_f1_thresh_max_outer_folds = data.frame() # best performances outer folds of best HP combinations for each outer fold
    rr_outer_list = list() # create list to save models of outer loop
    pred_outer_train_list = list() # create list to save predictions on train set of outer loop
    pred_outer_test_list = list() # create list to save predictions on test set of outer loop
    outer_train_task_list = list() # create list to save train task of outer loop
    outer_test_task_list = list() # create list to save test task of outer loop
    for (m in 1:5){
    
      cat("\n", paste0("seed : ", w, "; outer fold: ", m), "\n") # to keep track
      outer_fold_train = outer_resampling$train_set(m)
      outer_fold_test = outer_resampling$test_set(m)
      outer_train_task = tsk$clone(deep = TRUE)$filter(outer_fold_train)
      outer_test_task = tsk$clone(deep = TRUE)$filter(outer_fold_test)
      outer_train_task_list[[m]] = outer_train_task$row_roles$use
      outer_test_task_list[[m]] = outer_test_task$row_roles$use
      # outer_train_task$row_roles
      # outer_test_task$row_roles
      # tsk$row_roles

      # train logistic regression model on training set and predict on test set using probabilities
      lrn_logreg = lrn("classif.log_reg", predict_type = "prob")
      rr_logreg = lrn_logreg$train(outer_train_task)
      rr_logreg$model$predictions
      # predict on train and test set
      pred_train_bench = rr_logreg$predict(outer_train_task)
      pred_test_bench = rr_logreg$predict(outer_test_task)
  
      # perform decision threshold optimisation based on train sets
      prec_rec_f1_inner = data.frame()
      for (v in seq(0.0, 1.0, 0.001)){ # for each threshold separately
        pred_train_bench$set_threshold(v)
        df1 = perf_rec_prec_f1(pred_train_bench) # call custom function to calculate recall, precision and f1-score for each class independently
        df1$threshold = v
        prec_rec_f1_inner = dplyr::bind_rows(prec_rec_f1_inner, df1)
      }
      
      # calculate macro-f1-score for each threshold and select threshold with maximal f1 score
      prec_rec_f1_inner = prec_rec_f1_inner %>% dplyr::mutate(across(1:3, as.numeric))
      f1_macro_thresh = prec_rec_f1_inner %>%
        dplyr::group_by(threshold) %>%
        dplyr::summarize(macro_f1_train = mean(f1_score, na.rm = FALSE))
      f1_thresh_max_inner = f1_macro_thresh[which.max(f1_macro_thresh$macro_f1_train),] # best threshold
      
      # apply best threshold optimized in inner folds
      pred_train_bench$set_threshold(f1_thresh_max_inner$threshold)
      pred_test_bench$set_threshold(f1_thresh_max_inner$threshold)
      
      # call custom function to calculate recall, precision and f1-score for train and test set for each class independently
      perf_train_bench = perf_rec_prec_f1(pred_train_bench)
      perf_test_bench = perf_rec_prec_f1(pred_test_bench)
      perf_bench = dplyr::bind_cols(perf_train_bench[, 1:3], perf_test_bench)
      colnames(perf_bench) = c("recall_train", "precision_train", "f1_score_train", "recall_test", "precision_test", "f1_score_test", "class")
      perf_bench = perf_bench %>% dplyr::mutate(across(1:6, as.numeric))
      perf_bench$threshold = f1_thresh_max_inner$threshold
      perf_bench$outer_fold = m
      HP_f1_thresh_max_outer_folds = dplyr::bind_rows(HP_f1_thresh_max_outer_folds, perf_bench)
      
      rr_outer_list[[m]] = rr_logreg
      pred_outer_train_list[[m]] = pred_train_bench
      pred_outer_test_list[[m]] = pred_test_bench
        
    }
    
    # average measures over all classes for each outer fold
    HP_thresh_max_perf_outer_folds =
      HP_f1_thresh_max_outer_folds %>%
      dplyr::group_by(outer_fold) %>%
      dplyr::summarise(
        macro_rec_train_outer = mean(recall_train, na.rm = TRUE),
        macro_prec_train_outer = mean(precision_train, na.rm = TRUE),
        macro_f1_train_outer = mean(f1_score_train, na.rm = TRUE),
        macro_rec_test_outer = mean(recall_test, na.rm = TRUE),
        macro_prec_test_outer = mean(precision_test, na.rm = TRUE),
        macro_f1_test_outer = mean(f1_score_test, na.rm = TRUE),
        threshold = first(threshold))
    
    # average measures over all outer folds
    HP_thresh_max_perf_final =
      column_means <- sapply(HP_thresh_max_perf_outer_folds %>% dplyr::select(-c(outer_fold, threshold)), mean, na.rm = TRUE); HP_thresh_max_perf_final
    
    
    ### save important files
    saveRDS(HP_f1_thresh_max_outer_folds, file=paste(path_4_saving, "HP_f1_thresh_max_outer_folds_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(HP_thresh_max_perf_outer_folds, file=paste(path_4_saving, "HP_thresh_max_perf_outer_folds_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(HP_thresh_max_perf_final, file=paste(path_4_saving, "HP_thresh_max_perf_final_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(rr_outer_list, file=paste(path_4_saving, "rr_outer_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(pred_outer_train_list, file=paste(path_4_saving, "pred_outer_train_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(pred_outer_test_list, file=paste(path_4_saving, "pred_outer_test_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(outer_train_task_list, file=paste(path_4_saving, "outer_train_task_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(outer_test_task_list, file=paste(path_4_saving, "outer_test_task_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    
        
    #---------------------------------------------------------------------------
    ### 2. train model on entire dataset using the best HP configuration
    
    cat("\n", "final model training", "\n")
    lrner_best = lrn("classif.log_reg", predict_sets = c("train","test"), predict_type="prob")
    train_full = lrner_best$train(tsk)
    
    
    ### save important files
    saveRDS(train_full, file=paste(path_4_saving,"full_model_", compound_name,"_conc_group_", w,".RDS", sep="")) # final model
    
    #---------------------------------------------------------------------------
    

 } # loop for different seeds
} # loop for different compounds


