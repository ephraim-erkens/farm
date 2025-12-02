#### script for random forest training using mlr3 ranger

# packages for general data handling 
library(tidyverse)
library(dplyr)
library(stringi)
# mlr3 packages
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3verse)
library(ranger) # specific learner (efficient random forest)
# packages for parallel computing
library(future)
library(future.apply)


#---------------------------------------------------------------------------


#### load custom performance measures for evaluation
# custom measure class for F1-score in multiclass classification
MeasureF1Multiclass = readRDS(file = "../custom_performance_measures/MeasureF1Multiclass.RDS")
# Register the measure
f1_multiclass_measure <- MeasureF1Multiclass$new() # Instantiate the custom measure
mlr3::mlr_measures$add("f1_multiclass_measure", MeasureF1Multiclass) # add the new measure to the dictionary
# load function to calculate precision, recall and f1 score per class for both mlr3 prediction objects and prepared dataframes
perf_rec_prec_f1 = readRDS(file = "../custom_performance_measures/perf_rec_prec_f1.RDS")


#-------------------------------------------------------------------------------


#-------------------------------------------
### hyperparameter settings
SplitRule = "gini"
mtry = c(1, 2, 3, 4, 5, 6)
mns = c(10)
num_trees = c(1000)
#-------------------------------------------


# set working directory (adapt to own directory!)
### 
setwd(".../training/RF/")
###

# create output folder if not existent
path_4_saving0 = paste("RF_output_data/")
dir.create(path_4_saving0)

# read compound names for which models are build individually
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

### start of main loop to choose compound
for(k in seq_along(compound_names)) {

  # select and read compound file (has been shuffled before to make sure data frame has no underlying structure)
  compound_name = compound_names[k] 
  df_compound_k <- read.csv(paste0("../input_data/input_data_", compound_name, ".csv"))
  cat("\n", k, "\n", compound_name, "\n", paste("number of sites:", nrow(df_compound_k)), "\n") # to keep track
  
  # create output folder for compound if not existent
  path_4_saving = paste0(path_4_saving0, compound_name, "/")
  dir.create(path_4_saving)


  #-----------------------------------------------------------------------------
  ### start of random forest training
  
  
  # define and assign task for data set
  tsk = as_task_classif(conc_group ~ ., data = df_compound_k)
  tsk$feature_names
  # save tasks
  saveRDS(tsk, file = paste(path_4_saving, "task_", compound_name, "_conc_group.RDS", sep = ""))
  
  ### loop over different seeds to create several reproducible data set splits 
  seeds = c(123, 234, 345, 456, 567, 678, 789, 987, 876, 765)
  
  for (w in 1:length(seeds)){
    
    seed = seeds[w]
    set.seed(seed)
    
    ### set up outer cross-validation splits
    outer_resampling = rsmp("cv", folds = 5)
    outer_resampling$instantiate(tsk); outer_resampling
    
    #---------------------------------------------------------------------------
    ### 1. perform nested cross-validation
    # after tuning the hyperparameters in the inner loop, train the model on the respective outer train set using the best hyperparameters, respectively
    
    ### loop over different outer folds of nested cross-validation
    HP_f1_thresh_max_inner_folds = data.frame() # best performances inner folds for all HP combinations
    HP_f1_thresh_max_outer_folds = data.frame() # best performances outer folds of best HP combinations for each outer fold
    rr_outer_list = list() # create list to save models of outer loop
    pred_outer_train_list = list() # create list to save predictions on train set of outer loop
    pred_outer_test_list = list() # create list to save predictions on test set of outer loop
    outer_train_task_list = list() # create list to save train task of outer loop
    outer_test_task_list = list() # create list to save test task of outer loop
    for (m in 1:5){
    
      cat("\n", paste0("outer fold: ", m), "\n") # to keep track
      outer_fold_train = outer_resampling$train_set(m)
      outer_fold_test = outer_resampling$test_set(m)
      outer_train_task = tsk$clone(deep = TRUE)$filter(outer_fold_train)
      outer_test_task = tsk$clone(deep = TRUE)$filter(outer_fold_test)
      outer_train_task_list[[m]] = outer_train_task$row_roles$use
      outer_test_task_list[[m]] = outer_test_task$row_roles$use
      # outer_train_task$row_roles
      # outer_test_task$row_roles
      # tsk$row_roles
    
      
      # loop over different HP combinations for tuning
      for (r in 1:length(mtry)){
        for (s in 1:length(mns)){
          for (t in 1:length(num_trees)){
    
          HP_mtry = mtry[r]
          HP_mns = mns[s]
          HP_num_trees = num_trees[t]
          
          cat("\n", paste0("mtry = ", HP_mtry, "; mns = ", HP_mns, "; num_tree = ", HP_num_trees, "\n"))
          
          lrner = lrn("classif.ranger", predict_sets=c("train","test"), importance ="impurity", predict_type="prob",
                      splitrule = SplitRule, num.trees = HP_num_trees, seed = seed,
                      mtry = HP_mtry, min.node.size = HP_mns, oob.error = TRUE, replace = TRUE, keep.inbag = TRUE)
          
          future::plan("cluster",workers=50, gc=T)
          rr_inner = lrner$train(outer_train_task)
          rr_inner$model$inbag.counts
    
          response = rep(NA, length(rr_inner$model$predictions[,1]))
          df1_test = as.data.frame(bind_cols(outer_train_task$data(cols = outer_train_task$target_names),
                                             response, outer_train_task$row_ids, rr_inner$model$predictions))
          names(df1_test) = c("truth", "response", "row_ids", paste0("prob", 1:(length(df1_test)-3)))
        
        ### calculate precision, recall and f1 score for different prob thresholds
        prob_mat_test = as.matrix(df1_test[,c(4,5)])
        colnames(prob_mat_test) = c("1", "2")
        # write the bound test sets of all inner folds as prediction to easily optimize threshold and access performances
        new_prediction = PredictionClassif$new(
          task = outer_train_task,
          response = df1_test$response,  
          truth = df1_test$truth,              
          prob = prob_mat_test)
    
        # perform decision threshold optimisation based on train sets
        prec_rec_f1_inner = data.frame()
        for (v in seq(0.0, 1.0, 0.001)){ # for each threshold separately
          new_prediction$set_threshold(v)
          df1 = perf_rec_prec_f1(new_prediction) # call custom function to calculate recall, precision and f1-score for each class independently
          df1$threshold = v
          prec_rec_f1_inner = dplyr::bind_rows(prec_rec_f1_inner, df1)
        }
        
        # calculate macro-f1-score for each threshold and select threshold with maximal f1 score
        prec_rec_f1_inner = prec_rec_f1_inner %>% dplyr::mutate(across(1:3, as.numeric))
        f1_macro_thresh = prec_rec_f1_inner %>%
          dplyr::group_by(threshold) %>%
          dplyr::summarize(macro_f1_test = mean(f1_score, na.rm = FALSE))
        f1_thresh_max_inner = f1_macro_thresh[which.max(f1_macro_thresh$macro_f1_test),] # best threshold
        # apply optimal threshold to inner train predictions to calculate f1 train score
        pred_inner_train = rr_inner$predict(outer_train_task)
        pred_inner_train$set_threshold(f1_macro_thresh[which.max(f1_macro_thresh$macro_f1_test),]$threshold)
        df1_train = pred_inner_train$score(msr("f1_multiclass_measure")) ;df1_train
        
        f1_thresh_max_inner$macro_f1_train = df1_train
        f1_thresh_max_inner$mtry = HP_mtry
        f1_thresh_max_inner$mns = HP_mns
        f1_thresh_max_inner$num_trees = HP_num_trees
        f1_thresh_max_inner$outer_fold = m
    
        HP_f1_thresh_max_inner_folds = dplyr::bind_rows(HP_f1_thresh_max_inner_folds, f1_thresh_max_inner)
            
          } # loops for different HP combinations
        }
      }
    
      
      ### retrain best model of outer fold on training data of outer fold and obtain unbiased performances estimates by predicting on the outer test sets, respectively
      # extract HP combination of best performing model of the respective outer fold
      best_HP_f1_thresh_max_outer_fold = HP_f1_thresh_max_inner_folds %>% filter(outer_fold == m)
      best_HP_f1_thresh_max_outer_fold = best_HP_f1_thresh_max_outer_fold[which.max(best_HP_f1_thresh_max_outer_fold$macro_f1_test),]
      lrner_best_outer_fold = lrn("classif.ranger", predict_sets=c("train","test"), importance ="impurity", predict_type="prob",
                                  splitrule = SplitRule, num.trees = best_HP_f1_thresh_max_outer_fold$num_trees, seed = seed,
                                  mtry = best_HP_f1_thresh_max_outer_fold$mtry,min.node.size = best_HP_f1_thresh_max_outer_fold$mns,
                                  oob.error = TRUE, replace = TRUE, keep.inbag = TRUE)
      future::plan("cluster", workers = 50, gc = T)
      rr_outer = lrner_best_outer_fold$train(outer_train_task)
      # predict on train and test set
      pred_outer_train = rr_outer$predict(outer_train_task)
      pred_outer_test = rr_outer$predict(outer_test_task)
      # apply best threshold optimized in inner folds
      pred_outer_train$set_threshold(best_HP_f1_thresh_max_outer_fold$threshold)
      pred_outer_test$set_threshold(best_HP_f1_thresh_max_outer_fold$threshold)
    
      # call custom function to calculate recall, precision and f1-score for train and test set for each class independently
      perf_outer_train = perf_rec_prec_f1(pred_outer_train) ;perf_outer_train
      perf_outer_test = perf_rec_prec_f1(pred_outer_test) ;perf_outer_test
      perf_outer = dplyr::bind_cols(perf_outer_train[, 1:3], perf_outer_test)
      colnames(perf_outer) = c("recall_train", "precision_train", "f1_score_train", "recall_test", "precision_test", "f1_score_test", "class")
      perf_outer = perf_outer %>% dplyr::mutate(across(1:6, as.numeric))
      perf_outer = perf_outer %>% mutate(mtry = best_HP_f1_thresh_max_outer_fold$mtry,
                                         mns = best_HP_f1_thresh_max_outer_fold$mns,
                                         num_trees = best_HP_f1_thresh_max_outer_fold$num_trees,
                                         threshold = best_HP_f1_thresh_max_outer_fold$threshold)
      perf_outer$outer_fold = m
      HP_f1_thresh_max_outer_folds = dplyr::bind_rows(HP_f1_thresh_max_outer_folds, perf_outer)
      
      # save models and predictions with best HP combination trained on the respective outer folds
      rr_outer_list[[m]] = rr_outer
      pred_outer_train_list[[m]] = pred_outer_train
      pred_outer_test_list[[m]] = pred_outer_test
    
    } # loop for outer folds
    
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
        mtry = first(mtry),
        mns = first(mns),
        num_trees = first(num_trees),
        threshold = first(threshold))
    
    # average measures over all outer folds
    HP_thresh_max_perf_final =
      column_means <- sapply(HP_thresh_max_perf_outer_folds %>% select(-c(outer_fold, mtry, mns, num_trees, threshold)),
                             mean, na.rm = TRUE); HP_thresh_max_perf_final
    
    
    ### save important files
    saveRDS(HP_f1_thresh_max_inner_folds, file=paste(path_4_saving, "HP_f1_thresh_max_inner_folds_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(HP_f1_thresh_max_outer_folds, file=paste(path_4_saving, "HP_f1_thresh_max_outer_folds_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(HP_thresh_max_perf_outer_folds, file=paste(path_4_saving, "HP_thresh_max_perf_outer_folds_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(HP_thresh_max_perf_final, file=paste(path_4_saving, "HP_thresh_max_perf_final_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(rr_outer_list, file=paste(path_4_saving, "rr_outer_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(pred_outer_train_list, file=paste(path_4_saving, "pred_outer_train_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(pred_outer_test_list, file=paste(path_4_saving, "pred_outer_test_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(outer_train_task_list, file=paste(path_4_saving, "outer_train_task_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(outer_test_task_list, file=paste(path_4_saving, "outer_test_task_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    
    
    #---------------------------------------------------------------------------
    ### 2. perform single cross-validation
    # based on entire data set to determine best hyperparameter configuration to be used for a single final model
    
    set.seed(seed)
    single_cv = rsmp("cv", folds = 5)
    single_cv$instantiate(tsk); single_cv
    cat("\n", "single cv", "\n")
    
    HP_f1_thresh_max_single_folds = data.frame() # best performances single cv of best HP combinations for each outer fold
    for (r in 1:length(mtry)){
      for (s in 1:length(mns)){
        for (t in 1:length(num_trees)){
    
        HP_mtry = mtry[r]
        HP_mns = mns[s]
        HP_num_trees = num_trees[t]
    
        cat("\n", paste0("mtry = ", HP_mtry, "; mns = ", HP_mns, "; num_tree = ", HP_num_trees, "\n"))
    
        lrner_single = lrn("classif.ranger", predict_sets=c("train","test"), importance ="impurity", predict_type="prob",
                           splitrule = SplitRule, num.trees = HP_num_trees, seed = seed,
                    mtry = HP_mtry, min.node.size = HP_mns, oob.error = TRUE, replace = TRUE, keep.inbag = TRUE)
    
        future::plan("cluster",workers=50, gc=T)
        rr_single = lrner$train(tsk)
    
        response = rep(NA, length(rr_single$model$predictions[,1]))
        df1_test = as.data.frame(bind_cols(tsk$data(cols = tsk$target_names), response, tsk$row_ids, rr_single$model$predictions))
        names(df1_test) = c("truth", "response", "row_ids", paste0("prob", 1:(length(df1_test)-3)))
    
        ### calculate precision, recall and f1 score for different prob thresholds - only for binary!
        prob_mat_test = as.matrix(df1_test[,c(4,5)])
        colnames(prob_mat_test) = c("1", "2"); prob_mat_test
        # write the bound test sets of all inner folds as prediction to easily optimize threshold and access performances
        new_prediction = PredictionClassif$new(
          task = tsk,
          response = df1_test$response,
          truth = df1_test$truth,
          prob = prob_mat_test)
    
        prec_rec_f1_single = data.frame()
        for (v in seq(0.0, 1.0, 0.001)){ # for each threshold separately
          new_prediction$set_threshold(v)
          df1 = perf_rec_prec_f1(new_prediction) # call custom function to calculate recall, precision and f1-score for each class independently
          df1$threshold = v
          prec_rec_f1_single = dplyr::bind_rows(prec_rec_f1_single, df1)
        }
    
        # calculate macro-f1-score for each threshold and select threshold with maximal f1 score
        prec_rec_f1_single = prec_rec_f1_single %>% dplyr::mutate(across(1:3, as.numeric))
        f1_macro_thresh2 = prec_rec_f1_single %>%
          dplyr::group_by(threshold) %>%
          dplyr::summarize(macro_f1_test = mean(f1_score, na.rm = FALSE))
        f1_thresh_max_single = f1_macro_thresh2[which.max(f1_macro_thresh2$macro_f1_test),] # best threshold
        # apply optimal threshold to inner train predictions to calculate f1 train score
        pred_single_train = rr_single$predict(tsk)
        pred_single_train$set_threshold(f1_macro_thresh2[which.max(f1_macro_thresh2$macro_f1_test),]$threshold)
        df1_train = pred_single_train$score(msr("f1_multiclass_measure"))
    
        f1_thresh_max_single$macro_f1_train = df1_train
        f1_thresh_max_single$mtry = HP_mtry
        f1_thresh_max_single$mns = HP_mns
        f1_thresh_max_single$num_trees = HP_num_trees
    
        HP_f1_thresh_max_single_folds = dplyr::bind_rows(HP_f1_thresh_max_single_folds, f1_thresh_max_single)
    
    
    
          } # loops for different HP combinations
        }
      }
    
    # extract HP combination of best performing model of single cv
    best_HP_f1_thresh_max_single_cv = HP_f1_thresh_max_single_folds[which.max(HP_f1_thresh_max_single_folds$macro_f1_test),]
    
    
    #---------------------------------------------------------------------------
    ### 3. final model training
    # train a model with best configuration from single cross-validation on entire data set
    
    cat("\n", "final model training", "\n")
    lrner_best = lrn("classif.ranger", predict_sets=c("train","test"), importance ="impurity", predict_type="prob",
                     splitrule = SplitRule, num.trees = best_HP_f1_thresh_max_single_cv$num_trees, seed = seed,
                                mtry = best_HP_f1_thresh_max_single_cv$mtry, min.node.size = best_HP_f1_thresh_max_single_cv$mns,
                     oob.error = TRUE, replace = TRUE, keep.inbag = TRUE)
    future::plan("cluster",workers=50, gc=T)
    train_full = lrner_best$train(tsk)
    
    
    ### save important files
    saveRDS(HP_f1_thresh_max_single_folds, file=paste(path_4_saving, "HP_f1_thresh_max_single_folds_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(best_HP_f1_thresh_max_single_cv, file=paste(path_4_saving, "best_HP_f1_thresh_max_single_cv_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    saveRDS(train_full, file=paste(path_4_saving,"full_model_", compound_name,"_conc_group_", w,".RDS", sep="")) # final model
    
    
    #---------------------------------------------------------------------------
    
    
    ### train again using corrected gini impurity feature importance
    # only use HPs of best models for each outer fold, respectively, to reduce computational costs
    rr_outer_gini_cor_list = list() # create list to save models of outer loop for gini corrected
    cat("\n", "gini impurity corrected", "\n")
    for (m in 1:5){
    
      cat("\n", paste0("outer fold: ", m), "\n")
      outer_fold_train = outer_resampling$train_set(m)
      outer_fold_test = outer_resampling$test_set(m)
      outer_train_task = tsk$clone(deep = TRUE)$filter(outer_fold_train)
      outer_test_task = tsk$clone(deep = TRUE)$filter(outer_fold_test)
    
      ### train best model of outer fold on entire training data of outer fold evaluating feature importance by corrected gini impurity
      # extract HP combination of best performing model of the respective outer fold and assign it to new learner
      lrner_best_outer_fold_gini_cor = lrn("classif.ranger", predict_sets=c("train","test"), importance ="impurity_corrected", predict_type="prob",
                                           splitrule = SplitRule, num.trees = HP_thresh_max_perf_outer_folds[m,]$num_trees, seed = seed,
                                  mtry = HP_thresh_max_perf_outer_folds[m,]$mtry, min.node.size = HP_thresh_max_perf_outer_folds[m,]$mns,
                                  oob.error = TRUE, replace = TRUE, keep.inbag = TRUE)
      future::plan("cluster", workers = 50, gc = T)
      rr_outer_gini_cor = lrner_best_outer_fold_gini_cor$train(outer_train_task)
      # save models to later extract corrected gini impurity
      rr_outer_gini_cor_list[[m]] = rr_outer_gini_cor
    
    }
    saveRDS(rr_outer_gini_cor_list, file=paste(path_4_saving, "rr_outer_gini_cor_list_", compound_name, "_conc_group_", w, ".RDS", sep = ""))
    
    
  
  } # loop for different seeds
} # loop for different compounds


