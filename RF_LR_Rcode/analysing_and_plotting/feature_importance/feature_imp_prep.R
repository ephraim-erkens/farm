#### Skript to extract, calculate and process the feature importances and save them as files,
# that will be needed to execute the script feat_imp_plots.R

library(tidyverse)
library(dplyr)
library(mlr3measures)
library(mlr3)
library(patchwork)


setwd("V:/PROJEKTE/FARM/Dokumentation/Datenübergabe_UBA/ML_Code/code_paper_submission/training/RF/")
###

######################################################
### specify folder names

folder_name = "RF_output_data" # choose folder name of models


path0 = paste(folder_name,"/",  sep="")
path_input = paste("../input_data/")


### 1
##retrieve compound group names ----
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

colors_feature_categories = readRDS( file="../colors_feature_categories_ENG.RDS")

colors_feature_categories= colors_feature_categories %>% select(-c(abbr, kultur_deutsch, kultur_englisch, kultur_sentinel_id)) %>% distinct()

target_variable = "conc_group"


#### load custom performance measures to be evaluated ----
# custom measure class for F1-score in multiclass classification
 MeasureF1Multiclass = readRDS(file = "../custom_performance_measures/MeasureF1Multiclass.RDS")
  # Register the measure
f1_multiclass_measure <- MeasureF1Multiclass$new() # Instantiate the custom measure
mlr3::mlr_measures$add("f1_multiclass_measure", MeasureF1Multiclass) # add the new measure to the dictionary

######################################################


# Loop for feature interaction strength calculation
for(k in seq_along(compound_names)) {

  compound_name = compound_names[k] 
  cat("\n", "\n", compound_name, "\n", k, "\n") # to keep track
  path = paste(path0, compound_name, "/", sep="")

### loop over different seeds

for (w in c(1:10)){
#w = 0
  cat("\n", paste("model seed =", w)) # for orientation

### manually calculate permutation importance applying the optimal threshold as well as extract gini impurity from outer model folds
# randomly shuffle single variables and evaluate drop in performance (F1-score)
# make predictions for each learner and use the corresponding test set for evaluation
# set vector of seeds for repetitive permutation and later averaging for each fold, respectively
if (k == 8 & folder_name=="LR_output_data") { # TFA and LR
  seeds = c(123, 234, 345, 456, 567, 987, 765, 543, 432, 321)
} else {
  seeds = c(123, 234, 345, 456, 567, 678, 789, 987, 876, 765)
}
# read list of outer learners, predictions as well as tasks/instances
rr_outer_list = readRDS(file=paste(path, "rr_outer_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))
if(folder_name =="RF_output_data"){
rr_outer_gini_cor_list = readRDS(file=paste(path, "rr_outer_gini_cor_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = "")) # model run to calculate corrected gini impurity
}
pred_outer_test_list1 = readRDS(file=paste(path, "pred_outer_test_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))
outer_test_task_list = readRDS(file=paste(path, "outer_test_task_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))
# read performances, optimized HP combinations and thresholds for each fold, respectively
HP_thresh_max_perf_outer_folds1 = readRDS(file=paste(path, "HP_thresh_max_perf_outer_folds_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))

feat_imp_permut_thresh = list() # create empty list to store permutation importances for each fold
feat_imp_gini_thresh = list() # create empty list to store gini impurtity importances for each fold
feat_imp_gini_cor_thresh = list() # create empty list to store corrected gini impurtity importances for each fold
f1_loss_seeds_folds = data.frame()# create data frame to safe permutation importances for all features, seeds and folds
for (j in 1:length(rr_outer_list)){ # for each learner
  
  cat("\n", paste(compound_name, "fold =", j), "\n") # for orientation
  rr_outer_j = rr_outer_list[[j]]
  # compute f1 scores for each fold from outer test set predictions (baseline performance)
  pred_outer_test = pred_outer_test_list1[[j]]
  f1_score_thresh = pred_outer_test$score(msr("f1_multiclass_measure"));f1_score_thresh
  outer_test_inst = outer_test_task_list[[j]]
  
  # load preprocessed dataset for compound and make sure target variable is positioned as last column
  df_data_final = readRDS(file=paste(path_input,"df_data_final_",compound_name,"_", target_variable, ".RDS", sep=""))
  df_data_final = df_data_final %>% select(-all_of(target_variable), everything(), all_of(target_variable))
  # calculate drops in performance after randomly shuffling each variable and predicting again on test set
  f1_loss_seeds = data.frame(matrix(NA, nrow = length(seeds), ncol = ncol(df_data_final) - 1)) # create data frame to safe permutation importances of all features and seeds per fold
  names(f1_loss_seeds) = colnames(df_data_final)[-length(colnames(df_data_final))] # remove last column name as it is the target variable
  f1_loss_seeds$seed = seeds
  for (q in 1:(length(df_data_final) - 1)){  # for each variable
    
    cat(paste0("\t", colnames(df_data_final)[q]))
    ### perform several permutations for later averaging for each fold, respectively
    for (p in 1:length(seeds)){ # for each seed
      
      df_data_perm = df_data_final
      # randomly shuffle the values in column q
      set.seed(seeds[p])
      df_data_perm[, q] = sample(df_data_perm[, q])
      
      # define task for new dataset
      tsk_perm = as_task_classif(conc_group ~ ., data = df_data_perm)
      # predict on test set applying new dataset to learner of corresponding fold (model should NOT be trained again with train set of new dataset!)
      pred_perm = rr_outer_j$predict(tsk_perm, row_ids = outer_test_inst)
      # check number of classes for setting threshold correctly
      if (length(unique(df_data_perm$conc_group)) == 2) {
        pred_perm$set_threshold(HP_thresh_max_perf_outer_folds1$threshold[[j]])
        # df1_perm = as.data.frame(bind_cols(pred_perm$truth, pred_perm$response, pred_perm$row_ids, pred_perm$prob))
        # names(df1_perm) = c("truth", "response", "row_ids", paste0("prob", 1:(length(df1_perm)-3)))
      } else if (length(unique(df_data_perm$conc_group)) == 3) {
        pred_perm$set_threshold(c("1" = HP_thresh_max_perf_outer_folds1$threshold1[[j]], "2" = HP_thresh_max_perf_outer_folds1$threshold2[[j]], "3" = HP_thresh_max_perf_outer_folds1$threshold3[[j]]))
      }
      
      f1_score_thresh_perm = pred_perm$score(msr("f1_multiclass_measure")); f1_score_thresh_perm
      f1_loss = f1_score_thresh - f1_score_thresh_perm; f1_loss
      f1_loss_seeds[p, q] = f1_loss
      f1_loss_seeds$outer_fold = j
      
    } # loop seeds
    
  } # loop variables
  f1_loss_seeds_folds = dplyr::bind_rows(f1_loss_seeds_folds, f1_loss_seeds)
  
  if(folder_name =="RF_output_data"){
  feat_imp_gini_thresh[[j]] = rr_outer_j$model$variable.importance # extract gini impurity folds
  feat_imp_gini_cor_thresh[[j]] = rr_outer_gini_cor_list[[j]]$model$variable.importance # extract corrected gini impurity folds
  }
} # loop outer folds

# average permutation importances per fold over all seeds
f1_mean_loss_folds = f1_loss_seeds_folds %>%
  dplyr::group_by(outer_fold) %>%
  dplyr::summarise(across(-seed, mean))
# Convert rows into a list of named vectors
feat_imp_permut_thresh = split(f1_mean_loss_folds [, -1], seq(nrow(f1_mean_loss_folds)))
feat_imp_permut_thresh = lapply(feat_imp_permut_thresh, function(x) setNames(as.vector(t(x)), names(f1_mean_loss_folds[, -1])))

saveRDS(feat_imp_permut_thresh, file = paste(path, "feat_imp_permutthresh_", compound_name,"_", target_variable, "_", w,".RDS", sep="")) # safe list for each compound
saveRDS(f1_loss_seeds_folds, file = paste(path, "feat_imp_permutthreshseeds_", compound_name,"_", target_variable, "_", w,".RDS", sep="")) # safe list for each compound
if(folder_name =="RF_output_data"){
saveRDS(feat_imp_gini_thresh, file = paste(path, "feat_imp_ginithresh_", compound_name,"_", target_variable, "_", w,".RDS", sep="")) # safe list for each compound
saveRDS(feat_imp_gini_cor_thresh, file = paste(path, "feat_imp_gini_corthresh_", compound_name,"_", target_variable, "_", w,".RDS", sep="")) # safe list for each compound
}
} # loop seeds

} # loop compounds
