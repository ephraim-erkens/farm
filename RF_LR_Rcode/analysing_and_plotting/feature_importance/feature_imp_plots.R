### Script to analyse and plot the feature importance ranking plots -----
### 0. Load all performance results computed in script feature_imp_prep.R

# load libraries for data handling + plotting
library(dplyr)
library(ggplot2)
library(mlr3)
library(mlr3measures)
library(tidyverse)
library(RColorBrewer)

#------------------------------------------------------------------------


# set working directory (adapt to own directory!)
### 
setwd("V:/PROJEKTE/FARM/Dokumentation/Datenübergabe_UBA/ML_Code/code_paper_submission/training/RF/")
###

######################################################
### specify folder names 

folder_name = "RF_output_data" # choose folder name of models

path0 = paste(folder_name,"/",  sep="")
path_plots= ("plots/")
dir.create(path_plots)
path_plots = paste(path_plots, "feature_importance/", sep="") #
dir.create(path_plots)


######################################################

# read compound names for which models are build individually
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

target_variable = "conc_group"

## custom plotting colour code
colors_feature_categories = readRDS(file = "../colors_feature_categories_ENG.RDS")


######################################################


target_variable = "conc_group"

#### Preprocessing for plotting the ranks of the different feature importance measures for all compounds (script: plots_feature_imp_ranks_nested_cv_FJ_plain.R)
### prepare dataframes for imp, gini & gini_cor for all compounds

#--------------------
# considered compounds and feature measures for adaption
feats <- c(1)

#--------------------

# create empty dataframe of right format for later binding in for-loop
feat_imp <- function() {
  feat_imp = data.frame(
    feature = character(),
    imp = numeric(),
    targets = character(),
    compound_name = character()
  )
}

feat_imp_permut_thresh_compound_df_long <- feat_imp()
feat_imp_gini_compound_df_long <- feat_imp()
feat_imp_gini_cor_compound_df_long <- feat_imp()


### apply function defined above for each feature importance method and save dataframes

### permutation thresholding

  ### start of main loop to choose compound
  for(k in seq_along(compound_names)) {
    compound_name = compound_names[k] 
    compound_name_k = compound_names[k] 
    cat("\n", compound_name, "\n", k, "\n") # to keep track
    
    compound_name_title = compound_names_key %>% filter(lawa_name == compound_name) %>% select(engl_name)
    compound_name_title  = compound_name_title$engl_name
    
  file_names = list.files(paste(path0, compound_name, "/", sep = ""), pattern = paste("^", "feat_imp_permutthresh_", sep = ""))
  if (length(file_names) > 0){
    for (i in 1:length(file_names)) {
      df1 = as.data.frame(readRDS(paste(path0, compound_name, "/", file_names[i], sep = "")))
      names(df1) = c("imp_fold1", "imp_fold2", "imp_fold3", "imp_fold4", "imp_fold5")
      df1 <- rownames_to_column(df1, var = "feature")
      df1$targets = gsub(paste(c(paste("feat_imp_permut_thresh_", compound_name, "_", sep = ""), ".RDS"), collapse = "|"), "", file_names)[i]
      df1$compound_name <- compound_name
      df1$mean_imp_folds = apply(df1[, c("imp_fold1", "imp_fold2", "imp_fold3", "imp_fold4", "imp_fold5")], 1, mean)
      df1_long <- df1 %>% pivot_longer(cols = starts_with("imp_fold"), names_to = "Iteration", values_to = "imp") # convert to long format
      feat_imp_permut_thresh_compound_df_long = dplyr::bind_rows(feat_imp_permut_thresh_compound_df_long, df1_long)
    }
  }
}



### gini
for(k in seq_along(compound_names)) {
  
  compound_name = compound_names[k] 
  compound_name_k = compound_names[k] 
  file_names = list.files(paste(path0, compound_name, "/", sep = ""), pattern = paste("^", "feat_imp_ginithresh_", sep = ""))
  if (length(file_names) > 0){
    for (i in 1:length(file_names)) {
      df1 = as.data.frame(readRDS(paste(path0, compound_name, "/", file_names[i], sep = "")))
      names(df1) = c("imp_fold1", "imp_fold2", "imp_fold3", "imp_fold4", "imp_fold5")
      df1 <- rownames_to_column(df1, var = "feature")
      df1$targets = gsub(paste(c(paste("feat_imp_gini_", compound_name, "_", sep = ""), ".RDS"), collapse = "|"), "", file_names)[i]
      df1$compound_name <- compound_name
      df1$mean_imp_folds = apply(df1[, c("imp_fold1", "imp_fold2", "imp_fold3", "imp_fold4", "imp_fold5")], 1, mean)
      df1_long <- df1 %>% pivot_longer(cols = starts_with("imp_fold"), names_to = "Iteration", values_to = "imp") # convert to long format
      feat_imp_gini_compound_df_long = dplyr::bind_rows(feat_imp_gini_compound_df_long, df1_long)
    }
  }
}

### gini corrected
for(k in seq_along(compound_names)) {
  
  compound_name = compound_names[k] 
  compound_name_k = compound_names[k] 
  file_names = list.files(paste(path0, compound_name, "/", sep = ""), pattern = paste("^", "feat_imp_gini_corthresh_", sep = ""))
  if (length(file_names) > 0){
    for (i in 1:length(file_names)) {
      df1 = as.data.frame(readRDS(paste(path0, compound_name, "/", file_names[i], sep = "")))
      names(df1) = c("imp_fold1", "imp_fold2", "imp_fold3", "imp_fold4", "imp_fold5")
      df1 <- rownames_to_column(df1, var = "feature")
      df1$targets = gsub(paste(c(paste("feat_imp_ginicor_", compound_name, "_", sep = ""), ".RDS"), collapse = "|"), "", file_names)[i]
      df1$compound_name <- compound_name
      df1$mean_imp_folds = apply(df1[, c("imp_fold1", "imp_fold2", "imp_fold3", "imp_fold4", "imp_fold5")], 1, mean)
      df1_long <- df1 %>% pivot_longer(cols = starts_with("imp_fold"), names_to = "Iteration", values_to = "imp") # convert to long format
      feat_imp_gini_cor_compound_df_long = dplyr::bind_rows(feat_imp_gini_cor_compound_df_long, df1_long)
    }
  }
}

# final saving of dataframes
# if-loops as feature importance method do not exist for all models (especially gini and gini corrected)
if (nrow(feat_imp_permut_thresh_compound_df_long) > 1) {
  saveRDS(feat_imp_permut_thresh_compound_df_long, file = paste(path_plots, "feat_imp_permut_thresh_compound_df_long.RDS", sep=""))
}

if (nrow(feat_imp_gini_compound_df_long) > 1) {
  saveRDS(feat_imp_gini_compound_df_long, file = paste(path_plots, "feat_imp_gini_compound_df_long.RDS", sep=""))
}
if (nrow(feat_imp_gini_cor_compound_df_long) > 1) {
  saveRDS(feat_imp_gini_cor_compound_df_long, file = paste(path_plots, "feat_imp_gini_cor_compound_df_long.RDS", sep=""))
}

#------------------------------------------------------------------------
#------------------------------------------------------------------------
### 1.1 Plotting the ranks of the different feature importance measures for all compounds ----

library(rcartocolor)
# display_carto_all(colorblind_friendly = TRUE)
col_set = carto_pal(12, "Safe")
col_set2 = carto_pal(12, "Teal")

### define colors

colors_feature_categories = colors_feature_categories %>% distinct()
cols = c("meteorology" = col_set2[2], "groundwater chemistry" = col_set[3], "soil" = col_set[10], "hydrogeology" = col_set[12], "topography" = col_set[2], 
         "land use" = col_set[4], "water body" = col_set[5], "site" = col_set[6], "mohp" = "grey80", "water budget" = col_set[1])

colors_feature_categories$short_name_feature = ifelse(grepl("other agr", colors_feature_categories$short_name_feature, fixed = TRUE), "other agric. areas", colors_feature_categories$short_name_feature)
colors_feature_categories$short_name_feature = ifelse(grepl("other non-", colors_feature_categories$short_name_feature, fixed = TRUE), "other non-agric. areas", colors_feature_categories$short_name_feature)

#-------------------------------
#### Loops to generate plots for permutation each compound k and each target i ----
feat_imp_plot <- function(df_name, df_method, path_dir_method_1, path_dir_method_2, x_name_method, title) {
  
  # if-loops as feature importance method do not exist for all models (especially gini and gini corrected)
  if (file.exists(paste(path_plots, df_name,".RDS", sep=""))) {

    # read feature importance data in long format prepared in script feature_imp_prep.R and join colors to data.frame by matching "feature"
    df_method = readRDS(file = paste(path_plots, df_name,".RDS", sep=""))
    df_method = left_join(df_method, colors_feature_categories)
    
    filtered_df = df_method %>% ungroup()
    
    # define xlim based on minimum and maximum values of averaged folds across all compounds for better comparison (including standard deviation)
    filtered_df_imp_mean_sd = filtered_df %>%
      group_by(short_name_feature, compound_name) %>%
      summarize(mean_imp = mean(mean_imp_folds), sd_imp_low = mean(mean_imp_folds) - sd(mean_imp_folds), sd_imp_high = mean(mean_imp_folds) + sd(mean_imp_folds))
    axis = c(min(filtered_df_imp_mean_sd$sd_imp_low), max(filtered_df_imp_mean_sd$sd_imp_high))
    
    
    
    for(k in seq_along(compound_names)) {
      compound_name = compound_names[k] 
      compound_name_k = compound_names[k] 
      print(k)

      for(i in feats){
        
        ### try to put variability of folds into a single plot by horizontal boxplots
        filtered_df = df_method %>% ungroup() %>% filter(compound_name == compound_name_k)
     
        filtered_df_mean_all = filtered_df %>%
          group_by(short_name_feature) %>%
          summarize(imp = mean(imp)) %>% arrange(desc(imp)) %>% arrange(imp)
        filtered_df_all =  filtered_df %>%
          filter(short_name_feature %in% filtered_df_mean_all$short_name_feature)
        filtered_df_all$short_name_feature = factor(filtered_df_all$short_name_feature, levels = filtered_df_mean_all$short_name_feature)
        dir.create(paste0(path_plots, path_dir_method_1))
        dir.create(paste0(path_plots, path_dir_method_1, path_dir_method_2))
       
        # filter for the 20 features with largest mean
        filtered_df_mean = filtered_df %>%
          group_by(short_name_feature) %>%
          summarize(imp = mean(imp)) %>% arrange(desc(imp)) %>% slice_head(n = 20) %>% arrange(imp)
        # only use Top20 for plotting
        filtered_df_20 =  filtered_df %>%
          filter(short_name_feature %in% filtered_df_mean$short_name_feature)
        # arrange order
        filtered_df_20$short_name_feature = factor(filtered_df_20$short_name_feature, levels = filtered_df_mean$short_name_feature)
      
  
        # calculate mean and standard deviation of all seeds per feature (measures the variability between the seeds)
        filtered_df_20_summary_folds = filtered_df_20 %>%
          dplyr::group_by(short_name_feature, category) %>%
          dplyr::summarize(mean_imp_folds = mean(imp),
                    sd_imp_folds = sd(imp))        
   
        
        #-------------------------------
        
        
        ### do the same plot for different seeds
        # therefor make sure each mean of a seed i only considered once
        # remove columns for folds and delete duplicates
        filtered_df_all_seeds = filtered_df_all[, !(names(filtered_df_20) %in% c("imp", "Iteration"))]
        filtered_df_all_seeds = filtered_df_all_seeds %>% distinct()
        filtered_df_20_seeds = filtered_df_20[, !(names(filtered_df_20) %in% c("imp", "Iteration"))]
        filtered_df_20_seeds = filtered_df_20_seeds %>% distinct()
        dir.create(paste0(path_plots, path_dir_method_1, path_dir_method_2, "plots_seeds/"))
     
        # calculate mean and standard deviation of all seeds per feature (measures the variability between the seeds)
        filtered_df_all_summary = filtered_df_all_seeds %>%
          group_by(short_name_feature, category) %>%
          summarize(mean_imp_seeds = mean(mean_imp_folds),
                    sd_imp_seeds = sd(mean_imp_folds))
        
        # Plot the barplot with error bars
        bxplt_seeds2 = 
          ggplot(filtered_df_all_summary, aes(y = short_name_feature, x = mean_imp_seeds)) +
          geom_bar(aes(fill = category), stat = "identity", color = "black") +
          geom_errorbar(aes(xmin = mean_imp_seeds - sd_imp_seeds,
                            xmax = mean_imp_seeds + sd_imp_seeds), width = 0.2, color = "black") +
          scale_fill_manual(values = cols) +
          labs(fill = "feature category", y = "", x = x_name_method, title = compound_name_title) + #subtitle = title
          geom_text(aes(label = short_name_feature, x = - 0.12 * (axis[2] - axis[1])), hjust = 1, size = 4) +
          coord_cartesian(xlim = axis, clip = "off") +
          theme_bw() + 
          theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 12),
                axis.text = element_text(size = 12), axis.title = element_text(size = 12),
                axis.text.y = element_blank(),    
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 170, unit = "pt"),
                legend.text = element_text(size = 11),        # Increase legend text size
                legend.title = element_text(size = 12),      # Increase legend title size
                legend.key.size = unit(0.7, "cm"))
        ggsave(file = paste(path_plots, path_dir_method_1, path_dir_method_2, "plots_seeds/", "all_", compound_name_k, "_", target_variable, ".png", sep = ""), bxplt_seeds2, width = 8, height = 10)
        
       
        # calculate mean and standard deviation of all seeds per feature (measures the variability between the seeds)
        filtered_df_20_summary = filtered_df_20_seeds %>%
          group_by(short_name_feature, category) %>%
          summarize(mean_imp_seeds = mean(mean_imp_folds),
                    sd_imp_seeds = sd(mean_imp_folds))
        
        # Plot the barplot with error bars
        bxplt_seeds4 = 
          ggplot(filtered_df_20_summary, aes(y = short_name_feature, x = mean_imp_seeds)) +
          geom_bar(aes(fill = category), stat = "identity", color = "black") +
          geom_errorbar(aes(xmin = mean_imp_seeds - sd_imp_seeds,
                            xmax = mean_imp_seeds + sd_imp_seeds), width = 0.2, color = "black") +
          scale_fill_manual(values = cols) +
          labs(fill = "feature category", y = "", x = x_name_method, title = compound_name_title) + #subtitle = title
          geom_text(aes(label = short_name_feature, x = - 0.12 * (axis[2] - axis[1])), hjust = 1, size = 4) +
          coord_cartesian(xlim = axis, clip = "off") +
          theme_bw() + 
          theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 12),
                axis.text = element_text(size = 12), axis.title = element_text(size = 12),
                axis.text.y = element_blank(),    
                plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 170, unit = "pt"),
                legend.text = element_text(size = 11),        # Increase legend text size
                legend.title = element_text(size = 12),      # Increase legend title size
                legend.key.size = unit(0.7, "cm"))
        ggsave(file = paste(path_plots, path_dir_method_1, path_dir_method_2, "plots_seeds/", "top20_", compound_name_k, "_", target_variable, ".png", sep = ""), bxplt_seeds4, width = 8, height = 10)
        saveRDS(bxplt_seeds4,file = paste(path_plots, path_dir_method_1, path_dir_method_2, "plots_seeds/", "top20_", compound_name_k, "_", target_variable, "_PLOT.RDS", sep = ""))
        

      }
      }
  }
}

# call function for each method
feat_imp_plot("feat_imp_permut_thresh_compound_df_long", feat_imp_permut_thresh_compound_df_long, "", "permutation_thresh/", "F1-Score Permutation Importance", "Permutation Importance")
feat_imp_plot("feat_imp_gini_compound_df_long", feat_imp_gini_compound_df_long, "", "gini/", "Gini Importance", "Gini Importance")
feat_imp_plot("feat_imp_gini_cor_compound_df_long", feat_imp_gini_cor_compound_df_long, "", "gini_cor/", "Gini Importance Corrected", "Gini Importance Corrected")



# combine individual compound plots into a grid
list.plots = list.files(paste(path_plots,"permutation_thresh/plots_seeds/", sep=""), pattern="PLOT.RDS", full.names = T) #permutation_thresh

list.plots
k_10_plot =readRDS(list.plots[1])
k_26_plot =readRDS(list.plots[2])
k_27_plot =readRDS(list.plots[3])
k_28_plot =readRDS(list.plots[4])
k_30_plot =readRDS(list.plots[5])
k_31_plot =readRDS(list.plots[6])
k_38_plot =readRDS(list.plots[7])
k_40_plot =readRDS(list.plots[8])

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


k_10_plot = k_10_plot + ggtitle(firstup(compound_names_key$engl_name[1]))
k_26_plot = k_26_plot + ggtitle(firstup(compound_names_key$engl_name[2]))
k_27_plot = k_27_plot + ggtitle(firstup(compound_names_key$engl_name[3]))
k_28_plot = k_28_plot + ggtitle(firstup(compound_names_key$engl_name[4]))
k_30_plot = k_30_plot + ggtitle(firstup(compound_names_key$engl_name[5]))
k_31_plot = k_31_plot + ggtitle(firstup(compound_names_key$engl_name[6]))
k_38_plot = k_38_plot + ggtitle(firstup(compound_names_key$engl_name[7]))
k_40_plot = k_40_plot + ggtitle(firstup(compound_names_key$engl_name[8]))

r = ggarrange(k_10_plot, k_28_plot, k_26_plot,  k_27_plot,  k_30_plot,  k_31_plot, k_38_plot,k_40_plot,
              ncol=2, nrow=4, common.legend = T, legend="bottom",
              widths = c(rep(1,7),2), heights = c(rep(1,7),2))

r

ggsave(file = paste(path_plots,"permutation_thresh/plots_seeds/plots_permutation_importance.png", sep = ""), r, width = 12, height = 15)



