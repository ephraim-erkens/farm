#### script for feature effect (partial dependence plots) calculations and creation of plots


library(mlr3)
library(mlr3learners)
library(mlr3measures)
library(rpart)
library(tidyverse)
library(future)
library(future.apply)
library(iml)
library(rcartocolor)
library(ggtext)
library(patchwork)
library(ggplot2)


setwd("V:/PROJEKTE/FARM/Dokumentation/Datenübergabe_UBA/ML_Code/code_paper_submission/training/RF/")
###

######################################################
### specify folder names 

folder_name = "RF_output_data" # choose folder name of models


path0 = paste(folder_name,"/",  sep="")
path_plots= ("plots/")
dir.create(path_plots)
path_plots = paste(path_plots, "feature_effects/", sep="") #
dir.create(path_plots)

dir.create(paste(path_plots, "pdp/", sep=""), showWarnings = F)
dir.create(paste(path_plots, "pdp_ice/", sep=""), showWarnings = F)
dir.create(paste(path_plots, "ale/", sep=""), showWarnings = F)


### 1
##retrieve compound group names ----
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

colors_feature_categories = readRDS( file="../colors_feature_categories_ENG.RDS")

colors_feature_categories= colors_feature_categories %>% select(-c(abbr, kultur_deutsch, kultur_englisch, kultur_sentinel_id)) %>% distinct()


target_variable = "conc_group"


# load color blind friendly color ramp
#display_carto_all(colorblind_friendly = TRUE)
col_set = carto_pal(12, "Safe")
col_set2 = carto_pal(12, "Teal")

# Loop for feature effects calculation (caculations may take long)
for(k in seq_along(compound_names)) {
  
  compound_name = compound_names[k] 
  print(compound_name)
  path = paste(path0, compound_name, "/", sep="")
  path_plots1 = paste(path_plots, compound_name, "/", sep="")
  
  # read in the tsk
  tsk = readRDS(file=paste(path,"task_",compound_name, "_", target_variable,".RDS", sep=""))
  # features in train data
  data_x = tsk$data(cols = tsk$feature_names)
  # target in train data
  data_y = tsk$data(cols = tsk$target_names)

  ## take top 20 features from permutation importance ---
  feat_imp_permut_compound_df_long = readRDS(file = paste(path_plots, "feat_imp_permut_thresh_compound_df_long.RDS", sep = ""))

  compound_name_k = filter_prio[k] ###

  feat_imp_permut_top20 = feat_imp_permut_compound_df_long %>%
    dplyr::filter(compound_name == compound_name_k) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(mean_imp_folds = mean(imp)) %>% # Use na.rm to ignore NAs
    arrange(desc(mean_imp_folds)) %>%
    slice_head(n = 20)

  ### Feature effects of top20 permutation features ----
  pdp_data_list = list()
  for (w in c(1:10)){

    model = readRDS(file=paste(path, "full_model_", compound_name,"_", target_variable, "_", w,".RDS", sep=""))
    set.seed(123)
    predictor = Predictor$new(model,  data = data_x, y = data_y) #tsk,
    plan(sequential) # Disable parallel processing
    effs_pdp = FeatureEffects$new(predictor, features = feat_imp_permut_top20$feature, method="pdp")

    pdp_data_list[[w]] = effs_pdp$results

  }

  # save pdp effects for each compound list seperately
  saveRDS(pdp_data_list, file = paste(path, "pdp_data_list_", compound_name, "_", target_variable, ".RDS", sep = ""))
}

# start new loop to separate preprocessing from plotting
for(k in seq_along(compound_names)) {
  k=1
  compound_name = compound_names[k] 
  print(compound_name)
  path = paste(path0, compound_name, "/", sep="")

  pdp_data_list = readRDS(file = paste(path, "pdp_data_list_", compound_name, "_", target_variable, ".RDS", sep = "")) 
  df_final = readRDS(file = paste(path, "df_data_final_", compound_name, "_", target_variable, ".RDS", sep = ""))
  df_final$corg_gehalt_depth_weighted_mean =   df_final$corg_gehalt_depth_weighted_mean/(1.72^2)
  # bind respective dataframes stored in different sublists of a list together in long format
  pdp_data_list_long = lapply(1:length(pdp_data_list[[1]]), function(i) {
    Reduce(bind_rows, lapply(pdp_data_list, `[[`, i))  # Bind respective data frames
  })
  
  ### merge with other data frame defining feature categories and colors
  cols = c("meteorology" = col_set2[2], "groundwater chemistry" = col_set[3], "soil" = col_set[10], "hydrogeology" = col_set[12], "topography" = col_set[2], 
           "land use" = col_set[4], "water body" = col_set[5], "site" = col_set[6], "mohp" = "grey80", "water budget" = col_set[1])
  
   cols = setNames(paste0(cols, "80"), names(cols)) # add transparency
  
  # Apply left_join to each data frame in the list
  colors_feature_categories = colors_feature_categories %>%
    dplyr::rename(.feature = feature)
  pdp_data_list_long = lapply(pdp_data_list_long, function(df) {
    dplyr::left_join(df, colors_feature_categories, by = ".feature")
  })
  
  # loop through all data frames in the list to find the maximum value in the specified column across all data frames
  max_value = max(sapply(pdp_data_list_long, function(df) max(df[[".value"]], na.rm = TRUE)), na.rm = TRUE)
  min_value = min(sapply(pdp_data_list_long, function(df) min(df[[".value"]], na.rm = TRUE)), na.rm = TRUE)
  
  pdp_plot_list = list()
  # Loop through each data frame in the list
  for (i in seq_along(pdp_data_list_long)) {
    
    pdp_feat = pdp_data_list_long[[i]]
    if(pdp_feat$.feature[1] == "corg_gehalt_depth_weighted_mean"){
      pdp_feat$.borders= pdp_feat$.borders/(1.72^2)
    }
    # calculate mean and standard deviation of all seeds per feature (measures the variability between the seeds)
    pdp_feat_summary = pdp_feat %>%
      group_by(.borders, .class, .feature, short_name_feature, category, unit) %>%
      dplyr::summarize(mean_value = mean(.value),
                       sd_value = sd(.value))
    
    # convert units for land use features (m2 to ha)
    if(unique(pdp_feat_summary$category) == "land use"){
      pdp_feat_summary$.borders = pdp_feat_summary$.borders / 10000
      df_final[[unique(pdp_feat_summary$.feature)]] = df_final[[unique(pdp_feat_summary$.feature)]] / 10000
    }
    # convert units for groundwater chemistry features (ug to mg)
    if(unique(pdp_feat_summary$category) == "gw chemistry"){
      pdp_feat_summary$.borders = pdp_feat_summary$.borders / 1000
      df_final[[unique(pdp_feat_summary$.feature)]] = df_final[[unique(pdp_feat_summary$.feature)]] / 1000
    }
    # convert units for LP_, DSD_ and SD_ distance_m features (m to km)
    if(unique(pdp_feat_summary$.feature) %in% c("LP_1", "LP_2", "SD_1", "SD_2", "DSD_1", "DSD_2", "distance_m")){
      pdp_feat_summary$.borders = pdp_feat_summary$.borders / 1000
      df_final[[unique(pdp_feat_summary$.feature)]] = df_final[[unique(pdp_feat_summary$.feature)]] / 1000
    }
    
    # convert units for LP_, DSD_ and SD_ distance_m features (m to km)
    # if(unique(pdp_feat_summary$.feature) %in% c("HA")){
    #   pdp_feat_summary$.borders = ifelse(pdp_feat_summary$.borders=="K/P", "F/P", pdp_feat_summary$.borders)
    #   pdp_feat_summary$.borders = ifelse(pdp_feat_summary$.borders=="K/Ka", "F/Ka", pdp_feat_summary$.borders)
    #   df_final$HA = ifelse(df_final$HA =="K/P", "F/P",df_final$HA )
    #   df_final$HA = ifelse(df_final$HA =="K/Ka", "F/Ka",df_final$HA )
    #   
    # }
    
    if(is.numeric(pdp_feat_summary$.borders) == TRUE){
      
      # create a separate data frame for the rug plot with the same names as in pdp_feat_summary to subset data for the class, respectively
      df_rug = data.frame(cbind(df_final[[unique(pdp_feat_summary$.feature)]], df_final$conc_group))
      names(df_rug) = c(unique(pdp_feat_summary$.feature), ".class")
      df_rug_median = df_rug %>%
        group_by(.class) %>%
        summarise(median_value = median(!!sym(unique(pdp_feat_summary$.feature))))

      pdp_plot =
        ggplot() +
        geom_ribbon(pdp_feat_summary, mapping = aes(x = .borders, ymin = mean_value - sd_value, ymax = mean_value + sd_value), alpha = 0.6, fill = "lightblue") +
        geom_line(pdp_feat_summary, mapping = aes(x = .borders, y = mean_value), stat = "identity", color = "darkblue", lwd = 0.8) +
        geom_rug(data = df_rug, aes(x = !!sym(unique(pdp_feat_summary$.feature)), y = 0), sides = "b", color = "darkblue", size = 0.3, alpha = 0.05, inherit.aes = FALSE) +
        #geom_rug(data = data.frame(x = df_final[[unique(pdp_feat_summary$.feature)]]), aes(x = x, y = 0), sides = "b", color = "darkblue", size = 0.3, alpha = 0.05, inherit.aes = FALSE) +
        geom_rug(data = df_rug_median, aes(x = median_value, y = 0), sides = "b", color = "red", size = 1, inherit.aes = FALSE) +
        facet_wrap(~ .class) +
        scale_y_continuous(limits = c(min_value, max_value)) + 
        labs(y = "probability", x = paste0(pdp_feat_summary$short_name_feature, " [", unique(pdp_feat_summary$unit), "]")) +
        theme_bw() +
        theme(strip.background = element_rect(fill = cols[unique(pdp_feat_summary$category)]),
              strip.text = element_markdown(face = "bold", color = "black", size = 11),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 12))
      
    } else if(is.factor(pdp_feat_summary$.borders) == TRUE) {
      
      # define bar colors dependent on how many feature categories exist
      if (length(unique(pdp_feat_summary$.borders)) == 3){
        col_safe = col_set[c(7, 8, 9)]
        } else if (length(unique(pdp_feat_summary$.borders)) == 4){
        col_safe = col_set[c(1, 5, 3, 4)]
        } else {
          col_safe = col_set
        }
      
      # create a separate data frame for the jitter plot with the same names as in pdp_feat_summary to subset data for the class, respectively
      df_jitter = data.frame(cbind(df_final[[unique(pdp_feat_summary$.feature)]], df_final$conc_group))
      names(df_jitter) = c(unique(pdp_feat_summary$.feature), ".class")
      
      pdp_plot =
        ggplot(pdp_feat_summary, aes(x = .borders, y = mean_value)) +
        geom_bar(aes(fill = .borders), stat = "identity", alpha = 0.7) +
        scale_fill_manual(values = col_safe) +
        geom_errorbar(aes(ymin = mean_value - sd_value,
                          ymax = mean_value + sd_value), width = 0.3, color = "black") +
        # geom_jitter(data = df_jitter, aes(x = !!sym(unique(pdp_feat_summary$.feature)), y = min_value - 0.02), height = 0, color = "darkblue", alpha = 0.025, shape = 3, size = 2.5) +
        geom_jitter(data = data.frame(x = df_final[[unique(pdp_feat_summary$.feature)]]), aes(x = x, y = min_value - 0.02), height = 0, color = "darkblue", alpha = 0.025, shape = 3, size = 2.5) +
        facet_wrap(~ .class) + 
        #scale_y_continuous(limits = c(0, max_value)) +
        coord_cartesian(ylim = c(min_value, max_value)) +
        labs(y = "probability", x = paste0(pdp_feat_summary$short_name_feature, " [", unique(pdp_feat_summary$unit), "]")) +
        theme_bw() +
        theme(legend.position = "none",
              strip.background = element_rect(fill = cols[unique(pdp_feat_summary$category)]),
              strip.text = element_markdown(face = "bold", color = "black", size = 11),
              axis.text = element_text(size = 11),
              axis.title = element_text(size = 12))

    }
    pdp_plot_list[[i]] = pdp_plot
  }
  
  # adjust layout based on whether 2 or 3 class problem

    if (length(pdp_plot_list) > 10) { # reduce the length of the list to 15
      pdp_plot_list = pdp_plot_list[1:10]
    }
  pdp_plot_wrap = wrap_plots(pdp_plot_list, ncol = 5)
  plot_height = ceiling(length(pdp_plot_list) / 5) * 10 #12

  # save plots while making sure single plots always have same height
  ggsave(paste(path_plots,"pdp/", compound_name, "_pdp_top20permut_", target_variable,".png", sep=""), pdp_plot_wrap, units = "cm", dpi = 300, height = plot_height, width = 40) 
  
  

}
