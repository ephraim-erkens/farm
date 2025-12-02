#### 

library(tidyverse)
library(dplyr)
library(sf)
library(ggplot2)
library(svglite)
library(rcartocolor)
library(ggpubr)
library(ggExtra)

library(mlr3)
library(patchwork)

library(RColorBrewer)
library(ggpubr)

# # Load and prepare data including information about censoring -----


# set working directory (adapt to own directory!)
### 
setwd("V:/PROJEKTE/FARM/Dokumentation/Datenübergabe_UBA/ML_Code/code_paper_submission/training/RF/")
###

######################################################
### specify folder names and plot titles

folder_name = "RF_output_data" # choose folder name of models


path0 = paste(folder_name,"/",  sep="")
path_plots= ("plots/")
dir.create(path_plots)
path_plots = paste(path_plots, "residuals/", sep="") #
dir.create(path_plots)



### 1
##retrieve compound group names ----
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

colors_feature_categories = readRDS( file="../colors_feature_categories_ENG.RDS")

colors_feature_categories= colors_feature_categories %>% select(-c(abbr, kultur_deutsch, kultur_englisch, kultur_sentinel_id)) %>% distinct()

colors_feature_categories$variable = colors_feature_categories$feature
colors_feature_categories$feature = NULL


target_variable = "conc_group"


col_set3 = brewer.pal(n = 10, name = "Paired") # color ramp
col_set4 = brewer.pal(n = 11, name = "RdYlBu") # color ramp
col_set5 = brewer.pal(n = 8, name = "Set1") # color ramp
col_set6 = brewer.pal(n = 9, name = "YlOrBr") # color ramp
col_set7 = brewer.pal(n = 8, name = "PiYG") # color ramp


### Which compounds are used for which sentinel_id? # 
### for information about compound application
compounds_crops = readRDS("../compounds_crops.RDS")

# create dataframe with main landuses for most important compounds
main_landuses <- compounds_crops %>%
  filter(lawa_name != lag(lawa_name, default = first(lawa_name)))
main_landuses = rbind(compounds_crops[1, ], main_landuses) # otherwise first row not added
main_landuses[main_landuses$lawa_name == "s-metolachlor-metabolit cga 357704", 2] = "Weizen" # fill NA vlaue
main_landuses[main_landuses$lawa_name == "s-metolachlor-metabolit cga 357704", 3] = "corn" # fill NA vlaue
main_landuses$sentinel_area = NA
main_landuses[main_landuses$kultur_filter_englisch == "wheat", 4] =  1 # fill NA vlaue
main_landuses[main_landuses$kultur_filter_englisch == "sugarbeet", 4] =  14 # fill NA vlaue
main_landuses[main_landuses$kultur_filter_englisch == "rape seed", 4] =  12 # fill NA vlaue
main_landuses[main_landuses$kultur_filter_englisch == "corn", 4] =  10 # fill NA vlaue


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

for(k in seq_along(compound_names)) {


  compound_name = compound_names[k] 
  
  path = paste(path0,compound_name, "/", sep="")
  
  if (compound_name == "s-metolachlor-metabolit noa 413173"){
    compound_name_v2= "noa 413173" } else {
      compound_name_v2 = compound_name
    }
  print(compound_name) # to keep track
  print(k) # to keep track
  sent_id = main_landuses %>% filter(lawa_name == compound_name )
  

  hp_list =list()
  for(i in 1:10){
    hp = readRDS(file=paste(path, "HP_thresh_max_perf_final_", compound_name,"_", "conc_group","_",i,".RDS", sep=""))
    hp = as.data.frame(hp)
    hp$measure = row.names(hp)
    row.names(hp) = NULL
    hp$seed_nr = i
    hp_list[[i]] = hp
  }
  
  hp_df =bind_rows(hp_list)
  head(hp_df )
  max_hp= hp_df %>% filter(measure == "macro_f1_test_outer")
  max_hp= hp_df %>% filter(hp== max(hp, na.rm = TRUE))
  w = unique(max_hp$seed_nr)
  
  ### read in thresholds for this seed number

  
  rr = readRDS(file=paste(path, "rr_outer_list_", compound_name,"_", target_variable,"_",w,".RDS", sep=""))
  pred_outer_test_list = readRDS(file=paste(path, "pred_outer_test_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))
  outer_test_task_list = readRDS(file=paste(path, "outer_test_task_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))
  
  tsk = readRDS(file=paste(path, "task_", compound_name,"_", target_variable,".RDS", sep=""))
  df_data_final = readRDS(file=paste(path,"df_data_final_", compound_name,"_", target_variable,".RDS", sep=""))
  df_data_final$corg_gehalt_depth_weighted_mean =   df_data_final$corg_gehalt_depth_weighted_mean/(1.72^2)
  train_full = readRDS(file=paste(path,"full_model_", compound_name,"_", target_variable,"_", w,".RDS", sep=""))
  
  # convert area from m2 to % within buffer
  df_data = df_data_final %>%
    mutate(across(c(matches("area")), ~ . / 10000)) #(pi*1000^2/100)))
  df_data$lawa_name = compound_name
  ### 1
  ### read in predictions
  # stacking of classification results of different folds for more meaningfull analysis of accuracy
  test_perf_folds_df = data.frame(
    truth = numeric(),
    response = numeric(),
    row_ids = numeric(),
    prob1 = numeric(),
    prob2 = numeric(), #,
  prob3 = numeric())
  # extract splits during cv
  # inst = rr$resampling$instance
  test_perf_folds = test_perf_folds_df
  # write ground truth and predictions for all folds in long format
  for (j in 1:length(rr)){

    outer_test_inst = outer_test_task_list[[j]]
    
    rr_outer_test = rr[[j]]
    prediction = pred_outer_test_list[[j]] #rr_outer_test$predict(rr_outer_test$task, row_ids =  outer_test_inst)
    df1 = as.data.frame(bind_cols(prediction$truth, prediction$response, prediction$row_ids, prediction$prob))
    names(df1) = c("truth", "response", "row_ids", paste0("prob", 1:(length(df1)-3)))
    df1$truth = as.numeric(df1$truth)
    df1$response = as.numeric(df1$response)
    test_perf_folds = dplyr::bind_rows(test_perf_folds, df1) 
  }
  test_perf_folds = test_perf_folds %>%
    mutate_at(vars("truth", "response"), as.numeric)
  test_perf_folds$truth_minus_response = test_perf_folds$truth - test_perf_folds$response
  test_perf_folds = test_perf_folds %>%
    mutate_at(vars("truth", "response", "truth_minus_response"), as.factor)
  df_data$row_ids = 1:nrow(df_data)
  ## add the prediction to the original dataframe (which includes the site id etc)
  data_pred_df =full_join(df_data, test_perf_folds, by=c("row_ids") )
  data_pred_df$Makroporen_rounded = as.factor(data_pred_df$Makroporen_rounded)
  
  
  best_thresh = readRDS(file=paste(path, "best_HP_f1_thresh_max_single_cv_", compound_name,"_", "conc_group","_",w,".RDS", sep=""))
  


  exclude_columns <- c(
    "prob1", "prob2", "prob3", 
    "prob_thresh_ratio_1", "prob_thresh_ratio_2", "prob_thresh_ratio_3", 
    "ratio_1_0.8", "ratio_2_0.8", "ratio_3_0.8"
  )
  
  df_long <- data_pred_df %>%
    ungroup() %>%
    filter(lawa_name == compound_name) %>%
    pivot_longer(
      cols = where(is.numeric) & setdiff(names(data_pred_df), exclude_columns), # Dynamically exclude columns
      values_to = "value"
    )
  

  head(df_long)
  numeric_cols <- data_pred_df %>%
    select_if(~ is.numeric(.)) %>%
    colnames()
  


  df_long$variable = df_long$name
  df_long$name= NULL
  df_long = left_join(df_long, colors_feature_categories)
  

  df_long = left_join(df_long,  compound_names_key )
  
  
  
  vars  = numeric_cols
  
  
  colors_3_classes = c("#A50026", "#F46D43", "#FFFFBF", "#74ADD1", "#313695")
  # colors_binary = c( "#F46D43", "#FFFFBF", "#74ADD1")
  colors_binary = c( "#A50026", "#FFFFBF", "#313695")
  

  
  library(viridis)
  colors_3_classes = viridis_pal(option = "D")(5) # Option D for high contrast
  colors_binary = viridis_pal(option = "D")(3) # Option D for high contrast
  
  
    fill_colors = colors_binary
    fill_colors <- c("-1" = fill_colors[1], 
                     "0" = fill_colors[2], 
                     "1"  = fill_colors[3],
                     "All" = "grey")

  
  legend_pos="NONE" # bottom
  
 
  custom_labels <- c(
    "-1" = "overestimation",
    "0" = "correct",
    "1" = "underestimation",
    "All" = "all") 
    vars <- setdiff(vars, c("x", "y" , "row_ids", "prob1", "prob2", "prob3", "prob_thresh_ratio_1", "prob_thresh_ratio_2", "prob_thresh_ratio_3","ratio_1_0.8","ratio_2_0.8", "ratio_3_0.8"))


  for(var in vars){

    df_new = df_long %>% filter(variable %in% var, lawa_name == compound_name)  
    df_new_density <- df_new %>%
      group_by(truth_minus_response) %>%
      do({
        # Calculate the density for each group
        dens <- density(.$value)
        data.frame(x = dens$x, y = dens$y)
      })
    
    # Step 2: Find the maximum density value (width) for each group
    df_new_density_summary <- df_new_density %>%
      group_by(truth_minus_response) %>%
      summarize(max_density = max(y)) %>%
      arrange(max_density)  # Arrange by max density
    
    # Step 3: Reorder the factor levels based on max_density (widest first)
    df_new$truth_minus_response <- factor(df_new$truth_minus_response, levels = df_new_density_summary$truth_minus_response)
    
    df_new$alpha_values = ifelse(df_new$truth_minus_response =="0", 1, NA)
    df_new$alpha_values = ifelse(df_new$truth_minus_response =="1", 1,  df_new$alpha_values)
    df_new$alpha_values = ifelse(df_new$truth_minus_response =="-1", 1,  df_new$alpha_values)
    
    data_summary <- df_new  %>% group_by(truth_minus_response) %>% summarize(count = n())
    custom_labels_add =  custom_labels
    df_names = data.frame(matrix(NA, ncol=2, nrow=4))
    names(df_names ) = c("truth_minus_response", "name")
    
    df_names$truth_minus_response= names(custom_labels)
    df_names$name =  custom_labels
    
    data_summary = left_join(data_summary, df_names)
    custom_labels_add =      custom_labels
    names(custom_labels_add) = data_summary$truth_minus_response
    custom_labels_add[1] <- paste(data_summary$name[1],"\n n = ",data_summary$count[1], sep="" )
    custom_labels_add[2] <- paste(data_summary$name[2], "\n n = ",data_summary$count[2], sep="" )
    custom_labels_add[3] <- paste(data_summary$name[3], "\n n = ",data_summary$count[3], sep="" )
    custom_labels_add
    df_new<- df_new %>%
      mutate(facet_group = as.factor(truth_minus_response)) #%>%
      #bind_rows(df_new %>% mutate(facet_group = "All"))  # Add "All" category
    df_new$facet_group = factor( df_new$facet_group, levels=c("-1", "0", "1"))
  
      plot1 = ggplot(data=df_new, aes(x=1, y=value)) +  
       geom_violin(aes(  fill=facet_group, alpha= alpha_values), color="black", scale="count", position = "identity") + ylab("") + xlab("") + #,
       # facet_wrap(~truth_minus_response) +#variable #lawa_name
       theme_minimal() + theme(axis.text.x = element_blank(),text = element_text(size = 22),
                               strip.text = element_text(size = 22),axis.text = element_text(size = 22),
                               panel.grid.major.y = element_line(color = "grey", linewidth = 0.3),
                               panel.grid.minor.y = element_line(color = "grey", linewidth = 0.3), legend.position = legend_pos) +
        scale_alpha(guide = "none") +
         facet_wrap(~ facet_group, nrow = 1,  labeller = as_labeller(custom_labels_add)) +
    
        scale_fill_manual(values=fill_colors,  labels = as_labeller(custom_labels_add)) + 
        labs(title = firstup(unique(df_new$engl_name)) , x="", fill = "", y = paste(unique(df_new$short_name_feature)," [",unique(df_new$unit),"]", sep=""))
  

      saveRDS(plot1,file=paste(path_plots,compound_name_v2,"_", var, ".RDS", sep=""))
      
    ggsave(plot1,file=paste(path_plots,compound_name_v2,"_", var, ".png", sep=""), width=10, height=10, dpi = 300)
    
  }  
  
}  


# vars
# compound_names_key
mtz_rape = readRDS(file=paste(path_plots,"metazachlor esa","_","area_sentinel_id_12" , ".RDS", sep=""))
dc_sugar = readRDS(file=paste(path_plots,"desphenyl-chloridazon","_","area_sentinel_id_14" , ".RDS", sep=""))
smoc_maize = readRDS(file=paste(path_plots,"metolachlor esa","_","area_sentinel_id_10" , ".RDS", sep=""))
smoc_well = readRDS(file=paste(path_plots,"metolachlor-ca","_","filter_ok_unter_gok" , ".RDS", sep=""))
mtz_well = readRDS(file=paste(path_plots,"metazachlorsäure","_","filter_ok_unter_gok"  , ".RDS", sep=""))

combined_plot = ggarrange(mtz_rape, dc_sugar, smoc_maize, smoc_well, mtz_well,
              ncol=3, nrow=2, labels=c("a)", "b)", "c)", "d)", "e)"),
              font.label = list(size = 22, color = "black", face="plain"))#, #common.legend = T, legend="bottom",
             # widths = c(rep(1,7),2), heights = c(rep(1,7),2))

combined_plot
ggsave(combined_plot,file=paste(path_plots, "fig07.png", sep=""), width=24, height=15, dpi = 300)

