#### script for feature effect plots

# library(tidyverse)
# library(dplyr)

library(mlr3)
library(mlr3learners)
library(mlr3measures)
# library(mlr3tuning)
# library(mlr3viz)
# library(mlr3mbo)
library(rpart)
library(tidyverse)
# library(mlr3verse)
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
### specify folder names and plot titles

folder_name = "RF_output_data" # choose folder name of models


path0 = paste(folder_name,"/",  sep="")
path_plots= ("plots/")
dir.create(path_plots)
path_plots = paste(path_plots, "/feature_interactions/", sep="") #
dir.create(path_plots)



### 1
##retrieve compound group names ----
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

colors_feature_categories = readRDS( file="../colors_feature_categories_ENG.RDS")

colors_feature_categories= colors_feature_categories %>% select(-c(abbr, kultur_deutsch, kultur_englisch, kultur_sentinel_id)) %>% distinct()
# 
# colors_feature_categories$variable = colors_feature_categories$feature
# colors_feature_categories$feature = NULL


target_variable = "conc_group"



# load color blind friendly color ramp
#display_carto_all(colorblind_friendly = TRUE)
col_set = carto_pal(12, "Safe")
col_set2 = carto_pal(12, "Teal")


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
  feat_imp_permut_compound_df_long = readRDS(file = paste(path_plots1, "feat_imp_permut_thresh_compound_df_long.RDS", sep = ""))
  
  compound_name_k = compound_names[k] 
  
  feat_imp_permut_top20 = feat_imp_permut_compound_df_long %>%
    dplyr::filter(compound_name == compound_name_k) %>%
    dplyr::group_by(feature) %>%
    dplyr::summarize(mean_imp_folds = mean(imp)) %>% # Use na.rm to ignore NAs
    arrange(desc(mean_imp_folds)) %>%
    slice_head(n = 20)
  
  ### Feature effects of top20 permutation features ----
 interact_data_list = list()
 
  for (w in c(1:10)){
print(paste("fold w=", w))
    model = readRDS(file=paste(path, "full_model_", compound_name,"_", target_variable, "_", w,".RDS", sep=""))
    set.seed(123)
    predictor = Predictor$new(model, data = data_x, y = data_y) #tsk,
    plan(sequential) # Disable parallel processing
    # effs_pdp = FeatureEffects$new(predictor, features = feat_imp_permut_top20$feature, method="pdp")
    interact <- Interaction$new(predictor, grid.size = 15) #features = feat_imp_permut_top20$feature

   interact_data_list[[w]] = interact$results 
   
 
  }
  
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
  
  model = readRDS(file=paste(path, "full_model_", compound_name,"_", target_variable, "_", ww,".RDS", sep=""))
  set.seed(123)
  predictor = Predictor$new(model, data = data_x, y = data_y) #tsk,
  plan(sequential) 
  
  interact_feat_list = list()
  for(l in 1:length(feat_imp_permut_top20$feature)){
    print(paste("feat l = ", l))
    interact_feat <- Interaction$new(predictor, feature= feat_imp_permut_top20$feature[l], grid.size = 15)
    interact_feat_list[[l]] = interact_feat$results
  }
  interact_feat_df = bind_rows(interact_feat_list) 
  
  

  # save pdp effects for each compound list seperately
  saveRDS(interact_data_list, file = paste(path, "interaction_data_list_", compound_name, "_", target_variable, ".RDS", sep = ""))
  saveRDS(interact_feat_df, file = paste(path, "interaction_eachfeat_", compound_name, "_", target_variable, ".RDS", sep = ""))
  
  }

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


# start new loop to separate preprocessing from plotting
for(k in seq_along(compound_names)) {

  compound_name = compound_names[k] 

  compound_name_title = compound_names_key %>% filter(lawa_name == compound_name) %>% select(engl_name)
  compound_name_title  = compound_name_title$engl_name #engl_name
  print(compound_name)
  path = paste(path0, compound_name, "/", sep="")
  
 interact_data_list = readRDS(file = paste(path, "interaction_data_list_", compound_name, "_", target_variable, ".RDS", sep = "")) 
 interaction_eachfeat = readRDS(file = paste(path, "interaction_eachfeat_", compound_name, "_", target_variable, ".RDS", sep = "")) 
 
 df_final = readRDS(file = paste(path, "df_data_final_", compound_name, "_", target_variable, ".RDS", sep = ""))
  
 
 interact_data_df =  bind_rows(interact_data_list)
 
 interact_data_stat = interact_data_df  %>% group_by(.feature) %>% summarise(mean_int=mean(.interaction),
                                                                                     sd_int=sd(.interaction)) %>% ungroup()

 ### merge with other data frame defining feature categories and colors

 cols = c("meteorology" = col_set2[2], "groundwater chemistry" = col_set[3], "soil" = col_set[10], "hydrogeology" = col_set[12], "topography" = col_set[2],
          "land use" = col_set[4], "water body" = col_set[5], "site" = col_set[6], "mohp" = "grey80", "water budget" = col_set[1])

 colors_feature_categories$short_name_feature = ifelse(grepl("other agr", colors_feature_categories$short_name_feature, fixed = TRUE), "other agric. areas", colors_feature_categories$short_name_feature)
 colors_feature_categories$short_name_feature = ifelse(grepl("other non-", colors_feature_categories$short_name_feature, fixed = TRUE), "other non-agric. areas", colors_feature_categories$short_name_feature)

 # colors_feature_categories = readRDS(file = "V:/PROJEKTE/FARM/Daten/Datenanalyse/ML/ml_code/output/useful/colors_feature_categories_DE_update.RDS")
 # colors_feature_categories = colors_feature_categories %>% distinct()
 # cols = c("Meteorologie" = col_set2[2], "GW Chemie" = col_set[3], "Boden" = col_set[10], "Hydrogeologie" = col_set[12], "Topographie" = col_set[2], 
 #          "Landnutzung" = col_set[4], "Gewässer" = col_set[5], "Messstelle" = col_set[6], "mohp" = "grey80", "Wasserhaushalt" = col_set[1])
 # 
 # colors_feature_categories$short_name_feature = ifelse(grepl("Sonstige land", colors_feature_categories$short_name_feature, fixed = TRUE), "Sonstige LaWi-Fläche", colors_feature_categories$short_name_feature)
 # colors_feature_categories$short_name_feature = ifelse(grepl("Sonstige nicht", colors_feature_categories$short_name_feature, fixed = TRUE), "Sonstige nicht-LaWi-Fläche", colors_feature_categories$short_name_feature)
 # 
 # cols = setNames(paste0(cols, "80"), names(cols)) # add transparency
 
 # Apply left_join to each data frame in the list
names(colors_feature_categories) 
interact_data_stat$feature = interact_data_stat$.feature
interact_data_stat$.feature= NULL
names(interact_data_stat)
colors_feature_categories = colors_feature_categories %>% select(c(feature, short_name_feature, category))
interact_data_stat = left_join(interact_data_stat, colors_feature_categories)
interact_data_stat = interact_data_stat%>% distinct()
axis = c(min(interact_data_stat$sd_int), max(c(interact_data_stat$sd_int+interact_data_stat$mean_int)))

overall_int_plot = 
   ggplot(interact_data_stat, aes(y = reorder(short_name_feature, mean_int), x = mean_int)) +
   geom_bar(aes(fill = category), stat = "identity", color = "black") +
   geom_errorbar(aes(xmin = mean_int - sd_int,
                     xmax = mean_int + sd_int), width = 0.2, color = "black") +
   scale_fill_manual(values = cols) + #facet_grid(~.class, scales = "free", space = "free") +
   labs(fill = "feature category", y = "", x = "Overall feature interaction strength",  title = firstup(compound_name_title ) )+
   geom_text(aes(label = short_name_feature, x = - 0.12 * (axis[2] - axis[1])), hjust = 1, size = 4) +
   coord_cartesian(xlim = axis, clip = "on") +
   theme_bw() + 
   theme(plot.title = element_text(face = "bold", size = 14), 
         plot.subtitle = element_text(size = 12),
         axis.text = element_text(size = 12), 
         axis.title = element_text(size = 12),
         # plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 170, unit = "pt"),
         legend.text = element_text(size = 12),
         legend.title = element_text(size = 12),
         legend.key.size = unit(0.7, "cm"),
         axis.text.y = element_text(size = 12),  # Enable y-axis text
         strip.placement = "outside",
         strip.text = element_text(size = 12),
         axis.text.y.right = element_blank(),  # Remove right-side y-axis labels
         axis.ticks.y.right = element_blank()) 


overall_int_plot
saveRDS(overall_int_plot, paste(path_plots, compound_name, "_interact_", target_variable,"_PLOT.RDS", sep="")) 

ggsave(paste(path_plots, compound_name, "_interact_", target_variable,".png", sep=""), overall_int_plot, dpi = 300, height =8 , width = 12) 
  
}


list.plots = list.files(path_plots, pattern="_PLOT.RDS", full.names = T) #permutation_thresh


list.plots
k_10_plot =readRDS(list.plots[1])
k_26_plot =readRDS(list.plots[2])
k_27_plot =readRDS(list.plots[3])
k_28_plot =readRDS(list.plots[4])
k_30_plot =readRDS(list.plots[5])
k_31_plot =readRDS(list.plots[6])
k_38_plot =readRDS(list.plots[7])
k_40_plot =readRDS(list.plots[8])

# 10, 26, 27, 28, 30, 31, 38, 40
axis=c(0,0.6)
k_10_plot = k_10_plot +  coord_cartesian(xlim = axis, clip = "on") #+ ggtitle(compound_names_key$engl_name[1])
k_26_plot = k_26_plot + coord_cartesian(xlim = axis, clip = "on") #+ggtitle(compound_names_key$engl_name[2])
k_27_plot = k_27_plot + coord_cartesian(xlim = axis, clip = "on") #+ggtitle(compound_names_key$engl_name[3])
k_28_plot = k_28_plot + coord_cartesian(xlim = axis, clip = "on") #ggtitle(compound_names_key$engl_name[4])
k_30_plot = k_30_plot + coord_cartesian(xlim = axis, clip = "on") #ggtitle(compound_names_key$engl_name[5])
k_31_plot = k_31_plot + coord_cartesian(xlim = axis, clip = "on") #ggtitle(compound_names_key$engl_name[6])
k_38_plot = k_38_plot + coord_cartesian(xlim = axis, clip = "on") #ggtitle(compound_names_key$engl_name[7])
k_40_plot = k_40_plot + coord_cartesian(xlim = axis, clip = "on") #ggtitle(compound_names_key$engl_name[8])

r = ggpubr::ggarrange(k_10_plot, k_28_plot, k_26_plot,  k_27_plot,  k_30_plot,  k_31_plot, k_38_plot,k_40_plot,
              ncol=2, nrow=4, common.legend = T, legend="bottom",
               widths = c(rep(1,7),2), heights = c(rep(1,7),2))

r
# 
ggsave(file = paste(path_plots,"fig06.png", sep = ""), r, width = 12, height = 15)
# ggsave(file = paste("V:/PROJEKTE/FARM/Paper/Paper_ML/plots/plots_overall_feat_interaction.png", sep = ""), r, width = 25, height = 25)


# ggsave(file = paste("V:/PROJEKTE/FARM/Dokumentation/Projektbericht/Endbericht/Endbericht/plots/plots_overall_feat_interaction.png", sep = ""), r, width = 25, height = 25)

# 
# for(k in c( 10, 26, 27, 28, 30, 31, 38)){ # 3 classes: c(10, 26, 27, 28, 30, 31, 38, 40) ; # 2 classes: c(6, 10, 14, 23, 26, 27, 28, 30, 31, 36, 37, 38, 40)
#   # for(k in c( 40)){ # 3 classes: c(10, 26, 27, 28, 30, 31, 38, 40) ; # 2 classes: c(6, 10, 14, 23, 26, 27, 28, 30, 31, 36, 37, 38, 40)
#   
#   
#   
#   compound_name = filter_prio[k] 
#   compound_name_title = compound_names_key %>% filter(lawa_name == compound_name) %>% select(dt_name)
#   compound_name_title  = compound_name_title$dt_name
#   print(compound_name)
#   path = paste(path0, compound_name, "/", sep="")
#   
# 
#   interaction_eachfeat = readRDS(file = paste(path, "interaction_eachfeat_", compound_name, "_", target_variable, ".RDS", sep = "")) 
#   
#   df_final = readRDS(file = paste(path, "df_data_final_", compound_name, "_", target_variable, ".RDS", sep = ""))
# feats = names(df_final)
# feats = setdiff(feats ,c("conc_group"))
# 
# for(i in 1:length(feats)){
# one_feat = interaction_eachfeat %>% 
#   filter(.class=="1", stringr::str_detect(.feature, paste("^", feats[i], sep=""))) %>%
#   mutate(.feature = str_remove(.feature, "^[^:]+:"), ref_feat = feats[i])   %>% distinct()
# 
# names(colors_feature_categories) 
# one_feat$feature = one_feat$.feature
# one_feat$.feature= NULL
# colors_feature_categories = colors_feature_categories %>% select(c(feature, short_name_feature, category))
# one_feat = left_join(one_feat, colors_feature_categories)
# one_feat = one_feat%>% distinct()
# 
# ref_feat =unique(one_feat$ref_feat)
# ref_feat_title = colors_feature_categories %>% filter(feature == unique(one_feat$ref_feat))
# ref_feat_title  = ref_feat_title$short_name_feature
# axis = c(0,max(interaction_eachfeat$.interaction)+0.05)
# 
# 
# feat_int_plot =
#   ggplot(one_feat, aes(y = reorder(short_name_feature, .interaction), x = .interaction)) +
#   geom_bar(aes(fill = category), stat = "identity", color = "black") +
#   # geom_errorbar(aes(xmin = mean_int - sd_int,
#                     # xmax = mean_int + sd_int), width = 0.2, color = "black") +
#   scale_fill_manual(values = cols) + #facet_grid(~.class, scales = "free", space = "free") +
#   labs(fill = "feature category", y = "", x = "Specific feature interaction strength", subtitle = ref_feat_title, title =   compound_name_title ) +
#   geom_text(aes(label = short_name_feature, x = - 0.12 * (axis[2] - axis[1])), hjust = 1, size = 4) +
#   coord_cartesian(xlim = axis, clip = "on") +
#   theme_bw() + 
#   theme(plot.title = element_text(face = "bold", size = 14), 
#         plot.subtitle = element_text(size = 15),
#         axis.text = element_text(size = 15), 
#         axis.title = element_text(size = 15),
#         plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 170, unit = "pt"),
#         legend.text = element_text(size = 11),
#         legend.title = element_text(size = 15),
#         legend.key.size = unit(0.7, "cm"),
#         axis.text.y = element_text(size = 15),  # Enable y-axis text
#         strip.placement = "outside",
#         strip.text = element_text(size = 15),
#         axis.text.y.right = element_blank(),  # Remove right-side y-axis labels
#         axis.ticks.y.right = element_blank()) 
# 
# dir.create(paste(path_plots0,"feature_interactions/",compound_name,"/", sep=""))
# saveRDS(feat_int_plot, paste(path_plots0,"feature_interactions/",compound_name,"/", compound_name, "_interact_", ref_feat,"_PLOT.RDS", sep="")) 
# 
# ggsave(paste(path_plots0,"feature_interactions/", compound_name,"/", compound_name, "_interact_", ref_feat_title,".png", sep=""), feat_int_plot, dpi = 300, height =10 , width = 15) 
# 
# }
# 
#   
#   interact_data_df =  bind_rows(interact_data_list)
#   
#   interact_data_stat = interact_data_df  %>% group_by(.feature) %>% summarise(mean_int=mean(.interaction),
#                                                                                       sd_int=sd(.interaction)) %>% ungroup()
#   
#   ### merge with other data frame defining feature categories and colors
#   # colors_feature_categories = readRDS(file = "V:/PROJEKTE/FARM/Daten/Datenanalyse/ML/ml_code/output/useful/colors_feature_categories_ENG_update.RDS")
#   # unique(colors_feature_categories$category)
#   # 
#   # cols = c("meteorology" = col_set2[2], "groundwater chemistry" = col_set[3], "soil" = col_set[10], "hydrogeology" = col_set[12], "topography" = col_set[2], 
#   #          "land use" = col_set[4], "water body" = col_set[5], "site" = col_set[6], "mohp" = "grey80", "water budget" = col_set[1])
#   # 
#   colors_feature_categories = readRDS(file = "V:/PROJEKTE/FARM/Daten/Datenanalyse/ML/ml_code/output/useful/colors_feature_categories_DE_update.RDS")
#   colors_feature_categories = colors_feature_categories %>% distinct()
#   cols = c("Meteorologie" = col_set2[2], "GW Chemie" = col_set[3], "Boden" = col_set[10], "Hydrogeologie" = col_set[12], "Topographie" = col_set[2], 
#            "Landnutzung" = col_set[4], "Gewässer" = col_set[5], "Messstelle" = col_set[6], "mohp" = "grey80", "Wasserhaushalt" = col_set[1])
#   
#   colors_feature_categories$short_name_feature = ifelse(grepl("Sonstige land", colors_feature_categories$short_name_feature, fixed = TRUE), "Sonstige LaWi-Fläche", colors_feature_categories$short_name_feature)
#   colors_feature_categories$short_name_feature = ifelse(grepl("Sonstige nicht", colors_feature_categories$short_name_feature, fixed = TRUE), "Sonstige nicht-LaWi-Fläche", colors_feature_categories$short_name_feature)
#   
#   # cols = setNames(paste0(cols, "80"), names(cols)) # add transparency
#   
#   # Apply left_join to each data frame in the list
#   names(colors_feature_categories) 
#   interact_data_stat$feature = interact_data_stat$.feature
#   interact_data_stat$.feature= NULL
#   colors_feature_categories = colors_feature_categories %>% select(c(feature, short_name_feature, category))
#   interact_data_stat = left_join(interact_data_stat, colors_feature_categories)
#   interact_data_stat = interact_data_stat%>% distinct()
#   axis = c(min(interact_data_stat$sd_int), max(c(interact_data_stat$sd_int+interact_data_stat$mean_int)))
#   
#   overall_int_plot = 
#     ggplot(interact_data_stat, aes(y = reorder(short_name_feature, mean_int), x = mean_int)) +
#     geom_bar(aes(fill = category), stat = "identity", color = "black") +
#     geom_errorbar(aes(xmin = mean_int - sd_int,
#                       xmax = mean_int + sd_int), width = 0.2, color = "black") +
#     scale_fill_manual(values = cols) + #facet_grid(~.class, scales = "free", space = "free") +
#     labs(fill = "feature category", y = "", x = "Overall feature interaction strength",  title =   compound_name_title ) +
#     geom_text(aes(label = short_name_feature, x = - 0.12 * (axis[2] - axis[1])), hjust = 1, size = 4) +
#     coord_cartesian(xlim = axis, clip = "on") +
#     theme_bw() + 
#     theme(plot.title = element_text(face = "bold", size = 14), 
#           plot.subtitle = element_text(size = 15),
#           axis.text = element_text(size = 15), 
#           axis.title = element_text(size = 15),
#           plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 170, unit = "pt"),
#           legend.text = element_text(size = 11),
#           legend.title = element_text(size = 15),
#           legend.key.size = unit(0.7, "cm"),
#           axis.text.y = element_text(size = 15),  # Enable y-axis text
#           strip.placement = "outside",
#           strip.text = element_text(size = 15),
#           axis.text.y.right = element_blank(),  # Remove right-side y-axis labels
#           axis.ticks.y.right = element_blank()) 
#   
#   saveRDS(overall_int_plot, paste(path_plots0,"feature_interactions/overall/", compound_name, "_interact_", target_variable,"_PLOT.RDS", sep="")) 
#   
#   ggsave(paste(path_plots0,"feature_interactions/overall/", compound_name, "_interact_", target_variable,".png", sep=""), overall_int_plot, dpi = 300, height =8 , width = 12) 
#   
# }
# 
