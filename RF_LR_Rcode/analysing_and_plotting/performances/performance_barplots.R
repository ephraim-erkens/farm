#### Skript for plotting Randomf-Forest performances (F1-score, precision, recall) and class distribution as barplots per compound
#### Paper: Figure 3
library(tidyverse)
library(dplyr)
library(sf)
library(ggplot2)
library(svglite)
library(rcartocolor)
library(ggpubr)
library(ggExtra)
library(mlr3measures)
library(mlr3)
library(patchwork)
# library(mlr3viz)
library(RColorBrewer)
library(precrec)
library(yardstick)
library(ggpattern)
library(pROC)
library(reshape2)

# set working directory (adapt to own directory!)
### 
setwd("V:/PROJEKTE/FARM/Dokumentation/Datenübergabe_UBA/ML_Code/code_paper_submission/training/RF/")
###
                  
######################################################
### specify folder names

folder_name = "RF_output_data" # choose folder name of models

######################################################

#### load custom performance measures to be evaluated ----
# custom measure class for F1-score in multiclass classification
MeasureF1Multiclass = readRDS(file = "../custom_performance_measures/MeasureF1Multiclass.RDS")
# Register the measure
f1_multiclass_measure <- MeasureF1Multiclass$new() # Instantiate the custom measure
mlr3::mlr_measures$add("f1_multiclass_measure", MeasureF1Multiclass) # add the new measure to the dictionary
# load function to calculate precision, recall and f1 score per class for both mlr3 prediction objects and prepared dataframes
perf_rec_prec_f1 = readRDS(file = "../custom_performance_measures/perf_rec_prec_f1.RDS")


# read compound names for which models are build individually
compound_names_key = readRDS(file = "../compound_names_key.RDS")
compound_names = compound_names_key$lawa_name

target_variable = "conc_group"


path0 = paste(folder_name,"/",  sep="")
path_plots= ("plots/")
dir.create(path_plots)
path_plots = paste(path_plots, "performances/", sep="") #
dir.create(path_plots)

######################################################

### start of main loop to choose compound
  for(k in seq_along(compound_names)) {
    
    
compound_name = compound_names[k] 
cat("\n", compound_name, "\n", k, "\n") # to keep track

compound_name_title = compound_names_key %>% filter(lawa_name == compound_name) %>% select(engl_name)
compound_name_title  = compound_name_title$engl_name


path = paste(path0, "/", compound_name, "/", sep="") # todo



cm_df_long1 = data.frame()
HP_f1_thresh_max_outer_folds1_long = data.frame()

### loop over different seeds
# seeds = c() # 123, 234, 345, 456, 567, 678, 789, 987, 876, 765 #seed = seeds[w]
for (w in c(1:10)){
# w = 1
cat("\n", paste("model seed =", w), "\n") # for orientation
  
# read list of predictions on test data of outer folds
pred_outer_test_list1 = readRDS(file=paste(path, "pred_outer_test_list_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))

# read performances, optimized HP combinations and thresholds for each fold, respectively
HP_f1_thresh_max_outer_folds1 = readRDS(file=paste(path, "HP_f1_thresh_max_outer_folds_", compound_name, "_", target_variable, "_", w, ".RDS", sep = ""))
HP_f1_thresh_max_outer_folds1$seed = w
HP_f1_thresh_max_outer_folds1_long = dplyr::bind_rows(HP_f1_thresh_max_outer_folds1_long, HP_f1_thresh_max_outer_folds1)

# read predictions of all outer folds and write them into long format
pred_outer_folds1 = data.frame()
for (m in 1:length(pred_outer_test_list1)){

  pred_outer_test = pred_outer_test_list1[[m]]
  df1 = as.data.frame(bind_cols(pred_outer_test$truth, pred_outer_test$response, pred_outer_test$row_ids, pred_outer_test$prob))
  names(df1) = c("truth", "response", "row_ids", paste0("prob", 1:(length(df1)-3)))
  df1$outer_fold = m
  pred_outer_folds1 = dplyr::bind_rows(pred_outer_folds1, df1)
}
cm1 = table(Truth = pred_outer_folds1$truth, Predicted = pred_outer_folds1$response)
cm_df1 = melt(cm1)
cm_df1$seed = w
cm_df_long1 = dplyr::bind_rows(cm_df_long1, cm_df1)


} # loop seeds



cm_df_summary1 =
  cm_df_long1 %>%
  dplyr::group_by(Truth, Predicted) %>%
  dplyr::summarise(
    mean_count = mean(value),
    sd_count = sd(value)) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(Truth, desc(Predicted)) %>% # sort for calculating cumulative sum for recall
  dplyr::group_by(Truth) %>% 
  dplyr::mutate(mean_count_cum_rec = cumsum(mean_count)) %>% # cumulatively add values for predicted classes to plot error bar
  dplyr::ungroup() %>%
  dplyr::arrange(Predicted, desc(Truth)) %>% # sort for calculating cumulative sum for precision
  dplyr::group_by(Predicted) %>% 
  dplyr::mutate(mean_count_cum_prec = cumsum(mean_count)) # cumulatively add values for predicted classes to plot error bar
cm_df_summary1$Truth = as.factor(cm_df_summary1$Truth)
cm_df_summary1$Predicted = as.factor(cm_df_summary1$Predicted)


# summarize performances by averaging measures over different seeds and outer folds for each class
HP_thresh_max_perf_classes1 =
  HP_f1_thresh_max_outer_folds1_long %>%
  dplyr::group_by(class) %>%
  dplyr::summarise(
    macro_rec_test_class = mean(recall_test, na.rm = TRUE),
    macro_prec_test_class = mean(precision_test, na.rm = TRUE),
    macro_f1_test_class = mean(f1_score_test, na.rm = TRUE))
HP_thresh_max_perf_classes1$class = as.numeric(HP_thresh_max_perf_classes1$class)

### prepare plots

### colorblind friendly colormaps
# display_carto_all(colorblind_friendly = TRUE)
colors_safe = carto_pal(n = 12, name = "Safe")
colors_teal = carto_pal(n = 12, name = "Teal")

if(compound_name != "trifluoressigsaeure"){
  col_set1 = colors_safe[c(2, 9)]
  class_ranges1 = c("1" = expression("1 (x" <= "0.05 " * mu * "g/L)"), "2" = expression("2 (0.05 " * mu * "g/L < x)"))
  } else {
  col_set1 = colors_safe[c(1, 5)] #TFA
  class_ranges1 = c("1" = expression("1 (x" <= "0.05 " * mu * "g/L)"), "2" = expression("2 (0.05 " * mu * "g/L < x)"))
  }
    
x_label = 1.45
# define custom settings applied to all plots
custom_settings1 =    
  list(
    #scale_fill_manual(values=c("#69b3a2", "#404080", "red", "blue", "orange")) +
    scale_fill_manual(values = col_set1, labels = class_ranges1), # col_set3[c(2, 4, 6, 8, 10, 1, 3)]
    theme_bw(),
    theme(text = element_text(size = 14), legend.text = element_text(hjust = 0), plot.subtitle = element_text(size = 14))
  )


firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
### visualisation of recall (height of columns = ground truth)
plot_pred_rec_thresh1 =
  ggplot(cm_df_summary1, aes(x = Truth, y = mean_count)) +
    geom_bar(aes(fill = Predicted), alpha = 0.8, stat = "identity") +
    geom_errorbar(aes(ymin = mean_count_cum_rec - sd_count,
                      ymax = mean_count_cum_rec + sd_count), width = 0.2, color = "black") +
    custom_settings1 +
    scale_x_discrete(name = "concentration class") +
    scale_y_continuous(name = 'number of sites', limits = c(0,5000)) +
    labs(fill = "predicted concentration class", subtitle = paste("macro f1-score =", round(mean(HP_thresh_max_perf_classes1$macro_f1_test_class), 3))) + #, "\nDurchschnitt ROC AUC =", "\nDurchschnitt PRC AUC ="
    ggtitle(firstup(compound_name_title)) + # paste("Accuracy = ", round(df_scores_compounds[k, 2], 3), "\nBalanced accuracy = ", round(balanced_accuracy, 3), "\nMacro-averaged F2 = ", round(macro_f_beta_2, 3), "\nMacro-averaged AUC = ", round(macro_auc, 3))
    geom_text(data = HP_thresh_max_perf_classes1, aes(y = class, x = class, label = round(macro_rec_test_class, 3)), vjust = 1.8, colour = 'royalblue', size = 4) +
    geom_text(data = HP_thresh_max_perf_classes1, aes(y = class, x = class, label = round(macro_prec_test_class, 3)), vjust = - 49, colour = 'green4', size = 4) + #- 53.2
    annotate("text" , x = x_label, y = 0, label = "recall", color = "royalblue", size = 4.3, hjust = 4.55, vjust = 1.6) + 
    annotate("text" , x = x_label, y = 0, label = "precision", color = "green4", size = 4.3, hjust = 3.21, vjust = - 46.5) +
    coord_cartesian(clip = "off") +  # Allow annotations outside plot boundaries
    theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 14),
          axis.text = element_text(size = 14), axis.title = element_text(size = 14),
          plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 15, unit = "pt"))

q = plot_pred_rec_thresh1 + plot_layout(guides = "collect") 

assign(paste("k_",k,"_plot",sep=""), q)

ggsave(file = paste(path_plots,  compound_name, "_classification_plots_stacked_folds.png", sep = ""), q, width = 10, height = 7.95)


} # loop compounds


#### combine all plots #####

k_1_plot_new = k_1_plot + theme(legend.position = "none") 
k_2_plot_new = k_2_plot + theme(legend.position = "none") 
k_3_plot_new = k_3_plot + theme(legend.position = "none") 
k_4_plot_new = k_4_plot + theme(legend.position = "none")
k_5_plot_new = k_5_plot + theme(legend.position = "none")
k_6_plot_new = k_6_plot + theme(legend.position = "none")
k_7_plot_new = k_7_plot + theme(legend.position = "bottom")
k_8_plot_new = k_8_plot +  theme(legend.position = "bottom")

### all 8 metabolites plotted together
r = ggarrange(k_1_plot_new, k_2_plot_new, k_3_plot_new,  k_4_plot_new,  k_5_plot_new,  k_6_plot_new, k_7_plot_new,k_8_plot_new,
              ncol=4, nrow=2, common.legend = T, legend="bottom")

r 
ggsave(file = paste(path_plots,  "plots_performance_all_compounds.png", sep = ""), r, width = 15, height = 15)

### Only TFA
r2 = ggarrange(k_8_plot_new,
              ncol=4, nrow=2, common.legend = T, legend="bottom")

ggsave(file = paste(path_plots,  "plots_performance_tfa.png", sep = ""), r2, width = 15, height = 15)

