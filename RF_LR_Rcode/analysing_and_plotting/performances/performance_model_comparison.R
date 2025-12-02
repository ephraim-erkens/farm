### script for plotting performance comparison across model types and metabolites

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(ggpattern)
library(ggrepel)
library(patchwork)


# set working directory (adapt to own directory!)
###
wd = ""
setwd(paste0(wd, "analysing_and_plotting/Performances/"))
###

# read compound names for which models are build individually
path_training = "../../../training/"
compound_names_key = readRDS(file = paste0(path_training, "compound_names_key.RDS"))
compound_names = compound_names_key$lawa_name


################################################################################
### extract performances of FNN models

df_comps_folds_FNN = data.frame()
for(k in seq_along(compound_names)){
  
  compound_name = compound_names[k]
  path_input_FNN = paste0(wd, "training/FNN/FNN_output_data/", compound_name, "/")
  performances_FNN = read_csv(file = paste0(path_input_FNN, "performances.csv"))
  df_comps_folds_FNN = rbind(df_comps_folds_FNN, performances_FNN)
}

  # df1 = df_comps_folds_FNN %>%
  #   filter(test_set == "outer")

  df_comps_final_FNN = df_comps_folds_FNN %>% 
    dplyr::group_by(compound_name) %>%
    dplyr::summarise(
      f1_test_sd = sd(f1_tuned, na.rm = TRUE),
      f1_test = mean(f1_tuned, na.rm = TRUE)) %>%
    dplyr::ungroup()
  df_comps_final_FNN$f1_train = NA
  df_comps_final_FNN$f1_train_sd = NA
  df_comps_final_FNN$model = "FNN"
  
  df_comps_final_FNN$model = factor(df_comps_final_FNN$model)
  
  # join english names from lookup table
  df_comps_final_FNN <- df_comps_final_FNN %>%
    left_join(compound_names_key, by = c("compound_name" = "lawa_name")) %>%
    select(-c(compound_name, dt_name))
  

################################################################################
### extract performances of RF and LR benchmark models

folder_names = c("RF_output_data", "LR_output_data") # model_sentinel_binary_logreg_paper
model_name = c("RF", "LR")
df_comps_folds_bench = data.frame()
for (f in 1:length(folder_names)){
  f_folder = folder_names[f]
  for (w in c(1:10)){ # for each seed
    for(k in seq_along(compound_names)){ 
      compound_name = compound_names[k]
      path_input = paste0(wd, "training/", model_name[f], "/", f_folder, "/", compound_name, "/")
      HP_thresh_max_perf_outer_folds = readRDS(file = paste0(path_input, "HP_thresh_max_perf_outer_folds_", compound_name, "_conc_group_", w, ".RDS"))
      perf_df = as.data.frame(bind_cols(HP_thresh_max_perf_outer_folds[1], HP_thresh_max_perf_outer_folds[4], HP_thresh_max_perf_outer_folds[7]))
      names(perf_df) = c("fold", "f1_train", "f1_test")
      perf_df$compound_name = compound_name
      perf_df$model = model_name[f]
      perf_df$seed = w
      df_comps_folds_bench = dplyr::bind_rows(df_comps_folds_bench, perf_df)
    }
  }
}
  
df_comps_folds_bench$fold = as.factor(df_comps_folds_bench$fold)
# join english names from lookup table
df_comps_folds_bench <- df_comps_folds_bench %>%
  left_join(compound_names_key, by = c("compound_name" = "lawa_name"))

# calculate means of train and test performances
df_comps_folds_bench_seed =  df_comps_folds_bench %>% 
  dplyr::group_by(engl_name, model, seed) %>%
  dplyr::summarise(
    f1_train = mean(f1_train, na.rm = TRUE),
    f1_test = mean(f1_test, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()
# calculate mean and standard deviation over all seeds
df_comps_final_bench =  df_comps_folds_bench_seed %>% 
  dplyr::group_by(engl_name, model) %>%
  dplyr::summarise(
    f1_train_sd = sd(f1_train, na.rm = TRUE),
    f1_train = mean(f1_train, na.rm = TRUE),
    f1_test_sd = sd(f1_test, na.rm = TRUE),
    f1_test = mean(f1_test, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()
df_comps_final_bench$model = factor(df_comps_final_bench$model, levels = c(model[1], model[2]))

################################################################################
### join, process and plot performances

comp_final_all = rbind(df_comps_final_bench, df_comps_final_FNN)
comp_final_long_all_mean = comp_final_all %>%
  select(-c(f1_train_sd, f1_test_sd)) %>%
  pivot_longer(cols = c("f1_train", "f1_test"), names_to = "set", values_to = "f1_mean") # convert to long format
comp_final_long_all_sd = comp_final_all %>% 
  select(-c(f1_train, f1_test)) 
colnames(comp_final_long_all_sd) <- c("engl_name", "model", "f1_train", "f1_test")
comp_final_long_all_sd = comp_final_long_all_sd %>%
  pivot_longer(cols = c("f1_train", "f1_test"), names_to = "set", values_to = "f1_sd")
comp_final_long_all = left_join(comp_final_long_all_mean, comp_final_long_all_sd, by = c("engl_name", "model", "set"))

comp_final_long_all$model = factor(comp_final_long_all$model, levels = c("RF", "FNN", "LR"))


comp_final_long_test <- comp_final_long_all %>%
  filter(set == "f1_test")
comp_final_long_test$engl_name <- factor(comp_final_long_test$engl_name, levels = c("desphenyl chloridazone", "methyl desphenyl chloridazone", "metazachlor esa" ,
                                                                                                    "metazachlor oxa", "metolachlor esa", "metolachlor oxa", "metolachlor noa 413173", "trifluoroethanoic acid"))

model_comp =
ggplot(comp_final_long_test, aes(x = model, y = f1_mean, shape = set), size = 2) +
  geom_point(aes(color = set), shape = 15, size = 3, col = "orange") +
  scale_y_continuous(limits = c(0.6, 0.8)) +
  geom_errorbar(data = comp_final_long_test, aes(ymin = f1_mean - f1_sd,
                                                ymax = f1_mean + f1_sd), width = 0.08, color = "black") +
  facet_wrap(~ engl_name, nrow = 2) +
  labs(x = "model type", y = "macro-F1 score", color = NULL) +
  theme_bw() +
  geom_text_repel(data = comp_final_long_test, aes(label = round(f1_mean, 3)), hjust = -0.5, vjust = 0.3, size = 3, segment.color = NA) + # Add text labels next to points
  theme(text = element_text(size = 12), legend.text = element_text(size = 12, hjust = 0), plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 10))
ggsave(file = paste0("model_comp.png"), model_comp, width = 10, height = 6, dpi = 600)



