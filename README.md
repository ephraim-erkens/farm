
This repository accompanies the manuscript “Concentrations of Plant Protection Product Metabolites in Groundwater in Germany: Assessing the Influence of Site Characteristics using Random Forests and Feedforward Neural Network” by Cooke, A.-K., Mohanty, S., Joger, F., Erkens, E. and Broda, S.

In this repository, training code of the Random Forest (RF), Logistic Regression (LR) and Feedforward Neural Network (FNN) models are provided. Input data and the final trained models of  RF, LR and FNN models are published on Zenodo ([https://doi.org/10.5281/zenodo.17710126](https://doi.org/10.5281/zenodo.17710126)). On Zenodo, a readme-file provides more information on the deposited model files.


## Abstract

The requirements for a good chemical status of groundwater defined by the European Water Framework Directive are mainly breached  due to nitrate, agricultural pesticides and their metabolites. Understanding the drivers that control groundwater contamination is vital to protect groundwater resources. In this study, we used Random Forest, Feedforward Neural Network and Logistic Regression models to identify the most important features for the prediction of plant protection product metabolite concentrations in groundwater. We trained our models individually for eight frequently detected metabolites based on an extensive database of publicly available groundwater monitoring data comprising active substances and metabolites of plant protection products in Germany.  The most important explanatory variables were the main crops for which the respective parent compounds of the metabolites are authorised. Deeper well screens were associated with lower metabolite concentrations. Soil sand content and long-term groundwater recharge exhibited comparatively high feature importance, although their influence varied across metabolite models and they may serve as spatial proxies for agricultural areas. Aquifer property features showed the lowest feature importance.
The models classified concentration classes with macro-averaged F1-scores ranging from 0.626 to 0.762. Random Forest models consistently performed better than the Logistic Regression benchmark models, while the Feedforward Neural Network models achieved only marginally higher performances than the benchmarks on average. This study demonstrates the benefit of machine learning models to foster the understanding of groundwater vulnerability, despite the challenges associated with strongly and multiply censored compound concentration data. Our results indicate a pronounced link between pesticide application volumes and concentrations of their metabolites in groundwater. This highlights the need for spatially and temporally explicit data on actual pesticide application volumes to fully capture groundwater vulnerability.


## Set up

- Add environment details. 
- Add requirements file. 
- separate folders for 3 models
- Add a data folder (?)

## Code details

### Feedforward Neural Network

FNN models were written in Python. The files are located in the folder `FNN_Python_code`:

- `neural_networks.py` contains functions to build, train and evaluate a Feedforward Neural Network
- `FNN_nestedCV.ipynb` was used to select compounds and features, to format the model input data and to perform a nested cross-validation using the functions provided in `neural_networks.py`


### Random Forest and Logistic Regression

RF and LR models were written in R. The files are located in the folder `RF_LR_Rcode`.

#### 1. R-Training code for Radom Forest and Logistic Regression:

The training code for Random Forest and Logistic Regression models is provided in the following folders and files:

- `..\training\LR\LR_model_training.R`
- `..\training\RF\RF_model_training.R`

Training further requires the file `..\compound_names_key.RDS` and the functions defined in the .R files provided in the subfolder …`\training\custom_performance_measures\`.

#### 2. Input-Data

`..\training\input_data\`

This folder contains the training input-data for each of the 8 metabolites used in the training scripts mentioned above. Data is accessed as .RDS files in the code, for better viewing .csv files have been provided, additionally.


The tables have the following columns. Each row represents one groundwater measuring site.

- conc_group: concentration class of metabolites (1,2)
- filter_ok_unter_gok: well screen depth in m below surface
-   HA: cavity type of aquifer (P : porous, K/Ka: fractured/karstic, K: fractured, P/K: porous/fractured)
-   kf_bez: aquifer hydraulic conductivity class (1,2,3)
-   gwn_mean: mean annual groundwater recharge rate in mm
-   sand_depth_weighted_mean: soil sand content (averaged over 1 m soil depth) in %
-   corg_gehalt_depth_weighted_mean: soil organic carbon content (averaged over 1 m soil depth) in %
-   Makroporen_rounded: macropore classes (1,2,3)
-   area_sentinel_id_XX : area in m² of crop type, “XX” as specified in the following table:

| Crop Class                                               | Kultur Sentinel ID   |
| -------------------------------------------------------- | -------------------- |
| Winter wheat                                             | 1                    |
| Winter barley                                            | 2                    |
| Winter rye                                               | 3                    |
| Other winter cereals                                     | 4                    |
| Spring barley                                            | 5                    |
| Spring oat                                               | 6                    |
| Other spring cereals                                     | 7                    |
| Grassland                                                | 8                    |
| Legumes                                                  | 9                    |
| Grain maize and silage maize                             | 10                   |
| Rape seed                                                | 12                   |
| Sun flower                                               | 13                   |
| Sugar beet                                               | 14                   |
| Potato                                                   | 15                   |
| Vegetable                                                | 16                   |
| Wineyard                                                 | 17                   |
| Fruit                                                    | 18                   |
| Hops                                                     | 19                   |
| Small woody features                                     | 20                   |
| Other agricultural land with potential pesticide use     | 21                   |
| Fodder crops                                             | 22                   |
| Industrial crops                                         | 23                   |
| Spring wheat                                             | 24                   |
| Other agricultural land without pesticide use            | 25                   |
| Other non-agricultural land                              | 26                   |
| Other non-agricultural land with potential pesticide use | 27                   |

  
  

Further details on the source and pre-processing of the features are given in the manuscript.

#### 3. Random Forest and Logistic Regression trained models

The following folders have to downloaded from Zenodo (https://doi.org/10.5281/zenodo.17710126). The folders only contain a link to the repository. There these folders contain the trained models, and several files required for further calculations, display of the feature importances etc. The plotting scripts need to access some of these files, some are generated by the plotting scripts, but have been provided here for completeness.

`.\LR\LR_output_data\`

`.\RF\RF_output_data\`

They are further divided into subfolders for each of the 8 metabolites.

As an example, the folder of the Random Forest results for the compound desphenyl chloridazone can be found in:

`..\training\RF\RF_output_data\desphenyl-chloridazon`.

Each of these subfolders contain the following file types (.RDS) 10 times each, respectively (for each random seed, see training code for specifications):

**Trained model:**
- `full_model_desphenyl-chloridazon_conc_group_1.RDS`

**F1-Score and decision threshold-based hyper parameter tuning results:**
- `HP_f1_thresh_max_inner_folds_desphenyl-chloridazon_conc_group_1.RDS`
- `HP_thresh_max_perf_outer_folds_desphenyl-chloridazon_conc_group_1.RDS`
- `HP_thresh_max_perf_final_desphenyl-chloridazon_conc_group_1.RDS`
- `best_HP_f1_thresh_max_single_cv_desphenyl-chloridazon_conc_group_1.RDS`

**Tasks for test and train, respectively, of outer folds of cross validation:**
- `outer_train_task_list_desphenyl-chloridazon_conc_group_1.RDS`
- `outer_test_task_list_desphenyl-chloridazon_conc_group_1.RDS`

**Prediction on outer folds training and test folds of cross-validation:**
- `pred_outer_train_list_desphenyl-chloridazon_conc_group_1.RDS`
- `pred_outer_test_list_desphenyl-chloridazon_conc_group_1.RDS`

**Further training files:**
- `rr_outer_list_desphenyl-chloridazon_conc_group_1.RDS`
- `rr_outer_gini_cor_list_desphenyl-chloridazon_conc_group_1.RDS`

**mlr3-Classifiation Task:**
- `task_desphenyl-chloridazon_conc_group.RDS`

**Data-frame with inputfeatures**
- `df_data_final_desphenyl-chloridazon_conc_group.RDS`

**Feature Importances:**

**Gini impurity feature importance**
- `feat_imp_ginithresh_desphenyl-chloridazon_conc_group_1.RDS`

**bias-correct Gini Impurity feature importance**
- `feat_imp_gini_corthresh_desphenyl-chloridazon_conc_group_1.RDS`

**Permutation feature importance**
- `feat_imp_permutthresh_desphenyl-chloridazon_conc_group_1.RDS`

**Permutation feature importance per seeds**
- `feat_imp_permutthreshseeds_desphenyl-chloridazon_conc_group_1.RDS`

**Feature Effects (partial dependences)**
- `pdp_data_list_desphenyl-chloridazon_conc_group.RDS`

**Feature Interactions:**
- `interaction_data_list_desphenyl-chloridazon_conc_group.RDS`
- `interaction_eachfeat_desphenyl-chloridazon_conc_group.RDS`

#### 4. Plotting Scripts

Plotting requires the files `..\compound_crops.RDS` and `..\colors_feature_categories_ENG.RDS`, and access to the model-out folders (e.g. `..\RF\RF_output_data\` )

The Scripts are provided in the folder `..\analysing_and_plotting`, divided into subfolders:

- Feature_Effects (script for calculation and plotting of partial dependence plots)
- Feature_Importance (scripts for calculation and plotting of feature importances)
- Feature_Interactions (script for calculation and plotting of feature interactions)
- Performances (script for plotting of performance barplots per model and cross-model comparison plot)
- Residual_analysis (script for plotting of residuals against feature values)
- R-environment

Specifics of the used R-version and library information are provided in the folder:
`..\renv\`
as well as by the files
`..\.Rhistory`, `.Rprofile`, and `renv.lock`