
# Concentrations of Plant Protection Product Metabolites in Groundwater in Germany: Assessing the Influence of Site Characteristics using Random Forests and Feedforward Neural Networks


## Abstract

The requirements for a good chemical status of groundwater defined by the European Water Framework Directive are mainly breached  due to nitrate, agricultural pesticides and their metabolites. 
In earlier work, an extensive database of publicly available groundwater monitoring data on active substances and metabolites of plant protection products in Germany has been gathered and harmonised. In this study, we analysed data for eight frequently detected metabolites and trained Random Forest and Feedforward Neural Network models on temporally aggregated data, with Logistic Regression models serving as a baseline. The models classified concentration classes with macro-F1 scores ranging from 0.626 to 0.762. Random Forest models consistently outperformed the Logistic Regression benchmark models, while the Feedforward Neural Network models achieved only marginally higher performances than the benchmarks on average. The most important explanatory variables were the main crops for which the respective parent compounds of the metabolites are authorised. Deeper well screens were associated with lower metabolite concentrations. Aquifer property features showed the lowest feature importance. Soil sand content and long-term groundwater recharge exhibited comparatively high feature importance, although their influence varied across metabolite models and they may serve as spatial proxies for agricultural areas. This study highlights the need for spatially and temporally explicit data on actual pesticide application volumes to foster the understanding of groundwater vulnerability. Despite the challenges associated with strongly and multiply censored compound concentration data, this study generated promising results while critically discussing the implications for interpretation.

## Set up

- Add environment details. 
- Add requirements file. 
- separate folders for 3 models
- Add a data folder (?)

## Code details

- `neural_networks.py` contains functions to build, train and evaluate a Feedforward Neural Network
- `FNN_nestedCV.ipynb` was used to select compounds and features, to format the model input data and to perform a nested cross-validation using the functions provided in `neural_networks.py`
