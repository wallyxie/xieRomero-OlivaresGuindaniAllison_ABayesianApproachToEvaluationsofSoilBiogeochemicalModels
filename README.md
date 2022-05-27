# Scripts and code corresponding to "A Bayesian Approach to Evaluation of Soil Biogeochemical Models
## by Hua Xie, Adriana Romero-Olivares, Michele Guindani, and Steven D. Allison

Included in this repository is Stan and R code for evaluating goodness of fit (by leave-one-out cross-validation (LOO-CV) and widely applicable information criterion (WAIC)) of the CON and AWB soil biogeochemical ODE models from Allison et al. 2010 to soil microbe respiration data. If the R code, Stan code, and soil respiration data set ("adriana_data_s_min.csv") are all downloaded to the same directory, the code should run on a cluster after being in a job script, given that the required R packages have been installed and can be loaded.
