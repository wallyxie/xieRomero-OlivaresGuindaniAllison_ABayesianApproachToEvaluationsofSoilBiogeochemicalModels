library(ggplot2,quietly=True)
library(rstan)
library(loo)
library(bayesplot)
library(coda)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

filename = "AWB_adriana_pools5b"

#############
##FUNCTIONS##
#############

#Save summary, sample parameters, and posteriors
fit_summary = function(stan_fit, stan_fit_ex, S_0, D_0, M_0, E_0, filename) {
    #Posteriors
    write.csv(stan_fit_ex, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "posteriors", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".csv", sep = "_"))
    #Summary
    stan_fit_sum <- summary(stan_fit)$summary
    stan_fit_sum_write <- stan_fit_sum[c("V_ref", "V_U_ref", "Ea_V", "Ea_VU", "Ea_K", "Ea_KU", "E_C_ref", "m_t", "a_MS", "sigma"),
                                 c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
    write.csv(stan_fit_sum_write, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "summary", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".csv", sep = "_"))
    #Sampler parameters
    sampler_params <- get_sampler_params(stan_fit, inc_warmup = FALSE) #Access sampler values
    write.csv(sampler_params, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "sampler_params", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".csv", sep = "_"))
}

#IC calculations
calc_ic = function(stan_fit, stan_fit_ex, S_0, D_0, M_0, E_0, filename) {
    fit_log_lik = extract_log_lik(stan_fit)
    rel_n_eff <- relative_eff(exp(fit_log_lik))
    stan_fit_waic = waic(fit_log_lik, r_eff = rel_n_eff, cores = 2)
    stan_fit_loo = loo(fit_log_lik, r_eff = rel_n_eff, cores = 2)
    LPML = sum(log(1 / colMeans(stan_fit_ex$CPOi)))
    waic = stan_fit_waic$waic
    p_waic = stan_fit_waic$estimates["p_waic",]
    loo = stan_fit_loo$loo
    p_loo = stan_fit_loo$estimates["p_loo",]
    sink(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "ic", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".txt", sep = "_"))    
    cat("WAIC = ", waic, "\nLOO = ", loo, "\nLPML = ", LPML, "\np_waic = ", p_waic, "\np_loo = ", p_loo)
    sink()    
    return(c(waic, loo, LPML, p_waic, p_loo))
}

#Convert Stan object into coda object
stan2coda <- function(fit) {
    mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

#Diagnostics including posterior, traceplot, autocorrelation, and Rhat. Specific to AWB because of parameters. Cannot copy-pate to other models
bayes_diagnostics = function(stan_fit, S_0, D_0, M_0, E_0, filename) {
    #Posterior credible areas
    stan_fit.array <- as.array(stan_fit)
    CI_plot1 <- mcmc_areas(stan_fit.array, pars = c("Ea_V", "Ea_VU", "Ea_K", "Ea_KU"), prob = 0.8, prob_outer = 0.95) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "CI_Ea", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = CI_plot1)
    CI_plot2 <- mcmc_areas(stan_fit.array, pars = c("V_ref", "E_C_ref", "a_MS"), prob = 0.8, prob_outer = 0.95) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "CI_V_ref", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = CI_plot2)
    CI_plot3 <- mcmc_areas(stan_fit.array, pars = c("V_U_ref"), prob = 0.8, prob_outer = 0.95) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "CI_V_U_ref", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = CI_plot3)
    CI_plot4 <- mcmc_areas(stan_fit.array, pars = c("m_t"), prob = 0.8, prob_outer = 0.95) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "CI_m_t", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = CI_plot4)
    #Traceplot
    traceplot1 <- mcmc_trace(stan_fit.array, pars = c("Ea_V", "Ea_VU", "Ea_K", "Ea_KU")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "trace_Ea", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = traceplot1)
    traceplot2 <- mcmc_trace(stan_fit.array, pars = c("V_ref", "E_C_ref", "a_MS")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "trace_V_ref", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = traceplot2)
    traceplot3 <- mcmc_trace(stan_fit.array, pars = c("V_U_ref")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "trace_V_U_ref", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = traceplot3)
    traceplot4 <- mcmc_areas(stan_fit.array, pars = c("m_t")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "trace_m_t", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = traceplot4)
    #Autocorrelation
    acf_plot1 <- mcmc_acf(stan_fit.array, pars = c("Ea_V", "Ea_VU", "Ea_K", "Ea_KU")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "acf_Ea", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = acf_plot1)
    acf_plot2 <- mcmc_acf(stan_fit.array, pars = c("V_ref", "E_C_ref", "a_MS")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "acf_V_ref", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = acf_plot2)
    acf_plot3 <- mcmc_acf(stan_fit.array, pars = c("V_U_ref")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "acf_V_U_ref", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = acf_plot3)
    acf_plot4 <- mcmc_acf(stan_fit.array, pars = c("m_t")) + yaxis_text()
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "acf_m_t", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = acf_plot4)
    #Rhat plot
    stan_fit_coda <- stan2coda(stan_fit)
    pdf(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "Rhat", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"))
    gelman.plot(stan_fit_coda[,1:9])
    dev.off()
}

#Pearson's correlation calculation
calc_r <- function(model_vector, data_vector, S_0, D_0, M_0, E_0, filename) {
    sst_mean <- mean(data_vector)
    sst_vector <- (data_vector - sst_mean) ^ 2
    sst <- sum(sst_vector)
    ssr_vector <- (data_vector - model_vector) ^ 2
    ssr = sum(ssr_vector)
    r_sq = 1 - (ssr / sst)
    cat(r_sq, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "r_sq", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".txt", sep = "_"))
    return(r_sq)
}

#Plot posterior fit and forward simulation prediction
plot_fits <- function(stan_fit_ex, N_t, N_p, obs_times, pred_times, data_vector, S_0, D_0, M_0, E_0, filename) {
    sigma_hat <- median(stan_fit_ex$sigma)
    CO2_flux_ratios_hat_median <- rep(NA, N_t) #Pre-allocate vector for median model fit
    CO2_flux_ratios_pred_median <- rep(NA, N_p)
    sigma_hat <- median(stan_fit_ex$sigma)    
    for (t in 1:N_t) {
    CO2_flux_ratios_hat_median[t] <- median(stan_fit_ex$CO2_flux_ratios_hat_vector[,t])
    }
    for (t in 1:N_p) {
    CO2_flux_ratios_pred_median[t] <- median(stan_fit_ex$CO2_flux_ratios_new_vector[,t])
    }
    #Pearson's correlation calculation
    r_sq = calc_r(CO2_flux_ratios_hat_median, data_vector, S_0, D_0, M_0, E_0, filename)
    #Plotting posterior fit
    df_post <- data.frame(list(obs_times = obs_times / (24 * 365), CO2_flux_ratios_vector = data_vector, CO2_flux_ratios_hat_median = CO2_flux_ratios_hat_median))
    post_plot <- ggplot(df_post, aes(x = obs_times)) +
    geom_ribbon(aes(ymin = CO2_flux_ratios_hat_median - 1.96 * sigma_hat, ymax = CO2_flux_ratios_hat_median + 1.96 * sigma_hat), fill = "bisque") +
    geom_line(aes(y = CO2_flux_ratios_hat_median), colour = "black", size = 2) +
    geom_point(aes(y = CO2_flux_ratios_vector), shape = 21, colour = "black", fill = "white", size = 3) +
    geom_errorbar(interval_95pct, width = 0, color = "blue") +    
    labs(x = "Duration (days)",
         y = "Response Ratio") +
    ggtitle("AWB Response Ratio vs. Days")
    post_plot <- post_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "post", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = post_plot)    
    #Plotting predictive fit
    df_pre <- data.frame(list(pred_times = pred_times / (24 * 365), CO2_flux_ratios_pred_median = CO2_flux_ratios_pred_median))
    pre_plot <- ggplot(df_pre, aes(x = pred_times)) +
    geom_ribbon(aes(ymin = CO2_flux_ratios_pred_median - 1.96 * sigma_hat, ymax = CO2_flux_ratios_pred_median + 1.96 * sigma_hat), fill = "bisque") +
    geom_line(aes(y = CO2_flux_ratios_pred_median), colour = "black", size = 2) +
    labs(x = "Duration (years)",
         y = "Response Ratio") +
    ggtitle("AWB Predictive Response Ratio vs. Year")
    pre_plot <- pre_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%m"), filename, "pred", "S", S_0, "D", D_0, "M", M_0, "E", E_0, ".pdf", sep = "_"), plot = pre_plot)   
}

########
##DATA##
########

df = read.csv("adriana_data_s_min.csv", header=TRUE) #Read in data
hour_index_list <- df$Hours
CO2_flux_ratios_vector <- df$Response_Ratio
sd <- df$Standard_Deviation
interval_95pct <- aes(ymin = CO2_flux_ratios_vector - 1.96 * sd,
                      ymax = CO2_flux_ratios_vector + 1.96 * sd)

N_t <- length(hour_index_list)
t0 <- 0
ts_p <- seq(5000, 500000, by = 5000)
N_p <- length(ts_p)
mean_SOC_input <- 0.0009
mean_DOC_input <- 0.0001
T_ref <- 283.15
R_g <- 0.008314

S_01 <- 50; S_02 <- 75; S_03 <- 100; S_04 <- 125; S_05 <- 150; S_06 <- 175; S_07 <- 200;
D_0 <- 0.2
M_0 <- 2
E_0 <- 0.1
r_L <- 0.0005
r_E <- r_L * E_0 / M_0

adriana_datTest <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = 90, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat1 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_01, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat2 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_02, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat3 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_03, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat4 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_04, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat5 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_05, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat6 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_06, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

adriana_dat7 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_07, D_0 = D_0, M_0 = M_0, E_0 = E_0, r_L = r_L, r_E = r_E)

file_path <- "AWB_adriana_pools5.stan" #Read in Stan model code. Stan file must be in same directory.
lines <- readLines(file_path, encoding = "ASCII")
for (n in 1:length(lines)) cat(lines[n],'\n')

#############
##EXECUTION##
#############

AWB_fitTest <- stan("AWB_adriana_pools5.stan", data = adriana_datTest, iter = 500, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fitTest_ex = extract(AWB_fitTest)
AWB_fitTest_ic = calc_ic(AWB_fitTest, AWB_fitTest_ex, 90, D_0, M_0, E_0, filename)
fit_summary(AWB_fitTest, AWB_fitTest_ex, 90, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fitTest, 90, D_0, M_0, E_0, filename)
plot_fits(AWB_fitTest_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, 90, D_0, M_0, E_0, filename)

AWB_fit1 <- stan("AWB_adriana_pools5.stan", data = adriana_dat1, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit1_ex = extract(AWB_fit1)
AWB_fit1_ic = calc_ic(AWB_fit1, AWB_fit1_ex, S_01, D_0, M_0, E_0, filename)
fit_summary(AWB_fit1, AWB_fit1_ex, S_01, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit1, S_01, D_0, M_0, E_0, filename)
plot_fits(AWB_fit1_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_01, D_0, M_0, E_0, filename)

AWB_fit2 <- stan("AWB_adriana_pools5.stan", data = adriana_dat2, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit2_ex = extract(AWB_fit2)
AWB_fit2_ic = calc_ic(AWB_fit2, AWB_fit2_ex, S_02, D_0, M_0, E_0, filename)
fit_summary(AWB_fit2, AWB_fit2_ex, S_02, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit2, S_02, D_0, M_0, E_0, filename)
plot_fits(AWB_fit2_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_02, D_0, M_0, E_0, filename)

AWB_fit3 <- stan("AWB_adriana_pools5.stan", data = adriana_dat3, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit3_ex = extract(AWB_fit3)
AWB_fit3_ic = calc_ic(AWB_fit3, AWB_fit3_ex, S_03, D_0, M_0, E_0, filename)
fit_summary(AWB_fit3, AWB_fit3_ex, S_03, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit3, S_03, D_0, M_0, E_0, filename)
plot_fits(AWB_fit3_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_03, D_0, M_0, E_0, filename)

AWB_fit4 <- stan("AWB_adriana_pools5.stan", data = adriana_dat4, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit4_ex = extract(AWB_fit4)
AWB_fit4_ic = calc_ic(AWB_fit4, AWB_fit4_ex, S_04, D_0, M_0, E_0, filename)
fit_summary(AWB_fit4, AWB_fit4_ex, S_04, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit4, S_04, D_0, M_0, E_0, filename)
plot_fits(AWB_fit4_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_04, D_0, M_0, E_0, filename)

AWB_fit5 <- stan("AWB_adriana_pools5.stan", data = adriana_dat5, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit5_ex = extract(AWB_fit5)
AWB_fit5_ic = calc_ic(AWB_fit5, AWB_fit5_ex, S_05, D_0, M_0, E_0, filename)
fit_summary(AWB_fit5, AWB_fit5_ex, S_05, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit5, S_05, D_0, M_0, E_0, filename)
plot_fits(AWB_fit5_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_05, D_0, M_0, E_0, filename)

AWB_fit6 <- stan("AWB_adriana_pools5.stan", data = adriana_dat6, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit6_ex = extract(AWB_fit6)
AWB_fit6_ic = calc_ic(AWB_fit6, AWB_fit6_ex, S_06, D_0, M_0, E_0, filename)
fit_summary(AWB_fit6, AWB_fit6_ex, S_06, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit6, S_06, D_0, M_0, E_0, filename)
plot_fits(AWB_fit6_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_06, D_0, M_0, E_0, filename)

AWB_fit7 <- stan("AWB_adriana_pools5.stan", data = adriana_dat7, iter = 20000, refresh = 1, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
AWB_fit7_ex = extract(AWB_fit7)
AWB_fit7_ic = calc_ic(AWB_fit7, AWB_fit7_ex, S_07, D_0, M_0, E_0, filename)
fit_summary(AWB_fit7, AWB_fit7_ex, S_07, D_0, M_0, E_0, filename)

bayes_diagnostics(AWB_fit7, S_07, D_0, M_0, E_0, filename)
plot_fits(AWB_fit7_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_07, D_0, M_0, E_0, filename)
