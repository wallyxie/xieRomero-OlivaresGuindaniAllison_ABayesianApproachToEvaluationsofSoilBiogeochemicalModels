library(ggplot2,quietly=True)
library(rstan)
library(loo)
library(bayesplot)
library(coda)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

filename = "CON_adriana_pools5i_vary_mic"
color_scheme_set("viridisA")

#############
##FUNCTIONS##
#############

#Save summary, sample parameters, and posteriors
fit_summary <- function(stan_fit, stan_fit_ex, S_0, D_0, M_0, filename) {
    #Posteriors
    write.csv(stan_fit_ex, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "posteriors", "S", S_0, "D", D_0, "M", M_0, ".csv", sep = "_"))
    #Summary
    stan_fit_sum <- summary(stan_fit)$summary
    stan_fit_sum_write <- stan_fit_sum[c("Ea_S", "Ea_D", "Ea_M", "a_DS", "a_SD", "a_M", "a_MS", "sigma"), c("mean", "sd", "2.5%", "50%", "97.5%", "n_eff", "Rhat")]
    write.csv(stan_fit_sum_write, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "summary", "S", S_0, "D", D_0, "M", M_0, ".csv", sep = "_"))
    #Divergent Transitions
    div_trans <- sum(nuts_params(stan_fit, pars = "divergent__")$Value)
    write.csv(div_trans, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "div_trans", "S", S_0, "D", D_0, "M", M_0, ".csv", sep = "_"))
    #Sampler parameters
    sampler_params <- get_sampler_params(stan_fit, inc_warmup = FALSE) #Access sampler values
    write.csv(sampler_params, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "sampler_params", "S", S_0, "D", D_0, "M", M_0, ".csv", sep = "_"))
}

#IC calculations
calc_ic <- function(stan_fit, stan_fit_ex, S_0, D_0, M_0, filename) {
    fit_log_lik <- extract_log_lik(stan_fit)
    #r_eff <- relative_eff(exp(fit_log_lik), cores = 4)
    stan_fit_waic <- waic(fit_log_lik, cores = 4)
    stan_fit_loo <- loo(fit_log_lik, cores = 4)
    LPML <- sum(log(1 / colMeans(stan_fit_ex$CPOi)))
    waic <- stan_fit_waic$estimates["waic",]
    p_waic <- stan_fit_waic$estimates["p_waic",]
    loo <- stan_fit_loo$estimates["looic",]
    p_loo <- stan_fit_loo$estimates["p_loo",]
    sink(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "loo_print", "S", S_0, "D", D_0, "M", M_0, ".txt", sep = "_"))
    print(stan_fit_loo)
    sink()
    sink(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "pareto-k", "S", S_0, "D", D_0, "M", M_0, ".txt", sep = "_"))
    print(pareto_k_table(stan_fit_loo))
    sink()
    sink(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "ic", "S", S_0, "D", D_0, "M", M_0, ".txt", sep = "_"))    
    cat("WAIC = ", waic, "\nLOO = ", loo, "\nLPML = ", LPML, "\np_waic = ", p_waic, "\np_loo = ", p_loo)
    sink()    
    return(c(waic, loo, LPML, p_waic, p_loo))
}

#Convert Stan object into coda object
stan2coda <- function(fit) {
    mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
}

#Diagnostics including posterior, traceplot, autocorrelation, and Rhat. Specific to CON because of parameters. Cannot copy-pate to other models
bayes_diagnostics <- function(stan_fit, S_0, D_0, M_0, filename) {
    stan_fit.array <- as.array(stan_fit)
    stan_fit.np <- nuts_params(stan_fit)
    stan_fit.lp <- log_posterior(stan_fit)
    #Divergent transitions plots
    mcmc_ND <- mcmc_nuts_divergence(stan_fit.np, stan_fit.lp)
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "mcmc_ND", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = mcmc_ND)
    mcmc_PC <- mcmc_parcoord(stan_fit.array, np = stan_fit.np)
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "mcmc_PC", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = mcmc_PC)
    #MCMC NUTS energy
    mcmc_ne <- mcmc_nuts_energy(stan_fit.np)
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "mcmc_ne", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = mcmc_ne)   
    #Posterior credible areas
    CI_plot1 <- mcmc_areas(stan_fit.array, pars = c("Ea_S", "Ea_D", "Ea_M"), prob = 0.8, prob_outer = 0.95) + yaxis_text() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "CI_Ea", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = CI_plot1)
    CI_plot2 <- mcmc_areas(stan_fit.array, pars = c("a_DS", "a_SD", "a_M", "a_MS"), prob = 0.8, prob_outer = 0.95) + yaxis_text() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "CI_a", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = CI_plot2)
    #Pairs - Does not work with too many iterations
    # pairs_plot1 <- mcmc_pairs(stan_fit.array, diag_fun = "dens")
    # ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "pairs", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = pairs_plot1)
    # ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "pairs", "S", S_0, "D", D_0, "M", M_0, ".png", sep = "_"), plot = pairs_plot1)
    #Traceplot
    # traceplot1 <- rstan::traceplot(stan_fit, pars = c("Ea_S", "Ea_D", "Ea_M"))
    # ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "trace_Ea", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = traceplot1)
    # traceplot2 <- rstan::traceplot(stan_fit, pars = c("a_DS", "a_SD", "a_M", "a_MS"))
    # ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "trace_a", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = traceplot2)
    trace_plot1 <- mcmc_trace(stan_fit.array, pars = c("Ea_S", "Ea_D", "Ea_M"), np = stan_fit.np, facet_args = list(ncol = 1, strip.position = "left")) + xlab("Post-warmup iteration") + yaxis_text() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "trace_Ea", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = trace_plot1)
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "trace_Ea", "S", S_0, "D", D_0, "M", M_0, ".png", sep = "_"), plot = trace_plot1)
    trace_plot2 <- mcmc_trace(stan_fit.array, pars = c("a_DS", "a_SD", "a_M", "a_MS"), np = stan_fit.np, facet_args = list(ncol = 1, strip.position = "left")) + xlab("Post-warmup iteration") + yaxis_text() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "trace_a", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = trace_plot2)
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "trace_a", "S", S_0, "D", D_0, "M", M_0, ".png", sep = "_"), plot = trace_plot2)
    #Autocorrelation
    acf_plot1 <- mcmc_acf(stan_fit.array, pars = c("Ea_S", "Ea_D", "Ea_M")) + yaxis_text() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "acf_Ea", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = acf_plot1)
    acf_plot2 <- mcmc_acf(stan_fit.array, pars = c("a_DS", "a_SD", "a_M", "a_MS")) + yaxis_text() + theme(axis.text = element_text(size = 16), axis.title = element_text(size = 20))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "acf_a", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = acf_plot2)
    #Rhat plot
    stan_fit_coda <- stan2coda(stan_fit)
    pdf(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "Rhat", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"))
    gelman.plot(stan_fit_coda[,1:8])
    dev.off()
}

#Pearson's correlation calculation
calc_r <- function(model_vector, data_vector, S_0, D_0, M_0, filename) {
    sst_mean <- mean(data_vector)
    sst_vector <- (data_vector - sst_mean) ^ 2
    sst <- sum(sst_vector)
    ssr_vector <- (data_vector - model_vector) ^ 2
    ssr <- sum(ssr_vector)
    r_sq <- 1 - (ssr / sst)
    cat(r_sq, file = paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "r_sq", "S", S_0, "D", D_0, "M", M_0, ".txt", sep = "_"))
    return(r_sq)
}

#Plot posterior fit and forward simulation prediction
plot_fits <- function(stan_fit_ex, N_t, N_p, obs_times, pred_times, data_vector, S_0, D_0, M_0, filename) {
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
    r_sq <- calc_r(CO2_flux_ratios_hat_median, data_vector, S_0, D_0, M_0, filename)
    #Plotting posterior fit
    df_post <- data.frame(list(obs_times = obs_times / (24 * 365), CO2_flux_ratios_vector = data_vector, CO2_flux_ratios_hat_median = CO2_flux_ratios_hat_median))
    post_plot <- ggplot(df_post, aes(x = obs_times)) +
    geom_ribbon(aes(ymin = CO2_flux_ratios_hat_median - 1.96 * sigma_hat, ymax = CO2_flux_ratios_hat_median + 1.96 * sigma_hat), fill = "bisque") +
    geom_line(aes(y = CO2_flux_ratios_hat_median), colour = "black", size = 2) +
    geom_point(aes(y = CO2_flux_ratios_vector), shape = 21, colour = "black", fill = "white", size = 3) +
    geom_errorbar(interval_95pct, width = 0, color = "blue") +    
    labs(x = "Duration (years)",
         y = "Response Ratio") +
    ggtitle("CON Response Ratio vs. Year")
    post_plot <- post_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16), axis.title = element_text(size = 20)) + coord_cartesian(xlim = c(0, 12.5), ylim = c(0.5, 2.0))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "post", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = post_plot)    
    #Plotting predictive fit
    df_pre <- data.frame(list(pred_times = pred_times / (24 * 365), CO2_flux_ratios_pred_median = CO2_flux_ratios_pred_median))
    pre_plot <- ggplot(df_pre, aes(x = pred_times)) +
    geom_ribbon(aes(ymin = CO2_flux_ratios_pred_median - 1.96 * sigma_hat, ymax = CO2_flux_ratios_pred_median + 1.96 * sigma_hat), fill = "bisque") +
    geom_line(aes(y = CO2_flux_ratios_pred_median), colour = "black", size = 2) +
    labs(x = "Duration (years)",
         y = "Response Ratio") +
    ggtitle("CON Predictive Response Ratio vs. Year")
    pre_plot <- pre_plot + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                   panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(size = 16), axis.title = element_text(size = 20)) + coord_cartesian(xlim = c(0, 60), ylim = c(0.5, 2.0))
    ggsave(paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "pred", "S", S_0, "D", D_0, "M", M_0, ".pdf", sep = "_"), plot = pre_plot)   
}

########
##DATA##
########

df <- read.csv("adriana_data_s_min.csv",header=TRUE) #Read in data
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

S_0 <- 100;
D_0 <- 0.2
M_01 <- 1; M_02 <- 2; M_03 <- 3; M_04 <- 4; M_05 <- 5; M_06 <- 6; M_08 <- 8;
u_M <- 0.002
k_M_ref1 <- u_M * D_0 / M_01
k_M_ref2 <- u_M * D_0 / M_02
k_M_ref3 <- u_M * D_0 / M_03
k_M_ref4 <- u_M * D_0 / M_04
k_M_ref5 <- u_M * D_0 / M_05
k_M_ref6 <- u_M * D_0 / M_06
k_M_ref8 <- u_M * D_0 / M_08

CON_dat1 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_01, u_M = u_M, k_M_ref = k_M_ref1)

CON_dat2 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_02, u_M = u_M, k_M_ref = k_M_ref2)

CON_dat3 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_03, u_M = u_M, k_M_ref = k_M_ref3)

CON_dat4 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_04, u_M = u_M, k_M_ref = k_M_ref4)

CON_dat5 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_05, u_M = u_M, k_M_ref = k_M_ref5)

CON_dat6 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_06, u_M = u_M, k_M_ref = k_M_ref6)

CON_dat8 <- list(N_t = N_t, N_p = N_p, t0 = t0, ts = hour_index_list, ts_p = ts_p, CO2_flux_ratios_vector = CO2_flux_ratios_vector, SOC_input = mean_SOC_input, DOC_input = mean_DOC_input, T_ref = T_ref, R_g = R_g, S_0 = S_0, D_0 = D_0, M_0 = M_08, u_M = u_M, k_M_ref = k_M_ref8)

file_path <- "CON_adriana_pools5i.stan" #Read in Stan model code. Stan file must be in same directory.
lines <- readLines(file_path, encoding = "ASCII")
for (n in 1:length(lines)) cat(lines[n],'\n')

#############
##EXECUTION##
#############

CON_fit1 <- stan("CON_adriana_pools5i.stan", data = CON_dat1,  iter = 35000, warmup = 10000, refresh = 10, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 12))
saveRDS(CON_fit1, paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "stanfit", "S", S_0, "D", D_0, "M", M_01, ".rds", sep = "_"))
CON_fit1_ex <- extract(CON_fit1)
CON_fit1_ic <- calc_ic(CON_fit1, CON_fit1_ex, S_0, D_0, M_01, filename)
fit_summary(CON_fit1, CON_fit1_ex, S_0, D_0, M_01, filename)
bayes_diagnostics(CON_fit1, S_0, D_0, M_01, filename)
plot_fits(CON_fit1_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_0, D_0, M_01, filename)

CON_fit3 <- stan("CON_adriana_pools5i.stan", data = CON_dat3,  iter = 35000, warmup = 10000, refresh = 10, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 12))
saveRDS(CON_fit3, paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "stanfit", "S", S_0, "D", D_0, "M", M_03, ".rds", sep = "_"))
CON_fit3_ex <- extract(CON_fit3)
CON_fit3_ic <- calc_ic(CON_fit3, CON_fit3_ex, S_0, D_0, M_03, filename)
fit_summary(CON_fit3, CON_fit3_ex, S_0, D_0, M_03, filename)
bayes_diagnostics(CON_fit3, S_0, D_0, M_03, filename)
plot_fits(CON_fit3_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_0, D_0, M_03, filename)

CON_fit4 <- stan("CON_adriana_pools5i.stan", data = CON_dat4,  iter = 35000, warmup = 10000, refresh = 10, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 12))
saveRDS(CON_fit4, paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "stanfit", "S", S_0, "D", D_0, "M", M_04, ".rds", sep = "_"))
CON_fit4_ex <- extract(CON_fit4)
CON_fit4_ic <- calc_ic(CON_fit4, CON_fit4_ex, S_0, D_0, M_04, filename)
fit_summary(CON_fit4, CON_fit4_ex, S_0, D_0, M_04, filename)
bayes_diagnostics(CON_fit4, S_0, D_0, M_04, filename)
plot_fits(CON_fit4_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_0, D_0, M_04, filename)

CON_fit5 <- stan("CON_adriana_pools5i.stan", data = CON_dat5,  iter = 35000, warmup = 10000, refresh = 10, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 12))
saveRDS(CON_fit5, paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "stanfit", "S", S_0, "D", D_0, "M", M_05, ".rds", sep = "_"))
CON_fit5_ex <- extract(CON_fit5)
CON_fit5_ic <- calc_ic(CON_fit5, CON_fit5_ex, S_0, D_0, M_05, filename)
fit_summary(CON_fit5, CON_fit5_ex, S_0, D_0, M_05, filename)
bayes_diagnostics(CON_fit5, S_0, D_0, M_05, filename)
plot_fits(CON_fit5_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_0, D_0, M_05, filename)

CON_fit6 <- stan("CON_adriana_pools5i.stan", data = CON_dat6,  iter = 35000, warmup = 10000, refresh = 10, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 12))
saveRDS(CON_fit6, paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "stanfit", "S", S_0, "D", D_0, "M", M_06, ".rds", sep = "_"))
CON_fit6_ex <- extract(CON_fit6)
CON_fit6_ic <- calc_ic(CON_fit6, CON_fit6_ex, S_0, D_0, M_06, filename)
fit_summary(CON_fit6, CON_fit6_ex, S_0, D_0, M_06, filename)
bayes_diagnostics(CON_fit6, S_0, D_0, M_06, filename)
plot_fits(CON_fit6_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_0, D_0, M_06, filename)

CON_fit8 <- stan("CON_adriana_pools5i.stan", data = CON_dat8,  iter = 35000, warmup = 10000, refresh = 10, chains = 4, seed = 1234, open_progress = "False", control = list(adapt_delta = 0.95, stepsize = 0.1, max_treedepth = 12))
saveRDS(CON_fit8, paste(format(Sys.time(),"%Y_%m_%d_%H_%M"), filename, "stanfit", "S", S_0, "D", D_0, "M", M_08, ".rds", sep = "_"))
CON_fit8_ex <- extract(CON_fit8)
CON_fit8_ic <- calc_ic(CON_fit8, CON_fit8_ex, S_0, D_0, M_08, filename)
fit_summary(CON_fit8, CON_fit8_ex, S_0, D_0, M_08, filename)
bayes_diagnostics(CON_fit8, S_0, D_0, M_08, filename)
plot_fits(CON_fit8_ex, N_t = N_t, N_p = N_p, obs_times = hour_index_list, pred_times = ts_p, data_vector = CO2_flux_ratios_vector, S_0, D_0, M_08, filename)
