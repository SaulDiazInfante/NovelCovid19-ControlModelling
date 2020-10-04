library(lubridate)
mcmcm_stain_summary <- function(fit = nuts_fit){
    # Set initial values:
    parameters = c("theta", "R_0", "y_init")
    nuts_fit_summary <- summary(nuts_fit, pars = parameters)$summary
    df_pars <- data.frame(nuts_fit_summary)
    names(df_pars)[4] <-"2.5%"
    names(df_pars)[5] <-"25%"
    names(df_pars)[6] <-"50%"
    names(df_pars)[7] <-"75%"
    names(df_pars)[8] <-"97.5%"
    colnames(df_pars) <- c('mean', 'se_mean','sd',
                            '2.5%', '25%','50%',
                            '75%' ,'97.5%', 'n_eff', 'Rhat')
    row.names(df_pars) <- c("beta_s", "beta_a", "kappa",
                            "p", "delta_H", 
                            "mu_{I_S}", "mu_H", 
                            "gamma_S", "gamma_A", "gamma_H",
                            "R_zero",
                            "S_0", "E_0", "I_{S_0}", "I_{A_0}", 
                            "H_0", "R_0", "D_0", "ICS_0")
    print(df_pars)
    prefix_time <- Sys.time()
    prefix_time <- paste(date(prefix_time),
                    hour(prefix_time),
                    minute(prefix_time), sep="_")
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/COVID19-VACINATION/r_sources"
    sub_path_3 <- "mcmc_parameter_estimation/UNISON-UADY-Project/estimation"
    file_name <- paste("parameters", prefix_time,".csv", sep="")
    path <- paste(sub_path_1,
                    sub_path_2,
                    sub_path_3,
                    file_name,
                    sep = "/")
    path_last <- paste(sub_path_1,
                        sub_path_2,
                        sub_path_3, 
                        "last_est.csv",
                            sep = "/")
    write.csv(df_pars, path)
    write.csv(df_pars, path_last)
}
