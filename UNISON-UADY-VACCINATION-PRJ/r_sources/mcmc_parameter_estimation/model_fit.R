library(lubridate)
librry(ggthemr)
model_fit <- function(fit = nuts_fit, dates = onset){
    posts <-  rstan::extract(nuts_fit)
    mod_diagnostics  <- rstan::get_sampler_params(nuts_fit)
    color_scheme_set("viridisE")
    fit_CIS <- posts$y_hat[, , 8]
    fit_SIR <- fit_CIS
    median_I = apply(fit_SIR, 2, median)
    low_I = apply(fit_SIR, 2, quantile, probs = c(0.025))
    high_I = apply(fit_SIR, 2, quantile, probs = c(0.975))
    df_sample_N = data.frame(cum_cases, dates)
    df_fit_CIS = data.frame(median_I, low_I, high_I, dates)
    #
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/"
    sub_path_3 <- "NovelCovid19-ControlModellingGitHub"
    sub_path_4 <- "UNISON-UADY-VACCINATION-PRJ"
    sub_path_5 <- "r_sources/mcmc_parameter_estimation/data"
    file_name <- "data.Rda"
    data_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, sub_path_4, sub_path_5,
            file_name, sep = "/")
    save(df_sample_N, file =data_path)
        file_name <- "df_I_det_Poiss.Rda"
    data_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, sub_path_4, 
            file_name, sep = "/")
    save(df_fit_CIS, file = data_path)
    #
    ggthemr_reset()
    ggthemr('greyscale')
    plt <- ggplot(df_sample_N,
                    aes(x = as.Date(FECHA_SINTOMAS), y = cum_cases)) +
            geom_ribbon(aes(x = as.Date(FECHA_SINTOMAS),
                            ymin = low_I,
                            ymax = high_I),
                        fill = c_light_highlight,
                        alpha = 0.6) +
            geom_line(data = df_fit_CIS,
                        aes(x = as.Date(FECHA_SINTOMAS),
                            y = median_I,
                            color = "Median"),
                            size = 1.3) +
            geom_point(shape = 19,
                        size = 3,
                        (aes(color = "Data"))) +
            scale_colour_manual(name = 'dust',
                                values = c('Data' = 'black',
                                           'Median' = c_dark)) +
            guides(colour = 
                            guide_legend(override.aes = 
                                            list(shape = c(16, NA),
                                                    linetype = c(0, 1)))) +
            labs(x = "Time (days)",
                    y = "Cumulative Infected Cases") +
            # scale_x_continuous(limits = c(0, 22)) +
            scale_y_continuous(limits = c(0, 1000),
                       breaks = seq(from = 0, to = 1000, by = 200)) #+
            #theme_bw() + theme(text = element_text(size = 20))
    #
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/NovelCovid19-ControlModellingGitHub"
    sub_path_3 <- "UNISON-UADY-VACCINATION-PRJ"
    sub_path_4 <- "r_sources/mcmc_parameter_estimation/plots"
    file_name <- "cdmx_CIs_data_begining_fit.pdf"
    plot_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, sub_path_4,
            file_name, sep = "/")
    ggsave(plot_path)
}
