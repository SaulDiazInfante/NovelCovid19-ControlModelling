library(bayesplot)
library(ggplot2)
library(GGally)
library(dplyr)

mcmcm_post_analysis <- function(nuts_sample = nuts_fit){
    #' Create and save post analysis plots
    #' @param nuts_sample a list with the output of the sampling method
    #' @return trace, pairs, and density plots for evaluate the quality
    #' of the estimation
    #'
    posts <-  rstan::extract(nuts_fit)
    mod_diagnostics  <- rstan::get_sampler_params(nuts_fit)
    #
    # Check for divergent transitions
    #
    rstan::check_divergences(nuts_fit)
    posterior <- as.array(nuts_fit)
    color_scheme_set("viridisE")
    #
    # Markov chain traceplots
    #
    pars <- c("theta[1]",
            "theta[2]",
            "theta[3]",
            "theta[4]",
            "theta[5]",
            "theta[6]",
            "theta[7]",
            "theta[8]",
            "theta[9]",
            "theta[10]",
             "y_init[2]",
             "y_init[3]"
            )
    p1 <- mcmc_trace_highlight(posterior, pars = "lp__")
    p2 <- mcmc_trace_highlight(posterior, pars = "R_0")
    p3 <- mcmc_trace(posterior, pars = pars)
    #
#
    p4 <- mcmc_pairs(posterior,
                        pars = pars,
                        off_diag_args = list(size = 0.75),
                        np_style = pairs_style_np(
                                            div_color = "red",
                                            div_shape = 4,
                                            div_size = 1,
                                            div_alpha = 1,
                                            td_color = "yellow2",
                                            td_shape = 3,
                                            td_size = 1,
                                            td_alpha = 0.5))
#
    p5 <- stan_dens(nuts_fit, pars = pars, separate_chains = TRUE)
    p6 <- mcmc_intervals(posterior, pars = pars)
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/COVID19-VACINATION/r_sources"
    sub_path_3 <- "mcmc_parameter_estimation/UNISON-UADY-Project/plots"
    file_name_1 <- "trace_posterior.pdf.pdf"
    file_name_2 <- "trace_r_zero.pdf"
    file_name_3 <- "trace_parameters.pdf"
    file_name_4 <- "pairs.pdf"
    file_name_5 <- "posterior_densities.pdf"
    file_name_6 <- "intervals.pdf"
#    
    plot_path_1 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_1, sep = "/")
    plot_path_2 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_2, sep = "/")
    plot_path_3 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_3, sep = "/")
    plot_path_4 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_4, sep = "/")
    plot_path_5 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_5, sep = "/")
    plot_path_6 <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name_6, sep = "/")
    ggsave(plot_path_1, plot = p1)
    ggsave(plot_path_2, plot = p2)
    ggsave(plot_path_3 , plot = p3)
    ggsave(plot_path_4,
            plot = p4,
            width = 11,
            height = 8,
            dpi = 100,
            units = "in")
    ggsave(plot_path_5, plot = p5)
    ggsave(plot_path_6, plot = p6)
}
