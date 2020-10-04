library(ggplot2)
library(dplyr)
library(rstan)
library(gridExtra)
library(outbreaks)
library(bayesplot)
library(data.table)
library(knitr)
library(kableExtra)
load_data <- function(path="/data", location="cdmx"){
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/COVID19-VACINATION/r_sources"
    sub_path_3 <- "mcmc_parameter_estimation/UNISON-UADY-Project/data"
    file_name <- "cdmx_prevalence_data.csv"
    data_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/")
    covid19_data <- fread(data_path,
                                    select = c("FECHA_SINTOMAS",
                                              "i_s",
                                              "cumulative_i_s"))
    covid19_data <- data.frame(covid19_data)
#
    reference_date <- as.Date('2020-03-10')
    final_date_sample <- as.Date('2020-03-30')
    data_star_dynamics <- covid19_data %>%
        filter(as.Date(FECHA_SINTOMAS) >= reference_date &
                   as.Date(FECHA_SINTOMAS) <= final_date_sample)
    head(data_star_dynamics)
    data_plot <- ggplot(data = data_star_dynamics,
                        aes(x = FECHA_SINTOMAS, cumulative_i_s)) +
                    geom_bar(stat="identity", width = 0.05) +
                    geom_point() +
                    theme(axis.text.x = element_text(angle = 90)) +
                    ggtitle("CDMX")
    #
    sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
    sub_path_2 <- "NovelCovid19-ControlModelling/COVID19-VACINATION/r_sources"
    sub_path_3 <- "mcmc_parameter_estimation/UNISON-UADY-Project/plots"
    file_name <- "cdmx_input_data.pdf"
    plot_path <- 
        paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/")
    ggsave(plot_path)
    onset <- data_star_dynamics %>%
        select(FECHA_SINTOMAS)
    cum_cases <- data_star_dynamics %>%
        select(cumulative_i_s)

    cases <- data_star_dynamics %>%
        select(i_s)
    # cum_cases <- unlist(cum_cases, use.names = FALSE)
    names(cum_cases)[1] <-"cum_cases"
    data <- list(onset, cum_cases, data_star_dynamics)
    return (data)
}
