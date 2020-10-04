library(tidyverse)
library(ggplot2)
# Load the flu da1aset of reported cases
sub_path_1 <-"/home/saul/sauld@cimat.mx/UNISON/Articles"
sub_path_2 <-"NovelCovid-19/datos_abiertos_covid19"
file_name <- "200728COVID19MEXICO.csv"
path <-paste(sub_path_1, sub_path_2, file_name, sep="/")
raw_data <- read.csv(path,
                     stringsAsFactors = TRUE,
                     header = TRUE,
                     sep = ",",
                     dec = ".")
df <- data.frame(raw_data)
working_data <- df %>%
    select(ID_REGISTRO,
           SEXO,
           ENTIDAD_NAC,
           ENTIDAD_UM,
           ENTIDAD_RES,
           MUNICIPIO_RES,
           TIPO_PACIENTE,
           FECHA_SINTOMAS,
           FECHA_DEF,
           EDAD,
           NEUMONIA,
           DIABETES,
           EPOC,
           ASMA,
           HIPERTENSION,
           CARDIOVASCULAR,
           OBESIDAD,
           TABAQUISMO,
           RESULTADO
    )
# identifiers
# CDMX: 09 DF
# SINALOA identifier: 25 SL
#   Culiacán: 006
# SONORA identifier: 26 SR
# Hermosillo clave_municipio: 030
#
reference_date <- as.Date('2020-03-01')
final_date_sample <- as.Date('2020-06-30')
raw_data_cdmx <- working_data %>%
    filter(as.numeric(ENTIDAD_RES) == 09,
               as.numeric(RESULTADO) == 1 &
               as.Date(FECHA_SINTOMAS) >= reference_date &
               as.Date(FECHA_SINTOMAS) <= final_date_sample)

prevalence_date_data_cdmx <- raw_data_cdmx %>%
    select(FECHA_SINTOMAS) %>%
    arrange(FECHA_SINTOMAS)
#
prevalence_date_data_cdmx <- raw_data_cdmx %>%
    group_by(FECHA_SINTOMAS)  %>%
    summarise(
        i_s = n()
    )
prevalence_date_data_cdmx <- prevalence_date_data_cdmx %>%
    mutate(cumulative_i_s = cumsum( i_s ))
#
sub_path_1 <- "/home/saul/sauld@cimat.mx/UNISON/Articles/NovelCovid-19"
sub_path_2 <- "NovelCovid19-ControlModelling/COVID19-VACINATION/r_sources"
sub_path_3 <- "mcmc_parameter_estimation/UNISON-UADY-Project/data/"
file_name <- "cdmx_prevalence_data.csv"
#
path <- paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/")
write.csv(prevalence_date_data_cdmx, path)
#
pl_0 <- ggplot(prevalence_date_data_cdmx,
                aes(x = as.Date(FECHA_SINTOMAS), y = i_s)) +
                geom_point(shape = 1) +
                geom_line() + 
                scale_x_date(
                                breaks = 
                                    function(x) 
                                        seq.Date(from = min(x),
                                                    to = max(x),
                                                    by = "1 weeks"),
                                minor_breaks = 
                                    function(x) 
                                        seq.Date(from = min(x),
                                                    to = max(x),
                                                    by = "2 weeks")
                ) +
                theme(axis.text.x = 
                        element_text(angle = 90, vjust = 0.5, hjust=1))

integer_offset_cdmx <-
    as.numeric(prevalence_date_data_cdmx$FECHA_SINTOMAS[1])
prevalence_data_cdmx <- prevalence_date_data_cdmx %>%
    mutate(integer_time = as.numeric(FECHA_SINTOMAS) - integer_offset_cdmx)

pl_1 <- ggplot(data=prevalence_data_cdmx,
                aes(x = integer_time, y = i_s)) +
                geom_point(shape = 1) +
                geom_line() +
                geom_point(data = prevalence_data_cdmx,
                shape = 3,
                color = "darkred",
                aes(x = integer_time, y=i_s))
write.csv(prevalence_date_data_cdmx, path)
pl_0
file_name <- "new cases_data_cdmx.pdf"
path <- paste(sub_path_1, sub_path_2, sub_path_3, file_name, sep = "/")
ggsave(path, plot=pl_0)