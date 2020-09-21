library(tidyverse)
library(ggplot2)
setwd("~/Insync/saul.diazinfante@unison.mx/OneDrive Biz/UNISON/ARTICLES/Covid19/COVID19-Sonora/r_sources")
# Load the flu da1aset of reported cases
raw_data <- read.csv(file = "./data/200706COVID19MEXICO.csv",
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
final_date_sample <- as.Date('2020-03-23')
raw_data_culiacan <- working_data %>%
    filter(as.numeric(ENTIDAD_RES) == 25,
           as.numeric(MUNICIPIO_RES) == 6 &
           as.numeric(ENTIDAD_UM) == 25 &
           as.numeric(RESULTADO) == 1 &
           as.Date(FECHA_SINTOMAS) >= reference_date &
           as.Date(FECHA_SINTOMAS) <= final_date_sample)
N <- 962871
prevalence_date_data_culiacan <- raw_data_culiacan %>%
    select(FECHA_SINTOMAS) %>%
    arrange(FECHA_SINTOMAS)
#
prevalence_date_data_culiacan <- raw_data_culiacan %>%
    group_by(FECHA_SINTOMAS)  %>%
    summarise(
        i_s = n()
    )
prevalence_date_data_culiacan <- prevalence_date_data_culiacan %>%
    mutate(cumulative_i_s = cumsum(i_s))
#
#
raw_data_hermosillo <- filter(working_data, MUNICIPIO_RES == 30)
prevalence_date_data_hermosillo <- raw_data_hermosillo %>%
    group_by(FECHA_SINTOMAS) %>%
    summarise(
        i_s = n()
    )
#
write.csv(prevalence_date_data_culiacan, './data/culiacan_prevalence_data.csv')
write.csv(prevalence_date_data_hermosillo,
          './data/hermosillo_prevalence_data.csv')
#
pl_0 <- ggplot(prevalence_date_data_culiacan, aes(x = FECHA_SINTOMAS, y = cumulative_i_s)) +
    geom_point(shape = 1) +
    geom_line()
    # geom_point(data = prevalence_date_data_hermosillo,
    #           shape = 3,
    #           color="darkred",
    #           aes(x=FECHA_SINTOMAS, y=i_s))
    # geom_line(data = prevalence_data_hermosillo)

integer_offset_culiacan <-
    as.numeric(prevalence_date_data_culiacan$FECHA_SINTOMAS[1])
prevalence_data_culiacan <- prevalence_date_data_culiacan %>%
    mutate(integer_time = as.numeric(FECHA_SINTOMAS) - integer_offset_culiacan)

integer_offset_hermosillo <-
    as.numeric(prevalence_date_data_hermosillo$FECHA_SINTOMAS[1])
prevalence_data_hermosillo <- prevalence_date_data_hermosillo %>%
    mutate(integer_time =
               as.numeric(FECHA_SINTOMAS) - integer_offset_hermosillo) %>%
    select(as.Date(FECHA_SINTOMAS) > reference_date )

pl_1 <- ggplot(data=prevalence_data_culiacan, aes(x = integer_time, y = i_s)) +
    geom_point(shape = 1) +
    geom_line() +
    geom_point(data = prevalence_data_hermosillo,
               shape = 3,
               color = "darkred",
               aes(x = integer_time, y=i_s))
write.csv(prevalence_data_culiacan, './data/culiacan_prevalence_data.csv')
write.csv(prevalence_data_hermosillo,
          './data/hermosillo_prevalence_data.csv')
pl_0
