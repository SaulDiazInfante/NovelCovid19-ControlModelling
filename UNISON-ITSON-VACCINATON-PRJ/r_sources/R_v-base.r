library(plotly)
#
#
epsilon_bar   <- 0.3905
beta_a  <- 0.6441 * epsilon_bar
beta_s  <- 0.9303 * epsilon_bar
p  <-  0.1213
alpha_s  <-  1/10.81
alpha_a  <-  1/5.97
mu      <- 0.00003653
delta_e   <- 1/5.1
delta_r   <- 1/365
delta_v   <- 1/180
#
R0 <- 
    (p * delta_e * beta_s) / 
    (
        (mu + delta_e) * (mu + alpha_s)
    ) + ((1 - p) * delta_e * beta_a) / 
    (
        (mu + alpha_a) * (mu + delta_e)
    )
epsilon <- seq(0,1,0.001)
lambda_v   <- seq(0,0.009,0.000009)
# r_v generation data
Rv <- matrix(data = 0.0, 
nrow = length(epsilon), ncol = length(lambda_v))
#
for (i in seq(1,length(epsilon),1)){
    for (j in seq(1,length(lambda_v),1)){
        reduction_factor <- 
            1.0 -  
                (epsilon[i] * lambda_v[j]) / (mu + delta_v + lambda_v[j])
        Rv[i,j] <- R0 * reduction_factor
    #  
  }
}
#
# Figures R0 
f <- list(
  family = "Arial",
  size = 14
)
x_a <- list(
  title = "Vaccine Efficacy",
  titlefont = f
)
y_a <- list(
  title = "Vaccination Rate",
  titlefont = f
)
#
data_1  <- data.frame(epsilon, lambda_v, Rv)
#
#Contour figures
fig1 <- plot_ly(showscale = TRUE)
fig1 <- fig1 %>% 
    add_contour(
        data_1,
        x = ~epsilon,
        y = ~lambda_v,
        z = ~Rv,
        contours = 
            list(
                showlabels = TRUE,
                labelfont = list(size = 15, color = 'white')
            ),
        showlegend = TRUE
    )
fig1 <- fig1 %>% 
    layout(xaxis = x_a, yaxis = y_a)
fig1 <- fig1 %>% 
    add_lines(
        x = ~epsilon,
        y = 0.000611352,
        line = list(color = "lightgray", width = 2),
        showlegend = FALSE
    )
fig1 <- fig1 %>%
    add_lines(x = ~epsilon, 
        y = 0.0071,
        line = list(color = "white", width = 2),
        showlegend = FALSE
    )
fig1 <- fig1 %>% 
    add_lines(
        x = 0.8, 
        y = ~lambda_v,
        line = list(color = "gray", width = 2),
        showlegend = FALSE
    )
# fig1
orca(fig1, "R0_contour.png")
