# density parameters
library(ggplot2)
library(plotly)

mu_h <- function(){
    # We suppose mu_. as a exponential r.v.
    lambda <- 25
    x <- seq(from = 0, to = 1, by=.01)
    mu_h_ <- dexp(x, rate = lambda)
    df <- data.frame(x=x, mu_h = mu_h_)
    p <- ggplot(data=df,
                aes(x = x,
                    y = mu_h)) + 
        geom_line(colour = "#000000") +
        geom_area(alpha = 0.6, fill = "lightgray")
        
    fig <- ggplotly(p)
    fig
}

mu_h()
gamma_s <- function(){
    # We suppose mu_. as a exponential r.v.
    alpha <- 10
    theta <- 50
    beta <- 1/theta
     
    x <- seq(from = 0, to = 1, by=.001)
    gamma_s_ <- dgamma(x, shape = alpha, 
                        scale= beta)
    df <- data.frame(x=x, gamma_s = gamma_s_)
    p <- ggplot(data=df,
                aes(x = x,
                    y = gamma_s)) + 
        geom_line(colour = "#000000") +
        geom_area(alpha = 0.6, fill = "lightgray")
        
    fig <- ggplotly(p)
    fig
}
gamma_s()a
