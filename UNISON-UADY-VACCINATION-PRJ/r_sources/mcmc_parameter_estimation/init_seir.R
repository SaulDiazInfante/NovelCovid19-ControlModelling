init_seirh_norm <- function(){
    list(theta = c(beta_s = rnorm(1, mean = 1.0, sd = 0.1),
                 beta_a = rnorm(1, mean = 1.0, sd = 0.1),
                 kappa = rgamma(1, 10, 40),
                 p = runif(1, 0.1, 0.5),
                 delta_h = rgamma(1, 10, 100),
                 mu_i_s = rexp(1, 0.03),
                 mu_h = rexp(1, 0.04),
                 gamma_s = rgamma(1, 10, 100),
                 gamma_a =rgamma(1, 10, 50),
                 gamma_h =rgamma(1, 10, 20)
                ),
            E0 = runif(1,  210 / pop, 560 / pop),
            IA0 = runif(1, 70 / pop, 630 / pop)
    )
}

init_seirh <- function(){
    list(theta = c(
                    beta_s = rnorm(1, mean = 1.0, sd = 0.1),
                    beta_a = rnorm(1, mean = 1.0, sd = 0.1),
                    kappa = rgamma(1, 10, 40),
                    p = runif(1, 0.1, 0.5),
                    delta_h = rgamma(1, 10, 100),
                    mu_i_s = rexp(1, 0.03),
                    mu_h = rexp(1, 0.04),
                    gamma_s = rgamma(1, 10, 100),
                    gamma_a =rgamma(1, 10, 50),
                    gamma_h =rgamma(1, 10, 20)
                ),
            E0 = runif(1,  74, 2100),
            IA0 = runif(1, 74, 2100)
    )
}
