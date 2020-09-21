rm(list = ls())
require(deSolve)
library(Rtwalk)
library(fields)
library(rgl)
library(nlme)
library(caTools)

#### Ejemplo de densidades a priori para los parámetros ####
x11();
par(mfrow = c(1, 2), mar = c(4, 5, 4, 1))
curve(dgamma(x, 2, 2),
      xlab = expression(beta),
      ylab = expression( pi(beta)),
      lwd = 2,
      xlim = c(0,3),
      cex.lab = 1.5,
      cex.main = 1.5,
      main = "Gam(2,2)",
      col = "navyblue")

curve(dgamma(x,1.5,1.5),
      xlab = expression(gamma),
      ylab = expression(pi(gamma)),
      lwd = 2,
      xlim = c(0,3),
      cex.lab = 1.5,
      cex.main = 1.5,
      main = "Gam(1.5,1.5)",
      col = "navyblue")

#### Solver del modelo SIR determinista ####
X_theta <- function(theta, t = tiempos){
    SIRmod <- function(t, x, theta) {
        with(as.list(c(theta, x)),
             {
                 ds <- -bet * s * i
                 di <- bet * s * i - gam * i
                 dr <- gam * i
                 res <- c(ds, di, dr)
                 list(res) }
        ) }
    ## Solver
    out <- lsoda(X_ini, t, SIRmod, theta)
    out[which(out[, 3] < 0), 3] <- 0
    return(out)
}

#### Integral numérica de la tasa de nuevos infectados entre dos tiempos
#### s * i
spori <- function(bet0, gam0, t){
    ifelse(t == 0,
           si <- as.numeric(X_ini["s"] * X_ini["i"]),
           si <- prod(X_theta(c(bet = bet0, gam = gam0), c(0,t))[-1,2:3]))
    si
}
#### integral_ ti ^ tf N * beta * si ####

parampois <- function(bet0, gam0, ti, tf){
    rejilla <- seq(ti, tf, length = 15)
    nbetasi <- N * bet0 * Vectorize(spori)(bet0, gam0, rejilla)
    return(trapz(rejilla, nbetasi))
}
#### Logverosimilitud ####

logver <- function(bet0, gam0, dat = y){
    lv <- vector()
    for(i in 1:(length(tiempos) - 1))
        lv[i] <- dpois(y[i],
                       parampois(bet0, gam0, tiempos[i], tiempos[i + 1]),
                      log = T)
    return(sum(lv))
}

#### A prioris ####
#### Poner una Gamma(2,2) para beta y Gamma(1.5,1.5) para gamma ####

logaprioris <- function(bet0,
                        gam0,
                        id,
                        a_b = 2,
                        b_b = 2,
                        a_g = 1.5,
                        b_g = 1.5){
    if (id == 1)
        logap <-
            sum(c(dgamma(bet0, a_b, b_b, log = T),
                  dgamma(gam0, a_g, b_g, log = T))
            )
    else logap <- sum(c(log(1/4) +
                        log(bet0 > 0 & bet0 < 4),
                        log(1/4) + log(gam0 > 0 & gam0 < 4)))
    return(logap)
}

#### Log posterior ####
logpost <- function(bet0, gam0, id){
    theta <- c(bet0, gam0)
    lp <- logver(theta[1], theta[2]) + logaprioris(theta[1], theta[2], id = id)
    if (is.nan(lp))
        lp <- -Inf
    return(suppressWarnings(lp))
}
#### Posterior ####
post <- function(bet0, gam0, id){
    ps <- suppressWarnings(exp( logpost(bet0, gam0, id = id)))
    ifelse(is.nan(ps) || is.infinite(ps),
           out <- 0,
           out <- ps
    )
    return(out)
}
# Para colores
add.alpha <- function(col, alpha=1){
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb) / 255, 2,
          function(x)
              rgb(x[1], x[2], x[3], alpha = alpha))
}

#### Integrated Autocorrelation Time ####
IAT <- function(matsims, from=1, to=dim(matsims)[1]) {
    dim = dim(matsims)[2]
    # par elige los parámetros relevantes a considerar
    ##-lag/log(GetAutoCorr( info, lag, par=par, from=from, to=to))
    IATs <- c()
    for (k in 1:dim) {
        dt <-  as.matrix(matsims[from:to, k], ncol = 1)
        n <- to - from
        mu <- mean(dt)  ### with its mean and variance
        s2 <- var(dt)
        ### The maximum lag is half the sample size
        maxlag <- max( 3, floor(n/2))
        #### The gammas are sums of two consecutive autocovariances
        Ga <- rep(0, 2)  ## two consecutive gammas
        lg <- 0
        Ga[1] <- s2  #sum( (dt[1:(n-lg)]-mu)*(dt[(lg+1):n]-mu) )/n
        lg <- 1
        Ga[1] <- Ga[1] +
            sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu) ) / n
        m <- 1
        lg <- 2 * m
        Ga[2] <-
            sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu)) / (n - lg)
        lg <- 2 * m + 1
        Ga[2] <- Ga[2] +
            sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu) ) / (n - lg)
        IAT <- Ga[1] / s2
        ### Add the autocorrelations
        ### RULE: while Gamma stays positive and decreasing
        while ((Ga[2] > 0.0) & (Ga[2] < Ga[1])) {
            m <- m + 1
            if (2 * m + 1 > maxlag) {
                cat("Not enough data, maxlag=", maxlag, "\n")
                break
            }
            Ga[1] <- Ga[2]

            lg <- 2 * m
            Ga[2] <- sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu))/n
            lg <- 2 * m + 1
            Ga[2] <- Ga[2] +
                sum((dt[1:(n - lg)] - mu) * (dt[(lg + 1):n] - mu) ) / n
            IAT <- IAT + Ga[1] / s2
        }
        IAT <- -1 + 2 * IAT   ##Calculates the IAT from the gammas
        #cat("IAT: IAT=", IAT, ", last lag=", 2*m+1, ", last Gamma=", Ga[1], "\n")

        IATs[k] <- IAT }
    return(mean(IATs))
}
#### Valores iniciales para poder simular de la densidad posterior.####
# Valores reales de los parámetros
betreal = 0.55;
gamreal = 0.25;

#### Secuencias de los parámetros a evaluar ####
seqbeta <- seq(0.1, 2, length = 50)
seqgamma <- seq(0.1, 2, length = 50)
seqbeta_ <- seq(0.45, 0.68, length = 100)
seqgamma_ <- seq(0.1, 0.4, length = 100)

# Tamaño de la cadena a simular
nsim = 10;
# Tiempos considerados para las simulaciones
tiempos = seq(0, 30, length = 31);
# Vector de valores iniciales para el sistema de ecuaciones diferenciales
X_ini = c(s = 0.95, i = 0.05, r = 0);
# Tamaño total de la población
N = 500
# Simular datos observados
y <- vector()
set.seed(8)
for (i in 1:(length(tiempos) - 1))
    y[i] <- rpois(1,
                  parampois(betreal, gamreal, tiempos[i], tiempos[i + 1])
            )
#### Gráfica de la logverosimilitud de los datos ####
lv <- outer(seqbeta_, seqgamma_, Vectorize(logver))
lv[which(is.infinite(lv), arr.ind = T)] <- min(lv) - 500
image.plot(seqbeta_, seqgamma_, lv)
points(betreal, gamreal, pch = 16, col = "white")
#### Log a prioris 1 ####
lap1 <- outer(seqbeta, seqgamma, Vectorize(logaprioris), id = 1)
image.plot(seqbeta, seqgamma, lap1)
points(betreal, gamreal, pch = 16, col = "white")
#### Log a prioris 2 ####
lap2 <- outer(seqbeta, seqgamma, Vectorize(logaprioris), id = 2)
lap2[which(is.infinite(lap2), arr.ind = T)] <- 0
image.plot(seqbeta, seqgamma, lap2)
contour(seqbeta, seqgamma, lv, add = T)
points(betreal, gamreal, pch = 16, col = "white")

# Constante de integración 1
K1 <- integrate( function(y) {
    sapply(y, function(y) {
        integrate(function(x) Vectorize(post)(x, y, id = 1), 0, 5)$value
    })
    }, 0, 5)$value

# Logposterior 1
lp1 <- outer(seqbeta_, seqgamma_, Vectorize(logpost), id = 1)
image.plot(seqbeta_, seqgamma_, lp1)

#### A prioris 1 ####
par(mfrow = c(1, 2), mar = c(5, 5, 3, 2))
curve(dgamma(x, 2, 2),
      xlab = expression(beta),
      ylab = expression(pi(beta)),
      cex.lab = 1.5,
      xlim = c(0, 2),
      lwd = 2,
      main = "Gam(2,2)")
####
curve(dgamma(x, 1.5, 1.5),
      xlab = expression(gamma),
      ylab = expression(pi(gamma)),
      cex.lab = 1.5,
      xlim = c(0, 2),
      lwd = 2,
      main = "Gam(1.5,1.5)")

# Posterior 1
post1 <- outer(seqbeta_, seqgamma_, Vectorize(post), id = 1)
image.plot(seqbeta_,
           seqgamma_,
           post1 / K1, xlab = expression(beta),
           ylab = expression(gamma))

#### Constante de integración 2 ####
K2 <- integrate( function(y) {
    sapply(y, function(y) {
        integrate(function(x) Vectorize(post)(x, y, id = 2), 0, 5)$value
    }) }, 0, 5)$value

#### Log posterior 2 ####
lp2 <- outer(seqbeta_, seqgamma_, Vectorize(logpost), id = 2)
image.plot(seqbeta_, seqgamma_, lp2)

# Posterior 2
post2 <- outer(seqbeta_, seqgamma_, Vectorize(post), id = 2)
image.plot(seqbeta_, seqgamma_, post2 / K2)
#### Simulaciones de la posterior 1 ####
#### Simular un punto inicial ####
x0 = c(rexp(1), rexp(1))
#### Soporte 1 ####
sup1 <- function(theta){
    if (theta[1] > 0 & theta[2] > 0) return(TRUE)
    else(return(FALSE))
}
#### -logposterior 1 ####
menoslogpost1 <- function(theta)
    return(-logpost(theta[1], theta[2], id = 1))

# t-walk con aprioris Gamma
simtwalk1 <- Runtwalk(10000,
                      dim = 2,
                      Obj = menoslogpost1,
                      Supp = sup1,
                      x0 = x0,
                      xp0 = c(rexp(1),rexp(1)))
# Para ver convergencia
x11()
plot(0:100,
     -simtwalk1$Us[1:101],
     type = "l",
     xlab = "Iteración",
     ylab = "Logposterior")

burnin1 <- 30
abline(v = burnin1, col = "navyblue", lwd = 2, lty = 2)

# Para ver dependencia
Ana(info = simtwalk1)
IAT(simtwalk1$output)
tau1 <- ceiling(IAT(simtwalk1$output))
(tamrealmues1 <- round((nsim - burnin1) / tau1, 0))
#415
cadsinrez <- simtwalk1$output[burnin1:nsim, ]
cadsinrez <- cadsinrez[seq(1, dim(cadsinrez)[1], tau1),]

# Contornos y puntos simulados
x11();
par(mar = c(4, 5, 1, 1))
image.plot(seqbeta_,
           seqgamma_,
           post1 / K1,
           xlim = c(0.48, 0.68),
           ylim = c(0.1, 0.4),
           xlab = expression(beta),
           ylab = expression(gamma),
           cex.lab = 2)
points(cadsinrez,
       pch = 16,
       cex = 0.6)
points(betreal,
       gamreal,
       pch = 16,
       col = "white",
       cex = 1)

x11()
par(mfrow = c(1, 2), mar = c(4, 5, 2, 2))
hist(cadsinrez[,1],
     xlab = expression(beta),
     freq = F,
     ylab = "Densidad",
     col = add.alpha("gray",0.4),
     border = "gray",
     main = "",
     cex.lab = 1.5)
####
abline(v = betreal,
       lwd = 2,
       col = "navyblue",
       lty = 2)
####
hist(cadsinrez[,2],
     xlab = expression(gamma),
     freq = F,
     ylab = "Densidad",
     col = add.alpha("gray", 0.4), border = "gray", main ="", , cex.lab = 1.5)
abline(v = gamreal, lwd = 2, col = "navyblue", lty = 2, ylab = "Densidad")

require(xtable)
xtable(apply(cadsinrez, 2, quantile, prob = c(0.025, 0.5, 0.975)),
       digits = c(0, 4, 4))
# Simulaciones de la posterior 2
# Simular un punto inicial
x02 = c(rexp(1), rexp(1))
# Soporte 2
sup2 <- function(theta){
    if (theta[1] > 0 & theta[2] > 0 & theta[1] < 4 & theta[2] < 4)
        return(TRUE)
    else(return(FALSE))
}

# -logposterior 2
menoslogpost2 <- function(theta)
    return(-logpost(theta[1], theta[2], id = 2))

# Simulaciones con t-walk (a prioris uniformes)
simtwalk2 <- Runtwalk(10000, dim = 2, Obj = menoslogpost2,
                    Supp = sup2, x0 = x02, xp0 = c(rexp(1),rexp(1)))



# Para ver convergencia
x11()
plot(0:(nsim/4), simtwalk2$Us[1:(nsim / 4 + 1)], type = "l")
burnin2 <- 200
abline(v = burnin2, col = "tomato", lwd = 2)
# Para ver dependencia
Ana(info = simtwalk2)
tau2 <- 42
(tamrealmues1 <- round((nsim - burnin2) / tau2, 0))
cadsinrez2 <- simtwalk2$output[burnin2:nsim, ]
cadsinrez2 <- cadsinrez2[seq(1, dim(cadsinrez2)[1], tau2),]

# Contornos y puntos simulados
x11()
image.plot(seqbeta_,
           seqgamma_,
           post2 / K2,
           xlim = c(0.45, 0.63),
           ylim = c(0.1, 0.4))
points(cadsinrez2, pch = 16, cex = 0.6)
points(betreal, gamreal, pch = 16, col = "white", cex = 1)

x11()
par(mfrow = c(1, 2),mar = c(4, 4, 4, 2))
hist(cadsinrez2[, 1],
     xlab = expression(beta),
     freq = F,
     col = add.alpha("darksalmon", 0.4),
     border = "darksalmon", main = "")
abline(v = betreal, lwd = 2, col = "DeepPink4",lty = 2)
hist(cadsinrez2[, 2],
     xlab = expression(gamma),
     freq = F,
     col = add.alpha("darksalmon", 0.4), border = "darksalmon",main = "")
abline(v = gamreal,
       lwd = 2,
       col = "DeepPink4",
       lty = 2)
