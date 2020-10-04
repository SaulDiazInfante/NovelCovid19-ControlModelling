divergence_analysis <- function(fit = nuts_fit){
    #
    #
    c_light <- c("#DCBCBC")
    c_light_highlight <- c("#C79999")
    c_mid <- c("#B97C7C")
    c_mid_highlight <- c("#A25050")
    c_dark <- c("#8F2727")
    c_dark_highlight <- c("#7C0000")
    divergent <- get_sampler_params(nuts_fit,
                                    inc_warmup=FALSE)[[1]][,'divergent__']
    sum(divergent)
    sum(divergent) / 20000
    params_cp <- as.data.frame(extract(nuts_fit, permuted=FALSE))
    names(params_cp) <- gsub("chain:1.", "", names(params_cp), fixed = TRUE)
    names(params_cp) <- gsub("[", ".", names(params_cp), fixed = TRUE)
    names(params_cp) <- gsub("]", "", names(params_cp), fixed = TRUE)
    params_cp$iter <- 1:20000
    par(mar = c(4, 4, 0.5, 0.5))
    plot(params_cp$iter,
        log(params_cp[["chain:5.theta.1"]]),
        col=c_dark,
        pch=16,
        cex=0.8,
        xlab="Iteration",
        ylab="log(tau)", ylim=c(-6, 4))
    params_cp$divergent <- divergent
    div_params_cp <- params_cp[params_cp$divergent == 1,]
    nondiv_params_cp <- params_cp[params_cp$divergent == 0,]
    par(mar = c(4, 4, 0.5, 0.5))
    plot(nondiv_params_cp[["theta.1"]],
            log(nondiv_params_cp[["theta.3"]]),
            col=c_dark,
            pch=16,
            cex=0.8,
            xlab="theta.2",
            ylab="log(beta_a)",
            xlim=c(0, 20),
            ylim=c(-10, 5))
    points(div_params_cp[["theta.1"]],
            log(div_params_cp[["theta.3"]]),
            col = rgb(0, 1, 0, 0.2),
            pch=14,
            cex=0.8)
}
