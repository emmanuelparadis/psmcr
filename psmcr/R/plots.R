## plots.R (2019-05-06)

## Plots Output of PSMC

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `psmcr'.
## See the file ../COPYING for licensing issues.

plot.psmc <- function(x, type = "s",
                      xlab = if (scaled) "Scaled time" else "Time",
                      ylab = if (scaled) expression(Theta) else "N",
                      show.present = TRUE, mutation.rate = 1e-8,
                      scaled = FALSE, bin.size = 1, ...)
{
    RS <- x$RS[x$RS[, "iter"] == x$niters, ]
    theta0 <- x$parameters[nrow(x$parameters), "theta0"]
    if (scaled) {
        xx <- RS[, "t_k"] / (theta0 / bin.size)
        yy <- theta0 * RS[, "lambda_k"] / bin.size
    } else {
        N0 <- theta0/(4 * mutation.rate / bin.size)
        xx <- 2 * N0 * RS[, "t_k"]
        yy <- N0 * RS[, "lambda_k"]
    }
    boot <- x$bootstrap
    withbootstrap <- !is.null(boot)
    yl <- max(yy)
    if (withbootstrap) {
        THETA0 <- rep(boot$theta0, each = x$n)
        if (scaled) {
            Tk.boot <- boot$tk / THETA0
            Nk.boot <- boot$lk * THETA0
        } else {
            N0.boot <- boot$theta0/(4 * mutation.rate)
            Tk.boot <- 2 * N0.boot * boot$tk
            Nk.boot <- N0.boot * boot$lk
        }
        yl <- max(yl, Nk.boot)
    }
    plot(xx, yy, type = type, lwd = 3, xlab = xlab, ylab = ylab,
         ylim = c(0, yl), ...)
    if (show.present) mtext("Present", 1, 2.5, at = 0, font = 3)
    if (withbootstrap)
        for (j in 1:ncol(Tk.boot))
            points(Tk.boot[, j], Nk.boot[, j], type = type, col = rgb(.5, .5, .5, .3))
}

