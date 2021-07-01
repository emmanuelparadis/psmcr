## plots.R (2021-07-01)

## Plots Output of PSMC

## Copyright 2019-2021 Emmanuel Paradis

## This file is part of the R-package `psmcr'.
## See the file ../COPYING for licensing issues.

.getXYplot.psmc <- function(x, mutation.rate, g, scaled, bin.size)
{
    RS <- x$RS[x$RS[, "iter"] == x$niters, ]
    theta0 <- x$parameters[nrow(x$parameters), "theta0"]
    if (scaled) {
        ##xx <- RS[, "t_k"] / (theta0 / bin.size)
        xx <- RS[, "t_k"] / (theta0 * bin.size)
        yy <- theta0 * RS[, "lambda_k"] / bin.size
    } else {
        ##N0 <- theta0/(4 * mutation.rate / bin.size)
        denom <- 4 * mutation.rate * g * bin.size
        N0 <- theta0 / denom
        xx <- 2 * N0 * RS[, "t_k"]
        yy <- N0 * RS[, "lambda_k"]
    }
    list(xx = xx, yy = yy)
}

.getBootstrap.psmc <- function(x, mutation.rate, g, scaled, bin.size)
{
    boot <- x$bootstrap
    THETA0 <- rep(boot$theta0, each = x$n)
    if (scaled) {
        Tk.boot <- boot$tk / (THETA0 * bin.size)
        Nk.boot <- THETA0 * boot$lk / bin.size
    } else {
        denom <- 4 * mutation.rate * g * bin.size
        N0.boot <- boot$theta0 / denom
        Tk.boot <- 2 * N0.boot * boot$tk
        Nk.boot <- N0.boot * boot$lk
    }
    list(Tk.boot = Tk.boot, Nk.boot = Nk.boot)
}

plot.psmc <- function(x, type = "s", xlim = NULL, ylim = NULL, col = "grey",
                      xlab = if (scaled) "Scaled time" else "Time",
                      ylab = if (scaled) expression(Theta) else "N",
                      show.present = TRUE, mutation.rate = 1e-8, g = 1,
                      scaled = FALSE, bin.size = 100, ...)
{
    xy <- .getXYplot.psmc(x, mutation.rate = mutation.rate, g = g,
                          scaled = scaled, bin.size = bin.size)
    xx <- xy$xx
    yy <- xy$yy
    if (is.null(ylim)) yl <- max(yy)
    withbootstrap <- !is.null(x$bootstrap)
    if (withbootstrap) {
        obj <- .getBootstrap.psmc(x, mutation.rate = mutation.rate, g = g,
                                  scaled = scaled, bin.size = bin.size)
        Tk.boot <- obj$Tk.boot
        Nk.boot <- obj$Nk.boot
        if (is.null(ylim)) yl <- max(yl, Nk.boot)
    }
    if (is.null(xlim)) xlim <- range(xx)
    if (is.null(ylim)) ylim <- c(0, yl)
    plot(xx, yy, type = type, lwd = 3, col = col, xlab = xlab, ylab = ylab,
         xlim = xlim, ylim = ylim, ...)
    if (show.present) mtext("Present", 1, 2.5, at = 0, font = 3)
    if (withbootstrap) {
        col <- col2rgb(col)
        col <- rgb(col[1L], col[2L], col[3L], .3, maxColorValue = 255)
        for (j in 1:ncol(Tk.boot))
            points(Tk.boot[, j], Nk.boot[, j], type = type, col = col)
    }
}

lines.psmc <- function(x, type = "s", mutation.rate = 1e-8, g = 1,
                       col = "blue", scaled = FALSE, bin.size = 100, ...)
{
    xy <- .getXYplot.psmc(x, mutation.rate = mutation.rate, g = g,
                          scaled = scaled, bin.size = bin.size)
    xx <- xy$xx
    yy <- xy$yy
    boot <- x$bootstrap
    withbootstrap <- !is.null(boot)
    if (withbootstrap) {
        obj <- .getBootstrap.psmc(x, mutation.rate = mutation.rate, g = g,
                                  scaled = scaled, bin.size = bin.size)
        Tk.boot <- obj$Tk.boot
        Nk.boot <- obj$Nk.boot
    }
    lines(xx, yy, type = type, lwd = 3, col = col, ...)
    if (withbootstrap) {
        col <- col2rgb(col)
        col <- rgb(col[1L], col[2L], col[3L], .3, maxColorValue = 255)
        for (j in 1:ncol(Tk.boot))
            points(Tk.boot[, j], Nk.boot[, j], type = type, col = col)
    }
}
