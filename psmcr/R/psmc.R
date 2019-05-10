## psmc.R (2019-05-06)

##   PSMC

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `psmcr'.
## See the file ../COPYING for licensing issues.

psmc <- function(x, parapattern = "4+5*3+4", maxt = 15, niters = 30,
                 trratio = 4, B = 0, trunksize = 5e5, decoding = FALSE,
                 quiet = FALSE, raw.output = FALSE, mc.cores = 1)
{
    if (decoding) {
        decoding <- FALSE
        warning("'decoding' is not yet available")
    }
    if (is.matrix(x)) x <- as.list(x)
    if (is.null(names(x))) stop("sequences have no labels")
    fl <- tempfile()
    flcpg <- paste(tempdir(), "CpG_counts", sep = "/")
    para <- list(parapattern, maxt, as.integer(niters),
                 trratio, as.integer(decoding),
                 as.integer(quiet), 0L, fl, flcpg)
    o <- .Call(Rpsmc_C, x, para) # o should be NULL
    out <- scan(fl, what = "", sep = "\n", quiet = TRUE)

    ## number of intervals:
    nintervs <- eval(parse(text = parapattern))

    ##if (B && length(x) == 1) {
    ##    B <- 0
    ##    warning("Only one chromosome: bootstrap not performed.")
    ##}

    if (B) {
        ## split the chromosomes:
        xboot <- list()
        trunksize <- as.integer(trunksize)
        if (any(trunksize > lengths(x)))
            warning("some sequences shorter than 'trunksize'")
        for (i in seq_along(x)) {
            z <- x[[i]]
            lz <- length(z)
            a <- 1; b <- trunksize
            repeat {
                if (b > lz) b <- lz
                if (a >= lz) break
                xboot <- c(xboot, list(z[a:b]))
                a <- b + 1L
                b <- b + trunksize
            }
        }
        llx <- ceiling(lengths(x)/trunksize)
        names(xboot) <- paste(unlist(mapply(rep, names(x), each = llx)),
                              unlist(mapply(":", 1, llx)),
                              sep = "_segment")

        theta0 <- numeric(B)
        tk <- lk <- matrix(0, nintervs, B)
        para[[6]] <- 1L
        para[[7]] <- 1L
        keep <- (nintervs * (niters - 1) + 1):(nintervs * niters)
        if (mc.cores == 1) {
            for (i in 1:B) {
                if (!quiet) cat("\rBootstrapping ", i, "/", B, "...", sep = "")
                unlink(fl)
                o <- .Call(Rpsmc_C, xboot, para)
                bout <- scan(fl, what = "", sep = "\n", quiet = TRUE)
                s <- grep("^RS", bout, value = TRUE)[keep]
                s <- gsub("^RS\t", "", s)
                s <- strsplit(s, "\t")
                s <- matrix(as.numeric(unlist(s)), ncol = 6, byrow = TRUE)
                tk[, i] <- s[, 2]
                lk[, i] <- s[, 3]
                ## get theta0
                theta <- grep("^TR", out, value = TRUE)
                theta <- gsub("^TR\t", "", theta)
                theta <- as.numeric(gsub("\t.*$", "", theta))
                theta0[i] <- theta[length(theta)]
                BOOT <- list(theta0 = theta0, tk = tk, lk = lk)
            }
        } else {
            foo <- function(i) {
                para[[8]] <- fl <- tempfile()
                o <- .Call(Rpsmc_C, xboot, para)
                bout <- scan(fl, what = "", sep = "\n", quiet = TRUE)
                s <- grep("^RS", bout, value = TRUE)[keep]
                s <- gsub("^RS\t", "", s)
                s <- strsplit(s, "\t")
                s <- matrix(as.numeric(unlist(s)), ncol = 6, byrow = TRUE)
                tk <- s[, 2]
                lk <- s[, 3]
                ## get theta0
                theta <- grep("^TR", out, value = TRUE)
                theta <- gsub("^TR\t", "", theta)
                theta <- as.numeric(gsub("\t.*$", "", theta))
                theta0 <- theta[length(theta)]
                unlink(fl)
                list(theta0, tk, lk)
            }
            if (!quiet) cat("Running parallel bootstraps...")
            BOOT <- mclapply(1:B, foo, mc.cores = mc.cores)
            theta0 <- sapply(BOOT, "[[", i = 1)
            tk <- sapply(BOOT, "[[", i = 2)
            lk <- sapply(BOOT, "[[", i = 3)
            BOOT <- list(theta0 = theta0, tk = tk, lk = lk)
        }
        if (!quiet) cat(" Done.\n")
    }

    ## return raw output if asked for it:
    if (raw.output) return(out)

    res <- list()
    res$niters <- niters
    res$n <- nintervs
    res$n_free_lambdas <- as.integer(gsub("^.*, n_free_lambdas:", "", out[13]))

    ## get the log-likelihoods:
    s <- grep("^LK", out, value = TRUE)
    res$logLik <- as.numeric(gsub("^LK\t", "", s))

    ## get Q differences:
    s <- grep("^QD", out, value = TRUE)
    s <- gsub("^QD\t", "", s)
    s <- unlist(strsplit(s, " -> "))
    s <- matrix(as.numeric(s), ncol = 2, byrow = TRUE)
    colnames(s) <- c("before", "after")
    res$EMQ <- s

    ## get RI (relative information or KL distance):
    s <- grep("^RI", out, value = TRUE)
    res$RI <- as.numeric(gsub("^RI\t", "", s))

    ## get C_pi and n_recomb
    s <- grep("^MM\tC_pi: ", out, value = TRUE)
    s <- gsub("^MM\tC_pi: ", "", s)
    s <- unlist(strsplit(s, ", n_recomb: "))
    s <- matrix(as.numeric(s), ncol = 2, byrow = TRUE)
    res$Cpi <- s[, 1]
    res$Nrecomb <- s[, 2]

    ## get theta, rho, and max T:
    s <- grep("^TR", out, value = TRUE)
    s <- gsub("^TR\t", "", s)
    s <- as.numeric(unlist(strsplit(s, "\t")))
    s <- matrix(s, ncol = 2, byrow = TRUE)
    ## append max T:
    maxt <- grep("^MT", out, value = TRUE)
    maxt <- as.numeric(gsub("^MT\t", "", maxt))
    s <- cbind(s, maxt)
    colnames(s) <- c("theta0", "rho0", "maxT")
    res$parameters <- s

    s <- grep("^RS", out, value = TRUE)
    s <- gsub("^RS\t", "", s)
    s <- strsplit(s, "\t")
    s <- matrix(as.numeric(unlist(s)), ncol = 6, byrow = TRUE)
    colnames(s) <- c("k", "t_k", "lambda_k", "pi_k", "sum_A_kl", "A_kk")
    iter <- rep(0:niters, each = nintervs)
    s <- cbind(s, iter = iter)
    res$RS <- s

    s <- grep("^TC", out, value = TRUE)
    s <- gsub("^TC\t", "", s)
    s <- strsplit(s, "\t")
    res$TC <- matrix(as.numeric(unlist(s)), ncol = 4, byrow = TRUE)

    if (decoding) {
        s <- grep("^DC", out, value = TRUE)
        s <- gsub("^DC\t", "", s)
        s <- strsplit(s, "\t")
        ns <- length(s)
        s <- as.data.frame(matrix(unlist(s), ns, 6, byrow = TRUE),
                           stringsAsFactors = FALSE)
        for (j in 2:6) s[[j]] <- as.numeric(s[[j]])
        names(s) <- c("Chromosome", "begin", "end", "best-k", "t_k+Delta_k", "max-prob")
        res$decoding <- s
    }
    if (B) res$bootstrap <- BOOT
    class(res) <- "psmc"
    res
}

print.psmc <- function(x, ...)
{
    cat("*** Pairwise sequential Markovian coalescent ***\n\n")
    cat("Number of time intervals: ", x$n, "\n", sep = "")
    cat("Number of free parameters: ", x$n_free_lambdas, "\n", sep = "")
}

logLik.psmc <- function(object, ...)
{
    object$logLik
}
