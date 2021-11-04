## import.psmc.R (2021-11-04)

##   Import PSMC Files

## Copyright 2021 Emmanuel Paradis

## This file is part of the R-package `psmcr'.
## See the file ../COPYING for licensing issues.

.extract_psmc <- function(out, res, decoding)
{
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
    iter <- rep(0:res$niters, each = res$n)
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
    res
}

import.psmc <- function(file)
{
    out <- scan(file, what = "", sep = "\n", quiet = TRUE)

    niters <- as.integer(gsub("^.*n_iterations:|,.*$", "", out[14]))
    parapattern <- gsub("^.*pattern:|,.*$", "", out[13])
    nintervs <- eval(parse(text = parapattern))
    decoding <- as.integer(gsub("^.*is_decoding:", "", out[15]))

    res <- list()
    res$niters <- niters
    res$n <- nintervs
    res <- .extract_psmc(out, res, decoding)
    class(res) <- "psmc"
    res
}
