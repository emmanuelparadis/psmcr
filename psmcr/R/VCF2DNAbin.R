## vcf2DNAbin.R (2019-11-15)

## Conversion

## Copyright 2019 Emmanuel Paradis

## This file is part of the R-package `psmcr'.
## See the file ../COPYING for licensing issues.

## different from pegas: all files are handled by readBin()
.VCFconnection <- function(file)
{
    file <- path.expand(file)
    remote <- if (length(grep("^(ht|f)tp(s|):", file))) TRUE else FALSE
    GZ <- if (length(grep("\\.gz$", file))) TRUE else FALSE
    if (GZ) {
        file <- if (remote) url(file) else gzfile(file)
        file <- gzcon(file)
    } else {
        stop("the VCF file must be compressed with GZ (*.vcf.gz)")
    }
    x <- readChar(file, 16L, TRUE)
    if (!identical(x, "##fileformat=VCF"))
        stop("file apparently not in VCF format")
    file
}

.VCFheader <- function(file)
{
    f <- .VCFconnection(file)
    x <- readBin(f, "raw", 1e9)
    CHROM <- charToRaw("#CHROM")
    i <- 1L
    while (!identical(x[i + 0:5], CHROM)) i <- i + 1L
    rawToChar(x[1:(i - 1)])
}

VCF2DNAbin <- function(file, refgenome = NULL, individual = 1, quiet = FALSE)
{
    chunck.size <- 1e9
    LF <- charToRaw("\n")
    hdr <- .VCFheader(file)
    sizeofheader <- nchar(hdr) + 1L
    hdr <- unlist(strsplit(hdr, "\n"))
    if (is.null(refgenome)) {
        refgenome <- grep("^##reference=", hdr, value = TRUE)
        if (nchar(refgenome) == 0)
            stop("no reference genome in VCF file")
        refgenome <- gsub("^##reference=", "", refgenome)
        refgenome <- read.FASTA(refgenome)
    } else {
        if (!inherits(refgenome, "DNAbin")) {
            if (!is.character(refgenome))
                stop("'refgenome' should be either the name of a FASTA file or a DNAbin object")
            refgenome <- read.FASTA(refgenome)
        }
    }
    f <- .VCFconnection(file)
    ## GZ <- if (inherits(f, "connection")) TRUE else FALSE
    ## the VCF file must be gzipped, so the above test is skipped
    open(f)
    trail <- raw()
    vol <- 0L
    res <- refgenome
    repeat {
        Y <- readBin(f, "raw", chunck.size)
        if (length(trail)) Y <- c(trail, Y)
        nY <- length(Y)
        if (!nY) break
        i <- nY
        vol <- vol + nY
        if (!quiet) cat("\rScanning", vol/1e6, "MB")
        while (Y[i] != LF) i <- i - 1L
        if (i < nY) {
            trail <- Y[(i + 1):nY]
            Y <- Y[1:i]
        } else trail <- raw()
        res <- .Call(vcf2DNAbinFromRaw, Y, res, sizeofheader,
                     as.integer(individual))
        sizeofheader <- 0L
    }
    close(f)
    if (!quiet) cat("\nDone.\n")
    res
}

seqBinning <- function(x, bin.size = 100)
{
    if (!inherits(x, "DNAbin"))
        stop("'x' must be of class DNAbin")
    res <- .Call("seqBinning_DNAbin", x, as.integer(bin.size))
    names(res) <- names(x)
    class(res) <- "DNAbin"
    res
}
