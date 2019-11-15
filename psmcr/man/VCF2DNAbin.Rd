\name{VCF2DNAbin}
\alias{VCF2DNAbin}
\title{Conversion}
\description{
  This function outputs the consensus genome of an individual from its
  diploid genotypes stored in a VCF file and a reference genome.
}
\usage{
VCF2DNAbin(file, refgenome = NULL, individual = 1, quiet = FALSE)
}
\arguments{
  \item{file}{a character string giving the VCF file.}
  \item{refgenome}{the reference genome; see details.}
  \item{individual}{an integer value giving the individual to be
    analysed.}
  \item{quiet}{a logical value specifying whether to print the progress
    of the calculations.}
}
\details{
  The genotypes in the VCF file \emph{must all be} diploid. See the
  package \pkg{pegas} for various functions to handle and check VCF
  files.

  The reference genome can be given in different ways. By default, a
  file name (or URL) is looked for in the VCF file. Alternatively,
  \code{refgenome} can be either a \code{"DNAbin"} object or the name of
  a FASTA file.
}
\value{
  an object of class \code{"DNAbin"}.
}
\author{Emmanuel Paradis}
\keyword{manip}
