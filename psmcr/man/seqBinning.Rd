\name{seqBinning}
\alias{seqBinning}
\title{Binning DNA Sequences}
\description{
  This function bins a set of DNA sequences in order to code SNPs for
  input to \code{\link{psmc}}.
}
\usage{
seqBinning(x, bin.size = 100)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"}.}
  \item{bin.size}{an integer giving the size (or width) of the bins.}
}
\details{
  The main argument \code{x} must be the consensus sequence of an
  individual, so that the ambiguous bases (R, Y, \dots) code for
  heterozygous sites. A consensus sequence can be produced with
  \code{\link{vcf2DNAbin}}.

  The sequences are split into bins of size \code{bin.size}: if there is
  at least one heterozygous site within a bin, then the output sequence
  has K, otherwise T. The sequence is to be used as input to
  \code{\link{psmc}}.
}
\value{
  an object of class \code{"DNAbin"}.
}
\author{Emmanuel Paradis}
\seealso{\code{\link{psmc}}, \code{\link{vcf2DNAbin}}}
\keyword{manip}
