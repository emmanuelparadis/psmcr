\name{psmc}
\alias{psmc}
\alias{print.psmc}
\alias{logLik.psmc}
\title{PSMC}
\description{
  This function fits the pairwise sequential Markovian coalescent (PSMC)
  model to a diploid genome sequence.
}
\usage{
psmc(x, parapattern = "4+5*3+4", maxt = 15, niters = 30, trratio = 4,
     B = 0, trunksize = 5e5, decoding = FALSE, quiet = FALSE,
     raw.output = FALSE, mc.cores = 1)
\method{print}{psmc}(x, \dots)
\method{logLik}{psmc}(object, \dots)
}
\arguments{
  \item{x}{an object of class \code{"DNAbin"} or of class \code{"psmc"}.}
  \item{parapattern}{a character string giving the parameterization (see
    details).}
  \item{maxt}{a numeric value giving the largest possible value for time
    to the most recent common ancestor (MRCA).}
  \item{niters}{an integer giving the number of iterations.}
  \item{trratio}{a numeric value giving the ratio
    \eqn{\theta/\rho}{theta/rho}.}
  \item{B}{the number of bootstrap replicates. By default, no bootstrap
    is done.}
  \item{trunksize}{the length of the segments generated under the
    boostrap (see details).}
  \item{decoding}{a logical value; whether to perform decoding.}
  \item{quiet}{a logical value specifying whether to print the progress
    of the calculations.}
  \item{raw.output}{a logical value; if \code{TRUE} the function returns
    the raw output produced by the C code.}
  \item{mc.cores}{the number of cores (CPUs) to be used for
    bootstrapping.}
  \item{object}{an object of class \code{"psmc"}.}
  \item{\dots}{arguments passed between methods.}
}
\details{
  This is a port to R of a program written in C by Heng Li.

  The default value of \code{parapatt} is only indicative and must be
  adapated. The present value codes for 23 time atomic intervals and 7
  free parameters: the first parameter spans the first 4 atomic time
  intervals, each of the next 5 parameters spans 3 intervals, and the
  7th one spans 4 intervals.

  The options \code{parapatt} and \code{maxt} are manually chosen such
  that after 20 rounds of iterations, at least ~10 recombinations are
  inferred to occur in the intervals each parameter spans. Impropriate
  settings may lead to overfitting. For example, \code{parapatt = "4+25*2+4+6", maxt = 15} have been shown to be suitable for modern humans.

  \code{psmc} infers the scaled mutation rate, the recombination rate and
  the free population size parameters. All these parameters are scaled
  to \eqn{2N_0}{2N0}.

  To perform bootstrapping, long chromosome sequences are first split
  into shorter segments, \code{psmc} will then randomly sample with
  replacement from these segments.
}
\value{
  a list with the class \code{"psmc"}, or a vector of mode character if
  \code{raw.output = TRUE}.
}
\references{
  Li, H. and Durbin, R. (2011) Inference of human population history
  from individual whole-genome sequences. \emph{Nature}, \bold{475},
  493--U84.
}
\author{Original C code by Heng Li; port to R by Emmanuel Paradis.}
\seealso{\code{\link{plot.psmc}}}
\keyword{model}
