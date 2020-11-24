\name{plot.psmc}
\alias{plot.psmc}
\title{Plotting Results of PSMC}
\description{
  Graphical display of the results of PSMC fitting.
}
\usage{
\method{plot}{psmc}(x, type = "s", xlab = if (scaled) "Scaled time" else "Time",
     ylab = if (scaled) expression(Theta) else "N",
     show.present = TRUE, mutation.rate = 1e-8,
     scaled = FALSE, bin.size = 100, \dots)
}
\arguments{
  \item{x}{an object of class \code{"psmc"}.}
  \item{type}{a character string giving the type of plot.}
  \item{xlab}{the label under the horizontal axis.}
  \item{ylab}{the label beside the vertical axis.}
  \item{show.present}{a logical value: whether to indicate the present
    in italics (near the zero coordinate on the x-axis).}
  \item{mutation.rate}{the mutation rate for calculating population
    sizes and times (ignored if \code{scaled = TRUE}; see details).}
  \item{scaled}{a logical value. By default, the absolute times and
    population sizes are calculated assuming a value of mutation rate
    given by the previous argument. If \code{scaled = TRUE}, these
    values are scaled (see details).}
  \item{bin.size}{an integer giving the value of the bin size used in
    \code{\link{seqBinning}()}; the default value is the same than in
    this function.}
  \item{\dots}{further arguments passed to \code{plot}.}
}
\details{
  The PSMC output is scaled to \eqn{2N_0}{2N0}. There are two ways
  of rescaling the time and the population size more
  meaningfully. Firstly, suppose we know the per-site, per-generation
  mutation rate \eqn{\mu}{mu}, we can compute \eqn{N_0} as:

  \deqn{N_0 = \theta_0 / (4\mu s)}{N_0 = theta_0 / (4 mu s)}

  where \eqn{\theta_0}{theta0} is the estimated value from
  \code{\link{psmc}}, and \eqn{s} is the bin size we use for generating
  the PSMC input. Knowing \eqn{N_0}, we can scale time to generations
  and relative population size to effective size by:

  \deqn{T_k = 2N_0 t_k}
  \deqn{N_k = N_0 \lambda_k}

  where \eqn{t_k} and \eqn{\lambda_k}{lambda_k} are output by
  \code{\link{psmc}}.

  A problem with the above strategy is that we do not know a definite
  answer of \eqn{\mu}{mu} and in fact it varies with regions and
  mutation types. An alternative way is to use per-site pairwise
  sequence divergence to represent time:

  \deqn{d_k = 2\mu T_k = t_k \theta_0  s}{d_k = 2 mu T_k = t_k theta_0  s}

  and use scaled mutation rate to represent population size:

  \deqn{\theta_k = 4N_k \mu = \lambda_k \theta_0  s}{theta_k = 4N_k mu = lambda_k * theta_0 * s}

  This second method is used if \code{scaled = TRUE}.
}
\value{NULL.}
\references{
  Li, H. and Durbin, R. (2011) Inference of human population history
  from individual whole-genome sequences. \emph{Nature}, \bold{475},
  493--U84.
}
\author{Heng Li and Emmanuel Paradis}
\keyword{hplot}
