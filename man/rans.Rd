\name{rans}
\alias{rans}
\title{Re-analysis of Nested Sampling (NS) output(s).}
\description{
\code{rans} allows a re-analysis of one or more NS outputs for the fixed and variable tree topology cases.  In case of multiple NS outputs, this information is combined as one NS run.  The NS analysis is redone "\code{S}" times using the stochastic approach in order to estimate the prior mass (see Maturana R. et al. (2017) for more details).
}
\usage{
rans(R, S = 100, n_ns_outputs = 1, act_plot = TRUE)
}
\arguments{
\item{R}{NS outputs obtained from \code{ns} or \code{vns} functions, see details.}
\item{S}{Number of the re-analyses.}
\item{n_ns_outputs}{Number of NS runs.}
\item{act_plot}{Plots the marginal likelihood distribution.}
}
\details{

The NS outputs obtained from \code{ns} or \code{vns} functions are used in \code{rans} as inputs.  For instance, in case of one NS output, that is \code{R = vns(...)}, the re-analysis can be performed as \code{rans(R)}.  In case of two NS outputs, that is \code{R1 = vns(...)} and \code{R2 = vns(...)}, these outputs must be stored in a list as \code{R = list(R1 = R1, R2 = R2)} in order to perform the re-analysis as \code{rans(R, n_ns_outputs = 2)}.

\code{ns} and \code{vns} use the deterministic approach in order to estimate the cumulative prior probabilities.  \code{rans} uses the stochastic approach using the Beta distribution for the cumulative prior mass estimation.  Thus, at each re-analysis (out of \code{S}) different estimates are obtained allowing to estimate the uncertainty of the estimates.

}
\value{
\item{sim_mu_est}{95\% confident intervals for the means for the continuous parameters.}
\item{sim_sd_est}{95\% confident intervals for the standard deviations for the continuous parameters.}
\item{sim_med_est}{95\% confident intervals for the medians for the continuous parameters.}
\item{vecZ_est}{95\% confident intervals for the marginal likelihood.}
\item{vecH_est}{95\% confident intervals for the information value.}
\item{vecZ}{\code{S} marginal likelihood estimates.}
\item{vecH}{\code{S} information value estimates.}
\item{sim_mu}{Means from the NS posterior samples obtained at each re-analysis.  }
\item{sim_sd}{Standard deviations from the NS posterior samples obtained at each re-analysis.}
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\references{

Maturana R., P., Brewer, B.J., Klaere, S., and Bouckaert, R. 2017. Model selection and parameter inference in phylogenetics using nested sampling. ArXiv preprint arXiv:1703.05471v2.

Skilling, J., 2006. Nested sampling for general Bayesian computation. \emph{Bayesian
Analysis} \bold{1}(4), 833--860.

}
\examples{
data(rbcl)

R1 = vns(data, N = 1, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
         step = 500, max.iter = 130, end = Inf, act_plot = TRUE,
         model = "GTR_G")

R2 = vns(data, N = 1, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
         step = 500, max.iter = 130, end = Inf, act_plot = TRUE,
         model = "GTR_G")

# Re-analysis of a NS output
rans(R1, S = 100)

# Re-analysis of 2 NS outputs
rans(list(R1 = R1, R2 = R2), S = 100, n_ns_outputs = 2)

}
