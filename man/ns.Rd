\name{ns}
\alias{ns}
\title{Marginal likelihood estimation using Nested Sampling (NS)}
\description{
\code{ns} computes the marginal likelihood for a fixed tree topology using nested sampling and also provides posterior samples.
}
\usage{
ns(data, tree, N = 5, k = 4, a = 3, b = 0.2, al = 1, bl = 1,
   step = 50, max.iter = 2000, end = 2, act_plot = TRUE,
   model)
}
\arguments{
  \item{data}{An alignment, object of class \code{phyDat}.}
  \item{tree}{A phylogenetic \code{tree}, object of class phylo.}
  \item{N}{Number of active points used by the NS algorithm.}
  \item{k}{Number of intervals of the discrete gamma distribution.}
  \item{a}{Shape parameter for the prior over the branch lengths, see details.}
  \item{b}{Rate parameter for the prior over the branch lengths, see details.}
  \item{al}{Shape parameter for the prior over the gamma shape parameter, see details.}
  \item{bl}{Scale parameterfor the prior over the gamma shape parameter, see details.}
  \item{step}{Thining factor used in Metropolis algorithm.}
  \item{max.iter}{Maximum number of iterations.}
  \item{end}{Factor in the termination criterion, see details.}
  \item{act_plot}{Allows to produce relevant plots, including the tree with the posterior means as branch lengths.}
  \item{model}{Allows to choose a nucleotide model, see details.}
  \item{}{}
}
\details{
Specific priors have been assigned over the parameters, but they offer some degree of flexibility.  For the branch lengths, the following hierarchical structure is considered:

 \deqn{t_{i} | \mu ~ Exp(1/\mu), for i= 1,...,n}

 \deqn{\mu ~ Inverse-Gamma(a, b)}

with
\deqn{f(t_{i} | \mu) = \frac{1}{\mu} e^{ -t_{i} / \mu}}
\deqn{f(\mu) = \frac{b^a}{\Gamma(a)} \mu^{-a-1} e^{- b / \mu}.}

Depending on the chosen model (see below), the following priors are considered:  for the relative rates \eqn{q_{j} | \phi ~ Exp(1/\phi)} --(\eqn{j=1} and \eqn{j=5} for the HKY and GTR models, respectively)-- with \eqn{\phi ~ Exp(1)}; a \eqn{Dirichlet(1,1,1,1)} is used for the joint prior of the four nucleotide frequencies; and \eqn{Gamma(al, bl)} for the gamma shape parameter, density function given by
\deqn{f(\lambda) = \frac{1}{\Gamma(al) bl^al} \lambda^{al-1} e^{-\lambda / bl}.}

Via \code{model} the following models can be specified: \code{JC} (Jukes-Cantor), \code{JC_G} (Jukes-Cantor + Gamma), \code{HKY} (Hasegawa-Kishino), \code{HKY_G} (Hasegawa-Kishino + Gamma), \code{GTR} (general time reversible), and \code{GTR_G} (general time reversible + Gamma).

Skilling(2006) proposed "\code{end} * \eqn{H * N} << \code{iter}" as a termination criterion.  In practice, \code{end}\eqn{ = 2} is commonly used.

}
\value{
\code{ns} returns a list with the following elements:
\item{logZ}{Log-marginal likelihood estimate.}
\item{sd_logZ}{Standar deviation estimate of \code{logZ}.}
\item{H}{}
\item{iter}{Number of iterations.}
\item{time}{Elapsed time in seconds.}
\item{sampled_par}{Posterior samples drawn from NS.}
The next elements can be used in further analysis:
\item{logX}{log-cummulative prior mass}
\item{logLd}{}
\item{Lw.z}{}
\item{seq_lz}{Evolution of \code{logZ} through the iterations.}
\item{dTheta}{Discarded points.}

}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\references{

Maturana R., P., Brewer, B.J., and Klaere, S., 2017. Model selection and parameter inference in phylogenetics using nested sampling. ArXiv preprint arXiv:1703.05471v1.

Skilling, J., 2006. Nested sampling for general Bayesian computation. \emph{Bayesian
Analysis} \bold{1}(4), 833--860.

}
\examples{
data(rbcl)

library(ape) # read.tree function
tree <- read.tree(text = "(((((((Avena_sativa_L15300,Iris_unguicularis_AJ309693),Nicotiana_tabacum_J01450),Picea_pungens_AF456382  ),(Asplenium_nidus_U05907,Osmunda_cinnamomea_D14882)),Sphagnum_palustre_L13485),Bazzania_trilobata_L11056),Conocephalum_conicum_U87067,Chara_connivens_L13476);")

R = ns(data, tree, N = 1, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
       step = 50, max.iter = 100, end = Inf, act_plot = TRUE,
       model = "GTR_G")

R = ns(data, tree, N = 2, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
       step = 50, max.iter = 2000, end = 2, act_plot = TRUE,
       model = "GTR_G")
}
