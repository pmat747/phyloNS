\name{nis}
\alias{nis}
\title{Nested Importance Sampling (NIS) analysis - Fixed tree topology case.}
\description{
\code{nis} computes the marginal likelihood for a fixed tree topology using nested importance sampling. It requires posterior samples to calibrate the reference distributions.
}
\usage{
nis(data, tree, N = 1, k = 4, a = 3, b = 0.2, al = NULL, bl = NULL,
    step = 50, max.iter = 2000, tol = 1e-07, model, post_sample = NULL)
}
\arguments{
  \item{data}{An alignment, object of class \code{phyDat}.}
  \item{tree}{A phylogenetic \code{tree}, object of class \code{phylo}.}
  \item{N}{Number of active points used by the NS algorithm.}
  \item{k}{Number of intervals of the discrete gamma distribution.}
  \item{a}{Shape parameter for the prior over the branch lengths, see \code{ns} function for details.}
  \item{b}{Rate parameter for the prior over the branch lengths, see \code{ns} function for details.}
  \item{al}{Shape parameter for the prior over the gamma shape parameter, see \code{ns} function for details.}
  \item{bl}{Scale parameterfor the prior over the gamma shape parameter, see \code{ns} function for details.}
  \item{step}{Thining factor used in Metropolis algorithm.}
  \item{max.iter}{Maximum number of iterations.}
  \item{tol}{Tolerance value in the relative error.}
  \item{model}{Allows to choose a nucleotide model, see \code{ns} function for details.}
  \item{post_sample}{Posterior samples.}

}
\details{
The reference distribution is defined by using the prior calibrated by posterior samples, method proposed by Fan et al. (2011).  The marginal posterior sample means (\eqn{\hat{\mu}}) and variances (\eqn{\hat{\sigma}^2}) are matched to the respective prior parameters in order to calibrate the reference distributions (see Fan et al. (2011) for details).

}
\value{
\code{nis} returns a list with the following elements:
\item{logZ}{Log-marginal likelihood estimate.}
\item{H}{Amount of information.}
\item{logPsLd}{Sequence of pseudo-likelihood values}
\item{logPsPd}{Sequence of reference distribution values}
\item{logLd}{Sequence of likelihood values.}
\item{dTheta}{Discarded points.}
\item{seq_lz}{Evolution of \code{logZ} through the iterations.}
\item{time}{Elapsed time in seconds.}
\item{info}{Specifications of the process.}
}
\author{
Patricio Maturana Russel \email{p.russel@auckland.ac.nz}
}
\references{

Fan, Y., Wu, R., Chen, M.-H., Kuo, L., Lewis, P. O., 2011. Choosing among partition models in Bayesian phylogenetics. \emph{Mol. Biol. Evol.} \bold{28}(1), 523--532.

Maturana R., P., Brewer, B.J., and Klaere, S., 2017. Model selection and parameter inference in phylogenetics using nested sampling. ArXiv preprint arXiv:1703.05471v1.

Skilling, J., 2006. Nested sampling for general Bayesian computation. \emph{Bayesian
Analysis} \bold{1}(4), 833--860.

}
\examples{
data(rbcl)

library(ape) # read.tree function
tree <- read.tree(text = "(((((((Avena_sativa_L15300,Iris_unguicularis_AJ309693),Nicotiana_tabacum_J01450),Picea_pungens_AF456382  ),(Asplenium_nidus_U05907,Osmunda_cinnamomea_D14882)),Sphagnum_palustre_L13485),Bazzania_trilobata_L11056),Conocephalum_conicum_U87067,Chara_connivens_L13476);")

# Only to get samples from the posterior
R = ns(data, tree, N = 10, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
       step = 100, max.iter = 800, end = Inf, act_plot = TRUE,
       model = "GTR_G")
# NIS
R1 = nis(data, tree, N = 1, k = 4, a = 3, b = 0.2, al = 10, bl = 0.026,
         step = 50, max.iter = 100, tol = 1e-10,
         model = "GTR_G", post_sample = R$sampled_par)
}
