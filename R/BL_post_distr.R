#'Computes the Black-Litterman posterior distribution.
#'
#'@description BL_post_distr computes posterior distribution in the Black-Litterman model starting from arbitrary prior distribution
#' given as a discrete time series \code{dat} and using \code{views_distr} -- submitted by the user distribution of views.
#'
#'@usage BL_post_distr (dat, returns_freq, prior_type = c("elliptic", NULL), market_portfolio,
#' SR, P, q, tau, risk = c("CVAR", "DCVAR", "LSAD", "MAD"),  alpha = NULL,
#' views_distr, views_cov_matrix_type = c("diag", "full"), cov_matrix = NULL)
#'
#'@param dat Time series of returns data; dat = cbind(rr, pk), where \eqn{rr} is an array (time series) of market asset returns,
#' for \eqn{n} returns and \eqn{k} assets it is an array with \eqn{\dim(rr) = (n, k)},
#' \eqn{pk} is a vector of length \eqn{n} containing probabilities of returns.
#'@param returns_freq Frequency of data in time series \code{dat}; given as a number of data rows corresponding to the period of 1 year,
#'i.e. 52 for weekly data or 12 for monthly data.
#'@param prior_type Type of distribution in time series \code{dat}; can be "elliptic" -- \eqn{rr} is distributed according
#'to (any) elliptical distribution, NULL -- \eqn{rr} is distributed according to any non-elliptical distribution.
#'@param market_portfolio Market portfolio -- benchmark (equilibrium) portfolio (for details see Palczewski&Palczewski).
#'@param SR Benchmark Sharpe ratio.
#'@param P "Pick" matrix in the Black-Litterman model (see Palczewski&Palczewski).
#'@param q Vector of investor's views on future returns in the Black-Litterman model (see Palczewski&Palczewski).
#'@param tau Confidence parameter in the Black-Litterman model.
#'@param risk Risk measure chosen for optimization; one of "CVAR", "DCVAR", "LSAD", "MAD", where
#' "CVAR" – denotes Conditional Value-at-Risk (CVaR),
#' "DCVAR" – denotes deviation CVaR,
#' "LSAD" – denotes Lower Semi Absolute Deviation,
#' "MAD" – denotes Mean Absolute Deviation.
#'@param alpha Value of alpha quantile in the definition of risk measures CVAR and DCVAR. Can be any number when risk measure is parameter free.
#'@param views_distr Distribution of views. An external function submitted by the user which computes densities of the distribution of views in given data points.
#' It is assumed implicitly that this distribution is an elliptical distribution (but any other distribution type can be used
#'provided calling to this function will preserve described below structure).
#' Call to that function has to be of the following form
#' \code{FUN(x,q,covmat,COF = NULL)}, where \code{x} is a data points matrix which collects in rows the coordinates of the points in which density is  computed,
#'\code{q} is a vector of investor's views,
#' \code{covmat} is covariance matrix of the distribution and \code{COF} is a vector of additional parameters characterizing the distribution (if needed).
#'@param views_cov_matrix_type Type of the covariance matrix of the distribution of views; can be:
#'"diag" --   diagonal part of the covariance matrix is used;
#'"full" -- the complete covariance matrix is used;
#'(for details see Palczewski&Palczewski).
#'@param cov_matrix Covariance matrix used for computation of market expected return (\code{RM}) from the formula
#' \code{RM =  SR * sqrt( t(w_m) * cov_matrix * w_m)} where \code{w_m} is market portfolio
#' and  \code{SR} -- benchmark Sharpe ratio.
#'When \code{cov_matrix} = NULL covariance matrix is computed   from matrix \eqn{rr} in data set \code{dat}.
#'@return
#' \tabular{llll}{
#'\code{post_distr}  \tab a time series of data for posterior distribution; for a time series of length \eqn{n} and \eqn{k} assets \cr
#'
#'\code{} \tab it is a matrix \eqn{(n, k+1)}, where columns (1:k) contain return vectors and the last column \cr
#'
#'\code{} \tab probabilities of returns.
#'}
#'
#'@examples
#'library(mvtnorm)
#'k = 3 
#'num =100
#'dat <-  cbind(rmvnorm (n=num, mean = rep(0,k), sigma=diag(k)), matrix(1/num,num,1)) 
#'# a data sample with num rows and (k+1) columns for k assets;
#'returns_freq = 52 # we assume that data frequency is 1 week
#'w_m <- rep(1/k,k) # benchmark portfolio, a vector of length k,
#'SR = 0.5 # Sharpe ratio
#'Pe <- diag(k) # we assume that views are "absolute views"
#'qe <- rep(0.05, k) # user's opinions on future returns (views)
#'tau = 0.02
#'BL_post_distr(dat, returns_freq, NULL, w_m, SR, Pe, qe, tau, risk = "MAD", alpha = 0,
#' views_distr = observ_normal, "diag", cov_matrix = NULL)
#'
#' 
#'
#'
#'@references Palczewski, J., Palczewski, A., Black-Litterman Model for Continuous Distributions (2016). Available at SSRN: https://ssrn.com/abstract=2744621.  
#'
#'@export

BL_post_distr <- function(dat, returns_freq, prior_type = c("elliptic", NULL),  market_portfolio, SR, P, q, tau, risk = c("CVAR", "DCVAR", "LSAD", "MAD"),
alpha = NULL, views_distr, views_cov_matrix_type = c("diag", "full"), cov_matrix = NULL)
{
 FUN = match.fun(views_distr, descend = FALSE)
 PARAM = views_cov_matrix_type
 risk <- match.arg(risk)
 prior_type <- match.arg(prior_type)
 w_m = market_portfolio/sum(market_portfolio)
 if (!is.null(cov_matrix)){
   sample_cov = cov_matrix
 }
 else {
   disc_data = .discrete_variance (returns_freq, dat)
   sample_cov = disc_data$variance
 }
 nVar <- length(w_m)
 if( !all.equal( dim( sample_cov ), c( nVar, nVar ) ) == TRUE ) {
   stop( paste( "Number of asset for which return data are provided  ",
                "must be the same as length of portfolio weights.\n" ) )
 }
 if (is.null(prior_type)) {
   RM = SR * sqrt( t(w_m) %*% sample_cov %*% w_m)
   eq_mu = equilibrium_mean (dat, w_m, RM/returns_freq, risk, alpha)}
 else {
   Agamma = SR / sqrt( t(w_m) %*% sample_cov %*% w_m)
   eq_mu = .equilibrium_mean_elliptic (sample_cov/returns_freq, w_m, Agamma[1,1])
 }

 tmp = .post_distr_new (dat, t(eq_mu$market_returns), q/returns_freq, P, sample_cov/returns_freq, tau, FUN, PARAM)

 return(list(post_distr =tmp))
}
