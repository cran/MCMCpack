#' Stochastic search variable selection for quantile regression
#'
#' This function uses stochastic search to select promising regression models
#' at a fixed quantile \eqn{\tau}.  Indicator variables
#' \eqn{\gamma} are used to represent whether a predictor is included in
#' the model or not.  The user supplies the data and the prior distribution on
#' the model size.  A list is returned containing the posterior sample of
#' \eqn{\gamma} and the associated regression parameters
#' \eqn{\beta}.
#'
#'
#' \code{SSVSquantreg} implements stochastic search variable selection
#' over a set of potential predictors to obtain promising models.  The
#' models considered take the following form:
#'
#' \deqn{Q_{\tau}(y_i|x_{i\gamma}) = x_{i\gamma} ' \beta_{\gamma},}
#'
#' where \eqn{Q_{\tau}(y_i|x_{i\gamma})} denotes the conditional
#' \eqn{\tau}th quantile of \eqn{y_i} given \eqn{x_{i\gamma}},
#' \eqn{x_{i\gamma}} denotes \eqn{x_i} with those predictors
#' \eqn{x_{ij}} for which \eqn{\gamma_j=0} removed and
#' \eqn{\beta_{\gamma}} denotes the model specific regression
#' parameters.
#'
#' The likelihood is formed based on the assumption of independent
#' asymmetric Laplace distributions on the \eqn{y_i} with skewness
#' parameter \eqn{\tau} and location parameters \eqn{ x_{i\gamma} '
#' \beta_{\gamma}}. This assumption ensures that the likelihood
#' function is maximised by the \eqn{\tau}th conditional quantile of
#' the response variable.
#'
#' The prior on each \eqn{\beta_j} is
#'
#' \deqn{(1-\gamma_j)\delta_0+\gamma_j\mbox{Cauchy}(0,1),}
#'
#' where \eqn{\delta_0} denotes a degenerate distribution with all
#' mass at 0.  A standard Cauchy distribution is chosen conditional on
#' \eqn{\gamma_j=1}.  This allows for a wider range of nonzero
#' values of \eqn{\beta_j} than a standard Normal distribution,
#' improving the robustness of the method.  Each of the indicator variables
#' \eqn{\gamma_j} is independently assigned a Bernoulli prior, with
#' prior probability of inclusion \eqn{\pi_0}.  This in turn is assigned
#' a beta distribution, resulting in a beta-binomial prior on the model size.
#' The user can supply the hyperparameters for the beta distribution.  Starting
#' values are randomly generated from the prior distribution.
#'
#' It is recommended to standardise any non-binary predictors in order to
#' compare these predictors on the same scale.  This can be achieved using the
#' \code{scale} function.
#'
#' If it is certain that a predictor should be included, all predictors
#' specified are brought to the first positions for computational convenience.
#' The regression parameters associated with these predictors are given
#' independent improper priors. Users may notice a small speed advantage if
#' they specify the predictors that they feel certain should appear in the
#' model, particularly for large models with a large number of observations.
#'
#' @param formula Model formula.
#'
#' @param data Data frame.
#'
#' @param tau The quantile of interest. Must be between 0 and 1. The default
#' value of 0.5 corresponds to median regression model selection.
#'
#' @param include The predictor(s) that should definitely appear in the model.
#' Can be specified by name, or their position in the formula (taking into
#' account the intercept).
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of MCMC iterations after burnin.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' MCMC iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen.  If \code{verbose} is greater than 0 the
#' iteration number, the most recently sampled values of \eqn{\gamma}
#' and the associated values of \eqn{\beta} are printed to the screen
#' every \code{verbose}th iteration.
#'
#' @param seed The seed for the random number generator. If NA, the Mersenne
#' Twister generator is used with default seed 12345; if an integer is passed
#' it is used to seed the Mersenne twister. The default value for this argument
#' is a random integer between 1 and 1,000,000.  This default value ensures
#' that if the function is used again with a different value of
#' \eqn{\tau}, it is extremely unlikely that the seed will be identical.
#' The user can also pass a list of length two to use the L'Ecuyer random
#' number generator, which is suitable for parallel computation.  The first
#' element of the list is the L'Ecuyer seed, which is a vector of length six or
#' NA (if NA a default seed of \code{rep(12345,6)} is used).  The second
#' element of list is a positive substream number. See the MCMCpack
#' specification for more details.
#'
#' @param pi0a0,pi0b0 Hyperparameters of the beta prior on \eqn{\pi_0},
#' the prior probability of including a predictor. Default values of (1,1) are
#' equivalent to a uniform distribution.
#'
#' @param ... Further arguments
#'
#' @return A list containing:
#'
#' \item{gamma}{The posterior sample of \eqn{\gamma}. This has
#' associated summary and plot methods.}
#'
#' \item{beta}{The posterior sample of the associated regression parameters
#' \eqn{\beta}. This can be analysed with functions from the coda
#' package.}
#'
#' @export
#'
#' @author Craig Reed
#'
#' @seealso \code{\link[MCMCpack]{MCMCquantreg}},
#' \code{\link[MCMCpack]{summary.qrssvs}}, \code{\link[MCMCpack]{plot.qrssvs}},
#' \code{\link[MCMCpack]{mptable}}, \code{\link[MCMCpack]{topmodels}},
#' \code{\link[base]{scale}}, \code{\link[quantreg]{rq}}
#'
#' @references Craig Reed, David B. Dunson and Keming Yu. 2010. "Bayesian
#' Variable Selection for Quantile Regression" Technical Report.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.2.} \url{http://scythe.lsa.umich.edu}.
#'
#' Keming Yu and Jin Zhang. 2005. "A Three Parameter Asymmetric Laplace
#' Distribution and it's extensions." \emph{Communications in Statistics -
#' Theory and Methods}, 34, 1867-1879.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' @keywords models
#'
#' @examples
#'
#' \dontrun{
#'
#' set.seed(1)
#' epsilon<-rnorm(100)
#' set.seed(2)
#' x<-matrix(rnorm(1000),100,10)
#' y<-x[,1]+x[,10]+epsilon
#' qrssvs<-SSVSquantreg(y~x)
#' model.50pc<-SSVSquantreg(y~x)
#' model.90pc<-SSVSquantreg(y~x,tau=0.9)
#' summary(model.50pc) ## Intercept not in median probability model
#' summary(model.90pc) ## Intercept appears in median probability model
#' }
#'
"SSVSquantreg" <-
  function(formula, data=NULL, tau=0.5, include=NULL, burnin=1000, mcmc = 10000, thin=1,
           verbose = 0, seed = sample(1:1000000,1), pi0a0=1, pi0b0=1,
           ...) {

    ## checks
    check.offset(list(...))
    if (pi0a0<=0 || pi0b0<=0){
	stop("Parameters pi0a0 and pi0b0 must be positive. \nPlease respecify and call again.\n")
	}
    cl <- match.call()
    if (tau<=0 || tau>=1){
	stop("tau must be in (0,1). \nPlease respecify and call again.\n")
	}


    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    ## form response and model matrices
    holder <- parse.formula(formula, data=data)
    Y <- holder[[1]]
    X <- holder[[2]]
    xnames <- holder[[3]]
    K <- ncol(X)  # number of covariates
    q <- length(include) #number of covariates that are pre-specified to appear in the model

    if (!is.null(include)){
	if(is.character(include)){
	  if(!all(include%in%xnames)){
		include.positions<-NA
		}
          else{
		include.positions<-match(include, xnames)
		}
	  }
	else{
	  if(max(include)>length(xnames)){
		include.positions<-NA
		}
	  else{
	  	include.positions<-include
	  	}
	  }
	if (any(is.na(include.positions))){
	  stop("One or more covariates to be included are not present in the design matrix\n or one or more indices are out of range. Please respecify and call again.\n")
	  }

    ## Bring those covariates that are pre-specified to appear in the model to the first positions in the X matrix
	X <- cbind(X[,include.positions], X[,-include.positions])
	xnames <- c(xnames[include.positions], xnames[-include.positions])
	}

    ## define holder for posterior sample
    sample <- matrix(data=0, mcmc/thin, 2*K)
    posterior <- NULL

    ## call C++ code to draw sample
    auto.Scythe.call(output.object="posterior", cc.fun.name="cSSVSquantreg",
                     sample.nonconst=sample, tau=as.double(tau), Y=Y, X=X, q=as.integer(q),
                     burnin=as.integer(burnin), mcmc=as.integer(mcmc), thin=as.integer(thin),
                     lecuyer=as.integer(lecuyer),
                     seedarray=as.integer(seed.array),
                     lecuyerstream=as.integer(lecuyer.stream),
                     verbose=as.integer(verbose),
                     pi0a0 = as.double(pi0a0), pi0b0=as.double(pi0b0),
		     package="MCMCpack")


output <- form.mcmc.object(posterior,                                                names=rep(xnames, times=2),
                               title="SSVSquantreg Posterior Sample",
                               y=Y, call=cl
                               )

gammasample<-output[,1:K]
attr(gammasample, "tau")<-tau
attr(gammasample, "xnames")<-xnames
attr(gammasample, "class")<-"qrssvs"
betasample<-output[,-(1:K)]
attr(betasample,"tau")<-tau
attr(gammasample, "call") <- attr(betasample, "call") <- cl
return(list(gamma=gammasample,beta=betasample))
}
