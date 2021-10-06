## sample from the posterior distribution of a one-dimensional item
## response theory model in R using linked C++ code in Scythe.
##
## ADM and KQ 1/23/2003
## updated extensively ADM & KQ 7/28/2004
## store.ability arg added KQ 1/27/2006
## Hierarchical Subject Parameters MJM 2007-11-07
## Parameter Expansion 2008-11-18
## Grouped Subject Parameters (Wdata) started MJM,YS 2008-11

#' Markov Chain Monte Carlo for Hierarchical One Dimensional Item Response
#' Theory Model, Covariates Predicting Latent Ideal Point (Ability)
#'
#' This function generates a sample from the posterior distribution of a one
#' dimensional item response theory (IRT) model, with multivariate Normal
#' priors on the item parameters, and a Normal-Inverse Gamma hierarchical prior
#' on subject ideal points (abilities).  The user supplies item-response data,
#' subject covariates, and priors. Note that this identification strategy
#' obviates the constraints used on theta in \code{\link[MCMCpack]{MCMCirt1d}}.
#' A sample from the posterior distribution is returned as an mcmc object,
#' which can be subsequently analyzed with functions provided in the coda
#' package.
#'
#' If you are interested in fitting K-dimensional item response theory models,
#' or would rather identify the model by placing constraints on the item
#' parameters, please see \code{\link[MCMCpack]{MCMCirtKd}}.
#'
#' \code{MCMCirtHier1d} simulates from the posterior distribution using
#' standard Gibbs sampling using data augmentation (a Normal draw for the
#' subject abilities, a multivariate Normal draw for (second-level) subject
#' ability predictors, an Inverse-Gamma draw for the (second-level) variance of
#' subject abilities, a multivariate Normal draw for the item parameters, and a
#' truncated Normal draw for the latent utilities). The simulation proper is
#' done in compiled C++ code to maximize efficiency.  Please consult the coda
#' documentation for a comprehensive list of functions that can be used to
#' analyze the posterior sample.
#'
#' The model takes the following form.  We assume that each subject has an
#' subject ability (ideal point) denoted \eqn{\theta_j} and that each
#' item has a difficulty parameter \eqn{a_i} and discrimination parameter
#' \eqn{b_i}.  The observed choice by subject \eqn{j} on item
#' \eqn{i} is the observed data matrix which is \eqn{(I \times J)}.
#' We assume that the choice is dictated by an unobserved utility:
#'
#' \deqn{z_{i,j} = -\alpha_i + \beta_i \theta_j + \varepsilon_{i,j}}
#'
#' Where the errors are assumed to be distributed standard Normal.
#' This constitutes the measurement or level-1 model. The subject abilities
#' (ideal points) are modeled by a second level Normal linear predictor for
#' subject covariates \code{Xjdata}, with common variance
#' \eqn{\sigma^2}. The parameters of interest are the subject
#' abilities (ideal points), item parameters, and second-level coefficients.
#'
#' We assume the following priors.  For the subject abilities (ideal points):
#'
#' \deqn{\theta_j \sim \mathcal{N}(\mu_{\theta} ,T_{0}^{-1})}
#'
#' For the item parameters, the prior is:
#'
#' \deqn{\left[a_i, b_i \right]' \sim \mathcal{N}_2 (ab_{0},AB_{0}^{-1})}
#'
#' The model is identified by the proper priors on the item parameters and
#' constraints placed on the ability parameters.
#'
#' As is the case with all measurement models, make sure that you have plenty
#' of free memory, especially when storing the item parameters.
#'
#' @param datamatrix The matrix of data.  Must be 0, 1, or missing values.  The
#' rows of \code{datamatrix} correspond to subjects and the columns correspond
#' to items.
#'
#' @param Xjdata A \code{data.frame} containing second-level predictor
#' covariates for ideal points \eqn{\theta}. Predictors are modeled as a
#' linear regression on the mean vector of \eqn{\theta}; the posterior
#' sample contains regression coefficients \eqn{\beta} and common
#' variance \eqn{\sigma^2}. See Rivers (2003) for a thorough
#' discussion of identification of IRT models.
#'
#' @param burnin The number of burn-in iterations for the sampler.
#'
#' @param mcmc The number of Gibbs iterations for the sampler.
#'
#' @param thin The thinning interval used in the simulation.  The number of
#' Gibbs iterations must be divisible by this value.
#'
#' @param verbose A switch which determines whether or not the progress of the
#' sampler is printed to the screen. If \code{verbose} is greater than 0 then
#' every \code{verbose}th iteration will be printed to the screen.
#'
#' @param seed The seed for the random number generator.  If NA, the Mersenne
#' Twister generator is used with default seed 12345; if an integer is passed
#' it is used to seed the Mersenne twister.  The user can also pass a list of
#' length two to use the L'Ecuyer random number generator, which is suitable
#' for parallel computation.  The first element of the list is the L'Ecuyer
#' seed, which is a vector of length six or NA (if NA a default seed of
#' \code{rep(12345,6)} is used).  The second element of list is a positive
#' substream number. See the MCMCpack specification for more details.
#'
#' @param theta.start The starting values for the subject abilities (ideal
#' points). This can either be a scalar or a column vector with dimension equal
#' to the number of voters.  If this takes a scalar value, then that value will
#' serve as the starting value for all of the thetas.  The default value of NA
#' will choose the starting values based on an eigenvalue-eigenvector
#' decomposition of the agreement score matrix formed from the
#' \code{datamatrix}.
#'
#' @param a.start The starting values for the \eqn{a} difficulty parameters.
#' This can either be a scalar or a column vector with dimension equal to the
#' number of items.  If this takes a scalar value, then that value will serve
#' as the starting value for all \eqn{a}.  The default value of NA will set
#' the starting values based on a series of probit regressions that condition
#' on the starting values of theta.
#'
#' @param b.start The starting values for the \eqn{b} discrimination
#' parameters. This can either be a scalar or a column vector with dimension
#' equal to the number of items.  If this takes a scalar value, then that value
#' will serve as the starting value for all \eqn{b}.  The default value of
#' NA will set the starting values based on a series of probit regressions that
#' condition on the starting values of theta.
#'
#' @param beta.start The starting values for the \eqn{\beta} regression
#' coefficients that predict the means of ideal points \eqn{\theta}.
#' This can either be a scalar or a column vector with length equal to the
#' number of covariates. If this takes a scalar value, then that value will
#' serve as the starting value for all of the betas.  The default value of NA
#' will set the starting values based on a linear regression of the covariates
#' on (either provided or generated) \code{theta.start}.
#'
#' @param b0 The prior mean of \eqn{\beta}. Can be either a scalar or a
#' vector of length equal to the number of subject covariates. If a scalar all
#' means with be set to the passed value.
#'
#' @param B0 The prior precision of \eqn{\beta}.  This can either be a
#' scalar or a square matrix with dimensions equal to the number of betas.  If
#' this takes a scalar value, then that value times an identity matrix serves
#' as the prior precision of beta. A default proper but diffuse value of .01
#' ensures finite marginal likelihood for model comparison.  A value of 0 is
#' equivalent to an improper uniform prior for beta.
#'
#' @param c0 \eqn{c_0/2} is the shape parameter for the inverse Gamma
#'   prior on \eqn{\sigma^2} (the variance of \eqn{\theta}). The
#'   amount of information in the inverse Gamma prior is something
#'   like that from \eqn{c_0} pseudo-observations.
#'
#' @param d0 \eqn{d_0/2} is the scale parameter for the inverse Gamma
#' prior on \eqn{\sigma^2} (the variance of \eqn{\theta}). In
#' constructing the inverse Gamma prior, \eqn{d_0} acts like the sum of
#' squared errors from the \eqn{c_0} pseudo-observations.
#'
#' @param ab0 The prior mean of \code{(a, b)}. Can be either a scalar or a
#' 2-vector. If a scalar both means will be set to the passed value. The prior
#' mean is assumed to be the same across all items.
#'
#' @param AB0 The prior precision of \code{(a, b)}.This can either be ascalar
#' or a 2 by 2 matrix. If this takes a scalar value, then that value times an
#' identity matrix serves as the prior precision. The prior precision is
#' assumed to be the same across all items.
#'
#' @param store.item A switch that determines whether or not to store the item
#' parameters for posterior analysis.  \emph{NOTE: In situations with many
#' items storing the item parameters takes an enormous amount of memory, so
#' \code{store.item} should only be \code{TRUE} if the chain is thinned
#' heavily, or for applications with a small number of items}.  By default, the
#' item parameters are not stored.
#'
#' @param store.ability A switch that determines whether or not to store the
#' ability parameters for posterior analysis.  \emph{NOTE: In situations with
#' many individuals storing the ability parameters takes an enormous amount of
#' memory, so \code{store.ability} should only be \code{TRUE} if the chain is
#' thinned heavily, or for applications with a small number of individuals}.
#' By default, ability parameters are stored.
#'
#' @param drop.constant.items A switch that determines whether or not items
#' that have no variation should be deleted before fitting the model. Default =
#' TRUE.
#'
#' @param marginal.likelihood Should the marginal likelihood of the
#' second-level model on ideal points be calculated using the method of Chib
#' (1995)? It is stored as an attribute of the posterior \code{mcmc} object and
#' suitable for comparison using \code{\link[MCMCpack]{BayesFactor}}.
#'
#' @param px Use Parameter Expansion to reduce autocorrelation in the chain?
#' PX introduces an unidentified parameter \eqn{alpha} for the residual
#' variance in the latent data (Liu and Wu 1999). Default = TRUE
#'
#' @param px_a0 Prior shape parameter for the inverse-gamma distribution on
#' \eqn{alpha}, the residual variance of the latent data. Default=10.
#'
#' @param px_b0 Prior scale parameter for the inverse-gamma distribution on
#' \eqn{alpha}, the residual variance of the latent data. Default = 10
#'
#' @param ... further arguments to be passed
#'
#' @return An \code{mcmc} object that contains the sample from the posterior
#' distribution. This object can be summarized by functions provided by the
#' coda package. If \code{marginal.likelihood = "Chib95"} the object will have
#' attribute \code{logmarglike}.
#'
#' @export
#'
#' @author Michael Malecki, \email{mike@@crunch.io},
#' \url{https://github.com/malecki}.
#'
#' @seealso \code{\link[coda]{plot.mcmc}},\code{\link[coda]{summary.mcmc}},
#' \code{\link[MCMCpack]{MCMCirtKd}}
#'
#' @references James H. Albert. 1992. ``Bayesian Estimation of Normal Ogive
#' Item Response Curves Using Gibbs Sampling." \emph{Journal of Educational
#' Statistics}.  17: 251--269.
#'
#' Joshua Clinton, Simon Jackman, and Douglas Rivers. 2004. ``The Statistical
#' Analysis of Roll Call Data."  \emph{American Political Science Review} 98:
#' 355--370.
#'
#' Valen E. Johnson and James H. Albert. 1999. ``Ordinal Data Modeling."
#' Springer: New York.
#'
#' Liu, Jun S. and Ying Nian Wu. 1999. ``Parameter Expansion for Data
#' Augmentation.'' \emph{Journal of the American Statistical Association} 94:
#' 1264--1274.
#'
#' Andrew D. Martin, Kevin M. Quinn, and Jong Hee Park. 2011.  ``MCMCpack:
#' Markov Chain Monte Carlo in R.'', \emph{Journal of Statistical Software}.
#' 42(9): 1-21.  \doi{10.18637/jss.v042.i09}.
#'
#' Daniel Pemstein, Kevin M. Quinn, and Andrew D. Martin.  2007.  \emph{Scythe
#' Statistical Library 1.0.} \url{http://scythe.lsa.umich.edu}.
#'
#' Martyn Plummer, Nicky Best, Kate Cowles, and Karen Vines. 2006.  ``Output
#' Analysis and Diagnostics for MCMC (CODA)'', \emph{R News}. 6(1): 7-11.
#' \url{https://CRAN.R-project.org/doc/Rnews/Rnews_2006-1.pdf}.
#'
#' Douglas Rivers.  2004.  ``Identification of Multidimensional Item-Response
#' Models."  Stanford University, typescript.
#'
#' @keywords models
#'
#' @examples
#'
#'    \dontrun{
#' data(SupremeCourt)
#'
#' Xjdata <- data.frame(presparty= c(1,1,0,1,1,1,1,0,0),
#'                      sex= c(0,0,1,0,0,0,0,1,0))
#'
#' ## Parameter Expansion reduces autocorrelation.
#'   posterior1 <- MCMCirtHier1d(t(SupremeCourt),
#'                    burnin=50000, mcmc=10000, thin=20,
#'                    verbose=10000,
#'                    Xjdata=Xjdata,
#'                    marginal.likelihood="Chib95",
#' 		   px=TRUE)
#'
#' ## But, you can always turn it off.
#'   posterior2 <- MCMCirtHier1d(t(SupremeCourt),
#'                    burnin=50000, mcmc=10000, thin=20,
#'                    verbose=10000,
#'                    Xjdata=Xjdata,
#'                    #marginal.likelihood="Chib95",
#' 		   px=FALSE)
#' ## Note that the hierarchical model has greater autocorrelation than
#' ## the naive IRT model.
#'   posterior0 <- MCMCirt1d(t(SupremeCourt),
#'                         theta.constraints=list(Scalia="+", Ginsburg="-"),
#'                         B0.alpha=.2, B0.beta=.2,
#'                         burnin=50000, mcmc=100000, thin=100, verbose=10000,
#'                         store.item=FALSE)
#'
#' ## Randomly 10% Missing -- this affects the expansion parameter, increasing
#' ## the variance of the (unidentified) latent parameter alpha.
#'
#'    scMiss <- SupremeCourt
#'    scMiss[matrix(as.logical(rbinom(nrow(SupremeCourt)*ncol(SupremeCourt), 1, .1)),
#'       dim(SupremeCourt))] <- NA
#'
#'    posterior1.miss <- MCMCirtHier1d(t(scMiss),
#'                    burnin=80000, mcmc=10000, thin=20,
#'                    verbose=10000,
#'                    Xjdata=Xjdata,
#'                    marginal.likelihood="Chib95",
#' 		   px=TRUE)
#'
#'    }
#'
"MCMCirtHier1d" <-
  function(datamatrix, Xjdata,
           burnin = 1000, mcmc = 20000, thin=1,
           verbose = 0, seed = NA,
           theta.start = NA, a.start = NA, b.start = NA,
           beta.start=NA, b0=0, B0=.01, c0=.001, d0=.001,
           ab0=0, AB0=.25, store.item = FALSE, store.ability=TRUE,
           drop.constant.items=TRUE,
           marginal.likelihood=c("none","Chib95"), px=TRUE,
           px_a0 = 10, px_b0=10,
           ... ) {

    ## checks
    check.offset(list(...))
    check.mcmc.parameters(burnin, mcmc, thin)

    ## check vote matrix and convert to work with C++ code
    if (drop.constant.items==TRUE){
      x.col.var <- apply(datamatrix, 2, var, na.rm=TRUE)
      keep.inds <- x.col.var>0
      keep.inds[is.na(keep.inds)] <- FALSE
      datamatrix <- datamatrix[,keep.inds]
    }
    datamatrix <- as.matrix(datamatrix)
    K <- ncol(datamatrix)   # cases, bills, items, etc
    J <- nrow(datamatrix)   # justices, legislators, subjects, etc
    L <- ncol(Xjdata)       # predictors on theta Xj

    if(sum(datamatrix==1 | datamatrix==0 | is.na(datamatrix)) != (J*K)) {
      cat("Error: Data matrix contains elements other than 0, 1 or NA.\n")
      stop("Please check data and call ", calling.function(), " again.\n",
              call.=FALSE)
    }
    datamatrix[is.na(datamatrix)] <- 9
    item.names <- colnames(as.data.frame(datamatrix))
    subject.names <- rownames(as.data.frame(datamatrix))
    beta.names <- c(names(as.data.frame(Xjdata)),"sigmasq")

    ## names
    item.names <- colnames(datamatrix)
    if (is.null(item.names)){
      item.names <- paste("item", 1:K, sep="")
    }

    ## check Xj matrix and set up betastart
    Xjdata <- as.matrix(Xjdata)
    if(nrow(Xjdata) != nrow(datamatrix)) {
      cat("Error: subject covariates not of same length as datamatrix\n")
      stop("Please check data and try ",calling.function()," again.\n",call.=FALSE)
    }


    ## prior for (alpha, beta)
    holder <- form.mvn.prior(ab0, AB0, 2)
    ab0 <- holder[[1]]
    AB0 <- holder[[2]]

    ## starting values for theta error checking
    ## could use factor.score.start.check EXCEPT
    ## We have done away with eq and ineq constraints.
    if (max(is.na(theta.start))==1) {
      theta.start <- factor.score.eigen.start(agree.mat(datamatrix), 1)
    }
    else if(is.numeric(theta.start) & length(theta.start) == J ) {
      theta.start <- theta.start * matrix(1, J, 1)
    }
    else {
      cat("Inappropriate value of theta.start passed.\n")
      stop("Please respecify and call", calling.function(), " again.\n",
           call.=FALSE)
    }
    ## starting values for (a, b)
    ab.starts <- matrix(NA, K, 2)
    for (i in 1:K){
      local.y <-  datamatrix[,i]
      local.y[local.y==9] <- NA
      if (var(na.omit(local.y))==0){
        ab.starts[i,] <- c(0,10)
      }
      else {
        ab.starts[i,] <- coef(suppressWarnings(glm(local.y~theta.start,
                                                   family=binomial(link="probit"),
                                                   control=glm.control(
                                                     maxit=8, epsilon=1e-3)
                                                   )))
      }
    }
    ab.starts[,1] <- -1 * ab.starts[,1] # make this into a difficulty param

    ## starting values for a and b error checking
    if (is.na(a.start)) {
      a.start <- ab.starts[,1]
    }
    else if(is.null(dim(a.start))) {
      a.start <- a.start * matrix(1,K,1)
    }
    else if((dim(a.start)[1] != K) || (dim(a.start)[2] != 1)) {
      cat("Error: Starting value for a not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
           call.=FALSE)
    }
    if (is.na(b.start)) {
      b.start <- ab.starts[,2]
    }
    else if(is.null(dim(b.start))) {
      b.start <- b.start * matrix(1,K,1)
    }
    else if((dim(b.start)[1] != K) || (dim(b.start)[2] != 1)) {
      cat("Error: Starting value for b not conformable.\n")
       stop("Please respecify and call ", calling.function(), " again.\n",
              call.=FALSE)
    }

    cat("Generating starting values (glm) for hierarchical parameters:\n")
    ## starting values are regression of theta.start on Xj
    ## or passed vector, or same value all beta
    if (max(is.na(beta.start))==1) {         # beta.start NA
      beta.start <- coef(suppressWarnings(glm.fit(Xjdata,theta.start)))
    } else if ( length(beta.start) == L ) {  # beta.start vector
      beta.start <- matrix(beta.start,L,1)
    } else if ( length(beta.start) == 1 ) {  # beta.start scalar
      beta.start <- beta.start * matrix(1,L,1)
    } else {
      cat("Error: Starting value for beta not conformable.\n")
      stop("Please respecify and call ", calling.function(), " again.\n",
           call.=FALSE)
    }
    print(beta.start)
    ## prior for beta
    holder <- form.mvn.prior(b0, B0, L)
    b0 <- holder[[1]]
    B0 <- holder[[2]]
    check.ig.prior(c0, d0)

        ## get marginal likelihood argument
    marginal.likelihood  <- match.arg(marginal.likelihood)
    B0.eigenvalues <- eigen(B0)$values
    if (min(B0.eigenvalues) < 0){
      stop("B0 is not positive semi-definite.\nPlease respecify and call again.\n")
    }
    if (isTRUE(all.equal(min(B0.eigenvalues), 0))){
      if (marginal.likelihood != "none"){
        warning("Cannot calculate marginal likelihood with improper prior\n")
        marginal.likelihood <- "none"
      }
    }
    logmarglike <- NULL
    chib <- 0
    if (marginal.likelihood == "Chib95"){
      chib <- 1
    }

    cat("setting up posterior holder\n" )

    ## define holder for posterior sample
    if(store.item == FALSE & store.ability == TRUE) {
      sample <- matrix(data=0, mcmc/thin, J+L+1)
    }
    else if (store.item == TRUE & store.ability == FALSE){
      sample <- matrix(data=0, mcmc/thin, L+1 + 2*K)
    }
    else if (store.item == TRUE & store.ability == TRUE){
      sample <- matrix(data=0, mcmc/thin, L+1 + J + 2 * K)
    }
    else{
      stop("Either store.item or store.ability should be true.\n")
    }

    ## seeds
    seeds <- form.seeds(seed)
    lecuyer <- seeds[[1]]
    seed.array <- seeds[[2]]
    lecuyer.stream <- seeds[[3]]

    # call C++ code to draw sample
    posterior <- .C("cMCMCirtHier1d",
                    sampledata = as.double(sample),
                    samplerow = as.integer(nrow(sample)),
                    samplecol = as.integer(ncol(sample)),
                    Xdata = as.integer(datamatrix),
                    Xrow = as.integer(nrow(datamatrix)),
                    Xcol = as.integer(ncol(datamatrix)),
                    burnin = as.integer(burnin),
                    mcmc = as.integer(mcmc),
                    thin = as.integer(thin),
                    lecuyer = as.integer(lecuyer),
                    seedarray = as.integer(seed.array),
                    lecuyerstream = as.integer(lecuyer.stream),
                    verbose = as.integer(verbose),
                    thetastartdata = as.double(theta.start),
                    thetastartrow = as.integer(nrow(theta.start)),
                    thetastartcol = as.integer(ncol(theta.start)),
                    astartdata = as.double(a.start),
                    astartrow = as.integer(length(a.start)),
                    astartcol = as.integer(1),
                    bstartdata = as.double(b.start),
                    bstartrow = as.integer(length(b.start)),
                    bstartcol = as.integer(1),
                    ab0data = as.double(ab0),
                    ab0row = as.integer(nrow(ab0)),
                    ab0col = as.integer(ncol(ab0)),
                    AB0data = as.double(AB0),
                    AB0row = as.integer(nrow(AB0)),
                    AB0col = as.integer(ncol(AB0)),
                    Xjdata = as.double(Xjdata),
                    Xjrow = as.integer(nrow(Xjdata)),
                    Xjcol = as.integer(ncol(Xjdata)),
                    betastartdata = as.double(beta.start),
                    betastartrow = as.integer(length(beta.start)),
                    betastartcol = as.integer(1),
                    b0data = as.double(b0),
                    b0row = as.integer(length(b0)),
                    b0col = as.integer(1),
                    B0data = as.double(B0),
                    B0row = as.integer(nrow(B0)),
                    B0col = as.integer(ncol(B0)),
                    c0 = as.double(c0),
                    d0 = as.double(d0),
                    storei = as.integer(store.item),
                    storea = as.integer(store.ability),
                    logmarglikeholder = as.double(0.0),
                    chib = as.integer(chib),
                    px= as.integer(px),
                    px_a0 = as.double(px_a0),
                    px_b0 = as.double(px_b0),
                    PACKAGE="MCMCpack"
                  )

    beta.names <- paste("beta.",beta.names,sep="")
    theta.names <- paste("theta.", subject.names, sep = "")
    alpha.beta.names <- paste(rep(c("a.","b."), K),
                              rep(item.names, each = 2),
                              sep = "")

    # put together matrix and build MCMC object to return
    sample <- matrix(posterior$sampledata, posterior$samplerow,
                     posterior$samplecol,
                     byrow=FALSE)

    output <- mcmc(data=sample, start=burnin+1, end=burnin+mcmc, thin=thin)
    if (marginal.likelihood == "Chib95"){
      logmarglike <- posterior$logmarglikeholder
    }

    names <- NULL
    if(store.ability == TRUE) {
      names <- c(names, theta.names)
    }
    if (store.item == TRUE){
      names <- c(names, alpha.beta.names)
    }
    names <- c(names,beta.names)

    try( varnames(output) <- names)
    attr(output,"title") <-
      "MCMCirtHier1d Posterior Sample"
    attr(output,"logmarglike") <- posterior$logmarglikeholder
    return(output)

  }
