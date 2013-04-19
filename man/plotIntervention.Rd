\name{plotIntervention}
\alias{plotIntervention}
\title{Plot of intervention analysis}
\description{Plot the results of changepoint internvetion analysis}

\usage{
    plotIntervention(mcmcout, forward = TRUE, start = 1,
    alpha = 0.05, col = "red", main="", ylab="", xlab="")
}

\arguments{

\item{mcmcout}{The \code{mcmcout} object containing the posterior density sample from changepoint intervention analysis.}

\item{forward}{If TRUE, draw the result of forward prediction. If FALSE, 
draw the result of backward prediction.} 

\item{start}{The time of the first observation to be shown in the time series plot.}

\item{alpha}{The size of the prediction error band in the plot. By default, alpha = 0.05, which means the 95 percent Bayesian credible interval.}

\item{col}{The color of predicted data.}

\item{main}{Title of the plot} 

\item{ylab}{Label for the y-axis.}

\item{xlab}{Label for the x-axis.}

}

\references{
 Jong Hee Park. 2012. "A Change-point Approach to Intervention Analysis Using Bayesian Inference"
 Presented at the 2012 Annual Meeting of Korean Statistical Society.    
}

\seealso{\code{\link{MCMCintervention}}}
