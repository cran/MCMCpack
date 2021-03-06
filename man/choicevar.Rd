% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCmnl.R
\name{choicevar}
\alias{choicevar}
\title{Handle Choice-Specific Covariates in Multinomial Choice Models}
\usage{
choicevar(var, varname, choicelevel)
}
\arguments{
\item{var}{The is the name of the variable in the dataframe.}

\item{varname}{The name of the new variable to be created.}

\item{choicelevel}{The level of \code{y} that the variable corresponds to.}
}
\value{
The new variable used by the \code{MCMCmnl()} function.
}
\description{
This function handles choice-specific covariates in multinomial choice
models.  See the example for an example of useage.
}
\seealso{
\code{\link{MCMCmnl}}
}
\keyword{manip}
