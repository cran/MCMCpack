\name{choicevar}
\alias{choicevar}

\title{Handle Choice-Specific Covariates in Multinomial Choice Models}
\description{
This function handles choice-specific covariates in multinomial choice models.  
See the example for an example of useage.
}
\usage{
choicevar(var, varname, choicelevel)
}

\arguments{
  \item{var}{The is the name of the variable in the dataframe.}
  \item{varname}{The name of the new variable to be created.}
  \item{choicelevel}{The level of \code{y} that the variable
  corresponds to.}
}

\value{
  The new variable used by the \code{MCMCmnl()} function.
}

\seealso{\code{\link{MCMCmnl}}}
\keyword{manip}
