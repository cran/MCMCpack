\name{PErisk}
\alias{PErisk}
\docType{data}
\title{Political Economic Risk Data from 62 Countries in 1987}
\description{
Political Economic Risk Data from 62 Countries in 1987.
}
\usage{data(CountryRisk)}
\format{
  A data frame with 62 observations on the following 9 variables. All
  data points are from 1987. See Quinn (2004) for more details. 
  \describe{
    \item{country}{a factor with levels \code{Argentina} \code{Australia} \code{Austria} \code{Bangladesh} \code{Belgium} \code{Bolivia} \code{Botswana} \code{Brazil} \code{Burma} \code{Cameroon} \code{Canada} \code{Chile} \code{Colombia} \code{Congo-Kinshasa} \code{Costa Rica} \code{Cote d'Ivoire} \code{Denmark} \code{Dominican Republic} \code{Ecuador} \code{Finland} \code{Gambia, The} \code{Ghana} \code{Greece} \code{Hungary} \code{India} \code{Indonesia} \code{Iran} \code{Ireland} \code{Israel} \code{Italy} \code{Japan} \code{Kenya} \code{Korea, South} \code{Malawi} \code{Malaysia} \code{Mexico} \code{Morocco} \code{New Zealand} \code{Nigeria} \code{Norway} \code{Papua New Guinea} \code{Paraguay} \code{Philippines} \code{Poland} \code{Portugal} \code{Sierra Leone} \code{Singapore} \code{South Africa} \code{Spain} \code{Sri Lanka} \code{Sweden} \code{Switzerland} \code{Syria} \code{Thailand} \code{Togo} \code{Tunisia} \code{Turkey} \code{United Kingdom} \code{Uruguay} \code{Venezuela} \code{Zambia} \code{Zimbabwe}}
    \item{courts}{an ordered factor with levels \code{0} <
      \code{1}.\code{courts} is an indicator of whether the country in
      question is judged to have an independent judiciary. From Henisz
      (2002).}
    \item{barb2}{a numeric vector giving the natural log of the black
      market premium in each country. The black market premium is coded
      as the black market exchange rate (local currency per dollar)
      divided by the official exchange  rate minus 1. From Marshall,
      Gurr, and Harff (2002). }
    \item{prsexp2}{an ordered factor with levels \code{0} < \code{1} <
      \code{2} < \code{3} < \code{4} < \code{5}, giving the lack of
      expropriation risk. From Marshall, Gurr, and Harff (2002).}
    \item{prscorr2}{an ordered factor with levels \code{0} < \code{1} <
      \code{2} < \code{3} < \code{4} < \code{5}, measuring the lack of
      corruption. From Marshall, Gurr, and Harff (2002).}
    \item{gdpw2}{a numeric vector giving the natural log of real GDP per
      worker in 1985 international prices. From Alvarez et al. (1999).}
  }
}
\source{
  Mike Alvarez, Jose Antonio Cheibub, Fernando Limongi, and Adam
  Przeworski. 1999. ``ACLP Political and Economic Database.''
  \url{http://www.ssc.upenn.edu/~cheibub/data/}.

  Witold J. Henisz. 2002. ``The Political Constraint Index (POLCON)
  Dataset.''
  \url{http://www-management.wharton.upenn.edu/henisz/POLCON/ContactInfo.html}.

  Monty G. Marshall, Ted Robert Gurr, and Barbara Harff. 2002. ``State
  Failure Task Force Problem Set.''
  \url{http://www.cidcm.umd.edu/inscr/stfail/index.htm}. 
}
\references{
  Kevin M. Quinn. 2004. ``Bayesian Factor Analysis for Mixed Ordinal
  and Continuous Response.'' \emph{Political Analyis}. 12: 338-353.
}
\keyword{datasets}
