#' Dutch Voting Behavior in 1989
#'
#' Dutch Voting Behavior in 1989.
#'
#' @name Nethvote
#'
#' @docType data
#'
#' @format A data frame with 1754 observations and 11 variables from the 1989
#' Dutch Parliamentary Election Study (Anker and Oppenhuis, 1993). Each
#' observation is a survey respondent.  These data are a subset of one of five
#' multiply imputed datasets used in Quinn and Martin (2002). For more
#' information see Quinn and Martin (2002).
#' \describe{
#'   \item{vote}{A factor giving the self-reported vote choice of each respondent.
#'     The levels are CDA (Christen Democratisch Appel), D66 (Democraten 66), Pvda
#'     (Partij van de Arbeid), and VVD (Volkspartij voor Vrijheid en Democratie).}
#'   \item{distD66}{A numeric variable giving the squared ideological distance
#'     between the respondent and the D66. Larger values indicate ideological
#'     dissimilarity between the respondent and the party.}
#'   \item{distPvdA}{A numeric variable giving the squared ideological distance
#'     between the respondent and the PvdA. Larger values indicate ideological
#'     dissimilarity between the respondent and the party.}
#'   \item{distVVD}{A numeric variable giving the squared ideological distance between
#'     the respondent and the VVD. Larger values indicate ideological dissimilarity
#'     between the respondent and the party.}
#'   \item{distCDA}{A numeric variable giving the squared ideological
#'     distance between the respondent and the CDA. Larger values indicate
#'     ideological dissimilarity between the respondent and the party.}
#'   \item{relig}{An indicator variable equal to 0 if the respondent is not
#'     religious and 1 if the respondent is religious.}
#'   \item{class}{Social class of respondent. 0 is the lowest social class,
#'     4 is the highest social class.}
#'   \item{income}{Income of respondent. 0 is lowest and 6 is highest.}
#'   \item{educ}{Education of respondent. 0 is lowest and 4 is highest.}
#'   \item{age}{Age category of respondent. 0 is lowest and 12 is highest.}
#'   \item{urban}{Indicator variable equal to 0 if the respondent is not a
#'     resident of an urban area and 1 if the respondent is a resident of an urban
#'     area.}
#' }
#'
#' @references Kevin M. Quinn and Andrew D. Martin. 2002. ``An Integrated
#' Computational Model of Multiparty Electoral Competition.'' \emph{Statistical
#' Science}. 17: 405-419.
#'
#' @source H. Anker and E.V. Oppenhuis. 1993. ``Dutch Parliamentary Election
#' Study.'' (computer file). Dutch Electoral Research Foundation and
#' Netherlands Central Bureau of Statistics, Amsterdam.
#'
#' @keywords datasets
NULL

#' Political Economic Risk Data from 62 Countries in 1987
#'
#' Political Economic Risk Data from 62 Countries in 1987.
#'
#' @name PErisk
#'
#' @docType data
#'
#' @format A data frame with 62 observations on the following 9 variables. All
#' data points are from 1987. See Quinn (2004) for more details.
#' \describe{
#'   \item{country}{a factor with levels \code{Argentina} through
#'     \code{Zimbabwe}} \item{courts}{an ordered factor with levels \code{0} <
#'     \code{1}.\code{courts} is an indicator of whether the country in question is
#'     judged to have an independent judiciary. From Henisz (2002).}
#'   \item{barb2}{a numeric vector giving the natural log of the black market
#'     premium in each country. The black market premium is coded as the black market
#'     exchange rate (local currency per dollar) divided by the official exchange rate
#'     minus 1. From Marshall, Gurr, and Harff (2002). }
#'   \item{prsexp2}{an ordered factor
#'     with levels \code{0} < \code{1} < \code{2} < \code{3} < \code{4} < \code{5},
#'     giving the lack of expropriation risk. From Marshall, Gurr, and Harff
#'     (2002).}
#'   \item{prscorr2}{an ordered factor with levels \code{0} < \code{1} <
#'     \code{2} < \code{3} < \code{4} < \code{5}, measuring the lack of corruption.
#'     From Marshall, Gurr, and Harff (2002).}
#'   \item{gdpw2}{a numeric vector giving the natural log of real GDP per worker in
#'    1985 international prices. From Alvarez et al. (1999).}
#' }
#'
#' @references Kevin M. Quinn. 2004. ``Bayesian Factor Analysis for Mixed
#' Ordinal and Continuous Response.'' \emph{Political Analyis}. 12: 338-353.
#'
#' @source Mike Alvarez, Jose Antonio Cheibub, Fernando Limongi, and Adam
#' Przeworski. 1999. ``ACLP Political and Economic Database.''
#'
#' Witold J. Henisz. 2002. ``The Political Constraint Index (POLCON) Dataset.''
#'
#' Monty G. Marshall, Ted Robert Gurr, and Barbara Harff. 2002. ``State Failure
#' Task Force Problem Set.''
#'
#' @keywords datasets
NULL


#' U.S. Supreme Court Vote Matrix, Rehnquist Court (1994-2004)
#'
#' This dataframe contains a matrix of votes cast by U.S. Supreme Court
#' justices by all cases in the 1994-2004 terms.
#'
#' @name Rehnquist
#'
#' @docType data
#'
#' @format The dataframe has contains data for justices Rehnquist, Stevens,
#' O'Connor, Scalia, Kennedy, Souter, Thomas, Ginsburg, and Breyer for the
#' 1994-2004 terms of the U.S. Supreme Court.  The dataframe also contains the
#' term of the case, and a time variable that counts from term 1 to 11.  The
#' votes are coded liberal (1) and conservative (0) using the protocol of
#' Spaeth (2003).  The unit of analysis is the case citation (ANALU=0).  We are
#' concerned with formally decided cases issued with written opinions, after
#' full oral argument and cases decided by an equally divided vote
#' (DECTYPE=1,5,6,7).
#'
#' @source Harold J. Spaeth. 2005. \emph{Original United States Supreme Court
#' Database: 1953-2004 Terms.}
#'
#' @keywords datasets
NULL


#' 106th U.S. Senate Roll Call Vote Matrix
#'
#' This dataframe contains a matrix of votes cast by U.S. Senators in the 106th
#' Congress.
#'
#' @name Senate
#'
#' @docType data
#'
#' @format The dataframe contains roll call data for all Senators in the 106th
#' Senate.  The first column (id) is the ICPSR member ID number, the second
#' column (statecode) is the ICPSR state code, the third column (party) is the
#' member's state name, and the fourth column (member) is the member's name.
#' This is followed by all roll call votes (including unanimous ones) in the
#' 106th.  Nay votes are coded 0, yea votes are coded 1, and NAs are missing
#' votes.
#'
#' @source Keith Poole. 2005. \emph{106th Roll Call Vote Data}.
#'
#' @keywords datasets
NULL

#' U.S. Supreme Court Vote Matrix
#'
#' This dataframe contains a matrix votes cast by U.S. Supreme Court justices
#' in all cases in the 2000 term.
#'
#' @name SupremeCourt
#'
#' @docType data
#'
#' @format The dataframe has contains data for justices Rehnquist, Stevens,
#' O'Connor, Scalia, Kennedy, Souter, Thomas, Ginsburg, and Breyer for the 2000
#' term of the U.S. Supreme Court.  It contains data from 43 non-unanimous
#' cases. The votes are coded liberal (1) and conservative (0) using the
#' protocol of Spaeth (2003).  The unit of analysis is the case citation
#' (ANALU=0).  We are concerned with formally decided cases issued with written
#' opinions, after full oral argument and cases decided by an equally divided
#' vote (DECTYPE=1,5,6,7).
#'
#' @source Harold J. Spaeth. 2005. \emph{Original United States Supreme Court
#' Database: 1953-2004 Terms.} \url{http://supremecourtdatabase.org}.
#'
#' @keywords datasets
NULL

#' Euro 2016 data
#'
#' Data on head-to-head outcomes from the 2016 UEFA European Football
#' Championship.
#'
#' @name Euro2016
#'
#' @docType data
#'
#' @format This dataframe contains all of the head-to-head results from
#' Euro 2016. This includes results from both the group stage and the
#' knock-out rounds.
#' \describe{
#'   \item{dummy.rater}{An artificial "dummy" rater equal to 1 for all
#' matches. Included so that \code{Euro2016} can be used directly with
#' \code{MCMCpack}'s models for pairwise comparisons.}
#'   \item{team1}{The home team}
#'   \item{team2}{The away team }
#'   \item{winner}{The winner of the match. \code{NA} if a draw.}
#' }
#' 
#'
#' 
#' @source \url{https://en.wikipedia.org/wiki/UEFA_Euro_2016}
#' 
#' @keywords datasets
NULL


