\name{nordr}
\title{Nonlinear Ordinal Regression}
\usage{
nordr(y, distribution="proportional", mu, linear=NULL, pmu, 
	pintercept, wt=NULL, envir=sys.frame(sys.parent()),
	print.level=0, ndigit=10, gradtol=0.00001,
	steptol=0.00001, fscale=1, iterlim=100, typsiz=abs(p),
	stepmax=10*sqrt(p\%*\%p))
}
\arguments{
\item{y}{A vector of ordinal responses, integers numbered from one to
the maximum value.}
\item{distribution}{The ordinal distribution: proportional odds,
continuation ratio, or adjacent categories.}
\item{mu}{User-specified function of \code{pmu}, and possibly \code{linear},
giving the logistic regression equation. This must contain the first
intercept. It may contain a linear part as the second argument to the
function. It may also be a formula beginning with ~, specifying a
logistic regression function for the location parameter, either a
linear one using the Wilkinson and Rogers notation or a general
function with named unknown parameters. If none is supplied, the
location is taken to be constant unless the linear argument is given.}
\item{linear}{A formula beginning with ~, specifying the linear part
of the logistic regression function.}
\item{pmu}{Vector of initial estimates for the regression parameters,
including the first intercept. If \code{mu} is a formula with unknown
parameters, their estimates must be supplied either in their order of
appearance in the expression or in a named list.}
\item{pintercept}{Vector of initial estimates for the contrasts with
the first intercept parameter (difference in intercept for successive
categories): two less than the number of different ordinal values.}
\item{wt}{Weight vector for use with contingency tables.}
\item{envir}{Environment in which model formulae are to be
interpreted.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{nordr} fits arbitrary nonlinear regression functions (with
logistic link) to ordinal response data by proportional odds,
continuation ratio, or adjacent categories.

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)
}
\value{
A list of class nordr is returned.
The printed output includes the -log likelihood (not the deviance),
the corresponding AIC, the maximum likelihood estimates, standard
errors, and correlations. A list is returned that contains all of the
relevant information calculated, including error codes.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{finterp}}, \code{\link{fmr}, \code{\link{glm}},
\code{\link{glmm}}, \code{\link{gnlmm}}, \code{\link{gnlr}},
\code{\link{gnlr3}}, \code{\link{nlr}}.
}
\examples{
# McCullagh (1980) JRSS B42, 109-142
# Tonsil size: 2x3 contingency table
y <- c(1:3,1:3)
carrier <- c(rep(0,3),rep(1,3))
carrierf <- gl(2,3,6)
wt <- c(19,29,24,497,560,269)
pmu <- c(-1,0.5)
mu <- function(p) c(rep(p[1],3),rep(p[1]+p[2],3))
# proportional odds
# with mean function
nordr(y, dist="prop", mu=mu, pmu=pmu, wt=wt, pintercept=1.5)
# using Wilkinson and Rogers notation
nordr(y, dist="prop", mu=~carrierf, pmu=pmu, wt=wt, pintercept=1.5)
# using formula with unknowns
nordr(y, dist="prop", mu=~b0+b1*carrier, pmu=pmu, wt=wt, pintercept=1.5)
# continuation ratio
nordr(y, dist="cont", mu=mu, pmu=pmu, wt=wt, pintercept=1.5)
# adjacent categories
nordr(y, dist="adj", mu=~carrierf, pmu=pmu, wt=wt, pintercept=1.5)
}
