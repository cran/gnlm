\name{gnlr3}
\title{Fit a Generalized Nonlinear Regression Model for Three
Parameter Distributions}
\usage{
gnlr3(y, dist="normal", mu=NULL, shape=NULL, family=NULL,
	linear=NULL, pmu=NULL, pshape=NULL, pfamily=NULL, censor=F,
	exact=F, wt=1, delta=1, common=F, envir=sys.frame(sys.parent()),
	print.level=0, typsiz=abs(p), ndigit=10, gradtol=0.00001,
	stepmax=10*sqrt(p\%*\%p), steptol=0.00001, iterlim=100, fscale=1)
}
\alias{gnlr3}
\alias{coefficients.gnlr3}
\alias{deviance.gnlr3}
\alias{df.residual.gnlr3}
\alias{fitted.values.gnlr3}
\alias{residuals.gnlr3}
\alias{weights.gnlr3}
\arguments{
\item{y}{The response vector for uncensored data, two columns for
censored data, with the second being the censoring indicator (1:
uncensored, 0: right censored, -1: left censored.), or an object of
class, response (created by \code{\link{restovec}}) or repeated
(created by \code{\link{rmna}}).}
\item{distribution}{Either a character string containing the name of
the distribution or a function giving the -log likelihood and calling
the location and shape functions.}
\item{mu}{A user-specified function of \code{pmu}, and possibly
\code{linear}, giving the regression equation for the location. This
may contain a linear part as the second argument to the function.
It may also be a formula beginning with ~, specifying either a linear
regression function for the location parameter in the Wilkinson and
Rogers notation or a general function with named unknown
parameters. If none is supplied, the location is taken to be constant
unless the linear argument is given.}
\item{shape}{A user-specified function of \code{pshape}, and possibly
\code{linear}, giving the regression equation for the dispersion or shape
parameter. This may contain a linear part as the second argument
to the function. It may also be a formula beginning with ~, specifying
either a linear regression function for the shape parameter in the
Wilkinson and Rogers notation or a general function with named unknown
parameters. If none is supplied, this parameter is taken to be
constant unless the linear argument is given. This parameter is the
logarithm of the usual one.}
\item{family}{A user-specified function of \code{pfamily}, and possibly
\code{linear}, for the regression equation of the third (family)
parameter of the distribution. This may contain a linear part that is
the second argument to the function. It may also be a formula
beginning with ~, specifying either a linear regression function for
the family parameter in the Wilkinson and Rogers notation or a general
function with named unknown parameters. If neither is supplied, this
parameter is taken to be constant unless the linear argument is
given. In most cases, this parameter is the logarithm of the usual one.}
\item{linear}{A formula beginning with ~, specifying the linear part
of the regression function for the location parameters or list of
three such expressions for the location, shape, and/or family parameters.}
\item{pmu}{Vector of initial estimates for the location parameters.
If \code{mu} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pshape}{Vector of initial estimates for the shape parameters.
If \code{shape} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pfamily}{Vector of initial estimates for the family parameters.
If \code{family} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{exact}{If TRUE, fits the exact likelihood function for
continuous data by integration over intervals of observation,
i.e. interval censoring.}
\item{wt}{Weight vector.}
\item{delta}{Scalar or vector giving the unit of measurement (always
one for discrete data) for each response value, set to unity by
default - for example, if a response is measured to two decimals,
delta=0.01. If the response is transformed, this must be multiplied by
the Jacobian. The transformation cannot contain unknown
parameters. For example, with a log transformation,
\code{delta=1/y}. (The delta values for the censored response are
ignored.)}
\item{common}{If TRUE, at least two of \code{mu}, \code{shape}, and
\code{family} must both be functions with, as argument, a vector of
parameters having some or all elements in common between them so that
indexing is in common between them; all parameter estimates must be
supplied in \code{pmu}. If FALSE, parameters are distinct between the
two functions and indexing starts at one in each function.}
\item{envir}{Environment in which model formulae are to be
interpreted or a data object of class, repeated, tccov, or tvcov.
If \code{y} has class \code{repeated}, it is used as the
environment.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{gnlr3} fits user specified nonlinear regression equations to one,
two, or all three parameters of three parameter distributions
(Box-Cox transformed normal, generalized inverse Gauss, generalized
logistic, Hjorth, generalized gamma, Burr, generalized Weibull,
power exponential, Student t, and generalized extreme value).

Nonlinear regression models can be supplied as formulae where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)
}
\value{
A list of class gnlr is returned.
The printed output includes the -log likelihood (not the deviance),
the corresponding AIC, the maximum likelihood estimates, standard
errors, and correlations. A list is returned that contains all of the
relevant information calculated, including error codes.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{fmr}}, \code{\link{finterp}}, \code{\link{glm}},
\code{\link{gnlr}}, \code{\link{lm}}.
}
\examples{
y <- rgamma(20,2,1)
sex <- c(rep(0,10),rep(1,10))
sexf <- gl(2,10)
age <- rpois(20,10)
# linear regression with the generalized gamma distribution
mu <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr3(y, dist="gamma", mu=mu, pmu=rep(1,3), pshape=0, pfamily=0)
# or equivalently
gnlr3(y, dist="gamma", mu=~sexf+age, pmu=rep(1,3),
	pshape=0, pfamily=0)
# or
gnlr3(y, dist="gamma", linear=~sex+age, pmu=rep(1,3),
	pshape=0, pfamily=0)
# or
gnlr3(y, dist="gamma", mu=~b0+b1*sex+b2*age,
	pmu=list(b0=1,b1=1,b2=1), pshape=0, pfamily=0)
#
# nonlinear regression with generalized gamma distribution
mu <- function(p, linear) p[1]+exp(linear)
gnlr3(y, dist="gamma", mu=mu, linear=~sex+age, pmu=rep(1,4),
	pshape=0, pfamily=0)
# or equivalently
gnlr3(y, dist="gamma", mu=~b4+exp(b0+b1*sex+b2*age),
	pmu=list(b0=1,b1=1,b2=1,b4=1), pshape=0, pfamily=0)
#
# include regression for the shape parameter with same mu function
shape <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr3(y, dist="gamma", mu=mu, linear=~sexf+age, shape=shape,
	pmu=rep(1,4), pshape=rep(0,3), pfamily=0)
# or equivalently
gnlr3(y, dist="gamma", mu=mu, linear=list(~sexf+age,~sex+age,NULL),
	pmu=rep(1,4), pshape=rep(0,3), pfamily=0)
# or
gnlr3(y, dist="gamma", mu=mu, linear=~sexf+age,
	shape=~c0+c1*sex+c2*age, pmu=rep(1,4),
	pshape=list(c0=0,c1=0,c2=0), pfamily=0)
# include regression for the family parameter with same mu
# and shape functions
family <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr3(y, dist="gamma", mu=mu, linear=~sexf+age, shape=shape,
	family=shape, pmu=rep(1,4), pshape=rep(0,3), pfamily=rep(0,3))
# or equivalently
gnlr3(y, dist="gamma", mu=mu,
	linear=list(~sex+age,~sex+age,~sex+age), pmu=rep(1,4),
	pshape=rep(0,3), pfamily=rep(0,3))
# or
gnlr3(y, dist="gamma", mu=~b4+exp(b0+b1*sex+b2*age),
	shape=~c0+c1*sex+c2*age, family=~d0+d1*sex+d2*age,
	pmu=list(b0=1,b1=1,b2=1,b4=1), pshape=list(c0=0,c1=0,c2=0),
	pfamily=list(d0=0,d1=0,d2=0))
}
