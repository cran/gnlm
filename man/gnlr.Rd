\name{gnlr}
\title{Fit a Generalized Nonlinear Regression Model}
\usage{
gnlr(y, distribution="normal", mu=NULL, shape=NULL, linear=NULL,
	pmu=NULL, pshape=NULL, exact=F, wt=1, delta=1, shfn=F,
	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p\%*\%p),
	steptol=0.00001, iterlim=100, fscale=1)
}
\alias{gnlr}
\alias{coefficients.gnlr}
\alias{deviance.gnlr}
\alias{df.residual.gnlr}
\alias{fitted.values.gnlr}
\alias{residuals.gnlr}
\alias{weights.gnlr}
\arguments{
\item{y}{A response vector for uncensored data, a two column matrix
for binomial data or censored data, with the second column being the
censoring indicator (1: uncensored, 0: right censored, -1: left
censored), or an object of class, response (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).}
\item{distribution}{Either a character string containing the name of the
distribution or a function giving the -log likelihood and calling the
location and shape functions.}
\item{mu}{A user-specified function of \code{pmu}, and possibly
\code{linear}, giving the regression equation for the location. This
may contain a linear part as the second argument to the function. It
may also be a formula beginning with ~, specifying either
a linear regression function for the location parameter in the
Wilkinson and Rogers notation or a general function with named unknown
parameters. If none is supplied, the location is taken to be constant
unless the linear argument is given.}
\item{shape}{A user-specified function of \code{pshape}, and possibly
\code{linear} and/or \code{mu}, giving the regression equation for the
dispersion or shape parameter. This may contain a linear part as the
second argument to the function and the location function as last
argument (in which case \code{shfn} must be set to TRUE). It may also
be a formula beginning with ~, specifying either a linear regression
function for the shape parameter in the Wilkinson and Rogers notation
or a general function with named unknown parameters. If none is
supplied, this parameter is taken to be constant unless the linear
argument is given. This parameter is the logarithm of the usual one.}
\item{linear}{A formula beginning with ~, specifying the linear part
of the regression function for the location parameter or list of two
such expressions for the location and/or shape parameters.}
\item{pmu}{Vector of initial estimates for the location parameters.
If \code{mu} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pshape}{Vector of initial estimates for the shape parameters.
If \code{shape} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{exact}{If TRUE, fits the exact likelihood function for
continuous data by integration over intervals of observation,
i.e. interval censoring.}
\item{wt}{Weight vector.}
\item{delta}{Scalar or vector giving the unit of measurement (always
one for discrete data) for each response value, set to unity by
default. For example, if a response is measured to two decimals,
delta=0.01. If the response is transformed, this must be multiplied by
the Jacobian. The transformation cannot contain unknown
parameters. For example, with a log transformation,
\code{delta=1/y}. (The delta values for the censored response are
ignored.)}
\item{shfn}{If true, the supplied shape function depends on the
location (function). The name of this location function must be the
last argument of the shape function.}
\item{envir}{Environment in which model formulae are to be
interpreted.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{gnlr} fits user-specified nonlinear regression equations to one
or both parameters of the common one and two parameter distributions
(binomial, beta binomial, double binomial, mult(iplicative) binomial,
Poisson, negative binomial, double Poisson, mult(iplicative) Poisson,
gamma count, Consul generalized Poisson, logarithmic series,
geometric, normal, inverse Gauss, logistic, exponential, gamma,
Weibull, extreme value, Cauchy, Pareto, Laplace, and Levy; all but the
binomial-based distributions may be right and/or left censored). A
user-specified -log likelihood can also be supplied for the distribution.

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
\code{\link{lm}}, \code{\link{finterp}}, \code{\link{fmr}},
\code{\link{glm}}, \code{\link{glmm}}, \code{\link{gnlmm}},
\code{\link{gnlr3}}, \code{\link{fmr}, \code{\link{nlr}}.
}
\examples{
y <- rgamma(10,2,5)
sex <- c(rep(0,5),rep(1,5))
sexf <- gl(2,5)
age <- rpois(10,10)
# linear regression with inverse Gauss distribution
mu <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr(y, dist="inverse Gauss", mu=mu, pmu=rep(1,3), pshape=1)
# or equivalently
gnlr(y, dist="inverse Gauss", mu=~sexf+age, pmu=rep(1,3), pshape=1)
# or
gnlr(y, dist="inverse Gauss", linear=~sex+age, pmu=rep(1,3), pshape=1)
# or
gnlr(y, dist="inverse Gauss", mu=~b0+b1*sex+b2*age,
	pmu=list(b0=1,b1=1,b2=1), pshape=1)
#
# nonlinear regression with inverse Gauss distribution
mu <- function(p, linear) p[4]+exp(linear)
gnlr(y, dist="inverse Gauss", mu=mu, linear=~sex+age, pmu=rep(1,4),
	pshape=1)
# or equivalently
gnlr(y, dist="inverse Gauss", mu=~b4+exp(b0+b1*sex+b2*age),
	pmu=list(b0=1,b1=1,b2=1,b4=1), pshape=1)
# one explicit parameter in mu, three in linear, one for shape
#
# include regression for the shape parameter with same mu function
shape <- function(p) p[1]+p[2]*sex+p[3]*age
gnlr(y, dist="inverse Gauss", mu=mu, linear=~sex+age, shape=shape,
	pmu=rep(1,4), pshape=rep(1,3))
# or equivalently
gnlr(y, dist="inverse Gauss", mu=mu, linear=~sexf+age,
	shape=~sexf+age, pmu=rep(1,4), pshape=rep(1,3))
# or
gnlr(y, dist="inverse Gauss", mu=mu, linear=list(~sex+age,~sex+age),
	pmu=rep(1,4),pshape=rep(1,3))
# or
gnlr(y, dist="inverse Gauss", mu=mu, linear=~sex+age,
	shape=~c0+c1*sex+c2*age, pmu=rep(1,4),
	pshape=list(c0=1,c1=1,c2=1))
# shape as a function of the mean
shape <- function(p, mu) p[1]+p[2]*sex+p[3]*mu
gnlr(y, dist="inverse Gauss", mu=~age, shape=shape, pmu=rep(1,2),
	pshape=rep(1,3), shfn=T)
}
