\name{fmr}
\title{Fit a Generalized Nonlinear Regression Model with a Two
or Three Point Mixture}
\usage{
fmr(y, distribution="normal", mu=NULL, mix=NULL, linear=NULL, pmu=NULL,
	pshape=NULL, pmix=NULL, censor="right", exact=F, wt=1, delta=1,
	Common=F, envir=sys.frame(sys.parent()), print.level=0,
	typsiz=abs(p), ndigit=10, gradtol=0.00001,
	stepmax=10*sqrt(p\%*\%p), steptol=0.00001, iterlim=100, fscale=1)
}
\alias{fmr}
\alias{coefficients.fmr}
\alias{deviance.fmr}
\alias{df.residual.fmr}
\alias{fitted.values.fmr}
\alias{residuals.fmr}
\alias{weights.fmr}
\arguments{
\item{y}{A response vector for uncensored data, a two column matrix
for binomial data or censored data, with the second column being the
censoring indicator (1: uncensored, 0: right censored, -1: left
censored), or an object of class, response (created by
\code{\link{restovec}}) or repeated (created by \code{\link{rmna}}).}
\item{distribution}{Either a character string containing the name of the
distribution or a function giving the -log likelihood and calling the
location and mixture functions.}
\item{mu}{A user-specified function of \code{pmu}, and possibly
\code{linear}, giving the regression equation for the location. This
may contain a linear part as the second argument to the function. It
may also be a formula beginning with ~, specifying either a linear
regression function for the location parameter in the Wilkinson and
Rogers notation or a general function with named unknown
parameters. If none is supplied, the location is taken to be constant
unless the linear argument is given.}
\item{mix}{A user-specified function of \code{pmix}, and possibly
\code{linear}, giving the regression equation for the mixture
parameter. This may contain a linear part as the second argument to
the function. It may also be a formula beginning with ~, specifying
either a linear regression function for the mixture parameter in the
Wilkinson and Rogers notation or a general function with named unknown
parameters. If none is supplied, this parameter is taken to be
constant. This parameter is the logit of the mixture probability.}
\item{linear}{A formula beginning with ~, or list of two
such expressions, specifying the linear part of the regression
function for the location or location and mixture parameters.}
\item{pmu}{Vector of initial estimates for the location parameters.
If \code{mu} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{pshape}{An initial estimate for the shape parameter.}
\item{pmix}{Vector of initial estimates for the mixture parameters.
If \code{mix} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{censor}{\code{right}, \code{left}, or \code{both} indicating
where the mixing distribution is placed. \code{both} is only possible
for binomial data.}
\item{exact}{If TRUE, fits the exact likelihood function for
continuous data by integration over intervals of observation,
i.e. interval censoring.}
\item{wt}{Weight vector.}
\item{delta}{Scalar or vector giving the unit of measurement (always
one for discrete data) for each response value, set to unity by
default - for example, if a response is measured to two decimals,
delta=0.01. If the response is transformed, this must be multiplied by
the Jacobian. The transformation cannot contain unknown parameters.
For example, with a log transformation, \code{delta=1/y}.}
\item{common}{If TRUE, \code{mu} and \code{mix} must both be
functions with, as argument, a vector of parameters having some or all
elements in common between them so that indexing is in common
between them; all parameter estimates must be supplied in \code{pmu}.
If FALSE, parameters are distinct between the two functions and
indexing starts at one in each function.}
\item{envir}{Environment in which model formulae are to be
interpreted or a data object of class, repeated, tccov, or tvcov.
If \code{y} has class \code{repeated}, it is used as the
environment.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{fmr} fits user specified nonlinear regression equations to the
location parameter of the common one and two parameter distributions
(binomial, beta binomial, double binomial, multiplicative binomial,
Poisson, negative binomial, double Poisson, multiplicative Poisson,
gamma count, Consul, geometric, normal, inverse Gauss, logistic,
exponential, gamma, Weibull, extreme value, Pareto, Cauchy, Student t,
Laplace, and Levy). For the Poisson and negative binomial, the mixture
involves the zero category. For the (beta) binomial, it involves the
two extreme categories. For all other distributions, it involves
either left or right censored individuals. A user-specified -log
likelihood can also be supplied for the distribution.

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
\code{\link{finterp}}, \code{\link{glm}}, \code{\link{gnlr}},
\code{\link{gnlr3}}, \code{\link{lm}}.
}
\examples{
y <- cbind(rweibull(20,2,5),rbinom(20,1,0.7))
sex <- c(rep(0,10),rep(1,10))
sexf <- gl(2,10)
age <- rpois(20,10)
# linear regression with Weibull distribution with a point mass
# for right censored individuals
mu <- function(p) p[1]+p[2]*sex+p[3]*age
fmr(y, dist="Weibull", mu=mu, pmu=rep(1,3), pmix=1, pshape=1)
# or equivalently
fmr(y, dist="Weibull", mu=~sexf+age, pmu=rep(1,3), pmix=1, pshape=1)
# or
fmr(y, dist="Weibull", linear=~sex+age, pmu=rep(1,3), pmix=1, pshape=1)
# or
fmr(y, dist="Weibull", mu=~b0+b1*sex+b2*age, pmu=list(b0=1,b1=1,b2=1),
	pmix=1, pshape=1)
#
# nonlinear regression with Weibull distribution
mu <- function(p, linear) p[1]*exp(linear)
fmr(y, dist="Weibull", mu=mu, linear=~sex+age, pmu=rep(1,4),
	pmix=1, pshape=1)
# or equivalently
fmr(y, dist="Weibull", mu=~b4*exp(b0+b1*sex+b2*age),
	pmu=list(b0=1,b1=1,b2=1,b4=1), pmix=1, pshape=1)
#
# include logistic regression for the mixture parameter
mix <- function(p) p[1]+p[2]*sex
fmr(y, dist="Weibull", mu=~age, mix=mix, pmu=rep(1,2),
	pmix=rep(1,2), pshape=1)
# or equivalently
fmr(y, dist="Weibull", linear=list(~age,~sex), pmu=rep(1,2),
	pmix=rep(1,2), pshape=1)
# or
fmr(y, dist="Weibull", mu=~b0+b1*age, mix=~c0+c1*sex,
	pmu=list(b0=1,b1=1), pmix=list(c0=1,c1=1), pshape=1)
}
