\name{nlr}
\title{Nonlinear Regression}
\usage{
nlr(y, mu=NULL, pmu=NULL, distribution="normal", wt=1, delta=1,
	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(pmu),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(pmu\%*\%pmu),
	steptol=0.00001, iterlim=100, fscale=1)
}
\alias{nlr}
\alias{coefficients.nlr}
\alias{deviance.nlr}
\alias{df.residual.nlr}
\alias{fitted.values.nlr}
\alias{residuals.nlr}
\alias{weights.nlr}
\arguments{
\item{y}{The response vector.}
\item{mu}{A function of \code{p} giving the regression equation for
the mean or a formula beginning with ~, specifying either a linear
regression function in the Wilkinson and Rogers notation or a general
nonlinear function with named unknown parameters.}
\item{pmu}{Vector of initial estimates of the parameters.
If \code{mu} is a formula with unknown parameters, their estimates
must be supplied either in their order of appearance in the expression
or in a named list.}
\item{distribution}{The distribution to be used: normal, gamma, or
inverse Gauss.}
\item{wt}{Weight vector.}
\item{delta}{Scalar or vector giving the unit of measurement for each
response value, set to unity by default. For example, if a response is
measured to two decimals, delta=0.01. If the response is transformed,
this must be multiplied by the Jacobian. For example, with a log
transformation, \code{delta=1/y}.}
\item{envir}{Environment in which model formulae are to be
interpreted.}
\item{others}{Arguments controlling \code{\link{nlm}}.}
}
\description{
\code{nlr} fits a user-specified nonlinear regression equation by
least squares (normal) or its generalization for the gamma and inverse
Gauss distributions.
}
\value{
A list of class nlr is returned.
The printed output includes the -log likelihood (not the deviance),
the corresponding AIC, the parameter estimates, standard
errors, and correlations. A list is returned that contains all of the
relevant information calculated, including error codes.

A nonlinear regression model can be supplied as a formula where
parameters are unknowns. Factor variables cannot be used and
parameters must be scalars. (See \code{\link{finterp}}.)
}
\author{J.K. Lindsey}
\seealso{
\code{\link{finterp}}, \code{\link{fmr}}, \code{\link{glm}},
\code{\link{glmm}}, \code{\link{gnlmm}}, \code{\link{gnlr}},
\code{\link{gnlr3}}, \code{\link{lm}}.
}
\examples{
x <- c(3,5,0,0,0,3,2,2,2,7,4,0,0,2,2,2,0,1,3,4)
y <- c(-26,-191,2,1,10,-25,9,0,-8,-598,-80,3,-4,-7,-10,-22,-5,-8,-52,-84)
# 2*x^2-2*x^3+rnorm(20,0,10)
# linear regression
mu1 <- function(p) p[1]+p[2]*x
summary(lm(y~x))
nlr(y, mu=mu1, pmu=c(3,2))
# or
nlr(y, mu=~x, pmu=c(3,2))
# or
nlr(y, mu=~b0+b1*x, pmu=c(3,2))
# nonlinear regression
mu2 <- function(p) p[1]*x^p[2]+p[3]*x^p[4]
nlr(y, mu=mu2, pmu=c(1,1,-1,2))
# or
nlr(y, mu=~b0*x^b1+c0*x^c1, pmu=list(b0=1,b1=1,c0=-1,c1=2))
# with gamma distribution
y <- rgamma(20,2,0.2+2*exp(0.1*x))
nlr(y, dist="gamma", mu=~b0+c0*exp(c1*x),
	pmu=list(b0=0.2,c0=3,c1=0.2))
}
