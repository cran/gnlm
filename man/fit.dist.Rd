\name{fit.dist}
\title{Fit a Probability Distribution to Frequency Data}
\usage{
fit.dist(y, ni, dist, breaks=F, delta=1, censor=F, exact=T, plot=F,
add=F, main, xlab, ...)
}
\alias{fit.dist}
\arguments{
\item{y}{Vector of observations.}
\item{ni}{Corresponding vector of frequencies.}
\item{dist}{Character string specifying the distribution.}
\item{breaks}{If TRUE, \code{y} contains breaks between categories
instead of mid-points.}
\item{delta}{Scalar or vector giving the unit of measurement (always
one for discrete data) for each response value, set to unity by
default. For example, if a response is measured to two decimals,
delta=0.01.}
\item{censor}{If TRUE, the last category is right censored.}
\item{exact}{If FALSE, uses the approximations for certain
distributions in Lindsey (1995).}
\item{plot}{If TRUE, plots the histogram of observed frequencies and
the frequency polygon of fitted values.}
\item{add}{If TRUE, adds a new frequency polygon of fitted values
without replotting the histogram.}
\item{others}{Plotting control options.}
}
\description{
\code{fit.dist} fits the distributions in Chapter 4 of Lindsey (1995):
binomial, Poisson, negative binomial, geometric, zeta, normal, log
normal, inverse Gauss, logistic, exponential, Pareto, gamma, and Weibull.
}
\references{
Lindsey, J.K. (1995) Introductory Statistics: A Modelling Approach.
Oxford: Oxford University Press. 
}
\author{J.K. Lindsey}
\examples{
f <- c(215, 1485, 5331, 10649, 14959, 11929, 6678, 2092, 342)
y <- seq(0,8)
fit.dist(y, f, "binomial", plot=T, xlab="Number of males",
	main="Distribution of males in families of 8 children")
#
f <- c(1,1,6,3,4,3,9,6,5,16,4,11,6,11,3,4,5,6,4,4,5,1,1,4,1,2,
	0,2,0,0,1)
y <- seq(1100,4100,by=100)
fit.dist(y, f, "normal", delta=100, plot=T,
	xlab="Monthly salary (dollars)",
	main="Distribution of women mathematicians' salaries")
fit.dist(y, f, "log normal", delta=100, plot=T, add=T, lty=3)
fit.dist(y, f, "logistic", delta=100, exact=F, plot=T, add=T,
	lty=2)
}
