\name{rs3}
\title{Fit a Three-factor Box-Tidwell Nonlinear Response Surface}
\usage{
rs3(y, x1, x2, x3, power=c(1,1,1), weight=rep(1,length(y)),
	family=normal, iterlim=20)
}
\arguments{
\item{y}{Response variable}
\item{x1}{First covariate}
\item{x2}{Second covariate}
\item{x3}{Third covariate}
\item{power}{Initial estimates of the three power transformations}
\item{weight}{Weight vector}
\item{family}{glm family}
\item{iterlim}{Iteration limit}
}
\description{
\code{rs3} fits a three-covariate power-transformed response surface by
iterating the function, \code{\link{glm}}.}
\value{
A list of class, rs, is returned containing the model and the power
estimates.
}
\author{J.K. Lindsey}
\seealso{
\code{\link{lm}}, \code{\link{glm}}, \code{\link{gnlr}}, \code{\link{gnlr3}},
\code{\link{fmr}}, \code{\link{rs2}}.
}
\examples{
x1 <- rep(1:4,5)
x2 <- rep(1:5,rep(4,5))
x3 <- c(rep(1:3,6),1,2)
y <- rpois(20,1+2*sqrt(x1)+3*log(x2)+1/x3+4*x1+log(x2)^2+1/x3^2+
	2*sqrt(x1)*log(x2)+sqrt(x1)/x3+log(x2)/x3)
rs3(y, x1, x2, x3, family=poisson)
}
