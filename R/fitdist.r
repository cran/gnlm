#
#  gnlm : A Library of Special Functions for Nonlinear Regression
#  Copyright (C) 1998 J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     fit.dist(y, ni, dist="normal", breaks=F, delta=1, censor=F,
#	exact=T, plot=F, add=F, main, xlab, ...)
#
#  DESCRIPTION
#
#    A function to fit common distributions to frequency data and to
# plot the resulting curve along with the histogram

fit.dist <- function(y, ni, dist="normal", breaks=F, delta=1, censor=F,
	exact=T, plot=F, add=F, main, xlab, ...)
{
	if(!missing(dist))dist <- match.arg(dist,c("binomial","Poisson",
		"negative binomial","geometric","zeta","normal","log normal",
		"inverse Gauss","logistic","exponential","Pareto","gamma",
		"Weibull"))
	if(missing(main))main <- paste("Histogram of", deparse(substitute(y)))
	if(missing(xlab))xlab <- paste(deparse(substitute(y)))
	n <- length(ni)
	if (length(delta)==1) delta <- rep(delta,n)
	if (breaks){
		if(length(y)!=n+1)
			stop("Breaks vector must be one longer than frequency vector")
		yi <- (y[1:n]+y[2:(n+1)])/2
		delta <- diff(y)}
	else yi <- y
	pi.hat <- ni/sum(ni)
        ybar <- weighted.mean(yi,ni)
	s2 <- weighted.mean((yi-ybar)^2,ni)
	if (dist=="binomial"){
		m <- length(yi)-1
		nu <- ybar/m
		pi.tilde <- gamma(m+1)/gamma(yi+1)/gamma(m-yi+1)*
			nu^yi*(1-nu)^(m-yi)
		param <- nu
		names(param) <- "nu.hat"
		p <- 1}
	else if (dist=="Poisson"){
		pi.tilde <- exp(-ybar)*ybar^yi/gamma(yi+1)
		param <- ybar
		names(param) <- "mu.hat"
		p <- 1}
	else if (dist=="geometric"){
		nu <- 1/(1+ybar)
		pi.tilde <- nu*(1-nu)^yi
		param <- nu
		names(param) <- "nu.hat"
		p <- 1}
	else if (dist=="negative binomial"){
		nu <- ybar/s2
		gam <- ybar^2/(s2-ybar)
		if(exact){
			nu <- log(nu/(1-nu))
			fcn <- function(p)
				-sum(ni*(lgamma(yi+p[2])-lgamma(p[2])+
				p[2]*p[1]-(yi+p[2])*log(1+exp(p[1]))))
#				p[2]*log(p[1])+yi*log(1-p[1])))
			z <- nlm(fcn, p=c(nu, gam), stepmax=sqrt(nu^2+gam^2)/2)
			nu <- exp(z$estimate[1])
			nu <- nu/(1+nu)
			gam <- z$estimate[2]}
		pi.tilde <- gamma(yi+gam)/gamma(yi+1)/gamma(gam)*
			nu^gam*(1-nu)^yi
		param <- c(nu, gam)
		names(param) <- c("nu.hat","gamma.hat")
		p <- 2}
	else if (dist=="zeta"){
		pi.tilde <- 1/yi
		nu <- sum(pi.tilde)
		pi.tilde <- pi.tilde/nu
		rho <- round(pi.hat[1]/pi.tilde[1]+0.1)
		if(exact){
			fcn <- function(p) {
				if(censor) const <- sum(yi^(-p[1]))
				else const <- sum(seq(1,30)^(-p[1]))
				sum(ni*(p[1]*log(yi)+log(const)))}
			z <- nlm(fcn, p=rho, stepmax=1)
			rho <- z$estimate[1]}
		pi.tilde <- yi^(-rho)
		nu <- sum(pi.tilde)
		pi.tilde <- pi.tilde/nu
		param <- rho
		names(param) <- "rho.hat"
		p <- 1}
        else if (dist=="normal"){
		mu.hat <- ybar
		sigma2.hat <- s2
		pi.tilde <- exp(-(yi-mu.hat)^2/(2*sigma2.hat))/
			sqrt(2*pi*sigma2.hat)
		param <- c(mu.hat,sigma2.hat)
		names(param) <- c("mu.hat","sigma2.hat")
		p <- 2}
        else if (dist=="log normal"){
		mu.hat <- weighted.mean(log(yi),ni)
		sigma2.hat <- weighted.mean((log(yi)-mu.hat)^2,ni)
		pi.tilde <- exp(-(log(yi)-mu.hat)^2/(2*sigma2.hat))/
			(yi*sqrt(2*pi*sigma2.hat))
		param <- c(mu.hat,sigma2.hat)
		names(param) <- c("mu.hat","sigma2.hat")
		p <- 2}
        else if (dist=="inverse Gauss"){
		mu.hat <- ybar
		sigma2.hat <- weighted.mean(1/yi,ni)-(1/ybar)
		pi.tilde <- exp(-(yi-mu.hat)^2/(2*yi*sigma2.hat*mu.hat^2))/
			sqrt(2*pi*yi^3*sigma2.hat)
		param <- c(mu.hat,sigma2.hat)
		names(param) <- c("mu.hat","sigma2.hat")
		p <- 2}
	else if (dist=="logistic"){
		mu.hat <- ybar
		sigma <- sqrt(s2)
		if(exact){
			fcn <- function(p)
				sum(ni*(pi*(yi-p[1])/p[2]/sqrt(3)+log(p[2])
				+2*log(1+exp(-pi*(yi-p[1])/p[2]/sqrt(3)))))
			z <- nlm(fcn, p=c(mu.hat,sigma), stepmax=10)
			mu.hat <- z$estimate[1]
			sigma <- z$estimate[2]}
		pi.tilde <- exp(-pi*(yi-mu.hat)/(sigma*sqrt(3)))
		pi.tilde <- pi*pi.tilde/(sigma*sqrt(3)*(1+pi.tilde)^2)
		param <- c(mu.hat,sigma)
		names(param) <- c("mu.hat","sigma.hat")
		p <- 2}
        else if (dist=="exponential"){
		pi.tilde <- (1/ybar)*exp(-yi/ybar)
		param <- ybar
		names(param) <- "mu.hat"
		p <- 1}
        else if (dist=="Pareto"){
		delta.hat <- yi[1]-delta[1]/2
		alpha.hat <- sum(ni)/sum(ni*log(yi/delta.hat))
		pi.tilde <- alpha.hat*delta.hat^alpha.hat/(yi^(alpha.hat+1))
		param <- c(alpha.hat,delta.hat)
		names(param) <- c("alpha.hat","delta.hat")
		p <- 2}
        else if (dist=="gamma"){
		alpha.hat <- ybar^2/s2
		mu.hat <- ybar
		if(exact){
			fcn <- function(p)
				-sum(ni*(p[2]*log(p[2])
				+(p[2]-1)*log(yi)-p[2]*log(p[1])
				-p[2]*yi/p[1]-lgamma(p[2])))
			z <- nlm(fcn, p=c(mu.hat,alpha.hat), stepmax=10)
			mu.hat <- z$estimate[1]
			alpha.hat <- z$estimate[2]}
		pi.tilde <- (alpha.hat^alpha.hat)*(yi^(alpha.hat-1))/
			(mu.hat^alpha.hat)*exp(-alpha.hat*yi/mu.hat)/
			gamma(alpha.hat)
		param <- c(alpha.hat,mu.hat)
		names(param) <- c("alpha.hat","mu.hat")
		p <- 2}
        else if (dist=="Weibull"){
		tamp <- ybar^2/(s2+ybar^2)
		Alpha.Weibull.fn <- function(y){
			alpha.trans.fn <- function(al)
				gamma(1+1/al)*gamma(1+1/al)/gamma(1+2/al)
			tol <- 0.001
			al.start <- 0.0001
			al.end <- 50
			al.mid <- 0.5*(al.start+al.end)
			y.tamp <- alpha.trans.fn(al.mid)
			while (abs(y.tamp-y)>tol){
				if ((y.tamp-y)>0) al.end <- al.mid
				else al.start <- al.mid
				al.mid <- 0.5*(al.start+al.end)
				y.tamp <- alpha.trans.fn(al.mid)}
			al.mid}
		alpha.hat <- Alpha.Weibull.fn(tamp)
		mu.hat <- weighted.mean(yi^alpha.hat,ni)
		if(exact){
			fcn <- function(p)
				-sum(ni*(log(p[2])+(p[2]-1)*log(yi)
				-log(p[1])-yi^p[2]/p[1]))
			z <- nlm(fcn, p=c(mu.hat,alpha.hat), stepmax=10)
			mu.hat <- z$estimate[1]
			alpha.hat <- z$estimate[2]}
		pi.tilde <- alpha.hat*yi^(alpha.hat-1)/mu.hat*
			exp(-yi^alpha.hat/mu.hat)
		param <- c(alpha.hat,mu.hat)
		names(param) <- c("alpha.hat","mu.hat")
		p <- 2}
	pi.tilde <- pi.tilde*delta
        if (censor)
		pi.tilde[length(pi.tilde)] <- 1-sum(pi.tilde[1:(length(pi.tilde)-1)])
	dev.comp <- rep(0,length(pi.tilde))
	dev.comp[ni>0] <- -2*ni[ni>0]*log(pi.tilde[ni>0]/pi.hat[ni>0])
	deviance <- sum(dev.comp)
	AIC <- deviance+2*p
	resid <- (ni-sum(ni)*pi.tilde)/sqrt(sum(ni)*pi.tilde)
	result.output <- c(ybar,s2,param,deviance,AIC)
	names(result.output) <- c("mean","variance",names(param),"deviance","AIC")
	cat(dist," distribution,","  n = ",sum(ni),"\n",sep="")
	print(result.output)
        cat("\n")
	if(censor)nn <- n-1
	else nn <- n
	if (plot){
		if (!add){
			cum.histo (c(yi-delta/2, yi[n]+delta[n]/2), ni, prob=T,
				ylab="Probability", main, xlab, ...)
			lines(yi[1:nn], pi.tilde[1:nn]/delta[1:nn])}
		else lines(yi[1:nn], pi.tilde[1:nn]/delta[1:nn], add=T, ...)}
	if(length(unique(delta))==1)
		df <- data.frame(yi,ni,pi.hat,pi.tilde,dev.comp,resid)
	else
		df <- data.frame(yi,ni,delta,pi.hat,pi.tilde,dev.comp,resid)
#	list(parameters=result.output,df=df)
	df}

cum.histo <- function (breaks, freq, prob = F,
	main = paste("Histogram of", deparse(substitute(breaks))),
	xlab = deparse(substitute(breaks)), ylab,
        xlim = range(breaks), ...)
{
	if (prob) {
                freq <- freq/(sum(freq) * diff(breaks))
                if (missing(ylab)) 
                        ylab <- "Relative Frequency"
        }
        else if (missing(ylab)) 
                ylab <- "Frequency"
	plot(breaks, c(freq, 0), type = "n", main = main,
		xlab = xlab, ylab = ylab, ...)
	rect(breaks[-length(breaks)], 0, breaks[-1], freq, border = par("fg"))}
