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
#     nlr(y, mu=NULL, pmu=NULL, distribution="normal", wt=1, delta=1,
#	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(pmu),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(pmu%*%pmu),
#	steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models for distributions
# in the exponential family.

nlr <- function(y, mu=NULL, pmu=NULL, distribution="normal", wt=1, delta=1,
	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(pmu),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(pmu%*%pmu),
	steptol=0.00001, iterlim=100, fscale=1){
call <- sys.call()
if(!missing(distribution))distribution <- match.arg(distribution,c("normal","inverse Gauss","gamma"))
if(missing(pmu))stop("Initial parameter estimates must be supplied")
np <- length(pmu)
if(!is.vector(y,mode="numeric"))stop("y must be a vector")
if(any(is.na(y)))stop("NAs in y - use rmna")
n <- length(y)
mu2 <- NULL
respenv <- inherits(y,"repeated")
envname <- if(respenv)paste(deparse(substitute(y)))
	else NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")){
	if(is.null(envname))envname <- paste(deparse(substitute(envir)))
	if(inherits(mu,"formula")){
		mu2 <- if(respenv)finterp(mu,envir=y,name=envname)
			else finterp(mu,envir=envir,name=envname)
		class(mu) <- c(class(mu),type)}
	else if(is.function(mu)){
		tmp <- parse(text=paste(deparse(mu))[-1])
		class(mu) <- if(respenv||inherits(envir,"repeated"))"repeated"
			else "tccov"
		mu <- if(respenv)fnenvir(mu,envir=y,name=envname)
			else fnenvir(mu,envir=envir,name=envname)
		mu2 <- mu
		if(respenv)attr(mu2,"model") <- tmp}}
if(inherits(mu,"formula")){
	mu1 <- mu
	mu <- if(respenv)finterp(mu,envir=y,name=envname)
		else finterp(mu,envir=envir,name=envname)
	npt1 <- length(attr(mu,"parameters"))
	if(is.matrix(attr(mu,"model"))){
		if(all(dim(attr(mu,"model"))==1)){
			tmp <- attributes(mu)
			mu <- function(p) p[1]*rep(1,n)
			attributes(mu) <- tmp}}
	else {
		if(np!=npt1){
			cat("\nParameters are ")
			cat(attr(mu,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}}
if(is.null(attributes(mu))){
	attributes(mu) <- if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,envir=y))
			else attributes(fnenvir(mu,envir=envir))}
		else attributes(mu)}
nlp <- if(is.function(mu1))length(attr(mu,"parameters"))
       else npt1
if(nlp!=np)stop(paste("pmu should have",nlp,"initial estimates"))
if(missing(mu)||!is.function(mu))stop("A mean function or formula must be supplied")
fn <- switch(distribution,
normal=function(p) sum(wt*(y-mu(p))^2),
gamma=function(p) -sum(wt*(log(y/mu(p))-(y-mu(p))/mu(p))),
"inverse Gauss"=function(p) sum(wt*((y-mu(p))^2)/(y*mu(p)^2)))
if(fscale==1)fscale <- fn(pmu)
if(is.na(fn(pmu)))
	stop("Non-numerical function value: probably invalid initial values")
z0 <- nlm(fn, p=pmu, hessian=T, print.level=print.level, typsiz=typsiz,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
	steptol=steptol, iterlim=iterlim, fscale=fscale)
if(length(delta)==1)delta <- rep(delta,length(y))
if(length(wt)==1)wt <- rep(wt,length(y))
n <- sum(wt)
disp <- z0$minimum/n
p <- z0$estimate
switch(distribution,
normal=maxlike <- length(y)*(log(2*pi*disp)+1)/2,
gamma=maxlike <- (sum(wt*(y/mu(p)+log(mu(p))-log(y)))+n*log(disp))/
	disp+n*lgamma(1/disp)+sum(log(y)*wt),
"inverse Gauss"=maxlike <- (sum(wt)*(log(disp*2*pi)+1)+3*sum(log(y)*wt))/2)
maxlike <- maxlike-sum(log(delta))
fitted.values <-  as.vector(mu(z0$estimate))
residuals <-  y-fitted.values
if(np==1)cov <- 1/z0$hessian
else {
	a <- qr(z0$hessian)
	if(a$rank==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)}
cov <- 2*cov*z0$minimum/sum(wt)
se <- sqrt(diag(cov))
z1 <- list(
	call=call,
	distribution=distribution,
	delta=delta,
	mu=mu,
	prior.weights=wt,
	maxlike=maxlike,
	dispersion=disp,
	fitted.values=fitted.values,
	residuals=residuals,
	aic=maxlike+np+1,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	np=np,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "nlr"
return(z1)}

residuals.nlr <- function(z) z$residuals
fitted.values.nlr <- function(z) z$fitted.values
coefficients.nlr <- function(z) z$coefficients
weights.nlr <- function(z) z$prior.weights
df.residual.nlr <- function(z) z$df
deviance.nlr <- function(z) 2*z$maxlike

print.nlr <- function(z) {
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat(z$distribution,"distribution\n\n")
	cat("Mean function:\n")
	if(!is.null(attr(z$mu,"formula")))cat(deparse(attr(z$mu,"formula")),sep="\n")
	else if(!is.null(attr(z$mu,"model"))){
		t <- deparse(attr(z$mu,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cat("\n-Log likelihood   ",z$maxlike,"\n")
	cat("Degrees of freedom",z$df,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n\n")
	cat("Mean parameters:\n")
	coef.table <- cbind(z$coefficients[1:z$np], z$se[1:z$np])
	if(inherits(z$mu,"formulafn"))
		cname <- if(is.matrix(attr(z$mu,"model")))
				colnames(attr(z$mu,"model"))
			else attr(z$mu,"parameters")
	else cname <- seq(1,z$np)
	dimnames(coef.table) <- list(cname, c("estimate", "se"))
	print.default(coef.table, digits=4, print.gap=2)
	cat("\nDispersion estimate:",z$dispersion,"\n")
	if(z$np>1){
		cat("\nCorrelations:\n")
		dimnames(z$corr) <- list(seq(1,z$np),seq(1,z$np))
		print.default(z$corr, digits=4)}
	invisible(z)}
