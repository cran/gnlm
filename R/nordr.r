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
#     nordr(y, distribution="proportional", mu, linear=NULL, pmu, 
#	pintercept, wt=NULL, envir=sys.frame(sys.parent()), print.level=0,
#	ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
#	typsiz=abs(p), stepmax=10*sqrt(p%*%p))
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models for ordinal responses.

nordr <- function(y, distribution="proportional", mu, linear=NULL, pmu, 
	pintercept, wt=NULL, envir=sys.frame(sys.parent()), print.level=0,
	ndigit=10, gradtol=0.00001, steptol=0.00001, fscale=1, iterlim=100,
	typsiz=abs(p), stepmax=10*sqrt(p%*%p)){
lf <- function(p){
	g <- exp(mu1(p[1:npl])+block%*%p[npl1:np])
	g <- g/(1+g)
	if(mdl==1){
		g <- c(g,ext)
		g <- g[1:nlen]/g[nrows1:nlenr]
		g <- ifelse(g>=1,0.99,g)}
	-sum(pwt*(resp*log(g)+(1-resp)*log(1-g)))}
lf3 <- function(p){
	mu <- mu1(p[1:npl])
	g <- exp(mu*(y-1)+resp%*%p[npl1:np])/
	exp(mu%o%(0:my)+matrix(rep(cumsum(c(0,0,p[npl1:np])),nrows),ncol=my+1,byrow=T))%*%ext
	-sum(pwt*log(g))}
call <- sys.call()
tmp <- c("proportional odds","continuation ratio","adjacent categories")
mdl <- match(distribution <- match.arg(distribution,tmp),tmp)
if(!is.vector(y,mode="numeric"))stop("y must be a vector")
if(min(y)!=1)stop("ordinal values must start at 1")
else if(any(y!=trunc(y)))stop("ordinal values must be integers")
else my <- max(y)-1
nrows <- length(y)
nrows1 <- nrows+1
nlen <- my*nrows
nlenr <- nlen+nrows
npl <- length(pmu)
npl1 <- npl+1
if(missing(pintercept)||length(pintercept)!=my-1)
	stop(paste(my-1,"initial values of intercept parameters must be supplied"))
if(inherits(mu,"formula"))linear <- mu
lin1a <- mu2 <- name <- NULL
if(inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(inherits(envir,"repeated"))"repeated"
		else "tccov"
	if(inherits(linear,"formula")){
		lin1a <- finterp(linear)
		class(linear) <- c(class(linear),type)}
	name <- paste(deparse(substitute(envir)))
	if(is.function(mu)){
		mu2 <- mu
		attributes(mu2) <- attributes(fnenvir(mu))
		class(mu) <- type
		mu <- fnenvir(mu,envir=envir,name=name)}}
if(inherits(linear,"formula")){
	mu1 <- finterp(linear,envir=envir,name=name)
	npt1 <- length(attr(mu1,"parameters"))
	if(is.matrix(attr(mu1,"model"))){
		if(all(dim(attr(mu1,"model"))==1)){
			if(is.function(mu)){
				dm1 <- attr(mu1,"model")
				mu1 <- function(p) mu(p,p[1]*rep(1,n))}
			else {
				tmp <- attributes(mu1)
				mu1 <- function(p) p[1]*rep(1,n)
				attributes(mu1) <- tmp}}
		else {
			if(nrow(attr(mu1,"model"))!=nrows)stop("mu model matrix does not match number of response observations")
			if(is.function(mu)){
				dm1 <- attr(mu1,"model")
				linear <- mu1
				mu1 <- function(p) mu(p,dm1%*%p[1:npt1])}}}
	else {
		if(npl!=npt1){
			cat("\nParameters are ")
			cat(attr(mu1,"parameters"),"\n")
			stop(paste("pmu should have",npt1,"estimates"))}
		if(is.list(pmu)){
			if(!is.null(names(pmu))){
				o <- match(attr(mu1,"parameters"),names(pmu))
				pmu <- unlist(pmu)[o]
				if(sum(!is.na(o))!=length(pmu))stop("invalid estimates for mu - probably wrong names")}
			else pmu <- unlist(pmu)}}
	if(npl<npt1)stop("Not enough initial estimates for mu")}
else if(!is.function(mu)){
	mu1 <- function(p) p[1]*rep(1,nrows)
	npt1 <- 1}
else {
	mu1 <- mu
	if(length(mu1(pmu))==1)mu1 <- function(p) mu(p)*rep(1,nrows)}
if(is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn"))attributes(fnenvir(mu))
		else attributes(mu)}
		else attributes(fnenvir(mu1))}
nlp <- if(is.function(mu)){
		if(is.null(linear))length(attr(mu1,"parameters"))
		else length(attr(mu1,"parameters"))-1+npt1}
       else npt1
if(nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
if(any(is.na(mu1(pmu))))stop("The location model returns NAs: probably invalid initial values")
if(mdl==1)ext <- rep(1,nrows)
else if(mdl==3)ext <- rep(1,my+1)
if(mdl==3)resp <- NULL
else resp <- matrix(as.integer(y==1),ncol=1)
block <- NULL
pwt <- matrix(as.integer(y<3),ncol=1,nrow=nrows)
for(i in 2:my){
	resp <- cbind(resp,as.integer(y<=i))
	block <- cbind(block,as.integer(c(rep(0,nrows*(i-1)),
		rep(1,nrows),rep(0,nrows*(my-i)))))
	pwt <- cbind(pwt,as.integer(y<i+2))}
if(mdl!=1)resp <- 1-resp
if(mdl!=3){
	resp <- as.vector(resp)
	pwt <- as.vector(pwt)}
else pwt <- rep(1,length(y))
if(!is.null(wt)){
	if(!is.vector(wt,mode="numeric"))stop("wt must be a vector")
	else if(length(wt)!=nrows)stop(paste("wt must have length",nrows))
	if(mdl==3)pwt <- wt
	else pwt <- rep(wt,my)*pwt}
p <- c(pmu,pintercept)
np <- length(p)
if(mdl==3){
	if(fscale==1)fscale <- lf3(p)
	if(is.na(lf3(p)))stop("Likelihood returns NAs: probably invalid initial values")
	z <- nlm(lf3, p, hessian=T, print.level=print.level,
		typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
		steptol=steptol, iterlim=iterlim, fscale=fscale)}
else {
	if(fscale==1)fscale <- lf(p)
	if(is.na(lf(p)))stop("Likelihood returns NAs: probably invalid initial values")
	z <- nlm(lf, p, hessian=T, print.level=print.level,
		typsiz=typsiz, ndigit=ndigit, gradtol=gradtol, stepmax=stepmax,
		steptol=steptol, iterlim=iterlim, fscale=fscale)}
maxlike <- z$minimum
a <- qr(z$hessian)
if(a$rank==np)cov <- solve(z$hessian)
else cov <- matrix(NA,ncol=np,nrow=np)
se <- sqrt(diag(cov))
corr <- cov/(se%o%se)
dimnames(corr) <- list(1:np,1:np)
if(!is.null(mu2))mu1 <- mu2
if(!is.null(lin1a))linear <- lin1a
z1 <- list(
   call=call,
   distribution=distribution,
   wt=wt,
   maxlike=maxlike,
   aic=maxlike+np,
   mu=mu1,
   linear=linear,
   coefficients=z$estimate[1:npl],
   np=np,
   npl=npl1-1,
   nrows=nrows,
   intercept=z$estimate[npl1:np],
   cov=cov,
   corr=corr,
   se=se,
   iterations=z$iter,
   code=z$code)
class(z1) <- "nordr"
z1}

print.nordr <- function(z, digits = max(3, .Options$digits - 3)){
	m <- z$states
	cat("\nCall:",deparse(z$call),sep="\n")
	cat("\n")
	cat(z$distribution,"model\n\n")
	if(z$code>2)cat("Warning: no convergence - error",z$code,"\n\n")
	cat("-Log likelihood   ",z$maxlike,"\n")
	cat("AIC               ",z$aic,"\n")
	cat("Iterations        ",z$iterations,"\n")
	cat("\nLocation coefficients\n")
	if(inherits(z$mu,"formulafn")){
		cat("Location function:\n")
		if(!is.null(attr(z$mu,"formula")))cat(deparse(attr(z$mu,"formula")),sep="\n")
		else if(!is.null(attr(z$mu,"model"))){
			t <- deparse(attr(z$mu,"model"))
			t[1] <- sub("expression\\(","",t[1])
			t[length(t)] <- sub("\\)$","",t[length(t)])
			cat(t,sep="\n")}
		if(inherits(z$linear,"formulafn"))
			cat("Linear part: ",deparse(attr(z$linear,"formula")),sep="\n")}
	cname <- if(is.matrix(attr(z$mu,"model")))colnames(attr(z$mu,"model"))
		else if(length(grep("linear",attr(z$mu,"parameters")))>0)
		attr(z$mu,"parameters")[grep("\\[",attr(z$mu,"parameters"))]
		else attr(z$mu,"parameters")
	if(!is.null(z$linear)&&!is.null(attr(z$linear,"parameters")))
		cname <- c(colnames(attr(z$linear,"model")),cname)
	coef.table <- cbind(z$coef,z$se[1:z$npl])
	dimnames(coef.table) <- list(cname,c("estimate","s.e."))
	print.default(coef.table, digits=digits, print.gap=2)
	cat("\nIntercept coefficients\n")
	coef.table <- cbind(z$intercept,z$se[(z$npl+1):z$np])
	dimnames(coef.table) <- list(paste("b[",2:(z$np-z$npl+1),"]",sep=""),
			     c("estimate","s.e."))
	print.default(coef.table, digits=digits, print.gap=2)
	cat("\nCorrelation matrix\n")
	print.default(z$corr, digits=digits)
}
