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
#     fmr(y, distribution="normal", mu=NULL, mix=NULL, linear=NULL,
#	pmu=NULL, pmix=NULL, pshape=NULL, censor="right", exact=F,
#	wt=1, delta=1, common=F, envir=sys.frame(sys.parent()),
#	print.level=0, typsiz=abs(p), ndigit=10, gradtol=0.00001,
#	stepmax=10*sqrt(p%*%p), steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models with a variety of
# distributions and a mixture in the tail(s).

fmr <- function(y, distribution="normal", mu=NULL, mix=NULL, linear=NULL,
	pmu=NULL, pmix=NULL, pshape=NULL, censor="right", exact=F,
	wt=1, delta=1, common=F, envir=sys.frame(sys.parent()), print.level=0,
	typsiz=abs(p), ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p),
	steptol=0.00001, iterlim=100, fscale=1){

pinvgauss <- function(y,m,s){
	t <- y/m
	v <- sqrt(y*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}
plaplace <- function(y){
	t <- exp(-abs(y))/2
	ifelse(y<0,t,1-t)}
plevy <- function(y, m, s)
	.C("plevy",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(1),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		DUP=F)$res

call <- sys.call()
if(!missing(distribution)&&!is.function(distribution)){
	distribution <- match.arg(distribution,c("binomial","beta binomial",
	"double binomial","mult binomial","Poisson","negative binomial",
	"double Poisson","mult Poisson","gamma count","Consul","geometric",
	"normal","inverse Gauss","logistic","exponential","gamma","Weibull",
	"extreme value","Pareto","Cauchy","Student t","Laplace","Levy"))}
if(common){
	if(!is.function(mu))stop("with common parameters, mu must be a function")
	if(!is.function(mix))stop("with common parameters, mix must be a function")
	if(!is.null(linear))stop("linear cannot be used with common parameters")}
if(!missing(pmu))npl <- length(pmu)
else npl <- 0
if(!missing(pmix))npm <- length(pmix)
else npm <- 0
sht <- distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"
if(sht&&missing(pshape))
	stop("An estimate of the shape parameter must be given")
np <- npl+npm+sht
if(is.function(distribution)){
	fcn <- distribution
	distribution <- "own"}
respenv <- inherits(y,"repeated")
envname <- if(respenv)paste(deparse(substitute(y)))
	else NULL
lin1 <- lin2 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]}
else lin1 <- linear
if(inherits(mu,"formula"))lin1 <- mu
if(inherits(mix,"formula"))lin2 <- mix
lin1a <- lin2a <- mu2 <- mixt2 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(respenv||inherits(envir,"repeated"))"repeated"
		else "tccov"
	if(is.null(envname))envname <- paste(deparse(substitute(envir)))
	if(inherits(lin1,"formula")){
		if(is.function(mu)){
			lin1a <- if(respenv)finterp(lin1,envir=y,name=envname)
			else finterp(lin1,envir=envir,name=envname)}
		class(lin1) <- c(class(lin1),type)}
	if(inherits(lin2,"formula")){
		if(is.function(mix)){
			lin2a <- if(respenv)finterp(lin2,envir=y,name=envname)
			else finterp(lin2,envir=envir,name=envname)}
		class(lin2) <- c(class(lin2),type)}
	if(is.function(mu)){
		tmp <- parse(text=paste(deparse(mu))[-1])
		class(mu) <- type
		mu <- if(respenv)fnenvir(mu,envir=y,name=envname)
			else fnenvir(mu,envir=envir,name=envname)
		mu2 <- mu
		if(respenv)attr(mu2,"model") <- tmp}
	if(is.function(mix)){
		tmp <- parse(text=paste(deparse(mix))[-1])
		class(mix) <- type
		mix <- if(respenv)fnenvir(mix,envir=y,name=envname)
			else fnenvir(mix,envir=envir,name=envname)
		mixt2 <- mix
		if(respenv)attr(mixt2,"model") <- tmp}}
if(inherits(lin1,"formula")){
	mu1 <- if(respenv)finterp(lin1,envir=y,name=envname)
		else finterp(lin1,envir=envir,name=envname)
	npt1 <- length(attr(mu1,"parameters"))
	if(is.matrix(attr(mu1,"model"))){
		if(all(dim(attr(mu1,"model"))==1)){
			if(is.function(mu)){
				lin1 <- mu1
				mu1 <- function(p) mu(p,p[npl]*rep(1,n))}
			else {
				tmp <- attributes(mu1)
				mu1 <- function(p) p[1]*rep(1,n)
				attributes(mu1) <- tmp}}
		else {
			if(is.function(mu)){
				lf <- if(inherits(mu,"formulafn"))length(attr(mu,"parameters"))
					else length(if(respenv)attr(fnenvir(mu,envir=y),"parameters")
						else attr(fnenvir(mu,envir=envir),"parameters"))
				dm1 <- attr(mu1,"model")
				lin1 <- mu1
				mu1 <- function(p) mu(p,dm1%*%p[lf:(lf+npt1-1)])}}}
	else {
		if(is.function(mu)){
			warning("ignoring mu function\n")
			mu <- mu2 <- NULL}
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
	mu1 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else {
	mu1 <- mu
	if(length(mu1(pmu))==1)mu1 <- function(p) mu(p)*rep(1,n)}
if(is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn")){
			if(respenv)attributes(fnenvir(mu,envir=y))
			else attributes(fnenvir(mu,envir=envir))}
		else attributes(mu)}
		else {
			if(respenv)attributes(fnenvir(mu1,envir=y))
			else attributes(fnenvir(mu1,envir=envir))}}
nlp <- if(is.function(mu)){
		if(is.null(lin1))length(attr(mu1,"parameters"))
		else length(attr(mu1,"parameters"))-1+npt1}
       else npt1
if(!common&&nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
npl1 <- if(common) 1 else npl+1
if(inherits(lin2,"formula")){
	mixt1 <- if(respenv)finterp(lin2,envir=y,start=npl1,name=envname)
		else finterp(lin2,envir=envir,start=npl1,name=envname)
	npt2 <- length(attr(mixt1,"parameters"))
	if(is.matrix(attr(mixt1,"model"))){
		if(all(dim(attr(mixt1,"model"))==1)){
			if(is.function(mix)){
				lin2 <- mixt1
				mixt <- function(p) {
					mf <- mix(p[npl1:np],p[np]*rep(1,n))
					exp(mf)/(1+exp(mf))}}
			else {
				mixt <- function(p) {
					mf <- p[npl1]*rep(1,n)
					exp(mf)/(1+exp(mf))}
				mixt2 <- fnenvir(function(p) {
					mf <- p[1]*rep(1,n)
					exp(mf)/(1+exp(mf))})
				attributes(mixt) <- attributes(mixt1)}
			rm(mixt1)}
		else {
			dm2 <- attr(mixt1,"model")
			if(is.function(mix))mixt <- function(p) {
				lfm <- if(inherits(mix,"formulafn"))length(attr(mix,"parameters"))
					else length(if(respenv)attr(fnenvir(mix,envir=y),"parameters")
						else attr(fnenvir(mix,envir=envir),"parameters"))
				lin2 <- mixt1
				mf <- mix(p[npl1:np],dm2%*%p[(npl+lfm):np])
				exp(mf)/(1+exp(mf))}
			else {
				mixt <- function(p) {
					mf <- dm2%*%p[npl1:(npl1+npt2-1)]
					exp(mf)/(1+exp(mf))}
				attributes(mixt) <- attributes(mixt1)}}}
	else {
		if(is.function(mix)){
			warning("ignoring mix function\n")
			mix <- mixt2 <- NULL}
		if(npm!=npt2){
			cat("\nParameters are ")
			cat(attr(mixt1,"parameters"),"\n")
			stop(paste("pmix should have",npt2,"estimates"))}
		mixt <- function(p) {
			mf <- mixt1(p)
			exp(mf)/(1+exp(mf))}
		attributes(mixt) <- attributes(mixt1)
		if(is.list(pmix)){
			if(!is.null(names(pmix))){
				o <- match(attr(mixt,"parameters"),names(pmix))
				pmix <- unlist(pmix)[o]
				if(sum(!is.na(o))!=length(pmix))stop("invalid estimates for mix - probably wrong names")}
			else pmix <- unlist(pmix)}}}
else if(!is.function(mix)){
	mixt <- function(p) exp(p[npl1])/(1+exp(p[npl1]))*rep(1,n)
	mixt2 <- fnenvir(function(p) exp(p[1])/(1+exp(p[1]))*rep(1,n))
	npt2 <- 1}
else mixt <- function(p) {
	mf <- mix(p[npl1:np])
	exp(mf)/(1+exp(mf))}
if(is.null(attributes(mixt))){
	attributes(mixt) <- if(is.function(mix)){
		if(!inherits(mix,"formulafn")){
			if(respenv)attributes(fnenvir(mix,envir=y))
			else attributes(fnenvir(mix,envir=envir))}
		else attributes(mix)}
		else {
			if(respenv)attributes(fnenvir(mixt,envir=y))
			else attributes(fnenvir(mixt,envir=envir))}}
nlp <- if(is.function(mix)){
		if(is.null(lin2))length(attr(mixt,"parameters"))
		else length(attr(mixt,"parameters"))-1+npt2}
       else npt2
if(!common&&nlp!=npm)stop(paste("pshape should have",nlp,"initial estimates"))
if(common){
	nlp <- length(unique(c(attr(mu1,"parameters"),attr(mixt,"parameters"))))-shfn
	if(nlp!=npl)stop(paste("with a common parameter model, pmu should contain",nlp,"estimates"))}
p <- c(pmu,pmix,pshape)
if(respenv){
	if(inherits(envir,"repeated")&&(length(y$response$nobs)!=length(envir$response$nobs)||any(y$response$nobs!=envir$response$nobs)))stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&!is.na(y$response$wt))wt <- y$response$wt
	if(!is.null(y$response$delta))delta <- y$response$delta
	y <- response(y)}
else if(inherits(y,"response")){
	if(!is.null(y$wt)&&!is.na(y$wt))wt <- y$wt
	if(!is.null(y$delta))delta <- y$delta
	y <- response(y)}
if(any(is.na(y)))stop("NAs in y - use rmna")
if(distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="double Poisson"||distribution=="mult Poisson"||
	distribution=="gamma count"||distribution=="Consul"){
	if(!is.vector(y,mode="numeric"))stop("y must be a vector")
	n <- length(y)
	censor <- NULL
	cens <- ifelse(y==0,1,0)}
else {
	if(length(dim(y))!=2||ncol(y)!=2)
		stop(paste("Two column matrix required for response:",
		if(distribution=="binomial"||distribution=="beta binomial"||
			distribution=="double binomial"||
			distribution=="mult binomial")"successes and failures"
		else "times and censor indicator"))
	else {
		n <- nrow(y)
		if(distribution=="binomial"||distribution=="beta binomial"||
			distribution=="double binomial"||
			distribution=="mult binomial"){
			if(missing(censor))
				stop("Censoring must be left, right, or both")
			if(censor!="left"&&censor!="right"&&censor!="both")
				stop("Censoring must be left, right, or both")
			lcens <- ifelse((censor=="left"|censor=="both")&
				y[,1]==0,1,0)
			rcens <-ifelse((censor=="right"|censor=="both")&
				y[,2]==0,1,0)
			if(censor=="both"){
				lcens <- lcens/2
				rcens <- rcens/2}
			n <- nrow(y)
			nn <- y[,1]+y[,2]}
		else {
			if(any(delta<=0&y[,2]==1))
				stop("All deltas for uncensored data must be positive")
			else {
				delta <- ifelse(delta<=0,0.000001,delta)
				delta <- ifelse(y[,1]-delta/2<=0,delta-0.00001
				,delta)}
			y[,2] <- as.integer(y[,2])
			if(any(y[,2]!=-1&y[,2]!=0&y[,2]!=1))
				stop("Censor indicator must be -1, 0, or 1")
			if(censor!="left"&&censor!="right")
				stop("Censoring must be left or right")
			if(censor=="left"&!any(y[,2]==-1))
				stop("No left censored observations")
			if(censor=="right"&!any(y[,2]==0))
				stop("No right censored observations")
			cens <- as.integer(y[,2]==1)
			b <- as.integer((censor=="right"&y[,2]==0)|
				(censor=="left"&y[,2]==-1))
			r <- as.integer(censor=="left"&y[,2]==0)
			l <- as.integer(censor=="right"&y[,2]==-1)
			lc <- ifelse(censor=="left",1,0)
			rc <- ifelse(censor=="right",-1,1)}}
	if(distribution=="double Poisson"||distribution=="mult Poisson")
				my <- min(3*max(y),100)}
if((distribution!="normal"&&distribution!="logistic"&&distribution!="Cauchy"&&
	distribution!="Laplace"&&distribution!="Student t"&&
	distribution!="Poisson"&&distribution!="negative binomial"&&
	distribution!="Consul"&&distribution!="double Poisson"&&
	distribution!="mult Poisson"&&distribution!="gamma count"&&
	distribution!="binomial"&& distribution!="beta binomial"&&
	distribution!="double binomial"&&distribution!="mult binomial")&&
	(any(y[,1]<=0)))stop("All response values must be > 0")
else if((distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="gamma count"||distribution=="double Poisson"||
	distribution=="mult Poisson"||distribution=="Consul"||
	distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial")
	&&(any(y<0)))stop("All response values must be >= 0")
if(min(wt)<0)stop("All weights must be non-negative")
if(length(wt)==1)wt <- rep(wt,n)
if(length(delta)==1)delta <- rep(delta,n)
if(any(is.na(mu1(pmu))))stop("The location model returns NAs: probably invalid initial values")
if(distribution=="Levy"&&any(y[,1]<=mu1(p)))
	stop("location parameter must be strictly less than corresponding observation")
if(sht&&any(is.na((mixt(p)))))stop("The mix function returns NAs: probably invalid initial values")
ret <- switch(distribution,
	binomial={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			-wt*log((1-s)*(lcens+rcens)+s*m^y[,1]*(1-m)^y[,2])}
		const <- -wt*(lchoose(nn,y[,1]))},
	"beta binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			v <- exp(p[np])
			t <- v*m
			u <- v*(1-m)
			-wt*log((1-s)*(lcens+rcens)+s*beta(y[,1]+t,y[,2]+u)/
				beta(t,u))}
		const <- -wt*(lchoose(nn,y[,1]))},
	"double binomial"={
		fcn <- function(p) {
			-wt*log((1-s)*(lcens+rcens)+s*exp(.C("ddb",
				as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(p[np])),
				as.integer(n),as.double(wt),
				res=double(n),DUP=F)$res))}
		const <- 0},
	"mult binomial"={
		fcn <- function(p) {
			-wt*log((1-s)*(lcens+rcens)+s*exp(.C("dmb",
				as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(p[np])),
				as.integer(n),as.double(wt),
				res=double(n),DUP=F)$res))}
		const <- 0},
	Poisson={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			-wt*log((1-s)*cens+s*exp(-m)*m^y)}
		const <- wt*lgamma(y+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			t <- exp(p[np])
			-wt*log((1-s)*cens+s*gamma(y+t)/gamma(t)
				*t^t*m^y/(t+m)^(y+t))}
		const <- wt*lgamma(y+1)},
	"double Poisson"={
		fcn <- function(p) {
			-wt*log((1-s)*cens+s*exp(.C("ddp",as.integer(y),
				as.integer(my),as.double(mu1(p)),
				as.double(exp(p[np])),as.integer(length(y)),
				as.double(wt),res=double(length(y)),DUP=F)$res))}
		const <- 0},
	"mult Poisson"={
		fcn <- function(p) {
			-wt*log((1-s)*cens+s*exp(.C("dmp",as.integer(y),
				as.integer(my),as.double(mu1(p)),
				as.double(exp(p[np])),as.integer(length(y)),
				as.double(wt),res=double(length(y)),DUP=F)$res))}
		const <- 0},
	"gamma count"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			t <- exp(p[np])
			-wt*log((1-s)*cens+s*ifelse(y==0,1-pgamma(m*t,
				(y+1)*t,1),pgamma(m*t,y*t+(y==0),1)-
				pgamma(m*t,(y+1)*t,1)))}
		const <- 0},
	Consul={
		fcn <- function(p) {
			m <- mu1(p)
			s <- mixt(p)
			u <- exp(p[np])
			-wt*log((1-s)*cens+s*m*exp(-(m+y*(u-1))/u-y*p[np])*
				(m+y*(u-1))^(y-1))}
		const <- wt*lgamma(y+1)},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pn <- pnorm(y[,1],m,t)
				-wt*log(s*cens*(pnorm(y[,1]+delta/2,m,t)-
					pnorm(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pn-lc))*b
					+s*(r+pn*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				pn <- pnorm(y[,1],m,exp(p[np]/2))
				-wt*log(s*cens*exp((-(p[np]+(y[,1]-m)^2
					*exp(-p[np]))/2))+(1-cens)*
					((1+s*(rc*pn-lc))*b+s*(r+pn*(l-r))))}
			const <- wt*cens*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pit <- pinvgauss(y[,1],m,t)
				-wt*log(s*cens*(pinvgauss(y[,1]+delta/2,m,t)-
					pinvgauss(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pit-lc))*b
					+s*(r+pit*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pit <- pinvgauss(y[,1],m,t)
				-wt*log(s*cens*exp(-(p[np]+(y[,1]-m)^2/
					(y[,1]*t*m^2))/2)
					+(1-cens)*((1+s*(rc*pit-lc))*b
					+s*(r+pit*(l-r))))}
			const <- wt*cens*(log(2*pi*y[,1]^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])*sqrt(3)/pi
				pl <- plogis(y[,1],m,t)
				-wt*log(s*cens*(plogis(y[,1]+delta/2,m,t)-
					plogis(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])*sqrt(3)/pi
				y1 <- (y[,1]-m)/t
				pl <- plogis(y[,1],m,t)
				-wt*log(s*cens*exp(-y1-log(t)
					-2*log(1+exp(-y1)))
					+(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r))))}
			const <- -wt*cens*log(delta)}},
        "Student t"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ps <- pt(y[,1]-m,t)
				-wt*log(s*cens*(pt(y[,1]+delta/2-m,t)-
					pt(y[,1]-delta/2-m,t))
					+(1-cens)*((1+s*(rc*ps-lc))*b
					+s*(r+ps*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ps <- pt(y[,1]-m,t)
				-wt*log(s*cens*gamma((t+1)/2)/gamma(t/2)*
					exp(-p[np]/2-((t+1)/2)*
					log(1+(y[,1]-m)^2/t))
					+(1-cens)*((1+s*(rc*ps-lc))*b
					+s*(r+ps*(l-r))))}
			const <- wt*cens*(log(pi)/2-log(delta))}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pc <- pcauchy(y[,1],m,t)
				-wt*log(s*cens*(pcauchy(y[,1]+delta/2,m,t)-
					pcauchy(y[,1]-delta/2,m,t))
					+(1-cens)*((1+s*(rc*pc-lc))*b
					+s*(r+pc*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np]/2)
				pc <- pcauchy(y[,1],m,t)
				-wt*log(s*cens/(t*(1+(y[,1]-m)^2/t^2))
					+(1-cens)*((1+s*(rc*pc-lc))*b
					+s*(r+pc*(l-r))))}
			const <- -wt*cens*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plaplace((y[,1]-m)/t)
				-wt*log(s*cens*(plaplace((y[,1]+delta/2-m)/t)
					-plaplace((y[,1]-delta/2-m)/t))+
					(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plaplace((y[,1]-m)/t)
				-wt*log(s*cens*exp(-abs(y[,1]-m)/t-p[np])+
					(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r))))}
			const <- -wt*cens*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plevy(y[,1],m,t)
				-wt*log(s*cens*(plevy(y[,1]+delta/2,m,t)
					-plevy(y[,1]-delta/2,m,t))+
					(1-cens)*((1+s*(rc*pl-lc))*b
					+s*(r+pl*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pl <- plevy(y[,1],m,t)
				-wt*log(s*cens*sqrt(t/(2*pi))*log(y[,1]-m)^-1.5
					*exp(-t/(2*(y[,1]-m)))+(1-cens)*
					((1+s*(rc*pl-lc))*b+s*(r+pl*(l-r))))}
			const <- -wt*cens*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- mixt(p)
				u <- exp(p[np])
				t <- 1/(mu1(p)*u)
				pp <- 1-(1+y[,1]*t)^-u
				-wt*log(s*cens*((1+(y[,1]-delta/2)*t)^-u-
					(1+(y[,1]+delta/2)*t)^-u)
					+(1-cens)*((1+s*(rc*pp-lc))*b
					+s*(r+pp*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- mixt(p)
				u <- exp(p[np])
				t <- 1/(mu1(p)*u)
				pp <- 1-(1+y[,1]*t)^-u
				-wt*log(s*cens*u*t*(1+y[,1]*t)^(-(u+1))+
					(1-cens)*
					((1+s*(rc*pp-lc))*b+s*(r+pp*(l-r))))}
			const <- -wt*cens*log(delta)}},
	exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				u <- exp(-y[,1]/m)
				-wt*log(s*cens*(-exp(-(y[,1]+delta/2)/m)+
					exp(-(y[,1]-delta/2)/m))
					+(1-cens)*((1+s*(rc*(1-u)-lc))*b
					+s*(r+(1-u)*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				u <- exp(-y[,1]/m)
				-wt*log(s*cens*exp(-y[,1]/m)/m
					+(1-cens)*((1+s*(rc*(1-u)-lc))*b
					+s*(r+(1-u)*(l-r))))}
			const <- -wt*cens*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				u <- m/t
				pg <- pgamma(y[,1],t,u)
				-wt*log(s*cens*(pgamma(y[,1]+delta/2,t,u)-
					pgamma(y[,1]-delta/2,t,u))
					+(1-cens)*((1+s*(rc*pg-lc))*b
					+s*(r+pg*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				u <- m/t
				pg <- pgamma(y[,1],t,u)
				-wt*log(s*cens*y[,1]^(t-1)*exp(-y[,1]/u)/
					(u^t*gamma(t))
					+(1-cens)*((1+s*(rc*pg-lc))*b
					+s*(r+pg*(l-r))))}
			const <- -wt*cens*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				pw <- pweibull(y[,1],t,m)
				-wt*log(s*cens*(pweibull(y[,1]+delta/2,t,m)-
					pweibull(y[,1]-delta/2,t,m))
					+(1-cens)*((1+s*(rc*pw-lc))*b
					+s*(r+pw*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				v <- y[,1]/m
				u <- exp(-v^t)
				-wt*log(s*cens*t*v^(t-1)*u/m+
					(1-cens)*((1+s*(rc*(1-u)-lc))*b
					+s*(r+(1-u)*(l-r))))}
			const <- -wt*cens*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				ey <- exp(y[,1])
				pw <- pweibull(ey,t,m)
				-wt*log(s*cens*(pweibull(ey+ey*delta/2,
					t,m)-pweibull(ey-ey*delta/2,t,m))+
					(1-cens)*((1+s*(rc*pw-lc))*b
					+s*(r+pw*(l-r))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- mixt(p)
				t <- exp(p[np])
				v <- exp(y[,1])/m
				u <- exp(-v^t)
				-wt*log(s*cens*t*v^(t-1)*u/m+
					(1-cens)*((1+s*(rc*(1-u)-lc))*b
					+s*(r+(1-u)*(l-r))))}
			const <- -wt*cens*log(delta)}},
	own={const <- 0})
fn <- function(p) sum(fcn(p))
if(fscale==1)fscale <- fn(p)
if(is.na(fn(p)))stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(fn, p=p, hessian=T, print.level=print.level, typsiz=typsiz,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
z0$minimum <- z0$minimum+sum(const)
fitted.values <- if(distribution=="binomial"||distribution=="beta binomial"||
		distribution=="double binomial"||distribution=="mult binomial")
		as.vector((y[,1]+y[,2])*mu1(z0$estimate))
	else as.vector(mu1(z0$estimate))
residuals <- if(distribution!="Poisson"&&distribution!="negative binomial"&&
	distribution!="Consul"&&distribution!="double Poisson"&&
	distribution!="mult Poisson"&&distribution!="gamma count")
		y[,1]-fitted.values
	else y-fitted.values
if(np==1)cov <- 1/z0$hessian
else {
	a <- qr(z0$hessian)
	if(a$rank==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)}
se <- sqrt(diag(cov))
like.comp <- as.vector(fcn(z0$estimate)+const)
if(!is.null(mu2))mu1 <- mu2
if(!is.null(mixt2))mixt <- mixt2
if(!is.null(lin1a))lin1 <- lin1a
if(!is.null(lin2a))lin2 <- lin2a
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	mu=mu1,
	mix=mixt,
	linear=list(lin1,lin2),
	common=common,
	prior.weights=wt,
	censor=censor,
	maxlike=z0$minimum,
	fitted.values=fitted.values,
	residuals=residuals,
	like.comp=like.comp,
	aic=z0$minimum+np,
	df=sum(wt)-np,
	coefficients=z0$estimate,
	npl=npl,
	npm=npm,
	nps=as.numeric(sht),
	npf=0,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlr"
return(z1)}
