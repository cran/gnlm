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
#     gnlr(y, distribution="normal", mu=NULL, shape=NULL, linear=NULL,
#	pmu=NULL, pshape=NULL, exact=F, wt=1, delta=1, shfn=F,
#	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(p),
#	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
#	iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models with a variety of
# one and two parameter distributions.

gnlr <- function(y, distribution="normal", pmu=NULL, pshape=NULL, mu=NULL,
	shape=NULL, linear=NULL, exact=F, wt=1, delta=1, shfn=F,
	envir=sys.frame(sys.parent()), print.level=0, typsiz=abs(p),
	ndigit=10, gradtol=0.00001, stepmax=10*sqrt(p%*%p), steptol=0.00001,
	iterlim=100, fscale=1){

pinvgauss <- function(y,m,s){
	t <- y/m
	v <- sqrt(y*s)
	pnorm((t-1)/v)+exp(2/(m*s))*pnorm(-(t+1)/v)}
plaplace <- function(y,m,s){
	u <- (y-m)/s
	t <- exp(-abs(u))/2
	ifelse(u<0,t,1-t)}
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
	"double Poisson","mult Poisson","gamma count","Consul","logarithmic",
	"geometric","normal","inverse Gauss","logistic","exponential","gamma",
	"Weibull","extreme value","Pareto","Cauchy","Laplace","Levy"))}
shp <- distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"
if(!is.null(pmu))npl <- length(pmu)
else npl <- 0
if(!is.null(pshape))nps <- length(pshape)
else nps <- 0
np <- npl+nps
if(np<1)stop("At least one parameter must be estimated")
if(is.function(distribution)){
	fcn <- distribution
	distribution <- "own"}
if(inherits(y,"repeated")){
	if(inherits(envir,"repeated")&&(length(y$response$nobs)!=length(envir$response$nobs)||any(y$response$nobs!=envir$response$nobs)))stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&!is.na(y$response$wt))wt <- y$response$wt
	if(!is.null(y$response$delta))delta <- y$response$delta
	if(!is.null(y$response$censor)&&any(y$response$censor!=1))y <- cbind(y$response$y,y$response$censor)
	else if(!is.null(y$response$n))y <- cbind(y$response$y,y$response$n-y$response$y)
	else y <- y$response$y}
else if(inherits(y,"response")){
	if(!is.null(y$wt)&&!is.na(y$wt))wt <- y$wt
	if(!is.null(y$delta))delta <- y$delta
	if(!is.null(y$censor)&&any(y$censor!=1))y <- cbind(y$y,y$censor)
	else if(!is.null(y$n))y <- cbind(y$y,y$n-y$y)
	else y <- y$y}
if(distribution=="binomial"||distribution=="double binomial"||
	distribution=="beta binomial"||distribution=="mult binomial"){
	if(length(dim(y))!=2||ncol(y)!=2)
		stop(paste("Two column matrix required for response: successes and failures"))
	if(any(y<0))stop("All response values must be positive")
	n <- nrow(y)
	nn <- y[,1]+y[,2]
	censor <- F}
else {
	censor <- length(dim(y))==2&&ncol(y)==2
	if(censor&&all(y[,2]==1)){
		y <- y[,1]
		censor <- F}
	if(!censor){
		if(!is.vector(y,mode="numeric"))stop("y must be a vector")
		n <- length(y)
		if(distribution=="double Poisson"||distribution=="mult Poisson")
			my <- 3*max(y)}}
if((distribution=="inverse Gauss"||distribution=="exponential"||
	distribution=="gamma"||distribution=="Weibull"||
	distribution=="extreme value")&&((censor&&any(y[,1]<=0))||
	(!censor&&any(y<=0))))stop("All response values must be > 0")
if((distribution=="Poisson"||distribution=="negative binomial"||
	distribution=="gamma count"||distribution=="double Poisson"||
	distribution=="mult Poisson")&&(any(y<0)))
	stop("All response values must be >= 0")
if(distribution=="logarithmic"&&any(y[wt>0]<1))
	stop("All response values must be integers > 0")
if(censor){
	n <- nrow(y)
	y[,2] <- as.integer(y[,2])
	if(any(y[,2]!=-1&y[,2]!=0&y[,2]!=1))
		stop("Censor indicator must be -1s, 0s, and 1s")
	cc <- ifelse(y[,2]==1,1,0)
	rc <- ifelse(y[,2]==0,1,ifelse(y[,2]==-1,-1,0))
	lc <- ifelse(y[,2]==-1,0,1)
	if(any(delta<=0&y[,2]==1))
		stop("All deltas for uncensored data must be positive")
	else {
		delta <- ifelse(delta<=0,0.000001,delta)
		delta <- ifelse(y[,1]-delta/2<=0,delta-0.00001,delta)}}
else {
	if(min(delta)<=0)stop("All deltas for must be positive")}
if(length(wt)==1)wt <- rep(wt,n)
else if(length(wt)!=n)stop("wt must be the same length as the other variables")
if(min(wt)<0)stop("All weights must be non-negative")
if(length(delta)==1)delta <- rep(delta,n)
else if(length(delta)!=n)stop("delta must be the same length as the other variables")
lin1 <- lin2 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]}
else lin1 <- linear
if(inherits(mu,"formula"))lin1 <- mu
if(inherits(shape,"formula"))lin2 <- shape
lin1a <- lin2a <- mu2 <- sh2 <- name <- NULL
if(inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(inherits(envir,"repeated"))"repeated"
		else "tccov"
	if(inherits(lin1,"formula")){
		lin1a <- finterp(lin1)
		class(lin1) <- c(class(lin1),type)}
	if(inherits(lin2,"formula")){
		lin2a <- finterp(lin2)
		class(lin2) <- c(class(lin2),type)}
	name <- paste(deparse(substitute(envir)))
	if(is.function(mu)){
		mu2 <- mu
		attributes(mu2) <- attributes(fnenvir(mu))
		class(mu) <- type
		mu <- fnenvir(mu,envir=envir,name=name)}
	if(is.function(shape)){
		sh2 <- shape
		attributes(sh2) <- attributes(fnenvir(shape))
		class(shape) <- type
		shape <- fnenvir(shape,envir=envir,name=name)}}
if(inherits(lin1,"formula")){
	mu1 <- finterp(lin1,envir=envir,name=name)
	npt1 <- length(attr(mu1,"parameters"))
	if(is.matrix(attr(mu1,"model"))){
		if(all(dim(attr(mu1,"model"))==1)){
			if(is.function(mu)){
				lin1 <- mu1
				mu1 <- function(p) mu(p,p[1]*rep(1,n))}
			else {
				tmp <- attributes(mu1)
				mu1 <- function(p) p[1]*rep(1,n)
				attributes(mu1) <- tmp}}
		else {
			if(nrow(attr(mu1,"model"))!=n)stop("mu model matrix does not match number of response observations")
			if(is.function(mu)){
				dm1 <- attr(mu1,"model")
				lin1 <- mu1
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
			else pmu <- unlist(pmu)}}}
else if(!is.function(mu)){
	mu1 <- function(p) p[1]*rep(1,n)
	npt1 <- 1}
else {
	mu1 <- mu
	if(length(mu1(pmu))==1)mu1 <- function(p) mu(p)*rep(1,n)}
if(is.null(attributes(mu1))){
	attributes(mu1) <- if(is.function(mu)){
		if(!inherits(mu,"formulafn"))attributes(fnenvir(mu))
		else attributes(mu)}
		else attributes(fnenvir(mu1))}
nlp <- if(is.function(mu)){
		if(is.null(lin1))length(attr(mu1,"parameters"))
		else length(attr(mu1,"parameters"))-1+npt1}
       else npt1
if(nlp!=npl)stop(paste("pmu should have",nlp,"initial estimates"))
if(any(is.na(mu1(pmu))))stop("The location model returns NAs: probably invalid initial values")
npl1 <- npl+1
if(inherits(lin2,"formula")){
	sh1 <- finterp(lin2,envir=envir,start=npl1,name=name)
	npt2 <- length(attr(sh1,"parameters"))
	if(is.matrix(attr(sh1,"model"))){
		if(all(dim(attr(sh1,"model"))==1)){
			if(is.function(shape)){
				lin2 <- sh1
				sh1 <- if(shfn)function(p) shape(p[npl1:np],p[npl1]*rep(1,n), mu1(p))
					else function(p) shape(p[npl1:np],p[npl1])}
			else {
				tmp <- attributes(sh1)
				sh1 <- function(p) p[npl1]*rep(1,n)
				sh2 <- fnenvir(function(p) p[1]*rep(1,n))
				attributes(sh1) <- tmp}}
		else {
			if(nrow(attr(sh1,"model"))!=n)stop("shape model matrix does not match number of response observations")
			if(is.function(shape)){
				dm2 <- attr(sh1,"model")
				lin2 <- sh1
				sh1 <- if(shfn)function(p) shape(p[npl1:np], dm2%*%p[npl1:(npl1+npt2-1)], mu1(p))
					else function(p) shape(p[npl1:np],dm2%*%p[npl1:(npl1+npt2-1)])}}}
	else {
		if(nps!=npt2){
			cat("\nParameters are ")
			cat(attr(sh1,"parameters"),"\n")
			stop(paste("pshape should have",npt2,"estimates"))}
		if(is.list(pshape)){
			if(!is.null(names(pshape))){
				o <- match(attr(sh1,"parameters"),names(pshape))
				pshape <- unlist(pshape)[o]
				if(sum(!is.na(o))!=length(pshape))stop("invalid estimates for shape - probably wrong names")}
			else pshape <- unlist(pshape)}}}
else if(!is.function(shape)&&shp){
	sh1 <- function(p) p[npl1]*rep(1,n)
	sh2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt2 <- 1}
else {
	sh1 <- if(shfn)function(p) shape(p[npl1:np], mu1(p))
		else function(p) shape(p[npl1:np])}
if(shp){
	if(is.null(attributes(sh1))){
		attributes(sh1) <- if(is.function(shape)){
			if(!inherits(shape,"formulafn"))attributes(fnenvir(shape))
			else attributes(shape)}
			else attributes(fnenvir(sh1))}
	nlp <- if(is.function(shape)){
			if(is.null(lin2))length(attr(sh1,"parameters"))-shfn
			else length(attr(sh1,"parameters"))-1+npt2-shfn}
		else npt2
	if(nlp!=nps)stop(paste("pshape should have",nlp,"initial estimates"))}
p <- c(pmu,pshape)
if(distribution=="Levy"&&((!censor&&any(y<=mu1(p)))||(censor&&any(y[,1]<=mu1(p)))))
	stop("location parameter must be strictly less than corresponding observation")
if(distribution=="Pareto"&&exp(sh1(p))<=1)stop("shape parameters must be > 0")
if(distribution!="binomial"&&distribution!="Poisson"&&
	distribution!="exponential"&&distribution!="geometric"&&
	distribution!="logarithmic"&&any(is.na(sh1(p))))
	stop("The shape model returns NAs: probably invalid initial values")
if (!censor){
	ret <- switch(distribution,
	binomial={
		fcn <- function(p) {
			m <- mu1(p)
			-wt*(y[,1]*log(m)+y[,2]*log(1-m))}
		const <- -wt*lchoose(nn,y[,1])},
	"beta binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- exp(sh1(p))
			t <- s*m
			u <- s*(1-m)
			-wt*(lbeta(y[,1]+t,y[,2]+u)-lbeta(t,u))}
		const <- -wt*lchoose(nn,y[,1])},
	"double binomial"={
		fcn <- function(p) {
			-.C("ddb",as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(n),as.double(wt),res=double(n),DUP=F)$res}
		const <- 0},
	"mult binomial"={
		fcn <- function(p) {
			-.C("dmb",as.integer(y[,1]),as.integer(nn),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(n),as.double(wt),res=double(n),DUP=F)$res}
		const <- 0},
	Poisson={
		fcn <- function(p) {
			m <- mu1(p)
			-wt*(-m+y*log(m))}
		const <- wt*lgamma(y+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(lgamma(y+s)-lgamma(s)+s*t+y*log(m)
				-(y+s)*log(s+m))}
		const <- wt*lgamma(y+1)},
	"double Poisson"={
		fcn <- function(p) {
			-.C("ddp",as.integer(y),as.integer(my),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(length(y)),as.double(wt),
				res=double(length(y)),DUP=F)$res}
		const <- 0},
	"mult Poisson"={
		fcn <- function(p) {
			-.C("dmp",as.integer(y),as.integer(my),
				as.double(mu1(p)),as.double(exp(sh1(p))),
				as.integer(length(y)),as.double(wt),
				res=double(length(y)),DUP=F)$res}
		const <- 0},
	"gamma count"={
		fcn <- function(p) {
			m <- mu1(p)
			s <- exp(sh1(p))
			u <- m*s
			-wt*log(ifelse(y==0,1-pgamma(u,(y+1)*s,1),
				pgamma(u,y*s+(y==0),1)-
				pgamma(u,(y+1)*s,1)))}
		const <- 0},
	Consul={
		fcn <- function(p) {
			m <- mu1(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(log(m)-(m+y*(s-1))/s+(y-1)*log(m+y*(s-1))-y*t)}
		const <- wt*lgamma(y+1)},
	logarithmic={
		fcn <- function(p) {
			m <- exp(mu1(p))
			m <- m/(1+m)
			-wt*(y*log(m)-log(y)-log(-log(1-m)))}
		const <- 0},
	geometric={
		fcn <- function(p) {
			m <- mu1(p)
			-wt*(y*log(m)-(y+1)*log(1+m))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				-wt*log(pnorm(y+delta/2,m,s)
					-pnorm(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				wt*(t+(y-mu1(p))^2/exp(t))/2}
			const <- wt*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*log(pinvgauss(y+delta/2,m,s)-
					pinvgauss(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				wt*(t+(y-m)^2/(y*exp(t)*m^2))/2}
			const <- wt*(log(2*pi*y^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				-wt*log(plogis(y+delta/2,m,s)
					-plogis(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- (y-m)*pi/(exp(t)*sqrt(3))
				wt*(s+t+2*log(1+exp(-s)))}
			const <- -wt*(log(pi/sqrt(3))+log(delta))}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				-wt*log(pcauchy(y+delta/2,m,s)
					-pcauchy(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				wt*log(s*(1+((y-m)/s)^2))}
			const <- -wt*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*log(plaplace(y+delta/2,m,s)
					-plaplace(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				wt*(abs(y-mu1(p))/exp(t)+t)}
			const <- -wt*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*log(plevy(y+delta/2,m,s)
					-plevy(y-delta/2,m,s))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(0.5*log(s/(2*pi))-1.5*log(y-m)-
					s/(2*(y-m)))}
			const <- -wt*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*log((1+(y-delta/2)*t)^-s
					-(1+(y+delta/2)*t)^-s)}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*(log(s*t)-(s+1)*log(1+y*t))}
			const <- -wt*log(delta)}},
        exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				-wt*log(-exp(-(y+delta/2)/m)
					+exp(-(y-delta/2)/m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				wt*(log(m)+y/m)}
			const <- -wt*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				u <- m/s
				-wt*log(pgamma(y+delta/2,s,u)
					-pgamma(y-delta/2,s,u))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(s*(t-log(m)-y/m)+(s-1)*log(y)-lgamma(s))}
			const <- -wt*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*log(pweibull(y+delta/2,s,m)
					-pweibull(y-delta/2,s,m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(t+(s-1)*log(y)-s*log(m)-(y/m)^s)}
			const <- -wt*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				ey <- exp(y[,1])
				-wt*log(pweibull(ey+ey*delta/2,s,m)
					-pweibull(ey-ey*delta/2,s,m))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(t+s*y-s*log(m)-(exp(y)/m)^s)}
			const <- -wt*log(delta)}},
	own={ const <- 0})}
else {
	ret <- switch(distribution,
	Poisson={
		fcn <- function(p) {
			m <- mu1(p)
			-wt*(cc*(-m+y[,1]*log(m))+
				log(lc-rc*ppois(y[,1],m)))}
		const <- wt*cc*lgamma(y[,1]+1)},
	"negative binomial"={
		fcn <- function(p) {
			m <- mu1(p)
			t <- sh1(p)
			s <- exp(t)
			-wt*(cc*(lgamma(y[,1]+s)-lgamma(s)
				+s*t+y[,1]*log(m)-(y[,1]+s)*log(s+m))+
				log(lc-rc*pnbinom(y[,1],s,1/(1+m/s))))}
		const <- wt*cc*lgamma(y[,1]+1)},
	geometric={
		fcn <- function(p) {
			m <- mu1(p)
			-wt*(cc*(y[,1]*log(m)-(y[,1]+1)*log(1+m))+
				log(lc-rc*pgeom(y[,1],1/(1+m))))}
		const <- 0},
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				-wt*(cc*log(pnorm(y[,1]+delta/2,m,s)-
					pnorm(y[,1]-delta/2,m,s))
					+log(lc-rc*pnorm(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-(t+(y[,1]-m)^2/s)/2)+log(lc-rc
					*pnorm(y[,1],m,sqrt(s))))}
			const <- wt*cc*(log(2*pi)/2-log(delta))}},
        "inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(cc*log(pinvgauss(y[,1]+delta/2,m,s)-
					pinvgauss(y[,1]-delta/2,m,s))
					+log(lc-rc*pinvgauss(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-(t+(y[,1]-m)^2/(y[,1]*s*m^2))/2)+
					log(lc-rc*pinvgauss(y[,1],m,s)))}
			const <- wt*cc*(log(2*pi*y[,1]^3)/2-log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				-wt*(cc*log(plogis(y[,1]+delta/2,m,s)-
					plogis(y[,1]-delta/2,m,s))
					+log(lc-rc*plogis(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				y1 <- (y[,1]-m)/s
				-wt*(cc*(-y1-log(s)-2*log(1+exp(-y1)))
					+log(lc-rc*plogis(y[,1],m,s)))}
			const <- -wt*cc*log(delta)}},
	Cauchy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				-wt*(cc*log(pcauchy(y[,1]+delta/2,m,s)-
					pcauchy(y[,1]-delta/2,m,s))
					+log(lc-rc*pcauchy(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				-wt*(-cc*log(s*(1+((y[,1]-m)/s)^2))
					+log(lc-rc*pcauchy(y[,1],m,s)))}
			const <- -wt*cc*log(delta/pi)}},
        Laplace={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(cc*log(plaplace(y[,1]+delta/2,m,s)-
					plaplace(y[,1]-delta/2,m,s))
					+log(lc-rc*plaplace(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(-abs(y[,1]-m)/s-t)+log(lc-rc
					*plaplace(y[,1],m,s)))}
			const <- -wt*cc*log(delta/2)}},
        Levy={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(cc*log(plevy(y[,1]+delta/2,m,s)-
					plevy(y[,1]-delta/2,m,s))
					+log(lc-rc*plevy(y[,1],m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(0.5*log(s/(2*pi))-1.5*log(y[,1]-m)-
					s/(2*(y[,1]-m)))+log(lc-rc
					*plevy(y[,1],m,s)))}
			const <- -wt*cc*log(delta/2)}},
        Pareto={
		if(exact){
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*(cc*log((1+(y[,1]-delta/2)*t)^-s-
					(1+(y[,1]+delta/2)*t)^-s)
					+log(lc-rc*(-(1+(y[,1])*t)^-s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(sh1(p))
				t <- 1/(mu1(p)*(s-1))
				-wt*(cc*(log(s*t)-(s+1)*log(1+y[,1]*t))
					+log(lc-rc*(1-(1+y[,1]*t)^-s)))}
			const <- -wt*cc*log(delta)}},
	exponential={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				-wt*(cc*log(-exp(-(y[,1]+delta/2)/m)
					+exp(-(y[,1]-delta/2)/m))+
					log(lc-rc*(1-exp(-y[,1]/m))))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				-wt*(cc*(-log(m)-y[,1]/m)+log(lc-rc*
					(1-exp(-y[,1]/m))))}
			const <- -wt*cc*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				u <- m/s
				-wt*(cc*log(pgamma(y[,1]+delta/2,s,u)-
					pgamma(y[,1]-delta/2,s,u))
					+log(lc-rc*pgamma(y[,1],s,u)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(s*(t-log(m)-y[,1]/m)+(s-1)*log(y[,1])
					-lgamma(s))+log(lc-rc
					*pgamma(y[,1],s,m/s)))}
			const <- -wt*cc*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				-wt*(cc*log(pweibull(y[,1]+delta/2,s,m)-
					pweibull(y[,1]-delta/2,s,m))
					+log(lc-rc*pweibull(y[,1],s,m)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				-wt*(cc*(t+(s-1)*log(y[,1])-s*log(m)
					-(y[,1]/m)^s)+log(lc-rc*
					pweibull(y[,1],s,m)))}
			const <- -wt*cc*log(delta)}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				ey <- exp(y[,1])
				pw <- pweibull(ey-ey*delta/2,s,m)
				-wt*(cc*log(pweibull(ey+ey*delta/2,s,m)-
					pweibull(ey-ey*delta/2,s,m))
					+log(lc-rc*pweibull(ey,s,m)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				ey <- exp(y[,1])
				-wt*(cc*(t+s*y[,1]-s*log(m)-(ey/m)^s)+log(lc-
					rc*pweibull(ey,s,m)))}
			const <- -wt*cc*log(delta)}},
	own={const <- 0})}
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
residuals <- if(distribution=="binomial"||distribution=="beta binomial"||
	distribution=="double binomial"||distribution=="mult binomial"||censor)
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
if(!is.null(sh2))sh1 <- sh2
if(!is.null(lin1a))lin1 <- lin1a
if(!is.null(lin2a))lin2 <- lin2a
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	mu=mu1,
	shape=sh1,
	linear=list(lin1,lin2),
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
	npm=0,
	nps=nps,
	npf=0,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlr"
return(z1)}
