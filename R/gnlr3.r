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
#     gnlr3(y, distribution="normal", mu=NULL, shape=NULL,
#	family=NULL, linear=NULL, pmu=NULL, pshape=NULL,
#	pfamily=NULL, exact=F, wt=1, delta=1, envir=sys.frame(sys.parent()),
#	print.level=0,typsiz=abs(p), ndigit=10, gradtol=0.00001,
#	stepmax=10*sqrt(p%*%p), steptol=0.00001, iterlim=100, fscale=1)
#
#  DESCRIPTION
#
#    A function to fit nonlinear regression models with a variety of
# three parameter distributions.

require(rmutil)

gnlr3 <- function(y, distribution="normal", mu=NULL, shape=NULL,
	family=NULL, linear=NULL, pmu=NULL, pshape=NULL,
	pfamily=NULL, exact=F, wt=1, delta=1, envir=sys.frame(sys.parent()),
	print.level=0, typsiz=abs(p), ndigit=10, gradtol=0.00001,
	stepmax=10*sqrt(p%*%p), steptol=0.00001, iterlim=100, fscale=1){

pburr <- function(q, m, s, f) 1-(1+(q/m)^s/f)^-f
pglogis <- function(y, m, s, f) (1+exp(-sqrt(3)*(y-m)/(s*pi)))^-f
pgweibull <- function(y, s, m, f) (1-exp(-(y/m)^s))^f
phjorth <- function(y, m, s, f) 1-(1+s*y)^(-f/s)*exp(-(y/m)^2/2)
pginvgauss <- function(y, m, s, f)
	.C("pginvgauss",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(f),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		DUP=F)$res
ppowexp <- function(y, m, s, f){
	z <- .C("ppowexp",
		as.double(y),
		as.double(m),
		as.double(s),
		as.double(f),
		len=as.integer(n),
		eps=as.double(1.0e-6),
		pts=as.integer(5),
		max=as.integer(16),
		err=integer(1),
		res=double(n),
		DUP=F)$res
	ifelse(y-m>0,0.5+z,0.5-z)}

call <- sys.call()
if(!missing(distribution)&&!is.function(distribution)){
	distribution <- match.arg(distribution,c("normal","inverse Gauss",
	"logistic","Hjorth","gamma","Burr","Weibull","extreme value",
	"Student t","power exponential"))}
if(!is.null(pmu))npl <- length(pmu)
else npl <- 0
if(!is.null(pshape))nps <- length(pshape)
else nps <- 0
if(!is.null(pfamily))npf <- length(pfamily)
else npf <- 0
np <- npl+nps+npf
if(np<1)stop("At least one parameter must be estimated")
if(is.function(distribution)){
	fcn <- distribution
	distribution <- "own"}
if(inherits(y,"repeated")){
	if(inherits(envir,"repeated")&&(length(y$response$nobs)!=length(envir$response$nobs)||any(y$response$nobs!=envir$response$nobs)))stop("y and envir objects are incompatible")
	if(!is.null(y$response$wt)&&!is.na(y$response$wt))wt <- y$response$wt
	if(!is.null(y$response$delta))delta <- y$response$delta
	if(!is.null(y$response$censor)&&any(y$response$censor!=1))y <- cbind(y$response$y,y$response$censor)
	else y <- y$response$y}
else if(inherits(y,"response")){
	if(!is.null(y$wt)&&!is.na(y$wt))wt <- y$wt
	if(!is.null(y$delta))delta <- y$delta
	if(!is.null(y$censor)&&any(y$censor!=1))y <- cbind(y$y,y$censor)
	else y <- y$y}
censor <- length(dim(y))==2&&ncol(y)==2
if(censor&&all(y[,2]==1)){
	y <- y[,1]
	censor <- F}
if(censor){
	n <- nrow(y)
	y[,2] <- as.integer(y[,2])
	if(any(y[,2]!=-1&y[,2]!=0&y[,2]!=1))
		stop("Censor indicator must be -1s, 0s, and 1s")
	cc <- ifelse(y[,2]==1,1,0)
	rc <- ifelse(y[,2]==0,1,ifelse(y[,2]==-1,-1,0))
	lc <- ifelse(y[,2]==-1,0,1)
	if(delta<=0&y[,2]==1)
		stop("All deltas for uncensored data must be positive")
	else {
		delta <- ifelse(delta<=0,0.000001,delta)
		delta <- ifelse(y[,1]-delta/2<=0,delta-0.00001,delta)}}
else {
	if(!is.vector(y,mode="numeric"))stop("y must be a vector")
	n <- length(y)
	if(min(delta)<=0)stop("All deltas for must be positive")}
if((distribution!="logistic"&&distribution!="Student t"&&
	distribution!="power exponential")&&((censor&&any(y[,1]<=0))||
	(!censor&&any(y<=0))))stop("All response values must be > 0")
if(min(wt)<0)stop("All weights must be non-negative")
if(length(wt)==1)wt <- rep(wt,n)
if(length(delta)==1)delta <- rep(delta,n)
lin1 <- lin2 <- lin3 <- NULL
if(is.list(linear)){
	lin1 <- linear[[1]]
	lin2 <- linear[[2]]
	lin3 <- linear[[3]]}
else lin1 <- linear
if(inherits(mu,"formula"))lin1 <- mu
if(inherits(shape,"formula"))lin2 <- shape
if(inherits(family,"formula"))lin3 <- family
lin1a <- lin2a <- lin3a <- mu2 <- sh2 <- fa2 <- name <- NULL
if(inherits(envir,"repeated")||inherits(envir,"tccov")){
	type <- if(inherits(envir,"repeated"))"repeated"
		else "tccov"
	if(inherits(lin1,"formula")){
		lin1a <- finterp(lin1)
		class(lin1) <- c(class(lin1),type)}
	if(inherits(lin2,"formula")){
		lin2a <- finterp(lin2)
		class(lin2) <- c(class(lin2),type)}
	if(inherits(lin[[3]],"formula")){
		lin3a <- finterp(lin3)
		class(lin3) <- c(class(lin3),type)}
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
		shape <- fnenvir(shape,envir=envir,name=name)}
	if(is.function(shape)){
		fa2 <- shape
		attributes(fa2) <- attributes(fnenvir(shape))
		class(family) <- type
		family <- fnenvir(family,envir=envir,name=name)}}
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
				sh1 <- function(p) shape(p[npl1:np],p[npl1]*rep(1,n))}
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
				sh1 <- function(p) shape(p[npl1:np],dm2%*%p[npl1:(npl1+npt2-1)])}}}
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
else if(!is.function(shape)){
	sh1 <- function(p) p[npl1]*rep(1,n)
	sh2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt2 <- 1}
else sh1 <- function(p) shape(p[npl1:np])
if(is.null(attributes(sh1))){
	attributes(sh1) <- if(is.function(shape)){
		if(!inherits(shape,"formulafn"))attributes(fnenvir(shape))
		else attributes(shape)}
		else attributes(fnenvir(sh1))}
nlp <- if(is.function(shape)){
		if(is.null(lin2))length(attr(sh1,"parameters"))
		else length(attr(sh1,"parameters"))-1+npt2}
       else npt2
if(nlp!=nps)stop(paste("pshape should have",nlp,"initial estimates"))
np1 <- npl+nps
nps1 <- np1+1
if(inherits(lin3,"formula")){
	fa1 <- finterp(lin3,envir=envir,start=nps1,name=name)
	npt3 <- length(attr(fa1,"parameters"))
	if(is.matrix(attr(fa1,"model"))){
		if(all(dim(attr(fa1,"model"))==1)){
			if(is.function(family)){
				lin3 <- fa1
				fa1 <- function(p) family(p[nps1:np],p[nps1]*rep(1,n))}
			else {
				tmp <- attributes(fa1)
				fa1 <- function(p) p[nps1]*rep(1,n)
				fa2 <- fnenvir(function(p) p[1]*rep(1,n))
				attributes(fa1) <- tmp}}
		else {
			if(nrow(attr(fa1,"model"))!=n)stop("family model matrix does not match number of response observations")
			if(is.function(family)){
				dm3 <- attr(fa1,"model")
				lin3 <- fa1
				fa1 <- function(p) family(p[nps1:np],dm3%*%p[nps1:(nps1+npt3-1)])}}}
	else {
		if(npf!=npt3){
			cat("\nParameters are ")
			cat(attr(fa1,"parameters"),"\n")
			stop(paste("pfamily should have",npt3,"estimates"))}
		if(is.list(pfamily)){
			if(!is.null(names(pfamily))){
				o <- match(attr(fa1,"parameters"),names(pfamily))
				pfamily <- unlist(pfamily)[o]
				if(sum(!is.na(o))!=length(pfamily))stop("invalid estimates for family - probably wrong names")}
			else pfamily <- unlist(pfamily)}}
	if(npf<npt3)stop("Not enough initial estimates for family")}
else if(!is.function(family)){
	fa1 <- function(p) p[nps1]*rep(1,n)
	fa2 <- fnenvir(function(p) p[1]*rep(1,n))
	npt3 <- 1}
else fa1 <- function(p) family(p[nps1:np])
if(is.null(attributes(fa1))){
	attributes(fa1) <- if(is.function(family)){
		if(!inherits(family,"formulafn"))attributes(fnenvir(family))
		else attributes(family)}
		else attributes(fnenvir(fa1))}
nlp <- if(is.function(family)){
		if(is.null(lin3))length(attr(fa1,"parameters"))
		else length(attr(fa1,"parameters"))-1+npt3}
       else npt3
if(nlp!=npf)stop(paste("pfamily should have",nlp,"initial estimates"))
p <- c(pmu,pshape,pfamily)
if(any(is.na(sh1(p))))stop("The shape model returns NAs: probably invalid initial values")
if(any(is.na(fa1(p))))stop("The family model returns NAs: probably invalid initial values")
if (!censor){
	ret <- switch(distribution,
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				y <- y^f/f
				jy <- y^(2*f-1)*delta/(2*f)
				norm <- sign(f)*pnorm(0,m,s)
				ind <- f<0
				-wt*(log((pnorm(y+jy,m,s)-pnorm(y-jy,m,s)))
					-log(1-ind-norm))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				norm <- sign(f)*pnorm(0,m,s)
				ind <- f<0
				-wt*((f-1)*log(y)+log(dnorm(y^f/f,m,s))
					-log(1-ind-norm))}
			const <- -wt*log(delta)}},
	"power exponential"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*log(ppowexp(y+delta/2,m,s)
					-ppowexp(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- 0.5*sh1(p)
				f <- exp(fa1(p))
				b <- 1+1/(2*f)
				wt*(t+(abs(y-mu1(p))/exp(t))^(2*f)/2+
					lgamma(b)+b*log(2))}
			const <- -wt*log(delta)}},
	"inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*log(pginvgauss(y+delta/2,m,s,f)
					-pginvgauss(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(-f*log(m)+(f-1)*log(y)-
					log(2*besselK(1/(s*m),abs(f)))-
					(1/y+y/m^2)/(2*s))}
			const <- -wt*log(delta)}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pglogis(y+delta/2,m,s,f)
					-pglogis(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				t <- sh1(p)
				m <- (y-mu1(p))/exp(t)*sqrt(3)/pi
				wt*(-fa1(p)+m+t+(exp(fa1(p))+1)*
					log(1+exp(-m)))}
			const <- -wt*(log(delta*sqrt(3)/pi))}},
	Hjorth={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*log(phjorth(y+delta/2,m,s,f)-
					phjorth(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(-f*log(1+s*y)/s-(y/m)^2/2+
					log(y/m^2+f/(1+s*y)))}
			const <- -wt*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				u <- (m/s)^f
				-wt*log(pgamma((y+delta/2)^f,s,u)
					-pgamma((y-delta/2)^f,s,u))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				v <- s*f
				-wt*(v*(t-log(m))-(s*y/m)^f+u+(v-1)*log(y)
					-lgamma(s))}
			const <- -wt*log(delta)}},
	Burr={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pburr(y+delta/2,m,s,f)-
					pburr(y-delta/2,m,s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				y1 <- y/m
				-wt*(log(s/m)+(s-1)*log(y1)
					-(f+1)*log(1+y1^s/f))}
			const <- -wt*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*log(pgweibull(y+delta/2,s,m,f)
					-pgweibull(y-delta/2,s,m,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				y1 <- (y/m)^s
				-wt*(t+u+(s-1)*log(y)-s*log(m)+
					(f-1)*log(1-exp(-y1))-y1)}
			const <- -wt*log(delta)}},
        "Student t"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*log(pt((y+delta/2-m)/s,f)-
					pt((y-delta/2-m)/s,f))}
			const <- 0}
		else {
			fcn <- function(p) {
				s <- exp(0.5*sh1(p))
				-wt*log(dt((y-mu1(p))/s,exp(fa1(p)))/s)}
			const <- -wt*(log(delta))}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				y1 <- y^f/f
				ey <- exp(y1)
				jey <- y^(f-1)*ey*delta/2
				norm <- sign(f)*exp(-m^-s)
				ind <- f>0
				-wt*(log((pweibull(ey+jey,s,m)
					-pweibull(ey-jey,s,m))/
					(1-ind+norm)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				f <- fa1(p)
				y1 <- y^f/f
				norm <- sign(f)*exp(-m^-s)
				ind <- f>0
				-wt*(t+s*(y1-log(m))-(exp(y1)/m)^s
					+(f-1)*log(y)-
					log(1-ind+norm))}
			const <- -wt*log(delta)}},
	own={ const <- 0})}
else {
	ret <- switch(distribution,
	normal={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				y <- y[,1]^f/f
				jy <- y[,1]^(2*f-1)*delta/(2*f)
				norm <- sign(f)*pnorm(0,m,s)
				ind <- f<0
				-wt*(cc*log((pnorm(y+jy,m,s)-pnorm(y-jy,m,s)))+
					log(lc-rc*(pnorm(y,m,s)-(f>0)*norm)))/
					(1-ind-norm)}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(0.5*t)
				f <- fa1(p)
				norm <- sign(f)*pnorm(0,m,s)
				ind <- f<0
				-wt*(cc*(-(t+((y[,1]^f/f-m)/s)^2)/2+(f-1)*
					log(y[,1]))+log(lc-rc
					*(pnorm(y[,1]^f/f,m,s)
					-(f>0)*norm)))/(1-ind-norm)}
			const <- wt*cc*(log(2*pi)/2-log(delta))}},
	"power exponential"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- fa1(p)
				-wt*(cc*log(ppowexp(y[,1]+delta/2,m,s,f)-
					ppowexp(y[,1]-delta/2,m,s,f))
					+log(lc-rc*ppowexp(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- 0.5*sh1(p)
				s <- exp(t)
				f <- exp(fa1(p))
				b <- 1+1/(2*f)
				-wt*(cc*(-t-(abs(y[,1]-mu1(p))/s)^(2*f)/2-
					lgamma(b)-b*log(2))+log(lc-rc
					*ppowexp(y[,1],m,s,f)))}
			const <- -wt*cc*(log(delta))}},
	"inverse Gauss"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p)/2)
				f <- fa1(p)
				-wt*(cc*log(pginvgauss(y[,1]+delta/2,m,s,f)-
					pginvgauss((y[,1]-delta/2)^f/f,m,s))+
					log(lc-rc*pginvgauss(y[,1]^f/f,m,s)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(cc*(-f*log(m)+(f-1)*log(y[,1])-
					log(2*besselK(1/(s*m),abs(f)))-
					(1/y[,1]+y[,1]/m^2)/(2*s))+log(lc-rc
					*pginvgauss(y[,1],m,s,f)))}
			const <- -wt*cc*(log(delta))}},
	logistic={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				f <- exp(fa1(p))
				-wt*(cc*log(pglogis(y[,1]+delta/2,m,s,f)-
					pglogis(y[,1]-delta/2,m,s,f))
					+log(lc-rc*pglogis(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))*sqrt(3)/pi
				y1 <- (y[,1]-m)/s
				u <- fa1(p)
				f <- exp(u)
				-wt*(cc*(u-y1-log(s)-(f+1)*log(1+exp(-y1)))
					+log(lc-rc*pglogis(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}},
	Hjorth={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(cc*log(phjorth(y[,1]+delta/2,m,s,f)-
					phjorth(y[,1]-delta/2,m,s,f))
					+log(lc-rc*phjorth(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				-wt*(cc*(-f*log(1+s*y[,1])/s-(y[,1]/m)^2/2+
					log(y[,1]/m^2+f/(1+s*y[,1])))+
					log(lc-rc*phjorth(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}},
        gamma={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				u <- (m/s)^f
				-wt*(cc*log(pgamma((y[,1]+delta/2)^f,s,u)-
					pgamma((y[,1]-delta/2)^f,s,u))
					+log(lc-rc*pgamma(y[,1]^f,s,u)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				v <- s*f
				-wt*(cc*(v*(t-log(m))-(s*y[,1]/m)^f+u+(v-1)*
					log(y[,1])-lgamma(s))+log(lc-rc
					*pgamma(y[,1]^f,s,(m/s)^f)))}
			const <- -wt*cc*log(delta)}},
	Burr={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pburr(y[,1]+delta/2,m,s,f)-
					pburr(y[,1]-delta/2,m,s,f))
					+log(lc-rc*pburr(y[,1],m,s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				y1 <- y[,1]/m
				-wt*(cc*(log(s/m)+(s-1)*log(y1)
					-(f+1)*log(1+y1^s/f))+
					log(lc-rc*pburr(y[,1],m,s,f)))}
			const <- -wt*cc*log(delta)}},
        Weibull={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pgweibull(y[,1]+delta/2,s,m,f)-
					pgweibull(y[,1]-delta/2,s,m,f))
					+log(lc-rc*pgweibull(y[,1],s,m,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				u <- fa1(p)
				f <- exp(u)
				y1 <- (y[,1]/m)^s
				-wt*(cc*(t+u+(s-1)*log(y[,1])-s*log(m)+(f-1)
					*log(1-exp(-y1))-y1)+log(lc-rc*
					pgweibull(y[,1],s,m,f)))}
			const <- -wt*cc*log(delta)}},
        "Student t"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(pt((y[,1]+delta/2-m)/s,f)-
					pt((y[,1]-delta/2-m)/s,f))
					+log(lc-rc*pt((y[,1]-m)/s,f)))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(0.5*sh1(p))
				f <- exp(fa1(p))
				-wt*(cc*log(dt((y[,1]-m)/s,f)/s)
					+log(lc-rc*pt((y[,1]-m)/s,f)))}
			const <- -wt*cc*(log(delta))}},
        "extreme value"={
		if(exact){
			fcn <- function(p) {
				m <- mu1(p)
				s <- exp(sh1(p))
				f <- fa1(p)
				y1 <- y[,1]^f/f
				ey <- exp(y1)
				jey <- y[,1]^(f-1)*ey*delta/2
				norm <- sign(f)*exp(-m^-s)
				ind <- f>0
				-wt*(cc*log(pweibull(ey+jey,s,m)-
					pweibull(ey-jey,s,m))
					+log(lc-rc*(pweibull(ey,s,m)-ind+
					(f>0)*norm))-log(1-ind+norm))}
			const <- 0}
		else {
			fcn <- function(p) {
				m <- mu1(p)
				t <- sh1(p)
				s <- exp(t)
				f <- fa1(p)
				y1 <- y[,1]^f/f
				ey <- exp(y1)
				norm <- sign(f)*exp(-m^-s)
				ind <- f>0
				-wt*(cc*(t+s*(y1-log(m))-(ey/m)^s
					+(f-1)*log(y[,1]))+log(lc-rc*
					(pweibull(ey,s,m)-ind+(f>0)*norm))-
					log(1-ind+norm))}
			const <- -wt*cc*log(delta)}},
	own={const <- 0})}
fn <- function(p) sum(fcn(p))
if(fscale==1)fscale <- fn(p)
if(is.na(fn(p)))stop("Likelihood returns NAs: probably invalid initial values")
z0 <- nlm(fn, p=p, hessian=T, print.level=print.level, typsiz=typsiz,
	ndigit=ndigit, gradtol=gradtol, stepmax=stepmax, steptol=steptol,
	iterlim=iterlim, fscale=fscale)
z0$minimum <- z0$minimum+sum(const)
fitted.values <- as.vector(mu1(z0$estimate))
residuals <- y-fitted.values
if(np==1){
	cov <- 1/z0$hessian
	se <- as.vector(sqrt(cov))}
else {
	a <- qr(z0$hessian)
	if(a$rank==np)cov <- solve(z0$hessian)
	else cov <- matrix(NA,ncol=np,nrow=np)
	se <- sqrt(diag(cov))}
like.comp <- as.vector(fcn(z0$estimate)+const)
if(!is.null(mu2))mu1 <- mu2
if(!is.null(sh2))sh1 <- sh2
if(!is.null(fa2))fa1 <- fa2
if(!is.null(lin1a))lin1 <- lin1a
if(!is.null(lin2a))lin2 <- lin2a
if(!is.null(lin3a))lin3 <- lin3a
z1 <- list(
	call=call,
	delta=delta,
	distribution=distribution,
	likefn=fcn,
	mu=mu1,
	shape=sh1,
	family=fa1,
	linear=list(lin1,lin2,lin3),
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
	npf=npf,
	se=se,
	cov=cov,
	corr=cov/(se%o%se),
	gradient=z0$gradient,
	iterations=z0$iterations,
	code=z0$code)
class(z1) <- "gnlr"
return(z1)}
