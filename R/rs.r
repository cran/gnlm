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
#     rs2(y, x1, x2, power=c(1,1), weight=rep(1,length(x1)),
#	family=normal, iterlim=20)
#     rs3(y, x1, x2, x3, power=c(1,1,1), weight=rep(1,length(x1)),
#	family=normal, iterlim=20)
#
#  DESCRIPTION
#
#    A function to fit two- and three-covariate power-transformed
# response surface models for glms (Box-Tidwell transformation)

rs2 <- function(y, x1, x2, power=c(1,1), weight=rep(1,length(x1)),
	family=normal, iterlim=20){
	if(length(power)!=2)
		stop("Two estimates of power parameters must be supplied\n")
	if(any(c(x1,x2)<0))stop("All covariates must be non-negative")
	a <- power[1]
	b <- power[2]
	test <- T
	i <- 0
	while(test){
		xx1 <- x1^a
		xx2 <- x2^b
		u <- glm(y~xx1+xx2+I(xx1^2)+I(xx2^2)+xx1:xx2,family=family,
			weight=weight)
		z1 <- (u$coef[2]*xx1+2*u$coef[4]*xx1^2+u$coef[6]*xx1*xx2)*
			log(ifelse(x1==0,1,x1))
		z2 <- (u$coef[3]*xx2+2*u$coef[5]*xx2^2+u$coef[6]*xx1*xx2)*
			log(ifelse(x2==0,1,x2))
		if(any(is.na(c(z1,z2))))stop(paste("NAs in calculating estimates:",a,b))
		u <- glm(y~xx1+xx2+I(xx1^2)+I(xx2^2)+xx1:xx2+z1+z2,
			family=family,weight=weight)
		a <- a+u$coef[6]
		b <- b+u$coef[7]
		if(any(is.na(c(a,b))))stop(paste("NAs in calculating estimates:",a,b))
		i <- i+1
		test <- ((u$coef[6]^2>0.00001)||(u$coef[7]^2>0.00001))&&(i<iterlim)}
	z <- glm(y~xx1+xx2+I(xx1^2)+I(xx2^2)+xx1:xx2,family=family,
		weight=weight)
	z$df.residual <- z$df.residual-2
	z$aic <- z$aic+4
	z$powers <- c(a,b)
	z$iterations <- i
	class(z) <- c("rs",class(z))
	return(z)}

rs3 <- function(y, x1, x2, x3, power=c(1,1,1), weight=rep(1,length(x1)),
	family=normal, iterlim=20){
	if(length(power)!=3)
		stop("Three estimates of power parameters must be supplied\n")
	if(any(c(x1,x2,x3)<0))stop("All covariates must be non-negative")
	a <- power[1]
	b <- power[2]
	d <- power[3]
	test <- T
	i <- 0
	while(test){
		xx1 <- x1^a
		xx2 <- x2^b
		xx3 <- x3^d
		xx12 <- xx1*xx2
		xx13 <- xx1*xx3
		xx23 <- xx2*xx3
		u <- glm(y~xx1+xx2+xx3+I(xx1^2)+I(xx2^2)+I(xx3^2)+
			xx12+xx13+xx23,family=family,weight=weight)
		z1 <- (u$coef[2]*xx1+2*u$coef[5]*xx1^2+u$coef[8]*xx12+
			+u$coef[9]*xx13)*log(ifelse(x1==0,1,x1))
		z2 <- (u$coef[3]*xx2+2*u$coef[6]*xx2^2+u$coef[8]*xx12+
			u$coef[10]*xx23)*log(ifelse(x2==0,1,x2))
		z3 <- (u$coef[4]*xx2+2*u$coef[7]*xx2^2+u$coef[9]*xx13+
			u$coef[10]*xx23)*log(ifelse(x3==0,1,x3))
		if(any(is.na(c(z1,z2,z3))))stop(paste("NAs in calculating estimates:",a,b,d))
		u <- glm(y~xx1+xx2+xx3+I(xx1^2)+I(xx2^2)+I(xx3^2)+
			xx12+xx13+xx23+z1+z2+z3,family=family,weight=weight)
		a <- a+u$coef[11]
		b <- b+u$coef[12]
		d <- d+u$coef[13]
		if(any(is.na(c(a,b,d))))stop(paste("NAs in calculating estimates:",a,b,d))
		i <- i+1
		test <- ((u$coef[11]^2>0.00001)||(u$coef[12]^2>0.00001)||
			(u$coef[13]^2>0.00001))&&(i<iterlim)}
	z <- glm(y~xx1+xx2+xx3+I(xx1^2)+I(xx2^2)+I(xx3^2)+xx12+xx13+xx23,
		family=family,weight=weight)
	z$df.residual <- z$df.residual-3
	z$aic <- z$aic+6
	z$powers <- c(a,b,d)
	z$iterations <- i
	class(z) <- c("rs",class(z))
	return(z)}

print.rs <- function(z,...){
	cat("\nPowered transformed response surface\n\n")
	cat("Powers:",z$powers,"\n")
	cat("Iterations:",z$iterations,"\n")
	print.glm(z,...)}

print.summary.rs <- function(z,...){
	cat("\nPowered transformed response surface\n\n")
	cat("Powers:",z$powers,"\n")
	cat("Iterations:",z$iterations,"\n")
	print.summary.glm(z,...)}

summary.rs <- function(z,...){
	zz <- summary.glm(z,...)
	class(zz) <- c("summary.rs",class(zz))
	if(!is.null(z$powers))zz$powers <- z$powers
	if(!is.null(z$iterations))zz$iterations <- z$iterations
	zz}
