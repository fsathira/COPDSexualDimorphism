\name{COPDSexualDimorphism-package}
\alias{COPDSexualDimorphism-package}
\alias{COPDSexualDimorphism}
\docType{package}
\title{
Sexual Dimorphic and COPD Differential Analysis
}
\description{
Sexual dimoprhic and COPD differential (SDCD) analysis contrasts regression coefficients
from two stratified analysis. Stratification can be done in two ways: by COPD status or 
by sex. For COPD-stratified analysis, SDCD analysis contrasts sexual dimorphism between
cases and controls, while sex-stratified SDCD analsysis contrasts COPD differential 
expression pattern between males and females. The package is meant to be used in 
conjunction with the package limma.
}
\details{
\tabular{ll}{
Package: \tab COPDSexualDimorphism\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-09-02\cr
License: \tab LGPL-2.1\cr
}
Data and main methods to accompany "Integrative Genomics of Sexual Dimorphism in COPD." 
The main functions are \code{sdcd}, which contrasts regression coefficients from 
sex-stratified and COPD-stratified analyses.
}
\author{
J Fah Sathirapongsasuti

Maintainer: J Fah Sathirapongsasuti <fah@cs.stanford.edu>
}
\references{
Sathirapongsasuti JF, Glass K, Huttenhower C, Quackenbush J, DeMeo DL. Integrative
Genomics of Sexual Dimorphism in COPD. (In Prep).

www.lung-genomics.org
}
\keyword{ package }
\keyword{ COPD }
\keyword{ sexual dimorphism }
\keyword{ SDCD }
\examples{
data(lgrc.expr.meta)
data(lgrc.expr)
data(lgrc.genes)

library(limma)

## Sex-stratified
design.mtx = cbind(ctrl=1,
		copd=as.integer(grepl("COPD",colnames(expr))),
		age=expr.meta$age,
		pkyr=expr.meta$pkyrs)

good.idx = apply(design.mtx,1,function(x){!any(is.na(x))}) & (expr.meta$gender == "1-Male")
male.fit = lmFit(log(expr)[,good.idx], design.mtx[good.idx,])
male.fit = eBayes(male.fit)

good.idx = apply(design.mtx,1,function(x){!any(is.na(x))}) & (expr.meta$gender == "2-Female")
female.fit = lmFit(log(expr)[,good.idx], design.mtx[good.idx,])
female.fit = eBayes(female.fit)

male.female.copd.beta.diff.genes = sdcd(male.fit, female.fit, "copd", lgrc.genes, fdr.cutoff=0.25, file.prefix="male.female.copd", write.file=FALSE)

## COPD-stratified
design.mtx = cbind(ctrl=1,
		gender=expr.meta$gender,
		age=expr.meta$age,
		pkyr=expr.meta$pkyrs)

good.idx = apply(design.mtx,1,function(x){!any(is.na(x))}) & grepl("COPD",colnames(expr))
copd.fit = lmFit(log(expr)[,good.idx], design.mtx[good.idx,])
copd.fit = eBayes(copd.fit)

good.idx = apply(design.mtx,1,function(x){!any(is.na(x))}) & grepl("CTRL",colnames(expr))
ctrl.fit = lmFit(log(expr)[,good.idx], design.mtx[good.idx,])
ctrl.fit = eBayes(ctrl.fit)

copd.ctrl.gender.beta.diff.genes = sdcd(copd.fit, ctrl.fit, "gender", lgrc.genes, fdr.cutoff=0.25, file.prefix="copd.ctrl.gender", class.names=c("copd","ctrl"), write.file=FALSE)

## Combine
sdcd.genes = merge(copd.ctrl.gender.beta.diff.genes, male.female.copd.beta.diff.genes, by=setdiff(intersect(names(copd.ctrl.gender.beta.diff.genes), names(male.female.copd.beta.diff.genes)),c("beta.diff","beta.diff.pooled.sd")))
sdcd.genes = unique(sdcd.genes)
print(paste("There are", nrow(sdcd.genes),"SDCD genes"))
}
