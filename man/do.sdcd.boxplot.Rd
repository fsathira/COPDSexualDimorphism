\name{do.sdcd.boxplot}
\alias{do.sdcd.boxplot}
\alias{my.smart.plot}
\title{
Represent sexual dimorphic data as boxplots.
}
\description{
Stratify and plot data by gender and by disease status.
}
\usage{
do.sdcd.boxplot(marker, data, copd.bool, male.bool, symbol = marker, filename = paste(marker,".pdf",sep=""), take.log = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{marker}{
Character string representing the name of the marker (e.g. Ensembl gene, VMR) to plot.
}
  \item{data}{
Data to be plot. E.g. gene expression or percent methylation.
}
  \item{copd.bool}{
Boolean array corresponding to the COPD subject columns of \code{data}. 
}
  \item{male.bool}{
Boolean array corresponding to the male subject columns of \code{data}. 
}
  \item{symbol}{
Gene symbol associated with the marker.
}
  \item{filename}{
File name. If NA, display plot on screen.
}
  \item{take.log}{
Boolean. If true, plot in log scale.
}
}
\references{
Sathirapongsasuti JF, Glass K, Huttenhower C, Quackenbush J, DeMeo DL. Integrative
Genomics of Sexual Dimorphism in COPD. (In Prep).
}
\author{
J Fah Sathirapongsausti
}

\examples{
	data(lgrc.sdcd.genes)
	data(lgrc.expr.meta)
	data(lgrc.expr)
	copd.bool = (expr.meta$diagmaj == "2-COPD/Emphysema")
	male.bool = (expr.meta$gender == "1-Male")
	do.sdcd.boxplot("ENSG00000182472", expr, copd.bool, male.bool, symbol="CAPN12", filename=NA)
}

\keyword{ COPD }
\keyword{ sexual dimorphism }
\keyword{ boxplot }
