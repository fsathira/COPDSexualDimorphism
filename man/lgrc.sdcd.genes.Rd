\name{lgrc.sdcd.genes}
\alias{lgrc.sdcd.genes}
\alias{sdcd.genes}
\docType{data}
\title{
Sexual dimorphic and COPD differential genes
}
\description{
A list of 959 genes whose expressions were significantly different between
males and females and between COPD and controls. These were identified by
the function \code{sdcd}. 
}
\usage{data(lgrc.sdcd.genes)}
\format{
  The format is:
 num [1:14557, 1:229] 10.45 1.88 12 9.02 8.81 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:14557] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" ...
  ..$ : chr [1:229] "LT001098RU_COPD" "LT001796RU_CTRL" "LT005419RU_COPD" "LT007392RU_COPD" ...
}
\references{
Sathirapongsasuti JF, Glass K, Huttenhower C, Quackenbush J, DeMeo DL. Integrative
Genomics of Sexual Dimorphism in COPD. (In Prep).
}
\examples{
data(lgrc.sdcd.genes)
head(sdcd.genes)
}
\keyword{datasets}
