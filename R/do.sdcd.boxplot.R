do.sdcd.boxplot <-
function(marker, data, copd.bool, male.bool, symbol=marker, filename=paste(marker,".pdf",sep=""), take.log=FALSE) {
  `%+%` <- function(x,y) paste(x,y,sep="")
  require(beeswarm)
  if (!marker %in% row.names(data)) return()
  if (take.log) iexpr.c = log(data[marker,]) else iexpr.c = data[marker,]
  iexpr.list = list("Male COPD"=iexpr.c[copd.bool & male.bool],
                    "Male Control"=iexpr.c[!copd.bool & male.bool],
                    "Female COPD"=iexpr.c[copd.bool & !male.bool],
                    "Female Control"=iexpr.c[!copd.bool & !male.bool])
  
  if (!is.na(filename)) pdf(filename)
  boxplot(iexpr.list, main="Methylation of " %+% symbol, ylab="Percent methylation", border="black", boxwex=0.6)
  beeswarm(iexpr.list, col=c("cadetblue","cadetblue","pink","pink"), add=TRUE)
  if (!is.na(filename)) dev.off()
}
