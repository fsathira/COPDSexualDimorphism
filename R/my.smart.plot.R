my.smart.plot <-
function(x, y, colSet=c("Blues","Reds","Greens"), type=c("plot","points"), pch=20, ...) {
  require(RColorBrewer)
  colSet = colSet[1]
  type = type[1]
  shades = densCols(x, y, colramp=colorRampPalette(brewer.pal(9,colSet)[-(1:5)]))
  if (type == "plot") {
    plot(x, y, col=shades, pch=pch, ...)
  } else if (type == "points") {
    points(x, y, col=shades, pch=pch, ...)
  }
}
