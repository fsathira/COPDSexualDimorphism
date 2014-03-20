sdcd.core <-
function(male.fit, female.fit, coeff, stat=c("z","t")) {
  stat = stat[1]
  if (stat == "z") {
    get.pval = function(z, se, mydf) { 2*(1-pnorm(abs(z), mean=0, sd=se)) }
  } else if (stat == "t") {
    get.pval = function(z, se, mydf) { 2*(1-pt(abs(z)/se, mydf)) }
  }

  mysdcd = data.frame(beta.diff = male.fit$coefficients[,coeff] - female.fit$coefficients[,coeff])
  mysdcd$beta.diff.norm = mysdcd$beta.diff
  names(mysdcd$beta.diff.norm) = row.names(male.fit$coefficients)
  mysdcd$beta.diff.pooled.sd = pool.sd(male.fit, female.fit, coeff)
  mysdcd$beta.diff.pval = get.pval(mysdcd$beta.diff.norm, mysdcd$beta.diff.pooled.sd, male.fit$df.residual + female.fit$df.residual)
  mysdcd$beta.diff.pval.adj = p.adjust(mysdcd$beta.diff.pval,"BH")
  row.names(mysdcd) = row.names(male.fit$coefficients)
  
  return(mysdcd)
}
