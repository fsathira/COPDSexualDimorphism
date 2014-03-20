pool.sd <-
function(fit1, fit2, coef, use.df=c("residual","prior")) {
  use.df = use.df[1]
  fit1.beta.var = ((fit1$stdev.unscaled*fit1$sigma)[,coef])^2
  fit2.beta.var = ((fit2$stdev.unscaled*fit2$sigma)[,coef])^2
  if (use.df == "residual") {
    df1 = fit1$df.residual
    df2 = fit2$df.residual
  } else {
    df1 = fit1$df.prior
    df2 = fit2$df.prior
  }
  return(sqrt((df1*fit1.beta.var + df2*fit2.beta.var)/(df1 + df2)))
}
