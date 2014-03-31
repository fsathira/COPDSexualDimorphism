sdcd <-
function(male.fit, female.fit, coeff, genes, fdr.cutoff=0.25, stat=c("z","t"), file.prefix="male.female.copd", class.names=c("male","female"), write.file=TRUE) {
  `%+%` <- function(x,y) paste(x,y,sep="")
  mysdcd = sdcd.core(male.fit, female.fit, coeff, stat)
  
  male.female.copd.beta.diff.ensgenes = row.names(mysdcd)
  male.female.copd.beta.diff.genes.all = cbind(genes[male.female.copd.beta.diff.ensgenes,],
                                               male.beta=male.fit$coefficients[,coeff],
                                               female.beta=female.fit$coefficients[,coeff],
                                               male.sd=(male.fit$stdev.unscaled*male.fit$sigma)[,coeff],
                                               female.sd=(female.fit$stdev.unscaled*female.fit$sigma)[,coeff],
                                               male.t=male.fit$t[,coeff],
                                               female.t=female.fit$t[,coeff],
                                               male.p.value=male.fit$p.value[,coeff],
                                               female.p.value=female.fit$p.value[,coeff]
  )
  male.female.copd.beta.diff.genes.all$male.female.p = mysdcd$beta.diff.pval
  male.female.copd.beta.diff.genes.all$male.female.p.adj = mysdcd$beta.diff.pval.adj
  male.female.copd.beta.diff.genes.all$beta.diff = mysdcd$beta.diff
  male.female.copd.beta.diff.genes.all$beta.diff.pooled.sd = mysdcd$beta.diff.pooled.sd
  
  male.female.copd.beta.diff.idx = which(mysdcd$beta.diff.pval.adj < fdr.cutoff)
  male.female.copd.beta.diff.genes = male.female.copd.beta.diff.genes.all[male.female.copd.beta.diff.idx,]
  
  names(male.female.copd.beta.diff.genes) = sub("male", class.names[1], sub("female", class.names[2], names(male.female.copd.beta.diff.genes)))
  
  if (write.file) write.table(male.female.copd.beta.diff.genes, file=file.prefix %+% ".beta.diff.genes.withXY.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
  
  male.female.copd.beta.diff.genes = subset(male.female.copd.beta.diff.genes, !chromosome_name %in% c("X","Y"))
  
  if (write.file) {
    write.table(male.female.copd.beta.diff.genes, file=file.prefix %+% ".beta.diff.genes.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    write.table(cbind("chr" %+% male.female.copd.beta.diff.genes$chromosome_name, male.female.copd.beta.diff.genes[,c("start_position","end_position")]), file=file.prefix %+% ".beta.diff.genes.bed", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
  }
  
  print("Number of probes with sexual dimorphic differential expression: " %+% nrow(male.female.copd.beta.diff.genes))

  return(male.female.copd.beta.diff.genes)
}
