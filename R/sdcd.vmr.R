sdcd.vmr <-
function(male.fit, female.fit, coeff, genes, fdr.cutoff=0.25, stat=c("z","t"), annotate=FALSE, annotate.with=c("genes","NCBI"), file.prefix="male.female.copd", class.names=c("male","female"), write.file=TRUE) {
  `%+%` <- function(x,y) paste(x,y,sep="")
  mysdcd = sdcd.core(male.fit, female.fit, coeff, stat)
  
  male.female.copd.beta.diff.vmr = row.names(mysdcd)[mysdcd$beta.diff.pval.adj < fdr.cutoff]
  
  mf.vmr.df = data.frame(vmr=male.female.copd.beta.diff.vmr)
  mf.vmr.mtx = matrix(unlist(strsplit(as.character(male.female.copd.beta.diff.vmr),"_")),ncol=4,byrow=TRUE)
  mf.vmr.df$Chr = substring(mf.vmr.mtx[,2],4)
  mf.vmr.df$start = as.numeric(mf.vmr.mtx[,3])
  mf.vmr.df$end = as.numeric(mf.vmr.mtx[,4])
  mf.vmr.df$pos = as.integer((mf.vmr.df$start + mf.vmr.df$end)/2)
  mf.vmr.df$vmr = as.character(mf.vmr.df$vmr)
  
  if (annotate) {
    if (annotate.with == "NCBI") {
      require(NCBI2R) # there is a problem with hg18 vs hg19. Charm uses hg18 while NCBI2R uses hg19, so I have to convert everything to hg19 even before this analysis. 
      vmr.genes = GetNeighGenes(mf.vmr.df$Chr, mf.vmr.df$pos, FlankingDistance = 1e+04, web=FALSE, pbar=FALSE) # 10 kb based on Bell et al. Genome Biology 2011, 12:R10
      vmr.genes$chr = "chr" %+% vmr.genes$chr
    } else {
      vmr.genes = getNeighbors(mf.vmr.df, genes, dist=1e4)
    }
  } else {
    vmr.genes = rep(NA, nrow(mf.vmr.df))
  }
  
  male.female.copd.beta.diff.genes = cbind(mf.vmr.df, vmr.genes,
                                           male.beta=male.fit$coefficients[mf.vmr.df$vmr,coeff],
                                           female.beta=female.fit$coefficients[mf.vmr.df$vmr,coeff],
                                           male.sd=(male.fit$stdev.unscaled*male.fit$sigma)[mf.vmr.df$vmr,coeff],
                                           female.sd=(female.fit$stdev.unscaled*female.fit$sigma)[mf.vmr.df$vmr,coeff],
                                           male.t=male.fit$t[mf.vmr.df$vmr,coeff],
                                           female.t=female.fit$t[mf.vmr.df$vmr,coeff])
  male.female.copd.beta.diff.genes$male.female.p = mysdcd[mf.vmr.df$vmr,"beta.diff.pval"]
  male.female.copd.beta.diff.genes$male.female.p.adj = mysdcd[mf.vmr.df$vmr,"beta.diff.pval.adj"]
  male.female.copd.beta.diff.genes$beta.diff = mysdcd[mf.vmr.df$vmr,"beta.diff"]
  
  names(male.female.copd.beta.diff.genes) = sub("male", class.names[1], sub("female", class.names[2], names(male.female.copd.beta.diff.genes)))
  
  if (write.file) write.table(male.female.copd.beta.diff.genes, file=file.prefix %+% ".beta.diff.vmr.genes.withXY.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
  
  male.female.copd.beta.diff.genes = subset(male.female.copd.beta.diff.genes, !grepl("chr[YX]", vmr))
  
  if (write.file) {
    if (annotate) write.table(male.female.copd.beta.diff.genes[,c("chr","start","end","vmr")], file=file.prefix %+% ".beta.diff.vmr.bed", quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')
    write.table(male.female.copd.beta.diff.genes, file=file.prefix %+% ".beta.diff.vmr.genes.txt", quote=FALSE, col.names=TRUE, row.names=FALSE, sep='\t')
  }
  
  print("Number of probes with sexual dimorphic VMR: " %+% nrow(male.female.copd.beta.diff.genes))

  return(male.female.copd.beta.diff.genes)
}
