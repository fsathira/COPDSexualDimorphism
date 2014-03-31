getNeighbors <-
function(region, genes, dist) {
  # region has to have "chr", "start", "end", and "name"
  # genes has to have "chromosome_name", "start_position", "end_position", "ensembl_gene_id", and "strand"
  require(GenomicRanges)
  
  genes.bed = GRanges(seqnames=Rle(paste("chr",genes$chromosome_name,sep="")),
                         ranges=IRanges(start=genes$start_position, 
                                       end=genes$end_position, 
                                       names=genes$ensembl_gene_id),
                         strand=Rle(strand(genes$strand)))
  
  region.bed = GRanges(seqnames=Rle(paste("chr",region$Chr,sep="")),
                    ranges=IRanges(start=region$start,
                                  end=region$end,
                                  names=region$name))
  
  window = 2*dist
  genes.region.map = findOverlaps(resize(region.bed, window, fix="center"), genes.bed)
  
  region.genes = matrix(NA, nrow=nrow(region), ncol=3, dimnames=list(region$name,c("ensembl_gene_id","genesymbol","entrezgene")))
  for (i in 1:nrow(region)) {
    queryHit.idx = which(queryHits(genes.region.map) == i)
    if (length(queryHit.idx) > 0) {
      subjHit.idx = subjectHits(genes.region.map)[queryHit.idx]
      region.genes[i,"ensembl_gene_id"] = paste(genes$ensembl_gene_id[subjHit.idx], collapse=",")
      region.genes[i,"genesymbol"] = paste(genes$hgnc_symbol[subjHit.idx], collapse=",")
      region.genes[i,"entrezgene"] = paste(genes$entrezgene[subjHit.idx], collapse=",")
    }
  }
  region.genes = cbind(region.genes, chr=paste("chr",region$Chr,sep=""))
  
  return(region.genes)
}
