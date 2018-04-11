#Submit GFF table and generate 5/3 UTR annotations with CDS 
build.utr.annot<-function(gff.table){
  cds.mrna<-filter(gff.table, type %in% c('mRNA', 'CDS'))
  gene.ids<-cds.mrna$Dbxref
  gene.ids<-sub(".*, \"","",gene.ids)
  gene.ids<-sub("\".*","", gene.ids)
  cds.mrna.sh<-cds.mrna[, 1:10]
  cds.mrna.sh$ID<-gene.ids
  gene.ID<-unique(cds.mrna.sh$ID)
  utr.table<-data.frame(matrix(data=NA, ncol=10), stringsAsFactors = F)
  cds.mrna.sh$seqnames<-as.character(cds.mrna.sh$seqnames)
  cds.mrna.sh$strand<-as.character(cds.mrna.sh$strand)
  cds.mrna.sh$source<-as.character(cds.mrna.sh$source)
  cds.mrna.sh$type<-as.character(cds.mrna.sh$type)
  colnames(utr.table)<-colnames(cds.mrna.sh)

  for (gene in gene.ID){
    one.gene.tab<-filter(cds.mrna.sh, ID==gene)
    nrow.1.gene.tab<-nrow(one.gene.tab)
    utr.gene.tab<-data.frame(matrix(data=0, nrow=(nrow.1.gene.tab+1), ncol=10), stringsAsFactors = F)
    colnames(utr.gene.tab)<-colnames(cds.mrna.sh)
    k<-nrow(utr.gene.tab)
    if (one.gene.tab$strand[1]=='+'){
      utr.gene.tab[2:(k-1),]<-one.gene.tab[-1,]
      utr.gene.tab[1,]<-one.gene.tab[1,]
      utr.gene.tab$end[1]<-(one.gene.tab$start[2])
      utr.gene.tab$width[1]<-(utr.gene.tab$end[1]-utr.gene.tab$start[1])
      utr.gene.tab$type[1]<-'5\'-UTR'
      utr.gene.tab[k,]<-one.gene.tab[1,]
      utr.gene.tab$start[k]<-(one.gene.tab$end[nrow.1.gene.tab])
      utr.gene.tab$width[k]<-(utr.gene.tab$end[k]-utr.gene.tab$start[k])
      utr.gene.tab$type[k]<-'3\'-UTR'
      if (utr.gene.tab$width[1]>0) {
        utr.gene.tab$end[1]<-(utr.gene.tab$end[1]-1)}
      if (utr.gene.tab$width[k]>0){
        utr.gene.tab$start[k]<-(utr.gene.tab$start[k]+1)}
      utr.table<-rbind(utr.table, utr.gene.tab)}
    if (one.gene.tab$strand[1]=='-'){
      utr.gene.tab[2:(k-1),]<-one.gene.tab[-1,]
      utr.gene.tab[1,]<-one.gene.tab[1,]
      utr.gene.tab$start[1]<-(one.gene.tab$end[2]) #speify the end of 5UTR nearby TSS
      utr.gene.tab$width[1]<-(utr.gene.tab$end[1]-utr.gene.tab$start[1])
      utr.gene.tab$type[1]<-'5\'-UTR'
      utr.gene.tab[k,]<-one.gene.tab[1,]
      utr.gene.tab$end[k]<-(one.gene.tab$start[nrow.1.gene.tab]) #specify end of 3UTR nearby CDSs-end
      utr.gene.tab$width[k]<-(utr.gene.tab$end[k]-utr.gene.tab$start[k])
      utr.gene.tab$type[k]<-'3\'-UTR'
      if (utr.gene.tab$width[1]>0) {
        utr.gene.tab$start[1]<-(utr.gene.tab$start[1]+1)} #adjust not to overlap with gene boundaries
      if (utr.gene.tab$width[k]>0){
        utr.gene.tab$end[k]<-(utr.gene.tab$end[k]-1)}#adjust not to overlap with gene boundaries 
      utr.table<-rbind(utr.table, utr.gene.tab)}
    }
  
  utr.table
}
  
