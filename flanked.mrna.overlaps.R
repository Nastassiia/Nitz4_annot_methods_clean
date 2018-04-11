#Find overlaps between mrnas with 1000 flanking regions basing on gff annotation.
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(IRanges))

flanked.mrna.overlaps<-function(gff.table, flank){
  tab.mrna<-filter(gff.table, type=='mRNA')
  tab.mrna_<-tab.mrna
  colnames(tab.mrna_)<-paste0(colnames(tab.mrna),'_')
  tab.mrna_$start_<-((tab.mrna_$start_)-flank)
  tab.mrna_$end_<-((tab.mrna_$end_)+flank)
  chr.names<-unique(tab.mrna_$seqnames_) %>% as.character
  
  sumoverlap<-data.frame(matrix(ncol=2))
  colnames(sumoverlap)<-c('queryHits', 'subjectHits')
  for (nam in chr.names){
    chr.tab<-filter(tab.mrna_, seqnames_==nam) 
    prot.names<-unique(chr.tab$transcript_id_)
    for (tr.name in prot.names){
      prot.check<-filter(chr.tab, transcript_id_==tr.name)
      other.prot<-filter(chr.tab, transcript_id_!=tr.name)
      prot.check<-makeGRangesFromDataFrame(prot.check, start.field='start_', end.field='end_', seqnames.field='seqnames_', strand.field='strand_', keep.extra.columns = T)
      other.prot<-makeGRangesFromDataFrame(other.prot, start.field='start_', end.field='end_', seqnames.field='seqnames_', strand.field='strand_', keep.extra.columns = T)
      overlap<-findOverlaps(prot.check, other.prot) %>% as.data.frame
      if (nrow(overlap)>0) {sumoverlap<-rbind(sumoverlap, c(nam, tr.name))}
      sumoverlap<-rbind(sumoverlap,overlap)}
  }
  sumoverlap
}


