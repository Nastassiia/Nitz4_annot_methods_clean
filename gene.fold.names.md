---
title: "mito.genes.folder.names"
output: html_document
data: "23 Nov 2016"
---

```{r}
knitr::opts_chunk$set(echo = TRUE)
```

Remove '.fasta' part
```{r}
setwd('/home/nata/Desktop/Nitz4/sep.mito.aln/')
gene.list<-readLines('gene.folder.names.txt')
length(gene.list)
61
i<-1
for (i in 1:length(gene.list)){
  a<-gene.list[i]
  gene.list[i]<-gsub("\\.fasta", "",a)}
length(gene.list)
write(gene.list, file='gene.fold.short.nam.txt', ncolumns = 61)

```


Write code for blasting Nitz4 scaffolds against mitochondrial gene databases
```{r}
/share/apps/blast/ncbi-blast-2.3.0+/bin/makeblastdb -in /storage/ao006/Nitz4_project/raw_data/mito.gene.groups/'genenme'.fasta \
-dbtype nucl  -out /'genenme'/atp6.all.db -title paste0('Nitz4.db', 'genenme')

genenme<-gene.list[3]

for (i in 1:length(gene.list)){
  genename<-gene.list[i]
  a<-paste0('/share/apps/blast/ncbi-blast-2.3.0+/bin/makeblastdb -in /storage/ao006/Nitz4_project/raw_data/mito.gene.groups/', genename, '.fasta -dbtype nucl  -out ', genename,'/', genename,'.db')
  write(a, 'command.lines.txt', append = T)
}

for (i in 1:length(gene.list)){
  genename<-gene.list[i]
  a<-paste0('/share/apps/blast/ncbi-blast-2.3.0+/bin/tblastx -num_threads 16 -db /scratch/ao006/Nitz4_project/find.mitoch.reads.cont/mito.genes.dbs/', genename, '/', genename, '.db -evalue 1e-4 -query /scratch/ao006/Nitz4_project/RAY.outp.after.bowt.plast.filter/RAY_noplast_K45/nitz4.ray.noplast.K45/Contigs.fasta -out aln.', genename)
   write(a, 'comm.lines.tblastx.txt', append = T)}
```


The same for ORF genes
```{r}

genes.ORF<-readLines('gene.ORF.names')
length(genes.ORF)
30
i<-1
for (i in 1:length(genes.ORF)){
  a<-genes.ORF[i]
  genes.ORF[i]<-gsub("\\.fasta", "",a)}
length(genes.ORF)
write(genes.ORF, file='genes.ORF.short', ncolumns = 61)
```

```{r}

for (i in 1:length(genes.ORF)){
  genename<-genes.ORF[i]
  a<-paste0('/share/apps/blast/ncbi-blast-2.3.0+/bin/makeblastdb -in /storage/ao006/Nitz4_project/raw_data/mito.genes.ORF/', genename, '.fasta -dbtype nucl  -out ', genename,'/', genename,'.db')
  write(a, 'command.lines.ORF.txt', append = T)
}

for (i in 1:length(genes.ORF)){
  genename<-genes.ORF[i]
  a<-paste0('/share/apps/blast/ncbi-blast-2.3.0+/bin/tblastx -num_threads 16 -db /scratch/ao006/Nitz4_project/find.mitoch.reads.cont/mito.ORF.genes.dbs/', genename, '/', genename, '.db -evalue 1e-4 -query /scratch/ao006/Nitz4_project/RAY.outp.after.bowt.plast.filter/RAY_noplast_K45/nitz4.ray.noplast.K45/Contigs.fasta -out aln.', genename)
   write(a, 'comm.lines.ORF.tblastx.txt', append = T)}
```
