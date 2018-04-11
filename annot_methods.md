Code for the paper



- Correction of substitution errors


```
ace 30000000 /storage/ao006/Nitz4_project/raw_data/genome_reads/FCD20HPACXX-SZAIPI025576-84_L1_1.fq /scratch/ao006/Nitz4_project/ACE_output/unz_fw_nitz4_ACE.fq

ace 30000000 /storage/ao006/Nitz4_project/raw_data/genome_reads/FCD20HPACXX-SZAIPI025576-84_L1_2.fq /scratch/ao006/Nitz4_project/ACE_output/unz_rev__nitz4_ACE.fq
```



- Trimming adapters and low-quality read-ends


```
java -jar /share/apps/trimmomatic/Trimmomatic-0.32/trimmomatic-0.32.jar PE \
/scratch/ao006/Nitz4_project/ACE_output/unz_fw_nitz4_ACE.fq.gz /scratch/ao006/Nitz4_project/ACE_output/unz_rev__nitz4_ACE.fq.gz \
Nitz4_trimmed_fw.fq.gz Nitz4_trim_junk_fw.fq.gz Nitz4_trimmed_rev.fq.gz Nitz4_trim_junk_rev.fq.gz \
ILLUMINACLIP:/storage/ao006/nitzschia_sp_nitz4/illumina_adapter_database/TruSeq_adapters.fa:2:40:15 LEADING:2 TRAILING:2 \
SLIDINGWINDOW:4:2 MINLEN:30 TOPHRED64 2>>nitz4_std_eror.log
```


- Fastqc quality check


```
# for untrimmed reads
/share/apps/fastqc/FastQC_v0.11.5/fastqc /storage/ao006/Nitz4_project/ACE_output/unz_fw_nitz4_ACE.fq
/share/apps/fastqc/FastQC_v0.11.5/fastqc /storage/ao006/Nitz4_project/ACE_output/unz_rev__nitz4_ACE.fq

#for trimmed reads
/share/apps/fastqc/FastQC_v0.11.5/fastqc /storage/ao006/Nitz4_project/trimmomat_output/Nitz4_trimmed_fw.fq
/share/apps/fastqc/FastQC_v0.11.5/fastqc /storage/ao006/Nitz4_project/trimmomat_output/Nitz4_trimmed_rev.fq
```

- Assembly with range of K-mer lengths

```
# a є {15, 21,27,33,39,45,51,59,63}
mpirun -n 16 /share/apps/Ray/Ray-2.3.1/Ray -k a -amos -p /storage/ao006/Nitz4_project/trimmomat_output/Nitz4_trimmed_fw.fq \
/storage/ao006/Nitz4_project/trimmomat_output/Nitz4_trimmed_rev.fq -o nitz4.ray.K45
```
- R quality check
Here will be the graphs for K-mer statistics


- Get rid of plastid reads with bowtie2

```
#database from Nitzschia4 plastid genome
/share/apps/bowtie2/bowtie2-2.2.8/bin/bowtie2-build --seed 1 /storage/ao006/Nitz4_project/raw_data/plastom/nitzschia4_complete_plastom.fa nitz4.plast

#align reads to plastid-genome db

/share/apps/bowtie2/bowtie2-2.2.8/bin/bowtie2 --seed 1 -p 16 -x /storage/ao006/Nitz4_project/bowt2_output/bowt.DB.nitz4.plastom/nitz4.plast \
-1 /storage/ao006/Nitz4_project/trimmomat_output/Nitz4_trimmed_fw.fq \
-2 /storage/ao006/Nitz4_project/trimmomat_output/Nitz4_trimmed_rev.fq \
--phred64 --local --un-conc-gz Nitz4_noplastom_%.fq.gz \
--al-conc-gz Nitz4_plast_reads_%.fq.gz

```
- Assemble genome without plastid reads

```
mpirun -n 16 /share/apps/Ray/Ray-2.3.1/Ray -show-memory-usage -k 45 -amos -p /storage/ao006/Nitz4_project/bowt2_output/bowt_nitz4_aln/Nitz4_noplastom_1.fq  \
/storage/ao006/Nitz4_project/bowt2_output/bowt_nitz4_aln/Nitz4_noplastom_2.fq  -o nitz4.ray.noplast.K45
```

- Find mitochondrial scaffolds


Blast plastid free assembly to diatom mitochondrial genes from custom databases


[R code to find mitochondrial scaffolds ](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/gene.fold.names.md)




- Get rid of mitochondrial reads

```
#create database from mitochondrial contig (Scaffold-4000015)
/share/apps/bowtie2/bowtie2-2.2.8/bin/bowtie2-build --seed 1 /scratch/ao006/Nitz4_project/bowt.get.rid.mito.reads/cont4000015.fasta cont4.15

#align reads to mitochondrial contig
/share/apps/bowtie2/bowtie2-2.2.8/bin/bowtie2 --seed 1 -p 16 -x /scratch/ao006/Nitz4_project/bowt.get.rid.mito.reads/cont.4000015.DB/cont4.15 \
-1 /storage/ao006/Nitz4_project/bowt2_output/bowt_nitz4_aln/Nitz4_noplastom_1.fq.gz \
-2 /storage/ao006/Nitz4_project/bowt2_output/bowt_nitz4_aln/Nitz4_noplastom_2.fq.gz \
--phred64 --local --un-conc-gz Nitz4_no_organel_%.fq.gz \
--al-conc-gz Nitz4_mito_reads_%.fq.gz

```

- Assemble genome without organellar reads with range of K-mer length

```
# a є {21,27,33,39,45,51,57,63}
mpirun -n 16 /share/apps/Ray/Ray-2.3.1/Ray -k a -amos -show-memory-usage -p /storage/ao006/Nitz4_project/bowt.get.rid.mito.reads/nitz4_cont4.15_aln/Nitz4_no_organel_1.fq \
/storage/ao006/Nitz4_project/bowt.get.rid.mito.reads/nitz4_cont4.15_aln/Nitz4_no_organel_2.fq  -o nitz4.no.org.K45

```

#### Identification of non-target genomes with blobplots

[Blobplots for whole, plastid-free and organelle-free reads scaffolds](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/make_blobloplots.md)

#### Generate Phaeodactylum tricornutum chosen gene-set for retraining augustus gene-prediction parameters

[Generate training gene-set for augustus (no UTR training)](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/Phaeod.train.whole.Rmd)  

[Generate training gene-set for augustus (UTR training)](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/utr.training.Rmd)
[Check coding sequences for correct start codons](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/check_ATGs_correct.Rmd)
 - The main difference between gene sets for UTR and general parameters training is that for UTRs we do not need to filter for intron-less genes, so eventually we get more genes for training set.  

[Generate training gene-set for augustus UTR-parameters training](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/utr.training.Rmd)
  - Custome scripts used

[Buid.utr.annot.R](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/build.utr.annot.R)  
[flanked.mrna.overlaps.R](https://github.com/Nastassiia/Nitz4_annot_methods_clean/blob/master/flanked.mrna.overlaps.R)

#### Train augustus

```
# Train parameters, EXCEPT UTRs. Make gene bank file out of generated GFF. File *.fna -Phaeodactylum_tricornutum fasta from Genbank with locus sequences (named like >NC_011701.1).
gff2gbSmallDNA.pl exec.utr.cds.train.1.gff Ph.tricorn.short.names_genomic.fna 1000 exec.utr.cds.train.1.gb  

randomSplit.pl exec.utr.cds.train.1.gb 150 # it generates testing gene set of 150 random chosen genes. All other genes go to the training set named exec.utr.cds.train.1.gb.train.

new_species.pl --species=Phaeodactylum_tricornutum # generate files for training
etraining --species=Phaeodactylum_tricornutum exec.utr.cds.train.1.gb.train #initial evaluation of parameters

optimize_augustus.pl --cpus=4 --rounds=5 --species=Phaeodactylum_tricornutum exec.utr.cds.train.1.gb.train #actual trraining, with no UTRs

#utr parameters training
optimize_augustus.pl --species=Phaeodactylum_tricornutum --cpus=4 --rounds=5 ../utr.training/just.utr.training.gb --UTR=on --metapars=~/soft/augustus-3.2.2/config/species/Phaeodactylum_tricornutum/Phaeodactylum_tricornutum_metapars.utr.cfg --trainOnlyUtr=1

# test trained parameters on test gene set, generated from splitted Genbank file exec.utr.cds.train.1.gb
augustus --species=Phaeodactylum_tricornutum exec.utr.cds.train.1.gb.test

```
#### Annotation of Nitzschia4 nuclear genome with MAKER
