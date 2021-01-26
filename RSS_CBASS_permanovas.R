library(pairwiseAdonis)
library(compositions)
library(vegan)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/")

map=read.table("Input_files/RSS_metadata.txt", header = T, row.names = 1, sep = "\t")
asv=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,1:68]

#asv_clr=as.data.frame(t(apply(asv,2,clr)))
asv_trn=as.data.frame(t(sweep(asv,2,colSums(asv),"/")))
rownames(asv_trn)=gsub('\\.', '-', rownames(asv_trn))
asv_trn$groups=paste(map$Experiment, map$Time, map$Temperature, sep = "_")[match(rownames(asv_trn), rownames(map))]
asv_trn$Temperature=map$Temperature[match(rownames(asv_trn), rownames(map))]
asv_trn$Time=map$Time[match(rownames(asv_trn), rownames(map))]

## CB
asv.acu=subset(asv_trn, rownames(asv_trn) %like% "^CB")
adonis(asv.acu[,1:8660] ~ asv.acu$Temperature + asv.acu$Time)
acu.pai=pairwise.adonis(asv.acu[,1:8660], asv.acu$groups, p.adjust.m ='fdr', sim.method = 'bray')

## RSS
asv.chr=subset(asv_trn, rownames(asv_trn) %like% "^RSS")
adonis(asv.chr[,1:8660] ~ asv.chr$Temperature + asv.chr$Time)
chr.pai=pairwise.adonis(asv.chr[,1:8660], asv.chr$groups, p.adjust.m ='fdr', sim.method = 'bray')




