library(pairwiseAdonis)
library(compositions)
library(vegan)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/")

map=read.table("Input_files/RSS_metadata.txt", header = T, row.names = 1, sep = "\t")
asv=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,1:68]

#clr trasnformation
asv_clr=as.data.frame(t(apply(asv,2,clr)))
rownames(asv_clr)=gsub('\\.', '-', rownames(asv_clr))
asv_clr$groups=paste(map$Experiment, map$Time, map$Temperature, sep = "_")[match(rownames(asv_clr), rownames(map))]
asv_clr$Temperature=map$Temperature[match(rownames(asv_clr), rownames(map))]
asv_clr$Time=map$Time[match(rownames(asv_clr), rownames(map))]
asv_clr$Genotype=map$Genotype[match(rownames(asv_clr), rownames(map))]

## CB
asv.acu=subset(asv_clr, rownames(asv_clr) %like% "^CB")
adonis(asv.acu[,1:8660] ~ asv.acu$Temperature + asv.acu$Time +  asv.acu$Genotype, method = "euclidean")
pairwise.adonis(asv.acu[,1:8660], asv.acu$groups, p.adjust.m ='fdr', sim.method = 'euclidean')

## RSS
asv.chr=subset(asv_clr, rownames(asv_clr) %like% "^RSS")
adonis(asv.chr[,1:8660] ~ asv.chr$Temperature + asv.chr$Time + asv.chr$Genotype,  method = "euclidean")
pairwise.adonis(asv.chr[,1:8660], asv.chr$groups, p.adjust.m ='fdr', sim.method = 'euclidean')


