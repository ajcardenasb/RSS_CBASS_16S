setwd("~/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/")

library(phyloseq)
library(ggplot2)
library(reshape2)
library(microbiome)

asv=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,1:68]
map=read.table("Input_files/RSS_metadata.txt", header = T, row.names = 1, sep = "\t")
tax=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,70:76]
colnames(asv)=gsub("\\.", "-", colnames(asv))

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

##transform data and subset phyloseq objects
phy.t=microbiome::transform(phy.all, transform = "clr", target = "OTU", shift = 0, scale = 1)
cb=subset_samples(phy.t, Experiment=="CB")
rss=subset_samples(phy.t, Experiment=="RSS")

P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535") #27, 29, 32, 34
cb_PCOA_br = ordinate(cb, method = "PCoA", distance = "euclidean")
c.p=plot_ordination(cb,cb_PCOA_br, color = "Temperature", shape = "Time")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("CBASS") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + scale_shape_manual(values=c(0, 15))#+ scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6))
rss_PCOA_br = ordinate(rss, method = "PCoA", distance = "euclidean")
r.p=plot_ordination(rss,rss_PCOA_br, color = "Temperature",  shape = "Time")  + geom_point(size = 3, alpha = 1) + theme_bw()  + ggtitle("RSS") + theme(plot.title = element_text(hjust = 0.5)) + scale_colour_manual(values=P4) + scale_shape_manual(values=c(0, 15))#+ scale_shape_manual(values=c(15, 16, 17, 18, 25, 0, 1, 2, 5, 6))

pdf("CBASS_vs_RSS_ordination.pdf", onefile = TRUE, width=15,height=15)
gridExtra::grid.arrange(c.p, r.p,  ncol=2)
dev.off()

write.table(cb_PCOA_br$vectors[,c(1,2)], "output_plots/PCA_ASVs_CBASS.txt", quote= FALSE, row.names = TRUE, sep = "\t" )
write.table(rss_PCOA_br$vectors[,c(1,2)], "output_plots/PCA_ASVs_RSS.txt", quote= FALSE, row.names = TRUE, sep = "\t" )

