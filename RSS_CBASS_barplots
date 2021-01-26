
library(reshape2)
library(ggplot2)
library(scales)
library(dplyr)
library(gridExtra)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/")


#############################################################
#####Taxonomic profiles of the 20 most abundant families#####
#############################################################

asv=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)
map=read.table("Input_files/RSS_metadata.txt", header = T, row.names = 1, sep = "\t")
colnames(asv)=gsub("\\.", "-", colnames(asv))

names(asv)
asv.tax.ag=aggregate(asv[, 1:68], by = list(asv[, 74]), FUN =  sum) #define sample range and group factor
#topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[, 2:ncol(asv.tax.ag)]),decreasing = TRUE),][1:20,1]
topFamilies=asv.tax.ag[order(rowSums(asv.tax.ag[,2:ncol(asv.tax.ag) ]),decreasing = TRUE),][1:20,1] # top only in skeleton control samples 
fam.top=subset(asv.tax.ag, asv.tax.ag$Group.1 %in% topFamilies) 
fam.bot=subset(asv.tax.ag, !asv.tax.ag$Group.1 %in% topFamilies) 
fam.bot$Group.1=gsub(".*","zOthers", fam.bot$Group.1)
others=aggregate(fam.bot[, 2:ncol(fam.bot)], by = list(fam.bot[, 1]), FUN =  sum)
all.2 =rbind(fam.top, others)
all.l=melt(all.2, id.vars=c("Group.1"), variable.name = "Family", value.name = "Abundance")
colnames(all.l)=c("Family","Sample","Abundance")


## Add sample information
all.l$Experiment=map$Experiment[match(all.l$Sample, rownames(map))]
all.l$Time=map$Time[match(all.l$Sample, rownames(map))]
all.l$Temperature=map$Temperature[match(all.l$Sample, rownames(map))]
all.l$Genotype=map$Genotype[match(all.l$Sample, rownames(map))]

all.l.1=all.l %>% filter(Experiment == "CB") %>% group_by( Time, Temperature, Genotype, Family) %>% dplyr::summarise(Abundance=sum(Abundance))
all.l.2=all.l %>% filter(Experiment == "RSS") %>% group_by( Time, Temperature, Genotype,  Family) %>%  dplyr::summarise(Abundance=sum(Abundance))


## Plot
P21=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#C0C0C0")

cb=ggplot() +geom_bar(aes(y = Abundance, x = Genotype, fill = Family), data = all.l.1, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1), legend.position = "none",  plot.title = element_text(hjust = 0.5)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P21) + facet_grid( Time ~ Temperature) + labs(title="short-term heat stress") 
rs=ggplot() +geom_bar(aes(y = Abundance, x = Genotype, fill = Family), data = all.l.2, stat="identity", position = "fill")  +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "", x="Host colony") + scale_fill_manual(values=P21) + facet_grid(Time ~ Temperature) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + labs(title="long-term heat stress") 

pdf("/Users/anny/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/output_plots/CBASSvsRSS_barplots.pdf", onefile = TRUE, width=7,height=5, pointsize = 12)
grid.arrange(cb, rs, ncol = 2)
dev.off() 


pdf("/Users/anny/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/output_plots/labels_barplots.pdf", onefile = TRUE, width=7,height=5, pointsize = 12)
ggplot() +geom_bar(aes(y = Abundance, x = Genotype, fill = Family), data = all.l, stat="identity", position = "fill") +  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) + theme(axis.text.x=element_text(angle=90,hjust=1)) + labs( y= "Percentage of 16S rRNA sequences", x="Host colony") + scale_fill_manual(values=P21) + facet_grid( Time ~ Temperature)  + labs(title="CBASS")  +  theme( legend.key.size = unit(0.2, "cm"),legend.key.width = unit(0.2,"cm"), legend.position = 'right') + guides(fill=guide_legend(ncol=1))
dev.off() 

