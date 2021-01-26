library(vegan)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(GUniFrac)
setwd("~/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/")

map=read.table("Input_files/RSS_metadata.txt", header = T, row.names = 1, sep = "\t")
asv=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,1:68]
colnames(asv)=gsub("\\.", "-", colnames(asv))
map$Temperature=paste(map$Temperature, "ºC", sep = "")

###rarefying
cnts=t(asv[, 1:68])
min(rowSums(cnts)) # determine sample with lowest counts
asv.rar=Rarefy(cnts, 12485)$otu.tab.rff


###plots
alpha=as.data.frame(t(estimateR(asv.rar,  smallsample = TRUE)))
alpha$Shannon=diversity(asv.rar, index = "shannon")#$shannon
alpha$Evenness=(alpha$Shannon/log10(alpha$S.obs)) ### adding eveness.  J = H'/ln(S) where H' is Shannon Weiner diversity and S is the total number of species 
alpha$Temperature=map$Temperature[match(rownames(alpha),rownames(map))]
alpha$Time=map$Time[match(rownames(alpha),rownames(map))]
alpha$Experiment=map$Experiment[match(rownames(alpha),rownames(map))]

stats=alpha %>% group_by(Experiment, Temperature, Time) %>% summarise(averageObservedASVs=mean(S.obs),sdevObservedASVs=sd(S.obs) )
mean(stats$averageObservedASVs)
mean(stats$sdevObservedASVs)
# boxplots
P4=c("#2E33D1", "#FFEE32","#D37D47", "#F43535")
a=ggplot(alpha, aes(x=Temperature, y=S.obs, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Observed number of ASVs", x="", title = "Observed") + facet_grid(Experiment~Time)
b=ggplot(alpha, aes(x=Temperature, y=S.ACE, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Observed number of ASVs", x="", title = "ACE") + facet_grid(Experiment~Time)
c=ggplot(alpha, aes(x=Temperature, y=Shannon, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Shannon diversity", x="", title = "Shannon") + facet_grid(Experiment~Time)
d=ggplot(alpha, aes(x=Temperature, y=Evenness, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Evenness", x="", title = "Evenness")  + facet_grid(Experiment~Time)
grid.arrange(a,b,c,d, ncol=2)


##################################################
##################### Stats ######################
##################################################


shapiro.test(alpha$S.ACE) # p-value > 0.05 implying we can assume the normality.
shapiro.test(alpha$S.obs) # p-value > 0.05 implying we can assume the normality.
shapiro.test(alpha$Shannon) # p-value > 0.05 implying we can assume the normality.
shapiro.test(alpha$Evenness) # p-value > 0.05 implying we can assume the normality.

#ANOVAs
summary(aov(alpha$S.obs ~ alpha$Temperature))
summary(aov(alpha$S.ACE ~ alpha$Temperature))
summary(aov(alpha$Shannon ~ alpha$Temperature))
summary(aov(alpha$Evenness ~ alpha$Temperature))

#Pairwise t-test
pairwise.t.test(alpha$S.obs,alpha$Temperature, p.adj = "fdr")
#Calculate paiwise for ACE, Shannon and Evenness

## Plot with significance groups
#pdf("./outputs/ASVs_skeleton16S_Shannon.pdf", width=6.5,height=3, pointsize = 12)
ggplot(alpha, aes(x=Temperature, y=S.obs, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Observed number of ASVs", x="", title = "Observed") + annotate(geom="text", x=1, y=150, label= "A") + annotate(geom="text", x=2, y=185, label= "B")+ annotate(geom="text", x=3, y=295, label= "C") + annotate(geom="text", x=4, y=245, label= "BC")  ") +theme(legend.position = "none")      
ggplot(alpha, aes(x=Temperature, y=S.ACE, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "ACE estimated ASV richness", x="", title = "ACE") + annotate(geom="text", x=1, y=150, label= "?") + annotate(geom="text", x=2, y=185, label= "?")+ annotate(geom="text", x=3, y=295, label= "?") + annotate(geom="text", x=4, y=245, label= "?")  +theme(legend.position = "none")
ggplot(alpha, aes(x=Temperature, y=Shannon, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Shannon diversity", x="", title = "Shannon") + annotate(geom="text", x=1, y=3.4, label= "?") + annotate(geom="text", x=2, y=4.45, label= "?")+ annotate(geom="text", x=3, y=4.55, label= "?") + annotate(geom="text", x=4, y=4.65, label= "?") +theme(legend.position = "none")      
ggplot(alpha, aes(x=Temperature, y=Evenness, fill=Temperature)) + stat_boxplot(geom = "errorbar")  + geom_boxplot(alpha = 1) +  scale_fill_manual(values=P4)   +  theme_classic() + labs( y= "Evenness", x="", title = "Evenness") + annotate(geom="text", x=1, y=1.73, label= "?") + annotate(geom="text", x=2, y=1.97, label= "?")+ annotate(geom="text", x=3, y=1.9, label= "?") + annotate(geom="text", x=4, y=1.99, label= "?") +theme(legend.position = "none")     

#dev.off()
#grid.arrange(a,b,c,d, ncol=2)


