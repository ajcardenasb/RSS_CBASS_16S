library(ANCOMBC)
library(phyloseq)

setwd("~/Documents/Bioinformatics_scripts/R_scripts/RSS_16S/ASVs/")


asv=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,1:68]
map=read.table("Input_files/RSS_metadata.txt", header = T, row.names = 1, sep = "\t")
tax=read.table("Input_files/RSSvsCBASS_ASV_table_pooled.txt", header = TRUE, row.names = 1)[,70:76]
colnames(asv)=gsub("\\.", "-", colnames(asv))

#add lowest tax level
tax$label=ifelse(tax$Genus == "NA" | is.na(tax$Genus), as.character(paste("Unclassified", tax$Family, sep = "_")), as.character(tax$Genus))
tax$label=ifelse(tax$label == "Unclassified_NA", as.character(paste("Unclassified", tax$Family, sep = "_")), as.character(tax$label))
tax$label=ifelse(tax$label == "Unclassified_NA", as.character(paste("Unclassified", tax$Order, sep = "_")), as.character(tax$label))
tax$label=ifelse(tax$label == "Unclassified_NA", as.character(paste("Unclassified", tax$Class, sep = "_")), as.character(tax$label))
tax$label=ifelse(tax$label == "Unclassified_NA", as.character(paste("Unclassified", tax$Phylum, sep = "_")), as.character(tax$label))
tax$label=ifelse(tax$label == "Unclassified_NA", as.character(paste("Unclassified", tax$Kingdom, sep = "_")), as.character(tax$label))
unique(tax$label)

otu.t= otu_table(asv, taxa_are_rows=TRUE)
sam.t= sample_data(data.frame(map))
tax.t= tax_table(as.matrix(tax))

phy.all= phyloseq(otu.t, tax.t,  sam.t)

###########################################################################################################################################
###`Comparisons based on 4 points ##### T1 vs T2 acute, T1 vs T2 chronic, intercepts acute and chronic, comparint T1 across experimetns ###
###########################################################################################################################################

#######################################
##### 27-T1 vs 29|32|34 T1 acute #####
######################################
#27 T1 vs 29 T1
C1=subset_samples(phy.all, Temperature %in% c("27", "29") & Experiment == "CB" & Time == "T1")
res1=ancombc(phyloseq=C1,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res1_df=data.frame(Beta=res1[["res"]][["beta"]], Beta_se=res1[["res"]][["se"]], W=res1[["res"]][["W"]],pval=res1[["res"]][["p_val"]], qval=res1[["res"]][["q_val"]],DA=res1[["res"]][["diff_abn"]])
colnames(res1_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res1_sig=subset(res1_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res1_sig$Diff_more_abundant=ifelse( res1_sig$W < 0 , "27", "29")
res1_sig$ASV=rownames(res1_sig)
res1_sig$Taxa=tax$label[match(rownames(res1_sig),rownames(tax))]
res1_sig$Comparison="Comparison 1: CB-T1-27vs29"
message("Number of DA ASVs: ", nrow(res1_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res1_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res1_sig, Diff_more_abundant == "29" )))

#27 T1 vs 32 T1
C2=subset_samples(phy.all, Temperature %in% c("27", "32") & Experiment == "CB" & Time == "T1")
res2=ancombc(phyloseq=C2,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res2_df=data.frame(Beta=res2[["res"]][["beta"]], Beta_se=res2[["res"]][["se"]], W=res2[["res"]][["W"]],pval=res2[["res"]][["p_val"]], qval=res2[["res"]][["q_val"]],DA=res2[["res"]][["diff_abn"]])
colnames(res2_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res2_sig=subset(res2_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res2_sig$Diff_more_abundant=ifelse( res2_sig$W < 0 , "27", "32")
res2_sig$ASV=rownames(res2_sig)
res2_sig$Taxa=tax$label[match(rownames(res2_sig),rownames(tax))]
res2_sig$Comparison="Comparison 1: CB-T1-27vs32"
message("Number of DA ASVs: ", nrow(res2_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res2_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 32: ", nrow(subset(res2_sig, Diff_more_abundant == "32" )))

#27 T1 vs 34 T1
C3=subset_samples(phy.all, Temperature %in% c("27", "34") & Experiment == "CB" & Time == "T1")
res3=ancombc(phyloseq=C3,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res3_df=data.frame(Beta=res3[["res"]][["beta"]], Beta_se=res3[["res"]][["se"]], W=res3[["res"]][["W"]],pval=res3[["res"]][["p_val"]], qval=res3[["res"]][["q_val"]],DA=res3[["res"]][["diff_abn"]])
colnames(res3_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res3_sig=subset(res3_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res3_sig$Diff_more_abundant=ifelse( res3_sig$W < 0 , "27", "34")
res3_sig$ASV=rownames(res3_sig)
res3_sig$Taxa=tax$label[match(rownames(res3_sig),rownames(tax))]
res3_sig$Comparison="Comparison 1: CB-T1-27vs34"
message("Number of DA ASVs: ", nrow(res3_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res3_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 34: ", nrow(subset(res3_sig, Diff_more_abundant == "34" )))

#######################################
##### 27-T3 vs 29|32|34 T3 acute #####
######################################
#27 T3 vs 29 T3
C4=subset_samples(phy.all, Temperature %in% c("27", "29") & Experiment == "CB" & Time == "T3")
res4=ancombc(phyloseq=C4,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res4_df=data.frame(Beta=res4[["res"]][["beta"]], Beta_se=res4[["res"]][["se"]], W=res4[["res"]][["W"]],pval=res4[["res"]][["p_val"]], qval=res4[["res"]][["q_val"]],DA=res4[["res"]][["diff_abn"]])
colnames(res4_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res4_sig=subset(res4_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res4_sig$Diff_more_abundant=ifelse( res4_sig$W < 0 , "27", "29")
res4_sig$ASV=rownames(res4_sig)
res4_sig$Taxa=tax$label[match(rownames(res4_sig),rownames(tax))]
res4_sig$Comparison="Comparison 1: CB-T3-27vs29"
message("Number of DA ASVs: ", nrow(res4_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res4_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res4_sig, Diff_more_abundant == "29" )))

#27 T3 vs 32 T3
C5=subset_samples(phy.all, Temperature %in% c("27", "32") & Experiment == "CB" & Time == "T3")
res5=ancombc(phyloseq=C5,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res5_df=data.frame(Beta=res5[["res"]][["beta"]], Beta_se=res5[["res"]][["se"]], W=res5[["res"]][["W"]],pval=res5[["res"]][["p_val"]], qval=res5[["res"]][["q_val"]],DA=res5[["res"]][["diff_abn"]])
colnames(res5_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res5_sig=subset(res5_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res5_sig$Diff_more_abundant=ifelse( res5_sig$W < 0 , "27", "32")
res5_sig$ASV=rownames(res5_sig)
res5_sig$Taxa=tax$label[match(rownames(res5_sig),rownames(tax))]
res5_sig$Comparison="Comparison 1: CB-T3-27vs32"
message("Number of DA ASVs: ", nrow(res5_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res5_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 32: ", nrow(subset(res5_sig, Diff_more_abundant == "32" )))

#27 T3 vs 34 T3
C6=subset_samples(phy.all, Temperature %in% c("27", "34") & Experiment == "CB" & Time == "T3")
res6=ancombc(phyloseq=C6,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res6_df=data.frame(Beta=res6[["res"]][["beta"]], Beta_se=res6[["res"]][["se"]], W=res6[["res"]][["W"]],pval=res6[["res"]][["p_val"]], qval=res6[["res"]][["q_val"]],DA=res6[["res"]][["diff_abn"]])
colnames(res6_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res6_sig=subset(res6_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res6_sig$Diff_more_abundant=ifelse( res6_sig$W < 0 , "27", "34")
res6_sig$ASV=rownames(res6_sig)
res6_sig$Taxa=tax$label[match(rownames(res6_sig),rownames(tax))]
res6_sig$Comparison="Comparison 1: CB-T3-27vs34"
message("Number of DA ASVs: ", nrow(res6_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res6_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 34: ", nrow(subset(res6_sig, Diff_more_abundant == "34" )))

#######################################
##### 27-T1 vs 29|32|34 T1 chronic #####
#######################################

#27 T1 vs 29 T1
C7=subset_samples(phy.all, Temperature %in% c("27", "29") & Experiment == "RSS" & Time == "T1")
res7=ancombc(phyloseq=C7,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res7_df=data.frame(Beta=res7[["res"]][["beta"]], Beta_se=res7[["res"]][["se"]], W=res7[["res"]][["W"]],pval=res7[["res"]][["p_val"]], qval=res7[["res"]][["q_val"]],DA=res7[["res"]][["diff_abn"]])
colnames(res7_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res7_sig=subset(res7_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res7_sig$Diff_more_abundant=ifelse( res7_sig$W < 0 , "27", "29")
res7_sig$ASV=rownames(res7_sig)
res7_sig$Taxa=tax$label[match(rownames(res7_sig),rownames(tax))]
res7_sig$Comparison="Comparison 1: RSS-T1-27vs29"
message("Number of DA ASVs: ", nrow(res7_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res7_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res7_sig, Diff_more_abundant == "29" )))

#27 T1 vs 32 T1
C8=subset_samples(phy.all, Temperature %in% c("27", "32") & Experiment == "RSS" & Time == "T1")
res8=ancombc(phyloseq=C8,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res8_df=data.frame(Beta=res8[["res"]][["beta"]], Beta_se=res8[["res"]][["se"]], W=res8[["res"]][["W"]],pval=res8[["res"]][["p_val"]], qval=res8[["res"]][["q_val"]],DA=res8[["res"]][["diff_abn"]])
colnames(res8_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res8_sig=subset(res8_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res8_sig$Diff_more_abundant=ifelse( res8_sig$W < 0 , "27", "32")
res8_sig$ASV=rownames(res8_sig)
res8_sig$Taxa=tax$label[match(rownames(res8_sig),rownames(tax))]
res8_sig$Comparison="Comparison 1: RSS-T1-27vs32"
message("Number of DA ASVs: ", nrow(res8_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res8_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 32: ", nrow(subset(res8_sig, Diff_more_abundant == "32" )))

#27 T1 vs 34 T1
C9=subset_samples(phy.all, Temperature %in% c("27", "34") & Experiment == "RSS" & Time == "T1")
res9=ancombc(phyloseq=C9,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res9_df=data.frame(Beta=res9[["res"]][["beta"]], Beta_se=res9[["res"]][["se"]], W=res9[["res"]][["W"]],pval=res9[["res"]][["p_val"]], qval=res9[["res"]][["q_val"]],DA=res9[["res"]][["diff_abn"]])
colnames(res9_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res9_sig=subset(res9_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res9_sig$Diff_more_abundant=ifelse( res9_sig$W < 0 , "27", "34")
res9_sig$ASV=rownames(res9_sig)
res9_sig$Taxa=tax$label[match(rownames(res9_sig),rownames(tax))]
res9_sig$Comparison="Comparison 1: RSS-T1-27vs34"
message("Number of DA ASVs: ", nrow(res9_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res9_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 34: ", nrow(subset(res9_sig, Diff_more_abundant == "34" )))

#######################################
##### 27-T2 vs 29|32|34 T2 chronic #####
#######################################

#27 T2 vs 29 T2
C10=subset_samples(phy.all, Temperature %in% c("27", "29") & Experiment == "RSS" & Time == "T3")
res10=ancombc(phyloseq=C10,formula="Temperature",p_adj_method = "fdr",zero_cut = 0.9,lib_cut=1000,group = "Temperature",struc_zero =TRUE,neg_lb = FALSE,tol = 1e-05,max_iter = 100,conserve = F,alpha = 0.05,global = TRUE)
res10_df=data.frame(Beta=res10[["res"]][["beta"]], Beta_se=res10[["res"]][["se"]], W=res10[["res"]][["W"]],pval=res10[["res"]][["p_val"]], qval=res10[["res"]][["q_val"]],DA=res10[["res"]][["diff_abn"]])
colnames(res10_df)=c(	"Beta",	"se", "W",	"pval",	"qval", "Diff_abundant")
res10_sig=subset(res10_df, Diff_abundant == "TRUE")#[,c(1,2,6,8,10)]
res10_sig$Diff_more_abundant=ifelse( res10_sig$W < 0 , "27", "29")
res10_sig$ASV=rownames(res10_sig)
res10_sig$Taxa=tax$label[match(rownames(res10_sig),rownames(tax))]
res10_sig$Comparison="Comparison 1: RSS-T1-27vs29"
message("Number of DA ASVs: ", nrow(res10_sig), "\nNumber of DA ASVs enriched in 27: ", nrow(subset(res10_sig, Diff_more_abundant == "27" )), "\nNumber of DA ASVs enriched in 29: ", nrow(subset(res10_sig, Diff_more_abundant == "29" )))


ANCOMresults=rbind(res1_sig,res2_sig,res3_sig,res4_sig,res5_sig)
ANCOMresults$Taxa=paste(tax$Phylum,tax$Class,tax$Order,tax$Family,tax$Genus, sep= ";" )[match(ANCOMresults$ASV, rownames(tax))]
write.table(ANCOMresults,  "output_plots/ANCOMBC_ASVs_results.txt", sep = "\t", quote = F, row.names = T )

library(reshape2)
library(ggplot2)

ANCOMresults_plot=ANCOMresults %>% group_by(Diff_more_abundant, Comparison) %>% tally()
ANCOMresults_plot$n=ifelse(ANCOMresults_plot$Diff_more_abundant %in% c("IUI", "day", "full"), ANCOMresults_plot$n*-1, ANCOMresults_plot$n*1)

pdf("./outputs/ANCOMBC_DA_barplots.pdf", width=6,height=4, pointsize = 12)
ggplot(data=ANCOMresults_plot, aes(x=Comparison, y=n)) + geom_bar(stat="identity", position = "dodge")  + geom_text(aes(label=n), vjust=0.5, color="white", position = position_dodge(0.8), size=3) +  theme_classic() + theme(axis.text.x=element_text(angle=90,hjust=1)) 
dev.off()




########################
##### intersects  #####
########################

#################################################
##### point 1 = 27-T1 vs 29|32|34 T1 acute #####
#################################################

# comparing T1 and T2 in 29 
message("Present in T1 and T2 (not resilient) in 29ºC : ", length(intersect(rownames(res1_sig), rownames(res4_sig))), " ASVs")
ove1=intersect(rownames(res1_sig), rownames(res4_sig))
subset(res1_sig, rownames(res1_sig) %in% ove1) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T1 (resilient) in 29ºC : ", length(setdiff(rownames(res1_sig),rownames(res4_sig))), " ASVs")
ove2=setdiff(rownames(res1_sig),rownames(res4_sig))
subset(res1_sig, rownames(res1_sig) %in% ove2) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T2 (resilient) in 29ºC : ", length(setdiff(rownames(res4_sig),rownames(res1_sig))), " ASVs")
ove3=setdiff(rownames(res4_sig),rownames(res1_sig))
tmp=subset(res4_sig, rownames(res4_sig) %in% ove3) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))


# comparing T1 and T2 in 32
message("Present in T1 and T2 (not resilient) in 32ºC : ", length(intersect(rownames(res2_sig), rownames(res5_sig))), " ASVs")
ove4=intersect(rownames(res2_sig), rownames(res5_sig))
subset(res2_sig, rownames(res2_sig) %in% ove4) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T1 (resilient) in 32ºC : ", length(setdiff(rownames(res2_sig),rownames(res5_sig))), " ASVs")
ove5=setdiff(rownames(res2_sig),rownames(res5_sig))
subset(res2_sig, rownames(res2_sig) %in% ove5) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T2 (resilient) in 32ºC : ", length(setdiff(rownames(res5_sig),rownames(res2_sig))), " ASVs")
ove6=setdiff(rownames(res5_sig),rownames(res2_sig))
tmp=subset(res5_sig, rownames(res5_sig) %in% ove6) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

# comparing T1 and T2 in 34
message("Present in T1 and T2 (not resilient) in 34ºC : ", length(intersect(rownames(res3_sig), rownames(res6_sig))), " ASVs")
ove7=intersect(rownames(res3_sig), rownames(res6_sig))
subset(res3_sig, rownames(res3_sig) %in% ove7) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T1 (resilient) in 34ºC : ", length(setdiff(rownames(res3_sig),rownames(res6_sig))), " ASVs")
ove8=setdiff(rownames(res3_sig),rownames(res6_sig))
subset(res3_sig, rownames(res3_sig) %in% ove8) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T2 (resilient) in 34ºC : ", length(setdiff(rownames(res6_sig),rownames(res3_sig))), " ASVs")
ove9=setdiff(rownames(res6_sig),rownames(res3_sig))
tmp=subset(res6_sig, rownames(res6_sig) %in% ove9) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))


#################################################
##### point 2 = 27-T1 vs 29|32|34 T1 chronic #####
#################################################

# comparing T1 and T2 in 29 
message("Present in T1 and T2 (not resilient) in 29ºC : ", length(intersect(rownames(res7_sig), rownames(res10_sig))), " ASVs")
ove10=intersect(rownames(res7_sig), rownames(res10_sig))
subset(res7_sig, rownames(res7_sig) %in% ove10) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T1 (resilient) in 29ºC : ", length(setdiff(rownames(res7_sig),rownames(res10_sig))), " ASVs")
ove11=setdiff(rownames(res7_sig),rownames(res10_sig))
subset(res7_sig, rownames(res7_sig) %in% ove11) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T2 (resilient) in 29ºC : ", length(setdiff(rownames(res10_sig),rownames(res7_sig))), " ASVs")
ove12=setdiff(rownames(res10_sig),rownames(res7_sig))
tmp=subset(res10_sig, rownames(res10_sig) %in% ove12) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))



##############################################################################
##### point 3 = common DAF present only in T1 between acute and chronic #####
#############################################################################

### 29
message("Present only in T1 at 29ºC in acute: ", length(setdiff(rownames(res1_sig),rownames(res4_sig))), " OTUs")
ove13=setdiff(rownames(res1_sig),rownames(res4_sig))
subset(res1_sig, rownames(res1_sig) %in% ove13) %>% group_by(Taxa) %>%  tally()  %>% arrange(desc(n))

message("Present only in T1 at 29ºC in chronic: ", length(setdiff(rownames(res7_sig),rownames(res10_sig))), " OTUs")
ove14=setdiff(rownames(res7_sig),rownames(res10_sig))
subset(res7_sig, rownames(res7_sig) %in% ove14) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

length(intersect(ove13, ove14))
#intersect(setdiff(rownames(res1_sig),rownames(res4_sig)), setdiff(rownames(res7_sig),rownames(res10_sig)))



### 32 and 34 cant' be done as there are not T2 samples  in chronic

#######################################################################################
##### point 4 = DAF 27-T1 vs 29|32|34-T1 acute | DAF 27-T1 vs 29|32|34 T1 chronic #####
#######################################################################################


# common DAF between (27-T1 vs 29-T1) acute and (27-T1 vs 29-T1) chronic
message("Present in T1 acute and T1 chronic at 29ºC: ", length(intersect(rownames(res1_sig), rownames(res7_sig))), " OTUs")
ove15=intersect(rownames(res1_sig),rownames(res7_sig))
subset(res1_sig, rownames(res1_sig) %in% ove15) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

message("Present ony in T1 acute at 29ºC: ", length(setdiff(rownames(res1_sig),rownames(res7_sig))), " OTUs")
ove16=setdiff(rownames(res1_sig),rownames(res7_sig))
subset(res1_sig, rownames(res1_sig) %in% ove16) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

message("Present ony in T1 chronic at 29ºC: ", length(setdiff(rownames(res7_sig),rownames(res1_sig))), " OTUs")
ove17=setdiff(rownames(res7_sig),rownames(res1_sig))
subset(res7_sig, rownames(res7_sig) %in% ove17) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

# common DAF between (27-T1 vs 32-T1) acute and (27-T1 vs 32-T1) chronic
message("Present in T1 acute and T1 chronic at 32ºC: ", length(intersect(rownames(res2_sig), rownames(res8_sig))), " OTUs")
ove18=intersect(rownames(res2_sig),rownames(res8_sig))
subset(res2_sig, rownames(res2_sig) %in% ove18) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

message("Present ony in T1 acute at 32ºC: ", length(setdiff(rownames(res2_sig),rownames(res8_sig))), " OTUs")
ove19=setdiff(rownames(res2_sig),rownames(res8_sig))
subset(res2_sig, rownames(res2_sig) %in% ove19) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

message("Present ony in T1 chronic at 32ºC: ", length(setdiff(rownames(res8_sig),rownames(res2_sig))), " OTUs")
ove20=setdiff(rownames(res8_sig),rownames(res2_sig))
subset(res8_sig, rownames(res8_sig) %in% ove20) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

# common DAF between (27-T1 vs 34-T1) acute and (27-T1 vs 34-T1) chronic
message("Present in T1 acute and T1 chronic at 34ºC: ", length(intersect(rownames(res3_sig), rownames(res9_sig))), " OTUs")
ove18=intersect(rownames(res3_sig),rownames(res9_sig))
subset(res3_sig, rownames(res3_sig) %in% ove18) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

message("Present ony in T1 acute at 34ºC: ", length(setdiff(rownames(res3_sig),rownames(res9_sig))), " OTUs")
ove19=setdiff(rownames(res3_sig),rownames(res9_sig))
subset(res3_sig, rownames(res3_sig) %in% ove19) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

message("Present ony in T1 chronic at 34ºC: ", length(setdiff(rownames(res9_sig),rownames(res3_sig))), " OTUs")
ove20=setdiff(rownames(res9_sig),rownames(res3_sig))
subset(res9_sig, rownames(res9_sig) %in% ove20) %>% group_by(Taxa) %>%  tally() %>% arrange(desc(n))

