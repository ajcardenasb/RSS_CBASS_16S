# RSS_CBASS_16S
Scripts for cDNA-based 16S rRNA amplicon data analysis of Stylophora pistillata in response to short-term (18 hours) and long-term (14 days) heat stress from the paper: "Fast and pervasive transcriptomic resilience and acclimation of extremely heat tolerant coral holobionts from the northern Red Sea"

The script 'RSS_CBASS_dada2.R' runs DADA2 to infer ASVs
The script 'RSS_CBASS_barplots.R' plots the 20 most abundant bacterial families
The script 'RSS_CBASS_alphaDiv.R' calculates alpha diversity across treatments
The script 'RSS_CBASS_ordination.R' plots PCAs based on Euclidian distances of clr-transformed counts
The script 'RSS_CBASS_permanovas.R' tests ASV overall composition across treatments based on Euclidian distances of clr-trasnformed counts
The script 'RSS_CBASS_ANCOMBC.R' runs [ANCOMC](https://github.com/FrederickHuangLin/ANCOMBC) to identify differentially abundant ASVs across treatments

