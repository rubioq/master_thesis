###  TRANSCRIPTOMIC DIFFERENTIAL ANALYSIS  ###
###    AND GENE ONTOLOGY TERM ENRICHMENT   ###

# First, we source the script where custom functions are defined.
source("functions.R")

# We will use limma for differential analysis and ggplot2 to generate all plots.
library(limma)
library(ggplot2)

# We import the transcriptomic data into R
transc_data <- read.delim("../data/total_rna_data.csv",sep=";")
# and extract the transcriptomic data in a new matrix
# (columns 2 to 13 correspond to the degradome profiling, which
# we are not interested in)
transc_matrix <- transc_data[,14:25]
# We set the gene IDs as row names.
rownames(transc_matrix) <- transc_data[,1]
rownames(transc_data) <- transc_data[,1]
# We log transform the values in the matrix
transc_matrix <- log2(transc_matrix)
# and visualize the data in a boxplot
boxplot(transc_matrix) 
# As we can see, the data is already normalized

# Remove TSN data
transc_matrix <- transc_matrix[-c(31170,25396,25397),]

# PCA ----
library(FactoMineR)

# We perform principal component analysis using the data in transc_matrix
# and plot the corresponding scores plot.
pca.result <- PCA(t(transc_matrix),graph=FALSE, scale.unit=F)

pca.result$ind$coord <- as.data.frame(pca.result$ind$coord)
pca.result$ind$coord[,"Group"] <- rep(c("Col_CTL", "Ler_CTL", "tsnKO_CTL",
                                        "Col_Heat","Ler_Heat","tsnKO_Heat"),each=2)

ggplot(pca.result$ind$coord, aes(x=Dim.1, y=Dim.2)) + 
  geom_point(shape=21, aes(colour=Group),fill="black",size=2, stroke=2.5) +
  scale_color_manual( values=c("#e35a42","#d67c09","#f0c713","#9940b3","#32b85c","#4e6ccf") ) + 
  ylim (-100,100) + xlim(-120,120) +
  labs (title="PCA (transcriptomics)", x="PC1 (26.8%)", y="PC2 (19.4%)") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
# Together, principal components 1 and 2 explain 46.2 % of variability.
# We can see that samples belonging to the same group are very close
# to each other.

# Differential expression analysis ----

# We specify the experimental design (6 groups; 2 replicates each)
# In this study, instead of using the Columbia ecotype as control (wild-type),
# we also have Ler (Landsberg erecta) samples, because the genetic background 
# of tsnKO plants was still hybrid between the two ecotypes.
exp_design <- model.matrix(~ -1+factor(c(1,1,2,2,3,3,4,4,5,5,6,6)))
# In this script, we will often refer to a group using only two letters.
# The first one represents genotype (w/t) and the second one, treatment (c/h)
# wc: wildtype(col+ler)_ctl , wh: wildtype(col+ler)_heat , 
# tc: tsnko_ctl , th: tsnko_heat
colnames(exp_design) <- c("col_ctl","ler_ctl","tsnko_ctl","col_hs","ler_hs","tsnko_hs")

linearfit <- lmFit(transc_matrix, exp_design)

contrast_matrix <- makeContrasts((col_hs+ler_hs)/2-(col_ctl+ler_ctl)/2, # WT_CTL    vs WT_Heat
                                 tsnko_hs-tsnko_ctl,                    # tsnKO_CTL vs tsnKO_Heat
                                 (col_ctl+ler_ctl)/2-tsnko_ctl,         # WT_CTL    vs tsnKO_CTL
                                 (col_hs+ler_hs)/2-tsnko_hs,            # WT_Heat   vs tsnKO_Heat
                                 levels=colnames(exp_design))
# For comparisons with WT groups (Col and Ler), we need to specify in the model
# that they should be divided by two so that the comparison is possible with
# tsnko groups, which has only 2 replicates.

limma_results <- eBayes(contrasts.fit(linearfit, contrast_matrix))

# For each of the comparisons, we will extract differential analysis results
# and visualize the corresponding volcanoplot. The code used in each of the
# following sections is analogous, therefore only the first case will be
# commented.

# WT_CTL vs WT_Heat -------

# We store the results of the first contrast (coef=1) in a temporal variable.
temp_res <- topTable(limma_results, number=nrow(transc_matrix), coef=1, sort.by="logFC")
# The log fold-change and adjusted p values are saved in a data.frame, called
# trans_wc_wh (transcriptomics, wildtype_ctl vs wildtype_heat)
transc_wc_wh <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
# The rownames of this dataframe are gene IDs
rownames(transc_wc_wh) <- rownames(temp_res)
# We remove the temporal variable.
rm(temp_res)

# We add a column indicating whether a gene is activated (up),
# repressed (down) or not differentially expressed (ns).

for (i in 1:nrow(transc_wc_wh)){
  
  # If the change is significant (adjusted p < 0.05) and logfc > 1,
  # the gene is activated
  if (transc_wc_wh[i,"adj.p.val"]< 0.05 & transc_wc_wh[i,"logfc"]>1){
    transc_wc_wh[i,"significant"] <- "up" }
  
  # If the change is significant (adjusted p < 0.05) and logfc < -1,
  # the gene is repressed
  else if (transc_wc_wh[i,"adj.p.val"]< 0.05 & transc_wc_wh[i,"logfc"]<(-1)){
    transc_wc_wh[i,"significant"] <- "down" }
  
  # In any other case, we consider that the change is not significant.
  else {
    transc_wc_wh[i,"significant"] <- "ns" }
}

# Volcano plot
ggplot(transc_wc_wh, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#7f0fb8", "grey40", "#11ad3e") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="WT_CTL vs WT_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-12.5,12.5) + ylim(0,5) +
  # We will show the number of differentially expressed genes (up and down)
  annotate( geom="text", x=-11.25, y= 4.75, label = as.character(nrow(transc_wc_wh[transc_wc_wh$significant=="down",])), color="#7f0fb8", fontface="bold") +
  annotate( geom="text", x= 11.25, y= 4.75, label = as.character(nrow(transc_wc_wh[transc_wc_wh$significant=="up",])), color="#11ad3e", fontface="bold")

# tsnKO_CTl vs tsnKO_Heat -------

temp_res <- topTable(limma_results, number=nrow(transc_matrix), coef=2, sort.by="logFC")
transc_tc_th <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(transc_tc_th) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(transc_tc_th)){
  if (transc_tc_th[i,"adj.p.val"]< 0.05 & transc_tc_th[i,"logfc"]>1){
    transc_tc_th[i,"significant"] <- "up" }
  else if (transc_tc_th[i,"adj.p.val"]< 0.05 & transc_tc_th[i,"logfc"]<(-1)){
    transc_tc_th[i,"significant"] <- "down" }
  else {
    transc_tc_th[i,"significant"] <- "ns" }
}

ggplot(transc_tc_th, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#7f0fb8", "grey40", "#11ad3e") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="tsnKO_CTL vs tsnKO_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-12.5,12.5) + ylim(0,5) +
  annotate( geom="text", x=-11.25, y= 4.75, label = as.character(nrow(transc_tc_th[transc_tc_th$significant=="down",])), color="#7f0fb8", fontface="bold") +
  annotate( geom="text", x= 11.25, y= 4.75, label = as.character(nrow(transc_tc_th[transc_tc_th$significant=="up",])), color="#11ad3e", fontface="bold")

# WT_CTL vs tsnKO_CTL -------

temp_res <- topTable(limma_results, number=nrow(transc_matrix), coef=3, sort.by="logFC")
transc_wc_tc <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(transc_wc_tc) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(transc_wc_tc)){
  if (transc_wc_tc[i,"adj.p.val"]< 0.05 & transc_wc_tc[i,"logfc"]>1){
    transc_wc_tc[i,"significant"] <- "up" }
  else if (transc_wc_tc[i,"adj.p.val"]< 0.05 & transc_wc_tc[i,"logfc"]<(-1)){
    transc_wc_tc[i,"significant"] <- "down" }
  else {
    transc_wc_tc[i,"significant"] <- "ns" }
}

ggplot(transc_wc_tc, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#7f0fb8", "grey40", "#11ad3e") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="WT_CTL vs tsnKO_CTL") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-12.5,12.5) + ylim(0,5) +
  annotate( geom="text", x=-11.25, y= 4.75, label = as.character(nrow(transc_wc_tc[transc_wc_tc$significant=="down",])), color="#7f0fb8", fontface="bold") +
  annotate( geom="text", x= 11.25, y= 4.75, label = as.character(nrow(transc_wc_tc[transc_wc_tc$significant=="up",])), color="#11ad3e", fontface="bold")

# WT_Heat vs tsnKO_Heat -------

temp_res <- topTable(limma_results, number=nrow(transc_matrix), coef=4, sort.by="logFC")
transc_wh_th <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(transc_wh_th) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(transc_wh_th)){
  if (transc_wh_th[i,"adj.p.val"]< 0.05 & transc_wh_th[i,"logfc"]>1){
    transc_wh_th[i,"significant"] <- "up" }
  else if (transc_wh_th[i,"adj.p.val"]< 0.05 & transc_wh_th[i,"logfc"]<(-1)){
    transc_wh_th[i,"significant"] <- "down" }
  else {
    transc_wh_th[i,"significant"] <- "ns" }
}

ggplot(transc_wh_th, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#7f0fb8", "grey40", "#11ad3e") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="WT_Heat vs tsnKO_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-12.5,12.5) + ylim(0,5) +
  annotate( geom="text", x=-11.25, y= 4.75, label = as.character(nrow(transc_wh_th[transc_wh_th$significant=="down",])), color="#7f0fb8", fontface="bold") +
  annotate( geom="text", x= 11.25, y= 4.75, label = as.character(nrow(transc_wh_th[transc_wh_th$significant=="up",])), color="#11ad3e", fontface="bold")

# Export results ----

# We define an R list object in which all differentially expressed genes 
# lists are stored and export it. This variable contains only the gene IDs.
transc_diff<-list(
  wc_wh_up=rownames(transc_wc_wh[transc_wc_wh$significant=="up",]),wc_wh_down=rownames(transc_wc_wh[transc_wc_wh$significant=="down",]),
  tc_th_up=rownames(transc_tc_th[transc_tc_th$significant=="up",]),tc_th_down=rownames(transc_tc_th[transc_tc_th$significant=="down",]),
  wc_tc_up=rownames(transc_wc_tc[transc_wc_tc$significant=="up",]),wc_tc_down=rownames(transc_wc_tc[transc_wc_tc$significant=="down",]),
  wh_th_up=rownames(transc_wh_th[transc_wh_th$significant=="up",]),wh_th_down=rownames(transc_wh_th[transc_wh_th$significant=="down",]))

save(transc_diff,file="../results/diff_lists_transc.RData")

# The next three lines of code allow us to get the gene symbol (ex: MKK4 for
# mitogen-activated protein kinase kinase 4) of Arabidopsis genes.
library(org.At.tair.db)
getsymbols <- org.At.tairSYMBOL
getsymbols <- as.list(getsymbols[mappedkeys(getsymbols)])

# The export_genes_list writes tables containing all the results of our
# differential analysis: protein ID, -log10(adjusted p value), log fold-change,
# gene description and gene symbol
export_genes_list(vpdata=transc_wc_wh,path="../results/transc/transc_wt_ctl_wt_heat")
export_genes_list(vpdata=transc_tc_th,path="../results/transc/transc_tsnKO_ctl_tsnKO_heat")
export_genes_list(vpdata=transc_wc_tc,path="../results/transc/transc_wt_ctl_tsnko_ctl")
export_genes_list(vpdata=transc_wh_th,path="../results/transc/transc_wt_heat_tsnko_heat")

# Gene Ontology Over-Representation Analysis (ORA) -----

# We use the clusterProfiler library to perform enrichment analysis in ontologies
# Biological Process, Molecular Function and Cellular Component in the up- and
# down-regulated genes.
library(clusterProfiler)

# WT_CTL vs WT_Heat
bp_ora_wc_wh_up <- enrichGO(gene=get_proteins(transc_diff$wc_wh_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_wc_wh_up <- enrichGO(gene=get_proteins(transc_diff$wc_wh_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_wc_wh_up <- enrichGO(gene=get_proteins(transc_diff$wc_wh_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_wc_wh_down <- enrichGO(gene=get_proteins(transc_diff$wc_wh_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_wc_wh_down <- enrichGO(gene=get_proteins(transc_diff$wc_wh_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_wc_wh_down <- enrichGO(gene=get_proteins(transc_diff$wc_wh_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# tsnKO_CTL vs tsnKO_Heat
bp_ora_tc_th_up <- enrichGO(gene=get_proteins(transc_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_tc_th_up <- enrichGO(gene=get_proteins(transc_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_tc_th_up <- enrichGO(gene=get_proteins(transc_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_tc_th_down <- enrichGO(gene=get_proteins(transc_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_tc_th_down <- enrichGO(gene=get_proteins(transc_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_tc_th_down <- enrichGO(gene=get_proteins(transc_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# WT_CTL vs tsnKO_CTL
bp_ora_wc_tc_up <- enrichGO(gene=get_proteins(transc_diff$wc_tc_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_wc_tc_up <- enrichGO(gene=get_proteins(transc_diff$wc_tc_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_wc_tc_up <- enrichGO(gene=get_proteins(transc_diff$wc_tc_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_wc_tc_down <- enrichGO(gene=get_proteins(transc_diff$wc_tc_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_wc_tc_down <- enrichGO(gene=get_proteins(transc_diff$wc_tc_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_wc_tc_down <- enrichGO(gene=get_proteins(transc_diff$wc_tc_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# WT_Heat vs tsnKO_Heat
bp_ora_wh_th_up <- enrichGO(gene=get_proteins(transc_diff$wh_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_wh_th_up <- enrichGO(gene=get_proteins(transc_diff$wh_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_wh_th_up <- enrichGO(gene=get_proteins(transc_diff$wh_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_wh_th_down <- enrichGO(gene=get_proteins(transc_diff$wh_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_wh_th_down <- enrichGO(gene=get_proteins(transc_diff$wh_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_wh_th_down <- enrichGO(gene=get_proteins(transc_diff$wh_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# Because the enrichGO function takes some time, we will save all the ORA
# results in an .RData file to avoid running it again.
save(bp_ora_wc_wh_up,  mf_ora_wc_wh_up,  cc_ora_wc_wh_up,
     bp_ora_wc_wh_down,mf_ora_wc_wh_down,cc_ora_wc_wh_down,
     bp_ora_tc_th_up,  mf_ora_tc_th_up,  cc_ora_tc_th_up,
     bp_ora_tc_th_down,mf_ora_tc_th_down,cc_ora_tc_th_down,
     bp_ora_wc_tc_up,  mf_ora_wc_tc_up,  cc_ora_wc_tc_up,
     bp_ora_wc_tc_down,mf_ora_wc_tc_down,cc_ora_wc_tc_down,
     bp_ora_wh_th_up,  mf_ora_wh_th_up,  cc_ora_wh_th_up,
     bp_ora_wh_th_down,mf_ora_wh_th_down,cc_ora_wh_th_down,
     file="../results/enrichment_results_transc.RData")

# For each of the gene ontologies, results of the ORA will be summarized in one
# dot plot showing the enriched GO terms across all comparisons.
# Only the Biological Process dotplot is commented, since the code is analogous
# to that used for Cellular Component and Molecular Function plots.

# Biological Process (BP) dotplots ----

# Individual dotplots

# First, we take a look at the dot plots of each gene set, to ensure that no 
# important GO terms will be left behind in the final visualization.

options(enrichplot.colours=c("red","blue"))
dotplot( bp_ora_wc_wh_up,   title="WT_CTL vs WT_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_wc_wh_down, title="WT_CTL vs WT_Heat DOWN (BP)", showCategory=15)
dotplot( bp_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (BP)", showCategory=15)
dotplot( bp_ora_wc_tc_up,   title="WT_CTL vs tsnKO_CTL UP (BP)", showCategory=15)
dotplot( bp_ora_wc_tc_down, title="WT_CTL vs tsnKO_CT DOWN (BP)", showCategory=15)
dotplot( bp_ora_wh_th_up,   title="WT_Heat vs tsnKO_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_wh_th_down, title="WT_Heat vs tsnKO_Heat DOWN (BP)", showCategory=15)

# Since no significant enrichment is found in 5 of the 8 datasets, we will extract
# more than 5 terms of the first sets

# Dotplot with all comparisons

# We extract the IDs of the top 8 GO terms enriched in the gene sets
# WT UP, WT DOWN and tsnKO UP
dpgoids <- sort(unique(c(bp_ora_wc_wh_up@result$ID[1:8],bp_ora_wc_wh_down@result$ID[1:8],
                         bp_ora_tc_th_up@result$ID[1:8])))

# We still efine the eight gene sets we are studying
dpsets <- c("WT UP","WT DOWN","tsnKO UP","tsnKO DOWN",
            "WT/tsnKO CTL UP","WT/tsnKO CTL DOWN",
            "WT/tsnKO HS UP","WT/tsnKO HS DOWN")
dpsets <- factor(dpsets,levels=dpsets)

# We initialize two vectors: one for storing gene ratios and another for adjusted
# p values.
dpgeneratios <- c()
dppval <- c()

for (i.go in dpgoids){ # For each of the GO terms we will represent,
  
  # We extract the gene ratio
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_wc_wh_up@result[i.go,"GeneRatio"])))
  # and adjusted p value
  dppval <- c(dppval,bp_ora_wc_wh_up@result[i.go,"p.adjust"])
  # in each set
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_wc_wh_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_wc_wh_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_tc_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_tc_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_tc_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_tc_th_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_wc_tc_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_wc_tc_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_wc_tc_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_wc_tc_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_wh_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_wh_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_wh_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_wh_th_down@result[i.go,"p.adjust"])
  
}

# The dpgoids vector contains the IDs of GO terms (ex: GO:0010118). For the 
# visualization, we want to show the description of the GO term (ex: stomatal 
# movement), which we will store in the dpdescription vector.
dpdescription <- c()

for (i in 1:length(dpgoids)){ # For each GO term ID,
  # If the term is enriched in a specific set,
  if (!is.na(bp_ora_wc_wh_up@result[dpgoids[i],"Description"])){
    # we extract the description from the corresponding enrichResult object
    dpdescription[i] <- bp_ora_wc_wh_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(bp_ora_wc_wh_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- bp_ora_wc_wh_down@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(bp_ora_tc_th_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- bp_ora_tc_th_up@result[dpgoids[i],"Description"]
  }
}

# We need to extract descriptions from different comparisons because not all
# terms are enriched in the first set.

# We organize the data in a data frame that is compatible with ggplot2.
dpdata <- data.frame("GOs"=rep(dpdescription,each=length(dpsets)),
                     "Set"=rep(dpsets,length(dpdescription)),
                     "GeneRatio"=dpgeneratios,
                     "p.adjust"=dppval)

# In a dotplot, significance is represented by colours, mapped to the adjusted
# p values in logarithmic scale. We add a new column to the dataframe with
# -log10(adj p value).

for (i in 1:nrow(dpdata)){ # For each row (each term in each comparison)
  if (!is.na(dpdata[i,"p.adjust"])){ # If there exists an adjusted p value
    if (dpdata[i,"p.adjust"]>0.05){ # which is higher than 0.05 (not significant)
      dpdata[i,"log10p.adjust"] <- NA # we will add an NA to the -log10 column.
    }
    else { # If the enrichment is significant, we calculate the -log10 value.
      dpdata[i,"log10p.adjust"] <- log10(dpdata[i,"p.adjust"])
    }
  }
}


# Finally, we generate the plot, specifying that NA values should be ploted as
# transparent dots.
ggplot(data = dpdata, aes(x = Set, y = GOs, color = `log10p.adjust`, size = GeneRatio)) + 
  geom_point() + scale_color_gradient(low = "red", high = "blue" ,na.value="#00000000") +
  ylab("") + theme_gray() + xlab("") + ggtitle("Biological Process (BP)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"), axis.text.y=element_text(color="black"))

# Cellular Component (CC) dotplots ----

# Individual dotplots
# dotplot( cc_ora_wc_wh_up,   title="WT_CTL vs WT_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_wc_wh_down, title="WT_CTL vs WT_Heat DOWN (cc)", showCategory=15)
dotplot( cc_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (cc)", showCategory=15)
dotplot( cc_ora_wh_th_up,   title="WT_Heat vs tsnKO_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_wh_th_down, title="WT_Heat vs tsnKO_Heat DOWN (cc)", showCategory=15)

# Dotplot with all comparisons

dpgoids <- sort(unique(c(cc_ora_wc_wh_down@result$ID[1:8],
                         cc_ora_tc_th_up@result$ID[1:8],cc_ora_tc_th_down@result$ID[1:8])))

dpsets <- c("WT UP","WT DOWN","tsnKO UP","tsnKO DOWN",
            "WT/tsnKO CTL UP","WT/tsnKO CTL DOWN",
            "WT/tsnKO HS UP","WT/tsnKO HS DOWN")
dpsets <- factor(dpsets,levels=dpsets)

dpgeneratios <- c()
dppval <- c()

for (i.go in dpgoids){
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_wc_wh_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_wc_wh_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_wc_wh_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_wc_wh_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_tc_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_tc_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_tc_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_tc_th_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_wc_tc_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_wc_tc_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_wc_tc_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_wc_tc_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_wh_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_wh_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_wh_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_wh_th_down@result[i.go,"p.adjust"])
  
}

dpdescription <- c()

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_wc_wh_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_wc_wh_down@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_wc_wh_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_wc_wh_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_tc_th_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_tc_th_down@result[dpgoids[i],"Description"]
  }
}

dpdata <- data.frame("GOs"=rep(dpdescription,each=length(dpsets)),
                     "Set"=rep(dpsets,length(dpdescription)),
                     "GeneRatio"=dpgeneratios,
                     "p.adjust"=dppval)

for (i in 1:nrow(dpdata)){
  if (!is.na(dpdata[i,"p.adjust"])){
    if (dpdata[i,"p.adjust"]>0.05){
      dpdata[i,"log10p.adjust"] <- NA
    }
    else {
      dpdata[i,"log10p.adjust"] <- log10(dpdata[i,"p.adjust"])
    }
  }
}

ggplot(data = dpdata, aes(x = Set, y = GOs, color = `log10p.adjust`, size = GeneRatio)) + 
  geom_point() + scale_color_gradient(low = "red", high = "blue" ,na.value="#00000000") +
  ylab("") + theme_gray() + xlab("") + ggtitle("Cellular Component (CC)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"), axis.text.y=element_text(color="black"))


# Molecular Function (MF) dotplots ----

# Individual dotplots
dotplot( mf_ora_wc_wh_up,   title="WT_CTL vs WT_Heat UP (MF)", showCategory=15)
dotplot( mf_ora_wc_wh_down, title="WT_CTL vs WT_Heat DOWN (MF)", showCategory=15)
dotplot( mf_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (MF)", showCategory=15)
dotplot( mf_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (MF)", showCategory=15)
dotplot( mf_ora_wh_th_up,   title="WT_Heat vs tsnKO_Heat UP (MF)", showCategory=15)
dotplot( mf_ora_wh_th_down, title="WT_Heat vs tsnKO_Heat DOWN (MF)", showCategory=15)

# Dotplot with all comparisons

dpgoids <- sort(unique(c(mf_ora_wc_wh_up@result$ID[1:5],mf_ora_wc_wh_down@result$ID[1:5],
                         mf_ora_tc_th_up@result$ID[1:5],mf_ora_tc_th_down@result$ID[1:5],
                         mf_ora_wc_tc_up@result$ID[1])))

dpsets <- c("WT UP","WT DOWN","tsnKO UP","tsnKO DOWN",
            "WT/tsnKO CTL UP","WT/tsnKO CTL DOWN",
            "WT/tsnKO HS UP","WT/tsnKO HS DOWN")
dpsets <- factor(dpsets,levels=dpsets)

dpgeneratios <- c()
dppval <- c()

for (i.go in dpgoids){
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_wc_wh_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_wc_wh_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_wc_wh_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_wc_wh_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_tc_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_tc_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_tc_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_tc_th_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_wc_tc_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_wc_tc_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_wc_tc_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_wc_tc_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_wh_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_wh_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_wh_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_wh_th_down@result[i.go,"p.adjust"])
  
}

dpdescription <- c()

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_wc_wh_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_wc_wh_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_wc_wh_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_wc_wh_down@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_tc_th_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_tc_th_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_tc_th_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_tc_th_down@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_wc_tc_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_wc_tc_up@result[dpgoids[i],"Description"]
  }
}

dpdata <- data.frame("GOs"=rep(dpdescription,each=length(dpsets)),
                     "Set"=rep(dpsets,length(dpdescription)),
                     "GeneRatio"=dpgeneratios,
                     "p.adjust"=dppval)

for (i in 1:nrow(dpdata)){
  if (!is.na(dpdata[i,"p.adjust"])){
    if (dpdata[i,"p.adjust"]>0.05){
      dpdata[i,"log10p.adjust"] <- NA
    }
    else {
      dpdata[i,"log10p.adjust"] <- log10(dpdata[i,"p.adjust"])
    }
  }
}

ggplot(data = dpdata, aes(x = Set, y = GOs, color = `log10p.adjust`, size = GeneRatio)) + 
  geom_point() + scale_color_gradient(low = "red", high = "blue" ,na.value="#00000000") +
  ylab("") + theme_gray() + xlab("") + ggtitle("Molecular Function (MF)") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"), axis.text.y=element_text(color="black"))
