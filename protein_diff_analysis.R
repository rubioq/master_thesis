###  PROTEOMIC DIFFERENTIAL ANALYSIS  ###
### AND GENE ONTOLOGY TERM ENRICHMENT ###

# First, we source the script where custom functions are defined.
source("functions.R")

# We will use limma for differential analysis and ggplot2 to generate all plots.
library(limma)
library(ggplot2)

# We load the preprocessed proteomic data, exported from the 
# data_preprocessing_protein.R script.
load(file="../data/protein_data.RData")

# We specify the experimental design (4 groups; 4 replicates each)
exp_design <- model.matrix(~ -1+factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)))
# In this script, we will often refer to a group using only two letters.
# The first one represents genotype (c/t) and the second one, treatment (c/h)
# cc: col_ctl , ch: col_heat , tc: tsnko_ctl , th: tsnko_heat
colnames(exp_design) <- c("cc","ch","tc","th")

linearfit <- lmFit(protein_normalized, exp_design)

contrast_matrix <- makeContrasts(ch-cc, # Col_CTL   vs  Col_Heat
                                 th-tc, # tsnKO_CTL vs  tsnKO_Heat
                                 tc-cc, # Col_CTL   vs  tsnKO_CTL
                                 th-ch, # Col_Heat  vs  tsnKO_Heat
                                 levels=c("cc","ch","tc","th"))

limma_results <- eBayes(contrasts.fit(linearfit, contrast_matrix))

# For each of the comparisons, we will extract differential analysis results
# and visualize the corresponding volcanoplot. The code used in each of the
# following sections is analogous, therefore only the first case will be
# commented.

# Col_CTl vs Col_Heat -------

# We store the results of the first contrast (coef=1) in a temporal variable.
temp_res <- topTable(limma_results, number=nrow(protein_normalized), coef=1, sort.by="logFC")
# The log fold-change and adjusted p values are saved in a data.frame, called
# protein_cc_ch (proteomics, col_ctl vs col_heat)
protein_cc_ch <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
# The rownames of this dataframe are protein names.
rownames(protein_cc_ch) <- rownames(temp_res)
# We remove the temporal variable.
rm(temp_res)

# We add a column indicating whether a protein is up-regulated (up),
# down-regulated (down) or not differentially expressed (ns).

for (i in 1:nrow(protein_cc_ch)){
  
  # If the change is significant (adjusted p < 0.05) and logfc > 1,
  # the protein is up-regulated
  if (protein_cc_ch[i,"adj.p.val"]< 0.05 & protein_cc_ch[i,"logfc"]>1){
    protein_cc_ch[i,"significant"] <- "up" }
  
  # If the change is significant (adjusted p < 0.05) and logfc < -1,
  # the protein is down-regulated
  else if (protein_cc_ch[i,"adj.p.val"]< 0.05 & protein_cc_ch[i,"logfc"]<(-1)){
    protein_cc_ch[i,"significant"] <- "down" }
  
  # In any other case, we consider that the change is not significant.
  else {
    protein_cc_ch[i,"significant"] <- "ns" }
}

# Volcano plot
ggplot(protein_cc_ch, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="Col_CTL vs Col_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-10,10) + ylim(0,10) +
  # We will show the number of differentially expressed proteins (up and down)
  annotate( geom="text", x=-9, y= 9.5, label = as.character(nrow(protein_cc_ch[protein_cc_ch$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x= 9, y= 9.5, label = as.character(nrow(protein_cc_ch[protein_cc_ch$significant=="up",])), color="#e35023", fontface="bold")

# tsnKO_CTl vs tsnKO_Heat -------

temp_res <- topTable(limma_results, number=nrow(protein_normalized), coef=2, sort.by="logFC")
protein_tc_th <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(protein_tc_th) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(protein_tc_th)){
  if (protein_tc_th[i,"adj.p.val"]< 0.05 & protein_tc_th[i,"logfc"]>1){
    protein_tc_th[i,"significant"] <- "up" }
  else if (protein_tc_th[i,"adj.p.val"]< 0.05 & protein_tc_th[i,"logfc"]<(-1)){
    protein_tc_th[i,"significant"] <- "down" }
  else {
    protein_tc_th[i,"significant"] <- "ns" }
}

ggplot(protein_tc_th, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="tsnKO_CTL vs tsnKO_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-10,10) + ylim(0,10) +
  annotate( geom="text", x=-9, y= 9.5, label = as.character(nrow(protein_tc_th[protein_tc_th$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x= 9, y= 9.5, label = as.character(nrow(protein_tc_th[protein_tc_th$significant=="up",])), color="#e35023", fontface="bold") 

# Col_CTL vs tsnKO_CTL -------

temp_res <- topTable(limma_results, number=nrow(protein_normalized), coef=3, sort.by="logFC")
protein_cc_tc <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(protein_cc_tc) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(protein_cc_tc)){
  if (protein_cc_tc[i,"adj.p.val"]< 0.05 & protein_cc_tc[i,"logfc"]>1){
    protein_cc_tc[i,"significant"] <- "up" }
  else if (protein_cc_tc[i,"adj.p.val"]< 0.05 & protein_cc_tc[i,"logfc"]<(-1)){
    protein_cc_tc[i,"significant"] <- "down" }
  else {
    protein_cc_tc[i,"significant"] <- "ns" }
}

ggplot(protein_cc_tc, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="Col_CTL vs tsnKO_CTL") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-10,10) + ylim(0,10) +
  annotate( geom="text", x=-9, y= 9.5, label = as.character(nrow(protein_cc_tc[protein_cc_tc$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x= 9, y= 9.5, label = as.character(nrow(protein_cc_tc[protein_cc_tc$significant=="up",])), color="#e35023", fontface="bold")
  
# Col_Heat vs tsnKO_Heat -------

temp_res <- topTable(limma_results, number=nrow(protein_normalized), coef=4, sort.by="logFC")
protein_ch_th <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(protein_ch_th) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(protein_ch_th)){
  if (protein_ch_th[i,"adj.p.val"]< 0.05 & protein_ch_th[i,"logfc"]>1){
    protein_ch_th[i,"significant"] <- "up" }
  else if (protein_ch_th[i,"adj.p.val"]< 0.05 & protein_ch_th[i,"logfc"]<(-1)){
    protein_ch_th[i,"significant"] <- "down" }
  else {
    protein_ch_th[i,"significant"] <- "ns" }
}

ggplot(protein_ch_th, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="Col_Heat vs tsnKO_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-1),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 1),          linetype=2, linewidth=0.5) + 
  xlim(-10,10) + ylim(0,10) +
  annotate( geom="text", x=-9, y= 9.5, label = as.character(nrow(protein_ch_th[protein_ch_th$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x= 9, y= 9.5, label = as.character(nrow(protein_ch_th[protein_ch_th$significant=="up",])), color="#e35023", fontface="bold")
  
# Export results ----

# We define an R list object in which all differentially expressed proteins 
# lists are stored and export it. This variable contains only the protein IDs.
protein_diff<-list(
  cc_ch_up=rownames(protein_cc_ch[protein_cc_ch$significant=="up",]),cc_ch_down=rownames(protein_cc_ch[protein_cc_ch$significant=="down",]),
  tc_th_up=rownames(protein_tc_th[protein_tc_th$significant=="up",]),tc_th_down=rownames(protein_tc_th[protein_tc_th$significant=="down",]),
  cc_tc_up=rownames(protein_cc_tc[protein_cc_tc$significant=="up",]),cc_tc_down=rownames(protein_cc_tc[protein_cc_tc$significant=="down",]),
  ch_th_up=rownames(protein_ch_th[protein_ch_th$significant=="up",]),ch_th_down=rownames(protein_ch_th[protein_ch_th$significant=="down",]))

save(protein_diff,file="../results/diff_lists_protein.RData")

# The next three lines of code allow us to get the gene symbol (ex: MKK4 for
# mitogen-activated protein kinase kinase 4) of Arabidopsis genes.
library(org.At.tair.db)
getsymbols <- org.At.tairSYMBOL
getsymbols <- as.list(getsymbols[mappedkeys(getsymbols)])

# The export_proteins_list writes tables containing all the results of our
# differential analysis: protein ID, -log10(adjusted p value), log fold-change,
# gene description and gene symbol
export_proteins_list(vpdata=protein_cc_ch,path="../results/protein/protein_col_ctl_col_heat")
export_proteins_list(vpdata=protein_tc_th,path="../results/protein/protein_tsnKO_ctl_tsnKO_heat")
export_proteins_list(vpdata=protein_cc_tc,path="../results/protein/protein_col_ctl_tsnko_ctl")
export_proteins_list(vpdata=protein_ch_th,path="../results/protein/protein_col_heat_tsnko_heat")

# Gene Ontology Over-Representation Analysis (ORA) -----

# We use the clusterProfiler library to perform enrichment analysis in ontologies
# Biological Process, Molecular Function and Cellular Component in the up- and
# down-regulated proteins.
library(clusterProfiler)

# Col_CTL vs Col_Heat
bp_ora_cc_ch_up <- enrichGO(gene=get_proteins(protein_diff$cc_ch_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_cc_ch_up <- enrichGO(gene=get_proteins(protein_diff$cc_ch_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_cc_ch_up <- enrichGO(gene=get_proteins(protein_diff$cc_ch_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_cc_ch_down <- enrichGO(gene=get_proteins(protein_diff$cc_ch_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_cc_ch_down <- enrichGO(gene=get_proteins(protein_diff$cc_ch_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_cc_ch_down <- enrichGO(gene=get_proteins(protein_diff$cc_ch_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# tsnKO_CTL vs tsnKO_Heat
bp_ora_tc_th_up <- enrichGO(gene=get_proteins(protein_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_tc_th_up <- enrichGO(gene=get_proteins(protein_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_tc_th_up <- enrichGO(gene=get_proteins(protein_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_tc_th_down <- enrichGO(gene=get_proteins(protein_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_tc_th_down <- enrichGO(gene=get_proteins(protein_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_tc_th_down <- enrichGO(gene=get_proteins(protein_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# Col_CTL vs tsnKO_CTL
bp_ora_cc_tc_up <- enrichGO(gene=get_proteins(protein_diff$cc_tc_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_cc_tc_up <- enrichGO(gene=get_proteins(protein_diff$cc_tc_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_cc_tc_up <- enrichGO(gene=get_proteins(protein_diff$cc_tc_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_cc_tc_down <- enrichGO(gene=get_proteins(protein_diff$cc_tc_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_cc_tc_down <- enrichGO(gene=get_proteins(protein_diff$cc_tc_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_cc_tc_down <- enrichGO(gene=get_proteins(protein_diff$cc_tc_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# Col_Heat vs tsnKO_Heat
bp_ora_ch_th_up <- enrichGO(gene=get_proteins(protein_diff$ch_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_ch_th_up <- enrichGO(gene=get_proteins(protein_diff$ch_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_ch_th_up <- enrichGO(gene=get_proteins(protein_diff$ch_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_ch_th_down <- enrichGO(gene=get_proteins(protein_diff$ch_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_ch_th_down <- enrichGO(gene=get_proteins(protein_diff$ch_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_ch_th_down <- enrichGO(gene=get_proteins(protein_diff$ch_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# Because the enrichGO function takes some time, we will save all the ORA
# results in an .RData file to avoid running it again.
save(bp_ora_cc_ch_up,  mf_ora_cc_ch_up,  cc_ora_cc_ch_up,
     bp_ora_cc_ch_down,mf_ora_cc_ch_down,cc_ora_cc_ch_down,
     bp_ora_tc_th_up,  mf_ora_tc_th_up,  cc_ora_tc_th_up,
     bp_ora_tc_th_down,mf_ora_tc_th_down,cc_ora_tc_th_down,
     bp_ora_cc_tc_up,  mf_ora_cc_tc_up,  cc_ora_cc_tc_up,
     bp_ora_cc_tc_down,mf_ora_cc_tc_down,cc_ora_cc_tc_down,
     bp_ora_ch_th_up,  mf_ora_ch_th_up,  cc_ora_ch_th_up,
     bp_ora_ch_th_down,mf_ora_ch_th_down,cc_ora_ch_th_down,
     file="../results/enrichment_results_protein.RData")

# For each of the gene ontologies, results of the ORA will be summarized in one
# dot plot showing the enriched GO terms across all comparisons.
# Only the Biological Process dotplot is commented, since the code is analogous
# to that used for Cellular Component and Molecular Function plots.

# Biological Process (BP) dotplots ----

# First, we take a look at the dot plots of each gene set, to ensure that no 
# important GO terms will be left behind in the final visualization.

options(enrichplot.colours=c("red","blue"))
dotplot( bp_ora_cc_ch_up,   title="Col_CTL vs Col_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_cc_ch_down, title="Col_CTL vs Col_Heat DOWN (BP)", showCategory=15)
dotplot( bp_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (BP)", showCategory=15)
dotplot( bp_ora_cc_tc_up,   title="Col_CTL vs tsnKO_CTL UP (BP)", showCategory=15)
dotplot( bp_ora_cc_tc_down, title="Col_CTL vs tsnKO_CT DOWN (BP)", showCategory=15)
dotplot( bp_ora_ch_th_up,   title="Col_Heat vs tsnKO_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_ch_th_down, title="Col_Heat vs tsnKO_Heat DOWN (BP)", showCategory=15)

# We extract the IDs of the top 5 GO terms enriched in each gene set
dpgoids <- sort(unique(c(bp_ora_cc_ch_up@result$ID[1:5],bp_ora_cc_ch_down@result$ID[1:5],
                         bp_ora_tc_th_up@result$ID[1:5],bp_ora_tc_th_down@result$ID[1:5],
                         bp_ora_cc_tc_down@result$ID[1:5],
                         bp_ora_ch_th_up@result$ID[1:5],bp_ora_ch_th_down@result$ID[1:5])))

# We define the eight gene sets we are studying
dpsets <- c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN",
            "Col/tsnKO CTL UP","Col/tsnKO CTL DOWN",
            "Col/tsnKO HS UP","Col/tsnKO HS DOWN")
dpsets <- factor(dpsets,levels=dpsets)

# We initialize two vectors: one for storing gene ratios and another for adjusted
# p values.
dpgeneratios <- c()
dppval <- c()

for (i.go in dpgoids){ # For each of the GO terms we will represent,
  
  # We extract the gene ratio
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_cc_ch_up@result[i.go,"GeneRatio"])))
  # and adjusted p value
  dppval <- c(dppval,bp_ora_cc_ch_up@result[i.go,"p.adjust"])
  # in each set
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_cc_ch_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_cc_ch_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_tc_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_tc_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_tc_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_tc_th_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_cc_tc_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_cc_tc_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_cc_tc_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_cc_tc_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_ch_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_ch_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=bp_ora_ch_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,bp_ora_ch_th_down@result[i.go,"p.adjust"])
  
}

# The dpgoids vector contains the IDs of GO terms (ex: GO:0010118). For the 
# visualization, we want to show the description of the GO term (ex: stomatal 
# movement), which we will store in the dpdescription vector.
dpdescription <- c()

for (i in 1:length(dpgoids)){ # For each GO term ID,
  # If the term is enriched in a specific set,
  if (!is.na(bp_ora_cc_ch_up@result[dpgoids[i],"Description"])){
    # we extract the description from the corresponding enrichResult object
    dpdescription[i] <- bp_ora_cc_ch_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(bp_ora_cc_ch_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- bp_ora_cc_ch_down@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(bp_ora_tc_th_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- bp_ora_tc_th_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(bp_ora_tc_th_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- bp_ora_tc_th_down@result[dpgoids[i],"Description"]
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
# dotplot( cc_ora_cc_ch_up,   title="Col_CTL vs Col_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_cc_ch_down, title="Col_CTL vs Col_Heat DOWN (cc)", showCategory=15)
dotplot( cc_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (cc)", showCategory=15)
dotplot( cc_ora_ch_th_up,   title="Col_Heat vs tsnKO_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_ch_th_down, title="Col_Heat vs tsnKO_Heat DOWN (cc)", showCategory=15)

# Dotplot with all comparisons

dpgoids <- sort(unique(c(cc_ora_cc_ch_down@result$ID[1:5],
                         cc_ora_tc_th_up@result$ID[1:5],cc_ora_tc_th_down@result$ID[1:5],
                         cc_ora_cc_tc_up@result$ID[1:5],cc_ora_cc_tc_down@result$ID[1:5],
                         cc_ora_ch_th_up@result$ID[1:5],cc_ora_ch_th_down@result$ID[1:5])))

dpsets <- c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN",
            "Col/tsnKO CTL UP","Col/tsnKO CTL DOWN",
            "Col/tsnKO HS UP","Col/tsnKO HS DOWN")
dpsets <- factor(dpsets,levels=dpsets)

dpgeneratios <- c()
dppval <- c()

for (i.go in dpgoids){
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_cc_ch_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_cc_ch_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_cc_ch_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_cc_ch_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_tc_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_tc_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_tc_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_tc_th_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_cc_tc_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_cc_tc_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_cc_tc_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_cc_tc_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_ch_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_ch_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_ch_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_ch_th_down@result[i.go,"p.adjust"])
  
}

dpdescription <- c()

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_cc_ch_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_cc_ch_down@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_tc_th_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_tc_th_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_tc_th_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_tc_th_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_ch_th_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_ch_th_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_ch_th_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_ch_th_down@result[dpgoids[i],"Description"]
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
dotplot( mf_ora_cc_ch_up,   title="Col_CTL vs Col_Heat UP (MF)", showCategory=15)
dotplot( mf_ora_cc_ch_down, title="Col_CTL vs Col_Heat DOWN (MF)", showCategory=15)
dotplot( mf_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (MF)", showCategory=15)
dotplot( mf_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (MF)", showCategory=15)
dotplot( mf_ora_ch_th_up,   title="Col_Heat vs tsnKO_Heat UP (MF)", showCategory=15)
dotplot( mf_ora_ch_th_down, title="Col_Heat vs tsnKO_Heat DOWN (MF)", showCategory=15)

# Dotplot with all comparisons

dpgoids <- sort(unique(c(mf_ora_cc_ch_up@result$ID[1:5],mf_ora_cc_ch_down@result$ID[1:5],
                         mf_ora_tc_th_up@result$ID[1:5],mf_ora_tc_th_down@result$ID[1:5],
                         mf_ora_cc_tc_down@result$ID[1:5],
                         mf_ora_ch_th_up@result$ID[1:5],mf_ora_ch_th_down@result$ID[1:5])))

dpsets <- c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN",
            "Col/tsnKO CTL UP","Col/tsnKO CTL DOWN",
            "Col/tsnKO HS UP","Col/tsnKO HS DOWN")
dpsets <- factor(dpsets,levels=dpsets)

dpgeneratios <- c()
dppval <- c()

for (i.go in dpgoids){
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_cc_ch_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_cc_ch_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_cc_ch_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_cc_ch_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_tc_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_tc_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_tc_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_tc_th_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_cc_tc_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_cc_tc_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_cc_tc_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_cc_tc_down@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_ch_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_ch_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=mf_ora_ch_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,mf_ora_ch_th_down@result[i.go,"p.adjust"])
  
}

dpdescription <- c()

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_cc_ch_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_cc_ch_up@result[dpgoids[i],"Description"]
  }
}

for (i in 1:length(dpgoids)){
  if (!is.na(mf_ora_cc_ch_down@result[dpgoids[i],"Description"])){
    dpdescription[i] <- mf_ora_cc_ch_down@result[dpgoids[i],"Description"]
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
