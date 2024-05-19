### PHOSPHOPROTEOMIC DIFFERENTIAL ANALYSIS ###
###    AND GENE ONTOLOGY TERM ENRICHMENT   ###

# First, we source the script where custom functions are defined.
source("functions.R")

# We will use limma for differential analysis and ggplot2 to generate all plots.
library(limma)
library(ggplot2)

# We load the preprocessed phosphoproteomic data, exported from the 
# data_preprocessing_phospho.R script.
load(file="../data/phospho_data.RData")

# We specify the experimental design (4 groups; 4 replicates each)
exp_design <- model.matrix(~ -1+factor(c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4)))
# In this script, we will often refer to a group using only two letters.
# The first one represents genotype (c/t) and the second one, treatment (c/h)
# cc: col_ctl , ch: col_heat , tc: tsnko_ctl , th: tsnko_heat
colnames(exp_design) <- c("cc","ch","tc","th")

linearfit <- lmFit(phospho_normalized, exp_design)

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
temp_res <- topTable(limma_results, number=nrow(phospho_normalized), coef=1, sort.by="logFC")
# The log fold-change and adjusted p values are saved in a data.frame, called
# phospho_cc_ch (phosphoproteomics, col_ctl vs col_heat)
phospho_cc_ch <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
# The rownames of this dataframe are phosphosite names.
rownames(phospho_cc_ch) <- rownames(temp_res)
# We remove the temporal variable.
rm(temp_res)

# We add a column indicating whether a phosphosite is up-phosphorylated (up),
# down-phosphorylated (down) or not differentially phosphorylated (ns).

for (i in 1:nrow(phospho_cc_ch)){
  
  # If the change is significant (adjusted p < 0.05) and logfc > 3,
  # the phosphosite is up-phosphorylated.
  if (phospho_cc_ch[i,"adj.p.val"]< 0.05 & phospho_cc_ch[i,"logfc"]>3){
    phospho_cc_ch[i,"significant"] <- "up" }
  
  # If the change is significant (adjusted p < 0.05) and logfc < -3,
  # the phosphosite is down-phosphorylated.
  else if (phospho_cc_ch[i,"adj.p.val"]< 0.05 & phospho_cc_ch[i,"logfc"]<(-3)){
    phospho_cc_ch[i,"significant"] <- "down" }
  
  # In any other case, we consider that the change is not significant.
  else {
    phospho_cc_ch[i,"significant"] <- "ns" }
}

# Volcano plot
ggplot(phospho_cc_ch, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="Col_CTL vs Col_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-3),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 3),          linetype=2, linewidth=0.5) + 
  xlim(-17,17) + ylim(0,11) +
  # We will show the number of differentially phosphorylated sites and 
  # differentially phosphorylated phosphoproteins (up and down)
  # To calculate the number of proteins in which phosphosites are found, we use
  # substr(x,1,9) to extract the protein code and count how many are unique.
  # (length(unique(-)))
  annotate( geom="text", x=-15, y= 10.7, label = as.character(nrow(phospho_cc_ch[phospho_cc_ch$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x=-15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_cc_ch[phospho_cc_ch$significant=="down",]),1,9)))),color="#058dfc", size=3) +
  annotate( geom="text", x= 15, y= 10.7, label = as.character(nrow(phospho_cc_ch[phospho_cc_ch$significant=="up",])), color="#e35023", fontface="bold") +
  annotate( geom="text", x= 15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_cc_ch[phospho_cc_ch$significant=="up",]),1,9)))),color="#e35023", size=3)

# tsnKO_CTl vs tsnKO_Heat -------

temp_res <- topTable(limma_results, number=nrow(phospho_normalized), coef=2, sort.by="logFC")
phospho_tc_th <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(phospho_tc_th) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(phospho_tc_th)){
  if (phospho_tc_th[i,"adj.p.val"]< 0.05 & phospho_tc_th[i,"logfc"]>3){
    phospho_tc_th[i,"significant"] <- "up" }
  else if (phospho_tc_th[i,"adj.p.val"]< 0.05 & phospho_tc_th[i,"logfc"]<(-3)){
    phospho_tc_th[i,"significant"] <- "down" }
  else {
    phospho_tc_th[i,"significant"] <- "ns" }
}

ggplot(phospho_tc_th, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="tsnKO_CTL vs tsnKO_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-3),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 3),          linetype=2, linewidth=0.5) + 
  xlim(-17,17) + ylim(0,11) +
  annotate( geom="text", x=-15, y= 10.7, label = as.character(nrow(phospho_tc_th[phospho_tc_th$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x=-15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_tc_th[phospho_tc_th$significant=="down",]),1,9)))),color="#058dfc", size=3) +
  annotate( geom="text", x= 15, y= 10.7, label = as.character(nrow(phospho_tc_th[phospho_tc_th$significant=="up",])), color="#e35023", fontface="bold") +
  annotate( geom="text", x= 15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_tc_th[phospho_tc_th$significant=="up",]),1,9)))),color="#e35023", size=3)

# Col_CTL vs tsnKO_CTL -------

temp_res <- topTable(limma_results, number=nrow(phospho_normalized), coef=3, sort.by="logFC")
phospho_cc_tc <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(phospho_cc_tc) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(phospho_cc_tc)){
  if (phospho_cc_tc[i,"adj.p.val"]< 0.05 & phospho_cc_tc[i,"logfc"]>3){
    phospho_cc_tc[i,"significant"] <- "up" }
  else if (phospho_cc_tc[i,"adj.p.val"]< 0.05 & phospho_cc_tc[i,"logfc"]<(-3)){
    phospho_cc_tc[i,"significant"] <- "down" }
  else {
    phospho_cc_tc[i,"significant"] <- "ns" }
}

ggplot(phospho_cc_tc, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="Col_CTL vs tsnKO_CTL") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-3),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 3),          linetype=2, linewidth=0.5) + 
  xlim(-17,17) + ylim(0,11) +
  annotate( geom="text", x=-15, y= 10.7, label = as.character(nrow(phospho_cc_tc[phospho_cc_tc$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x=-15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_cc_tc[phospho_cc_tc$significant=="down",]),1,9)))),color="#058dfc", size=3) +
  annotate( geom="text", x= 15, y= 10.7, label = as.character(nrow(phospho_cc_tc[phospho_cc_tc$significant=="up",])), color="#e35023", fontface="bold") +
  annotate( geom="text", x= 15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_cc_tc[phospho_cc_tc$significant=="up",]),1,9)))),color="#e35023", size=3)

# Col_Heat vs tsnKO_Heat -------

temp_res <- topTable(limma_results, number=nrow(phospho_normalized), coef=4, sort.by="logFC")
phospho_ch_th <- data.frame(logfc=temp_res$logFC, adj.p.val=temp_res$adj.P.Val)
rownames(phospho_ch_th) <- rownames(temp_res)
rm(temp_res)

for (i in 1:nrow(phospho_ch_th)){
  if (phospho_ch_th[i,"adj.p.val"]< 0.05 & phospho_ch_th[i,"logfc"]>3){
    phospho_ch_th[i,"significant"] <- "up" }
  else if (phospho_ch_th[i,"adj.p.val"]< 0.05 & phospho_ch_th[i,"logfc"]<(-3)){
    phospho_ch_th[i,"significant"] <- "down" }
  else {
    phospho_ch_th[i,"significant"] <- "ns" }
}

ggplot(phospho_ch_th, aes(x=logfc, y=-log10(adj.p.val))) + 
  geom_point( aes(color=significant), shape=16, alpha=0.4) + 
  scale_color_manual( values=c("#058dfc", "grey40", "#e35023") ) +
  labs (x="\nlog2fc", y="-log10 (adjusted p value)\n", title="Col_Heat vs tsnKO_Heat") + 
  theme( plot.title=element_text( hjust=0.5, face="bold" ), legend.position="none") +
  geom_hline(color="grey20", aes(yintercept=-log10(0.05)),linetype=2, linewidth=0.5) + 
  geom_vline(color="grey20", aes(xintercept=-3),          linetype=2, linewidth=0.5) +  
  geom_vline(color="grey20", aes(xintercept= 3),          linetype=2, linewidth=0.5) + 
  xlim(-17,17) + ylim(0,11) +
  annotate( geom="text", x=-15, y= 10.7, label = as.character(nrow(phospho_ch_th[phospho_ch_th$significant=="down",])), color="#058dfc", fontface="bold") +
  annotate( geom="text", x=-15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_ch_th[phospho_ch_th$significant=="down",]),1,9)))),color="#058dfc", size=3) +
  annotate( geom="text", x= 15, y= 10.7, label = as.character(nrow(phospho_ch_th[phospho_ch_th$significant=="up",])), color="#e35023", fontface="bold") +
  annotate( geom="text", x= 15, y=  9.7, label = as.character(length(unique(substr(rownames(phospho_ch_th[phospho_ch_th$significant=="up",]),1,9)))),color="#e35023", size=3)

# Export results ----

# We define an R list object in which all differentially phosphorylated sites 
# lists are stored and export it. This variable contains only the phosphosites
# names.
phospho_diff<-list(
  cc_ch_up=rownames(phospho_cc_ch[phospho_cc_ch$significant=="up",]),cc_ch_down=rownames(phospho_cc_ch[phospho_cc_ch$significant=="down",]),
  tc_th_up=rownames(phospho_tc_th[phospho_tc_th$significant=="up",]),tc_th_down=rownames(phospho_tc_th[phospho_tc_th$significant=="down",]),
  cc_tc_up=rownames(phospho_cc_tc[phospho_cc_tc$significant=="up",]),cc_tc_down=rownames(phospho_cc_tc[phospho_cc_tc$significant=="down",]),
  ch_th_up=rownames(phospho_ch_th[phospho_ch_th$significant=="up",]),ch_th_down=rownames(phospho_ch_th[phospho_ch_th$significant=="down",]))

save(phospho_diff,file="../results/diff_lists_phospho.RData")

# The next three lines of code allow us to get the gene symbol (ex: MKK4 for
# mitogen-activated protein kinase kinase 4) of Arabidopsis genes.
library(org.At.tair.db)
getsymbols <- org.At.tairSYMBOL
getsymbols <- as.list(getsymbols[mappedkeys(getsymbols)])

# The export_phosphosites_list writes tables containing all the results of our
# differential analysis: protein ID, modified residue, multiplicity, -log10(adjusted 
# p value), log fold-change, gene description and gene symbol
export_phosphosites_list(vpdata=phospho_cc_ch,path="../results/phospho/phospho_col_ctl_col_heat")
export_phosphosites_list(vpdata=phospho_tc_th,path="../results/phospho/phospho_tsnKO_ctl_tsnKO_heat")
export_phosphosites_list(vpdata=phospho_cc_tc,path="../results/phospho/phospho_col_ctl_tsnko_ctl")
export_phosphosites_list(vpdata=phospho_ch_th,path="../results/phospho/phospho_col_heat_tsnko_heat")

# Gene Ontology Over-Representation Analysis (ORA) -----

# From this point on, analyses will not use the lists of the Col_CTL vs tsnKO_CTL
# comparison, since we detected a very small numbers of differentially 
# phosphorylated sites.

# We use the clusterProfiler library to perform enrichment analysis in ontologies
# Biological Process, Molecular Function and Cellular Component in the up- and
# down-phosphorylated proteins.
library(clusterProfiler)

# Col_CTL vs Col_Heat
bp_ora_cc_ch_up <- enrichGO(gene=get_proteins(phospho_diff$cc_ch_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_cc_ch_up <- enrichGO(gene=get_proteins(phospho_diff$cc_ch_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_cc_ch_up <- enrichGO(gene=get_proteins(phospho_diff$cc_ch_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_cc_ch_down <- enrichGO(gene=get_proteins(phospho_diff$cc_ch_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_cc_ch_down <- enrichGO(gene=get_proteins(phospho_diff$cc_ch_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_cc_ch_down <- enrichGO(gene=get_proteins(phospho_diff$cc_ch_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# tsnKO_CTL vs tsnKO_Heat
bp_ora_tc_th_up <- enrichGO(gene=get_proteins(phospho_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_tc_th_up <- enrichGO(gene=get_proteins(phospho_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_tc_th_up <- enrichGO(gene=get_proteins(phospho_diff$tc_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_tc_th_down <- enrichGO(gene=get_proteins(phospho_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_tc_th_down <- enrichGO(gene=get_proteins(phospho_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_tc_th_down <- enrichGO(gene=get_proteins(phospho_diff$tc_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# Col_Heat vs tsnKO_Heat
bp_ora_ch_th_up <- enrichGO(gene=get_proteins(phospho_diff$ch_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_ch_th_up <- enrichGO(gene=get_proteins(phospho_diff$ch_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_ch_th_up <- enrichGO(gene=get_proteins(phospho_diff$ch_th_up),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")
bp_ora_ch_th_down <- enrichGO(gene=get_proteins(phospho_diff$ch_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "BP")
mf_ora_ch_th_down <- enrichGO(gene=get_proteins(phospho_diff$ch_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "MF")
cc_ora_ch_th_down <- enrichGO(gene=get_proteins(phospho_diff$ch_th_down),OrgDb="org.At.tair.db",keyType = "TAIR",ont = "CC")

# Because the enrichGO function takes some time, we will save all the ORA
# results in an .RData file to avoid running it again.
save(bp_ora_cc_ch_up,  mf_ora_cc_ch_up,  cc_ora_cc_ch_up,
     bp_ora_cc_ch_down,mf_ora_cc_ch_down,cc_ora_cc_ch_down,
     bp_ora_tc_th_up,  mf_ora_tc_th_up,  cc_ora_tc_th_up,
     bp_ora_tc_th_down,mf_ora_tc_th_down,cc_ora_tc_th_down,
     bp_ora_ch_th_up,  mf_ora_ch_th_up,  cc_ora_ch_th_up,
     bp_ora_ch_th_down,mf_ora_ch_th_down,cc_ora_ch_th_down,
     file="../results/enrichment_results_phospho.RData")

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
dotplot( bp_ora_ch_th_up,   title="Col_Heat vs tsnKO_Heat UP (BP)", showCategory=15)
dotplot( bp_ora_ch_th_down, title="Col_Heat vs tsnKO_Heat DOWN (BP)", showCategory=15)

# We extract the IDs of the top 5 GO terms enriched in each gene set
dpgoids <- sort(unique(c(bp_ora_cc_ch_up@result$ID[1:5],bp_ora_cc_ch_down@result$ID[1:5],
                         bp_ora_tc_th_up@result$ID[1:5],bp_ora_tc_th_down@result$ID[1:5],
                         bp_ora_ch_th_up@result$ID[1:5],bp_ora_ch_th_down@result$ID[1:5])))
# We removed the terms "nuclear-transcribed mRNA catabolic process, 
# deadenylation-dependent decay" and "RNA splicing, via transesterification 
# reactions with bulged adenosine as nucleophile", since they are redundant
# with other terms that will be represented and make the final plot look worse.
dpgoids <- dpgoids[-c(1,3)]

# We define the six gene sets we are studying
dpsets <- c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN","Col/tsnKO UP","Col/tsnKO DOWN")
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
  # If the term is enriched in the Col UP set,
  if (!is.na(bp_ora_cc_ch_up@result[dpgoids[i],"Description"])){
    # we extract the description from the corresponding enrichResult object
    dpdescription[i] <- bp_ora_cc_ch_up@result[dpgoids[i],"Description"]
  }
}

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
dotplot( cc_ora_cc_ch_up,   title="Col_CTL vs Col_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_cc_ch_down, title="Col_CTL vs Col_Heat DOWN (cc)", showCategory=15)
dotplot( cc_ora_tc_th_up,   title="tsnKO_CTL vs tsnKO_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_tc_th_down, title="tsnKO_CTL vs tsnKO_Heat DOWN (cc)", showCategory=15)
dotplot( cc_ora_ch_th_up,   title="Col_Heat vs tsnKO_Heat UP (cc)", showCategory=15)
dotplot( cc_ora_ch_th_down, title="Col_Heat vs tsnKO_Heat DOWN (cc)", showCategory=15)

# Dotplot with all comparisons

dpgoids <- sort(unique(c(cc_ora_cc_ch_up@result$ID[1:5],cc_ora_cc_ch_down@result$ID[1:5],
                         cc_ora_tc_th_up@result$ID[1:5],cc_ora_tc_th_down@result$ID[1:5],
                         cc_ora_ch_th_up@result$ID[1:5],cc_ora_ch_th_down@result$ID[1:5])))

dpsets <- c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN","Col/tsnKO UP","Col/tsnKO DOWN")
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
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_ch_th_up@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_ch_th_up@result[i.go,"p.adjust"])
  
  dpgeneratios <- c(dpgeneratios,eval(parse(text=cc_ora_ch_th_down@result[i.go,"GeneRatio"])))
  dppval <- c(dppval,cc_ora_ch_th_down@result[i.go,"p.adjust"])
  
}

dpdescription <- c()

for (i in 1:length(dpgoids)){
  if (!is.na(cc_ora_cc_ch_up@result[dpgoids[i],"Description"])){
    dpdescription[i] <- cc_ora_cc_ch_up@result[dpgoids[i],"Description"]
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
# dotplot( mf_ora_ch_th_down, title="Col_Heat vs tsnKO_Heat DOWN (MF)", showCategory=15)

# Dotplot with all comparisons

dpgoids <- sort(unique(c(mf_ora_cc_ch_up@result$ID[1:5],mf_ora_cc_ch_down@result$ID[1:5],
                         mf_ora_tc_th_up@result$ID[1:5],mf_ora_tc_th_down@result$ID[1:5],
                         mf_ora_ch_th_up@result$ID[1:5],mf_ora_ch_th_down@result$ID[1:5])))
dpgoids <- dpgoids[-c(13,18)]

dpsets <- c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN","Col/tsnKO UP","Col/tsnKO DOWN")
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
