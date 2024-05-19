### PHOSPHOPROTEOMIC DOWNSTREAM ###
###          ANALYSIS           ###

# First, we source the script where custom functions are defined.
source("functions.R")

# In this script, we will make use of the following databases
library(org.At.tair.db)
library(GO.db)

# All plots will be generated using ggplot2
library(ggplot2)

# We load the results of our phosphoproteomic analysis, which we generated
# using data_preprocessing_phospho.R
load("../data/phospho_data.RData")
# or phospho_downstream_analysis.R
load("../results/diff_lists_phospho.RData")

# Cellular Component barplots ----

# The goal of this section is to represent the gene ratios of specific GO terms
# in each set. 

# The db_data variable contains the GO IDs associated with each protein ID.
db_data <- as.list(org.At.tairGO[mappedkeys(org.At.tairGO)])
# We also define the ancestor_data variable, which will allow us to extract
# the so-called ancestor terms of each term. For example, the GO term "nuclear
# body" is an ancestor of terms such as "Cajal body" and "nuclear speck". 
# Taking this into account is necessary because a protein might be annotated
# with the term "Cajal body" in the database and we need to take into account
# that it is indeed a "nuclear body".
ancestor_data <- as.list(GOCCANCESTOR)

# We define a vector containing the protein IDs of all the phosphoproteins 
# in our study.
my_proteins <- get_proteins(rownames(phospho_imputed))

# CC - Nucleolus ----

# We will analyze the GO terms "nucleolus", "nuclear body" and "cytoplasmic
# stress granule", using analogous code in each case. Therefore, only the
# first section of code is commented.

# Using the custom function go_term_present, we will create the 
# nucleolus_go_info variable. It consists of a named vector where names
# are protein IDs and the values are boolean (TRUE if the protein ID
# is associated with the nucleolus and FALSE if it is not).
# This vector contains information regarding all phosphoproteins detected
# in our assay. 
nucleolus_go_info <- go_term_present(my_proteins,"GO:0005730",db_data,ancestor_data)

# The barplot_data is a data frame, formatted in a ggplot2-friendly way, which
# contains gene ratios of the term "nucleolus" in the background (whole
# Arabidopsis proteome) and each differentially phosphorylated proteins set.
# The gene ratio of the background was extracted from enrichResult objects.
barplot_data <- data.frame(set = factor(c("Background",
                                        "Col UP","Col DOWN",
                                        "tsnKO UP","tsnKO DOWN",
                                        "Col/tsnKO UP", "Col/tsnKO DOWN"),
                                      levels=c("Background","Col UP","Col DOWN","tsnKO UP",
                                               "tsnKO DOWN","Col/tsnKO UP", "Col/tsnKO DOWN")),
                           
                           count = c(500/26290,
                                    sum(nucleolus_go_info[get_proteins(phospho_diff$cc_ch_up)])   / length(get_proteins(phospho_diff$cc_ch_up)),
                                    sum(nucleolus_go_info[get_proteins(phospho_diff$cc_ch_down)]) / length(get_proteins(phospho_diff$cc_ch_down)),
                                    sum(nucleolus_go_info[get_proteins(phospho_diff$tc_th_up)])   / length(get_proteins(phospho_diff$tc_th_up)),
                                    sum(nucleolus_go_info[get_proteins(phospho_diff$tc_th_down)]) / length(get_proteins(phospho_diff$tc_th_down)),
                                    sum(nucleolus_go_info[get_proteins(phospho_diff$ch_th_up)])   / length(get_proteins(phospho_diff$ch_th_up)),
                                    sum(nucleolus_go_info[get_proteins(phospho_diff$ch_th_down)]) / length(get_proteins(phospho_diff$ch_th_down))))

# We generate the barplot, where heights represent gene ratios
ggplot(data=barplot_data,aes(x=set,y=count)) + geom_bar(stat="identity",fill=c("#404040",rep(c("#ed7753","#5fbee3"),3))) +
  labs(title= "Nucleolus", y="Gene Ratio\n", x=NULL) +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + ylim(0,0.08) + 
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,color="black")) +
  # and add the gene counts as numbers
  annotate(geom = "text", x = 1, y = barplot_data[1,2]+0.006, label = "500", size=3.3) +
  # For each set, we simply show how many protein IDs have a "TRUE" value 
  annotate(geom = "text", x = 2, y = barplot_data[2,2]+0.006, label = as.character(sum(nucleolus_go_info[get_proteins(phospho_diff$cc_ch_up)])), size=3.3) + 
  annotate(geom = "text", x = 3, y = barplot_data[3,2]+0.006, label = as.character(sum(nucleolus_go_info[get_proteins(phospho_diff$cc_ch_down)])), size=3.3) +
  annotate(geom = "text", x = 4, y = barplot_data[4,2]+0.006, label = as.character(sum(nucleolus_go_info[get_proteins(phospho_diff$tc_th_up)])), size=3.3) +
  annotate(geom = "text", x = 5, y = barplot_data[5,2]+0.006, label = as.character(sum(nucleolus_go_info[get_proteins(phospho_diff$tc_th_down)])), size=3.3) +
  annotate(geom = "text", x = 6, y = barplot_data[6,2]+0.006, label = as.character(sum(nucleolus_go_info[get_proteins(phospho_diff$ch_th_up)])), size=3.3) +
  annotate(geom = "text", x = 7, y = barplot_data[7,2]+0.006, label = as.character(sum(nucleolus_go_info[get_proteins(phospho_diff$ch_th_down)])), size=3.3)

# CC - Nuclear body ----

nuclear_body_go_info <- go_term_present(my_proteins,"GO:0016604",db_data,ancestor_data)

barplot_data <- data.frame(set = factor(c("Background",
                                          "Col UP","Col DOWN",
                                          "tsnKO UP","tsnKO DOWN",
                                          "Col/tsnKO UP", "Col/tsnKO DOWN"),
                                        levels=c("Background","Col UP","Col DOWN","tsnKO UP",
                                                 "tsnKO DOWN","Col/tsnKO UP", "Col/tsnKO DOWN")),
                           
                           count = c(105/26290,
                                     sum(nuclear_body_go_info[get_proteins(phospho_diff$cc_ch_up)])   / length(get_proteins(phospho_diff$cc_ch_up)),
                                     sum(nuclear_body_go_info[get_proteins(phospho_diff$cc_ch_down)]) / length(get_proteins(phospho_diff$cc_ch_down)),
                                     sum(nuclear_body_go_info[get_proteins(phospho_diff$tc_th_up)])   / length(get_proteins(phospho_diff$tc_th_up)),
                                     sum(nuclear_body_go_info[get_proteins(phospho_diff$tc_th_down)]) / length(get_proteins(phospho_diff$tc_th_down)),
                                     sum(nuclear_body_go_info[get_proteins(phospho_diff$ch_th_up)])   / length(get_proteins(phospho_diff$ch_th_up)),
                                     sum(nuclear_body_go_info[get_proteins(phospho_diff$ch_th_down)]) / length(get_proteins(phospho_diff$ch_th_down))))

ggplot(data=barplot_data,aes(x=set,y=count)) + geom_bar(stat="identity",fill=c("#404040",rep(c("#ed7753","#5fbee3"),3))) +
  labs(title= "Nuclear body", y="Gene Ratio\n", x=NULL) +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + ylim(0,0.045) + 
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,color="black")) +
  annotate(geom = "text", x = 1, y = barplot_data[1,2]+0.0035, label = "105", size=3.3) +
  annotate(geom = "text", x = 2, y = barplot_data[2,2]+0.0035, label = as.character(sum(nuclear_body_go_info[get_proteins(phospho_diff$cc_ch_up)])), size=3.3) + 
  annotate(geom = "text", x = 3, y = barplot_data[3,2]+0.0035, label = as.character(sum(nuclear_body_go_info[get_proteins(phospho_diff$cc_ch_down)])), size=3.3) +
  annotate(geom = "text", x = 4, y = barplot_data[4,2]+0.0035, label = as.character(sum(nuclear_body_go_info[get_proteins(phospho_diff$tc_th_up)])), size=3.3) +
  annotate(geom = "text", x = 5, y = barplot_data[5,2]+0.0035, label = as.character(sum(nuclear_body_go_info[get_proteins(phospho_diff$tc_th_down)])), size=3.3) +
  annotate(geom = "text", x = 6, y = barplot_data[6,2]+0.0035, label = as.character(sum(nuclear_body_go_info[get_proteins(phospho_diff$ch_th_up)])), size=3.3) +
  annotate(geom = "text", x = 7, y = barplot_data[7,2]+0.0035, label = as.character(sum(nuclear_body_go_info[get_proteins(phospho_diff$ch_th_down)])), size=3.3)

# CC - Cytoplasmic tress granule ----

cyt_sg_go_info <- go_term_present(my_proteins,"GO:0010494",db_data,ancestor_data)

barplot_data <- data.frame(set = factor(c("Background",
                                          "Col UP","Col DOWN",
                                          "tsnKO UP","tsnKO DOWN",
                                          "Col/tsnKO UP", "Col/tsnKO DOWN"),
                                        levels=c("Background","Col UP","Col DOWN","tsnKO UP",
                                                 "tsnKO DOWN","Col/tsnKO UP", "Col/tsnKO DOWN")),
                           
                           count = c(152/26290,
                                     sum(cyt_sg_go_info[get_proteins(phospho_diff$cc_ch_up)])   / length(get_proteins(phospho_diff$cc_ch_up)),
                                     sum(cyt_sg_go_info[get_proteins(phospho_diff$cc_ch_down)]) / length(get_proteins(phospho_diff$cc_ch_down)),
                                     sum(cyt_sg_go_info[get_proteins(phospho_diff$tc_th_up)])   / length(get_proteins(phospho_diff$tc_th_up)),
                                     sum(cyt_sg_go_info[get_proteins(phospho_diff$tc_th_down)]) / length(get_proteins(phospho_diff$tc_th_down)),
                                     sum(cyt_sg_go_info[get_proteins(phospho_diff$ch_th_up)])   / length(get_proteins(phospho_diff$ch_th_up)),
                                     sum(cyt_sg_go_info[get_proteins(phospho_diff$ch_th_down)]) / length(get_proteins(phospho_diff$ch_th_down))))

ggplot(data=barplot_data,aes(x=set,y=count)) + geom_bar(stat="identity",fill=c("#404040",rep(c("#ed7753","#5fbee3"),3))) +
  labs(title= "Cytoplasmic stress granule", y="Gene Ratio\n", x=NULL) +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) + ylim(0,0.11) + 
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,color="black")) +
  annotate(geom = "text", x = 1, y = barplot_data[1,2]+0.008, label = "152", size=3.3) +
  annotate(geom = "text", x = 2, y = barplot_data[2,2]+0.008, label = as.character(sum(cyt_sg_go_info[get_proteins(phospho_diff$cc_ch_up)])), size=3.3) + 
  annotate(geom = "text", x = 3, y = barplot_data[3,2]+0.008, label = as.character(sum(cyt_sg_go_info[get_proteins(phospho_diff$cc_ch_down)])), size=3.3) +
  annotate(geom = "text", x = 4, y = barplot_data[4,2]+0.008, label = as.character(sum(cyt_sg_go_info[get_proteins(phospho_diff$tc_th_up)])), size=3.3) +
  annotate(geom = "text", x = 5, y = barplot_data[5,2]+0.008, label = as.character(sum(cyt_sg_go_info[get_proteins(phospho_diff$tc_th_down)])), size=3.3) +
  annotate(geom = "text", x = 6, y = barplot_data[6,2]+0.008, label = as.character(sum(cyt_sg_go_info[get_proteins(phospho_diff$ch_th_up)])), size=3.3) +
  annotate(geom = "text", x = 7, y = barplot_data[7,2]+0.008, label = as.character(sum(cyt_sg_go_info[get_proteins(phospho_diff$ch_th_down)])), size=3.3)

# LLPS analysis - PPID ----

# We obtain the in silico predictions of intrinsic disorder that were generated
# using IUPRED. Results are imported as a data frame where each row represents
# a protein and each column, the amino acid position. For example, column 75 of
# the row "AT1G01120.1" shows the score of the 75th amino acid in that protein.
iupred_df <- read.csv(file="../results/iupred_results.csv",header=FALSE,
                      col.names = paste0("V",0:5400), fill = TRUE, row.names = 1)

# Now, we will calculate the PPID (Predicted Percentage of Intrinsic Disorder)
# of all proteins in the dataframe.

# ppid_info is the variable where we will store PPID values.
ppid_info <- c()

for (i in 1:nrow(iupred_df)){ # For each protein,
  # We extract the information in that row (scores of all amino acids)
  i.vector <- as.numeric(as.vector(iupred_df[i,])) 
  # and remove NAs (positions where column index > protein length)
  i.vector <- i.vector[!is.na(i.vector)]
  # The PPID is calculated as the proportion of numbers that are above 0.5,
  # that is, we consider amino acids with score > 0.5 to be disordered and
  # obtain the proportion of how many residues are disordered.
  ppid_info[i] <- sum(i.vector>0.5)/length(i.vector)
  # We name each element with the corresponding protein ID
  names(ppid_info)[i] <- rownames(iupred_df)[i]
}
# Since we calculated proportions, we need to multiply the values by 100.
ppid_info <- ppid_info*100
# We remove i.vector and i, since we will not use them again.
rm(i.vector,i)

# We prepare the data to plot it using ggplot2, extracting the necessary PPID
# values from the ppid_info vector.
background_bp <- data.frame(set="Background",value=ppid_info)
cc_ch_up_bp <- data.frame(set="Col UP",value=ppid_info[get_protein_ids(phospho_diff$cc_ch_up)])
cc_ch_down_bp <- data.frame(set="Col DOWN",value=ppid_info[get_protein_ids(phospho_diff$cc_ch_down)])
tc_th_up_bp <- data.frame(set="tsnKO UP",value=ppid_info[get_protein_ids(phospho_diff$tc_th_up)])
tc_th_down_bp <- data.frame(set="tsnKO DOWN",value=ppid_info[get_protein_ids(phospho_diff$tc_th_down)])
ch_th_up_bp <- data.frame(set="Col/tsnKO UP",value=ppid_info[get_protein_ids(phospho_diff$ch_th_up)])
ch_th_down_bp <- data.frame(set="Col/tsnKO DOWN",value=ppid_info[get_protein_ids(phospho_diff$ch_th_down)])

boxplot_data <- rbind(background_bp,cc_ch_up_bp,cc_ch_down_bp,tc_th_up_bp,tc_th_down_bp,ch_th_up_bp,ch_th_down_bp)
boxplot_data$set <- factor(boxplot_data$set,levels=unique(boxplot_data$set))

# PPIDs in each set do not have a normal distribution
shapiro.test(background_bp$value) 
shapiro.test(cc_ch_up_bp$value) 
shapiro.test(cc_ch_down_bp$value) 
shapiro.test(tc_th_up_bp$value) 
shapiro.test(tc_th_down_bp$value) 
shapiro.test(ch_th_up_bp$value) 
shapiro.test(ch_th_down_bp$value) 

# We run non-parametric tests comparing the PPIDs in each set with the background.
wilcox.test(background_bp$value,cc_ch_up_bp$value) # ***
wilcox.test(background_bp$value,cc_ch_down_bp$value) # ***
wilcox.test(background_bp$value,tc_th_up_bp$value) # ***
wilcox.test(background_bp$value,tc_th_down_bp$value) # ***
wilcox.test(background_bp$value,ch_th_up_bp$value) # ***
wilcox.test(background_bp$value,ch_th_down_bp$value) # ***

# We remove the variables we do not need for plotting
rm(background_bp,cc_ch_up_bp,cc_ch_down_bp,tc_th_up_bp,tc_th_down_bp,ch_th_up_bp,ch_th_down_bp)

ggplot(boxplot_data, aes(x=set, y=value, fill=set)) +
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,color="black")) + 
  theme(legend.position = "none") + 
  geom_boxplot(width=0.75,color="#b30e0e", linewidth=0.6, fill="#dc7e86") +
  labs(y="% IDR", x=NULL, title="PPID") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  ylim(0,106) +
  annotate(geom = "text", x = 2, y = 105, label = "***",size=4) +
  annotate(geom = "text", x = 3, y = 105, label = "***",size=4) +
  annotate(geom = "text", x = 4, y = 105, label = "***",size=4) +
  annotate(geom = "text", x = 5, y = 105, label = "***",size=4) +
  annotate(geom = "text", x = 6, y = 105, label = "***",size=4) +
  annotate(geom = "text", x = 7, y = 105, label = "***",size=4)

# LLPS analysis - PSP Score ----

# We import the PSP (Phase Separating Protein) score for the whole Arabidopsis
# thaliana proteome. Data is imported as a data frame but then converted to
# a named vector (names are protein IDs), where values represent the PSP Score
# of each protein.
psp_df <- read.csv(file="../results/psp_score_results.csv", header= T)
psp_info <- psp_df$Score
names(psp_info) <- rownames(psp_df)
# Since the psp_info vector contains the same information as the psp_df, we 
# remove the latter.
rm(psp_df)

background_bp <- data.frame(set="Background",value=psp_info)

# In the psp_info variable, some values are shared among different proteins,
# which makes it difficult for us to extract PSP Scores of proteins if its
# corresponding entry includes more than one protein. To solve this, we will
# run the following code.

for (i in 1:length(psp_info)){ # For each protein,
  if (nchar(names(psp_info)[i])>13){ 
    # If the name contains over 13 characters, it means that several proteins
    # share the same PSP Score. If that is the case, we will separate them.
    ids <- strsplit(names(psp_info)[i],";")[[1]]
    for (j in 1:length(ids)){ 
      # and then for each ID, we will add a new element so that we can access
      # the PSP Score of any protein.
      psp_info <- c(psp_info,unname(psp_info[i]))
      names(psp_info)[length(psp_info)] <- ids[j]
    }
  }
}

# We prepare the data to plot it using ggplot2, extracting the necessary PSP 
# Scores from the psp_info vector.
cc_ch_up_bp <- data.frame(set="Col UP",value=psp_info[get_protein_ids(phospho_diff$cc_ch_up)])
cc_ch_down_bp <- data.frame(set="Col DOWN",value=psp_info[get_protein_ids(phospho_diff$cc_ch_down)])
tc_th_up_bp <- data.frame(set="tsnKO UP",value=psp_info[get_protein_ids(phospho_diff$tc_th_up)])
tc_th_down_bp <- data.frame(set="tsnKO DOWN",value=psp_info[get_protein_ids(phospho_diff$tc_th_down)])
ch_th_up_bp <- data.frame(set="Col/tsnKO UP",value=psp_info[get_protein_ids(phospho_diff$ch_th_up)])
ch_th_down_bp <- data.frame(set="Col/tsnKO DOWN",value=psp_info[get_protein_ids(phospho_diff$ch_th_down)])

boxplot_data <- rbind(background_bp,cc_ch_up_bp,cc_ch_down_bp,tc_th_up_bp,tc_th_down_bp,ch_th_up_bp,ch_th_down_bp)
boxplot_data$set <- factor(boxplot_data$set,levels=unique(boxplot_data$set))

# PSP Scores in each set do not have a normal distribution
shapiro.test(background_bp$value) 
shapiro.test(cc_ch_up_bp$value) 
shapiro.test(cc_ch_down_bp$value) 
shapiro.test(tc_th_up_bp$value) 
shapiro.test(tc_th_down_bp$value) 
shapiro.test(ch_th_up_bp$value) 
shapiro.test(ch_th_down_bp$value) 

# We run non-parametric tests comparing the PSP SCores in each set with the 
# background.
wilcox.test(background_bp$value,cc_ch_up_bp$value) # ***
wilcox.test(background_bp$value,cc_ch_down_bp$value) # ***
wilcox.test(background_bp$value,tc_th_up_bp$value) # ***
wilcox.test(background_bp$value,tc_th_down_bp$value) # ***
wilcox.test(background_bp$value,ch_th_up_bp$value) # ***
wilcox.test(background_bp$value,ch_th_down_bp$value) # ***
rm(background_bp,cc_ch_up_bp,cc_ch_down_bp,tc_th_up_bp,tc_th_down_bp,ch_th_up_bp,ch_th_down_bp)

ggplot(boxplot_data, aes(x=set, y=value, fill=set)) +
  theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,color="black")) + 
  theme(legend.position = "none") + 
  geom_boxplot(width=0.75,color="#b30e0e", linewidth=0.6, fill="#dc7e86") +
  labs(y="PSP Score", x=NULL, title="PSP Score") + 
  theme(plot.title = element_text(hjust = 0.5, face="bold"))+ 
  ylim(0,1.06) +
  annotate(geom = "text", x = 2, y = 1.05, label = "***",size=4) +
  annotate(geom = "text", x = 3, y = 1.05, label = "***",size=4) +
  annotate(geom = "text", x = 4, y = 1.05, label = "***",size=4) +
  annotate(geom = "text", x = 5, y = 1.05, label = "***",size=4) +
  annotate(geom = "text", x = 6, y = 1.05, label = "***",size=4) +
  annotate(geom = "text", x = 7, y = 1.05, label = "***",size=4)

# LLPS analysis - Phosphosites in IDRs ----

# For each set, we will use the custom function "phosphosites_in_idr" to determine
# whether a phosphosite is in a disordered region (IUPRED score > 0.5) or not.
barplot_vector <- 100 * c(as.vector(prop.table(table(phosphosites_in_idr(phospho_diff$cc_ch_up)))),
                          as.vector(prop.table(table(phosphosites_in_idr(phospho_diff$cc_ch_down)))),
                          as.vector(prop.table(table(phosphosites_in_idr(phospho_diff$tc_th_up)))),
                          as.vector(prop.table(table(phosphosites_in_idr(phospho_diff$tc_th_down)))),
                          as.vector(prop.table(table(phosphosites_in_idr(phospho_diff$ch_th_up)))),
                          as.vector(prop.table(table(phosphosites_in_idr(phospho_diff$ch_th_down)))))

barplot_sets <- rep(c("Col UP","Col DOWN","tsnKO UP","tsnKO DOWN",
                      "Col/tsnKO UP", "Col/tsnKO DOWN"), each=4)
# We remove one Col/tsnKO UP because it has only three of the four possibilities
# defined ("Shares protein with another phosphosite" = 0)
barplot_sets <- barplot_sets[-20]
barplot_sets <- factor(barplot_sets,levels=unique(barplot_sets))

barplot_disorder <- rep(c("Disordered","Ordered","Shares IDR with another phosphosite","Shares protein with another phosphosite"),6)
# Again, this is to adjust for the fact that there are 0 phosphosites in Col/tsnKO
# UP in the group "Shares protein with another phosphosite"
barplot_disorder <- barplot_disorder[-20]
barplot_disorder <- factor(barplot_disorder,levels=c("Shares IDR with another phosphosite",
                                                     "Shares protein with another phosphosite",
                                                     "Disordered",
                                                     "Ordered"))
# We include all the data in a variable that is compatible with ggplot2 and 
# remove the ones we will not need anymore.
barplot_data <- data.frame(set=barplot_sets,
                           disorder=barplot_disorder,
                           percentage=barplot_vector)
rm(barplot_sets,barplot_disorder,barplot_vector)

ggplot(data=barplot_data, aes(x=set, y=percentage, fill=disorder)) + geom_bar(stat="identity", width=0.7) +
  scale_fill_manual(values=c("#EE4C5A","#F87B56","#F3AD6A","#0D5861")) + 
  labs(title= "Phosphosites in IDRs", y="% phosphosites", x=NULL) +
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, color = "black"), axis.text.y=element_text(color="black"))

# Venn diagrams - SG proteins ----

library(ggvenn)

# The sg_pool_hs.txt file contains a list of putative stress granule proteins,
# since they interacted with TSN2, RBP47 or RBGD2/4 in Heat Stress conditions
# (see Solis-Miranda et al., 2023; doi: 10.1093/plcell/koad127)
sg_proteins <- read.table("../data/sg_pool_hs.txt")[,1]

# We define a list with the sets we are interested in representing: stress
# granules proteins and the differentially phosphorylated proteins we found.
venn_data <- list(`SG proteins`=sg_proteins,
                  `Col UP`=get_proteins(phospho_diff$cc_ch_up),
                  `Col DOWN`=get_proteins(phospho_diff$cc_ch_down),
                  `tsnKO UP`=get_proteins(phospho_diff$tc_th_up),
                  `tsnKO DOWN`=get_proteins(phospho_diff$tc_th_down),
                  `Col/tsnKO UP`=get_proteins(phospho_diff$ch_th_up),
                  `Col/tsnKO DOWN`=get_proteins(phospho_diff$ch_th_down))

venn_palette <- c("#e35023","#148733","#058dfc")

ggvenn(venn_data, columns = c("Col UP","SG proteins","Col DOWN"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.6, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("tsnKO UP","SG proteins","tsnKO DOWN"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.6, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("Col/tsnKO UP","SG proteins","Col/tsnKO DOWN"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.6, set_name_color = venn_palette, text_size=6)

# Clustering - Heatmap ----

# In this section, we will represent the phosphorylation levels of phosphosites
# in a heatmap, but we are only interested in phosphosites that, according to
# our analysis, are differentially phosphorylated in at least one comparison.

# The diffpps vector contains all the identified DPSs
diffpps <- unique(c(phospho_diff[[1]],phospho_diff[[2]],phospho_diff[[3]],phospho_diff[[4]],
                    phospho_diff[[5]],phospho_diff[[6]],phospho_diff[[7]],phospho_diff[[8]]))

# To make different phosphosites comparable among them, we need to calculate
# their z-scores.
z_scores <- scale(t(phospho_normalized[diffpps,]))

col <- colorRampPalette(c("blue", "white", "red"))(255)
length(diffpps) # A total of 2608 phosphosites will be visualized.

# We rename the rows, since the z_scores matrix retains the names originally
# in phospho_normalized (ex: Intensity.Col_CTL_1 instead of just Col_CTL_1)
rownames(z_scores) <- c("Col_CTL_1","Col_CTL_2","Col_CTL_3","Col_CTL_4",
                        "Col_Heat_1","Col_Heat_2","Col_Heat_3","Col_Heat_4",
                        "tsnKO_CTL_1","tsnKO_CTL_2","tsnKO_CTL_3","tsnKO_CTL_4",
                        "tsnKO_Heat_1","tsnKO_Heat_2","tsnKO_Heat_3","tsnKO_Heat_4")
myheatmap <- heatmap(t(z_scores), col=col, scale="none",keep.dendro=TRUE)

# Using the results generated, we will group the phosphosites into 6clusters
clusters <- cutree(as.hclust(myheatmap$Rowv), k=6)  # k = 6 clusters

# cluster_ids is a list with 6 elements, where each element is the vector of
# phosphosites names in each cluster
cluster_ids <- list()
for (i in 1:max(clusters)){
  cluster_ids[[i]] <- names(clusters[clusters==i])
}

cluster_palette <- c("#32b85c","#e35a42","#d67c09","#9940b3","#4e6ccf","#f0c713")

# Using the rowcols vector will allow us to represent the clusters in the heatmap
rowcols <- rep(NA,length(clusters))
for (i in 1:length(clusters)){
  rowcols[i] <- cluster_palette[clusters[i]]
}

heatmap(t(z_scores), col=col, scale="none",keep.dendro=TRUE, RowSideColors = rowcols)
rm(myheatmap)

# Clustering - Boxplots ----

# We reorder the clusters so that they appear from A to E when represented
# in the heatmap
clusters <- list(clusterA=cluster_ids[[2]],
                 clusterB=cluster_ids[[3]],
                 clusterC=cluster_ids[[6]],
                 clusterD=cluster_ids[[1]],
                 clusterE=cluster_ids[[5]],
                 clusterF=cluster_ids[[4]])

boxplot_colors <- c("#e35a42","#d67c09","#f0c713","#32b85c","#4e6ccf","#9940b3")
cluster_letters <- c("A","B","C","D","E","F")

for (i in 1:length(clusters)){ # For each of the 6 clusters,
  
  # we calculate the mean z-score of the corresponding phosphosites
  boxplot_data <- rbind(data.frame(cluster=paste0(cluster_letters[i]," (",length(clusters[[i]]),")"),set="Col_CTL",value=unname(colMeans(z_scores[1:4,clusters[[i]]]))),
                        data.frame(cluster=paste0(cluster_letters[i]," (",length(clusters[[i]]),")"),set="tsnKO_CTL",value=unname(colMeans(z_scores[9:12,clusters[[i]]]))),
                        data.frame(cluster=paste0(cluster_letters[i]," (",length(clusters[[i]]),")"),set="Col_Heat",value=unname(colMeans(z_scores[5:8,clusters[[i]]]))),
                        data.frame(cluster=paste0(cluster_letters[i]," (",length(clusters[[i]]),")"),set="tsnKO_Heat",value=unname(colMeans(z_scores[13:16,clusters[[i]]]))))
  boxplot_data$set <- factor(boxplot_data$set,levels=c("Col_CTL","tsnKO_CTL","Col_Heat","tsnKO_Heat"))
  
  # and generate a boxplot using ggplot2
  print(ggplot(boxplot_data, aes(x=set, y=value, fill=set)) + geom_jitter(color=boxplot_colors[i],position=position_jitter(0.3)) +
    theme(legend.position = "none", axis.text.x=element_text(color="black"), axis.text.y=element_text(color="black")) + 
    geom_boxplot(width=0.7, fill=NA,color="black", linewidth=0.6) +
    labs(title= boxplot_data[1,"cluster"], y="mean Z-score", x=NULL) + 
    theme(plot.title = element_text(hjust = 0.5, face="bold")))
   
}

# The next three lines of code allow us to get the gene symbol (ex: MKK4 for
# mitogen-activated protein kinase kinase 4) of Arabidopsis genes.
library(org.At.tair.db)
getsymbols <- org.At.tairSYMBOL
getsymbols <- as.list(getsymbols[mappedkeys(getsymbols)])

# Using the custom export_cluster_list, we generate txt files containing the
# lists of phosphosites in each cluster.
export_cluster_list(clusters$clusterA,"../results/phospho/cluster_a.txt")
export_cluster_list(clusters$clusterB,"../results/phospho/cluster_b.txt")
export_cluster_list(clusters$clusterC,"../results/phospho/cluster_c.txt")
export_cluster_list(clusters$clusterD,"../results/phospho/cluster_d.txt")
export_cluster_list(clusters$clusterE,"../results/phospho/cluster_e.txt")
export_cluster_list(clusters$clusterF,"../results/phospho/cluster_f.txt")

# Clustering - Upset plot ----

library(UpSetR)

# To study the overlap between sets, w first need to get the protein IDs from
# the phosphosites names.
upset_data <- list(clusterA=get_proteins(clusters$clusterA),
                   clusterB=get_proteins(clusters$clusterB),
                   clusterC=get_proteins(clusters$clusterC),
                   clusterD=get_proteins(clusters$clusterD),
                   clusterE=get_proteins(clusters$clusterE),
                   clusterF=get_proteins(clusters$clusterF))

# We take a first look at how the groups overlap 
upset(fromList(upset_data),order.by="freq",nsets=6,nintersects=10)

# and then define which overlaps will be plotted
upset(fromList(upset_data),order.by = "freq",
      keep.order = T, sets=names(upset_data),
      sets.bar.color=rev(boxplot_colors),text.scale = 2,
      intersections =
        list(
          list("clusterF","clusterE","clusterD","clusterC","clusterB","clusterA"),
          list("clusterE","clusterF"),
          list("clusterA","clusterB"),
          list("clusterE","clusterA"),
          list("clusterA","clusterF"),
          list("clusterF","clusterB"),
          list("clusterB","clusterE"),
          list("clusterF","clusterD"),
          list("clusterD","clusterA"),
          list("clusterA","clusterF","clusterE"),
          list("clusterB","clusterF","clusterE"),
          list("clusterC","clusterB")        )) 

# Clustering - Motif analysis ----

# For motif analysis, we will extract the sequence windows (that is, the amino 
# acid sequences around each phosphosite) in clusters C and D, which are the ones
# that behave most differently in wild-type and tsnKO seedlings.

# We import the protein sequences of all Arabidopsis proteins, which were 
# downloaded from the TAIR database. aa_seqs is a dataframe with only one
# column, which has the amino acid sequences, and the name of each row
# corresponds to the protein ID.
aa_seqs <- read.delim("../data/arabidopsis_aa_seqs.csv",sep=",",row.names = 1)

seq_windows_C <- c()

for (i in 1:length(clusters$clusterC)){ # For each element in the cluster,
  
  # For example, for phosphosite AT2G41900.1_510_1
  # we extract the protein ID, which is the text before the first underscore
  # ex: AT2G41900.1
  protein <- strsplit(clusters$clusterC[i],"_")[[1]][1]
  # We extract the position of the phosphosite within the protein, which
  # is indicated under the first underscore (ex: 510)
  position <- as.numeric(strsplit(clusters$clusterC[i],"_")[[1]][2])
  # We extract the surrounding amino acids (sequence window of size 15)
  seq_windows_C <- c(seq_windows_C,substr(aa_seqs[protein,"Protein.Sequence"],position-7,position+7))
}

# We use similar code to extract the sequence windows corresponding to cluster D
# phosphosites.
seq_windows_D <- c()

for (i in 1:length(clusters$clusterC)){
  protein <- strsplit(clusters$clusterD[i],"_")[[1]][1]
  position <- as.numeric(strsplit(clusters$clusterC[i],"_")[[1]][2])
  seq_windows_D <- c(seq_windows_C,substr(aa_seqs[protein,"Protein.Sequence"],position-7,position+7))
}

write.table(seq_windows_C,quote=F,row.names=F,col.names=F,file="../results/seq_windows_clusterC.txt")
write.table(seq_windows_D,quote=F,row.names=F,col.names=F,file="../results/seq_windows_clusterD.txt")

# Venn diagrams - MPK6 substrates ----

library(ggvenn)

# The putative_mpk6_substrates.txt file contains a list of MPK6 substrates 
# identified by Wang et al. (2020a); doi: 10.1073/pnas. 1919901117
mpk6_substrates <- read.table("../data/putative_mpk6_substrates.txt")[,1]

# We define a list with the sets we are interested in representing: MPK6
# substrates and the differentially phosphorylated proteins we found.
venn_data <- list(`MPK6 substrates`=mpk6_substrates,
                  `Col UP`=get_proteins(phospho_diff$cc_ch_up),
                  `Col DOWN`=get_proteins(phospho_diff$cc_ch_down),
                  `tsnKO UP`=get_proteins(phospho_diff$tc_th_up),
                  `tsnKO DOWN`=get_proteins(phospho_diff$tc_th_down),
                  `Col/tsnKO UP`=get_proteins(phospho_diff$ch_th_up),
                  `Col/tsnKO DOWN`=get_proteins(phospho_diff$ch_th_down))

venn_palette <- c("#e35023","#148733","#058dfc")

ggvenn(venn_data, columns = c("Col UP","MPK6 substrates","Col DOWN"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.6, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("tsnKO UP","MPK6 substrates","tsnKO DOWN"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.6, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("Col/tsnKO UP","MPK6 substrates","Col/tsnKO DOWN"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.6, set_name_color = venn_palette, text_size=6)
