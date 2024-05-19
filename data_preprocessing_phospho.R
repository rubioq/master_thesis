### PHOSPHOPROTEOMIC DATA PREPROCESSING ###

# First, we source the script where custom functions are defined.
source("functions.R")

library(preprocessCore)
library(matrixStats)

# We load the processed proteomic data, which will be used in this script in
# a normalization step.
load("../data/protein_data.RData")

# For imputation, values will be generated in a pseudorandom manner. To allow 
# for reproducibility of our analysis, we set a seed.
set.seed(367684115)

phospho_raw <- read.delim("../data/Phospho (STY)Sites.txt", header=TRUE)
# The imported data frame consists of 24908 rows
# Our experimental design is identical to that of the proteomic data
# For each of the following 4 groups, 4 replicates were analyzed:
#   - Col_CTL (Columbia seedlings in control conditions)
#   - tsnKO_CTL (tsn1 tsn2 knockout mutant seedlings in control conditions)
#   - Col_Heat (Columbia seedlings in heat stress conditions)
#   - tsnKO_Heat (tsn1 tsn2 knockout mutant seedlings in heat stress conditions)
# Heat stressed Arabidopsis seedlings were incubated at 42ºC for 1 hour.

# Using a custom function, we will discard rows corresponding to phosphosites
# detected as reverse sequences or potential contaminants. Numbers next to each
# command show how many rows are retained after the correspondent filter.
phospho_raw <- filter_categorical(phospho_raw,"Reverse") # 24660
phospho_raw <- filter_categorical(phospho_raw,"Potential.contaminant") # 24637

# We calculate the number of proteins in which the identified phosphosites are
# located. Protein IDs are extracted from the fasta headers column (found before
# the vertical bar). 
length(unique(sapply(strsplit(phospho_raw$Fasta.headers,"\\|"), `[[`, 1))) 
# 6060 proteins

# We also obtain the number of peptides.
length(unique(phospho_raw$Sequence.window)) 
# 24631 peptides

# Phosphosites are often divided into different groups depending on the 
# localization probability, i.e. how confident we are that a specific
# residue is modified. 
classes_table <- table(cut(phospho_raw$Localization.prob,breaks = c(0,0.25,0.5,0.75,1)))
classes_table <- classes_table[-1] # Remove phosphosites with very low probabilities.
prop.table(classes_table)*100 # Calculate percentages multiplying proportions by 100.
# We show the results in a pie plot.
pie(classes_table, labels = c("16.10% (3894)", "15.62% (3777)", "68.29% (16517)"), 
    border=NA, col=c("#FFFFBB","#D9F0A3","#78B679"),init.angle=180,radius=1,cex=1.3)
# To remove phosphosites with localization probability under 0.75, we use
# a custom function.
phospho_raw <- filter_class_i(phospho_raw) 
# The resulting dataframe has 16517 rows (=number of class I phosphosites)

# We extract the names of intensity columns. 
intensity_cols <- get_phospho_intensity_cols(phospho_raw)
# Then, we apply the custom log2_transformation function, which log transforms
# the columns specified in the second argument.
phospho_log <- log2_transformation(phospho_raw,intensity_cols)

# In the phospho_log dataframe, there are three intensities in each sample for
# each row. This corresponds to phosphopeptides that are phosphorylated in one,
# two or three positions. The number of modifications within one peptide is
# called "multiplicity". We want to expand our data frame so that these three
# cases are in three different rows.
phospho_expanded <- expand_site_table(phospho_log,intensity_cols)

# The add_phosphosite_name creates a new column with an identifier for each
# phosphosite, which is a string consisting of the protein ID, the modified
# residue position and the multiplicity. For example, AT1G01050.2_24_1 
# represents the phosphorylation of protein AT1G01050.2 in the 24th residue,
# with a multiplicity of 1 (no other residue modified in the same peptide)
phospho_expanded <- add_phosphosite_name(phospho_expanded)

# We filter out rows with a high number of missing values. Our filter_missing 
# function will only keep phosphosites that have a valid value in three samples
# of at least one group. 
samples <- extract_sample_names(intensity_cols)
phospho_expanded <- filter_missing(phospho_expanded,samples) 
# 10737 rows are kept after this filter.

phospho_imputed <- impute_missing_values(phospho_expanded,samples)
# The impute_missing_values performs data imputation of missing values using the
# Missing Not At Random (MNAR) method. Imputed values are drawn from a normal
# distribution with centered around mean-1.8*sd and with a standard deviation of
# 0.3*sd, where mean and sd represent the mean and standard deviation of the 
# intensities in each sample.
rownames(phospho_imputed) <- phospho_imputed$Phosphosite.name
# Phosphosite names are used as the row names of our data.

# To ensure that values are representative of the phosphorylation status and 
# not a consequence of the corresponding protein level, we normalize 
# phosphoproteomic data in relation to proteomic intensities.
# We use the medians of intensities for this normalization.
protein_medians <- get_proteome_medians(protein_intensities)
phospho_corrected <- normalization_to_proteome(phospho_imputed,protein_medians,samples)

# We quantile normalize data and name columns and rows accordingly
phospho_normalized <- normalize.quantiles(as.matrix(phospho_corrected[,samples]))
colnames(phospho_normalized) <- samples
rownames(phospho_normalized) <- phospho_corrected$Phosphosite.name

# We are not interested in comparing the phosphorylation status of TSN residues.
# Because mutants do not express TSN, normalization to protein intensities results
# in falsely high phosphorylation values in tsnKO mutants. We remove the 
# phosphosites found in TSN1 and TSN2.
tsn_phospho_data <- phospho_normalized[substr(rownames(phospho_normalized),1,9)=="AT5G61780" | substr(rownames(phospho_normalized),1,9)=="AT5G07350",]
phospho_normalized <- phospho_normalized[!(substr(rownames(phospho_normalized),1,9)=="AT5G61780" | substr(rownames(phospho_normalized),1,9)=="AT5G07350"),]

# We export the data we will use in differential and downstream analyses in a
# .RData file.
save(phospho_log,phospho_imputed,phospho_normalized,tsn_phospho_data,file="../data/phospho_data.RData")

# Pie plot showing the frequency of each modified amino acid (Serine, Threonine
# or Tyrosine)
prop.table(table(phospho_imputed[,"Amino.acid"]))*100
pie(table(phospho_imputed[,"Amino.acid"]) , labels = c("87.43% (9387)","12.37% (1328)","0.20% (22)"), border=NA, 
    col=c("#e7793d","#f0ad5f","#e42626"),main="Phosphorylated residue",radius=1,cex=1.3)

# PRINCIPAL COMPONENT ANALYSIS
library(FactoMineR)
library(ggplot2)

# We perform principal component analysis using the phospho_normalized matrix
# and plot the corresponding scores plot.
pca.result <- PCA(t(phospho_normalized),graph=FALSE, scale.unit=F)

pca.result$ind$coord <- as.data.frame(pca.result$ind$coord)
pca.result$ind$coord[,"Group"] <- rep(c("Col_CTL","Col_Heat","tsnKO_CTL","tsnKO_Heat"),each=4)

ggplot(pca.result$ind$coord, aes(x=Dim.1, y=Dim.2)) + 
  geom_point(shape=21, aes(colour=Group),fill="black",size=2, stroke=2.5) +
  scale_color_manual( values=c("#e35a42","#d67c09","#32b85c","#4e6ccf") ) + 
  ylim (-150,200) + xlim(-200,200) +
  labs (title="PCA (phosphoproteomics)", x="PC1 (30.5%)", y="PC2 (10.2%)") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
# Together, principal components 1 and 2 explain 40.7 % of variability.
# The four groups are clearly separated and heat stress is the main source
# of variability.