### PROTEOMIC DATA PREPROCESSING ###

# First, we source the script where custom functions are defined.
source("functions.R")

library(preprocessCore)

# Later on, values will be generated in a pseudorandom manner. To allow for
# reproducibility of our analysis, we set a seed.
set.seed(-611543874)

protein_raw <- read.delim("../data/proteinGroups_proteome_2KO.txt", header=TRUE)
# The imported data frame consists of 8184 rows (proteins)
# Our experimental design includes 4 groups (4 replicates each): 
#   - Col_CTL (Columbia seedlings in control conditions)
#   - tsnKO_CTL (tsn1 tsn2 knockout mutant seedlings in control conditions)
#   - Col_Heat (Columbia seedlings in heat stress conditions)
#   - tsnKO_Heat (tsn1 tsn2 knockout mutant seedlings in heat stress conditions)
# Heat stressed Arabidopsis seedlings were incubated at 42ºC for 1 hour.

# Using a custom function, we will discard rows corresponding to proteins only
# identified by site, as well as those flagged as reverse sequences or
# potential contaminants. Numbers next to each command show how many proteins
# are retained after the correspondent filter.
protein_raw <- filter_categorical(protein_raw,"Only.identified.by.site") #7930
protein_raw <- filter_categorical(protein_raw,"Reverse") #7852
protein_raw <- filter_categorical(protein_raw,"Potential.contaminant") #7842

# First, we extract the names of columns with label-free quantification (LFQ)
# intensity values. All of them start with "LFQ".
intensity_cols <- colnames(protein_raw)[startsWith(colnames(protein_raw),"LFQ")]
# Then, we apply the custom log2_transformation function, which log transforms
# the columns specified in the second argument.
protein_log <- log2_transformation(protein_raw,intensity_cols)

# After log-transformation, we filter out rows with a high number of missing 
# values. Our filter_missing function will discard sequences if it does not 
# have three valid values in at least one of the groups.
protein_log <- filter_missing(protein_log,intensity_cols)
# 6599 pass the filter applied.
rownames(protein_log) <- protein_log$Majority.protein.IDs
# We name the rows of the protein_log matrix to make data accession easier.

protein_imputed <- impute_missing_values(protein_log,intensity_cols)
# The impute_missing_values performs data imputation of missing values using the
# Missing Not At Random (MNAR) method. Imputed values are drawn from a normal
# distribution with centered around mean-1.8*sd and with a standard deviation of
# 0.3*sd, where mean and sd represent the mean and standard deviation of the 
# intensities in each sample.

# After imputation, we quantile normalize data and name columns and rows 
# accordingly.
protein_normalized <- as.data.frame(normalize.quantiles(as.matrix(protein_imputed[,intensity_cols])))
colnames(protein_normalized) <- intensity_cols
rownames(protein_normalized) <- protein_imputed$Majority.protein.IDs

# Additionally, we generate the protein_intensities matrix, in which log2
# transformation is reversed. These values will be used in the processing
# of phosphoproteomic data.
protein_intensities <- revert_log2_transformation(protein_normalized,intensity_cols)

# We remove rows 5165 and 6356, which correspond to TSN1 and TSN2 (knocked-out
# proteins in the mutant).
protein_normalized <- protein_normalized[-c(5165,6356),]

# We export the data we will use in subsequent analyses, saving the following
# data frames in an .RData file
save(protein_log,protein_normalized,protein_intensities,file="../data/protein_data.RData")

# PRINCIPAL COMPONENT ANALYSIS
library(FactoMineR)
library(ggplot2)

# We perform principal component analysis using the protein_normalized matrix
# and plot the corresponding scores plot.
pca.result <- PCA(t(protein_normalized),graph=FALSE, scale.unit=F)

pca.result$ind$coord <- as.data.frame(pca.result$ind$coord)
pca.result$ind$coord[,"Group"] <- rep(c("Col_CTL","Col_Heat","tsnKO_CTL","tsnKO_Heat"),each=4)

ggplot(pca.result$ind$coord, aes(x=Dim.1, y=Dim.2)) + 
  geom_point(shape=21, aes(colour=Group),fill="black",size=2, stroke=2.5) +
  scale_color_manual( values=c("#e35a42","#d67c09","#32b85c","#4e6ccf") ) + 
  ylim (-50,30) + xlim(-30,50) +
  labs (title="PCA (proteomics)", x="PC1 (22.3%)", y="PC2 (16.1%)") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))
# Together, principal components 1 and 2 explain 38.4 % of variability.
# The four groups are clearly separated, though Col_Heat samples show
# the least reproducibility.