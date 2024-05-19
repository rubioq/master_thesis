### MULTI-OMIC VENN DIAGRAMS ###

# We first load the lists of 
load("../results/diff_lists_phospho.RData") # differentially phosphorylated sites,
load("../results/diff_lists_protein.RData") # differentially expressed proteins and
load("../results/diff_lists_transc.RData") # differentially expressed genes
# Also, we source the functions.R script, in which custom functions are defined.
source("functions.R")
# We will use ggvenn for visualization
library(ggvenn)

# We will unite all lists in a single object, so we add a prefix indicating
# the omic to which each list refers.
names(phospho_diff) <- paste0("phospho_",names(phospho_diff))
names(protein_diff) <- paste0("protein_",names(protein_diff))
names(transc_diff) <- paste0("transc_",names(transc_diff))

# In each list, we need to extract the protein names.
for (i in 1:length(phospho_diff)){
  phospho_diff[[i]] <- get_proteins(phospho_diff[[i]])
  protein_diff[[i]] <- get_proteins(protein_diff[[i]])
  transc_diff[[i]] <- get_proteins(transc_diff[[i]])
}

# Now that we have made the changes described above, we can create a single list
# with all the information we need.
venn_data <- c(phospho_diff,protein_diff,transc_diff)

venn_palette <- c("#4548D0","#2F6129","#CB181D")

# We generate the Venn diagrams of interest.
ggvenn(venn_data, columns = c("phospho_cc_ch_up","protein_cc_ch_up","transc_wc_wh_up"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.5, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("phospho_cc_ch_down","protein_cc_ch_down","transc_wc_wh_down"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.5, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("phospho_tc_th_up","protein_tc_th_up","transc_tc_th_up"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.5, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("phospho_tc_th_down","protein_tc_th_down","transc_tc_th_down"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.5, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("phospho_ch_th_up","protein_ch_th_up","transc_wh_th_up"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.5, set_name_color = venn_palette, text_size=6)

ggvenn(venn_data, columns = c("phospho_ch_th_down","protein_ch_th_down","transc_wh_th_down"),stroke_size = 0.3,
       fill_color = venn_palette, stroke_color="white", show_percentage=F,
       fill_alpha=0.5, set_name_color = venn_palette, text_size=6)
