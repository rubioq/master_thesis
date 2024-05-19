### CUSTOM FUNCTIONS USED FOR THE ANALYSIS ###

# In this script, we define functions that will be useful for the analysis of
# phosphoproteomic, proteomic and/or transcriptomic data.

filter_categorical <- function(dataframe,column,remove_plus=TRUE){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   column: name of the categorical (boolean) column of interest
  #   remove_plus: if TRUE, rows where the column contains "+" will be removed
  #                else, only rows where the column contains "+" will be kept
  
  # Depending on remove_plus, we obtain a vector with the positions we wish to
  # keep in the output of this function
              if (remove_plus){ positions_keep <- which(dataframe[,column] != "+")  }
              else { positions_keep <- which(dataframe[,column] == "+") }
  
  # Using positions_keep, we can directly return the filtered dataframe
              return(dataframe[positions_keep,])
            }

filter_class_i <- function(dataframe,column="Localization.prob"){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   column: name of the column containing localization probabilities
  
  # In the positions_keep vector, we store the indexes of rows where the 
  # localization probability is greater than 0.75 (class I phosphosites)
  positions_keep <- which(dataframe[,column] > 0.75)
  
  # Using positions_keep, we can directly return the filtered dataframe
  return (dataframe[positions_keep,])
}

get_phospho_intensity_cols <- function(dataframe){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object  
  
  # This function simply returns the names of intensity columns, all of which 
  # start with "Intensity." and end with three underscores followed by the
  # multiplicity (1, 2 or 3)
  result <- c(
    colnames(dataframe)[startsWith(colnames(dataframe),"Intensity.")&endsWith(colnames(dataframe),"___1")],
    colnames(dataframe)[startsWith(colnames(dataframe),"Intensity.")&endsWith(colnames(dataframe),"___2")],
    colnames(dataframe)[startsWith(colnames(dataframe),"Intensity.")&endsWith(colnames(dataframe),"___3")]
    )
  
  return(result)
}

extract_sample_names <- function(phospho_intensities){
  
  # INPUT
  # phospho_intensities: a vector of MaxQuant-style column names
  
  result <- phospho_intensities[1:((length(phospho_intensities))/3)]
  for (i in 1:length(result)){
    # Using the substr function, we remove the ___X suffix
    result[i] <- substr(result[i],1,nchar(result[i])-4) 
  }
  return(result)
}

log2_transformation <- function(dataframe,columns){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   columns: a vector containing the names of columns on which the logarithmic
  #            transformation will be performed

              for (j in columns){ # For each column,
                # we apply log2 to the whole column
                dataframe[,j] <- log2(dataframe[,j])
              }
  
  
              for (j in columns){ # In every column
                for (i in 1:nrow(dataframe)){ # In every row (i.e., in each cell)
                  if (dataframe[i,j] == -Inf){ 
                    # If the value is -Inf (intensity was 0),
                    # we replace the cell with a NA (not assigned)
                    dataframe[i,j] <- NA
                  }
                }
              }
  
            # The output of the function is the same dataframe used as input,
            # except now the intensities have been log transformed.
              return(dataframe)
            }

revert_log2_transformation <- function(dataframe,columns){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   columns: a vector containing the names of columns on which the logarithmic
  #            transformation will be reversed
  
  for (j in columns){ # For each column,
    # we revert the log2 transformation
    dataframe[,j] <- 2^(dataframe[,j])
  }
  return(dataframe)
}

expand_site_table <- function(dataframe,intensity_columns,columns_keep=c("Amino.acid","Charge","Reverse","Potential.contaminant","Localization.prob","PEP","Score","Score.for.localization","Proteins","Positions.within.proteins","Leading.proteins","Protein","Sequence.window","Fasta.headers")) {

  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   intensity_columns: a vector containing the names of intensity columns
  #   columns_keep: a vector with the names of additional columns that will
  #                 be kept in the output
  
  # This function has the same purpose as the "Expand site table" option in
  # Perseus. Basically, multiplicities of each phosphosite are restructured
  # so that each one is its own row.
  
  # First, we extract sample names (taking the first column names and removing
  # the ___X suffix)
  samples <- intensity_columns[1:((length(intensity_columns))/3)]
  for (i in 1:length(samples)){
    samples[i] <- substr(samples[i],1,nchar(samples[i])-4) 
  }
  
  # We initialize the result (output) of this function as a dataframe.
  result <- data.frame()

  for (n in 1:nrow(dataframe)){ # For each row,
    for (multiplicity in 1:3){ # we repeat it three times (one per multiplicity)
      # First, we will add the columns that are not intensities (columns_keep)
      result <- rbind(result,dataframe[n,columns_keep])
    }
  }
  
  # We initialize a second dataframe where intensities will be temporarily stored
  # before binding it with the result dataframe.
  intensities <- data.frame()
  
  for (i in 1:nrow(dataframe)){ # Again, for each row
    for (j in samples) { # each column
      for (multiplicity in 1:3){ # and each multiplicity
        # we get the index of the corresponding row in the result dataframe,
        # taking into consideration the desired final structure
        currenti <- 1+(i-1)*3+(multiplicity-1) 
        # We pass the intensity value into the intensities dataframe
        intensities[currenti,j] <- dataframe[i,paste0(j,"___",multiplicity)]
      }
    }
  }
  
  # Finally, we get the final result by joining the two dataframes generated.
  result <- cbind(intensities,result)
  # and add an additional column, called "Multiplicity"
  result[,"Multiplicity"] <- rep(1:3,nrow(dataframe))
  return(result)
}

filter_missing <- function(dataframe,columns){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   columns: vector containing the names of intensity columns
  
  # It is worth mentioning that, unlike other functions in this script, this one
  # is only suitable for our analysis, since the information regarding the 
  # experimental design is not received as an input.
  
  # We will store the row indexes that pass the filter in positions_keep 
    positions_keep <- c()
    
    # For each row (phosphosite), 
    for (i in 1:nrow(dataframe)){
      
      if(sum(c(sum(is.na(dataframe[i,columns[1:4]]))<=1), # If in at least one of the samples,
             c(sum(is.na(dataframe[i,columns[5:8]]))<=1), # the number of NAs is 1 or less
             c(sum(is.na(dataframe[i,columns[9:12]]))<=1), # (at least 3 out of 4 valid values,)
             c(sum(is.na(dataframe[i,columns[13:16]]))<=1))>0){
        # Then that phosphosite passes the filter and its index is added to positions_keep
        positions_keep <- c(positions_keep,i)
      }
    }
    
    return(dataframe[positions_keep,])
}

impute_missing_values <- function(dataframe,columns,offset=-1.8,width=0.3){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   columns: vector containing the names of intensity columns
  #   offset: determines how much imputed values will be dragged to the lower
  #           tail of the data
  #   width: determines the standard deviation of the imputed values
  
            for (j in columns){ # For each column (sample),
              
              # we fix the mean of imputed values, using the following formula
              mean_norm <- mean(dataframe[,j],na.rm=TRUE)+offset*sd(dataframe[,j],na.rm=TRUE)
              # and the corresponding standard deviation
              sd_norm <- width*sd(dataframe[,j],na.rm=TRUE)
              
              for (i in 1:nrow(dataframe)){ # For each row
                if (is.na(dataframe[i,j])){ # if a value is NA, we impute it
                  # by drawing a random value from a normal distribution
                  dataframe[i,j] <- rnorm(1,mean=mean_norm,sd=sd_norm)
                }
              }
            }
          return(dataframe)
}

get_proteome_medians <- function(dataframe,nreplicates = c(4,4,4,4),
                                 conditions=c("Col_CTL","Col_Heat","tsnKO_CTL","tsnKO_Heat")){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 
  #   nreplicates: numeric vector defining how many replicates in each group
  #   conditions: character vector of group names
  
  rows <- rownames(dataframe) # We extract row names (proteins) from the dataframe
  dataframe <- as.matrix(dataframe) # We convert the dataframe to a matrix
  nconditions <- length(conditions) # We get the number of conditions (groups)
  
  # We initialize the result dataframe with NA values
  result <- data.frame(NAs=rep(NA,nrow(dataframe)))
  
  # We will start at column 1
  firstcol <- 1
  
  for (condition in 1:length(conditions)){ # For each condition,
    # we calculate the index of the last column corresponding to that group
    lastcol <- firstcol+nreplicates[condition]-1 
    # and get the corresponding submatrix, where only the data of one group
    # is present
    submatrix <- dataframe[,firstcol:lastcol]
    # We calculate the medians of each row and add it to the result dataframe
    result <- cbind(result,rowMedians(submatrix))
    # We get the index of the first column for the next iteration
    # (first column of the next group)
    firstcol <- lastcol+1 
  }
  
  # We remove the first column (NA) of the result dataframe
  result <- result[,-1]
  
  # We name columns and rows accordingly before returning the result
  colnames(result) <- paste0("Median LFQ intensity ",conditions)
  rownames(result) <- rows
  return (result)
}

normalization_to_proteome <- function(phospho_df,protein_medians,columns,protein_column="Protein"){

  # INPUT
  #   phospho_df: MaxQuant output, imported as a dataframe object (phosphoproteomic data)
  #   protein_medians: result of the get_proteome_medians function
  #   columns: a vector containing the names of intensity columns
  #   protein_column: name of the column containing the protein IDs (in phospho_df)
  
  # First, we get the intensities values (phosphoproteomic data)
  phospho_df <- revert_log2_transformation(phospho_df,columns) 
  
  # In protein_rows, we store the names of proteins, extracted from protein_medians
  protein_rows <- rownames(protein_medians)
  
  # We get all protein IDs in a list, where each element contains all of the
  # protein IDs that are present in a row in protein_medians.
  protein_list <- list()
  for (i in protein_rows){
    protein_list <- c(protein_list,strsplit(i,";"))
  }
  
  for (iphospho in 1:nrow(phospho_df)){ # For each phosphosite,
    # we try to find whether there is proteomic_data regarding the corresponding protein
    for (iprotein in 1:nrow(protein_medians)){ 
      if (sum(phospho_df[iphospho,"Protein"]==protein_list[[iprotein]])==1){ # If there is,
        
        # We perform the normalization in the four groups
        # Normalized intensity = Phospho intensity * (Phospho intensity / Median protein intensity)
        
        phospho_df[iphospho,"Intensity.Col_CTL_1"] <- phospho_df[iphospho,"Intensity.Col_CTL_1"]*(phospho_df[iphospho,"Intensity.Col_CTL_1"]/protein_medians[iprotein,"Median LFQ intensity Col_CTL"])
        phospho_df[iphospho,"Intensity.Col_CTL_2"] <- phospho_df[iphospho,"Intensity.Col_CTL_2"]*(phospho_df[iphospho,"Intensity.Col_CTL_2"]/protein_medians[iprotein,"Median LFQ intensity Col_CTL"])
        phospho_df[iphospho,"Intensity.Col_CTL_3"] <- phospho_df[iphospho,"Intensity.Col_CTL_3"]*(phospho_df[iphospho,"Intensity.Col_CTL_3"]/protein_medians[iprotein,"Median LFQ intensity Col_CTL"])
        phospho_df[iphospho,"Intensity.Col_CTL_4"] <- phospho_df[iphospho,"Intensity.Col_CTL_4"]*(phospho_df[iphospho,"Intensity.Col_CTL_4"]/protein_medians[iprotein,"Median LFQ intensity Col_CTL"])
        
        phospho_df[iphospho,"Intensity.Col_Heat_1"] <- phospho_df[iphospho,"Intensity.Col_Heat_1"]*(phospho_df[iphospho,"Intensity.Col_Heat_1"]/protein_medians[iprotein,"Median LFQ intensity Col_Heat"])
        phospho_df[iphospho,"Intensity.Col_Heat_2"] <- phospho_df[iphospho,"Intensity.Col_Heat_2"]*(phospho_df[iphospho,"Intensity.Col_Heat_2"]/protein_medians[iprotein,"Median LFQ intensity Col_Heat"])
        phospho_df[iphospho,"Intensity.Col_Heat_3"] <- phospho_df[iphospho,"Intensity.Col_Heat_3"]*(phospho_df[iphospho,"Intensity.Col_Heat_3"]/protein_medians[iprotein,"Median LFQ intensity Col_Heat"])
        phospho_df[iphospho,"Intensity.Col_Heat_4"] <- phospho_df[iphospho,"Intensity.Col_Heat_4"]*(phospho_df[iphospho,"Intensity.Col_Heat_4"]/protein_medians[iprotein,"Median LFQ intensity Col_Heat"])
        
        phospho_df[iphospho,"Intensity.tsnKO_CTL_1"] <- phospho_df[iphospho,"Intensity.tsnKO_CTL_1"]*(phospho_df[iphospho,"Intensity.tsnKO_CTL_1"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_CTL"])
        phospho_df[iphospho,"Intensity.tsnKO_CTL_2"] <- phospho_df[iphospho,"Intensity.tsnKO_CTL_2"]*(phospho_df[iphospho,"Intensity.tsnKO_CTL_2"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_CTL"])
        phospho_df[iphospho,"Intensity.tsnKO_CTL_3"] <- phospho_df[iphospho,"Intensity.tsnKO_CTL_3"]*(phospho_df[iphospho,"Intensity.tsnKO_CTL_3"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_CTL"])
        phospho_df[iphospho,"Intensity.tsnKO_CTL_4"] <- phospho_df[iphospho,"Intensity.tsnKO_CTL_4"]*(phospho_df[iphospho,"Intensity.tsnKO_CTL_4"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_CTL"])
        
        phospho_df[iphospho,"Intensity.tsnKO_Heat_1"] <- phospho_df[iphospho,"Intensity.tsnKO_Heat_1"]*(phospho_df[iphospho,"Intensity.tsnKO_Heat_1"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_Heat"])
        phospho_df[iphospho,"Intensity.tsnKO_Heat_2"] <- phospho_df[iphospho,"Intensity.tsnKO_Heat_2"]*(phospho_df[iphospho,"Intensity.tsnKO_Heat_2"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_Heat"])
        phospho_df[iphospho,"Intensity.tsnKO_Heat_3"] <- phospho_df[iphospho,"Intensity.tsnKO_Heat_3"]*(phospho_df[iphospho,"Intensity.tsnKO_Heat_3"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_Heat"])
        phospho_df[iphospho,"Intensity.tsnKO_Heat_4"] <- phospho_df[iphospho,"Intensity.tsnKO_Heat_4"]*(phospho_df[iphospho,"Intensity.tsnKO_Heat_4"]/protein_medians[iprotein,"Median LFQ intensity tsnKO_Heat"])
        break 
        }
    }
  }
  
  # We log transform the intensities again before returning the result
  phospho_df <- log2_transformation(phospho_df,columns) 
  return(phospho_df)
  
}

add_phosphosite_name <- function(dataframe){
  
  # INPUT
  #   dataframe: MaxQuant output, imported as a dataframe object 

  # This function will add a new column, Phosphosite.name, which is a code
  # showing protein ID, modified residue position and multiplicity
  # (ex: AT5G58950.1_124_1)
  
  # For each row,
  for (i in 1:nrow(dataframe)){
    # We extract the position of the modified amino acid
    position <- strsplit(dataframe[i,"Positions.within.proteins"],"[;]")[[1]][1]
    # and the protein ID
    protein <- strsplit(dataframe[i,"Fasta.headers"]," ")[[1]][1]
    # and generate the phosphosite name
    dataframe[i,"Phosphosite.name"] <- paste0(protein,"_",position,"_",dataframe[i,"Multiplicity"])
  }
  
  # The result is the same dataframe with a new additional column
  return(dataframe)
}

go_term_present <- function(protein.ids,go.term,db.info,ancestor.info){
  
  # INPUT
  #   protein.ids: vector of protein IDs
  #   go.term: GO term ID of interest. This function will check whether
  #            this term is present in the protein IDs
  #   db.info and ancestor.info are objects from data base R packages
  
  # We initialize the result vector
  result <- c()
  
  for (i in protein.ids){ # For each protein ID,
    # We get the list of GO terms annotated to that protein
    go.list <- names(sapply(db.info[i][[1]],`[[`,1))
    for (j in go.list){ # We add the ancestors of the terms
      go.list <- c(go.list,ancestor.info[j][[1]])
    }
    # before checking whether the GO term is present in the annotation of that protein
    result[i] <- (go.term %in% go.list)
  }
  
  # The result is a named logical vector, where names are proteins and TRUE/FALSE
  # indicate whether the GO term is present in the annotation of the protein
  return(result)
}

get_protein_ids <- function(phosphosites){
  
  # INPUT
  #   phosphosites: character vector of phosphosite names 
  
  # This function simply extracts the protein IDs from phosphosite names,
  # that is, the code before the first underscore (ex: AT1G35580.2)
  return(unique(sapply(strsplit(phosphosites,"_"), `[[`, 1)))
  
}

get_proteins <- function(phosphosites){
  # INPUT
  #   phosphosites: character vector of phosphosite names 
  
  # This function simply extracts the proteins from phosphosite names,
  # that is, the code before the first dot (ex: AT1G35580)
  return(unique(substr(phosphosites,1,9)))
  
}

export_phosphosites_list <- function(vpdata,path){
  
  # INPUT
  #   vpdata: dataframe containing volcano plot data generated in the 
  #   phospho_diff_analysis.
  #   path: output path where txt files will be saved
  
  # We extract the phosphosite names of features that are upregulated or downregulated
  upregulated <- vpdata[vpdata$significant=="up",]
  downregulated <- vpdata[vpdata$significant=="down",]
  
  # We extract the information we will export for upregulated phosphosites:
  upregulated$`-log.adj.p` <- -log10(upregulated$adj.p.val) # -log10 (adjusted p value)
  upregulated <- upregulated[,c(1,4)] # We only keep the first four columns of the dataframe
  upregulated$`Protein ID` <- sapply(strsplit(rownames(upregulated),"_"), `[[`, 1) # protein IDs
  upregulated$Position <- sapply(strsplit(rownames(upregulated),"_"), `[[`, 2) # modified residue position
  upregulated$Multiplicity <- sapply(strsplit(rownames(upregulated),"_"), `[[`, 3) # multiplicity
  upregulated$`Fasta headers` <- phospho_imputed[rownames(upregulated),"Fasta.headers"] # Fasta headers, from which
  upregulated$`Fasta headers` <- sapply(strsplit(upregulated$`Fasta headers`,"\\|"), `[[`, 2) # we remove the protein IDs
  
  # The following for loop will add gene symbols to the Fasta headers, to allow
  # for easier identification of interesting proteins
  for (i in 1:nrow(upregulated)){
    symbols <- paste(getsymbols[substr(upregulated[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      upregulated[i,"Fasta headers"] <- paste0(upregulated[i,"Fasta headers"],"(",symbols,")")
    }
  }
  
  # We perform the same process for downregulated phosphosites
  downregulated$`-log.adj.p` <- -log10(downregulated$adj.p.val)
  downregulated <- downregulated[,c(1,4)]
  downregulated$`Protein ID` <- sapply(strsplit(rownames(downregulated),"_"), `[[`, 1)
  downregulated$Position <- sapply(strsplit(rownames(downregulated),"_"), `[[`, 2)
  downregulated$Multiplicity <- sapply(strsplit(rownames(downregulated),"_"), `[[`, 3)
  downregulated$`Fasta headers` <- phospho_imputed[rownames(downregulated),"Fasta.headers"]
  downregulated$`Fasta headers` <- sapply(strsplit(downregulated$`Fasta headers`,"\\|"), `[[`, 2)
  
  for (i in 1:nrow(downregulated)){
    symbols <- paste(getsymbols[substr(downregulated[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      downregulated[i,"Fasta headers"] <- paste0(downregulated[i,"Fasta headers"],"(",symbols,")")
    }
  }
  
  # and export both tables as txt files
  write.table(upregulated,file=paste0(path,"_up.txt"),row.names=F,quote=F, sep="\t", dec=",")
  write.table(downregulated,file=paste0(path,"_down.txt"),row.names=F,quote=F, sep="\t", dec=",")
  
}

export_proteins_list <- function(vpdata,path){
  
  # INPUT
  #   vpdata: dataframe containing volcano plot data generated in the 
  #   protein_diff_analysis.
  #   path: output path where txt files will be saved
  
  # This function is very similar to export_phosphosites_list, but less information
  # needs to be extracted because the exported data is simply a list of proteins.
  
  upregulated <- vpdata[vpdata$significant=="up",]
  downregulated <- vpdata[vpdata$significant=="down",]
  
  upregulated$`-log.adj.p` <- -log10(upregulated$adj.p.val)
  upregulated <- upregulated[,c(1,4)]
  upregulated$`Protein ID` <- rownames(upregulated)
  upregulated$`Fasta headers` <- protein_log[rownames(upregulated),"Fasta.headers"]
  upregulated$`Fasta headers` <- sapply(strsplit(upregulated$`Fasta headers`,"\\|"), `[[`, 2)
  
  for (i in 1:nrow(upregulated)){
    symbols <- paste(getsymbols[substr(upregulated[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      upregulated[i,"Fasta headers"] <- paste0(upregulated[i,"Fasta headers"],"(",symbols,")")
    }
  }
  
  downregulated$`-log.adj.p` <- -log10(downregulated$adj.p.val)
  downregulated <- downregulated[,c(1,4)]
  downregulated$`Protein ID` <- rownames(downregulated)
  downregulated$`Fasta headers` <- protein_log[rownames(downregulated),"Fasta.headers"]
  downregulated$`Fasta headers` <- sapply(strsplit(downregulated$`Fasta headers`,"\\|"), `[[`, 2)
  
  for (i in 1:nrow(downregulated)){
    symbols <- paste(getsymbols[substr(downregulated[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      downregulated[i,"Fasta headers"] <- paste0(downregulated[i,"Fasta headers"],"(",symbols,")")
    }
  }
  
  upregulated$`Protein ID` <- substr(upregulated$`Protein ID`,1,9)
  downregulated$`Protein ID` <- substr(downregulated$`Protein ID`,1,9)
  
  write.table(upregulated,file=paste0(path,"_up.txt"),row.names=F,quote=F, sep="\t", dec=",")
  write.table(downregulated,file=paste0(path,"_down.txt"),row.names=F,quote=F, sep="\t", dec=",")
  
}

phosphosites_in_idr <- function(pps_list){
  
  # INPUT
  #   pps_list: character vector of phosphosite names to be analyzed
  
  # Note: this function uses objects that are generated in the 
  # phospho_downstream_analysis.R script, even though they are not 
  # passed to the function as arguments
  
  # We intialize the result vector
  result <- c()
  
  for (i in 1:length(pps_list)){ # For each phosphosite,
    
    protein.id <- strsplit(pps_list[i],"_")[[1]][1] # we extract the protein ID
    position <- as.numeric(strsplit(pps_list[i],"_")[[1]][2]) # and the modified residue position
    
    # Using the IUPRED results (iupred_df), we check whether that position is in
    # an IDR or not
      
      # Positions with scores below 0.5 are considered ordered
      if (iupred_df[protein.id,position]<=0.5){result[pps_list[i]] <- "Ordered"} 
      else {
        
        # If it is a disordered residue, we additionally check for other residues within the protein
        result[pps_list[i]] <- "Disordered"
        
        # We get a list of the other phosphorylated positions that are found in the same protein
        closest.positions <- sort(as.numeric(sapply(strsplit(pps_list[sapply(strsplit(pps_list,"_"), `[[`, 1)==protein.id],"_"), `[[`, 2)))
        closest.positions <- closest.positions[c(which(closest.positions==position)-1,which(closest.positions==position)+1)]
        closest.positions <- closest.positions[!is.na(closest.positions)]
        
        if (length(closest.positions)>1) {
          # If there was more than one residue in the protein, we annotate accordingly
        result[pps_list[i]] <- "Shares protein with another phosphosite"
        
        # Finally, we check if those positions share an IDR with the current residue
        for (j.position in closest.positions){
          if (position != j.position){
            if(sum(iupred_df[protein.id,position:j.position]>0.5)/(abs(position-j.position+1))>0.9){
              result[pps_list[i]] <- "Shares IDR with another phosphosite"
            }
          }
        }
      }
    }
  }
  return(result)
}

export_genes_list <- function(vpdata,path){
  
  # INPUT
  #   vpdata: dataframe containing volcano plot data generated in the 
  #   transc_diff_analysis.
  #   path: output path where txt files will be saved
  
  # This function is very similar to export_phosphosites_list, but is modified
  # appropriately for exporting transcriptomic results
  
  upregulated <- vpdata[vpdata$significant=="up",]
  downregulated <- vpdata[vpdata$significant=="down",]
  
  upregulated$`-log.adj.p` <- -log10(upregulated$adj.p.val)
  upregulated <- upregulated[,c(1,4)]
  upregulated$`Protein ID` <- rownames(upregulated)
  upregulated$`Description` <- transc_data[upregulated$`Protein ID`,"Description"]
  
  for (i in 1:nrow(upregulated)){
    symbols <- paste(getsymbols[substr(upregulated[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      upregulated[i,"Description"] <- paste0(upregulated[i,"Description"],"(",symbols,")")
    }
  }
  
  downregulated$`-log.adj.p` <- -log10(downregulated$adj.p.val)
  downregulated <- downregulated[,c(1,4)]
  downregulated$`Protein ID` <- rownames(downregulated)
  downregulated$`Description` <- transc_data[downregulated$`Protein ID`,"Description"]
  
  for (i in 1:nrow(downregulated)){
    symbols <- paste(getsymbols[substr(downregulated[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      downregulated[i,"Description"] <- paste0(downregulated[i,"Description"]," (",symbols,")")
    }
  }
  
  upregulated$`Protein ID` <- substr(upregulated$`Protein ID`,1,9)
  downregulated$`Protein ID` <- substr(downregulated$`Protein ID`,1,9)
  
  write.table(upregulated,file=paste0(path,"_up.txt"),row.names=F,quote=F, sep="\t", dec=",")
  write.table(downregulated,file=paste0(path,"_down.txt"),row.names=F,quote=F, sep="\t", dec=",")
  
}

export_cluster_list <- function(phosphosites,path){
  
  # INPUT
  #   phosphosites: character vector containing phosphosite names
  #   path: output path where the txt file will be saved
  
  # This function is very similar to export_phosphosites_list, but only contains
  # the information regarding the phosphosites in each cluster, such as the 
  # protein ID, modified residue position, multiplicity, fasta headers and
  # gene symbols (no intensity values)
  
  result <- data.frame(phosphosites)
  
  result$`Protein ID` <- sapply(strsplit(phosphosites,"_"), `[[`, 1)
  result$Position <- sapply(strsplit(phosphosites,"_"), `[[`, 2)
  result$Multiplicity <- sapply(strsplit(phosphosites,"_"), `[[`, 3)
  result$`Fasta headers` <- phospho_imputed[phosphosites,"Fasta.headers"]
  result$`Fasta headers` <- sapply(strsplit(result$`Fasta headers`,"\\|"), `[[`, 2)
  
  result <- result[,-1]
  
  for (i in 1:nrow(result)){
    symbols <- paste(getsymbols[substr(result[i,"Protein ID"],1,9)][[1]],collapse=", ")
    if (symbols!=""){
      result[i,"Fasta headers"] <- paste0(result[i,"Fasta headers"],"(",symbols,")")
    }
  }
  
  write.table(result,file=path,row.names=F,quote=F, sep="\t", dec=",")
  
}
