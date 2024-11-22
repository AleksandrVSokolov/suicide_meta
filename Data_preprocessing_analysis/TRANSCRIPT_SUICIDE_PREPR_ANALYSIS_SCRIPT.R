setwd("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis")

# Setting options
getOption("scipen") # Default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)
options(show.error.messages = TRUE)

################### Package import ###################
library(R.utils)
library(fun)
library(stringr)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(grid)
library(gdata)
library(RColorBrewer)
library(networkD3)
library(webshot)
library(htmlwidgets)
library(magrittr)
library(igraph)
library(visNetwork)
library(data.table)
library(XML)
library(rvest)
library(RCurl)
library(HGNChelper)
library(stringi)
library(httr)
library(lubridate)
library(rjson)
library(rtracklayer)
library(rstudioapi)
library(tidyr)
library(Gviz)
library(limma)
library(FactoMineR)
library(ggthemes)
library(igraph)
library(RSelenium)
library(lumi)
library(outliers)
library(svglite)
library(scatterplot3d)
library(sva)
library(jsonlite)
library(ggrepel)
library(parallel)
library(bacon)
library(gridExtra)
library(ggplotify)
library(HGNChelper)
library(jetset)
library(GEOquery)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(chromoMap)
library(RIdeogram)
library(ggVennDiagram)
library(seqinr)
library(Biostrings)
library(dbparser)
library(LDlinkR)
library(parallel)
library(HardyWeinberg)
library(sva)
library(MatrixEQTL)
library(MASS)
library(oligo)
library(pd.huex.1.0.st.v2)
library(biomaRt)
library(affydata)
library(SRAdb)
library(rentrez)
library(org.Hs.eg.db)
library(edgeR)
library(Matrix)
library(Seurat)

################### Defining functions ###################

# NOT IN operator
'%!in%' = function(x,y){!('%in%'(x,y))}

# A function to convert list to data frame
list_to_df = function(data_list){
  if (length(data_list) > 1){
    data_list = do.call(rbind, data_list)
  } else {
    data_list = data_list[[1]]
  }
  return(data_list)
}

# A function to replace multiple patterns by multiple replacements in a string
multiple_stri_replacer = function(string, pattern_vector, replacement_vector){
  
  # Pattern_vector and replacement_vector should have the same length
  for (i in 1:length(pattern_vector)){
    string = stri_replace_all_fixed(str = string, pattern = pattern_vector[i], replacement = replacement_vector[i])
  }
  return(string)
}

# A function to read text files fast; uses data.table::fread
smart_fread = function(x, ...){
  x = as.data.frame(fread(x, nThread = 10, ...))
  if ("V1" %in% colnames(x)){
    rownames(x) = x$V1
    x$V1 = NULL
  }
  return(x)
}

# A function to detect at least one pattern in a string
multiple_stri_detector = function(string, pattern_vector){
  output_list = list()
  for (i in 1:length(pattern_vector)){
    output_list[[i]] = stri_detect_fixed(str = string, pattern = pattern_vector[i])
  }
  output_list = do.call(rbind, output_list)
  apply(output_list, 2, any)
}

# A function to expand a data frame where several columns contain condensed cells with a specified separator
multiple_expander = function(df, cols_to_expand, pattern){
  #
  orig_colnames = colnames(df)
  df_modif = df[, cols_to_expand, drop = FALSE]
  df_const = df[,-cols_to_expand, drop = FALSE]
  orig_colnames_modif = colnames(df_modif)
  
  # Running expansion
  df_list = list()
  for (i in 1:nrow(df_const)){
    print(i)
    curr_df_const = df_const[i,, drop = FALSE]
    curr_df_modif = df_modif[i,, drop = FALSE]
    curr_df_modif = apply(curr_df_modif, 2, function(x) unlist(stri_split_fixed(x, pattern = pattern)), simplify = FALSE)
    
    if (length(cols_to_expand) > 1){
      curr_df_modif = do.call(cbind, curr_df_modif)
    } else {
      curr_df_modif = unlist(curr_df_modif)
    }
    
    if (is.matrix(curr_df_modif)){
      for (b in 1:nrow(curr_df_modif)){
        curr_df_const[b, ] = curr_df_const[1, ]
      }
    } else {
      for (b in 1:length(curr_df_modif)){
        curr_df_const[b, ] = curr_df_const[1, ]
      }
    }
    
    curr_df_combined = cbind(curr_df_const, curr_df_modif)
    
    if (length(cols_to_expand) == 1){
      colnames(curr_df_combined)[ncol(curr_df_combined)] = orig_colnames[cols_to_expand]
    }
    df_list[[i]] = curr_df_combined
  }
  
  if (length(df_list) > 1){
    df_list = do.call(rbind, df_list)
  } else {
    df_list = df_list[[1]]
  }
  
  df_list = df_list[, orig_colnames]
  return(df_list)
}

fix_columns = function(dataset){
  
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = ":ch1", "")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = " ", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = ".", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = ",", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = ":", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = "(", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = ")", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = "-", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = "=", "_")
  colnames(dataset) = stri_replace_all_fixed(colnames(dataset), pattern = ":", "_")
  colnames(dataset) = toupper(colnames(dataset))
  return(dataset)
  
}



# A function to check gene symbols
# Importing NIH dataset
Homo_Sapiens_Gene_info_NIH = smart_fread("/home/aleksandr/Desktop/WORK/UCSC_ID_MAP/Homo_sapiens.gene_info")
# https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/ (replace with an appropriate path)
Homo_Sapiens_Gene_info_NIH_expanded = multiple_expander(df = Homo_Sapiens_Gene_info_NIH, cols_to_expand = 5, pattern = "|")
check_gene_symbol_NIH = function(PRF_gene_symbols, PRF_ref_NIH_expanded, PRF_replace_NA_with_old = FALSE){
  PRF_gene_symbols_check = lapply(PRF_gene_symbols, function(x){
    if (x %in% PRF_ref_NIH_expanded$Symbol_from_nomenclature_authority){
      Curr_gene = x
      Approved = TRUE
      Suggested.Symbol = x
    } else if (x %in% PRF_ref_NIH_expanded$Symbol){
      RRF_df = PRF_ref_NIH_expanded[PRF_ref_NIH_expanded$Symbol == x,]
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = RRF_df[,"Symbol_from_nomenclature_authority"]
      Suggested.Symbol = unique(Suggested.Symbol)
      Suggested.Symbol = Suggested.Symbol[Suggested.Symbol != "-"]
      if (length(Suggested.Symbol) >= 1){
        Suggested.Symbol = Suggested.Symbol[1]
        if (Suggested.Symbol == x){
          Approved = TRUE
        }
      } else {
        Suggested.Symbol = RRF_df[,"Symbol"]
        Suggested.Symbol = unique(Suggested.Symbol)
        Suggested.Symbol = Suggested.Symbol[1]
      }
    } else if (x %in% PRF_ref_NIH_expanded$Synonyms){
      RRF_df = PRF_ref_NIH_expanded[PRF_ref_NIH_expanded$Synonyms == x,]
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = RRF_df[,"Symbol_from_nomenclature_authority"]
      Suggested.Symbol = unique(Suggested.Symbol)
      Suggested.Symbol = Suggested.Symbol[Suggested.Symbol != "-"]
      if (length(Suggested.Symbol) >= 1){
        Suggested.Symbol = Suggested.Symbol[1]
      } else {
        Suggested.Symbol = RRF_df[,"Symbol"]
        Suggested.Symbol = unique(Suggested.Symbol)
        Suggested.Symbol = Suggested.Symbol[1]
      }
    } else {
      Curr_gene = x
      Approved = FALSE
      Suggested.Symbol = NA
      if (PRF_replace_NA_with_old){
        Suggested.Symbol = x
      }
    }
    Dataset = data.frame(x = Curr_gene, Approved = Approved, Suggested.Symbol = Suggested.Symbol)
  })
  PRF_gene_symbols_check = list_to_df(PRF_gene_symbols_check)
  return(PRF_gene_symbols_check)
}


map2color=function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

# A function to process a variable from GSE dataset (in the form of "Sex: Male" or "Age: 20", etc)
preprocess_GSE_var = function(x){
  x = sapply(x, function(z){
    z = stri_split_fixed(z, pattern = ":")
    z = unlist(z)
    z = z[2]
    z = str_trim(z)
    return(z)
  })
  return(x)
}

# A function to process GSE variable names
preprocess_GSE_var_names = function(x){
  x = sapply(x, function(z){
    z = stri_split_fixed(z, pattern = ":")
    z = unlist(z)
    z = z[1]
    z = str_trim(z)
    return(z)
  })
  return(x)
}

# A function to get value table from GEO series file (sometimes it is different from the standard import)
geo_data_table_gsm_extract = function(GEO_ID){
  Data_GEO = getGEOfile(
    GEO = GEO_ID,
    destdir = getwd(), amount = "full"
  )
  gunzip(Data_GEO)
  Data_GEO = stri_replace_all_fixed(Data_GEO, pattern = ".gz", replacement = "") #useful!
  Data_full = parseGEO(Data_GEO)
  file.remove(Data_GEO)
  GSMS = Data_full@gsms
  Full_GSMS = lapply(GSMS, function(x) x@dataTable@table)
  Full_Columns = lapply(GSMS, function(x) x@dataTable@columns)
  Full_Columns = Full_Columns[[1]]
  Full_GSMS = do.call(cbind, Full_GSMS)
  Output = list()
  Output[[1]] = Full_Columns
  Output[[2]] = Full_GSMS
  return(Output)
}

# A function to get phenotype table from GEO series file (sometimes it is different from the standard import)
geo_data_table_gsm_extract_phenotypes = function(GEO_ID){
  Data_GEO = getGEOfile(
    GEO = GEO_ID,
    destdir = getwd(), amount = "full"
  )
  gunzip(Data_GEO, remove = TRUE, overwrite = TRUE)
  Data_GEO = stri_replace_all_fixed(Data_GEO, pattern = ".gz", replacement = "")
  Data_full = parseGEO(Data_GEO)
  GSMS = Data_full@gsms
  file.remove(Data_GEO)
  Full_GSMS = lapply(GSMS, function(x) x@header)
  Colnames_GSMS = names(Full_GSMS[[1]])
  Full_GSMS_DFs = lapply(Full_GSMS, function(x){
    x = lapply(x, function(j){
      if (length(j) > 1){
        j = paste0(j, collapse = "___")
      }
      return(j)
    })
    x = data.frame(x)
    char_vector = x$characteristics_ch1
    char_vector = unlist(stri_split_fixed(char_vector, pattern =  "___"))
    vars_GSM = preprocess_GSE_var(char_vector)
    names_GSM = preprocess_GSE_var_names(char_vector)
    vars_GSM = as.data.frame(vars_GSM)
    vars_GSM = t(vars_GSM)
    vars_GSM = as.data.frame(vars_GSM)
    colnames(vars_GSM) = names_GSM
    x$characteristics_ch1 = NULL
    x$data_row_count = NULL
    x = x[1,]
    x = cbind(x, vars_GSM)
    rownames(x) = NULL
    return(x)
  })
  Data_pheno_GSMS = as.data.frame(do.call(rbind, Full_GSMS_DFs))
  return(Data_pheno_GSMS)
}


generate_palette = function(variable, random = TRUE) {
  # Determine the number of unique categories
  num_categories = length(unique(variable))
  
  print(num_categories)
  # Choose a palette type based on the number of categories
  if (num_categories > 9) {
    # For more than 9 categories, use the Set3 palette but recycle colors if necessary
    colors = colorRampPalette(brewer.pal(9, "Spectral"))(num_categories)
  } else {
    # Use distinct colors from Set1, Set2, or Set3 as appropriate
    if (num_categories <= 3) {
      colors = brewer.pal(num_categories, "Dark2")
    } else if (num_categories <= 8) {
      colors = brewer.pal(num_categories, "Dark2")
    } else {
      colors = brewer.pal(num_categories, "Dark2")
    }
  }
  
  # Return the color palette
  if (random){
    colors = colors[sample(1:length(colors), size = length(colors))]
    print(colors)
  }
  names(colors) = unique(variable)
  
  return(colors)
}

characterize_categorical_variable = function(df, Category_1, Category_2, Variable_name, P.val.valid = TRUE, keep_missing = TRUE){
  Total_variable = df[, Variable_name]
  
  if (!is.ordered(Variable_name)){
    df[, Variable_name] = ordered(Total_variable, levels = sort(unique(Total_variable)))
    Total_variable = df[, Variable_name]
  }
  
  Contingency_table = table(df[,Variable_name], df[,"variable_to_split"])
  
  if (keep_missing){
    Contingency_table = table(df[,Variable_name], df[,"variable_to_split"], exclude = NULL)
  }
  
  Contingency_table_percents_1 = (Contingency_table[,Category_1]/(sum(Contingency_table[,Category_1])) * 100) %>% round(., digits = 1)
  Contingency_table_percents_2 = (Contingency_table[,Category_2]/(sum(Contingency_table[,Category_2])) * 100) %>% round(., digits = 1)
  Report_1 = paste0(row.names(Contingency_table), ": ")
  Report_1 = paste0(Report_1, Contingency_table[,Category_1], " (",Contingency_table_percents_1, "%)", collapse = "\n")
  Report_2 = paste0(row.names(Contingency_table), ": ")
  Report_2 = paste0(Report_2, Contingency_table[,Category_2], " (",Contingency_table_percents_2, "%)", collapse = "\n")
  
  Contingency_table_pval = table(df[,Variable_name], df[,"variable_to_split"])
  if (P.val.valid){
    tryCatch({
      Var.pval <<- chisq.test(Contingency_table_pval)
      Var.pval <<- Var.pval$p.value}, warning = function(w){
        message = paste0(Variable_name, " has small counts in some groups, using P.val from Fisher's Exact Test")
        writeLines(message)
        Var.pval <<- fisher.test(Contingency_table_pval)
        Var.pval <<- Var.pval$p.value
      })
  } else {
    Var.pval = "Not valid"
  }
  
  Output_list = list()
  Output_list[[1]] = Report_1
  Output_list[[2]] = Report_2
  Output_list[[3]] = Var.pval
  return(Output_list)
}



characterize_dataset_generelized_two_subgroups = function(dataset,
                                                          study_char,
                                                          contrast_col_number,
                                                          contrast_vector,
                                                          participants_col_number,
                                                          model_covariates_col_vector,
                                                          columns_to_characterise_vector,
                                                          Remove_NA_predictors = FALSE,
                                                          drop_P = TRUE,
                                                          simplif_P = 3){
  
  # Running Raw Statistics
  total_number = nrow(dataset)
  total_number_message = paste0("Initial dataset includes ", total_number, " participants")
  writeLines(total_number_message)
  first_row_output = data.frame(Name = total_number_message,Category.1 = "", Category.2 = "", P.value = NA)
  
  # Removing missing data
  colnames(dataset)[participants_col_number] = "Participant_ID_FUN"
  indeces = c(participants_col_number, contrast_col_number, model_covariates_col_vector)
  indeces = unique(indeces)
  Model_df = dataset[,indeces]
  if (Remove_NA_predictors){
    
    writeLines("Removing participants with missing data")
    Participants_missing_data = apply(Model_df, 1, function(x){
      x = as.character(x)
      
      if (any(is.na(x))){
        return(TRUE)
      } else{
        return(FALSE)
      }
      
    })
    Participants_missing_data = Model_df$Participant_ID_FUN[Participants_missing_data]
    
    if (length(Participants_missing_data) < 1){
      excluded_number = 0
      
    } else {
      
      excluded_number = length(Participants_missing_data)
    }
    
    Participants_excluded_df = dataset[dataset$Participant_ID_FUN %in% Participants_missing_data, ]
    dataset = dataset[dataset$Participant_ID_FUN %!in% Participants_missing_data, ]
    excluded_message = paste0("Participants with missing data excluded: ", excluded_number,"\n",
                              "Resulting number of participants: ", nrow(dataset))
    second_row_output = data.frame(Name = excluded_message,Category.1 = "", Category.2 = "", P.value = NA)
    writeLines(excluded_message)
    
  } else {
    
    writeLines("Exclusion of participants was not performed")
    second_row_output = data.frame(Name = "Exclusion of participants was not performed", Category.1 = "", Category.2 = "", P.value = NA)
    Participants_excluded_df = NA
  }
  
  main_indeces = c(contrast_col_number, columns_to_characterise_vector)
  main_indeces = unique(main_indeces)
  Main_dataset = list()
  
  
  # Making categories
  Initial_names = colnames(dataset)[main_indeces]
  colnames(dataset)[contrast_col_number] = "variable_to_split"
  levels_var = dataset[,"variable_to_split"]
  levels_var = levels(levels_var)
  names(levels_var) = NULL
  Category_1 = levels_var[1]
  Category_2 = levels_var[2]
  
  for (i in 1:length(main_indeces)){
    curr_variable = dataset[, main_indeces[i]]
    if (is.factor(curr_variable)){
      Characterised_var = characterize_categorical_variable(df = dataset,
                                                            Category_1 = Category_1, 
                                                            Category_2 = Category_2, 
                                                            Variable_name = colnames(dataset)[main_indeces[i]], 
                                                            keep_missing = TRUE)
    } else if (is.character(curr_variable)){
      stop("All columns should be either a Factor or Numeric")
    } else if (is.numeric(curr_variable)){
      Characterised_var = characterize_numeric_variable(df = dataset, 
                                                        Category_1 = Category_1, 
                                                        Category_2 = Category_2, 
                                                        Variable_name = colnames(dataset)[main_indeces[i]], 
                                                        Mention_NAs = TRUE)
    }
    Current_DF = data.frame(Name = Initial_names[i], 
                            Category.1 = Characterised_var[[1]], 
                            Category.2 = Characterised_var[[2]], 
                            P.value = Characterised_var[[3]])
    Main_dataset[[i]] = Current_DF
  }
  Main_dataset = list_to_df(data_list = Main_dataset)
  Main_dataset$P.value = round(Main_dataset$P.value, digits = simplif_P)
  
  # Compiling full dataframe
  header = data.frame(Name = study_char, Category.1 = "", Category.2 = "", P.value = NA)
  Full_table = rbind(header, first_row_output, second_row_output, Main_dataset)
  
  if(drop_P){
    Full_table$P.value = NULL
  }
  
  output = list()
  output$Table = Full_table
  output$Excluded = Participants_excluded_df
  return(output)
}


characterize_numeric_variable = function(df, Category_1, Category_2, Variable_name, Mention_NAs = TRUE){
  Total_variable = df[, Variable_name]
  Variable_1 = df[df$variable_to_split == Category_1, Variable_name]
  Variable_1_report = describe_vector_numeric(Variable_1, Mention_NAs = Mention_NAs)
  Variable_2 = df[df$variable_to_split == Category_2, Variable_name]
  Variable_2_report = describe_vector_numeric(Variable_2, Mention_NAs = Mention_NAs)
  Shapiro_normality_check = shapiro.test(Total_variable)
  Shapiro_normality_check = Shapiro_normality_check$p.value
  if (Shapiro_normality_check < 0.05){
    message = paste0(Variable_name, " is not normally distributed. Using P-val from Mann Whitney U Test")
    writeLines(message)
    formula_test = paste0(Variable_name, "~ variable_to_split")
    Total_pval = wilcox.test(formula = as.formula(formula_test), data = df)
    Total_pval = Total_pval$p.value
  } else {
    message = paste0(Variable_name, " is normally distributed. Using P-val from T Test")
    writeLines(message)
    formula_test = paste0(Variable_name, "~ variable_to_split")
    Var_test_check = var.test(formula	= as.formula(formula_test), data = df)
    Var_test_check = Var_test_check$p.value 
    if (Var_test_check >= 0.05){
      Total_pval = t.test(formula = as.formula(formula_test), data = df, var.equal = TRUE)
    } else  {
      Total_pval = t.test(formula = as.formula(formula_test), data = df, var.equal = FALSE)
    }
    Total_pval = Total_pval$p.value
  }
  Output_list = list()
  Output_list[[1]] = Variable_1_report
  Output_list[[2]] = Variable_2_report
  Output_list[[3]] = Total_pval
  return(Output_list)
}



################### Array-based cohorts ###################


################### GSE208338 ###################
geo_record_GSE208338 = getGEO("GSE208338")
geo_record_GSE208338 = geo_record_GSE208338[[1]]
geo_pheno_GSE208338 = pData(geo_record_GSE208338)
geo_pheno_GSE208338 = fix_columns(geo_pheno_GSE208338)

# Variables used in the paper https://www.medrxiv.org/content/10.1101/2022.08.09.22278128v2.full.pdf
"
For the differential expression analyses, we collapsed all cases (schizophrenia, BD, MDD) and limma version
3.42.2 package in R 3.6.1 was used. We calculated the effect of case-control status on the difference in gene
expression controlling for Age, Sex, pH, PMI, RIN, RIN2
, Suicide and Type of Death (Natural, Violent and
Non-Violent) as well as the first four dimensions of the genotype defined ancestry (Dim1-Dim4) and one
found surrogate variable (SV1), which was not correlated to any other covariate of the model (Figure S2C,E).

"
geo_pheno_GSE208338$AGE
geo_pheno_GSE208338$SEX
geo_pheno_GSE208338$PH
geo_pheno_GSE208338$PMI
geo_pheno_GSE208338$RIN
geo_pheno_GSE208338$C1
geo_pheno_GSE208338$C2
geo_pheno_GSE208338$C3
geo_pheno_GSE208338$C4
geo_pheno_GSE208338$DIAGNOSIS
geo_pheno_GSE208338$SV1
geo_pheno_GSE208338$SUICIDE
geo_pheno_GSE208338$TYPE_OF_DEATH

"
SVA version 3.34.0 93,94 package in R 3.6.1 was used for batch correction of known and hidden batches. Batch
correction was conducted at probeset level before annotation and summarization to different genetic levels.
We first removed the five known batches with SVA’s ComBat function and then applied surrogate variable
(SV) analysis to remove hidden effects. One SV was found, which was the only one explaining high variance
in the expression dataset 
"

# This batch covariate is not shown in the phenotype -> we need to use preprocessed data
"Filtered and batch corrected RMA signal estimates from oligo 1.50.0 package in R 3.6.1"
expression_GSE208338 = exprs(geo_record_GSE208338)


# phenotype curation
GSE208338_pheno_curated = geo_pheno_GSE208338

GSE208338_pheno_curated$AGE = as.numeric(GSE208338_pheno_curated$AGE)
GSE208338_pheno_curated$SEX = ifelse(GSE208338_pheno_curated$SEX == "M", "Male", "Female")
GSE208338_pheno_curated$SEX = factor(GSE208338_pheno_curated$SEX, levels = c("Female", "Male"))

GSE208338_pheno_curated$PH = as.numeric(GSE208338_pheno_curated$PH)
GSE208338_pheno_curated$PMI = as.numeric(GSE208338_pheno_curated$PMI)
GSE208338_pheno_curated$RIN = as.numeric(GSE208338_pheno_curated$RIN)
GSE208338_pheno_curated$C1 = as.numeric(GSE208338_pheno_curated$C1)
GSE208338_pheno_curated$C2 = as.numeric(GSE208338_pheno_curated$C2)
GSE208338_pheno_curated$C3 = as.numeric(GSE208338_pheno_curated$C3)
GSE208338_pheno_curated$C4 = as.numeric(GSE208338_pheno_curated$C4)
GSE208338_pheno_curated$DIAGNOSIS = as.factor(GSE208338_pheno_curated$DIAGNOSIS)
GSE208338_pheno_curated$SV1 = as.numeric(GSE208338_pheno_curated$SV1)
GSE208338_pheno_curated$SUICIDE = ifelse(GSE208338_pheno_curated$SUICIDE == "Y", "Suicide", "Control")
GSE208338_pheno_curated$SUICIDE = factor(GSE208338_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
GSE208338_pheno_curated$TYPE_OF_DEATH = factor(GSE208338_pheno_curated$TYPE_OF_DEATH, levels = c("Natural", "Non-violent", "Violent"))
table(GSE208338_pheno_curated$SUICIDE, GSE208338_pheno_curated$TYPE_OF_DEATH)

"          Natural Non-violent Violent
  Control      95           3      12
  Suicide       0          22      37
"

'
 1) "natural"
for all natural causes, e.g. ischemic heart disease or pneumonia, 2) "violent" for violent suicides or accidents,
such as hanging, drowning or a car accident, and 3) "non-violent" for all types of poisoning.
'
# TYPE_OF_DEATH is included as Non-violent Violent separates both suicide individuals and controls

### check matching between expression and phenotypes
all(colnames(expression_GSE208338) == GSE208338_pheno_curated$GEO_ACCESSION) # TRUE
GSE208338_probes = featureData(geo_record_GSE208338)
GSE208338_probes = GSE208338_probes@data

probes_huex_1 = smart_fread("GSE248260/GPL5188-122.txt")

#  Use BioMart to get the gene symbols (performed in python) affy_huex_1_0_st_v2
# http://grch37.ensembl.org/biomart

# mapping probe IDs to gene symbols
mapped_ensembl_huex_1 = as.data.frame(fread("mapping_huex_enseml.txt", nThread = 10, sep ="\t"))
GSE208338_probes_ids = GSE208338_probes$ID

GSE208338_probes_gene_symbols = mclapply(GSE208338_probes_ids, function(x){
  
  if (x %!in% mapped_ensembl_huex_1$V4){
    return(c(NA, NA))
  }
  
  gene_names = mapped_ensembl_huex_1[mapped_ensembl_huex_1$V4 == x, "V5"]
  gene_names = str_trim(gene_names)
  gene_names = unique(gene_names)
  gene_names = gene_names[gene_names!=""]
  gene_names = gene_names[!is.na(gene_names)]
  gene_names = paste0(gene_names, collapse = ";")
  
  gene_names_non_HGNC = mapped_ensembl_huex_1[mapped_ensembl_huex_1$V4 == x, "V2"]
  gene_names_non_HGNC = str_trim(gene_names_non_HGNC)
  gene_names_non_HGNC = unique(gene_names_non_HGNC)
  gene_names_non_HGNC = gene_names_non_HGNC[gene_names_non_HGNC!=""]
  gene_names_non_HGNC = gene_names_non_HGNC[!is.na(gene_names_non_HGNC)]
  gene_names_non_HGNC = paste0(gene_names_non_HGNC, collapse = ";")
  
  return(c(gene_names, gene_names_non_HGNC))
}, mc.cores = 6)


GSE208338_probes$Gene_symbol = sapply(GSE208338_probes_gene_symbols, function(x){
  return(x[1])
})

GSE208338_probes$Gene_symbol_non_hgnc = sapply(GSE208338_probes_gene_symbols, function(x){
  return(x[2])
})

# Row match check
all(rownames(expression_GSE208338) == GSE208338_probes$ID) # TRUE
rownames(GSE208338_probes) = GSE208338_probes$ID
all(colnames(expression_GSE208338) == GSE208338_pheno_curated$GEO_ACCESSION) # TRUE

# Differential expression
Design.matrix = model.matrix(~ SUICIDE + AGE + SEX + PH + PMI + RIN + C1 + C2 + C3 + C4 + DIAGNOSIS + SV1 + TYPE_OF_DEATH, data = GSE208338_pheno_curated)
fit = lmFit(expression_GSE208338, Design.matrix)
fitE = eBayes(fit)
GSE208338_Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE208338_Top_table$ID = rownames(GSE208338_Top_table)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE208338_Top_table$ID]
GSE208338_Top_table$SE = SE

# Annotating results
GSE208338_probes = GSE208338_probes[GSE208338_Top_table$ID, ]
GSE208338_Top_table$Gene_symbol = GSE208338_probes$Gene_symbol
GSE208338_Top_table$Gene_symbol_non_hgnc = GSE208338_probes$Gene_symbol_non_hgnc
GSE208338_Top_table$Tissue = GSE208338_pheno_curated$TISSUE[1]
GSE208338_Top_table$Tissue_type = "Brain"
GSE208338_Top_table$Technology = "Array"

# Saving output
dir.create("GSE208338_results")
write.csv(expression_GSE208338, "GSE208338_results/GSE208338_expression.csv")
write.csv(GSE208338_Top_table, "GSE208338_results/GSE208338_Top_table.csv")
write.csv(GSE208338_probes, "GSE208338_results/GSE208338_probes.csv")
write.csv(GSE208338_pheno_curated, "GSE208338_results/GSE208338_pheno_curated.csv")
rm(list = ls(pattern = "GSE208338"))
gc()


# analysis (no covariates)
GSE208338_expression = smart_fread("GSE208338_results/GSE208338_expression.csv")
GSE208338_probes = smart_fread("GSE208338_results/GSE208338_probes.csv")
GSE208338_pheno_curated = smart_fread("GSE208338_results/GSE208338_pheno_curated.csv")

GSE208338_probes = GSE208338_probes[rownames(GSE208338_expression),]
all(rownames(GSE208338_expression) == rownames(GSE208338_probes)) # TRUE
all(colnames(GSE208338_expression) == GSE208338_pheno_curated$GEO_ACCESSION) # TRUE

GSE208338_pheno_curated$SUICIDE = factor(GSE208338_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))

# Differential expression (no covar)
Design.matrix = model.matrix(~ SUICIDE, data = GSE208338_pheno_curated)
fit = lmFit(GSE208338_expression, Design.matrix)
fitE = eBayes(fit)
GSE208338_Top_table_no_covar = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE208338_Top_table_no_covar$ID = rownames(GSE208338_Top_table_no_covar)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE208338_Top_table_no_covar$ID]
GSE208338_Top_table_no_covar$SE = SE

# Annotating results
GSE208338_probes = GSE208338_probes[GSE208338_Top_table_no_covar$ID, ]
GSE208338_Top_table_no_covar$Gene_symbol = GSE208338_probes$Gene_symbol
GSE208338_Top_table_no_covar$Gene_symbol_non_hgnc = GSE208338_Top_table_no_covar$Gene_symbol_non_hgnc
GSE208338_Top_table_no_covar$Tissue = GSE208338_pheno_curated$TISSUE[1]
GSE208338_Top_table_no_covar$Tissue_type = "Brain"
GSE208338_Top_table_no_covar$Technology = "Array"

write.csv(GSE208338_Top_table_no_covar, "GSE208338_results/GSE208338_Top_table_no_covar.csv")
rm(list = ls(pattern = "GSE208338"))
gc()

################### GSE5388 ###################
GSE5388_geo_record = getGEO("GSE5388")
GSE5388_geo_record = GSE5388_geo_record[[1]]
GSE5388_pheno = pData(GSE5388_geo_record)
GSE5388_pheno = fix_columns(GSE5388_pheno)

# Variables used in the paper https://www.nature.com/articles/4001875

"
Isolated total RNA was carried through the Affymetrix preparation protocol19 and hybridized to HG-U133A GeneChips (Affymetrix). 
Affymetrix Microarray Suite 5 (MAS5) was used for image processing and data acquisition and the raw data were then subjected 
to our stringent quality control (QC) procedures.18, 20 For the DLPFC, 30 bipolar and 31 control samples were included in the 
final analysis. For the OFC, 10 bipolar and 11 control samples were included in the final analysis. Expression measures were
computed using the robust multi-chip average (RMA) method.21 Finally, the gene expression matrix was filtered to exclude 
control probesets and to include only probesets that had a ‘class A’ assignment, that is, at least nine of its 11 probes
overlapped the target transcript. This filtering step reduced the gene expression matrix to 19 537 probesets.

"
# Separate cohorts were used for the two different brain regions.

"
Quantitative variables that were significantly different between the disease and control groups and showed correlation with 
gene expression were considered as potential confounds. Drug treatment, brain pH and post-mortem interval (PMI) met these criteria 
in at least one brain region (Table 1). PMI, in contrast to pH, does not influence mRNA levels greatly25, 26 and was not included 
in the ANCOVA models to maximize the degrees of freedom.

"

# Age, Gender, PH, PMI; treatment should be also used
table(GSE5388_pheno$FLUPHENAZINE_MG__EQUIVALENTS, useNA = "always") # 31 is NA
table(GSE5388_pheno$ELECTROCONVULSIVE_THERAPY) # 31 is NA
table(GSE5388_pheno$VALPROATE_TREATMENT) # 31 is NA
table(GSE5388_pheno$DISEASE_STATUS)
"Bipolar disorder  Healthy control 
              30               31 "
GSE5388_pheno$TISSUE = "DLPFC"

# BASED on TABLE 1 -> NA is 0 and indicates NO

GSE5388_pheno$FLUPHENAZINE_MG__EQUIVALENTS_CURATED = ifelse(GSE5388_pheno$FLUPHENAZINE_MG__EQUIVALENTS == "N/A",
                                                            "0", GSE5388_pheno$FLUPHENAZINE_MG__EQUIVALENTS)

GSE5388_pheno$ELECTROCONVULSIVE_THERAPY_CURATED = sapply(GSE5388_pheno$ELECTROCONVULSIVE_THERAPY, function(x){
  
  if (x == "Yes (60 treatments)"){
    return("Yes")
  }
  
  if (x == "N/A"){
    return("No")
  }
  
  return(x)
  
})

GSE5388_pheno$VALPROATE_TREATMENT_CURATED = ifelse(GSE5388_pheno$VALPROATE_TREATMENT == "N/A",
                                                            "No", GSE5388_pheno$VALPROATE_TREATMENT)

GSE5388_pheno$LITHIUM_TREATMENT_CURATED = ifelse(GSE5388_pheno$LITHIUM_TREATMENT == "N/A",
                                                 "No", GSE5388_pheno$LITHIUM_TREATMENT)

table(GSE5388_pheno$FLUPHENAZINE_MG__EQUIVALENTS_CURATED)
# GSM123189 FLUPHENAZINE_MG__EQUIVALENTS_CURATED -> unknown
# GSM123189 ELECTROCONVULSIVE_THERAPY -> yes (1 out of 2)
table(GSE5388_pheno$VALPROATE_TREATMENT_CURATED)

# GSM123189 is removed during analysis as it has no data on fluphenazine, and ELECTROCONVULSIVE_THERAPY is not included (1 level)
GSE5388_pheno_curated = GSE5388_pheno
# Age, Gender, PH, PMI; treatment should be also used

GSE5388_pheno_curated$AGE__YEARS_ = as.numeric(GSE5388_pheno_curated$AGE__YEARS_)
GSE5388_pheno_curated$GENDER = factor(GSE5388_pheno_curated$GENDER, levels = c("Female", "Male"))
GSE5388_pheno_curated$BRAIN_PH = as.numeric(GSE5388_pheno_curated$BRAIN_PH)
GSE5388_pheno_curated$POST_MORTEM_INTERVAL__HOURS_ = as.numeric(GSE5388_pheno_curated$POST_MORTEM_INTERVAL__HOURS_)
GSE5388_pheno_curated$SUICIDE = factor(GSE5388_pheno_curated$SUICIDE, levels = c("No", "Yes"))
GSE5388_pheno_curated$FLUPHENAZINE_MG__EQUIVALENTS_CURATED = as.numeric(GSE5388_pheno_curated$FLUPHENAZINE_MG__EQUIVALENTS_CURATED) # 1 NA is added for unknown
GSE5388_pheno_curated$VALPROATE_TREATMENT_CURATED = factor(GSE5388_pheno_curated$VALPROATE_TREATMENT_CURATED, levels = c("No", "Yes"))
GSE5388_pheno_curated$LITHIUM_TREATMENT_CURATED = factor(GSE5388_pheno_curated$LITHIUM_TREATMENT_CURATED, levels = c("No", "Yes"))
GSE5388_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED = stri_replace_all_fixed(GSE5388_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE, 
                                                                                    pattern = "0 (no use) to 6 (heavy use)): ",
                                                                                    replacement = "")
GSE5388_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED = as.numeric(GSE5388_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED) # 1 NA is added for unknown
GSE5388_pheno_curated$DISEASE_STATUS = factor(GSE5388_pheno_curated$DISEASE_STATUS, levels = c("Healthy control", "Bipolar disorder"))

GSE5388_expression = exprs(GSE5388_geo_record)
GSE5388_probes = featureData(GSE5388_geo_record)
GSE5388_probes = GSE5388_probes@data

# annotation from biomart is incorrect here, example: 211600_at -> using annotation from GEO/manufacturer (it is the same)
# this is supported by UCSC and https://www.thermofisher.com/order/catalog/product/900469?SID=srch-srp-900469 supp files

probes_affy_hg_u133a_2 = smart_fread("affy_hg_u133a_2_annot.txt", skip=16)
affy_hg_u133a_genes = unlist(stri_split_fixed(probes_affy_hg_u133a_2$`Gene Symbol`, pattern = " /// "))
affy_hg_u133a_genes = unique(affy_hg_u133a_genes)
affy_hg_u133a_genes = affy_hg_u133a_genes[affy_hg_u133a_genes != ""]
affy_hg_u133a_genes = affy_hg_u133a_genes[!is.na(affy_hg_u133a_genes)]
affy_hg_u133a_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = affy_hg_u133a_genes, 
                                                  PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                  PRF_replace_NA_with_old = TRUE)
affy_hg_u133a_genes_check_dict = affy_hg_u133a_genes_check$Suggested.Symbol
names(affy_hg_u133a_genes_check_dict) = affy_hg_u133a_genes_check$x


probes_affy_hg_u133a_2$Gene_symbol = sapply(probes_affy_hg_u133a_2$`Gene Symbol`, function(x){
  x = unlist(stri_split_fixed(x, pattern = " /// "))
  x = unique(x)
  x = x[x != ""]
  x = x[!is.na(x)]
  
  if (length(x) < 1){
    return(NA)
  }
  x = sapply(x, function(z) affy_hg_u133a_genes_check_dict[z])
  x = paste0(x, collapse = ";")
  return(x)
})
write.csv(probes_affy_hg_u133a_2, "probes_affy_hg_u133a_2.csv")

GSE5388_probes$Gene_symbol = sapply(GSE5388_probes$ID, function(x){
  x = probes_affy_hg_u133a_2[probes_affy_hg_u133a_2$ID == x, "Gene_symbol"]
  return(x)
})

# getting expression values
# Getting supplementary files
getGEOSuppFiles("GSE5388")
files = list.files("GSE5388")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE5388", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE5388_CEL"))

# Inspecting and preparing CEL files
files = list.files("GSE5388_CEL")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE5388_CEL", "/", files)
# All files are fixed
GSE5388_rawdata = affy::ReadAffy(filenames = files)

# QC
affy::hist(GSE5388_rawdata)
affy::boxplot(GSE5388_rawdata) # To preview boxplots

GSE5388_expression = affy::rma(GSE5388_rawdata)

affy::boxplot(GSE5388_expression)
affy::hist(GSE5388_expression)

# everything is normalized
degradation = affy::AffyRNAdeg(GSE5388_rawdata)
colors_deg = map2color(degradation$slope,rainbow(200),limits=c(1,10))
plotAffyRNAdeg(degradation, transform = "shift.scale", cols = colors_deg)


degradation$sample.names[which(colors_deg == "#00A3FF")] # "GSM123233.cel.gz"
degradation$slope[which(colors_deg == "#00A3FF")] #  6.07362
GSE5388_pheno_curated[GSE5388_pheno_curated$GEO_ACCESSION == "GSM123233", ] # it is a control
# 1 sample seems to be odd but was kept in the original manuscript -> we remove it as it is control anyways
GSE5388_pheno_curated_2 = GSE5388_pheno_curated[GSE5388_pheno_curated$GEO_ACCESSION != "GSM123233", ]

colnames(GSE5388_expression) = stri_replace_all_fixed(colnames(GSE5388_expression), pattern = ".cel.gz", replacement = "")
all(rownames(GSE5388_expression) == GSE5388_probes$ID) # FALSE
all(rownames(GSE5388_expression) %in% GSE5388_probes$ID) # TRUE
all(GSE5388_pheno_curated_2$GEO_ACCESSION %in% colnames(GSE5388_expression)) # TRUE

# matching order
GSE5388_expression = GSE5388_expression[GSE5388_probes$ID,]
GSE5388_expression = GSE5388_expression[,GSE5388_pheno_curated_2$GEO_ACCESSION]
dim(GSE5388_expression) # 60 samples, 22283 probes

# check
all(rownames(GSE5388_expression) == GSE5388_probes$ID) # TRUE
all(colnames(GSE5388_expression) == GSE5388_pheno_curated_2$GEO_ACCESSION) # TRUE

# Differential expression (no covar)
Design.matrix = model.matrix(~ SUICIDE, data = GSE5388_pheno_curated_2)
fit = lmFit(GSE5388_expression, Design.matrix)
fitE = eBayes(fit)
GSE5388_Top_table_no_covar = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE5388_Top_table_no_covar$ID = rownames(GSE5388_Top_table_no_covar)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE5388_Top_table_no_covar$ID]
GSE5388_Top_table_no_covar$SE = SE

# Annotating results
GSE5388_probes = GSE5388_probes[GSE5388_Top_table_no_covar$ID, ]
GSE5388_Top_table_no_covar$Gene_symbol = GSE5388_probes$Gene_symbol
GSE5388_Top_table_no_covar$Tissue = GSE5388_pheno_curated_2$TISSUE[1]
GSE5388_Top_table_no_covar$Tissue_type = "Brain"
GSE5388_Top_table_no_covar$Technology = "Array"

# Differential expression (with covar)
table(GSE5388_pheno_curated_2$DRUG_ABUSE__RATINGS_SCALE, GSE5388_pheno_curated_2$SUICIDE) # 29 are unknown
table(GSE5388_pheno_curated_2$ALCOHOL_ABUSE__RATINGS_SCALE, GSE5388_pheno_curated_2$SUICIDE) # 1 is unknown


# drug abuse is ommitted due to high number of missing values
# alcohol abuse is curated 
table(GSE5388_pheno_curated_2$DISEASE_STATUS, GSE5388_pheno_curated_2$SUICIDE)
"No Yes
  Healthy control  30   0
  Bipolar disorder 18  12
"

Design.matrix = model.matrix(~ SUICIDE + DISEASE_STATUS + GENDER + AGE__YEARS_ + BRAIN_PH + POST_MORTEM_INTERVAL__HOURS_ + FLUPHENAZINE_MG__EQUIVALENTS_CURATED +
                               VALPROATE_TREATMENT_CURATED + LITHIUM_TREATMENT_CURATED + ALCOHOL_ABUSE__RATINGS_SCALE_CURATED, data = GSE5388_pheno_curated_2)
# 58 participants! 
GSE5388_pheno_curated_2$GEO_ACCESSION[GSE5388_pheno_curated_2$GEO_ACCESSION %!in% rownames(Design.matrix)] # "GSM123186" "GSM123189" are excluded

fit = lmFit(GSE5388_expression[,rownames(Design.matrix)], Design.matrix)
fitE = eBayes(fit)
GSE5388_Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE5388_Top_table$ID = rownames(GSE5388_Top_table)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE5388_Top_table$ID]
GSE5388_Top_table$SE = SE

# Annotating results
GSE5388_probes = GSE5388_probes[GSE5388_Top_table$ID, ]
GSE5388_Top_table$Gene_symbol = GSE5388_probes$Gene_symbol
GSE5388_Top_table$Tissue = GSE5388_pheno_curated_2$TISSUE[1]
GSE5388_Top_table$Tissue_type = "Brain"
GSE5388_Top_table$Technology = "Array"

# Saving outputs
dir.create("GSE5388_results")
write.csv(GSE5388_expression, "GSE5388_results/GSE5388_expression.csv")
write.csv(GSE5388_probes, "GSE5388_results/GSE5388_probes.csv")
write.csv(GSE5388_pheno_curated_2, "GSE5388_results/GSE5388_pheno_curated.csv") # "GSM123186" "GSM123189" are excluded (phenotypes); GSM123233 is excluded (degradation)
write.csv(GSE5388_Top_table_no_covar, "GSE5388_results/GSE5388_Top_table_no_covar.csv")
write.csv(GSE5388_Top_table, "GSE5388_results/GSE5388_Top_table.csv")
rm(list = ls(pattern = "GSE5388"))
gc()


################### GSE5389 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5389

# Variables used in the paper https://www.nature.com/articles/4001875
# Separate cohorts were used for the two different brain regions.
"Separate cohorts were used for the two different brain regions" # samples are independent and could be viewed as "true" replicates

GSE5389_geo_record = getGEO("GSE5389")
GSE5389_geo_record = GSE5389_geo_record[[1]]
GSE5389_pheno = pData(GSE5389_geo_record)
GSE5389_pheno = fix_columns(GSE5389_pheno)

# Age, Gender, PH, PMI; treatment should be also used but many are missing...
table(GSE5389_pheno$FLUPHENAZINE_MG__EQUIVALENTS, useNA = "always") # 11 is NA
table(GSE5389_pheno$ELECTROCONVULSIVE_THERAPY) # 11 is NA
table(GSE5389_pheno$VALPROATE_TREATMENT) # 11 is NA

GSE5389_pheno$TISSUE = "Orbitofrontal cortex"

# BASED on TABLE 1 -> NA is 0 and indicates NO
# BASED on TABLE 1 -> NA is 0 and indicates NO

GSE5389_pheno$FLUPHENAZINE_MG__EQUIVALENTS_CURATED = ifelse(GSE5389_pheno$FLUPHENAZINE_MG__EQUIVALENTS == "N/A",
                                                            "0", GSE5389_pheno$FLUPHENAZINE_MG__EQUIVALENTS)

GSE5389_pheno$ELECTROCONVULSIVE_THERAPY_CURATED = sapply(GSE5389_pheno$ELECTROCONVULSIVE_THERAPY, function(x){
  
  if (x == "Yes (5 treatments)"){
    return("Yes")
  }
  
  if (x == "N/A"){
    return("No")
  }
  
  return(x)
  
})

GSE5389_pheno$VALPROATE_TREATMENT_CURATED = ifelse(GSE5389_pheno$VALPROATE_TREATMENT == "N/A",
                                                   "No", GSE5389_pheno$VALPROATE_TREATMENT)

GSE5389_pheno$LITHIUM_TREATMENT_CURATED = ifelse(GSE5389_pheno$LITHIUM_TREATMENT == "N/A",
                                                 "No", GSE5389_pheno$LITHIUM_TREATMENT)

table(GSE5389_pheno$FLUPHENAZINE_MG__EQUIVALENTS_CURATED)
table(GSE5389_pheno$ELECTROCONVULSIVE_THERAPY_CURATED)
table(GSE5389_pheno$VALPROATE_TREATMENT_CURATED)
table(GSE5389_pheno$LITHIUM_TREATMENT_CURATED)
# all phenotypes are OK!

GSE5389_pheno_curated = GSE5389_pheno
# Age, Gender, PH, PMI; treatment should be also used but many are missing...

GSE5389_pheno_curated$AGE__YEARS_ = as.numeric(GSE5389_pheno_curated$AGE__YEARS_)
GSE5389_pheno_curated$GENDER = factor(GSE5389_pheno_curated$GENDER, levels = c("Female", "Male"))
GSE5389_pheno_curated$BRAIN_PH = as.numeric(GSE5389_pheno_curated$BRAIN_PH)
GSE5389_pheno_curated$POST_MORTEM_INTERVAL__HOURS_ = as.numeric(GSE5389_pheno_curated$POST_MORTEM_INTERVAL__HOURS_)
GSE5389_pheno_curated$SUICIDE = factor(GSE5389_pheno_curated$SUICIDE, levels = c("No", "Yes"))
GSE5389_pheno_curated$FLUPHENAZINE_MG__EQUIVALENTS_CURATED = as.numeric(GSE5389_pheno_curated$FLUPHENAZINE_MG__EQUIVALENTS_CURATED)
GSE5389_pheno_curated$VALPROATE_TREATMENT_CURATED = factor(GSE5389_pheno_curated$VALPROATE_TREATMENT_CURATED, levels = c("No", "Yes"))
GSE5389_pheno_curated$LITHIUM_TREATMENT_CURATED = factor(GSE5389_pheno_curated$LITHIUM_TREATMENT_CURATED, levels = c("No", "Yes"))
GSE5389_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED = stri_replace_all_fixed(GSE5389_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE, 
                                                                                    pattern = "0 (no use) to 6 (heavy use)): ",
                                                                                    replacement = "")
GSE5389_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED = as.numeric(GSE5389_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED)

GSE5389_pheno_curated$DRUG_ABUSE__RATINGS_SCALE_CURATED = stri_replace_all_fixed(GSE5389_pheno_curated$DRUG_ABUSE__RATINGS_SCALE, 
                                                                                    pattern = "0 (no use) to 5 (heavy use)): ",
                                                                                    replacement = "")
GSE5389_pheno_curated$DRUG_ABUSE__RATINGS_SCALE_CURATED = as.numeric(GSE5389_pheno_curated$DRUG_ABUSE__RATINGS_SCALE_CURATED)
GSE5389_pheno_curated$DISEASE_STATUS = factor(GSE5389_pheno_curated$DISEASE_STATUS, levels = c("Healthy control", "Bipolar disorder"))

GSE5389_expression = exprs(GSE5389_geo_record)
GSE5389_probes = featureData(GSE5389_geo_record)
GSE5389_probes = GSE5389_probes@data

# annotation from biomart is incorrect here, example: 211600_at -> using annotation from GEO/manufacturer (it is the same)
# this is supported by UCSC and https://www.thermofisher.com/order/catalog/product/900469?SID=srch-srp-900469 supp files

# reading the same annotation as in the previous cohort
probes_affy_hg_u133a_2 = read.csv("probes_affy_hg_u133a_2.csv")
probes_affy_hg_u133a_2$X = NULL
rownames(probes_affy_hg_u133a_2) = probes_affy_hg_u133a_2$ID

GSE5389_probes$Gene_symbol = sapply(GSE5389_probes$ID, function(x){
  x = probes_affy_hg_u133a_2[probes_affy_hg_u133a_2$ID == x, "Gene_symbol"]
  return(x)
})


# getting expression values
# Getting supplementary files
getGEOSuppFiles("GSE5389")
files = list.files("GSE5389")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE5389", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE5389_CEL"))

# Inspecting and preparing CEL files
files = list.files("GSE5389_CEL")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE5389_CEL", "/", files)
# All files are fixed
GSE5389_rawdata = affy::ReadAffy(filenames = files)

# QC
affy::hist(GSE5389_rawdata) # To preview density
affy::boxplot(GSE5389_rawdata) # To preview boxplots

GSE5389_expression = affy::rma(GSE5389_rawdata)

affy::boxplot(GSE5389_expression)
affy::hist(GSE5389_expression)

# everything is normalized and seems OK!
degradation = affy::AffyRNAdeg(GSE5389_rawdata)
colors_deg = map2color(degradation$slope,rainbow(200),limits=c(1,10))
plotAffyRNAdeg(degradation, transform = "shift.scale", cols = colors_deg)

colnames(GSE5389_expression) = stri_replace_all_fixed(colnames(GSE5389_expression), pattern = ".cel.gz", replacement = "")
all(rownames(GSE5389_expression) == GSE5389_probes$ID) # FALSE
all(rownames(GSE5389_expression) %in% GSE5389_probes$ID) # TRUE
all(GSE5389_pheno_curated$GEO_ACCESSION %in% colnames(GSE5389_expression)) # TRUE

# matching order
GSE5389_expression = GSE5389_expression[GSE5389_probes$ID,]
GSE5389_expression = GSE5389_expression[,GSE5389_pheno_curated$GEO_ACCESSION]
dim(GSE5389_expression) # 21 samples, 22283 probes

# check
all(rownames(GSE5389_expression) == GSE5389_probes$ID) # TRUE
all(colnames(GSE5389_expression) == GSE5389_pheno_curated$GEO_ACCESSION) # TRUE


# Differential expression (no covar)
Design.matrix = model.matrix(~ SUICIDE, data = GSE5389_pheno_curated)
fit = lmFit(GSE5389_expression, Design.matrix)
fitE = eBayes(fit)
GSE5389_Top_table_no_covar = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE5389_Top_table_no_covar$ID = rownames(GSE5389_Top_table_no_covar)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE5389_Top_table_no_covar$ID]
GSE5389_Top_table_no_covar$SE = SE

# Annotating results
GSE5389_probes = GSE5389_probes[GSE5389_Top_table_no_covar$ID, ]
GSE5389_Top_table_no_covar$Gene_symbol = GSE5389_probes$Gene_symbol
GSE5389_Top_table_no_covar$Tissue = GSE5389_pheno_curated$TISSUE[1]
GSE5389_Top_table_no_covar$Tissue_type = "Brain"
GSE5389_Top_table_no_covar$Technology = "Array"


# Differential expression (with covar)
# drug abuse is included
# alcohol abuse is included 
# ECT is included
table(GSE5389_pheno_curated$DISEASE_STATUS, GSE5389_pheno_curated$SUICIDE)
"
                   No Yes
  Healthy control  11   0
  Bipolar disorder  4   6
"
Design.matrix = model.matrix(~ SUICIDE + DISEASE_STATUS + GENDER + AGE__YEARS_ + BRAIN_PH + POST_MORTEM_INTERVAL__HOURS_ + FLUPHENAZINE_MG__EQUIVALENTS_CURATED +
                               VALPROATE_TREATMENT_CURATED + LITHIUM_TREATMENT_CURATED + ELECTROCONVULSIVE_THERAPY_CURATED +
                               ALCOHOL_ABUSE__RATINGS_SCALE_CURATED + DRUG_ABUSE__RATINGS_SCALE_CURATED, data = GSE5389_pheno_curated)

fit = lmFit(GSE5389_expression, Design.matrix)
fitE = eBayes(fit)
GSE5389_Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE5389_Top_table$ID = rownames(GSE5389_Top_table)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE5389_Top_table$ID]
GSE5389_Top_table$SE = SE

# Annotating results
GSE5389_probes = GSE5389_probes[GSE5389_Top_table$ID, ]
GSE5389_Top_table$Gene_symbol = GSE5389_probes$Gene_symbol
GSE5389_Top_table$Tissue = GSE5389_pheno_curated$TISSUE[1]
GSE5389_Top_table$Tissue_type = "Brain"
GSE5389_Top_table$Technology = "Array"

# Saving outputs
dir.create("GSE5389_results")
write.csv(GSE5389_expression, "GSE5389_results/GSE5389_expression.csv")
write.csv(GSE5389_probes, "GSE5389_results/GSE5389_probes.csv")
write.csv(GSE5389_pheno_curated, "GSE5389_results/GSE5389_pheno_curated.csv")
write.csv(GSE5389_Top_table_no_covar, "GSE5389_results/GSE5389_Top_table_no_covar.csv")
write.csv(GSE5389_Top_table, "GSE5389_results/GSE5389_Top_table.csv")
rm(list = ls(pattern = "GSE5389"))
gc()

################### GSE66937 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE66937
# no original paper and has technical replicates

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8458545/ this is original paper

# affy_hta_2_0 probe mapping with http://grch37.ensembl.org/biomart

# Variables used in the paper (nothing found)
GSE66937_geo_record = getGEO("GSE66937")
GSE66937_geo_record = GSE66937_geo_record[[1]]
GSE66937_pheno = pData(GSE66937_geo_record)
GSE66937_pheno = fix_columns(GSE66937_pheno)


GSE66937_pheno$BRAIN_REGION
"         amygdala       hippocampus prefrontal cortex          thalamus 
               17                16                17                17
"
GSE66937_pheno$GROUP # "suicide victim" "control person"
GSE66937_pheno$RIN

GSE66937_pheno_curated = GSE66937_pheno

# inspecting technical replicates (will be removed before analysis, used for better normalization)
table(stri_detect_fixed(GSE66937_pheno_curated$TITLE, pattern = "technical replicate"))
"
FALSE  TRUE 
   59     8 
"


GSE66937_pheno_curated$SUICIDE = ifelse(GSE66937_pheno_curated$GROUP == "suicide victim", "Suicide", "Control")
GSE66937_pheno_curated$SUICIDE = factor(GSE66937_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
GSE66937_pheno_curated$RIN = as.numeric(GSE66937_pheno_curated$RIN)

# Getting expression data
# Getting supplementary files
getGEOSuppFiles("GSE66937")
files = list.files("GSE66937")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE66937", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE66937_CEL"))

# Inspecting and preparing CEL files
files = list.files("GSE66937_CEL")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE66937_CEL", "/", files)

# files_2 = files[stri_detect_fixed(files, pattern = "Thalamus")]

# weird distibution regardless of tissue
# All files are fixed
GSE66937_rawdata = oligo::read.celfiles(files, pkgname = "pd.hta.2.0")

# QC
oligo::hist(x = GSE66937_rawdata, col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
oligo::boxplot(GSE66937_rawdata, col = darkColors(16), target = "probeset")

# sample names
GSE66937_samples = GSE66937_rawdata@phenoData@data
GSE66937_samples$GEO_ACCESSION = sapply(rownames(GSE66937_samples), function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
cases = which(GSE66937_samples$GEO_ACCESSION %in% GSE66937_pheno_curated[GSE66937_pheno_curated$SUICIDE == "Suicide", "GEO_ACCESSION"])
controls = which(GSE66937_samples$GEO_ACCESSION %in% GSE66937_pheno_curated[GSE66937_pheno_curated$SUICIDE != "Suicide", "GEO_ACCESSION"])

oligo::hist(x = GSE66937_rawdata[,cases], col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
oligo::hist(x = GSE66937_rawdata[,controls], col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
# difference is not between cases and controls ...

# it seems like batch from 2014 has different distribution...
GSE66937_pheno_curated$BATCH = sapply(GSE66937_pheno_curated$HYBRIDIZIATION_DATE, function(x){
  if (stri_detect_fixed(x, pattern = "2013")){
    return("2013")
  }
  
  if (stri_detect_fixed(x, pattern = "2014")){
    return("2014")
  }
  
})
table(GSE66937_pheno_curated$BATCH)
"
2013 2014 
  40   27 
"

GSE66937_2013 = which(GSE66937_samples$GEO_ACCESSION %in% GSE66937_pheno_curated[GSE66937_pheno_curated$BATCH == "2013", "GEO_ACCESSION"])
GSE66937_2014 = which(GSE66937_samples$GEO_ACCESSION %in% GSE66937_pheno_curated[GSE66937_pheno_curated$BATCH == "2014", "GEO_ACCESSION"])
oligo::hist(x = GSE66937_rawdata[,GSE66937_2013], col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
oligo::hist(x = GSE66937_rawdata[,GSE66937_2014], col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")

# samples from 2014 have completely different distribution. Reason: RNA degradation?
dir.create("GSE66937_results")

png(filename = "GSE66937_results/batch_2013.png", width = 1920, height = 1080, units = "px")
oligo::hist(x = GSE66937_rawdata[,GSE66937_2013], col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
dev.off()

png(filename = "GSE66937_results/batch_2014.png", width = 1920, height = 1080, units = "px")
oligo::hist(x = GSE66937_rawdata[,GSE66937_2014], col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
dev.off()

mean(GSE66937_pheno_curated[GSE66937_pheno_curated$BATCH == "2013", "RIN"]) # 6.6375
mean(GSE66937_pheno_curated[GSE66937_pheno_curated$BATCH == "2014", "RIN"]) # 7.307407
GSE66937_pheno_curated$BATCH = factor(GSE66937_pheno_curated$BATCH, levels = c("2013", "2014"))


GSE66937_expression = oligo::rma(GSE66937_rawdata)

oligo::hist(x = GSE66937_expression, transfo="identity", col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
oligo::boxplot(x=GSE66937_expression,  transfo="identity", col = darkColors(16), target = "probeset")
# it seems normalized but now follows distribution from 2014!

table(GSE66937_pheno_curated$BATCH, GSE66937_pheno_curated$BRAIN_REGION)

"
amygdala hippocampus prefrontal cortex thalamus
  2013       10          10                10       10
  2014        7           6                 7        7

"
table(GSE66937_pheno_curated$HYBRIDIZIATION_DATE, GSE66937_pheno_curated$BRAIN_REGION)

"             amygdala hippocampus prefrontal cortex thalamus
  2013-11-28        5           5                 5        5
  2013-12-03        5           5                 5        5
  2014-01-12        4           3                 4        4
  2014-01-17        3           3                 3        3

"

table(GSE66937_pheno_curated$HYBRIDIZIATION_DATE, GSE66937_pheno_curated$BRAIN_REGION, GSE66937_pheno_curated$SUICIDE)

"
, ,  = Control

             amygdala hippocampus prefrontal cortex thalamus
  2013-11-28        2           2                 2        2
  2013-12-03        3           3                 3        3
  2014-01-12        1           1                 1        1
  2014-01-17        1           1                 1        1

, ,  = Suicide

            
             amygdala hippocampus prefrontal cortex thalamus
  2013-11-28        3           3                 3        3
  2013-12-03        2           2                 2        2
  2014-01-12        3           2                 3        3
  
"

# affy_hta_2_0 probe mapping with http://grch37.ensembl.org/biomart

GSE66937_probes = featureData(GSE66937_geo_record)
GSE66937_probes = GSE66937_probes@data

#  Use BioMart to get the gene symbols (performed in python) affy_hta_2_0
# http://grch37.ensembl.org/biomart

# mapping probe IDs to gene symbols
probes_affy_hta2 = as.data.frame(fread("mapping_affy_hta_2_0_enseml.txt", nThread = 10, sep ="\t", header = FALSE))
all(GSE66937_probes$ID %in% probes_affy_hta2$V4) # FALSE
all(probes_affy_hta2$V4 %in% GSE66937_probes$ID) # FALSE

table(GSE66937_probes$ID %in% probes_affy_hta2$V4) # FALSE 70523
table(probes_affy_hta2$V4 %in% GSE66937_probes$ID)  # FALSE 117287

GSE66937_probes$CURATED_ID = stri_replace_all_fixed(GSE66937_probes$ID, pattern = ".hg.1", replacement = ".hg")

any(stri_detect_fixed(GSE66937_probes$CURATED_ID, pattern = ".hg.1")) # FALSE

table(GSE66937_probes$CURATED_ID %in% probes_affy_hta2$V4) # FALSE 19296 TRUE 51227
table(probes_affy_hta2$V4 %in% GSE66937_probes$CURATED_ID)  # FALSE 41898 TRUE 75389

GSE66937_mismatchd_IDs = GSE66937_probes[GSE66937_probes$CURATED_ID %!in% probes_affy_hta2$V4, ]
table(GSE66937_mismatchd_IDs$`locus type`)
"
    Coding   control->affx->bac_spike        control->affx->ercc  control->affx->ercc->step control->affx->polya_spike  control->bgp->antigenomic      main///normgene->exon 
      5402                          4                         92                         63                          4                         23                       1465 
 NonCoding             normgene->exon           normgene->intron 
     10899                        698                        646 
"

GSE66937_matched_IDs = GSE66937_probes[GSE66937_probes$CURATED_ID %in% probes_affy_hta2$V4, ]
table(GSE66937_matched_IDs$`locus type`)
"
   Coding NonCoding 
    39297     11930 
"
GSE66937_probes_ids = GSE66937_probes$CURATED_ID

GSE66937_probes_gene_symbols = mclapply(GSE66937_probes_ids, function(x){
  
  if (x %!in% probes_affy_hta2$V4){
    return(c(NA, NA))
  }
  
  gene_names = probes_affy_hta2[probes_affy_hta2$V4 == x, "V5"]
  gene_names = str_trim(gene_names)
  gene_names = unique(gene_names)
  gene_names = gene_names[gene_names!=""]
  gene_names = gene_names[!is.na(gene_names)]
  gene_names = paste0(gene_names, collapse = ";")
  
  gene_names_non_HGNC = probes_affy_hta2[probes_affy_hta2$V4 == x, "V2"]
  gene_names_non_HGNC = str_trim(gene_names_non_HGNC)
  gene_names_non_HGNC = unique(gene_names_non_HGNC)
  gene_names_non_HGNC = gene_names_non_HGNC[gene_names_non_HGNC!=""]
  gene_names_non_HGNC = gene_names_non_HGNC[!is.na(gene_names_non_HGNC)]
  gene_names_non_HGNC = paste0(gene_names_non_HGNC, collapse = ";")
  
  return(c(gene_names, gene_names_non_HGNC))
}, mc.cores = 6)


GSE66937_probes$Gene_symbol = sapply(GSE66937_probes_gene_symbols, function(x){
  return(x[1])
})

GSE66937_probes$Gene_symbol_non_hgnc = sapply(GSE66937_probes_gene_symbols, function(x){
  return(x[2])
})

colnames(GSE66937_expression) = sapply(colnames(GSE66937_expression), function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})
all(rownames(GSE66937_expression) == GSE66937_probes$ID) # FALSE
all(rownames(GSE66937_expression) %in% GSE66937_probes$ID) # TRUE
all(GSE66937_pheno_curated$GEO_ACCESSION %in% colnames(GSE66937_expression)) # TRUE

# matching order
GSE66937_expression = GSE66937_expression[GSE66937_probes$ID,]
GSE66937_expression = GSE66937_expression[,GSE66937_pheno_curated$GEO_ACCESSION]
dim(GSE66937_expression) # 67 samples, 70523 probes

# check
all(rownames(GSE66937_expression) == GSE66937_probes$ID) # TRUE
all(colnames(GSE66937_expression) == GSE66937_pheno_curated$GEO_ACCESSION) # TRUE

# removing technical replicates!
table(stri_detect_fixed(GSE66937_pheno_curated$TITLE, pattern = "technical replicate"))
"
FALSE  TRUE 
   59     8 
"
GSE66937_pheno_curated = GSE66937_pheno_curated[!stri_detect_fixed(GSE66937_pheno_curated$TITLE, pattern =  "technical replicate"), ]
table(stri_detect_fixed(GSE66937_pheno_curated$TITLE, pattern = "technical replicate"))

"
FALSE 
   59 
"
GSE66937_expression = GSE66937_expression[,GSE66937_pheno_curated$GEO_ACCESSION]
# check
all(rownames(GSE66937_expression) == GSE66937_probes$ID) # TRUE
all(colnames(GSE66937_expression) == GSE66937_pheno_curated$GEO_ACCESSION) # TRUE

# EXPRSSION MUST BE ANALYZED IN A LOOP PER TISSUE!

# Differential expression (no covar)

GSE66937_analysis_list_no_covar = list()

for (x in 1:length(unique(GSE66937_pheno_curated$BRAIN_REGION))){
  
  TMP_tissue = unique(GSE66937_pheno_curated$BRAIN_REGION)[x]
  writeLines(paste0("Working on: ", TMP_tissue))
  
  GSE66937_TMP_DF = GSE66937_pheno_curated[GSE66937_pheno_curated$BRAIN_REGION == TMP_tissue, ]
  
  Design.matrix = model.matrix(~ SUICIDE, data = GSE66937_TMP_DF)
  fit = lmFit(GSE66937_expression[, rownames(Design.matrix)], Design.matrix)
  fitE = eBayes(fit)
  GSE66937_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE66937_Top_table_no_covar_TMP$ID = rownames(GSE66937_Top_table_no_covar_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE66937_Top_table_no_covar_TMP$ID]
  GSE66937_Top_table_no_covar_TMP$SE = SE
  
  # Annotating results
  GSE66937_probes_TMP = GSE66937_probes
  GSE66937_probes_TMP = GSE66937_probes_TMP[GSE66937_Top_table_no_covar_TMP$ID, ]
  GSE66937_Top_table_no_covar_TMP$Gene_symbol = GSE66937_probes_TMP$Gene_symbol
  GSE66937_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = GSE66937_probes_TMP$Gene_symbol_non_hgnc
  GSE66937_Top_table_no_covar_TMP$Tissue = TMP_tissue
  GSE66937_Top_table_no_covar_TMP$Tissue_type = "Brain"
  GSE66937_Top_table_no_covar_TMP$Technology = "Array"
  
  GSE66937_analysis_list_no_covar[[x]] = GSE66937_Top_table_no_covar_TMP
  
}
lapply(GSE66937_analysis_list_no_covar, head)
GSE66937_Top_table_no_covar = do.call(rbind, GSE66937_analysis_list_no_covar)
rownames(GSE66937_Top_table_no_covar) = NULL

# Differential expression (with covar)

GSE66937_analysis_list = list()

for (x in 1:length(unique(GSE66937_pheno_curated$BRAIN_REGION))){
  
  TMP_tissue = unique(GSE66937_pheno_curated$BRAIN_REGION)[x]
  writeLines(paste0("Working on: ", TMP_tissue))
  
  GSE66937_TMP_DF = GSE66937_pheno_curated[GSE66937_pheno_curated$BRAIN_REGION == TMP_tissue, ]
  
  Design.matrix = model.matrix(~ SUICIDE + RIN + BATCH, data = GSE66937_TMP_DF)
  fit = lmFit(GSE66937_expression[, rownames(Design.matrix)], Design.matrix)
  fitE = eBayes(fit)
  GSE66937_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE66937_Top_table_TMP$ID = rownames(GSE66937_Top_table_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE66937_Top_table_TMP$ID]
  GSE66937_Top_table_TMP$SE = SE
  
  # Annotating results
  GSE66937_probes_TMP = GSE66937_probes
  GSE66937_probes_TMP = GSE66937_probes_TMP[GSE66937_Top_table_TMP$ID, ]
  GSE66937_Top_table_TMP$Gene_symbol = GSE66937_probes_TMP$Gene_symbol
  GSE66937_Top_table_TMP$Gene_symbol_non_hgnc = GSE66937_probes_TMP$Gene_symbol_non_hgnc
  GSE66937_Top_table_TMP$Tissue = TMP_tissue
  GSE66937_Top_table_TMP$Tissue_type = "Brain"
  GSE66937_Top_table_TMP$Technology = "Array"
  
  GSE66937_analysis_list[[x]] = GSE66937_Top_table_TMP
  
}
lapply(GSE66937_analysis_list, head)
GSE66937_Top_table = do.call(rbind, GSE66937_analysis_list)
rownames(GSE66937_Top_table) = NULL

# Saving outputs
write.csv(GSE66937_expression, "GSE66937_results/GSE66937_expression.csv")
write.csv(GSE66937_probes, "GSE66937_results/GSE66937_probes.csv")
write.csv(GSE66937_pheno_curated, "GSE66937_results/GSE66937_pheno_curated.csv")
write.csv(GSE66937_Top_table_no_covar, "GSE66937_results/GSE66937_Top_table_no_covar.csv")
write.csv(GSE66937_Top_table, "GSE66937_results/GSE66937_Top_table.csv")
rm(list = ls(pattern = "GSE66937"))
gc()


################### GSE199536 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE199536

# Variables used in the paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9134578/

# Probe mapping with affy_hugene_2_0_st_v1 in python

"For this study, Hb tissues were prepared from 10 suicide subjects diagnosed with major depression 
and 10 psychiatrically healthy control subjects. All subjects were male Caucasians, and groups 
(suicides and controls) were matched for age, pH, and postmortem interval (PMI). Hb was carefully 
dissected, as described in our previous study [1]. Briefly, the human brain was divided into left 
and right hemispheres, and the meninges of the hemispheres were carefully removed. The trunk from 
the hemisphere was separated with a scalpel between the mammillary body and the superior colliculus. 
The pH of the sample was measured using the cerebellum. To make slabs, the brain hemispheres were 
placed with the median face down on a cutting plate and cut coronally into 18–20 pieces. Since Hb 
typically exists in the 11th or 12th slabs, Hb was gently removed from the frozen slab using a burin 
with a round-shaped tip."

"
Genes exhibiting > 1.2-fold changes at p < 0.05, as determined by Student’s t-test, were regarded 
as differentially expressed genes (DEGs). The datasets used and analyzed in the present study are 
available from the corresponding authors upon reasonable request. Raw microarray data 
have been submitted to the Gene Expression Omnibus (GEO) repository

"

"Associations between DEGs and psychiatric diseases were identified using the PsyGeNET database [2]. 
Enrichment for each psychiatric disease was analyzed with the PsyGeNET2r package. The evidence index
of PsyGeNET is expressed as the number of pieces of evidence supporting the existence of gene-disease 
associations divided by the total number of pieces of evidence. 

"

# Data seems to be paired! However, the text claims that "groups were matched"
# Extra information from the paper: 
"
All subjects were male Caucasians, and groups (suicides and controls, n = 10 for each group) were matched 
for age, postmortem interval (PMI), and pH; Age: 44.20 ± 3.35 years (mean ± SE) for controls and 40.50 ± 4.35 
for suicides; PMI: 26.45 ± 4.95 h for controls and 36.10 ± 4.32 for suicides; pH: 6.671 ± 0.047 for controls and 
6.427 ± 0.095 for suicides (see Additional file 1 for experimental procedures and Additional file 2: Table S1 
for more detailed subject information).

"
# Exctracted phenotypes from Data S1
GSE199536_DATA_S1 = read.xlsx("13041_2022_934_MOESM2_ESM.xlsx", sheet = 1, startRow = 2)
GSE199536_DATA_S1 = fix_columns(GSE199536_DATA_S1)
GSE199536_DATA_S1$GROUP = NULL
GSE199536_DATA_S1 = GSE199536_DATA_S1[!is.na(GSE199536_DATA_S1$SUBJECT_ID), ]


GSE199536_geo_record = getGEO("GSE199536")
GSE199536_geo_record = GSE199536_geo_record[[1]]
GSE199536_pheno = pData(GSE199536_geo_record)
GSE199536_pheno = fix_columns(GSE199536_pheno)

GSE199536_pheno_curated = GSE199536_pheno

GSE199536_pheno_curated$SUICIDE = ifelse(GSE199536_pheno_curated$GROUP == "suicide", "Suicide", "Control")
GSE199536_pheno_curated$SUICIDE = factor(GSE199536_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))

GSE199536_pheno_curated$ID = stri_replace_all_fixed(GSE199536_pheno_curated$TITLE, pattern = "Habenula_", replacement = "")
all(GSE199536_pheno_curated$ID == GSE199536_DATA_S1$SUBJECT_ID) # TRUE
GSE199536_pheno_curated = cbind(GSE199536_pheno_curated, GSE199536_DATA_S1)
GSE199536_pheno_curated$AGE__YR_
GSE199536_pheno_curated$PMI__H_A
GSE199536_pheno_curated$PH

# Getting expression data
# Getting supplementary files
getGEOSuppFiles("GSE199536")
files = list.files("GSE199536", full.names = TRUE)
files = files[stri_detect_fixed(files, pattern = ".tar")]
lapply(files, function(x) untar(x, exdir = "GSE199536_CEL"))

# Inspecting and preparing CEL files
files = list.files("GSE199536_CEL",  full.names = TRUE)
files = files[stri_detect_fixed(files, pattern = ".gz")]

GSE199536_rawdata = oligo::read.celfiles(files, pkgname = "pd.hugene.2.0.st")

# QC
oligo::hist(x = GSE199536_rawdata, col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
oligo::boxplot(GSE199536_rawdata, col = darkColors(16), target = "probeset")
# not identical but OK!

GSE199536_expression = oligo::rma(GSE199536_rawdata)

oligo::hist(x = GSE199536_expression, transfo="identity", col = darkColors(16), lty = 1, xlim = c(6, 12), target = "probeset")
oligo::boxplot(x=GSE199536_expression,  transfo="identity", col = darkColors(16), target = "probeset")
# perfect!

#  Use BioMart to get the gene symbols (performed in python) affy_hugene_2_0_st_v1
# http://grch37.ensembl.org/biomart

# mapping probe IDs to gene symbols
probes_affy_hugene_2_0_st_v1 = as.data.frame(fread("mapping_affy_hugene_2_0_st_v1_enseml.txt", nThread = 10, sep ="\t", header = FALSE))
annotation_affy_hugene_2_0_st_v1 =  as.data.frame(fread("HuGene-2_0-st-v1.na36.hg19.transcript.csv", nThread = 10))

all(GSE199536_probes$ID %in% annotation_affy_hugene_2_0_st_v1$probeset_id) # TRUE
all(GSE199536_probes$ID %in% annotation_affy_hugene_2_0_st_v1$transcript_cluster_id) # TRUE
all(annotation_affy_hugene_2_0_st_v1$transcript_cluster_id == annotation_affy_hugene_2_0_st_v1$probeset_id) # TRUE

all(GSE199536_probes$ID %in% probes_affy_hugene_2_0_st_v1$V4) # FALSE
all(probes_affy_hugene_2_0_st_v1$V4 %in% GSE199536_probes$ID) # FALSE

GSE199536_probes = annotation_affy_hugene_2_0_st_v1

table(GSE199536_probes$probeset_id %in% probes_affy_hugene_2_0_st_v1$V4) # FALSE 13342 TRUE 40275
table(probes_affy_hugene_2_0_st_v1$V4 %in% GSE199536_probes$probeset_id)  # FALSE 38378 TRUE 59090

GSE199536_mismatchd_IDs = GSE199536_probes[GSE199536_probes$probeset_id %!in% probes_affy_hugene_2_0_st_v1$V4, ]

GSE199536_probes_ids = GSE199536_probes$probeset_id

GSE199536_probes_gene_symbols = mclapply(GSE199536_probes_ids, function(x){
  
  if (x %!in% probes_affy_hugene_2_0_st_v1$V4){
    return(c(NA, NA))
  }
  
  gene_names = probes_affy_hugene_2_0_st_v1[probes_affy_hugene_2_0_st_v1$V4 == x, "V5"]
  gene_names = str_trim(gene_names)
  gene_names = unique(gene_names)
  gene_names = gene_names[gene_names!=""]
  gene_names = gene_names[!is.na(gene_names)]
  gene_names = paste0(gene_names, collapse = ";")
  
  gene_names_non_HGNC = probes_affy_hugene_2_0_st_v1[probes_affy_hugene_2_0_st_v1$V4 == x, "V2"]
  gene_names_non_HGNC = str_trim(gene_names_non_HGNC)
  gene_names_non_HGNC = unique(gene_names_non_HGNC)
  gene_names_non_HGNC = gene_names_non_HGNC[gene_names_non_HGNC!=""]
  gene_names_non_HGNC = gene_names_non_HGNC[!is.na(gene_names_non_HGNC)]
  gene_names_non_HGNC = paste0(gene_names_non_HGNC, collapse = ";")
  
  return(c(gene_names, gene_names_non_HGNC))
}, mc.cores = 6)

GSE199536_probes$Gene_symbol = sapply(GSE199536_probes_gene_symbols, function(x){
  return(x[1])
})

GSE199536_probes$Gene_symbol_non_hgnc = sapply(GSE199536_probes_gene_symbols, function(x){
  return(x[2])
})

colnames(GSE199536_expression) = sapply(colnames(GSE199536_expression), function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[1]
  return(x)
})

all(rownames(GSE199536_expression) == GSE199536_probes$probeset_id) # FALSE
all(rownames(GSE199536_expression) %in% GSE199536_probes$probeset_id) # TRUE
all(GSE199536_pheno_curated$GEO_ACCESSION %in% colnames(GSE199536_expression)) # TRUE
rownames(GSE199536_probes) = GSE199536_probes$probeset_id

# matching order
GSE199536_probes = GSE199536_probes[rownames(GSE199536_expression),]
GSE199536_expression = GSE199536_expression[,GSE199536_pheno_curated$GEO_ACCESSION]
dim(GSE199536_expression) # 20 samples, 53617 probes

# check
all(rownames(GSE199536_expression) == GSE199536_probes$probeset_id) # TRUE
all(colnames(GSE199536_expression) == GSE199536_pheno_curated$GEO_ACCESSION) # TRUE


# Differential expression (no covar)
Design.matrix = model.matrix(~ SUICIDE, data = GSE199536_pheno_curated)
fit = lmFit(GSE199536_expression, Design.matrix)
fitE = eBayes(fit)
GSE199536_Top_table_no_covar = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE199536_Top_table_no_covar$ID = rownames(GSE199536_Top_table_no_covar)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE199536_Top_table_no_covar$ID]
GSE199536_Top_table_no_covar$SE = SE

# Annotating results
GSE199536_probes = GSE199536_probes[GSE199536_Top_table_no_covar$ID, ]
GSE199536_Top_table_no_covar$Gene_symbol = GSE199536_probes$Gene_symbol
GSE199536_Top_table_no_covar$Gene_symbol_non_hgnc = GSE199536_probes$Gene_symbol_non_hgnc
GSE199536_Top_table_no_covar$Tissue = GSE199536_pheno_curated$TISSUE[1]
GSE199536_Top_table_no_covar$Tissue_type = "Brain"
GSE199536_Top_table_no_covar$Technology = "Array"


# Differential expression (with covar)
Design.matrix = model.matrix(~ SUICIDE + AGE__YR_ + PMI__H_A + PH, data = GSE199536_pheno_curated)

fit = lmFit(GSE199536_expression, Design.matrix)
fitE = eBayes(fit)
GSE199536_Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE199536_Top_table$ID = rownames(GSE199536_Top_table)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE199536_Top_table$ID]
GSE199536_Top_table$SE = SE

# Annotating results
GSE199536_probes = GSE199536_probes[GSE199536_Top_table$ID, ]
GSE199536_Top_table$Gene_symbol = GSE199536_probes$Gene_symbol
GSE199536_Top_table$Gene_symbol_non_hgnc = GSE199536_probes$Gene_symbol_non_hgnc
GSE199536_Top_table$Tissue = GSE199536_pheno_curated$TISSUE[1]
GSE199536_Top_table$Tissue_type = "Brain"
GSE199536_Top_table$Technology = "Array"

# Saving outputs
dir.create("GSE199536_results")
write.csv(GSE199536_expression, "GSE199536_results/GSE199536_expression.csv")
write.csv(GSE199536_probes, "GSE199536_results/GSE199536_probes.csv")
write.csv(GSE199536_pheno_curated, "GSE199536_results/GSE199536_pheno_curated.csv")
write.csv(GSE199536_Top_table_no_covar, "GSE199536_results/GSE199536_Top_table_no_covar.csv")
write.csv(GSE199536_Top_table, "GSE199536_results/GSE199536_Top_table.csv")

rm(list = ls(pattern = "GSE199536"))
gc()



################### GSE92538 U133A ###################

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92538

# paper: https://pubmed.ncbi.nlm.nih.gov/30016334/
# paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3683716/

GSE92538_geo_record = getGEO("GSE92538")
GSE92538_geo_record = GSE92538_geo_record[[1]]
GSE92538_pheno = pData(GSE92538_geo_record)
GSE92538_pheno = fix_columns(GSE92538_pheno)

# GSE92538_pheno is incomplete

# importing correct phenotypes
GSE92538_pheno_2 = geo_data_table_gsm_extract_phenotypes("GSE92538")
GSE92538_pheno_2 = fix_columns(GSE92538_pheno_2)
GSE92538_pheno_2$GENDER = ifelse(GSE92538_pheno_2$GENDER == "M", "Male", "Female")
GSE92538_pheno_2$GENDER = factor(GSE92538_pheno_2$GENDER, levels = c("Female", "Male"))
GSE92538_pheno_2$AGE = as.numeric(GSE92538_pheno_2$AGE)
GSE92538_pheno_2$POST_MORTEM_INTERVAL = as.numeric(GSE92538_pheno_2$POST_MORTEM_INTERVAL)
GSE92538_pheno_2$QC_BATCH = factor(GSE92538_pheno_2$QC_BATCH)
GSE92538_pheno_2$RACE = factor(GSE92538_pheno_2$RACE, levels = c("Caucasian", "Asian", "Other", "Hispanic", "African American"))
GSE92538_pheno_2$SUICIDE__1_YES_ = ifelse(GSE92538_pheno_2$SUICIDE__1_YES_ == "0", "Control", "Suicide")
GSE92538_pheno_2$SUICIDE__1_YES_ = factor(GSE92538_pheno_2$SUICIDE__1_YES_, levels = c("Control", "Suicide"))
GSE92538_pheno_2$TISSUE_PH__CEREBELLUM_ = as.numeric(GSE92538_pheno_2$TISSUE_PH__CEREBELLUM_) # produced 4 NA

table(GSE92538_pheno_2$PLATFORM_ID)
# GPL10526 HG-U133 Plus 2 128 samples
# GPL17027 U133A Array  235 samples

GSE92538_pheno_2_PLUS2 = GSE92538_pheno_2[GSE92538_pheno_2$PLATFORM_ID == "GPL10526", ]
GSE92538_pheno_2_PLUS2
table(GSE92538_pheno_2_PLUS2$SUICIDE__1_YES_)
"Control Suicide 
     95      33 "

length(unique(GSE92538_pheno_2_PLUS2$SUBJECT_ID)) # 75

GSE92538_pheno_2_U133A  = GSE92538_pheno_2[GSE92538_pheno_2$PLATFORM_ID == "GPL17027", ]
GSE92538_pheno_2_U133A
table(GSE92538_pheno_2_U133A$SUICIDE__1_YES_)
"Control Suicide 
     192      43 "
length(unique(GSE92538_pheno_2_U133A$SUBJECT_ID)) # 105 -> Higher power!
table(GSE92538_pheno_2_U133A$RACE)

# double array sample inspection
intersect(unique(GSE92538_pheno_2_U133A$SUBJECT_ID), unique(GSE92538_pheno_2_PLUS2$SUBJECT_ID))
# [1] "101078" "102118" "103311" "104126" "104259" "106229" "106589" "109779"
length(intersect(unique(GSE92538_pheno_2_U133A$SUBJECT_ID), unique(GSE92538_pheno_2_PLUS2$SUBJECT_ID))) # 8

GSE92538_dual_array_samples = intersect(unique(GSE92538_pheno_2_U133A$SUBJECT_ID), unique(GSE92538_pheno_2_PLUS2$SUBJECT_ID))

# There are multiple samples analyzed in the several labs
# Before proceeding to analysis we need to select unique samples after normalization
# There are 8 samples analyzed on both platforms -> should be removed from the larger sample

GSE92538_pheno_curated = GSE92538_pheno_2_U133A

GSE92538_replicated_samples = table(GSE92538_pheno_curated$SUBJECT_ID)
GSE92538_replicated_samples = as.data.frame(GSE92538_replicated_samples)
GSE92538_replicated_samples = GSE92538_replicated_samples[GSE92538_replicated_samples$Freq > 1, ]
GSE92538_replicated_samples = GSE92538_replicated_samples$Var1
GSE92538_replicated_samples = as.character(GSE92538_replicated_samples)


# reading the same annotation as in the previous u133a cohorts
probes_affy_hg_u133a_2 = read.csv("probes_affy_hg_u133a_2.csv")
probes_affy_hg_u133a_2$X = NULL
rownames(probes_affy_hg_u133a_2) = probes_affy_hg_u133a_2$ID
GSE92538_probes = probes_affy_hg_u133a_2


# Getting expression data
# Getting supplementary files
getGEOSuppFiles("GSE92538")
files = list.files("GSE92538")
files = files[stri_detect_fixed(files, pattern = ".tar")]
files = paste0("GSE92538", "/", files)
lapply(files, function(x) untar(x, exdir = "GSE92538_CEL"))

# Inspecting and preparing CEL files
files = list.files("GSE92538_CEL")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE92538_CEL", "/", files)
files = files[stri_detect_fixed(files, pattern = "GSM")] # 363 files -> original phenotypes are incorrect
files = files[stri_detect_fixed(files, pattern = "133A")] # 235 files from U133A

# All files are fixed
GSE92538_rawdata = affy::ReadAffy(filenames = files)

dir.create("GSE92538_U133A_results")
png(filename = "GSE92538_U133A_results/box_raw.png", width = 1920, height = 1080, units = "px")
affy::boxplot(GSE92538_rawdata)
dev.off()
png(filename = "GSE92538_U133A_results/dens_raw.png", width = 1920, height = 1080, units = "px")
affy::hist(GSE92538_rawdata)
dev.off()


GSE92538_rawdata_exprs = exprs(GSE92538_rawdata)
GSE92538_rawdata_meds = apply(GSE92538_rawdata_exprs, 2, median)
GSE92538_rawdata_meds = log2(GSE92538_rawdata_meds)
# it seems like values between 6.6 and 10 should be used...

samples_too_low = GSE92538_rawdata_meds[GSE92538_rawdata_meds < 6.6]
samples_too_low = names(samples_too_low)
samples_too_low = unlist(stri_split_fixed(samples_too_low, pattern = "_"))
samples_too_low = samples_too_low[stri_detect_fixed(samples_too_low, pattern = "GSM")]

df_too_low = GSE92538_pheno_curated[GSE92538_pheno_curated$GEO_ACCESSION %in% samples_too_low, ]
table(df_too_low$QC_BATCH) # mostly batch 8, some are from 9 and 1
table(df_too_low$SITE_OF_PROCESSING) # UC_Irvine: 18, UC_Davis: 1, U_Michigan: 4
table(df_too_low$SUBJECT_ID) # no duplicates

all(df_too_low$SUBJECT_ID %in% GSE92538_replicated_samples) # TRUE -> we can discard them!


samples_too_high = GSE92538_rawdata_meds[GSE92538_rawdata_meds > 10 ]
samples_too_high = names(samples_too_high)
samples_too_high = unlist(stri_split_fixed(samples_too_high, pattern = "_"))
samples_too_high = samples_too_high[stri_detect_fixed(samples_too_high, pattern = "GSM")]

df_too_high = GSE92538_pheno_curated[GSE92538_pheno_curated$GEO_ACCESSION %in% samples_too_high, ]
table(df_too_high$QC_BATCH) # batch 11
table(df_too_high$SITE_OF_PROCESSING) # UU_Michigan
df_too_high$SUBJECT_ID %in% GSE92538_replicated_samples # TRUE -> we can discard it!


GSE92538_bad_samples = c(df_too_low$GEO_ACCESSION, df_too_high$GEO_ACCESSION)
GSE92538_pheno_curated_2 = GSE92538_pheno_curated[GSE92538_pheno_curated$GEO_ACCESSION %!in% GSE92538_bad_samples, ]
table(GSE92538_pheno_curated_2$QC_BATCH)


files = files[stri_detect_fixed(files, pattern = "133A")] # 235 files from U133A

# All files are fixed
fixed_files = files
fixed_files_idex = sapply(fixed_files, function(x){
  
  x = unlist(stri_split_fixed(x, pattern = "/"))
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[stri_detect_fixed(x, pattern = "GSM")]
  
  if (x %in% GSE92538_pheno_curated_2$GEO_ACCESSION){
    return(TRUE)
  } else {
    return(FALSE)
  }
  
})
fixed_files = fixed_files[fixed_files_idex]

GSE92538_rawdata = affy::ReadAffy(filenames = fixed_files)

png(filename = "GSE92538_U133A_results/box_raw_selected.png", width = 1920, height = 1080, units = "px")
affy::boxplot(GSE92538_rawdata)
dev.off()
png(filename = "GSE92538_U133A_results/dens_raw_selected.png", width = 1920, height = 1080, units = "px")
affy::hist(GSE92538_rawdata)
dev.off()

GSE92538_expression = affy::rma(GSE92538_rawdata)

png(filename = "GSE92538_U133A_results/box_normal.png", width = 1920, height = 1080, units = "px")
affy::boxplot(GSE92538_expression)
dev.off()

png(filename = "GSE92538_U133A_results/dens_normal.png", width = 1920, height = 1080, units = "px")
affy::hist(GSE92538_expression)
dev.off()

rownames(GSE92538_expression)
# everything is normalized and seems OK!

colnames(GSE92538_expression) = sapply(colnames(GSE92538_expression), function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[stri_detect_fixed(x, pattern = "GSM")]
  return(x)
})
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # FALSE
all(rownames(GSE92538_expression) %in% GSE92538_probes$ID) # TRUE
all(GSE92538_pheno_curated_2$GEO_ACCESSION %in% colnames(GSE92538_expression)) # TRUE

# matching order
GSE92538_expression = GSE92538_expression[GSE92538_probes$ID,]
GSE92538_expression = GSE92538_expression[,GSE92538_pheno_curated_2$GEO_ACCESSION]
dim(GSE92538_expression) # 211 samples, 22283 probes

# check
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # TRUE
all(colnames(GSE92538_expression) == GSE92538_pheno_curated_2$GEO_ACCESSION) # TRUE


# removing duplicated measurements!
duplicated(GSE92538_pheno_curated_2$SUBJECT_ID)
table(duplicated(GSE92538_pheno_curated_2$SUBJECT_ID)) # 106 duplicates
table(GSE92538_pheno_2_PLUS2$SUBJECT_ID %in% GSE92538_pheno_curated_2$SUBJECT_ID) # 113 participants are missing
table(GSE92538_pheno_curated_2$QC_BATCH)


GSE92538_pheno_curated_3 = GSE92538_pheno_curated_2
GSE92538_pheno_curated_3 = dplyr::arrange(GSE92538_pheno_curated_3, SUBJECT_ID)
duplicated(GSE92538_pheno_curated_3$SUBJECT_ID)

GSE92538_pheno_curated_3 = distinct(GSE92538_pheno_curated_3, SUBJECT_ID, .keep_all = TRUE)
table(GSE92538_pheno_curated_3$QC_BATCH)
GSE92538_expression = GSE92538_expression[,GSE92538_pheno_curated_3$GEO_ACCESSION]
dim(GSE92538_expression) # 105 samples, 22283 probes

# check
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # TRUE
all(colnames(GSE92538_expression) == GSE92538_pheno_curated_3$GEO_ACCESSION) # TRUE

# Removing dual-array samples
GSE92538_pheno_curated_3 = GSE92538_pheno_curated_3[GSE92538_pheno_curated_3$SUBJECT_ID %!in% GSE92538_dual_array_samples,]
GSE92538_expression = GSE92538_expression[,GSE92538_pheno_curated_3$GEO_ACCESSION]

# check
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # TRUE
all(colnames(GSE92538_expression) == GSE92538_pheno_curated_3$GEO_ACCESSION) # TRUE
dim(GSE92538_expression)
"Features  Samples 
   22283       97 "
dim(GSE92538_pheno_curated_3)
"[1] 97 43"

# Differential expression (no covar)
Design.matrix = model.matrix(~ SUICIDE__1_YES_, data = GSE92538_pheno_curated_3)
fit = lmFit(GSE92538_expression, Design.matrix)
fitE = eBayes(fit)
GSE92538_U133A_Top_table_no_covar = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE92538_U133A_Top_table_no_covar$ID = rownames(GSE92538_U133A_Top_table_no_covar)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE92538_U133A_Top_table_no_covar$ID]
GSE92538_U133A_Top_table_no_covar$SE = SE

# Annotating results
GSE92538_probes = GSE92538_probes[GSE92538_U133A_Top_table_no_covar$ID, ]
GSE92538_U133A_Top_table_no_covar$Gene_symbol = GSE92538_probes$Gene_symbol
GSE92538_U133A_Top_table_no_covar$Tissue = GSE92538_pheno_curated_3$TISSUE[1]
GSE92538_U133A_Top_table_no_covar$Tissue_type = "Brain"
GSE92538_U133A_Top_table_no_covar$Technology = "Array"



# Differential expression (with covar)
GSE92538_pheno_curated_3$QC_BATCH = drop.levels(GSE92538_pheno_curated_3$QC_BATCH, reorder = FALSE)
table(GSE92538_pheno_curated_3$AGONAL_FACTOR, GSE92538_pheno_curated_3$SUICIDE__1_YES_)
" Control Suicide
  0      35      14
  1      33       2
  2       6       1
  3       5       0"
GSE92538_pheno_curated_3$AGONAL_FACTOR = ifelse(GSE92538_pheno_curated_3$AGONAL_FACTOR=="NA", NA, GSE92538_pheno_curated_3$AGONAL_FACTOR)

GSE92538_pheno_curated_3$AGONAL_FACTOR_binary = ifelse(GSE92538_pheno_curated_3$AGONAL_FACTOR=="0", "AFS=0", "AFS>=1")
GSE92538_pheno_curated_3$AGONAL_FACTOR_binary = factor(GSE92538_pheno_curated_3$AGONAL_FACTOR_binary, levels = c("AFS=0", "AFS>=1"))

table(GSE92538_pheno_curated_3$AGONAL_FACTOR_binary, GSE92538_pheno_curated_3$SUICIDE__1_YES_)

"         Control Suicide
  AFS=0       35      14
  AFS>=1      45       3
"

table(GSE92538_pheno_curated_3$DIAGNOSIS, GSE92538_pheno_curated_3$SUICIDE__1_YES_)
"                              Control Suicide
  Bipolar Disorder               12       7
  Control                        51       0
  Major Depressive Disorder      13       9
  Schizophrenia                   4       1
"
GSE92538_pheno_curated_3$DIAGNOSIS = factor(GSE92538_pheno_curated_3$DIAGNOSIS, levels = c("Control", 
                                                                                           "Major Depressive Disorder",
                                                                                           "Bipolar Disorder", 
                                                                                           "Schizophrenia"))

# 13 participants are excluded due to missing tissue pH or AGONAL_FACTOR_binary
Design.matrix = model.matrix(~ SUICIDE__1_YES_ + GENDER + AGE + DIAGNOSIS + AGONAL_FACTOR_binary + POST_MORTEM_INTERVAL + QC_BATCH + RACE + TISSUE_PH__CEREBELLUM_, data = GSE92538_pheno_curated_3)
GSE92538_pheno_curated_3_TMP = GSE92538_pheno_curated_3[GSE92538_pheno_curated_3$GEO_ACCESSION %in% row.names(Design.matrix), ]
table(GSE92538_pheno_curated_3_TMP$SUICIDE__1_YES_)
"
Control Suicide 
     67      17
"
fit = lmFit(GSE92538_expression[,rownames(Design.matrix)], Design.matrix)
fitE = eBayes(fit)
GSE92538_U133A_Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE92538_U133A_Top_table$ID = rownames(GSE92538_U133A_Top_table)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE92538_U133A_Top_table$ID]
GSE92538_U133A_Top_table$SE = SE

# Annotating results
GSE92538_probes = GSE92538_probes[GSE92538_U133A_Top_table$ID, ]
GSE92538_U133A_Top_table$Gene_symbol = GSE92538_probes$Gene_symbol
GSE92538_U133A_Top_table$Tissue = GSE92538_pheno_curated_3$TISSUE[1]
GSE92538_U133A_Top_table$Tissue_type = "Brain"
GSE92538_U133A_Top_table$Technology = "Array"

write.csv(GSE92538_expression, "GSE92538_U133A_results/GSE92538_expression.csv")
write.csv(GSE92538_probes, "GSE92538_U133A_results/GSE92538_probes.csv")
write.csv(GSE92538_pheno_curated_3, "GSE92538_U133A_results/GSE92538_pheno_curated.csv")
write.csv(GSE92538_U133A_Top_table_no_covar, "GSE92538_U133A_results/GSE92538_U133A_Top_table_no_covar.csv")
write.csv(GSE92538_U133A_Top_table, "GSE92538_U133A_results/GSE92538_U133A_Top_table.csv")

rm(list = ls(pattern = "GSE92538"))
gc()


################### GSE92538 PLUS2 ###################

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92538

# paper: https://pubmed.ncbi.nlm.nih.gov/30016334/


GSE92538_geo_record = getGEO("GSE92538")
GSE92538_geo_record = GSE92538_geo_record[[1]]
GSE92538_pheno = pData(GSE92538_geo_record)
GSE92538_pheno = fix_columns(GSE92538_pheno)

# GSE92538_pheno is incomplete

# importing correct phenotypes
GSE92538_pheno_2 = geo_data_table_gsm_extract_phenotypes("GSE92538")
GSE92538_pheno_2 = fix_columns(GSE92538_pheno_2)
GSE92538_pheno_2$GENDER = ifelse(GSE92538_pheno_2$GENDER == "M", "Male", "Female")
GSE92538_pheno_2$GENDER = factor(GSE92538_pheno_2$GENDER, levels = c("Female", "Male"))
GSE92538_pheno_2$AGE = as.numeric(GSE92538_pheno_2$AGE)
GSE92538_pheno_2$POST_MORTEM_INTERVAL = as.numeric(GSE92538_pheno_2$POST_MORTEM_INTERVAL)
GSE92538_pheno_2$QC_BATCH = factor(GSE92538_pheno_2$QC_BATCH)
GSE92538_pheno_2$RACE = factor(GSE92538_pheno_2$RACE, levels = c("Caucasian", "Asian", "Other", "Hispanic", "African American"))
GSE92538_pheno_2$SUICIDE__1_YES_ = ifelse(GSE92538_pheno_2$SUICIDE__1_YES_ == "0", "Control", "Suicide")
GSE92538_pheno_2$SUICIDE__1_YES_ = factor(GSE92538_pheno_2$SUICIDE__1_YES_, levels = c("Control", "Suicide"))
GSE92538_pheno_2$TISSUE_PH__CEREBELLUM_ = as.numeric(GSE92538_pheno_2$TISSUE_PH__CEREBELLUM_) # produced 4 NA

table(GSE92538_pheno_2$PLATFORM_ID)
# GPL10526 HG-U133 Plus 2 128 samples
# GPL17027 U133A Array  235 samples

GSE92538_pheno_2_PLUS2 = GSE92538_pheno_2[GSE92538_pheno_2$PLATFORM_ID == "GPL10526", ]
GSE92538_pheno_2_PLUS2
table(GSE92538_pheno_2_PLUS2$SUICIDE__1_YES_)
"Control Suicide 
     95      33 "

length(unique(GSE92538_pheno_2_PLUS2$SUBJECT_ID)) # 75


# There are multiple samples analyzed in the several labs
# Before proceeding to analysis we need to select unique samples after normalization

GSE92538_pheno_curated = GSE92538_pheno_2_PLUS2

GSE92538_replicated_samples = table(GSE92538_pheno_curated$SUBJECT_ID)
GSE92538_replicated_samples = as.data.frame(GSE92538_replicated_samples)
GSE92538_replicated_samples = GSE92538_replicated_samples[GSE92538_replicated_samples$Freq > 1, ]
GSE92538_replicated_samples = GSE92538_replicated_samples$Var1
GSE92538_replicated_samples = as.character(GSE92538_replicated_samples)


# getting annotation

probes_affy_hg_u133_plus2 = smart_fread("HG_U133_Plus_2.txt", skip=16)
affy_hg_u133_plus2_genes = unlist(stri_split_fixed(probes_affy_hg_u133_plus2$`Gene Symbol`, pattern = " /// "))
affy_hg_u133_plus2_genes = unique(affy_hg_u133_plus2_genes)
affy_hg_u133_plus2_genes = affy_hg_u133_plus2_genes[affy_hg_u133_plus2_genes != ""]
affy_hg_u133_plus2_genes = affy_hg_u133_plus2_genes[!is.na(affy_hg_u133_plus2_genes)]
affy_hg_u133_plus2_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = affy_hg_u133_plus2_genes, 
                                                  PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                  PRF_replace_NA_with_old = TRUE)
affy_hg_u133_plus2_genes_check_dict = affy_hg_u133_plus2_genes_check$Suggested.Symbol
names(affy_hg_u133_plus2_genes_check_dict) = affy_hg_u133_plus2_genes_check$x

probes_affy_hg_u133_plus2$Gene_symbol = sapply(probes_affy_hg_u133_plus2$`Gene Symbol`, function(x){
  x = unlist(stri_split_fixed(x, pattern = " /// "))
  x = unique(x)
  x = x[x != ""]
  x = x[!is.na(x)]
  
  if (length(x) < 1){
    return(NA)
  }
  x = sapply(x, function(z) affy_hg_u133_plus2_genes_check_dict[z])
  
  x = paste0(x, collapse = ";")
  return(x)
})

write.csv(probes_affy_hg_u133_plus2, "probes_affy_hg_u133_plus2.csv")

GSE92538_probes  = probes_affy_hg_u133_plus2
rownames(GSE92538_probes) = GSE92538_probes$ID

# Inspecting and preparing CEL files
files = list.files("GSE92538_CEL")
files = files[stri_detect_fixed(files, pattern = ".gz")]
files = paste0("GSE92538_CEL", "/", files)
files = files[stri_detect_fixed(files, pattern = "GSM")] # 363 files -> original phenotypes are incorrect
files = files[stri_detect_fixed(files, pattern = "133P")] # 128 files from U133 PLUS2

# All files are fixed
GSE92538_rawdata = affy::ReadAffy(filenames = files)

dir.create("GSE92538_U133_PLUS2_results")
png(filename = "GSE92538_U133_PLUS2_results/box_raw.png", width = 1920, height = 1080, units = "px")
affy::boxplot(GSE92538_rawdata)
dev.off()
png(filename = "GSE92538_U133_PLUS2_results/dens_raw.png", width = 1920, height = 1080, units = "px")
affy::hist(GSE92538_rawdata)
dev.off()
# quite messy but seems like acceptable...

GSE92538_expression = affy::rma(GSE92538_rawdata)

png(filename = "GSE92538_U133_PLUS2_results/box_normal.png", width = 1920, height = 1080, units = "px")
affy::boxplot(GSE92538_expression)
dev.off()

png(filename = "GSE92538_U133_PLUS2_results/dens_normal.png", width = 1920, height = 1080, units = "px")
affy::hist(GSE92538_expression)
dev.off()
# normalization seems good

colnames(GSE92538_expression) = sapply(colnames(GSE92538_expression), function(x){
  x = unlist(stri_split_fixed(x, pattern = "_"))
  x = x[stri_detect_fixed(x, pattern = "GSM")]
  return(x)
})
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # FALSE
all(rownames(GSE92538_expression) %in% GSE92538_probes$ID) # TRUE
all(GSE92538_pheno_curated$GEO_ACCESSION %in% colnames(GSE92538_expression)) # TRUE

# matching order
GSE92538_expression = GSE92538_expression[GSE92538_probes$ID,]
GSE92538_expression = GSE92538_expression[,GSE92538_pheno_curated$GEO_ACCESSION]
dim(GSE92538_expression) # 128 samples, 54675 probes

# check
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # TRUE
all(colnames(GSE92538_expression) == GSE92538_pheno_curated$GEO_ACCESSION) # TRUE

# removing duplicated measurements!
duplicated(GSE92538_pheno_curated$SUBJECT_ID)
table(duplicated(GSE92538_pheno_curated$SUBJECT_ID)) # 53 duplicates
table(GSE92538_pheno_curated$QC_BATCH)
table(is.na(GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_)) # 4 missing

GSE92538_pheno_curated_3 = GSE92538_pheno_curated
GSE92538_pheno_curated_3 = dplyr::arrange(GSE92538_pheno_curated_3, SUBJECT_ID)
duplicated(GSE92538_pheno_curated_3$SUBJECT_ID)

GSE92538_pheno_curated_3 = distinct(GSE92538_pheno_curated_3, SUBJECT_ID, .keep_all = TRUE)
table(GSE92538_pheno_curated_3$QC_BATCH)
GSE92538_expression = GSE92538_expression[,GSE92538_pheno_curated_3$GEO_ACCESSION]
dim(GSE92538_expression) # 75 samples, 54675 probes
table(is.na(GSE92538_pheno_curated_3$TISSUE_PH__CEREBELLUM_)) # 2 missing

# check
all(rownames(GSE92538_expression) == GSE92538_probes$ID) # TRUE
all(colnames(GSE92538_expression) == GSE92538_pheno_curated_3$GEO_ACCESSION) # TRUE


# Differential expression (no covar)
Design.matrix = model.matrix(~ SUICIDE__1_YES_, data = GSE92538_pheno_curated_3)
fit = lmFit(GSE92538_expression, Design.matrix)
fitE = eBayes(fit)
GSE92538_U133_PLUS2_Top_table_no_covar = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE92538_U133_PLUS2_Top_table_no_covar$ID = rownames(GSE92538_U133_PLUS2_Top_table_no_covar)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE92538_U133_PLUS2_Top_table_no_covar$ID]
GSE92538_U133_PLUS2_Top_table_no_covar$SE = SE

# Annotating results
GSE92538_probes = GSE92538_probes[GSE92538_U133_PLUS2_Top_table_no_covar$ID, ]
GSE92538_U133_PLUS2_Top_table_no_covar$Gene_symbol = GSE92538_probes$Gene_symbol
GSE92538_U133_PLUS2_Top_table_no_covar$Tissue = GSE92538_pheno_curated_3$TISSUE[1]
GSE92538_U133_PLUS2_Top_table_no_covar$Tissue_type = "Brain"
GSE92538_U133_PLUS2_Top_table_no_covar$Technology = "Array"


# Differential expression (with covar)
GSE92538_pheno_curated_3$QC_BATCH = drop.levels(GSE92538_pheno_curated_3$QC_BATCH, reorder = FALSE)
GSE92538_pheno_curated_3$RACE = drop.levels(GSE92538_pheno_curated_3$RACE, reorder = FALSE)

table(GSE92538_pheno_curated_3$AGONAL_FACTOR, GSE92538_pheno_curated_3$SUICIDE__1_YES_)
"   Control Suicide
  0      54      20
  2       1       0
"
GSE92538_pheno_curated_3$AGONAL_FACTOR_binary = ifelse(GSE92538_pheno_curated_3$AGONAL_FACTOR=="0", "AFS=0", "AFS>=1")
GSE92538_pheno_curated_3$AGONAL_FACTOR_binary = factor(GSE92538_pheno_curated_3$AGONAL_FACTOR_binary, levels = c("AFS=0", "AFS>=1"))

table(GSE92538_pheno_curated_3$DIAGNOSIS, GSE92538_pheno_curated_3$SUICIDE__1_YES_)
" 
                            Control Suicide
  Bipolar Disorder                4       3
  Control                        32       0
  Major Depressive Disorder       6      13
  Schizophrenia                  13       4
"
GSE92538_pheno_curated_3$DIAGNOSIS = factor(GSE92538_pheno_curated_3$DIAGNOSIS, levels = c("Control", 
                                                                                           "Major Depressive Disorder",
                                                                                           "Bipolar Disorder", 
                                                                                           "Schizophrenia"))

table(GSE92538_pheno_curated_3$AGONAL_FACTOR_binary)

# only one participant has value aboe 0 in AGONAL_FACTOR_binary -> skipped
Design.matrix = model.matrix(~ SUICIDE__1_YES_ + GENDER + AGE + DIAGNOSIS + POST_MORTEM_INTERVAL + QC_BATCH + RACE + TISSUE_PH__CEREBELLUM_, data = GSE92538_pheno_curated_3)
GSE92538_pheno_curated_3_TMP = GSE92538_pheno_curated_3[GSE92538_pheno_curated_3$GEO_ACCESSION %in% row.names(Design.matrix), ]
table(GSE92538_pheno_curated_3_TMP$SUICIDE__1_YES_)
"
Control Suicide 
     54      19
"
# 2 participants are dropped ue to missing TISSUE_PH__CEREBELLUM_

fit = lmFit(GSE92538_expression[,rownames(Design.matrix)], Design.matrix)
fitE = eBayes(fit)
GSE92538_U133_PLUS2_Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE92538_U133_PLUS2_Top_table$ID = rownames(GSE92538_U133_PLUS2_Top_table)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE92538_U133_PLUS2_Top_table$ID]
GSE92538_U133_PLUS2_Top_table$SE = SE

# Annotating results
GSE92538_probes = GSE92538_probes[GSE92538_U133_PLUS2_Top_table$ID, ]
GSE92538_U133_PLUS2_Top_table$Gene_symbol = GSE92538_probes$Gene_symbol
GSE92538_U133_PLUS2_Top_table$Tissue = GSE92538_pheno_curated_3$TISSUE[1]
GSE92538_U133_PLUS2_Top_table$Tissue_type = "Brain"
GSE92538_U133_PLUS2_Top_table$Technology = "Array"


write.csv(GSE92538_expression, "GSE92538_U133_PLUS2_results/GSE92538_expression.csv")
write.csv(GSE92538_probes, "GSE92538_U133_PLUS2_results/GSE92538_probes.csv")
write.csv(GSE92538_pheno_curated_3, "GSE92538_U133_PLUS2_results/GSE92538_pheno_curated.csv")
write.csv(GSE92538_U133_PLUS2_Top_table_no_covar, "GSE92538_U133_PLUS2_results/GSE92538_U133_PLUS2_Top_table_no_covar.csv")
write.csv(GSE92538_U133_PLUS2_Top_table, "GSE92538_U133_PLUS2_results/GSE92538_U133_PLUS2_Top_table.csv")

rm(list = ls(pattern = "GSE92538"))
gc()



################### RNA Seq cohorts ###################

################### Preparing genome index for STAR ###################

# This file shows commands that were used to perform RNAseq aanalyses with STAR
"

STAR documentation is available here : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/ and here: https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html

and here: https://github.com/alexdobin/STAR

"

# some extra infor on paired-end data
"
# The data was just paired-end instead of single-end, so I had to fastq-dump -split-files the fastq files before I aligned them to the genome. After that, my mapped reads shot up to ~91% as well.
https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump

"

## Making genome index.

# We need to download genome file and genome annotation file.


# TERMINAL
"

conda activate tf-py38

conda install -c bioconda star

cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/genome_files_hg19

wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa.gz (genome file)
wget https://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz (annotation file)

gunzip *.gz

"

# making genome index with STAR
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/

mkdir genome_index_STAR_hg19

STAR --runThreadN 10 --runMode genomeGenerate --genomeDir genome_index_STAR_hg19/ --genomeFastaFiles genome_files_hg19/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# installing SRA-Toolkit and SRA tools

"
conda activate tf-py38
cd /home/aleksandr/Downloads
chmod +x setup-apt.sh
sudo sh setup-apt.sh
source /etc/profile.d/sra-tools.sh
"

# Example on SRR5961796 (paired-end seq)
"
mkdir GSE102556_fastq
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE102556_fastq
source /etc/profile.d/sra-tools.sh
fasterq-dump --split-files --threads 4 SRR5961796
"
# produced 2 files (13.5 GB each)

# --split-files and split-3 produce mostly the same results


# example with prefetch 
# MUCH FASTER

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE102556_pref
source /etc/profile.d/sra-tools.sh
prefetch SRR5961796 -O GSE102556_pref

"

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE102556_pref
source /etc/profile.d/sra-tools.sh
fasterq-dump --threads 9 SRR5961796

"

# Qality control with FastQC

"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR5961796_QC
FastQC/fastqc GSE102556_pref/SRR5961796_1.fastq -o SRR5961796_QC/

"


# Alignment to genome. 83.01% uniquely mapped, % of reads mapped to multiple loci |	10.78%; % of reads unmapped: too short |	6.01%


"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE102556_mapped

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE102556_pref/SRR5961796_1.fastq GSE102556_pref/SRR5961796_2.fastq \
--outFileNamePrefix GSE102556_mapped/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# Qality control with FastQC post-alignment (duolicates are gone!)

"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
FastQC/fastqc GSE102556_mapped/Aligned.sortedByCoord.out.bam -o SRR5961796_QC/

"

####
####
#### 
#### Alternatives 

# trimming with Trimmomatic (specified adapters) -> VERY SLOW!

"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE102556_trimmed
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE GSE102556_pref/SRR5961796_1.fastq GSE102556_pref/SRR5961796_2.fastq \
GSE102556_trimmed/output_forward_paired.fq.gz GSE102556_trimmed/output_forward_unpaired.fq.gz \
GSE102556_trimmed/output_reverse_paired.fq.gz GSE102556_trimmed/output_reverse_unpaired.fq.gz \
ILLUMINACLIP:Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36 -threads 9

"

# no specified adapter
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE102556_trimmed
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE GSE102556_pref/SRR5961796_1.fastq GSE102556_pref/SRR5961796_2.fastq \
GSE102556_trimmed/output_forward_paired.fq.gz GSE102556_trimmed/output_forward_unpaired.fq.gz \
GSE102556_trimmed/output_reverse_paired.fq.gz GSE102556_trimmed/output_reverse_unpaired.fq.gz \
LEADING:3 TRAILING:3 MINLEN:36

"

# trimming with fastp
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE102556_trimmed
conda activate trimming
fastp --thread 8 -i GSE102556_pref/SRR5961796_1.fastq -I GSE102556_pref/SRR5961796_2.fastq -o GSE102556_trimmed/out.R1.fq.gz -O GSE102556_trimmed/out.R2.fq.gz

"

# STAR after trimming -> only minor improvements
# skipping trimming is OK...

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE102556_mapped_trimmed

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE102556_trimmed/out.R1.fq.gz GSE102556_trimmed/out.R2.fq.gz --readFilesCommand gunzip -c \
--outFileNamePrefix GSE102556_mapped_trimmed/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# counting reads for paired-end
# installing: conda install -c bioconda subread

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE102556_counts

bin/featureCounts -T 10 -p --countReadPairs -a genome_files_hg19/Homo_sapiens.GRCh37.87.gtf \
-t exon -g gene_id -o GSE102556_counts/counts.txt GSE102556_mapped/Aligned.sortedByCoord.out.bam

"
# 52% for gene and 35% for exon


# STAR Alignment to genome and counting (does not match subread feature counts)
# This option almost matches exon counting with featureCounts


"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE102556_mapped_2

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE102556_fastq/SRR5961796_1.fastq GSE102556_fastq/SRR5961796_2.fastq \
--outFileNamePrefix GSE102556_mapped_2/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

################### Triiming notes ###################
# some notes on trimming
# https://dnatech.genomecenter.ucdavis.edu/faqs/when-should-i-trim-my-illumina-reads-and-how-should-i-do-it/#:~:text=For%20counting%20applications%20such%20as,pseudo%2Daligners%20should%20be%20used.
# https://www.biostars.org/p/474013/
# https://github.com/alexdobin/STAR/issues/242
# https://www.biostars.org/p/298619/#298642
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7671312/


################### GSE102556 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE102556

# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5734943/

"
Any reads that fell in multiple genes were excluded from the analysis. Threshold for filtering 
out low expressed genes was set to >5 reads in at least 80% of the samples as previously described123.
"

"
 In the first stage, linear models implemented in the ‘limma’ package124 of Bioconductor125 were used to 
 compute the variance of gene expression across all groups. Gene expression was transformed and normalized 
 using voom in the limma package. In humans, models were adjusted for age, RIN, alcohol abuse, and medication 
 status. 
 
"

"
TPM (transcripts per million)
FPKM (fragments per kilobase of transcript per million fragments mapped)

"


"Samples were barcoded for multiplexing and sequenced at 50 bp paired-end on Illumina HiSeq2500. 
Samples were pooled eight per lane and sequenced twice at a depth of 50 million reads per sample."

# FPKM values can't be used with voom: https://support.bioconductor.org/p/56275/

"
Dear Jon,

No, it is absolutely not ok to input FPKM values to voom. The results will be nonsense. 
Same goes for edgeR and DESeq, for the reasons explained in the documentation for those packages.

Note that a matrix of FPKM is not a matrix of counts.

It seems to me that you have three options in decreasing order of desirability:

1) Get the actual integer counts from which the FPKM were computed and do a proper analysis, for example using voom or edgeR.

2) Get the gene lengths and library sizes used to compute the FPKM and convert the FPKM back to counts.

3) If FPKM is really all you have, then convert the values to a log2 scale (y = log2(FPKM+0.1) say) and do an 
ordinary limma analysis as you would for microarray data, using eBayes() with trend=TRUE. Do not use voom,
do not use edgeR, do not use DESeq. (Do not pass go and do not collect $200.) This isn't 100% ideal, but is 
probably the best analysis available. You make this method somewhat better by using arrayWeights() as well which, 
in this context, will attempt to estimate the library sizes the FPKMs were originally computed from. Nevertheless,
the mean-variance trend estimated by limma from the logFPKMs will never be as smooth or as informative as the trend 
that would have been estimated had you had the real counts.

The third option is similar to the limma-trend analysis described in the limma preprint, except that it is applied to 
the logFPKM instead of logCPM. Statistically this will not perform as well as it would applied to the logCPM.

Best wishes
Gordon


"
## Makes sense to use STAR aligner based on many papers and tutorials
## Data has different provided formats -> we should obtain the same format for analysis to be valid...
GSE102556_geo_record = getGEO("GSE102556")
GSE102556_geo_record = GSE102556_geo_record[[1]]
GSE102556_pheno = pData(GSE102556_geo_record)
GSE102556_pheno = fix_columns(GSE102556_pheno)


GSE102556_pheno_curated = GSE102556_pheno
table(GSE102556_pheno_curated$ORGANISM_CH1) # 263 human

GSE102556_pheno_curated$AGE
GSE102556_pheno_curated$ALCOOL
GSE102556_pheno_curated$CAUSE_OF_DEATH
GSE102556_pheno_curated$DRUG_TYPE
GSE102556_pheno_curated$DRUGS
GSE102556_pheno_curated$GENDER
GSE102556_pheno_curated$MEDICATION
GSE102556_pheno_curated$MEDICATION_TYPE
GSE102556_pheno_curated$PH
GSE102556_pheno_curated$PHENOTYPE
GSE102556_pheno_curated$PMI
GSE102556_pheno_curated$RIN
GSE102556_pheno_curated$SMOKING


# quantification
table(GSE102556_pheno_curated$ALCOOL)

"NA  no  No yes 
116  97  23  27 "

table(GSE102556_pheno_curated$CAUSE_OF_DEATH)

"Accident  Natural  Suicide 
      21       39      203 "

table(GSE102556_pheno_curated$DRUG_TYPE)
"cannabis       NA       no 
       5      116      142 "


table(GSE102556_pheno_curated$DRUGS)
" NA  no  No yes 
116 119  23   5 "


table(GSE102556_pheno_curated$MEDICATION)

"NA  no yes 
 75  99  89"

table(GSE102556_pheno_curated$MEDICATION_TYPE)
" AD anticonvulsivant               AP          Lithium               NA               no 
  55                5               11               12               81               99 "

table(GSE102556_pheno_curated$PHENOTYPE)

"CTRL  MDD 
 122  141 "

table(GSE102556_pheno_curated$SMOKING)
"           NA            no   yes (heavy) yes(moderate) 
           77            79           101             6 "

table(GSE102556_pheno_curated$CAUSE_OF_DEATH, GSE102556_pheno_curated$PHENOTYPE)
"           CTRL MDD
  Accident   21   0
  Natural    39   0
  Suicide    62 141"
# suicide will be converted to binary

"In humans, models were adjusted for age, RIN, alcohol abuse, and medication status. Among the extensive 
information collected for the individuals in this study, these variables were selected for adjustment based 
on a combination of domain knowledge and variance analysis of the RNAseq data. Eigen-R2126 was used to estimate 
the amount of variance in RNAseq data explained by each variable. The estimate for each variable is similar to taking 
the average of the correlations between the variable and the expression values for each gene. Correlation averages are 
vulnerable to technical artifacts such as stochastic noise for genes with little or no expression values so Eigen-R2 uses 
principal component analysis to reduce the contribution of these and other problematic genes. Using Eigen-R2 we found evidence 
that phenotype (PC4: p<0.01), sex (PC3: p<1.0e-4) and brain regions (PC1: p<5.0e-36) explained the most variance of any variable. 
Differential gene expression was assessed through a generalized linear model implemented in limma, with phenotype (MDD vs CTRL) 
and sex (male and female) as main factors for every brain region. Our analyses were adjusted for age (PC2: p<0.01), RIN (PC1: 
p<1.0e-3) and alcohol abuse (PC2; p<0.05) because these contributed to the variance according to the Eigen-R2 analyses, and we 
chose to adjust for medication status because medication is well known to change gene expression in the brain75."

# They only used phenotype, sex, age, RIN, alcohol, medication status
# I would also add PH and PMI


# curation
GSE102556_pheno_curated$ALCOOL = toupper(GSE102556_pheno_curated$ALCOOL)
GSE102556_pheno_curated$CAUSE_OF_DEATH = toupper(GSE102556_pheno_curated$CAUSE_OF_DEATH)
GSE102556_pheno_curated$DRUG_TYPE  = toupper(GSE102556_pheno_curated$DRUG_TYPE)
GSE102556_pheno_curated$DRUGS  = toupper(GSE102556_pheno_curated$DRUGS)
GSE102556_pheno_curated$GENDER  = toupper(GSE102556_pheno_curated$GENDER)
GSE102556_pheno_curated$MEDICATION  = toupper(GSE102556_pheno_curated$MEDICATION)
GSE102556_pheno_curated$MEDICATION_TYPE  = toupper(GSE102556_pheno_curated$MEDICATION_TYPE)
GSE102556_pheno_curated$PHENOTYPE  = toupper(GSE102556_pheno_curated$PHENOTYPE)
GSE102556_pheno_curated$SMOKING  = toupper(GSE102556_pheno_curated$SMOKING)

GSE102556_pheno_curated$RELATION_1
GSE102556_pheno_curated$SRA_file = stri_replace_all_fixed(GSE102556_pheno_curated$RELATION_1, 
                                                          pattern = "SRA: https://www.ncbi.nlm.nih.gov/sra?term=", 
                                                          replacement = "")

GSE102556_SRA_run_table = smart_fread("GSE102556/SraRunTable.txt")
all(GSE102556_pheno_curated$SRA_file %in% GSE102556_SRA_run_table$Experiment) # TRUE
GSE102556_SRA_run_table = GSE102556_SRA_run_table[GSE102556_SRA_run_table$Experiment %in% GSE102556_pheno_curated$SRA_file, ]
GSE102556_SRA_run_table = GSE102556_SRA_run_table[GSE102556_SRA_run_table$`GEO_Accession (exp)` %in% GSE102556_pheno_curated$GEO_ACCESSION,]

any(duplicated(GSE102556_SRA_run_table$Experiment)) # TRUE
any(duplicated(GSE102556_SRA_run_table$`GEO_Accession (exp)`)) # TRUE

GSE102556_dupl_SRA = GSE102556_SRA_run_table[duplicated(GSE102556_SRA_run_table$Experiment), "Experiment"]
GSE102556_SRA_run_table[GSE102556_SRA_run_table$Experiment %in% GSE102556_dupl_SRA, ] # complete duplicates but runs (SRR) are different as well as file types
GSE102556_SRA_run_table = GSE102556_SRA_run_table[stri_detect_fixed(GSE102556_SRA_run_table$`DATASTORE filetype`, pattern = "fastq"), ]

# check 
all(GSE102556_SRA_run_table$`GEO_Accession (exp)` %in% GSE102556_pheno_curated$GEO_ACCESSION) # TRUE
all(GSE102556_pheno_curated$GEO_ACCESSION %in% GSE102556_SRA_run_table$`GEO_Accession (exp)`) # TRUE

GSE102556_SRA_run_table$Run

# Test on SRR5961796
# write run table to a file test_id.txt
GSE102556_SRA_run_table$Run[1] #SRR5961796
writeLines(GSE102556_SRA_run_table$Run[1:3], "test_id.txt")

any(duplicated(GSE102556_SRA_run_table$Run)) # FALSE
writeLines(GSE102556_SRA_run_table$Run, "GSE102556_SRR_ids.txt")
GSE102556_pheno_curated$RUN_ID = sapply(GSE102556_pheno_curated$GEO_ACCESSION, function(x){
  x = GSE102556_SRA_run_table[GSE102556_SRA_run_table$`GEO_Accession (exp)` == x, "Run"]
  return(x)
})


# Obtaining counts for ALL runs (Terminal). It takes approx. 10 minutes per sample
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTAR.py GSE102556_SRR_ids.txt PAIRED GSE102556_

"

# Inspecting runs
GSE102556_run_logs = list.files("GSE102556_OUTPUT", pattern = "logs", full.names = TRUE)
GSE102556_mapped_percent = sapply(GSE102556_run_logs, function(x){
  x = readLines(x)
  x = x[10]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads % |\t", replacement = "")
  x = stri_replace_all_fixed(str = x, pattern = "%", replacement = "")
  return(x)
})
GSE102556_mapped_percent = as.numeric(GSE102556_mapped_percent)
min(GSE102556_mapped_percent) # 26.7
max(GSE102556_mapped_percent) # 89.1
GSE102556_mapped_percent_df = data.frame(file = GSE102556_run_logs, percent = GSE102556_mapped_percent)
# GSE102556_OUTPUT/SRR5961809_star_logs.txt 66 % of reads are too short -> exclude ?


# preparing counts data
GSE102556_runs = list.files("GSE102556_OUTPUT", pattern = "counts", full.names = TRUE)
GSE102556_count_df = lapply(GSE102556_runs, function(x){
  x = smart_fread(x)
  return(x)
})
names(GSE102556_count_df) = GSE102556_runs
reference = GSE102556_count_df[[1]]$gene_ID
GSE102556_compar_to_ref = sapply(GSE102556_count_df, function(x){
  compar = reference == x$gene_ID
  compar = all(compar)
  return(compar)
})
all(GSE102556_compar_to_ref) # all row names match


### checking strandness (package developer)
"Hi Wendy,

I believe for Illumina Tru-seq stranded protocol you need to use the 4th column.
However, the general rule is to compare the total counts over genes in the 3rd and 4th column, and select the column which has much larger total counts.
If the counts in two columns are not very different, the protocol is unstranded and you need to use the 2nd column.
Yet even an easier method is to look at the N_noFeature line for the 3rd and 4th column and pick the column with the lowest count.

Cheers
Alex"

GSE102556_participant_1 = GSE102556_count_df[[1]]
GSE102556_participant_1 = GSE102556_participant_1[-(1:4), ]
GSE102556_participant_1$gene_ID = NULL
GSE102556_participant_1 = GSE102556_participant_1[rowSums(GSE102556_participant_1) > 0, ]
table(GSE102556_participant_1$Strand_1 > GSE102556_participant_1$Strand_2)
table(GSE102556_participant_1$Strand_1 < GSE102556_participant_1$Strand_2)
# Strand 1 is systematically bigger than 2
sum(GSE102556_participant_1$Strand_1) # 24927202
sum(GSE102556_participant_1$Strand_2) # 1235431
sum(GSE102556_participant_1$Strand_1)/sum(GSE102556_participant_1$Strand_2) # 20.17693 20x difference

GSE102556_count = lapply(GSE102556_count_df, function(x){
  x = x[,"Strand_1"]
  return(x)
})
GSE102556_count = do.call(cbind,GSE102556_count)
rownames(GSE102556_count) = reference
colnames(GSE102556_count) = sapply(colnames(GSE102556_count), function(x){
  x = stri_replace_all_fixed(x, pattern = "GSE102556_OUTPUT/", replacement = "")
  x = stri_replace_all_fixed(x, pattern = "_star_counts.csv", replacement = "")
  return(x)
})
GSE102556_count = as.data.frame(GSE102556_count)

# check
all(colnames(GSE102556_count) == GSE102556_pheno_curated$RUN_ID) # FALSE
all(colnames(GSE102556_count) %in% GSE102556_pheno_curated$RUN_ID) # TRUE

GSE102556_count = GSE102556_count[,GSE102556_pheno_curated$RUN_ID]
all(colnames(GSE102556_count) == GSE102556_count$RUN_ID) # TRUE

# filtering bad probes
GSE102556_pheno_curated_2 = GSE102556_pheno_curated[GSE102556_pheno_curated$RUN_ID != "SRR5961809",]
GSE102556_count = GSE102556_count[,GSE102556_pheno_curated_2$RUN_ID]

# Pheno curation
# Initial paper used phenotype, sex, age, RIN, alcohol, medication status
# I would also add PH and PMI
# smoking and drugs were not used due to high missing counts
GSE102556_pheno_curated_2$AGE
GSE102556_pheno_curated_2$ALCOOL
GSE102556_pheno_curated_2$CAUSE_OF_DEATH
GSE102556_pheno_curated_2$GENDER
GSE102556_pheno_curated_2$MEDICATION
GSE102556_pheno_curated_2$PH
GSE102556_pheno_curated_2$PHENOTYPE
GSE102556_pheno_curated_2$PMI
GSE102556_pheno_curated_2$RIN

GSE102556_pheno_curated_2$AGE = as.numeric(GSE102556_pheno_curated_2$AGE)
GSE102556_pheno_curated_2$ALCOOL = sapply(GSE102556_pheno_curated_2$ALCOOL, function(x){
  
  if (x == "NA"){return(NA)}
  if (x == "YES"){return("YES")}
  if (x == "NO"){return("NO")}
  
})
GSE102556_pheno_curated_2$ALCOOL = factor(GSE102556_pheno_curated_2$ALCOOL, levels = c("NO", "YES"))
GSE102556_pheno_curated_2$SUICIDE = ifelse(GSE102556_pheno_curated_2$CAUSE_OF_DEATH == "SUICIDE", "SUICIDE", "CONTROL")
table(GSE102556_pheno_curated_2$SUICIDE, GSE102556_pheno_curated_2$TISSUE)
"         Anterior Insula (aINS) Cingulate gyrus 25 (Cg25) Dorsolateral prefrontal cortex (dlPFC; BA8/9) Nucleus Accumbens (Nac) Orbitofrontal (OFC; BA11) Subiculum (Sub)
  CONTROL                     11                         8                                            11                      11                        11               8
  SUICIDE                     37                        20                                            37                      37                        36              35"
GSE102556_pheno_curated_2$SUICIDE = factor(GSE102556_pheno_curated_2$SUICIDE, levels = c("CONTROL", "SUICIDE"))
GSE102556_pheno_curated_2$GENDER = factor(GSE102556_pheno_curated_2$GENDER, levels = c("FEMALE", "MALE"))
GSE102556_pheno_curated_2$MEDICATION = sapply(GSE102556_pheno_curated_2$MEDICATION, function(x){
  
  if (x == "NA"){return(NA)}
  if (x == "YES"){return("YES")}
  if (x == "NO"){return("NO")}
  
})
GSE102556_pheno_curated_2$MEDICATION = factor(GSE102556_pheno_curated_2$MEDICATION, levels = c("NO", "YES"))
GSE102556_pheno_curated_2$PH = as.numeric(GSE102556_pheno_curated_2$PH)
GSE102556_pheno_curated_2$PMI = as.numeric(GSE102556_pheno_curated_2$PMI)
GSE102556_pheno_curated_2$RIN = as.numeric(GSE102556_pheno_curated_2$RIN)
GSE102556_pheno_curated_2$PHENOTYPE = factor(GSE102556_pheno_curated_2$PHENOTYPE, levels = c("CTRL", "MDD"))

# Mapping probes
# Mapping RNA seq to ENSEMBL
ENSEML_DF = read.csv("mapping_enseml_all.txt", sep = "\t", header = FALSE)
colnames(ENSEML_DF) = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', "hgnc_symbol")
ENSEML_DF$external_gene_name = ifelse(ENSEML_DF$external_gene_name == "", NA, ENSEML_DF$external_gene_name)
ENSEML_DF$hgnc_symbol = ifelse(ENSEML_DF$hgnc_symbol == "", NA, ENSEML_DF$hgnc_symbol)
ENSEML_DF$NIH_check = mapply(function(x, y){
  
  if (!is.na(y)){
    # Y is available (HGNC symbol)
    gene_symbol_NIH = check_gene_symbol_NIH(PRF_gene_symbols = y, 
                                            PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                            PRF_replace_NA_with_old = TRUE)
    gene_symbol_NIH = gene_symbol_NIH$Suggested.Symbol
    
  } else {
    # Y is NOT available (HGNC symbol) -> using other symbol
    
    gene_symbol_NIH = check_gene_symbol_NIH(PRF_gene_symbols = x, 
                                            PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                            PRF_replace_NA_with_old = TRUE)
    gene_symbol_NIH = gene_symbol_NIH$Suggested.Symbol
    
  }
  
  return(gene_symbol_NIH)
  
}, ENSEML_DF$external_gene_name, ENSEML_DF$hgnc_symbol)
write.csv(ENSEML_DF, "ENSEML_DF.csv")
ENSEML_DF_dipl = ENSEML_DF[ENSEML_DF$ensembl_gene_id %in% ENSEML_DF$ensembl_gene_id[duplicated(ENSEML_DF$ensembl_gene_id)], ]


RNAseq_genes = rownames(GSE102556_count)[-(1:4)]
RNAseq_gene_symbols = mapIds(org.Hs.eg.db, keys = RNAseq_genes, keytype = "ENSEMBL", column = "SYMBOL")
RNAseq_gene_probes = data.frame(ID = RNAseq_genes, Base_symbol = RNAseq_gene_symbols)
RNAseq_gene_probes$Gene_symbol = sapply(RNAseq_gene_probes$ID, function(x){
  value = ENSEML_DF[ENSEML_DF$ensembl_gene_id == x, "NIH_check"]
  value = as.character(value)
  value = unique(value)
  value = paste0(value, collapse = ";")
  return(value)
})

rownames(RNAseq_gene_probes) = RNAseq_gene_probes$ID
write.csv(RNAseq_gene_probes, "RNAseq_gene_probes.csv")

# Analysis without covariates per tissue in a loop

"Threshold for filtering out low expressed genes was set to >5 reads in at least 80% of the samples as previously described"

GSE102556_tissues = unique(GSE102556_pheno_curated_2$TISSUE)

GSE102556_analysis_list_no_covar = list()

for (x in 1:length(GSE102556_tissues)){
  
  
  GSE102556_count_tmp = GSE102556_count[-(1:4), ]
  GSE102556_TMP_tissue = GSE102556_tissues[x]
  
  print(paste0("Working on: ", GSE102556_TMP_tissue))
  
  GSE102556_pheno_curated_2_tmp = GSE102556_pheno_curated_2[GSE102556_pheno_curated_2$TISSUE == GSE102556_TMP_tissue,]
  GSE102556_count_tmp = GSE102556_count_tmp[,GSE102556_pheno_curated_2_tmp$RUN_ID]
  
  # 80% filter
  #GSE102556_TMP_transcr_sums = rowSums(GSE102556_count_tmp)
  #quantile(GSE102556_TMP_transcr_sums, probs = seq(from=0.1, to=1, by=0.05))
  #GSE102556_TMP_transcr_sums_more_than_5 = apply(GSE102556_count_tmp, 2, function(x) as.numeric(x > 5))  # FALSE become 0, and TRUE becomes 1
  #GSE102556_TMP_row_selection = which(rowSums(GSE102556_count_tmp) >= 0.8 * ncol(GSE102556_TMP_transcr_sums_more_than_5))
  #GSE102556_count_tmp = GSE102556_count_tmp[GSE102556_TMP_row_selection,]
  
  # Design matrix
  tmp_design = model.matrix(~ SUICIDE, data = GSE102556_pheno_curated_2_tmp)
  
  # egeR filter
  GSE102556_tmp_dge = DGEList(counts=GSE102556_count_tmp)
  keep = filterByExpr(GSE102556_tmp_dge, tmp_design)
  GSE102556_tmp_dge = GSE102556_tmp_dge[keep,,keep.lib.sizes=FALSE]
  GSE102556_tmp_dge = calcNormFactors(GSE102556_tmp_dge)
  GSE102556_tmp_voom = voom(GSE102556_tmp_dge, tmp_design, plot=TRUE)
  
  # limma 
  fit = lmFit(GSE102556_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  GSE102556_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE102556_Top_table_no_covar_TMP$ID = rownames(GSE102556_Top_table_no_covar_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE102556_Top_table_no_covar_TMP$ID]
  GSE102556_Top_table_no_covar_TMP$SE = SE
  all(names(SE) == GSE102556_Top_table_no_covar_TMP$ID) # TRUE
  
  # Annotating results
  GSE102556_probes_TMP = RNAseq_gene_probes
  GSE102556_probes_TMP = GSE102556_probes_TMP[GSE102556_Top_table_no_covar_TMP$ID, ]
  GSE102556_Top_table_no_covar_TMP$Gene_symbol = GSE102556_probes_TMP$Gene_symbol
  GSE102556_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
  GSE102556_Top_table_no_covar_TMP$Tissue = GSE102556_TMP_tissue
  GSE102556_Top_table_no_covar_TMP$Tissue_type = "Brain"
  GSE102556_Top_table_no_covar_TMP$Technology = "RNA-seq"
  
  GSE102556_analysis_list_no_covar[[x]] = GSE102556_Top_table_no_covar_TMP
  
}

GSE102556_Top_table_no_covar = do.call(rbind, GSE102556_analysis_list_no_covar)
rownames(GSE102556_Top_table_no_covar) = NULL

# analysis with covariates in a loop
GSE102556_tissues = unique(GSE102556_pheno_curated_2$TISSUE)

GSE102556_analysis_list = list()

for (x in 1:length(GSE102556_tissues)){
  
  
  GSE102556_count_tmp = GSE102556_count[-(1:4), ]
  GSE102556_TMP_tissue = GSE102556_tissues[x]
  
  print(paste0("Working on: ", GSE102556_TMP_tissue))
  
  GSE102556_pheno_curated_2_tmp = GSE102556_pheno_curated_2[GSE102556_pheno_curated_2$TISSUE == GSE102556_TMP_tissue,]
  GSE102556_count_tmp = GSE102556_count_tmp[,GSE102556_pheno_curated_2_tmp$RUN_ID]
  
  # 80% filter (too conservative)
  #GSE102556_TMP_transcr_sums = rowSums(GSE102556_count_tmp)
  #quantile(GSE102556_TMP_transcr_sums, probs = seq(from=0.1, to=1, by=0.05))
  #GSE102556_TMP_transcr_sums_more_than_5 = apply(GSE102556_count_tmp, 2, function(x) as.numeric(x > 5))  # FALSE become 0, and TRUE becomes 1
  #GSE102556_TMP_row_selection = which(rowSums(GSE102556_count_tmp) >= 0.8 * ncol(GSE102556_TMP_transcr_sums_more_than_5))
  #GSE102556_count_tmp = GSE102556_count_tmp[GSE102556_TMP_row_selection,]
  
  # Design matrix
  tmp_design = model.matrix(~ SUICIDE + PHENOTYPE + GENDER + AGE + ALCOOL + MEDICATION + RIN + PMI + PH, data = GSE102556_pheno_curated_2_tmp)
  GSE102556_pheno_curated_2_tmp_small = GSE102556_pheno_curated_2_tmp[GSE102556_pheno_curated_2_tmp$GEO_ACCESSION %in% rownames(tmp_design),]
  table(GSE102556_pheno_curated_2_tmp_small$SUICIDE, GSE102556_pheno_curated_2_tmp_small$ALCOOL)
  table(GSE102556_pheno_curated_2_tmp_small$SUICIDE, GSE102556_pheno_curated_2_tmp_small$MEDICATION)
  GSE102556_count_tmp = GSE102556_count_tmp[,GSE102556_pheno_curated_2_tmp_small$RUN_ID]
  
  # egeR filter
  GSE102556_tmp_dge = DGEList(counts=GSE102556_count_tmp)
  keep = filterByExpr(GSE102556_tmp_dge, tmp_design)
  GSE102556_tmp_dge = GSE102556_tmp_dge[keep,,keep.lib.sizes=FALSE]
  GSE102556_tmp_dge = calcNormFactors(GSE102556_tmp_dge)
  GSE102556_tmp_voom = voom(GSE102556_tmp_dge, tmp_design, plot=TRUE)
  
  # limma 
  fit = lmFit(GSE102556_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  GSE102556_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE102556_Top_table_TMP$ID = rownames(GSE102556_Top_table_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE102556_Top_table_TMP$ID]
  GSE102556_Top_table_TMP$SE = SE
  all(names(SE) == GSE102556_Top_table_TMP$ID) # TRUE
  
  # Annotating results
  GSE102556_probes_TMP = RNAseq_gene_probes
  GSE102556_probes_TMP = GSE102556_probes_TMP[GSE102556_Top_table_TMP$ID, ]
  GSE102556_Top_table_TMP$Gene_symbol = GSE102556_probes_TMP$Gene_symbol
  GSE102556_Top_table_TMP$Gene_symbol_non_hgnc = NA
  GSE102556_Top_table_TMP$Tissue = GSE102556_TMP_tissue
  GSE102556_Top_table_TMP$Tissue_type = "Brain"
  GSE102556_Top_table_TMP$Technology = "RNA-seq"
  
  GSE102556_analysis_list[[x]] = GSE102556_Top_table_TMP
}

GSE102556_Top_table = do.call(rbind, GSE102556_analysis_list)
rownames(GSE102556_Top_table) = NULL

# saving results
dir.create("GSE102556_results")

write.csv(GSE102556_count, "GSE102556_results/GSE102556_expression.csv")
write.csv(RNAseq_gene_probes, "GSE102556_results/GSE102556_probes.csv")
write.csv(GSE102556_pheno_curated_2, "GSE102556_results/GSE102556_pheno_curated.csv")
write.csv(GSE102556_Top_table_no_covar, "GSE102556_results/GSE102556_Top_table_no_covar.csv")
write.csv(GSE102556_Top_table, "GSE102556_results/GSE102556_Top_table.csv")
rm(list = ls(pattern = "GSE102556"))
gc()


################### GSE243356 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243356

# paper: https://pubmed.ncbi.nlm.nih.gov/37938766/

"Raw RNA-sequencing reads were first trimmed using program Trim Galore (v0.6) to remove low quality bases 
and adapter sequences [30]. Next, trimmed reads were aligned to human reference genome (hg38) using STAR (2.7.3a) [31]. 
Count tables of all samples were then imported into edgeR (v2.34.1) [32]. Normalization was carried out using “calcNormFactors” 
function in edgeR which utilizes TMM method (trimmed mean of M-values) to estimate scale factors between samples [32]. 
Genes that had less than 0.5 count per million (cpm) normalized counts and presented in less than 20 (30% of all the samples) 
samples were excluded from the downstream analysis. "

GSE243356_geo_record = getGEO("GSE243356")
GSE243356_geo_record = GSE243356_geo_record[[1]]
GSE243356_pheno = pData(GSE243356_geo_record)
GSE243356_pheno = fix_columns(GSE243356_pheno)


# curating pheno
GSE243356_pheno_curated = GSE243356_pheno
table(GSE243356_pheno_curated$ORGANISM_CH1) # 61 human
GSE243356_pheno_curated$GROUP = factor(GSE243356_pheno_curated$GROUP, levels = c("healthy", "suicide"))


# geting SRA runs
GSE243356_SRA_run_table = smart_fread("GSE243356_SraRunTable.txt")
all(GSE243356_pheno_curated$GEO_ACCESSION %in% GSE243356_SRA_run_table$`Sample Name`) # TRUE
all(GSE243356_SRA_run_table$`Sample Name` %in% GSE243356_pheno_curated$GEO_ACCESSION) # TRUE
GSE243356_pheno_curated$RUN_ID = sapply(GSE243356_pheno_curated$GEO_ACCESSION, function(x){
  x = GSE243356_SRA_run_table[GSE243356_SRA_run_table$`Sample Name` == x, "Run"]
  return(x)
})
GSE243356_pheno_curated$TISSUE = "temporal cortex (BA20 and BA36)"

"Libraries were prepared by the Van Andel Genomics Core from 500 ng of total RNA using the KAPA RNA HyperPrep Kit (Kapa Biosystems, Wilmington, MA USA).
Ribosomal RNA material was reduced using the QIAseq FastSelect –rRNA HMR Kit to remove cytoplasmic and mitochondrial rRNAs and QIAseq FastSelect –5S/16S/23S Kit 
to remove pan bacterial rRNAs (Qiagen, Germantown, MD, USA). RNA was sheared to 300-400 bp. Prior to PCR amplification, cDNA fragments were ligated to adapters for 
Illumina TruSeq UD Indexed adapters (Illumina Inc, San Diego CA, USA). Quality and quantity of the finished libraries were assessed using a combination of Agilent
DNA High Sensitivity chip (Agilent Technologies, Inc.), QuantiFluor® dsDNA System (Promega Corp., Madison, WI, USA), and Kapa Illumina Library Quantification qPCR 
assays (Kapa Biosystems). Barcodes that are unique to each sample (from each patient) were added to each DNA fragment; and resulting DNA libraries were pooled 
(combined) in the sequencing run. The sequencing was performed on an Illumina NovaSeq6000 targeting depth of 100 M per sample with 100 bp, paired-end configuration. 
Base calling was done by Illumina Real Time Analysis and output of NextSeq Control Software was demultiplexed and converted to
Fastq format with Illumina Bcl2fastq v1.9.0."

####
####
GSE243356_pheno_curated$RUN_ID[1] # SRR26073596
# PAIRED LAYOUT

# prefetch download
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE243356_pref
source /etc/profile.d/sra-tools.sh
prefetch SRR26073596 -O GSE243356_pref

"

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE243356_pref
source /etc/profile.d/sra-tools.sh
fasterq-dump --threads 9 SRR26073596

"

# FastQC
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR5961796_QC
FastQC/fastqc GSE243356_pref/SRR26073596_1.fastq -o SRR5961796_QC/

"

# STAR (no trimming) Uniquely mapped reads number |	59697851 Uniquely mapped reads % |	89.68%  % of reads mapped to multiple loci |	5.92% % of reads unmapped: too short |	4.07%
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE243356_mapped

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE243356_pref/SRR26073596_1.fastq GSE243356_pref/SRR26073596_2.fastq \
--outFileNamePrefix GSE243356_mapped/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# trimming with fastp
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE243356_trimmed
conda activate trimming
fastp --thread 8 -i GSE243356_pref/SRR26073596_1.fastq -I GSE243356_pref/SRR26073596_2.fastq -o GSE243356_trimmed/out.R1.fq.gz -O GSE243356_trimmed/out.R2.fq.gz

"

# STAR after trimming -> only minor improvements Uniquely mapped reads number |	59618764 Uniquely mapped reads % |	90.11% % of reads mapped to multiple loci |	5.95% % of reads unmapped: too short |	3.61%
# skipping trimming is OK as difference is minor

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE243356_mapped_trimmed

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE243356_trimmed/out.R1.fq.gz GSE243356_trimmed/out.R2.fq.gz --readFilesCommand gunzip -c \
--outFileNamePrefix GSE243356_mapped_trimmed/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

####
####


any(duplicated(GSE243356_pheno_curated$RUN_ID)) # FALSE
writeLines(GSE243356_pheno_curated$RUN_ID, "GSE243356_SRR_ids.txt")

# Obtaining counts for ALL runs (Terminal). It takes approx. 10 minutes per sample
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTAR.py GSE243356_SRR_ids.txt PAIRED GSE243356_

"

# Inspecting runs
GSE243356_run_logs = list.files("GSE243356_OUTPUT", pattern = "logs", full.names = TRUE)
GSE243356_mapped_percent = sapply(GSE243356_run_logs, function(x){
  x = readLines(x)
  x = x[10]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads % |\t", replacement = "")
  x = stri_replace_all_fixed(str = x, pattern = "%", replacement = "")
  return(x)
})
GSE243356_mapped_percent = as.numeric(GSE243356_mapped_percent)
min(GSE243356_mapped_percent) # 81.49 %
max(GSE243356_mapped_percent) # 91.21 %
GSE243356_mapped_percent_df = data.frame(file = GSE243356_run_logs, percent = GSE243356_mapped_percent)
# all samples have relatively high percent of mapped reads

# preparing counts data
GSE243356_runs = list.files("GSE243356_OUTPUT", pattern = "counts", full.names = TRUE)
GSE243356_count_df = lapply(GSE243356_runs, function(x){
  x = smart_fread(x)
  return(x)
})
names(GSE243356_count_df) = GSE243356_runs
reference = GSE243356_count_df[[1]]$gene_ID
GSE243356_compar_to_ref = sapply(GSE243356_count_df, function(x){
  compar = reference == x$gene_ID
  compar = all(compar)
  return(compar)
})
all(GSE243356_compar_to_ref) # all row names match

# inspecting stranding

GSE243356_participant_1 = GSE243356_count_df[[1]]
GSE243356_participant_1 = GSE243356_participant_1[-(1:4), ]
GSE243356_participant_1$gene_ID = NULL
GSE243356_participant_1 = GSE243356_participant_1[rowSums(GSE243356_participant_1) > 0, ]
table(GSE243356_participant_1$Strand_1 > GSE243356_participant_1$Strand_2)
table(GSE243356_participant_1$Strand_1 < GSE243356_participant_1$Strand_2)
# Strand 2 is systematically bigger than 1
sum(GSE243356_participant_1$Strand_1) # 1783255
sum(GSE243356_participant_1$Strand_2) # 31940659
sum(GSE243356_participant_1$Strand_1)/sum(GSE243356_participant_1$Strand_2) # 0.05583025
sum(GSE243356_participant_1$Strand_2)/sum(GSE243356_participant_1$Strand_1) # 17.91144 

GSE243356_count = lapply(GSE243356_count_df, function(x){
  x = x[,"Strand_2"]
  return(x)
})
GSE243356_count = do.call(cbind,GSE243356_count)
rownames(GSE243356_count) = reference
colnames(GSE243356_count) = sapply(colnames(GSE243356_count), function(x){
  x = stri_replace_all_fixed(x, pattern = "GSE243356_OUTPUT/", replacement = "")
  x = stri_replace_all_fixed(x, pattern = "_star_counts.csv", replacement = "")
  return(x)
})
GSE243356_count = as.data.frame(GSE243356_count)

# check
all(colnames(GSE243356_count) == GSE243356_pheno_curated$RUN_ID) # FALSE
all(colnames(GSE243356_count) %in% GSE243356_pheno_curated$RUN_ID) # TRUE

GSE243356_count = GSE243356_count[,GSE243356_pheno_curated$RUN_ID]
all(colnames(GSE243356_count) == GSE243356_count$RUN_ID) # TRUE

# covariates are not available...

# analysis without covariates

GSE243356_count_tmp = GSE243356_count[-(1:4), ]

# 80% filter (too conservative)
# GSE243356_TMP_transcr_sums = rowSums(GSE243356_count_tmp)
# quantile(GSE243356_TMP_transcr_sums, probs = seq(from=0.1, to=1, by=0.05))
# GSE243356_TMP_transcr_sums_more_than_5 = apply(GSE243356_count_tmp, 2, function(x) as.numeric(x > 5))  # FALSE become 0, and TRUE becomes 1
# GSE243356_TMP_row_selection = which(rowSums(GSE243356_count_tmp) >= 0.8 * ncol(GSE243356_TMP_transcr_sums_more_than_5))
# GSE243356_count_tmp = GSE243356_count_tmp[GSE243356_TMP_row_selection,]

# Design matrix
tmp_design = model.matrix(~ GROUP, data = GSE243356_pheno_curated)

# egeR filter
GSE243356_tmp_dge = DGEList(counts=GSE243356_count_tmp)
keep = filterByExpr(GSE243356_tmp_dge, tmp_design)
GSE243356_tmp_dge = GSE243356_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE243356_tmp_dge = calcNormFactors(GSE243356_tmp_dge)
GSE243356_tmp_voom = voom(GSE243356_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE243356_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE243356_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE243356_Top_table_no_covar_TMP$ID = rownames(GSE243356_Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE243356_Top_table_no_covar_TMP$ID]
GSE243356_Top_table_no_covar_TMP$SE = SE
all(names(SE) == GSE243356_Top_table_no_covar_TMP$ID) # TRUE

# Annotating results
GSE243356_probes_TMP = RNAseq_gene_probes
GSE243356_probes_TMP = GSE243356_probes_TMP[GSE243356_Top_table_no_covar_TMP$ID, ]
GSE243356_Top_table_no_covar_TMP$Gene_symbol = GSE243356_probes_TMP$Gene_symbol
GSE243356_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
GSE243356_Top_table_no_covar_TMP$Tissue = GSE243356_pheno_curated$TISSUE[1]
GSE243356_Top_table_no_covar_TMP$Tissue_type = "Brain"
GSE243356_Top_table_no_covar_TMP$Technology = "RNA-seq"

# Analysis with covariates is the same as there are no other phenotypes
GSE243356_Top_table = GSE243356_Top_table_no_covar_TMP
GSE243356_Top_table_no_covar = GSE243356_Top_table_no_covar_TMP


dir.create("GSE243356_results")
write.csv(GSE243356_count, "GSE243356_results/GSE243356_expression.csv")
write.csv(RNAseq_gene_probes, "GSE243356_results/GSE243356_probes.csv")
write.csv(GSE243356_pheno_curated, "GSE243356_results/GSE243356_pheno_curated.csv")
write.csv(GSE243356_Top_table_no_covar, "GSE243356_results/GSE243356_Top_table_no_covar.csv")
write.csv(GSE243356_Top_table, "GSE243356_results/GSE243356_Top_table.csv")

rm(list = ls(pattern = "GSE243356"))
gc()

################### GSE248260 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE248260
# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11189724/


"For live samples, PAXgene blood tubes were thawed quickly and brought to room temperature prior to total RNA 
isolation according to the manufacturer’s instructions, and Globin mRNA was removed using GLOBINclear (Invitrogen/Ambion). 
For postmortem human brain specimens, total RNA was isolated from ~30–50 mg white matter tissue using RNeasy Lipid Tissue
Mini Kit (Qiagen, #74804) according to the manufacturer’s instructions. RNA quality was measured via BioAnalyzer and RNA integrity Number (RIN) 
for live and white matter postmortem samples is reported in Table 1 and Table S2, respectively. All RNA samples were subject to Ribo-zero depletion, 
and libraries were prepared using the Illumina Truseq library preparation kit and sequenced on the Illumina HiSeq 2500 (2 × 50), generating a mean of 48
million reads per sample. Reads were aligned to hg19 using STAR aligner and annotated to transcripts using Gencode v18 annotation."

GSE248260_geo_record = getGEO("GSE248260")
GSE248260_geo_record = GSE248260_geo_record[[1]]
GSE248260_pheno = pData(GSE248260_geo_record)
GSE248260_pheno = fix_columns(GSE248260_pheno)


# curating pheno
GSE248260_pheno_curated = GSE248260_pheno
table(GSE248260_pheno_curated$ORGANISM_CH1) # 24 human

GSE248260_pheno_curated$AGE = as.numeric(GSE248260_pheno_curated$AGE)
GSE248260_pheno_curated$DIAGNOSIS = factor(GSE248260_pheno_curated$DIAGNOSIS, levels = c("Normal", "MDD"))
GSE248260_pheno_curated$SEX = ifelse(GSE248260_pheno_curated$SEX == "F", "Female", "Male")
GSE248260_pheno_curated$SEX = factor(GSE248260_pheno_curated$SEX, levels = c("Female", "Male"))
GSE248260_pheno_curated$SUICIDE = ifelse(GSE248260_pheno_curated$SUICIDE == "Y", "Suicide", "Control")
GSE248260_pheno_curated$SUICIDE = factor(GSE248260_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))

# getting_extra_pheno
GSE248260_pheno_extra = read.xlsx("GSE248260_files_pheno/41380_2024_2420_MOESM3_ESM.xlsx", sheet = 1)
all(GSE248260_pheno_extra$SampleID %in% GSE248260_pheno_curated$TITLE) # TRUE
any(duplicated(GSE248260_pheno_extra$SampleID)) # FALSE
rownames(GSE248260_pheno_extra) = GSE248260_pheno_extra$SampleID
GSE248260_pheno_extra = GSE248260_pheno_extra[GSE248260_pheno_curated$TITLE,]
GSE248260_pheno_extra_2 = GSE248260_pheno_extra[,c("RIN", "PMI", "pH", 
                                                   "Brain_tox_antidepressant",
                                                   "Brain_tox_antipsychotic",
                                                   "Brain_tox_benzodiazepine",
                                                   "Brain_tox_opioid")]
GSE248260_pheno_curated_2 = cbind(GSE248260_pheno_curated, GSE248260_pheno_extra_2)
GSE248260_pheno_curated_2$RIN
GSE248260_pheno_curated_2$PMI
GSE248260_pheno_curated_2$pH
GSE248260_pheno_curated_2$Drug_intake = mapply(function(a,b,c,d){
  combined_vector = c(a, b, c, d)
  
  if (any(combined_vector == "Yes")){
    
    return("Drugs")
    
  } else {
    return("Not reported")
  }
  
}, GSE248260_pheno_curated_2$Brain_tox_antidepressant, GSE248260_pheno_curated_2$Brain_tox_antipsychotic, 
GSE248260_pheno_curated_2$Brain_tox_benzodiazepine, GSE248260_pheno_curated_2$Brain_tox_opioid)
table(GSE248260_pheno_curated_2$Drug_intake, GSE248260_pheno_curated_2$SUICIDE)
"               Control Suicide
  Not reported       7       6
  Drugs              2       9
"
table(GSE248260_pheno_curated_2$Drug_intake, GSE248260_pheno_curated_2$DIAGNOSIS)
table(GSE248260_pheno_curated_2$SUICIDE, GSE248260_pheno_curated_2$DIAGNOSIS) # Diagnosis is 100% related to suicide
"
          Normal MDD
  Control      9   0
  Suicide      0  15
"
GSE248260_pheno_curated_2$Drug_intake = factor(GSE248260_pheno_curated_2$Drug_intake, levels = c("Not reported", "Drugs"))

# geting SRA runs
GSE248260_SRA_run_table = smart_fread("GSE248260_SraRunTable.txt")
all(GSE248260_pheno_curated_2$GEO_ACCESSION %in% GSE248260_SRA_run_table$`Sample Name`) # TRUE
all(GSE248260_SRA_run_table$`Sample Name` %in% GSE248260_pheno_curated_2$GEO_ACCESSION) # TRUE
GSE248260_pheno_curated_2$RUN_ID = sapply(GSE248260_pheno_curated_2$GEO_ACCESSION, function(x){
  x = GSE248260_SRA_run_table[GSE248260_SRA_run_table$`Sample Name` == x, "Run"]
  return(x)
})

GSE248260_pheno_curated_2$RUN_ID[1] # SRR26894803
GSE248260_pheno_curated_2$TITLE[1] # 104A
writeLines(GSE248260_pheno_curated_2$RUN_ID, "GSE248260_SRR_ids.txt")


# Layout SINGLE (based on SRA)
# However, it looks like paired based on the reads per spot and data
# prefetch download

# Download of main file is usually finished after: 2024-08-08T10:23:12 prefetch.3.1.1: 1) 'SRR26894803' was downloaded successfully
# Termination at this step does not work as it continues download...

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE248260_pref
source /etc/profile.d/sra-tools.sh
prefetch SRR26894803 -O GSE248260_pref

"
# it can throw an error (sometimes):

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE248260_pref
source /etc/profile.d/sra-tools.sh
fasterq-dump --threads 9 SRR26894803

"

# inspection 
"
head SRR26894803_1.fastq
head SRR26894803_2.fastq
"
# FastQC
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR26894803_QC
FastQC/fastqc GSE248260_pref/SRR26894803_1.fastq -o SRR26894803_QC/

"
# STAR (no trimming)
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE248260_mapped

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE248260_pref/SRR26894803_1.fastq GSE248260_pref/SRR26894803_2.fastq \
--outFileNamePrefix GSE248260_mapped/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf
"
# Finished correctly but the number of reads is double?!

# Comparison with reported counts for 104A

# ENSG00000227232 We: 121 Paper: 190
# ENSG00000228794 We: 430 Paper: 815
# ENSG00000078808 We: 2123 Paper: 4197
# ENSG00000157933 We: 5502 Paper: 10921
# ENSG00000196581 We: 2337 Paper: 4612
# ENSG00000196586 We: 13741 Paper: 26125
# ENSG00000112706 We: 0 Paper: 0
# ENSG00000083123 We 1333 Paper: 2595
# ENSG00000112159 We 11590 Paper: 22651
# Paper reports double counts?!

# ENSG00000112297 We 150 Paper: 142
# ENSG00000197498 We 348 Paper: 659
# ENSG00000203778 We 1403 Paper: 356
# ENSG00000112769 We 906 Paper: 755
# ENSG00000155130 We 6122 Paper: 12233
# ENSG00000047932 We 2115 Paper 4089

# -> we preprocess them as paired-end as method description indicates around 50M reads (x2) -> likely paired

# Obtaining counts for ALL runs (Terminal). It takes approx. 25 minutes per sample

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARprefSEP.py GSE248260_SRR_ids.txt PAIRED GSE248260_

"

# Inspecting runs
GSE248260_run_logs = list.files("GSE248260_OUTPUT", pattern = "logs", full.names = TRUE)
GSE248260_mapped_percent = sapply(GSE248260_run_logs, function(x){
  x = readLines(x)
  x = x[10]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads % |\t", replacement = "")
  x = stri_replace_all_fixed(str = x, pattern = "%", replacement = "")
  return(x)
})
GSE248260_mapped_percent = as.numeric(GSE248260_mapped_percent)
min(GSE248260_mapped_percent) # 82.53 %
max(GSE248260_mapped_percent) # 87.88 %
GSE248260_mapped_percent_df = data.frame(file = GSE248260_run_logs, percent = GSE248260_mapped_percent)
# all samples have relatively high percent of mapped reads


# preparing counts data
GSE248260_runs = list.files("GSE248260_OUTPUT", pattern = "counts", full.names = TRUE)
GSE248260_count_df = lapply(GSE248260_runs, function(x){
  x = smart_fread(x)
  return(x)
})
names(GSE248260_count_df) = GSE248260_runs
reference = GSE248260_count_df[[1]]$gene_ID
GSE248260_compar_to_ref = sapply(GSE248260_count_df, function(x){
  compar = reference == x$gene_ID
  compar = all(compar)
  return(compar)
})
all(GSE248260_compar_to_ref) # all row names match


GSE248260_participant_1 = GSE248260_count_df[[1]]
GSE248260_participant_1 = GSE248260_participant_1[-(1:4), ]
GSE248260_participant_1$gene_ID = NULL
GSE248260_participant_1 = GSE248260_participant_1[rowSums(GSE248260_participant_1) > 0, ]
table(GSE248260_participant_1$Strand_1 > GSE248260_participant_1$Strand_2)
table(GSE248260_participant_1$Strand_1 < GSE248260_participant_1$Strand_2)
# Strand 2 is systematically bigger than 1
sum(GSE248260_participant_1$Strand_1) # 2529887
sum(GSE248260_participant_1$Strand_2) # 33957026
sum(GSE248260_participant_1$Strand_1)/sum(GSE248260_participant_1$Strand_2) # 0.07450261
sum(GSE248260_participant_1$Strand_2)/sum(GSE248260_participant_1$Strand_1) # 13.42235 

GSE248260_count = lapply(GSE248260_count_df, function(x){
  x = x[,"Strand_2"]
  return(x)
})
GSE248260_count = do.call(cbind,GSE248260_count)
rownames(GSE248260_count) = reference
colnames(GSE248260_count) = sapply(colnames(GSE248260_count), function(x){
  x = stri_replace_all_fixed(x, pattern = "GSE248260_OUTPUT/", replacement = "")
  x = stri_replace_all_fixed(x, pattern = "_star_counts.csv", replacement = "")
  return(x)
})
GSE248260_count = as.data.frame(GSE248260_count)

# check
all(colnames(GSE248260_count) == GSE248260_pheno_curated_2$RUN_ID) # FALSE
all(colnames(GSE248260_count) %in% GSE248260_pheno_curated_2$RUN_ID) # TRUE

GSE248260_count = GSE248260_count[,GSE248260_pheno_curated_2$RUN_ID]
all(colnames(GSE248260_count) == GSE248260_count$RUN_ID) # TRUE


# Last phenotype inspection before analysis
# Selected variables:  AGE PMI pH SUICIDE RIN Drug_intake 
GSE248260_pheno_curated_2$AGE
GSE248260_pheno_curated_2$SEX
GSE248260_pheno_curated_2$PMI
GSE248260_pheno_curated_2$pH
GSE248260_pheno_curated_2$SUICIDE
GSE248260_pheno_curated_2$TISSUE
GSE248260_pheno_curated_2$RIN
GSE248260_pheno_curated_2$Drug_intake

table(GSE248260_pheno_curated_2$DIAGNOSIS, GSE248260_pheno_curated_2$SUICIDE)

"
         Control Suicide
  Normal       9       0
  MDD          0      15
"

# DIAGNOSIS and SUICIDE ARE FULLY CORRELATED AND THUS WE CANT ADJUST FOR IT

# analysis without covariates

GSE248260_count_tmp = GSE248260_count[-(1:4), ]

# 80% filter (too conservative)
# GSE248260_TMP_transcr_sums = rowSums(GSE248260_count_tmp)
# quantile(GSE248260_TMP_transcr_sums, probs = seq(from=0.1, to=1, by=0.05))
# GSE248260_TMP_transcr_sums_more_than_5 = apply(GSE248260_count_tmp, 2, function(x) as.numeric(x > 5))  # FALSE become 0, and TRUE becomes 1
# GSE248260_TMP_row_selection = which(rowSums(GSE248260_count_tmp) >= 0.8 * ncol(GSE248260_TMP_transcr_sums_more_than_5))
# GSE248260_count_tmp = GSE248260_count_tmp[GSE248260_TMP_row_selection,]

# Design matrix
tmp_design = model.matrix(~ SUICIDE, data = GSE248260_pheno_curated_2)

# egeR filter
GSE248260_tmp_dge = DGEList(counts=GSE248260_count_tmp)
keep = filterByExpr(GSE248260_tmp_dge, tmp_design)
GSE248260_tmp_dge = GSE248260_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE248260_tmp_dge = calcNormFactors(GSE248260_tmp_dge)
GSE248260_tmp_voom = voom(GSE248260_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE248260_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE248260_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE248260_Top_table_no_covar_TMP$ID = rownames(GSE248260_Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE248260_Top_table_no_covar_TMP$ID]
GSE248260_Top_table_no_covar_TMP$SE = SE
all(names(SE) == GSE248260_Top_table_no_covar_TMP$ID) # TRUE

# Annotating results
GSE248260_probes_TMP = RNAseq_gene_probes
GSE248260_probes_TMP = GSE248260_probes_TMP[GSE248260_Top_table_no_covar_TMP$ID, ]
GSE248260_Top_table_no_covar_TMP$Gene_symbol = GSE248260_probes_TMP$Gene_symbol
GSE248260_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
GSE248260_Top_table_no_covar_TMP$Tissue = GSE248260_pheno_curated$TISSUE[1]
GSE248260_Top_table_no_covar_TMP$Tissue_type = "Brain"
GSE248260_Top_table_no_covar_TMP$Technology = "RNA-seq"


# analysis with covariates

GSE248260_count_tmp = GSE248260_count[-(1:4), ]

# 80% filter (too conservative)
# GSE248260_TMP_transcr_sums = rowSums(GSE248260_count_tmp)
# quantile(GSE248260_TMP_transcr_sums, probs = seq(from=0.1, to=1, by=0.05))
# GSE248260_TMP_transcr_sums_more_than_5 = apply(GSE248260_count_tmp, 2, function(x) as.numeric(x > 5))  # FALSE become 0, and TRUE becomes 1
# GSE248260_TMP_row_selection = which(rowSums(GSE248260_count_tmp) >= 0.8 * ncol(GSE248260_TMP_transcr_sums_more_than_5))
# GSE248260_count_tmp = GSE248260_count_tmp[GSE248260_TMP_row_selection,]

# Design matrix
tmp_design = model.matrix(~ SUICIDE + AGE + SEX + PMI + pH + RIN + Drug_intake, data = GSE248260_pheno_curated_2)

# egeR filter
GSE248260_tmp_dge = DGEList(counts=GSE248260_count_tmp)
keep = filterByExpr(GSE248260_tmp_dge, tmp_design)
GSE248260_tmp_dge = GSE248260_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE248260_tmp_dge = calcNormFactors(GSE248260_tmp_dge)
GSE248260_tmp_voom = voom(GSE248260_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE248260_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE248260_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE248260_Top_table_TMP$ID = rownames(GSE248260_Top_table_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE248260_Top_table_TMP$ID]
GSE248260_Top_table_TMP$SE = SE
all(names(SE) == GSE248260_Top_table_TMP$ID) # TRUE

# Annotating results
GSE248260_probes_TMP = RNAseq_gene_probes
GSE248260_probes_TMP = GSE248260_probes_TMP[GSE248260_Top_table_TMP$ID, ]
GSE248260_Top_table_TMP$Gene_symbol = GSE248260_probes_TMP$Gene_symbol
GSE248260_Top_table_TMP$Gene_symbol_non_hgnc = NA
GSE248260_Top_table_TMP$Tissue = GSE248260_pheno_curated$TISSUE[1]
GSE248260_Top_table_TMP$Tissue_type = "Brain"
GSE248260_Top_table_TMP$Technology = "RNA-seq"

# Saving
GSE248260_Top_table = GSE248260_Top_table_TMP
GSE248260_Top_table_no_covar = GSE248260_Top_table_no_covar_TMP


dir.create("GSE248260_results")
write.csv(GSE248260_count, "GSE248260_results/GSE248260_expression.csv")
write.csv(RNAseq_gene_probes, "GSE248260_results/GSE248260_probes.csv")
write.csv(GSE248260_pheno_curated_2, "GSE248260_results/GSE248260_pheno_curated.csv")
write.csv(GSE248260_Top_table_no_covar, "GSE248260_results/GSE248260_Top_table_no_covar.csv")
write.csv(GSE248260_Top_table, "GSE248260_results/GSE248260_Top_table.csv")

rm(list = ls(pattern = "GSE248260"))
gc()


################### GSE247998 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247998
# paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11189724/


"Demographics and relevant clinical measurements for the N = 100 live subjects are shown in Table 1, 
by recruitment group. Groups did not differ by age, sex, race, ethnicity, smoking status, and RNA quality,
but as expected, differed on depression and suicidal ideation scales."

# Extra phenotypes are not available and not provided as supplement

GSE247998_geo_record = getGEO("GSE247998")
GSE247998_geo_record = GSE247998_geo_record[[1]]
GSE247998_pheno = pData(GSE247998_geo_record)
GSE247998_pheno = fix_columns(GSE247998_pheno)


# curating pheno
GSE247998_pheno_curated = GSE247998_pheno
table(GSE247998_pheno_curated$ORGANISM_CH1) # 100 human

# variables
GSE247998_pheno_curated$AGE
GSE247998_pheno_curated$DIAGNOSIS
GSE247998_pheno_curated$SEX
GSE247998_pheno_curated$SI_GROUP # PATH OF CAUSAL PATHWAY?!
GSE247998_pheno_curated$SUICIDE_ATTEMPT
GSE247998_pheno_curated$TISSUE


GSE247998_pheno_curated$AGE = as.numeric(GSE247998_pheno_curated$AGE)
GSE247998_pheno_curated$DIAGNOSIS = factor(GSE247998_pheno_curated$DIAGNOSIS, levels = c("Healthy control", "MDD"))
GSE247998_pheno_curated$SEX = ifelse(GSE247998_pheno_curated$SEX == "F", "Female", "Male")
GSE247998_pheno_curated$SEX = factor(GSE247998_pheno_curated$SEX, levels = c("Female", "Male"))
GSE247998_pheno_curated$SUICIDE = ifelse(GSE247998_pheno_curated$SUICIDE_ATTEMPT == "Y", "Suicide", "Control")
GSE247998_pheno_curated$SUICIDE = factor(GSE247998_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
table(GSE247998_pheno_curated$SUICIDE, GSE247998_pheno_curated$DIAGNOSIS)
table(GSE247998_pheno_curated$SUICIDE, GSE247998_pheno_curated$SI_GROUP)


# geting SRA runs
GSE247998_SRA_run_table = smart_fread("GSE247998_SraRunTable.txt")
all(GSE247998_pheno_curated$GEO_ACCESSION %in% GSE247998_SRA_run_table$`Sample Name`) # TRUE
all(GSE247998_SRA_run_table$`Sample Name` %in% GSE247998_pheno_curated$GEO_ACCESSION) # TRUE
GSE247998_pheno_curated$RUN_ID = sapply(GSE247998_pheno_curated$GEO_ACCESSION, function(x){
  x = GSE247998_SRA_run_table[GSE247998_SRA_run_table$`Sample Name` == x, "Run"]
  return(x)
})

GSE247998_pheno_curated$RUN_ID[1] # SRR26854272
writeLines(GSE247998_pheno_curated$RUN_ID, "GSE247998_SRR_ids.txt")



# prefetch download
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE247998_pref
source /etc/profile.d/sra-tools.sh
prefetch SRR26854272 -O GSE247998_pref

"

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE247998_pref
source /etc/profile.d/sra-tools.sh
fasterq-dump --threads 9 SRR26854272

"

# FastQC
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR26854272_QC
FastQC/fastqc GSE247998_pref/SRR26854272_1.fastq -o SRR26854272_QC/

"

# STAR (no trimming) Uniquely mapped reads number 
# Uniquely mapped reads % |	86.31%   Average mapped length |	99.72 % of reads mapped to multiple loci |	9.43% % of reads mapped to too many loci |	0.09% % of reads unmapped: too short |	4.02%
# Number of input reads |	48869976
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE247998_mapped

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE247998_pref/SRR26854272_1.fastq GSE247998_pref/SRR26854272_2.fastq \
--outFileNamePrefix GSE247998_mapped/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# trimming with fastp
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE247998_trimmed
conda activate trimming
fastp --thread 8 -i GSE247998_pref/SRR26854272_1.fastq -I GSE247998_pref/SRR26854272_2.fastq -o GSE247998_trimmed/out.R1.fq.gz -O GSE247998_trimmed/out.R2.fq.gz

"

# STAR after trimming: niquely mapped reads % |	86.76% (almost no difference) Number of input reads |	48581503
# skipping trimming is OK as difference is minor

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE247998_mapped_trimmed

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE247998_trimmed/out.R1.fq.gz GSE247998_trimmed/out.R2.fq.gz --readFilesCommand gunzip -c \
--outFileNamePrefix GSE247998_mapped_trimmed/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# Obtaining counts for ALL runs (Terminal). It takes approx. 25 minutes per sample

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARprefSEP.py GSE247998_SRR_ids.txt PAIRED GSE247998_

"

# Pipeline has failed without raising an error -> rerun with more control for remaining samples
GSE247998_runs_batch_1 = list.files("GSE247998_OUTPUT", pattern = "counts")
GSE247998_runs_batch_1 = stri_replace_all_fixed(GSE247998_runs_batch_1, pattern = "_star_counts.csv", replacement = "")
GSE247998_runs_not_finished = GSE247998_pheno_curated$RUN_ID[GSE247998_pheno_curated$RUN_ID %!in% GSE247998_runs_batch_1]
GSE247998_runs_not_finished %in% GSE247998_runs_batch_1
writeLines(GSE247998_runs_not_finished, "GSE247998_SRR_ids_batch_2.txt")

GSE247998_runs_not_finished_2 = GSE247998_runs_not_finished[GSE247998_runs_not_finished != "SRR26854217"]
writeLines(GSE247998_runs_not_finished_2, "GSE247998_SRR_ids_batch_3.txt")
# SRR26854217 crashes every time -> need to obtain manually! possibly through ENA

# Obtaining counts for missing runs (Terminal). It takes approx. 25 minutes per sample
# Missing runs worked fined besides SRR26854217
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARControl.py GSE247998_SRR_ids_batch_3.txt PAIRED GSE247998_

"

# Obtaining counts for SRR26854217 (Terminal)

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR26854217_fastq
cd SRR26854217_fastq
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR268/017/SRR26854217/SRR26854217_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR268/017/SRR26854217/SRR26854217_2.fastq.gz

gunzip SRR26854217_1.fastq.gz
gunzip SRR26854217_2.fastq.gz

cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR26854217_STAR

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn SRR26854217_fastq/SRR26854217_1.fastq SRR26854217_fastq/SRR26854217_2.fastq \
--outFileNamePrefix SRR26854217_STAR/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"


# Inspecting runs
GSE247998_run_logs = list.files("GSE247998_OUTPUT", pattern = "logs", full.names = TRUE)
GSE247998_mapped_percent = sapply(GSE247998_run_logs, function(x){
  x = readLines(x)
  x = x[10]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads % |\t", replacement = "")
  x = stri_replace_all_fixed(str = x, pattern = "%", replacement = "")
  return(x)
})
GSE247998_mapped_percent = as.numeric(GSE247998_mapped_percent)
min(GSE247998_mapped_percent) # 77.52 %
max(GSE247998_mapped_percent) # 89.41 %
GSE247998_mapped_percent_df = data.frame(file = GSE247998_run_logs, percent = GSE247998_mapped_percent)
# all samples have relatively high percent of mapped reads


# preparing counts data
GSE247998_runs = list.files("GSE247998_OUTPUT", pattern = "counts", full.names = TRUE)
GSE247998_count_df = lapply(GSE247998_runs, function(x){
  x = smart_fread(x)
  return(x)
})
names(GSE247998_count_df) = GSE247998_runs
reference = GSE247998_count_df[[1]]$gene_ID
GSE247998_compar_to_ref = sapply(GSE247998_count_df, function(x){
  compar = reference == x$gene_ID
  compar = all(compar)
  return(compar)
})
all(GSE247998_compar_to_ref) # all row names match
all(reference == rownames(RNAseq_gene_probes$ID))

GSE247998_participant_1 = GSE247998_count_df[[1]]
GSE247998_participant_1 = GSE247998_participant_1[-(1:4), ]
GSE247998_participant_1$gene_ID = NULL
GSE247998_participant_1 = GSE247998_participant_1[rowSums(GSE247998_participant_1) > 0, ]
table(GSE247998_participant_1$Strand_1 > GSE247998_participant_1$Strand_2)
table(GSE247998_participant_1$Strand_1 < GSE247998_participant_1$Strand_2)
# Strand 2 is systematically bigger than 1
sum(GSE247998_participant_1$Strand_1) # 1206175
sum(GSE247998_participant_1$Strand_2) # 16040133
sum(GSE247998_participant_1$Strand_1)/sum(GSE247998_participant_1$Strand_2) # 0.07519732
sum(GSE247998_participant_1$Strand_2)/sum(GSE247998_participant_1$Strand_1) # 13.29835 



GSE247998_count = lapply(GSE247998_count_df, function(x){
  x = x[,"Strand_2"]
  return(x)
})
GSE247998_count = do.call(cbind,GSE247998_count)
rownames(GSE247998_count) = reference
colnames(GSE247998_count) = sapply(colnames(GSE247998_count), function(x){
  x = stri_replace_all_fixed(x, pattern = "GSE247998_OUTPUT/", replacement = "")
  x = stri_replace_all_fixed(x, pattern = "_star_counts.csv", replacement = "")
  return(x)
})
GSE247998_count = as.data.frame(GSE247998_count)

# check
all(colnames(GSE247998_count) == GSE247998_pheno_curated$RUN_ID) # FALSE
all(colnames(GSE247998_count) %in% GSE247998_pheno_curated$RUN_ID) # TRUE

GSE247998_count = GSE247998_count[,GSE247998_pheno_curated$RUN_ID]
all(colnames(GSE247998_count) == GSE247998_count$RUN_ID) # TRUE


# Last phenotype inspection before analysis
GSE247998_pheno_curated$AGE
GSE247998_pheno_curated$DIAGNOSIS
GSE247998_pheno_curated$SEX
GSE247998_pheno_curated$SUICIDE


# analysis without covariates

GSE247998_count_tmp = GSE247998_count[-(1:4), ]

# Design matrix
tmp_design = model.matrix(~ SUICIDE, data = GSE247998_pheno_curated)

# egeR filter
GSE247998_tmp_dge = DGEList(counts=GSE247998_count_tmp)
keep = filterByExpr(GSE247998_tmp_dge, tmp_design)
GSE247998_tmp_dge = GSE247998_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE247998_tmp_dge = calcNormFactors(GSE247998_tmp_dge)
GSE247998_tmp_voom = voom(GSE247998_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE247998_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE247998_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE247998_Top_table_no_covar_TMP$ID = rownames(GSE247998_Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE247998_Top_table_no_covar_TMP$ID]
GSE247998_Top_table_no_covar_TMP$SE = SE
all(names(SE) == GSE247998_Top_table_no_covar_TMP$ID) # TRUE

# Annotating results
GSE247998_probes_TMP = RNAseq_gene_probes
GSE247998_probes_TMP = GSE247998_probes_TMP[GSE247998_Top_table_no_covar_TMP$ID, ]
GSE247998_Top_table_no_covar_TMP$Gene_symbol = GSE247998_probes_TMP$Gene_symbol
GSE247998_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
GSE247998_Top_table_no_covar_TMP$Tissue = GSE247998_pheno_curated$TISSUE[1]
GSE247998_Top_table_no_covar_TMP$Tissue_type = "Blood"
GSE247998_Top_table_no_covar_TMP$Technology = "RNA-seq"


# analysis with covariates

GSE247998_count_tmp = GSE247998_count[-(1:4), ]

# Design matrix
tmp_design = model.matrix(~ SUICIDE + AGE + SEX + DIAGNOSIS, data = GSE247998_pheno_curated)

# egeR filter
GSE247998_tmp_dge = DGEList(counts=GSE247998_count_tmp)
keep = filterByExpr(GSE247998_tmp_dge, tmp_design)
GSE247998_tmp_dge = GSE247998_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE247998_tmp_dge = calcNormFactors(GSE247998_tmp_dge)
GSE247998_tmp_voom = voom(GSE247998_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE247998_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE247998_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE247998_Top_table_TMP$ID = rownames(GSE247998_Top_table_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE247998_Top_table_TMP$ID]
GSE247998_Top_table_TMP$SE = SE
all(names(SE) == GSE247998_Top_table_TMP$ID) # TRUE

# Annotating results
GSE247998_probes_TMP = RNAseq_gene_probes
GSE247998_probes_TMP = GSE247998_probes_TMP[GSE247998_Top_table_TMP$ID, ]
GSE247998_Top_table_TMP$Gene_symbol = GSE247998_probes_TMP$Gene_symbol
GSE247998_Top_table_TMP$Gene_symbol_non_hgnc = NA
GSE247998_Top_table_TMP$Tissue = GSE247998_pheno_curated$TISSUE[1]
GSE247998_Top_table_TMP$Tissue_type = "Blood"
GSE247998_Top_table_TMP$Technology = "RNA-seq"

# Saving
GSE247998_Top_table = GSE247998_Top_table_TMP
GSE247998_Top_table_no_covar = GSE247998_Top_table_no_covar_TMP

dir.create("GSE247998_results")
write.csv(GSE247998_count, "GSE247998_results/GSE247998_expression.csv")
write.csv(RNAseq_gene_probes, "GSE247998_results/GSE247998_probes.csv")
write.csv(GSE247998_pheno_curated, "GSE247998_results/GSE247998_pheno_curated.csv")
write.csv(GSE247998_Top_table_no_covar, "GSE247998_results/GSE247998_Top_table_no_covar.csv")
write.csv(GSE247998_Top_table, "GSE247998_results/GSE247998_Top_table.csv")

rm(list = ls(pattern = "GSE247998"))
gc()

################### GSE202537 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202537
# paper: https://pubmed.ncbi.nlm.nih.gov/36302706/

GSE202537_geo_record = getGEO("GSE202537")
GSE202537_geo_record = GSE202537_geo_record[[1]]
GSE202537_pheno = pData(GSE202537_geo_record)
GSE202537_pheno = fix_columns(GSE202537_pheno)


# curating pheno
GSE202537_pheno_curated = GSE202537_pheno
table(GSE202537_pheno_curated$ORGANISM_CH1) # 215 human

# pheno inspection 
GSE202537_pheno_curated$TISSUE # Caudate 72 Nac 71 Putamen 72
GSE202537_pheno_curated$AGE
GSE202537_pheno_curated$BMI
# time of death (TOD)
GSE202537_pheno_curated$CORRECTED_TOD
GSE202537_pheno_curated$DISEASE_STATE
table(GSE202537_pheno_curated$DISEASE_STATE) # control 108 bipolar 23 scz 84
GSE202537_pheno_curated$GENDER
GSE202537_pheno_curated$LIBRARY_SIZE # probably not needed as adjusted in limma
GSE202537_pheno_curated$MANNER_OF_DEATH
table(GSE202537_pheno_curated$MANNER_OF_DEATH) # Accidental 27 Natural 159 Suicide 23 Undetermined 6
GSE202537_pheno_curated$PH
GSE202537_pheno_curated$PMI
GSE202537_pheno_curated$RACE
GSE202537_pheno_curated$RIN
GSE202537_pheno_curated$TISSUESTORAGETIME
"
Brain specimens (n=72) were obtained, following consent from the next of kin, during autopsies 
conducted at the Allegheny County (Pittsburgh, PA; n=69) or the Davidson County (Nashville, TN; n=3) 
Medical Examiner’s Office. All procedures were approved by the University of Pittsburgh Institutional Review Board for 
Biomedical Research and Committee for Oversight of Research and Clinical Training Involving Decedents. An independent committee 
of experienced clinicians made consensus, lifetime DSM-IV diagnoses for each subject using the results of an expanded psychological 
autopsy, including structured interviews with family members and review of medical records, as well as toxicological and neuropathological 
reports (1). Unaffected comparison subjects had no known history of psychiatric or neurological disorders except for minor or in remission
psychiatric diagnoses in two subjects specified in the footnotes of Dataset S1. Subjects were included based on four criteria: 1) known time of 
death (TOD) within a 4-hour window; 2) age less than 65 years; 3) postmortem interval (PMI) less than 30 hours; and 4) died by accident, natural causes 
or suicide, suddenly and out of hospital, with no evidence of an agonal state. One subject (ID – 1180; BD with psychosis) had low expression data in the
NAc and therefore was removed from NAc-specific analyses. "


"Libraries were prepped for RNA-seq using the TruSeq Stranded Total RNA Sample Preparation Kit (Illumina). Paired-end dual-indexed sequencing (75 bp)
was performed using the NextSeq 500 platform (Illumina). Trimmomatic v0.38 (http://www.usadellab.org/cms/?page=trimmomatic) was used to filter poor
quality reads and trim poor-quality bases. FastQC v0.11.3 (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was performed to assess the quality 
of the data. Per base sequence quality was high (Quality score generally>30), indicating good data quality. HISAT2 (HISAT2v2.1.0) was used to align reads to 
the reference (Homo sapiens Ensembl GRCh38) using default parameters (3). The resulting bam files from HISAT2 were converted to expression count data using 
HTSeq (HTSeq v0.10.0, GTF file: Ensembl CRCh38) with default union mode (4). For the NAc and dorsal striatum sequencing runs (caudate and putamen samples run 
simultaneously), the average total paired end reads before alignment was approximately 47.0 (NAc), 47.0 (caudate), and 44.7 (putamen) million reads. The average 
mapped reads for the NAc, caudate, and putamen were 31.2, 34.4, and 32.8 million reads, respectively. RNA-seq count data were transformed to log2 continuous counts per
million (cpm) data using the cpm function of the Bioconductor edgeR package (5, 6). Transcripts were retained for analysis if log2(cpm) was greater than 1 in 50% or 
more of subjects. All Y-chromosome gene were eliminated from analysis. After filtering, 15,300 (NAc), 15,041 (caudate), and 14,866 (putamen) transcripts remained. Raw 
and processed RNA-seq data were deposited into the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) database under accession no. GSE202537."


# geting SRA runs
GSE202537_SRA_run_table = smart_fread("GSE202537_SraRunTable.txt")
all(GSE202537_pheno_curated$GEO_ACCESSION %in% GSE202537_SRA_run_table$`Sample Name`) # TRUE
all(GSE202537_SRA_run_table$`Sample Name` %in% GSE202537_pheno_curated$GEO_ACCESSION) # TRUE
GSE202537_pheno_curated$RUN_ID = sapply(GSE202537_pheno_curated$GEO_ACCESSION, function(x){
  x = GSE202537_SRA_run_table[GSE202537_SRA_run_table$`Sample Name` == x, "Run"]
  return(x)
})

GSE202537_pheno_curated$RUN_ID[1] # SRR19147648
writeLines(GSE202537_pheno_curated$RUN_ID, "GSE202537_SRR_ids.txt")

GSE202537_pheno_curated$TITLE[1] # "NAc, Caudate, and Putamen [10003N]"

# inspection for SRR19147648

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE202537_pref
source /etc/profile.d/sra-tools.sh
prefetch SRR19147648 -O GSE202537_pref

"
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE202537_pref
source /etc/profile.d/sra-tools.sh
fasterq-dump --threads 9 SRR19147648

"
# FastQC
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR19147648_QC
FastQC/fastqc GSE202537_pref/SRR19147648_1.fastq -o SRR19147648_QC/

"

# trimming with fastp
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir  GSE202537_trimmed
conda activate trimming
fastp --thread 8 -p --detect_adapter_for_pe -i  GSE202537_pref/SRR19147648_1.fastq -I GSE202537_pref/SRR19147648_2.fastq -o GSE202537_trimmed/out.R1.fq.gz -O GSE202537_trimmed/out.R2.fq.gz

"


# STAR (no trimming)  Uniquely mapped reads number |	38794275 Uniquely mapped reads % |	71.01%  Average mapped length |	149.18 % of reads mapped to multiple loci |	20.58%
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE202537_mapped

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE202537_pref/SRR19147648_1.fastq GSE202537_pref/SRR19147648_2.fastq \
--outFileNamePrefix GSE202537_mapped/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# STAR after trimming:  Uniquely mapped reads number |	38077139 Uniquely mapped reads % |	75.48%  % of reads mapped to multiple loci |	21.66%
# skipping trimming is OK as difference is minor

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE202537_mapped_trimmed

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE202537_trimmed/out.R1.fq.gz GSE202537_trimmed/out.R2.fq.gz --readFilesCommand gunzip -c \
--outFileNamePrefix GSE202537_mapped_trimmed/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# subreads for comparison score is twice as low compared to initial paper...
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE202537_counts

bin/featureCounts -T 10 -p --countReadPairs -a genome_files_hg19/Homo_sapiens.GRCh37.87.gtf \
-t gene -g gene_id -o GSE202537_counts/counts.txt GSE202537_mapped/Aligned.sortedByCoord.out.bam

"

# comparing scores for SRR19147648 (10003N or GSM6123790) using paper, SRAcounts, obtained counts

# SDF4, ENSG00000078808: paper 759, sra 446, we 401
# ACAP3 116983 ENSG00000131584  paper 2833, sra 1688, we 1529
# RER1 11079 ENSG00000157916 paper 1476, sra 792, we 781
# KIF1B 23095 ENSG00000054523 paper 29038, sra 16295, we 15257

# our scores are more consistent with SRA and are twice as low compared to paper...

# Our trimmed scores for the same genes: 390, 1483, 768, 14992 -> even lower
# Trimming does not improve things here -> probably SRA corrected files already
# We skip trimming in this data

# getting counts for all reads
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARControl.py GSE202537_SRR_ids.txt PAIRED GSE202537_

"
# SRR19147547 crashed on STAR?
# SRR19147547 inspection on https://trace.ncbi.nlm.nih.gov/Traces/?view=run_browser&page=100001&page_size=10&acc=SRR19147547&display=reads
# page 100001  in fact shows a lot of single reads? -> SRA corruption?
# ENA has the same issue and files for read are very small, whereas the bulk of SRR19147547 in unpaired...
# SRR19147547 should be removed

# SRR19147540 unpaired reads


# Rerun for remaining samples
GSE202537_runs_batch_1 = list.files("GSE202537_OUTPUT", pattern = "counts")
GSE202537_runs_batch_1 = stri_replace_all_fixed(GSE202537_runs_batch_1, pattern = "_star_counts.csv", replacement = "")
GSE202537_runs_not_finished = GSE202537_pheno_curated$RUN_ID[GSE202537_pheno_curated$RUN_ID %!in% GSE202537_runs_batch_1]
any(GSE202537_runs_not_finished %in% GSE202537_runs_batch_1) # FALSE
GSE202537_runs_not_finished = GSE202537_runs_not_finished[GSE202537_runs_not_finished != "SRR19147547"]
writeLines(GSE202537_runs_not_finished, "GSE202537_SRR_ids_batch_2.txt")

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARControl.py GSE202537_SRR_ids_batch_2.txt PAIRED GSE202537_

"
# Inspecting runs
GSE202537_run_logs = list.files("GSE202537_OUTPUT", pattern = "logs", full.names = TRUE)
GSE202537_mapped_percent = sapply(GSE202537_run_logs, function(x){
  x = readLines(x)
  x = x[10]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads % |\t", replacement = "")
  x = stri_replace_all_fixed(str = x, pattern = "%", replacement = "")
  return(x)
})
GSE202537_mapped_percent = as.numeric(GSE202537_mapped_percent)
min(GSE202537_mapped_percent) # 70.2 %
max(GSE202537_mapped_percent) # 92.48 %
GSE202537_mapped_percent_df = data.frame(file = GSE202537_run_logs, percent = GSE202537_mapped_percent)
# all samples have relatively high percent of mapped reads but some mapping may be corrupted due to many unpaired reads...


GSE202537_mapped_count = sapply(GSE202537_run_logs, function(x){
  x = readLines(x)
  x = x[9]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads number |\t", replacement = "")
  x = str_trim(x)
  return(x)
})
GSE202537_mapped_count = as.numeric(GSE202537_mapped_count)
min(GSE202537_mapped_count) # 10695974
max(GSE202537_mapped_count) # 60903495
GSE202537_mapped_percent_df$mapped_count = GSE202537_mapped_count
# SRR19147547 is removed
# SRR19147540 and SRR19147610 have only one third of counts
# as several reads in the cohort were unpaired for unknown reason -> suspect library contamination/scanner fail
# the samples are removed


# preparing counts data
GSE202537_runs = list.files("GSE202537_OUTPUT", pattern = "counts", full.names = TRUE)
GSE202537_count_df = lapply(GSE202537_runs, function(x){
  x = smart_fread(x)
  return(x)
})
names(GSE202537_count_df) = GSE202537_runs
reference = GSE202537_count_df[[1]]$gene_ID
GSE202537_compar_to_ref = sapply(GSE202537_count_df, function(x){
  compar = reference == x$gene_ID
  compar = all(compar)
  return(compar)
})
all(GSE202537_compar_to_ref) # all row names match
all(reference == rownames(RNAseq_gene_probes$ID))

GSE202537_participant_1 = GSE202537_count_df[[1]]
GSE202537_participant_1 = GSE202537_participant_1[-(1:4), ]
GSE202537_participant_1$gene_ID = NULL
GSE202537_participant_1 = GSE202537_participant_1[rowSums(GSE202537_participant_1) > 0, ]
table(GSE202537_participant_1$Strand_1 > GSE202537_participant_1$Strand_2) # TRUE 11111
table(GSE202537_participant_1$Strand_1 < GSE202537_participant_1$Strand_2) # TRUE 27170
# Strand 2 is systematically bigger than 1
sum(GSE202537_participant_1$Strand_1) # 1172968
sum(GSE202537_participant_1$Strand_2) # 20234154
sum(GSE202537_participant_1$Strand_1)/sum(GSE202537_participant_1$Strand_2) # 0.05796971
sum(GSE202537_participant_1$Strand_2)/sum(GSE202537_participant_1$Strand_1) # 17.25039 



GSE202537_count = lapply(GSE202537_count_df, function(x){
  x = x[,"Strand_2"]
  return(x)
})
GSE202537_count = do.call(cbind,GSE202537_count)
rownames(GSE202537_count) = reference
colnames(GSE202537_count) = sapply(colnames(GSE202537_count), function(x){
  x = stri_replace_all_fixed(x, pattern = "GSE202537_OUTPUT/", replacement = "")
  x = stri_replace_all_fixed(x, pattern = "_star_counts.csv", replacement = "")
  return(x)
})
GSE202537_count = as.data.frame(GSE202537_count)

# check
all(colnames(GSE202537_count) == GSE202537_pheno_curated$RUN_ID) # FALSE
all(colnames(GSE202537_count) %in% GSE202537_pheno_curated$RUN_ID) # TRUE

# removal of bad samples from pheno and counts
GSE202537_pheno_curated_2 = GSE202537_pheno_curated
GSE202537_pheno_curated_2 = GSE202537_pheno_curated_2[GSE202537_pheno_curated_2$RUN_ID %!in% c("SRR19147540","SRR19147610", "SRR19147547"), ]

GSE202537_count = GSE202537_count[,GSE202537_pheno_curated_2$RUN_ID]
all(colnames(GSE202537_count) == GSE202537_count$RUN_ID) # TRUE

# curating phenotypes
# pheno inspection 
GSE202537_pheno_curated_2$TISSUE 
table(GSE202537_pheno_curated_2$TISSUE) # Caudate 72 Nac 69 Putamen 71
GSE202537_pheno_curated_2$AGE
GSE202537_pheno_curated_2$BMI
# all missing BMIs are in controls 
# time of death (TOD)
GSE202537_pheno_curated_2$CORRECTED_TOD
GSE202537_pheno_curated_2$DISEASE_STATE
table(GSE202537_pheno_curated_2$DISEASE_STATE) # control 105 bipolar 23 scz 84
GSE202537_pheno_curated_2$GENDER
GSE202537_pheno_curated_2$LIBRARY_SIZE # probably not needed as adjusted in limma
GSE202537_pheno_curated_2$MANNER_OF_DEATH
table(GSE202537_pheno_curated_2$MANNER_OF_DEATH, GSE202537_pheno_curated_2$TISSUE) # Accidental 27 Natural 156 Suicide 23 Undetermined 6
"
               Caudate Nac Putamen
  Accidental         9   9       9
  Natural           53  51      52
  Suicide            8   7       8
  Undetermined       2   2       2
"
GSE202537_pheno_curated_2$PH
GSE202537_pheno_curated_2$PMI
GSE202537_pheno_curated_2$RACE
GSE202537_pheno_curated_2$RIN
GSE202537_pheno_curated_2$TISSUESTORAGETIME

# curation
GSE202537_pheno_curated_2$AGE = as.numeric(GSE202537_pheno_curated_2$AGE)
GSE202537_pheno_curated_2$BMI = as.numeric(GSE202537_pheno_curated_2$BMI) # produced NAs
GSE202537_pheno_curated_2$DISEASE_STATE = factor(GSE202537_pheno_curated_2$DISEASE_STATE, levels = c("match control",
                                                                                                     "psychosis_schizophrenia",
                                                                                                     "psychosis_bipolar"))
GSE202537_pheno_curated_2$GENDER = factor(GSE202537_pheno_curated_2$GENDER, levels = c("Female", "Male"))
GSE202537_pheno_curated_2$SUICIDE = ifelse(GSE202537_pheno_curated_2$MANNER_OF_DEATH == "Suicide", "Suicide", "Control")
GSE202537_pheno_curated_2$SUICIDE = factor(GSE202537_pheno_curated_2$SUICIDE, levels = c("Control", "Suicide"))
GSE202537_pheno_curated_2$PH = as.numeric(GSE202537_pheno_curated_2$PH)
GSE202537_pheno_curated_2$PMI = as.numeric(GSE202537_pheno_curated_2$PMI)
GSE202537_pheno_curated_2$RACE = factor(GSE202537_pheno_curated_2$RACE, levels = c("White", "Black"))
GSE202537_pheno_curated_2$RIN = as.numeric(GSE202537_pheno_curated_2$RIN)

# filtering of data
GSE202537_pheno_curated_2 = GSE202537_pheno_curated_2[GSE202537_pheno_curated_2$MANNER_OF_DEATH != "Undetermined",]
GSE202537_count = GSE202537_count[,GSE202537_pheno_curated_2$RUN_ID]
all(colnames(GSE202537_count) == GSE202537_count$RUN_ID) # TRUE



GSE202537_tissues = unique(GSE202537_pheno_curated_2$TISSUE)

GSE202537_analysis_list_no_covar = list()

for (x in 1:length(GSE202537_tissues)){
  
  
  GSE202537_count_tmp = GSE202537_count[-(1:4), ]
  GSE202537_TMP_tissue = GSE202537_tissues[x]
  
  print(paste0("Working on: ", GSE202537_TMP_tissue))
  
  GSE202537_pheno_curated_2_tmp = GSE202537_pheno_curated_2[GSE202537_pheno_curated_2$TISSUE == GSE202537_TMP_tissue,]
  GSE202537_count_tmp = GSE202537_count_tmp[,GSE202537_pheno_curated_2_tmp$RUN_ID]
  
  # Design matrix
  tmp_design = model.matrix(~ SUICIDE, data = GSE202537_pheno_curated_2_tmp)
  
  # egeR filter
  GSE202537_tmp_dge = DGEList(counts=GSE202537_count_tmp)
  keep = filterByExpr(GSE202537_tmp_dge, tmp_design)
  GSE202537_tmp_dge = GSE202537_tmp_dge[keep,,keep.lib.sizes=FALSE]
  GSE202537_tmp_dge = calcNormFactors(GSE202537_tmp_dge)
  GSE202537_tmp_voom = voom(GSE202537_tmp_dge, tmp_design, plot=TRUE)
  
  # limma 
  fit = lmFit(GSE202537_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  GSE202537_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE202537_Top_table_no_covar_TMP$ID = rownames(GSE202537_Top_table_no_covar_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE202537_Top_table_no_covar_TMP$ID]
  GSE202537_Top_table_no_covar_TMP$SE = SE
  all(names(SE) == GSE202537_Top_table_no_covar_TMP$ID) # TRUE
  
  # Annotating results
  GSE202537_probes_TMP = RNAseq_gene_probes
  GSE202537_probes_TMP = GSE202537_probes_TMP[GSE202537_Top_table_no_covar_TMP$ID, ]
  GSE202537_Top_table_no_covar_TMP$Gene_symbol = GSE202537_probes_TMP$Gene_symbol
  GSE202537_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
  GSE202537_Top_table_no_covar_TMP$Tissue = GSE202537_TMP_tissue
  GSE202537_Top_table_no_covar_TMP$Tissue_type = "Brain"
  GSE202537_Top_table_no_covar_TMP$Technology = "RNA-seq"
  
  GSE202537_analysis_list_no_covar[[x]] = GSE202537_Top_table_no_covar_TMP
  
}

GSE202537_Top_table_no_covar = do.call(rbind, GSE202537_analysis_list_no_covar)


# analysis with covariates in a loop
GSE202537_tissues = unique(GSE202537_pheno_curated_2$TISSUE)

# MANNER_OF_DEATH is not included as it is related to suicide group and Natural and Accidental are just subclasses of control

GSE202537_analysis_list = list()

for (x in 1:length(GSE202537_tissues)){
  
  GSE202537_count_tmp = GSE202537_count[-(1:4), ]
  GSE202537_TMP_tissue = GSE202537_tissues[x]
  
  print(paste0("Working on: ", GSE202537_TMP_tissue))
  
  GSE202537_pheno_curated_2_tmp = GSE202537_pheno_curated_2[GSE202537_pheno_curated_2$TISSUE == GSE202537_TMP_tissue,]
  GSE202537_count_tmp = GSE202537_count_tmp[,GSE202537_pheno_curated_2_tmp$RUN_ID]

  # Design matrix
  tmp_design = model.matrix(~ SUICIDE + AGE + BMI + DISEASE_STATE + GENDER + PH + PMI + RACE + RIN, data = GSE202537_pheno_curated_2_tmp)
  GSE202537_pheno_curated_2_tmp_small = GSE202537_pheno_curated_2_tmp[GSE202537_pheno_curated_2_tmp$GEO_ACCESSION %in% rownames(tmp_design),]
  GSE202537_count_tmp = GSE202537_count_tmp[,GSE202537_pheno_curated_2_tmp_small$RUN_ID]
  
  # egeR filter
  GSE202537_tmp_dge = DGEList(counts=GSE202537_count_tmp)
  keep = filterByExpr(GSE202537_tmp_dge, tmp_design)
  GSE202537_tmp_dge = GSE202537_tmp_dge[keep,,keep.lib.sizes=FALSE]
  GSE202537_tmp_dge = calcNormFactors(GSE202537_tmp_dge)
  GSE202537_tmp_voom = voom(GSE202537_tmp_dge, tmp_design, plot=TRUE)
  
  # limma 
  fit = lmFit(GSE202537_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  GSE202537_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE202537_Top_table_TMP$ID = rownames(GSE202537_Top_table_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE202537_Top_table_TMP$ID]
  GSE202537_Top_table_TMP$SE = SE
  all(names(SE) == GSE202537_Top_table_TMP$ID) # TRUE
  
  # Annotating results
  GSE202537_probes_TMP = RNAseq_gene_probes
  GSE202537_probes_TMP = GSE202537_probes_TMP[GSE202537_Top_table_TMP$ID, ]
  GSE202537_Top_table_TMP$Gene_symbol = GSE202537_probes_TMP$Gene_symbol
  GSE202537_Top_table_TMP$Gene_symbol_non_hgnc = NA
  GSE202537_Top_table_TMP$Tissue = GSE202537_TMP_tissue
  GSE202537_Top_table_TMP$Tissue_type = "Brain"
  GSE202537_Top_table_TMP$Technology = "RNA-seq"
  
  GSE202537_analysis_list[[x]] = GSE202537_Top_table_TMP
}

GSE202537_Top_table = do.call(rbind, GSE202537_analysis_list)
rownames(GSE202537_Top_table) = NULL

dir.create("GSE202537_results")
write.csv(GSE202537_count, "GSE202537_results/GSE202537_expression.csv")
write.csv(RNAseq_gene_probes, "GSE202537_results/GSE202537_probes.csv")
write.csv(GSE202537_pheno_curated_2, "GSE202537_results/GSE202537_pheno_curated.csv")
write.csv(GSE202537_Top_table_no_covar, "GSE202537_results/GSE202537_Top_table_no_covar.csv")
write.csv(GSE202537_Top_table, "GSE202537_results/GSE202537_Top_table.csv")

rm(list = ls(pattern = "GSE202537"))
gc()

################### GSE101521 ###################
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE101521
# paper: https://pubmed.ncbi.nlm.nih.gov/27528462/

GSE101521_geo_record = getGEO("GSE101521")
GSE101521_geo_record_1 = GSE101521_geo_record[[1]]
GSE101521_geo_record_2 = GSE101521_geo_record[[2]]

GSE101521_pheno_1 = pData(GSE101521_geo_record_1)
GSE101521_pheno_1 = fix_columns(GSE101521_pheno_1)

GSE101521_pheno_2 = pData(GSE101521_geo_record_2)
GSE101521_pheno_2 = fix_columns(GSE101521_pheno_2)


# curating pheno
table(GSE101521_pheno_1$ORGANISM_CH1) # 27 human
table(GSE101521_pheno_2$ORGANISM_CH1) # 59 human

table(GSE101521_pheno_1$DIAGNOSIS)
"
DSM-IV major depressive disorder non-suicides (MDD)   DSM-IV major depressive disorder suicides (MDD-S)                      non-psychiatric controls (CON) 
                                                  9                                                   9                                                   9 
"

table(GSE101521_pheno_2$DIAGNOSIS)
"
DSM-IV major depressive disorder non-suicides (MDD)   DSM-IV major depressive disorder suicides (MDD-S)                      non-psychiatric controls (CON) 
                                                  9                                                  21                                                  29 
"

# Pheno 2 makes sense to use

"
Fifty-nine clinical samples were obtained from the brain collection of The Division of Molecular Imaging 
and Neuropathology, at the New York State Psychiatric Institute and Columbia University. All procedures 
for brain collection and psychological autopsy were approved by the applicable Institutional Review Boards. 
Psychiatric diagnosis in the suicides and depressed individuals as well as absence of diagnoses in the 
controls was determined by the Structured Clinical Interview for DSM IV (SCID-I and II) as part of a 
validated psychological autopsy method described elsewhere (25). All subjects were selected because they
died suddenly to avoid metabolic complications related to agonal effects. All brains were free of gross 
neuropathology and had negative brain toxicology for psychotropic, illicit psychoactive drugs and
neurotoxic drugs. There were no diagnoses of Alcohol or Drug Use Disorders. Antemortem medication history
for three months ruled out recent exposure to psychotropic medication and confirmed results of comprehensive 
peripheral and brain toxicology. Brain samples were dissected from Brodmann Area 9 as previously reported (26) 
with an attempt to enrich for gray matter in when dissecting tissue from BA9. However, samples inevitably also
contained a small amount of white matter.

We applied whole-exome, gene-level analysis of count data using DESeq2 (27) to examine differential expression 
between 21 MDD-S (subjects with major depressive disorder and suicide), 9 MDD (subjects with MDD and no suicide), 
and 29 sudden death healthy control (CON) subjects with no MDD and no suicide (59 samples total). Small RNA 
analysis (miRNA differential expression) used an age and sex matched subset that included 9 MDD-S, 9 MDD and 9 
CON. The main results reported throughout the text derive from the full dataset of 59 samples, and adjust for 
age and gender. Unfortunately RIN scores from two samples (one MDD-S and one CON) were not estimated due to an
unexpected ribosomal ratio, low 18S or low RNA concentration. Given that RIN score captures the integrity of 
ribosomal RNAs but fails to measure total mRNA integrity, and is not necessarily a determining factor in generating 
good quality RNA sequencing data (personal communication with Peter Nagy, Director of Columbia University Genomics 
Core in the Department of Pathology), in the main text we present results in the full sample (N=59) without covarying 
for RIN score, and compared these results to an additional analysis that includes RIN score as a covariate in a subset 
of 57 samples (results listed in Supplement). See Supplementary Material for further discussion on accounting for the
possible effects of pH, PMI and RIN. The raw RNA-seq data will be deposited in publically available and structured 
repositories (i.e. Gene Expression Omnibus, see Supplementary Materials for more details).


"

"
Paired-end, strand specific sequencing for total RNA was performed on Illumina HiSeq 2500 with 100 bp read lengths, 
while single-end sequencing was performed for microRNA on Illumina MiSeq with 50 bp read lengths. For each clinical 
sample, raw RNA-seq reads were aligned and mapped to the Ensembl GRCh37 human reference genome and assembled using 
Tophat v2.0.9 resulting in BAM files for each of 59 samples. Between 17,000,000 and 57,000,000 reads were obtained 
for each sample, and ~75–90% of reads were successfully mapped to the genome for each sample. Of these reads, ~10–15% 
aligned to multiple genomic loci. Read statistics for each sample (subject) are listed in Supplementary Table 1. Note
that at the time of writing of this manuscript, the most recent GTF file compatible with Tophat was GRCh37. Hence 
some of the IDs listed in the main text tables have been deprecated and are therefore not included (i.e. they been 
replaced by one or more IDs) in the most recent (GRCh38) EnsEMBL database. The archived GRCh37 database can be 
accessed at http://grch37.ensembl.org/index.html.

"
GSE101521_pheno_curated = GSE101521_pheno_2

# variable curation
GSE101521_pheno_curated$AGE__YRS_
GSE101521_pheno_curated$BRAIN_PH
GSE101521_pheno_curated$DIAGNOSIS
GSE101521_pheno_curated$PMI
GSE101521_pheno_curated$RIN
GSE101521_pheno_curated$SEX


GSE101521_pheno_curated$AGE__YRS_ = as.numeric(GSE101521_pheno_curated$AGE__YRS_)
GSE101521_pheno_curated$BRAIN_PH = as.numeric(GSE101521_pheno_curated$BRAIN_PH)
GSE101521_pheno_curated$PMI = as.numeric(GSE101521_pheno_curated$PMI) # produced NAs
GSE101521_pheno_curated$RIN = as.numeric(GSE101521_pheno_curated$RIN) # produced NAs
GSE101521_pheno_curated$SEX = factor(GSE101521_pheno_curated$SEX, levels = c("Female", "Male"))
GSE101521_pheno_curated$SUICIDE = ifelse(GSE101521_pheno_curated$DIAGNOSIS == "DSM-IV major depressive disorder suicides (MDD-S)", "Suicide", "Non-suicide")
GSE101521_pheno_curated$SUICIDE = factor(GSE101521_pheno_curated$SUICIDE, levels = c("Non-suicide", "Suicide"))
GSE101521_pheno_curated$Depression = ifelse(stri_detect_fixed(GSE101521_pheno_curated$DIAGNOSIS, "MDD"), "MDD", "Control")
GSE101521_pheno_curated$Depression = factor(GSE101521_pheno_curated$Depression, levels = c("Control", "MDD"))

# geting SRA runs
GSE101521_SRA_run_table = smart_fread("GSE101521_SraRunTable.txt")
all(GSE101521_pheno_curated$GEO_ACCESSION %in% GSE101521_SRA_run_table$`Sample Name`) # TRUE
all(GSE101521_SRA_run_table$`Sample Name` %in% GSE101521_pheno_curated$GEO_ACCESSION) # FALSE as dataset does not include small RNAs
GSE101521_pheno_curated$RUN_ID = sapply(GSE101521_pheno_curated$GEO_ACCESSION, function(x){
  x = GSE101521_SRA_run_table[GSE101521_SRA_run_table$`Sample Name` == x, "Run"]
  return(x)
})

GSE101521_pheno_curated$RUN_ID[1] # SRR5831944
writeLines(GSE101521_pheno_curated$RUN_ID, "GSE101521_SRR_ids.txt")


"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir GSE101521_pref
source /etc/profile.d/sra-tools.sh
prefetch SRR5831944 -O GSE101521_pref

"

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE101521_pref
source /etc/profile.d/sra-tools.sh
fasterq-dump --threads 9 SRR5831944

"
# FastQC
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir SRR5831944_QC
FastQC/fastqc GSE101521_pref/SRR5831944_1.fastq -o SRR5831944_QC/

"

# trimming with fastp
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
mkdir  GSE101521_trimmed
conda activate trimming
fastp --thread 8 -p --detect_adapter_for_pe -i  GSE101521_pref/SRR5831944_1.fastq -I GSE101521_pref/SRR5831944_2.fastq -o GSE101521_trimmed/out.R1.fq.gz -O GSE101521_trimmed/out.R2.fq.gz

"


# STAR (no trimming) #  Uniquely mapped reads number |	15462942 # Uniquely mapped reads % |	90.11% # Number of reads mapped to multiple loci |	1001185
# % of reads mapped to multiple loci |	5.83%
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE101521_mapped

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE101521_pref/SRR5831944_1.fastq GSE101521_pref/SRR5831944_2.fastq \
--outFileNamePrefix GSE101521_mapped/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"

# STAR after trimming Uniquely mapped reads number |	15549603 Uniquely mapped reads % |	93.40% Number of reads mapped to multiple loci |	984834
# % of reads mapped to multiple loci |	5.92%
# quality seems to slightly improve 
# we keep without trimming for consistency as difference is not large

"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis

mkdir GSE101521_mapped_trimmed

STAR --runMode alignReads --runThreadN 10 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
--readFilesIn GSE101521_trimmed/out.R1.fq.gz GSE101521_trimmed/out.R2.fq.gz --readFilesCommand gunzip -c \
--outFileNamePrefix GSE101521_mapped_trimmed/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf

"
# Comparing counts with paper:

# V1 is SRR5831944	

# ENSG00000000003 paper 77.6591310902939 we 55
# ENSG00000004478 963.271914485376 we 884
# ENSG00000005436 428.618665825276 we 190
# ENSG00000007237 6759.33129451289 we 4843

# counts are different due to methods
"
For each BAM file resulting from Tophat2, read counts per gene were summarized using the ‘summarizeOverlaps’
function in the GenomicAlignments R library (29) and a transcript database derived from the GRCh37 human genome assembly.

"

# obtaining counts
"
conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARControl.py GSE101521_SRR_ids.txt PAIRED GSE101521_

"

# Inspecting runs
GSE101521_run_logs = list.files("GSE101521_OUTPUT", pattern = "logs", full.names = TRUE)
GSE101521_mapped_percent = sapply(GSE101521_run_logs, function(x){
  x = readLines(x)
  x = x[10]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads % |\t", replacement = "")
  x = stri_replace_all_fixed(str = x, pattern = "%", replacement = "")
  return(x)
})
GSE101521_mapped_percent = as.numeric(GSE101521_mapped_percent)
min(GSE101521_mapped_percent) # 62.39 %
max(GSE101521_mapped_percent) # 90.11 %
GSE101521_mapped_percent_df = data.frame(file = GSE101521_run_logs, percent = GSE101521_mapped_percent)
# SRR5831982 has 32% mapped to multiple loci
# SRR5831976 has 23.6% mapped to multiple loci
# SRR5831976 has 23.6% mapped to multiple loci
# SRR5831990 has 24.65% mapped to multiple loci


GSE101521_mapped_count = sapply(GSE101521_run_logs, function(x){
  x = readLines(x)
  x = x[9]
  x = stri_replace_all_fixed(str = x, pattern = "Uniquely mapped reads number |\t", replacement = "")
  x = str_trim(x)
  return(x)
})
GSE101521_mapped_count = as.numeric(GSE101521_mapped_count)
min(GSE101521_mapped_count) # 6721160
max(GSE101521_mapped_count) # 52841934
GSE101521_mapped_percent_df$mapped_count = GSE101521_mapped_count
6721160/mean(GSE101521_mapped_percent_df$mapped_count) # 0.255545
# /SRR5831961 has low count but this is due to overall low reads number
# since all reads in all samples were paired and the % of mapped reads is high 
# compared to initial reads number -> samples are kept

# preparing counts data
GSE101521_runs = list.files("GSE101521_OUTPUT", pattern = "counts", full.names = TRUE)
GSE101521_count_df = lapply(GSE101521_runs, function(x){
  x = smart_fread(x)
  return(x)
})
names(GSE101521_count_df) = GSE101521_runs
reference = GSE101521_count_df[[1]]$gene_ID
GSE101521_compar_to_ref = sapply(GSE101521_count_df, function(x){
  compar = reference == x$gene_ID
  compar = all(compar)
  return(compar)
})
all(GSE101521_compar_to_ref) # all row names match
all(reference == rownames(RNAseq_gene_probes$ID))

GSE101521_participant_1 = GSE101521_count_df[[1]]
GSE101521_participant_1 = GSE101521_participant_1[-(1:4), ]
GSE101521_participant_1$gene_ID = NULL
GSE101521_participant_1 = GSE101521_participant_1[rowSums(GSE101521_participant_1) > 0, ]
table(GSE101521_participant_1$Strand_1 > GSE101521_participant_1$Strand_2) # TRUE 9136
table(GSE101521_participant_1$Strand_1 < GSE101521_participant_1$Strand_2) # TRUE 24465
# Strand 2 is systematically bigger than 1
sum(GSE101521_participant_1$Strand_1) # 360847
sum(GSE101521_participant_1$Strand_2) # 7924989
sum(GSE101521_participant_1$Strand_1)/sum(GSE101521_participant_1$Strand_2) # 0.04553281
sum(GSE101521_participant_1$Strand_2)/sum(GSE101521_participant_1$Strand_1) # 21.96219 


GSE101521_count = lapply(GSE101521_count_df, function(x){
  x = x[,"Strand_2"]
  return(x)
})
GSE101521_count = do.call(cbind,GSE101521_count)
rownames(GSE101521_count) = reference
colnames(GSE101521_count) = sapply(colnames(GSE101521_count), function(x){
  x = stri_replace_all_fixed(x, pattern = "GSE101521_OUTPUT/", replacement = "")
  x = stri_replace_all_fixed(x, pattern = "_star_counts.csv", replacement = "")
  return(x)
})
GSE101521_count = as.data.frame(GSE101521_count)

# check
all(colnames(GSE101521_count) == GSE101521_pheno_curated$RUN_ID) # TRUE
all(colnames(GSE101521_count) %in% GSE101521_pheno_curated$RUN_ID) # TRUE


# variable curation
GSE101521_pheno_curated$AGE__YRS_
GSE101521_pheno_curated$BRAIN_PH
GSE101521_pheno_curated$SUICIDE
GSE101521_pheno_curated$PMI
GSE101521_pheno_curated$RIN
GSE101521_pheno_curated$SEX
GSE101521_pheno_curated$TISSUE
GSE101521_pheno_curated$Depression

# analysis without covariates

GSE101521_count_tmp = GSE101521_count[-(1:4), ]

# Design matrix
tmp_design = model.matrix(~ SUICIDE, data = GSE101521_pheno_curated)

# egeR filter
GSE101521_tmp_dge = DGEList(counts=GSE101521_count_tmp)
keep = filterByExpr(GSE101521_tmp_dge, tmp_design)
GSE101521_tmp_dge = GSE101521_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE101521_tmp_dge = calcNormFactors(GSE101521_tmp_dge)
GSE101521_tmp_voom = voom(GSE101521_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE101521_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE101521_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE101521_Top_table_no_covar_TMP$ID = rownames(GSE101521_Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE101521_Top_table_no_covar_TMP$ID]
GSE101521_Top_table_no_covar_TMP$SE = SE
all(names(SE) == GSE101521_Top_table_no_covar_TMP$ID) # TRUE

# Annotating results
GSE101521_probes_TMP = RNAseq_gene_probes
GSE101521_probes_TMP = GSE101521_probes_TMP[GSE101521_Top_table_no_covar_TMP$ID, ]
GSE101521_Top_table_no_covar_TMP$Gene_symbol = GSE101521_probes_TMP$Gene_symbol
GSE101521_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
GSE101521_Top_table_no_covar_TMP$Tissue = GSE101521_pheno_curated$TISSUE[1]
GSE101521_Top_table_no_covar_TMP$Tissue_type = "Brain"
GSE101521_Top_table_no_covar_TMP$Technology = "RNA-seq"


# analysis with covariates

GSE101521_count_tmp = GSE101521_count[-(1:4), ]

# Design matrix
GSE101521_pheno_curated_tmp = GSE101521_pheno_curated
rownames(GSE101521_pheno_curated_tmp) = GSE101521_pheno_curated_tmp$RUN_ID
tmp_design = model.matrix(~ SUICIDE + Depression + AGE__YRS_ + SEX + BRAIN_PH + PMI + RIN, data = GSE101521_pheno_curated_tmp)
GSE101521_count_tmp = GSE101521_count_tmp[, rownames(tmp_design)]

# egeR filter
GSE101521_tmp_dge = DGEList(counts=GSE101521_count_tmp)
keep = filterByExpr(GSE101521_tmp_dge, tmp_design)
GSE101521_tmp_dge = GSE101521_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE101521_tmp_dge = calcNormFactors(GSE101521_tmp_dge)
GSE101521_tmp_voom = voom(GSE101521_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE101521_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE101521_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE101521_Top_table_TMP$ID = rownames(GSE101521_Top_table_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE101521_Top_table_TMP$ID]
GSE101521_Top_table_TMP$SE = SE
all(names(SE) == GSE101521_Top_table_TMP$ID) # TRUE

# Annotating results
GSE101521_probes_TMP = RNAseq_gene_probes
GSE101521_probes_TMP = GSE101521_probes_TMP[GSE101521_Top_table_TMP$ID, ]
GSE101521_Top_table_TMP$Gene_symbol = GSE101521_probes_TMP$Gene_symbol
GSE101521_Top_table_TMP$Gene_symbol_non_hgnc = NA
GSE101521_Top_table_TMP$Tissue = GSE101521_pheno_curated_tmp$TISSUE[1]
GSE101521_Top_table_TMP$Tissue_type = "Brain"
GSE101521_Top_table_TMP$Technology = "RNA-seq"

# Saving
GSE101521_Top_table = GSE101521_Top_table_TMP
GSE101521_Top_table_no_covar = GSE101521_Top_table_no_covar_TMP

dir.create("GSE101521_results")
write.csv(GSE101521_count, "GSE101521_results/GSE101521_expression.csv")
write.csv(RNAseq_gene_probes, "GSE101521_results/GSE101521_probes.csv")
write.csv(GSE101521_pheno_curated, "GSE101521_results/GSE101521_pheno_curated.csv")
write.csv(GSE101521_Top_table_no_covar, "GSE101521_results/GSE101521_Top_table_no_covar.csv")
write.csv(GSE101521_Top_table, "GSE101521_results/GSE101521_Top_table.csv")

rm(list = ls(pattern = "GSE101521"))
gc()

################### Single-cell RNA Seq cohorts ###################

################### GSE144136 and GSE213982 ###################
# GSE144136 Male; GSE213982 Female
GSE213982_geo_record = getGEO("GSE213982")
GSE213982_geo_record = GSE213982_geo_record[[1]]
GSE213982_pheno = pData(GSE213982_geo_record)
GSE213982_pheno = fix_columns(GSE213982_pheno)


GSE144136_geo_record = getGEO("GSE144136")
GSE144136_geo_record = GSE144136_geo_record[[1]]
GSE144136_pheno = pData(GSE144136_geo_record)
GSE144136_pheno = fix_columns(GSE144136_pheno)

# supplementary files
getGEOSuppFiles("GSE213982")
files = list.files("GSE213982", full.names = TRUE)
sapply(files, gunzip)

# inspecting files
gene_rows = read.csv("GSE213982/GSE213982_combined_counts_matrix_genes_rows.csv")
gene_rows = gene_rows$x
length(gene_rows) # 36588
any(duplicated(gene_rows)) # FALSE

cells_colums =  read.csv("GSE213982/GSE213982_combined_counts_matrix_cells_columns.csv")
cells_colums = cells_colums$x
any(duplicated(cells_colums)) # FALSE

# 79058 + 81653 = 160711

"After filtering, we obtained 79,058 nuclei in the male cohort (43,347 from cases, 35,711 from controls) 
and 81,653 nuclei in the female cohort (49,926 from cases, 31,727 from controls). In the female cohort, 
after filtering, the median across samples of the median number of UMIs per cell and the median the number
of genes per cell were 2758.5 and 1711.5 respectively (Supplementary Data 1). In the males, the corresponding 
numbers were 2530.5 and 1638.25 respectively (Supplementary Data 1)."

# https://www.nature.com/articles/nmeth.4612 limma voom could be used for SC-RNA
# https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# MALE and FEMALE sample are obtained with different platforms
# GSE144136 Male Illumina HiSeq 4000 (Homo sapiens)
# GSE213982 Female Illumina NovaSeq 6000 (Homo sapiens
# Need to be analyzed separately

# phenotype curation
GSE213982_pheno$GROUP = factor(GSE213982_pheno$GROUP, levels = c("Control", "Case"))
GSE144136_pheno$GROUP = factor(GSE144136_pheno$GROUP, levels = c("Control", "Major Depressive Disorder (MDD)"))

GSE213982_pheno$PARTICIPANT = GSE213982_pheno$TITLE
GSE144136_pheno$PARTICIPANT = sapply(GSE144136_pheno$TITLE, function(x){
  x = unlist(stri_split_fixed(x, pattern = ":"))
  x = x[1]
  x = paste0("M", x)
  return(x)
})

GSE144136_GSE213982_count_matrix = readMM("GSE213982/GSE213982_combined_counts_matrix.mtx")
rownames(GSE144136_GSE213982_count_matrix) = gene_rows
colnames(GSE144136_GSE213982_count_matrix) = cells_colums
GSE144136_GSE213982_seurat = CreateSeuratObject(counts = GSE144136_GSE213982_count_matrix)


GSE144136_GSE213982_metadata = data.frame("sample_name" = colnames(GSE144136_GSE213982_seurat))

GSE144136_GSE213982_metadata$person = sapply(GSE144136_GSE213982_metadata$sample_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[1]
  return(x)
})


GSE144136_GSE213982_metadata$barcode = sapply(GSE144136_GSE213982_metadata$sample_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[2]
  return(x)
})

GSE144136_GSE213982_metadata$cell_overall = sapply(GSE144136_GSE213982_metadata$sample_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[3]
  return(x)
})

GSE144136_GSE213982_metadata$cell_cluster = sapply(GSE144136_GSE213982_metadata$sample_name, function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[4]
  return(x)
})

GSE144136_GSE213982_seurat = AddMetaData(object = GSE144136_GSE213982_seurat, metadata =  GSE144136_GSE213982_metadata)

GSE144136_GSE213982_cell_proprotions = table(GSE144136_GSE213982_metadata$person, GSE144136_GSE213982_metadata$cell_overall)
GSE144136_GSE213982_cell_proprotions = as.data.frame.matrix(GSE144136_GSE213982_cell_proprotions)
GSE144136_GSE213982_cell_proprotions = apply(GSE144136_GSE213982_cell_proprotions, 1, function(x){
  x = x/sum(x)
  return(x)
}, simplify = FALSE)
GSE144136_GSE213982_cell_proprotions = do.call(rbind, GSE144136_GSE213982_cell_proprotions)
GSE144136_GSE213982_cell_proprotions = as.data.frame(GSE144136_GSE213982_cell_proprotions)
GSE144136_GSE213982_cell_proprotions$PARTICIPANT = rownames(GSE144136_GSE213982_cell_proprotions)

apply(GSE144136_GSE213982_cell_proprotions[,1:8],2, mean)
"
Ast        End        ExN        InN        Mic        Mix        Oli        OPC 
0.08224416 0.02845957 0.48762526 0.19154091 0.02146172 0.02444783 0.11564245 0.04857810 
"
# MIC is the smallest on average


###
###
###
###
### summarization of counts
GSE144136_GSE213982_counts_mixed = AggregateExpression(GSE144136_GSE213982_seurat, group.by = c("person"), return.seurat = FALSE)
GSE144136_GSE213982_counts_mixed = GSE144136_GSE213982_counts_mixed$RNA
GSE144136_GSE213982_counts_mixed = as.matrix(GSE144136_GSE213982_counts_mixed)

GSE144136_GSE213982_counts_cell_type = AggregateExpression(GSE144136_GSE213982_seurat, group.by = c("person", "cell_overall"), return.seurat = FALSE)
GSE144136_GSE213982_counts_cell_type = GSE144136_GSE213982_counts_cell_type$RNA
GSE144136_GSE213982_counts_cell_type = as.matrix(GSE144136_GSE213982_counts_cell_type)

GSE144136_GSE213982_counts_cell_type_list = list()
cell_types = unique(GSE144136_GSE213982_metadata$cell_overall)

for (i in 1:length(cell_types)){
  current_cell = cell_types[i]
  GSE144136_GSE213982_current_counts = GSE144136_GSE213982_counts_cell_type[, stri_detect_fixed(colnames(GSE144136_GSE213982_counts_cell_type), pattern = current_cell)]
  colnames(GSE144136_GSE213982_current_counts) = stri_replace_all_fixed(colnames(GSE144136_GSE213982_current_counts), 
                                                                        pattern = paste0("_", current_cell), 
                                                                        replacement = "")
  GSE144136_GSE213982_counts_cell_type_list[[i]] = GSE144136_GSE213982_current_counts
}
names(GSE144136_GSE213982_counts_cell_type_list) = cell_types

GSE144136_GSE213982_probes = data.frame(init_genes = rownames(GSE144136_GSE213982_counts_mixed))

GSE144136_GSE213982_probes_check = check_gene_symbol_NIH(PRF_gene_symbols = GSE144136_GSE213982_probes$init_genes, 
                                                       PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                                       PRF_replace_NA_with_old = TRUE)
GSE144136_GSE213982_probes_check_dict = GSE144136_GSE213982_probes_check$Suggested.Symbol
names(GSE144136_GSE213982_probes_check_dict) = GSE144136_GSE213982_probes_check$x

GSE144136_GSE213982_probes$Gene_symbol = sapply(GSE144136_GSE213982_probes$init_genes, function(x){
  gene = GSE144136_GSE213982_probes_check_dict[x]
  return(gene)
})
rownames(GSE144136_GSE213982_probes) = GSE144136_GSE213982_probes$init_genes

################### GSE144136 ###################

GSE144136_pheno$TISSUE = "Dorsolateral prefrontal cortex (BA9)"
GSE144136_pheno$GROUP
GSE144136_pheno$PARTICIPANT
GSE144136_pheno$SEX


GSE144136_selection = base::intersect(GSE144136_pheno$PARTICIPANT, colnames(GSE144136_GSE213982_counts_mixed))
GSE144136_selection = GSE144136_selection[GSE144136_selection != "M24"]

GSE144136_pheno_curated = GSE144136_pheno
rownames(GSE144136_pheno_curated) = GSE144136_pheno_curated$PARTICIPANT
GSE144136_pheno_curated = GSE144136_pheno_curated[GSE144136_selection,]
GSE144136_cell_proprotions = GSE144136_GSE213982_cell_proprotions[GSE144136_selection,]
all(rownames(GSE144136_cell_proprotions) == GSE144136_cell_proprotions$PARTICIPANT) # TRUE
all(rownames(GSE144136_cell_proprotions) == GSE144136_pheno_curated$PARTICIPANT) # TRUE
GSE144136_cell_proprotions$PARTICIPANT = NULL
GSE144136_pheno_curated = cbind(GSE144136_pheno_curated, GSE144136_cell_proprotions)

# analysis without covariates (male)
GSE144136_count_tmp = GSE144136_GSE213982_counts_mixed[,GSE144136_selection]

# Design matrix
tmp_design = model.matrix(~ GROUP, data = GSE144136_pheno_curated)

# egeR filter
GSE144136_tmp_dge = DGEList(counts=GSE144136_count_tmp)
keep = filterByExpr(GSE144136_tmp_dge, tmp_design)
GSE144136_tmp_dge = GSE144136_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE144136_tmp_dge = calcNormFactors(GSE144136_tmp_dge)
GSE144136_tmp_voom = voom(GSE144136_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE144136_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE144136_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE144136_Top_table_no_covar_TMP$ID = rownames(GSE144136_Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE144136_Top_table_no_covar_TMP$ID]
GSE144136_Top_table_no_covar_TMP$SE = SE
all(names(SE) == GSE144136_Top_table_no_covar_TMP$ID) # TRUE

# Annotating results
GSE144136_probes_TMP = GSE144136_GSE213982_probes
GSE144136_probes_TMP = GSE144136_probes_TMP[GSE144136_Top_table_no_covar_TMP$ID, ]
GSE144136_Top_table_no_covar_TMP$Gene_symbol = GSE144136_probes_TMP$Gene_symbol
GSE144136_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
GSE144136_Top_table_no_covar_TMP$Tissue = GSE144136_pheno_curated$TISSUE[1]
GSE144136_Top_table_no_covar_TMP$Tissue_type = "Brain"
GSE144136_Top_table_no_covar_TMP$Technology = "scRNA-seq"


# analysis with covariates (male)
GSE144136_count_tmp = GSE144136_GSE213982_counts_mixed[,GSE144136_selection]

# Design matrix
tmp_design = model.matrix(~ GROUP + Ast + End + ExN + InN + Mix + Oli + OPC, data = GSE144136_pheno_curated)

# egeR filter
GSE144136_tmp_dge = DGEList(counts=GSE144136_count_tmp)
keep = filterByExpr(GSE144136_tmp_dge, tmp_design)
GSE144136_tmp_dge = GSE144136_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE144136_tmp_dge = calcNormFactors(GSE144136_tmp_dge)
GSE144136_tmp_voom = voom(GSE144136_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE144136_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE144136_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE144136_Top_table_TMP$ID = rownames(GSE144136_Top_table_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE144136_Top_table_TMP$ID]
GSE144136_Top_table_TMP$SE = SE
all(names(SE) == GSE144136_Top_table_TMP$ID) # TRUE

# Annotating results
GSE144136_probes_TMP = GSE144136_GSE213982_probes
GSE144136_probes_TMP = GSE144136_probes_TMP[GSE144136_Top_table_TMP$ID, ]
GSE144136_Top_table_TMP$Gene_symbol = GSE144136_probes_TMP$Gene_symbol
GSE144136_Top_table_TMP$Gene_symbol_non_hgnc = NA
GSE144136_Top_table_TMP$Tissue = GSE144136_pheno_curated$TISSUE[1]
GSE144136_Top_table_TMP$Tissue_type = "Brain"
GSE144136_Top_table_TMP$Technology = "scRNA-seq"

# analysis without covariates per cell
GSE144136_Top_table_no_covar_cell_types = list()

for (i in 1:length(cell_types)){
  
  current_cell = cell_types[i]
  GSE144136_count_tmp = GSE144136_GSE213982_counts_cell_type_list[current_cell][[1]]
  GSE144136_count_tmp = GSE144136_count_tmp[,GSE144136_selection]
  
  # Design matrix
  tmp_design = model.matrix(~ GROUP, data = GSE144136_pheno_curated)
  
  # egeR filter
  GSE144136_tmp_dge = DGEList(counts=GSE144136_count_tmp)
  keep = filterByExpr(GSE144136_tmp_dge, tmp_design)
  GSE144136_tmp_dge = GSE144136_tmp_dge[keep,,keep.lib.sizes=FALSE]
  GSE144136_tmp_dge = calcNormFactors(GSE144136_tmp_dge)
  GSE144136_tmp_voom = voom(GSE144136_tmp_dge, tmp_design, plot=TRUE)
  
  # limma 
  fit = lmFit(GSE144136_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  GSE144136_Top_table_no_covar__cell_type_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE144136_Top_table_no_covar__cell_type_TMP$ID = rownames(GSE144136_Top_table_no_covar__cell_type_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE144136_Top_table_no_covar__cell_type_TMP$ID]
  GSE144136_Top_table_no_covar__cell_type_TMP$SE = SE
  all(names(SE) == GSE144136_Top_table_no_covar__cell_type_TMP$ID) # TRUE
  
  # Annotating results
  GSE144136_probes_TMP = GSE144136_GSE213982_probes
  GSE144136_probes_TMP = GSE144136_probes_TMP[GSE144136_Top_table_no_covar__cell_type_TMP$ID, ]
  GSE144136_Top_table_no_covar__cell_type_TMP$Gene_symbol = GSE144136_probes_TMP$Gene_symbol
  GSE144136_Top_table_no_covar__cell_type_TMP$Gene_symbol_non_hgnc = NA
  GSE144136_Top_table_no_covar__cell_type_TMP$Tissue = paste0(GSE144136_pheno_curated$TISSUE[1], ": ", current_cell)
  
  GSE144136_Top_table_no_covar__cell_type_TMP$Tissue_type = "Brain"
  GSE144136_Top_table_no_covar__cell_type_TMP$Technology = "scRNA-seq"
  
  GSE144136_Top_table_no_covar_cell_types[[i]] = GSE144136_Top_table_no_covar__cell_type_TMP
  
}
GSE144136_Top_table_no_covar_cell_types = do.call(rbind, GSE144136_Top_table_no_covar_cell_types)


dir.create("GSE144136_results")
write.csv(GSE144136_GSE213982_counts_mixed, "GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
write.csv(GSE144136_GSE213982_counts_cell_type, "GSE144136_results/GSE144136_GSE213982_counts_cell_type.csv")
write.csv(GSE144136_pheno_curated, "GSE144136_results/GSE144136_pheno_curated.csv")
write.csv(GSE144136_Top_table_no_covar_TMP, "GSE144136_results/GSE144136_Top_table_no_covar.csv")
write.csv(GSE144136_Top_table_TMP, "GSE144136_results/GSE144136_Top_table.csv")
write.csv(GSE144136_Top_table_no_covar_cell_types, "GSE144136_results/GSE144136_Top_table_no_covar_cell_types.csv")

################### GSE213982 ###################

GSE213982_pheno$TISSUE = "Dorsolateral prefrontal cortex (BA9)"
GSE213982_pheno$GROUP
GSE213982_pheno$PARTICIPANT
GSE213982_pheno$SEX

GSE213982_selection = base::intersect(GSE213982_pheno$PARTICIPANT, colnames(GSE144136_GSE213982_counts_mixed))

GSE213982_pheno_curated = GSE213982_pheno
rownames(GSE213982_pheno_curated) = GSE213982_pheno_curated$PARTICIPANT
GSE213982_pheno_curated = GSE213982_pheno_curated[GSE213982_selection,]
GSE213982_cell_proprotions = GSE144136_GSE213982_cell_proprotions[GSE213982_selection,]
all(rownames(GSE213982_cell_proprotions) == GSE213982_cell_proprotions$PARTICIPANT) # TRUE
all(rownames(GSE213982_cell_proprotions) == GSE213982_pheno_curated$PARTICIPANT) # TRUE
GSE213982_cell_proprotions$PARTICIPANT = NULL
GSE213982_pheno_curated = cbind(GSE213982_pheno_curated, GSE213982_cell_proprotions)

# analysis without covariates (female)
GSE213982_count_tmp = GSE144136_GSE213982_counts_mixed[,GSE213982_selection]

# Design matrix
tmp_design = model.matrix(~ GROUP, data = GSE213982_pheno_curated)

# egeR filter
GSE213982_tmp_dge = DGEList(counts=GSE213982_count_tmp)
keep = filterByExpr(GSE213982_tmp_dge, tmp_design)
GSE213982_tmp_dge = GSE213982_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE213982_tmp_dge = calcNormFactors(GSE213982_tmp_dge)
GSE213982_tmp_voom = voom(GSE213982_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE213982_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE213982_Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE213982_Top_table_no_covar_TMP$ID = rownames(GSE213982_Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE213982_Top_table_no_covar_TMP$ID]
GSE213982_Top_table_no_covar_TMP$SE = SE
all(names(SE) == GSE213982_Top_table_no_covar_TMP$ID) # TRUE

# Annotating results
GSE213982_probes_TMP = GSE144136_GSE213982_probes
GSE213982_probes_TMP = GSE213982_probes_TMP[GSE213982_Top_table_no_covar_TMP$ID, ]
GSE213982_Top_table_no_covar_TMP$Gene_symbol = GSE213982_probes_TMP$Gene_symbol
GSE213982_Top_table_no_covar_TMP$Gene_symbol_non_hgnc = NA
GSE213982_Top_table_no_covar_TMP$Tissue = GSE213982_pheno_curated$TISSUE[1]
GSE213982_Top_table_no_covar_TMP$Tissue_type = "Brain"
GSE213982_Top_table_no_covar_TMP$Technology = "scRNA-seq"


# analysis with covariates (female)
GSE213982_count_tmp = GSE144136_GSE213982_counts_mixed[,GSE213982_selection]

# Design matrix
tmp_design = model.matrix(~ GROUP + Ast + End + ExN + InN + Mix + Oli + OPC, data = GSE213982_pheno_curated)

# egeR filter
GSE213982_tmp_dge = DGEList(counts=GSE213982_count_tmp)
keep = filterByExpr(GSE213982_tmp_dge, tmp_design)
GSE213982_tmp_dge = GSE213982_tmp_dge[keep,,keep.lib.sizes=FALSE]
GSE213982_tmp_dge = calcNormFactors(GSE213982_tmp_dge)
GSE213982_tmp_voom = voom(GSE213982_tmp_dge, tmp_design, plot=TRUE)

# limma 
fit = lmFit(GSE213982_tmp_voom, tmp_design)
fitE = eBayes(fit)
GSE213982_Top_table_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
GSE213982_Top_table_TMP$ID = rownames(GSE213982_Top_table_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[GSE213982_Top_table_TMP$ID]
GSE213982_Top_table_TMP$SE = SE
all(names(SE) == GSE213982_Top_table_TMP$ID) # TRUE

# Annotating results
GSE213982_probes_TMP = GSE144136_GSE213982_probes
GSE213982_probes_TMP = GSE213982_probes_TMP[GSE213982_Top_table_TMP$ID, ]
GSE213982_Top_table_TMP$Gene_symbol = GSE213982_probes_TMP$Gene_symbol
GSE213982_Top_table_TMP$Gene_symbol_non_hgnc = NA
GSE213982_Top_table_TMP$Tissue = GSE213982_pheno_curated$TISSUE[1]
GSE213982_Top_table_TMP$Tissue_type = "Brain"
GSE213982_Top_table_TMP$Technology = "scRNA-seq"


# analysis without covariates per cell (female)

GSE213982_Top_table_no_covar_cell_types = list()

for (i in 1:length(cell_types)){
  
  current_cell = cell_types[i]
  GSE213982_count_tmp = GSE144136_GSE213982_counts_cell_type_list[current_cell][[1]]
  GSE213982_count_tmp = GSE213982_count_tmp[,GSE213982_selection]
  
  # Design matrix
  tmp_design = model.matrix(~ GROUP, data = GSE213982_pheno_curated)
  
  # egeR filter
  GSE213982_tmp_dge = DGEList(counts=GSE213982_count_tmp)
  keep = filterByExpr(GSE213982_tmp_dge, tmp_design)
  GSE213982_tmp_dge = GSE213982_tmp_dge[keep,,keep.lib.sizes=FALSE]
  GSE213982_tmp_dge = calcNormFactors(GSE213982_tmp_dge)
  GSE213982_tmp_voom = voom(GSE213982_tmp_dge, tmp_design, plot=TRUE)
  
  # limma 
  fit = lmFit(GSE213982_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  GSE213982_Top_table_no_covar__cell_type_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  GSE213982_Top_table_no_covar__cell_type_TMP$ID = rownames(GSE213982_Top_table_no_covar__cell_type_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[GSE213982_Top_table_no_covar__cell_type_TMP$ID]
  GSE213982_Top_table_no_covar__cell_type_TMP$SE = SE
  all(names(SE) == GSE213982_Top_table_no_covar__cell_type_TMP$ID) # TRUE
  
  # Annotating results
  GSE213982_probes_TMP = GSE144136_GSE213982_probes
  GSE213982_probes_TMP = GSE213982_probes_TMP[GSE213982_Top_table_no_covar__cell_type_TMP$ID, ]
  GSE213982_Top_table_no_covar__cell_type_TMP$Gene_symbol = GSE213982_probes_TMP$Gene_symbol
  GSE213982_Top_table_no_covar__cell_type_TMP$Gene_symbol_non_hgnc = NA
  GSE213982_Top_table_no_covar__cell_type_TMP$Tissue = paste0(GSE213982_pheno_curated$TISSUE[1], ": ", current_cell)
  
  GSE213982_Top_table_no_covar__cell_type_TMP$Tissue_type = "Brain"
  GSE213982_Top_table_no_covar__cell_type_TMP$Technology = "scRNA-seq"
  
  GSE213982_Top_table_no_covar_cell_types[[i]] = GSE213982_Top_table_no_covar__cell_type_TMP
  
}
GSE213982_Top_table_no_covar_cell_types = do.call(rbind, GSE213982_Top_table_no_covar_cell_types)

dir.create("GSE213982_results")
write.csv(GSE144136_GSE213982_counts_mixed, "GSE213982_results/GSE144136_GSE213982_counts_mixed.csv")
write.csv(GSE144136_GSE213982_counts_cell_type, "GSE213982_results/GSE144136_GSE213982_counts_cell_type.csv")
write.csv(GSE213982_pheno_curated, "GSE213982_results/GSE213982_pheno_curated.csv")
write.csv(GSE213982_Top_table_no_covar_TMP, "GSE213982_results/GSE213982_Top_table_no_covar.csv")
write.csv(GSE213982_Top_table_TMP, "GSE213982_results/GSE213982_Top_table.csv")
write.csv(GSE213982_Top_table_no_covar_cell_types, "GSE213982_results/GSE213982_Top_table_no_covar_cell_types.csv")


rm(list = ls(pattern = "GSE144136"))
gc()
rm(list = ls(pattern = "GSE213982"))
gc()
