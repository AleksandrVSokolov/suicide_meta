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
    
    tryCatch(
      {
        Total_pval = wilcox.test(formula = as.formula(formula_test), data = df)
        Total_pval = Total_pval$p.value
      },
      error = function(cond) {
        message(conditionMessage(cond))
        message("Returning NA for p-value")
        Total_pval <<- NA
      }
    )
    
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

describe_vector_numeric = function(x, Mention_NAs, show_as_percent = FALSE, Describe_min_max = TRUE){
  
  if (all(is.na(x))){
    Missing_count = length(x[is.na(x)])
    Missing_val_percent = (Missing_count/length(x) * 100) %>% round(., digits = 2)
    Var_report = paste0("Missing val: ", Missing_count, " (", Missing_val_percent, "%)")
    return(Var_report)
  }
  
  Mean_var = mean(x, na.rm = TRUE) %>% round(., digits = 2)
  SD_var = sd(x, na.rm = TRUE) %>% round(., digits = 2)
  Var_report = paste0(Mean_var, " ± ", SD_var)
  
  if (is.na(SD_var)){
    Var_report = paste0("Single non-NA value")
  }
  
  if (show_as_percent & Var_report != "Single non-NA value"){
    Mean_var = Mean_var * 100/Mean_var
    SD_var = SD_var * 100/Mean_var
    Var_report = paste0(Mean_var, " ± ", SD_var, "%")
  }
  
  if (Mention_NAs){
    Missing_count = length(x[is.na(x)])
    Missing_val_percent = (Missing_count/length(x) * 100) %>% round(., digits = 2)
    if (Missing_count > 0){
      Var_report = paste0(Var_report, "\nMissing val: ", Missing_count, " (", Missing_val_percent, "%)")
    }
  }
  
  
  if (Describe_min_max){
    Min_var = min(x, na.rm = TRUE) %>% round(., digits = 2)
    Max_var = max(x, na.rm = TRUE) %>% round(., digits = 2)
    if (show_as_percent & !is.na(SD_var)){
      Min_var = Min_var * 100/Mean_var
      Max_var = Max_var * 100/Mean_var
      Var_report = paste0(Var_report, "\n", "Min: ", Min_var, "%, Max: ",Max_var, "%")
    } else {
      Var_report = paste0(Var_report, "\n", "Min: ", Min_var, ", Max: ",Max_var)
    }
  }
  return(Var_report)
}

characterize_dataset_generelized_two_subgroups = function(dataset,
                                                          study_char,
                                                          contrast_col_number,
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
      message = paste0("Bad index: ", main_indeces[i])
      print(message)
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


################### Arrays ###################
################### GSE208338 ###################
# paper https://www.sciencedirect.com/science/article/pii/S0924977X22003169?via%3Dihub
# preprint https://www.medrxiv.org/content/10.1101/2022.08.09.22278128v1.full.pdf
GSE208338_pheno_curated = read.csv("GSE208338_results/GSE208338_pheno_curated.csv")
rownames(GSE208338_pheno_curated) = GSE208338_pheno_curated$X
GSE208338_pheno_curated$X = NULL

# SUICIDE + AGE + SEX + PH + PMI + RIN + C1 + C2 + C3 + C4 + DIAGNOSIS + SV1 + TYPE_OF_DEATH
# Real demographics
GSE208338_pheno_curated$AGE = as.numeric(GSE208338_pheno_curated$AGE)
GSE208338_pheno_curated$SEX = factor(GSE208338_pheno_curated$SEX, levels = c("Female", "Male"))
GSE208338_pheno_curated$PH = as.numeric(GSE208338_pheno_curated$PH)
GSE208338_pheno_curated$PMI = as.numeric(GSE208338_pheno_curated$PMI)
GSE208338_pheno_curated$RIN = as.numeric(GSE208338_pheno_curated$RIN)
GSE208338_pheno_curated$DIAGNOSIS = as.factor(GSE208338_pheno_curated$DIAGNOSIS)
GSE208338_pheno_curated$SUICIDE = factor(GSE208338_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
GSE208338_pheno_curated$TYPE_OF_DEATH = factor(GSE208338_pheno_curated$TYPE_OF_DEATH, levels = c("Natural", "Non-violent", "Violent"))
GSE208338_pheno_curated$TISSUE = factor(GSE208338_pheno_curated$TISSUE, levels = GSE208338_pheno_curated$TISSUE[1])
# Technical covariate
GSE208338_pheno_curated$C1 = as.numeric(GSE208338_pheno_curated$C1)
GSE208338_pheno_curated$C2 = as.numeric(GSE208338_pheno_curated$C2)
GSE208338_pheno_curated$C3 = as.numeric(GSE208338_pheno_curated$C3)
GSE208338_pheno_curated$C4 = as.numeric(GSE208338_pheno_curated$C4)
GSE208338_pheno_curated$SV1 = as.numeric(GSE208338_pheno_curated$SV1)


which(colnames(GSE208338_pheno_curated) == "AGE") # 59
which(colnames(GSE208338_pheno_curated) == "SEX") # 77
which(colnames(GSE208338_pheno_curated) == "PH") # 71
which(colnames(GSE208338_pheno_curated) == "PMI") # 72
which(colnames(GSE208338_pheno_curated) == "RIN") # 73
which(colnames(GSE208338_pheno_curated) == "DIAGNOSIS") # 66
which(colnames(GSE208338_pheno_curated) == "TYPE_OF_DEATH") # 82
which(colnames(GSE208338_pheno_curated) == "SUICIDE") # 78
which(colnames(GSE208338_pheno_curated) == "GEO_ACCESSION") # 2
which(colnames(GSE208338_pheno_curated) == "TISSUE") # 2

GSE208338_charact_df = GSE208338_pheno_curated

GSE208338_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE208338_charact_df, 
                                                                   study_char = "GSE208338", 
                                                                   contrast_col_number = 78,
                                                                   participants_col_number = 2, 
                                                                   model_covariates_col_vector = c(77, 59, 66, 82, 71, 72, 73), 
                                                                   columns_to_characterise_vector = c(77, 59, 66, 82, 71, 72, 73, 81), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)
GSE208338_summary = GSE208338_summary$Table
GSE208338_summary = GSE208338_summary[-3,]
GSE208338_summary = apply(X = GSE208338_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "AGE",
                                                "SEX",
                                                "DIAGNOSIS",
                                                "TYPE_OF_DEATH",
                                                "SUICIDE",
                                                "TISSUE",
                                                "postmortem dorsolateral prefrontal cortex (DLPFC)"), 
                             replacement_vector = c("Dataset includes",
                                                    "Age",
                                                    "Sex",
                                                    "Diagnosis",
                                                    "Type of death",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Dorsolateral prefrontal cortex (DLPFC)"))
  
})

GSE208338_summary = as.data.frame(GSE208338_summary)
GSE208338_summary
openxlsx::write.xlsx(GSE208338_summary, file = "GSE208338_results/GSE208338_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE208338"))


################### GSE5388 ###################
GSE5388_pheno_curated = read.csv("GSE5388_results/GSE5388_pheno_curated.csv")
rownames(GSE5388_pheno_curated) = GSE5388_pheno_curated$X
GSE5388_pheno_curated$X = NULL


GSE5388_pheno_curated$AGE__YEARS_ = as.numeric(GSE5388_pheno_curated$AGE__YEARS_)
GSE5388_pheno_curated$GENDER = factor(GSE5388_pheno_curated$GENDER, levels = c("Female", "Male"))
GSE5388_pheno_curated$BRAIN_PH = as.numeric(GSE5388_pheno_curated$BRAIN_PH)
GSE5388_pheno_curated$POST_MORTEM_INTERVAL__HOURS_ = as.numeric(GSE5388_pheno_curated$POST_MORTEM_INTERVAL__HOURS_)
GSE5388_pheno_curated$SUICIDE = factor(GSE5388_pheno_curated$SUICIDE, levels = c("No", "Yes"))
GSE5388_pheno_curated$FLUPHENAZINE_MG__EQUIVALENTS_CURATED = as.numeric(GSE5388_pheno_curated$FLUPHENAZINE_MG__EQUIVALENTS_CURATED) # 1 NA is added for unknown
GSE5388_pheno_curated$VALPROATE_TREATMENT_CURATED = factor(GSE5388_pheno_curated$VALPROATE_TREATMENT_CURATED, levels = c("No", "Yes"))
GSE5388_pheno_curated$LITHIUM_TREATMENT_CURATED = factor(GSE5388_pheno_curated$LITHIUM_TREATMENT_CURATED, levels = c("No", "Yes"))
GSE5388_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED = as.numeric(GSE5388_pheno_curated$ALCOHOL_ABUSE__RATINGS_SCALE_CURATED) # 1 NA is added for unknown

GSE5388_pheno_curated$DRUG_ABUSE__RATINGS_SCALE_CURATED = stri_replace_all_fixed(GSE5388_pheno_curated$DRUG_ABUSE__RATINGS_SCALE, 
                                                                                 pattern = "0 (no use) to 5 (heavy use)): ",
                                                                                 replacement = "")
GSE5388_pheno_curated$DRUG_ABUSE__RATINGS_SCALE_CURATED = as.numeric(GSE5388_pheno_curated$DRUG_ABUSE__RATINGS_SCALE_CURATED)
# Tissue row
table(GSE5388_pheno_curated$TISSUE)
GSE5388_pheno_curated$TISSUE = factor(GSE5388_pheno_curated$TISSUE, levels = GSE5388_pheno_curated$TISSUE[1])

GSE5388_pheno_curated$ELECTROCONVULSIVE_THERAPY_CURATED = factor(GSE5388_pheno_curated$ELECTROCONVULSIVE_THERAPY_CURATED, levels = c("No", "Yes"))
GSE5388_pheno_curated$DISEASE_STATUS = factor(GSE5388_pheno_curated$DISEASE_STATUS, levels = c("Healthy control", "Bipolar disorder"))

# GSM123233 was excluded due to weird RNA degradation
"GSM123233" %in% GSE5388_pheno_curated$GEO_ACCESSION # FALSE
# drug abuse is ommitted due to high number of missing values

"
Design.matrix = model.matrix(~ SUICIDE + DISEASE_STATUS + GENDER + AGE__YEARS_ + BRAIN_PH + POST_MORTEM_INTERVAL__HOURS_ + FLUPHENAZINE_MG__EQUIVALENTS_CURATED +
                               VALPROATE_TREATMENT_CURATED + LITHIUM_TREATMENT_CURATED + ALCOHOL_ABUSE__RATINGS_SCALE_CURATED, data = GSE5388_pheno_curated_2)

"

# GSM123189 is removed during analysis as it has no data on fluphenazine, and ELECTROCONVULSIVE_THERAPY is not included (1 level)
# "GSM123186" "GSM123189" are excluded (phenotypes); GSM123233 is excluded (degradation)
"GSM123233" %in% GSE5388_pheno_curated$GEO_ACCESSION # FALSE -> already removed

which(colnames(GSE5388_pheno_curated) == "SUICIDE") # 59
which(colnames(GSE5388_pheno_curated) == "GENDER") # 55
which(colnames(GSE5388_pheno_curated) == "AGE__YEARS_") # 46
which(colnames(GSE5388_pheno_curated) == "FLUPHENAZINE_MG__EQUIVALENTS_CURATED") # 62
which(colnames(GSE5388_pheno_curated) == "VALPROATE_TREATMENT_CURATED") # 64
which(colnames(GSE5388_pheno_curated) == "LITHIUM_TREATMENT_CURATED") # 65
which(colnames(GSE5388_pheno_curated) == "ELECTROCONVULSIVE_THERAPY_CURATED") # 63 - NOT USED!
which(colnames(GSE5388_pheno_curated) == "ALCOHOL_ABUSE__RATINGS_SCALE_CURATED") # 66
which(colnames(GSE5388_pheno_curated) == "DRUG_ABUSE__RATINGS_SCALE_CURATED") # 67 - NOT USED!
which(colnames(GSE5388_pheno_curated) == "DISEASE_STATUS") # 50

which(colnames(GSE5388_pheno_curated) == "BRAIN_PH") # 49
which(colnames(GSE5388_pheno_curated) == "POST_MORTEM_INTERVAL__HOURS_") # 57
which(colnames(GSE5388_pheno_curated) == "TISSUE") # 61
which(colnames(GSE5388_pheno_curated) == "GEO_ACCESSION") # 2

GSE5388_charact_df = GSE5388_pheno_curated

GSE5388_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE5388_charact_df, 
                                                                 study_char = "GSE5388", 
                                                                 contrast_col_number = 59,
                                                                 participants_col_number = 2, 
                                                                 model_covariates_col_vector = c(55,46,50,62,64,65,66,49,57), 
                                                                 columns_to_characterise_vector = c(55,46,62,50,64,63, 65,66,67,49,57,61), 
                                                                 Remove_NA_predictors = FALSE, 
                                                                 drop_P = TRUE)
GSE5388_summary = GSE5388_summary$Table
GSE5388_summary = GSE5388_summary[-3,]
GSE5388_summary = apply(X = GSE5388_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "AGE__YEARS_",
                                                "FLUPHENAZINE_MG__EQUIVALENTS_CURATED",
                                                "VALPROATE_TREATMENT_CURATED",
                                                "LITHIUM_TREATMENT_CURATED",
                                                "ALCOHOL_ABUSE__RATINGS_SCALE_CURATED",
                                                "DRUG_ABUSE__RATINGS_SCALE_CURATED",
                                                "BRAIN_PH", 
                                                "POST_MORTEM_INTERVAL__HOURS_",
                                                "GENDER",
                                                "TISSUE",
                                                "DLPFC",
                                                "ELECTROCONVULSIVE_THERAPY_CURATED",
                                                "DISEASE_STATUS"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Age",
                                                    "Fluphenazine (mg. equiv.)",
                                                    "Valproate treatment",
                                                    "Lithium treatment",
                                                    "Alcohol abuse\nrating scale",
                                                    "Drug abuse\nrating scale\n(not used)",
                                                    "PH",
                                                    "PMI(h)",
                                                    "Sex",
                                                    "Tissue",
                                                    "Dorsolateral prefrontal cortex (DLPFC)",
                                                    "Electroconvulsive ther.\n(not used)",
                                                    "Disease status")
                             )
  
})
GSE5388_summary = as.data.frame(GSE5388_summary)
GSE5388_summary
openxlsx::write.xlsx(GSE5388_summary, file = "GSE5388_results/GSE5388_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE5388"))

################### GSE5389 ###################
GSE5389_pheno_curated = read.csv("GSE5389_results/GSE5389_pheno_curated.csv")
rownames(GSE5389_pheno_curated) = GSE5389_pheno_curated$X
GSE5389_pheno_curated$X = NULL

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
GSE5389_pheno_curated$ELECTROCONVULSIVE_THERAPY_CURATED = sapply(GSE5389_pheno_curated$ELECTROCONVULSIVE_THERAPY, function(x){
  
  if (x == "Yes (5 treatments)"){
    return("Yes")
  }
  
  if (x == "N/A"){
    return("No")
  }
  
  return(x)
  
})
GSE5389_pheno_curated$ELECTROCONVULSIVE_THERAPY_CURATED = factor(GSE5389_pheno_curated$ELECTROCONVULSIVE_THERAPY_CURATED, levels = c("No", "Yes"))
GSE5389_pheno_curated$DISEASE_STATUS = factor(GSE5389_pheno_curated$DISEASE_STATUS, levels = c("Healthy control", "Bipolar disorder"))

# Tissue row
table(GSE5389_pheno_curated$TISSUE)
GSE5389_pheno_curated$TISSUE = factor(GSE5389_pheno_curated$TISSUE, levels = GSE5389_pheno_curated$TISSUE[1])

"
Design.matrix = model.matrix(~ SUICIDE + DISEASE_STATUS + GENDER + AGE__YEARS_ + BRAIN_PH + POST_MORTEM_INTERVAL__HOURS_ + FLUPHENAZINE_MG__EQUIVALENTS_CURATED +
                               VALPROATE_TREATMENT_CURATED + LITHIUM_TREATMENT_CURATED + ELECTROCONVULSIVE_THERAPY_CURATED +
                               ALCOHOL_ABUSE__RATINGS_SCALE_CURATED + DRUG_ABUSE__RATINGS_SCALE_CURATED, data = GSE5389_pheno_curated)
"

which(colnames(GSE5389_pheno_curated) == "SUICIDE") # 59
which(colnames(GSE5389_pheno_curated) == "GENDER") # 55
which(colnames(GSE5389_pheno_curated) == "AGE__YEARS_") # 46
which(colnames(GSE5389_pheno_curated) == "FLUPHENAZINE_MG__EQUIVALENTS_CURATED") # 62
which(colnames(GSE5389_pheno_curated) == "VALPROATE_TREATMENT_CURATED") # 64
which(colnames(GSE5389_pheno_curated) == "LITHIUM_TREATMENT_CURATED") # 65
which(colnames(GSE5389_pheno_curated) == "ELECTROCONVULSIVE_THERAPY_CURATED") # 63
which(colnames(GSE5389_pheno_curated) == "ALCOHOL_ABUSE__RATINGS_SCALE_CURATED") # 66
which(colnames(GSE5389_pheno_curated) == "DRUG_ABUSE__RATINGS_SCALE_CURATED") # 67
which(colnames(GSE5389_pheno_curated) == "DISEASE_STATUS") # 50

which(colnames(GSE5389_pheno_curated) == "BRAIN_PH") # 49
which(colnames(GSE5389_pheno_curated) == "POST_MORTEM_INTERVAL__HOURS_") # 57
which(colnames(GSE5389_pheno_curated) == "TISSUE") # 61
which(colnames(GSE5389_pheno_curated) == "GEO_ACCESSION") # 2

GSE5389_charact_df = GSE5389_pheno_curated

GSE5389_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE5389_charact_df, 
                                                                 study_char = "GSE5389", 
                                                                 contrast_col_number = 59,
                                                                 participants_col_number = 2, 
                                                                 model_covariates_col_vector = c(55,46,50,62,63,64,65,66,67,49,57), 
                                                                 columns_to_characterise_vector = c(55,46,50,62,63,64,65,66,67,49,57), 
                                                                 Remove_NA_predictors = FALSE, 
                                                                 drop_P = TRUE)
GSE5389_summary = GSE5389_summary$Table
GSE5389_summary = GSE5389_summary[-3,]
GSE5389_summary
GSE5389_summary = apply(X = GSE5389_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "AGE__YEARS_",
                                                "FLUPHENAZINE_MG__EQUIVALENTS_CURATED",
                                                "VALPROATE_TREATMENT_CURATED",
                                                "LITHIUM_TREATMENT_CURATED",
                                                "ALCOHOL_ABUSE__RATINGS_SCALE_CURATED",
                                                "DRUG_ABUSE__RATINGS_SCALE_CURATED",
                                                "BRAIN_PH", 
                                                "POST_MORTEM_INTERVAL__HOURS_",
                                                "GENDER",
                                                "TISSUE",
                                                "DLPFC",
                                                "ELECTROCONVULSIVE_THERAPY_CURATED",
                                                "DISEASE_STATUS"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Age",
                                                    "Fluphenazine (mg. equiv.)",
                                                    "Valproate treatment",
                                                    "Lithium treatment",
                                                    "Alcohol abuse\nrating scale",
                                                    "Drug abuse\nrating scale",
                                                    "PH",
                                                    "PMI(h)",
                                                    "Sex",
                                                    "Tissue",
                                                    "Dorsolateral prefrontal cortex (DLPFC)",
                                                    "Electroconvulsive ther.",
                                                    "Disease status")
  )
  
})
GSE5389_summary = as.data.frame(GSE5389_summary)
GSE5389_summary
openxlsx::write.xlsx(GSE5389_summary, file = "GSE5389_results/GSE5389_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE5389"))

################### GSE66937 ###################
GSE66937_pheno_curated = read.csv("GSE66937_results/GSE66937_pheno_curated.csv")
rownames(GSE66937_pheno_curated) = GSE66937_pheno_curated$X
GSE66937_pheno_curated$X = NULL

GSE66937_pheno_curated$SUICIDE = factor(GSE66937_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
GSE66937_pheno_curated$RIN = as.numeric(GSE66937_pheno_curated$RIN)
GSE66937_pheno_curated$BATCH = factor(GSE66937_pheno_curated$BATCH, levels = c("2013", "2014"))

which(colnames(GSE66937_pheno_curated) == "SUICIDE") # 39
which(colnames(GSE66937_pheno_curated) == "RIN") # 38
which(colnames(GSE66937_pheno_curated) == "BATCH") # 40
which(colnames(GSE66937_pheno_curated) == "BRAIN_REGION") # 34
which(colnames(GSE66937_pheno_curated) == "GEO_ACCESSION") # 2

GSE66937_charact_df = GSE66937_pheno_curated

GSE66937_tissues = unique(GSE66937_charact_df$BRAIN_REGION)
GSE66937_summary = list()

for (i in 1:length(GSE66937_tissues)){
  GSE66937_tmp_df = GSE66937_charact_df[GSE66937_charact_df$BRAIN_REGION == GSE66937_tissues[i],]
  GSE66937_tmp_df$BRAIN_REGION = factor(GSE66937_tmp_df$BRAIN_REGION, levels = GSE66937_tissues[i])
  GSE66937_summary_tmp = characterize_dataset_generelized_two_subgroups(dataset = GSE66937_tmp_df, 
                                                                        study_char = paste0("GSE66937: ", GSE66937_tissues[i]), 
                                                                        contrast_col_number = 39,
                                                                        participants_col_number = 2, 
                                                                        model_covariates_col_vector = c(38,40), 
                                                                        columns_to_characterise_vector = c(38,40,34), 
                                                                        Remove_NA_predictors = FALSE, 
                                                                        drop_P = TRUE)
  GSE66937_summary_tmp = GSE66937_summary_tmp$Table
  GSE66937_summary_tmp = GSE66937_summary_tmp[-3,]
  GSE66937_summary[[i]] = GSE66937_summary_tmp
}

GSE66937_summary = do.call(rbind, GSE66937_summary)
GSE66937_summary = apply(X = GSE66937_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "BATCH",
                                                "2013",
                                                "2014",
                                                "BRAIN_REGION",
                                                "amygdala",
                                                "hippocampus",
                                                "thalamus",
                                                "prefrontal cortex"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Expression batch",
                                                    "Batch2013",
                                                    "Batch2014",
                                                    "Tissue",
                                                    "Amygdala",
                                                    "Hippocampus",
                                                    "Thalamus",
                                                    "Prefrontal cortex")
  )
  
})
GSE66937_summary = as.data.frame(GSE66937_summary)
GSE66937_summary
openxlsx::write.xlsx(GSE66937_summary, file = "GSE66937_results/GSE66937_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE66937"))

################### GSE199536 ###################
GSE199536_pheno_curated = read.csv("GSE199536_results/GSE199536_pheno_curated.csv")
rownames(GSE199536_pheno_curated) = GSE199536_pheno_curated$X
GSE199536_pheno_curated$X = NULL

# Design.matrix = model.matrix(~ SUICIDE + AGE__YR_ + PMI__H_A + PH, data = GSE199536_pheno_curated)
# ALL ARE MALE
GSE199536_pheno_curated$SUICIDE
GSE199536_pheno_curated$AGE__YR_
GSE199536_pheno_curated$PMI__H_A
GSE199536_pheno_curated$PH
GSE199536_pheno_curated$SUICIDE = factor(GSE199536_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
GSE199536_pheno_curated$TISSUE = "Habenula"
GSE199536_pheno_curated$TISSUE = factor(GSE199536_pheno_curated$TISSUE, levels = GSE199536_pheno_curated$TISSUE[1])
GSE199536_pheno_curated$SEX = "Male"
GSE199536_pheno_curated$SEX = factor(GSE199536_pheno_curated$SEX, levels = GSE199536_pheno_curated$SEX[1])

which(colnames(GSE199536_pheno_curated) == "SUICIDE") # 35
which(colnames(GSE199536_pheno_curated) == "SEX") # 41
which(colnames(GSE199536_pheno_curated) == "AGE__YR_") # 38
which(colnames(GSE199536_pheno_curated) == "PMI__H_A") # 39
which(colnames(GSE199536_pheno_curated) == "PH") # 40
which(colnames(GSE199536_pheno_curated) == "TISSUE") # 34
which(colnames(GSE199536_pheno_curated) == "GEO_ACCESSION") # 2

GSE199536_charact_df = GSE199536_pheno_curated

GSE199536_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE199536_charact_df, 
                                                                 study_char = "GSE199536", 
                                                                 contrast_col_number = 35,
                                                                 participants_col_number = 2, 
                                                                 model_covariates_col_vector = c(41,38,39,40), 
                                                                 columns_to_characterise_vector = c(41,38,39,40,34), 
                                                                 Remove_NA_predictors = FALSE, 
                                                                 drop_P = TRUE)
GSE199536_summary = GSE199536_summary$Table
GSE199536_summary = GSE199536_summary[-3,]
GSE199536_summary
GSE199536_summary = apply(X = GSE199536_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "AGE__YR_",
                                                "PMI__H_A",
                                                "TISSUE",
                                                "SEX"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Age",
                                                    "PMI(h)",
                                                    "Tissue",
                                                    "Sex")
  )
  
})
GSE199536_summary = as.data.frame(GSE199536_summary)
GSE199536_summary
openxlsx::write.xlsx(GSE199536_summary, file = "GSE199536_results/GSE199536_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE199536"))

################### GSE92538 U133A ###################
GSE92538_pheno_curated = read.csv("GSE92538_U133A_results/GSE92538_pheno_curated.csv")
rownames(GSE92538_pheno_curated) = GSE92538_pheno_curated$X
GSE92538_pheno_curated$X = NULL


GSE92538_pheno_curated$GENDER = factor(GSE92538_pheno_curated$GENDER, levels = c("Female", "Male"))
GSE92538_pheno_curated$AGE = as.numeric(GSE92538_pheno_curated$AGE)
GSE92538_pheno_curated$POST_MORTEM_INTERVAL = as.numeric(GSE92538_pheno_curated$POST_MORTEM_INTERVAL)
GSE92538_pheno_curated$QC_BATCH = paste0("Batch", GSE92538_pheno_curated$QC_BATCH)
GSE92538_pheno_curated$QC_BATCH = factor(GSE92538_pheno_curated$QC_BATCH)
GSE92538_pheno_curated$RACE = factor(GSE92538_pheno_curated$RACE, levels = c("Caucasian", "Asian", "Other", "Hispanic", "African American"))
GSE92538_pheno_curated$SUICIDE__1_YES_ = factor(GSE92538_pheno_curated$SUICIDE__1_YES_, levels = c("Control", "Suicide"))
GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_ = as.numeric(GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_) # produced 4 NA

GSE92538_pheno_curated$AGONAL_FACTOR = as.character(GSE92538_pheno_curated$AGONAL_FACTOR)
GSE92538_pheno_curated$AGONAL_FACTOR_binary = ifelse(GSE92538_pheno_curated$AGONAL_FACTOR=="0", "AFS=0", "AFS>=1")
GSE92538_pheno_curated$AGONAL_FACTOR_binary = factor(GSE92538_pheno_curated$AGONAL_FACTOR_binary, levels = c("AFS=0", "AFS>=1"))
GSE92538_pheno_curated$DIAGNOSIS = factor(GSE92538_pheno_curated$DIAGNOSIS, levels = c("Control", 
                                                                                           "Major Depressive Disorder",
                                                                                           "Bipolar Disorder", 
                                                                                           "Schizophrenia"))

table(GSE92538_pheno_curated$TISSUE)
"Dorsolateral Prefrontal Cortex 
                           97 "

GSE92538_pheno_curated$TISSUE = "Dorsolateral prefrontal cortex (DLPFC)"
GSE92538_pheno_curated$TISSUE = factor(GSE92538_pheno_curated$TISSUE, levels = GSE92538_pheno_curated$TISSUE[1])

GSE92538_pheno_curated$GENDER
GSE92538_pheno_curated$AGE
GSE92538_pheno_curated$POST_MORTEM_INTERVAL
GSE92538_pheno_curated$QC_BATCH
GSE92538_pheno_curated$RACE
GSE92538_pheno_curated$SUICIDE__1_YES_
GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_


which(colnames(GSE92538_pheno_curated) == "SUICIDE__1_YES_") # 41
which(colnames(GSE92538_pheno_curated) == "GENDER") # 37
which(colnames(GSE92538_pheno_curated) == "AGE") # 39
which(colnames(GSE92538_pheno_curated) == "DIAGNOSIS") # 33
which(colnames(GSE92538_pheno_curated) == "AGONAL_FACTOR_binary") # 44
which(colnames(GSE92538_pheno_curated) == "RACE") # 38
which(colnames(GSE92538_pheno_curated) == "POST_MORTEM_INTERVAL") # 40
which(colnames(GSE92538_pheno_curated) == "QC_BATCH") # 43
which(colnames(GSE92538_pheno_curated) == "TISSUE_PH__CEREBELLUM_") # 36
which(colnames(GSE92538_pheno_curated) == "TISSUE") # 42
which(colnames(GSE92538_pheno_curated) == "GEO_ACCESSION") # 2


GSE92538_charact_df = GSE92538_pheno_curated

GSE92538_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE92538_charact_df, 
                                                                 study_char = "GSE92538 U133A", 
                                                                 contrast_col_number = 41,
                                                                 participants_col_number = 2, 
                                                                 model_covariates_col_vector = c(37,39,33,44,38,40,43,36), 
                                                                 columns_to_characterise_vector = c(37,39,33,44,38,40,43,36,42), 
                                                                 Remove_NA_predictors = FALSE, 
                                                                 drop_P = TRUE)
GSE92538_summary = GSE92538_summary$Table
GSE92538_summary = GSE92538_summary[-3,]
GSE92538_summary
GSE92538_summary = apply(X = GSE92538_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE__1_YES_",
                                                "GENDER",
                                                "AGE",
                                                "RACE",
                                                "POST_MORTEM_INTERVAL",
                                                "QC_BATCH",
                                                "TISSUE_PH__CEREBELLUM_",
                                                "TISSUE",
                                                "DIAGNOSIS",
                                                "AGONAL_FACTOR_binary"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Sex",
                                                    "Age",
                                                    "Race",
                                                    "PMI",
                                                    "Expression batch",
                                                    "PH (cerebellum)",
                                                    "Tissue",
                                                    "Diagnosis",
                                                    "Agonal factor score")
  )
  
})
GSE92538_summary = as.data.frame(GSE92538_summary)
GSE92538_summary
openxlsx::write.xlsx(GSE92538_summary, file = "GSE92538_U133A_results/GSE92538_U133A_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE92538"))

################### GSE92538 PLUS2 ###################
GSE92538_pheno_curated = read.csv("GSE92538_U133_PLUS2_results/GSE92538_pheno_curated.csv")
rownames(GSE92538_pheno_curated) = GSE92538_pheno_curated$X
GSE92538_pheno_curated$X = NULL


GSE92538_pheno_curated$GENDER = factor(GSE92538_pheno_curated$GENDER, levels = c("Female", "Male"))
GSE92538_pheno_curated$AGE = as.numeric(GSE92538_pheno_curated$AGE)
GSE92538_pheno_curated$POST_MORTEM_INTERVAL = as.numeric(GSE92538_pheno_curated$POST_MORTEM_INTERVAL)
GSE92538_pheno_curated$QC_BATCH = paste0("Batch", GSE92538_pheno_curated$QC_BATCH)
GSE92538_pheno_curated$QC_BATCH = factor(GSE92538_pheno_curated$QC_BATCH)
GSE92538_pheno_curated$RACE = factor(GSE92538_pheno_curated$RACE, levels = c("Caucasian", "Asian", "Other", "Hispanic", "African American"))
GSE92538_pheno_curated$SUICIDE__1_YES_ = factor(GSE92538_pheno_curated$SUICIDE__1_YES_, levels = c("Control", "Suicide"))
GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_ = as.numeric(GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_)
GSE92538_pheno_curated$QC_BATCH = drop.levels(GSE92538_pheno_curated$QC_BATCH, reorder = FALSE)
GSE92538_pheno_curated$RACE = drop.levels(GSE92538_pheno_curated$RACE, reorder = FALSE)

GSE92538_pheno_curated$AGONAL_FACTOR = as.character(GSE92538_pheno_curated$AGONAL_FACTOR)
GSE92538_pheno_curated$AGONAL_FACTOR_binary = ifelse(GSE92538_pheno_curated$AGONAL_FACTOR=="0", "AFS=0", "AFS>=1")
GSE92538_pheno_curated$AGONAL_FACTOR_binary = factor(GSE92538_pheno_curated$AGONAL_FACTOR_binary, levels = c("AFS=0", "AFS>=1"))
GSE92538_pheno_curated$DIAGNOSIS = factor(GSE92538_pheno_curated$DIAGNOSIS, levels = c("Control", 
                                                                                       "Major Depressive Disorder",
                                                                                       "Bipolar Disorder", 
                                                                                       "Schizophrenia"))



table(GSE92538_pheno_curated$TISSUE)
"Dorsolateral Prefrontal Cortex 
                           75 "

GSE92538_pheno_curated$TISSUE = "Dorsolateral prefrontal cortex (DLPFC)"
GSE92538_pheno_curated$TISSUE = factor(GSE92538_pheno_curated$TISSUE, levels = GSE92538_pheno_curated$TISSUE[1])

GSE92538_pheno_curated$GENDER
GSE92538_pheno_curated$AGE
GSE92538_pheno_curated$POST_MORTEM_INTERVAL
GSE92538_pheno_curated$QC_BATCH
GSE92538_pheno_curated$RACE
GSE92538_pheno_curated$SUICIDE__1_YES_
GSE92538_pheno_curated$TISSUE_PH__CEREBELLUM_


which(colnames(GSE92538_pheno_curated) == "SUICIDE__1_YES_") # 41
which(colnames(GSE92538_pheno_curated) == "GENDER") # 37
which(colnames(GSE92538_pheno_curated) == "AGE") # 39
which(colnames(GSE92538_pheno_curated) == "DIAGNOSIS") # 33
which(colnames(GSE92538_pheno_curated) == "AGONAL_FACTOR_binary") # 44
which(colnames(GSE92538_pheno_curated) == "RACE") # 38
which(colnames(GSE92538_pheno_curated) == "POST_MORTEM_INTERVAL") # 40
which(colnames(GSE92538_pheno_curated) == "QC_BATCH") # 43
which(colnames(GSE92538_pheno_curated) == "TISSUE_PH__CEREBELLUM_") # 36
which(colnames(GSE92538_pheno_curated) == "TISSUE") # 42
which(colnames(GSE92538_pheno_curated) == "GEO_ACCESSION") # 14


GSE92538_charact_df = GSE92538_pheno_curated

GSE92538_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE92538_charact_df, 
                                                                  study_char = "GSE92538 PLUS2", 
                                                                  contrast_col_number = 41,
                                                                  participants_col_number = 14, 
                                                                  model_covariates_col_vector = c(37,39,33,44,38,40,43,36), 
                                                                  columns_to_characterise_vector = c(37,39,33,44,38,40,43,36,42), 
                                                                  Remove_NA_predictors = FALSE, 
                                                                  drop_P = TRUE)
GSE92538_summary = GSE92538_summary$Table
GSE92538_summary = GSE92538_summary[-3,]
GSE92538_summary
GSE92538_summary = apply(X = GSE92538_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE__1_YES_",
                                                "GENDER",
                                                "AGE",
                                                "RACE",
                                                "POST_MORTEM_INTERVAL",
                                                "QC_BATCH",
                                                "TISSUE_PH__CEREBELLUM_",
                                                "TISSUE",
                                                "DIAGNOSIS",
                                                "AGONAL_FACTOR_binary"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Sex",
                                                    "Age",
                                                    "Race",
                                                    "PMI",
                                                    "Expression batch",
                                                    "PH (cerebellum)",
                                                    "Tissue",
                                                    "Diagnosis",
                                                    "Agonal factor score")
  )
  
})
GSE92538_summary = as.data.frame(GSE92538_summary)
GSE92538_summary
openxlsx::write.xlsx(GSE92538_summary, file = "GSE92538_U133_PLUS2_results/GSE92538_PLUS2_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE92538"))

################### RNA Seq cohorts ###################

################### GSE102556 ###################
GSE102556_pheno_curated_2 = read.csv("GSE102556_results/GSE102556_pheno_curated.csv")
rownames(GSE102556_pheno_curated_2) = GSE102556_pheno_curated_2$X
GSE102556_pheno_curated_2$X = NULL


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
GSE102556_pheno_curated_2$ALCOOL = factor(GSE102556_pheno_curated_2$ALCOOL, levels = c("NO", "YES"))
table(GSE102556_pheno_curated_2$SUICIDE, GSE102556_pheno_curated_2$TISSUE)
"         Anterior Insula (aINS) Cingulate gyrus 25 (Cg25) Dorsolateral prefrontal cortex (dlPFC; BA8/9) Nucleus Accumbens (Nac) Orbitofrontal (OFC; BA11) Subiculum (Sub)
  CONTROL                     11                         8                                            11                      11                        11               8
  SUICIDE                     37                        20                                            37                      37                        36              35"
GSE102556_pheno_curated_2$SUICIDE = factor(GSE102556_pheno_curated_2$SUICIDE, levels = c("CONTROL", "SUICIDE"))
GSE102556_pheno_curated_2$GENDER = factor(GSE102556_pheno_curated_2$GENDER, levels = c("FEMALE", "MALE"))
GSE102556_pheno_curated_2$MEDICATION = factor(GSE102556_pheno_curated_2$MEDICATION, levels = c("NO", "YES"))
GSE102556_pheno_curated_2$PH = as.numeric(GSE102556_pheno_curated_2$PH)
GSE102556_pheno_curated_2$PMI = as.numeric(GSE102556_pheno_curated_2$PMI)
GSE102556_pheno_curated_2$RIN = as.numeric(GSE102556_pheno_curated_2$RIN)
GSE102556_pheno_curated_2$PHENOTYPE = factor(GSE102556_pheno_curated_2$PHENOTYPE, levels = c("CTRL", "MDD"))

# tmp_design = model.matrix(~ SUICIDE + PHENOTYPE + GENDER + AGE + ALCOOL + MEDICATION + RIN + PMI + PH, data = GSE102556_pheno_curated_2_tmp)

which(colnames(GSE102556_pheno_curated_2) == "SUICIDE") # 65
which(colnames(GSE102556_pheno_curated_2) == "GENDER") # 54
which(colnames(GSE102556_pheno_curated_2) == "AGE") # 49
which(colnames(GSE102556_pheno_curated_2) == "PHENOTYPE") # 58
which(colnames(GSE102556_pheno_curated_2) == "ALCOOL") # 50
which(colnames(GSE102556_pheno_curated_2) == "MEDICATION") # 56
which(colnames(GSE102556_pheno_curated_2) == "RIN") # 60
which(colnames(GSE102556_pheno_curated_2) == "PMI") # 59
which(colnames(GSE102556_pheno_curated_2) == "PH") # 57
which(colnames(GSE102556_pheno_curated_2) == "TISSUE") # 62
which(colnames(GSE102556_pheno_curated_2) == "GEO_ACCESSION") # 2


GSE102556_charact_df = GSE102556_pheno_curated_2

GSE102556_tissues = unique(GSE102556_charact_df$TISSUE)
GSE102556_summary = list()

for (i in 1:length(GSE102556_tissues)){
  GSE102556_tmp_df = GSE102556_charact_df[GSE102556_charact_df$TISSUE == GSE102556_tissues[i],]
  GSE102556_tmp_df$TISSUE = factor(GSE102556_tmp_df$TISSUE, levels = GSE102556_tissues[i])
  GSE102556_summary_tmp = characterize_dataset_generelized_two_subgroups(dataset = GSE102556_tmp_df, 
                                                                        study_char = paste0("GSE102556: ", GSE102556_tissues[i]), 
                                                                        contrast_col_number = 65,
                                                                        participants_col_number = 2, 
                                                                        model_covariates_col_vector = c(54,49,58,50,56,60,59,57), 
                                                                        columns_to_characterise_vector = c(54,49,58,50,56,60,59,57,62), 
                                                                        Remove_NA_predictors = FALSE, 
                                                                        drop_P = TRUE)
  GSE102556_summary_tmp = GSE102556_summary_tmp$Table
  GSE102556_summary_tmp = GSE102556_summary_tmp[-3,]
  GSE102556_summary[[i]] = GSE102556_summary_tmp
}

GSE102556_summary = do.call(rbind, GSE102556_summary)
GSE102556_summary

GSE102556_summary = apply(X = GSE102556_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "\nSUICIDE",
                                                "CONTROL",
                                                "SUICIDE",
                                                "GENDER",
                                                "AGE",
                                                "PHENOTYPE",
                                                "ALCOOL",
                                                "MEDICATION",
                                                "YES",
                                                "NO",
                                                "CTRL",
                                                "dlPFC",
                                                "FEMALE",
                                                "MALE",
                                                "TISSUE"), 
                             replacement_vector = c("Dataset includes",
                                                    "\nSuicide",
                                                    "Control",
                                                    "Suicide status",
                                                    "Sex",
                                                    "Age",
                                                    "Diagnosis",
                                                    "Alcohol\nintake",
                                                    "Medication\nintake",
                                                    "Yes",
                                                    "No",
                                                    "Control",
                                                    "DLPFC",
                                                    "Female",
                                                    "Male",
                                                    "Tissue")
  )
  
})
GSE102556_summary = as.data.frame(GSE102556_summary)
GSE102556_summary
openxlsx::write.xlsx(GSE102556_summary, file = "GSE102556_results/GSE102556_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE102556"))

################### GSE243356 ###################
GSE243356_pheno_curated = read.csv("GSE243356_results/GSE243356_pheno_curated.csv")
rownames(GSE243356_pheno_curated) = GSE243356_pheno_curated$X
GSE243356_pheno_curated$X = NULL

GSE243356_pheno_curated$GROUP = factor(GSE243356_pheno_curated$GROUP, levels = c("healthy", "suicide"))
GSE243356_pheno_curated$TISSUE = "Temporal cortex (BA20 and BA36)"
GSE243356_pheno_curated$TISSUE  = factor(GSE243356_pheno_curated$TISSUE , levels = GSE243356_pheno_curated$TISSUE[1])

which(colnames(GSE243356_pheno_curated) == "GROUP") # 35
which(colnames(GSE243356_pheno_curated) == "TISSUE") # 36
which(colnames(GSE243356_pheno_curated) == "GEO_ACCESSION") # 2


GSE243356_charact_df = GSE243356_pheno_curated

GSE243356_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE243356_charact_df, 
                                                                  study_char = "GSE243356", 
                                                                  contrast_col_number = 35,
                                                                  participants_col_number = 2, 
                                                                  model_covariates_col_vector = c(35), 
                                                                  columns_to_characterise_vector = c(35,36), 
                                                                  Remove_NA_predictors = FALSE, 
                                                                  drop_P = TRUE)
GSE243356_summary = GSE243356_summary$Table
GSE243356_summary = GSE243356_summary[-3,]
GSE243356_summary
GSE243356_summary = apply(X = GSE243356_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "GROUP",
                                                "healthy",
                                                "suicide",
                                                "TISSUE"), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Control",
                                                    "Suicide",
                                                    "Tissue")
  )
  
})
GSE243356_summary = as.data.frame(GSE243356_summary)
GSE243356_summary
openxlsx::write.xlsx(GSE243356_summary, file = "GSE243356_results/GSE243356_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE243356"))

################### GSE248260 ###################
GSE248260_pheno_curated = read.csv("GSE248260_results/GSE248260_pheno_curated.csv")
rownames(GSE248260_pheno_curated) = GSE248260_pheno_curated$X
GSE248260_pheno_curated$X = NULL

GSE248260_pheno_curated$AGE = as.numeric(GSE248260_pheno_curated$AGE)
GSE248260_pheno_curated$DIAGNOSIS = factor(GSE248260_pheno_curated$DIAGNOSIS, levels = c("Normal", "MDD"))
GSE248260_pheno_curated$SEX = factor(GSE248260_pheno_curated$SEX, levels = c("Female", "Male"))
GSE248260_pheno_curated$SUICIDE = factor(GSE248260_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
GSE248260_pheno_curated$Drug_intake = factor(GSE248260_pheno_curated$Drug_intake, levels = c("Not reported", "Drugs"))

table(GSE248260_pheno_curated$TISSUE)
"ventral white matter (BA 47) 
                          24 "
GSE248260_pheno_curated$TISSUE = "Ventral white matter (BA 47)"
GSE248260_pheno_curated$TISSUE = factor(GSE248260_pheno_curated$TISSUE, GSE248260_pheno_curated$TISSUE[1])
# tmp_design = model.matrix(~ SUICIDE + AGE + SEX + PMI + pH + RIN + Drug_intake, data = GSE248260_pheno_curated_2)

which(colnames(GSE248260_pheno_curated) == "SUICIDE") # 45
which(colnames(GSE248260_pheno_curated) == "SEX") # 44
which(colnames(GSE248260_pheno_curated) == "AGE") # 42
which(colnames(GSE248260_pheno_curated) == "Drug_intake") # 54
which(colnames(GSE248260_pheno_curated) == "PMI") # 48
which(colnames(GSE248260_pheno_curated) == "pH") # 49
which(colnames(GSE248260_pheno_curated) == "RIN") # 47
which(colnames(GSE248260_pheno_curated) == "TISSUE") # 46
which(colnames(GSE248260_pheno_curated) == "GEO_ACCESSION") # 2

GSE248260_charact_df = GSE248260_pheno_curated

GSE248260_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE248260_charact_df, 
                                                                   study_char = "GSE248260", 
                                                                   contrast_col_number = 45,
                                                                   participants_col_number = 2, 
                                                                   model_covariates_col_vector = c(44,42,54,48,49,47), 
                                                                   columns_to_characterise_vector = c(44,42,54,48,49,47,46), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)
GSE248260_summary = GSE248260_summary$Table
GSE248260_summary = GSE248260_summary[-3,]
GSE248260_summary
GSE248260_summary = apply(X = GSE248260_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "TISSUE",
                                                "SEX",
                                                "AGE",
                                                "Drug_intake",
                                                "pH"),
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Sex",
                                                    "Age",
                                                    "Drug\nintake",
                                                    "PH")
  )
  
})
GSE248260_summary = as.data.frame(GSE248260_summary)
GSE248260_summary
openxlsx::write.xlsx(GSE248260_summary, file = "GSE248260_results/GSE248260_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE248260"))

################### GSE247998 ###################

GSE247998_pheno_curated = read.csv("GSE247998_results/GSE247998_pheno_curated.csv")
rownames(GSE247998_pheno_curated) = GSE247998_pheno_curated$X
GSE247998_pheno_curated$X = NULL


GSE247998_pheno_curated$AGE = as.numeric(GSE247998_pheno_curated$AGE)
GSE247998_pheno_curated$DIAGNOSIS = factor(GSE247998_pheno_curated$DIAGNOSIS, levels = c("Healthy control", "MDD"))
GSE247998_pheno_curated$SEX = factor(GSE247998_pheno_curated$SEX, levels = c("Female", "Male"))
GSE247998_pheno_curated$SUICIDE = factor(GSE247998_pheno_curated$SUICIDE, levels = c("Control", "Suicide"))
table(GSE247998_pheno_curated$SUICIDE, GSE247998_pheno_curated$DIAGNOSIS)
table(GSE247998_pheno_curated$SUICIDE, GSE247998_pheno_curated$SI_GROUP)
GSE247998_pheno_curated$TISSUE = factor(GSE247998_pheno_curated$TISSUE, levels = GSE247998_pheno_curated$TISSUE[1])

which(colnames(GSE247998_pheno_curated) == "SUICIDE") # 49
which(colnames(GSE247998_pheno_curated) == "SEX") # 45
which(colnames(GSE247998_pheno_curated) == "AGE") # 43
which(colnames(GSE247998_pheno_curated) == "DIAGNOSIS") # 44
which(colnames(GSE247998_pheno_curated) == "TISSUE") # 48
which(colnames(GSE247998_pheno_curated) == "GEO_ACCESSION") # 2

# tmp_design = model.matrix(~ SUICIDE + AGE + SEX + DIAGNOSIS, data = GSE247998_pheno_curated)


GSE247998_charact_df = GSE247998_pheno_curated

GSE247998_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE247998_charact_df, 
                                                                   study_char = "GSE247998", 
                                                                   contrast_col_number = 49,
                                                                   participants_col_number = 2, 
                                                                   model_covariates_col_vector = c(45,43,44), 
                                                                   columns_to_characterise_vector = c(45,43,44,48), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)
GSE247998_summary = GSE247998_summary$Table
GSE247998_summary = GSE247998_summary[-3,]
GSE247998_summary
GSE247998_summary = apply(X = GSE247998_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "TISSUE",
                                                "SEX",
                                                "AGE",
                                                "DIAGNOSIS"),
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Sex",
                                                    "Age",
                                                    "Diagnosis")
  )
  
})
GSE247998_summary = as.data.frame(GSE247998_summary)
GSE247998_summary
openxlsx::write.xlsx(GSE247998_summary, file = "GSE247998_results/GSE247998_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE247998"))


################### GSE202537 ###################

GSE202537_pheno_curated_2 = read.csv("GSE202537_results/GSE202537_pheno_curated.csv")
rownames(GSE202537_pheno_curated_2) = GSE202537_pheno_curated_2$X
GSE202537_pheno_curated_2$X = NULL

GSE202537_pheno_curated_2$PH
GSE202537_pheno_curated_2$PMI
GSE202537_pheno_curated_2$RACE
GSE202537_pheno_curated_2$RIN
GSE202537_pheno_curated_2$TISSUESTORAGETIME

# curation
GSE202537_pheno_curated_2$AGE = as.numeric(GSE202537_pheno_curated_2$AGE)
GSE202537_pheno_curated_2$BMI = as.numeric(GSE202537_pheno_curated_2$BMI)
GSE202537_pheno_curated_2$DISEASE_STATE = factor(GSE202537_pheno_curated_2$DISEASE_STATE, levels = c("match control",
                                                                                                     "psychosis_schizophrenia",
                                                                                                     "psychosis_bipolar"))
GSE202537_pheno_curated_2$GENDER = factor(GSE202537_pheno_curated_2$GENDER, levels = c("Female", "Male"))
GSE202537_pheno_curated_2$SUICIDE = factor(GSE202537_pheno_curated_2$SUICIDE, levels = c("Control", "Suicide"))
GSE202537_pheno_curated_2$PH = as.numeric(GSE202537_pheno_curated_2$PH)
GSE202537_pheno_curated_2$PMI = as.numeric(GSE202537_pheno_curated_2$PMI)
GSE202537_pheno_curated_2$RACE = factor(GSE202537_pheno_curated_2$RACE, levels = c("White", "Black"))
GSE202537_pheno_curated_2$RIN = as.numeric(GSE202537_pheno_curated_2$RIN)

#  tmp_design = model.matrix(~ SUICIDE + AGE + BMI + DISEASE_STATE + GENDER + PH + PMI + RACE + RIN, data = GSE202537_pheno_curated_2_tmp)

which(colnames(GSE202537_pheno_curated_2) == "SUICIDE") # 69
which(colnames(GSE202537_pheno_curated_2) == "GENDER") # 57
which(colnames(GSE202537_pheno_curated_2) == "AGE") # 53
which(colnames(GSE202537_pheno_curated_2) == "BMI") # 54
which(colnames(GSE202537_pheno_curated_2) == "DISEASE_STATE") # 56
which(colnames(GSE202537_pheno_curated_2) == "RACE") # 63
which(colnames(GSE202537_pheno_curated_2) == "RIN") # 64
which(colnames(GSE202537_pheno_curated_2) == "PMI") # 62
which(colnames(GSE202537_pheno_curated_2) == "PH") # 61
which(colnames(GSE202537_pheno_curated_2) == "TISSUE") # 66
which(colnames(GSE202537_pheno_curated_2) == "GEO_ACCESSION") # 2

GSE202537_charact_df = GSE202537_pheno_curated_2

GSE202537_tissues = unique(GSE202537_charact_df$TISSUE)
GSE202537_summary = list()

for (i in 1:length(GSE202537_tissues)){
  GSE202537_tmp_df = GSE202537_charact_df[GSE202537_charact_df$TISSUE == GSE202537_tissues[i],]
  GSE202537_tmp_df$TISSUE = factor(GSE202537_tmp_df$TISSUE, levels = GSE202537_tissues[i])
  GSE202537_summary_tmp = characterize_dataset_generelized_two_subgroups(dataset = GSE202537_tmp_df, 
                                                                         study_char = paste0("GSE202537: ", GSE202537_tissues[i]), 
                                                                         contrast_col_number = 69,
                                                                         participants_col_number = 2, 
                                                                         model_covariates_col_vector = c(57,53,54,56,63,64,62,61), 
                                                                         columns_to_characterise_vector = c(57,53,54,56,63,64,62,61,66), 
                                                                         Remove_NA_predictors = FALSE, 
                                                                         drop_P = TRUE)
  GSE202537_summary_tmp = GSE202537_summary_tmp$Table
  GSE202537_summary_tmp = GSE202537_summary_tmp[-3,]
  GSE202537_summary[[i]] = GSE202537_summary_tmp
}

GSE202537_summary = do.call(rbind, GSE202537_summary)
GSE202537_summary
GSE202537_summary = apply(X = GSE202537_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "TISSUE",
                                                "GENDER",
                                                "AGE",
                                                "DISEASE_STATE",
                                                "match control",
                                                "psychosis_schizophrenia",
                                                "psychosis_bipolar",
                                                "RACE"),
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Sex",
                                                    "Age",
                                                    "Diagnosis",
                                                    "Control",
                                                    "Schizophrenia",
                                                    "Bipolar disorder",
                                                    "Race")
  )
  
})
GSE202537_summary = as.data.frame(GSE202537_summary)
GSE202537_summary
openxlsx::write.xlsx(GSE202537_summary, file = "GSE202537_results/GSE202537_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE202537"))

################### GSE101521 ###################

GSE101521_pheno_curated = read.csv("GSE101521_results/GSE101521_pheno_curated.csv")
rownames(GSE101521_pheno_curated) = GSE101521_pheno_curated$X
GSE101521_pheno_curated$X = NULL


GSE101521_pheno_curated$AGE__YRS_ = as.numeric(GSE101521_pheno_curated$AGE__YRS_)
GSE101521_pheno_curated$BRAIN_PH = as.numeric(GSE101521_pheno_curated$BRAIN_PH)
GSE101521_pheno_curated$PMI = as.numeric(GSE101521_pheno_curated$PMI) 
GSE101521_pheno_curated$RIN = as.numeric(GSE101521_pheno_curated$RIN)
GSE101521_pheno_curated$SEX = factor(GSE101521_pheno_curated$SEX, levels = c("Female", "Male"))
GSE101521_pheno_curated$SUICIDE = factor(GSE101521_pheno_curated$SUICIDE, levels = c("Non-suicide", "Suicide"))
GSE101521_pheno_curated$Depression = factor(GSE101521_pheno_curated$Depression, levels = c("Control", "MDD"))
GSE101521_pheno_curated$TISSUE = "Dorsolateral prefrontal cortex (DLPFC, BA9)"
GSE101521_pheno_curated$TISSUE = factor(GSE101521_pheno_curated$TISSUE, levels = GSE101521_pheno_curated$TISSUE[1])

# tmp_design = model.matrix(~ SUICIDE + Depression + AGE__YRS_ + SEX + BRAIN_PH + PMI + RIN, data = GSE101521_pheno_curated_tmp)


which(colnames(GSE101521_pheno_curated) == "SUICIDE") # 53
which(colnames(GSE101521_pheno_curated) == "SEX") # 51
which(colnames(GSE101521_pheno_curated) == "AGE__YRS_") # 46
which(colnames(GSE101521_pheno_curated) == "Depression") # 55
which(colnames(GSE101521_pheno_curated) == "BRAIN_PH") # 47
which(colnames(GSE101521_pheno_curated) == "PMI") # 49
which(colnames(GSE101521_pheno_curated) == "RIN") # 50
which(colnames(GSE101521_pheno_curated) == "TISSUE") # 52
which(colnames(GSE101521_pheno_curated) == "GEO_ACCESSION") # 2

GSE101521_charact_df = GSE101521_pheno_curated

GSE101521_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE101521_charact_df, 
                                                                   study_char = "GSE101521", 
                                                                   contrast_col_number = 53,
                                                                   participants_col_number = 2, 
                                                                   model_covariates_col_vector = c(51,46,55,47,49,50), 
                                                                   columns_to_characterise_vector = c(51,46,55,47,49,50,52), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)
GSE101521_summary = GSE101521_summary$Table
GSE101521_summary = GSE101521_summary[-3,]
GSE101521_summary
GSE101521_summary = apply(X = GSE101521_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUICIDE",
                                                "TISSUE",
                                                "SEX",
                                                "AGE__YRS_",
                                                "Depression",
                                                "BRAIN_PH",
                                                "Non-suicide"),
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Sex",
                                                    "Age",
                                                    "Diagnosis",
                                                    "PH",
                                                    "Control")
  )
  
})
GSE101521_summary = as.data.frame(GSE101521_summary)
GSE101521_summary
openxlsx::write.xlsx(GSE101521_summary, file = "GSE101521_results/GSE101521_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE101521"))


################### Single-cell RNA Seq cohorts ###################


################### GSE144136 ###################

GSE144136_pheno_curated = read.csv("GSE144136_results/GSE144136_pheno_curated.csv")
rownames(GSE144136_pheno_curated) = GSE144136_pheno_curated$X
GSE144136_pheno_curated$X = NULL

GSE144136_pheno_curated$GROUP = factor(GSE144136_pheno_curated$GROUP, levels = c("Control", "Major Depressive Disorder (MDD)"))
GSE144136_pheno_curated$SEX = factor(GSE144136_pheno_curated$SEX, levels = GSE144136_pheno_curated$SEX[1])
GSE144136_pheno_curated$TISSUE = "Dorsolateral prefrontal cortex (DLPFC, BA9)"
GSE144136_pheno_curated$TISSUE = factor(GSE144136_pheno_curated$TISSUE, levels = GSE144136_pheno_curated$TISSUE[1])


# tmp_design = model.matrix(~ GROUP + Ast + End + ExN + InN + Mix + Oli + OPC, data = GSE144136_pheno_curated)

which(colnames(GSE144136_pheno_curated) == "GROUP") # 42
which(colnames(GSE144136_pheno_curated) == "SEX") # 43
which(colnames(GSE144136_pheno_curated) == "Ast") # 46
which(colnames(GSE144136_pheno_curated) == "End") # 47
which(colnames(GSE144136_pheno_curated) == "ExN") # 48
which(colnames(GSE144136_pheno_curated) == "InN") # 49
which(colnames(GSE144136_pheno_curated) == "Mix") # 51
which(colnames(GSE144136_pheno_curated) == "Oli") # 52
which(colnames(GSE144136_pheno_curated) == "OPC") # 53
which(colnames(GSE144136_pheno_curated) == "TISSUE") # 44
which(colnames(GSE144136_pheno_curated) == "GEO_ACCESSION") # 2

GSE144136_charact_df = GSE144136_pheno_curated

GSE144136_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE144136_charact_df, 
                                                                   study_char = "GSE144136", 
                                                                   contrast_col_number = 42,
                                                                   participants_col_number = 2, 
                                                                   model_covariates_col_vector = c(46,47,48,49,51,52,53), 
                                                                   columns_to_characterise_vector = c(43,46,47,48,49,51,52,53), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)
GSE144136_summary = GSE144136_summary$Table
GSE144136_summary = GSE144136_summary[-3,]
GSE144136_summary
GSE144136_summary = apply(X = GSE144136_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "GROUP",
                                                "TISSUE",
                                                "SEX",
                                                "Major Depressive Disorder (MDD)",
                                                "Ast",
                                                "End",
                                                "ExN",
                                                "InN",
                                                "Mix",
                                                "Oli",
                                                "OPC"),
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Sex",
                                                    "MDD & Suicide",
                                                    "Proportion: Ast",
                                                    "Proportion: End",
                                                    "Proportion: ExN",
                                                    "Proportion: InN",
                                                    "Proportion: Mix",
                                                    "Proportion: Oli",
                                                    "Proportion: OPC")
  )
  
})
GSE144136_summary = as.data.frame(GSE144136_summary)
GSE144136_summary
openxlsx::write.xlsx(GSE144136_summary, file = "GSE144136_results/GSE144136_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE144136"))


################### GSE213982 ###################

GSE213982_pheno_curated = read.csv("GSE213982_results/GSE213982_pheno_curated.csv")
rownames(GSE213982_pheno_curated) = GSE213982_pheno_curated$X
GSE213982_pheno_curated$X = NULL

GSE213982_pheno_curated$GROUP = factor(GSE213982_pheno_curated$GROUP, levels = c("Control", "Case"))
GSE213982_pheno_curated$SEX = factor(GSE213982_pheno_curated$SEX, levels = GSE213982_pheno_curated$SEX[1])
GSE213982_pheno_curated$TISSUE = "Dorsolateral prefrontal cortex (DLPFC, BA9)"
GSE213982_pheno_curated$TISSUE = factor(GSE213982_pheno_curated$TISSUE, levels = GSE213982_pheno_curated$TISSUE[1])


# tmp_design = model.matrix(~ GROUP + Ast + End + ExN + InN + Mix + Oli + OPC, data = GSE213982_pheno_curated)

which(colnames(GSE213982_pheno_curated) == "GROUP") # 43
which(colnames(GSE213982_pheno_curated) == "SEX") # 44
which(colnames(GSE213982_pheno_curated) == "Ast") # 47
which(colnames(GSE213982_pheno_curated) == "End") # 48
which(colnames(GSE213982_pheno_curated) == "ExN") # 49
which(colnames(GSE213982_pheno_curated) == "InN") # 50
which(colnames(GSE213982_pheno_curated) == "Mix") # 52
which(colnames(GSE213982_pheno_curated) == "Oli") # 53
which(colnames(GSE213982_pheno_curated) == "OPC") # 54
which(colnames(GSE213982_pheno_curated) == "TISSUE") # 45
which(colnames(GSE213982_pheno_curated) == "GEO_ACCESSION") # 2

GSE213982_charact_df = GSE213982_pheno_curated

GSE213982_summary = characterize_dataset_generelized_two_subgroups(dataset = GSE213982_charact_df, 
                                                                   study_char = "GSE213982", 
                                                                   contrast_col_number = 43,
                                                                   participants_col_number = 2, 
                                                                   model_covariates_col_vector = c(47,48,49,50,52,53,54), 
                                                                   columns_to_characterise_vector = c(44,47,48,49,50,52,53,54), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)
GSE213982_summary = GSE213982_summary$Table
GSE213982_summary = GSE213982_summary[-3,]
GSE213982_summary
GSE213982_summary = apply(X = GSE213982_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "GROUP",
                                                "TISSUE",
                                                "SEX",
                                                "Case",
                                                "Ast",
                                                "End",
                                                "ExN",
                                                "InN",
                                                "Mix",
                                                "Oli",
                                                "OPC"),
                             replacement_vector = c("Dataset includes",
                                                    "Suicide status",
                                                    "Tissue",
                                                    "Sex",
                                                    "MDD & Suicide",
                                                    "Proportion: Ast",
                                                    "Proportion: End",
                                                    "Proportion: ExN",
                                                    "Proportion: InN",
                                                    "Proportion: Mix",
                                                    "Proportion: Oli",
                                                    "Proportion: OPC")
  )
  
})
GSE213982_summary = as.data.frame(GSE213982_summary)
GSE213982_summary
openxlsx::write.xlsx(GSE213982_summary, file = "GSE213982_results/GSE213982_demographics.xlsx", overwrite = TRUE)
rm(list = ls(pattern = "GSE213982"))