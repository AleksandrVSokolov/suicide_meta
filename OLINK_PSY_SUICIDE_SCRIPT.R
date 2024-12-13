setwd("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project")

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
library(metafor)
library(AnnotationDbi)
library(org.Hs.eg.db)

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

test_genes_OLINK_linear_models_generalized = function(PREFIX_Pheno_df,
                                                      PREFIX_Expression_matrix,
                                                      PREFIX_contrast_col_number, 
                                                      PREFIX_participants_col_number, 
                                                      PREFIX_contrast_vector, 
                                                      PREFIX_model.covariates,
                                                      PREFIX_Remove_NA_predictors = TRUE,
                                                      PREFIX_show_all_coef = FALSE,
                                                      PREFIX_round_output = TRUE){
  
  writeLines("Testing all genes")
  
  #Preparing data
  colnames(PREFIX_Pheno_df)[PREFIX_contrast_col_number] = "Case_Control"
  PREFIX_Pheno_df = PREFIX_Pheno_df[PREFIX_Pheno_df$Case_Control %in% PREFIX_contrast_vector, ]
  PREFIX_Pheno_df$Case_Control = factor(PREFIX_Pheno_df$Case_Control, levels = c(PREFIX_contrast_vector[1], PREFIX_contrast_vector[2]))
  
  colnames(PREFIX_Pheno_df)[PREFIX_participants_col_number] = "Participant_ID_FUN"
  
  # getting expression values
  PREFIX_Tested_expression = PREFIX_Expression_matrix
  PREFIX_Tested_expression = PREFIX_Tested_expression[,colnames(PREFIX_Tested_expression) %in% PREFIX_Pheno_df$Participant_ID_FUN]
  
  # preparing pheno df before model
  if (!is.null(PREFIX_model.covariates)){
    PREFIX_Pheno_df = PREFIX_Pheno_df[,c("Participant_ID_FUN", "Case_Control", PREFIX_model.covariates)]
  } else {
    PREFIX_Pheno_df = PREFIX_Pheno_df[,c("Participant_ID_FUN", "Case_Control")]
  }
  PREFIX_MODEL_LIST = list()
  
  for (PREFIX_i in 1:nrow(PREFIX_Tested_expression)){
    PREFIX_Curr_Test = rownames(PREFIX_Tested_expression)[PREFIX_i]
    PREFIX_Curr_values = PREFIX_Tested_expression[PREFIX_Curr_Test,]
    PREFIX_tmp_df = cbind(PREFIX_Pheno_df, t(PREFIX_Curr_values))
    
    if (PREFIX_Remove_NA_predictors){
      writeLines("Removing NA predictors")
      PREFIX_Participants_missing_data = apply(PREFIX_tmp_df, 1, function(x){
        x = as.character(x)
        if (any(is.na(x))){
          return(TRUE)
        } else{
          return(FALSE)
        }
      })
      PREFIX_Participants_missing_data = PREFIX_tmp_df$Participant_ID_FUN[PREFIX_Participants_missing_data]
    }
    
    PREFIX_tmp_df = PREFIX_tmp_df[PREFIX_tmp_df$Participant_ID_FUN %!in% PREFIX_Participants_missing_data,]
    
    if (!is.null(PREFIX_model.covariates)){
      PREFIX_Model.formula.string = paste0(PREFIX_model.covariates, collapse = " + ")
      PREFIX_Model.formula.string = paste0(PREFIX_Curr_Test, " ~ ", "Case_Control", " + ", PREFIX_Model.formula.string)
    } else {
      PREFIX_Model.formula.string = paste0(PREFIX_Curr_Test, " ~ ", "Case_Control")
    }
    
    PREFIX_model_fit = lm(formula = as.formula(PREFIX_Model.formula.string), data = PREFIX_tmp_df)
    if (PREFIX_show_all_coef){
      PREFIX_model_output = modify_coef_log_model(PREFIX_model_fit, model = "lin.model", rounding = FALSE)
    } else {
      PREFIX_model_output = modify_coef_log_model(PREFIX_model_fit, model = "lin.model", rounding = FALSE)
      PREFIX_model_output = PREFIX_model_output[2,]
    }
    PREFIX_model_output$Model.formula.string = PREFIX_Model.formula.string
    PREFIX_model_output$Gene = PREFIX_Curr_Test
    PREFIX_model_output$Contrast = paste0(PREFIX_contrast_vector, collapse = "/")
    PREFIX_model_output$Excluded_particip = paste0(PREFIX_Participants_missing_data, collapse = ";")
    PREFIX_model_output$Excluded_particip_number = length(PREFIX_Participants_missing_data)
    PREFIX_model_output$Df = PREFIX_model_fit$df.residual
    
    # Calculate logFC
    PREFIX_expression_values_controls = PREFIX_tmp_df[PREFIX_tmp_df$Case_Control == PREFIX_contrast_vector[1], ncol(PREFIX_tmp_df)]
    PREFIX_expression_values_cases = PREFIX_tmp_df[PREFIX_tmp_df$Case_Control == PREFIX_contrast_vector[2], ncol(PREFIX_tmp_df)]
    PREFIX_expression_values_controls_exp = 2^PREFIX_expression_values_controls
    PREFIX_expression_values_cases_exp = 2^PREFIX_expression_values_cases
    Mean_1 = mean(PREFIX_expression_values_controls_exp)
    Mean_2 = mean(PREFIX_expression_values_cases_exp)
    FC = Mean_2/Mean_1
    Log_FC = log2(FC)
    PREFIX_model_output$logFC = Log_FC
    
    PREFIX_model_output = PREFIX_model_output[,c("Model.formula.string", "Contrast", "Gene", "Excluded_particip", "Excluded_particip_number",
                                                 "Coef.", "logFC", "β","95% CI β","Std. Error β","T.value","P.value", "Df")]
    PREFIX_MODEL_LIST[[PREFIX_i]] = PREFIX_model_output
  }
  
  if (length(PREFIX_MODEL_LIST) > 1){
    PREFIX_MODEL_LIST = do.call(rbind, PREFIX_MODEL_LIST)
  } else {
    PREFIX_MODEL_LIST = PREFIX_MODEL_LIST[[1]]
  }
  
  if (PREFIX_show_all_coef == FALSE){
    PREFIX_MODEL_LIST$Adj.P.val =  p.adjust(p = PREFIX_MODEL_LIST$P.value, method = "bonferroni")
    PREFIX_MODEL_LIST = dplyr::arrange(PREFIX_MODEL_LIST, P.value)
  } else {
    PREFIX_MODEL_LIST$Adj.P.val = NA
  }
  colnames(PREFIX_MODEL_LIST) = stri_replace_all_fixed(colnames(PREFIX_MODEL_LIST), pattern = "PREFIX_", replacement = "")
  
  if (PREFIX_round_output){
    writeLines("Rounding Output")
    PREFIX_MODEL_LIST$β = as.numeric(PREFIX_MODEL_LIST$β)
    PREFIX_MODEL_LIST$`Std. Error β` = as.numeric(PREFIX_MODEL_LIST$`Std. Error β`)
    PREFIX_MODEL_LIST$P.value = as.numeric(PREFIX_MODEL_LIST$P.value)
    PREFIX_MODEL_LIST$`95% CI β` = sapply(PREFIX_MODEL_LIST$`95% CI β`, function(x){
      x = stri_split_fixed(x, pattern = " \u2013 ")
      x = unlist(x)
      x = x[x!=""]
      x = as.numeric(x)
      x = round(x, digits = 3)
      x = paste0(x, collapse = " \u2013 ")
    })
    PREFIX_MODEL_LIST = PREFIX_MODEL_LIST %>%
      mutate_if(is.numeric, function(x) (round(x, digits = 3)))
  }
  return(PREFIX_MODEL_LIST)
}


modify_coef_log_model = function(x, digits = 3, model = "log_regr", rounding = TRUE){
  if (model == "log_regr"){
    Summary_model = summary(x)[["coefficients"]]
    if (rounding){
      Summary_model[,-ncol(Summary_model)] = apply(Summary_model[,-ncol(Summary_model)], 2, function(x) round(x, digits = digits))
    }
    #Summary_model[,ncol(Summary_model)] = round(Summary_model[,ncol(Summary_model)], digits = 9)
    Conf_inter = confint(x, level = 0.95)
    if (rounding){
      Conf_inter = apply(Conf_inter, 2, function(x) round(x, digits = digits))
    }
    Conf_inter = mapply(Conf_inter[,1], Conf_inter[,2], FUN = function(x,y){paste0(x, " \u2013 " ,y)})
    Exponent_beta_coef = round(exp(Summary_model[,1]), digits = digits)
    Table_output = cbind(rownames(Summary_model), Summary_model[,1], Conf_inter, Exponent_beta_coef, Summary_model[,-1])
    colnames(Table_output) = c("Coef.", "\U03B2", "95% CI \U03B2", "exp(\U03B2)", "Std. Error \U03B2","Z.value", "P.value")
    rownames(Table_output) = NULL
    Table_output = as.data.frame(Table_output)
  } else {
    Summary_model = summary(x)[["coefficients"]]
    if (rounding){
      Summary_model[,-ncol(Summary_model)] = apply(Summary_model[,-ncol(Summary_model)], 2, function(x) round(x, digits = digits))
    }
    #Summary_model[,ncol(Summary_model)] = round(Summary_model[,ncol(Summary_model)], digits = 9)
    Conf_inter = confint(x, level = 0.95)
    if (rounding){
      Conf_inter = apply(Conf_inter, 2, function(x) round(x, digits = digits))
    }
    Conf_inter = mapply(Conf_inter[,1], Conf_inter[,2], FUN = function(x,y){paste0(x, " \u2013 " ,y)})
    Table_output = cbind(rownames(Summary_model), Summary_model[,1], Conf_inter, Summary_model[,-1])
    colnames(Table_output) = c("Coef.", "\U03B2", "95% CI \U03B2", "Std. Error \U03B2", "T.value", "P.value")
    rownames(Table_output) = NULL
    Table_output = as.data.frame(Table_output)
  }
  return(Table_output)
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


# SCAS score calculator
SCAS_calculator = function(df_SCAS_ordered, NA_param = FALSE){
  
  sum_function = function(x){
    x = sum(x, na.rm = NA_param)
    return(x)
  }
  
  Total_score = apply(df_SCAS_ordered[,-c(11, 17, 26, 31, 38, 43)], MARGIN = 1, sum_function)
  
  Social_phobia = apply(df_SCAS_ordered[,c(6, 7, 9, 10, 29, 35)], MARGIN = 1, sum_function)
  
  GAD = apply(df_SCAS_ordered[,c(1, 3, 4, 20, 22, 24)], MARGIN = 1, sum_function)
  
  FOPH = apply(df_SCAS_ordered[,c(2, 18, 23, 25, 33)], MARGIN = 1, sum_function)
  
  Panic_agoraphobia = apply(df_SCAS_ordered[,c(13, 21, 28, 30, 32, 34, 36, 37, 39)], MARGIN = 1, sum_function)
  
  OCD = apply(df_SCAS_ordered[,c(14, 19, 27, 40, 41, 42)], MARGIN = 1, sum_function)
  
  SAD = apply(df_SCAS_ordered[,c(5, 8, 12, 15, 16, 44)], MARGIN = 1, sum_function)
  
  Score_df = cbind(Total_score, Social_phobia, GAD, FOPH, Panic_agoraphobia, OCD, SAD)
  
  colnames(Score_df) = paste0("SCAS_", colnames(Score_df))
  
  Score_df = as.data.frame(Score_df)
  
  return(Score_df)
  
}

# plot forest for genes
plot_forest_meta_gene = function(gene_name, meta_df){
  
  tmp_df_genes = meta_df[meta_df$Corrected_symbol == gene_name,]
  
  if(nrow(tmp_df_genes)<1){
    writeLines("Gene is missing in the data")
    return(NULL)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  # Example Visualization
  weights = fmtx(weights(tmp_meta_model), digits=1)
  sav = forest(tmp_meta_model, slab = Study, ilab = weights)
  k = nrow(tmp_df_genes)
  colp <- "red"
  segments(coef(tmp_meta_model), 0, coef(tmp_meta_model), k, col=colp, lty="33", lwd=0.8)
  
  # Add text
  par(xpd=NA)
  par(cex=sav$cex, font=2)
  # Headers
  text(sav$xlim[1], k+2.5, pos=4, "Cohort")
  text(-2, k+2.5, pos=4, "Weight %")
  text(0, k+2.7, "Log2FC,\n(95% CI)")
  segments(sav$ilab.xpos[1]-0.22, k+2.5, sav$ilab.xpos[2]+0.13, k+2.8)
  text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")
  
  # Use a non-bold font for the rest of the text
  par(cex=sav$cex, font=1)
  text(sav$ilab.xpos[3], 0, "100.0")
  text(sav$xlim[1], -2, pos=4, bquote(paste("Test for heterogeneity: ",
                                            tau^2, "=", .(fmtx(tmp_meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(tmp_meta_model$QE, digits=2)),
                                            ", df=", .(tmp_meta_model$k - tmp_meta_model$p), ", ",
                                            .(fmtp(tmp_meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)))))
  title(paste0("Gene: ", gene_name))
}

plot_forest_simple_gene = function(gene_name, meta_df, study_column_char = "Study"){
  tmp_df_genes = meta_df[meta_df$Corrected_symbol == gene_name,]
  if(nrow(tmp_df_genes)<1){
    writeLines("Gene is missing in the data")
    return(NULL)
  }
  forest(tmp_df_genes$avgLog2FC, tmp_df_genes$maxSE**2, slab = tmp_df_genes[,study_column_char], xlab = paste0("Gene: ", gene_name))
}

# A function to make Venn diagrams from named lists (uses ggplot2 and ggVennDiagram)
make_Venn_digram_list = function(named_list, plot_full_path = NULL, ...){
  
  # Customizable Venn diagram
  venn = Venn(named_list)
  data = process_data(venn)
  data@region$full_lable = sapply(data@region$count, function(x){
    Number = x
    Percent = x/sum(data@region$count)
    Percent = Percent*100
    Percent = round(Percent, digits = 1)
    Label_full = paste0(Number, "\n","(", Percent, "%)")
    return(Label_full)
  })
  
  # to see available shapes plot_shapes()
  plot = ggplot() +
    
    # 1. region count layer
    geom_sf(aes(fill = id), data = venn_region(data), ...) +
    
    # 2. set edge layer
    geom_sf(color="black", size = 0.5, data = venn_setedge(data), show.legend = FALSE, ...) +
    
    # 3. set label layer
    geom_sf_text(aes(label = name), data = venn_setlabel(data), ...) +
    
    # 4. region label layer
    geom_sf_label(aes(label = full_lable), data = venn_region(data), ...) +
    scale_fill_brewer(...) +
    scale_x_continuous(expand = expansion(mult = .2)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(1,1,1,1, "cm")
    )
  print(plot)
  if (is.character(plot_full_path)){
    pdf(plot_full_path, width = 10, height = 10)
    print(plot)
    dev.off()
  }
}

run_enrichment_GO_KEGG_gene_set = function(genes, universe, categories_to_show = 30, folder, plot_name_pref){
  
  # genes should be submitted as Entrez IDs
  if (length(genes) < 10) {
    Minsize = length(genes)
  } else {
    Minsize = 10
  }
  
  # IDs should be characters
  genes = as.character(genes)
  universe = as.character(universe)
  
  # Running GO analysis
  GO_String = c("BP", "MF", "CC")
  GO_Enrichment = lapply(GO_String, function(x){
    result = clusterProfiler::enrichGO(gene = genes, OrgDb = 'org.Hs.eg.db', ont = x, universe = universe, minGSSize = Minsize)
    return(result)
  })
  KEGG_Enrichment = clusterProfiler::enrichKEGG(gene = genes,  organism = 'hsa' , universe = universe, minGSSize = Minsize)
  
  # Saving GO results
  dir.create(folder)
  combined_enrich = list()
  for (i in 1:length(GO_String)){
    filename = paste0(plot_name_pref, "_GO_", GO_String[i], ".xlsx")
    filename =  paste0(folder, "/", filename)
    Curr_GO = GO_Enrichment[[i]]
    Curr_GO_readable = setReadable(Curr_GO, 'org.Hs.eg.db', 'ENTREZID')
    Curr_GO_result = Curr_GO@result
    Curr_GO_result$geneSymbol = Curr_GO_readable@result$geneID
    combined_enrich[[i]] = Curr_GO_result
    openxlsx::write.xlsx(x = Curr_GO_result, file = filename, overwrite = TRUE)
  }
  
  # Saving KEGG results
  if (!is.null(KEGG_Enrichment)){
    KEGG_Enrichment_readable = setReadable(KEGG_Enrichment, 'org.Hs.eg.db', 'ENTREZID')
    KEGG_Enrichment_result = KEGG_Enrichment@result
    KEGG_Enrichment_result$geneSymbol = KEGG_Enrichment_readable@result$geneID
    combined_enrich[[4]] = KEGG_Enrichment_result
    openxlsx::write.xlsx(x = KEGG_Enrichment_result, file = paste0(folder, "/", plot_name_pref, "_KEGG.xlsx"), overwrite = TRUE)
  } else {
    combined_enrich[[4]] = NA
  }
  # Saving full results
  openxlsx::write.xlsx(x = combined_enrich, file = paste0(folder, "/", plot_name_pref, "_enrichment_combined_sheet.xlsx"), overwrite = TRUE)
  
  # Making plots
  GO_Enrichment_readable = lapply(GO_Enrichment, function(x) setReadable(x, 'org.Hs.eg.db', 'ENTREZID'))
  
  if (!is.null(KEGG_Enrichment)){
    GO_Enrichment_readable[[4]] = KEGG_Enrichment_readable
  } else {
    GO_Enrichment_readable[[4]] = NA
  }
  
  GO_String = c("BP", "MF", "CC", "KEGG")
  gc()
  
  for (i in 1:length(GO_String)){
    Curr_GO = GO_Enrichment_readable[[i]]
    
    if (is.na(Curr_GO)){
      writeLines(paste0(filename, " has NO results and analysis was not conducted"))
      next
    }
    
    filename = paste0(plot_name_pref, "_GO_", GO_String[i])
    filename =  paste0(folder, "/", filename)
    
    if (GO_String[i] == "KEGG"){
      # Figures for KEGG
      dimensions = c(19*nrow(summary(Curr_GO))/categories_to_show, 11*nrow(summary(Curr_GO))/categories_to_show)
      
      if (any(dimensions < 1)){
        dimensions = dimensions*10
      }
      
      if (nrow(summary(Curr_GO)) < 1){
        writeLines(paste0(filename, " has NO significant results"))
        next
      }
      
      # Barplot
      pdf(paste0(filename, "_bar.pdf")
          , width = 10, height = 10)
      print(barplot(Curr_GO, showCategory = categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_bar_full.pdf")
          , width = 10, height = 10*nrow(summary(Curr_GO))/categories_to_show)
      print(barplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Dotplot
      pdf(paste0(filename, "_dot.pdf")
          , width = 10, height = 10)
      print(dotplot(Curr_GO, showCategory= categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_dot_full.pdf")
          , width = 10, height = dimensions[2])
      print(dotplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Gene-Concept Network
      pdf(paste0(filename, "_gene_conc.pdf")
          , width = 19, height = 11)
      print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                     showCategory = categories_to_show, cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
      dev.off()
      
      tryCatch({
        pdf(paste0(filename, "_gene_conc_full.pdf")
            , width = dimensions[1], height = dimensions[2])
        print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                       showCategory = nrow(summary(Curr_GO)), cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
        dev.off()
      }, error = function(e) {writeLines(paste0("Full gene conc plot is not available for ", filename))})
      
    } else {
      # Figures for GO
      dimensions = c(19*nrow(summary(Curr_GO))/categories_to_show, 11*nrow(summary(Curr_GO))/categories_to_show)
      
      if (any(dimensions < 1)){
        dimensions = dimensions*10
      }
      
      if (nrow(summary(Curr_GO)) < 1){
        writeLines(paste0(filename, " has NO significant results"))
        next
      }
      
      # Barplot
      pdf(paste0(filename, "_bar.pdf")
          , width = 10, height = 10)
      print(barplot(Curr_GO, showCategory = categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_bar_full.pdf")
          , width = 10, height = dimensions[2])
      print(barplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Dotplot
      pdf(paste0(filename, "_dot.pdf")
          , width = 10, height = 10)
      print(dotplot(Curr_GO, showCategory= categories_to_show, font.size = 10))
      dev.off()
      
      pdf(paste0(filename, "_dot_full.pdf")
          , width = 10, height = dimensions[2])
      print(dotplot(Curr_GO, showCategory = nrow(summary(Curr_GO)), font.size = 10))
      dev.off()
      
      # Gene-Concept Network
      pdf(paste0(filename, "_gene_conc.pdf")
          , width = 19, height = 11)
      print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                     showCategory = categories_to_show, cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
      dev.off()
      
      tryCatch({
        pdf(paste0(filename, "_gene_conc_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(cnetplot(Curr_GO, circular = FALSE,  colorEdge = TRUE, color_category='firebrick', color_gene='steelblue', layout = "kk",
                       showCategory = nrow(summary(Curr_GO)), cex_gene = 0.7, cex_label_gene = 0.7, shadowtext = "gene", max.overlaps = 1000, force = 3, force_pull = 0.5, max.time = 2))
        dev.off()
      }, error = function(e) {writeLines(paste0("Full gene conc plot is not available for ", filename))})
      
      # GO induced graph
      tryCatch({
        pdf(paste0(filename, "_go_graph.pdf"), 
            width = 19, height = 11)
        print(goplot(Curr_GO, showCategory = categories_to_show))
        dev.off()
      }, error = function(e) {writeLines(paste0("Go graph is not available for ", filename))})
      
      tryCatch({
        pdf(paste0(filename, "_go_graph_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(goplot(Curr_GO, showCategory = nrow(summary(Curr_GO))))
        dev.off()
      }, error = function(e) {writeLines(paste0("Go graph (big) is not available for ", filename))})
      
      # Tree plot
      Curr_GO_pairwise = enrichplot::pairwise_termsim(Curr_GO)
      dimensions = c(19*nrow(Curr_GO_pairwise@termsim)/categories_to_show, 11*nrow(Curr_GO_pairwise@termsim)/categories_to_show)
      if (any(dimensions < 1)){
        dimensions = dimensions*10
      }
      
      tryCatch({
        pdf(paste0(filename, "_tree.pdf"), 
            width = 19, height = 11)
        print(treeplot(Curr_GO_pairwise, showCategory = categories_to_show, nCluster = 10))
        dev.off()
      }, error = function(e) {writeLines(paste0("Tree plot is not available for ", filename))})
      
      tryCatch({
        pdf(paste0(filename, "_tree_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(treeplot(Curr_GO_pairwise, showCategory = nrow(Curr_GO_pairwise@termsim), nCluster = 10))
        dev.off()
      }, error = function(e) {writeLines(paste0("Tree plot (big) is not available for ", filename))})
      
      
      # Enrichment map
      tryCatch({
        pdf(paste0(filename, "_enr_map.pdf"), 
            width = 19, height = 11)
        print(emapplot(Curr_GO_pairwise, showCategory = categories_to_show, layout = "fr", cex_label_category = 0.9, 
                       cex_line = 0.8, repel = TRUE))
        dev.off()
      }, error = function(e) {writeLines(paste0("Enrich. map is not available for ", filename))})
      
      tryCatch({
        pdf(paste0(filename, "_enr_map_full.pdf"), 
            width = dimensions[1], height = dimensions[2])
        print(emapplot(Curr_GO_pairwise, showCategory = nrow(Curr_GO_pairwise@termsim), layout = "fr", cex_label_category = 0.9, 
                       cex_line = 0.8, repel = TRUE))
        dev.off()
      }, error = function(e) {writeLines(paste0("Enrich. map is not available for ", filename))})
      
    }
  }
  
  # Preparing outputs to load into the environment
  GO_Enrichment[[4]] = KEGG_Enrichment
  
  if (is.null(KEGG_Enrichment)){
    GO_Enrichment[[4]] = NA
  }
  
  Output = list()
  Output[[1]] = GO_Enrichment
  Output[[2]] = GO_Enrichment_readable
  Output[[3]] = combined_enrich
  names(Output) = c("Enrichment", "Mapped Enrichment", "Combined Results")
  
  for (i in 1:length(Output)){
    names(Output[[i]]) = GO_String
  }
  
  return(Output)
}

run_enrichment_GO_KEGG_gene_set_only_tables = function(genes, universe, folder, plot_name_pref){
  
  # genes should be submitted as Entrez IDs
  if (length(genes) < 10) {
    Minsize = length(genes)
  } else {
    Minsize = 10
  }
  
  # IDs should be characters
  genes = as.character(genes)
  universe = as.character(universe)
  
  # Running GO analysis
  GO_String = c("BP", "MF", "CC")
  GO_Enrichment = lapply(GO_String, function(x){
    result = clusterProfiler::enrichGO(gene = genes, OrgDb = 'org.Hs.eg.db', ont = x, universe = universe, minGSSize = Minsize)
    return(result)
  })
  KEGG_Enrichment = clusterProfiler::enrichKEGG(gene = genes,  organism = 'hsa' , universe = universe, minGSSize = Minsize)
  
  # Saving GO results
  dir.create(folder)
  combined_enrich = list()
  for (i in 1:length(GO_String)){
    filename = paste0(plot_name_pref, "_GO_", GO_String[i], ".xlsx")
    filename =  paste0(folder, "/", filename)
    Curr_GO = GO_Enrichment[[i]]
    Curr_GO_readable = setReadable(Curr_GO, 'org.Hs.eg.db', 'ENTREZID')
    Curr_GO_result = Curr_GO@result
    Curr_GO_result$geneSymbol = Curr_GO_readable@result$geneID
    combined_enrich[[i]] = Curr_GO_result
    openxlsx::write.xlsx(x = Curr_GO_result, file = filename, overwrite = TRUE)
  }
  
  # Saving KEGG results
  if (!is.null(KEGG_Enrichment)){
    KEGG_Enrichment_readable = setReadable(KEGG_Enrichment, 'org.Hs.eg.db', 'ENTREZID')
    KEGG_Enrichment_result = KEGG_Enrichment@result
    KEGG_Enrichment_result$geneSymbol = KEGG_Enrichment_readable@result$geneID
    combined_enrich[[4]] = KEGG_Enrichment_result
    openxlsx::write.xlsx(x = KEGG_Enrichment_result, file = paste0(folder, "/", plot_name_pref, "_KEGG.xlsx"), overwrite = TRUE)
  } else {
    combined_enrich[[4]] = NA
  }
  # Saving full results
  openxlsx::write.xlsx(x = combined_enrich, file = paste0(folder, "/", plot_name_pref, "_enrichment_combined_sheet.xlsx"), overwrite = TRUE)
  
  GO_Enrichment_readable = lapply(GO_Enrichment, function(x) setReadable(x, 'org.Hs.eg.db', 'ENTREZID'))
  
  
  # Preparing outputs to load into the environment
  GO_Enrichment[[4]] = KEGG_Enrichment
  
  if (is.null(KEGG_Enrichment)){
    GO_Enrichment[[4]] = NA
  }
  
  Output = list()
  Output[[1]] = GO_Enrichment
  Output[[2]] = GO_Enrichment_readable
  Output[[3]] = combined_enrich
  names(Output) = c("Enrichment", "Mapped Enrichment", "Combined Results")
  
  for (i in 1:length(Output)){
    names(Output[[i]]) = GO_String
  }
  
  return(Output)
}

################### PSY Recall phenotype preprocessing ###################
DAWBA_DS = read.xlsx("DAWBA PSY June 2021.xlsx", sheet = 1)


Recall_2_Pheno = list.files()
Recall_2_Pheno = Recall_2_Pheno[stri_detect_fixed(Recall_2_Pheno, "RECALL_2")]
Recall_2_Pheno = lapply(Recall_2_Pheno, function(x) read.xlsx(x, sheet = 1, startRow = 1, detectDates = TRUE))
which(colnames(Recall_2_Pheno[[1]]) == "Medication")
which(colnames(Recall_2_Pheno[[2]]) == "Medication")
names(Recall_2_Pheno) = c("Recall_2A", "Recall_2B")
Columns_Recall_2A = c(1, 4, 5, 6, 7, 8:16, 17, 18,20, 21:24, 25:33, 111:132, 168:178, 179:226, 152:167) # HADS - 152:167
Columns_Recall_2B = c(1:30, 31:52, 53:63, 64:111)
Recall_2A_Pheno = Recall_2_Pheno[[1]]
Recall_2A_Pheno = Recall_2A_Pheno[,Columns_Recall_2A]
Recall_2B_Pheno = Recall_2_Pheno[[2]]
Recall_2B_Pheno = Recall_2B_Pheno[,Columns_Recall_2B]
Recall_2B_Pheno = Recall_2B_Pheno[!is.na(Recall_2B_Pheno$Code),]

############ Recall 2A
Recall_2A_Pheno_small = Recall_2A_Pheno
Columns_to_select_RC2A = c(1,2, 22:26, 30, 31:51, 52, 53:63, 64:108, 111, 112, 113)
Recall_2A_Pheno_small = Recall_2A_Pheno_small[,Columns_to_select_RC2A]
Recall_2A_Pheno_small$DAWBA_DEPBAND = sapply(Recall_2A_Pheno_small$DAWBA.ID, function(x){
  index = which(DAWBA_DS$dawbaID == x)
  output = DAWBA_DS$depband[index]
  if (length(output)<1){return(NA)}
  return(output)
})
SCAS_DF_RC2A = Recall_2A_Pheno_small[,42:86]
SCAS_DF_RC2A = apply(SCAS_DF_RC2A, 2, as.numeric)

# corrupt total score check
corrupt_total_score = which(apply(SCAS_DF_RC2A[,-c(11, 17, 26, 31, 38, 43)], 1, function(x) any(is.na(x))))
corrupt_total_score
# Rows 2 199 200 336 should be technically excluded

corrupt_total_score_item_count = apply(SCAS_DF_RC2A[,-c(11, 17, 26, 31, 38, 43)], 1, function(x) sum(is.na(x)))
corrupt_total_score_item_count[corrupt_total_score_item_count>0]
# All rows just miss 1 item
SCAS_scores_RC2A = SCAS_calculator(df_SCAS_ordered = SCAS_DF_RC2A, NA_param = TRUE)
SCAS_scores_RC2A[corrupt_total_score, ]
Recall_2A_Pheno_small = cbind(Recall_2A_Pheno_small, SCAS_scores_RC2A)
Recall_2A_Pheno_small$Code[2] # "PSY0021"

# item check
item_numbers = sapply(colnames(SCAS_DF_RC2A), function(x) unlist(stri_split_fixed(x, pattern = "."))[2])
names(item_numbers) = NULL
is.unsorted(as.numeric(item_numbers)) # FALSE


############ Recall 2B
Recall_2B_Pheno_small = Recall_2B_Pheno
Columns_to_select_RC2B = c(1, 2, 22:26, 30, 31:51, 52, 53:63, 64:108, 111)
Recall_2B_Pheno_small = Recall_2B_Pheno_small[,Columns_to_select_RC2B]
Recall_2B_Pheno_small$DAWBA_DEPBAND = sapply(Recall_2B_Pheno_small$DAWBA.ID, function(x){
  index = which(DAWBA_DS$dawbaID == x)
  output = DAWBA_DS$depband[index]
  if (length(output)<1){return(NA)}
  return(output)
})
SCAS_DF_RC2B = Recall_2B_Pheno_small[,42:86]
SCAS_DF_RC2B = apply(SCAS_DF_RC2B, 2, as.numeric)

# corrupt total score check
corrupt_total_score = which(apply(SCAS_DF_RC2B[,-c(11, 17, 26, 31, 38, 43)], 1, function(x) any(is.na(x))))
corrupt_total_score
# ALL rows are OK!

SCAS_scores_RC2B = SCAS_calculator(df_SCAS_ordered = Recall_2B_Pheno_small[,42:86])
Recall_2B_Pheno_small = cbind(Recall_2B_Pheno_small, SCAS_scores_RC2B)

# item check
item_numbers = sapply(colnames(SCAS_DF_RC2B), function(x) unlist(stri_split_fixed(x, pattern = "."))[2])
names(item_numbers) = NULL
is.unsorted(as.numeric(item_numbers)) # FALSE

############ Saving phenotypes
write.csv(Recall_2A_Pheno_small, "Recall_2A_Pheno_small.csv", row.names = FALSE)
write.csv(Recall_2B_Pheno_small, "Recall_2B_Pheno_small.csv", row.names = FALSE)

################### PSY Recall phenotype curation ###################
PSY_PHENO_RC2A = smart_fread("Recall_2A_Pheno_small.csv")
PSY_PHENO_RC2B = smart_fread("Recall_2B_Pheno_small.csv")
intersecting_cols_RC2 = intersect(colnames(PSY_PHENO_RC2A), colnames(PSY_PHENO_RC2B))
PSY_PHENO_RC_TOTAL = rbind(PSY_PHENO_RC2A[,intersecting_cols_RC2], PSY_PHENO_RC2B[,intersecting_cols_RC2])

# Curating recall
PSY_PHENO_RC_TOTAL$Code = toupper(PSY_PHENO_RC_TOTAL$Code)
PSY_PHENO_RC_TOTAL$`Gender.(0)male.(1)female` = sapply(PSY_PHENO_RC_TOTAL$`Gender.(0)male.(1)female`, function(x){
  if (x == 0){
    return("M")
  } else {
    return("W")
  }
})
colnames(PSY_PHENO_RC_TOTAL) = stri_replace_all_fixed(colnames(PSY_PHENO_RC_TOTAL), pattern = "Gender.(0)male.(1)female", replacement = "Gender")
############ Saving phenotypes
write.csv(PSY_PHENO_RC_TOTAL, "PSY_PHENO_RC_TOTAL.csv", row.names = FALSE)

################### Curating phenotype medication intake and removing duplicates ###################
ANXIETY_CURATED_PSY_PHENO_RC_TOTAL = openxlsx::read.xlsx("ANXIETY_PSY_PHENO_RC_TOTAL_CURATED.xlsx", sheet = 1)
CURATED_PSY_PHENO_RC_TOTAL = PSY_PHENO_RC_TOTAL
CURATED_PSY_PHENO_RC_TOTAL$Code = toupper(CURATED_PSY_PHENO_RC_TOTAL$Code)

any(duplicated(CURATED_PSY_PHENO_RC_TOTAL$Code)) # TRUE There are duplicates in phenotypes
CURATED_PSY_PHENO_RC_TOTAL$Code[duplicated(CURATED_PSY_PHENO_RC_TOTAL$Code)] # "PSY0662" "PSY0867"

duplicated_recall = CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code %in%  c("PSY0662", "PSY0867"),]
duplicated_recall_orig = Recall_2A_Pheno[Recall_2A_Pheno$Code %in%  c("PSY0662", "PSY0867"),]

# "PSY0662" should be female as indicated in screening "Collected data_Screening.xlsx", height=159
as.Date(43158, origin = "1899-12-30") # "2018-02-27", timestamp from PSY_recall_1_2_Biobank_masterfile - Updated_Aleksandr_Sokolov_Oct2021.xlsx
# "PSY0867" # 02/27/2018 is the correct date -> Age 16
duplicated_recall_selected = duplicated_recall[c(1,3),]

# curating recall dataset
CURATED_PSY_PHENO_RC_TOTAL = CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code %!in% duplicated_recall$Code,]
CURATED_PSY_PHENO_RC_TOTAL = rbind(CURATED_PSY_PHENO_RC_TOTAL, duplicated_recall_selected)

# curating recall anxiety dataset
duplicated_recall = ANXIETY_CURATED_PSY_PHENO_RC_TOTAL[ANXIETY_CURATED_PSY_PHENO_RC_TOTAL$Code %in%  c("PSY0662", "PSY0867"),]
duplicated_recall_selected = duplicated_recall[c(1,3),]
ANXIETY_CURATED_PSY_PHENO_RC_TOTAL = ANXIETY_CURATED_PSY_PHENO_RC_TOTAL[ANXIETY_CURATED_PSY_PHENO_RC_TOTAL$Code %!in% duplicated_recall$Code,]
ANXIETY_CURATED_PSY_PHENO_RC_TOTAL = rbind(ANXIETY_CURATED_PSY_PHENO_RC_TOTAL, duplicated_recall_selected)


# check
any(duplicated(CURATED_PSY_PHENO_RC_TOTAL$Code)) # FALSE
any(duplicated(ANXIETY_CURATED_PSY_PHENO_RC_TOTAL$Code)) # FALSE


CURATED_PSY_PHENO_RC_TOTAL$Antidepr_anxiety = sapply(CURATED_PSY_PHENO_RC_TOTAL$Code, function(x){
  if (x %in% ANXIETY_CURATED_PSY_PHENO_RC_TOTAL$Code){
    treatment_status = ANXIETY_CURATED_PSY_PHENO_RC_TOTAL[ANXIETY_CURATED_PSY_PHENO_RC_TOTAL$Code == x, "Antidepr_anxiety"]
    treatment_status = as.numeric(treatment_status)
  } else {
    treatment_status = NA
  }
  return(treatment_status)
})

# PSY0021
CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code == "PSY0021", "Antidepr_anxiety"] = 0


################### Initial Variable prep ###################
any(is.na(CURATED_PSY_PHENO_RC_TOTAL$SUAS.Total.score)) # FALSE


CURATED_PSY_PHENO_RC_TOTAL$Gender = factor(CURATED_PSY_PHENO_RC_TOTAL$Gender, levels = c("W", "M"))
CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary = ifelse(CURATED_PSY_PHENO_RC_TOTAL$SUAS.Total.score >= 31,
                                                          "Risk group", "Normal")
CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary = factor(CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary, levels = c("Normal", "Risk group"))

# Variables
CURATED_PSY_PHENO_RC_TOTAL$Gender
CURATED_PSY_PHENO_RC_TOTAL$Age
CURATED_PSY_PHENO_RC_TOTAL$BMI
CURATED_PSY_PHENO_RC_TOTAL$SUAS.Total.score
CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary


################### Importing proteomics data ###################
# Preprocess OLINK function
preprocess_OLINK_main = function(PRFX_OLINK_DS){
  # modifying original DS
  PRFX_idx = which(PRFX_OLINK_DS[,1] == "Panel")
  PRFX_OLINK_DS = PRFX_OLINK_DS[PRFX_idx:nrow(PRFX_OLINK_DS),]
  PRFX_OLINK_DS = as.data.frame(t(PRFX_OLINK_DS))
  colnames(PRFX_OLINK_DS) = PRFX_OLINK_DS[1,]
  PRFX_OLINK_DS = PRFX_OLINK_DS[-1,]
  rownames(PRFX_OLINK_DS) = PRFX_OLINK_DS$Assay
  
  #preparing data
  PRFX_probes = PRFX_OLINK_DS[,c(1:4, (ncol(PRFX_OLINK_DS)-2):ncol(PRFX_OLINK_DS))]
  PRFX_probes = PRFX_probes[-((nrow(PRFX_OLINK_DS)-1):nrow(PRFX_OLINK_DS)),]
  PRFX_probes$`Missing Data freq.` = as.numeric(PRFX_probes$`Missing Data freq.`)
  PRFX_probes$LOD = as.numeric(PRFX_probes$LOD)
  
  PRFX_sample_status = PRFX_OLINK_DS[(nrow(PRFX_OLINK_DS)-1):nrow(PRFX_OLINK_DS), ]
  PRFX_sample_status = PRFX_sample_status[,colnames(PRFX_sample_status) %!in% colnames(PRFX_probes)]
  PRFX_sample_status = as.data.frame(t(PRFX_sample_status))
  PRFX_sample_status$Sample = rownames(PRFX_sample_status)
  
  PRFX_expression = PRFX_OLINK_DS[, colnames(PRFX_OLINK_DS) %in% PRFX_sample_status$Sample]
  PRFX_expression = PRFX_expression[-((nrow(PRFX_OLINK_DS)-1):nrow(PRFX_OLINK_DS)),]
  PRFX_expression = as.matrix(PRFX_expression)
  PRFX_matrix_names = rownames(PRFX_expression)
  PRFX_expression = apply(PRFX_expression, 2, as.numeric)
  rownames(PRFX_expression) = PRFX_matrix_names
  
  #filtering probes and samples 
  PRFX_sample_status$QC_result = ifelse(PRFX_sample_status$`QC Warning` == "Pass", "Included", "Excluded")
  PRFX_expression_pass = PRFX_expression[, colnames(PRFX_expression) %in% PRFX_sample_status[PRFX_sample_status$`QC Warning` == "Pass", "Sample"]]
  
  PRFX_probes$Probe_QC = ifelse(PRFX_probes$`Missing Data freq.` > 0.2, "Excluded", "Pass")
  for (PRFX_I in 1 : nrow(PRFX_probes)){
    if (PRFX_probes$LOD[PRFX_I] > 2.5 & PRFX_probes$Probe_QC[PRFX_I] == "Excluded"){
      PRFX_probes$Probe_QC[PRFX_I] = "High LOD"
    }
  }
  PRFX_passed_probes = PRFX_probes[PRFX_probes$Probe_QC == "Pass", "Assay"]
  PRFX_expression_pass = PRFX_expression_pass[rownames(PRFX_expression_pass) %in% PRFX_passed_probes,]
  PRFX_expression_pass_strict = PRFX_expression_pass
  
  for (PRFX_I in 1 : nrow(PRFX_expression_pass_strict)){
    PRFX_Curr_Probe = rownames(PRFX_expression_pass_strict)[PRFX_I]
    PRFX_Curr_LOD = PRFX_probes[PRFX_probes$Assay == PRFX_Curr_Probe, "LOD"]
    PRFX_Curr_row = PRFX_expression_pass_strict[PRFX_I,]
    PRFX_Curr_row = ifelse(PRFX_Curr_row >= PRFX_Curr_LOD, PRFX_Curr_row, NA)
    PRFX_expression_pass_strict[PRFX_I,] = PRFX_Curr_row
  }
  
  PRFX_expression_pass = as.data.frame(PRFX_expression_pass)
  PRFX_expression_pass_strict = as.data.frame(PRFX_expression_pass_strict)
  PRFX_Output = list(Expression_matrix = PRFX_expression_pass,
                     Expression_matrix_strict = PRFX_expression_pass_strict,
                     Probes = PRFX_probes,
                     Samples = PRFX_sample_status)
  return(PRFX_Output)
}

# Batch 2021
OLINK_DS_2021 = openxlsx::read.xlsx("OLINK_DIANA/Proteomics in healthy adolescents_PSY_NPX_no duplicates.xlsx")
OLINK_DS_2021_prep = preprocess_OLINK_main(PRFX_OLINK_DS = OLINK_DS_2021)
View(OLINK_DS_2021_prep$Samples)
View(OLINK_DS_2021_prep$Probes)
table(OLINK_DS_2021_prep$Samples$QC_result) # 4 - excluded, 313 - included
OLINK_DS_2021_prep$Samples[OLINK_DS_2021_prep$Samples$QC_result == "Excluded", ]
table(OLINK_DS_2021_prep$Probes$Probe_QC)  # 43 - excluded, 49 - passed
OLINK_DS_2021_EXPR = OLINK_DS_2021_prep$Expression_matrix

# Batch 2022
OLINK_DS_2022 = openxlsx::read.xlsx("OLINK_AS/Proteomic analyses_Plasma_NPX.xlsx")
OLINK_DS_2022_prep = preprocess_OLINK_main(PRFX_OLINK_DS = OLINK_DS_2022)
View(OLINK_DS_2022_prep$Samples)
View(OLINK_DS_2022_prep$Probes)
table(OLINK_DS_2022_prep$Samples$QC_result) # 5 - excluded, 233 - included
OLINK_DS_2022_prep$Samples[OLINK_DS_2022_prep$Samples$QC_result == "Excluded", ]
table(OLINK_DS_2022_prep$Probes$Probe_QC) # 39 - excluded, 10 - High LOD, 43 - passed
OLINK_DS_2022_EXPR = OLINK_DS_2022_prep$Expression_matrix

# Check for repeated measurement
table(colnames(OLINK_DS_2022_EXPR) %in% colnames(OLINK_DS_2021_EXPR)) # 8 samples were analyzed the 2nd time
over_analyzed_samples = colnames(OLINK_DS_2022_EXPR)[colnames(OLINK_DS_2022_EXPR) %in% colnames(OLINK_DS_2021_EXPR)]
OLINK_DS_2022_EXPR = OLINK_DS_2022_EXPR[,colnames(OLINK_DS_2022_EXPR) %!in% over_analyzed_samples]

# Updating phenotype information
ALL_PSY_OLINK = c(OLINK_DS_2021_prep$Samples$Sample, OLINK_DS_2022_prep$Samples$Sample)
ALL_PSY_OLINK = unique(ALL_PSY_OLINK)

CURATED_PSY_PHENO_RC_TOTAL$OLINK_STATUS = NA
CURATED_PSY_PHENO_RC_TOTAL$OLINK_BATCH = NA
CURATED_PSY_PHENO_RC_TOTAL$OLINK_PLATE = NA

for (i in 1:length(ALL_PSY_OLINK)){
  
  PRF_sample = ALL_PSY_OLINK[i]
  PRF_name_frag = unlist(stri_split_fixed(PRF_sample, pattern = "_"))
  PRF_name_frag = unlist(stri_split_fixed(PRF_name_frag, pattern = "-"))
  
  # Meta data
  PRF_praticip = PRF_name_frag[1]
  PRF_visit = PRF_name_frag[2]
  
  # Get Batch info
  if (PRF_sample %in% OLINK_DS_2021_prep$Samples$Sample){
    PRF_OLINK_BATCH = "2021"
    PRF_SAMPLE_DF = OLINK_DS_2021_prep$Samples
  } else if (PRF_sample %in% OLINK_DS_2022_prep$Samples$Sample){
    PRF_OLINK_BATCH = "2022"
    PRF_SAMPLE_DF = OLINK_DS_2022_prep$Samples
  }
  
  # Get Plate Info and Status
  PRF_OLINK_STATUS = PRF_SAMPLE_DF[PRF_SAMPLE_DF$Sample == PRF_sample, "QC_result"]
  PRF_OLINK_PLATE = PRF_SAMPLE_DF[PRF_SAMPLE_DF$Sample == PRF_sample, "Plate ID"]
  
  if (stri_detect_fixed(PRF_visit, pattern = "2")){
    # Data from recall
    CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code == PRF_praticip, "OLINK_STATUS"] = PRF_OLINK_STATUS
    CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code == PRF_praticip, "OLINK_BATCH"] = PRF_OLINK_BATCH
    CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code == PRF_praticip, "OLINK_PLATE"] = PRF_OLINK_PLATE
  }
  
  # Cleaning
  rm(list = c("PRF_sample", 
              "PRF_name_frag",
              "PRF_praticip",
              "PRF_visit", 
              "PRF_SAMPLE_DF",
              "PRF_OLINK_STATUS", 
              "PRF_OLINK_BATCH", 
              "PRF_OLINK_PLATE"))
  
}

table(CURATED_PSY_PHENO_RC_TOTAL$OLINK_STATUS)
table(CURATED_PSY_PHENO_RC_TOTAL$OLINK_BATCH)
table(CURATED_PSY_PHENO_RC_TOTAL$OLINK_PLATE)

################### Participant Selection ###################
colnames(OLINK_DS_2022_EXPR) = toupper(colnames(OLINK_DS_2022_EXPR))
colnames(OLINK_DS_2021_EXPR) = toupper(colnames(OLINK_DS_2021_EXPR))
colnames(OLINK_DS_2022_EXPR) = stri_replace_all_fixed(colnames(OLINK_DS_2022_EXPR), pattern =  "-", replacement = "_")
colnames(OLINK_DS_2021_EXPR) = stri_replace_all_fixed(colnames(OLINK_DS_2021_EXPR), pattern =  "-", replacement = "_")

combined_cols = c(colnames(OLINK_DS_2022_EXPR), colnames(OLINK_DS_2021_EXPR))
combined_cols = toupper(combined_cols)
any(duplicated(combined_cols)) # FALSE -> no duplicated tubes
combined_cols = stri_replace_all_fixed(combined_cols, pattern =  "-", replacement = "_")

recall_samples = combined_cols[multiple_stri_detector(combined_cols, pattern_vector = c("_2_", "_2B_", "_2A_"))]
recall_particip = sapply(recall_samples, function(x) unlist(stri_split_fixed(x, pattern = "_"))[1])
PSY_PHENO_RECALL_OLINK_VALID = CURATED_PSY_PHENO_RC_TOTAL[CURATED_PSY_PHENO_RC_TOTAL$Code %in% recall_particip, ]

# Filtering Samples
OLINK_DS_2021_EXPR_KEEP = OLINK_DS_2021_EXPR[, colnames(OLINK_DS_2021_EXPR) %in% recall_samples]
OLINK_DS_2022_EXPR_KEEP = OLINK_DS_2022_EXPR[, colnames(OLINK_DS_2022_EXPR) %in% recall_samples]
all(recall_samples %in% c(colnames(OLINK_DS_2021_EXPR), colnames(OLINK_DS_2022_EXPR))) # TRUE

# Check for duplicated individuals
total_samples_keep_individ = sapply(recall_samples, function(x) unlist(stri_split_fixed(x, pattern = "_"))[1])
any(duplicated(total_samples_keep_individ)) # -> FALSE, no duplicated persons

# Filtering phenotypes
any(duplicated(PSY_PHENO_RECALL_OLINK_VALID$Code)) # FALSE There are NO duplicates in phenotypes
PSY_PHENO_RECALL_OLINK_FINAL_curated = PSY_PHENO_RECALL_OLINK_VALID

# QC
any(duplicated(PSY_PHENO_RECALL_OLINK_FINAL_curated$Code)) # FALSE

# Finding overlapping proteins
proteins_passed_2021 = rownames(OLINK_DS_2021_EXPR_KEEP)
proteins_passed_2022 = rownames(OLINK_DS_2022_EXPR_KEEP)
proteins_passed_both_batches = intersect(proteins_passed_2021, proteins_passed_2022)

# Obtaining combined proteomic dataset
OLINK_DS_2021_EXPR_KEEP_both = OLINK_DS_2021_EXPR_KEEP[proteins_passed_both_batches,]
OLINK_DS_2022_EXPR_KEEP_both = OLINK_DS_2022_EXPR_KEEP[proteins_passed_both_batches,]
all(rownames(OLINK_DS_2021_EXPR_KEEP_both) == rownames(OLINK_DS_2022_EXPR_KEEP_both))  # -> TRUE
OLINK_DS_BOTH_BATCHES_FINAL = cbind(OLINK_DS_2021_EXPR_KEEP_both, OLINK_DS_2022_EXPR_KEEP_both)

# Adding Sample column
PSY_PHENO_RECALL_OLINK_FINAL_curated$Matching_sample = sapply(PSY_PHENO_RECALL_OLINK_FINAL_curated$Code, function(x){
  sample = colnames(OLINK_DS_BOTH_BATCHES_FINAL)[stri_detect_fixed(str = colnames(OLINK_DS_BOTH_BATCHES_FINAL), pattern = x)]
  return(sample)
})

# QC 
any(duplicated(colnames(OLINK_DS_BOTH_BATCHES_FINAL))) # FALSE
any(duplicated(colnames(PSY_PHENO_RECALL_OLINK_FINAL_curated$Matching_sample))) # FALSE
# uncovered participant in phenotype
colnames(OLINK_DS_BOTH_BATCHES_FINAL)[colnames(OLINK_DS_BOTH_BATCHES_FINAL) %!in% PSY_PHENO_RECALL_OLINK_FINAL_curated$Matching_sample]
# all participants are covered!

PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE = as.factor(PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE)
table(PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE)

# Reordering proteomics
OLINK_DS_BOTH_BATCHES_FINAL = OLINK_DS_BOTH_BATCHES_FINAL[,PSY_PHENO_RECALL_OLINK_FINAL_curated$Matching_sample]
all(PSY_PHENO_RECALL_OLINK_FINAL_curated$Matching_sample == colnames(OLINK_DS_BOTH_BATCHES_FINAL)) # -> TRUE, order is matching

# Checking missing values in the data
corrupt_samples_proteomics = apply(OLINK_DS_BOTH_BATCHES_FINAL, 2, function(x) any(is.na(x)))
OLINK_DS_BOTH_BATCHES_FINAL = OLINK_DS_BOTH_BATCHES_FINAL[,!corrupt_samples_proteomics]
PSY_PHENO_RECALL_OLINK_FINAL_curated = PSY_PHENO_RECALL_OLINK_FINAL_curated[!corrupt_samples_proteomics, ]

all(colnames(OLINK_DS_BOTH_BATCHES_FINAL) == PSY_PHENO_RECALL_OLINK_FINAL_curated$Matching_sample) # TRUE

######## Saving processed data
write.csv(OLINK_DS_BOTH_BATCHES_FINAL, "OLINK_DS_BOTH_BATCHES_FINAL.csv", row.names = TRUE) # 404 samples, 43 proteins
write.csv(PSY_PHENO_RECALL_OLINK_FINAL_curated, "PSY_PHENO_RECALL_OLINK_FINAL_curated.csv", row.names = FALSE)


################### SUAS cut-off cpecs ###################

CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary = ifelse(CURATED_PSY_PHENO_RC_TOTAL$SUAS.Total.score >= 29,
                                                "Risk group", "Normal")
CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary = factor(CURATED_PSY_PHENO_RC_TOTAL$SUAS_binary, levels = c("Normal", "Risk group"))


################### Adjusting for batch covariates ###################

########### Raw Adjustment (no Covariates)
OLINK_DS_BOTH_BATCHES_FINAL_ADJUSTED = ComBat(OLINK_DS_BOTH_BATCHES_FINAL, batch = PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE)
OLINK_DS_BOTH_BATCHES_FINAL_ADJUSTED = as.data.frame(OLINK_DS_BOTH_BATCHES_FINAL_ADJUSTED)


########### Adjustment with covar (no antidepr)

# Batch effect
# Participants must be columns
model_batch = model.matrix(~ SUAS_binary + Gender + Age, data = PSY_PHENO_RECALL_OLINK_FINAL_curated)
OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ADJ = ComBat(dat=OLINK_DS_BOTH_BATCHES_FINAL, batch = PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE, mod=model_batch)
OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ADJ = as.data.frame(OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ADJ)

########### Adjustment with covar (including antidepr)

# Selecting participants that all have SUAS score and information on antidepressants
PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR = PSY_PHENO_RECALL_OLINK_FINAL_curated[!is.na(PSY_PHENO_RECALL_OLINK_FINAL_curated$Antidepr_anxiety),]
OLINK_DS_BOTH_BATCHES_FINAL_ANTIDEPR = OLINK_DS_BOTH_BATCHES_FINAL[,PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$Matching_sample]
all(PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$Matching_sample == colnames(OLINK_DS_BOTH_BATCHES_FINAL_ANTIDEPR)) # -> TRUE, order is matching
PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$Antidepr_anxiety = factor(PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$Antidepr_anxiety, levels = c(0,1))

# Batch effect
# Participants must be columns
PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$OLINK_PLATE = factor(PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$OLINK_PLATE)
model_batch = model.matrix(~ SUAS_binary + Antidepr_anxiety + Gender + Age, data = PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR)

test_table = model.matrix(~ SUAS_binary + Antidepr_anxiety + Gender + Age + OLINK_PLATE, data = PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR)
test_table = cor(test_table, method = "pearson", use = "complete.obs")

OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ANTIDEPR_ADJ = ComBat(dat=OLINK_DS_BOTH_BATCHES_FINAL_ANTIDEPR, batch = PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR$OLINK_PLATE, mod=model_batch)
OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ANTIDEPR_ADJ = as.data.frame(OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ANTIDEPR_ADJ)

################### Limma differential expression ###################
PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE = factor(PSY_PHENO_RECALL_OLINK_FINAL_curated$OLINK_PLATE)

# Analysis with limma
design.matrix = model.matrix(~ SUAS_binary + Gender + Age, data = PSY_PHENO_RECALL_OLINK_FINAL_curated)
fit = lmFit(OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ADJ, design.matrix)
fitE = eBayes(fit, robust = TRUE)
Top_table = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf)
Top_table$Gene = rownames(Top_table) # RBKS CCL27 AKT1S1

################### Limma differential expression (antidepress) ###################

# Analysis with limma
design.matrix = model.matrix(~ SUAS_binary + Antidepr_anxiety + Gender + Age, data = PSY_PHENO_RECALL_OLINK_FINAL_ANTIDEPR)
fit = lmFit(OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ANTIDEPR_ADJ, design.matrix)
fitE = eBayes(fit, robust = TRUE)
Top_table_antidepr = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf)
Top_table_antidepr$Gene = rownames(Top_table_antidepr) # PLA2G10 

# Boxplot for CCL27 in PSY
probe = "CCL27"
plot_df = PSY_PHENO_RECALL_OLINK_FINAL_curated

probe_vals = OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ADJ[probe,plot_df$Matching_sample]
all(colnames(probe_vals) == plot_df$Matching_sample)
probe_vals = as.numeric(probe_vals)
plot_df = cbind(plot_df, probe_vals)

plot = ggplot(data = plot_df, aes(x = SUAS_binary, y = probe_vals, fill = SUAS_binary)) +
  geom_boxplot(notchwidth = 0.25) +
  scale_fill_manual(values = c("#C44652","#EAC151")) +
  geom_violin(alpha = 0.25) +
  labs(y = "NPX: CCL27", x = "SUAS risk group", fill = "Disease status") +
  ggthemes::theme_clean(base_size = 14)
plot
ggsave(file = "CCL27_PSY_box.pdf", plot = plot, width = 15, height = 20, units = "cm")

# Boxplot for AKT1S1 in PSY
probe = "AKT1S1"
plot_df = PSY_PHENO_RECALL_OLINK_FINAL_curated

probe_vals = OLINK_DS_BOTH_BATCHES_FINAL_SUAS_ADJ[probe,plot_df$Matching_sample]
all(colnames(probe_vals) == plot_df$Matching_sample)
probe_vals = as.numeric(probe_vals)
plot_df = cbind(plot_df, probe_vals)

plot = ggplot(data = plot_df, aes(x = SUAS_binary, y = probe_vals, fill = SUAS_binary)) +
  geom_boxplot(notchwidth = 0.25) +
  scale_fill_manual(values = c("#C44652","#EAC151")) +
  geom_violin(alpha = 0.25) +
  labs(y = "NPX: AKT1S1", x = "SUAS risk group", fill = "Disease status") +
  ggthemes::theme_clean(base_size = 14)
plot
ggsave(file = "AKT1S1_PSY_box.pdf", plot = plot, width = 15, height = 20, units = "cm")

################### Saving PSY data ###################

PSY_results_limma = list(
  "PSY_no_antidepr" = Top_table,
  "PSY_antidep_adj" = Top_table_antidepr
  
)

names_for_list_PSY = c(
  "PSY_no_antidepr",
  "PSY_antidep_adj"
  )

wb = createWorkbook()

for (i in 1:2) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, names_for_list_PSY[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = names_for_list_PSY[i], PSY_results_limma[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = names_for_list_PSY[i], cols = 1:ncol(PSY_results_limma[[i]]), widths = "auto")
}

saveWorkbook(wb, "PSY_proteomics_SUAS.xlsx", overwrite = TRUE)

# Dataset caracterization

PSY_charact_df = PSY_PHENO_RECALL_OLINK_FINAL_curated

# variables
# SUAS_binary + Antidepr_anxiety + Gender + Age
which(colnames(PSY_charact_df) == "Code") # 1
which(colnames(PSY_charact_df) == "Gender") # 4
which(colnames(PSY_charact_df) == "Age") # 3
which(colnames(PSY_charact_df) == "SUAS_binary") # 93
which(colnames(PSY_charact_df) == "Antidepr_anxiety") # 92
PSY_charact_df$Antidepr_anxiety = ifelse(PSY_charact_df$Antidepr_anxiety == 0, "No treatment", "Treatment")
PSY_charact_df$Antidepr_anxiety = factor(PSY_charact_df$Antidepr_anxiety, levels = c("No treatment", "Treatment"))

PSY_summary = characterize_dataset_generelized_two_subgroups(dataset = PSY_charact_df, 
                                                                   study_char = "PSY", 
                                                                   contrast_col_number = 93,
                                                                   participants_col_number = 1, 
                                                                   model_covariates_col_vector = c(4,3,92), 
                                                                   columns_to_characterise_vector = c(93,4,3,92), 
                                                                   Remove_NA_predictors = FALSE, 
                                                                   drop_P = TRUE)

PSY_summary = PSY_summary$Table
PSY_summary = PSY_summary[-3,]
PSY_summary
PSY_summary = apply(X = PSY_summary, MARGIN = 2, function(x){
  x = multiple_stri_replacer(string = x,
                             pattern_vector = c("Initial dataset includes",
                                                "SUAS_binary",
                                                "Antidepr_anxiety",
                                                "W: ",
                                                "M: "), 
                             replacement_vector = c("Dataset includes",
                                                    "Suicide risk group (SUAS)",
                                                    "Treatment for depression or anxiety",
                                                    "Female: ",
                                                    "Male: "))
  
})

PSY_summary = as.data.frame(PSY_summary)
PSY_summary

wb = createWorkbook()

# Add a new worksheet with the sheet name
addWorksheet(wb, sheetName = "PSY_demographics")

# Write the data frame to the worksheet
writeData(wb, sheet = "PSY_demographics", PSY_summary)

# Adjust column widths to fit the text
setColWidths(wb, sheet = "PSY_demographics", cols = 1:ncol(PSY_summary), widths = "auto")

# Adjust row heights
setRowHeights(wb, sheet = "PSY_demographics", rows = 1:nrow(PSY_summary), heights = "auto")

saveWorkbook(wb, "PSY_demographics.xlsx", overwrite = TRUE)


################### Array Suicide datasets GEO ###################
# Jul 17 2024
# Organism - human
# Study type: ("suicide"[MeSH Terms] OR suicide[All Fields]) AND ("human"[Organism] AND ("Expression profiling by array"[Filter] OR "Expression profiling by genome tiling array"[Filter]))

# https://www.ncbi.nlm.nih.gov/gds/?term=suicide

## GSE146446
geo_test = getGEO("GSE146446")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
geo_test_pheno$`tissue:ch1` # only blood -> no data on suicide
# GSE146446 is only for antidepressant response

## GSE92538
# seems good
geo_test = getGEO("GSE92538")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
geo_test_pheno$`suicide (1=yes):ch1`
table(geo_test_pheno$`suicide (1=yes):ch1`) # 95 33 
# GSE92538 seems good

## GSE208338
geo_test = getGEO("GSE208338")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`suicide:ch1`, geo_test_pheno$`diagnosis:ch1`)
# GSE208338 this cohort is good!

## GSE5388 and GSE5389
geo_test = getGEO("GSE5388")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`Suicide:ch1`, geo_test_pheno$`Disease_status:ch1`)
# GSE5388 seems good

geo_test = getGEO("GSE5389")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`Suicide:ch1`, geo_test_pheno$`Disease_status:ch1`)
# GSE5389 seems good

## GSE66937
geo_test = getGEO("GSE66937")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`group:ch1`, geo_test_pheno$`brain region:ch1`)
# GSE66937 seems good

## GSE125681
geo_test = getGEO("GSE125681")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`suicide:ch1`)
# GSE125681 only suicide and there are no controls...


## GSE24095
geo_test = getGEO("GSE24095")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`Cause of death:ch1`)
table(geo_test_pheno$`Cause of death:ch2`)
# GSE24095 weird design -> better to avoid?
# Not a comparison of suicide VS control -> exclude

## GSE199536
geo_test = getGEO("GSE199536")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`group:ch1`)
# GSE199536 seems good

## GSE20568
geo_test = getGEO("GSE20568")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`cause of death:ch1`) # 1 suicide -> useless
# GSE20568 is useless

## GSE93112
geo_test = getGEO("GSE93112")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$`disease state:ch1`) 
table(geo_test_pheno$`response phenotype, alda scale:ch1`) 
# response to treatment Yes/No
# GSE93112 is useless

## GSE9058
geo_test = getGEO("GSE9058")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
# bad phenotypes and 1 suicide -> useless
# GSE9058 is useless

#### Array Cohorts to include: GSE208338, GSE5388, GSE5389, GSE66937, GSE199536, GSE92538


################### RNA Seq Suicide datasets GEO ###################

# ("suicide"[MeSH Terms] OR suicide[All Fields]) AND ("human"[Organism] AND "Expression profiling by high throughput sequencing"[Filter])

## GSE102556
geo_test = getGEO("GSE102556")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1) # 263 homo sapiens
table(geo_test_pheno$`tissue:ch1`, geo_test_pheno$`Cause of death:ch1`)
# GSE102556 is good

## GSE202537
geo_test = getGEO("GSE202537")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1) # 215 homo sapiens
table(geo_test_pheno$`tissue:ch1`, geo_test_pheno$`manner of death:ch1`)
# GSE202537 is good

## GSE243356
geo_test = getGEO("GSE243356")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1) # 61 homo sapiens
table(geo_test_pheno$`group:ch1`, geo_test_pheno$`tissue:ch1`)
# GSE243356 is good

## GSE247998
geo_test = getGEO("GSE247998")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 100 homo sapiens
table(geo_test_pheno$`suicide attempt:ch1`) 
# GSE247998 is good

## GSE101521
geo_test = getGEO("GSE101521")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 27 homo sapiens
table(geo_test_pheno$`diagnosis:ch1`) 
# GSE101521 is good

## GSE157197
geo_test = getGEO("GSE157197")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 45 homo sapiens
table(geo_test_pheno$`diagnosis:ch1`) 
# GSE157197 is useless

## GSE213982
geo_test = getGEO("GSE213982")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 36 homo sapiens
table(geo_test_pheno$`group:ch1`) 
# GSE213982 is good

## GSE144136
# GSE144136 is single cell -> skipped for now as difficult to compare to bulk

## GSE248260
geo_test = getGEO("GSE248260")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 24 homo sapiens
table(geo_test_pheno$`suicide:ch1`) 
# GSE248260 is good

## GSE12297
geo_test = getGEO("GSE12297")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 20 homo sapiens
# only 3 suicide
# GSE12297 is useless


## GSE193417
geo_test = getGEO("GSE193417")
geo_test = geo_test[[1]]
geo_test_pheno = pData(geo_test)
table(geo_test_pheno$organism_ch1)  # 12 homo sapiens
# only 3 suicide
# GSE193417 is useless

## GSE208438
# GSE208438 is on cell lines -> useless

# Single cell
# GSE144136 (SC) MDD vs control (suicide)
# GSE213982 (SC) MDD vs control (suicide)

#### RNA seq cohorts to include: GSE102556, GSE202537, GSE243356, GSE247998, GSE101521, GSE213982, GSE248260 

################### Selected cohorts ###################

#### Array Cohorts to include: GSE208338, GSE5388, GSE5389, GSE66937, GSE199536, GSE92538 (limma-trend (standard pipeline))
#### RNA seq cohorts to include: GSE102556, GSE202537, GSE243356, GSE247998, GSE101521, GSE213982, GSE144136, GSE248260 (limma-voom (RNAseq specific))

# genome version in records:
# GSE102556 hg19 n>240
# GSE243356 most likely hg19 n=61
# GSE248260 hg19 n=24 # Macedonian/New York State Psychiatric Institute Brain Collection
# GSE247998 hg19 n=100 # Study participants were adults 20–65 years old and were recruited in the New York metropolitan area by advertising, via the Columbia University Medical Center Portal 
# GSE202537 hg38 n=215
# GSE101521 hg19 n=86

# Single cell GSE144136
# Single cell GSE213982

# Criteria: More than 5 suicide cases
# Criteria: Individual signals/counts are available
# Criteria: Study/Methods allow comparison between suicide and control
# Criteria: Human samples only

# Single cell (separate)
# GSE144136 (SC) MDD vs control (suicide)
# GSE213982 (SC) MDD vs control (suicide)

#### Competitor: https://www.sciencedirect.com/science/article/pii/S0924977X21016631?via%3Dihub#sec0002
# competitor cohorts:
# GSE5389
# GSE5388
# GSE66937
# GSE102556
# GSE92538


#### Inspection of datasets from https://www.sciencedirect.com/science/article/pii/S0022395621000480?via%3Dihub (post-analysis)
# GSE21138 - no clear evidence for suicide in samples on GEO, paper has potentially matchable variables in the supplement https://www.sciencedirect.com/science/article/pii/S0006899308019410?via%3Dihub#app1
# GSE92538 - used by us
# GSE53987 - no clear evidence for suicide in samples on GEO, paper has NO variables (only grouped): https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0121744#pone-0121744-t001
# GSE54567 - no clear evidence for suicide in samples on GEO, paper does not show variables (shared https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090980)
# GSE54568 - no clear evidence for suicide in samples on GEO, paper does not show variables (shared https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090980)
# GSE5388 - used by us
# GSE5389 - used by us
# GSE54570 - no clear evidence for suicide in samples on GEO, paper does not show variables (shared https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0090980)
# SMRI 2 and SMRI 7 were not available as www.stanleyresearch.org was not responsive

# Conclusion -> GSE21138 could be theoretically useful due to 13 suicide cases however does not fit inclusion protocol based on GEO search
# Other datasets have no readily availbable infromation on suicide or already used by us

################### Importing DE summary stats ###################
summary_paths_no_covar = list.files(path = "Data_preprocessing_analysis", recursive = TRUE, full.names = TRUE)
summary_paths_no_covar = summary_paths_no_covar[stri_detect_fixed(summary_paths_no_covar, pattern = "no_covar")]
summary_paths_no_covar = summary_paths_no_covar[!stri_detect_fixed(summary_paths_no_covar, pattern = "_cell_types")]

summary_paths_with_covar = list.files(path = "Data_preprocessing_analysis", recursive = TRUE, full.names = TRUE)
summary_paths_with_covar = summary_paths_with_covar[stri_detect_fixed(summary_paths_with_covar, pattern = "Top_table")]
summary_paths_with_covar = summary_paths_with_covar[!stri_detect_fixed(summary_paths_with_covar, pattern = "no_covar")]

summary_paths_cell_types =  list.files(path = "Data_preprocessing_analysis", recursive = TRUE, full.names = TRUE)
summary_paths_cell_types = summary_paths_cell_types[stri_detect_fixed(summary_paths_cell_types, pattern = "_cell_types")]

data_list_no_covar = lapply(summary_paths_no_covar, function(x){
  study_name = unlist(stri_split_fixed(x, pattern = "/"))
  study_name = study_name[length(study_name)]
  study_name = stri_replace_all_fixed(study_name, pattern = "_Top_table_no_covar.csv", replacement = "")
  
  x = smart_fread(x)
  x$Study = study_name
  return(x)
})
reference_colnames = colnames(data_list_no_covar[[1]])
data_list_no_covar = lapply(data_list_no_covar, function(x){
  
  if ("Gene_symbol_non_hgnc" %!in% colnames(x)){
    x$Gene_symbol_non_hgnc = NA
  }
  
  x = x[,reference_colnames]
  return(x)
})


data_list_with_covar = lapply(summary_paths_with_covar, function(x){
  study_name = unlist(stri_split_fixed(x, pattern = "/"))
  study_name = study_name[length(study_name)]
  study_name = stri_replace_all_fixed(study_name, pattern = "_Top_table.csv", replacement = "")
  
  x = smart_fread(x)
  x$Study = study_name
  return(x)
})
reference_colnames = colnames(data_list_with_covar[[1]])
data_list_with_covar = lapply(data_list_with_covar, function(x){
  
  if ("Gene_symbol_non_hgnc" %!in% colnames(x)){
    x$Gene_symbol_non_hgnc = NA
  }
  
  x = x[,reference_colnames]
  return(x)
})

data_list_cell_types = lapply(summary_paths_cell_types, function(x){
  study_name = unlist(stri_split_fixed(x, pattern = "/"))
  study_name = study_name[length(study_name)]
  study_name = stri_replace_all_fixed(study_name, pattern = "_Top_table_no_covar_cell_types.csv", replacement = "")
  
  x = smart_fread(x)
  x$Study = study_name
  return(x)
})
data_list_cell_types = lapply(data_list_cell_types, function(x){
  
  if ("Gene_symbol_non_hgnc" %!in% colnames(x)){
    x$Gene_symbol_non_hgnc = NA
  }
  
  x = x[,reference_colnames]
  return(x)
})

################### Harmonizing gene names ###################

gene_names_no_covar = lapply(data_list_no_covar, function(x){
  
  gene_symbols = x$Gene_symbol
  gene_symbols_non_HGNC = x$Gene_symbol_non_hgnc
  
  gene_symbols_all = c(gene_symbols, gene_symbols_non_HGNC)
  gene_symbols_all = unlist(gene_symbols_all)
  gene_symbols_all = unique(gene_symbols_all)
  gene_symbols_all = unlist(stri_split_fixed(gene_symbols_all, pattern = ";"))
  gene_symbols_all = gene_symbols_all[!is.na(gene_symbols_all)]
  gene_symbols_all = gene_symbols_all[gene_symbols_all != ""]
  gene_symbols_all = gene_symbols_all[gene_symbols_all != " "]
  gene_symbols_all = gene_symbols_all[gene_symbols_all != "NA"]
  
  return(gene_symbols_all)
})
gene_names_no_covar = unlist(gene_names_no_covar)
gene_names_no_covar = unique(gene_names_no_covar)

gene_names_with_covar = lapply(data_list_with_covar, function(x){
  
  gene_symbols = x$Gene_symbol
  gene_symbols_non_HGNC = x$Gene_symbol_non_hgnc
  
  gene_symbols_all = c(gene_symbols, gene_symbols_non_HGNC)
  gene_symbols_all = unlist(gene_symbols_all)
  gene_symbols_all = unique(gene_symbols_all)
  gene_symbols_all = unlist(stri_split_fixed(gene_symbols_all, pattern = ";"))
  gene_symbols_all = gene_symbols_all[!is.na(gene_symbols_all)]
  gene_symbols_all = gene_symbols_all[gene_symbols_all != ""]
  gene_symbols_all = gene_symbols_all[gene_symbols_all != " "]
  gene_symbols_all = gene_symbols_all[gene_symbols_all != "NA"]
  
  return(gene_symbols_all)
})
gene_names_with_covar = unlist(gene_names_with_covar)
gene_names_with_covar = unique(gene_names_with_covar)

gene_names_full = c(gene_names_no_covar, gene_names_with_covar)
gene_names_full = unique(gene_names_full)
harmoniz_genes_check = check_gene_symbol_NIH(PRF_gene_symbols = gene_names_full, 
                                             PRF_ref_NIH_expanded = Homo_Sapiens_Gene_info_NIH_expanded,
                                             PRF_replace_NA_with_old = TRUE)
harmoniz_genes_check_dict = harmoniz_genes_check$Suggested.Symbol
names(harmoniz_genes_check_dict) = harmoniz_genes_check$x


combined_df_no_covar = do.call(rbind, data_list_no_covar)
combined_df_wth_covar = do.call(rbind, data_list_with_covar)
combined_df_cell_types = do.call(rbind, data_list_cell_types)

combined_df_no_covar$Corrected_symbol = mapply(function(x,y){
  combined_symbols = c(x, y)
  combined_symbols = combined_symbols[!is.na(combined_symbols)]
  
  if (length(combined_symbols)<1){
    return(NA)
  }
  
  combined_symbols = stri_split_fixed(combined_symbols, pattern = ";")
  combined_symbols = unlist(combined_symbols)
  combined_symbols = combined_symbols[combined_symbols!= ""]
  combined_symbols = combined_symbols[combined_symbols!= " "]
  combined_symbols = combined_symbols[combined_symbols!= "NA"]
  
  verified_symbols = sapply(combined_symbols, function(x) harmoniz_genes_check_dict[x])
  verified_symbols = unique(verified_symbols)
  verified_symbols = paste0(verified_symbols, collapse = ";")
  return(verified_symbols)
}, combined_df_no_covar$Gene_symbol, combined_df_no_covar$Gene_symbol_non_hgnc)
combined_df_wth_covar$Corrected_symbol = mapply(function(x,y){
  combined_symbols = c(x, y)
  combined_symbols = combined_symbols[!is.na(combined_symbols)]
  
  if (length(combined_symbols)<1){
    return(NA)
  }
  
  combined_symbols = stri_split_fixed(combined_symbols, pattern = ";")
  combined_symbols = unlist(combined_symbols)
  combined_symbols = combined_symbols[combined_symbols!= ""]
  combined_symbols = combined_symbols[combined_symbols!= " "]
  combined_symbols = combined_symbols[combined_symbols!= "NA"]
  
  verified_symbols = sapply(combined_symbols, function(x) harmoniz_genes_check_dict[x])
  verified_symbols = unique(verified_symbols)
  verified_symbols = paste0(verified_symbols, collapse = ";")
  return(verified_symbols)
}, combined_df_wth_covar$Gene_symbol, combined_df_wth_covar$Gene_symbol_non_hgnc)
combined_df_cell_types$Corrected_symbol = mapply(function(x,y){
  combined_symbols = c(x, y)
  combined_symbols = combined_symbols[!is.na(combined_symbols)]
  
  if (length(combined_symbols)<1){
    return(NA)
  }
  
  combined_symbols = stri_split_fixed(combined_symbols, pattern = ";")
  combined_symbols = unlist(combined_symbols)
  combined_symbols = combined_symbols[combined_symbols!= ""]
  combined_symbols = combined_symbols[combined_symbols!= " "]
  combined_symbols = combined_symbols[combined_symbols!= "NA"]
  
  verified_symbols = sapply(combined_symbols, function(x) harmoniz_genes_check_dict[x])
  verified_symbols = unique(verified_symbols)
  verified_symbols = paste0(verified_symbols, collapse = ";")
  return(verified_symbols)
}, combined_df_cell_types$Gene_symbol, combined_df_cell_types$Gene_symbol_non_hgnc)

################### Correcting cell type df ###################

combined_df_cell_types$Cell_type = stri_replace_all_fixed(combined_df_cell_types$Tissue, pattern = "Dorsolateral prefrontal cortex (BA9): ", replacement = "")
combined_df_cell_types$Study_Cell = paste0(combined_df_cell_types$Study, ": ", combined_df_cell_types$Cell_type)

table(combined_df_cell_types$Cell_type)
"
  Ast   End   ExN   InN   Mic   Mix   Oli   OPC 
 8510  6145 34875 26001  2093  1305  9025  6806 
"
# Ast -astrocyte
# End - endothelial
# ExN - excitatory neuron
# InN - inhibitory neuron
# Mic - microglia
# Mix - mixed expression profile
# Oli - oligodendrocytes
# OPC - oligodendrocyte precursor cells

combined_df_cell_types_meta = combined_df_cell_types %>%
  group_by(.,Study_Cell, Corrected_symbol)  %>%
  summarise(.,avgLog2FC = mean(logFC), 
            maxSE = max(SE)) 
combined_df_cell_types_meta = ungroup(combined_df_cell_types_meta)
combined_df_cell_types_meta = as.data.frame(combined_df_cell_types_meta)

################### Cohort stats on individual datasets ###################
length(unique(combined_df_no_covar$Study))
length(unique(combined_df_wth_covar$Study))

studies = unique(combined_df_no_covar$Study)

cohort_stats_vector_no_covariates = vector()
for (i in 1:length(studies)){
  
  separator_large = "\n\n"
  if (i == 1){
    
    string_1 = paste0("Total associations without covariates: ", nrow(combined_df_no_covar))
    string_2 = paste0("Total significant associations (nominal) without covariates: ", nrow(combined_df_no_covar[combined_df_no_covar$P.Value<=0.05,]))
    string_3 = paste0("Average absolute effect size (log2FC) without covariates: ", mean(abs(combined_df_no_covar$logFC)))
    string_4 = paste0("Average absolute effect size (log2FC) (nominal significance) without covariates: ", mean(abs(combined_df_no_covar[combined_df_no_covar$P.Value<=0.05,"logFC"])))
    cohort_stats_vector_no_covariates = c(string_1, string_2, string_3,  string_4, separator_large)
  }
  
  TMP_current_cohort = studies[i]
  TMP_current_dataset = combined_df_no_covar[combined_df_no_covar$Study == TMP_current_cohort,]
  TMP_study_type = unique(TMP_current_dataset$Tissue_type)
  TMP_tech_type = unique(TMP_current_dataset$Technology)
  
  main_header = paste(c(TMP_current_cohort, TMP_study_type, TMP_tech_type), collapse = "\t")
  
  aggregated_stats = TMP_current_dataset %>%
    group_by(.,Tissue) %>%
    summarise(
      total_associations = length(ID),
      unique_gene_symbols_p_lt_0_05 = n_distinct(Gene_symbol[P.Value <= 0.05]),
      unique_IDs_p_lt_0_05 = n_distinct(ID[P.Value <= 0.05]),
      unique_gene_symbols_adj_p_lt_0_05 = n_distinct(Gene_symbol[adj.P.Val <= 0.05]),
      unique_IDs_adj_p_lt_0_05 = n_distinct(ID[adj.P.Val <= 0.05]),
      mean_abs_log2FC = mean(abs(logFC)),
      mean_abs_log2FC_nomin_signif = mean(abs(logFC[P.Value <= 0.05]))
    )
  aggregated_stats = as.data.frame(aggregated_stats)
  header = paste(colnames(aggregated_stats), collapse = "\t")
  # char_aggregated_stats = apply(aggregated_stats, 1, function(x) paste(x, collapse = "\t"))
  
  char_aggregated_stats = knitr::kable(aggregated_stats)
  
  cohort_stats_vector_no_covariates = c(cohort_stats_vector_no_covariates, main_header, char_aggregated_stats, separator_large)
  if (i == length(studies)){
    cohort_stats_vector_no_covariates = c(cohort_stats_vector_no_covariates, separator_large)
  }
  
  
}
writeLines(text = cohort_stats_vector_no_covariates, con = "individual_DE_stats_no_covar.txt")

cohort_stats_vector_with_covariates = vector()
for (i in 1:length(studies)){
  
  separator_large = "\n\n"
  if (i == 1){
    
    string_1 = paste0("Total associations with covariates: ", nrow(combined_df_wth_covar))
    string_2 = paste0("Total significant associations (nominal) with covariates: ", nrow(combined_df_wth_covar[combined_df_wth_covar$P.Value<0.05,]))
    string_3 = paste0("Average absolute effect size (log2FC) with covariates: ", mean(abs(combined_df_wth_covar$logFC)))
    string_4 = paste0("Average absolute effect size (log2FC) (nominal significance) with covariates: ", mean(abs(combined_df_wth_covar[combined_df_wth_covar$P.Value<0.05,"logFC"])))
    cohort_stats_vector_with_covariates = c(string_1, string_2, string_3,  string_4, separator_large)
  }
  
  TMP_current_cohort = studies[i]
  TMP_current_dataset = combined_df_wth_covar[combined_df_wth_covar$Study == TMP_current_cohort,]
  TMP_study_type = unique(TMP_current_dataset$Tissue_type)
  TMP_tech_type = unique(TMP_current_dataset$Technology)
  
  main_header = paste(c(TMP_current_cohort, TMP_study_type, TMP_tech_type), collapse = "\t")
  
  aggregated_stats = TMP_current_dataset %>%
    group_by(.,Tissue) %>%
    summarise(
      total_associations = length(ID),
      unique_gene_symbols_p_lt_0_05 = n_distinct(Gene_symbol[P.Value < 0.05]),
      unique_IDs_p_lt_0_05 = n_distinct(ID[P.Value < 0.05]),
      unique_gene_symbols_adj_p_lt_0_05 = n_distinct(Gene_symbol[adj.P.Val < 0.05]),
      unique_IDs_adj_p_lt_0_05 = n_distinct(ID[adj.P.Val < 0.05]),
      mean_abs_log2FC = mean(abs(logFC)),
      mean_abs_log2FC_nomin_signif = mean(abs(logFC[P.Value < 0.05]))
    )
  aggregated_stats = as.data.frame(aggregated_stats)
  header = paste(colnames(aggregated_stats), collapse = "\t")
  # char_aggregated_stats = apply(aggregated_stats, 1, function(x) paste(x, collapse = "\t"))
  
  char_aggregated_stats = knitr::kable(aggregated_stats)
  
  cohort_stats_vector_with_covariates = c(cohort_stats_vector_with_covariates, main_header, char_aggregated_stats, separator_large)
  
  if (i == length(studies)){
    cohort_stats_vector_with_covariates = c(cohort_stats_vector_with_covariates, separator_large)
  }
}
writeLines(text = cohort_stats_vector_with_covariates, con = "individual_DE_stats_with_covar.txt")


################### Cohort inspection for specificity ###################
table(combined_df_no_covar$Tissue, combined_df_no_covar$Study)

# GSE101521 -> dorsal lateral prefrontal cortex (Brodmann Area 9) "with an attempt to enrich for gray matter" "samples inevitably also contained a small amount of white matter"
# GSE102556 -> many regions, no specifics on gray/white
# GSE144136 -> single cell, "Frozen gray matter samples were dissected from Brodmann area 9 (dlPFC)"
# GSE199536 -> habenula. "4.4 Habenula In 1872 he described a small mass of grey matter", paper does not specify
# GSE202537 -> several striatal regions, no specifics on gray/white matter content
# GSE208338 -> no original paper, GEO: "Total RNA was isolated from ∼100 mg frozen gray matter using 1.0 ml TRIzol reagent"
# GSE213982 -> "Frozen histological grade samples of gray and white matter were dissected from the dlPFC (Brodmann Area 9)" -> mixed matters
# GSE243356 -> "The temporal pole, Brodmann Area 20, gray and white matter were dissected for the assays" -> mix of matters
# Since the dissected brain gyri fold in three dimensions, it is impossible to exclude all white matter while still taking the full thickness (~ 3 mm) of the cortical"
# GSE247998 -> blood
# GSE248260 -> ventral white matter (BA 47) -> gray matter was not provided
# GSE5388 -> DLPFC (BA9) no specifics on gray/white matter content
# GSE5389 -> OFC (BA11) no specifics on gray/white matter content
# GSE66937 -> no original paper, no specifics on gray/white matter content, has technical replicates!
# GSE92538_U133A -> predominantly gray matter, has samples on 2 arrays
# GSE92538_U133_PLUS2 -> predominantly gray matter, has samples on 2 arrays

################### GSE92538 inspection ###################

GSE92538_U133A = read.csv("Data_preprocessing_analysis/GSE92538_U133A_results/GSE92538_pheno_curated.csv")
rownames(GSE92538_U133A) = GSE92538_U133A$X
GSE92538_U133A$X = NULL
nrow(GSE92538_U133A) # 97

GSE92538_U133_PLUS2 = read.csv("Data_preprocessing_analysis/GSE92538_U133_PLUS2_results/GSE92538_pheno_curated.csv")
rownames(GSE92538_U133_PLUS2) = GSE92538_U133_PLUS2$X
GSE92538_U133_PLUS2$X = NULL
nrow(GSE92538_U133_PLUS2) # 75
any(duplicated(GSE92538_U133A$SUBJECT_ID)) # FALSE
any(duplicated(GSE92538_U133_PLUS2$SUBJECT_ID)) # FALSE

length(intersect(GSE92538_U133A$SUBJECT_ID, GSE92538_U133_PLUS2$SUBJECT_ID))
# 0 samples are duplicated

################### Meta (no covar) all brain data ###################
combined_df_no_covar_meta_all = combined_df_no_covar
combined_df_no_covar_meta_all = combined_df_no_covar_meta_all[combined_df_no_covar_meta_all$Tissue_type == "Brain",]
blood_df_no_covar = combined_df_no_covar[combined_df_no_covar$Tissue_type == "Blood",]

table(combined_df_no_covar_meta_all$Tissue, combined_df_no_covar_meta_all$Study)
# GSE102556 many tissues
# GSE202537 many tissues
# GSE66937 many tissues

combined_df_no_covar_meta_all_1 = combined_df_no_covar_meta_all[combined_df_no_covar_meta_all$Study %!in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_no_covar_meta_all_2 = combined_df_no_covar_meta_all[combined_df_no_covar_meta_all$Study %in% c("GSE102556", "GSE202537", "GSE66937"), ]
# Caudate https://pubmed.ncbi.nlm.nih.gov/39164232/
# Nac https://pubmed.ncbi.nlm.nih.gov/38894648/ https://pubmed.ncbi.nlm.nih.gov/38965529/ https://pubmed.ncbi.nlm.nih.gov/39167467/
# GSE102556 -> Dorsolateral prefrontal cortex (dlPFC; BA8/9)
# GSE202537 -> Nac 22652  (more probes analyzed)
# GSE66937 -> prefrontal cortex
combined_df_no_covar_meta_all_2 = combined_df_no_covar_meta_all_2[combined_df_no_covar_meta_all_2$Tissue %in% c("Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                                                "Nac",
                                                                                                                "prefrontal cortex"),]
table(combined_df_no_covar_meta_all_2$Tissue)
"
Dorsolateral prefrontal cortex (dlPFC; BA8/9)                                           Nac                             prefrontal cortex 
                                        21365                                         22652                                         70523 
"
combined_df_no_covar_meta_full = rbind(combined_df_no_covar_meta_all_1, combined_df_no_covar_meta_all_2)
table(combined_df_no_covar_meta_full$Tissue, combined_df_no_covar_meta_full$Study)
combined_df_no_covar_meta_full_list = lapply(unique(combined_df_no_covar_meta_full$Study), function(x){
  df = combined_df_no_covar_meta_full[combined_df_no_covar_meta_full$Study == x,]
  return(df)
})

combined_df_no_covar_meta_full_list_reduced = lapply(combined_df_no_covar_meta_full_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_no_covar_meta_full_list_reduced = do.call(rbind, combined_df_no_covar_meta_full_list_reduced)
unique_genes_meta = unique(combined_df_no_covar_meta_full_list_reduced$Corrected_symbol)

meta_no_covar_all_brain = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_no_covar_meta_full_list_reduced[combined_df_no_covar_meta_full_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_no_covar[blood_df_no_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
meta_no_covar_all_brain = do.call(rbind, meta_no_covar_all_brain)
rownames(meta_no_covar_all_brain) = NULL
rm(list = ls(pattern = "tmp_"))
meta_no_covar_all_brain = dplyr::arrange(meta_no_covar_all_brain, meta_pval)
meta_no_covar_all_brain$meta_FDR = p.adjust(meta_no_covar_all_brain$meta_pval, method = "fdr")
meta_no_covar_all_brain = meta_no_covar_all_brain[,c("gene",
                                                     "detected_cohorts",
                                                     "cohorts_up",
                                                     "cohorts_down",
                                                     "commentary",
                                                     "meta_LFc",
                                                     "meta_se",
                                                     "meta_pval",
                                                     "meta_FDR",
                                                     "tau2",
                                                     "I2",
                                                     "H2" ,
                                                     "Q",
                                                     "Q.p",
                                                     "mean_blood_lfc",
                                                     "blood_dir",
                                                     "blood_signif",
                                                     "matching_with_blood")]

meta_no_covar_all_brain_signif = meta_no_covar_all_brain[meta_no_covar_all_brain$meta_pval < 0.05,]
meta_no_covar_all_brain_signif = meta_no_covar_all_brain_signif[!is.na(meta_no_covar_all_brain_signif$meta_pval),]
nrow(meta_no_covar_all_brain_signif)

# P2RY12
# https://pubmed.ncbi.nlm.nih.gov/37396924/
# https://www.ingentaconnect.com/content/ben/cnr/2022/00000019/00000003/art00004
# https://pubmed.ncbi.nlm.nih.gov/31739114/
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8908792/
# https://www.mdpi.com/1422-0067/25/11/5750
# https://psychiatryonline.org/doi/10.1176/appi.ajp.2021.22010026

plot_forest_simple_gene(gene_name = "P2RY12", meta_df = combined_df_cell_types_meta, study_column_char = "Study_Cell")
plot_forest_meta_gene(gene_name = "P2RY12", meta_df = combined_df_no_covar_meta_full_list_reduced)

# making plots
dir.create("forest_plots")
folder = "forest_plots/meta_no_covar_all_brain"
dir.create(folder)

for (i in 1:nrow(meta_no_covar_all_brain_signif)){
  
  
  gene = meta_no_covar_all_brain_signif$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_no_covar_meta_full_list_reduced))
  dev.off()
  
}

################### Meta (with covar) all brain data ###################
combined_df_with_covar_meta_all = combined_df_wth_covar
combined_df_with_covar_meta_all = combined_df_with_covar_meta_all[combined_df_with_covar_meta_all$Tissue_type == "Brain",]
blood_df_with_covar = combined_df_wth_covar[combined_df_wth_covar$Tissue_type == "Blood",]

table(combined_df_with_covar_meta_all$Tissue, combined_df_with_covar_meta_all$Study)

# GSE102556 many tissues
# GSE202537 many tissues
# GSE66937 many tissues

combined_df_with_covar_meta_all_1 = combined_df_with_covar_meta_all[combined_df_with_covar_meta_all$Study %!in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_with_covar_meta_all_2 = combined_df_with_covar_meta_all[combined_df_with_covar_meta_all$Study %in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_with_covar_meta_all_2 = combined_df_with_covar_meta_all_2[combined_df_with_covar_meta_all_2$Tissue %in% c("Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                                                "Nac",
                                                                                                                "prefrontal cortex"),]
table(combined_df_with_covar_meta_all_2$Tissue)
"
Dorsolateral prefrontal cortex (dlPFC; BA8/9)                                           Nac                             prefrontal cortex 
                                        23331                                         23919                                         70523 
"

combined_df_with_covar_meta_full = rbind(combined_df_with_covar_meta_all_1, combined_df_with_covar_meta_all_2)
table(combined_df_with_covar_meta_full$Tissue, combined_df_with_covar_meta_full$Study)
combined_df_with_covar_meta_full_list = lapply(unique(combined_df_with_covar_meta_full$Study), function(x){
  df = combined_df_with_covar_meta_full[combined_df_with_covar_meta_full$Study == x,]
  return(df)
})

combined_df_with_covar_meta_full_list_reduced = lapply(combined_df_with_covar_meta_full_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_with_covar_meta_full_list_reduced = do.call(rbind, combined_df_with_covar_meta_full_list_reduced)
unique_genes_meta = unique(combined_df_with_covar_meta_full_list_reduced$Corrected_symbol)

meta_with_covar_all_brain = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_with_covar_meta_full_list_reduced[combined_df_with_covar_meta_full_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_with_covar[blood_df_with_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
meta_with_covar_all_brain = do.call(rbind, meta_with_covar_all_brain)
rownames(meta_with_covar_all_brain) = NULL
rm(list = ls(pattern = "tmp_"))
meta_with_covar_all_brain = dplyr::arrange(meta_with_covar_all_brain, meta_pval)
meta_with_covar_all_brain$meta_FDR = p.adjust(meta_with_covar_all_brain$meta_pval, method = "fdr")
meta_with_covar_all_brain = meta_with_covar_all_brain[,c("gene",
                                                     "detected_cohorts",
                                                     "cohorts_up",
                                                     "cohorts_down",
                                                     "commentary",
                                                     "meta_LFc",
                                                     "meta_se",
                                                     "meta_pval",
                                                     "meta_FDR",
                                                     "tau2",
                                                     "I2",
                                                     "H2" ,
                                                     "Q",
                                                     "Q.p",
                                                     "mean_blood_lfc",
                                                     "blood_dir",
                                                     "blood_signif",
                                                     "matching_with_blood")]

meta_with_covar_all_brain_signif = meta_with_covar_all_brain[meta_with_covar_all_brain$meta_pval < 0.05,]

meta_with_covar_all_brain_signif = meta_with_covar_all_brain_signif[!is.na(meta_with_covar_all_brain_signif$meta_pval),]
nrow(meta_with_covar_all_brain_signif)

plot_forest_meta_gene(gene_name = "RB1", meta_df = combined_df_with_covar_meta_full_list_reduced)
plot_forest_meta_gene(gene_name = "P2RY12", meta_df = combined_df_with_covar_meta_full_list_reduced)

# Overlap full analysis
overlapping_genes_meta_full = intersect(meta_no_covar_all_brain_signif$gene, meta_with_covar_all_brain_signif$gene)
overlapping_genes_meta_full = unique(overlapping_genes_meta_full) 
length(overlapping_genes_meta_full) # 57 genes overlap

# making plots
dir.create("forest_plots")
folder = "forest_plots/meta_with_covar_all_brain"
dir.create(folder)

for (i in 1:nrow(meta_with_covar_all_brain_signif)){
  
  
  gene = meta_with_covar_all_brain_signif$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_with_covar_meta_full_list_reduced))
  dev.off()
  
}

################### Meta cortical regions (no covar) ###################
combined_df_no_covar_meta_cortical = combined_df_no_covar
combined_df_no_covar_meta_cortical = combined_df_no_covar_meta_cortical[combined_df_no_covar_meta_cortical$Tissue_type == "Brain",]
blood_df_no_covar = combined_df_no_covar[combined_df_no_covar$Tissue_type == "Blood",]

table(combined_df_no_covar_meta_cortical$Tissue, combined_df_no_covar_meta_cortical$Study)
# Tissues to include
# Cingulate gyrus 25 (Cg25) is also cortical but Dorsolateral prefrontal cortex (dlPFC; BA8/9)  was used instead
combined_df_no_covar_meta_cortical = combined_df_no_covar_meta_cortical[combined_df_no_covar_meta_cortical$Tissue %in% 
                                                                          c("DLPFC",
                                                                            "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                            "Dorsolateral Prefrontal Cortex",
                                                                            "Dorsolateral prefrontal cortex (BA9)",
                                                                            "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                            "Orbitofrontal cortex",
                                                                            "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                            "prefrontal cortex",
                                                                            "temporal cortex (BA20 and BA36)"),]
table(combined_df_no_covar_meta_cortical$Tissue, combined_df_no_covar_meta_cortical$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE243356 GSE5388 GSE5389 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0         0   22283       0        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     18677         0         0         0         0         0       0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0         0       0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     19271         0     19191         0       0       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     21365         0         0         0         0       0       0        0                   0              0
  Orbitofrontal cortex                                       0         0         0         0         0         0       0   22283        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0         0       0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0         0       0       0    70523                   0              0
  temporal cortex (BA20 and BA36)                            0         0         0         0         0     20955       0       0        0                   0              0

"
combined_df_no_covar_meta_cortical_list = lapply(unique(combined_df_no_covar_meta_cortical$Study), function(x){
  df = combined_df_no_covar_meta_cortical[combined_df_no_covar_meta_cortical$Study == x,]
  return(df)
})

combined_df_no_covar_meta_cortical_list_reduced = lapply(combined_df_no_covar_meta_cortical_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_no_covar_meta_cortical_list_reduced = do.call(rbind, combined_df_no_covar_meta_cortical_list_reduced)
unique_genes_meta = unique(combined_df_no_covar_meta_cortical_list_reduced$Corrected_symbol)

meta_no_covar_cortical = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_no_covar_meta_cortical_list_reduced[combined_df_no_covar_meta_cortical_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_no_covar[blood_df_no_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
meta_no_covar_cortical = do.call(rbind, meta_no_covar_cortical)
rownames(meta_no_covar_cortical) = NULL
rm(list = ls(pattern = "tmp_"))
meta_no_covar_cortical = dplyr::arrange(meta_no_covar_cortical, meta_pval)
meta_no_covar_cortical$meta_FDR = p.adjust(meta_no_covar_cortical$meta_pval, method = "fdr")
meta_no_covar_cortical = meta_no_covar_cortical[,c("gene",
                                                   "detected_cohorts",
                                                   "cohorts_up",
                                                   "cohorts_down",
                                                   "commentary",
                                                   "meta_LFc",
                                                   "meta_se",
                                                   "meta_pval",
                                                   "meta_FDR",
                                                   "tau2",
                                                   "I2",
                                                   "H2" ,
                                                   "Q",
                                                   "Q.p",
                                                   "mean_blood_lfc",
                                                   "blood_dir",
                                                   "blood_signif",
                                                   "matching_with_blood")]

meta_no_covar_cortical_signif = meta_no_covar_cortical[meta_no_covar_cortical$meta_pval < 0.05,]
meta_no_covar_cortical_signif = meta_no_covar_cortical_signif[!is.na(meta_no_covar_cortical_signif$meta_pval),]
meta_no_covar_cortical[meta_no_covar_cortical$gene == "P2RY12",]
nrow(meta_no_covar_cortical_signif)

# making plots
dir.create("forest_plots")
folder = "forest_plots/meta_no_covar_cortical"
dir.create(folder)

for (i in 1:nrow(meta_no_covar_cortical_signif)){
  
  
  gene = meta_no_covar_cortical_signif$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_no_covar_meta_cortical_list_reduced))
  dev.off()
  
}


################### Meta cortical regions (with covar) ###################
combined_df_with_covar_meta_cortical = combined_df_wth_covar
combined_df_with_covar_meta_cortical = combined_df_with_covar_meta_cortical[combined_df_with_covar_meta_cortical$Tissue_type == "Brain",]
blood_df_with_covar = combined_df_wth_covar[combined_df_wth_covar$Tissue_type == "Blood",]

# Tissues to include
combined_df_with_covar_meta_cortical = combined_df_with_covar_meta_cortical[combined_df_with_covar_meta_cortical$Tissue %in% 
                                                                          c("DLPFC",
                                                                            "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                            "Dorsolateral Prefrontal Cortex",
                                                                            "Dorsolateral prefrontal cortex (BA9)",
                                                                            "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                            "Orbitofrontal cortex",
                                                                            "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                            "prefrontal cortex",
                                                                            "temporal cortex (BA20 and BA36)"),]
table(combined_df_with_covar_meta_cortical$Tissue, combined_df_with_covar_meta_cortical$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE243356 GSE5388 GSE5389 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0         0   22283       0        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     21197         0         0         0         0         0       0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0         0       0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     22389         0     24182         0       0       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     23331         0         0         0         0       0       0        0                   0              0
  Orbitofrontal cortex                                       0         0         0         0         0         0       0   22283        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0         0       0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0         0       0       0    70523                   0              0
  temporal cortex (BA20 and BA36)                            0         0         0         0         0     20955       0       0        0                   0              0

"

combined_df_with_covar_meta_cortical_list = lapply(unique(combined_df_with_covar_meta_cortical$Study), function(x){
  df = combined_df_with_covar_meta_cortical[combined_df_with_covar_meta_cortical$Study == x,]
  return(df)
})

combined_df_with_covar_meta_cortical_list_reduced = lapply(combined_df_with_covar_meta_cortical_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_with_covar_meta_cortical_list_reduced = do.call(rbind, combined_df_with_covar_meta_cortical_list_reduced)
unique_genes_meta = unique(combined_df_with_covar_meta_cortical_list_reduced$Corrected_symbol)

meta_with_covar_cortical = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_with_covar_meta_cortical_list_reduced[combined_df_with_covar_meta_cortical_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_with_covar[blood_df_with_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
meta_with_covar_cortical = do.call(rbind, meta_with_covar_cortical)
rownames(meta_with_covar_cortical) = NULL
rm(list = ls(pattern = "tmp_"))
meta_with_covar_cortical = dplyr::arrange(meta_with_covar_cortical, meta_pval)
meta_with_covar_cortical$meta_FDR = p.adjust(meta_with_covar_cortical$meta_pval, method = "fdr")
meta_with_covar_cortical = meta_with_covar_cortical[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]
meta_with_covar_cortical_signif = meta_with_covar_cortical[meta_with_covar_cortical$meta_pval < 0.05,]
meta_with_covar_cortical_signif = meta_with_covar_cortical_signif[!is.na(meta_with_covar_cortical_signif$meta_pval),]
nrow(meta_with_covar_cortical_signif) # 205
meta_with_covar_cortical[meta_with_covar_cortical$gene == "P2RY12",]

# making plots
dir.create("forest_plots")
folder = "forest_plots/meta_with_covar_cortical"
dir.create(folder)

for (i in 1:nrow(meta_with_covar_cortical_signif)){
  
  
  gene = meta_with_covar_cortical_signif$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_with_covar_meta_cortical_list_reduced))
  dev.off()
  
}



################### Meta prefrontal regions (no covar) ###################
combined_df_no_covar_meta_prefrontal = combined_df_no_covar
combined_df_no_covar_meta_prefrontal = combined_df_no_covar_meta_prefrontal[combined_df_no_covar_meta_prefrontal$Tissue_type == "Brain",]
blood_df_no_covar = combined_df_no_covar[combined_df_no_covar$Tissue_type == "Blood",]

table(combined_df_no_covar_meta_prefrontal$Tissue, combined_df_no_covar_meta_prefrontal$Study)
# Tissues to include
combined_df_no_covar_meta_prefrontal = combined_df_no_covar_meta_prefrontal[combined_df_no_covar_meta_prefrontal$Tissue %in% 
                                                                          c("DLPFC",
                                                                            "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                            "Dorsolateral Prefrontal Cortex",
                                                                            "Dorsolateral prefrontal cortex (BA9)",
                                                                            "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                            "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                            "prefrontal cortex"),]
table(combined_df_no_covar_meta_prefrontal$Tissue, combined_df_no_covar_meta_prefrontal$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE5388 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0   22283        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     18677         0         0         0         0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     19271         0     19191       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     21365         0         0         0       0        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0       0    70523                   0              0
  
"
combined_df_no_covar_meta_prefrontal_list = lapply(unique(combined_df_no_covar_meta_prefrontal$Study), function(x){
  df = combined_df_no_covar_meta_prefrontal[combined_df_no_covar_meta_prefrontal$Study == x,]
  return(df)
})

combined_df_no_covar_meta_prefrontal_list_reduced = lapply(combined_df_no_covar_meta_prefrontal_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_no_covar_meta_prefrontal_list_reduced = do.call(rbind, combined_df_no_covar_meta_prefrontal_list_reduced)
unique_genes_meta = unique(combined_df_no_covar_meta_prefrontal_list_reduced$Corrected_symbol)

meta_no_covar_prefrontal = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_no_covar_meta_prefrontal_list_reduced[combined_df_no_covar_meta_prefrontal_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_no_covar[blood_df_no_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
meta_no_covar_prefrontal = do.call(rbind, meta_no_covar_prefrontal)
rownames(meta_no_covar_prefrontal) = NULL
rm(list = ls(pattern = "tmp_"))
meta_no_covar_prefrontal = dplyr::arrange(meta_no_covar_prefrontal, meta_pval)
meta_no_covar_prefrontal$meta_FDR = p.adjust(meta_no_covar_prefrontal$meta_pval, method = "fdr")
meta_no_covar_prefrontal = meta_no_covar_prefrontal[,c("gene",
                                                       "detected_cohorts",
                                                       "cohorts_up",
                                                       "cohorts_down",
                                                       "commentary",
                                                       "meta_LFc",
                                                       "meta_se",
                                                       "meta_pval",
                                                       "meta_FDR",
                                                       "tau2",
                                                       "I2",
                                                       "H2" ,
                                                       "Q",
                                                       "Q.p",
                                                       "mean_blood_lfc",
                                                       "blood_dir",
                                                       "blood_signif",
                                                       "matching_with_blood")]

meta_no_covar_prefrontal_signif = meta_no_covar_prefrontal[meta_no_covar_prefrontal$meta_pval < 0.05,]
meta_no_covar_prefrontal_signif = meta_no_covar_prefrontal_signif[!is.na(meta_no_covar_prefrontal_signif$meta_pval),]
nrow(meta_no_covar_prefrontal_signif)
meta_no_covar_prefrontal[meta_no_covar_prefrontal$gene == "P2RY12",]

# making plots
dir.create("forest_plots")
folder = "forest_plots/meta_no_covar_prefrontal"
dir.create(folder)

for (i in 1:nrow(meta_no_covar_prefrontal_signif)){
  
  
  gene = meta_no_covar_prefrontal_signif$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_no_covar_meta_prefrontal_list_reduced))
  dev.off()
  
}

################### Meta prefrontal regions (with covar) ###################
combined_df_with_covar_meta_prefrontal = combined_df_wth_covar
combined_df_with_covar_meta_prefrontal = combined_df_with_covar_meta_prefrontal[combined_df_with_covar_meta_prefrontal$Tissue_type == "Brain",]
blood_df_with_covar = combined_df_wth_covar[combined_df_wth_covar$Tissue_type == "Blood",]

# Tissues to include
combined_df_with_covar_meta_prefrontal = combined_df_with_covar_meta_prefrontal[combined_df_with_covar_meta_prefrontal$Tissue %in% 
                                                                              c("DLPFC",
                                                                                "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                                "Dorsolateral Prefrontal Cortex",
                                                                                "Dorsolateral prefrontal cortex (BA9)",
                                                                                "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                                "prefrontal cortex"),]
table(combined_df_with_covar_meta_prefrontal$Tissue, combined_df_with_covar_meta_prefrontal$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE5388 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0   22283        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     21197         0         0         0         0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     22389         0     24182       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     23331         0         0         0       0        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0       0    70523                   0              0
  
"

combined_df_with_covar_meta_prefrontal_list = lapply(unique(combined_df_with_covar_meta_prefrontal$Study), function(x){
  df = combined_df_with_covar_meta_prefrontal[combined_df_with_covar_meta_prefrontal$Study == x,]
  return(df)
})

combined_df_with_covar_meta_prefrontal_list_reduced = lapply(combined_df_with_covar_meta_prefrontal_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_with_covar_meta_prefrontal_list_reduced = do.call(rbind, combined_df_with_covar_meta_prefrontal_list_reduced)
unique_genes_meta = unique(combined_df_with_covar_meta_prefrontal_list_reduced$Corrected_symbol)

meta_with_covar_prefrontal = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_with_covar_meta_prefrontal_list_reduced[combined_df_with_covar_meta_prefrontal_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_with_covar[blood_df_with_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
meta_with_covar_prefrontal = do.call(rbind, meta_with_covar_prefrontal)
rownames(meta_with_covar_prefrontal) = NULL
rm(list = ls(pattern = "tmp_"))
meta_with_covar_prefrontal = dplyr::arrange(meta_with_covar_prefrontal, meta_pval)
meta_with_covar_prefrontal$meta_FDR = p.adjust(meta_with_covar_prefrontal$meta_pval, method = "fdr")
meta_with_covar_prefrontal = meta_with_covar_prefrontal[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]


meta_with_covar_prefrontal_signif = meta_with_covar_prefrontal[meta_with_covar_prefrontal$meta_pval < 0.05,]
meta_with_covar_prefrontal_signif = meta_with_covar_prefrontal_signif[!is.na(meta_with_covar_prefrontal_signif$meta_pval),]
nrow(meta_with_covar_prefrontal_signif)
meta_with_covar_prefrontal[meta_with_covar_prefrontal$gene == "P2RY12",]

# plots 
dir.create("forest_plots")
folder = "forest_plots/meta_with_covar_prefrontal"
dir.create(folder)

for (i in 1:nrow(meta_with_covar_prefrontal_signif)){
  
  
  gene = meta_with_covar_prefrontal_signif$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_with_covar_meta_prefrontal_list_reduced))
  dev.off()
  
}
################### Volcano plot for every meta-analysis ###################

plot_meta_volcano = function(meta_df, figure_path, title){
  
  meta_df = meta_df[!is.na(meta_df$meta_pval),]
  meta_df$adj.P.Val = p.adjust(meta_df$meta_pval, method = "fdr")
  
  Pval_treshold = meta_df[meta_df$adj.P.Val < 0.05,]
  Pval_treshold = max(Pval_treshold$meta_pval)
  Pval_treshold = -log10(Pval_treshold)
  if (is.na(Pval_treshold)){
    Pval_treshold = -log10(0.001)
  }
  PREFIX_logFC_threshold = 0.2
  
  meta_df$is.highlight = sapply(meta_df$meta_LFc, function(x){
    if (x > PREFIX_logFC_threshold){
      x = "Up"
    } else if (x < - PREFIX_logFC_threshold){
      x = "Down"
    } else {
      x = "None"
    }
  })
  meta_df = mutate(meta_df, is_annotate = ifelse(-log10(meta_pval) >= -log10(0.05) & is.highlight != "None", "yes", "no"))
  
  # Make the plot
  plot = ggplot(meta_df, aes(x=meta_LFc, y=-log10(meta_pval))) +
    
    # Show all points
    geom_point(aes(color= factor(is.highlight, levels = c("None", "Down", "Up"))), alpha=0.6, size=4) +
    scale_color_manual(values = c("grey", "skyblue", "red")) 
  
  # Add standard scale for plot
  plot = plot + scale_x_continuous(limits = c(-1.5,1.5))
  plot = plot + scale_y_continuous(limits = c(0,5))
  
  # Add title
  plot = plot + ggtitle(label = title)
  
  # Add pval line
  if (!is.null(Pval_treshold)){
    plot = plot + geom_hline(yintercept=Pval_treshold, linetype="dashed", 
                             color = "red", size=0.5)
  }
  
  # Add pval line 0.05
  plot = plot + geom_hline(yintercept= -log10(0.05), linetype="dashed", 
                           color = "blue", size=0.5)
  # Add logFC lines
  if (min(meta_df$meta_LFc) < -PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= -PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  if (max(meta_df$meta_LFc) > PREFIX_logFC_threshold){
    plot = plot + geom_vline(xintercept= PREFIX_logFC_threshold, linetype="dashed", 
                             color = "grey", size=0.5)
  }
  
  # Add label using ggrepel to avoid overlapping
  if (any(meta_df$is_annotate == "yes")){
    plot = plot + 
      geom_label_repel(data=subset(meta_df, is_annotate=="yes"), aes(label=gene), size=4, force = 10, 
                       max.overlaps = 50)
  } 
  
  plot = plot +
    labs(x = "meta-estimated Log2FC", y = "-log10 p-value") +
    # Custom the theme:
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"), # Modifying horizontal lines in the plot
      panel.background = element_blank(),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold")
    )
  ggsave(file = figure_path, plot = plot, width=2560, height=1440, units = "px", scale = 2)
}

# making volcanos in a loop
joined_list_all_metas = list(
  meta_no_covar_all_brain,
  meta_no_covar_cortical,
  meta_no_covar_prefrontal,
  meta_with_covar_all_brain,
  meta_with_covar_cortical,
  meta_with_covar_prefrontal
)
names_for_list_all = c("All brain", "Cortical", "DLPFC", "All brain (covar)", "Cortical (covar)", "DLPFC (covar)")
dir.create("volcano_plots")

for (i in 1:6){
  
  curr_df = joined_list_all_metas[[i]]
  plot_path = paste0("volcano_plots/", names_for_list_all[i], ".png")
  plot_meta_volcano(curr_df, figure_path = plot_path, title = names_for_list_all[i])
  
}

images_in_folder = list.files("volcano_plots", full.names = TRUE)
image_list  =  lapply(images_in_folder, png::readPNG)
image_grobs = lapply(image_list, grid::rasterGrob)
combined_file_name = paste0("Figure_S1_Volcano_combined_image.png")
height = nrow(image_list[[1]])
width = ncol(image_list[[1]])

spacer = rectGrob(gp=gpar(col=NA, fill=NA))
image_grobs_modif = list(
  
  image_grobs[[1]], image_grobs[[2]],
  spacer, spacer,
  image_grobs[[3]], image_grobs[[4]], 
  spacer, spacer,
  image_grobs[[5]], image_grobs[[6]]
)

row_heights = unit(c(1, 3, 1, 3, 1),   # Heights for each row
                    c("null", "inches", "null", "inches", "null"))

# generating PNG
png(filename = combined_file_name, units = "px", width = width*2, height = height*3)
grid::grid.newpage()
gridExtra::grid.arrange(grobs = image_grobs_modif, ncol = 2, heights = row_heights)
dev.off()


################### Sensitivity analysis ###################

################### REML Meta (no covar) all brain data ###################
combined_df_no_covar_meta_all = combined_df_no_covar
combined_df_no_covar_meta_all = combined_df_no_covar_meta_all[combined_df_no_covar_meta_all$Tissue_type == "Brain",]
blood_df_no_covar = combined_df_no_covar[combined_df_no_covar$Tissue_type == "Blood",]

table(combined_df_no_covar_meta_all$Tissue, combined_df_no_covar_meta_all$Study)
# GSE102556 many tissues
# GSE202537 many tissues
# GSE66937 many tissues

combined_df_no_covar_meta_all_1 = combined_df_no_covar_meta_all[combined_df_no_covar_meta_all$Study %!in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_no_covar_meta_all_2 = combined_df_no_covar_meta_all[combined_df_no_covar_meta_all$Study %in% c("GSE102556", "GSE202537", "GSE66937"), ]
# Caudate https://pubmed.ncbi.nlm.nih.gov/39164232/
# Nac https://pubmed.ncbi.nlm.nih.gov/38894648/ https://pubmed.ncbi.nlm.nih.gov/38965529/ https://pubmed.ncbi.nlm.nih.gov/39167467/
# GSE102556 -> Dorsolateral prefrontal cortex (dlPFC; BA8/9)
# GSE202537 -> Nac 22652  (more probes analyzed)
# GSE66937 -> prefrontal cortex
combined_df_no_covar_meta_all_2 = combined_df_no_covar_meta_all_2[combined_df_no_covar_meta_all_2$Tissue %in% c("Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                                                "Nac",
                                                                                                                "prefrontal cortex"),]
table(combined_df_no_covar_meta_all_2$Tissue)
"
Dorsolateral prefrontal cortex (dlPFC; BA8/9)                                           Nac                             prefrontal cortex 
                                        21365                                         22652                                         70523 
"
combined_df_no_covar_meta_full = rbind(combined_df_no_covar_meta_all_1, combined_df_no_covar_meta_all_2)
table(combined_df_no_covar_meta_full$Tissue, combined_df_no_covar_meta_full$Study)
combined_df_no_covar_meta_full_list = lapply(unique(combined_df_no_covar_meta_full$Study), function(x){
  df = combined_df_no_covar_meta_full[combined_df_no_covar_meta_full$Study == x,]
  return(df)
})

combined_df_no_covar_meta_full_list_reduced = lapply(combined_df_no_covar_meta_full_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_no_covar_meta_full_list_reduced = do.call(rbind, combined_df_no_covar_meta_full_list_reduced)
unique_genes_meta = unique(combined_df_no_covar_meta_full_list_reduced$Corrected_symbol)

REML_meta_no_covar_all_brain = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_no_covar_meta_full_list_reduced[combined_df_no_covar_meta_full_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_no_covar[blood_df_no_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "REML", 
                           control=list(stepadj=0.5, maxiter=10000),
                           test="knha",  weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
REML_meta_no_covar_all_brain = do.call(rbind, REML_meta_no_covar_all_brain)
rownames(REML_meta_no_covar_all_brain) = NULL
rm(list = ls(pattern = "tmp_"))
REML_meta_no_covar_all_brain = dplyr::arrange(REML_meta_no_covar_all_brain, meta_pval)
REML_meta_no_covar_all_brain$meta_FDR = p.adjust(REML_meta_no_covar_all_brain$meta_pval, method = "fdr")
REML_meta_no_covar_all_brain = REML_meta_no_covar_all_brain[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]

REML_meta_no_covar_all_brain_signif = REML_meta_no_covar_all_brain[REML_meta_no_covar_all_brain$meta_pval < 0.05,]
REML_meta_no_covar_all_brain_signif = REML_meta_no_covar_all_brain_signif[!is.na(REML_meta_no_covar_all_brain_signif$meta_pval),]
nrow(REML_meta_no_covar_all_brain_signif)


################### REML Meta (with covar) all brain data ###################
combined_df_with_covar_meta_all = combined_df_wth_covar
combined_df_with_covar_meta_all = combined_df_with_covar_meta_all[combined_df_with_covar_meta_all$Tissue_type == "Brain",]
blood_df_with_covar = combined_df_wth_covar[combined_df_wth_covar$Tissue_type == "Blood",]

table(combined_df_with_covar_meta_all$Tissue, combined_df_with_covar_meta_all$Study)

# GSE102556 many tissues
# GSE202537 many tissues
# GSE66937 many tissues

combined_df_with_covar_meta_all_1 = combined_df_with_covar_meta_all[combined_df_with_covar_meta_all$Study %!in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_with_covar_meta_all_2 = combined_df_with_covar_meta_all[combined_df_with_covar_meta_all$Study %in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_with_covar_meta_all_2 = combined_df_with_covar_meta_all_2[combined_df_with_covar_meta_all_2$Tissue %in% c("Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                                                      "Nac",
                                                                                                                      "prefrontal cortex"),]
table(combined_df_with_covar_meta_all_2$Tissue)
"
Dorsolateral prefrontal cortex (dlPFC; BA8/9)                                           Nac                             prefrontal cortex 
                                        23331                                         23919                                         70523 
"

combined_df_with_covar_meta_full = rbind(combined_df_with_covar_meta_all_1, combined_df_with_covar_meta_all_2)
table(combined_df_with_covar_meta_full$Tissue, combined_df_with_covar_meta_full$Study)
combined_df_with_covar_meta_full_list = lapply(unique(combined_df_with_covar_meta_full$Study), function(x){
  df = combined_df_with_covar_meta_full[combined_df_with_covar_meta_full$Study == x,]
  return(df)
})

combined_df_with_covar_meta_full_list_reduced = lapply(combined_df_with_covar_meta_full_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_with_covar_meta_full_list_reduced = do.call(rbind, combined_df_with_covar_meta_full_list_reduced)
unique_genes_meta = unique(combined_df_with_covar_meta_full_list_reduced$Corrected_symbol)

REML_meta_with_covar_all_brain = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_with_covar_meta_full_list_reduced[combined_df_with_covar_meta_full_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_with_covar[blood_df_with_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "REML", 
                           control=list(stepadj=0.5, maxiter=10000),
                           test="knha",  weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
REML_meta_with_covar_all_brain = do.call(rbind, REML_meta_with_covar_all_brain)
rownames(REML_meta_with_covar_all_brain) = NULL
rm(list = ls(pattern = "tmp_"))
REML_meta_with_covar_all_brain = dplyr::arrange(REML_meta_with_covar_all_brain, meta_pval)
REML_meta_with_covar_all_brain$meta_FDR = p.adjust(REML_meta_with_covar_all_brain$meta_pval, method = "fdr")
REML_meta_with_covar_all_brain = REML_meta_with_covar_all_brain[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]

REML_meta_with_covar_all_brain_signif = REML_meta_with_covar_all_brain[REML_meta_with_covar_all_brain$meta_pval < 0.05,]
REML_meta_with_covar_all_brain_signif = REML_meta_with_covar_all_brain_signif[!is.na(REML_meta_with_covar_all_brain_signif$meta_pval),]
nrow(REML_meta_with_covar_all_brain_signif)


################### REML Meta cortical regions (no covar) ###################
combined_df_no_covar_meta_cortical = combined_df_no_covar
combined_df_no_covar_meta_cortical = combined_df_no_covar_meta_cortical[combined_df_no_covar_meta_cortical$Tissue_type == "Brain",]
blood_df_no_covar = combined_df_no_covar[combined_df_no_covar$Tissue_type == "Blood",]

table(combined_df_no_covar_meta_cortical$Tissue, combined_df_no_covar_meta_cortical$Study)
# Tissues to include
# Cingulate gyrus 25 (Cg25) is also cortical but Dorsolateral prefrontal cortex (dlPFC; BA8/9)  was used instead
combined_df_no_covar_meta_cortical = combined_df_no_covar_meta_cortical[combined_df_no_covar_meta_cortical$Tissue %in% 
                                                                          c("DLPFC",
                                                                            "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                            "Dorsolateral Prefrontal Cortex",
                                                                            "Dorsolateral prefrontal cortex (BA9)",
                                                                            "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                            "Orbitofrontal cortex",
                                                                            "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                            "prefrontal cortex",
                                                                            "temporal cortex (BA20 and BA36)"),]
table(combined_df_no_covar_meta_cortical$Tissue, combined_df_no_covar_meta_cortical$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE243356 GSE5388 GSE5389 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0         0   22283       0        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     18677         0         0         0         0         0       0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0         0       0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     19271         0     19191         0       0       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     21365         0         0         0         0       0       0        0                   0              0
  Orbitofrontal cortex                                       0         0         0         0         0         0       0   22283        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0         0       0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0         0       0       0    70523                   0              0
  temporal cortex (BA20 and BA36)                            0         0         0         0         0     20955       0       0        0                   0              0

"
combined_df_no_covar_meta_cortical_list = lapply(unique(combined_df_no_covar_meta_cortical$Study), function(x){
  df = combined_df_no_covar_meta_cortical[combined_df_no_covar_meta_cortical$Study == x,]
  return(df)
})

combined_df_no_covar_meta_cortical_list_reduced = lapply(combined_df_no_covar_meta_cortical_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_no_covar_meta_cortical_list_reduced = do.call(rbind, combined_df_no_covar_meta_cortical_list_reduced)
unique_genes_meta = unique(combined_df_no_covar_meta_cortical_list_reduced$Corrected_symbol)

REML_meta_no_covar_cortical = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_no_covar_meta_cortical_list_reduced[combined_df_no_covar_meta_cortical_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_no_covar[blood_df_no_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "REML", 
                           control=list(stepadj=0.5, maxiter=10000),
                           test="knha",  weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
REML_meta_no_covar_cortical = do.call(rbind, REML_meta_no_covar_cortical)
rownames(REML_meta_no_covar_cortical) = NULL
rm(list = ls(pattern = "tmp_"))
REML_meta_no_covar_cortical = dplyr::arrange(REML_meta_no_covar_cortical, meta_pval)
REML_meta_no_covar_cortical$meta_FDR = p.adjust(REML_meta_no_covar_cortical$meta_pval, method = "fdr")
REML_meta_no_covar_cortical = REML_meta_no_covar_cortical[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]

REML_meta_no_covar_cortical_signif = REML_meta_no_covar_cortical[REML_meta_no_covar_cortical$meta_pval < 0.05,]
REML_meta_no_covar_cortical_signif = REML_meta_no_covar_cortical_signif[!is.na(REML_meta_no_covar_cortical_signif$meta_pval),]
nrow(REML_meta_no_covar_cortical_signif)


################### REML Meta cortical regions (with covar) ###################
combined_df_with_covar_meta_cortical = combined_df_wth_covar
combined_df_with_covar_meta_cortical = combined_df_with_covar_meta_cortical[combined_df_with_covar_meta_cortical$Tissue_type == "Brain",]
blood_df_with_covar = combined_df_wth_covar[combined_df_wth_covar$Tissue_type == "Blood",]

# Tissues to include
combined_df_with_covar_meta_cortical = combined_df_with_covar_meta_cortical[combined_df_with_covar_meta_cortical$Tissue %in% 
                                                                              c("DLPFC",
                                                                                "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                                "Dorsolateral Prefrontal Cortex",
                                                                                "Dorsolateral prefrontal cortex (BA9)",
                                                                                "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                "Orbitofrontal cortex",
                                                                                "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                                "prefrontal cortex",
                                                                                "temporal cortex (BA20 and BA36)"),]
table(combined_df_with_covar_meta_cortical$Tissue, combined_df_with_covar_meta_cortical$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE243356 GSE5388 GSE5389 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0         0   22283       0        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     21197         0         0         0         0         0       0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0         0       0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     22389         0     24182         0       0       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     23331         0         0         0         0       0       0        0                   0              0
  Orbitofrontal cortex                                       0         0         0         0         0         0       0   22283        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0         0       0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0         0       0       0    70523                   0              0
  temporal cortex (BA20 and BA36)                            0         0         0         0         0     20955       0       0        0                   0              0

"

combined_df_with_covar_meta_cortical_list = lapply(unique(combined_df_with_covar_meta_cortical$Study), function(x){
  df = combined_df_with_covar_meta_cortical[combined_df_with_covar_meta_cortical$Study == x,]
  return(df)
})

combined_df_with_covar_meta_cortical_list_reduced = lapply(combined_df_with_covar_meta_cortical_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_with_covar_meta_cortical_list_reduced = do.call(rbind, combined_df_with_covar_meta_cortical_list_reduced)
unique_genes_meta = unique(combined_df_with_covar_meta_cortical_list_reduced$Corrected_symbol)

REML_meta_with_covar_cortical = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_with_covar_meta_cortical_list_reduced[combined_df_with_covar_meta_cortical_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_with_covar[blood_df_with_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "REML", 
                           control=list(stepadj=0.5, maxiter=10000),
                           test="knha",  weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
REML_meta_with_covar_cortical = do.call(rbind, REML_meta_with_covar_cortical)
rownames(REML_meta_with_covar_cortical) = NULL
rm(list = ls(pattern = "tmp_"))
REML_meta_with_covar_cortical = dplyr::arrange(REML_meta_with_covar_cortical, meta_pval)
REML_meta_with_covar_cortical$meta_FDR = p.adjust(REML_meta_with_covar_cortical$meta_pval, method = "fdr")
REML_meta_with_covar_cortical = REML_meta_with_covar_cortical[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]

REML_meta_with_covar_cortical_signif = REML_meta_with_covar_cortical[REML_meta_with_covar_cortical$meta_pval < 0.05,]
REML_meta_with_covar_cortical_signif = REML_meta_with_covar_cortical_signif[!is.na(REML_meta_with_covar_cortical_signif$meta_pval),]
nrow(REML_meta_with_covar_cortical_signif)


################### REML Meta prefrontal regions (no covar) ###################
combined_df_no_covar_meta_prefrontal = combined_df_no_covar
combined_df_no_covar_meta_prefrontal = combined_df_no_covar_meta_prefrontal[combined_df_no_covar_meta_prefrontal$Tissue_type == "Brain",]
blood_df_no_covar = combined_df_no_covar[combined_df_no_covar$Tissue_type == "Blood",]

table(combined_df_no_covar_meta_prefrontal$Tissue, combined_df_no_covar_meta_prefrontal$Study)
# Tissues to include
combined_df_no_covar_meta_prefrontal = combined_df_no_covar_meta_prefrontal[combined_df_no_covar_meta_prefrontal$Tissue %in% 
                                                                              c("DLPFC",
                                                                                "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                                "Dorsolateral Prefrontal Cortex",
                                                                                "Dorsolateral prefrontal cortex (BA9)",
                                                                                "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                                "prefrontal cortex"),]
table(combined_df_no_covar_meta_prefrontal$Tissue, combined_df_no_covar_meta_prefrontal$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE5388 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0   22283        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     18677         0         0         0         0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     19271         0     19191       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     21365         0         0         0       0        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0       0    70523                   0              0
  
"
combined_df_no_covar_meta_prefrontal_list = lapply(unique(combined_df_no_covar_meta_prefrontal$Study), function(x){
  df = combined_df_no_covar_meta_prefrontal[combined_df_no_covar_meta_prefrontal$Study == x,]
  return(df)
})

combined_df_no_covar_meta_prefrontal_list_reduced = lapply(combined_df_no_covar_meta_prefrontal_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_no_covar_meta_prefrontal_list_reduced = do.call(rbind, combined_df_no_covar_meta_prefrontal_list_reduced)
unique_genes_meta = unique(combined_df_no_covar_meta_prefrontal_list_reduced$Corrected_symbol)

REML_meta_no_covar_prefrontal = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_no_covar_meta_prefrontal_list_reduced[combined_df_no_covar_meta_prefrontal_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_no_covar[blood_df_no_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "REML", 
                           control=list(stepadj=0.5, maxiter=10000),
                           test="knha",  weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
REML_meta_no_covar_prefrontal = do.call(rbind, REML_meta_no_covar_prefrontal)
rownames(REML_meta_no_covar_prefrontal) = NULL
rm(list = ls(pattern = "tmp_"))
REML_meta_no_covar_prefrontal = dplyr::arrange(REML_meta_no_covar_prefrontal, meta_pval)
REML_meta_no_covar_prefrontal$meta_FDR = p.adjust(REML_meta_no_covar_prefrontal$meta_pval, method = "fdr")
REML_meta_no_covar_prefrontal = REML_meta_no_covar_prefrontal[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]
REML_meta_no_covar_prefrontal_signif = REML_meta_no_covar_prefrontal[REML_meta_no_covar_prefrontal$meta_pval < 0.05,]
REML_meta_no_covar_prefrontal_signif = REML_meta_no_covar_prefrontal_signif[!is.na(REML_meta_no_covar_prefrontal_signif$meta_pval),]
nrow(REML_meta_no_covar_prefrontal_signif)


################### REML Meta prefrontal regions (with covar) ###################
combined_df_with_covar_meta_prefrontal = combined_df_wth_covar
combined_df_with_covar_meta_prefrontal = combined_df_with_covar_meta_prefrontal[combined_df_with_covar_meta_prefrontal$Tissue_type == "Brain",]
blood_df_with_covar = combined_df_wth_covar[combined_df_wth_covar$Tissue_type == "Blood",]

# Tissues to include
combined_df_with_covar_meta_prefrontal = combined_df_with_covar_meta_prefrontal[combined_df_with_covar_meta_prefrontal$Tissue %in% 
                                                                                  c("DLPFC",
                                                                                    "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                                    "Dorsolateral Prefrontal Cortex",
                                                                                    "Dorsolateral prefrontal cortex (BA9)",
                                                                                    "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                    "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                                    "prefrontal cortex"),]
table(combined_df_with_covar_meta_prefrontal$Tissue, combined_df_with_covar_meta_prefrontal$Study)

"
                                                     GSE101521 GSE102556 GSE144136 GSE208338 GSE213982 GSE5388 GSE66937 GSE92538_U133_PLUS2 GSE92538_U133A
  DLPFC                                                      0         0         0         0         0   22283        0                   0              0
  dorsal lateral prefrontal cortex (Brodmann Area 9)     21197         0         0         0         0       0        0                   0              0
  Dorsolateral Prefrontal Cortex                             0         0         0         0         0       0        0               54675          22283
  Dorsolateral prefrontal cortex (BA9)                       0         0     22389         0     24182       0        0                   0              0
  Dorsolateral prefrontal cortex (dlPFC; BA8/9)              0     23331         0         0         0       0        0                   0              0
  postmortem dorsolateral prefrontal cortex (DLPFC)          0         0         0    750011         0       0        0                   0              0
  prefrontal cortex                                          0         0         0         0         0       0    70523                   0              0
  
"

combined_df_with_covar_meta_prefrontal_list = lapply(unique(combined_df_with_covar_meta_prefrontal$Study), function(x){
  df = combined_df_with_covar_meta_prefrontal[combined_df_with_covar_meta_prefrontal$Study == x,]
  return(df)
})

combined_df_with_covar_meta_prefrontal_list_reduced = lapply(combined_df_with_covar_meta_prefrontal_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxSE = max(SE))
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  return(x)
})
combined_df_with_covar_meta_prefrontal_list_reduced = do.call(rbind, combined_df_with_covar_meta_prefrontal_list_reduced)
unique_genes_meta = unique(combined_df_with_covar_meta_prefrontal_list_reduced$Corrected_symbol)

REML_meta_with_covar_prefrontal = mclapply(unique_genes_meta, function(x){
  
  
  # meta df
  tmp_df_genes = combined_df_with_covar_meta_prefrontal_list_reduced[combined_df_with_covar_meta_prefrontal_list_reduced$Corrected_symbol == x,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  # blood df
  tmp_df_blood = blood_df_with_covar[blood_df_with_covar$Corrected_symbol == x, ]
  
  if (nrow(tmp_df_blood) > 0){
    
    tmp_blood_lfc = mean(tmp_df_blood$logFC)
    
    if (all(tmp_df_blood$logFC < 0)){
      tmp_lfc_blood_dir = "all negative"
    } else {
      tmp_lfc_blood_dir = "all positive"
    }
    
    if (all(tmp_df_blood$P.Value < 0.05)){
      tmp_blood_signif = "significant"
    } else {
      tmp_blood_signif = "mixed/non-significant"
    }
    
  } else {
    tmp_blood_lfc = NA
    tmp_lfc_blood_dir = "Not detected/analyzed"
    tmp_blood_signif = "Not detected/analyzed"
  }
  
  if (cohorts_total < 5){
    tmp_output_df = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts < 5: gene is excluded",
      meta_LFc = NA,
      meta_se = NA,
      meta_pval = NA,
      tau2 = NA,
      I2 = NA,
      H2 = NA,
      Q = NA,
      Q.p = NA,
      mean_blood_lfc = tmp_blood_lfc,
      blood_dir = tmp_lfc_blood_dir,
      blood_signif = tmp_blood_signif,
      matching_with_blood = NA)
    return(tmp_output_df)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "REML", 
                           control=list(stepadj=0.5, maxiter=10000),
                           test="knha",  weighted = TRUE)
  
  logFCcompar = c(tmp_blood_lfc, tmp_meta_model$b)
  
  if (any(is.na(logFCcompar))){
    matching = NA
  } else if (all(logFCcompar > 0)){
    matching = "Match"
  } else if (all(logFCcompar < 0)){
    matching = "Match"
  } else {
    matching = "Mismatch"
  }
  
  tmp_output_df = data.frame(
    gene = x,
    detected_cohorts = cohorts_total,
    cohorts_up = tmp_cohorts_up,
    cohorts_down = tmp_cohorts_down,
    commentary = "Number of cohorts >= 5: gene is analyzed",
    meta_LFc = tmp_meta_model$b,
    meta_se = tmp_meta_model$se,
    meta_pval = tmp_meta_model$pval,
    tau2 = tmp_meta_model$tau2,
    I2 = tmp_meta_model$I2,
    H2 = tmp_meta_model$H2,
    Q = tmp_meta_model$QE,
    Q.p = tmp_meta_model$QEp,
    mean_blood_lfc = tmp_blood_lfc,
    blood_dir = tmp_lfc_blood_dir,
    blood_signif = tmp_blood_signif,
    matching_with_blood = matching)
  
  return(tmp_output_df)
  
}, mc.cores = 9)
REML_meta_with_covar_prefrontal = do.call(rbind, REML_meta_with_covar_prefrontal)
rownames(REML_meta_with_covar_prefrontal) = NULL
rm(list = ls(pattern = "tmp_"))
REML_meta_with_covar_prefrontal = dplyr::arrange(REML_meta_with_covar_prefrontal, meta_pval)
REML_meta_with_covar_prefrontal$meta_FDR = p.adjust(REML_meta_with_covar_prefrontal$meta_pval, method = "fdr")
REML_meta_with_covar_prefrontal = REML_meta_with_covar_prefrontal[,c("gene",
                                                         "detected_cohorts",
                                                         "cohorts_up",
                                                         "cohorts_down",
                                                         "commentary",
                                                         "meta_LFc",
                                                         "meta_se",
                                                         "meta_pval",
                                                         "meta_FDR",
                                                         "tau2",
                                                         "I2",
                                                         "H2" ,
                                                         "Q",
                                                         "Q.p",
                                                         "mean_blood_lfc",
                                                         "blood_dir",
                                                         "blood_signif",
                                                         "matching_with_blood")]
REML_meta_with_covar_prefrontal_signif = REML_meta_with_covar_prefrontal[REML_meta_with_covar_prefrontal$meta_pval < 0.05,]
REML_meta_with_covar_prefrontal_signif = REML_meta_with_covar_prefrontal_signif[!is.na(REML_meta_with_covar_prefrontal_signif$meta_pval),]
nrow(REML_meta_with_covar_prefrontal_signif)

################### Saving sensitivity ###################
significant_genes_var = ls(pattern = "REML_meta_")
significant_genes_var = significant_genes_var[stri_detect_regex(significant_genes_var, pattern = "signif")]
significant_genes_var = significant_genes_var[!stri_detect_regex(significant_genes_var, pattern = "meta_list")]
significant_genes_counts = sapply(significant_genes_var, function(x) nrow(get(x)))

stats_meta = paste0("Signifi gene counts for: ", significant_genes_var)
stats_meta = paste0(stats_meta, " = ", significant_genes_counts)
writeLines(stats_meta, "stats_meta_REML.txt")

# saving files
REML_meta_list_no_covar_signif = list(
  "REML_meta_no_covar_all_brain_signif" = REML_meta_no_covar_all_brain_signif,
  "REML_meta_no_covar_cortical_signif" = REML_meta_no_covar_cortical_signif,
  "REML_meta_no_covar_prefrontal_signif" = REML_meta_no_covar_prefrontal_signif
)

REML_meta_list_with_covar_signif = list(
  "REML_meta_with_covar_all_brain_signif" = REML_meta_with_covar_all_brain_signif,
  "REML_meta_with_covar_cortical_signif" = REML_meta_with_covar_cortical_signif,
  "REML_meta_with_covar_prefrontal_signif" = REML_meta_with_covar_prefrontal_signif
)

REML_joined_list_significant = c(REML_meta_list_no_covar_signif, REML_meta_list_with_covar_signif)
REML_joined_list_significant = lapply(REML_joined_list_significant, function(x){
  x = dplyr::arrange(x, -meta_LFc)
  return(x)
})

names_for_list = c(
  "REML All brain",
  "REML Cortical regions",
  "REML DLPFC",
  "REML All brain (covar)",
  "REML Cortical regions (covar)",
  "REML DLPFC (covar)"
)
names(REML_joined_list_significant)

wb = createWorkbook()

for (i in 1:6) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, names_for_list[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = names_for_list[i], REML_joined_list_significant[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = names_for_list[i], cols = 1:ncol(REML_joined_list_significant[[i]]), widths = "auto")
}

saveWorkbook(wb, "Meta_suicide_significant_genes_REML.xlsx", overwrite = TRUE)

################### Other analyses ? ###################

table(combined_df_no_covar$Tissue, combined_df_no_covar$Study)

# Amygdala vs hippocampus vs thalamus
# Anterior Insula (aINS) vs Cingulate gyrus 25 (Cg25) vs  Nucleus Accumbens (Nac) vs Subiculum (Sub)

# Striatum: caudate nucleus and the putamen; The nucleus accumbens (NAc) is a major component of the ventral striatum 
# The subiculum (Latin for "support") is the most inferior component of the hippocampal formation
# The striatum sends projections directly to the medial habenula https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3894476/

################### Stats on Meta ###################
significant_genes_var = ls(pattern = "meta")
significant_genes_var = significant_genes_var[stri_detect_regex(significant_genes_var, pattern = "signif")]
significant_genes_var = significant_genes_var[!stri_detect_regex(significant_genes_var, pattern = "meta_list")]
significant_genes_counts = sapply(significant_genes_var, function(x) nrow(get(x)))

stats_meta = paste0("Signifi gene counts for: ", significant_genes_var)
stats_meta = paste0(stats_meta, " = ", significant_genes_counts)
writeLines(stats_meta, "stats_meta.txt")


################### Aggregating data ###################

meta_list_no_covar = list(
  "meta_no_covar_all_brain" = meta_no_covar_all_brain,
  "meta_no_covar_cortical" = meta_no_covar_cortical,
  "meta_no_covar_prefrontal" = meta_no_covar_prefrontal
)

meta_list_with_covar = list(
  "meta_with_covar_all_brain" = meta_with_covar_all_brain,
  "meta_with_covar_cortical" = meta_with_covar_cortical,
  "meta_with_covar_prefrontal" = meta_with_covar_prefrontal
)

meta_list_no_covar_signif = list(
  "meta_no_covar_all_brain_signif" = meta_no_covar_all_brain_signif,
  "meta_no_covar_cortical_signif" = meta_no_covar_cortical_signif,
  "meta_no_covar_prefrontal_signif" = meta_no_covar_prefrontal_signif
)

meta_list_with_covar_signif = list(
  "meta_with_covar_all_brain_signif" = meta_with_covar_all_brain_signif,
  "meta_with_covar_cortical_signif" = meta_with_covar_cortical_signif,
  "meta_with_covar_prefrontal_signif" = meta_with_covar_prefrontal_signif
)

################### Comparison with PSY SUAS ###################
suas_genes_PSY = c("PLA2G10", "RBKS", "CCL27", "ASRG1", "WWP2")
meta_no_covar_all_brain[meta_no_covar_all_brain$gene %in% suas_genes_PSY,]
meta_with_covar_all_brain[meta_with_covar_all_brain$gene %in% suas_genes_PSY,]
blood_df_no_covar[blood_df_no_covar$Corrected_symbol %in% suas_genes_PSY,]
blood_df_with_covar[blood_df_with_covar$Corrected_symbol %in% suas_genes_PSY,]
# None are validated in full analyses and in blood

################### Stats on matching with blood ###################

meta_list_signif_combined = list(
  "meta_no_covar_all_brain_signif" = meta_no_covar_all_brain_signif,
  "meta_no_covar_cortical_signif" = meta_no_covar_cortical_signif,
  "meta_no_covar_prefrontal_signif" = meta_no_covar_prefrontal_signif,
  "meta_with_covar_all_brain_signif" = meta_with_covar_all_brain_signif,
  "meta_with_covar_cortical_signif" = meta_with_covar_cortical_signif,
  "meta_with_covar_prefrontal_signif" = meta_with_covar_prefrontal_signif
)

tabulate_variable = function(factor_var, useNA = "always"){
  
  # Calculate counts for each level
  counts = table(factor_var,  useNA = useNA)
  
  # Calculate percentages
  percentages = round((counts / sum(counts)) * 100, 2)
  
  # Create the result as counts with percentages in parentheses
  result = paste0(counts, " (", percentages, "%)")
  
  # Name the results with factor levels
  names(result) = names(counts)
  
  result = paste0(names(result), ": ", result)
  
  result = paste0(result, collapse = "\n")
  
  return(result)
}

compare_with_blood = function(dataset, name){
  
  Data_name = name
  Blood_dir = paste0("Blood direction (probe-level): \n", tabulate_variable(dataset$blood_dir))
  Blood_signif = paste0("Blood significance (probe-level): \n", tabulate_variable(dataset$blood_signif))
  Matching_brain_blood = paste0("Directional match brain and blood (only for analyzed genes): \n", tabulate_variable(dataset$matching_with_blood, useNA="no"))
  
  Matched_genes = dataset[dataset$matching_with_blood == "Match",]
  Matched_genes = Matched_genes$gene
  Matched_genes = Matched_genes[!is.na(Matched_genes)]
  Matched_genes = paste0(Matched_genes, collapse = ";")
  Matched_genes = paste0("Matched genes: ", Matched_genes)
  
  Matched_genes_blood_signif = dataset[dataset$matching_with_blood == "Match", ]
  Matched_genes_blood_signif = Matched_genes_blood_signif[Matched_genes_blood_signif$blood_signif == "significant",]
  Matched_genes_blood_signif = Matched_genes_blood_signif$gene
  Matched_genes_blood_signif = Matched_genes_blood_signif[!is.na(Matched_genes_blood_signif)]
  Matched_genes_blood_signif = paste0(Matched_genes_blood_signif, collapse = ";")
  Matched_genes_blood_signif = paste0("Matched genes (signif in blood): ", Matched_genes_blood_signif)
  
  Mismatched_genes = dataset[dataset$matching_with_blood == "Mismatch",]
  Mismatched_genes = Mismatched_genes$gene
  Mismatched_genes = Mismatched_genes[!is.na(Mismatched_genes)]
  Mismatched_genes = paste0(Mismatched_genes, collapse = ";")
  Mismatched_genes = paste0("Mismatched genes: ", Mismatched_genes)
  
  separator_small = "\n\n"
  separator_large = "\n\n\n\n"
  output = c(Data_name, 
             separator_small,
             Blood_dir,
             separator_small,
             Blood_signif,
             separator_small,
             Matching_brain_blood,
             separator_small,
             Matched_genes,
             separator_small,
             Matched_genes_blood_signif,
             separator_small,
             Mismatched_genes,
             separator_large)
  output = paste0(output, collapse = "")
  return(output)
}

names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)"
)

Stats_matching_brain_blood_meta = vector()

for (i in 1:length(meta_list_signif_combined)){
  
  message = compare_with_blood(meta_list_signif_combined[[i]], names_for_list[i])
  Stats_matching_brain_blood_meta = c(Stats_matching_brain_blood_meta, message)
  
}
writeLines(text = Stats_matching_brain_blood_meta, con = "Stats_matching_brain_blood_meta.txt")


################### Venn diagrams ###################
# Without covariates
meta_gene_list_overlaps = lapply(meta_list_no_covar_signif, function(x){
  x = x$gene
  return(x)
})

names(meta_gene_list_overlaps) = c(
  "All brain",
  "Cortical regions",
  "DLPFC"
)

make_Venn_digram_list(named_list = meta_gene_list_overlaps, palette = 2, plot_full_path = "Venn_no_covar.pdf")

# With covariates
meta_gene_list_overlaps_covar = lapply(meta_list_with_covar_signif, function(x){
  x = x$gene
  return(x)
})

names(meta_gene_list_overlaps_covar) = c(
  "All brain",
  "Cortical regions",
  "DLPFC"
)

make_Venn_digram_list(named_list = meta_gene_list_overlaps_covar, palette = 7, plot_full_path = "Venn_with_covar.pdf")

################### Network visualization for meta ###################
joined_list = c(meta_list_no_covar_signif, meta_list_with_covar_signif)
joined_list = lapply(joined_list, function(x){
  x = x[abs(x$meta_LFc) >= 0.2,]
  return(x)
})

names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)"
)

meta_gene_joined_list = lapply(1:6, function(index){
  data_list = joined_list[[index]]
  name = names_for_list[index]
  output = data.frame(gene = data_list$gene, effect_size = data_list$meta_LFc, analysis = name)
  return(output)
})
meta_gene_joined_list = do.call(rbind, meta_gene_joined_list)
combined_graph_df = meta_gene_joined_list[,c("analysis", "gene")]
colnames(combined_graph_df) = c("Start", "Target")

# preparing data
network_data = graph_from_data_frame(d = combined_graph_df, directed = TRUE)
network_data = toVisNetworkData(network_data)

network_data$nodes$group = sapply(network_data$nodes$id, function(x){
if (stri_detect_fixed(x, pattern = "(covar)")){
  return("Analysis Covar")
}
if (x %in% names_for_list){
  return("Analysis")
} else {
  return("Gene")
}
})

network_data$edges$width = mapply(function(x,y){
  Curr_edge = x
  sub_df = meta_gene_joined_list[meta_gene_joined_list$analysis == x, ]
  beta = sub_df[sub_df$gene == y,"effect_size"]
  curr_parameter = beta*20
  return(curr_parameter)
}, network_data$edges$from, network_data$edges$to)

network_data$edges$color = sapply(network_data$edges$width, function(x){
  if (x < 0){
    return("#00FF00")
  } else {
    return("#AA4A44")
  }
})
network_data$edges$width = abs(network_data$edges$width)
net = visNetwork(nodes = network_data$nodes, edges = network_data$edges, height = "2000px", width = "2000px") %>%
  visIgraphLayout(layout = "layout_with_fr", physics = TRUE, randomSeed = 190) %>%
  visGroups(groupname = "Analysis", color = list("background" = "orange", border = "black")) %>% 
  visGroups(groupname = "Analysis Covar", color = list("background" = "blue", border = "black")) %>% 
  visGroups(groupname = "Gene", color = list("background" = "red", border = "black")) %>% 
  visNodes(size = 15, borderWidth = 1, font = list(size=30)) %>%
  visEdges(smooth = list("roundness" = 0.2, "type" = "diagonalCross"), arrows = list(to = list(enabled = TRUE, scaleFactor = 1)), color = "grey") %>%
  visPhysics(barnesHut = list(gravitationalConstant = -50000, centralGravity=0.001)) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = F)
visSave(net, file = "Meta_net.html")
Sys.setenv("OPENSSL_CONF"="/dev/null")
webshot("Meta_net.html", "Meta_net.png", vwidth = 2000, vheight = 2000, zoom = 1, delay = 10)


################### Enrichment analysis ###################
# Preparing Entrez dataset
# https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/, https://www.ncbi.nlm.nih.gov/gene/
ENTREZ_genes_Homo_Sapiens = smart_fread("/home/aleksandr/Desktop/WORK/GWAS_CAT_CHARACT/data") # Replace with another path if needed
ENTREZ_genes_Homo_Sapiens = ENTREZ_genes_Homo_Sapiens[ENTREZ_genes_Homo_Sapiens$`#tax_id` == "9606",]
ENTREZ_genes_Homo_Sapiens$Ensembl = sapply(ENTREZ_genes_Homo_Sapiens$dbXrefs, function(x){
  x = unlist(stri_split_fixed(x, pattern = "|"))
  x = x[stri_detect_fixed(str = x, pattern = "Ensembl")]
  x = stri_replace_all_fixed(str = x, pattern = "Ensembl:", replacement = "")
  names(x) = NULL
  
  if (length(x) > 1){
    x = paste0(x, collapse = ";")
  }
  
  if (length(x) > 0){
    return(x)
  } else {
    return(NA)
  }
  
})
ENTREZ_genes_Homo_Sapiens$chromosome = paste0("chr", ENTREZ_genes_Homo_Sapiens$chromosome)
ENTREZ_Meta_universe = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% harmoniz_genes_check$Suggested.Symbol){
    return(TRUE)
  }
  return(FALSE)
}),]

joined_list_enrichment = c(meta_list_no_covar_signif, meta_list_with_covar_signif)
names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)"
)
joined_list_enrichment = lapply(joined_list_enrichment, function(x){
  data = x$gene
  return(data)
})

# running enrichment
for (x in 1:6){
  print(paste0("Working on: ", names_for_list[x]))
  current_genes_tmp = joined_list_enrichment[[x]]
  ENTREZ_current_genes_tmp = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
    if (x %in% current_genes_tmp){
      return(TRUE)
    }
    return(FALSE)
  }),]
  
  # Enrichment
  ENRICHMENT_current_genes_tmp = run_enrichment_GO_KEGG_gene_set(genes = ENTREZ_current_genes_tmp$GeneID,
                                                                 universe = ENTREZ_Meta_universe$GeneID,
                                                                 folder = "Enrichment_folder",
                                                                 plot_name_pref = names_for_list[x])
  
}


enrichment_files = list.files(path = "Enrichment_folder", pattern = "_GO_BP", full.names = TRUE)
enrichment_files = enrichment_files[stri_detect_fixed(enrichment_files, pattern = ".xlsx")]
enrichment_file_names = multiple_stri_replacer(string = enrichment_files, pattern_vector = c("Enrichment_folder/", "_GO_BP.xlsx"), replacement_vector = c("",""))
enrichment_files_BP_list = lapply(enrichment_files, function(x){
  x = read.xlsx(x, sheet = 1)
  return(x)
  })
names(enrichment_files_BP_list) = enrichment_file_names

wb = createWorkbook()

for (sheet_name in names(enrichment_files_BP_list)) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, sheet_name)
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = sheet_name, enrichment_files_BP_list[[sheet_name]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = sheet_name, cols = 1:ncol(enrichment_files_BP_list[[sheet_name]]), widths = "auto")
}

saveWorkbook(wb, "Enrichment_BP.xlsx", overwrite = TRUE)

# resaving some images as PNG
library("animation")
im.convert("Enrichment_folder/DLPFC_GO_BP_bar.pdf", convert = "convert",output = "DLPFC_GO_BP_bar.png", extra.opts="-density 300x300 -units pixelsperinch ")
im.convert("Enrichment_folder/DLPFC_GO_CC_bar.pdf", convert = "convert",output = "DLPFC_GO_CC_bar.png", extra.opts="-density 300x300 -units pixelsperinch ")

################### AI-based classification and stats ###################
joined_list_significant = c(meta_list_no_covar_signif, meta_list_with_covar_signif)
joined_list_significant = lapply(joined_list_significant, function(x){
  x = dplyr::arrange(x, -meta_LFc)
  return(x)
})
unique_significant_genes = lapply(joined_list_significant, function(x){
  data = x$gene
  return(data)
})
unique_significant_genes = unlist(unique_significant_genes)
unique_significant_genes = unique(unique_significant_genes)
unique_significant_genes_concat = paste0(unique_significant_genes, collapse = ";")

split_into_chunks = function(vec, n) {
  split(vec, ceiling(seq_along(vec) / n))
}

unique_significant_genes_chunks = split_into_chunks(unique_significant_genes, 300)
unique_significant_genes_chunks = lapply(unique_significant_genes_chunks, function(x) paste0(x, collapse = ";"))
unique_significant_genes_chunks = unlist(unique_significant_genes_chunks)
unique_significant_genes_chunks

# Classification based on ChatGPT o1-preview

prompt_chat =
"
You are given the list of human gene names separated by semicolon ; . You need to do the following:

1. Classify the genes using these classes: 
'Receptor(non-immune)', 
'Immune receptor', 
'Cytokine', 
'Ligand', 
'Enzyme (non-kinase)',
'Transcription factor',
'Kinase',
'Structural protein', 
'Non-coding RNA',
'Pseudogene'

2. Only use scientific literature (papers) or scientific databases to perform classification

3. Give 1 the most suitable class for every gene

4. DON'T SKIP ANY GENES and provide classification for all of them

5. Provide output as CSV in one string where each row is as follows:

gene,class;\n

note there should be a semicolon ; after every gene,class pair

The list of genes: 

"
prompt_chat_full = paste0(prompt_chat, "\n", unique_significant_genes_chunks)

# Importing data

GPT01_classification = readLines("GPT01_classification.txt")
GPT01_classification = str_trim(GPT01_classification)
GPT01_classification = unlist(stri_split_fixed(GPT01_classification, pattern = ";"))
GPT01_classification = sapply(GPT01_classification, function(x){
  x = unlist(stri_split_fixed(x, pattern = ","))
  return(x)
})
GPT01_classification = do.call(rbind, GPT01_classification)
rownames(GPT01_classification) = NULL
colnames(GPT01_classification) = c("Gene", "GPT01_class")
GPT01_classification = as.data.frame(GPT01_classification)
GPT01_classification$Gene = str_trim(GPT01_classification$Gene)
GPT01_classification$GPT01_class = str_trim(GPT01_classification$GPT01_class)
table(unique_significant_genes %in% GPT01_classification$Gene) # TRUE 1403
GPT01_classification = GPT01_classification[GPT01_classification$Gene != "",]

# Extra mapping
GPT01_classification_unmapped = unique_significant_genes[unique_significant_genes %!in% GPT01_classification$Gene]
GPT01_classification_unmapped = paste0(GPT01_classification_unmapped, collapse = ";")
prompt_chat_remapping = paste0(prompt_chat, "\n", GPT01_classification_unmapped)
writeLines(prompt_chat_remapping)

# Mapping to classes and generating stats
joined_list_significant = c(meta_list_no_covar_signif, meta_list_with_covar_signif)
joined_list_significant = lapply(joined_list_significant, function(x){
  x = dplyr::arrange(x, -meta_LFc)
  x$GPT01_class = sapply(x$gene, function(symbol){
    value = GPT01_classification[GPT01_classification$Gene == symbol, "GPT01_class"]
    return(value)
  })
  x$GPT01_class = factor(x$GPT01_class)
  return(x)
})


joined_list_significant_high_LFc = lapply(joined_list_significant, function(x){
  x = dplyr::arrange(x, -meta_LFc)
  x = x[abs(x$meta_LFc) >= 0.2,]
  return(x)
})

names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)"
)

joined_list_significant_high_LFc = lapply(1:6, function(x){
  df = joined_list_significant_high_LFc[[x]]
  df$analysis = names_for_list[x]
  return(df)
})

stats_high_LFc = do.call(rbind, joined_list_significant_high_LFc)
stats_high_LFc = as.data.frame(table(stats_high_LFc$GPT01_class, stats_high_LFc$analysis))
colnames(stats_high_LFc) = c("Category", "Analysis", "Freq")
stats_high_LFc$Category = stri_replace_all_fixed(stats_high_LFc$Category, pattern = "Receptor(non-immune)", replacement = "Receptor (non-immune)")

plot = ggplot(data = stats_high_LFc, aes(x=Analysis, y = Freq, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, col = "black", size = 0.25) +
  labs(y = "Gene count", x = "Analysis type", fill = "Gene category (GPTo1 preview)", title = "Gene categories with abs(logFC)>=0.2") +
  scale_fill_brewer(palette = "Set3") +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.8), vjust = -0.5, size = 3, family = "Times") +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    panel.grid.major.y = element_line(color = "black", linewidth = 0.25, linetype = 2),
    axis.title = element_text(family = "Times"),
    legend.title = element_text(family = "Times"),
    plot.title = element_text(family = "Times", face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, family = "Times", color = "#5d5d5d"),
    legend.text = element_text(family = "Times", color = "#5d5d5d")
  )
  
ggsave(file = "stats_high_LFc.pdf", plot = plot, width = 20, height = 15, units = "cm")


################### Saving processed files ###################
names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)"
)
names(joined_list_significant)

wb = createWorkbook()

for (i in 1:6) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, names_for_list[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = names_for_list[i], joined_list_significant[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = names_for_list[i], cols = 1:ncol(joined_list_significant[[i]]), widths = "auto")
}

saveWorkbook(wb, "Meta_suicide_significant_genes.xlsx", overwrite = TRUE)

################### Saving full files ###################
joined_list_meta_all_data = c(meta_list_no_covar, meta_list_with_covar)
names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)"
)
names(joined_list_meta_all_data)

dir.create("full_meta_results")

for (i in 1:6) {
  
  path_to_save = paste0("full_meta_results/",names_for_list[i], ".csv")
  
  # Save data as CSV
  write.csv(joined_list_meta_all_data[[i]], file = path_to_save, row.names = FALSE)
}
