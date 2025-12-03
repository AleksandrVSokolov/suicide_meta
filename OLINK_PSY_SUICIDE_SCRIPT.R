# The script is intended viewed in R studio
# Blocks are titled
# "=" symbol was used as an assignment operator
# Data processing and DE in initial cohorts is located in another file "Data_preprocessing_analysis/TRANSCRIPT_SUICIDE_PREPR_ANALYSIS_SCRIPT.R"
# Calculation of cohort-level moderators is located in "Moderator_calculation.R"
# Cell expression deconvolution tests are located in "Cell_expression_imputation_tests.R"
# Cell expression deconvolution final runs and analysis are located in "Cell_expression_imputation_analysis.R"


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
library(RobustRankAggreg)

################### Defining basic functions ###################

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
  difference_coords = sav$xlim[2] - sav$xlim[1]
  # Headers
  text(sav$xlim[1], k+2.5, pos=4, "Cohort")
  text(sav$xlim[1]+difference_coords*0.25, k+2.5, pos=4, "Weight %")
  text(0, k+2.7, "Log2FC,\n(95% CI)")
  segments(sav$ilab.xpos[1]-0.22, k+2.5, sav$ilab.xpos[2]+0.13, k+2.8)
  text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")
  
  # Use a non-bold font for the rest of the text
  # Y coord for heterogen
  y_low = sav$ylim[1]
  y_high = sav$ylim[2]
  par(cex=sav$cex, font=1)
  text(sav$ilab.xpos[3], 0, "100.0")
  text(sav$xlim[1], y_low, pos=4, bquote(paste("Test for heterogeneity: ",
                                            I^2, "=", .(fmtx(tmp_meta_model$I2, digits=2)), ", ",
                                            tau^2, "=", .(fmtx(tmp_meta_model$tau2, digits=2)), ", ",
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
make_Venn_digram_list = function(named_list, plot_full_path = NULL, label_text_size=4, ...){
  
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
    geom_sf_text(aes(label = name), data = venn_setlabel(data), size = label_text_size, ...) +
    
    # 4. region label layer
    geom_sf_label(aes(label = full_lable), data = venn_region(data), size = label_text_size, ...) +
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
    
    if (stri_detect_fixed(plot_full_path, ".pdf")){
      pdf(plot_full_path, width = 10, height = 10)
      print(plot)
      dev.off()
    }
    
    if (stri_detect_fixed(plot_full_path, ".png")){
      png(plot_full_path, width = 1024, height = 1024)
      print(plot)
      dev.off()
    }
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
# Other datasets have no readily available infromation on suicide or already used by us

################### Importing DE summary stats ###################
summary_paths_no_covar = list.files(path = "Data_preprocessing_analysis", recursive = TRUE, full.names = TRUE)
summary_paths_no_covar = summary_paths_no_covar[stri_detect_fixed(summary_paths_no_covar, pattern = "no_covar")]
summary_paths_no_covar = summary_paths_no_covar[!stri_detect_fixed(summary_paths_no_covar, pattern = "_SV")]
summary_paths_no_covar = summary_paths_no_covar[!stri_detect_fixed(summary_paths_no_covar, pattern = "_cell_types")]

summary_paths_with_covar = list.files(path = "Data_preprocessing_analysis", recursive = TRUE, full.names = TRUE)
summary_paths_with_covar = summary_paths_with_covar[stri_detect_fixed(summary_paths_with_covar, pattern = "Top_table")]
summary_paths_with_covar = summary_paths_with_covar[!stri_detect_fixed(summary_paths_with_covar, pattern = "_SV")]
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

# SV datasets
summary_paths_SV = list.files(path = "Data_preprocessing_analysis", recursive = TRUE, full.names = TRUE)
summary_paths_SV = summary_paths_SV[stri_detect_fixed(summary_paths_SV, pattern = "Top_table")]
summary_paths_SV = summary_paths_SV[stri_detect_fixed(summary_paths_SV, pattern = "_SV")]

data_list_SV = lapply(summary_paths_SV, function(x){
  study_name = unlist(stri_split_fixed(x, pattern = "/"))
  study_name = study_name[length(study_name)]
  study_name = stri_replace_all_fixed(study_name, pattern = "_Top_table_SV.csv", replacement = "")
  
  x = smart_fread(x)
  x$Study = study_name
  return(x)
})
reference_colnames = colnames(data_list_SV[[1]])
data_list_SV = lapply(data_list_SV, function(x){
  
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

gene_names_from_sv = lapply(data_list_SV, function(x){
  
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
gene_names_from_sv = unlist(gene_names_from_sv)
gene_names_from_sv = unique(gene_names_from_sv)

gene_names_full = c(gene_names_no_covar, gene_names_with_covar, gene_names_from_sv)
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

# Harmonizing symbols for SV dataset
combined_df_SV = do.call(rbind, data_list_SV)
dim(combined_df_no_covar)[1] == dim(combined_df_SV)[1] # TRUE -> All genes are already covered in harmoniz_genes_check_dict
combined_df_SV$Corrected_symbol = mapply(function(x,y){
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
}, combined_df_SV$Gene_symbol, combined_df_SV$Gene_symbol_non_hgnc)

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
length(unique(combined_df_SV$Study))

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


cohort_stats_vector_SV = vector()
for (i in 1:length(studies)){
  
  separator_large = "\n\n"
  if (i == 1){
    
    string_1 = paste0("Total associations using surrogate varaibles (SVs): ", nrow(combined_df_SV))
    string_2 = paste0("Total significant associations (nominal) with SVs: ", nrow(combined_df_SV[combined_df_SV$P.Value<0.05,]))
    string_3 = paste0("Average absolute effect size (log2FC) with SVs: ", mean(abs(combined_df_SV$logFC)))
    string_4 = paste0("Average absolute effect size (log2FC) (nominal significance) with SVs: ", mean(abs(combined_df_SV[combined_df_SV$P.Value<0.05,"logFC"])))
    cohort_stats_vector_SV = c(string_1, string_2, string_3,  string_4, separator_large)
  }
  
  TMP_current_cohort = studies[i]
  TMP_current_dataset = combined_df_SV[combined_df_SV$Study == TMP_current_cohort,]
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
  
  cohort_stats_vector_SV = c(cohort_stats_vector_SV, main_header, char_aggregated_stats, separator_large)
  
  if (i == length(studies)){
    cohort_stats_vector_SV = c(cohort_stats_vector_SV, separator_large)
  }
}
writeLines(text = cohort_stats_vector_SV, con = "individual_DE_stats_SV.txt")

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


################### RRA (all brain, no covar, test) ###################
# This is a non-paramentic approach to aggregate ranked lists
combined_df_no_covar_RRA_all = combined_df_no_covar
combined_df_no_covar_RRA_all = combined_df_no_covar_RRA_all[combined_df_no_covar_RRA_all$Tissue_type == "Brain",]
table(combined_df_no_covar_RRA_all$Tissue, combined_df_no_covar_RRA_all$Study)
# GSE102556 many tissues
# GSE202537 many tissues
# GSE66937 many tissues

combined_df_no_covar_RRA_all_1 = combined_df_no_covar_RRA_all[combined_df_no_covar_RRA_all$Study %!in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_no_covar_RRA_all_2 = combined_df_no_covar_RRA_all[combined_df_no_covar_RRA_all$Study %in% c("GSE102556", "GSE202537", "GSE66937"), ]
# Caudate https://pubmed.ncbi.nlm.nih.gov/39164232/
# Nac https://pubmed.ncbi.nlm.nih.gov/38894648/ https://pubmed.ncbi.nlm.nih.gov/38965529/ https://pubmed.ncbi.nlm.nih.gov/39167467/
# GSE102556 -> Dorsolateral prefrontal cortex (dlPFC; BA8/9)
# GSE202537 -> Nac 22652  (more probes analyzed)
# GSE66937 -> prefrontal cortex
combined_df_no_covar_RRA_all_2 = combined_df_no_covar_RRA_all_2[combined_df_no_covar_RRA_all_2$Tissue %in% c("Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                                             "Nac",
                                                                                                             "prefrontal cortex"),]
table(combined_df_no_covar_RRA_all_2$Tissue)
"
Dorsolateral prefrontal cortex (dlPFC; BA8/9)                                           Nac                             prefrontal cortex 
                                        21365                                         22652                                         70523 
"
combined_df_no_covar_RRA_all = rbind(combined_df_no_covar_RRA_all_1, combined_df_no_covar_RRA_all_2)
table(combined_df_no_covar_RRA_all$Tissue, combined_df_no_covar_RRA_all$Study)
combined_df_no_covar_RRA_all_list = lapply(unique(combined_df_no_covar_RRA_all$Study), function(x){
  df = combined_df_no_covar_RRA_all[combined_df_no_covar_RRA_all$Study == x,]
  return(df)
})

combined_df_no_covar_RRA_all_list_reduced = lapply(combined_df_no_covar_RRA_all_list, function(x){
  
  orig_df = x
  
  x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
  x = x[!is.na(x$Corrected_symbol),]
  x = x[x$Corrected_symbol != "",]
  x = x %>% 
    dplyr::group_by(., Corrected_symbol)%>% 
    dplyr::summarise(., 
                     avgLog2FC = mean(logFC), 
                     maxP = max(P.Value)) # Most conservative
  
  x = dplyr::ungroup(x)
  x$Tissue = unique(orig_df$Tissue)
  x$Tissue_type = unique(orig_df$Tissue_type)
  x$Technology = unique(orig_df$Technology)
  x$Study = unique(orig_df$Study)
  x = as.data.frame(x)
  x = dplyr::arrange(x, maxP)
  return(x)
})
combined_df_no_covar_RRA_all_list_reduced = do.call(rbind, combined_df_no_covar_RRA_all_list_reduced)
total_number = length(unique(combined_df_no_covar_RRA_all_list_reduced$Corrected_symbol))

combined_df_no_covar_RRA_all_list_reduced_positive = combined_df_no_covar_RRA_all_list_reduced[combined_df_no_covar_RRA_all_list_reduced$avgLog2FC > 0,]
combined_df_no_covar_RRA_all_list_reduced_negative = combined_df_no_covar_RRA_all_list_reduced[combined_df_no_covar_RRA_all_list_reduced$avgLog2FC < 0,]


list_sizes = combined_df_no_covar_RRA_all_list_reduced_positive %>%
  group_by(., Study) %>%
  summarize(N=length(Corrected_symbol)) %>%
  ungroup(.)
mean(list_sizes$N)

list_sizes = combined_df_no_covar_RRA_all_list_reduced_negative %>%
  group_by(., Study) %>%
  summarize(N=length(Corrected_symbol)) %>%
  ungroup(.)
mean(list_sizes$N)

studies = unique(combined_df_no_covar_RRA_all_list_reduced_positive$Study)
no_covar_RRA_gene_lists_positive = lapply(studies, function(x){
  tmp_df = combined_df_no_covar_RRA_all_list_reduced_positive[combined_df_no_covar_RRA_all_list_reduced_positive$Study == x,]
  tmp_df = dplyr::arrange(tmp_df, maxP)
  ranked_list = tmp_df$Corrected_symbol
  return(ranked_list)
})
no_covar_RRA_gene_lists_positive_top = lapply(studies, function(x){
  tmp_df = combined_df_no_covar_RRA_all_list_reduced_positive[combined_df_no_covar_RRA_all_list_reduced_positive$Study == x,]
  tmp_df = dplyr::arrange(tmp_df, maxP)
  tmp_df = tmp_df[tmp_df$maxP<0.05,]
  ranked_list = tmp_df$Corrected_symbol
  return(ranked_list)
})

studies = unique(combined_df_no_covar_RRA_all_list_reduced_negative$Study)
no_covar_RRA_gene_lists_negative = lapply(studies, function(x){
  tmp_df = combined_df_no_covar_RRA_all_list_reduced_negative[combined_df_no_covar_RRA_all_list_reduced_negative$Study == x,]
  tmp_df = dplyr::arrange(tmp_df, maxP)
  ranked_list = tmp_df$Corrected_symbol
  return(ranked_list)
})
no_covar_RRA_gene_lists_negative_top = lapply(studies, function(x){
  tmp_df = combined_df_no_covar_RRA_all_list_reduced_negative[combined_df_no_covar_RRA_all_list_reduced_negative$Study == x,]
  tmp_df = dplyr::arrange(tmp_df, maxP)
  tmp_df = tmp_df[tmp_df$maxP<0.05,]
  ranked_list = tmp_df$Corrected_symbol
  return(ranked_list)
})

RRA_no_covar_all_brain_positive = aggregateRanks(glist = no_covar_RRA_gene_lists_positive, full = TRUE)
RRA_no_covar_all_brain_positive$P.adjusted = p.adjust(RRA_no_covar_all_brain_positive$Score, method = "fdr")
RRA_no_covar_all_brain_positive$N.cohorts = sapply(RRA_no_covar_all_brain_positive$Name, function(x){
  internal_count = sapply(no_covar_RRA_gene_lists_positive, function(sublist) x %in% sublist)
  internal_count = unlist(internal_count)
  internal_count = as.numeric(internal_count)
  internal_count = sum(internal_count)
  return(internal_count)
})

RRA_no_covar_all_brain_negative = aggregateRanks(glist = no_covar_RRA_gene_lists_negative, full = TRUE)
RRA_no_covar_all_brain_negative$P.adjusted = p.adjust(RRA_no_covar_all_brain_negative$Score, method = "fdr")
RRA_no_covar_all_brain_negative$N.cohorts = sapply(RRA_no_covar_all_brain_negative$Name, function(x){
  internal_count = sapply(no_covar_RRA_gene_lists_positive, function(sublist) x %in% sublist)
  internal_count = unlist(internal_count)
  internal_count = as.numeric(internal_count)
  internal_count = sum(internal_count)
  return(internal_count)
})

# Subset interesting genes
candidates = c("P2RY12", "AC159540.14", "CX3CR1", "DEPP1", "FOS", "GPR34", "KLF4", "PMP2", "SOX9")
RRA_no_covar_all_brain_positive[RRA_no_covar_all_brain_positive$Name %in% candidates,]
RRA_no_covar_all_brain_negative[RRA_no_covar_all_brain_negative$Name %in% candidates,]

# ---> Add RRA in the meta-calculation as columns

add_RRA_to_meta_df = function(meta_df, initial_study_df_list){
  
  # Reduce df for RRA (need p-values)
  print("Reducing df")
  initial_study_df_list_reduced = lapply(initial_study_df_list, function(x){
    orig_df = x
    x = x[!stri_detect_fixed(x$Corrected_symbol, pattern = ";"),]
    x = x[!is.na(x$Corrected_symbol),]
    x = x[x$Corrected_symbol != "",]
    x = x %>% 
      dplyr::group_by(., Corrected_symbol)%>% 
      dplyr::summarise(., 
                       avgLog2FC = mean(logFC), 
                       maxP = max(P.Value)) # Most conservative
    x = dplyr::ungroup(x)
    x$Tissue = unique(orig_df$Tissue)
    x$Tissue_type = unique(orig_df$Tissue_type)
    x$Technology = unique(orig_df$Technology)
    x$Study = unique(orig_df$Study)
    x = as.data.frame(x)
    x = dplyr::arrange(x, maxP)
    return(x)
  })
  
  initial_study_df_list_reduced = do.call(rbind, initial_study_df_list_reduced)
  RRA_all_list_reduced_positive = initial_study_df_list_reduced[initial_study_df_list_reduced$avgLog2FC > 0,]
  RRA_all_list_reduced_negative = initial_study_df_list_reduced[initial_study_df_list_reduced$avgLog2FC < 0,]
  
  # Prepare ranked lists
  studies = unique(RRA_all_list_reduced_positive$Study)
  RRA_gene_lists_positive = lapply(studies, function(x){
    tmp_df = RRA_all_list_reduced_positive[RRA_all_list_reduced_positive$Study == x,]
    tmp_df = dplyr::arrange(tmp_df, maxP)
    ranked_list = tmp_df$Corrected_symbol
    return(ranked_list)
  })
  RRA_gene_lists_negative = lapply(studies, function(x){
    tmp_df = RRA_all_list_reduced_negative[RRA_all_list_reduced_negative$Study == x,]
    tmp_df = dplyr::arrange(tmp_df, maxP)
    ranked_list = tmp_df$Corrected_symbol
    return(ranked_list)
  })
  
  # Run for positive log2FC
  print("Running RRA")
  RRA_all_brain_positive = aggregateRanks(glist = RRA_gene_lists_positive, full = TRUE)
  RRA_all_brain_positive$P.adjusted.RRA = p.adjust(RRA_all_brain_positive$Score, method = "fdr")
  RRA_all_brain_positive$N.cohorts.RRA = sapply(RRA_all_brain_positive$Name, function(x){
    internal_count = sapply(RRA_gene_lists_positive, function(sublist) x %in% sublist)
    internal_count = unlist(internal_count)
    internal_count = as.numeric(internal_count)
    internal_count = sum(internal_count)
    return(internal_count)
  })
  colnames(RRA_all_brain_positive) = c("Name_RRA", "P_RRA", "FDR_RRA", "N_cohorts_RRA")
  
  # Run for negative log2FC
  RRA_all_brain_negative = aggregateRanks(glist = RRA_gene_lists_negative, full = TRUE)
  RRA_all_brain_negative$P.adjusted.RRA = p.adjust(RRA_all_brain_negative$Score, method = "fdr")
  RRA_all_brain_negative$N.cohorts.RRA = sapply(RRA_all_brain_negative$Name, function(x){
    internal_count = sapply(RRA_gene_lists_negative, function(sublist) x %in% sublist)
    internal_count = unlist(internal_count)
    internal_count = as.numeric(internal_count)
    internal_count = sum(internal_count)
    return(internal_count)
  })
  colnames(RRA_all_brain_negative) = c("Name_RRA", "P_RRA", "FDR_RRA", "N_cohorts_RRA")
  
  # Attach stats to meta_df
  print("Preparing attachment")
  attachment = mclapply(meta_df$gene, function(x){
    
    selected_subset = meta_df[meta_df$gene == x,]
    selected_subset = selected_subset[1,] # always ensure 1st row
    
    meta_LFc = selected_subset$meta_LFc[1]
    cohorts_up = selected_subset$cohorts_up[1]
    cohorts_down = selected_subset$cohorts_down[1]
    
    
    if (is.na(meta_LFc)){
      # Number of cohorts less than 5
      # Consider cohort counts only
      
      if (cohorts_up >= cohorts_down){
        # Positive analysis
        P_RRA = RRA_all_brain_positive[RRA_all_brain_positive$Name_RRA == x, "P_RRA"]
        FDR_RRA = RRA_all_brain_positive[RRA_all_brain_positive$Name_RRA == x, "FDR_RRA"]
        N_cohorts_RRA = RRA_all_brain_positive[RRA_all_brain_positive$Name_RRA == x, "N_cohorts_RRA"]
      } else {
        # Negative analysis
        P_RRA = RRA_all_brain_negative[RRA_all_brain_negative$Name_RRA == x, "P_RRA"]
        FDR_RRA = RRA_all_brain_negative[RRA_all_brain_negative$Name_RRA == x, "FDR_RRA"]
        N_cohorts_RRA = RRA_all_brain_negative[RRA_all_brain_negative$Name_RRA == x, "N_cohorts_RRA"]
      }
    } else {
      
      if (meta_LFc >= 0){
        # Positive analysis
        P_RRA = RRA_all_brain_positive[RRA_all_brain_positive$Name_RRA == x, "P_RRA"]
        FDR_RRA = RRA_all_brain_positive[RRA_all_brain_positive$Name_RRA == x, "FDR_RRA"]
        N_cohorts_RRA = RRA_all_brain_positive[RRA_all_brain_positive$Name_RRA == x, "N_cohorts_RRA"]
      } else {
        # Negative analysis
        P_RRA = RRA_all_brain_negative[RRA_all_brain_negative$Name_RRA == x, "P_RRA"]
        FDR_RRA = RRA_all_brain_negative[RRA_all_brain_negative$Name_RRA == x, "FDR_RRA"]
        N_cohorts_RRA = RRA_all_brain_negative[RRA_all_brain_negative$Name_RRA == x, "N_cohorts_RRA"]
      }
    }
    
    att_df_row = data.frame(
      P_RRA = P_RRA,
      FDR_RRA = FDR_RRA,
      N_cohorts_RRA = N_cohorts_RRA
    )
    
    return(att_df_row)
  }, mc.cores = 9)
  
  attachment = do.call(rbind, attachment)
  meta_df = cbind(meta_df, attachment)
  return(meta_df)
}

################### Meta-functions ###################

run_meta_SJ = function(initial_study_list, blood_df_to_compare){
  
  
  # Reducing df for meta-analysis
  print("Aggregating Log2FC and SE per gene per study")
  
  initial_study_list_reduced = lapply(initial_study_list, function(x){
    
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
  
  initial_study_list_reduced = do.call(rbind, initial_study_list_reduced)
  unique_genes_meta = unique(initial_study_list_reduced$Corrected_symbol)
  
  print("Total genes: ")
  print(length(unique_genes_meta))
  
  print("Running meta-analysis per gene")
  
  meta_analysis_list = mclapply(unique_genes_meta, function(x){
    
    
    # meta df
    tmp_df_genes = initial_study_list_reduced[initial_study_list_reduced$Corrected_symbol == x,]
    tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
    tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
    cohorts_total = nrow(tmp_df_genes)
    
    # blood df
    tmp_df_blood = blood_df_to_compare[blood_df_to_compare$Corrected_symbol == x, ]
    
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
                             method = "SJ", 
                             weighted = TRUE)
    
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
  
  meta_analysis_df = do.call(rbind, meta_analysis_list)
  rownames(meta_analysis_df) = NULL
  rm(list = ls(pattern = "tmp_"))
  meta_analysis_df = dplyr::arrange(meta_analysis_df, meta_pval)
  meta_analysis_df$meta_FDR = p.adjust(meta_analysis_df$meta_pval, method = "fdr")
  meta_analysis_df = meta_analysis_df[,c("gene",
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
  
  out_list = list()
  out_list[[1]] = meta_analysis_df
  out_list[[2]] = initial_study_list_reduced
  names(out_list) = c("meta_analysis_df", "initial_study_list_reduced")
  
  return(out_list)
}

run_single_gene_meta_SJ = function(gene,
                                   initial_study_list){
  
  
  # Reducing df for meta-analysis
  print("Aggregating Log2FC and SE per gene per study")
  
  initial_study_list_reduced = lapply(initial_study_list, function(x){
    
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
  
  initial_study_list_reduced = do.call(rbind, initial_study_list_reduced)
  unique_genes_meta = unique(initial_study_list_reduced$Corrected_symbol)
  
  print("Running meta-analysis for gene")
    
    
  # meta df
  tmp_df_genes = initial_study_list_reduced[initial_study_list_reduced$Corrected_symbol == gene,]
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  
  if (cohorts_total < 5){
    print("N cohorts < 5")
  } else {
    
    tmp_meta_model = rma.uni(yi = avgLog2FC, 
                             vi = maxSE^2, 
                             data = tmp_df_genes, 
                             method = "SJ", 
                             weighted = TRUE)
    
    print(summary(tmp_meta_model))
  }
}

run_meta_REML = function(initial_study_list, blood_df_to_compare){
  
  
  # Reducing df for meta-analysis
  print("Aggregating Log2FC and SE per gene per study")
  
  initial_study_list_reduced = lapply(initial_study_list, function(x){
    
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
  
  initial_study_list_reduced = do.call(rbind, initial_study_list_reduced)
  unique_genes_meta = unique(initial_study_list_reduced$Corrected_symbol)
  
  print("Running meta-analysis per gene")
  
  meta_analysis_list = mclapply(unique_genes_meta, function(x){
    
    
    # meta df
    tmp_df_genes = initial_study_list_reduced[initial_study_list_reduced$Corrected_symbol == x,]
    tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
    tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
    cohorts_total = nrow(tmp_df_genes)
    
    # blood df
    tmp_df_blood = blood_df_to_compare[blood_df_to_compare$Corrected_symbol == x, ]
    
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
  
  meta_analysis_df = do.call(rbind, meta_analysis_list)
  rownames(meta_analysis_df) = NULL
  rm(list = ls(pattern = "tmp_"))
  meta_analysis_df = dplyr::arrange(meta_analysis_df, meta_pval)
  meta_analysis_df$meta_FDR = p.adjust(meta_analysis_df$meta_pval, method = "fdr")
  meta_analysis_df = meta_analysis_df[,c("gene",
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
  
  out_list = list()
  out_list[[1]] = meta_analysis_df
  out_list[[2]] = initial_study_list_reduced
  names(out_list) = c("meta_analysis_df", "initial_study_list_reduced")
  
  return(out_list)
}

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

# Run meta-analysis
analysis_container_no_covar_all_brain = run_meta_SJ(initial_study_list = combined_df_no_covar_meta_full_list,
                                                    blood_df_to_compare = blood_df_no_covar)
meta_no_covar_all_brain = analysis_container_no_covar_all_brain[["meta_analysis_df"]]
combined_df_no_covar_meta_full_list_reduced = analysis_container_no_covar_all_brain[["initial_study_list_reduced"]]

meta_no_covar_all_brain = add_RRA_to_meta_df(meta_no_covar_all_brain, combined_df_no_covar_meta_full_list)
non_na_substet = meta_no_covar_all_brain[!is.na(meta_no_covar_all_brain$meta_LFc), ]
cor.test(non_na_substet$meta_pval, non_na_substet$P_RRA, method = "spearman") # 0.1461325
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

# Heterogeneity summary (28644 were excluded from meta due to Ncohorts <5)
summary(meta_no_covar_all_brain$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.059  34.431  49.216  48.750  63.528  98.581   28644 
"

summary(meta_no_covar_all_brain_signif$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
 0.8334 20.4189 32.9392 34.0485 45.3801 78.5286
"

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

# Run meta-analysis
analysis_container_with_covar_all_brain = run_meta_SJ(initial_study_list = combined_df_with_covar_meta_full_list,
                                                    blood_df_to_compare = blood_df_with_covar)



meta_with_covar_all_brain = analysis_container_with_covar_all_brain[["meta_analysis_df"]]
combined_df_with_covar_meta_full_list_reduced = analysis_container_with_covar_all_brain[["initial_study_list_reduced"]]

meta_with_covar_all_brain = add_RRA_to_meta_df(meta_with_covar_all_brain, combined_df_with_covar_meta_full_list)
meta_with_covar_all_brain_signif = meta_with_covar_all_brain[meta_with_covar_all_brain$meta_pval < 0.05,]
meta_with_covar_all_brain_signif = meta_with_covar_all_brain_signif[!is.na(meta_with_covar_all_brain_signif$meta_pval),]
nrow(meta_with_covar_all_brain_signif)

plot_forest_meta_gene(gene_name = "RB1", meta_df = combined_df_with_covar_meta_full_list_reduced)
plot_forest_meta_gene(gene_name = "P2RY12", meta_df = combined_df_with_covar_meta_full_list_reduced)

# Overlap full analysis
overlapping_genes_meta_full = intersect(meta_no_covar_all_brain_signif$gene, meta_with_covar_all_brain_signif$gene)
overlapping_genes_meta_full = unique(overlapping_genes_meta_full) 
length(overlapping_genes_meta_full) # 57 genes overlap
overlapping_subset_allbrain_ncv = meta_no_covar_all_brain_signif[meta_no_covar_all_brain_signif$gene %in% overlapping_genes_meta_full, ]
overlapping_subset_allbrain_ncv = dplyr::arrange(overlapping_subset_allbrain_ncv, gene)
overlapping_subset_allbrain_cv = meta_with_covar_all_brain_signif[meta_with_covar_all_brain_signif$gene %in% overlapping_genes_meta_full, ]
overlapping_subset_allbrain_cv = dplyr::arrange(overlapping_subset_allbrain_cv, gene)

overlapping_subset_allbrain_stat = cbind(
  overlapping_subset_allbrain_ncv$gene,
  overlapping_subset_allbrain_ncv$meta_LFc,
  overlapping_subset_allbrain_cv$gene,
  overlapping_subset_allbrain_cv$meta_LFc
)
overlapping_subset_allbrain_stat = as.data.frame(overlapping_subset_allbrain_stat)
overlapping_subset_allbrain_stat$V2 = as.numeric(overlapping_subset_allbrain_stat$V2)
overlapping_subset_allbrain_stat$V4 = as.numeric(overlapping_subset_allbrain_stat$V4)
overlapping_subset_allbrain_stat$mean = mapply(function(x,y) mean(x,y),overlapping_subset_allbrain_stat$V2, overlapping_subset_allbrain_stat$V4)

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

# Heterogeneity estimates
summary(meta_with_covar_all_brain$I2)

"   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.001  25.237  39.893  41.247  56.093  99.471   29942 
"
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

# Run meta-analysis
analysis_container_no_covar_cortical = run_meta_SJ(initial_study_list = combined_df_no_covar_meta_cortical_list,
                                                      blood_df_to_compare = blood_df_no_covar)



meta_no_covar_cortical = analysis_container_no_covar_cortical[["meta_analysis_df"]]
combined_df_no_covar_meta_cortical_list_reduced = analysis_container_no_covar_cortical[["initial_study_list_reduced"]]


meta_no_covar_cortical = add_RRA_to_meta_df(meta_no_covar_cortical, combined_df_no_covar_meta_cortical_list)
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

# Heterogeneity summary (31325 were excluded from meta due to Ncohorts <5)
summary(meta_no_covar_cortical$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.021  31.841  48.919  48.347  65.182  98.987   31325 
"
summary(meta_no_covar_cortical_signif$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.215  22.009  34.357  35.388  48.552  83.960 
"

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

# Run meta-analysis
analysis_container_with_covar_cortical = run_meta_SJ(initial_study_list = combined_df_with_covar_meta_cortical_list,
                                                   blood_df_to_compare = blood_df_with_covar)



meta_with_covar_cortical = analysis_container_with_covar_cortical[["meta_analysis_df"]]
combined_df_with_covar_meta_cortical_list_reduced = analysis_container_with_covar_cortical[["initial_study_list_reduced"]]


meta_with_covar_cortical = add_RRA_to_meta_df(meta_with_covar_cortical, combined_df_with_covar_meta_cortical_list)
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

# Heterogeneity estimates
summary(meta_with_covar_cortical$I2)

"     Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.03   22.37   37.90   39.75   55.92   99.57   33360 
"


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

# Run meta-analysis
analysis_container_no_covar_prefrontal = run_meta_SJ(initial_study_list = combined_df_no_covar_meta_prefrontal_list,
                                                     blood_df_to_compare = blood_df_no_covar)



meta_no_covar_prefrontal = analysis_container_no_covar_prefrontal[["meta_analysis_df"]]
combined_df_no_covar_meta_prefrontal_list_reduced = analysis_container_no_covar_prefrontal[["initial_study_list_reduced"]]


meta_no_covar_prefrontal = add_RRA_to_meta_df(meta_no_covar_prefrontal, combined_df_no_covar_meta_prefrontal_list)
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

# Heterogeneity summary (33247 were excluded from meta due to Ncohorts <5)
summary(meta_no_covar_prefrontal$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00   30.27   48.05   47.93   65.83   99.22   33247 
"

summary(meta_no_covar_prefrontal_signif$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   0.28   21.44   35.29   37.44   50.67   90.73 
"

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

# Run meta-analysis
analysis_container_with_covar_prefrontal = run_meta_SJ(initial_study_list = combined_df_with_covar_meta_prefrontal_list,
                                                     blood_df_to_compare = blood_df_with_covar)



meta_with_covar_prefrontal = analysis_container_with_covar_prefrontal[["meta_analysis_df"]]
combined_df_with_covar_meta_prefrontal_list_reduced = analysis_container_with_covar_prefrontal[["initial_study_list_reduced"]]


meta_with_covar_prefrontal = add_RRA_to_meta_df(meta_with_covar_prefrontal, combined_df_with_covar_meta_prefrontal_list)
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

# Heterogeneity estimates
summary(meta_with_covar_prefrontal$I2)

"       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.03   19.87   35.90   38.37   55.02   99.65   35343 
"

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

# Run meta-analysis
REML_analysis_container_no_covar_all_brain = run_meta_REML(initial_study_list = combined_df_no_covar_meta_full_list,
                                                    blood_df_to_compare = blood_df_no_covar)
REML_meta_no_covar_all_brain = REML_analysis_container_no_covar_all_brain[["meta_analysis_df"]]

REML_meta_no_covar_all_brain = add_RRA_to_meta_df(REML_meta_no_covar_all_brain, combined_df_no_covar_meta_full_list)
REML_meta_no_covar_all_brain_signif = REML_meta_no_covar_all_brain[REML_meta_no_covar_all_brain$meta_pval < 0.05,]
REML_meta_no_covar_all_brain_signif = REML_meta_no_covar_all_brain_signif[!is.na(REML_meta_no_covar_all_brain_signif$meta_pval),]
nrow(REML_meta_no_covar_all_brain_signif)

# Heterogeneity
summary(REML_meta_no_covar_all_brain$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00    0.00   12.71   20.95   39.31   94.73   28644 
"

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

# Run meta-analysis
REML_analysis_container_with_covar_all_brain = run_meta_REML(initial_study_list = combined_df_with_covar_meta_full_list,
                                                           blood_df_to_compare = blood_df_with_covar)
REML_meta_with_covar_all_brain = REML_analysis_container_with_covar_all_brain[["meta_analysis_df"]]

REML_meta_with_covar_all_brain = add_RRA_to_meta_df(REML_meta_with_covar_all_brain, combined_df_with_covar_meta_full_list)
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

# Run meta-analysis
REML_analysis_container_no_covar_cortical = run_meta_REML(initial_study_list = combined_df_no_covar_meta_cortical_list,
                                                             blood_df_to_compare = blood_df_no_covar)
REML_meta_no_covar_cortical = REML_analysis_container_no_covar_cortical[["meta_analysis_df"]]

REML_meta_no_covar_cortical = add_RRA_to_meta_df(REML_meta_no_covar_cortical, combined_df_no_covar_meta_cortical_list)
REML_meta_no_covar_cortical_signif = REML_meta_no_covar_cortical[REML_meta_no_covar_cortical$meta_pval < 0.05,]
REML_meta_no_covar_cortical_signif = REML_meta_no_covar_cortical_signif[!is.na(REML_meta_no_covar_cortical_signif$meta_pval),]
nrow(REML_meta_no_covar_cortical_signif)

# Heterogeneity
summary(REML_meta_no_covar_cortical$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00    0.00   12.86   22.62   42.70   95.93   31325 
"

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

# Run meta-analysis
REML_analysis_container_with_covar_cortical = run_meta_REML(initial_study_list = combined_df_with_covar_meta_cortical_list,
                                                          blood_df_to_compare = blood_df_with_covar)
REML_meta_with_covar_cortical = REML_analysis_container_with_covar_cortical[["meta_analysis_df"]]

REML_meta_with_covar_cortical = add_RRA_to_meta_df(REML_meta_with_covar_cortical, combined_df_with_covar_meta_cortical_list)
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

# Run meta-analysis
REML_analysis_container_no_covar_prefrontal = run_meta_REML(initial_study_list = combined_df_no_covar_meta_prefrontal_list,
                                                            blood_df_to_compare = blood_df_no_covar)
REML_meta_no_covar_prefrontal = REML_analysis_container_no_covar_prefrontal[["meta_analysis_df"]]

REML_meta_no_covar_prefrontal = add_RRA_to_meta_df(REML_meta_no_covar_prefrontal, combined_df_no_covar_meta_prefrontal_list)
REML_meta_no_covar_prefrontal_signif = REML_meta_no_covar_prefrontal[REML_meta_no_covar_prefrontal$meta_pval < 0.05,]
REML_meta_no_covar_prefrontal_signif = REML_meta_no_covar_prefrontal_signif[!is.na(REML_meta_no_covar_prefrontal_signif$meta_pval),]
nrow(REML_meta_no_covar_prefrontal_signif)

# Heterogeneity
summary(REML_meta_no_covar_prefrontal$I2)
"
  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.00    0.00   14.48   24.13   45.66   97.75   33247 
"


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

# Run meta-analysis
REML_analysis_container_with_covar_prefrontal = run_meta_REML(initial_study_list = combined_df_with_covar_meta_prefrontal_list,
                                                            blood_df_to_compare = blood_df_with_covar)
REML_meta_with_covar_prefrontal = REML_analysis_container_with_covar_prefrontal[["meta_analysis_df"]]

REML_meta_with_covar_prefrontal = add_RRA_to_meta_df(REML_meta_with_covar_prefrontal, combined_df_with_covar_meta_prefrontal_list)
REML_meta_with_covar_prefrontal_signif = REML_meta_with_covar_prefrontal[REML_meta_with_covar_prefrontal$meta_pval < 0.05,]
REML_meta_with_covar_prefrontal_signif = REML_meta_with_covar_prefrontal_signif[!is.na(REML_meta_with_covar_prefrontal_signif$meta_pval),]
nrow(REML_meta_with_covar_prefrontal_signif)

################### Saving REML sensitivity ###################
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

################### SV meta (all brain) ###################
combined_df_SV_meta_all = combined_df_SV
combined_df_SV_meta_all = combined_df_SV_meta_all[combined_df_SV_meta_all$Tissue_type == "Brain",]
blood_df_SV = combined_df_SV[combined_df_SV$Tissue_type == "Blood",]

table(combined_df_SV_meta_all$Tissue, combined_df_SV_meta_all$Study)
# GSE102556 many tissues
# GSE202537 many tissues
# GSE66937 many tissues

combined_df_SV_meta_all_1 = combined_df_SV_meta_all[combined_df_SV_meta_all$Study %!in% c("GSE102556", "GSE202537", "GSE66937"), ]
combined_df_SV_meta_all_2 = combined_df_SV_meta_all[combined_df_SV_meta_all$Study %in% c("GSE102556", "GSE202537", "GSE66937"), ]
# Caudate https://pubmed.ncbi.nlm.nih.gov/39164232/
# Nac https://pubmed.ncbi.nlm.nih.gov/38894648/ https://pubmed.ncbi.nlm.nih.gov/38965529/ https://pubmed.ncbi.nlm.nih.gov/39167467/
# GSE102556 -> Dorsolateral prefrontal cortex (dlPFC; BA8/9)
# GSE202537 -> Nac 22652  (more probes analyzed)
# GSE66937 -> prefrontal cortex
combined_df_SV_meta_all_2 = combined_df_SV_meta_all_2[combined_df_SV_meta_all_2$Tissue %in% c("Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                                                                "Nac",
                                                                                                                "prefrontal cortex"),]
table(combined_df_SV_meta_all_2$Tissue)
"
Dorsolateral prefrontal cortex (dlPFC; BA8/9)                                           Nac                             prefrontal cortex 
                                        21365                                         22652                                         70523 
"
combined_df_SV_meta_full = rbind(combined_df_SV_meta_all_1, combined_df_SV_meta_all_2)
table(combined_df_SV_meta_full$Tissue, combined_df_SV_meta_full$Study)
combined_df_SV_meta_full_list = lapply(unique(combined_df_SV_meta_full$Study), function(x){
  df = combined_df_SV_meta_full[combined_df_SV_meta_full$Study == x,]
  return(df)
})

# Run meta-analysis
analysis_container_SV_all_brain = run_meta_SJ(initial_study_list = combined_df_SV_meta_full_list,
                                                    blood_df_to_compare = blood_df_SV)
meta_SV_all_brain = analysis_container_SV_all_brain[["meta_analysis_df"]]
combined_df_SV_meta_full_list_reduced = analysis_container_SV_all_brain[["initial_study_list_reduced"]]

meta_SV_all_brain = add_RRA_to_meta_df(meta_SV_all_brain, combined_df_SV_meta_full_list)
meta_SV_all_brain_signif = meta_SV_all_brain[meta_SV_all_brain$meta_pval < 0.05,]
meta_SV_all_brain_signif = meta_SV_all_brain_signif[!is.na(meta_SV_all_brain_signif$meta_pval),]
meta_SV_all_brain_signif = arrange(meta_SV_all_brain_signif, -meta_LFc)
nrow(meta_SV_all_brain_signif)

# Heterogeneity
summary(meta_SV_all_brain$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.031  25.926  40.445  41.201  55.704  98.730   28644 
"

# plots 
dir.create("forest_plots")
folder = "forest_plots/meta_SV_all_brain"
dir.create(folder)

for (i in 1:nrow(meta_SV_all_brain_signif)){
  
  
  gene = meta_SV_all_brain$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_SV_meta_full_list_reduced))
  dev.off()
  
}

################### SV meta (cortical) ###################
combined_df_SV_meta_cortical = combined_df_SV
combined_df_SV_meta_cortical = combined_df_SV_meta_cortical[combined_df_SV_meta_cortical$Tissue_type == "Brain",]
blood_df_SV = combined_df_SV[combined_df_SV$Tissue_type == "Blood",]

table(combined_df_SV_meta_cortical$Tissue, combined_df_SV_meta_cortical$Study)

combined_df_SV_meta_cortical = combined_df_SV_meta_cortical[combined_df_SV_meta_cortical$Tissue %in% 
                                                                          c("DLPFC",
                                                                            "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                            "Dorsolateral Prefrontal Cortex",
                                                                            "Dorsolateral prefrontal cortex (BA9)",
                                                                            "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                            "Orbitofrontal cortex",
                                                                            "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                            "prefrontal cortex",
                                                                            "temporal cortex (BA20 and BA36)"),]
table(combined_df_SV_meta_cortical$Tissue, combined_df_SV_meta_cortical$Study)

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

combined_df_SV_meta_cortical_list = lapply(unique(combined_df_SV_meta_cortical$Study), function(x){
  df = combined_df_SV_meta_cortical[combined_df_SV_meta_cortical$Study == x,]
  return(df)
})

# Run meta-analysis
analysis_container_SV_cortical = run_meta_SJ(initial_study_list = combined_df_SV_meta_cortical_list,
                                              blood_df_to_compare = blood_df_SV)
meta_SV_cortical = analysis_container_SV_cortical[["meta_analysis_df"]]
combined_df_SV_meta_cortical_list_reduced = analysis_container_SV_cortical[["initial_study_list_reduced"]]

meta_SV_cortical = add_RRA_to_meta_df(meta_SV_cortical, combined_df_SV_meta_cortical_list)
meta_SV_cortical_signif = meta_SV_cortical[meta_SV_cortical$meta_pval < 0.05,]
meta_SV_cortical_signif = meta_SV_cortical_signif[!is.na(meta_SV_cortical_signif$meta_pval),]
meta_SV_cortical_signif = arrange(meta_SV_cortical_signif, -meta_LFc)
nrow(meta_SV_cortical_signif)

summary(meta_SV_cortical$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
  0.033  22.011  37.720  39.168  54.900  97.410   31325 
"

dir.create("forest_plots")
folder = "forest_plots/meta_SV_cortical"
dir.create(folder)

for (i in 1:nrow(meta_SV_cortical_signif)){
  
  
  gene = meta_SV_cortical$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_SV_meta_cortical_list_reduced))
  dev.off()
  
}


################### SV meta (prefrontal) ###################
combined_df_SV_meta_prefrontal = combined_df_SV
combined_df_SV_meta_prefrontal = combined_df_SV_meta_prefrontal[combined_df_SV_meta_prefrontal$Tissue_type == "Brain",]
blood_df_SV = combined_df_SV[combined_df_SV$Tissue_type == "Blood",]

table(combined_df_SV_meta_prefrontal$Tissue, combined_df_SV_meta_prefrontal$Study)

combined_df_SV_meta_prefrontal = combined_df_SV_meta_prefrontal[combined_df_SV_meta_prefrontal$Tissue %in% 
                                                                c("DLPFC",
                                                                  "dorsal lateral prefrontal cortex (Brodmann Area 9)",
                                                                  "Dorsolateral Prefrontal Cortex",
                                                                  "Dorsolateral prefrontal cortex (BA9)",
                                                                  "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
                                                                  "postmortem dorsolateral prefrontal cortex (DLPFC)",
                                                                  "prefrontal cortex"),]
table(combined_df_SV_meta_prefrontal$Tissue, combined_df_SV_meta_prefrontal$Study)
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


combined_df_SV_meta_prefrontal_list = lapply(unique(combined_df_SV_meta_prefrontal$Study), function(x){
  df = combined_df_SV_meta_prefrontal[combined_df_SV_meta_prefrontal$Study == x,]
  return(df)
})

# Run meta-analysis
analysis_container_SV_prefrontal = run_meta_SJ(initial_study_list = combined_df_SV_meta_prefrontal_list,
                                             blood_df_to_compare = blood_df_SV)
meta_SV_prefrontal = analysis_container_SV_prefrontal[["meta_analysis_df"]]
combined_df_SV_meta_prefrontal_list_reduced = analysis_container_SV_prefrontal[["initial_study_list_reduced"]]

meta_SV_prefrontal = add_RRA_to_meta_df(meta_SV_prefrontal, combined_df_SV_meta_prefrontal_list)
meta_SV_prefrontal_signif = meta_SV_prefrontal[meta_SV_prefrontal$meta_pval < 0.05,]
meta_SV_prefrontal_signif = meta_SV_prefrontal_signif[!is.na(meta_SV_prefrontal_signif$meta_pval),]
meta_SV_prefrontal_signif = arrange(meta_SV_prefrontal_signif, -meta_LFc)
nrow(meta_SV_prefrontal_signif)

summary(meta_SV_prefrontal$I2)
"
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
   0.02   18.40   33.37   35.54   50.55   98.13   33247 
"

dir.create("forest_plots")
folder = "forest_plots/meta_SV_prefrontal"
dir.create(folder)

for (i in 1:nrow(meta_SV_prefrontal_signif)){
  
  
  gene = meta_SV_prefrontal$gene[i]
  plot_path =  paste0(folder, "/", gene, "_forest.pdf")
  
  pdf(file = plot_path, width = 10, height = 7)
  print(plot_forest_meta_gene(gene_name = gene, meta_df = combined_df_SV_meta_prefrontal_list_reduced))
  dev.off()
  
}

################### Saving SV sensitivity ###################

# SV numbers
# GSE208338 -> 1 suggested
# GSE5388 -> 5 suggested and 5 used
# GSE5389 -> 0 suggested
# GSE66937 -> DLPFC 0 suggested
# GSE199536 -> 5 suggested and 5 used
# GSE92538_U133A -> 3 suggested
# GSE92538_U133_PLUS2 -> 2 suggested
# GSE102556 -> very large number -> 5 used
# GSE243356 -> very large number (59) -> 5 used
# GSE248260 -> very large number (22) -> 5 used
# GSE247998 -> very large number (98) -> 5 used
# GSE202537 -> very large number (>60) -> 5 used
# GSE101521 -> very large number (57) -> 5 used
# GSE144136 -> very large number (30) -> 5 used
# GSE213982 -> very large number (34) -> 5 used

significant_genes_var = ls(pattern = "meta_SV")
significant_genes_var = significant_genes_var[stri_detect_regex(significant_genes_var, pattern = "signif")]
significant_genes_var = significant_genes_var[!stri_detect_regex(significant_genes_var, pattern = "meta_list")]
significant_genes_counts = sapply(significant_genes_var, function(x) nrow(get(x)))

stats_meta = paste0("Signifi gene counts for: ", significant_genes_var)
stats_meta = paste0(stats_meta, " = ", significant_genes_counts)
writeLines(stats_meta, "stats_meta_SV.txt")

# saving files
meta_list_SV_signif = list(
  "meta_SV_all_brain_signif" = meta_SV_all_brain_signif,
  "meta_SV_cortical_signif" = meta_SV_cortical_signif,
  "meta_SV_prefrontal_signif" = meta_SV_prefrontal_signif
)

names_for_list = c(
  "Meta all brain (SVs)",
  "Meta cortical (SVs)",
  "Meta prefrontal cortex (SVs)"
)

wb = createWorkbook()

for (i in 1:3) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, names_for_list[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = names_for_list[i], meta_list_SV_signif[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = names_for_list[i], cols = 1:ncol(meta_list_SV_signif[[i]]), widths = "auto")
}

saveWorkbook(wb, "Meta_suicide_significant_genes_SV.xlsx", overwrite = TRUE)


## Save volcano plots
# making volcanos in a loop
joined_list_all_metas_SV = list(
  meta_SV_all_brain,
  meta_SV_cortical,
  meta_SV_prefrontal
)
names_for_list_all = c("All brain (SVs)", "Cortical (SVs)", "DLPFC (SVs)")

for (i in 1:3){
  
  curr_df = joined_list_all_metas_SV[[i]]
  plot_path = paste0("volcano_plots/", names_for_list_all[i], ".png")
  plot_meta_volcano(curr_df, figure_path = plot_path, title = names_for_list_all[i])
  
}

################### Venn diagrams for SVs ###################

meta_list_SV_signif = list(
  meta_SV_all_brain_signif,
  meta_SV_cortical_signif,
  meta_SV_prefrontal_signif
)

meta_gene_list_overlaps_SV = lapply(meta_list_SV_signif, function(x){
  x = x$gene
  return(x)
})

names(meta_gene_list_overlaps_SV) = c(
  "All brain",
  "Cortical regions",
  "DLPFC"
)

make_Venn_digram_list(named_list = meta_gene_list_overlaps_SV, palette = 9, plot_full_path = "Venn_SV.pdf")


dir.create("analyses_comparison_covar_SV_nothing_venn")
meta_gene_list_SVs = list(
  meta_no_covar_all_brain_signif,
  meta_SV_all_brain_signif,
  meta_with_covar_all_brain_signif
)
meta_gene_list_SVs = lapply(meta_gene_list_SVs, function(x){
  x = x$gene
  return(x)
})

names(meta_gene_list_SVs) = c(
  "All brain (no covar)",
  "All brain (SVs)",
  "All brain (covar)"
)

make_Venn_digram_list(named_list = meta_gene_list_SVs,
                      palette = 5,
                      label_text_size = 7,
                      plot_full_path = "analyses_comparison_covar_SV_nothing_venn/1_Venn_SVs_vs_other_all.png")


meta_gene_list_SVs = list(
  meta_no_covar_cortical_signif,
  meta_SV_cortical_signif,
  meta_with_covar_cortical_signif
)
meta_gene_list_SVs = lapply(meta_gene_list_SVs, function(x){
  x = x$gene
  return(x)
})

names(meta_gene_list_SVs) = c(
  "Cortical (no covar)",
  "Cortical (SVs)",
  "Cortical (covar)"
)

make_Venn_digram_list(named_list = meta_gene_list_SVs, 
                      label_text_size = 7,
                      palette = 4, plot_full_path = "analyses_comparison_covar_SV_nothing_venn/2_Venn_SVs_vs_other_cortical.png")

meta_gene_list_SVs = list(
  meta_no_covar_prefrontal_signif,
  meta_SV_prefrontal_signif,
  meta_with_covar_prefrontal_signif
)
meta_gene_list_SVs = lapply(meta_gene_list_SVs, function(x){
  x = x$gene
  return(x)
})

names(meta_gene_list_SVs) = c(
  "Prefrontal (no covar)",
  "Prefrontal (SVs)",
  "Prefrontal (covar)"
)

make_Venn_digram_list(named_list = meta_gene_list_SVs, 
                      label_text_size = 7,
                      palette = 3, plot_full_path = "analyses_comparison_covar_SV_nothing_venn/3_Venn_SVs_vs_other_prefrontal.png")


# Define cropping function
crop_image <- function(img, row_range, col_range) {
  # img: 3D array (height, width, channels)
  img[row_range, col_range, , drop = FALSE]
}

images_in_folder = list.files("analyses_comparison_covar_SV_nothing_venn", full.names = TRUE)
images_in_folder = images_in_folder[!grepl("combined_image", images_in_folder)]

image_list  =  lapply(images_in_folder, png::readPNG)

# Define cropping region
row_range = 50:950
col_range = 50:950

cropped_list = lapply(image_list, crop_image, row_range, col_range)
image_grobs = lapply(cropped_list, rasterGrob)

# Combine
height = length(row_range)
width = length(col_range)
combined_file_name = "analyses_comparison_covar_SV_nothing_venn/Venn_DE_comparison_combined_image.png"

png(filename = combined_file_name, width = width*3, height = height, units = "px")
grid.arrange(grobs = image_grobs[1:3], ncol = 3)
dev.off()

################### Moderator analysis for RE with SJ ###################
cohort_moderators_data = openxlsx::read.xlsx("Cohort_moderators.xlsx")
colnames(cohort_moderators_data)[1] = "Study"
cohort_moderators_data = cohort_moderators_data[cohort_moderators_data$Study != "GSE66937", ]

cohort_moderators_data$Platform_binary = factor(cohort_moderators_data$Platform_binary, levels = c( "Array", "RNAseq"))
cohort_moderators_data$Primary_tissue_group = factor(cohort_moderators_data$Primary_tissue_group, 
                                                     levels = c( "Cortical", "Non_cortical")) # Better to stratify


# Moderators to include
cohort_moderators_data$Platform_binary
cohort_moderators_data$Mean_PMI_hours
cohort_moderators_data$Percent_male
cohort_moderators_data$Percent_Depr
cohort_moderators_data$Percent_BD
cohort_moderators_data$Percent_SCZ

# Check
all(cohort_moderators_data$Study %in% combined_df_no_covar$Study) # TRUE

run_meta_SJ_moderated = function(initial_study_list, 
                                 moderator_df,
                                 moderator_names){
  
  
  # Reducing df for meta-analysis
  print("Aggregating Log2FC and SE per gene per study")
  initial_study_list_reduced = lapply(initial_study_list, function(x){
    
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
  initial_study_list_reduced = do.call(rbind, initial_study_list_reduced)
  unique_genes_meta = unique(initial_study_list_reduced$Corrected_symbol)
  
  print("Running meta-analysis per gene with moderators")
  meta_analysis_list = mclapply(unique_genes_meta, function(x){
    
    
    # meta df
    tmp_df_genes = initial_study_list_reduced[initial_study_list_reduced$Corrected_symbol == x,]
    
    # Adding moderators
    tmp_df_genes = inner_join(tmp_df_genes, moderator_df)
    tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
    tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
    cohorts_total = nrow(tmp_df_genes)
    
    if (cohorts_total < 10){
      return(NA)
    }
    
    moderator_string = paste0(moderator_names, collapse = " + ")
    moderator_string = paste0("~ ", moderator_string)
    
    
    tmp_meta_model = rma.uni(yi = avgLog2FC, 
                             vi = maxSE^2, 
                             data = tmp_df_genes, 
                             method = "SJ",
                             mods = as.formula(moderator_string),
                             weighted = TRUE)
    
    tmp_meta_model_estimates = tmp_meta_model$b
    tmp_meta_model_estimates = t(tmp_meta_model_estimates)
    tmp_meta_model_estimates = as.data.frame(tmp_meta_model_estimates)
    
    tmp_meta_model_pvals = tmp_meta_model$pval
    tmp_meta_model_pvals = as.data.frame(tmp_meta_model_pvals)
    tmp_meta_model_pvals = t(tmp_meta_model_pvals)
    
    colnames(tmp_meta_model_pvals) = paste0("pval_", colnames(tmp_meta_model_estimates))
    colnames(tmp_meta_model_estimates) = paste0("beta_", colnames(tmp_meta_model_estimates))
    
    
    leading_row = data.frame(
      gene = x,
      detected_cohorts = cohorts_total,
      cohorts_up = tmp_cohorts_up,
      cohorts_down = tmp_cohorts_down,
      commentary = "Number of cohorts >= 10: gene is analyzed")
    
    estimate_row = cbind(tmp_meta_model_estimates, tmp_meta_model_pvals)
    full_row = cbind(leading_row, estimate_row)
    
    full_row$tau2 = tmp_meta_model$tau2
    full_row$I2 = tmp_meta_model$I2
    full_row$H2 = tmp_meta_model$H2
    full_row$Q = tmp_meta_model$QE
    full_row$Q.p = tmp_meta_model$QEp
    full_row$QM = tmp_meta_model$QM
    full_row$QM.p = tmp_meta_model$QMp
    full_row$R2 =  tmp_meta_model$R2
    
    return(full_row)
    
  }, mc.cores = 9)
  
  #print(meta_analysis_list)
  
  meta_analysis_list = meta_analysis_list[sapply(meta_analysis_list, is.data.frame)]
  meta_analysis_df = do.call(rbind, meta_analysis_list)
  rownames(meta_analysis_df) = NULL
  rm(list = ls(pattern = "tmp_"))
  out_list = list()
  out_list[[1]] = meta_analysis_df
  out_list[[2]] = initial_study_list_reduced
  names(out_list) = c("meta_analysis_df", "initial_study_list_reduced")
  
  return(out_list)
}


run_single_gene_meta_SJ_moderated = function(gene,
                                      initial_study_list, 
                                      moderator_df,
                                      moderator_names){
  
  
  # Reducing df for meta-analysis
  print("Aggregating Log2FC and SE per gene per study")
  initial_study_list_reduced = lapply(initial_study_list, function(x){
    
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
  initial_study_list_reduced = do.call(rbind, initial_study_list_reduced)
  unique_genes_meta = unique(initial_study_list_reduced$Corrected_symbol)
  
  print("Running meta-analysis for gene with moderators")
    
  # meta df
  tmp_df_genes = initial_study_list_reduced[initial_study_list_reduced$Corrected_symbol == gene,]
  
  # Adding moderators
  tmp_df_genes = inner_join(tmp_df_genes, moderator_df)
  tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
  tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
  cohorts_total = nrow(tmp_df_genes)
  
  if (cohorts_total < 10){
    return(NA)
  }
  
  moderator_string = paste0(moderator_names, collapse = " + ")
  moderator_string = paste0("~ ", moderator_string)
  
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, 
                           vi = maxSE^2, 
                           data = tmp_df_genes, 
                           method = "SJ",
                           mods = as.formula(moderator_string),
                           weighted = TRUE)
  
  print(summary(tmp_meta_model))
  
}


# Moderated all brain (no covar)
MODERATED_analysis_all_brain_no_covar = run_meta_SJ_moderated(
  initial_study_list = combined_df_no_covar_meta_full_list,
  moderator_df = cohort_moderators_data,
  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                      "Percent_male", "Percent_Depr", 
                      "Percent_BD", "Percent_SCZ")
)
MODERATED_analysis_all_brain_no_covar = MODERATED_analysis_all_brain_no_covar[[1]]

# Moderated cortical (no covar)
MODERATED_analysis_cortical_no_covar = run_meta_SJ_moderated(
  initial_study_list = combined_df_no_covar_meta_cortical_list,
  moderator_df = cohort_moderators_data,
  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                      "Percent_male", "Percent_Depr", 
                      "Percent_BD", "Percent_SCZ")
)
MODERATED_analysis_cortical_no_covar = MODERATED_analysis_cortical_no_covar[[1]]


# Moderated all brain (with covar)
MODERATED_analysis_all_brain_with_covar = run_meta_SJ_moderated(
  initial_study_list = combined_df_with_covar_meta_full_list,
  moderator_df = cohort_moderators_data,
  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                      "Percent_male", "Percent_Depr", 
                      "Percent_BD", "Percent_SCZ")
)
MODERATED_analysis_all_brain_with_covar = MODERATED_analysis_all_brain_with_covar[[1]]

# Moderated cortical (with covar)
MODERATED_analysis_cortical_with_covar = run_meta_SJ_moderated(
  initial_study_list = combined_df_with_covar_meta_cortical_list,
  moderator_df = cohort_moderators_data,
  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                      "Percent_male", "Percent_Depr", 
                      "Percent_BD", "Percent_SCZ")
)
MODERATED_analysis_cortical_with_covar = MODERATED_analysis_cortical_with_covar[[1]]


# Moderated all brain (SV)
MODERATED_analysis_all_brain_SV = run_meta_SJ_moderated(
  initial_study_list = combined_df_SV_meta_full_list,
  moderator_df = cohort_moderators_data,
  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                      "Percent_male", "Percent_Depr", 
                      "Percent_BD", "Percent_SCZ")
)
MODERATED_analysis_all_brain_SV = MODERATED_analysis_all_brain_SV[[1]]

# Moderated cortical (SV)
MODERATED_analysis_cortical_SV = run_meta_SJ_moderated(
  initial_study_list = combined_df_SV_meta_cortical_list,
  moderator_df = cohort_moderators_data,
  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                      "Percent_male", "Percent_Depr", 
                      "Percent_BD", "Percent_SCZ")
)
MODERATED_analysis_cortical_SV = MODERATED_analysis_cortical_SV[[1]]

# Cortical correlations
cor.test(MODERATED_analysis_cortical_no_covar$beta_Percent_BD, 
         MODERATED_analysis_cortical_no_covar$beta_Percent_Depr,
         method = "spearman") # 0.8972178; p-value < 2.2e-16

cor.test(MODERATED_analysis_cortical_no_covar$pval_Percent_BD, 
         MODERATED_analysis_cortical_no_covar$pval_Percent_Depr,
         method = "spearman") # 0.7594179; p-value < 2.2e-16

cor.test(MODERATED_analysis_cortical_no_covar$pval_Percent_BD, 
         MODERATED_analysis_cortical_no_covar$pval_Percent_SCZ,
         method = "spearman") # 0.9225474; p-value < 2.2e-16

cor.test(MODERATED_analysis_cortical_with_covar$pval_Percent_BD, 
         MODERATED_analysis_cortical_with_covar$pval_Percent_Depr,
         method = "spearman") # 0.9211632; p-value < 2.2e-16

cor.test(MODERATED_analysis_cortical_with_covar$pval_Percent_BD, 
         MODERATED_analysis_cortical_with_covar$pval_Percent_SCZ,
         method = "spearman") # 0.9226585; p-value < 2.2e-16

cor.test(MODERATED_analysis_cortical_SV$pval_Percent_BD, 
         MODERATED_analysis_cortical_SV$pval_Percent_Depr,
         method = "spearman") # 0.803596; p-value < 2.2e-16

# All brain correlations
cor.test(MODERATED_analysis_all_brain_no_covar$beta_Percent_BD, 
         MODERATED_analysis_all_brain_no_covar$beta_Percent_Depr,
         method = "spearman") # 0.6962317 ; p-value < 2.2e-16

cor.test(MODERATED_analysis_all_brain_no_covar$pval_Percent_BD, 
         MODERATED_analysis_all_brain_no_covar$pval_Percent_Depr,
         method = "spearman") # 0.6735521 ; p-value < 2.2e-16

cor.test(MODERATED_analysis_all_brain_no_covar$pval_Percent_BD, 
         MODERATED_analysis_all_brain_no_covar$pval_Percent_SCZ,
         method = "spearman") # 0.6448914; p-value < 2.2e-16

cor.test(MODERATED_analysis_all_brain_with_covar$pval_Percent_BD, 
         MODERATED_analysis_all_brain_with_covar$pval_Percent_Depr,
         method = "spearman") # 0.6958896; p-value < 2.2e-16

cor.test(MODERATED_analysis_all_brain_with_covar$pval_Percent_BD, 
         MODERATED_analysis_all_brain_with_covar$pval_Percent_SCZ,
         method = "spearman") # 0.574108; p-value < 2.2e-16

cor.test(MODERATED_analysis_all_brain_SV$pval_Percent_BD, 
         MODERATED_analysis_all_brain_SV$pval_Percent_Depr,
         method = "spearman") # 0.6541848; p-value < 2.2e-16

# Plots
dir.create("moderator_correlations")
plot_moderator_correlation = function(moderated_meta_df,
                                      variable_1,
                                      variable_2,
                                      title,
                                      label_1, 
                                      label_2,
                                      path_to_save,
                                      color_dots){
  
  
  
  # Rename variables
  var_num_1 = which(colnames(moderated_meta_df) == variable_1)
  var_num_2 = which(colnames(moderated_meta_df) == variable_2)
  
  colnames(moderated_meta_df)[var_num_1] = "Corr_variable_1"
  colnames(moderated_meta_df)[var_num_2] = "Corr_variable_2"
  
  
  # Run correlations
  n_genes = nrow(moderated_meta_df)
  cor_spearman = cor.test(moderated_meta_df$Corr_variable_1, moderated_meta_df$Corr_variable_2, method = "spearman")$estimate
  cor_pearson = cor.test(moderated_meta_df$Corr_variable_1, moderated_meta_df$Corr_variable_2, method = "pearson")$estimate
  cor_spearman = round(cor_spearman, digits = 2)
  cor_pearson = round(cor_pearson, digits = 2)
  
  
  caption = paste0("Genes: ", n_genes, "\n",
                   "R (Pearson): ", cor_pearson, "\n",
                   "R (Spearman): ", cor_spearman)
  
  # Generate plot
  meta_cor_plot = ggplot(data = moderated_meta_df, aes(x = Corr_variable_1, y = Corr_variable_2)) +
    geom_point(alpha = .5, col = color_dots, size = 1.5) +
    labs(title=title,
         x = label_1, 
         y = label_2,
         caption=caption) +
    
    # Custom the theme:
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_line(size = 0.1, linetype = 2, color =  "black"),
      panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
      axis.line = element_line(size = 1, linetype = 1, color =  "black"),
      panel.background = element_blank(),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 12),
    )
  ggsave(file = path_to_save, plot = meta_cor_plot, width=3000, height=2560, units = "px", scale = 1)
}

plot_moderator_correlation(moderated_meta_df=MODERATED_analysis_all_brain_no_covar,
                           variable_1 = "beta_Percent_Depr",
                           variable_2 = "beta_Percent_BD",
                           title = expression("Correlations of " * beta ~ " % Depression and " * beta ~ " % BD in all brain"),
                           label_1 =  expression(beta ~ " % depression"),
                           label_2 =  expression(beta ~ " % BD"),
                           path_to_save = "moderator_correlations/depr_vs_BD_all_brain.png",
                           color_dots = "darkgreen")

plot_moderator_correlation(moderated_meta_df=MODERATED_analysis_cortical_no_covar,
                           variable_1 = "beta_Percent_Depr",
                           variable_2 = "beta_Percent_BD",
                           title = expression("Correlations of " * beta ~ " % Depression and " * beta ~ " % BD in cortical"),
                           label_1 =  expression(beta ~ " % depression"),
                           label_2 =  expression(beta ~ " % BD"),
                           path_to_save = "moderator_correlations/depr_vs_BD_cortical.png",
                           color_dots = "#AA4A44")

plot_moderator_correlation(moderated_meta_df=MODERATED_analysis_cortical_no_covar,
                           variable_1 = "beta_Percent_Depr",
                           variable_2 = "beta_Percent_SCZ",
                           title = expression("Correlations of " * beta ~ " % Depression and " * beta ~ " % SCZ in cortical"),
                           label_1 =  expression(beta ~ " % depression"),
                           label_2 =  expression(beta ~ " % SCZ"),
                           path_to_save = "moderator_correlations/depr_vs_SCZ_cortical.png",
                           color_dots = "#00008B")


#### Generate combined image for moderators
images_in_folder = list.files("moderator_correlations", full.names = TRUE)
image_list  =  lapply(images_in_folder, png::readPNG)
image_grobs = lapply(image_list, grid::rasterGrob)
combined_file_name = paste0("Figure_S_Moderator_correlations.png")
height = nrow(image_list[[1]])
width = ncol(image_list[[1]])

spacer = rectGrob(gp=gpar(col=NA, fill=NA))
image_grobs_modif = list(
  
  image_grobs[[1]], image_grobs[[2]], image_grobs[[3]]
)


# generating PNG
png(filename = combined_file_name, units = "px", width = width*3, height = height*1)
grid.arrange(grobs = image_grobs[1:3], ncol = 3)
dev.off()


# Generate moderator tables for selected genes
list_of_moderated_analyses = list(
  MODERATED_analysis_all_brain_no_covar,
  MODERATED_analysis_cortical_no_covar,
  MODERATED_analysis_all_brain_with_covar,
  MODERATED_analysis_cortical_with_covar,
  MODERATED_analysis_all_brain_SV,
  MODERATED_analysis_cortical_SV
)

list_of_significant_genes_from_ma = list(
  meta_no_covar_all_brain_signif,
  meta_no_covar_cortical_signif,
  meta_with_covar_all_brain_signif,
  meta_with_covar_cortical_signif,
  meta_SV_all_brain_signif,
  meta_SV_cortical_signif
)

filtered_moderator_results = mapply(function(mod_df, meta_df){
  
  meta_df = dplyr::arrange(meta_df, -meta_LFc)
  meta_df = meta_df[meta_df$gene %in% mod_df$gene, ]
  
  rownames(mod_df) = mod_df$gene
  mod_df = mod_df[meta_df$gene, ]
  
  att_1 = meta_df[, c("gene", "meta_LFc", "meta_pval")]
  colnames(att_1) = c("gene", "init_meta_logFC", "init_meta_pval")
  att_2 = mod_df
  att_2$gene = NULL
  
  comb_df = cbind(att_1, att_2)
  
  return(comb_df)
  
}, list_of_moderated_analyses, list_of_significant_genes_from_ma, SIMPLIFY = FALSE)

names(filtered_moderator_results) = c(
  "Mod. for all brain (no covar) p (meta) <0.05",
  "Mod. for cortical (no covar) p (meta) <0.05",
  "Mod. for all brain (with covar) p (meta) <0.05",
  "Mod. for cortical (with covar) p (meta) <0.05",
  "Mod. for all brain (SV) p (meta) <0.05",
  "Mod. for cortical (SV) p (meta) <0.05"
)


run_single_gene_meta_SJ_moderated(gene="P2RY12",
                                  initial_study_list = combined_df_no_covar_meta_full_list,
                                  moderator_df = cohort_moderators_data,
                                  moderator_names = c("Platform_binary", "Mean_PMI_hours",
                                               "Percent_male", "Percent_Depr", 
                                               "Percent_BD", "Percent_SCZ"))
run_single_gene_meta_SJ(gene="P2RY12",
                        initial_study_list = combined_df_no_covar_meta_full_list)

# Run explanations for models
explain_moderators = function(moder_df, df_name){
  
  genes_one_signif = moder_df %>%
    filter(., pval_Platform_binaryRNAseq < 0.05 | pval_Mean_PMI_hours < 0.05 | pval_Percent_male < 0.05 | pval_Percent_Depr < 0.05 | pval_Percent_BD < 0.05 | pval_Percent_SCZ < 0.05) %>%
    pull(gene) %>%
    length(.) %>%
    paste0("N genes where at least one moderator was significant (p<0.05):", .)
  
  
  genes_signif_platform = moder_df %>%
    filter(., pval_Platform_binaryRNAseq < 0.05) %>%
    arrange(., pval_Platform_binaryRNAseq) %>%
    mutate(gene_with_effect = paste0(gene, " (", round(beta_Platform_binaryRNAseq, 2), ")")) %>%
    pull(gene_with_effect) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) genes where effect size was related to platform (p<0.05):", .)
  
  genes_signif_PMI = moder_df %>%
    filter(., pval_Mean_PMI_hours < 0.05) %>%
    arrange(., pval_Mean_PMI_hours) %>%
    mutate(gene_with_effect = paste0(gene, " (", round(beta_Mean_PMI_hours, 2), ")")) %>%
    pull(gene_with_effect) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) genes where effect size was related to mean PMI (p<0.05):", .)
  
  genes_signif_percent_male = moder_df %>%
    filter(., pval_Percent_male < 0.05) %>%
    arrange(., pval_Percent_male) %>%
    mutate(gene_with_effect = paste0(gene, " (", round(beta_Percent_male, 2), ")")) %>%
    pull(gene_with_effect) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) genes where effect size was related to % of male (p<0.05):", .) 
  
  genes_signif_percent_depr = moder_df %>%
    filter(., pval_Percent_Depr < 0.05) %>%
    arrange(., pval_Percent_Depr) %>%
    mutate(gene_with_effect = paste0(gene, " (", round(beta_Percent_Depr, 2), ")")) %>%
    pull(gene_with_effect) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) genes where effect size was related to % of depression (p<0.05):", .) 
  
  genes_signif_percent_BD = moder_df %>%
    filter(., pval_Percent_BD < 0.05) %>%
    arrange(., pval_Percent_BD) %>%
    mutate(gene_with_effect = paste0(gene, " (", round(beta_Percent_BD, 2), ")")) %>%
    pull(gene_with_effect) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) genes where effect size was related to % of BD (p<0.05):", .) 
  
  genes_signif_percent_SCZ = moder_df %>%
    filter(., pval_Percent_SCZ < 0.05) %>%
    arrange(., pval_Percent_SCZ) %>%
    mutate(gene_with_effect = paste0(gene, " (", round(beta_Percent_SCZ, 2), ")")) %>%
    pull(gene_with_effect) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) genes where effect size was related to % of SCZ (p<0.05):", .) 
  
  average_R2 = moder_df %>%
    pull(R2) %>%
    mean(.) %>%
    paste0("Mean R2 for analysis subset: ", .) 
  
  separator_small = "\n\n"
  separator_large = "\n\n\n\n"
  
  
  output = c(df_name,
             separator_small,
             genes_one_signif,
             separator_small,
             genes_signif_platform,  
             separator_small,
             genes_signif_PMI,
             separator_small,
             genes_signif_percent_male,
             separator_small,
             genes_signif_percent_depr,
             separator_small,
             genes_signif_percent_BD,
             separator_small,
             genes_signif_percent_SCZ,
             separator_small,
             average_R2,
             separator_large)
  output = paste0(output, collapse = "")
  return(output)
}

stats_moderator_analysis = vector()

for (i in 1:length(filtered_moderator_results)){
  message = explain_moderators(filtered_moderator_results[[i]], names(filtered_moderator_results)[i])
  stats_moderator_analysis = c(stats_moderator_analysis, message)
}

writeLines(text = stats_moderator_analysis, con = "Stats_moderator_analysis.txt")

# Save moderators as Excel
short_names = c(
  "Mod. all (no cov) p(met)<0.05",
  "Mod. cort (no cov) p(met)<0.05",
  "Mod. all (cov) p(met)<0.05",
  "Mod. cort (cov) p(met)<0.05",
  "Mod. all (SV) p(met)<0.05",
  "Mod. cortical (SV) p(met)<0.05"
)

wb = createWorkbook()
for (i in 1:6) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, short_names[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = short_names[i], filtered_moderator_results[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = short_names[i], cols = 1:ncol(filtered_moderator_results[[i]]), widths = "auto")
}

saveWorkbook(wb, "Moderator_analysis_results.xlsx", overwrite = TRUE)

################### Renaming function for meta-analyses ###################

rename_function = function(x){
  template = ""
  
  if (stri_detect_fixed(x, "all_brain")){
    template = paste0(template, "All brain ")
  } else if (stri_detect_fixed(x, "cortical")){
    template = paste0(template, "Cortical ")
  } else if (stri_detect_fixed(x, "prefrontal")){
    template = paste0(template, "Prefrontal ")
  }
  
  if (stri_detect_fixed(x, "no_covar")){
    template = paste0(template, "(no covariates, ")
  } else if (stri_detect_fixed(x, "with_covar")){
    template = paste0(template, "(covariates,  ")
  } else if (stri_detect_fixed(x, "SV")){
    template = paste0(template, "(SVs,  ")
  }
  
  if (stri_detect_fixed(x, "REML")){
    template = paste0(template, "REML, ")
  } else {
    template = paste0(template, "SJ,  ")
  }
  
  if (stri_detect_fixed(x, "MODERATED")){
    template = paste0(template, "moderators, ")
  }
  
  if (stri_detect_fixed(x, "signif")){
    template = paste0(template, "p<0.05)")
  } else {
    template = paste0(template, "all probes)")
  }
  return(template)
}

################### Other Sensitivity analyses? ###################

table(combined_df_no_covar$Tissue, combined_df_no_covar$Study)

# Amygdala vs hippocampus vs thalamus
# Anterior Insula (aINS) vs Cingulate gyrus 25 (Cg25) vs  Nucleus Accumbens (Nac) vs Subiculum (Sub)

# Striatum: caudate nucleus and the putamen; The nucleus accumbens (NAc) is a major component of the ventral striatum 
# The subiculum (Latin for "support") is the most inferior component of the hippocampal formation
# The striatum sends projections directly to the medial habenula https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3894476/


################### Correlations of Meta-effects ###################
dir.create("meta_effect_correlations")

meta_analysis_lists = ls(pattern = "meta")
meta_analysis_lists = meta_analysis_lists[stri_detect_regex(meta_analysis_lists, pattern = "all_brain")]
meta_analysis_lists = meta_analysis_lists[!stri_detect_regex(meta_analysis_lists, pattern = "signif")]
meta_analysis_lists = meta_analysis_lists[!stri_detect_regex(meta_analysis_lists, pattern = "REML")]

meta_analysis_lists_datasets = lapply(meta_analysis_lists, function(x) get(x))

sapply(meta_analysis_lists_datasets, function(x){
  x = x[!is.na(x$meta_LFc), ]
  rows = nrow(x)
  return(rows)
})

names_for_list = c(
  "Meta all brain (no covar.)",
  "Meta all brain (SVs)",
  "Meta all brain (covar.)"
)

plot_meta_correlation = function(meta_df_1,
                                 meta_df_2, 
                                 label_1, 
                                 label_2,
                                 path_to_save,
                                 color_dots){
  
  # Modify labels
  label_1 = paste0("Log2FC ", label_1)
  label_2 = paste0("Log2FC ", label_2)
  
  # Select genes where effects are calculated
  meta_df_1 = meta_df_1[!is.na(meta_df_1$meta_LFc), ]
  meta_df_2 = meta_df_2[!is.na(meta_df_2$meta_LFc), ]
  
  # Subset by overlapping genes
  overlapping_genes = intersect(meta_df_1$gene, meta_df_2$gene)
  meta_df_1 = meta_df_1[meta_df_1$gene %in% overlapping_genes, ]
  meta_df_2 = meta_df_2[meta_df_2$gene %in% overlapping_genes, ]
  
  # Arrange
  meta_df_1 = dplyr::arrange(meta_df_1, gene)
  meta_df_2 = dplyr::arrange(meta_df_2, gene)
  
  # Check gene match
  print(all(meta_df_1$gene == meta_df_2$gene))
  
  tmp_df = cbind(meta_df_1$meta_LFc, meta_df_2$meta_LFc)
  colnames(tmp_df) = c("lFC_1", "lFC_2")
  
  len_overlapping_genes = length(overlapping_genes)
  cor_spearman = cor.test(meta_df_1$meta_LFc, meta_df_2$meta_LFc, method = "spearman")$estimate
  cor_pearson = cor.test(meta_df_1$meta_LFc, meta_df_2$meta_LFc, method = "pearson")$estimate
  cor_spearman = round(cor_spearman, digits = 2)
  cor_pearson = round(cor_pearson, digits = 2)
  
  
  caption = paste0("Genes in both: ", len_overlapping_genes, "\n",
                   "R (Pearson): ", cor_pearson, "\n",
                   "R (Spearman): ", cor_spearman)
  
  # Generate plot
  meta_cor_plot = ggplot(data = tmp_df, aes(x = lFC_1, y = lFC_2)) +
    geom_point(alpha = .5, col = color_dots, size = 1.5) +
    labs(x = label_1, 
         y = label_2,
         caption=caption) +
    scale_x_continuous(limits = c(-1.5, 1.5)) +
    scale_y_continuous(limits = c(-1.5, 1.5)) +
    # Custom the theme:
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_line(size = 0.1, linetype = 2, color =  "black"),
      panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
      axis.line = element_line(size = 1, linetype = 1, color =  "black"),
      panel.background = element_blank(),
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.caption = element_text(size = 12),
    )
  ggsave(file = path_to_save, plot = meta_cor_plot, width=2560, height=2560, units = "px", scale = 1)
}

plot_meta_correlation(meta_analysis_lists_datasets[[1]],
                      meta_analysis_lists_datasets[[2]],
                      names_for_list[1],
                      names_for_list[2],
                      path_to_save = "meta_effect_correlations/F1.png",
                      color_dots = "darkgreen")

plot_meta_correlation(meta_analysis_lists_datasets[[1]],
                      meta_analysis_lists_datasets[[3]],
                      names_for_list[1],
                      names_for_list[3],
                      path_to_save = "meta_effect_correlations/F2.png",
                      color_dots = "#AA4A44")

plot_meta_correlation(meta_analysis_lists_datasets[[2]],
                      meta_analysis_lists_datasets[[3]],
                      names_for_list[2],
                      names_for_list[3],
                      path_to_save = "meta_effect_correlations/F3.png",
                      color_dots = "#00008B")


# Tissue correlations
meta_analysis_lists_tissue = ls(pattern = "meta")
meta_analysis_lists_tissue = meta_analysis_lists_tissue[stri_detect_regex(meta_analysis_lists_tissue, pattern = "no_covar")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[stri_detect_regex(meta_analysis_lists_tissue, pattern = "signif")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "REML")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "_list")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "combined_df_")]

meta_analysis_lists_datasets = lapply(meta_analysis_lists_tissue, function(x) get(x))
names_for_list = c(
  "Meta all brain (no covar., p<0.05)",
  "Meta cortical (no covar., p<0.05)",
  "Meta prefrontal (no covar., p<0.05)"
)

plot_meta_correlation(meta_analysis_lists_datasets[[1]],
                      meta_analysis_lists_datasets[[2]],
                      names_for_list[1],
                      names_for_list[2],
                      path_to_save = "meta_effect_correlations/F4.png",
                      color_dots = "darkred")

plot_meta_correlation(meta_analysis_lists_datasets[[1]],
                      meta_analysis_lists_datasets[[3]],
                      names_for_list[1],
                      names_for_list[3],
                      path_to_save = "meta_effect_correlations/F5.png",
                      color_dots = "darkred")

plot_meta_correlation(meta_analysis_lists_datasets[[2]],
                      meta_analysis_lists_datasets[[3]],
                      names_for_list[2],
                      names_for_list[3],
                      path_to_save = "meta_effect_correlations/F6.png",
                      color_dots = "darkred")

# Tissue correlations all genes
meta_analysis_lists_tissue = ls(pattern = "meta")
meta_analysis_lists_tissue = meta_analysis_lists_tissue[stri_detect_regex(meta_analysis_lists_tissue, pattern = "no_covar")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "signif")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "REML")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "_list")]
meta_analysis_lists_tissue = meta_analysis_lists_tissue[!stri_detect_regex(meta_analysis_lists_tissue, pattern = "combined_df_")]

meta_analysis_lists_datasets = lapply(meta_analysis_lists_tissue, function(x) get(x))
names_for_list = c(
  "Meta all brain (no covar., all genes)",
  "Meta cortical (no covar., all genes)",
  "Meta prefrontal (no covar., all genes)"
)

plot_meta_correlation(meta_analysis_lists_datasets[[1]],
                      meta_analysis_lists_datasets[[2]],
                      names_for_list[1],
                      names_for_list[2],
                      path_to_save = "meta_effect_correlations/F7.png",
                      color_dots = "darkblue")

plot_meta_correlation(meta_analysis_lists_datasets[[1]],
                      meta_analysis_lists_datasets[[3]],
                      names_for_list[1],
                      names_for_list[3],
                      path_to_save = "meta_effect_correlations/F8.png",
                      color_dots = "darkblue")

plot_meta_correlation(meta_analysis_lists_datasets[[2]],
                      meta_analysis_lists_datasets[[3]],
                      names_for_list[2],
                      names_for_list[3],
                      path_to_save = "meta_effect_correlations/F9.png",
                      color_dots = "darkblue")



images_in_folder = list.files("meta_effect_correlations", full.names = TRUE)
image_list  =  lapply(images_in_folder, png::readPNG)
image_grobs = lapply(image_list, grid::rasterGrob)
combined_file_name = paste0("Figure__Meta_effect_correlations.png")
height = nrow(image_list[[1]])
width = ncol(image_list[[1]])

spacer = rectGrob(gp=gpar(col=NA, fill=NA))
image_grobs_modif = list(
  
  image_grobs[[1]], image_grobs[[2]], image_grobs[[3]],
  spacer, spacer,spacer,
  image_grobs[[4]], image_grobs[[5]], image_grobs[[6]]
)

row_heights = unit(c(1, 3, 1, 3),   # Heights for each row
                   c("null", "inches", "null", "inches"))

# generating PNG
png(filename = combined_file_name, units = "px", width = width*3, height = height*2)
grid::grid.newpage()
gridExtra::grid.arrange(grobs = image_grobs_modif, ncol = 3, heights = row_heights)
dev.off()

#### Generate image for all genes
images_in_folder = list.files("meta_effect_correlations", full.names = TRUE)
images_in_folder = images_in_folder[7:9]
image_list  =  lapply(images_in_folder, png::readPNG)
image_grobs = lapply(image_list, grid::rasterGrob)
combined_file_name = paste0("Figure_S_Meta_effect_correlations_all.png")
height = nrow(image_list[[1]])
width = ncol(image_list[[1]])

spacer = rectGrob(gp=gpar(col=NA, fill=NA))
image_grobs_modif = list(
  
  image_grobs[[1]], image_grobs[[2]], image_grobs[[3]]
)


# generating PNG
png(filename = combined_file_name, units = "px", width = width*3, height = height*1)
grid.arrange(grobs = image_grobs[1:3], ncol = 3)
dev.off()


################### Extended Comparison of Meta-genes ###################

analysis_pairs = list(
  c("meta_no_covar_all_brain", "meta_with_covar_all_brain"),
  c("meta_no_covar_cortical", "meta_with_covar_cortical"),
  c("meta_no_covar_prefrontal", "meta_with_covar_prefrontal"),
  c("meta_no_covar_all_brain", "meta_SV_all_brain"),
  c("meta_no_covar_cortical", "meta_SV_cortical"),
  c("meta_no_covar_prefrontal", "meta_SV_prefrontal"),
  c("meta_with_covar_all_brain", "meta_SV_all_brain"),
  c("meta_with_covar_cortical", "meta_SV_cortical"),
  c("meta_with_covar_prefrontal", "meta_SV_prefrontal")
)

compare_analyses_results = function(pair){
  
  data_full_1 =  get(pair[1])
  data_full_2 =  get(pair[2])
  
  # Remove rows without estimated log2FC
  data_full_1 = data_full_1[!is.na(data_full_1$meta_LFc),]
  data_full_2 = data_full_2[!is.na(data_full_2$meta_LFc),]
  
  name1 = rename_function(pair[1])
  name2 = rename_function(pair[2])
  
  Data_name = paste0("######### Summary of comparison ", name1, " VS ",  name2," #########")
  
  overlapping_genes_meta_all_included = intersect(data_full_1$gene, data_full_2$gene)
  overlapping_genes_meta_all_included = unique(overlapping_genes_meta_all_included) 
  n_genes_overlap_all_included = length(overlapping_genes_meta_all_included)
  n_genes_overlap_all_included = paste0("Total number of shared analyzed genes in both strategies: ", n_genes_overlap_all_included)
  
  overlapping_subset_all_inc_1 = data_full_1[data_full_1$gene %in% overlapping_genes_meta_all_included, ]
  overlapping_subset_all_inc_1 = dplyr::arrange(overlapping_subset_all_inc_1, gene)
  overlapping_subset_all_inc_2 = data_full_2[data_full_2$gene %in% overlapping_genes_meta_all_included, ]
  overlapping_subset_all_inc_2 = dplyr::arrange(overlapping_subset_all_inc_2, gene)
  
  overlapping_subset_all_inc_stat = cbind(
    overlapping_subset_all_inc_1$gene,
    overlapping_subset_all_inc_1$meta_LFc,
    overlapping_subset_all_inc_2$gene,
    overlapping_subset_all_inc_2$meta_LFc
  )
  
  overlapping_subset_all_inc_stat = as.data.frame(overlapping_subset_all_inc_stat)
  overlapping_subset_all_inc_stat$V2 = as.numeric(overlapping_subset_all_inc_stat$V2)
  overlapping_subset_all_inc_stat$V4 = as.numeric(overlapping_subset_all_inc_stat$V4)
  overlapping_subset_all_inc_stat$mean = mapply(function(x,y) mean(x,y),overlapping_subset_all_inc_stat$V2, overlapping_subset_all_inc_stat$V4)
  overlapping_subset_all_inc_stat$dir_match = mapply(function(x,y){
    if (x<0 & y<0){
      return(TRUE)
    }
    if (x>0 & y>0){
      return(TRUE)
    }
    return(FALSE)
    
  },overlapping_subset_all_inc_stat$V2, overlapping_subset_all_inc_stat$V4)
  
  # Correlations
  cor_spearman = cor.test(overlapping_subset_all_inc_stat$V2, overlapping_subset_all_inc_stat$V4, method = "spearman")$estimate
  cor_pearson = cor.test(overlapping_subset_all_inc_stat$V2, overlapping_subset_all_inc_stat$V4, method = "pearson")$estimate
  cor_spearman = round(cor_spearman, digits = 2)
  cor_pearson = round(cor_pearson, digits = 2)
  
  cor_spearman = paste0("Spearman correlation of all log2FC: ", cor_spearman)
  cor_pearson = paste0("Pearson correlation of all log2FC: ", cor_pearson)
  
  
  percent_match_all_inc = table(overlapping_subset_all_inc_stat$dir_match)
  percent_match_all_inc = percent_match_all_inc["TRUE"]/sum(percent_match_all_inc)
  percent_match_all_inc = percent_match_all_inc * 100
  percent_match_all_inc = round(percent_match_all_inc, digits=2)
  percent_match_all_inc = paste0("Percent of genes with matching direction (All analyzed): ", percent_match_all_inc, "%")
  
  n_genes_overla_up_all_inc = overlapping_subset_all_inc_stat %>%
    filter(., dir_match, mean>0) %>%
    pull(V1) %>%
    length() %>%
    paste0("N genes with direction match and log2FC>0 (All analyzed): ", .)
  
  n_genes_overla_down_all_inc = overlapping_subset_all_inc_stat %>%
    filter(., dir_match, mean<0) %>%
    pull(V1) %>%
    length() %>%
    paste0("N genes with direction match and log2FC<0 (All analyzed): ", .)
  
  
  # Compatibility
  dat1 =  data_full_1[data_full_1$meta_pval < 0.05,]
  dat2 =  data_full_2[data_full_2$meta_pval < 0.05,]
  dat1 =  dat1[!is.na(dat1$meta_pval),]
  dat2 =  dat2[!is.na(dat2$meta_pval),]
  
  overlapping_genes_meta_full = intersect(dat1$gene, dat2$gene)
  overlapping_genes_meta_full = unique(overlapping_genes_meta_full) 
  n_genes_overlap = length(overlapping_genes_meta_full)
  n_genes_overlap = paste0("Total number of overlapping genes (p<0.05): ", n_genes_overlap)
  
  overlapping_subset_dat1 = dat1[dat1$gene %in% overlapping_genes_meta_full, ]
  overlapping_subset_dat1 = dplyr::arrange(overlapping_subset_dat1, gene)
  overlapping_subset_dat2 = dat2[dat2$gene %in% overlapping_genes_meta_full, ]
  overlapping_subset_dat2 = dplyr::arrange(overlapping_subset_dat2, gene)
  
  overlapping_subset_stat = cbind(
    overlapping_subset_dat1$gene,
    overlapping_subset_dat1$meta_LFc,
    overlapping_subset_dat2$gene,
    overlapping_subset_dat2$meta_LFc
  )
  
  overlapping_subset_stat = as.data.frame(overlapping_subset_stat)
  overlapping_subset_stat$V2 = as.numeric(overlapping_subset_stat$V2)
  overlapping_subset_stat$V4 = as.numeric(overlapping_subset_stat$V4)
  overlapping_subset_stat$mean = mapply(function(x,y) mean(x,y),overlapping_subset_stat$V2, overlapping_subset_stat$V4)
  overlapping_subset_stat$dir_match = mapply(function(x,y){
    if (x<0 & y<0){
      return(TRUE)
    }
    if (x>0 & y>0){
      return(TRUE)
    }
    return(FALSE)
    
  },overlapping_subset_stat$V2, overlapping_subset_stat$V4)
  
  n_genes_overla_up = overlapping_subset_stat %>%
    filter(., dir_match, mean>0) %>%
    pull(V1) %>%
    length() %>%
    paste0("N genes with direction match and log2FC>0 (p<0.05): ", .)
  
  n_genes_overla_down = overlapping_subset_stat %>%
    filter(., dir_match, mean<0) %>%
    pull(V1) %>%
    length() %>%
    paste0("N genes with direction match and log2FC<0 (p<0.05): ", .)
  
  top_genes_RE_nom_sig_up = overlapping_subset_stat %>%
    filter(., dir_match, mean>=0.1) %>%
    arrange(., -mean) %>%
    pull(V1) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) upregulated genes p<0.05 and log2FC>=0.1 sorted by log2FC:\n", .)
  
  top_genes_RE_nom_sig_down = overlapping_subset_stat %>%
    filter(., dir_match, mean<=-0.1) %>%
    arrange(., mean) %>%
    pull(V1) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) downregulated genes p<0.05 and log2FC<=-0.1 sorted by log2FC:\n", .)
  
  separator_small = "\n"
  separator_large = "\n\n\n\n"
  output = c(Data_name, 
             separator_small,
             n_genes_overlap_all_included,
             separator_small,
             cor_spearman,
             separator_small,
             cor_pearson,
             separator_small,
             percent_match_all_inc,
             separator_small,
             n_genes_overla_up_all_inc,
             separator_small,
             n_genes_overla_down_all_inc,
             separator_small,
             n_genes_overlap,
             separator_small,
             n_genes_overla_up,
             separator_small,
             n_genes_overla_down,
             separator_small,
             top_genes_RE_nom_sig_up,
             separator_small,
             top_genes_RE_nom_sig_down,
             separator_large)
  output = paste0(output, collapse = "")
  return(output)
  
}

stats_analyses_comparisons = vector()

for (i in 1:length(meta_list_signif_combined)){
  
  message = compare_analyses_results(analysis_pairs[[i]])
  stats_analyses_comparisons = c(stats_analyses_comparisons, message)
  
}
writeLines(text = stats_analyses_comparisons, con = "Data_S_stats_analyses_comparisons.txt")

################### Stats on Meta ###################
significant_genes_var = ls(pattern = "meta")
significant_genes_var = significant_genes_var[stri_detect_regex(significant_genes_var, pattern = "signif")]
significant_genes_var = significant_genes_var[!stri_detect_regex(significant_genes_var, pattern = "meta_list")]
#significant_genes_var = significant_genes_var[!stri_detect_regex(significant_genes_var, pattern = "SV")]
#significant_genes_var = significant_genes_var[!stri_detect_regex(significant_genes_var, pattern = "REML")]
significant_genes_counts = sapply(significant_genes_var, function(x) nrow(get(x)))

stats_meta = paste0("Signifi gene counts for: ", significant_genes_var)
stats_meta = paste0(stats_meta, " = ", significant_genes_counts)
writeLines(stats_meta, "stats_meta.txt")

################### Funnel plots for Meta ###################

plot_funnel_meta_gene = function(gene_name, meta_df){
  
  tmp_df_genes = meta_df[meta_df$Corrected_symbol == gene_name,]
  
  if(nrow(tmp_df_genes)<1){
    writeLines("Gene is missing in the data")
    return(NULL)
  }
  
  tmp_meta_model = rma.uni(yi = avgLog2FC, vi = maxSE^2, data = tmp_df_genes, method = "SJ", weighted = TRUE)
  
  # Example Visualization
  weights = fmtx(weights(tmp_meta_model), digits=1)
  
  funnel(tmp_meta_model)
  title(paste0("Gene: ", gene_name))
  
}

plot_funnel_global = function(meta_df){
  funnel(x = meta_df$avgLog2FC, vi = meta_df$maxSE**2,
         cex.axis = 2,   # axis tick labels
         cex.lab = 2)  # axis labels)
  title("Global funnel plot for all genes",cex.main = 3)
}

png(file = "Fig_S_Global_funnel_plot.png", width = 2560, height = 1440, units="px")
print(plot_funnel_global(combined_df_no_covar_meta_full_list_reduced))
dev.off()


# making plots
dir.create("funnel_plots")

directory_names = c(
  "meta_no_covar_all_brain",
  "meta_no_covar_cortical",
  "meta_no_covar_prefrontal",
  "meta_with_covar_all_brain",
  "meta_with_covar_cortical",
  "meta_with_covar_prefrontal",
  "meta_SV_all_brain",
  "meta_SV_cortical",
  "meta_SV_prefrontal"
)

mt_reduced_lists = list(
  combined_df_no_covar_meta_full_list_reduced,
  combined_df_no_covar_meta_cortical_list_reduced,
  combined_df_no_covar_meta_prefrontal_list_reduced,
  combined_df_with_covar_meta_full_list_reduced,
  combined_df_with_covar_meta_cortical_list_reduced,
  combined_df_with_covar_meta_prefrontal_list_reduced,
  combined_df_SV_meta_full_list_reduced,
  combined_df_SV_meta_cortical_list_reduced,
  combined_df_SV_meta_prefrontal_list_reduced
)

mt_signif_lists = list(
  meta_no_covar_all_brain_signif,
  meta_no_covar_cortical_signif,
  meta_no_covar_prefrontal_signif,
  meta_with_covar_all_brain_signif,
  meta_with_covar_cortical_signif,
  meta_with_covar_prefrontal_signif,
  meta_SV_all_brain_signif,
  meta_SV_cortical_signif,
  meta_SV_prefrontal_signif
)

for (i in 1:length(directory_names)){
  
  selected_sig_df = mt_signif_lists[[i]]
  selected_reduced_df = mt_reduced_lists[[i]]
  
  folder = paste0("funnel_plots/", directory_names[i])
  
  dir.create(folder)
  
  for (sub_index in 1:nrow(selected_sig_df)){
    
    
    gene = selected_sig_df$gene[sub_index]
    plot_path =  paste0(folder, "/", gene, "_funnel.pdf")
    
    pdf(file = plot_path, width = 10, height = 7)
    print(plot_funnel_meta_gene(gene_name = gene, meta_df = selected_reduced_df))
    dev.off()
    
  }
  
}


################### Heterogeneity stats on Meta ###################
het_dfs_meta = c(
  "meta_no_covar_all_brain",
  "meta_no_covar_all_brain_signif",
  "meta_no_covar_cortical",
  "meta_no_covar_cortical_signif",
  "meta_no_covar_prefrontal",
  "meta_no_covar_prefrontal_signif",
  "meta_with_covar_all_brain",
  "meta_with_covar_all_brain_signif",
  "meta_with_covar_cortical",
  "meta_with_covar_cortical_signif",
  "meta_with_covar_prefrontal",
  "meta_with_covar_prefrontal_signif",
  "meta_SV_all_brain",
  "meta_SV_all_brain_signif",
  "meta_SV_cortical",
  "meta_SV_cortical_signif",
  "meta_SV_prefrontal",
  "meta_SV_prefrontal_signif",
  "REML_meta_no_covar_all_brain",
  "REML_meta_no_covar_all_brain_signif",
  "REML_meta_no_covar_cortical",
  "REML_meta_no_covar_cortical_signif",
  "REML_meta_no_covar_prefrontal",
  "REML_meta_no_covar_prefrontal_signif",
  "REML_meta_with_covar_all_brain",
  "REML_meta_with_covar_all_brain_signif",
  "REML_meta_with_covar_cortical",
  "REML_meta_with_covar_cortical_signif",
  "REML_meta_with_covar_prefrontal",
  "REML_meta_with_covar_prefrontal_signif"
)


label_vector = vector()

for (i in het_dfs_meta){
  label_vector[i] = rename_function(i)
}

het_dfs_meta = lapply(het_dfs_meta, function(x){
  
  df = get(x)
  stats_on_meta = summary(df$I2)
  names_of_stats = names(stats_on_meta)
  stats_on_meta = as.numeric(stats_on_meta)
  names(stats_on_meta) = names_of_stats
  
  stats_on_meta_2 = summary(df$tau2)
  names_of_stats_2 = names(stats_on_meta_2)
  stats_on_meta_2 = as.numeric(stats_on_meta_2)
  names(stats_on_meta_2) = names_of_stats_2
  
  out_df = data.frame(
    name=label_vector[x],
    minimal_I2 = stats_on_meta["Min."],
    median_I2 = stats_on_meta["Median"],
    maximal_I2 = stats_on_meta["Max."],
    minimal_tau2 = stats_on_meta_2["Min."],
    median_tau2 = stats_on_meta_2["Median"],
    maximal_tau2 = stats_on_meta_2["Max."]
    
  )
  return(out_df)
})

het_dfs_meta = do.call(rbind, het_dfs_meta)
colnames(het_dfs_meta) = c(
  "Analysis",
  "Min. I2",
  "Med. I2",
  "Max. I2",
  "Min. tau2",
  "Med. tau2",
  "Max. tau2"
)

wb = createWorkbook()
# Add a new worksheet with the sheet name
addWorksheet(wb, sheetName = "Heterogeneity stats")
# Write the data frame to the worksheet
writeData(wb, sheet = "Heterogeneity stats", het_dfs_meta)
# Adjust column widths to fit the text
setColWidths(wb, sheet = "Heterogeneity stats", cols = 1:ncol(het_dfs_meta), widths = "auto")
# Adjust row heights
setRowHeights(wb, sheet = "Heterogeneity stats", rows = 1:nrow(het_dfs_meta), heights = "auto")
# Save
saveWorkbook(wb, "Heterogeneity stats.xlsx", overwrite = TRUE)

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

meta_list_SV = list(
  "meta_SV_all_brain" = meta_SV_all_brain,
  "meta_SV_cortical" = meta_SV_cortical,
  "meta_SV_prefrontal" = meta_SV_prefrontal
)

meta_list_SV_signif = list(
  "meta_SV_all_brain_signif" = meta_SV_all_brain_signif,
  "meta_SV_cortical_signif" = meta_SV_cortical_signif,
  "meta_SV_prefrontal_signif" = meta_SV_prefrontal_signif
)

################### Comparison with PSY SUAS ###################
suas_genes_PSY = c("PLA2G10", "RBKS", "CCL27", "ASRG1", "WWP2")
meta_no_covar_all_brain[meta_no_covar_all_brain$gene %in% suas_genes_PSY,]
meta_with_covar_all_brain[meta_with_covar_all_brain$gene %in% suas_genes_PSY,]
meta_SV_all_brain[meta_SV_all_brain$gene %in% suas_genes_PSY,]
blood_df_no_covar[blood_df_no_covar$Corrected_symbol %in% suas_genes_PSY,]
blood_df_with_covar[blood_df_with_covar$Corrected_symbol %in% suas_genes_PSY,]

meta_no_covar_cortical[meta_no_covar_cortical$gene %in% suas_genes_PSY,]
meta_with_covar_cortical[meta_with_covar_cortical$gene %in% suas_genes_PSY,]
meta_SV_cortical[meta_SV_cortical$gene %in% suas_genes_PSY,]

meta_no_covar_prefrontal[meta_no_covar_prefrontal$gene %in% suas_genes_PSY,]
meta_with_covar_prefrontal[meta_with_covar_prefrontal$gene %in% suas_genes_PSY,]
meta_SV_prefrontal[meta_SV_prefrontal$gene %in% suas_genes_PSY,]

################### Stats on matching with blood ###################

meta_list_signif_combined = list(
  "meta_no_covar_all_brain_signif" = meta_no_covar_all_brain_signif,
  "meta_no_covar_cortical_signif" = meta_no_covar_cortical_signif,
  "meta_no_covar_prefrontal_signif" = meta_no_covar_prefrontal_signif,
  "meta_with_covar_all_brain_signif" = meta_with_covar_all_brain_signif,
  "meta_with_covar_cortical_signif" = meta_with_covar_cortical_signif,
  "meta_with_covar_prefrontal_signif" = meta_with_covar_prefrontal_signif,
  "meta_SV_all_brain_signif" = meta_SV_all_brain_signif,
  "meta_SV_cortical_signif" = meta_SV_cortical_signif,
  "meta_SV_prefrontal_signif" = meta_SV_prefrontal_signif
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

names_for_list = vector()

for (i in names(meta_list_signif_combined)){
  names_for_list[i] = rename_function(i)
}

Stats_matching_brain_blood_meta = vector()

for (i in 1:length(meta_list_signif_combined)){
  
  message = compare_with_blood(meta_list_signif_combined[[i]], names_for_list[i])
  Stats_matching_brain_blood_meta = c(Stats_matching_brain_blood_meta, message)
  
}
writeLines(text = Stats_matching_brain_blood_meta, con = "Stats_matching_brain_blood_meta.txt")

# Figure

blood_cor_plot = ggplot(data = meta_no_covar_all_brain, aes(x = meta_LFc, y = mean_blood_lfc)) +
  geom_point(alpha = .5, col = "#8B0000", size = 2) +
  labs(x = "meta-estimated Log2FC (all brain)", y = "average Log2FC in GSE247998 (blood)") +
  # Custom the theme:
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_line(size = 0.1, linetype = 2, color =  "black"),
    panel.grid.major.y = element_line(size = 0.1, linetype = 2, color =  "black"),
    axis.line = element_line(size = 1, linetype = 1, color =  "black"),
    panel.background = element_blank(),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  )

ggsave(file = "Fig_S2.png", plot = blood_cor_plot, width=2560, height=1440, units = "px", scale = 1)
nrow(meta_no_covar_all_brain) - 35219 # 15859 transcripts plotted

################### RRA and Random effects stats ###################
meta_list_RRA_stats_combined = list(
  "meta_no_covar_all_brain" = meta_no_covar_all_brain,
  "meta_no_covar_cortical" = meta_no_covar_cortical,
  "meta_no_covar_prefrontal" = meta_no_covar_prefrontal,
  "meta_with_covar_all_brain" = meta_with_covar_all_brain,
  "meta_with_covar_cortical" = meta_with_covar_cortical,
  "meta_with_covar_prefrontal" = meta_with_covar_prefrontal,
  "meta_SV_all_brain" = meta_SV_all_brain,
  "meta_SV_cortical" = meta_SV_cortical,
  "meta_SV_prefrontal" = meta_SV_prefrontal
)

names_for_list = vector()
for (i in names(meta_list_RRA_stats_combined)){
  names_for_list[i] = rename_function(i)
}

compare_with_RRA = function(dataset, name){
  
  Data_name = name
  Data_name = paste0("######### Summary of RRA vs Random Effects Meta on ", Data_name, " #########")
  
  # RE subset
  signif_subset_nominal_RE = dataset[dataset$meta_pval < 0.05,]
  signif_subset_nominal_RE = signif_subset_nominal_RE[!is.na(signif_subset_nominal_RE$meta_pval),]
  signif_subset_nominal_RE = dplyr::arrange(signif_subset_nominal_RE, -meta_LFc)
  
  signif_subset_FDR_RE = dataset[dataset$meta_FDR < 0.05,]
  signif_subset_FDR_RE = signif_subset_FDR_RE[!is.na(signif_subset_FDR_RE$meta_pval),]
  signif_subset_FDR_RE = dplyr::arrange(signif_subset_FDR_RE, -meta_LFc)
  
  # RRA subset
  signif_subset_nominal_RRA = dataset[dataset$P_RRA < 0.05,]
  signif_subset_FDR_RRA = dataset[dataset$FDR_RRA < 0.05,]
  signif_subset_nominal_RRA = dplyr::arrange(signif_subset_nominal_RRA, P_RRA)
  signif_subset_FDR_RRA = dplyr::arrange(signif_subset_FDR_RRA, P_RRA)
  
  # Overlap
  signif_subset_nominal_overlap = signif_subset_nominal_RE[signif_subset_nominal_RE$P_RRA < 0.05,]
  
  # RE
  n_genes_RE_nom_sig_up = signif_subset_nominal_RE %>%
    filter(., meta_LFc > 0) %>%
    pull(gene) %>%
    length() %>%
    paste0("N genes RE p<0.05, log2FC>0: ", .)
  
  n_genes_RE_nom_sig_down = signif_subset_nominal_RE %>%
    filter(., meta_LFc < 0) %>%
    pull(gene) %>%
    length() %>%
    paste0("N genes RE p<0.05, log2FC<0: ", .)
  
  n_genes_RE_FDR = signif_subset_FDR_RE %>%
    filter(., meta_LFc < 0) %>%
    pull(gene) %>%
    length() %>%
    paste0("N genes RE FDR<0.05: ", .)
  
  
  top_genes_RE_nom_sig_up = signif_subset_nominal_RE %>%
    filter(., meta_LFc > 0) %>%
    arrange(., -meta_LFc) %>%
    pull(gene) %>%
    head(., 10) %>%
    paste0(collapse="; ")  %>%
    paste0("Top 10 (or less) upregulated genes RE p<0.05 sorted by log2FC:\n", .)
  
  top_genes_RE_nom_sig_down = signif_subset_nominal_RE %>%
    filter(., meta_LFc < 0) %>%
    arrange(., meta_LFc) %>%
    pull(gene) %>%
    head(., 10) %>%
    paste0(collapse="; ") %>%
    paste0("Top 10 (or less) downregulated genes RE p<0.05 sorted by log2FC:\n", .)
  
  
  # RRA
  n_genes_RRA_nom_sig = signif_subset_nominal_RRA %>%
    pull(gene) %>%
    length()  %>%
    paste0("N genes RRA p<0.05: ", .)
  n_genes_RRA_FDR = signif_subset_FDR_RRA %>%
    pull(gene) %>%
    length() %>%
    paste0("N genes RRA FDR<0.05: ", .)
  
  top_genes_RRA_nom_sig = signif_subset_nominal_RRA %>%
    arrange(., P_RRA) %>%
    pull(gene) %>%
    head(., 10) %>%
    paste0(collapse="; ") %>%
    paste0("Top 10 (or less) genes RRA p<0.05 sorted by p-value:\n", .)
  
  top_genes_RRA_FDR = signif_subset_FDR_RRA %>%
    arrange(., P_RRA) %>%
    pull(gene) %>%
    head(., 10) %>%
    paste0(collapse="; ") %>%
    paste0("Top 10 (or less) genes RRA FDR<0.05 sorted by p-value:\n", .)
  
  # Overlap
  n_overlaping_nom_sig = signif_subset_nominal_overlap %>%
    pull(gene) %>%
    length()  %>%
    paste0("N genes in overlap of RE and RRA, p<0.05: ", .)
  
  overlapping_genes_RE_RRA = signif_subset_nominal_overlap %>%
    arrange(., P_RRA) %>%
    pull(gene) %>%
    paste0(collapse="; ")%>%
    paste0("Overlapping genes between RE and RRA, p<0.05:\n", .)
  
  # Correlations
  non_na_substet = dataset[!is.na(meta_no_covar_all_brain$meta_LFc), ]
  p_corr_all = cor.test(non_na_substet$meta_pval,
                        non_na_substet$P_RRA, method = "spearman")$estimate
  p_corr_all = paste0("Spearman correlation for p-values in RE vs RRA (all genes): ", p_corr_all)
  
  p_corr_overlap = cor.test(signif_subset_nominal_overlap$meta_pval,
                            signif_subset_nominal_overlap$P_RRA, method = "spearman")$estimate
  p_corr_overlap = paste0("Spearman correlation for p-values in RE vs RRA (overlap): ", p_corr_overlap)
  
  
  separator_small = "\n"
  separator_large = "\n\n\n\n"
  output = c(Data_name, 
             separator_small,
             n_genes_RE_nom_sig_up,
             separator_small,
             n_genes_RE_nom_sig_down,
             separator_small,
             top_genes_RE_nom_sig_up,
             separator_small,
             top_genes_RE_nom_sig_down,
             separator_small,
             n_genes_RRA_nom_sig,
             separator_small,
             n_genes_RRA_FDR,
             separator_small,
             top_genes_RRA_nom_sig,
             separator_small,
             top_genes_RRA_FDR,
             separator_small,
             n_overlaping_nom_sig,
             separator_small,
             overlapping_genes_RE_RRA,
             separator_small,
             p_corr_all,
             separator_small,
             p_corr_overlap,
             separator_large)
  output = paste0(output, collapse = "")
  return(output)
}

Stats_matching_RE_RRA = vector()

for (i in 1:length(meta_list_signif_combined)){
  
  message = compare_with_RRA(meta_list_RRA_stats_combined[[i]], names_for_list[i])
  Stats_matching_RE_RRA = c(Stats_matching_RE_RRA, message)
  
}
writeLines(text = Stats_matching_RE_RRA, con = "Stats_RE_vs_RRA.txt")

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
joined_list = c(meta_list_no_covar_signif, meta_list_with_covar_signif, meta_list_SV_signif)
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
  "DLPFC (covar)",
  "All brain (SVs)",
  "Cortical regions (SVs)",
  "DLPFC (SVs)"
)

meta_gene_joined_list = lapply(1:9, function(index){
  data_list = joined_list[[index]]
  name = names_for_list[index]
  output = data.frame(gene = data_list$gene, effect_size = data_list$meta_LFc, analysis = name)
  return(output)
})
meta_gene_joined_list = do.call(rbind, meta_gene_joined_list)
combined_graph_df = meta_gene_joined_list[,c("analysis", "gene")]
colnames(combined_graph_df) = c("Start", "Target")

table_graph_targets = table(combined_graph_df$Target)
table_graph_targets = as.data.frame(table_graph_targets)
table_graph_targets = dplyr::arrange(table_graph_targets,-Freq)

# preparing data
network_data = graph_from_data_frame(d = combined_graph_df, directed = TRUE)
network_data = toVisNetworkData(network_data)

network_data$nodes$group = sapply(network_data$nodes$id, function(x){
if (stri_detect_fixed(x, pattern = "(covar)")){
  return("Analysis Covar")
}
if (stri_detect_fixed(x, pattern = "(SVs)")){
  return("Analysis SVs")
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
  visGroups(groupname = "Analysis", color = list("background" = "orange", border = "black"), size=20) %>% 
  visGroups(groupname = "Analysis Covar", color = list("background" = "blue", border = "black"), size=20) %>%
  visGroups(groupname = "Analysis SVs", color = list("background" = "skyblue", border = "black"), size=20) %>% 
  visGroups(groupname = "Gene", color = list("background" = "red", border = "black")) %>% 
  visNodes(size = 15, borderWidth = 1, font = list(size=30)) %>%
  visEdges(smooth = list("roundness" = 0.2, "type" = "diagonalCross"), arrows = list(to = list(enabled = TRUE, scaleFactor = 1)), color = "grey") %>%
  visPhysics(barnesHut = list(gravitationalConstant = -60000, centralGravity=0.001)) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = F)
visSave(net, file = "Meta_net.html")
Sys.setenv("OPENSSL_CONF"="/dev/null")
webshot("Meta_net.html", "Meta_net_2.png", vwidth = 2000, vheight = 2000, zoom = 1, delay = 10)

# Stats on graph
combined_graph_df_stats_base = table(combined_graph_df$Target)
combined_graph_df_stats_base = as.data.frame(combined_graph_df_stats_base)
combined_graph_df_stats_base = dplyr::arrange(combined_graph_df_stats_base, -Freq)

################### RRA-RE Network visualization for meta ###################
joined_list = c(meta_list_no_covar_signif, meta_list_with_covar_signif, meta_list_SV_signif)
joined_list = lapply(joined_list, function(x){
  x = x[x$P_RRA < 0.05,]
  x = x[abs(x$meta_LFc) >= 0.1,]
  return(x)
})

names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)",
  "All brain (SVs)",
  "Cortical regions (SVs)",
  "DLPFC (SVs)"
)

meta_gene_joined_list = lapply(1:9, function(index){
  data_list = joined_list[[index]]
  name = names_for_list[index]
  
  if (nrow(data_list) < 1){
    return(NA)
  }
  output = data.frame(gene = data_list$gene, effect_size = data_list$meta_LFc, analysis = name)
  return(output)
})
meta_gene_joined_list = meta_gene_joined_list[sapply(meta_gene_joined_list, is.data.frame)]
meta_gene_joined_list = do.call(rbind, meta_gene_joined_list)
combined_graph_df = meta_gene_joined_list[,c("analysis", "gene")]
colnames(combined_graph_df) = c("Start", "Target")

table_graph_targets = table(combined_graph_df$Target)
table_graph_targets = as.data.frame(table_graph_targets)
table_graph_targets = dplyr::arrange(table_graph_targets,-Freq)

head(table_graph_targets, 20)

# preparing data
network_data = graph_from_data_frame(d = combined_graph_df, directed = TRUE)
network_data = toVisNetworkData(network_data)

network_data$nodes$group = sapply(network_data$nodes$id, function(x){
  if (stri_detect_fixed(x, pattern = "(covar)")){
    return("Analysis Covar")
  }
  if (stri_detect_fixed(x, pattern = "(SVs)")){
    return("Analysis SVs")
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
  visGroups(groupname = "Analysis", color = list("background" = "orange", border = "black"), size=20) %>% 
  visGroups(groupname = "Analysis Covar", color = list("background" = "blue", border = "black"), size=20) %>%
  visGroups(groupname = "Analysis SVs", color = list("background" = "skyblue", border = "black"), size=20) %>% 
  visGroups(groupname = "Gene", color = list("background" = "red", border = "black")) %>% 
  visNodes(size = 15, borderWidth = 1, font = list(size=30)) %>%
  visEdges(smooth = list("roundness" = 0.2, "type" = "diagonalCross"), arrows = list(to = list(enabled = TRUE, scaleFactor = 1)), color = "grey") %>%
  visPhysics(barnesHut = list(gravitationalConstant = -60000, centralGravity=0.001)) %>%
  visOptions(highlightNearest = list(enabled = T, hover = T), 
             nodesIdSelection = F)

visSave(net, file = "Meta_net_RRA_RE.html")
Sys.setenv("OPENSSL_CONF"="/dev/null")
webshot("Meta_net_RRA_RE.html", "Meta_net_RRA_RE.png", vwidth = 2000, vheight = 2000, zoom = 1, delay = 10)

# Stats on graph
combined_graph_df_stats = table(combined_graph_df$Target)
combined_graph_df_stats = as.data.frame(combined_graph_df_stats)
combined_graph_df_stats = dplyr::arrange(combined_graph_df_stats, -Freq)

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

joined_list_enrichment = c(meta_list_no_covar_signif, meta_list_with_covar_signif, meta_list_SV_signif)
names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)",
  "All brain (SVs)",
  "Cortical regions (SVs)",
  "DLPFC (SVs)"
)
joined_list_enrichment = lapply(joined_list_enrichment, function(x){
  data = x$gene
  return(data)
})

# running enrichment
for (x in 1:length(joined_list_enrichment)){
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
joined_list_significant = c(meta_list_no_covar_signif, meta_list_with_covar_signif, meta_list_SV_signif)
joined_list_significant = lapply(joined_list_significant, function(x){
  x = dplyr::arrange(x, -meta_LFc)
  return(x)
})

names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)",
  "All brain (SVs)",
  "Cortical regions (SVs)",
  "DLPFC (SVs)"
)

joined_list_significant_high_LFc = lapply(1:length(joined_list_significant), function(x){
  df = joined_list_significant[[x]]
  df = df[abs(df$meta_LFc) >= 0.2,]
  df$analysis = names_for_list[x]
  return(df)
})

unique_significant_genes = lapply(joined_list_significant_high_LFc, function(x){
  data = x$gene
  return(data)
})
unique_significant_genes = unlist(unique_significant_genes)
unique_significant_genes = unique(unique_significant_genes)

old_unique_genes = read.csv("OLD/GPT_5_classified_symbols__fixed.csv")
unique_significant_genes[unique_significant_genes %!in% old_unique_genes$gene_symbol]

writeLines(unique_significant_genes, "unique_significant_genes_for_AI.txt")
gpt_classified_genes = read.csv("GPT_5_classified_symbols__fixed.csv")


stats_high_LFc = do.call(rbind, joined_list_significant_high_LFc)
stats_high_LFc$GPT_5_class = sapply(stats_high_LFc$gene, function(x){
  class =gpt_classified_genes[gpt_classified_genes$gene_symbol == x, "gene_class"]
  return(class)
})


stats_high_LFc = as.data.frame(table(stats_high_LFc$GPT_5_class, stats_high_LFc$analysis))

colnames(stats_high_LFc) = c("Category", "Analysis", "Freq")

plot = ggplot(data = stats_high_LFc, aes(x=Analysis, y = Freq, fill = Category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6, col = "black", size = 0.25) +
  labs(y = "Gene count", x = "Analysis type", fill = "Gene category (GPT-5)", title = "Gene categories with abs(logFC)>=0.2") +
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
  
ggsave(file = "stats_high_LFc.pdf", plot = plot, width = 22, height = 12, units = "cm")


################### Saving processed files ###################
joined_list_significant = c(meta_list_no_covar_signif, meta_list_with_covar_signif, meta_list_SV_signif)
joined_list_significant = lapply(joined_list_significant, function(x){
  x = dplyr::arrange(x, -meta_LFc)
  return(x)
})

names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)",
  "All brain (SVs)",
  "Cortical regions (SVs)",
  "DLPFC (SVs)"
)

wb = createWorkbook()

for (i in 1:length(joined_list_significant)) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, names_for_list[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = names_for_list[i], joined_list_significant[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = names_for_list[i], cols = 1:ncol(joined_list_significant[[i]]), widths = "auto")
}

saveWorkbook(wb, "Meta_suicide_significant_genes.xlsx", overwrite = TRUE)

################### Saving full files ###################
joined_list_meta_all_data = c(meta_list_no_covar, meta_list_with_covar, meta_list_SV)
names_for_list = c(
  "All brain",
  "Cortical regions",
  "DLPFC",
  "All brain (covar)",
  "Cortical regions (covar)",
  "DLPFC (covar)",
  "All brain (SVs)",
  "Cortical regions (SVs)",
  "DLPFC (SVs)"
)
names(joined_list_meta_all_data)

dir.create("full_meta_results")

for (i in 1:length(joined_list_meta_all_data)) {
  
  path_to_save = paste0("full_meta_results/",names_for_list[i], ".csv")
  
  # Save data as CSV
  write.csv(joined_list_meta_all_data[[i]], file = path_to_save, row.names = FALSE)
}

################### Inspection of some genes ###################
# BDNF
lapply(joined_list_meta_all_data, function(x) x[x$gene == "BDNF",])
plot_forest_meta_gene(gene_name = "BDNF", meta_df = combined_df_no_covar_meta_full_list_reduced)
plot_forest_meta_gene(gene_name = "BDNF", meta_df = combined_df_no_covar_meta_cortical_list_reduced)
plot_forest_meta_gene(gene_name = "BDNF", meta_df = combined_df_no_covar_meta_prefrontal_list_reduced)

plot_forest_meta_gene(gene_name = "BDNF", meta_df = combined_df_with_covar_meta_full_list_reduced)
plot_forest_meta_gene(gene_name = "BDNF", meta_df = combined_df_with_covar_meta_cortical_list_reduced)
plot_forest_meta_gene(gene_name = "BDNF", meta_df = combined_df_with_covar_meta_prefrontal_list_reduced)

# VEGF
lapply(joined_list_meta_all_data, function(x) x[x$gene == "VEGFA",])
plot_forest_meta_gene(gene_name = "VEGFA", meta_df = combined_df_no_covar_meta_full_list_reduced)
plot_forest_meta_gene(gene_name = "VEGFA", meta_df = combined_df_no_covar_meta_cortical_list_reduced)
plot_forest_meta_gene(gene_name = "VEGFA", meta_df = combined_df_no_covar_meta_prefrontal_list_reduced)

plot_forest_meta_gene(gene_name = "VEGFA", meta_df = combined_df_with_covar_meta_full_list_reduced)
plot_forest_meta_gene(gene_name = "VEGFA", meta_df = combined_df_with_covar_meta_cortical_list_reduced)
plot_forest_meta_gene(gene_name = "VEGFA", meta_df = combined_df_with_covar_meta_prefrontal_list_reduced)


################### Mini stats ###################
a = c(49.2, 48.9, 48)
b = c(12.7, 12.8, 14.5)
mean(b-a)



