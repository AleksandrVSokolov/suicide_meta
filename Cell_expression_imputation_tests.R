# Comments
# This code is intended to be viewed in RStudio and contains appropriate headings
# "=" symbol was used as an assignment operator
# Code for CIBERSORTx runs is presented as character strings within the file
# Real token for CIBERSORTx docker container is replaced with <token_from_cibersortx_website>
# This file contains only test runs

# Setting options
getOption("scipen") # Default number notation is 0
options(scipen=999)
options(stringsAsFactors = FALSE)
options(show.error.messages = TRUE)
setwd("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project")


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
library(Matrix)
library(Seurat)
library(magick)

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


################### Cohort overview ###################

#### GSE102556
# DLPFC, no cell sorting reported

#### GSE243356
# Temporal cortex, no cell sorting reported

#### GSE248260
# Ventral white matter, no cell sorting reported

#### GSE202537
# Nac, no cell sorting reported

#### GSE101521
# DLPFC, Attempted to enrich gray matter. No cell sorting reported

#### GSE144136
# DLPFC, single-cell seq. Gray matter samples

#### GSE213982
# DLPFC, single-cell seq. Gray matter samples



################### CIBERSORTx file requirements ###################

# We need -> Impute Gene Expression, High-Resolution Mode (Tutorial 5) -> 

# Before this -> Tutorial 1 - "Build a Signature Matrix File from Single-Cell RNA Sequencing Data"


################### Base SC Matrix File preparation ###################

# Read phenotypes
ref_cell_pheno_GSE144136 = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE144136_results/GSE144136_pheno_curated.csv")
ref_cell_pheno_GSE213982 = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE213982_results/GSE213982_pheno_curated.csv")

# Read annotation files
TMP_sig_matrix_gene_rows = read.csv("Data_preprocessing_analysis/GSE213982/GSE213982_combined_counts_matrix_genes_rows.csv")
TMP_sig_matrix_gene_rows = TMP_sig_matrix_gene_rows$x
length(TMP_sig_matrix_gene_rows) # 36588
any(duplicated(TMP_sig_matrix_gene_rows)) # FALSE

TMP_sig_matrix_cells_colums =  read.csv("Data_preprocessing_analysis/GSE213982/GSE213982_combined_counts_matrix_cells_columns.csv")
TMP_sig_matrix_cells_colums = TMP_sig_matrix_cells_colums$x
any(duplicated(TMP_sig_matrix_cells_colums)) # FALSE

GSE144136_GSE213982_count_matrix = readMM("Data_preprocessing_analysis/GSE213982/GSE213982_combined_counts_matrix.mtx")
rownames(GSE144136_GSE213982_count_matrix) = TMP_sig_matrix_gene_rows
colnames(GSE144136_GSE213982_count_matrix) = TMP_sig_matrix_cells_colums

GSE144136_GSE213982_metadata = data.frame("sample_name" = TMP_sig_matrix_cells_colums)

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
table(GSE144136_GSE213982_metadata$person, GSE144136_GSE213982_metadata$cell_overall)

################### Base Signature matrix File preparation ###################

# In total we have 15+17+18+18=68 people
# We can pick 10 males and 10 females as "Validation set" for matrix
# CIBERSORTx should work better if number of samples is way higher than cell types
length(unique(GSE144136_GSE213982_metadata$person)) #72 transcriptomes in total (4 samples are excluded in DE)

female_vector = GSE144136_GSE213982_metadata$person
female_vector = female_vector[stri_detect_fixed(female_vector, "F")]
female_vector = unique(female_vector)
set.seed(1234)
selected_validation_females = sample(female_vector, size=10, replace=FALSE)
set.seed(NULL)

male_vector = GSE144136_GSE213982_metadata$person
male_vector = male_vector[stri_detect_fixed(male_vector, "M")]
male_vector = unique(male_vector)
set.seed(1234)
selected_validation_males = sample(male_vector, size=10, replace=FALSE)
set.seed(NULL)
selected_validation_samples = c(selected_validation_females, selected_validation_males)

# selected_validation_samples
# "F34" "F23" "F29" "F17" "F13" "F8"  "F12" "F10" "F15" "F7"  "M34" "M23" "M32" "M29" "M13" "M2"  "M22" "M17" "M5"  "M14"
# "F19" "F5" are not used for validation matrix!

GSE144136_GSE213982_metadata$is_validation = GSE144136_GSE213982_metadata$person %in% selected_validation_samples


# Pseudobulked validation counts
validation_bulk_counts = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
rownames(validation_bulk_counts) = validation_bulk_counts$X
validation_bulk_counts$X = NULL
validation_bulk_counts = validation_bulk_counts[,colnames(validation_bulk_counts) %in% selected_validation_samples]

validation_cell_prop_1 = ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %in% selected_validation_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
validation_cell_prop_2 = ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %in% selected_validation_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
validation_cell_full = rbind(validation_cell_prop_1, validation_cell_prop_2)
validation_bulk_counts = validation_bulk_counts[,validation_cell_full$PARTICIPANT]

validation_bulk_counts_df = as.data.frame(validation_bulk_counts)
genes_column = data.frame("Gene" = rownames(validation_bulk_counts_df))
validation_bulk_counts_df = cbind(genes_column, validation_bulk_counts_df)
rownames(validation_bulk_counts_df) = NULL

fwrite(validation_bulk_counts_df,
       file = "validation_bulk_counts_df.txt",
       sep = "\t",
       row.names = FALSE,
       col.names = TRUE,
       quote = FALSE)


########## Subsetting matrix and exporting
# Males and Females come from slightly different platforms
# However, counts were obtained from the same alignment run
# Probably safe to use as mixed sample and would make sense to make robust matrix

GSE144136_GSE213982_count_matrix_processed = GSE144136_GSE213982_count_matrix
colnames(GSE144136_GSE213982_count_matrix_processed) = GSE144136_GSE213982_metadata$cell_overall

# Creating dense matrix
dense_count_matrix_processed = as.matrix(GSE144136_GSE213982_count_matrix_processed)
colnames(dense_count_matrix_processed)[1:10]
rownames(dense_count_matrix_processed)[1:10]

all(colnames(GSE144136_GSE213982_count_matrix_processed) == colnames(dense_count_matrix_processed))
all(rownames(GSE144136_GSE213982_count_matrix_processed) == rownames(dense_count_matrix_processed))

# Sampling 10000 cells from non-validation positions!
# This is hard cap on the website and maybe the bests option https://pmc.ncbi.nlm.nih.gov/articles/PMC7695353/#S28
# Req. 12 CIBERSORTx Gene Expression Analysis Mode works best when the number of samples to deconvolve is much larger than the number of cell types in the signature matrix. 

allowwed_indices = 1:ncol(dense_count_matrix_processed)
allowwed_indices = allowwed_indices[!GSE144136_GSE213982_metadata$is_validation]
set.seed(1234)
sampled_cell_idx = sample(allowwed_indices, size = 10000, replace = FALSE)
set.seed(NULL)

sampled_metadata_sig_matrix = GSE144136_GSE213982_metadata[sampled_cell_idx,]

table(sampled_metadata_sig_matrix$is_validation)
table(sampled_metadata_sig_matrix$person)
table(sampled_metadata_sig_matrix$cell_overall)
table(sampled_metadata_sig_matrix$person, 
      sampled_metadata_sig_matrix$cell_overall)

length(unique(sampled_metadata_sig_matrix$person)) # 52

length(unique(GSE144136_GSE213982_metadata[GSE144136_GSE213982_metadata$is_validation, "person"])) # 20
unique(GSE144136_GSE213982_metadata[GSE144136_GSE213982_metadata$is_validation, "person"])
# "F10" "F12" "F13" "F15" "F17" "F23" "F29" "F34" "F7"  "F8"  "M13" "M14" "M17" "M2"  "M22" "M23" "M29" "M32" "M34" "M5" 

dense_count_matrix_processed = dense_count_matrix_processed[,sampled_cell_idx]
all(colnames(dense_count_matrix_processed) == sampled_metadata_sig_matrix$cell_overall) # TRUE
dense_count_matrix_processed_df = as.data.frame(dense_count_matrix_processed)
genes_column = data.frame("Gene" = rownames(dense_count_matrix_processed_df))
dense_count_matrix_processed_df = cbind(genes_column, dense_count_matrix_processed_df)
colnames(dense_count_matrix_processed_df)[1:10]
rownames(dense_count_matrix_processed_df) = NULL

dim(dense_count_matrix_processed_df) # 36588 10001

# Saving as TSV in TXT file
fwrite(dense_count_matrix_processed_df,
            file = "ref_count_matrix_dense.txt", 
            sep = "\t",
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

# In the paper https://pubmed.ncbi.nlm.nih.gov/31960376/
"
While increasing the number of cells per phenotype and the number of biological
replicates can improve the quality of the signature matrix (up to ~20 cells per
phenotype and up to 2-3 donor samples), simulation experiments and empirical
observations suggest that as few as 3 cells per phenotype and as few as one
donor sample can still generate reliable results [23]. For this reason, and to limit
the amount of space and time necessary to run a CIBERSORTx job, we
recommend that users restrict the number of cells to at most 5,000 when
uploading the scRNA-seq reference profile to the CIBERSORTx website
(currently a strict upper limit of 10,000 cells is allowed). 
"


################### Signature matrix data (male only) ###################
# We limit to male dataset (GSE144136) used one chemistry
# We will use only 3 people with largest number of cells. Myaybe excessive amount of cells introduces too much noise for signature matrix
GSE144136_count_matrix_processed = GSE144136_GSE213982_count_matrix
GSE144136_GSE213982_metadata$male_selector = stri_detect_fixed(GSE144136_GSE213982_metadata$person, "M")

GSE144136_metadata = GSE144136_GSE213982_metadata[GSE144136_GSE213982_metadata$male_selector, ]
GSE144136_count_matrix_processed = GSE144136_count_matrix_processed[,GSE144136_metadata$sample_name]

# Creating dense matrix
dense_count_matrix_processed = as.matrix(GSE144136_count_matrix_processed)
colnames(dense_count_matrix_processed)[1:10]
rownames(dense_count_matrix_processed)[1:10]

# Selecting most prevalent people
GSE144136_metadata_person_stat = as.data.frame(table(GSE144136_metadata$person))
GSE144136_metadata_person_stat = GSE144136_metadata_person_stat[GSE144136_metadata_person_stat$Var1 %!in% selected_validation_samples, ]

# M20, M11, M15
table(GSE144136_metadata$person, GSE144136_metadata$cell_overall)
"
M11    622   27 2081  618   49   22  365  117
M20    292   32 1979  728  101   84  966  248
M15    419   39 1938  776   16   27  328  272
"

GSE144136_metadata$large_samples = GSE144136_metadata$person %in% c("M11", "M20", "M15")
sampled_metadata_sig_matrix_male = GSE144136_metadata[GSE144136_metadata$large_samples,]

table(sampled_metadata_sig_matrix_male$is_validation)
table(sampled_metadata_sig_matrix_male$person)
table(sampled_metadata_sig_matrix_male$cell_overall)
table(sampled_metadata_sig_matrix_male$person, 
      sampled_metadata_sig_matrix_male$cell_overall)

dense_count_matrix_processed = dense_count_matrix_processed[,GSE144136_metadata$large_samples]
colnames(dense_count_matrix_processed) = sampled_metadata_sig_matrix_male$cell_overall

dense_count_matrix_processed_df = as.data.frame(dense_count_matrix_processed)
genes_column = data.frame("Gene" = rownames(dense_count_matrix_processed_df))
dense_count_matrix_processed_df = cbind(genes_column, dense_count_matrix_processed_df)
colnames(dense_count_matrix_processed_df)[1:10]
rownames(dense_count_matrix_processed_df) = NULL


# Saving as TSV in TXT file
fwrite(dense_count_matrix_processed_df,
       file = "ref_count_matrix_dense_male.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Removing large files
rm(GSE144136_GSE213982_count_matrix_processed)
rm(GSE144136_count_matrix_processed)
rm(dense_count_matrix_processed)
rm(dense_count_matrix_processed_df)
gc()


################### Signature matrix data 30 000 cells ###################
GSE144136_GSE213982_count_matrix_processed = GSE144136_GSE213982_count_matrix
colnames(GSE144136_GSE213982_count_matrix_processed) = GSE144136_GSE213982_metadata$cell_overall

# Creating dense matrix
dense_count_matrix_processed = as.matrix(GSE144136_GSE213982_count_matrix_processed)
colnames(dense_count_matrix_processed)[1:10]
rownames(dense_count_matrix_processed)[1:10]

all(colnames(GSE144136_GSE213982_count_matrix_processed) == colnames(dense_count_matrix_processed))
all(rownames(GSE144136_GSE213982_count_matrix_processed) == rownames(dense_count_matrix_processed))

allowwed_indices = 1:ncol(dense_count_matrix_processed)
allowwed_indices = allowwed_indices[!GSE144136_GSE213982_metadata$is_validation]
set.seed(1234)
sampled_cell_idx = sample(allowwed_indices, size = 30000, replace = FALSE)
set.seed(NULL)

sampled_metadata_sig_matrix = GSE144136_GSE213982_metadata[sampled_cell_idx,]

table(sampled_metadata_sig_matrix$is_validation)
table(sampled_metadata_sig_matrix$person)
table(sampled_metadata_sig_matrix$cell_overall)
table(sampled_metadata_sig_matrix$person, 
      sampled_metadata_sig_matrix$cell_overall)


dense_count_matrix_processed = dense_count_matrix_processed[,sampled_cell_idx]
all(colnames(dense_count_matrix_processed) == sampled_metadata_sig_matrix$cell_overall) # TRUE
dense_count_matrix_processed_df = as.data.frame(dense_count_matrix_processed)
genes_column = data.frame("Gene" = rownames(dense_count_matrix_processed_df))
dense_count_matrix_processed_df = cbind(genes_column, dense_count_matrix_processed_df)
colnames(dense_count_matrix_processed_df)[1:10]
rownames(dense_count_matrix_processed_df) = NULL

dim(dense_count_matrix_processed_df) # 36588 10001

# Saving as TSV in TXT file
fwrite(dense_count_matrix_processed_df,
       file = "ref_count_matrix_dense_30K.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)


################### Sampled signature matrix ###################
# We can pick 5 males and 5 females
# Take approx 2000 cells of each type -> signature matrix
# Everything else -> Goes to validation
# We will set fraction to 0
GSE144136_GSE213982_metadata_sampled_matrix = GSE144136_GSE213982_metadata
GSE144136_GSE213982_metadata_sampled_matrix$is_validation = NULL

female_vector = GSE144136_GSE213982_metadata$person
female_vector = female_vector[stri_detect_fixed(female_vector, "F")]
female_vector = unique(female_vector)
set.seed(12345) # To try another seed
selected_signature_females = sample(female_vector, size=5, replace=FALSE)
set.seed(NULL)

male_vector = GSE144136_GSE213982_metadata$person
male_vector = male_vector[stri_detect_fixed(male_vector, "M")]
male_vector = unique(male_vector)
set.seed(12345) # To try another seed
selected_signature_males = sample(male_vector, size=5, replace=FALSE)
set.seed(NULL)
selected_signature_total = c(selected_signature_females, selected_signature_males)

# Inspection of cells VS participant
sampled_signature_metadata = GSE144136_GSE213982_metadata_sampled_matrix[GSE144136_GSE213982_metadata_sampled_matrix$person %in% selected_signature_total,]
table(sampled_signature_metadata$person, sampled_signature_metadata$cell_overall)

"       Ast  End  ExN  InN  Mic  Mix  Oli  OPC
  F21   85   97  530  265   65   37   54   82
  F23  142   74  436  183   32   27   32   46
  F30  113   46  288  160    8   58   17   58
  F32  207   44  671  285    6    6   51   81
  F34  165  103 1046  328  125   29  536  151
  M21   74    9 1265  364   11   66  237  149
  M23   14    2 3439  491    3   35   72    3
  M30   29   11 1377  567    4   85   70   94
  M32   88   25 1438  667   41   17  138  152
  M34   29    3 2165  783   38   18   95   61
"

# Solution: Ast, End, Mic, Mix, Oli, OPC -> TAKE all cells
# Solution: ExN, InN, -> Sample 2000 cells
tmp_ExN_indeces = which(sampled_signature_metadata$cell_overall == "ExN")
set.seed(12345)
tmp_selected_ExN_indeces = sample(tmp_ExN_indeces, size = 2000, replace = FALSE)
set.seed(NULL)

tmp_InN_indeces = which(sampled_signature_metadata$cell_overall == "InN")
set.seed(12345)
tmp_selected_InN_indeces = sample(tmp_InN_indeces, size = 2000, replace = FALSE)
set.seed(NULL)

sampled_signature_metadata$is_sampled = NA
for (i in 1:nrow(sampled_signature_metadata)){
  
  if (sampled_signature_metadata$cell_overall[i] %!in% c("ExN", "InN")){
    sampled_signature_metadata$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_ExN_indeces){
    sampled_signature_metadata$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_InN_indeces){
    sampled_signature_metadata$is_sampled[i] = TRUE
  } else {
    sampled_signature_metadata$is_sampled[i] = FALSE
  }
}

table(sampled_signature_metadata$is_sampled, sampled_signature_metadata$cell_overall)

"          Ast   End   ExN   InN   Mic   Mix   Oli   OPC
  FALSE     0     0 10655  2093     0     0     0     0
  TRUE    946   414  2000  2000   333   378  1302   877
"
sampled_signature_picked_cells = sampled_signature_metadata[sampled_signature_metadata$is_sampled, "sample_name"]

sampled_signature_matrix_processed = GSE144136_GSE213982_count_matrix
sampled_signature_matrix_processed = sampled_signature_matrix_processed[,sampled_signature_picked_cells]
dim(sampled_signature_matrix_processed) # 36588  8250

sampled_signature_matrix_processed = as.matrix(sampled_signature_matrix_processed)
dim(sampled_signature_matrix_processed) # 36588  8250
colnames(sampled_signature_matrix_processed)[1:10]
rownames(sampled_signature_matrix_processed)[1:10]

replacement_colnames = sapply(colnames(sampled_signature_matrix_processed), function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[3]
  return(x)
})
"replacement_colnames
 Ast  End  ExN  InN  Mic  Mix  Oli  OPC 
 946  414 2000 2000  333  378 1302  877 
"
colnames(sampled_signature_matrix_processed) = replacement_colnames
sampled_signature_matrix_processed = as.data.frame(sampled_signature_matrix_processed)
genes_column = data.frame("Gene" = rownames(sampled_signature_matrix_processed))
sampled_signature_matrix_processed = cbind(genes_column, sampled_signature_matrix_processed)
colnames(sampled_signature_matrix_processed)[1:10]
rownames(sampled_signature_matrix_processed) = NULL

# Saving as TSV in TXT file
fwrite(sampled_signature_matrix_processed,
       file = "ref_count_matrix_sampled.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

rm(sampled_signature_matrix_processed)
gc()

################### LSSMS Large Sampled Signature Matrix Simplif ###################
# We can pick 10 males and 10 females
# Everything else -> Goes to validation
# We will set fraction to 0.5
GSE144136_GSE213982_LSSMS = GSE144136_GSE213982_metadata
GSE144136_GSE213982_LSSMS$is_validation = NULL

female_vector = GSE144136_GSE213982_LSSMS$person
female_vector = female_vector[stri_detect_fixed(female_vector, "F")]
female_vector = unique(female_vector)
set.seed(12345) # To try another seed
LSSMS_females = sample(female_vector, size=10, replace=FALSE)
set.seed(NULL)

male_vector = GSE144136_GSE213982_LSSMS$person
male_vector = male_vector[stri_detect_fixed(male_vector, "M")]
male_vector = unique(male_vector)
set.seed(12345) # To try another seed
LSSMS_males = sample(male_vector, size=10, replace=FALSE)
set.seed(NULL)
LSSMS_total = c(LSSMS_females, LSSMS_males)

LSSMS_total_pheno_1 = ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %in% LSSMS_total, ]
LSSMS_total_pheno_2 = ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %in% LSSMS_total, ]

# Inspection of cells VS participant
LSSMS_metadata = GSE144136_GSE213982_LSSMS[GSE144136_GSE213982_LSSMS$person %in% LSSMS_total,]
table(LSSMS_metadata$person, LSSMS_metadata$cell_overall)

# Solution: Ast, End, Mic, Mix, Oli, OPC -> TAKE all cells
# Solution: ExN, InN, -> Sample 6000 cells
tmp_ExN_indeces = which(LSSMS_metadata$cell_overall == "ExN")
set.seed(12345)
tmp_selected_ExN_indeces = sample(tmp_ExN_indeces, size = 6000, replace = FALSE)
set.seed(NULL)

tmp_InN_indeces = which(LSSMS_metadata$cell_overall == "InN")
set.seed(12345)
tmp_selected_InN_indeces = sample(tmp_InN_indeces, size = 6000, replace = FALSE)
set.seed(NULL)

LSSMS_metadata$is_sampled = NA
for (i in 1:nrow(LSSMS_metadata)){
  
  if (LSSMS_metadata$cell_overall[i] %!in% c("ExN", "InN")){
    LSSMS_metadata$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_ExN_indeces){
    LSSMS_metadata$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_InN_indeces){
    LSSMS_metadata$is_sampled[i] = TRUE
  } else {
    LSSMS_metadata$is_sampled[i] = FALSE
  }
}

table(LSSMS_metadata$is_sampled, LSSMS_metadata$cell_overall)

"         Ast   End   ExN   InN   Mic   Mix   Oli   OPC
  FALSE     0     0 19268  2245     0     0     0     0
  TRUE   3616  1107  6000  6000   827  1410  5285  2228
"
LSSMS_picked_cells = LSSMS_metadata[LSSMS_metadata$is_sampled, "sample_name"]

LSSMS_processed = GSE144136_GSE213982_count_matrix
LSSMS_processed = LSSMS_processed[,LSSMS_picked_cells]
dim(LSSMS_processed) # 36588 26473

LSSMS_processed = as.matrix(LSSMS_processed)
dim(LSSMS_processed) # 36588 26473
colnames(LSSMS_processed)[1:10]
rownames(LSSMS_processed)[1:10]

replacement_colnames = sapply(colnames(LSSMS_processed), function(x){
  
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[3]
  
  # Ast Mic Oli OPC -> Glia
  # ExN InN -> Neuronal
  # End Mix -> Other cells
  
  if (x %in% c("Ast", "Mic", "Oli", "OPC")){
    return("Glia")
  }
  
  if (x %in% c("ExN", "InN")){
    return("Neuronal")
  }
  
  return("Other")
})

table(replacement_colnames)
"
replacement_colnames
    Glia Neuronal    Other 
   11956    12000     2517 
"
colnames(LSSMS_processed) = replacement_colnames
LSSMS_processed = as.data.frame(LSSMS_processed)
genes_column = data.frame("Gene" = rownames(LSSMS_processed))
LSSMS_processed = cbind(genes_column, LSSMS_processed)
colnames(LSSMS_processed)[1:10]
LSSMS_processed$Gene[1:10]
rownames(LSSMS_processed) = NULL

# Saving as TSV in TXT file
fwrite(LSSMS_processed,
       file = "ref_count_matrix_LSSMS.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

rm(LSSMS_processed)
gc()

################### Pseudobulking for LSSMS ###################
GSE144136_GSE213982_seurat = CreateSeuratObject(counts = GSE144136_GSE213982_count_matrix)
GSE144136_GSE213982_metadata_init = GSE144136_GSE213982_metadata
GSE144136_GSE213982_metadata_init$is_validation = NULL
GSE144136_GSE213982_metadata_init$cell_simplif = sapply(GSE144136_GSE213982_metadata_init$cell_overall, function(x){
  if (x %in% c("Ast", "Mic", "Oli", "OPC")){
    return("Glia")
  }
  
  if (x %in% c("ExN", "InN")){
    return("Neuronal")
  }
  
  return("Other")
})
table(GSE144136_GSE213982_metadata_init$cell_simplif)
GSE144136_GSE213982_seurat = AddMetaData(object = GSE144136_GSE213982_seurat, metadata =  GSE144136_GSE213982_metadata_init)

GSE144136_GSE213982_cell_proprotions_simplif = table(GSE144136_GSE213982_metadata_init$person, GSE144136_GSE213982_metadata_init$cell_simplif)
GSE144136_GSE213982_cell_proprotions_simplif = as.data.frame.matrix(GSE144136_GSE213982_cell_proprotions_simplif)
GSE144136_GSE213982_cell_proprotions_simplif = apply(GSE144136_GSE213982_cell_proprotions_simplif, 1, function(x){
  x = x/sum(x)
  return(x)
}, simplify = FALSE)
GSE144136_GSE213982_cell_proprotions_simplif = do.call(rbind, GSE144136_GSE213982_cell_proprotions_simplif)
GSE144136_GSE213982_cell_proprotions_simplif = as.data.frame(GSE144136_GSE213982_cell_proprotions_simplif)
GSE144136_GSE213982_cell_proprotions_simplif$PARTICIPANT = rownames(GSE144136_GSE213982_cell_proprotions_simplif)

GSE144136_GSE213982_counts_simplif_cell = AggregateExpression(GSE144136_GSE213982_seurat, group.by = c("person", "cell_simplif"), return.seurat = FALSE)
GSE144136_GSE213982_counts_simplif_cell = GSE144136_GSE213982_counts_simplif_cell$RNA
GSE144136_GSE213982_counts_simplif_cell = as.matrix(GSE144136_GSE213982_counts_simplif_cell)

colnames(GSE144136_GSE213982_counts_simplif_cell)[1:10]

rm(GSE144136_GSE213982_seurat)
gc()

################### File cleaning ###################
rm(GSE144136_GSE213982_count_matrix_processed)
rm(dense_count_matrix_processed)
rm(dense_count_matrix_processed_df)
gc()

################### Build signature matrix with docker (fraction 0) ###################

# Making directories
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/output", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/output", recursive = TRUE)

# Moving reference expression files
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/input/ref_count_matrix_dense.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_dense_male.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/input/ref_count_matrix_dense_male.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


"Average gene expression threshold (in log2 space) for cells with the same identity/phenotype showing evidence of expression (default = 0.75).
Although appropriate for plate-based approaches (e.g., SmartSeq2), this threshold may be too high for single cell experiments generated using droplet-based platforms,
such as 10x Chromium or DropSeq, which generally capture a much smaller number of genes (e.g., <1,500). For the latter case, we recommend reducing this parameter to 0.50 
or even 0. Otherwise, the sparsity of the data may yield too few genes for creating a reliable signature matrix."


"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense.txt \
--fraction 0

"

"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense_male.txt \
--fraction 0

"

################### Comments on cell-specif expression imputation ###################
"S-mode (single cell mode) batch correction is tailored for single cell-derived signature matrices generated from droplet-based or UMI-based platforms, including 10x Chromium."
# We need to select S-mode
# Disable quantile normalization (disabling is recommended for RNA-Seq data)

# First time
"
docker pull cibersortx/hires
"
# Signature matrix is created from Chromium kits -> droplet-based -> S-mode

################### Rinning high-res of cell-specif expression (docker) ###################

# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample/output_files", recursive = TRUE)

full_valid_bulk_counts = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
rownames(full_valid_bulk_counts) = full_valid_bulk_counts$X
full_valid_bulk_counts$X = NULL
full_valid_bulk_counts = full_valid_bulk_counts[,colnames(full_valid_bulk_counts) %in% selected_validation_samples]
genes_column = data.frame("Gene" = rownames(full_valid_bulk_counts))
full_valid_bulk_counts = cbind(genes_column, full_valid_bulk_counts)
rownames(full_valid_bulk_counts) = NULL

fwrite(full_valid_bulk_counts,
            file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files/validation_pseudobulked_counts_full.txt", 
            sep = "\t",
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)


validation_cell_specif_counts = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_cell_type.csv")
rownames(validation_cell_specif_counts) = validation_cell_specif_counts$X
validation_cell_specif_counts$X = NULL
validation_cell_specif_counts_selection_index = sapply(colnames(validation_cell_specif_counts), function(x){
  participant = unlist(stri_split_fixed(x, pattern = "_"))[1]
  if (participant %in% selected_validation_samples){
    return(TRUE)
  }
  return(FALSE)
})
validation_cell_specif_counts = validation_cell_specif_counts[,validation_cell_specif_counts_selection_index]

# Place files to the directory
genes_column = data.frame("Gene" = rownames(validation_cell_specif_counts))
validation_cell_specif_counts = cbind(genes_column, validation_cell_specif_counts)
rownames(validation_cell_specif_counts) = NULL
fwrite(validation_cell_specif_counts,
            file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files/validation_cell_specif_counts.txt", 
            sep = "\t",
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  \
--threads 8
'



combine_pngs <- function(input_dir,
                         output_file,
                         ncol = NULL,
                         nrow = NULL,
                         bg = "white",
                         border_px = 0,
                         resize_to = NULL) {
  stopifnot(dir.exists(input_dir))
  
  # 1) Collect PNGs
  files <- list.files(input_dir, pattern = "(?i)\\.png$", full.names = TRUE)
  if (!length(files)) stop("No PNG files found in ", input_dir)
  
  # 2) Read
  imgs <- image_read(files)
  
  # 3) Grid layout
  n <- length(imgs)
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n))
    nrow <- ceiling(n / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n / ncol)
  }
  
  # 4) Tile size
  if (is.null(resize_to)) {
    info <- image_info(imgs)
    tgt_w <- min(info$width)
    tgt_h <- min(info$height)
  } else {
    tgt_w <- resize_to[1]; tgt_h <- resize_to[2]
  }
  resize_geom <- sprintf("%dx%d!", tgt_w, tgt_h)
  
  # 5) Resize + montage
  imgs_resized <- image_resize(imgs, resize_geom)
  tile_spec <- sprintf("%dx%d", ncol, nrow)
  gap_geom  <- sprintf("+%d+%d", border_px, border_px)
  
  montage <- image_montage(
    imgs_resized,
    tile     = tile_spec,
    geometry = gap_geom,
    bg       = bg
  )
  montage <- image_background(montage, bg, flatten = TRUE)
  
  # 6) Write with 300 DPI metadata
  image_write(montage, path = output_file, format = "png", density = "300x300")
  
  invisible(montage)
}

################### Evaluating fractions (docker) ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)
validation_cell_full_ordered = validation_cell_full
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]

plot_cell_prop_correlation = function(estimated_proportions,
                                      real_proportions,
                                      cell_type,
                                      path_to_save,
                                      title_prefix,
                                      color_dots){
  
  
  # Run correlations
  cor_spearman = cor.test(estimated_proportions, real_proportions, method = "spearman")$estimate
  cor_pearson = cor.test(estimated_proportions, real_proportions, method = "pearson")$estimate
  cor_spearman = round(cor_spearman, digits = 2)
  cor_pearson = round(cor_pearson, digits = 2)
  
  caption = paste0("R (Pearson): ", cor_pearson, "\n",
                   "R (Spearman): ", cor_spearman)
  
  tmp_df = data.frame(estimated_proportions = estimated_proportions, real_proportions = real_proportions)
  
  # Generate plot
  meta_cor_plot = ggplot(data = tmp_df, aes(x = real_proportions, y = estimated_proportions)) +
    geom_point(alpha = .5, col = color_dots, size = 1.5) +
    geom_smooth(method = "lm") + 
    labs(title=paste0(title_prefix, cell_type),
         x = "Real proportions", 
         y = "Estimated proportions",
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


dir.create("cell_type_imputation_performance/imputation_correlations_docker")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "10K cells (52 donors), 0 frac: ",
                             color_dots = "blue")
  
})
combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker/combined_img.png",
  ncol=3,
  border_px = 10
)

################### Evaluating high-res expression (docker) ###################

# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample/output_files"
val_out_files = list.files(output_folder, recursive = TRUE, full.names = TRUE)
val_out_files = val_out_files[stri_detect_fixed(val_out_files, pattern = "CIBERSORTxHiRes")]
val_out_files = val_out_files[!stri_detect_fixed(val_out_files, pattern = "Heatmap")]
loaded_imputed_val_counts = lapply(val_out_files, function(x){
  x = fread(x)
  x = as.data.frame(x)
  rownames(x) = x$GeneSymbol
  x$GeneSymbol = NULL
  return(x)
})

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window10.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

plot_gene_expression_imput_corr = function(cell_type,
                                           loaded_imputed_list,
                                           real_counts_df,
                                           selected_validation_samples,
                                           path_to_save,
                                           title_prefix,
                                           color_dots
                                           ){
  
  if ("Gene" %in% colnames(real_counts_df)){
    rownames(real_counts_df) = real_counts_df$Gene
  }
  
  loaded_imputed_list = loaded_imputed_list[[cell_type]]
  loaded_imputed_list = loaded_imputed_list[, selected_validation_samples]
  
  real_counts_df = real_counts_df[,stri_detect_fixed(colnames(real_counts_df), pattern = cell_type)]
  colnames(real_counts_df) = stri_replace_all_fixed(colnames(real_counts_df), pattern = paste0("_", cell_type), replacement = "")
  real_counts_df = real_counts_df[, selected_validation_samples]
  
  # compare names
  names_df = data.frame(init_names = rownames(real_counts_df), new_names = rownames(loaded_imputed_list))
  names_df$matching = names_df$init_names == names_df$new_names
  print(table(names_df$matching))
  
  rownames(loaded_imputed_list) = rownames(real_counts_df)
  gene_selector = apply(loaded_imputed_list, 1, function(x) any(is.na(x)))
  loaded_imputed_list = loaded_imputed_list[!gene_selector,]
  real_counts_df = real_counts_df[rownames(loaded_imputed_list),]
  
  # Sample 5000 genes
  gene_sample = sample(1:nrow(real_counts_df), size = 5000, replace = FALSE)
  sampled_loaded_imputed_list = loaded_imputed_list[gene_sample,]
  sampled_real_counts_df = real_counts_df[gene_sample,]
  
  imputed_vector = as.vector(t(sampled_loaded_imputed_list))
  real_vector = as.vector(t(sampled_real_counts_df))
  
  tmp_df = data.frame(imputed_vector = imputed_vector, real_vector = real_vector)
  
  # Run correlations
  cor_spearman = cor.test(imputed_vector, real_vector, method = "spearman")$estimate
  cor_pearson = cor.test(imputed_vector, real_vector, method = "pearson")$estimate
  cor_spearman = round(cor_spearman, digits = 2)
  cor_pearson = round(cor_pearson, digits = 2)
  
  caption = paste0("R (Pearson): ", cor_pearson, "\n",
                   "R (Spearman): ", cor_spearman)
  
  
  # Generate plot
  cell_cor_plot = ggplot(data = tmp_df, aes(x = real_vector, y = imputed_vector)) +
    geom_point(alpha = .5, col = color_dots, size = 1.5) +
    geom_smooth(method = "lm") + 
    labs(title=paste0(title_prefix, cell_type),
         x = "Real counts", 
         y = "Estimated counts",
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
  
  ggsave(file = path_to_save, plot = cell_cor_plot, width=3000, height=2560, units = "px", scale = 1)
  
}


dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker/", x, ".png")
  plot_gene_expression_imput_corr(cell_type = x,
                             loaded_imputed_list = loaded_imputed_val_counts,
                             real_counts_df = validation_cell_specif_counts,
                             selected_validation_samples = selected_validation_samples,
                             path_to_save = path,
                             title_prefix = "10K cells (52 donors), 0 frac: ",
                             color_dots = "blue")
  
})

# Ast failed to be imputed !
# End failed to be imputed !
# ExN failed to be imputed ! (Duplicated values across genes)
# InN failed to be imputed !
# Mic failed to be imputed !
# Mix failed to be imputed ! (Duplicated values across genes)
# Oli failed to be imputed !
# OPC failed to be imputed ! (Duplicated values across genes)

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker/combined_img.png",
  ncol=3,
  border_px = 10
)

################### MO Rinning high-res of cell-specif expression (docker) ###################
# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/output_files", recursive = TRUE)

# Place files to the directory
fwrite(full_valid_bulk_counts,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/output/CIBERSORTx_ref_count_matrix_dense_male_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_male_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/input/ref_count_matrix_dense_male.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  \
--threads 8

'

################### MO Evaluating fractions (docker) ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)
validation_cell_full_ordered = validation_cell_full
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]


dir.create("cell_type_imputation_performance/imputation_correlations_docker_MO")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_MO/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~12K cells (3 male donors), 0 frac: ",
                             color_dots = "red")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_MO",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_MO/combined_img.png",
  ncol=3,
  border_px = 10
)

################### MO Evaluating high-res expression (docker) ###################

# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/output_files"
val_out_files = list.files(output_folder, recursive = TRUE, full.names = TRUE)
val_out_files = val_out_files[stri_detect_fixed(val_out_files, pattern = "CIBERSORTxHiRes")]
val_out_files = val_out_files[!stri_detect_fixed(val_out_files, pattern = "Heatmap")]
loaded_imputed_val_counts = lapply(val_out_files, function(x){
  x = fread(x)
  x = as.data.frame(x)
  rownames(x) = x$GeneSymbol
  x$GeneSymbol = NULL
  return(x)
})

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window10.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed


dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_MO")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_MO/", x, ".png")
  plot_gene_expression_imput_corr(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = validation_cell_specif_counts,
                                  selected_validation_samples = selected_validation_samples,
                                  path_to_save = path,
                                  title_prefix = "~12K cells (3 male donors), 0 frac: ",
                                  color_dots = "red")
  
})


combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_MO",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_MO/combined_img.png",
  ncol=3,
  border_px = 10
)

################### Large valid MO High-res running ###################
# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/output_files", recursive = TRUE)

# Place files to the directory
source_samples = c("M20", "M11", "M15")
full_valid_bulk_counts_large = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
rownames(full_valid_bulk_counts_large) = full_valid_bulk_counts_large$X
full_valid_bulk_counts_large$X = NULL
full_valid_bulk_counts_large = full_valid_bulk_counts_large[,colnames(full_valid_bulk_counts_large) %!in% source_samples]
genes_column = data.frame("Gene" = rownames(full_valid_bulk_counts_large))
full_valid_bulk_counts_large = cbind(genes_column, full_valid_bulk_counts_large)
rownames(full_valid_bulk_counts_large) = NULL

fwrite(full_valid_bulk_counts_large,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/output/CIBERSORTx_ref_count_matrix_dense_male_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_male_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_MO/input/ref_count_matrix_dense_male.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt

'


################### Large valid MO Evaluating fractions (docker) ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)


# Exploring reference fractions
source_samples = c("M20", "M11", "M15")
tmp_cell_1 = ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %!in% source_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
tmp_cell_2 = ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %!in% source_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
tmp_cell_full = rbind(tmp_cell_1, tmp_cell_2)

validation_cell_full_ordered = tmp_cell_full
imputed_fractions_full = imputed_fractions_full[imputed_fractions_full$Mixture %in% validation_cell_full_ordered$PARTICIPANT,]
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]


dir.create("cell_type_imputation_performance/imputation_correlations_docker_MO_LV")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_MO_LV/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~12K cells (3 male donors), 0 frac, large validation: ",
                             color_dots = "green")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_MO_LV",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_MO_LV/combined_img.png",
  ncol=3,
  border_px = 10
)

################### Large valid MO Evaluating high-res expression (docker) ###################

# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/output_files"
val_out_files = list.files(output_folder, recursive = TRUE, full.names = TRUE)
val_out_files = val_out_files[stri_detect_fixed(val_out_files, pattern = "CIBERSORTxHiRes")]
val_out_files = val_out_files[!stri_detect_fixed(val_out_files, pattern = "Heatmap")]
loaded_imputed_val_counts = lapply(val_out_files, function(x){
  x = fread(x)
  x = as.data.frame(x)
  rownames(x) = x$GeneSymbol
  x$GeneSymbol = NULL
  return(x)
})

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_MO_large_valid/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window32.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

# select appropriate counts
source_samples = c("M20", "M11", "M15")
tmp_counts = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_cell_type.csv")
rownames(tmp_counts) = tmp_counts$X
tmp_counts$X = NULL
tmp_counts_selection_index = sapply(colnames(tmp_counts), function(x){
  participant = unlist(stri_split_fixed(x, pattern = "_"))[1]
  if (participant %!in% source_samples){
    return(TRUE)
  }
  return(FALSE)
})
tmp_counts = tmp_counts[,tmp_counts_selection_index]

required_tmp_samples = colnames(tmp_counts)
required_tmp_samples = sapply(required_tmp_samples, function(x){
  x = stri_split_fixed(x, pattern = "_")
  x = unlist(x)
  x = x[1]
  return(x)
})
required_tmp_samples = unique(required_tmp_samples)

dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_MO_LV")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_MO_LV/", x, ".png")
  plot_gene_expression_imput_corr(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = tmp_counts,
                                  selected_validation_samples = required_tmp_samples,
                                  path_to_save = path,
                                  title_prefix = "~12K cells (3 male donors), 0 frac, large validation: ",
                                  color_dots = "green")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_MO_LV",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_MO_LV/combined_img.png",
  ncol=3,
  border_px = 10
)

################### Build signature matrix with docker (30K, 0.5 fraction) ###################

# Making directories
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/output", recursive = TRUE)

# Moving reference expression files
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_dense_30K.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/input/ref_count_matrix_dense.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense.txt \
--fraction 0.5

"

################### 30K, 0.5 Rinning high-res of cell-specif expression (docker) ###################
# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/output_files", recursive = TRUE)

# Place files to the directory
fwrite(full_valid_bulk_counts,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_30K/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  \
--threads 8

'
################### 30K, 0.5 Evaluating fractions (docker) ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)
validation_cell_full_ordered = validation_cell_full
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]


dir.create("cell_type_imputation_performance/imputation_correlations_docker_30K")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_30K/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "30K cells (52 donors), 0.5 frac: ",
                             color_dots = "purple")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_30K",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_30K/combined_img.png",
  ncol=3,
  border_px = 10
)

################### 30K, 0.5 Evaluating high-res expression (docker) ###################

# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/output_files"
val_out_files = list.files(output_folder, recursive = TRUE, full.names = TRUE)
val_out_files = val_out_files[stri_detect_fixed(val_out_files, pattern = "CIBERSORTxHiRes")]
val_out_files = val_out_files[!stri_detect_fixed(val_out_files, pattern = "Heatmap")]
loaded_imputed_val_counts = lapply(val_out_files, function(x){
  x = fread(x)
  x = as.data.frame(x)
  rownames(x) = x$GeneSymbol
  x$GeneSymbol = NULL
  return(x)
})

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_30K/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window10.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed


dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_30K")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_30K/", x, ".png")
  plot_gene_expression_imput_corr(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = validation_cell_specif_counts,
                                  selected_validation_samples = selected_validation_samples,
                                  path_to_save = path,
                                  title_prefix = "30K cells (52 donors), 0.5 frac: ",
                                  color_dots = "purple")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_30K",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_30K/combined_img.png",
  ncol=3,
  border_px = 10
)




################### Comments on cell-specif expression imputation ###################
# None of the suggested matrix strategies imputed cell-type specific expression in high resolution mode with good accuracy
# Cell proportion estimations are relatively accurate
# First -> Evaluate on large validation using 3 donor matrix -> OK expression for ExN
# It potentially makes sense to run cell-prop-adjusted analysis
# It potentially makes sense to reduce cell-type specificity to major cell types
# It potentially makes sense to run analysis in group mode per class -> compare cases vs controls



################### Build sampled signature matrix  ###################

# Making directories
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/output", recursive = TRUE)

# Moving reference expression files
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_sampled.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/input/ref_count_matrix_dense.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense.txt \
--fraction 0
"

################### Sampled signature running high res  ###################

# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/output_files", recursive = TRUE)

# Place files to the directory
source_samples = selected_signature_total
full_valid_bulk_counts_ssig = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
rownames(full_valid_bulk_counts_ssig) = full_valid_bulk_counts_ssig$X
full_valid_bulk_counts_ssig$X = NULL
full_valid_bulk_counts_ssig = full_valid_bulk_counts_ssig[,colnames(full_valid_bulk_counts_ssig) %!in% source_samples]
genes_column = data.frame("Gene" = rownames(full_valid_bulk_counts_ssig))
full_valid_bulk_counts_ssig = cbind(genes_column, full_valid_bulk_counts_ssig)
rownames(full_valid_bulk_counts_ssig) = NULL

fwrite(full_valid_bulk_counts_ssig,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  
'

################### Sampled signature Evaluating fractions (docker) ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)


# Exploring reference fractions
source_samples = selected_signature_total
tmp_cell_1 = ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %!in% source_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
tmp_cell_2 = ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %!in% source_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
tmp_cell_full = rbind(tmp_cell_1, tmp_cell_2)

validation_cell_full_ordered = tmp_cell_full
imputed_fractions_full = imputed_fractions_full[imputed_fractions_full$Mixture %in% validation_cell_full_ordered$PARTICIPANT,]
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]


dir.create("cell_type_imputation_performance/imputation_correlations_docker_sSig")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_sSig/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~8K balanced cells (10 donors), 0 frac: ",
                             color_dots = "orange")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_sSig",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_sSig/combined_img.png",
  ncol=3,
  border_px = 10
)


################### Sampled signature Evaluating High resolution (docker) ##################
# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/output_files"
val_out_files = list.files(output_folder, recursive = TRUE, full.names = TRUE)
val_out_files = val_out_files[stri_detect_fixed(val_out_files, pattern = "CIBERSORTxHiRes")]
val_out_files = val_out_files[!stri_detect_fixed(val_out_files, pattern = "Heatmap")]
loaded_imputed_val_counts = lapply(val_out_files, function(x){
  x = fread(x)
  x = as.data.frame(x)
  rownames(x) = x$GeneSymbol
  x$GeneSymbol = NULL
  return(x)
})

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window31.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

# select appropriate counts
source_samples = selected_signature_total
tmp_counts = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_cell_type.csv")
rownames(tmp_counts) = tmp_counts$X
tmp_counts$X = NULL
tmp_counts_selection_index = sapply(colnames(tmp_counts), function(x){
  participant = unlist(stri_split_fixed(x, pattern = "_"))[1]
  if (participant %!in% source_samples){
    return(TRUE)
  }
  return(FALSE)
})
tmp_counts = tmp_counts[,tmp_counts_selection_index]

required_tmp_samples = colnames(tmp_counts)
required_tmp_samples = sapply(required_tmp_samples, function(x){
  x = stri_split_fixed(x, pattern = "_")
  x = unlist(x)
  x = x[1]
  return(x)
})
required_tmp_samples = unique(required_tmp_samples)

dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_sSig")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_sSig/", x, ".png")
  plot_gene_expression_imput_corr(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = tmp_counts,
                                  selected_validation_samples = required_tmp_samples,
                                  path_to_save = path,
                                  title_prefix = "~8K balanced cells (10 donors), 0 frac: ",
                                  color_dots = "orange")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_sSig",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_sSig/combined_img.png",
  ncol=3,
  border_px = 10
)


################### Build LSSMS  ###################

# Making directories
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/output", recursive = TRUE)

# Moving reference expression files
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_LSSMS.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/input/ref_count_matrix_dense.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense.txt \
--fraction 0.5
"

################### Run LSSMS  ###################

# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/output_files", recursive = TRUE)

# Place files to the directory
source_samples = LSSMS_total
full_valid_bulk_counts_LSSMS = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
rownames(full_valid_bulk_counts_LSSMS) = full_valid_bulk_counts_LSSMS$X
full_valid_bulk_counts_LSSMS$X = NULL
full_valid_bulk_counts_LSSMS = full_valid_bulk_counts_LSSMS[,colnames(full_valid_bulk_counts_LSSMS) %!in% source_samples]
genes_column = data.frame("Gene" = rownames(full_valid_bulk_counts_LSSMS))
full_valid_bulk_counts_LSSMS = cbind(genes_column, full_valid_bulk_counts_LSSMS)
rownames(full_valid_bulk_counts_LSSMS) = NULL

fwrite(full_valid_bulk_counts_LSSMS,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  
'

################### Evaluate fractions LSSMS  ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)


# Exploring reference fractions
source_samples = LSSMS_total
tmp_cell_full = GSE144136_GSE213982_cell_proprotions_simplif[GSE144136_GSE213982_cell_proprotions_simplif$PARTICIPANT %!in% source_samples, ]

validation_cell_full_ordered = tmp_cell_full
imputed_fractions_full = imputed_fractions_full[imputed_fractions_full$Mixture %in% validation_cell_full_ordered$PARTICIPANT,]
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]

all(imputed_fractions_full$Mixture == validation_cell_full_ordered$PARTICIPANT) # TRUE

dir.create("cell_type_imputation_performance/imputation_correlations_docker_LSSMS")
cell_types = c("Glia", "Neuronal", "Other")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_LSSMS/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~26K LSSMS (20 donors), 0 frac: ",
                             color_dots = "red")
  
})


combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_LSSMS",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_LSSMS/combined_img.png",
  ncol=3,
  border_px = 10
)

################### Evaluate high resolution LSSMS  ###################

# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/output_files"
val_out_files = list.files(output_folder, recursive = TRUE, full.names = TRUE)
val_out_files = val_out_files[stri_detect_fixed(val_out_files, pattern = "CIBERSORTxHiRes")]
val_out_files = val_out_files[!stri_detect_fixed(val_out_files, pattern = "Heatmap")]
loaded_imputed_val_counts = lapply(val_out_files, function(x){
  x = fread(x)
  x = as.data.frame(x)
  rownames(x) = x$GeneSymbol
  x$GeneSymbol = NULL
  return(x)
})

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window12.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

# Select appropriate pseudobulked counts 
source_samples = LSSMS_total
tmp_counts = GSE144136_GSE213982_counts_simplif_cell
rownames(tmp_counts)[1:10]
colnames(tmp_counts)[1:10]

tmp_counts_selection_index = sapply(colnames(tmp_counts), function(x){
  participant = unlist(stri_split_fixed(x, pattern = "_"))[1]
  if (participant %!in% source_samples){
    return(TRUE)
  }
  return(FALSE)
})
tmp_counts = tmp_counts[,tmp_counts_selection_index]

required_tmp_samples = colnames(tmp_counts)
required_tmp_samples = sapply(required_tmp_samples, function(x){
  x = stri_split_fixed(x, pattern = "_")
  x = unlist(x)
  x = x[1]
  return(x)
})
required_tmp_samples = unique(required_tmp_samples)
required_tmp_samples = required_tmp_samples[required_tmp_samples %in% colnames(loaded_imputed_val_counts[[1]])]

dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS")
cell_types = c("Glia", "Neuronal", "Other")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS/", x, ".png")
  plot_gene_expression_imput_corr(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = tmp_counts,
                                  selected_validation_samples = required_tmp_samples,
                                  path_to_save = path,
                                  title_prefix ="~26K LSSMS (20 donors), 0 frac: ",
                                  color_dots = "red")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS/combined_img.png",
  ncol=3,
  border_px = 10
)


# Move all files in one place for combining
merged_figures = list.files("cell_type_imputation_performance", full.names = TRUE, recursive = TRUE, pattern = "combined_img")
merged_figures = merged_figures[stri_detect_fixed(merged_figures, pattern = "highres")]

for (i in 1:length(merged_figures)){
  
  names_split = unlist(stri_split_fixed(merged_figures[i], "/"))
  sel_name = names_split[2]
  sel_name = stri_replace_all_fixed(sel_name, pattern = "imputation_correlations_highres_docker", replacement = "")
  
  if (sel_name == ""){
    sel_name = "10K_unbalanced_signature"
  }
  if (sel_name == "_30K"){
    sel_name = "30K_unbalanced_signature"
  }
  if (sel_name == "_LSSMS"){
    sel_name = "LSSMS"
  }
  if (sel_name == "_MO_LV"){
    sel_name = "3males_large_val"
  }
  if (sel_name == "_MO"){
    sel_name = "3males_orig_val"
  }
  if (sel_name == "_sSig"){
    sel_name = "sampled_signature"
  }
  
  # rewrite file
  file_path = paste0("cell_type_imputation_performance/Merged_files_highres_performance/", sel_name, ".png")
  
  file.copy(merged_figures[i],
            file_path,
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
  
}