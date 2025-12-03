
# Comments
# This code is intended to be viewed in RStudio and contains appropriate headings
# "=" symbol was used as an assignment operator
# Code for CIBERSORTx runs is presented as character strings " <command code> " within the file
# Real token for CIBERSORTx docker container is replaced with <token_from_cibersortx_website> in associted character strings




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
library(edgeR)

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


plot_gene_expression_imput_corr_harmon = function(cell_type,
                                           loaded_imputed_list,
                                           real_counts_df,
                                           selected_validation_samples,
                                           path_to_save,
                                           title_prefix,
                                           color_dots
){
  
  
  
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


################### Comments on cell-specif expression imputation ###################
# None of the suggested matrix strategies imputed cell-type specific expression in high resolution mode with good accuracy
# Bulking approximated expression of glia of around 60% (pearson) based on validation samples
# Neuronal correlation was around 0.85% 
# Solutions: 
# 1) Run high resolution and only use sampled signature -> ExN
# 2) Run high resolution and use LSSMS -> Neuronal and Glia

#### We need to prepare matrix in a harmonized way!


################### Overview of cohorts ###################
# GSE102556 11+37 = 48 samples
# GSE243356 32+29 = 61 samples
# GSE248260 9+15 = 24 samples
# GSE202537 60+7 = 67 samples
# GSE101521 38+21 = 59 samples
# GSE144136 15+17 = 32 samples (already single cell)
# GSE213982 18+17= 36 samples (already single cell)


# -> GSE248260 is the most risky cohort due to low sample count!
# Need to check if everything works on this cohort


################### Gene name inspection ###################
# We also need to harmonize gene symbols with signature matrix and reference profile
TMP_sig_matrix_gene_rows_fixed = check_gene_symbol_NIH(PRF_gene_symbols=TMP_sig_matrix_gene_rows, 
                                                       PRF_ref_NIH_expanded=Homo_Sapiens_Gene_info_NIH_expanded, 
                                                       PRF_replace_NA_with_old=TRUE)

table(TMP_sig_matrix_gene_rows_fixed$Approved)



"
FALSE  TRUE 
13090 23498 

"
# Select any processed validation dataset
cibersortx_names = smart_fread("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS/output_files/CIBERSORTxHiRes_NA_Glia_Window12.txt")
cibersortx_names = cibersortx_names$GeneSymbol

# 
cibersortx_names_df = data.frame(intial_gene_names = TMP_sig_matrix_gene_rows, cibersortx_names=cibersortx_names)
table(cibersortx_names_df$intial_gene_names == cibersortx_names_df$cibersortx_names) # 2150 genes are not matching
all(TMP_sig_matrix_gene_rows_fixed$x == cibersortx_names_df$intial_gene_names) #TRUE

name_harmoniz_df = cbind(TMP_sig_matrix_gene_rows_fixed, cibersortx_names = cibersortx_names_df$cibersortx_names)
name_harmoniz_df$cybersort_change = name_harmoniz_df$x != name_harmoniz_df$cibersortx_names
table(name_harmoniz_df$cybersort_change)

# Inspect signature matrices
LSSMS_signature = smart_fread("CIBERSORTx_main/singature_matrix_build_LSSMS/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt")
LSSMS_signature$NAME[LSSMS_signature$NAME %in% name_harmoniz_df[name_harmoniz_df$cybersort_change, "cibersortx_names"]]

# QUITE MANY NAMES ARE MODIFIED -> WE NEED TO ASSEMBLE 2 last matrices using corrected gene symbols or fix gene symbols or fix matrices
# Easiest is to assemble new matrices




################### Base SC File preparation ###################

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

################### Harmonize SC data ###################
TMP_sig_matrix_gene_rows_gene_check = check_gene_symbol_NIH(PRF_gene_symbols=TMP_sig_matrix_gene_rows, 
                                                            PRF_ref_NIH_expanded=Homo_Sapiens_Gene_info_NIH_expanded, 
                                                            PRF_replace_NA_with_old=TRUE)


STAR_gene_names = smart_fread("Data_preprocessing_analysis/RNAseq_gene_probes.csv", header = TRUE)

STAR_gene_names_check = check_gene_symbol_NIH(PRF_gene_symbols=STAR_gene_names$Gene_symbol, 
                                              PRF_ref_NIH_expanded=Homo_Sapiens_Gene_info_NIH_expanded, 
                                              PRF_replace_NA_with_old=TRUE)


# Finding duplicates
duplicates_sc = TMP_sig_matrix_gene_rows_gene_check$Suggested.Symbol[duplicated(TMP_sig_matrix_gene_rows_gene_check$Suggested.Symbol)]
duplicates_star = STAR_gene_names_check$Suggested.Symbol[duplicated(STAR_gene_names_check$Suggested.Symbol)]
duplicted_all = unlist(c(duplicates_sc, duplicates_star))
duplicted_all = unique(duplicted_all)

# Finding intersect
intersecting_symbols = intersect(TMP_sig_matrix_gene_rows_gene_check$Suggested.Symbol, STAR_gene_names_check$Suggested.Symbol)
intersecting_symbols = unique(intersecting_symbols)

all(intersecting_symbols %in% TMP_sig_matrix_gene_rows_gene_check$Suggested.Symbol) # TRUE
all(intersecting_symbols %in% STAR_gene_names_check$Suggested.Symbol) # TRUE

# Labelling indeces
TMP_sig_matrix_gene_rows_gene_check$is_good_symbol = sapply(TMP_sig_matrix_gene_rows_gene_check$Suggested.Symbol, function(x){
  
  if (x %in% duplicted_all){
    return(FALSE)
  }
  
  if (x %!in% intersecting_symbols){
    return(FALSE)
  }
  
  return(TRUE)
  
})
table(TMP_sig_matrix_gene_rows_gene_check$is_good_symbol)

"
FALSE  TRUE 
13786 22802 
"

GSE144136_GSE213982_count_matrix_harmonized = GSE144136_GSE213982_count_matrix
all(rownames(GSE144136_GSE213982_count_matrix_harmonized) == TMP_sig_matrix_gene_rows_gene_check$x) # TRUE
GSE144136_GSE213982_count_matrix_harmonized = GSE144136_GSE213982_count_matrix_harmonized[TMP_sig_matrix_gene_rows_gene_check$is_good_symbol,]
dim(GSE144136_GSE213982_count_matrix_harmonized) # 22802 160711

Harmonized_matrix_gene_check = TMP_sig_matrix_gene_rows_gene_check[TMP_sig_matrix_gene_rows_gene_check$x %in% rownames(GSE144136_GSE213982_count_matrix_harmonized),]
all(Harmonized_matrix_gene_check$x == rownames(GSE144136_GSE213982_count_matrix_harmonized)) # TRUE
table(Harmonized_matrix_gene_check$is_good_symbol)

# Make adjustments for CIBERSORTx
# CIBERSORTx likes to change gene symbols so we have to fix it
Harmonized_matrix_gene_check$name_for_cibersort = stri_replace_all_fixed(Harmonized_matrix_gene_check$Suggested.Symbol, 
                                                                         pattern = " ",
                                                                         replacement = "_")
Harmonized_matrix_gene_check$name_for_cibersort = stri_replace_all_fixed(Harmonized_matrix_gene_check$name_for_cibersort, 
                                                                         pattern = "-",
                                                                         replacement = "__")
Harmonized_matrix_gene_check$name_for_cibersort = toupper(Harmonized_matrix_gene_check$name_for_cibersort)
any(duplicated(Harmonized_matrix_gene_check$name_for_cibersort)) # FALSE

sum(stri_count_fixed(Harmonized_matrix_gene_check$name_for_cibersort, pattern = "-")) # 0
sum(stri_count_fixed(Harmonized_matrix_gene_check$name_for_cibersort, pattern = " ")) # 0
sum(stri_count_fixed(Harmonized_matrix_gene_check$name_for_cibersort, pattern = ";")) # 0
sum(stri_count_fixed(Harmonized_matrix_gene_check$name_for_cibersort, pattern = ":")) # 0
sum(stri_count_fixed(Harmonized_matrix_gene_check$name_for_cibersort, pattern = "#")) # 0
sum(stri_count_fixed(Harmonized_matrix_gene_check$name_for_cibersort, pattern = "@")) # 0

# Saving names
all(rownames(GSE144136_GSE213982_count_matrix_harmonized) == Harmonized_matrix_gene_check$x) # TRUE
rownames(GSE144136_GSE213982_count_matrix_harmonized) = Harmonized_matrix_gene_check$name_for_cibersort
rownames(GSE144136_GSE213982_count_matrix_harmonized)[1:30]

################### Normalize to CPM or TPM? ###################

# Human GRCh38 (GENCODE v32/Ensembl98)
# /home/aleksandr/Downloads/Homo_sapiens.GRCh38.98.gtf
"
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project
python gtftools.py -l gene_length.txt Homo_sapiens.GRCh38.98.gtf
python gtftools.py -s isoform.txt Homo_sapiens.GRCh38.98.gtf
"

# Processed files
GRCh38_lengths = smart_fread("gene_length.txt")
GRCh38_names = fread("isoform.txt")
table(GRCh38_lengths$gene %in% GRCh38_names$V6)
"TRUE 
60527"
GRCh38_lengths$symbol = sapply(GRCh38_lengths$gene, function(x){
  out = GRCh38_names[GRCh38_names$V6 == x, "V7"]
  out = unlist(out)
  out = unique(out)
  out = out[1]
  return(out)
})
GRCh38_lengths = dplyr::arrange(GRCh38_lengths, symbol)

# Comment
# Limma requires raw counts to work with RNAseq -> It's probably better to keep it in such way
# Cibersort paper:
"
CIBERSORTx will automatically normalize the input data such that the sum of all normalized reads is the same for each transcriptome. 
If a gene length-normalized expression matrix is provided (e.g., RPKM), then the signature matrix will be adjusted to TPM (transcripts per million).
If a count matrix is provided, the signature matrix will be normalized to CPM (counts per million).
"


################### Sampled signature FINAL ###################
# We can pick 5 males and 5 females
# Take approx 2000 cells of each type -> signature matrix
# Everything else -> Goes to validation
# We will set fraction to 0
GSE144136_GSE213982_metadata_sampled_matrix_final = GSE144136_GSE213982_metadata
GSE144136_GSE213982_metadata_sampled_matrix_final$is_validation = NULL

female_vector = GSE144136_GSE213982_metadata$person
female_vector = female_vector[stri_detect_fixed(female_vector, "F")]
female_vector = unique(female_vector)
set.seed(12345678) # # This seed gives several cases and controls
selected_signature_F_females = sample(female_vector, size=5, replace=FALSE)
set.seed(NULL)

# Inspect phenotypes:
ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %in% selected_signature_F_females, ] # 2 cases, 3 controls

male_vector = GSE144136_GSE213982_metadata$person
male_vector = male_vector[stri_detect_fixed(male_vector, "M")]
male_vector = unique(male_vector)
set.seed(1234567) # This seed gives several cases and controls
selected_signature_F_males = sample(male_vector, size=5, replace=FALSE)
set.seed(NULL)
# Inspect phenotypes:
ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %in% selected_signature_F_males, ] # 3 cases, 2 controls

# Final subset
selected_signature_F_total = c(selected_signature_F_females, selected_signature_F_males)

# Inspection of cells VS participant
sampled_signature_F_metadata = GSE144136_GSE213982_metadata_sampled_matrix_final[GSE144136_GSE213982_metadata_sampled_matrix_final$person %in% selected_signature_F_total,]
table(sampled_signature_F_metadata$person, sampled_signature_F_metadata$cell_overall)

"     Ast  End  ExN  InN  Mic  Mix  Oli  OPC
  F12  212   81 1007  265   49   20  245   69
  F25  271  155  996  358   93  123  515  245
  F33  374  228 1774  505  174   67  741  247
  F34  165  103 1046  328  125   29  536  151
  F37  257   98  132  111   22   23  186  125
  M1    16    9 1860  707    8   55   70  138
  M13  263   11 1118  413   10   83  178  106
  M21   74    9 1265  364   11   66  237  149
  M33  284   11  770  886   51    8  204  149
  M34   29    3 2165  783   38   18   95   61
"

# Solution: Ast, End, Mic, Mix, Oli, OPC -> TAKE all cells
# Solution: ExN, InN, -> Sample 2000 cells
tmp_ExN_indeces = which(sampled_signature_F_metadata$cell_overall == "ExN")
set.seed(12345)
tmp_selected_ExN_indeces = sample(tmp_ExN_indeces, size = 2000, replace = FALSE)
set.seed(NULL)

tmp_InN_indeces = which(sampled_signature_F_metadata$cell_overall == "InN")
set.seed(12345)
tmp_selected_InN_indeces = sample(tmp_InN_indeces, size = 2000, replace = FALSE)
set.seed(NULL)

sampled_signature_F_metadata$is_sampled = NA
for (i in 1:nrow(sampled_signature_F_metadata)){
  
  if (sampled_signature_F_metadata$cell_overall[i] %!in% c("ExN", "InN")){
    sampled_signature_F_metadata$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_ExN_indeces){
    sampled_signature_F_metadata$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_InN_indeces){
    sampled_signature_F_metadata$is_sampled[i] = TRUE
  } else {
    sampled_signature_F_metadata$is_sampled[i] = FALSE
  }
}

table(sampled_signature_F_metadata$is_sampled, sampled_signature_F_metadata$cell_overall)

"         Ast   End   ExN   InN   Mic   Mix   Oli   OPC
  FALSE     0     0 10133  2720     0     0     0     0
  TRUE   1945   708  2000  2000   581   492  3007  1440
"
sampled_signature_F_picked_cells = sampled_signature_F_metadata[sampled_signature_F_metadata$is_sampled, "sample_name"]

sampled_signature_matrix_processed_F = GSE144136_GSE213982_count_matrix_harmonized
sampled_signature_matrix_processed_F = sampled_signature_matrix_processed_F[,sampled_signature_F_picked_cells]
dim(sampled_signature_matrix_processed_F) # 22802 12173

sampled_signature_matrix_processed_F = as.matrix(sampled_signature_matrix_processed_F)
colnames(sampled_signature_matrix_processed_F)[1:10]
rownames(sampled_signature_matrix_processed_F)[1:10]

replacement_colnames = sapply(colnames(sampled_signature_matrix_processed_F), function(x){
  x = unlist(stri_split_fixed(x, pattern = "."))
  x = x[3]
  return(x)
})
table(replacement_colnames)

"replacement_colnames
 Ast  End  ExN  InN  Mic  Mix  Oli  OPC 
1945  708 2000 2000  581  492 3007 1440 
"
colnames(sampled_signature_matrix_processed_F) = replacement_colnames
sampled_signature_matrix_processed_F = as.data.frame(sampled_signature_matrix_processed_F)
genes_column = data.frame("Gene" = rownames(sampled_signature_matrix_processed_F))
sampled_signature_matrix_processed_F = cbind(genes_column, sampled_signature_matrix_processed_F)
colnames(sampled_signature_matrix_processed_F)[1:10]
rownames(sampled_signature_matrix_processed_F) = NULL

# Saving as TSV in TXT file
fwrite(sampled_signature_matrix_processed_F,
       file = "ref_count_matrix_sampled_harmonized.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

rm(sampled_signature_matrix_processed_F)
gc()


################### Build sampled signature matrix  ###################

# Making directories
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/output", recursive = TRUE)

# Moving reference expression files
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_sampled_harmonized.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/input/ref_count_matrix_dense.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense.txt \
--fraction 0
"

################### Prepare harmonized full validation counts ###################
full_valid_bulk_counts_harmon = read.csv("Data_preprocessing_analysis/GSE144136_results/GSE144136_GSE213982_counts_mixed.csv")
rownames(full_valid_bulk_counts_harmon) = full_valid_bulk_counts_harmon$X
full_valid_bulk_counts_harmon$X = NULL
genes_column = data.frame("Gene" = rownames(full_valid_bulk_counts_harmon))
full_valid_bulk_counts_harmon = cbind(genes_column, full_valid_bulk_counts_harmon)
rownames(full_valid_bulk_counts_harmon) = NULL
dim(full_valid_bulk_counts_harmon) #36588    73

# Haramonize symbols
full_valid_bulk_counts_harmon_gene_check = check_gene_symbol_NIH(PRF_gene_symbols=full_valid_bulk_counts_harmon$Gene, 
                                                               PRF_ref_NIH_expanded=Homo_Sapiens_Gene_info_NIH_expanded, 
                                                               PRF_replace_NA_with_old=TRUE)

full_valid_bulk_counts_harmon_gene_check$name_for_cibersort = stri_replace_all_fixed(full_valid_bulk_counts_harmon_gene_check$Suggested.Symbol, 
                                                                                   pattern = " ",
                                                                                   replacement = "_")
full_valid_bulk_counts_harmon_gene_check$name_for_cibersort = stri_replace_all_fixed(full_valid_bulk_counts_harmon_gene_check$name_for_cibersort, 
                                                                                   pattern = "-",
                                                                                   replacement = "__")
full_valid_bulk_counts_harmon_gene_check$name_for_cibersort = toupper(full_valid_bulk_counts_harmon_gene_check$name_for_cibersort)

full_valid_bulk_counts_harmon$Gene = full_valid_bulk_counts_harmon_gene_check$name_for_cibersort




################### Run sampled signature matrix harmonized in high resolution ###################
# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/output_files", recursive = TRUE)


# Select required subset
source_samples = selected_signature_F_total
full_valid_bulk_counts_ssig = full_valid_bulk_counts_harmon[,colnames(full_valid_bulk_counts_harmon) %!in% source_samples]


fwrite(full_valid_bulk_counts_ssig,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_sampled_harmon/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  
'

#

################### Evaluating fractions sSig ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)


# Exploring reference fractions
source_samples = selected_signature_F_total
tmp_cell_1 = ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %!in% source_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
tmp_cell_2 = ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %!in% source_samples, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
tmp_cell_full = rbind(tmp_cell_1, tmp_cell_2)

validation_cell_full_ordered = tmp_cell_full
imputed_fractions_full = imputed_fractions_full[imputed_fractions_full$Mixture %in% validation_cell_full_ordered$PARTICIPANT,]
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]


dir.create("cell_type_imputation_performance/imputation_correlations_docker_sSig_harmon")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_sSig_harmon/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~12K bal. cells harmon. genes (10 donors), 0 frac: ",
                             color_dots = "#8f1402")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_sSig_harmon",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_sSig_harmon/combined_img.png",
  ncol=3,
  border_px = 10
)

################### High resolution sSig evaluation ###################

# Loading imputed counts
output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/output_files"
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

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_sSig_harmon/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window31.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

# Check naming
table(rownames(loaded_imputed_val_counts[[1]]) == full_valid_bulk_counts_harmon_gene_check$name_for_cibersort)
# 21 names in validation are duplicated
"
FALSE  TRUE 
   21 36567 
"

# select appropriate counts
source_samples = selected_signature_F_total
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

all(rownames(tmp_counts) == TMP_sig_matrix_gene_rows_gene_check$x)
tmp_counts_names_df = TMP_sig_matrix_gene_rows_gene_check
tmp_counts_names_df$name_for_cibersort = stri_replace_all_fixed(tmp_counts_names_df$Suggested.Symbol, 
                                                                                     pattern = " ",
                                                                                     replacement = "_")
tmp_counts_names_df$name_for_cibersort = stri_replace_all_fixed(tmp_counts_names_df$name_for_cibersort, 
                                                                                     pattern = "-",
                                                                                     replacement = "__")
tmp_counts_names_df$name_for_cibersort = toupper(tmp_counts_names_df$name_for_cibersort)



dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_sSig_harmon")
cell_types = c("Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_sSig_harmon/", x, ".png")
  plot_gene_expression_imput_corr_harmon(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = tmp_counts,
                                  selected_validation_samples = required_tmp_samples,
                                  path_to_save = path,
                                  title_prefix = "~12K bal. cells harmon. genes (10 donors), 0 frac: ",
                                  color_dots = "#8f1402")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_sSig_harmon",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_sSig_harmon/combined_img.png",
  ncol=3,
  border_px = 10
)



################### LSSMS harmonized ###################
GSE144136_GSE213982_LSSMS_hr = GSE144136_GSE213982_metadata
GSE144136_GSE213982_LSSMS_hr$is_validation = NULL

female_vector = GSE144136_GSE213982_LSSMS_hr$person
female_vector = female_vector[stri_detect_fixed(female_vector, "F")]
female_vector = unique(female_vector)
set.seed(123456789) # This seed gives 5 cases and controls
LSSMS_females_hr = sample(female_vector, size=10, replace=FALSE)
set.seed(NULL)
# Inspect phenotypes:
ref_cell_pheno_GSE213982[ref_cell_pheno_GSE213982$X %in% LSSMS_females_hr, ] # 5 cases, 5 controls

male_vector = GSE144136_GSE213982_LSSMS_hr$person
male_vector = male_vector[stri_detect_fixed(male_vector, "M")]
male_vector = unique(male_vector)
set.seed(123456) # This seed gives 5 cases and controls
LSSMS_males_hr = sample(male_vector, size=10, replace=FALSE)
set.seed(NULL)
# Inspect phenotypes:
ref_cell_pheno_GSE144136[ref_cell_pheno_GSE144136$X %in% LSSMS_males_hr, ] # 5 cases, 5 controls


LSSMS_total_hr = c(LSSMS_females_hr, LSSMS_males_hr)
LSSMS_metadata_hr = GSE144136_GSE213982_LSSMS_hr[GSE144136_GSE213982_LSSMS_hr$person %in% LSSMS_total_hr,]


# Solution: Ast, End, Mic, Mix, Oli, OPC -> TAKE all cells
# Solution: ExN, InN, -> Sample 6000 cells
tmp_ExN_indeces = which(LSSMS_metadata_hr$cell_overall == "ExN")
set.seed(12345)
tmp_selected_ExN_indeces = sample(tmp_ExN_indeces, size = 6000, replace = FALSE)
set.seed(NULL)

tmp_InN_indeces = which(LSSMS_metadata_hr$cell_overall == "InN")
set.seed(12345)
tmp_selected_InN_indeces = sample(tmp_InN_indeces, size = 6000, replace = FALSE)
set.seed(NULL)

LSSMS_metadata_hr$is_sampled = NA
for (i in 1:nrow(LSSMS_metadata_hr)){
  
  if (LSSMS_metadata_hr$cell_overall[i] %!in% c("ExN", "InN")){
    LSSMS_metadata_hr$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_ExN_indeces){
    LSSMS_metadata_hr$is_sampled[i] = TRUE
  } else if (i %in% tmp_selected_InN_indeces){
    LSSMS_metadata_hr$is_sampled[i] = TRUE
  } else {
    LSSMS_metadata_hr$is_sampled[i] = FALSE
  }
}

table(LSSMS_metadata_hr$is_sampled, LSSMS_metadata_hr$cell_overall)

"
         Ast   End   ExN   InN   Mic   Mix   Oli   OPC
  FALSE     0     0 16491  1943     0     0     0     0
  TRUE   4129   805  6000  6000   987   723  6548  1937
"
LSSMS_picked_cells_hr = LSSMS_metadata_hr[LSSMS_metadata_hr$is_sampled, "sample_name"]

LSSMS_processed_hr = GSE144136_GSE213982_count_matrix_harmonized
LSSMS_processed_hr = LSSMS_processed_hr[,LSSMS_picked_cells_hr]
dim(LSSMS_processed_hr) # 22802 27129

LSSMS_processed_hr = as.matrix(LSSMS_processed_hr)
colnames(LSSMS_processed_hr)[1:10]
rownames(LSSMS_processed_hr)[1:10]

replacement_colnames = sapply(colnames(LSSMS_processed_hr), function(x){
  
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
   13601    12000     1528 
"



colnames(LSSMS_processed_hr) = replacement_colnames
LSSMS_processed_hr = as.data.frame(LSSMS_processed_hr)
genes_column = data.frame("Gene" = rownames(LSSMS_processed_hr))
LSSMS_processed_hr = cbind(genes_column, LSSMS_processed_hr)
colnames(LSSMS_processed_hr)[1:10]
LSSMS_processed_hr$Gene[1:10]
rownames(LSSMS_processed_hr) = NULL

all(LSSMS_processed_hr$Gene == Harmonized_matrix_gene_check$name_for_cibersort) # TRUE

# Saving as TSV in TXT file
fwrite(LSSMS_processed_hr,
       file = "ref_count_matrix_LSSMS_harmon.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

rm(LSSMS_processed_hr)
gc()


################### Build LSSMS harmonized ###################

# Making directories
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/input", recursive = TRUE)
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/output", recursive = TRUE)

# Moving reference expression files
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/ref_count_matrix_LSSMS_harmon.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/input/ref_count_matrix_dense.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


"
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/input:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/output:/src/outdir \
cibersortx/fractions \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--single_cell TRUE \
--refsample ref_count_matrix_dense.txt \
--fraction 0.5
"
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


################### Run harmonized LSSMS  ###################

# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/output_files", recursive = TRUE)

# Place files to the directory
source_samples = LSSMS_total_hr
full_valid_bulk_counts_LSSMS_hr = full_valid_bulk_counts_harmon[,colnames(full_valid_bulk_counts_harmon) %!in% source_samples]

fwrite(full_valid_bulk_counts_LSSMS_hr,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  
'


################### Evaluate fractions LSSMS hr  ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)


# Exploring reference fractions
source_samples = LSSMS_total_hr
tmp_cell_full = GSE144136_GSE213982_cell_proprotions_simplif[GSE144136_GSE213982_cell_proprotions_simplif$PARTICIPANT %!in% source_samples, ]

validation_cell_full_ordered = tmp_cell_full
imputed_fractions_full = imputed_fractions_full[imputed_fractions_full$Mixture %in% validation_cell_full_ordered$PARTICIPANT,]
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]

all(imputed_fractions_full$Mixture == validation_cell_full_ordered$PARTICIPANT) # TRUE

dir.create("cell_type_imputation_performance/imputation_correlations_docker_LSSMS_harmon")
cell_types = c("Glia", "Neuronal", "Other")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_LSSMS_harmon/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~27K LSSMS harmon. (20 donors), 0 frac: ",
                             color_dots = "#00008B")
  
})


combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_LSSMS_harmon",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_LSSMS_harmon/combined_img.png",
  ncol=3,
  border_px = 10
)


################### Evaluate high res. LSSMS hr  ###################

output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/output_files"
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

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window12.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

# Select appropriate pseudobulked counts 
source_samples = LSSMS_total_hr
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


dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_harmon")
cell_types = c("Glia", "Neuronal", "Other")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_harmon/", x, ".png")
  plot_gene_expression_imput_corr_harmon(cell_type = x,
                                  loaded_imputed_list = loaded_imputed_val_counts,
                                  real_counts_df = tmp_counts,
                                  selected_validation_samples = required_tmp_samples,
                                  path_to_save = path,
                                  title_prefix = "~27K LSSMS harmon. (20 donors), 0 frac: ",
                                  color_dots = "#00008B")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_harmon",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_harmon/combined_img.png",
  ncol=2,
  border_px = 10
)



# Move all files in one place for combining
merged_figures = list.files("cell_type_imputation_performance", full.names = TRUE, recursive = TRUE, pattern = "combined_img")
merged_figures = merged_figures[stri_detect_fixed(merged_figures, pattern = "highres")]
merged_figures = merged_figures[stri_detect_fixed(merged_figures, pattern = "harmon")]

for (i in 1:length(merged_figures)){
  
  names_split = unlist(stri_split_fixed(merged_figures[i], "/"))
  sel_name = names_split[2]
  sel_name = stri_replace_all_fixed(sel_name, pattern = "imputation_correlations_highres_docker", replacement = "")
  
  if (sel_name == "_LSSMS_harmon"){
    sel_name = "LSSMS_harmon"
  }
  
  if (sel_name == "_sSig_harmon"){
    sel_name = "sampled_signature_harmon"
  }
  
  # rewrite file
  file_path = paste0("cell_type_imputation_performance/Merged_Highres_performance_harmon/", sel_name, ".png")
  
  file.copy(merged_figures[i],
            file_path,
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
  
}


################### Run CPM LSSMS   ###################

# Directory creation
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/input_files", recursive = TRUE)
dir.create("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/output_files", recursive = TRUE)

# Place files to the directory
source_samples = LSSMS_total_hr
full_valid_bulk_counts_LSSMS_hr = full_valid_bulk_counts_harmon[,colnames(full_valid_bulk_counts_harmon) %!in% source_samples]
full_valid_bulk_counts_LSSMS_hr_dge = full_valid_bulk_counts_LSSMS_hr
full_valid_bulk_counts_LSSMS_hr_dge$Gene = NULL

full_valid_bulk_counts_LSSMS_hr_dge = DGEList(counts = full_valid_bulk_counts_LSSMS_hr_dge)
full_valid_bulk_counts_LSSMS_hr_dge = calcNormFactors(full_valid_bulk_counts_LSSMS_hr_dge)
full_valid_bulk_counts_LSSMS_hr_norm = cpm(full_valid_bulk_counts_LSSMS_hr_dge)
full_valid_bulk_counts_LSSMS_hr_norm = as.data.frame(full_valid_bulk_counts_LSSMS_hr_norm)
full_valid_bulk_counts_LSSMS_hr_norm = cbind(Gene = full_valid_bulk_counts_LSSMS_hr$Gene, full_valid_bulk_counts_LSSMS_hr_norm)

fwrite(full_valid_bulk_counts_LSSMS_hr_norm,
       file = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/input_files/validation_pseudobulked_counts_full.txt", 
       sep = "\t",
       row.names = FALSE, 
       col.names = TRUE, 
       quote = FALSE)

# Place basic files in the directory
# Signature matrix file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/output/CIBERSORTx_ref_count_matrix_dense_inferred_phenoclasses.CIBERSORTx_ref_count_matrix_dense_inferred_refsample.bm.K999.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/input_files/custom_signature_mtx.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# Source GEP file
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/output/CIBERSORTx_cell_type_sourceGEP.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/input_files/source_GEP.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

# Reference profile
file.copy("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/CIBERSORTx_main/singature_matrix_build_LSSMS_harmon/input/ref_count_matrix_dense.txt",
          "/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/input_files/reference_profile.txt",
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)


# Calculation of cell-specific expression profile for validation sample
# Run in the terminal
'
docker run \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/input_files:/src/data \
-v /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/output_files:/src/outdir \
cibersortx/hires \
--username aleksandr.sokolov@uu.se  \
--token <token_from_cibersortx_website>  \
--mixture validation_pseudobulked_counts_full.txt  \
--sigmatrix custom_signature_mtx.txt  \
--QN FALSE  \
--rmbatchSmode TRUE \
--refsample reference_profile.txt  
'

################### Evaluate fractions CPM LSSMS  ###################

# Loading imputed fractions
imputed_fractions_full = fread("Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/output_files/CIBERSORTxGEP_NA_Fractions-Adjusted.txt")
imputed_fractions_full = as.data.frame(imputed_fractions_full)


# Exploring reference fractions
source_samples = LSSMS_total_hr
tmp_cell_full = GSE144136_GSE213982_cell_proprotions_simplif[GSE144136_GSE213982_cell_proprotions_simplif$PARTICIPANT %!in% source_samples, ]

validation_cell_full_ordered = tmp_cell_full
imputed_fractions_full = imputed_fractions_full[imputed_fractions_full$Mixture %in% validation_cell_full_ordered$PARTICIPANT,]
rownames(validation_cell_full_ordered) = validation_cell_full_ordered$PARTICIPANT
validation_cell_full_ordered = validation_cell_full_ordered[imputed_fractions_full$Mixture,]

all(imputed_fractions_full$Mixture == validation_cell_full_ordered$PARTICIPANT) # TRUE

dir.create("cell_type_imputation_performance/imputation_correlations_docker_LSSMS_CPM")
cell_types = c("Glia", "Neuronal", "Other")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_docker_LSSMS_CPM/", x, ".png")
  plot_cell_prop_correlation(estimated_proportions = imputed_fractions_full[,x],
                             real_proportions = validation_cell_full_ordered[,x],
                             cell_type = x,
                             path_to_save = path,
                             title_prefix = "~27K LSSMS harmon. CPM (20 donors), 0 frac: ",
                             color_dots = "#00008B")
  
})


combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_docker_LSSMS_CPM",
  output_file = "cell_type_imputation_performance/imputation_correlations_docker_LSSMS_CPM/combined_img.png",
  ncol=3,
  border_px = 10
)

# Fractions appear to be completely the same!

################### Evaluate high res. CPM LSSMS  ###################

output_folder = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/output_files"
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

names_for_imputed = stri_replace_all_fixed(val_out_files, pattern = "Data_preprocessing_analysis/cell_expression_imputation/validation_sample_LSSMS_harmon_CPM/output_files/CIBERSORTxHiRes_NA_", replacement = "")
names_for_imputed = stri_replace_all_fixed(names_for_imputed, pattern = "_Window12.txt", replacement = "")
names(loaded_imputed_val_counts) = names_for_imputed

# Select appropriate pseudobulked counts 
source_samples = LSSMS_total_hr
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

tmp_counts_dge = DGEList(counts = tmp_counts)
tmp_counts_dge = calcNormFactors(tmp_counts_dge)
tmp_counts_norm = cpm(tmp_counts_dge)
tmp_counts_norm = as.data.frame(tmp_counts_norm)
tmp_counts = tmp_counts_norm

required_tmp_samples = colnames(tmp_counts)
required_tmp_samples = sapply(required_tmp_samples, function(x){
  x = stri_split_fixed(x, pattern = "_")
  x = unlist(x)
  x = x[1]
  return(x)
})
required_tmp_samples = unique(required_tmp_samples)
required_tmp_samples = required_tmp_samples[required_tmp_samples %in% colnames(loaded_imputed_val_counts[[1]])]


dir.create("cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_CPM")
cell_types = c("Glia", "Neuronal", "Other")
lapply(cell_types, function(x){
  path=paste0("cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_CPM/", x, ".png")
  plot_gene_expression_imput_corr_harmon(cell_type = x,
                                         loaded_imputed_list = loaded_imputed_val_counts,
                                         real_counts_df = tmp_counts,
                                         selected_validation_samples = required_tmp_samples,
                                         path_to_save = path,
                                         title_prefix = "~27K LSSMS harmon. CPM (20 donors), 0 frac: ",
                                         color_dots = "#00008B")
  
})

combine_pngs(
  input_dir = "cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_CPM",
  output_file = "cell_type_imputation_performance/imputation_correlations_highres_docker_LSSMS_CPM/combined_img.png",
  ncol=2,
  border_px = 10
)


################### Cohort-specific imputations ###################
dir.create("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts")

# Gene mapping file
RNAseq_gene_probes = as.data.frame(fread("Data_preprocessing_analysis/RNAseq_gene_probes.csv", header = TRUE))
RNAseq_gene_probes$V1 = NULL
rownames(RNAseq_gene_probes) = RNAseq_gene_probes$ID
RNAseq_gene_probes_genes_check = check_gene_symbol_NIH(PRF_gene_symbols=RNAseq_gene_probes$Gene_symbol, 
                                                       PRF_ref_NIH_expanded=Homo_Sapiens_Gene_info_NIH_expanded, 
                                                       PRF_replace_NA_with_old=TRUE)
RNAseq_gene_probes$Suggested.Symbol = RNAseq_gene_probes_genes_check$Suggested.Symbol

RNAseq_gene_probes$name_for_cibersort = stri_replace_all_fixed(RNAseq_gene_probes$Suggested.Symbol, 
                                                                pattern = " ",
                                                                replacement = "_")
RNAseq_gene_probes$name_for_cibersort = stri_replace_all_fixed(RNAseq_gene_probes$name_for_cibersort, 
                                                                pattern = "-",
                                                                replacement = "__")
RNAseq_gene_probes$name_for_cibersort = toupper(RNAseq_gene_probes$name_for_cibersort)

any(duplicated(RNAseq_gene_probes$name_for_cibersort)) # There are duplicated symbols
RNAseq_gene_probes$name_for_cibersort[duplicated(RNAseq_gene_probes$name_for_cibersort)] %>% unique(.) %>% length(.) # 444 symbols
# To avoid imputation problems genes with multiple transcripts were discarded from imputation


count_preparator = function(count_df){
  
  # Exclude first rows
  count_df = count_df[-(1:4), ]
  
  print(dim(count_df))
  
  # Obtain appropriate gene symbols
  tmp_RNAseq_gene_probes = RNAseq_gene_probes
  
  # Printout check
  print(paste0("Total IDs: ", length(unique(count_df$V1))))
  
  print(paste0("All count IDs are in the reference IDs: ", all(count_df$V1 %in% tmp_RNAseq_gene_probes$ID)))
  tmp_RNAseq_gene_probes = tmp_RNAseq_gene_probes[count_df$V1,]
  
  # Remove duplicated expression
  duplicated_symbols = tmp_RNAseq_gene_probes$name_for_cibersort[duplicated(tmp_RNAseq_gene_probes$name_for_cibersort)]
  tmp_RNAseq_gene_probes = tmp_RNAseq_gene_probes[tmp_RNAseq_gene_probes$name_for_cibersort %!in% duplicated_symbols, ]
  count_df = count_df[count_df$V1 %in% tmp_RNAseq_gene_probes$ID,]
  print(dim(count_df))
  
  required_gene_symbols = tmp_RNAseq_gene_probes$name_for_cibersort
  
  count_df$V1 = required_gene_symbols
  
  colnames(count_df)[1] = "Gene"
  
  # Drop genes mapped to multiple IDs
  count_df = count_df[!stri_detect_fixed(count_df$Gene, pattern = ";"),]
  
  return(count_df)
  
}

place_source_highres_files_to_dir = function(fixed_counts,
                                             drirectory_source, 
                                             directory_target){
  
  fwrite(fixed_counts,
         file = paste0(directory_target, "/input/bulk_counts.txt"), 
         sep = "\t",
         row.names = FALSE, 
         col.names = TRUE, 
         quote = FALSE)
  
  # Get signature matrix
  drirectory_source_out = paste0(drirectory_source, "/output")
  sig_matrix_file = list.files(drirectory_source_out, full.names = TRUE, pattern = "_inferred_refsample.bm.K999.txt")
  file.copy(sig_matrix_file,
            paste0(directory_target, "/input/custom_signature_mtx.txt"),
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
  
  # Get reference profile
  drirectory_source_in = paste0(drirectory_source, "/input")
  ref_profile_file = list.files(drirectory_source_in, full.names = TRUE, pattern = "ref_count")
  file.copy(ref_profile_file,
            paste0(directory_target, "/input/reference_profile.txt"),
            overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
}

################### Processing of all single cell cohort with CIBERSORTx ###################

# GSE102556 11+37 = 48 samples - several tissues
# GSE243356 32+29 = 61 samples - one tissue
# GSE248260 9+15 = 24 samples - one tissue
# GSE202537 60+7 = 67 samples - 3 tissues
# GSE101521 38+21 = 59 samples one tissue
# GSE144136 15+17 = 32 samples (already single cell) - done separately (one tissue)
# GSE213982 18+17= 36 samples (already single cell) - done separately (one tissue)

impute_expession_cibersort = function(cohort_id,
                                      source_matrix_folder,
                                      result_folder_suffix){
  
  
  writeLines(paste0("Working on dataset: ", cohort_id, " Source matrix: ", source_matrix_folder))
  
  # Performs analysis with LSSMS and Exn
  # Target folders are placed in "home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts"
  
  pheno_path = paste0("Data_preprocessing_analysis/", cohort_id, "_results/", cohort_id, "_pheno_curated.csv")
  cohort_pheno = read.csv(pheno_path)
  
  # Filter required expression
  selected_tissues = c(
    "Dorsolateral prefrontal cortex (dlPFC; BA8/9)",
    "Nac"
  )
  
  if (length(unique(cohort_pheno$TISSUE)) > 1) {
    writeLines("Detected more than one tissue")
    cohort_pheno = cohort_pheno[cohort_pheno$TISSUE %in% selected_tissues, ]
    writeLines(paste0("Selected tissue: ", unique(cohort_pheno$TISSUE)))
  }
  
  count_path = paste0("Data_preprocessing_analysis/", cohort_id, "_results/", cohort_id, "_expression.csv")
  cohort_counts = as.data.frame(fread(count_path))
  cohort_counts = cohort_counts[,colnames(cohort_counts) %in% c("V1", cohort_pheno$RUN_ID), ]
  cohort_counts_fixed = count_preparator(cohort_counts)
  
  # Making all directories
  result_path_base = paste0("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/",cohort_id,result_folder_suffix)
  result_path_input = paste0("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/",cohort_id,result_folder_suffix,"/input")
  result_path_output = paste0("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/",cohort_id,result_folder_suffix,"/output")
  dir.create(result_path_input, recursive = TRUE)
  dir.create(result_path_output, recursive = TRUE)
  
  
  # Moving files to folder
  place_source_highres_files_to_dir(
    fixed_counts = cohort_counts_fixed,
    drirectory_source = source_matrix_folder,
    directory_target = result_path_base
  )
  
  # docker command
  system2(
    "docker", 
    args = c(
      "run",
      "-v", paste0(result_path_input, ":/src/data"),
      "-v", paste0(result_path_output, ":/src/outdir"),
      "cibersortx/hires",
      "--username", "aleksandr.sokolov@uu.se",
      "--token", "<token_from_cibersortx_website>",
      "--mixture", "bulk_counts.txt",
      "--sigmatrix", "custom_signature_mtx.txt",
      "--QN", "FALSE",
      "--rmbatchSmode", "TRUE",
      "--refsample", "reference_profile.txt"
    )
  )
}


cohorts_to_deconvolute = c("GSE248260", "GSE102556","GSE202537", "GSE243356", "GSE101521")

for (i in 1:length(cohorts_to_deconvolute)){
  
  impute_expession_cibersort(
    cohort_id = cohorts_to_deconvolute[i],
    source_matrix_folder = "CIBERSORTx_main/singature_matrix_build_LSSMS_harmon",
    result_folder_suffix = "_LSSMS"
  )
  
  impute_expession_cibersort(
    cohort_id = cohorts_to_deconvolute[i],
    source_matrix_folder = "CIBERSORTx_main/singature_matrix_build_sampled_harmon",
    result_folder_suffix = "_ExN"
  )
  
}


################### Test voom of output (raw initial counts) ###################
sample_out_matrix = smart_fread("Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/GSE202537_LSSMS/output/CIBERSORTxHiRes_NA_Neuronal_Window12.txt")
sample_pheno = smart_fread("Data_preprocessing_analysis/GSE202537_results/GSE202537_pheno_curated.csv")
sample_pheno = sample_pheno[sample_pheno$TISSUE == "Nac", ]
rownames(sample_out_matrix) = sample_out_matrix$GeneSymbol
sample_out_matrix$GeneSymbol = NULL
# = sample_out_matrix[!apply(sample_out_matrix, 1, function(x) length(unique(x)) == 1), ]
sample_out_matrix = sample_out_matrix[!apply(sample_out_matrix,1, function(x) any(is.na(x))), ]
sample_out_matrix = sample_out_matrix[apply(sample_out_matrix, 1, function(x){
  freq_table = table(x)
  max_prop = max(freq_table) / length(x)   # proportion of the most frequent value
  (1 - max_prop) > 0.8
}), ]

all(sample_pheno$RUN_ID == colnames(sample_out_matrix)) # TRUE

# Design matrix
sample_pheno$SUICIDE = factor(sample_pheno$SUICIDE, levels = c("Control", "Suicide"))

tmp_design = model.matrix(~ SUICIDE, data = sample_pheno)

# egeR filter
test_tmp_dge = DGEList(counts=sample_out_matrix)
keep = filterByExpr(test_tmp_dge, tmp_design)
test_tmp_dge = test_tmp_dge[keep,,keep.lib.sizes=FALSE]
test_tmp_dge = calcNormFactors(test_tmp_dge)
test_tmp_voom = voom(test_tmp_dge, tmp_design, plot=TRUE)

# DE analysis
fit = lmFit(test_tmp_voom, tmp_design)
fitE = eBayes(fit)
Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
Top_table_no_covar_TMP$ID = rownames(Top_table_no_covar_TMP)
# *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975

SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
SE = SE[,2]
SE = SE[Top_table_no_covar_TMP$ID]
Top_table_no_covar_TMP$SE = SE
all(names(SE) == Top_table_no_covar_TMP$ID) # TRUE


################### Run DE across cohorts in a loop (raw initial counts) ###################

cohorts_to_deconvolute = c("GSE248260", "GSE102556", "GSE202537", "GSE243356", "GSE101521")
dir.create("DEs_cell_specif")

run_DE_cell_spec_cohort = function(cohort_id, cell_suffix, type = "LSSMS"){
  
  m = paste0("Working on: ", cohort_id)
  print(m)
  
  if (type != "LSSMS"){
    imputation_folder = stri_replace_all_fixed("Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/SELECTED_ID_CELL_TYPE/output/", 
                                               pattern = "SELECTED_ID",
                                               replacement = cohort_id)
    imputation_folder = stri_replace_all_fixed(imputation_folder, 
                                               pattern = "CELL_TYPE",
                                               replacement = cell_suffix)
  } else {
    
    imputation_folder = stri_replace_all_fixed("Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/SELECTED_ID_LSSMS/output/", 
                                               pattern = "SELECTED_ID",
                                               replacement = cohort_id)
  }
  
  print(imputation_folder)
  
  out_files = list.files(imputation_folder, recursive = TRUE, full.names = TRUE, pattern = ".txt")
  out_files = out_files[stri_detect_fixed(out_files, pattern = "CIBERSORTxHiRes")]
  out_files = out_files[stri_detect_fixed(out_files, pattern = paste0(cell_suffix, "_Wind"))]
  
  print(out_files)
  
  phenotype_path = stri_replace_all_fixed("Data_preprocessing_analysis/SELECTED_ID_results/SELECTED_ID_pheno_curated.csv", 
                                          pattern = "SELECTED_ID",
                                          replacement = cohort_id)
  
  imported_pheno = smart_fread(phenotype_path)
  
  
  # Different blocks depending on dataset
  if (cohort_id == "GSE248260"){
    imported_pheno$SUICIDE = factor(imported_pheno$SUICIDE, levels = c("Control", "Suicide"))
    
  } else if (cohort_id == "GSE102556"){
    imported_pheno = imported_pheno[imported_pheno$TISSUE == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)", ]
    imported_pheno$SUICIDE == factor(imported_pheno$SUICIDE, levels = c("CONTROL", "SUICIDE"))
    
  } else if (cohort_id == "GSE202537"){
    imported_pheno = imported_pheno[imported_pheno$TISSUE == "Nac", ]
    imported_pheno$SUICIDE = factor(imported_pheno$SUICIDE, levels = c("Control", "Suicide"))
    
  } else if (cohort_id == "GSE243356"){
    imported_pheno$SUICIDE = factor(imported_pheno$GROUP, levels = c("healthy", "suicide"))
    
  } else if (cohort_id == "GSE101521"){
    imported_pheno$SUICIDE = factor(imported_pheno$SUICIDE, levels = c("Non-suicide", "Suicide"))
  }
  
  # Load counts
  imputed_counts = smart_fread(out_files)
  rownames(imputed_counts) = imputed_counts$GeneSymbol
  imputed_counts$GeneSymbol = NULL
  
  # Order check
  print(paste0("All participants match with RUN_ID: ", all(imported_pheno$RUN_ID == colnames(imputed_counts))))
  
  
  # DE part
  imputed_counts = imputed_counts[!apply(imputed_counts,1, function(x) any(is.na(x))), ]
  imputed_counts = imputed_counts[apply(imputed_counts, 1, function(x){
    freq_table = table(x)
    max_prop = max(freq_table) / length(x)   # proportion of the most frequent value
    (1 - max_prop) > 0.8
  }), ]
  
  if (nrow(imputed_counts) < 1){
    print("Bad imputation result -> Exit")
    return(FALSE)
  }
  
  all(imported_pheno$RUN_ID == colnames(imputed_counts)) # TRUE
  
  # Design matrix
  tmp_design = model.matrix(~ SUICIDE, data = imported_pheno)
  
  # egeR filter
  test_tmp_dge = DGEList(counts=imputed_counts)
  keep = filterByExpr(test_tmp_dge, tmp_design)
  test_tmp_dge = test_tmp_dge[keep,,keep.lib.sizes=FALSE]
  test_tmp_dge = calcNormFactors(test_tmp_dge)
  test_tmp_voom = voom(test_tmp_dge, tmp_design, plot=TRUE)
  
  # DE analysis
  fit = lmFit(test_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  Top_table_no_covar_TMP$ID = rownames(Top_table_no_covar_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[Top_table_no_covar_TMP$ID]
  Top_table_no_covar_TMP$SE = SE
  all(names(SE) == Top_table_no_covar_TMP$ID) # TRUE
  
  Top_table_no_covar_TMP$Gene_symbol = sapply(Top_table_no_covar_TMP$ID, function(x){
    out = RNAseq_gene_probes[RNAseq_gene_probes$name_for_cibersort == x, "Suggested.Symbol"]
    out = unlist(out)
    out = unique(out)
    out = out[1]
    return(out)
  })
  Top_table_no_covar_TMP$Tissue = imported_pheno$TISSUE[1]
  Top_table_no_covar_TMP$Tissue_type = "Brain"
  Top_table_no_covar_TMP$Cell_group = cell_suffix
  Top_table_no_covar_TMP$Technology = "RNA-seq"
  
  saving_path = paste0("DEs_cell_specif/", cohort_id, "_DE_", cell_suffix,".csv")
  write.csv(Top_table_no_covar_TMP, saving_path)
  
}

for (i in 1:length(cohorts_to_deconvolute)){
  
  run_DE_cell_spec_cohort(
    cohort_id = cohorts_to_deconvolute[i],
    cell_suffix = "Neuronal"
  )
  
  run_DE_cell_spec_cohort(
    cohort_id = cohorts_to_deconvolute[i],
    cell_suffix = "Glia"
  )
  
  run_DE_cell_spec_cohort(
    cohort_id = cohorts_to_deconvolute[i],
    cell_suffix = "ExN",
    type = "ExN"
  )
}

################### Count fractions in a loop across cohorts and plot them ###################

count_frac_cell_spec_cohort = function(cohort_id, cell_suffix, type = "LSSMS"){
  
  m = paste0("Working on: ", cohort_id)
  print(m)
  
  if (cohort_id == "GSE144136" & type == "LSSMS"){
    pheno = ref_cell_pheno_GSE144136
    cell_fractions = GSE144136_GSE213982_cell_proprotions_simplif[GSE144136_GSE213982_cell_proprotions_simplif$PARTICIPANT %in% pheno$PARTICIPANT,]
    cell_fractions$cohort_id = cohort_id
    return(cell_fractions)
    
  } else if (cohort_id == "GSE213982" & type == "LSSMS"){
    pheno = ref_cell_pheno_GSE213982
    cell_fractions = GSE144136_GSE213982_cell_proprotions_simplif[GSE144136_GSE213982_cell_proprotions_simplif$PARTICIPANT %in% pheno$PARTICIPANT,]
    cell_fractions$cohort_id = cohort_id
    return(cell_fractions)
    
  }  else if (cohort_id %in% c("GSE144136", "GSE213982") & type != "LSSMS"){
    
    if (cohort_id == "GSE144136"){
      pheno = ref_cell_pheno_GSE144136
    } else {
      pheno = ref_cell_pheno_GSE213982
    }
    
    cell_fractions = pheno[, c("PARTICIPANT","Ast","End" , "ExN","InN","Mic","Mix","Oli","OPC")]
    cell_fractions$cohort_id = cohort_id
    return(cell_fractions)
  }
  
  if (type != "LSSMS"){
    imputation_folder = stri_replace_all_fixed("Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/SELECTED_ID_CELL_TYPE/output/", 
                                               pattern = "SELECTED_ID",
                                               replacement = cohort_id)
    imputation_folder = stri_replace_all_fixed(imputation_folder, 
                                               pattern = "CELL_TYPE",
                                               replacement = cell_suffix)
  } else {
    
    imputation_folder = stri_replace_all_fixed("Data_preprocessing_analysis/cell_expression_imputation/actual_cohorts/SELECTED_ID_LSSMS/output/", 
                                               pattern = "SELECTED_ID",
                                               replacement = cohort_id)
  }
  
  print(imputation_folder)
  
  out_files = list.files(imputation_folder, recursive = TRUE, full.names = TRUE, pattern = "Fractions-Adjusted.txt")
  
  loaded_fractions = smart_fread(out_files)
  
  loaded_fractions$cohort_id = cohort_id
  
  colnames(loaded_fractions)[1] = "PARTICIPANT"
  
  return(loaded_fractions)
  
  
}
cohorts_to_deconvolute = c("GSE248260", "GSE102556", "GSE202537", "GSE243356", "GSE101521", "GSE144136", "GSE213982")


dir.create("cell_fraction_plots")
cell_subsets = c("ExN", "Neuronal", "Glia")
colors = c("skyblue", "darkblue", "orange")

for (i in 1:length(cell_subsets)){
  
  if (cell_subsets[i] %in% c("Neuronal", "Glia")){
    curr_type = "LSSMS"
  } else {
    curr_type = cell_subsets[i]
  }
  
  cell_fractions_est = lapply(cohorts_to_deconvolute, function(x){
    df = count_frac_cell_spec_cohort(cohort_id = x, cell_suffix=cell_subsets[i], type=curr_type)
    return(df)
  })
  overlapping_names = sapply(cell_fractions_est, colnames)
  overlapping_names = Reduce(intersect, overlapping_names)
  
  cell_fractions_est = lapply(cell_fractions_est, function(x){
    df = x[, overlapping_names]
    return(df)
  })
  cell_fractions_est = do.call(rbind, cell_fractions_est)
  
  colnames(cell_fractions_est) = stri_replace_all_fixed(colnames(cell_fractions_est), pattern = cell_subsets[i], replacement = "CURRENT_CELL")
  
  plot = ggplot(cell_fractions_est, aes(x = cohort_id, y = CURRENT_CELL)) +
    geom_jitter(col = colors[i]) +
    labs(x = "Cohort",
         y = cell_subsets[i]) +
    theme_bw(base_size = 15)
  
  png(file = paste0("cell_fraction_plots/", cell_subsets[i],".png"), width = 8, height = 4, units = "in", res = 300)
  print(plot)
  dev.off()
}

combine_pngs(
  input_dir = "cell_fraction_plots",
  output_file = "cell_fraction_plots/combined_img.png",
  ncol=1,
  border_px = 10
)

################### Run DE across PSEUDOBULK cohorts in a loop ###################

cohorts_to_deconvolute = c("GSE144136", "GSE213982")

run_DE_sc_cohorts = function(cohort_id, cell_suffix){
  
  
  m = paste0("Working on: ", cohort_id)
  print(m)
  
  base_df_path = paste0("Data_preprocessing_analysis/", cohort_id, "_results")
  
  if (cell_suffix %in% c("Ast","End", "ExN","InN","Mic","Mix","Oli","OPC")){
    
    m = paste0("Working on cell-specific subtupe: ", cell_suffix)
    print(m)
    count_file = list.files(base_df_path, full.names = TRUE, pattern = "_counts_cell_type.csv")
    tmp_counts = smart_fread(count_file)
    
  } else {
    m = paste0("Working on aggregated subtupe: ", cell_suffix)
    print(m)
    # We are reusing simplified pseudobulk from LSSMS
    tmp_counts = GSE144136_GSE213982_counts_simplif_cell
  }
  
  # Importing pheno
  pheno_file = list.files(base_df_path, full.names = TRUE, pattern = "_pheno_curated")
  imported_pheno = smart_fread(pheno_file)
  
  
  # Different blocks depending on dataset
  if (cohort_id == "GSE144136"){
    imported_pheno$SUICIDE = factor(imported_pheno$GROUP, levels = c("Control", "Major Depressive Disorder (MDD)"))
  } else if (cohort_id == "GSE213982"){
    imported_pheno$SUICIDE = factor(imported_pheno$GROUP, levels = c("Control", "Case"))
  } else {
    print("Wrog cohort!")
    return(FALSE)
  }
  
  # Order check
  tmp_counts_selection_index = sapply(colnames(tmp_counts), function(x){
    participant = unlist(stri_split_fixed(x, pattern = "_"))[1]
    cell_type =  unlist(stri_split_fixed(x, pattern = "_"))[2]
    
    if (participant %in% imported_pheno$PARTICIPANT & cell_type == cell_suffix){
      return(TRUE)
    } else {
      return(FALSE)
    }
      
  })
  tmp_counts = tmp_counts[,tmp_counts_selection_index]
  colnames(tmp_counts) = stri_replace_all_fixed(colnames(tmp_counts),
                                                pattern = paste0("_", cell_suffix),
                                                replacement = "")
  
  
  m = paste0("All count participants are in pheno: ", all(colnames(tmp_counts) %in% imported_pheno$PARTICIPANT))
  print(m)
  tmp_counts = tmp_counts[,imported_pheno$PARTICIPANT]
  
  # Design matrix
  tmp_design = model.matrix(~ SUICIDE, data = imported_pheno)
  
  # egeR filter
  test_tmp_dge = DGEList(counts=tmp_counts)
  keep = filterByExpr(test_tmp_dge, tmp_design)
  test_tmp_dge = test_tmp_dge[keep,,keep.lib.sizes=FALSE]
  test_tmp_dge = calcNormFactors(test_tmp_dge)
  test_tmp_voom = voom(test_tmp_dge, tmp_design, plot=TRUE)
  
  # DE analysis
  fit = lmFit(test_tmp_voom, tmp_design)
  fitE = eBayes(fit)
  Top_table_no_covar_TMP = limma::topTable(fit = fitE, coef = 2, adjust.method = "fdr", number = Inf, confint = TRUE)
  Top_table_no_covar_TMP$ID = rownames(Top_table_no_covar_TMP)
  # *qt(alpha, df=fitE$df.total) calculate T critical value for confint alpha = 0.975
  
  SE = sqrt(fitE$s2.post) * fit$stdev.unscaled
  SE = SE[,2]
  SE = SE[Top_table_no_covar_TMP$ID]
  Top_table_no_covar_TMP$SE = SE
  all(names(SE) == Top_table_no_covar_TMP$ID) # TRUE
  
  Top_table_no_covar_TMP$Gene_symbol = sapply(Top_table_no_covar_TMP$ID, function(x){
    out = TMP_sig_matrix_gene_rows_gene_check[TMP_sig_matrix_gene_rows_gene_check$x == x, "Suggested.Symbol"]
    out = unlist(out)
    out = unique(out)
    out = out[1]
    return(out)
  })
  Top_table_no_covar_TMP$Tissue = imported_pheno$TISSUE[1]
  Top_table_no_covar_TMP$Tissue_type = "Brain"
  Top_table_no_covar_TMP$Cell_group = cell_suffix
  Top_table_no_covar_TMP$Technology = "scRNA-seq"
  
  # Remove duplicated symbols
  duplicated_genes = Top_table_no_covar_TMP$Gene_symbol[duplicated(Top_table_no_covar_TMP$Gene_symbol)]
  Top_table_no_covar_TMP = Top_table_no_covar_TMP[Top_table_no_covar_TMP$Gene_symbol %!in% duplicated_genes,]
  
  saving_path = paste0("DEs_cell_specif/", cohort_id, "_DE_", cell_suffix,".csv")
  write.csv(Top_table_no_covar_TMP, saving_path)
  
}

for (i in 1:length(cohorts_to_deconvolute)){
  
  run_DE_sc_cohorts(
    cohort_id = cohorts_to_deconvolute[i],
    cell_suffix = "Neuronal"
  )
  
  run_DE_sc_cohorts(
    cohort_id = cohorts_to_deconvolute[i],
    cell_suffix = "Glia"
  )
  
  run_DE_sc_cohorts(
    cohort_id = cohorts_to_deconvolute[i],
    cell_suffix = "ExN"
  )
}


################### Run RE meta across selected datasets ###################

plot_forest_meta_gene = function(gene_name, meta_df){
  
  tmp_df_genes = meta_df[meta_df$Gene_symbol == gene_name,]
  
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

run_SJ_meta_datasets = function(
    data_folder,
    cell_suffix,
    cohort_subset,
    enforce_reference_signif = FALSE
  ){
  
  # select only DEs
  detected_files = list.files(data_folder, full.names = TRUE, pattern = "_DE_")
  
  # filter based on cell suffix
  detected_files = detected_files[stri_detect_fixed(detected_files, pattern = cell_suffix)]
  
  # filter based on cohorts
  detected_files = detected_files[sapply(detected_files, function(x){
  
    x = unlist(stri_split_fixed(x, pattern = .Platform$file.sep))
    x = x[length(x)]
    cohort_id = unlist(stri_split_fixed(x, pattern = "_"))[1]
    
    if (cohort_id %in% cohort_subset){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })]
  
  # Reading 
  input_list = lapply(detected_files, function(x){
    
    path = unlist(stri_split_fixed(x, pattern = .Platform$file.sep))
    path = path[length(path)]
    cohort_id = unlist(stri_split_fixed(path, pattern = "_"))[1]
    
    x = smart_fread(x)
    x$Study = cohort_id
    return(x)
  })
  
  # Insuring removal of duplicates
  input_list = lapply(input_list, function(x){
    duplicated_names = x$Gene_symbol[duplicated(x$Gene_symbol)]
    x = x[x$Gene_symbol %!in% duplicated_names, ]
    return(x)
  })
  
  # Summarising DF
  print("Aggregating Log2FC and SE per gene per study")
  
  initial_study_list_reduced = lapply(input_list, function(x){
    
    orig_df = x
    
    x = x[!stri_detect_fixed(x$Gene_symbol, pattern = ";"),]
    x = x[!is.na(x$Gene_symbol),]
    x = x[x$Gene_symbol != "",]
    x = x %>% 
      dplyr::group_by(., Gene_symbol)%>% 
      dplyr::summarise(., 
                       avgLog2FC = mean(logFC), 
                       maxSE = max(SE),
                       p_val_init = max(P.Value))
    
    x = dplyr::ungroup(x)
    x$Tissue = unique(orig_df$Tissue)
    x$Tissue_type = unique(orig_df$Tissue_type)
    x$Technology = unique(orig_df$Technology)
    x$Cell_group = unique(orig_df$Cell_group)
    x$Study = unique(orig_df$Study)
    x = as.data.frame(x)
    return(x)
  })
  
  
  
  initial_study_list_reduced = do.call(rbind, initial_study_list_reduced)
  
  if (enforce_reference_signif){
    
    signif_reference_subset = initial_study_list_reduced[initial_study_list_reduced$p_val_init < 0.05,]
    signif_reference_subset = signif_reference_subset[signif_reference_subset$Technology == "scRNA-seq",]
    signif_genes = signif_reference_subset$Gene_symbol
    signif_genes = unique(signif_genes)
    initial_study_list_reduced = initial_study_list_reduced[initial_study_list_reduced$Gene_symbol %in% signif_genes,]
  }
  
  unique_genes_meta = unique(initial_study_list_reduced$Gene_symbol)
  
  
  
  print("Running meta per gene with multiple cores")
  
  meta_analysis_list = mclapply(unique_genes_meta, function(x){
    
    
    # meta df
    tmp_df_genes = initial_study_list_reduced[initial_study_list_reduced$Gene_symbol == x,]
    tmp_cohorts_up = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC > 0])
    tmp_cohorts_down = length(tmp_df_genes$avgLog2FC[tmp_df_genes$avgLog2FC < 0])
    cohorts_total = nrow(tmp_df_genes)
    
    if (cohorts_total < 5){
      tmp_output_df = data.frame(
        gene = x,
        cell_group = cell_suffix,
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
        Q.p = NA)
      return(tmp_output_df)
    }
    
    tmp_meta_model = rma.uni(yi = avgLog2FC, 
                             vi = maxSE^2, 
                             data = tmp_df_genes, 
                             method = "SJ", 
                             weighted = TRUE)
    
    tmp_output_df = data.frame(
      gene = x,
      cell_group = cell_suffix,
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
      Q.p = tmp_meta_model$QEp)
    
    return(tmp_output_df)
    
  }, mc.cores = 9)
  
  meta_analysis_df = do.call(rbind, meta_analysis_list)
  rownames(meta_analysis_df) = NULL
  rm(list = ls(pattern = "tmp_"))
  meta_analysis_df = dplyr::arrange(meta_analysis_df, meta_pval)
  meta_analysis_df$meta_FDR = p.adjust(meta_analysis_df$meta_pval, method = "fdr")
  meta_analysis_df = meta_analysis_df[,c("gene",
                                         "cell_group",
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
                                         "Q.p")]
  
  out_list = list()
  out_list[[1]] = meta_analysis_df
  out_list[[2]] = initial_study_list_reduced
  names(out_list) = c("meta_analysis_df", "initial_study_list_reduced")
  
  return(out_list)
  
}


cell_subsets = c("ExN", "Neuronal", "Glia")
stacked_reduced_lists = list()
stacked_META_lists = list()

for (i in 1:length(cell_subsets)){
  
  ME_output = run_SJ_meta_datasets(data_folder = "DEs_cell_specif",
                                   cell_suffix = cell_subsets[i],
                                   cohort_subset = c("GSE248260", "GSE102556", "GSE202537", "GSE243356", "GSE101521","GSE144136", "GSE213982"),
                                   enforce_reference_signif = TRUE)
  ME_output_meta = ME_output[[1]]
  ME_output_meta_signif = ME_output_meta[!is.na(ME_output_meta$meta_pval),]
  ME_output_meta_signif = ME_output_meta_signif[ME_output_meta_signif$meta_pval < 0.05,]
  ME_output_meta_signif = dplyr::arrange(ME_output_meta_signif, -meta_LFc)
  
  # Saving things
  stacked_reduced_lists[[cell_subsets[i]]] = ME_output[[2]]
  stacked_META_lists[[cell_subsets[i]]] = ME_output_meta_signif
}


View(stacked_META_lists[["ExN"]])
View(stacked_META_lists[["Neuronal"]])
View(stacked_META_lists[["Glia"]])

plot_gene_helper = function(gene_name, reduced_lists, subset){
  meta_df = reduced_lists[[subset]]
  print(plot_forest_meta_gene(gene_name, meta_df))
}
plot_gene_helper("CHPF2", stacked_reduced_lists, "Neuronal")

### Add significance of normal meta-analyses
file_path = "Meta_suicide_significant_genes.xlsx"
sheet_names = getSheetNames(file_path)
last_meta_results = lapply(sheet_names, function(x){
  file = openxlsx::read.xlsx(file_path, sheet = x)
  file$Analysis_name = x
  return(file)
})

stacked_META_lists_upd = stacked_META_lists
stacked_META_lists_upd = lapply(stacked_META_lists_upd, function(x){
  
  genes_to_select = x$gene
  genes_to_select = unique(genes_to_select)
  
  supported_analyses = sapply(genes_to_select, function(gene){
    
    message_vector = vector()
    
    for (i in 1:length(last_meta_results)){
      
      curr_meta = last_meta_results[[i]]
      curr_meta_analysis = unique(curr_meta$Analysis_name)
      
      if (gene %in% curr_meta$gene){
        message_vector = c(message_vector, curr_meta_analysis)
      }
    }
    
    if (length(message_vector) >= 1){
      message_vector = paste0(message_vector, collapse = "; ")
      message_vector = paste0("Supporting analyses: ", message_vector)
    } else {
      message_vector = "None"
    }
    return(message_vector)
  })
  x$supported_analyses = supported_analyses
  return(x)
})

View(stacked_META_lists_upd[["ExN"]])
View(stacked_META_lists_upd[["Neuronal"]])

################### Enrichment analysis for datasets ###################

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


all_gene_symbols = c(TMP_sig_matrix_gene_rows_gene_check$Suggested.Symbol, STAR_gene_names_check$Suggested.Symbol)
all_gene_symbols = unique(all_gene_symbols)

ENTREZ_Meta_universe = ENTREZ_genes_Homo_Sapiens[sapply(ENTREZ_genes_Homo_Sapiens$Symbol, function(x){
  if (x %in% all_gene_symbols){
    return(TRUE)
  }
  return(FALSE)
}),]

joined_list_enrichment = stacked_META_lists_upd
joined_list_enrichment = lapply(joined_list_enrichment, function(x){
  data = x$gene
  return(data)
})

# running enrichment
for (x in 1:length(joined_list_enrichment)){
  print(paste0("Working on: ", names(joined_list_enrichment)[x]))
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
                                                                 folder = "Enrichment_folder_deconvoluted",
                                                                 plot_name_pref = names(joined_list_enrichment)[x])
  
}


################### Saving processed files ###################

selected_enrichments_files = list.files("Enrichment_folder_deconvoluted", full.names = TRUE, pattern = "_KEGG.xlsx")
selected_enrichments_files_names = sapply(selected_enrichments_files, function(x){
  x = stri_split_fixed(x, pattern = "/") %>% unlist(.)
  x = x[length(x)]
  x = stri_split_fixed(x, pattern = "_") %>% unlist(.)
  x = x[1]
  x = paste0("Enrich. KEGG ", x)
  return(x)
})

loaded_enrichments = lapply(selected_enrichments_files, read.xlsx)
names(loaded_enrichments) = selected_enrichments_files_names

combined_results = c(stacked_META_lists_upd, loaded_enrichments)

wb = createWorkbook()

for (i in 1:6) {
  # Add a new worksheet with the sheet name
  addWorksheet(wb, names(combined_results)[i])
  
  # Write the data frame to the worksheet
  writeData(wb, sheet = names(combined_results)[i], combined_results[[i]])
  
  # Adjust column widths to fit the text
  setColWidths(wb, sheet = names(combined_results)[i], cols = 1:ncol(combined_results[[i]]), widths = "auto")
}

saveWorkbook(wb, "cell_deconvoluted_analysis.xlsx", overwrite = TRUE)
