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

################### Moderator thoughts ###################
# https://pubmed.ncbi.nlm.nih.gov/11836738/
# How many we can include? The paper https://onlinelibrary.wiley.com/doi/epdf/10.1002/sim.1040  suggests around 5-10 "cases" per variable!


################### PMI calculations ###################
# GSE208338
GSE208338_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE208338_results/GSE208338_pheno_curated.csv")
mean(GSE208338_pheno_tmp$PMI) # 41.95562 -> 42

# GSE5388
GSE5388_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE5388_results/GSE5388_pheno_curated.csv")
mean(GSE5388_pheno_tmp$POST_MORTEM_INTERVAL__HOURS_) # 33.2 -> 33.2

# GSE5389
GSE5389_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE5389_results/GSE5389_pheno_curated.csv")
mean(GSE5389_pheno_tmp$POST_MORTEM_INTERVAL__HOURS_) # 28.19048 -> 28.2

# GSE66937
GSE66937_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE66937_results/GSE66937_pheno_curated.csv")
length(unique(GSE66937_pheno_tmp$INDIVIDUAL))
# PMI is not reported in the dataset!
# We calculate PMI based on Table 1 in https://pmc.ncbi.nlm.nih.gov/articles/PMC8458545/#Sec2
# It is approximate due to 20 sample vs available 15
GSE66937_pheno_tmp_scraped = read_html("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/HTMLs/Identification of transcriptome alterations in the prefrontal cortex, hippocampus, amygdala and hippocampus of suicide victims - PMC.html")
GSE66937_pheno_tmp_scraped = GSE66937_pheno_tmp_scraped %>% html_elements("table") %>% html_table(fill = TRUE)
GSE66937_pheno_tmp_scraped = GSE66937_pheno_tmp_scraped[[2]]
mean(GSE66937_pheno_tmp_scraped$`PMI (h)`, na.rm=TRUE) # 31.55556 -> 31.6

# GSE199536
GSE199536_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE199536_results/GSE199536_pheno_curated.csv")
mean(GSE199536_pheno_tmp$PMI__H_A) # 31.275 -> 31.3

# GSE92538_U133A
GSE92538_U133A_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE92538_U133A_results/GSE92538_pheno_curated.csv")
mean(GSE92538_U133A_pheno_tmp$POST_MORTEM_INTERVAL)  # 24.89433 -> 24.9

# GSE92538_U133_PLUS2
GSE92538_U133_PLUS2_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE92538_U133_PLUS2_results/GSE92538_pheno_curated.csv")
mean(GSE92538_U133_PLUS2_pheno_tmp$POST_MORTEM_INTERVAL)  # 23.294 -> 23.3

# GSE102556
GSE102556_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE102556_results/GSE102556_pheno_curated.csv")
GSE102556_pheno_tmp = GSE102556_pheno_tmp[GSE102556_pheno_tmp$TISSUE == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)", ]
mean(GSE102556_pheno_tmp$PMI)  # 27.09375 -> 27.1

# GSE243356
GSE243356_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE243356_results/GSE243356_pheno_curated.csv")
# PMI is not reported in the dataset!
# Paper https://www.nature.com/articles/s41380-023-02311-9#Sec2 reports median PMI of 15

# GSE248260
GSE248260_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE248260_results/GSE248260_pheno_curated.csv")
mean(GSE248260_pheno_tmp$PMI) # 14.91667 -> 14.9

# GSE202537
GSE202537_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE202537_results/GSE202537_pheno_curated.csv")
# Nucleus accumbens (Nac)
GSE202537_pheno_tmp = GSE202537_pheno_tmp[GSE202537_pheno_tmp$TISSUE == "Nac",]
mean(GSE202537_pheno_tmp$PMI) # 17.4791 -> 17.5

# GSE101521
GSE101521_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE101521_results/GSE101521_pheno_curated.csv")
mean(GSE101521_pheno_tmp$PMI, na.rm=TRUE) # 14.48276 -> 14.5

# GSE144136
GSE144136_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE144136_results/GSE144136_pheno_curated.csv")
# PMI is not reported in the dataset!
# https://www.nature.com/articles/s41593-020-0621-y#Tab1 we recalculate it from table 1 presenting averaged data per group
mean(c(34.01, 41.69)) # 37.85 -> 37.9, Ncases=Ncontrols

# Better to recalculate from GSE213982 https://pmc.ncbi.nlm.nih.gov/articles/PMC10203145/#Tab1 as we used data from there
(34.01*16 + 41.69*17)/(17+16) # 37.96636 -> 38

# GSE213982
GSE213982_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE213982_results/GSE213982_pheno_curated.csv")
# PMI is not reported in the dataset!
# https://pmc.ncbi.nlm.nih.gov/articles/PMC10203145/#Tab1 we recalculate it from table 1 presenting averaged data per group
(41.49*20 + 30.27*18)/(20+18) # 36.17526 -> 36.2


################### Sex ratio calculations (Male/Female) or percent of male ###################

# Sex ratio is not possible since we have several cohorts where either only male or female are present
# Ratio calculations in such instances produce Inf for male-only datasetes

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



# GSE208338
GSE208338_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE208338_results/GSE208338_pheno_curated.csv")
tabulate_variable(GSE208338_pheno_tmp$SEX) # "Female: 46 (27.22%)\nMale: 123 (72.78%)\nNA: 0 (0%)"
123/46 # 2.673913 -> 2.7

# GSE5388
GSE5388_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE5388_results/GSE5388_pheno_curated.csv")
tabulate_variable(GSE5388_pheno_tmp$GENDER) # "Female: 20 (33.33%)\nMale: 40 (66.67%)\nNA: 0 (0%)"
40/20 # 2

# GSE5389
GSE5389_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE5389_results/GSE5389_pheno_curated.csv")
tabulate_variable(GSE5389_pheno_tmp$GENDER)
12/9 # 1.333333

# GSE66937
GSE66937_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE66937_results/GSE66937_pheno_curated.csv")
length(unique(GSE66937_pheno_tmp$INDIVIDUAL))
# Sex is not reported in the dataset!
# We calculate Sex based on Table 1 in https://pmc.ncbi.nlm.nih.gov/articles/PMC8458545/#Sec2
# It is approximate due to 20 sample vs available 15
GSE66937_pheno_tmp_scraped = read_html("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/HTMLs/Identification of transcriptome alterations in the prefrontal cortex, hippocampus, amygdala and hippocampus of suicide victims - PMC.html")
GSE66937_pheno_tmp_scraped = GSE66937_pheno_tmp_scraped %>% html_elements("table") %>% html_table(fill = TRUE)
GSE66937_pheno_tmp_scraped = GSE66937_pheno_tmp_scraped[[2]]
tabulate_variable(GSE66937_pheno_tmp_scraped$Gender) # "F: 9 (45%)\nM: 11 (55%)\nNA: 0 (0%)"
11/9 # 1.222222

# GSE199536
GSE199536_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE199536_results/GSE199536_pheno_curated.csv")
tabulate_variable(GSE199536_pheno_tmp$GENDER) # "male: 20 (100%)\nNA: 0 (0%)" Ratio can't be calculated

# GSE92538_U133A
GSE92538_U133A_pheno_tmp = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE92538_U133A_results/GSE92538_pheno_curated.csv")
tabulate_variable(GSE92538_U133A_pheno_tmp$GENDER) # "Female: 36 (37.11%)\nMale: 61 (62.89%)\nNA: 0 (0%)"
61/36 # 1.694444 -> 1.7

# GSE92538_U133_PLUS2
GSE92538_U133_PLUS2_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE92538_U133_PLUS2_results/GSE92538_pheno_curated.csv")
tabulate_variable(GSE92538_U133_PLUS2_pheno_tmp$GENDER) # "Female: 14 (18.67%)\nMale: 61 (81.33%)\nNA: 0 (0%)"
61/14

# GSE102556
GSE102556_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE102556_results/GSE102556_pheno_curated.csv")
GSE102556_pheno_tmp = GSE102556_pheno_tmp[GSE102556_pheno_tmp$TISSUE == "Dorsolateral prefrontal cortex (dlPFC; BA8/9)", ]
tabulate_variable(GSE102556_pheno_tmp$GENDER) # "FEMALE: 22 (45.83%)\nMALE: 26 (54.17%)\nNA: 0 (0%)"
26/22 # 1.181818 -> 1.2

# GSE243356
GSE243356_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE243356_results/GSE243356_pheno_curated.csv")
# Sex is not reported in the dataset!
# Paper https://www.nature.com/articles/s41380-023-02311-9#Sec2 reports 34% of women -> men is 66%
66/34 # 1.941176 -> 1.9

# GSE248260
GSE248260_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE248260_results/GSE248260_pheno_curated.csv")
tabulate_variable(GSE248260_pheno_tmp$SEX) # "Female: 8 (33.33%)\nMale: 16 (66.67%)\nNA: 0 (0%)"
16/8 # 2

# GSE202537
GSE202537_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE202537_results/GSE202537_pheno_curated.csv")
# Nucleus accumbens (Nac)
GSE202537_pheno_tmp = GSE202537_pheno_tmp[GSE202537_pheno_tmp$TISSUE == "Nac",]
tabulate_variable(GSE202537_pheno_tmp$GENDER) # "Female: 19 (28.36%)\nMale: 48 (71.64%)\nNA: 0 (0%)"
48/19 # 2.526316 -> 2.5

# GSE101521
GSE101521_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE101521_results/GSE101521_pheno_curated.csv")
tabulate_variable(GSE101521_pheno_tmp$SEX) # "Female: 17 (28.81%)\nMale: 42 (71.19%)\nNA: 0 (0%)"
42/17 # 2.470588 -> 2.5

# GSE144136
GSE144136_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE144136_results/GSE144136_pheno_curated.csv")
tabulate_variable(GSE144136_pheno_tmp$SEX) #  "Male: 32 (100%)\nNA: 0 (0%)"


# GSE213982
GSE213982_pheno_tmp  = read.csv("/home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis/GSE213982_results/GSE213982_pheno_curated.csv")
tabulate_variable(GSE213982_pheno_tmp$SEX) #  "Female: 36 (100%)\nNA: 0 (0%)"
