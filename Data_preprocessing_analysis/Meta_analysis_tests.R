
############ Testing package ############

dat <- escalc(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg)
dat

dat.long <- to.long(measure="OR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat.bcg, append=FALSE)


rma(yi, vi, data=dat)
rma.mv(yi, vi, random = ~ 1 | trial, data=dat)


levels(dat.long$group) <- c("exp", "con")
dat.long$group <- relevel(dat.long$group, ref="con")

dat.long <- escalc(measure="PLO", xi=out1, mi=out2, data=dat.long)
dat.long


############ Testing rma.mv ############

## rma.mv should be used
dat$trial2 = dat$trial + sample(c(0,1), size = 13, prob = c(.5,.5), replace = TRUE)
model = rma.mv(yi, vi, random = ~ 1 | trial2, data=dat)

# looks like working...

# 1 RMA MV with study effect
# 1 Obtain average log2FC and average SE

test_data = smart_fread("Data_preprocessing_analysis/Eug_combined_top_tables_with_gse.csv")
test_data_genes_table = as.data.frame(table(test_data$Gene_Name))

############ Simple rma.uni (incorrect) ############
test_data_ABCF1 = test_data[test_data$Gene_Name %in% "ABCF1", ]
test_data_ABCF1$X = NULL
test_data_ABCF1 = distinct(test_data_ABCF1, .keep_all = TRUE)
col_vector = generate_palette(test_data_ABCF1$GSE)
test_data_ABCF1$color = sapply(test_data_ABCF1$GSE, function(x){return(col_vector[x])})

meta_model = rma.uni(logFC, SE^2, data=test_data_ABCF1)
meta_model$ci.lb
meta_model$ci.ub
# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights, colout = test_data_ABCF1$color)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

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
                                          tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                          chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                          ", df=", .(meta_model$k - meta_model$p), ", ",
                                          .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)))))
title(paste0("Gene: ", "ABCF1"))


############ Testing rma.mv on test data ############

test_data_ABCF1 = test_data[test_data$Gene_Name %in% "ABCF1", ]
test_data_ABCF1$X = NULL
test_data_ABCF1 = distinct(test_data_ABCF1, .keep_all = TRUE)
col_vector = generate_palette(test_data_ABCF1$GSE)
test_data_ABCF1$color = sapply(test_data_ABCF1$GSE, function(x){return(col_vector[x])})

meta_model = rma.mv(logFC, SE^2, random = ~ 1 | GSE, data=test_data_ABCF1)
meta_model$ci.lb
meta_model$ci.ub
# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights, colout = test_data_ABCF1$color)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

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
                                            tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                            ", df=", .(meta_model$k - meta_model$p), ", ",
                                            .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)))))
title(paste0("Gene: ", "ABCF1"))

# kind of working

############ Testing rma.uni on averaged data ############

compute_mean_SE = function(x){
  
  if (length(x) == 1){
    return(x)
  }
  
  pooled_SE = sqrt(sum(x**2)/length(x))
  return(pooled_SE)
}


GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = compute_mean_SE(SE))
# 0.20840432

meta_model = rma.uni(meanlog2FC, meanSE^2, data=test_data_ABCF1_aver)
meta_model$ci.lb
meta_model$ci.ub

# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

# Add text
par(xpd=NA)
par(cex=sav$cex, font=2)
# Headers
text(sav$xlim[1], k+2.5, pos=4, "Cohort")
text(-0.58, k+2.5, pos=4, "Weight %")
text(0, k+2.7, "Log2FC,\n(95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")

# Use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)
text(sav$ilab.xpos[3], 0, "100.0")
text(sav$xlim[1], -0.5, pos=4, bquote(paste("Test for heterogeneity: ",
                                            tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                            ", df=", .(meta_model$k - meta_model$p), ", ",
                                            .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)), "; ",
                                            I^2, "=", .(round(meta_model$I2)), "%")))
title(paste0("Gene: ", "ABCF1"))



############ Testing rma.uni on averaged data (approach 2) ############
# Approach 2

compute_mean_SE2 = function(x){
  
  if (length(x) == 1){
    return(x)
  }
  
  pooled_SE = x**2
  pooled_SE = 1/pooled_SE
  pooled_SE = sum(pooled_SE)
  pooled_SE = sqrt(1/pooled_SE)
  
  return(pooled_SE)
}


GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = compute_mean_SE2(SE))
# 0.03091282


meta_model = rma.uni(meanlog2FC, meanSE^2, data=test_data_ABCF1_aver)
meta_model$ci.lb
meta_model$ci.ub

# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

# Add text
par(xpd=NA)
par(cex=sav$cex, font=2)
# Headers
text(sav$xlim[1], k+2.5, pos=4, "Cohort")
text(-0.58, k+2.5, pos=4, "Weight %")
text(0, k+2.7, "Log2FC,\n(95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")

# Use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)
text(sav$ilab.xpos[3], 0, "100.0")
text(sav$xlim[1], -0.5, pos=4, bquote(paste("Test for heterogeneity: ",
                                            tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                            ", df=", .(meta_model$k - meta_model$p), ", ",
                                            .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)), "; ",
                                            I^2, "=", .(round(meta_model$I2)), "%")))
title(paste0("Gene: ", "ABCF1"))


############ Approach 3 ############

compute_mean_SE = function(x){
  
  if (length(x) == 1){
    return(x)
  }
  
  pooled_SE = sqrt(sum(x**2))/length(x)
  return(pooled_SE)
}


GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = compute_mean_SE(SE))

# 0.03684103

meta_model = rma.uni(meanlog2FC, meanSE^2, data=test_data_ABCF1_aver)
meta_model$ci.lb
meta_model$ci.ub

# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

# Add text
par(xpd=NA)
par(cex=sav$cex, font=2)
# Headers
text(sav$xlim[1], k+2.5, pos=4, "Cohort")
text(-0.58, k+2.5, pos=4, "Weight %")
text(0, k+2.7, "Log2FC,\n(95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")

# Use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)
text(sav$ilab.xpos[3], 0, "100.0")
text(sav$xlim[1], -0.5, pos=4, bquote(paste("Test for heterogeneity: ",
                                            tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                            ", df=", .(meta_model$k - meta_model$p), ", ",
                                            .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)), "; ",
                                            I^2, "=", .(round(meta_model$I2)), "%")))
title(paste0("Gene: ", "ABCF1"))

############ Approach 4 ############

compute_mean_SE = function(x){
  
  if (length(x) == 1){
    return(x)
  }
  
  pooled_SE = sum(x**2)
  pooled_SE = pooled_SE/length(x)**2
  pooled_SE = sqrt(pooled_SE)

  return(pooled_SE)
}


GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = compute_mean_SE(SE))

# 0.03684103


############ Approach 5 ############

# Not tested

GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = mean(SE))

# 0.19813110

############ Approach 6 (Lars, max SE) ############

test_data = smart_fread("Data_preprocessing_analysis/Eug_combined_top_tables_with_gse.csv")
test_data_genes_table = as.data.frame(table(test_data$Gene_Name))

test_data_ABCF1 = test_data[test_data$Gene_Name %in% "ABCF1", ]
test_data_ABCF1$X = NULL
test_data_ABCF1 = distinct(test_data_ABCF1, .keep_all = TRUE)

GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = max(SE))

# 0.03684103

meta_model = rma.uni(meanlog2FC, meanSE^2, method = "REML", control=list(maxiter=10000),data=test_data_ABCF1_aver)
meta_model$ci.lb
meta_model$ci.ub

# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

# Add text
par(xpd=NA)
par(cex=sav$cex, font=2)
# Headers
text(sav$xlim[1], k+2.5, pos=4, "Cohort")
text(-0.58, k+2.5, pos=4, "Weight %")
text(0, k+2.7, "Log2FC,\n(95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")

# Use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)
text(sav$ilab.xpos[3], 0, "100.0")
text(sav$xlim[1], -0.5, pos=4, bquote(paste("Test for heterogeneity: ",
                                            tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                            ", df=", .(meta_model$k - meta_model$p), ", ",
                                            .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)), "; ",
                                            I^2, "=", .(round(meta_model$I2)), "%")))
title(paste0("Gene: ", "ABCF1"))



############ Approach 6 (Lars, Fixed effect) ############

test_data = smart_fread("Eug_combined_top_tables_with_gse.csv")
test_data_genes_table = as.data.frame(table(test_data$Gene_Name))

test_data_ABCF1 = test_data[test_data$Gene_Name %in% "ABCF1", ]
test_data_ABCF1$X = NULL
test_data_ABCF1 = distinct(test_data_ABCF1, .keep_all = TRUE)


compute_mean_SE = function(x){
  
  if (length(x) == 1){
    return(x)
  }
  
  pooled_SE = sqrt(sum(x**2))/length(x)
  return(pooled_SE)
}


GSE9385_test = test_data_ABCF1[test_data_ABCF1$GSE == "GSE9385", ]
test_data_ABCF1_aver = test_data_ABCF1 %>%
  group_by(GSE)  %>%
  summarise(meanlog2FC = mean(logFC), meanSE = max(SE))

# 0.03684103

meta_model = rma.uni(meanlog2FC, meanSE^2, method = "EE", data=test_data_ABCF1_aver)
meta_model$ci.lb
meta_model$ci.ub

# Example Visualization
weights = fmtx(weights(meta_model), digits=1)
sav = forest(meta_model, slab = GSE, ilab = weights)
k = nrow(test_data_ABCF1)
colp <- "red"
segments(coef(meta_model), 0, coef(meta_model), k, col=colp, lty="33", lwd=0.8)

# Add text
par(xpd=NA)
par(cex=sav$cex, font=2)
# Headers
text(sav$xlim[1], k+2.5, pos=4, "Cohort")
text(-0.58, k+2.5, pos=4, "Weight %")
text(0, k+2.7, "Log2FC,\n(95% CI)")
segments(sav$ilab.xpos[1]-0.22, k+2.8, sav$ilab.xpos[2]+0.13, k+2.8)
text(sav$xlim[2]-0.10, k+2.7, "Log2FC\n(95% CI)")

# Use a non-bold font for the rest of the text
par(cex=sav$cex, font=1)
text(sav$ilab.xpos[3], 0, "100.0")
text(sav$xlim[1], -0.5, pos=4, bquote(paste("Test for heterogeneity: ",
                                            tau^2, "=", .(fmtx(meta_model$tau2, digits=2)), "; ",
                                            chi^2, "=", .(fmtx(meta_model$QE, digits=2)),
                                            ", df=", .(meta_model$k - meta_model$p), ", ",
                                            .(fmtp(meta_model$QEp, digits=2, pname="P", add0=TRUE, equal=TRUE)), "; ",
                                            I^2, "=", .(round(meta_model$I2)), "%")))
title(paste0("Gene: ", "ABCF1"))


