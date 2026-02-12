# Input education data, call analysis functions
library(ggplot2)
library(lme4)
library(openxlsx)
library(tidyverse)
library(FactoMineR)
library(corrplot)
library(MASS)
library(pheatmap)
library(lattice)
library(viridis)

source("class I to SPs.R")

# load("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/R data/edu exps pls BLCL data.Rdata")
# load("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/R data/BLCL EU AA functional haps.Rdata")
load("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/R data/basic BLCL and monocyte data.Rdata")
# load("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/R data/lm3C.Rdata")
names(data363) <- gsub(".", "_", names(data363), fixed = T)
# for Mathias 
# rm(binders_lm, allsps, edu_230, edu.48pbmc, hla_spc, pbmc_hla, pbmc_sps)
# rm(data723)
# save.image("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/R data/edu exp data.Rdata")
# vlexp <- readRDS("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/R data/mVLexp.rds")

##########################################################################################
##########################################################################################
# Figure 2 analysis, 3 E to T ratios, CD107a % vs. Sp score
##########################################################################################
##########################################################################################
# names(data363)
# str(data363)
# data363$N_of_E0103 <- as.numeric(data363$N_of_E0103)
# lm363_ctE <- lm(NKG2A ~ SP_1A + SP_2A + SP_6B + SP_1C + SP_2C + N_of_E0103, data = data363)
# summary(lm363_ctE)

et3_47_1 <- read.xlsx("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/exp data/education data/47 PBMCs exp/47 PBMCs.xlsx")
et3_47 <- read.xlsx("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/exp data/education data/47 PBMCs exp/47 PBMCs2.xlsx")
names(et3_47) <- gsub(".", "_", names(et3_47), fixed = T)


et3_47 <- spvars(et3_47, binders)
et3_47$predNKG <- predict(lm363, et3_47)



#### plots alt predictions ##############
plot(et3_47$predNKG, et3_47$`CD107a_3:1`, xlab = "Prediction from SP score, BL", ylab = "CD107a response 10-1 ET", pch = 21, bg = "red")
abline(reg = lm(et3_47$`CD107a_3:1` ~ et3_47$predNKG), col = "black", lwd = 1.9)
abline(reg = lm(et3_47$`CD107a_3:1` ~ et3_47$predNKG), col = "red", lwd = 1.5)
# text(10, 75, "R2 = 0.42\nP = 1e-6")
plot(et3_47$predNKG, et3_47$`CD107a_10:1`, xlab = "Prediction from SP score, BL", ylab = "CD107a response 10-1 ET", pch = 21, bg = "green")
abline(reg = lm(et3_47$`CD107a_10:1` ~ et3_47$predNKG), col = "black", lwd = 1.9)
abline(reg = lm(et3_47$`CD107a_10:1` ~ et3_47$predNKG), col = "green", lwd = 1.5)
text(10, 75, "R2 = 0.42\nP = 1e-6")
plot(et3_47$predNKG, et3_47$`CD107a_30:1`, xlab = "Prediction from SP score, BL", ylab = "CD107a response 10-1 ET", pch = 21, bg = "blue")
abline(reg = lm(et3_47$`CD107a_30:1` ~ et3_47$predNKG), col = "black", lwd = 1.9)
abline(reg = lm(et3_47$`CD107a_30:1` ~ et3_47$predNKG), col = "blue", lwd = 1.5)
# text(10, 75, "R2 = 0.42\nP = 1e-6")

plot(et3_47$predNKG, et3_47$`CD107a_3:1`, ylim = c(15, 95), xlab = "SP score", ylab = "CD107a+ cells, (%)", pch = 21, bg = "red")
abline(reg = lm(et3_47$`CD107a_3:1` ~ et3_47$predNKG), col = "black", lwd = 1.9)
abline(reg = lm(et3_47$`CD107a_3:1` ~ et3_47$predNKG), col = "red", lwd = 1.5)
points(et3_47$predNKG, et3_47$`CD107a_10:1`,  pch = 21, bg = "green")
abline(reg = lm(et3_47$`CD107a_10:1` ~ et3_47$predNKG), col = "black", lwd = 1.9)
abline(reg = lm(et3_47$`CD107a_10:1` ~ et3_47$predNKG), col = "green", lwd = 1.5)
points(et3_47$predNKG, et3_47$`CD107a_30:1`,  pch = 21, bg = "blue")
abline(reg = lm(et3_47$`CD107a_30:1` ~ et3_47$predNKG), col = "black", lwd = 1.9)
abline(reg = lm(et3_47$`CD107a_30:1` ~ et3_47$predNKG), col = "blue", lwd = 1.5)



#####################################################################################
# for PBMC SP score
# generate new score adding each non-functional SP (need data363)
# calculate the effect on education, for each of the 3 E to T ratios 
# (change code by hand for the three)
#####################################################################################
nonbinders <- allSPs[!(allSPs %in% binders)]
data363 <- spvars(data363, nonbinders)
et3_47 <- spvars(et3_47, nonbinders)
names(data363)
names(et3_47)
c("`CD107a_3:1`", "`CD107a_10:1`", "`CD107a_30:1`")
setwd("~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/writing/response to reviewers 2/more sps, 3C")
predsp <- c("predNKG", "SP_3A")
pred3 <- data.frame(row.names = predsp)
for (add_sp in nonbinders)
{
  model <- paste(bindermodel, "+", add_sp)
  nblm <- lm(as.formula(model), data = data363)
  et3_47$pred_nb <- predict(nblm, et3_47)
  rslt <- summary(lm(`CD107a_10:1` ~ pred_nb, data = et3_47))
  pred3[add_sp, "R2"] <- signif(rslt$r.squared, 2)
  pred3[add_sp, "P"] <- signif(rslt$coefficients[2,4], 2)
}
rslt <- summary(lm(`CD107a_10:1` ~ predNKG, data = et3_47))
pred3["predNKG", "R2"] <- signif(rslt$r.squared, 2)
pred3["predNKG", "P"] <- signif(rslt$coefficients[2,4], 2)

write.csv(pred3, file = paste("BLCL predicted for E to T ", "CD107a_10:1", ".csv", sep = ""))



# et3_47$predNKG_ctE <- predict(lm363_ctE, et3_47)
# et3_47$E3D12 <- et3_47$`3D12.MFI.on.PBMCs`
# et3_47$predNKG_3d12 <- predict(lm363_3d12, et3_47)


#########################################################
######## reproducing plot 10_3_25 ################
names(et3_47)
et3_47$E3D12 <- et3_47$`3D12_MFI_on_PBMCs` #. `3D12_MFI_on_PBMCs` holds original 3D12 value
et3_47$predNKG_3d12_new <- predict(lm363_3d12_new, et3_47)
et3_47$predNKG_3d12_old <- predict(lm363_3d12_old, et3_47)
hist(et3_47$E3D12)
hist(data363$E3D12)
######## rescaling et3_47$E3D12 
mean_blcl <- mean(data363$E3D12, na.rm = T)
std_blcl <- sd(data363$E3D12, na.rm = T)
mean_et3 <- mean(et3_47$E3D12, na.rm = T)
std_et3 <- sd(et3_47$E3D12, na.rm = T)
# et3_47$E3D12rscld <- et3_47$E3D12 - mean_et3
# et3_47$E3D12rscld <- et3_47$E3D12rscld*std_blcl/std_et3
# et3_47$E3D12rscld <- et3_47$E3D12rscld + mean_blcl
# but this rescaling gives weak association, why?
et3_47$E3D12rscld <- et3_47$`3D12_MFI_on_PBMCs` + mean_blcl - mean_et3

hist(et3_47$E3D12rscld, xlim = c(0, 50))
hist(data363$E3D12, xlim = c(0, 50))
# for prediction name of E level must be E3D12 to match BLCL vbls
et3_47$E3D12 <- et3_47$E3D12rscld
et3_47$pred_nkg_scld <- predict(lm363_3d12_old, et3_47)


plot(et3_47$pred_nkg_scld, et3_47$`CD107a_10:1`, xlab = "Prediction from SP score + 3D12", ylab = "CD107a response 10-1 ET")
summary(lm(et3_47$`CD107a_10:1` ~ et3_47$pred_nkg_scld))
abline(reg = lm(et3_47$`CD107a_10:1` ~ et3_47$pred_nkg_scld))
text(8.5, 79.5, expression(R^2 * "= 0.45"))
text(8.7, 76, "p = 2.8e-7")

plot(et3_47$predNKG, et3_47$`CD107a_10:1`, xlab = "Prediction from SP score", ylab = "CD107a response 10-1 ET")
summary(lm(et3_47$`CD107a_10:1` ~ et3_47$predNKG))
abline(reg = lm(et3_47$`CD107a_10:1` ~ et3_47$predNKG))
text(9.5, 79.5, expression(R^2 * "= 0.42"))
text(9.47, 76, "p = 1e-6")




############# Lin's 47, more prediction from %Bright etc, and %NK as covar ###########

et3_47$bright_in_pbmcs <- et3_47$`%NK`*et3_47$`%Bright`/100
et3_47$dim_in_pbmcs <- et3_47$`%NK`*(100 - et3_47$`%Bright`)/100

summary(lm(`CD107a_10:1` ~ `%NK`, data = et3_47))
summary(lm(`CD107a_30:1` ~ `%NK`, data = et3_47))

summary(lm(`CD107a_10:1` ~ `%Bright`, data = et3_47))
summary(lm(`CD107a_30:1` ~ `%Bright`, data = et3_47))

summary(lm(`CD107a_10:1` ~ bright_in_pbmcs, data = et3_47))
summary(lm(`CD107a_30:1` ~ bright_in_pbmcs, data = et3_47))

summary(lm(`CD107a_10:1` ~ dim_in_pbmcs, data = et3_47))
summary(lm(`CD107a_30:1` ~ dim_in_pbmcs, data = et3_47))

summary(lm(`CD107a_3:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_3:1` ~ predNKG + `%NK`, data = et3_47))
summary(lm(`CD107a_10:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_10:1` ~ predNKG + `%NK`, data = et3_47))
summary(lm(`CD107a_30:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_30:1` ~ predNKG + `%NK`, data = et3_47))

summary(lm(`CD107a_10:1` ~ bright_in_pbmcs, data = et3_47))
summary(lm(`CD107a_10:1` ~ dim_in_pbmcs, data = et3_47))


# does E expression correlate with N of E0103?
# E expression as a covariate; PBMCs and NKs

summary(lm(`3D12.MFI.on.PBMCs` ~ N.of.E0103, data = et3_47))
summary(lm(`3D12.MFI.on.NK.cells` ~ N.of.E0103, data = et3_47))

summary(lm(`CD107a_10:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_10:1` ~ predNKG + `3D12.MFI.on.PBMCs`, data = et3_47))
summary(lm(`CD107a_10:1` ~ predNKG + `3D12.MFI.on.NK.cells`, data = et3_47))

summary(lm(`CD107a_30:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_30:1` ~ predNKG + `3D12.MFI.on.PBMCs`, data = et3_47))
summary(lm(`CD107a_30:1` ~ predNKG + `3D12.MFI.on.NK.cells`, data = et3_47))


summary(lm(`CD107a_3:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_3:1` ~ predNKG + `3D12.MFI.on.PBMCs`, data = et3_47))

summary(lm(`CD107a_3:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_3:1` ~ predNKG + `3D12.MFI.on.NK.cells`, data = et3_47))

summary(lm(`CD107a_10:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_10:1` ~ predNKG_3d12, data = et3_47))

plot(et3_47$predNKG_3d12, et3_47$`CD107a_10:1`)

CMV_pls <- subset(et3_47, CMV == 1)
CMV_mns <- subset(et3_47, CMV == 0)
et3_47$CMV_status <- as.character(et3_47$CMV)
plot(et3_47$predNKG, et3_47$`CD107a_10:1`, xlab = "Prediction from SP score", ylab = "CD107a, ET 10 to 1")
points(CMV_pls$predNKG, CMV_pls$`CD107a_10:1`, col = "red", pch = 16)
points(CMV_mns$predNKG, CMV_mns$`CD107a_10:1`, col = "blue", pch = 16)
abline(reg = lm(`CD107a_10:1` ~ predNKG, data = et3_47))
abline(reg = lm(`CD107a_10:1` ~ predNKG, data = CMV_pls), col = "red")
abline(reg = lm(`CD107a_10:1` ~ predNKG, data = CMV_mns), col = "blue")
text(10, 78, "R2 = 0.49", col = "red")
text(10, 73, "R2 = 0.17", col = "blue")

summary(lm(`CD107a_10:1` ~ predNKG, data = et3_47))
summary(lm(`CD107a_10:1` ~ predNKG, data = CMV_pls))
summary(lm(`CD107a_10:1` ~ predNKG, data = CMV_mns))
summary(lm(`CD107a_30:1` ~ predNKG, data = CMV_pls))
summary(lm(`CD107a_30:1` ~ predNKG, data = CMV_mns))
summary(lm(`CD107a_10:1` ~ predNKG*CMV, data = et3_47))
summary(lm(`CD107a_10:1` ~ CMV, data = et3_47))
summary(lm(CMV ~ predNKG, data = et3_47))


delr2 <- 0.4931 - 0.1715
delr2 <- 0.3715 - 0.2634
table(et3_47$CMV)
frac <- 25/46
nperm <- 10000
ct_moresig <- 0
for (i in 1:nperm)
{
  et3_47$rand_cmv <- rbinom(47, 1, frac)
  CMV_pls <- subset(et3_47, rand_cmv == 1)
  CMV_mns <- subset(et3_47, rand_cmv == 0)
  sumpls <- summary(lm(`CD107a_30:1` ~ predNKG, data = CMV_pls))
  summns <- summary(lm(`CD107a_30:1` ~ predNKG, data = CMV_mns))
  diff <- sumpls$r.squared - summns$r.squared
  if (abs(diff) > abs(delr2)) ct_moresig <- ct_moresig + 1
}
ct_moresig/nperm




str(sumpls)

# p = 0.08 for 10:1; = 0.06 for 30:1



ggplot(subset(et3_47, !is.na(CMV)), aes(x = CMV_status, y = `CD107a_10:1`)) +
  geom_jitter(width = 0.1)
t.test(CMV_pls$`CD107a_10:1`, CMV_mns$`CD107a_10:1`)

ggplot(subset(et3_47, !is.na(CMV)), aes(x = CMV_status, y = `CD107a_30:1`)) +
  geom_jitter(width = 0.1)
t.test(CMV_pls$`CD107a_30:1`, CMV_mns$`CD107a_30:1`)




##########################################################################################
# auto selection
##########################################################################################
outcome <- "pctCD107a"
vars <- binders[c(3,1,2,4,5)]
# vars <- c("pctNK", vars)
modelstart <- paste(outcome ," ~ ", vars[1])
upperText <- lowerText <- paste(" ~ ", vars[1])
for (i in 2:length(vars))
{
  upperText <- paste(upperText, " + ", vars[i], sep = "")
}
modelstart <- as.formula(modelstart)
uppersignal <- as.formula(upperText)
lowersignal <- as.formula(lowerText)
setk <- 2
startmdl <- lm(modelstart,  data = edu_pbmc2)
bestAIC <- stepAIC(startmdl, scope = list(upper = uppersignal, lower = lowersignal), trace = 1, k = setk)
AICout <- summary(bestAIC)
AICout
##########################################################################################
# all SP interactions
##########################################################################################
# allinter <- "(SP_1A + SP_3A + SP_4A + SP_2A + SP_5A + SP_6B + SP_3B + SP_1B + SP_2B + SP_5B + SP_4B + SP_7B + SP_3C + SP_1C + SP_2C + SP_4C + SP_5C + SP_45C):(SP_1A + SP_3A + SP_4A + SP_2A + SP_5A + SP_6B + SP_3B + SP_1B + SP_2B + SP_5B + SP_4B + SP_7B + SP_3C + SP_1C + SP_2C + SP_4C + SP_5C + SP_45C)"
# upperText <- paste(upperText, "+", allinter, sep = " ")
##########################################################################################
# simple run; with or without interactions
##########################################################################################
modelstart <- as.formula(modelstart)
uppersignal <- as.formula(upperText)
lowersignal <- as.formula(lowerText)
setk <- 2
startmdl <- lm(modelstart,  data = edu_pbmc2)
bestAIC <- stepAIC(startmdl, scope = list(upper = uppersignal, lower = lowersignal), trace = 1, k = setk)
AICout <- summary(bestAIC)
AICout


#########################################################################################################
# exp 012 mixed model
#########################################################################################################

summary(lmer(pctCD107a ~ predNKG2A + (1|PID), data = edu.exp012))
summary(lmer(pctCD107a ~ predNKG2A + (1|exp), data = edu.exp012))
summary(lmer(pctCD107a ~ predNKG2A + (1|PID) + (1|exp), data = edu.exp012))

summary(lm(pctCD107a ~ predNKG2A, data = edu.exp012))
edu.exp01 <- subset(edu.exp012, exp != 2)
summary(lm(pctCD107a ~ predNKG2A, data = edu.exp01))
summary(lm(pctCD107a ~ predNKG2A, data = edu_pbmc2))





#########################################################################################################
# factors vs factors
#########################################################################################################

lm1 <- summary(lm(CD107a_pls ~ D14_purity, data = edu.exp1))
lm2 <- summary(lm(CD107a_pls ~ D14_purity, data = edu.exp2))

lm1 <- summary(lm(D14_purity ~ pctCD16_pls, data = edu.exp1))
lm2 <- summary(lm(D14_purity ~ pctCD16_pls, data = edu.exp2))





samesubs <- merge(edu.exp1, edu.exp2, by = "PID")
plot(samesubs$CD107a_pls.x, samesubs$CD107a_pls.y)
abline(lm(CD107a_pls.y ~ CD107a_pls.x, data = samesubs))
summary(lm(CD107a_pls.y ~ CD107a_pls.x, data = samesubs))
text(30,80, "R2 = 0.18")
plot(samesubs$cd107_corNK.x, samesubs$cd107_corNK.y)
abline(lm(cd107_corNK.y ~ cd107_corNK.x, data = samesubs))
summary(lm(cd107_corNK.y ~ cd107_corNK.x, data = samesubs))
text(57,80, "R2 = 0.10")
plot(samesubs$NKxCD16.x, samesubs$NKxCD16.y)
abline(lm(NKxCD16.y ~ cd107_corNK.x, data = samesubs))
summary(lm(cd107_corNK.y ~ cd107_corNK.x, data = samesubs))
text(57,80, "R2 = 0.10")


#########################################################################################################
# pca and correlation
#########################################################################################################
exp1mat <- as.matrix(edu.exp1in[,c(9,12:17)])
pca.exp1mat <- PCA(exp1mat)
exp2mat <- as.matrix(edu.exp2in[,c(9:17)])
pca.exp1mat <- PCA(exp2mat)

cor_matrix1 <- cor(exp1mat)
corrplot(cor_matrix, method = "color")
cor_matrix2 <- cor(exp2mat)
corrplot(cor_matrix, method = "color")
pheatmap(cor_matrix2)
wéq3


#########################################################################################################
# compare repeated subs, get BLCL prediction, basic analysis
#########################################################################################################
names(edu.exp1)
# compare data with same subs now
# bring in covariates
plot(samesubs$CD107a_pls.x, samesubs$CD107a_pls.y, xlab = "exp1 %CD107a", ylab = "exp2 %CD107a", main = "Exp 2 vs. Exp 1, 26 shared samples")
text(30, 80, "R2 = 0.18")
summary(lm(samesubs$CD107a_pls.x ~ samesubs$CD107a_pls.y))

source("~/Documents/analysis/R/R code and data/HLA R/HLA R 2019 on/signal peptides/SP functions/class I to SPs.R")
edu.exp1 <- spvars(edu.exp1, binders)
edu.exp1$predNKG2A <- predict(lm363, edu.exp1)
# Lin's vs. my BLCL score
plot(edu.exp1$predNKG2A, edu.exp1$BLCL363_score, xlab = "my predicted NKG2A recognition", ylab = "Lin's BLCL363 score", pch = 21,  cex = 2, cex.axis = 1.5, cex.lab = 1.4, mar = c(5, 8, 4, 2),  cex.main = 2, bg = "gray")

edu.exp2 <- spvars(edu.exp2, binders)
edu.exp2$predNKG2A <- predict(lm363, edu.exp2)
plot(edu.exp2$predNKG2A, edu.exp2$BLCL363_score, xlab = "my predicted NKG2A recognition", ylab = "Lin's BLCL363 score", pch = 21,  cex = 2, cex.axis = 1.5, cex.lab = 1.4, mar = c(5, 8, 4, 2),  cex.main = 2, bg = "gray")

