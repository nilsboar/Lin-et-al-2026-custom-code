library(MASS)

uc_europe <- readRDS("~/Documents/analysis/HLA/HLA 2019 on/UC/Europe.rds")
names(uc_europe)

table(uc_europe$SPsignal)
table(uc_europe$ethnicity)
names(uc_europe)[1:36]
table(substr(uc_europe$ID, 1, 2))
table(substr(overlaps$FID2, 1, 2))
table(nchar(uc_europe$ID))
overlaps <- read.table("~/Documents/analysis/HLA/HLA 2019 on/UC/list_ukibdgc_gwas1_gwas2_samples_in_ukb.tsv", header = T)
table(nchar(overlaps$FID2))
intersect(overlaps$FID2, uc_europe$ID)
write.csv(uc_europe$ID, file =  "~/Documents/analysis/HLA/HLA 2019 on/UC/IBDGC IDs.csv")


uc_bbk <- readRDS("~/Documents/analysis/HLA/HLA 2019 on/UC/UKB_UC.rds")
names(uc_bbk)
names(uc_bbk)[2] <- "sex"
for (i in 4:13) names(uc_bbk)[i] <- paste("PC", i-3, sep = "")
table(uc_bbk$SPsignal)
table(uc_bbk$ethnicity)
names(uc_bbk)[50:395]
cc_bbk <- complete.cases(uc_bbk[,c(50, 53:395)])
table(cc_bbk)
uc_data_uk <- uc_bbk[cc_bbk, ]
save(uc_data_uk, file = "~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/UC/complete.cases.IIDGBC.Rdata")
saveRDS(uc_data_uk, "~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/UC/complete.cases.UKB.rds")
table(uc_data_uk$case)
uc_data <- uc_data_uk

names(uc_europe)[38:358]
cc_eur <- complete.cases(uc_europe[,c(26:358)])
table(cc_eur)
uc_data_eur <- uc_europe[cc_eur, ]
saveRDS(uc_data_eur, "~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/UC/complete.cases.IIDGBC.rds")
table(uc_data_eur$case)

# 53:395 for biobank
uc_data <- uc_europe # glmnet doesn't need complete cases?  
sigallele <- vector(mode = "character")
sigalleledata <- data.frame()
# sigallele[1] <- "SPsignal"
j <- 1
for (i in 38:358)
# for (i in 53:395)
{
  if (dim(table(uc_data[,i])) > 1)
  {
    # smry <- summary(glm(uc_data[,"case"] ~ uc_data[,i], family = "binomial"))
    smry <- summary(glm(uc_data[,"case"] ~ uc_data[,i] + PC1 + PC2 + PC3 + PC4, data = uc_data, family = "binomial"))
    # print(c(names(uc_data)[i], round(smry$coefficients[2,1],3), signif(smry$coefficients[2,4],3)), sum(uc_data[,i]))
    if (smry$coefficients[2,4] < 0.01)
    {
      if (!(is.numeric(uc_data[,i])))
      {
        sigallele[j] <- names(uc_data)[i]
        print(c(names(uc_data)[i], round(smry$coefficients[2,1],3), signif(smry$coefficients[2,4],3)))
        # sigalleledata[j,1] <- as.data.frame(t(smry$coefficients[2,1]))
        # sigalleledata[j,2] <- as.data.frame(t(smry$coefficients[2,2]))
        # sigalleledata[j,3] <- as.data.frame(t(smry$coefficients[2,4]))
        # sigalleledata[j,4] <- names(uc_data)[i]
        j <- j+1
      }
      else if (sum(uc_data[,i], na.rm = T) > 19 & smry$coefficients[2,4] < 0.01)
      {
        sigallele[j] <- names(uc_data)[i]
        print(c(names(uc_data)[i], round(smry$coefficients[2,1],3), signif(smry$coefficients[2,4],3), sum(uc_data[,i], na.rm = T)))
        # sigalleledata[j,1] <- as.data.frame(t(smry$coefficients[2,1]))
        # sigalleledata[j,2] <- as.data.frame(t(smry$coefficients[2,2]))
        # sigalleledata[j,3] <- as.data.frame(t(smry$coefficients[2,4]))
        # sigalleledata[j,4] <- names(uc_data)[i]
        j <- j+1
      }
    }
  }
}
# write.csv(sigalleledata, file = "~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/UC/chosen eur.csv")
# sigallele <- sigallele[sample.int(length(sigallele))]
#####################################################################################################################################
# stepAIC
#####################################################################################################################################
# factors <- names(uc_europe)[37:358] # should be factors with p < 0.1. Or more stringent
# modelstart <- "case ~ PC1 + PC2 + PC3 + PC4 + sex + SPsignal"
modelstart <- "case ~ PC1 + PC2 + PC3 + PC4 + SPsignal"
# modelstart <- "case ~ PC1 + PC2 + PC3 + PC4 + sex + SPsignal"
lowerText <- " ~ PC1 + PC2 + PC3 + PC4" # or lower text is empty (just ~ ?)
# lowerText <- " ~ PC1 + PC2 + PC3 + PC4 + sex " # or lower text is empty (just ~ ?)
upperText <- " ~ PC1 + PC2 + PC3 + PC4 + SPsignal"
# upperText <- " ~ PC1 + PC2 + PC3 + PC4 + sex + SPsignal" # if SPsignal in sigalleles, otherwise add it
for (var in sigallele)
{
  upperText <- paste(upperText, " + ", var, sep = "")
}
# fullmodel <- as.formula(fullmodel)
# print(summary(glm(fullmodel, data = uc_data, family = "binomial")))
modelstart <- as.formula(modelstart)
uppersignal <- as.formula(upperText)
lowersignal <- as.formula(lowerText)
setk <- 2
startmdl <- glm(modelstart,  data = uc_data, family = "binomial")
startmdl <- glm(modelstart,  data = uc_data, family = binomial)
summary(startmdl)
bestAIC <- stepAIC(startmdl, scope = list(upper = uppersignal, lower = lowersignal), trace = 1, k = setk)
# bestAIC <- stepAIC(startmdl, scope = uppersignal, trace = 1, k = setk)
# str(bestAIC)
AICout <- summary(bestAIC)
print(AICout)
# write.csv(AICout$coefficients, file = "~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/UC/chosen UKB no sex.csv")

"DPB10601" %in% names(uc_data)
"DPB10601" %in% sigallele
#####################################################################################################################################
# bootstrap stepAIC
#####################################################################################################################################
nsubs <- nrow(uc_data)
nboot <- 100
table(uc_data$case)
bootsps_nosex <- matrix(nrow = nboot, ncol = 2, dimnames = list(NULL, c("coeff", "p")))
for (I in 1:nboot)
{
  t <- Sys.time()
  print(t)
  set.seed(as.integer(Sys.time()))
  boot_uc <- uc_data[sample.int(nsubs, replace = T),]
  startmdl <- glm(modelstart,  data = boot_uc, family = binomial)
  summary(startmdl)
  bestAIC <- stepAIC(startmdl, scope = list(upper = uppersignal, lower = lowersignal), trace = 0, k = setk)
  # bestAIC <- stepAIC(startmdl, scope = uppersignal, trace = 1, k = setk)
  # str(bestAIC)
  AICout <- summary(bestAIC)
  # print(AICout)
  # browser()
  if ("SPsignal" %in% rownames(AICout$coefficients))
  {
    bootsps_nosex[I, "coeff"] <- AICout$coefficients["SPsignal", 1]
    bootsps_nosex[I, "p"] <- AICout$coefficients["SPsignal", 4]
    print(c(I, AICout$coefficients["SPsignal", 4]))
  }
  else print(c(I, "SPsignal not chosen"))
  t <- Sys.time()
  print(t)
}
write.csv(bootsps_nosex, file = "~/Documents/analysis/HLA/HLA 2019 on/nk and ligands/UC/output/bootstrap UKBB no sex.csv")

#####################################################################################################################################
# glmnet
#####################################################################################################################################
# library(glmnet)
# # matrix of variables:
# pcs <- c("PC1", "PC2", "PC3", "PC4")
# varmat <- as.matrix(uc_europe[,c(pcs, sigallele)])
# outcome <- uc_europe[,"case"]
# 
# netout <- glmnet(varmat, outcome, "binomial")
# outcoeff <- netout[["beta"]]
# 
# str(outcoeff)
# outcoeff <- as.matrix(outcoeff)
# 













