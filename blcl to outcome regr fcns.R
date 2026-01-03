require(MASS)
# Create dataframe with data for each signal peptide genotype
spgt <- function(gtdata, sps, morevars = NULL, datacols = NULL) # subset return to needed data
{
  # browser()
  # vars <- c("SP_1A", "SP_2A", "SP_6B", "SP_1C", "SP_2C", morevars)
  vars <- c(sps, morevars)
  # vars <- binders
  bgt_gtdata <- unique(gtdata[,vars])
  # for (i in 1:nrow(bgt_gtdata)) row.names(bgt_gtdata)
  row.names(bgt_gtdata) <- 1:nrow(bgt_gtdata)
  spmatches <- list()
  for (i in 1:nrow(bgt_gtdata))
  {
    spmatches_gtdata <- gtdata
    # browser()
    gtname <- ""
    for (sp in vars) # "sp" may not be an sp
    {
      spmatches_gtdata <- subset(spmatches_gtdata, get(sp) == bgt_gtdata[i, sp])
      if(is.na(bgt_gtdata[i, sp])) browser()
      if (bgt_gtdata[i, sp] > 0) 
      {
        for (j in 1:bgt_gtdata[i, sp]) gtname <- paste(gtname, "_", substr(sp,4,5), sep = "")
      }
    }
    bgt_gtdata[i,"name"] <- gtname
    bgt_gtdata[i,"count"] <- nrow(spmatches_gtdata)
    # browser()
    for (col in datacols)
    {
      meanname <- paste("mean_", col, sep = "")
      variancename <- paste("var_", col, sep = "")
      bgt_gtdata[i,meanname] <- mean(spmatches_gtdata[,col], na.rm = T)
      bgt_gtdata[i,variancename] <- var(spmatches_gtdata[,col], na.rm = T)
    }
    # output for brights, dims
    # bgt_gtdata[i,"avall"] <- mean(spmatches_gtdata$all)
    # bgt_gtdata[i,"avbright"] <- mean(spmatches_gtdata$bright)
    # bgt_gtdata[i,"avdim"] <- mean(spmatches_gtdata$dim)
    # bgt_gtdata[i,"weight"] <- sum(spmatches_gtdata$weight)
    # bgt_gtdata[i,"ctB38"] <- nrow(subset(spmatches_gtdata, B3801 == 1))
    # if (nrow(spmatches_gtdata > 5)) browser()
    # print(nrow(spmatches_gtdata))
    spmatches[[i]] <- spmatches_gtdata
  }
  # return(list(bgt_gtdata, spmatches)) #n.b. RETURNING LIST
  return(bgt_gtdata)
}

predict_nkg_auto <- function(dataset, outcome, fixed_set, auto_set)
{
  start_text <- paste(outcome," ~ ", fixed_set[1])
  lowertext <- paste(" ~ ", fixed_set[1])
  uppertext <- paste(" ~ ", fixed_set[1])
  for (i in 2:length(fixed_set))
  {
    start_text <- paste(start_text, " + ", fixed_set[i], sep = "")
    lowertext <- paste(lowertext, " + ", fixed_set[i], sep = "")
    uppertext <- paste(uppertext, " + ", fixed_set[i], sep = "")
  }
  for (i in 1:length(auto_set))
  {
    uppertext <- paste(uppertext, " + ", auto_set[i], sep = "")
  }
  modelstart <- as.formula(start_text)
  uppersignal <- as.formula(uppertext)
  lowersignal <- as.formula(lowertext)
  setk <- 2
  startmdl <- lm(modelstart,  data = dataset)
  AICout <- stepAIC(startmdl, scope = list(upper = uppersignal, lower = lowersignal), trace = 1, k = setk)
  # browser()
  chosen_mdl <- AICout$call
  lm_out_binders_inter <- eval(chosen_mdl)
  # browser()
  return(list(startmdl, lm_out_binders_inter))
}


predict_outcome <- function(outcome_data, outcome, spset, covars_blcls = NULL, covar_outcome = NULL, blcl_data = NULL, blcl_response = NULL, coeff_frame = NULL, outcome_vs_sps = T, plottitle = NULL) 
# "blcl" stands for any data giving nkg2a response to SPs
# add weights.  need to edit function for multiple outcome covars.
{
  #######################################################################################################################
  # direct regression with signal peptides
  #######################################################################################################################
  # browser()
  # Regression predicting outcome from set of SPs carried
  # Distinct from predicting outcome from NKG2A score for SPs carried, determined from e.g. BLCL NKG2A response
  if (outcome_vs_sps)
  {
    model <- paste(outcome," ~ ", spset[1])
    for (i in 2:length(spset))
    {
      model <- paste(model, " + ", spset[i], sep = "")
    }
    if(is.null(coeff_frame)) # if this is inside a loop, coefficients frame is passed and updated
    {
      outcome_lm <- summary(lm(as.formula(model), data = outcome_data))
      coeff_frame <- as.data.frame(outcome_lm$coefficients)
      coeff_frame$run <- "edu_regr"
      coeff_frame$SP <- row.names(coeff_frame)
      coeff_frame <- coeff_frame[,1:6]
    }
    if (is.null(blcl_data)) return(list(coeff_frame, outcome_lm)) # just testing prediction of outcome from binders 
  }
  #######################################################################################################################
  # append BLCL coefficients for binders
  #######################################################################################################################
  blcl_lm <- blcl_regr(blcl_data, blcl_response, spset, covars_blcls)
  more_coeffs <- as.data.frame(summary(blcl_lm)$coefficients)
  more_coeffs$run <- blcl_response #  need real groupname
  more_coeffs$SP <- row.names(more_coeffs)
  if(is.null(coeff_frame)) coeff_frame <- more_coeffs
  else coeff_frame <- rbind(coeff_frame, more_coeffs) 
  #######################################################################################################################
  # predict edu from blcl coeffs
  #######################################################################################################################
  outcome_data$predict_cd69 <- predict(blcl_lm, outcome_data)
  # browser()
  # "edu" should be replaced by something more generic?
  if (is.null(covar_outcome))
  {
    edu_lm <- lm(outcome_data[,outcome] ~ outcome_data$predict_cd69, weights = outcome_data$weight)
    edu_smry <- summary(edu_lm)
    plottitle <- paste (plottitle, "R2 =", signif(edu_smry$r.squared, 3))
    print("test, plotting")
    plot(outcome_data$predict_cd69, outcome_data[,outcome], xlab = "Prediction from BLCL data", ylab = "%CD107a positive", main = plottitle, pch = 21, bg = "gray")
    abline(reg = edu_lm)
  } else 
  {
    edu_lm <- lm(outcome_data[,outcome] ~ outcome_data$predict_cd69 + outcome_data[,covar_outcome], weights = outcome_data$weight)
    edu_smry <- summary(edu_lm)
    plottitle <- paste (plottitle, "R2 =", signif(edu_smry$r.squared, 3))
    plot(outcome_data$predict_cd69, outcome_data[,outcome], xlab = "Prediction from BLCL data", ylab = "%CD107a positive", main = plottitle, pch = 21, bf = "gray")
    abline(reg = edu_lm)
  }
  return(list(blcl_lm, outcome_lm, edu_lm, outcome_data, coeff_frame))
}

blcl_regr <- function(blcl_data, blcl_response, spset, covars = NULL)
{
  # browser()
  model <- paste(blcl_response ," ~ ", spset[1])
  for (i in 2:length(spset))
  {
    model <- paste(model, " + ", spset[i], sep = "")
  }
  if (!is.null(covars))
  {
    for (i in 1:length(covars)) model <- paste(model, " + ", covars[i], sep = "")
    maintitle <- paste("BLCL CD69 response vs. BLCL prediction\n from SPs and ", covars[1], ", R2 = ", sep = "")
  }
  else maintitle <- paste("BLCL CD69 response vs. BLCL prediction\n from SPs, R2 =")
  blcl_lm <- lm(as.formula(model), data = blcl_data)
  # summary(blcl_lm)
  # FOLLOWING gives a plot of self prediction of SP + for BLCLs
  # blcl_data$predictNKG2A <- predict(blcl_lm, blcl_data)
  # maintitle <- paste(maintitle, signif(summary(blcl_lm)$r.squared, 3))
  # # browser()
  # plot(blcl_data$predictNKG2A, blcl_data[,blcl_response], xlab = "self prediction from BLCL data", ylab = "NKG2A CD69 response", main = maintitle)
  return(blcl_lm)
}



