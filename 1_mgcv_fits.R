##--------------------------------
## Fit trend and covariate models
## CM: 5/12/22
##--------------------------------
library(mgcv)
library(ggplot2); theme_set(theme_bw())

load("all_lcdat.RData")

##-----------------
## COUNT DATA FITS 
##-----------------
lakes <- c("Feeagh", "Furnace", "Bunaveela", "BOH")

## containers
count_fits <- list()
effects_pred <- NULL
year_pred <- NULL

## set the basis functions for GAMs
bs_year <- "cr"
bs_other <- "tp"

for(lake in lakes){
    print(lake)
    ## remove October sampling, which was out of the sampling season 
    sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
    sub_dat <- droplevels(sub_dat)
    ## sum to zero contrasts - doesn't matter currently
    contrasts(sub_dat$fSite) <- contr.sum
    if(lake == "Feeagh"){
        form <- as.formula(count ~
                               s(Year, m = 1, bs = bs_year) +
                               s(doy, m = 1, bs = bs_other) +
                               s(fSite, bs = "re") +
                               s(trap_depth, m = 1, bs = bs_other) +
                               s(trap_gradient, m = 1, bs = bs_other) +
                               s(trap_number, by = survey, m = 1, bs = bs_other) +
                               s(fchain, bs = "re") +
                               offset(log(Effort)))
    }
    if(lake == "BOH"){
        form <- as.formula(count ~ 
                               s(Year, m = 1, bs = bs_year) +
                               s(doy, k = 5, m = 1, bs = bs_other) +
                               s(fSite, bs = "re", k= 3) +
                               s(trap_depth, m = 1, bs = bs_other) +
                               s(trap_gradient, k = 5, m = 1, bs = bs_other) +
                               s(trap_number, m = 1, bs = bs_other) +
                               s(fchain, bs = "re") +
                               offset(log(Effort)))
    }
    if(lake %in% c("Furnace", "Bunaveela")){
        form <- as.formula(count ~ 
                               s(Year, m = 1, bs = bs_year) +
                               s(doy, k = 5, m = 1, bs = bs_other) +
                               s(fSite, bs = "re", k= 3) +
                               s(trap_depth, k= 5, m = 1, bs = bs_other) +
                               s(trap_gradient, k= 5, m = 1, bs = bs_other) +
                               s(trap_number, k= 5, m = 1, bs = bs_other) +
                               s(fchain, bs = "re") +
                               offset(log(Effort)))
    }
    ## fit the model
    f0 <- gam(form,
              select = TRUE,
              method = "REML",
              family = nb(),
              data = sub_dat)
    ##gam.check(f0)
    ##
    count_fits[[lake]] <- f0
    ## effect predictions
    m <- 100
    pred_df0 <- data.frame(
        Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
        doy = seq(min(sub_dat$doy), max(sub_dat$doy), length = m),
        trap_depth = seq(min(sub_dat$trap_depth), max(sub_dat$trap_depth), length = m),
        trap_gradient = seq(min(sub_dat$trap_gradient), max(sub_dat$trap_gradient), length = m),
        trap_number = rep(1:20, 5),
        survey = "Russell",
        fchain = unique(sub_dat$fchain)[1],
        fSite = unique(sub_dat$fSite)[1],
        Effort = 1)
    if(lake == "Feeagh"){
        ## two survey trap number effects in Feeagh
        tmp <- pred_df0
        tmp <- subset(tmp, trap_number <= 10)
        tmp$survey <- "IFI"
        pred_df0 <- rbind(pred_df0, tmp)
    }    
    ## predictions for non random effects
    pred0 <- predict(f0, newdata = pred_df0, type = "terms", se.fit = TRUE)
    zs <- colnames(pred0$fit)
    ## for Feeagh only
    zs <- zs[zs != "s(trap_number):surveyIFI"]
    lake_effects <- NULL
    idx <- which(pred_df0$survey == "Russell")
    for(z in zs){
        v <- gsub("(s\\(|\\))", "", z)
        if(v == "trap_number:surveyRussell"){
            v <- "trap_number"
        }
        tmp <- data.frame(Lake = lake,
                          variable = v,
                          survey = "Russell",
                          x = pred_df0[idx,v],
                          yhat = pred0$fit[idx, z],
                          ylwr = pred0$fit[idx, z] - 2 * pred0$se.fit[idx, z],
                          yupr = pred0$fit[idx, z] + 2 * pred0$se.fit[idx, z]
                          )
        tmp <- unique(tmp)
        lake_effects <- rbind(lake_effects, tmp)
        rm(tmp)
    }
    ## add in IFI on Feeagh
    if(lake == "Feeagh"){
        z <- "s(trap_number):surveyIFI"
        idx <- which(pred_df0$survey == "IFI")
        v <- "trap_number"
        tmp <- data.frame(Lake = lake,
                          variable = v,
                          survey = "IFI",
                          x = pred_df0[idx,v],
                          yhat = pred0$fit[idx, z],
                          ylwr = pred0$fit[idx, z] - 2 * pred0$se.fit[idx, z],
                          yupr = pred0$fit[idx, z] + 2 * pred0$se.fit[idx, z]
                          )
        tmp <- unique(tmp)
        lake_effects <- rbind(lake_effects, tmp)
    }
    ## predictions for site and chain
    pred_df1 <- expand.grid(
        Year = sub_dat$Year[1],
        doy = mean(sub_dat$doy),
        trap_depth = mean(sub_dat$trap_depth),
        trap_gradient = mean(sub_dat$trap_depth),
        trap_number = 10,
        survey = "Russell",
        fchain = unique(sub_dat$fchain),
        fSite = unique(sub_dat$fSite),
        Effort = 1)
    ## predictions for random effects
    pred1 <- predict(f0, newdata = pred_df1, type = "terms", se.fit = TRUE)
    zs <- c("s(fSite)", "s(fchain)")
    random_effects <- NULL
    ## 
    for(z in zs){
        v <- gsub("(s\\(|\\))", "", z)
        tmp <- data.frame(Lake = lake,
                          variable = v,
                          survey = "Russell",
                          x = pred_df1[,v],
                          yhat = pred1$fit[, z],
                          ylwr = pred1$fit[, z] - 2 * pred1$se.fit[, z],
                          yupr = pred1$fit[, z] + 2 * pred1$se.fit[, z]
                          )
        tmp <- unique(tmp)
        random_effects <- rbind(random_effects, tmp)
        rm(tmp)
    }
    all_effects <- rbind(lake_effects, random_effects)
    effects_pred<- rbind(effects_pred, all_effects)
    ##------------------------
    ## GET YEARLY PREDICTIONS
    ##------------------------
    ## here we set non-year continuous covariates to their mean and exclude random effects
    pred_df2 <- data.frame(Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
                           doy = mean(sub_dat$doy, na.rm = TRUE),
                           trap_depth = mean(sub_dat$trap_depth, na.rm = TRUE),
                           trap_gradient = mean(sub_dat$trap_gradient, na.rm = TRUE),
                           trap_number = 10,
                           survey = "Russell",
                           fchain = unique(sub_dat$fchain)[1],
                           fSite = unique(sub_dat$fSite)[1],
                           Effort = 2) ## 2 codends is a fyke net
    ##
    pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                     exclude = c("s(fchain)", "s(fSite)"))
    ##
    pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
    pred_year$yhat <- exp(pred2$fit)
    pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
    pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
    year_pred <- rbind(year_pred, pred_year)
}

save(year_pred, file = "year_pred.RData")
save(effects_pred, file = "effects_pred.RData")
save(count_fits, file = "count_fits.RData")


##-----------------
## WEIGHT DATA FITS 
##-----------------
load("wdat.RData")

lakes <- c("Feeagh", "Furnace", "Bunaveela", "BOH")

weight_fits <- list()
weight_effects_pred <- NULL
weight_year_pred <- NULL

for(lake in lakes){
    print(lake)
    sub_dat <- subset(wdat, Lake == lake & !is.na(wt))
    sub_dat <- droplevels(sub_dat)
    ## sum to zero contrasts - doesn't matter currently
    contrasts(sub_dat$fSite) <- contr.sum
    if(lake != "Furnace"){
        form <- as.formula(wt ~
                               s(Year, bs = bs_year) +
                               s(doy, bs = bs_other) +
                               s(fSite, bs = "re") +
                               offset(log(Effort))
                           )
    }else{
        form <- as.formula(wt ~
                               s(Year, bs = bs_year) +
                               s(doy, bs = bs_other, k = 5) +
                               s(fSite, bs = "re") +
                               offset(log(Effort))
                           )
    }
    ##
    f0 <- gam(form,
              select = TRUE,
              ##method = "ML",
              family = tw(),
              data = sub_dat)
    ###
    weight_fits[[lake]] <- f0
    ## effect predictions 
    m <- 100
    pred_df0 <- data.frame(
        Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
        doy = seq(min(sub_dat$doy), max(sub_dat$doy), length = m),
        fSite = unique(sub_dat$fSite)[1],
        Effort = 1) ## setting this to 1 net-night
    ## predictions for non random effects
    pred0 <- predict(f0, newdata = pred_df0, type = "terms", se.fit = TRUE)
    zs <- colnames(pred0$fit)
    ## for Feeagh only
    ##zs <- zs[zs != "s(trap_number):surveyIFI"]
    lake_effects <- NULL
    for(z in zs){
        v <- gsub("(s\\(|\\))", "", z)
        if(v == "trap_number:surveyRussell"){
            v <- "trap_number"
        }
        tmp <- data.frame(Lake = lake,
                          variable = v,
                          x = pred_df0[,v],
                          yhat = pred0$fit[, z],
                          ylwr = pred0$fit[, z] - 2 * pred0$se.fit[, z],
                          yupr = pred0$fit[, z] + 2 * pred0$se.fit[, z]
                          )
        tmp <- unique(tmp)
        lake_effects <- rbind(lake_effects, tmp)
        rm(tmp)
    }
    ## predictions for site
    pred_df1 <- expand.grid(
        Year = sub_dat$Year[1],
        doy = mean(sub_dat$doy),
        Effort = 1, ## setting this to one net
        fSite = unique(sub_dat$fSite))
    ## predictions for random effects
    pred1 <- predict(f0, newdata = pred_df1, type = "terms", se.fit = TRUE)
    zs <- c("s(fSite)")
    random_effects <- NULL
    for(z in zs){
        v <- gsub("(s\\(|\\))", "", z)
        tmp <- data.frame(Lake = lake,
                          variable = v,
                          x = pred_df1[,v],
                          yhat = pred1$fit[, z],
                          ylwr = pred1$fit[, z] - 2 * pred1$se.fit[, z],
                          yupr = pred1$fit[, z] + 2 * pred1$se.fit[, z]
                          )
        tmp <- unique(tmp)
        random_effects <- rbind(random_effects, tmp)
        rm(tmp)
    }
    all_effects <- rbind(lake_effects, random_effects)
    weight_effects_pred<- rbind(weight_effects_pred, all_effects)
    ##------------------------
    ## GET YEARLY PREDICTIONS
    ##------------------------
    pred_df2 <- data.frame(Year = seq(min(sub_dat$Year), max(sub_dat$Year), length = m),
                          doy = 200,
                          fSite = unique(sub_dat$fSite)[1],
                          Effort = 1) ## setting this to one net
    ##
    pred2 <- predict(f0, newdata = pred_df2, se.fit = TRUE,
                     ##exclude = c("s(fSite)", "s(doy)"))
                     exclude = c("s(fSite)"))
    ##
    pred_year <- data.frame(Lake = lake, Year = pred_df2$Year)
    pred_year$yhat <- exp(pred2$fit)
    pred_year$lwr <- exp(pred2$fit - 2 * pred2$se.fit)
    pred_year$upr <- exp(pred2$fit + 2 * pred2$se.fit)
    weight_year_pred <- rbind(weight_year_pred, pred_year)
}

save(weight_year_pred, file = "weight_year_pred.RData")
save(weight_effects_pred, file = "weight_effects_pred.RData")
save(weight_fits, file = "weight_fits.RData")

