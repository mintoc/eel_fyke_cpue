##-----------------------------
## Residual diagnostics
## CM: 30/01/23
##
##-----------------------------

library(mgcv)
library(DHARMa) 
load("count_fits.RData")
load("all_lcdat.RData")

##------------------------
## COUNT MODEL RESIDUALS - uses negative binomial distribution
##------------------------
lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

today <- format(Sys.time(), "_%d_%m_%Y")

cex_lab <- 0.9

png(paste0("../tex/figures/Fig_S0_residual_diagnostics", today, ".png"), height = 8, width = 7, units = "in", res = 400)
par(mfrow = c(4, 3), mar = c(2, 3, 1, 1), oma = c(2, 1, 1, 1))
set.seed(101)
for(lake in lakes){
    isb <- lake == "Bunaveela"
    f0 <- count_fits[[lake]]
    ## Calculate DHARMa randomised quantile residuals
    simulationOutput <- simulateResiduals(fittedModel = f0, plot = F)
    resids <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
    qqnorm(resids, bty = "l", main = "", col = "slategrey"); qqline(resids)
    legend("topleft", legend = lake, bty = "n", cex = 1.1)
    if(isb){mtext(side = 1, text = "Theoretical quantile", line = 2.5, cex = cex_lab)}
    ##
    xlim <- c(-1, 1) * max(abs(resids))
    hist(resids, probability = TRUE, bty = "l", xlim = xlim, breaks = 30, border = "lightgrey", main = "")
    curve(dnorm(x), col = "red", add = TRUE, n = 1e3)
    abline(v = 0, lty = 2)
    if(isb){mtext(side = 1, text = "Quantile residual", line = 2.5, cex = cex_lab)}
    ##
    plot(predict(f0), resids, bty = "l", col = "slategrey")
    abline(h = 0, lty = 2)
    if(isb){mtext(side = 1, text = "ln(Fitted value)", line = 2.5, cex = cex_lab)}
    if(isb){mtext(side = 2, text = "Sample quantile", line = -0.5, outer = TRUE, cex = cex_lab)}
    if(isb){mtext(side = 2, text = "Probability", line = -18, outer = TRUE, cex = cex_lab)}
    if(isb){mtext(side = 2, text = "Quantile residual", line = -35, outer = TRUE, cex = cex_lab)}
}
dev.off()


## residual autocorrelation across the fyke chain
acf_df <- NULL

##png(paste0("../tex/figures/Fig_S2_residual_autocorrelation_", today, ".png"), height = 8, width = 7, units = "in", res = 400)
## autocorrelation
set.seed(101)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
for(lake in lakes){
    isb <- lake == "Bunaveela"
    f0 <- count_fits[[lake]]
    simulationOutput <- simulateResiduals(fittedModel = f0, plot = F)
    resids <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
    sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
    sub_dat <- droplevels(sub_dat)
    ##
    chains <- levels(sub_dat$fchain)
    sub_dat$resid <- resids
    ar1_vec <- sapply(chains, function(z){
        tmp <- subset(sub_dat, chain == z)
        tmp <- tmp[order(tmp$trap_number),]
        acf(tmp$resid, plot = FALSE, na.action = na.omit)$acf[2,1,1]
    })
    tmp <- data.frame(Lake = lake, fchain = chains, acf = ar1_vec)
    acf_df <- rbind(acf_df, tmp)
    ## what would a random distibution look like?
    ##rand_vec <- sapply(1:length(chains), function(z){
    rand_vec <- sapply(1:1e5 , function(z){
         y <- rnorm(20)
         ## first order unbiased - true value of mean used
         ###1/19 * sum(y[1:19] * y[2:20])
         acf(y, plot = FALSE)$acf[2,1,1]
    })
    hist(ar1_vec, breaks = seq(-1, 1, by = 0.05), xlim = c(-1, 1), border = "lightgrey", col = "grey", probability = TRUE, ylim = c(0, 3.5), main = "", xlab = "", ylab = "")
    n <- 19 ## pairs
    ##curve((1 - x^2)^((n-4)/2) / beta(a = 1/2, b = (n-2)/2), col = "blue", add = TRUE)
    lines(density(rand_vec), col = "black")
    legend("topleft", legend = lake, bty = "n")
    abline(v = 0, lty = 2)
    ##lines(density(rand_vec), col = "red")
}
mtext(side = 1, line = 0.5, text = "Trap autocorrelation", outer = TRUE)
mtext(side = 2, line = 0.5, text = "Density", outer = TRUE)
##legend("topright", legend = c("Sample acf", "Random acf", "Unbiased acf"), lty = c(NA, 1, 1), col = c("grey", "black", "blue"), pch = c(15, NA, NA), bty = "n")
legend("topright", legend = c("Observed acf", "Random acf"), lty = c(NA, 1), col = c("grey", "black"), pch = c(15, NA), bty = "n")
##dev.off()

acf_df2 <- merge(unique(lcdat[, c("Year", "Date", "fSite", "fchain")]), acf_df)

## critical bands
cb <- acf_df2[, c("Lake", "fSite")]
cb$n <- 20
cb$n[cb$fSite == "IFI"] <- 10
cb$crit <- with(cb, 1.96 / (sqrt(n -1)))

png(paste0("../tex/figures/Figure_SX_acf_time_", today, ".png"), height = 9, width = 8, units = "in", res = 400)
ggplot(acf_df2, aes(x = Date, y = acf)) +
    geom_point() +
    facet_wrap(~ Lake + fSite) +
    geom_hline(yintercept = 0, lty = 1) +
    geom_hline(data = cb, aes(yintercept = c(1, -1) * crit), lty = 2) +
    xlab("Date") +
    ylab("First-order autocorrelation of residuals along a chain")
dev.off()

## residuals vs month

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

p_tab <- NULL

lcdat$fMonth <- factor(lcdat$Month, levels = c("Apr", "May", "Jun", "Jul", "Aug", "Sep"))

png(paste0("../tex/figures/Fig_S3_residual_vs_month_", today, ".png"), height = 8, width = 7, units = "in", res = 400)
par(mfrow = c(2, 2), mar = c(2, 2, 1, 1), oma = c(2, 2, 1, 1))
for(lake in lakes){
    simulationOutput <- simulateResiduals(fittedModel = count_fits[[lake]], plot = F)
    res <- residuals(simulationOutput, quantileFunction = qnorm)
    sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
    sub_dat <- droplevels(sub_dat)    
    boxplot(res ~ sub_dat$fMonth, notch = TRUE, bty = "l")
    legend("topleft", legend = lake, bty = "n")
    abline(h = 0, lty = 1)
}
mtext(side = 1, text = "Month", line = 0, outer = TRUE)
mtext(side = 2, text = "Residual", line = 0, outer = TRUE)
dev.off()

##------------------------
## WEIGHT MODEL RESIDUALS - uses Tweedie distribution
##------------------------
load("weight_fits.RData")
load("wdat.RData")

png(paste0("../tex/figures/Fig_S0_residual_diagnostics_weights", today, ".png"), height = 8, width = 7, units = "in", res = 400)
set.seed(101)
par(mfrow = c(4, 3), mar = c(2, 3, 1, 1), oma = c(2, 1, 1, 1))
for(lake in lakes){
    isb <- lake == "Bunaveela"
    f0 <- weight_fits[[lake]]
    ##
    simulationOutput <- simulateResiduals(fittedModel = f0, plot = F)
    resids <- residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-5,5))
    qqnorm(resids, bty = "l", main = "", col = "slategrey"); qqline(resids)
    legend("topleft", legend = lake, bty = "n", cex = 1.1)
    if(isb){mtext(side = 1, text = "Theoretical quantile", line = 2.5, cex = cex_lab)}
    ##
    xlim <- c(-1, 1) * max(abs(resids))
    hist(resids, probability = TRUE, bty = "l", xlim = xlim, breaks = 30, border = "lightgrey", main = "")
    curve(dnorm(x), col = "red", add = TRUE, n = 1e3)
    abline(v = 0, lty = 2)
    if(isb){mtext(side = 1, text = "Quantile residual", line = 2.5, cex = cex_lab)}
    ##
    plot(predict(f0), resids, bty = "l", col = "slategrey")
    abline(h = 0, lty = 2)
    if(isb){mtext(side = 1, text = "ln(Fitted value)", line = 2.5, cex = cex_lab)}
    if(isb){mtext(side = 2, text = "Sample quantile", line = -0.5, outer = TRUE, cex = cex_lab)}
    if(isb){mtext(side = 2, text = "Probability", line = -18, outer = TRUE, cex = cex_lab)}
    if(isb){mtext(side = 2, text = "Quantile residual", line = -35, outer = TRUE, cex = cex_lab)}
}
dev.off()
