##--------------------------------
## Plots for the paper
## CM: 16/2/23
##--------------------------------
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

load("../../data/all_lcdat.RData")

## basic data plots
p0 <- ggplot(lcdat, aes(x = jitter(Year), y = count)) +
    geom_point() +
    facet_wrap(~ Lake, ncol = 1) +
    ylab("Count per trap") +
    theme(strip.background = element_rect(fill = "grey"))

pdf("raw_data.pdf", height = 8, width = 10)
print(p0)
dev.off()


site_df <- as.data.frame(with(subset(lcdat, survey != "IFI"), table(Year, Site, Lake)))

site_df <- subset(site_df, Freq > 0)

site_df$Year <- as.numeric(as.character(site_df$Year))

library(viridis)
p1 <- ggplot(site_df, aes(x = Year, y = Site, fill = Freq)) +
    geom_tile(color="white", size=0.1) +
    ##facet_wrap(~ Lake, scales = "free")
    facet_wrap(~ Lake, scales="free_y", ncol = 1) +
    scale_fill_viridis(name="# Codends") +
    theme(strip.background = element_rect(fill = "grey"),
          legend.position = "none")

p2 <- ggplot(lcdat, aes(x = Year, y = trap_depth)) +
    geom_point(alpha = 0.1) +
    facet_wrap(~ Lake, ncol = 1) +
    ylab("Trap depth") +
    theme(strip.background = element_rect(fill = "grey"))

p3 <- ggplot(lcdat, aes(x = Year, y = trap_gradient)) +
    geom_point(alpha = 0.1) +
    facet_wrap(~ Lake, ncol = 1) +
    ylab("Trap gradient") +
    theme(strip.background = element_rect(fill = "grey"))

today <- format(Sys.time(), "_%d_%m_%Y")

png(paste0("../tex/figures/", "Variable_plots", today, ".png"), height = 6, width = 10, units = "in", res = 400)
grid.arrange(p0, p1, p2, p3, nrow = 1)
dev.off()

##---------------
## RESULTS PLOTS
##---------------
## covariate effects

load("../../data/effects_pred.RData")
load("../../data/count_fits.RData")

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

effects_pred$Lake <- factor(effects_pred$Lake, levels = lakes)

pvalue_df <- NULL

for(lake in lakes){
    tmp <- summary(count_fits[[lake]])
    tmp <- tmp$s.table
    colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
    df <- as.data.frame(tmp)
    df$var <- rownames(df)
    rownames(df) <- NULL
    df$Lake <- lake
    pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)


vars <- unique(effects_pred$variable)

vars <- vars[!vars %in% c("fSite", "fchain")]

for(v in vars){
    sub <- subset(effects_pred, variable == v)
    sub$x <- as.numeric(sub$x)
    ##sub$Lake <- factor(sub$Lake, levels = c("BOH", "Furnace", "Feeagh", "Bunaveela"))
    var_lab <- stringr::str_to_title(gsub("\\_", " ", v))
    ##
    pv <- pvalue_df[grep(v, pvalue_df$var),]
    pv$survey <- "Russell"
    pv$survey[grep("surveyIFI", pv$var)] <- "IFI"
    ##if(length() > 0){
    ##    pv <- pv[-grep("surveyIFI", pv$var), ] ## double survey trap number effect for IFI 
    ##}
    if(v == "trap_number"){
        tmpF <- subset(pv, Lake == "Feeagh")
        p_string <- paste(tmpF$p_value2, collapse = ", ") ## IFI first
        tmp <- subset(pv, Lake != "Feeagh")
        tmp0 <- subset(tmpF, var == "s(trap_number):surveyRussell")
        tmp0$p_value2 <- p_string
        pv <- rbind(tmp, tmp0)
    }
    p <-
    ggplot(sub, aes(x = x, y = yhat, group = survey)) +
        geom_ribbon(aes(ymin = ylwr, ymax = yupr), fill = "grey", alpha = 0.7) +
        geom_line(colour = "black", lwd = 0.5) +
        facet_wrap(~Lake, nrow = 1) +
        ##xlab(var_lab) +
        xlab("") +
        ylab(paste(var_lab, "effect")) +
        geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
                  hjust = 1, vjust  = 1.5, size = 3) +
        theme(plot.margin = unit(c(0, 1, 0, 1), "lines"), legend.position = "none") +
        scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))
    if(v != "Year"){
        p <- p + theme(##strip.background = element_rect("grey"),
                     strip.text.x = element_blank() , 
                     strip.background = element_blank())
    }
    assign(paste0("p_", v), p)
}

## chain effects
sub <- subset(effects_pred, variable == "fchain")
var_lab <- "Chain effect"
pv <- pvalue_df[grep("chain", pvalue_df$var),]

p_chain <-
    ggplot(sub, aes(x = exp(yhat))) +
    geom_histogram(fill = "slategrey", colour = "grey") +
    facet_wrap(~Lake, nrow = 1) +
    xlab("") +
    ylab("Chain effect") +
    theme(plot.margin = unit(c(0, 1, 0, 1), "lines"),
          strip.text.x = element_blank() , 
          strip.background = element_blank(),
          legend.position = "none")  +
    geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
              hjust = 1.2, vjust = 2, size = 3) +
    scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))


## site effects
sub <- subset(effects_pred, variable == "fSite")
sub$x[sub$x == "Back Weir"] <- "BW"
sub$x[sub$x == "Front Weir"] <- "FW"
sub$x[sub$x == "North"] <- "N"
sub$x[sub$x == "South"] <- "S"
sub$x[sub$x == "Weir"] <- "W"

var_lab <- "Site effect"

pv <- pvalue_df[grep("Site", pvalue_df$var),]

p_site <-
    ggplot(sub, aes(y = yhat, x = x)) +
    geom_point() +
    geom_segment(aes(y = ylwr, yend = yupr, xend = x)) +
    geom_hline(yintercept = 0, lty = 2) +
    facet_wrap(~Lake, nrow = 1, scales = "free_x") +
    xlab("") +
    ylab("Site effect") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.margin = unit(c(0, 1, 0, 1), "lines"),
          strip.text.x = element_blank() , 
          strip.background = element_blank(),
          legend.position = "none"          
          )  +
    geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
              hjust = 1.3, vjust = 2, size = 3) +
    scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))

png(paste0("../tex/figures/", "Count_effect_plots", today, ".png"), height = 10, width = 7, units = "in", res = 400)
grid.arrange(p_Year, p_doy, p_trap_depth, p_trap_gradient, p_trap_number, p_chain, p_site, ncol = 1)
dev.off()

##-----------------------
## AVERAGED YEARLY TREND
##-----------------------
load("../../data/year_pred.RData")

year_pred$Lake <- factor(year_pred$Lake, levels = lakes)
var_lab <- "Year effect"

lcdat$Lake <- factor(lcdat$Lake, levels = lakes)

mean_count <- aggregate(count ~ Year + Lake, mean, data = lcdat)

## get the count per pair
pair_tab <- table(lcdat$fchain_fpair)
pair_remove <- names(pair_tab)[pair_tab < 2]
## 2 unpaired codends out of 5302 pairs

pair_count <- aggregate(count ~ Year + Lake + fchain_fpair, sum, data = lcdat)

py0 <-
    ggplot(pair_count, aes(x = Year, y = count)) +
    geom_boxplot(aes(y = count, group = Year)) +
    stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 0.3) +
    facet_wrap(~Lake, scales = "free", ncol = 4) +
    geom_ribbon(data = year_pred, aes(ymin = lwr, ymax = upr, y = yhat), fill = "grey", alpha = 0.4) +
    geom_line(data = year_pred, aes(y = yhat), colour = "black", lwd = 0.5) +
    ##scale_y_continuous(trans = "sqrt") +
    xlab("Year") +
    ylab("Count per fyke net (two codends)") +
    theme(strip.background = element_rect(fill = "grey"))

py1 <-
    ggplot(year_pred, aes(x = Year, y = yhat)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey", alpha = 0.4) +
    geom_line(colour = "black", lwd = 0.5) +
    facet_wrap(~Lake, scales = "free", ncol = 4) +
    ##geom_point(data = mean_count, aes(y = count)) +
    xlab("Year") +
    ylab("Standardised count per fyke net") +
    theme(strip.background = element_rect(fill = "grey"))


## percentage change between the start and the end of the time series
library(MASS) ## for random multivariate normal simulation

percent_df <- NULL

for(lake in lakes){
    f0 <- count_fits[[lake]]
    sub_dat <- subset(lcdat, Lake == lake & Month != "Oct" & !is.na(count))
    sub_dat <- droplevels(sub_dat)
    pred_df2 <- data.frame(Year = c(min(sub_dat$Year), max(sub_dat$Year)),
                           doy = mean(sub_dat$doy, na.rm = TRUE),
                           trap_depth = mean(sub_dat$trap_depth, na.rm = TRUE),
                           trap_gradient = mean(sub_dat$trap_gradient, na.rm = TRUE),
                           trap_number = 10,
                           survey = "Russell",
                           fchain = unique(sub_dat$fchain)[1],
                           ##fpair = unique(sub_dat$fpair)[1],
                           fSite = unique(sub_dat$fSite)[1],
                           Effort = 2) ## per net = 2 codends
    ##
    if(lake == "Feeagh"){
        Xp <- predict(f0, newdata = pred_df2, se.fit = TRUE, type = "lpmatrix", exclude = c("s(fchain)", "s(fSite)"))
    }else{
        Xp <- predict(f0, newdata = pred_df2, se.fit = TRUE, type = "lpmatrix", exclude = c("s(fchain)", "s(fSite)"))
    }
    ## draw from posterior - see: https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
    br <- mvrnorm(10000,coef(f0), vcov(f0))
    lnN <- t(Xp %*% t(br))
    ## percentage decline
    percent_decline <- 100 * (1 - exp(lnN[, 2] - lnN[, 1]))
    df <- data.frame(Lake = lake, pd = percent_decline)
    percent_df <- rbind(percent_df, df)
    rm(df)
}

percent_df$Lake <- factor(percent_df$Lake, levels = c("BOH", "Furnace", "Feeagh", "Bunaveela"))

agg0 <- aggregate(pd ~ Lake, FUN = mean, data = percent_df)
agg1 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.025, data = percent_df)
names(agg1)[2] <- c("lwr")
agg2 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.975, data = percent_df)
names(agg2)[2] <- c("upr")

summary_df <- merge(merge(agg0, agg1), agg2)

write.csv(summary_df, file = "../data/count_percent_decline.csv", row.names = FALSE)

py2 <- ggplot(percent_df, aes(x = pd)) +
    geom_density(fill = "grey") +
    facet_wrap(~ Lake, ncol = 4, scales = "free_y") +
    ylab("Density") +
    geom_vline(data = summary_df, aes(xintercept = pd)) +
    geom_vline(data = summary_df, aes(xintercept = lwr), lty = 3) +
    geom_vline(data = summary_df, aes(xintercept = upr), lty = 3) +
    geom_vline(xintercept = 0, lty = 4) +
    xlab("Percentage decline in standardised count trend 1987-2022")

png(paste0("../tex/figures/", "Count_percent_decline_plots", today, ".png"), height = 8, width = 10, units = "in", res = 400)
grid.arrange(py0, py1, py2, nrow = 3)
dev.off()

##---------------
## WEIGHTS PLOTS
##---------------
load("../../data/weight_effects_pred.RData")
load("../../data/weight_fits.RData")
load("../../data/wdat.RData")

wdat <- subset(wdat, Month != "Oct" & !is.na(wt))

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

weight_effects_pred$Lake <- factor(weight_effects_pred$Lake, levels = lakes)

pvalue_df <- NULL

for(lake in lakes){
    tmp <- summary(weight_fits[[lake]])
    tmp <- tmp$s.table
    colnames(tmp)[colnames(tmp) == "p-value"] <- "p_value"
    df <- as.data.frame(tmp)
    df$var <- rownames(df)
    rownames(df) <- NULL
    df$Lake <- lake
    pvalue_df <- rbind(pvalue_df, df)
}

pvalue_df$p_value2 <- round(pvalue_df$p_value, 3)
pvalue_df$p_value2[pvalue_df$p_value < 0.001] <- "<0.001"

pvalue_df$Lake <- factor(pvalue_df$Lake, levels = lakes)

vars <- unique(weight_effects_pred$variable)

vars <- vars[!vars %in% c("fSite", "fchain")]

for(v in vars){
    sub <- subset(weight_effects_pred, variable == v)
    sub$x <- as.numeric(sub$x)
    var_lab <- stringr::str_to_title(gsub("\\_", " ", v))
    ##
    pv <- pvalue_df[grep(v, pvalue_df$var),]
    if(length(grep("surveyIFI", pv$var)) > 0){
        pv <- pv[-grep("surveyIFI", pv$var), ] ## double survey trap number effect for IFI 
    }
    p <- ggplot(sub, aes(x = x, y = yhat)) +
        geom_ribbon(aes(ymin = ylwr, ymax = yupr), fill = "grey", alpha = 0.7) +
        geom_line(colour = "black", lwd = 0.5) +
        facet_wrap(~Lake, nrow = 1) +
        ##xlab(var_lab) +
        xlab("") +
        ylab(paste(var_lab, "effect")) +
        geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
                  hjust = 1, vjust  = 1.5, size = 3) +
        theme(plot.margin = unit(c(0, 1, 0, 1), "lines"), legend.position = "none") +
        scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))
    if(v != "Year"){
        p <- p + theme(##strip.background = element_rect("grey"),
                     strip.text.x = element_blank() , 
                     strip.background = element_blank())
    }
    assign(paste0("p_", v), p)
}


## site effects
sub <- subset(weight_effects_pred, variable == "fSite")
sub$x[sub$x == "Back Weir"] <- "BW"
sub$x[sub$x == "Front Weir"] <- "FW"
sub$x[sub$x == "North"] <- "N"
sub$x[sub$x == "South"] <- "S"
sub$x[sub$x == "Weir"] <- "W"

var_lab <- "Site effect"

pv <- pvalue_df[grep("Site", pvalue_df$var),]

p_site <-
    ggplot(sub, aes(y = yhat, x = x)) +
    geom_point() +
    geom_segment(aes(y = ylwr, yend = yupr, xend = x)) +
    geom_hline(yintercept = 0, lty = 2) +
    facet_wrap(~Lake, nrow = 1, scales = "free_x") +
    xlab("") +
    ylab("Site effect") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          plot.margin = unit(c(0, 1, 0, 1), "lines"),
          strip.text.x = element_blank() , 
          strip.background = element_blank(),
          legend.position = "none"          
          )  +
    geom_text(data = pv, aes(x = Inf, y = Inf, label = p_value2, colour = p_value < 0.05),
              hjust = 1.3, vjust = 2, size = 3) +
    scale_colour_manual(values = c("black", "darkgrey"), breaks = c("TRUE", "FALSE"))

png(paste0("../tex/figures/", "Weight_effect_plots", today, ".png"), height = 7, width = 7, units = "in", res = 400)
grid.arrange(p_Year, p_doy, p_site, ncol = 1)
dev.off()

##----------------------
## AVERAGE YEARLY TREND
##----------------------
load("../../data/weight_year_pred.RData")
load("../../data/wdat.RData")

weight_year_pred$Lake <- factor(weight_year_pred$Lake, levels = lakes)

wdat$Lake <- factor(wdat$Lake, levels = lakes)

py0 <-
    ggplot(wdat, aes(x = Year, y = wt/Effort)) +
    geom_boxplot(aes(group = Year)) +
    stat_summary(fun.data = "mean_cl_boot", colour = "red", size = 0.3) +
    facet_wrap(~Lake, scales = "free", ncol = 4) +
    geom_ribbon(data = weight_year_pred, aes(ymin = lwr, ymax = upr, y = yhat), fill = "grey", alpha = 0.4) +
    geom_line(data = weight_year_pred, aes(y = yhat), colour = "black", lwd = 0.5) +
    ##scale_y_continuous(trans = "sqrt") +
    xlab("Year") +
    ylab("Mass per fyke net (kg)") +
    theme(strip.background = element_rect(fill = "grey"))

py1 <-
    ggplot(weight_year_pred, aes(x = Year, y = yhat)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "grey", alpha = 0.4) +
    geom_line(colour = "black", lwd = 0.5) +
    facet_wrap(~Lake, scales = "free", ncol = 4) +
    ##geom_point(data = mean_count, aes(y = count)) +
    xlab("Year") +
    ylab("Standardised mass per fyke net (kg)") +
    theme(strip.background = element_rect(fill = "grey"))

## percentage change between the start and the end of the time series
weight_percent_df <- NULL

for(lake in lakes){
    f0 <- weight_fits[[lake]]
    sub_dat <- subset(wdat, Lake == lake & Month != "Oct" & !is.na(wt))
    sub_dat <- droplevels(sub_dat)
    pred_df2 <- data.frame(Year = c(min(sub_dat$Year), max(sub_dat$Year)),
                           doy = 200,
                           fSite = unique(sub_dat$fSite)[1],
                           Effort = 1)
    ##
    Xp <- predict(f0, newdata = pred_df2, se.fit = TRUE, type = "lpmatrix", exclude = c("s(fSite)"))
    ## draw from posterior - see: https://www.maths.ed.ac.uk/~swood34/mgcv/check-select.pdf
    br <- mvrnorm(10000,coef(f0), vcov(f0))
    lnN <- t(Xp %*% t(br))
    ## percentage decline
    percent_decline <- 100 * (1 - exp(lnN[, 2] - lnN[, 1]))
    df <- data.frame(Lake = lake, pd = percent_decline)
    weight_percent_df <- rbind(weight_percent_df, df)
    rm(df)
}

weight_percent_df$Lake <- factor(weight_percent_df$Lake, levels = lakes)

agg0 <- aggregate(pd ~ Lake, FUN = mean, data = weight_percent_df)
agg1 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.025, data = weight_percent_df)
names(agg1)[2] <- c("lwr")
agg2 <- aggregate(pd ~ Lake, FUN = quantile, probs = 0.975, data = weight_percent_df)
names(agg2)[2] <- c("upr")

weight_summary_df <- merge(merge(agg0, agg1), agg2)

py2 <- ggplot(weight_percent_df, aes(x = pd)) +
    geom_density(fill = "grey") +
    facet_wrap(~ Lake, nrow = 1, scales = "free_y") +
    ylab("Density") +
    geom_vline(data = weight_summary_df, aes(xintercept = pd)) +
    geom_vline(data = weight_summary_df, aes(xintercept = lwr), lty = 3) +
    geom_vline(data = weight_summary_df, aes(xintercept = upr), lty = 3) +
        xlab("Percentage decline in standardised mass trend 1987-2022")


png(paste0("../tex/figures/", "Weight_percent_decline_plots", today, ".png"), height = 8, width = 10, units = "in", res = 400)
grid.arrange(py0, py1, py2, ncol = 1)
dev.off()

write.csv(weight_summary_df, file = "../data/weight_percent_decline.csv", row.names = FALSE)
