##-------------------------
## Some additional plots
## trends and otter guards
## CM: 11/10/2023
##-------------------------
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

load("../../data/all_lcdat.RData")

load("../../data/year_pred.RData")

lakes <- c("BOH", "Furnace", "Feeagh", "Bunaveela")

year_pred$Lake <- factor(year_pred$Lake, levels = lakes)
var_lab <- "Year effect"

lcdat$Lake <- factor(lcdat$Lake, levels = lakes)

mean_count <- aggregate(count ~ Year + Lake, mean, data = lcdat)

## get the count per pair
pair_tab <- table(lcdat$fchain_fpair)
pair_remove <- names(pair_tab)[pair_tab < 2]
## 2 unpaired codends out of 5302 pairs

pair_count <- aggregate(count ~ Year + Lake + fchain_fpair, sum, data = lcdat)

pogc <-
    ggplot(pair_count, aes(x = Year, y = count, colour = Year >= 2015)) +
    ##geom_boxplot(aes(y = count, group = Year)) +
    stat_summary(fun.data = "mean_cl_boot", size = 0.3) +
    facet_wrap(~Lake, scales = "free", ncol = 2) +
    geom_vline(xintercept = 2014.5, lty = 2) +
    scale_colour_manual("Otter guard", values = c("slategrey", "blue")) +
    xlab("Year") +
    ylab("Count per fyke net (two codends)") +
    theme(strip.background = element_rect(fill = "grey"))


## weight
load("../../data/wdat.RData")

wdat <- subset(wdat, Month != "Oct" & !is.na(wt))
load("../../data/weight_year_pred.RData")

weight_year_pred$Lake <- factor(weight_year_pred$Lake, levels = lakes)

wdat$Lake <- factor(wdat$Lake, levels = lakes)

pogw <-
    ggplot(wdat, aes(x = Year, y = wt/Effort, colour = Year >= 2015)) +
    ##geom_boxplot(aes(group = Year)) +
    stat_summary(fun.data = "mean_cl_boot", size = 0.3) +
    facet_wrap(~Lake, scales = "free", ncol = 2) +
    geom_vline(xintercept = 2014.5, lty = 2) +
    scale_colour_manual("Otter guard", values = c("slategrey", "blue")) +
    xlab("Year") +
    ylab("Mass per fyke net (kg)") +
    theme(strip.background = element_rect(fill = "grey"))

pdf("../../tex/figures/otter_guard_visual.pdf", height = 7, width = 9)
print(pogc)
print(pogw)
dev.off()
