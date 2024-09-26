# ~Library ####
# Load packages
library(AICcmodavg) # for aictab()
library(beeswarm)   # for beeswarm()
library(car)        # for leveneTest()
library(DHARMa)     # for model validation
library(FSA)        # for dunnTest ()
library(glmmTMB)    # for glmmTMB()
library(lattice)    # for dotplot(); xyplot()
library(lme4)       # for lmer()
library(MuMIn)      # for model.sel()
library(nlme)       # for gls(); lme()
library(pscl)
library(rcompanion) # for cldList()
library(usdm)       # for vifstep()

# 1) Dataset variables ####
## a. Response Variable:
# "*Chitons*": abundance of fluorescent chitons, i.e. *Ischnoplax pectinata* (G. B. Sowerby II, 1840).

## b. Explanatory Variable:
# "*Site*": study sites.
# "*Samp.Time*": sampling time (min).
# "*Temp*": water temperature (°C).
# "*Wt.level*": water level in relation to the base of the boulder (cm).
# "*Weight*": boulder weight (Kg).
# "*Boulder.area*": boulder surface area (cm^2^).
# "*Exp.side.area*": exposed surface side area (cm^2^).
# "*Flu.cover*": coverage of fluorescent substrates is the sum of '*Live.cca*' and '*Live.pey*' coverage (%).
# "*Live.cca*": live crustose coralline algae, i.e. *Lithothamnium* sp. (%).
# "*Live.pey*": live red soft crust algae, i.e. *Peyssonnelia* sp. (%).
# "*Unflu.cover*": coverage of unfluorescent substrates is the sum of '*Dead.red.crust*' and '*Unflu.others*' coverage (%).
# "*Dead.red.crust*": dead red crusts algae, i.e. *Lithothamnium* sp. and *Peyssonnelia* sp. together. (%).
# "*Unflu.other*": coverage of "bare rock" and microbial mats (%).
# "*Asc.cover*": coverage of ascidians (%).
# "*Bryo.cover*": coverage of bryozoans (%).
# "*Spong.cover*": coverage of sponges (%).
# "*Others*": coverage of holes and unidentifiable substrates (%).

# 2) Data exploration ####

# Dataset: "*chitons.data.csv*"
side <- read.csv("chitons.data.csv", sep=";",dec=".", header=T)
side$fSite <- factor(side$Site)
str(side)

## a. Missing data ?
summary(side) # NO

## b. Balanced sampling ?
tapply(side$Chitons,side$fSite,length) # YES

## c. Outliers Y & X ?

# Boxplot
boxplot(decostand(side[,2:17],
                  method="standardize"), ylab="Standardize value")

# Cleveland dotplot
dotplot(as.matrix(side[,2:17]), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data")

## d. Normality Y ?

# Histogram:
hist(side$Chitons, prob = T, breaks=15, ylim=c(0,1),
     xlab="Abundance of chitons", main="Histogram")
lines(density(side$Chitons), col="blue")

# QQ-plot:
qqnorm(side$Chitons)
qqline(side$Chitons, col="blue", lwd=2)

# Shapiro-Wilk's test:
shapiro.test(side$Chitons) # H0 rejected

## e. Homogeneity Y?

boxplot(Chitons~fSite, data=side) # I don't think so

variance <- tapply(side$Chitons,side$fSite,var)
variance

variance [2]/ variance [1]  # Buzios is 66x bigger than BF

# Teste de Bartlett:
bartlett.test(side$Chitons ~ side$fSite) # H0 rejected (Heteroscedastic)

# Teste de Fligner-Killeen:
fligner.test(side$Chitons ~ side$fSite) # H0 rejected (Heteroscedastic)

leveneTest(Chitons ~ fSite, data=side, center=mean) # H0 rejected (Heteroscedastic)

## f. Zero Trouble Y ?

# Frequency plot
table(side$Chitons)
barplot(table(side$Chitons), ylim=c(0,100),
        ylab="Frequency", xlab="Observed values")

# How many zeros?
sum(side$Chitons==0)/length(side$Chitons)*100 # 63.33%

# Possible inflation, but not guaranteed (See: Warton, 2005)
# Warton, D. I. (2005). Many zeros does not mean zero inflation: comparing the goodness-of-fit of parametric models to multivariate abundance data. Environmetrics 16(3), 275-289

## g. Collinearity X ?

# *pairs* fuctions:
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method = "kendall")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}

# Scatterplots correlations
pairs(side[,3:17],panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Collinear Variables:

# Weight ~ Boulder.area (0.68)
# Weight ~ Exp.side.area (0.56)
# Boulder.area ~ Exp.side.area (0.47)
# Flu.cover ~ Live.cca (0.72)
# Flu.cover ~ Unflu.cover (-0.63)
# Live.cca ~ Unflu.cover (-0.56)
# Unflu.cover ~ Dead.red.crust (0.39)


# Variance Inflation Factor
X<-side[,3:17]
vifstep(X, th=3,keep = "Flu.cover")

## h. Relationships Y & X

# Scatterplots
pairs(side[,c(2:10)],
      panel=panel.smooth,
      diag.panel=panel.hist,
      lower.panel=panel.cor
)

pairs(side[,c(2,11:17)],
      panel=panel.smooth,
      diag.panel=panel.hist,
      lower.panel=panel.cor
)


# Chitons ~ Flu.cover | Site
xyplot(Chitons ~ Flu.cover | Site, data=side, type=c("p","r"))

# Chitons ~ Wt.level | Site
xyplot(Chitons ~ Wt.level | Site, data=side, type=c("p","r"))

# Chitons ~ Boulder.area | Site
xyplot(Chitons ~ Boulder.area | Site, data=side, type=c("p","r"))

# Chitons ~ Exp.side.area | Site
xyplot(Chitons ~ Exp.side.area | Site, data=side, type=c("p","r"))

# Chitons ~ Weight | Site
xyplot(Chitons ~ Weight | Site, data=side, type=c("p","r"))

# Chitons ~ Samp.time | Site
xyplot(Chitons ~ Samp.time | Site, data=side, type=c("p","r"))

# Chitons ~ Live.cca | Site
xyplot(Chitons ~ Live.cca | Site, data=side, type=c("p","r"))

# Chitons ~ Live.pey | Site
xyplot(Chitons ~ Live.pey | Site, data=side, type=c("p","r"))

# Chitons ~ Unflu.cover | Site
xyplot(Chitons ~ Unflu.cover | Site, data=side, type=c("p","r"))

# Chitons ~ Dead.red.crust | Site
xyplot(Chitons ~ Dead.red.crust | Site, data=side, type=c("p","r"))

# Chitons ~ Asc.cover | Site
xyplot(Chitons ~ Asc.cover | Site, data=side, type=c("p","r"))

# Chitons ~ Bryo.cover | Site
xyplot(Chitons ~ Bryo.cover | Site, data=side, type=c("p","r"))

# Chitons ~ Spong.cover | Site
xyplot(Chitons ~ Spong.cover | Site, data=side, type=c("p","r"))

# Chitons ~ Others | Site
xyplot(Chitons ~ Others | Site, data=side, type=c("p","r"))

## i. Interactions
# Chitons ~ Flu.cover | Boulder.area * Wt.level
par(mfrow=c(1,1), mar = c(5,4,2,1))
coplot(Chitons ~ Flu.cover | Boulder.area*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Fluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Flu.cover | Exp.side.area * Wt.level
coplot(Chitons ~ Flu.cover | Exp.side.area*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Fluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Flu.cover | Weight * Wt.level
coplot(Chitons ~ Flu.cover | Weight*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Fluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Flu.cover | Bryo.cover * Weight
coplot(Chitons ~ Flu.cover | Bryo.cover*Weight, data=side,
       ylab = "Abundance",
       xlab = "Fluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Flu.cover | Bryo.cover * Wt.level
coplot(Chitons ~ Flu.cover | Bryo.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Fluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Flu.cover | Asc.cover * Exp.side.area
coplot(Chitons ~ Flu.cover | Asc.cover*Exp.side.area, data=side,
       ylab = "Abundance",
       xlab = "Fluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Flu.cover | Asc.cover * Wt.level
coplot(Chitons ~ Flu.cover | Asc.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Fluorescence cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.cca | Exp.side.area * Wt.level
coplot(Chitons ~ Live.cca | Exp.side.area*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "CCA cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.cca | Weight * Wt.level
coplot(Chitons ~ Live.cca | Weight*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "CCA cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.cca | Bryo.cover * Weight
coplot(Chitons ~ Live.cca | Bryo.cover*Weight, data=side,
       ylab = "Abundance",
       xlab = "CCA cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.cca | Bryo.cover * Wt.level
coplot(Chitons ~ Live.cca | Bryo.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "CCA cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.cca | Asc.cover * Exp.side.area
coplot(Chitons ~ Live.cca | Asc.cover*Exp.side.area, data=side,
       ylab = "Abundance",
       xlab = "CCA cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.cca | Asc.cover * Wt.level
coplot(Chitons ~ Live.cca | Asc.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "CCA cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.pey | Exp.side.area * Wt.level
coplot(Chitons ~ Live.pey | Exp.side.area*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Peyssonnelia cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.pey | Weight * Wt.level
coplot(Chitons ~ Live.pey | Weight*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Peyssonnelia cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.pey | Bryo.cover * Weight

coplot(Chitons ~ Live.pey | Bryo.cover*Weight, data=side,
       ylab = "Abundance",
       xlab = "Peyssonnelia cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) }

# Chitons ~ Live.pey | Bryo.cover * Wt.level
coplot(Chitons ~ Live.pey | Bryo.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Peyssonnelia cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.pey | Asc.cover * Exp.side.area
coplot(Chitons ~ Live.pey | Asc.cover*Exp.side.area, data=side,
       ylab = "Abundance",
       xlab = "Peyssonnelia cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Live.pey | Asc.cover * Wt.level
coplot(Chitons ~ Live.pey | Asc.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Peyssonnelia cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Unflu.cover | Exp.side.area * Wt.level
coplot(Chitons ~ Unflu.cover | Exp.side.area*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Unfluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Unflu.cover | Weight * Wt.level
coplot(Chitons ~ Unflu.cover | Weight*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Unfluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Unflu.cover | Bryo.cover * Weight
coplot(Chitons ~ Unflu.cover | Bryo.cover*Weight, data=side,
       ylab = "Abundance",
       xlab = "Unfluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Unflu.cover | Bryo.cover * Wt.level
coplot(Chitons ~ Unflu.cover | Bryo.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Unfluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Unflu.cover | Asc.cover * Exp.side.area
coplot(Chitons ~ Unflu.cover | Asc.cover*Exp.side.area, data=side,
       ylab = "Abundance",
       xlab = "Unfluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Unflu.cover | Asc.cover * Wt.level
coplot(Chitons ~ Unflu.cover | Asc.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Unfluorescent cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Dead.red.crust | Exp.side.area * Wt.level
coplot(Chitons ~ Dead.red.crust | Exp.side.area*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Dead red crust algae  cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Dead.red.crust | Weight * Wt.level
coplot(Chitons ~ Dead.red.crust | Weight*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Dead red crust algae  cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Dead.red.crust | Bryo.cover * Weight
coplot(Chitons ~ Dead.red.crust | Bryo.cover*Weight, data=side,
       ylab = "Abundance",
       xlab = "Dead red crust algae cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Dead.red.crust | Bryo.cover * Wt.level
coplot(Chitons ~ Dead.red.crust | Bryo.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Dead red crust algae cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Dead.red.crust | Asc.cover * Exp.side.area
coplot(Chitons ~ Dead.red.crust | Asc.cover*Exp.side.area, data=side,
       ylab = "Abundance",
       xlab = "Dead red crust algae cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# Chitons ~ Dead.red.crust | Asc.cover * Wt.level
coplot(Chitons ~ Dead.red.crust | Asc.cover*Wt.level, data=side,
       ylab = "Abundance",
       xlab = "Dead red crust algae cover (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

# 3) Variable statistics ####

# Dataset: "*chitons.data.csv*"
side <- read.csv("chitons.data.csv", sep=";",dec=".", header=T)
side$fSite <- factor(side$Site)
str(side)

## a. Chitons ####

# Total chiton abundance?

sum(side$Chitons) # 137

# Max abundance in a boulder?
max(side$Chitons) # 15

# Percentage of occupied boulders?
boulders <- length(side$Chiton)
occupied <- length(side$Chiton[side$Chiton!=0])
percentage <- (occupied*100)/boulders
percentage # 36.6 % occupied boulders

# Encounter rate
enc.rate <- sum(side$Chitons)/sum(side$Samp.time)
enc.rate # 0.15 chitons/min
enc.rate*60 # 9.35 chitons/hour

# Density
area <- (sum(side$Boulder.area)/100)
chitons <- sum(side$Chitons)
density <- chitons/area
density  # 0.427 individuals/m^2

# Visualizing
par(mfrow=c(1,1), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(side$Chitons~side$fSite, ylab="Abundance",xlab="")

# Maximum:
tapply(side$Chitons,side$Site,max)

# Average:
tapply(side$Chitons,side$Site,mean)

# Standard deviations:
tapply(side$Chitons,side$Site,sd)

# Kruskal-Wallis
kruskal.test(Chitons~fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Chitons ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

## b. Samp.time ####

# Visualizing:
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Samp.time, prob = T, ylim=c(0,0.5),breaks = 10, xlab="Sampling time (min)", main = "Histogram")
lines(density(side$Samp.time), col="blue")
dotchart(side$Samp.time, xlab = "Sampling time (min)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Samp.time)

# Standard Deviation
sd(side$Samp.time)

# Total sampling time
sum(side$Samp.time)

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(side$Samp.time~side$fSite, ylab="Samping time (min)", xlab="")
dotchart(side$Samp.time, xlab = "Sampling time (min)", group=side$fSite)

# Total Sampling time:
tapply(side$Samp.time,side$Site,sum)

# Average Sampling time:
tapply(side$Samp.time,side$Site,mean)

# Standard deviations:
tapply(side$Samp.time,side$Site,sd)

# Formal test:
kruskal.test(Samp.time ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Samp.time ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## c. Temp ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Temp, prob = T, ylim=c(0,0.5),breaks = 10, xlab="Water temperature (°C)", main = "Histogram")
lines(density(side$Temp), col="blue")
dotchart(side$Temp, xlab = "Water temperature (°C)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Temp)

# Standard Deviation
sd(side$Temp)

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(side$Temp~side$fSite, ylab="Water temperature (°C)", xlab="")
dotchart(side$Temp, xlab = "Water temperature (°C)", group=side$fSite)

# Maximum temperature
tapply(side$Temp,side$Site,max)

# Minimum temperature
tapply(side$Temp,side$Site,min)

# Average temperature:
tapply(side$Temp,side$Site,mean)

# Standard deviations:
tapply(side$Temp,side$Site,sd)

# Formal test:
kruskal.test(Temp ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Temp ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## d. Wt.level ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Wt.level, prob = T, ylim=c(0,0.05), breaks = 10, xlab="Water level (mm)", main = "Histogram")
lines(density(side$Wt.level), col="blue")
dotchart(side$Wt.level, xlab = "Water level (mm)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Wt.level)

# Standard Deviation
sd(side$Wt.level)

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(side$Wt.level~side$fSite, ylab="Water level (mm)",xlab="")
dotchart(side$Wt.level, xlab = "Water level (mm)",group=side$fSite)

# Maximum Wt.level:
tapply(side$Wt.level,side$Site,max)

# Minimum Wt.level:
tapply(side$Wt.level,side$Site,min)

# Average Wt.level:
tapply(side$Wt.level,side$Site,mean)
```
# Standard deviations:
tapply(side$Wt.level,side$Site,sd)

# Formal test:
kruskal.test(Wt.level ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Wt.level ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## e. Weight ####

# Visualizing:
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Weight, prob = T, ylim=c(0,0.2), breaks = 10, xlab="Boulder's weight (Kg)", main = "Histogram")
lines(density(side$Weight), col="blue")
dotchart(side$Weight, xlab = "Boulder's weight (Kg)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Weight)

# Standard Deviation
sd(side$Weight)

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Weight~fSite, data=side, ylab="Boulders weight (Kg)",xlab="")
dotchart(side$Weight, xlab = "Boulders weight (Kg)", group=side$fSite)

# Maximum Weight:
tapply(side$Weight,side$Site,max)

# Minimum Weight:
tapply(side$Weight,side$Site,min)

# Average Weight:
tapply(side$Weight,side$Site,mean)

# Standard deviations:
tapply(side$Weight,side$Site,sd)

# Formal test:
kruskal.test(Weight ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Weight ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## f. Boulder.area ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Boulder.area, prob = T, ylim=c(0,0.01),breaks = 15, xlab="Boulder surface area (mm)", main = "Histogram")
lines(density(side$Boulder.area), col="blue")
dotchart(side$Boulder.area, xlab = "Boulder surface area (mm)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Boulder.area)

# Standard Deviation
sd(side$Boulder.area)

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Boulder.area~fSite, data=side, ylab="Boulder surface area (mm)",xlab="")
dotchart(side$Boulder.area, xlab = "Boulder surface area (mm)", group=side$fSite)
```
# Total boulder area:
tapply(side$Boulder.area,side$Site,sum)

# Maximum Boulder.area:
tapply(side$Boulder.area,side$Site,max)

# Minimum Boulder.area:
tapply(side$Boulder.area,side$Site,min)

# Average Boulder.area:
tapply(side$Boulder.area,side$Site,mean)

# Standard deviations:
tapply(side$Boulder.area,side$Site,sd)

# Kruskal-Wallis
kruskal.test(Boulder.area~fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Boulder.area ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## g. Exp.side.area ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Exp.side.area, prob = T, ylim=c(0,0.01),breaks = 10, xlab="Exposed side surface (cm²)", main = "Histogram")
lines(density(side$Exp.side.area), col="blue")
dotchart(side$Exp.side.area, xlab = "Exposed side surface (cm²)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Exp.side.area)

# Standard Deviation
sd(side$Exp.side.area)

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Exp.side.area~fSite, data=side, ylab="Exposed side surface (cm²)",xlab="")
dotchart(side$Exp.side.area, xlab = "Exposed side surface (cm²)", group=side$fSite)

# Total Exposed side surface area:
tapply(side$Exp.side.area,side$Site,sum)

# Maximum Exp.side.area:
tapply(side$Exp.side.area,side$Site,max)

# Minimum Exp.side.area:
tapply(side$Exp.side.area,side$Site,min)

# Average Exp.side.area:
tapply(side$Exp.side.area,side$Site,mean)

# Standard deviations:
tapply(side$Exp.side.area,side$Site,sd)

# Formal test:
kruskal.test(Exp.side.area ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Exp.side.area ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## h. Flu.cover ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Flu.cover, prob = T, ylim=c(0,0.05),breaks = 15, xlab="Fluorescent substrate (%)", main = "Histogram")
lines(density(side$Flu.cover), col="blue")
dotchart(side$Flu.cover, xlab = "Fluorescent substrate (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Flu.cover)

# Standard Deviation
sd(side$Flu.cover)

# Percentage of boulders with red crust algae on the sides?*
boulders <- length(side$Site)
occupied <- length(side$Flu.cover[side$Flu.cover!=0])
percentage <- (occupied*100)/boulders
percentage # Red crust algae is present on the side faces of all boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Flu.cover~fSite, data=side, ylab="Fluorescent substrate (%)",xlab="")

# Maximum:
tapply(side$Flu.cover,side$Site,max)

# Minimum:
tapply(side$Flu.cover,side$Site,min)

# Average fluorescent substrate cover:
tapply(side$Flu.cover,side$Site,mean)

# Standard deviations:
tapply(side$Flu.cover,side$Site,sd)

# Formal test:
kruskal.test(Flu.cover ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Flu.cover ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

## i. Live.cca ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Live.cca, prob = T, ylim=c(0,0.05),breaks = 15, xlab="CCA cover (%)", main = "Histogram")
lines(density(side$Live.cca), col="blue")
dotchart(side$Live.cca, xlab = "CCA cover (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Live.cca)

# Standard Deviation
sd(side$Live.cca)

# Percentage of boulders with live Crustose Coralline Algae on the side face?*

boulders <- length(side$Site)
occupied <- length(side$Live.cca[side$Live.cca!=0])
percentage <- (occupied*100)/boulders
percentage # CCA is present on the side face of all boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Live.cca~fSite, data=side, ylab="CCA cover (%)",xlab="")
dotchart(side$Live.cca, xlab = "CCA cover (%)", group=side$fSite)

# Maximum:
tapply(side$Live.cca,side$Site,max)

# Minimum:
tapply(side$Live.cca,side$Site,min)

# Average:
tapply(side$Live.cca,side$Site,mean)

# Standard deviations:
tapply(side$Live.cca,side$Site,sd)

# Formal test:
kruskal.test(Live.cca ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Live.cca ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## j. Live.pey ####

par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Live.pey, prob = T, ylim=c(0,0.2),breaks = 15, xlab="Peyssonnelia cover (%)", main = "Histogram")
lines(density(side$Live.pey), col="blue")
dotchart(side$Live.pey, xlab = "Peyssonnelia cover (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Live.pey)

# Standard Deviation
sd(side$Live.pey)

# Percentage of boulders with Peyssonnelia crust algae on the sides?*
boulders <- length(side$Site)
occupied <- length(side$Live.pey[side$Live.pey!=0])
percentage <- (occupied*100)/boulders
percentage # Peyssonnelia crusts are present on the side faces of 66% of the boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Live.pey~fSite, data=side, ylab="Peyssonnelia cover (%)",xlab="")
dotchart(side$Live.pey, xlab = "Peyssonnelia cover (%)", group=side$fSite)

# Maximum:
tapply(side$Live.pey,side$Site,max)

# Non-zero minimum :
fill_Live.pey <- side[side$Live.pey !=0,]
tapply(fill_Live.pey$Live.pey,fill_Live.pey$Site,min)

# Average:
tapply(side$Live.pey,side$Site,mean)

# Standard deviations:
tapply(side$Live.pey,side$Site,sd)

# Formal test:
kruskal.test(Live.pey ~ fSite, data=side)

# Dunn's Kruskal-Wallis Multiple Comparisons
DT = dunnTest(Live.pey ~ fSite,
              data=side,
              method="bonferroni")   # Bonferroni correction
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## k. Unflu.cover ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Unflu.cover, prob = T, ylim=c(0,0.05),breaks = 15, xlab="Unfluorescent substrates cover (%)", main = "Histogram")
lines(density(side$Unflu.cover), col="blue")
dotchart(side$Unflu.cover, xlab = "Unfluorescent substrates cover (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Unflu.cover)

# Standard Deviation
sd(side$Unflu.cover)

# Percentage of boulders with red crust algae on the sides?*
boulders <- length(side$Site)
occupied <- length(side$Unflu.cover[side$Unflu.cover!=0])
percentage <- (occupied*100)/boulders
percentage # Unfluorescent substrates are present on the side face of all boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Unflu.cover~fSite, data=side, ylab="Unfluorescent substrates (%)",xlab="")
dotchart(side$Unflu.cover, xlab = "Unfluorescent substrates (%)", group=side$fSite)

# Maximum:
tapply(side$Unflu.cover,side$Site,max)

# Minimum:
tapply(side$Unflu.cover,side$Site,min)

# Average:
tapply(side$Unflu.cover,side$Site,mean)

# Standard deviations:
tapply(side$Unflu.cover,side$Site,sd)

# Formal test:
kruskal.test(Unflu.cover ~ fSite, data=side)

# Comparisons between reefs:
DT = dunnTest(Unflu.cover ~ fSite,
              data=side,
              method="bonferroni")
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## l. Dead.red.crust ####

# Vizualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Dead.red.crust, prob = T, ylim=c(0,0.1), xlab="Dead red crust cover (%)", main = "Histogram")
lines(density(side$Dead.red.crust), col="blue")
dotchart(side$Dead.red.crust, xlab = "Dead red crust cover (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Dead.red.crust)

# Standard Deviation
sd(side$Dead.red.crust)

# Percentage of boulders with dead red crust algae on the sides?*
boulders <- length(side$Site)
occupied <- length(side$Dead.red.crust[side$Dead.red.crust!=0])
percentage <- (occupied*100)/boulders
percentage # Dead red crust is present on the side faces of 95% of the boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Dead.red.crust~fSite, data=side, ylab="Dead red crust (%)",xlab="")
dotchart(side$Dead.red.crust, xlab = "Dead red crust (%)", group=side$fSite)

# Maximum:
tapply(side$Dead.red.crust,side$Site,max)

# Minimum:
tapply(side$Dead.red.crust,side$Site,min)

# Average:
tapply(side$Dead.red.crust,side$Site,mean)

# Standard deviations:
tapply(side$Dead.red.crust,side$Site,sd)

# Formal test:
kruskal.test(Dead.red.crust ~ fSite, data=side)

# Comparisons between reefs:
DT = dunnTest(Dead.red.crust ~ fSite,
              data=side,
              method="bonferroni")
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## m. Asc.cover ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Asc.cover, prob = T, ylim=c(0,0.2), breaks = 10, xlab="Ascidian cover (%)", main = "Histogram")
lines(density(side$Asc.cover), col="blue")
dotchart(side$Asc.cover, xlab = "Ascidian cover(%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Asc.cover)

# Standard Deviation
sd(side$Asc.cover)

# Percentage of boulders with ascidians on the sides?*
boulders <- length(side$Site)
occupied <- length(side$Asc.cover[side$Asc.cover!=0])
percentage <- (occupied*100)/boulders
percentage # Ascidians are present on the side faces of 78.33% of the boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Asc.cover~fSite, data=side, ylab="Ascidian cover (%)",xlab="")
dotchart(side$Asc.cover, xlab = "Ascidian cover (%)", group=side$fSite)

# Maximum:
tapply(side$Asc.cover,side$Site,max)

# Non-zero minimum :
fill_Asc.cover <- side[side$Asc.cover !=0,]
tapply(fill_Asc.cover$Asc.cover,fill_Asc.cover$Site,min)

# Average:
tapply(side$Asc.cover,side$Site,mean)

# Standard deviations:
tapply(side$Asc.cover,side$Site,sd)

# Formal test:
kruskal.test(Asc.cover ~ fSite, data=side)

# Comparisons between reefs:
DT = dunnTest(Asc.cover ~ fSite,
              data=side,
              method="bonferroni")
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## n. Bryo.cover #####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Bryo.cover, prob = T, ylim=c(0,0.5), breaks = 15, xlab="Bryozoan cover (%)", main = "Histogram")
lines(density(side$Bryo.cover), col="blue")
dotchart(side$Bryo.cover, xlab = "Bryozoan cover (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Bryo.cover)

# Standard Deviation
sd(side$Bryo.cover)

# Percentage of boulders with bryozoans on the sides:
boulders <- length(side$Site)
occupied <- length(side$Bryo.cover[side$Bryo.cover!=0])
percentage <- (occupied*100)/boulders
percentage # Bryozoans are present on the side faces of 26.66% of the boulders.

# Percentage of zeros:
sum(side$Bryo.cover==0)/length(side$Bryo.cover)*100 # 73.33%

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Bryo.cover~fSite, data=side, ylab="Bryozoan cover (%)",xlab="")
dotchart(side$Bryo.cover, xlab = "Bryozoan cover (%)", group=side$fSite)

# Maximum:
tapply(side$Bryo.cover,side$Site,max)

# Non-zero minimum :
fill_Bryo.cover <- side[side$Bryo.cover !=0,]
tapply(fill_Bryo.cover$Bryo.cover,fill_Bryo.cover$Site,min)

# Average:
tapply(side$Bryo.cover,side$Site,mean)

# Standard deviations:
tapply(side$Bryo.cover,side$Site,sd)

# Formal test:
kruskal.test(Bryo.cover ~ fSite, data=side)

# Comparisons between reefs:
DT = dunnTest(Bryo.cover ~ fSite,
              data=side,
              method="bonferroni")
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## o. Spong.cover ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Spong.cover, prob = T, ylim=c(0,0.8), breaks = 15, xlab="Sponge cover (%)", main = "Histogram")
lines(density(side$Spong.cover), col="blue")
dotchart(side$Spong.cover, xlab = "Sponge cover (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Spong.cover)

# Standard Deviation
sd(side$Spong.cover)

# Percentage of boulders with sponges on the sides?
boulders <- length(side$Site)
occupied <- length(side$Spong.cover[side$Spong.cover!=0])
percentage <- (occupied*100)/boulders
percentage # Sponges are present on the side faces of 57.5% of the boulders.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Spong.cover~fSite, data=side, ylab="Sponges cover (%)",xlab="")
dotchart(side$Spong.cover, xlab = "Sponges cover (%)", group=side$fSite)

# Maximum:
tapply(side$Spong.cover,side$Site,max)

# Non-zero minimum :
fill_Spong.cover <- side[side$Spong.cover !=0,]
tapply(fill_Spong.cover$Spong.cover,fill_Spong.cover$Site,min)

# Average:
tapply(side$Spong.cover,side$Site,mean)

# Standard deviations:
tapply(side$Spong.cover,side$Site,sd)

# Formal test:
kruskal.test(Spong.cover ~ fSite, data=side)

# Comparisons between reefs:
DT = dunnTest(Spong.cover ~ fSite,
              data=side,
              method="bonferroni")
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


## p. Others ####

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1))
hist(side$Others, prob = T, ylim=c(0,0.5),breaks = 10, xlab="Others coverage (%)", main = "Histogram")
lines(density(side$Others), col="blue")
dotchart(side$Others, xlab = "Others coverage (%)",
         ylab = "Order of the data", main = "Dotplot")

# Summary
summary(side$Others)

# Standard Deviation
sd(side$Others)

# Percentage of boulders with red crust algae on the sides?
boulders <- length(side$Site)
occupied <- length(side$Others[side$Others!=0])
percentage <- (occupied*100)/boulders
percentage # Others coverage are present in 65%.

# Visualizing
par(mfrow=c(1,2), mar = c(5,4,2,1), cex.axis=0.6)
boxplot(Others~fSite, data=side, ylab="Others coverage (%)",xlab="")
dotchart(side$Others, xlab = "Others coverage (%)", group=side$fSite)

# Maximum:
tapply(side$Others,side$Site,max)

# Non-zero minimum :
fill_Others <- side[side$Others !=0,]
tapply(fill_Others$Others,fill_Others$Site,min)

# Average:
tapply(side$Others,side$Site,mean)

# Standard deviations:
tapply(side$Others,side$Site,sd)

# Formal test:
kruskal.test(Others ~ fSite, data=side)

# Comparisons between reefs:
DT = dunnTest(Others ~ fSite,
              data=side,
              method="bonferroni")
PT = DT$res
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


# 4) Question 1 ####

# Objective: Investigate whether the abundance of fluorescent chitons ('*Chitons*')
# is influenced by the fluorescent substrate covering the sides of the boulders they inhabit.

# Hypothesis: The abundance of fluorescent chitons increases with greater
# fluorescent substrate coverage on the sides of boulders.

# Response: '*Chitons*'.
# Predictors: '*Flu.cover*', '*Asc.cover*', '*Bryo.cover*', and '*Spong.cover*'.

# Dataset: "*chitons.data.csv*"
side <- read.csv("chitons.data.csv", sep=";",dec=".", header=T)
side$fSite <- factor(side$Site)
str(side)

# Variance Inflation Factor
Z <- side[,c("Flu.cover","Asc.cover","Bryo.cover","Spong.cover")]
vifcor(Z)

# MODEL BUILDING

# STEP 1: Fitting a 'loaded' mean structure model.

# Adding full interactions:
M1 <- gls(
  Chitons ~ Flu.cover * Asc.cover * Bryo.cover * Spong.cover,
  method = "REML",
  data = side
)

summary(M1)

# STEP 2: Including a structure for the random effects.

# Adding random effect:
M1a <- lme(
  Chitons ~ Flu.cover * Asc.cover * Bryo.cover * Spong.cover,
  random = ~ 1 | fSite,
  method = "REML",
  data = side
)

summary(M1a)

# Akaike information criterion (AIC):
AIC(M1,M1a)

# Bayesian Information Criteria (BIC)
BIC(M1,M1a)

# Model selection:
model.sel(M1,M1a)

# Anova test:*
anova(M1,M1a)

#Log likelihood ratio test:
-2*(-337.913-(-324.673)) # Zuur et al. (2009)

# Correction:
0.5*(1-pchisq(26.48,1)) # Verbeke and Molenberghs (2000)

# STEP 3: Finding the optimal fixed structure

# Full
M1a.1 <- lme(
  Chitons ~  Flu.cover * Asc.cover * Bryo.cover * Spong.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

#Adding three-way interactions
M1a.2 <- lme(
  Chitons ~  Flu.cover * Bryo.cover +
   Flu.cover * Spong.cover +
   Flu.cover * Asc.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

# Add two-way interactions
M1a.3 <- lme(
  Chitons ~
    Flu.cover * Asc.cover +
    Flu.cover * Bryo.cover +
    Flu.cover * Spong.cover +
    Asc.cover * Bryo.cover +
    Asc.cover * Spong.cover +
    Bryo.cover * Spong.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

# without interactions
M1a.4 <- lme(
  Chitons ~  Flu.cover + Asc.cover + Bryo.cover + Spong.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

# Model selection:
model.sel(M1a.1,M1a.2,M1a.3,M1a.4)
anova(M1a.1,M1a.4)

#Log likelihood ratio test:
-2*(-254.5977-(-255.8333)) # Zuur et al. (2009)

# Correction:
0.5*(1-pchisq(-2.4712,1)) # Verbeke and Molenberghs (2000)

summary(M1a.4)

# STEP 4: Reduce the model by removing nonsignificant fixed effects.

# Original model:
M1a.4 <- lme(
  Chitons ~  Flu.cover + Asc.cover + Bryo.cover + Spong.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

drop1(M1a.4, test="Chi")

# Removing 'Asc.cover':
M1a.4a <- lme(
  Chitons ~  Flu.cover + Bryo.cover + Spong.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

drop1(M1a.4a, test="Chi")

# Removing 'Spong.cover':
M1a.4b <- lme(
  Chitons ~ Flu.cover + Bryo.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

drop1(M1a.4b, test="Chi")

# Removing 'Bryo.cover':
M1a.4c <- lme(
  Chitons ~  Flu.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

model.sel(M1a.4,M1a.4a,M1a.4b,M1a.4c)
anova(M1a.4,M1a.4c)

#Log likelihood ratio test
-2*(-255.8333-(-256.6482)) # Zuur et al. (2009)

# Correction:
0.5*(1-pchisq(-1.6298,1)) # Verbeke and Molenberghs (2000)

summary(M1a.4c)

# STEP 5: Adjusting to REML and model validation.

M1a.4c_reml <- lmer(
  Chitons ~ Flu.cover + (1|fSite),
  data=side
)

summary(M1a.4c_reml)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = M1a.4c_reml,
                            n = 1000)
plot(simRes)

# STEP 6: Generalizing to Poisson and reducing the model

P_M1a.4 <- glmmTMB(
  Chitons ~  Flu.cover + Asc.cover + Bryo.cover + Spong.cover + (1|fSite),
  family = poisson,
  data = side
)

summary(P_M1a.4)
drop1(P_M1a.4, test="Chi") # Asc.cover

P_M1a.4a <- glmmTMB(
  Chitons ~  Flu.cover + Bryo.cover + Spong.cover + (1|fSite),
  family = poisson,
  data = side
)

summary(P_M1a.4a)
drop1(P_M1a.4a, test="Chi") # Spong.cover

P_M1a.4b <- glmmTMB(
  Chitons ~  Flu.cover + Bryo.cover + (1|fSite),
  family = poisson,
  data = side
)

summary(P_M1a.4b)

model.sel(P_M1a.4,P_M1a.4a,P_M1a.4b)
anova(P_M1a.4,P_M1a.4b)

#Log likelihood ratio test
-2*(-161.37-(-162.97)) # Zuur et al. (2009)

# Correction:
0.5*(1-pchisq(-3.2,1)) # Verbeke and Molenberghs (2000)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = P_M1a.4b,
                            n = 1000)
plot(simRes)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,2))
plotResiduals(simRes, side$Flu.cover)
plotResiduals(simRes, side$Bryo.cover)

# STEP 7: Incorporing zero inflation

P_M1a.4b.1<- glmmTMB(
  Chitons ~  Flu.cover + Bryo.cover  + (1|fSite),
  family  = poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = P_M1a.4b.1,
                            n = 1000)
plot(simRes)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,2))
plotResiduals(simRes, side$Flu.cover)
plotResiduals(simRes, side$Bryo.cover)

# STEP 7: Final adjustments

P_M1a.4b.2 <- glmmTMB(
  Chitons ~  Flu.cover + sqrt(Bryo.cover) + (1|fSite),
  family  = poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

summary(P_M1a.4b.2)

P_M1a.4a.1 <- glmmTMB(
  Chitons ~  Flu.cover + Bryo.cover + Spong.cover + (1|fSite),
  family = poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

summary(P_M1a.4b.2)

model.sel(P_M1a.4a.1,P_M1a.4b.2)

AIC(P_M1a.4a.1,P_M1a.4b.2)
AICc(P_M1a.4a.1,P_M1a.4b.2)
BIC(P_M1a.4a.1,P_M1a.4b.2)

anova(P_M1a.4a.1,P_M1a.4b.2)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = P_M1a.4b.2,
                            n = 1000)
plot(simRes)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,2))
plotResiduals(simRes, side$Flu.cover)
plotResiduals(simRes, side$Bryo.cover)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = P_M1a.4a.1,
                            n = 1000)
plot(simRes)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,3))
plotResiduals(simRes, side$Flu.cover)
plotResiduals(simRes, side$Bryo.cover)
plotResiduals(simRes, side$Spong.cover)

# MODEL VALIDATION
# Final model: P_M1a.4b.2

P_M1a.4b.2 <- glmmTMB(
  Chitons ~  Flu.cover + sqrt(Bryo.cover) + (1|fSite),
  family  = poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

summary(P_M1a.4b.2)

# Calculating randomized quantile residuals.
simRes <- simulateResiduals(fittedModel = P_M1a.4b.2, n=1000,
                            plot=T)

# Formal tests for over/underdispersion.
testDispersion(simRes)

# Formal tests for zero-inflation.
testZeroInflation(simRes, plot=T)

# Outlier test.
testOutliers(simulationOutput = simRes)

# Detecting missing predictors or wrong functional assumptions.
par(mfrow = c(1,2))
plotResiduals(simRes, side$Flu.cover)
plotResiduals(simRes, side$Bryo.cover)

# 5) Question 2 ####

# Objective: Investigate whether the effect of the fluorescent substrate cover
# on the abundance of fluorescent chiton species can be translated into a
# specie-specific association to the type of red crust algae on the sides of boulders.

# Hypothesis: The fluorescent chiton-substrate association is independent of
# species-specific relationships.

# Response: '*Chitons*'.
# Predictors: '*Live.cca*', '*Live.pey*','*Asc.cover*', '*Bryo.cover*' and '*Spong.cover*'.

# Dataset: "*chitons.data.csv*"
side <- read.csv("chitons.data.csv", sep=";",dec=".", header=T)
side$fSite <- factor(side$Site)
str(side)

# Variance Inflation Factor:
Z <- side[,c("Live.cca","Live.pey","Asc.cover","Bryo.cover","Spong.cover")]
vifcor(Z)

# MODEL BUILDING

# STEP 1: Fitting a 'loaded' mean structure model.
M2 <- gls(
  Chitons ~ Live.cca + Live.pey + Asc.cover + Spong.cover + Bryo.cover,
  method = "REML",
  data = side
)

# STEP 2: Selecting a structure for the random effects.
M2a <- lme(
  Chitons ~ Live.cca + Live.pey + Asc.cover + Spong.cover + Bryo.cover,
  random = ~ 1 | fSite,
  method = "REML",
  data = side
)

# Akaike information criterion:
AIC(M2,M2a)

# Bayesian Information Criteria:
BIC(M2,M2a)

# Anova test:
anova(M2,M2a) # M2a

# Log likelihood ratio test:
-2*(-284.6273-(-269.0895)) # Zuur et al. (2009)

# Correction:
0.5*(1-pchisq(31.0756,1)) # Verbeke and Molenberghs (2000)

# STEP 3: Finding the optimal fixed structure.

# Original model:
M2a <- lme(
  Chitons ~ Live.cca + Live.pey + Asc.cover + Spong.cover + Bryo.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

# Interaction structure:
M2b <- lme(
  Chitons ~ Live.cca*Live.pey + Asc.cover + Spong.cover + Bryo.cover,
  random = ~ 1 | fSite,
  method = "ML",
  data = side
)

# Akaike information criterion:
AIC(M2a,M2b)

# Bayesian Information Criteria:
BIC(M2a,M2b)

# Second-order Akaike Information Criterion
AICc(M2a,M2b)

# Anova test:
anova(M2a,M2b)

# Log likelihood ratio test:
-2*(-255.6756-(-255.6334)) # Zuur et al. (2009)

# Correction:
0.5*(1-pchisq(0.0844,1)) # Verbeke and Molenberghs (2000)

# STEP 4: Generalizing to Poisson and Reducing the Model

P_M2a <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Asc.cover + Spong.cover + Bryo.cover + (1|fSite),
  family=poisson,
  data = side
)

summary(P_M2a)
drop1(P_M2a, test = "Chi")

# Removing Asc.cover:
P_M2b <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Spong.cover + Bryo.cover + (1|fSite),
  family=poisson,
  data = side
)

summary(P_M2b)
drop1(P_M2b, test = "Chi")

# Removing Spong.cover:
P_M2c <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Bryo.cover + (1|fSite),
  family=poisson,
  data = side
)

summary(P_M2c)
drop1(P_M2c, test = "Chi")

model.sel(P_M2a,P_M2b,P_M2c)

anova(P_M2a,P_M2b)
anova(P_M2a,P_M2c)

# STEP 5: Calculating randomized quantile residuals

# A) Model P_M2a:
simRes <- simulateResiduals(fittedModel = P_M2a,
                            n = 1000, plot = TRUE)

# B) Model P_M2b:
simRes <- simulateResiduals(fittedModel = P_M2b,
                            n = 1000, plot = TRUE)

# C)  Model P_M2c:
simRes <- simulateResiduals(fittedModel = P_M2c,
                            n = 1000, plot = TRUE)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(2,2))
plotResiduals(simRes, side$Live.cca)
plotResiduals(simRes, side$Bryo.cover)
plotResiduals(simRes, side$Live.pey)

# STEP 6: Zero inflation correction in Model P_M2c

# Original model:
P_M2c <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Bryo.cover + (1|fSite),
  family=poisson,
  data = side
)

# Zi structure 1:
P_M2c.1 <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Bryo.cover + (1|fSite),
  family=poisson,
  ziformula = ~ Live.pey + Bryo.cover,
  data = side
)

# Zi structure 2:
P_M2c.2 <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Bryo.cover + (1|fSite),
  family=poisson,
  ziformula = ~ Live.pey,
  data = side
)

# Zi structure 3:
P_M2c.3 <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Bryo.cover + (1|fSite),
  family=poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

model.sel(P_M2c,P_M2c.1,P_M2c.2,P_M2c.3)
anova(P_M2c,P_M2c.3)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = P_M2c.3,
                            n = 1000, plot = TRUE)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,3))
plotResiduals(simRes, side$Live.cca)
plotResiduals(simRes, side$Bryo.cover)
plotResiduals(simRes, side$Live.pey)

# STEP 7: Transforming predictors
P_M2c.3a <- glmmTMB(
  Chitons ~ Live.cca + sqrt(Live.pey) + sqrt(Bryo.cover) + (1|fSite),
  family=poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

P_M2c.3b <- glmmTMB(
  Chitons ~ Live.cca + sqrt(Live.pey) + Bryo.cover + (1|fSite),
  family=poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

P_M2c.3c <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + sqrt(Bryo.cover) + (1|fSite),
  family=poisson,
  ziformula = ~ Bryo.cover,
  data = side
)

model.sel(P_M2c.3, P_M2c.3a,P_M2c.3b,P_M2c.3c)
anova(P_M2c.3,P_M2c.3c)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = P_M2c.3c,
                            n = 1000, plot = TRUE)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,3))
plotResiduals(simRes, side$Live.cca)
plotResiduals(simRes, side$Bryo.cover)
plotResiduals(simRes, side$Live.pey)

# STEP 8: Switching to Negative Binomial

# Original fixed structure:
NB_M2 <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Asc.cover + Spong.cover + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

summary(NB_M2)
drop1(NB_M2, test = "Chi")

# Removing Asc.cover:
NB_M2a <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Spong.cover + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

summary(NB_M2a)
drop1(NB_M2a, test = "Chi")

# Removing Spong.cover:
NB_M2b <- glmmTMB(
  Chitons ~ Live.cca + Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

summary(NB_M2b)
drop1(NB_M2b, test = "Chi")

# Removing Live.cca:
NB_M2c <- glmmTMB(
  Chitons ~ Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

summary(NB_M2c)
drop1(NB_M2c, test = "Chi")

# Removing Live.pey:
NB_M2d <- glmmTMB(
  Chitons ~ Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

summary(NB_M2d)

model.sel(NB_M2,NB_M2a,NB_M2b, NB_M2c, NB_M2d)
anova(NB_M2, NB_M2c)
anova(NB_M2, NB_M2d)

# Model Validation:

# A) Model NB_M2c:
simRes <- simulateResiduals(fittedModel = NB_M2c,n = 1000, plot=T)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,2))
plotResiduals(simRes, side$Bryo.cover)
plotResiduals(simRes, side$Live.pey)

# B) Model NB_M2d:
simRes <- simulateResiduals(fittedModel = NB_M2d,n = 1000, plot=T)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,1))
plotResiduals(simRes, side$Bryo.cover)


# Step 9: Introducing polynomial structure in NB_M2d
# Original model:
NB_M2d <- glmmTMB(
  Chitons ~ Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

# Modified model:
NB_M2d.1 <- glmmTMB(
  Chitons ~ I(Bryo.cover^2) + (1|fSite),
  family=nbinom1,
  data = side
)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = NB_M2d.1,n = 1000, plot=T)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,1))
plotResiduals(simRes, side$Bryo.cover)

# Model selection
model.sel(NB_M2d, NB_M2d.1)
anova(NB_M2d, NB_M2d.1)

# Step 10: Proceeding with Model NB_M2c

# Introducing zero inflation:

# Original model
NB_M2c <- glmmTMB(
  Chitons ~ Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

# ZI structure 1:
NB_M2c.1 <- glmmTMB(
  Chitons ~ Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  ziformula = ~ Live.pey + Bryo.cover,
  data = side
)

# ZI structure 2:
NB_M2c.2 <- glmmTMB(
  Chitons ~ Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  ziformula = ~ Live.pey,
  data = side
)

# ZI structure 3:
NB_M2c.3 <- glmmTMB(
  Chitons ~ Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  ziformula = ~ Bryo.cover,
  data = side
)

model.sel(NB_M2c,NB_M2c.1,NB_M2c.2,NB_M2c.3)
anova(NB_M2c, NB_M2c.2)

# Calculating randomized quantile residuals:
simRes <- simulateResiduals(fittedModel = NB_M2c.2,
                            n = 1000, plot = TRUE)

# Square root transformation:

# Original model:
NB_M2c <- glmmTMB(
  Chitons ~ Live.pey + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

# sqrt structure 1:
NB_M2c.1 <- glmmTMB(
  Chitons ~ sqrt(Live.pey) + sqrt(Bryo.cover) + (1|fSite),
  family=nbinom1,
  data = side
)

# sqrt structure 2:
NB_M2c.2 <- glmmTMB(
  Chitons ~ sqrt(Live.pey) + Bryo.cover + (1|fSite),
  family=nbinom1,
  data = side
)

# sqrt structure 3:
NB_M2c.3 <- glmmTMB(
  Chitons ~ Live.pey +sqrt(Bryo.cover) + (1|fSite),
  family=nbinom1,
  data = side
)

model.sel(NB_M2c,NB_M2c.1,NB_M2c.2,NB_M2c.3)
anova(NB_M2c, NB_M2c.3)
anova(NB_M2c, NB_M2c.1)

# Model Validation

# A. Model NB_M2c.3
simRes <- simulateResiduals(fittedModel = NB_M2c.3,
                            n = 1000, plot = TRUE)

# B. Model NB_M2c.1
simRes <- simulateResiduals(fittedModel = NB_M2c.1,
                            n = 1000, plot = TRUE)

# Detecting missing predictors or wrong functional assumptions:
par(mfrow = c(1,2))
plotResiduals(simRes, side$Bryo.cover)
plotResiduals(simRes, side$Live.pey)

# Anova test
anova(NB_M2c.1, NB_M2d)  # NB_M2d

# Our final model
summary(NB_M2d)
