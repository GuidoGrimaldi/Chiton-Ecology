#++++++++++++++++++++++++#
# Doutorado - UFSC
# Chapter 3: Ecology of chitons
#
# Script for: "chitons.eco.data"
#++++++++++++++++++++++++#


# Data Raw ####
usethis::use_data_raw("chitons.eco.data.xlsx")

library(magrittr)

data <- "chitons.eco.data.xlsx" %>%
  readxl::read_xlsx() %>%


# Packages ####
install.packages("MuMIn")
library(DHARMa) # To validate model
library(glmmTMB) # To fit GLMM model zero-inflated
library(car) # For leveneTest and vif model
library(lattice)
library(AICcmodavg)
library(beeswarm) # For exploration graphs
library(MuMIn)

#install.packages("ggplot2")
#install.packages("rlang")
#install.packages("ggeffects")
#install.packages("effects")
library(ggplot2)
library(ggeffects)


#install.packages("glmtoolbox")
library(glmtoolbox)

adjR2(mod1c, digits = 4, verbose = TRUE)

install.packages("HH")
library(HH)
vif(data)

# Importing data ####

data <- read.csv("chitons.eco.data.csv", sep=";",dec=".", header=T)
str(data)
data$fSite <- factor(data$Site)

# Variables ####

# A. Response Variable:
# + "Chiton_F": abundance of fluorescent chitons (Ischonoplax pectinata)

# B. Explanatory Variable:
# + "Site": study sites (factor)
# + "Samp.Time": sampling time (min)
# + "Temp": water temperature (°C)
# + "Wt.level": water level in relation to the base of the boulder (cm)
# + "Weight": boulder weight (Kg)
# + "Area.total": boulder surface area (cm2)
# + "Ara.lateral": exposed side surface area (cm2)
# + "Cov.Flu": coverage of living fluorescent substrates on the side of the boulder (%)
# + "Cov.DeathCCAPey": coverage of no living fluorescent substrates (%)
# + "Cov.Others": coverage of other substrates, e.g. bare rock, 'mucus', others (%)
# + "Cov.Asc": ascidians cover (%)
# + "Cov.Bryo": briozoans cover (%)
# + "Cov.Spong": sponges cover (%)

summary(data)

# Data exploration ####

# Missing data ?
summary(data) # NO

# Balanced sampling ?
tapply(data$Chiton_F,data$fSite,length) # YES

# Percentage of occupied boulders?

boulders <- length(data$Chiton_F)
occupied <- length(data$Chiton_F[data$Chiton_F!=0])
percentage <- (occupied*100)/boulders
percentage # 36.6 % occupied boulders

# Chapman 2005 encontrou valor similar (30% de ocupacao)

# Total chiton abundance:
sum(data$Chiton_F) # 137

# Chiton abundance/ reef:
tapply(data$Chiton_F,data$fSite,sum) # Buzios eh a maior abundancia

# Zero trouble Y ?
table(data$Chiton_F)
plot(table(data$Chiton_F))
barplot(table(data$Chiton_F)) # Parece bastante zeros!

# Podemos saber a proporcao de zeros aplicando o seguinte calculo:
sum(data$Chiton_F==0)/length(data$Chiton_F)*100

# Proporcao dos dados que sao iguais a zero eh 63.33%
# Pode ser possível o uso de um modelos inflacionado por zeros.


# Vamos olhar em um histograma com a distribuicao de probabilidade
hist(data$Chiton_F, prob = T, breaks=15)
rug(data$Chiton_F)
lines(density(data$Chiton_F), col="blue")

# Nossa variavel resposta apresenta uma distribuicao unimodal nao normal, ja esperada para
# dados discretos (contagem).
# Mas  vamos examinar...

# Normality Y ?

# Qual seria a distribuicao normal teorica a partir da media e desvio padrao da
# abundancia observada ?
curve(dnorm(x, mean=mean(data$Chiton_F), sd=sd(data$Chiton_F)), add=T, col="red")

# Possivelmente os dados estao inflados por zeros!! Mas nao nescessariamente.
# Ver: Warton, D. I. (2005). Many zeros does not mean zero inflation: comparing the goodness-of-fit of parametric models to multivariate abundance data. Environmetrics 16(3), 275-289

# Vamos ver em um QQ-plot:
qqnorm(data$Chiton_F)
qqline(data$Chiton_F, col="blue", lwd=2)

# Noossa!! Nem um pouco normal!

# Vamos tentar ajustar a variancia a uma distribuicao de Poisson com lambda em
# referencia a media

# QQ-plot com distribuicao Poisson
library(car)
qqPlot(data$Chiton_F, distribution="pois", lambda=mean(data$Chiton_F))

# Os dados caem fora do limite da distribuicao de Poisson.
# E bem possivel que um glm Poisson, nao se ajuste bem aos dados.
# Provavelmente uma Binomial Negativa ressolva

# Apenas para conferir, segue o teste de normalidade...

# Shapiro test:
# A hipotese nula testada eh de que os dados da variavel possuem uma distribuicao
# normal.
# Dessa forma, um p significativo (< 0.05) indicara uma rejeicao da hipotese nula,
# aceitando-se a hipotese alternativa de nao normalidade dos dados.

# Vejamos:
shapiro.test(data$Chiton_F) # NAO NORMAL

# Conforme vizualizados, nossa variavel resposta nao possui uma distribuicao normal,
# portanto ha de se esperar que seus erros tambem nao sigam uma distribuicao normal

# Homogeneity test:
variancia <- tapply(data$Chiton_F,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [1]  # A variancia de Buzios eh 66x maior que BF

# Teste de Bartlett:
bartlett.test(data$Chiton_F ~ data$fSite) # H0 aceita (Homocedastica)

# Teste de Fligner-Killeen:
fligner.test(data$Chiton_F ~ data$fSite) # H0 aceita (Homocedastica)

library(car)
leveneTest(Chiton_F ~ fSite, data=data, center=mean) # H0 aceita (Homocedastica)

# Ha diferenca entre praias?

kruskal.test(Chiton_F ~ fSite,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:
library(FSA)
DT = dunnTest(Chiton_F ~ fSite,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:
PT = DT$res
PT
library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

# De fato, Buzio eh diferente das demais prais



# ~Sampling Time ####
tapply(data$Samp.time,data$f.Site,mean)
tapply(data$Samp.time,data$f.Site,sd)

# Normality test:
shapiro.test(data$Samp.time) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Samp.time,data$f.Site,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [3]  # Buzios eh 2x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Samp.time ~ data$f.Site) # H0 aceita (Homocedastica)

# Teste de Fligner-Killeen:
fligner.test(data$Samp.time ~ data$Site) # H0 aceita (Homocedastica)

library(car)
leveneTest(Samp.time ~ f.Site, data=data, center=mean) # H0 aceita (Homocedastica)

# Kruaskal-Wallis test:

kruskal.test(Samp.time ~ f.Site,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:
library(FSA)
DT = dunnTest(Samp.time ~ f.Site,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:
PT = DT$res
PT
library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)


# ~Temperature ####

## Global Water Temperature
summary(data$Temp)
mean(data$Temp)
sd(data$Temp)

  ## Water Temperature by Reef
tapply(data$Temp,data$Site,mean)
tapply(data$Temp,data$Site,sd)

boxplot(Temp~Site, data = data)


# Normality test:
shapiro.test(data$Temp) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Temp,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [1]/ variancia [2]  # BF eh 3.66x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Temp ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Temp ~ data$fSite) # H0 rejeitada (Heterocedastico)

library(car)
leveneTest(Temp ~ fSite, data=data, center=mean) # H0 rejeitada (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(Temp ~ Site,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(Temp ~ fSite,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

# ~Wt.level ####

mean(data$Wt.level)
sd(data$Wt.level)

min(data$Wt.level)

tapply(data$Wt.level,data$fSite,mean)
tapply(data$Wt.level,data$fSite,sd)

# Normality test:
shapiro.test(data$Wt.level) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Wt.level,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [4]  # Buzios eh 20x maior...provavelmente heterocedastico

# Teste de Bartlett:
bartlett.test(data$Wt.level ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Wt.level ~ data$fSite) # H0 rejeitada (Heterocedastico)

library(car)
leveneTest(Wt.level ~ fSite, data=data, center=mean) # H0 rejeitada (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(Wt.level ~ Site,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(Wt.level ~ fSite,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res
PT

library(rcompanion)
cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

# ~Weight ####
mean(data$Weight)
sd(data$Weight)
min(data$Weight)

tapply(data$Weight,data$fSite,mean)
tapply(data$Weight,data$fSite,sd)

# Normality test:
shapiro.test(data$Weight) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Weight,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [4]  # Buzios eh 1.83.66x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Weight ~ data$fSite) # H0 Aceitada (Homocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Weight ~ data$fSite) # H0 Aceitada (Homocedastico)

library(car)
leveneTest(Weight ~ fSite, data=data, center=mean) # H0 Aceitada (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Weight ~ Site,
             data = data) # H0 aceita, nao diferenca entre praias

# ~Area Total ####

tapply(data$Area.total,data$fSite,mean)
tapply(data$Area.total,data$fSite,sd)

# Normality test:
shapiro.test(data$Area.total) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Area.total,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [1]  # Buzios eh 3x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Area.total ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Area.total ~ data$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Area.total ~ fSite, data=data, center=mean) # H0 aceita (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Area.total ~ Site,
             data = data) # H0 aceita, nao ha diferenca entre praias


# ~Exposed.area ####
tapply(data$Expo.area,data$fSite,mean)
tapply(data$Expo.area,data$fSite,sd)

tapply(data$Area.total,data$fSite,mean)
tapply(data$Area.total,data$fSite,sd)

# Normality test:
shapiro.test(data$Expo.area) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Expo.area,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [1]  # Buzios eh 3x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Expo.area ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Expo.area ~ data$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Expo.area ~ fSite, data=data, center=mean) # H0 aceita (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Expo.area ~ Site,
             data = data) # H0 aceita, nao ha diferenca entre praias

# ~Cov.Flu ####
tapply(data$Cov.Flu,data$fSite,mean)
tapply(data$Cov.Flu,data$fSite,sd)

# Normality test:
shapiro.test(data$Cov.Flu) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Cov.Flu,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [3]/ variancia [2]  # Pitangui eh 2.35x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Cov.Flu ~ data$fSite) # H0 aceita (Homocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Cov.Flu ~ data$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Cov.Flu ~ fSite, data=data, center=mean) # H0 aceita (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Cov.Flu ~ Site,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(Cov.Flu ~ fSite,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



# ~Cov.Asc ####
tapply(data$Cov.Asc,data$fSite,mean)
tapply(data$Cov.Asc,data$fSite,sd)

# Normality test:
shapiro.test(data$Cov.Asc) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Cov.Asc,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [3]/ variancia [1]  # Pitangui eh 11x maior que BF

# Teste de Bartlett:
bartlett.test(data$Cov.Asc ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Cov.Asc ~ data$fSite) # H0 rejeitada (Heterocedastico)

library(car)
leveneTest(Cov.Asc ~ fSite, data=data, center=mean) # H0 rejeitada (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(Cov.Asc ~ Site,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(Cov.Asc ~ fSite,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)



# ~Cov.Bryo ####
tapply(data$Cov.Bryo,data$fSite,mean)
tapply(data$Cov.Bryo,data$fSite,sd)

# Normality test:
shapiro.test(data$Cov.Bryo) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Cov.Bryo,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [4]/ variancia [2]  # Santa Rita eh 13x maior que Buzios

# Teste de Bartlett:
bartlett.test(data$Cov.Bryo ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Cov.Bryo ~ data$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Cov.Bryo ~ fSite, data=data, center=mean) # H0 Rejeitada (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(Cov.Bryo ~ Site,
             data = data) # H0 aceita, nao ha diferenca entre praias


# ~Cov.Spong ####
tapply(data$Cov.Spong,data$fSite,mean)
tapply(data$Cov.Spong,data$fSite,sd)

# Normality test:
shapiro.test(data$Cov.Spong) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Cov.Spong,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [3]/ variancia [1]  # Pitangui eh 28x maior que BF

# Teste de Bartlett:
bartlett.test(data$Cov.Spong ~ data$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(data$Cov.Spong ~ data$fSite) # H0 rejeitada (Heterocedastico)

library(car)
leveneTest(Cov.Spong ~ fSite, data=data, center=mean) # H0 rejeitada (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(Cov.Spong ~ Site,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(Cov.Spong ~ fSite,
              data=data,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)




# 5. Outliers Y & X ?

# A ferramenta grafica tipica utilizada para deteccao de outliers eh o boxplot.

# Vamos analisar por partes:
# Quitons:
boxplot(data$Chiton_F)

# Coberturas:
boxplot(data[,c(8:13)])

# Abiotics:
boxplot(data[,c(4:7)])

# All:
boxplot(data[,c(2:13)])

head(data)
# Vamos examinar um pouco mais atraves do 'cleveland dotplot':
x11()
par(mfrow= c (4,4))

library(lattice)
Z <- cbind(data$Chiton_F,
           data$Samp.time,
           data$Temp,
           data$Wt.level,
           data$Weight,
           data$Expo.area)

colnames(Z) <- c("Chiton F",
                 "Samp.time",
                 "Temp",
                 "Wt.level",
                 "Weight",
                 "Expo.area")

dotplot(as.matrix(Z), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data")

# Uhh...parece haver um outlier em Samp.time

Z <- cbind(data$Cov.Flu,
           data$Cov.DeathCCAPey,
           data$Cov.Others,
           data$Cov.Asc,
           data$Cov.Bryo,
           data$Cov.Spong)

colnames(Z) <- c("Cov.Fluo",
                 "Cov.DeathCCAPey",
                 "Cov.Others",
                 "Cov.Asci",
                 "Cov.Bryo",
                 "Cov.Spong")

dotplot(as.matrix(Z), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data")


# Boxplots condicionados por praia
boxplot(Weight ~ Site,data=data)
boxplot(Weight ~ Site, varwidth = TRUE, data=data)
boxplot(Temp ~ Site, varwidth = TRUE, data=data)
boxplot(Wt.level ~ Site, varwidth = TRUE, data=data)
boxplot(Expo.area ~ Site, varwidth = TRUE, data=data)
boxplot(Samp.time ~ Site, varwidth = TRUE, data=data)
boxplot(Cov.Flu ~ Site, varwidth = TRUE, data=data)
boxplot(Cov.DeathCCAPey ~ Site, varwidth = TRUE, data=data)
boxplot(Cov.Asc ~ Site, varwidth = TRUE, data=data)
boxplot(Cov.Spong ~ Site, varwidth = TRUE, data=data)
boxplot(Cov.Bryo ~ Site, varwidth = TRUE, data=data)
boxplot(Rug.base ~ Site, varwidth = TRUE, data=data)
tapply(data$Rug.base,data$f.Site,sum)
boxplot(Brittle ~ Site, varwidth = TRUE, data=data)
tapply(data$Brittle,data$f.Site,sum)

library(beeswarm)
beeswarm::beeswarm(Weight ~ Site,data=data)
beeswarm::beeswarm(Temp ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Wt.level ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Expo.area ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Samp.time ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Cov.Flu ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Cov.DeathCCAPey ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Cov.Asc ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Cov.Spong ~ Site,data=data, col="red", pch=16, method="swarm")
beeswarm::beeswarm(Cov.Bryo ~ Site,data=data, col="red", pch=16, method="swarm")


stripchart(Weight ~ Site,data=data, method="stack")
stripchart(Temp ~ Site,data=data, method="stack")
stripchart(Wt.level ~ Site,data=data, method="stack")
stripchart(Expo.area ~ Site,data=data, method="stack")
stripchart(Samp.time ~ Site,data=data, method="stack")
stripchart(Cov.Flu ~ Site,data=data, method="stack")
stripchart(Cov.DeathCCAPey ~ Site,data=data, method="stack")
stripchart(Cov.Asc ~ Site,data=data, method="stack")
stripchart(Cov.Spong ~ Site,data=data, method="stack")
stripchart(Cov.Bryo ~ Site,data=data, method="stack")

# Podemos observar a dispersao dos dados para cada praia
# Ou com a funcao do pacote 'lattice' conseguimos ver a dispersao por praia
library(lattice)
xyplot(Chiton_F ~ Weight | Site, data=data)
xyplot(Chiton_F ~ Weight | Site, data=data, type=c("p","r")) # com linha de regressao
xyplot(Chiton_F ~ Temp | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Wt.level | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Expo.area | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Samp.time | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Cov.Flu | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Cov.DeathCCAPey | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Cov.Asc | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Cov.Spong | Site, data=data, type=c("p","r"))
xyplot(Chiton_F ~ Cov.Bryo | Site, data=data, type=c("p","r"))


# 6. Homogeneity Y ?

# Vamos verificar mais uma das premissas para a regressao linear.
# Para nossa unica variavel resposta.

variancia <- tapply(data$Chiton_F,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [1]

# Nossa maior variancia excede a menor em apenas 66.356 vezes.

# Os dads nao parecem ser homocedasticos...vamos realizar os testes:

# Para testar a hipotese nula (H0) de igualdade das variancias foram propostos
# alguns testes, chamados de Testes de Homocedasticidade das variancias. Um deles
# eh o Teste de Bartlett, um dos mais usados, mas que tende a mascarar diferencas
# que existem quando a curtose eh negativa, e encontrar diferencas que nao
# existem quando a curtose eh positiva.

# Teste de Bartlett:
bartlett.test(data$Chiton_F ~ data$fSite)

# Podemos observar que a hipotese nula foi rejeitada, ou seja, as variancias sao
# heterocedastica.

# Devido as ressalvas do Teste de Bartlett, alguns autores
# recomendam outros testes a serem aplicados para verificar a homocedasticidade
# das variancias.

# o Teste de Fligner-Killeen que nao eh tao sensivel a outliers como o teste de
# Bartlett.

# Teste de Fligner-Killeen:
fligner.test(data$Chiton_F ~ data$fSite)

# Rejeitamos a H0 com o Teste de Fligner-Killen. Ou seja as variancias sao
# heterocedasticas.

# Outro teste sugerido nos livros de estatistica eh o Teste de Levene como
# substituto ao Teste de Bartlett.

# No R, para usarmos o Teste de Levene precisamos instalar o pacote 'car'.

#install.packages("car")
library(car)
leveneTest(Chiton_F ~ fSite, data=data) # H0 rejeitada, dados heterocedasticos

# Usei a opcao 'center=median' para dar mais robustei a analise
# O help da funcao leveneTest() explica que o teste usando a mediana (default da
# funcao) torna-se mais robusto).

leveneTest(Chiton_F ~ fSite, data=data, center=mean)

# mesmo resultado usando a media. HETEROCEDASTICIDADE

# 7. Relationship X & Y ?

# Vamos analisar por partes:

# Var.abioticas primeiro
head(data)
pairs(data[,c(2:7)])

# Funcao encontrada no help da funcao 'pairs' para plotar histogramas
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# Agora plotamos com histograma nos paineis diagonais
pairs(data[,c(2:7)], diag.panel=panel.hist)

# Outra funcao para fornecer os coeficientes de correlacao
# com a fonte proporcional ao indice
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}

# Agora plotamos com indices de correlacao e histogramas
pairs(data[,c(2:7)], diag.panel=panel.hist, lower.panel=panel.cor)

# Incluimos tambem a linha de suavizacao
pairs(data[,c(2:7)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Temos uma Correlacao Positiva entre Weight X Expo.area (0.54)
# Temos uma Correlacao Positiva entre Weight X Samp.time (0.50)
# Temos uma Correlacao Positiva fraca entre Samp.time X Expo.area (0.38)
# Temos uma Correlacao Negativa fraca entre Temp X Wt.level (-0.32)


# Vamos ver as coberturas laterais...

pairs(data[,c(2,8:13)])

# Funcao encontrada no help da funcao 'pairs' para plotar histogramas
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

# Agora plotamos com histograma nos paineis diagonais
pairs(data[,c(2,8:13)], diag.panel=panel.hist)

# Outra funcao para fornecer os coeficientes de correlacao
# com a fonte proporcional ao indice
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}

# Agora plotamos com indices de correlacao e histogramas
pairs(data[,c(2,8:13)], diag.panel=panel.hist, lower.panel=panel.cor)

# Incluimos tambem a linha de suavizacao
pairs(data[,c(2,8:13)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Correlacao entre: Cov.Flu x Cov.DeathCCAPey (-0.36)
# Correlacao entre: Cov.Flu x Cov.Others (-0.71)
# Correlacao entre: Cov.Flu x Cov.Spong (-0.39)

# Vamos em ver todas:
x11()
pairs(data[,c(3:13)])

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

pairs(data[,c(3:13)], diag.panel=panel.hist)

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1)
}

pairs(data[,c(3:13)], diag.panel=panel.hist, lower.panel=panel.cor)

pairs(data[,c(3:13)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Correlacao entre: Temp x Cov.Spong (-0.32)
# Correlacao entre: Wt.level x Cov.Others (-0.44)
# Correlacao entre: Weight x Cov.Spong (0.32)

# Vamos determinar a signifiancia das correlacoes...

# Para isso precisamos testar primero as normalidades

shapiro.test(data$Expo.area) # Nao normal
shapiro.test(data$Weight)   # Nao normal
shapiro.test(data$Wt.level) # Nao normal
shapiro.test(data$Temp)   # Nao normal
shapiro.test(data$Samp.time) # Nao normal
shapiro.test(data$Cov.Flu) # Nao normal
shapiro.test(data$Cov.DeathCCAPey) # Nao normal
shapiro.test(data$Cov.Spong) # Nao normal
shapiro.test(data$Cov.Others) # Nao normal

# portanto....Correlacao de Spearman:

# 1. Weight X Expo.area (0.54):
cor(data$Expo.area, data$Weight)
cor.test(data$Expo.area, data$Weight, method= "spearman") # Significante
plot(data$Expo.area, data$Weight)

# 2. Weight X Samp.time (0.50)
cor.test(data$Weight, data$Samp.time, method= "spearman") # Significante
plot(data$Weight, data$Samp.time)

# 3. Samp.time X Expo.area (0.38)
cor.test(data$Expo.area, data$Samp.time, method= "spearman") # Significante
plot(data$Expo.area, data$Samp.time)

# 4. Temp X Wt.level (-0.32)
cor.test(data$Temp, data$Wt.level, method= "spearman") # Significante
plot(data$Temp, data$Wt.level)

# 5. Cov.Flu x Cov.DeathCCAPey (-0.36)
cor.test(data$Cov.Flu, data$Cov.DeathCCAPey, method= "spearman") # Significante
plot(data$Cov.Flu, data$Cov.DeathCCAPey)

# 6. Cov.Flu x Cov.Others (-0.71)
cor.test(data$Cov.Flu, data$Cov.Other, method= "spearman") # Significante
plot(data$Cov.Flu, data$Cov.Other)

# 7. Cov.Flu x Cov.Spong (-0.39)
cor.test(data$Cov.Flu, data$Cov.Spong, method= "spearman") # Significante
plot(data$Cov.Flu, data$Cov.Spong)

# 8. Temp x Cov.Spong (-0.32)
cor.test(data$Temp, data$Cov.Spong, method= "spearman") # Significante
plot(data$Temp, data$Cov.Spong)

# 9. Wt.level x Cov.Others (-0.44)
cor.test(data$Wt.level, data$Cov.Spong, method= "spearman") # NAO
plot(data$Wt.level, data$Cov.Spong)

# 10. Weight x Cov.Spong (0.32)
cor.test(data$Weight, data$Cov.Spong, method= "spearman") # Significante
plot(data$Weight, data$Cov.Spong)

# Discritivo
summary(data)



aggregate(data$Samp.time, list(data$Site), FUN=mean)

anova(aov(Samp.time ~ Site, data=data))
TukeyHSD(aov(Samp.time ~ Site, data=data),ordered=T)

# ~ Diagnostico ####

# - Nao ha dados faltantes
# - Zero inflado!
# - Nao normal
# - Heterocedasticos
# - Variaveis X colineares:
#   + Cov.Flu x Cov.DeathCCAPey (-0.36)
#   + Cov.Flu x Cov.Others (-0. 71)
#   + Cov.Flu x Cov.Spong (-0.39)
#   + Weight x Expo.area (0.54)
#   + Temp. X Wt.level (-0.32)


# PCA ####

env <- read.csv("pca.csv", sep=";",dec=".", header=T)
head(env)

library('vegan')
env.z <- decostand(env[,-8],method="standardize")

x11()
par(mfrow=c(2,1))
boxplot(env[,-8], col="bisque", main="Boxplot sem transformacao")
boxplot(env.z, col="bisque", main="Boxplot com estandardizacao")

euc <- vegdist(env.z, method="euc")
head(as.matrix(euc))
euc.s <- head(as.matrix(1-euc))

pca <- rda(env[,-8], scale=TRUE)
pca
summary(pca)

summary(pca)$species
summary(pca)$sites
pca$CA$eig

sum(pca$CA$eig)

explic <- (pca$CA$eig)/(sum(pca$CA$eig))*100
explic

## Seleção do númro e eixos
### Kaiser-Guttman

ev <- pca$CA$eig
ev

mean(ev)
ev[ev > mean(ev)]

### Broken stick

n <- length(ev)
bs.ev <- bstick(n,n)
bs.ev

### Visualizando

x11()
par(mfrow=c(2,1))
barplot(ev, main="Eigenvalues - Kaiser-Guttman criterion",
        col="bisque", las=2, ylim=c(0,120))
abline(h=mean(ev), col="red") #average eigenvalue
legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
barplot(t(cbind(ev,bs.ev)),beside=TRUE, col=c("bisque",2),
        main="Eigenvalues - Broken stick model", las=2, ylim=c(0,120))
legend("topright", c("eigenvalue", "Broken stick model"), pch=15,
       col=c("bisque",2), bty="n")

## Grafico PCA
x11()
par(mfrow=c(1,2))
biplot(pca,scaling=1,main="PCA - scaling 1")
biplot(env.std.pca,main="PCA - scaling 2")

site.sc <- scores(pca,display="wa",choices=1:2)
env.sc <- scores(pca,display="sp",choices=1:2)

site.sc <- as.data.frame(site.sc)
env.sc <- as.data.frame(env.sc)

site <- env[,8]
site.sc <- cbind(site.sc,site)
head(site.sc)

x11()
par(mar=c(5,5,2,2))
plot(pca,type="n",cex.lab=1.4,font.lab=2,cex.axis=1.2,
     xlab=paste("PC 1 - ",round(explic[1],1),"%",sep=""),
     ylab=paste("PC 2 - ",round(explic[2],1),"%",sep=""))
axis(1,lwd=2,labels=F)
axis(2,lwd=2,labels=F)
box(lwd=2)
points(site.sc$PC1[(site.sc$site=="Buzios")],
       site.sc$PC2[(site.sc$site=="Buzios")],
       pch=1,col="red",lwd=2,cex=1.6)

points(site.sc$PC1[(site.sc$site=="Pitangui")],
       site.sc$PC2[(site.sc$site=="Pitangui")],
       pch=2,col="darkorange",lwd=2,cex=1.6)

points(site.sc$PC1[(site.sc$site=="Baia Formosa")],
       site.sc$PC2[(site.sc$site=="Baia Formosa")],
       pch=3,col="forestgreen",lwd=2,cex=1.6)

points(site.sc$PC1[(site.sc$site=="Santa Rita")],
       site.sc$PC2[(site.sc$site=="Santa Rita")],
       pch=4,col="blue",lwd=2,cex=1.6)

arrows(x0=0,y0=0,x1=env.sc$PC1,y1=env.sc$PC2,length=0.15,angle=20,
       col="royalblue")

text(env.sc$PC1[(env.sc$PC1<0)&(env.sc$PC2>0)],
     env.sc$PC2[(env.sc$PC1<0)&(env.sc$PC2>0)],
     labels=rownames(env.sc)[(env.sc$PC1<0)&(env.sc$PC2>0)],
     col="royalblue",pos=3,offset=0.3)
text(env.sc$PC1[(env.sc$PC1<0)&(env.sc$PC2<0)],
     env.sc$PC2[(env.sc$PC1<0)&(env.sc$PC2<0)],
     labels=rownames(env.sc)[(env.sc$PC1<0)&(env.sc$PC2<0)],
     col="royalblue",pos=1,offset=0.3)
text(env.sc$PC1[env.sc$PC1>0],env.sc$PC2[env.sc$PC1>0],
     labels=rownames(env.sc)[env.sc$PC1>0],pos=4,offset=0.3,
     col="royalblue")

legend("topright",legend=c("Buzios","Pitangui","Baia Formosa","Santa Rita"),
       pch=c(1,2,3,4),col=c("red","darkorange","forestgreen","blue"),
       bty="n",pt.cex=1.6,cex=1.2,pt.lwd=2)
legend("topleft",legend=c("Buzios","Pitangui","Baia Formosa","Santa Rita"),
       pch=c(1,2),col="black",bty="n",pt.cex=1.6,cex=1.2,pt.lwd=2)


env.std=decostand(env[,-8], method="standardize")
env.std.pca <- rda(env.std)
env.std.pca




# MODELOS
## At boulder side ####

# O objetivo eh investigar se a abundancia de quitons fluorescentes ('Chiton_F') eh influenciada por fatores ambientais e bióticos
# presentes/ atuantes na escala espacial da lateral da rocha no quais os quítons habitam.

# Hipótese: eh que a cobertura de substratos fluorescentes (cca+peyssonnelia) será um fator importante para explica a distribuição spacial dos quítons fluorescentes.

data <- read.csv("rockside.csv", sep=";",dec=".", header=T)
str(data)
data$fSite <- factor(data$Site)

## 1° Formula ####
# Considerando apenas os parametros ambientais e cobertura do substrato:

library(pscl)
f1 <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.DeathCCAPey + Cov.Asc + Cov.Bryo + Cov.Spong + (1|fSite))

# GLMM Poisson:

# Modelo completo
library(glmmTMB)
mod1 <- glmmTMB(f1,family = poisson, data=data)
summary(mod1)

drop1(mod1, test="Chi")

# Retirando o 'Cov.Asc' do Mod1, obtemos o modelo mod1a:
f1a <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.DeathCCAPey + Cov.Bryo + Cov.Spong + (1|fSite))
mod1a <- glmmTMB(f1a,family = poisson, data=data)
summary(mod1a)
drop1(mod1a, test="Chi")

# Retirando o 'Cov.DeathCCAPey' do Mod1a, obtemos o modelo Mod1b:
f1b <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo + Cov.Spong + (1|fSite))
mod1b <- glmmTMB(f1b,family = poisson, data=data)
summary(mod1b)
drop1(mod1b, test="Chi")

# Retirando o 'Cov.Spong' do Mod1b, obtemos o modelo Mod1c:
f1c <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo + (1|fSite))
mod1c <- glmmTMB(f1c,family = poisson, data=data)
summary(mod1c)

# Vamos comparar os modelos:

library(MuMIn)
model.sel(mod1, mod1a, mod1b, mod1c)

# vendo a diferenca entre os modelos pela ANOVA:
anova(mod1, mod1a, mod1b, mod1c)

# Podemos considerar tanto o mod1c cquanto o mod1b, mas vamos coniderar o mod1c por ter o maior peso na explicação os modelos e ser o mais parcimonioso

r.squaredGLMM(mod1c)

install.packages("AICcmodavg")
library(AICcmodavg)

#setup a subset of models of Table 1
Cand.models <- list(mod1, mod1a, mod1b, mod1c)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)


#Vamos analisar o vif do modelo:
# Calculo da inflacao da variancia para fatores lineares vif - fator de inflacao - menor que 3

library(car)
vif(mod1c)

# Validacao:
install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod1c)
plot(simulationOutput)

# Teste para inflacao de zero
testZeroInflation(simulationOutput)

# b. Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,2))
plotResiduals(simulationOutput, data$Cov.Flu)
plotResiduals(simulationOutput, data$Area.lateral)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod1c,
                            n = 1000)
plot(simRes) # good


## 2° Formula ####
# Considerando a presença de esécies de outros taxons na lateral da rocha:

f2 <- formula(Chiton_F ~ Bryo + Asci + Crabs + Worms + Spong +
                Urchins + Zoanthus + Hermit + Fissu + (1|fSite))


library(tidyverse)
library(vegan) # decostand() function

# Tranformando os dados de abundancias de coespecies em presenca-ausencia:
data2 <- decostand(data[,15:29], "pa") %>% cbind(data[,1:14])
head(data2)
data2$fSite <- factor(data$Site)

# Modelo completo:

mod2 <- glmmTMB(f2,family = poisson, data=data2)
summary(mod2)

drop1(mod2, test="Chi")

# Retirando 'Spong' chegamos ao mod2a:
f2a <- formula(Chiton_F ~ Bryo + Asci + Crabs + Worms +
                Urchins + Zoanthus + Hermit + Fissu + (1|fSite))
mod2a <- glmmTMB(f2a,family = poisson, data=data2)
summary(mod2a)
drop1(mod2a, test="Chi")

# Retirando 'Worms' chegamos ao mod2b
f2b <- formula(Chiton_F ~ Bryo + Asci + Crabs +
                 Urchins + Zoanthus + Hermit + Fissu + (1|fSite))
mod2b <- glmmTMB(f2b,family = poisson, data=data2)
summary(mod2b)
drop1(mod2b, test="Chi")

# Retirando 'Hermit' chegamos ao mod2c:
f2c <- formula(Chiton_F ~ Bryo + Asci + Crabs +
                 Urchins + Zoanthus  + Fissu + (1|fSite))
mod2c <- glmmTMB(f2c,family = poisson, data=data2)
summary(mod2c)
drop1(mod2c, test="Chi")

# Retirando 'Asci' chegamos ao mod2d:
f2d <- formula(Chiton_F ~ Bryo + Crabs + Urchins + Zoanthus  + Fissu + (1|fSite))
mod2d <- glmmTMB(f2d,family = poisson, data=data2)
summary(mod2d)
drop1(mod2d, test="Chi")

# Retirando 'Bryo' chegamos ao mod2e:
f2e <- formula(Chiton_F ~ Crabs + Urchins + Zoanthus  + Fissu + (1|fSite))
mod2e <- glmmTMB(f2e,family = poisson, data=data2)
summary(mod2e)
drop1(mod2e, test="Chi")

# Retirando 'Zoanthus' chegamos ao mod2f:
f2f <- formula(Chiton_F ~ Crabs + Urchins + Fissu + (1|fSite))
mod2f <- glmmTMB(f2f,family = poisson, data=data2)
summary(mod2f)
drop1(mod2f, test="Chi")

# Retirando 'Zoanthus' chegamos ao mod2g:
f2g <- formula(Chiton_F ~ Crabs + Fissu + (1|fSite))
mod2g <- glmmTMB(f2g,family = poisson, data=data2)
summary(mod2g)

model.sel(mod2, mod2a, mod2b, mod2c,mod2d,mod2e,mod2f,mod2g)

# vendo a diferenca entre os modelos pela ANOVA:
anova(mod2, mod2a, mod2b, mod2c,mod2d,mod2e,mod2f,mod2g)

simulationOutput <- simulateResiduals(fittedModel = mod2d)
plot(simulationOutput)

# Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)

simRes <- simulateResiduals(fittedModel = mod2d,
                            n = 1000)
plot(simRes) # bad!

# Apresentou problema no ajuste dos valores observados vs esperados pelo resíduos simulados
# o Komogorov-Sminoff test falhou em aceitar a hipotese nulo (KS tes signifiativo, p=00.3)

# Ir para distribuição Binomial Negativa do mod2d:

f2d <- formula(Chiton_F ~ Bryo + Crabs + Urchins + Zoanthus  + Fissu + (1|fSite))
mod2dNB <- glmmTMB(f2d,family = nbinom1, data=data2)
summary(mod2dNB)

simulationOutput <- simulateResiduals(fittedModel = mod2dNB)
plot(simulationOutput)

simRes <- simulateResiduals(fittedModel = mod2dNB,
                            n = 1000)
plot(simRes) # Good!

# Modelo Ajustado!!

# Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)


## 3° Formula ####
# Considerando a abundância das espécies que co-ocorrem na lateral da rocha:

f3 <- formula(Chiton_F ~ Bryo + Asci + Crabs + Worms + Spong +
                Urchins + Zoanthus + Hermit + Fissu + (1|fSite))

# Modelo completo:

mod3 <- glmmTMB(f3,family = poisson, data=data)
summary(mod3)

drop1(mod3, test="Chi")

# Retirando 'Hermit' chegamos ao mod3a:
f3a <- formula(Chiton_F ~ Bryo + Asci + Crabs + Worms + Spong +
                 Urchins + Zoanthus + Fissu + (1|fSite))
mod3a <- glmmTMB(f3a,family = poisson, data=data)
summary(mod3a)
drop1(mod3a, test="Chi")

# Retirando 'Worms' chegamos ao mod3b
f3b <- formula(Chiton_F ~ Bryo + Asci + Crabs + Spong +
               Urchins + Zoanthus + Fissu + (1|fSite))
mod3b <- glmmTMB(f3b,family = poisson, data=data)
summary(mod3b)
drop1(mod3b, test="Chi")

# Retirando 'Spong' chegamos ao mod3c:
f3c <- formula(Chiton_F ~ Bryo + Asci + Crabs +
                 Urchins + Zoanthus + Fissu + (1|fSite))
mod3c <- glmmTMB(f3c,family = poisson, data=data)
summary(mod3c)
drop1(mod3c, test="Chi")

# Retirando 'Fissu' chegamos ao mod3d:
f3d <- formula(Chiton_F ~ Bryo + Asci + Crabs +
                 Urchins + Zoanthus + (1|fSite))
mod3d <- glmmTMB(f3d,family = poisson, data=data)
summary(mod3d)
drop1(mod3d, test="Chi")

# Retirando 'Zoanthus' chegamos ao mod3e:
f3e <- formula(Chiton_F ~ Bryo + Asci + Crabs +
                 Urchins + (1|fSite))
mod3e <- glmmTMB(f3e,family = poisson, data=data)
summary(mod3e)
drop1(mod3e, test="Chi")

model.sel(mod3, mod3a, mod3b, mod3c,mod3d,mod3e)

simulationOutput <- simulateResiduals(fittedModel = mod3d)
plot(simulationOutput)

# Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)

simRes <- simulateResiduals(fittedModel = mod3d,
                            n = 1000)
plot(simRes) # bad!

# Apresentou problema no ajuste dos valores observados vs esperados pelo resíduos simulados
# o Komogorov-Sminoff test falhou em aceitar a hipotese nulo (KS tes signifiativo, p=00.3)

# Ajustar o mod3d para Binomial Negativa:

f3d <- formula(Chiton_F ~ Bryo + Asci + Crabs +
                 Urchins + Zoanthus + (1|fSite))
mod3dNB <- glmmTMB(f3d,family = nbinom1, data=data)
summary(mod3dNB)

simulationOutput <- simulateResiduals(fittedModel = mod3dNB)
plot(simulationOutput)

# Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)

simRes <- simulateResiduals(fittedModel = mod3dNB,
                            n = 1000)
plot(simRes) # bad!

# Modelo não ajustado


## Modelo Lateral ####

# Vou criar juntando as variaveis eligidas por cada modelo:

f4 <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo + Bryo + Asci +
                Crabs + Urchins + Zoanthus + Fissu + (1|fSite))

# Modelo completo:

mod4 <- glmmTMB(f4,family = poisson, data=data) # Estou considerando abundancias
summary(mod4)
drop1(mod4, test="Chi")

# Tirando 'Asci' chegamos ao mod4a:

f4a <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo + Bryo +
                Crabs + Urchins + Zoanthus + Fissu + (1|fSite))
mod4a <- glmmTMB(f4a,family = poisson, data=data) # Estou considerando abundancias
summary(mod4a)
drop1(mod4a, test="Chi")

# Tirando 'Bryo' chegamos ao mod4b:

f4b <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo  +
                 Crabs + Urchins + Zoanthus + Fissu + (1|fSite))
mod4b <- glmmTMB(f4b,family = poisson, data=data) # Estou considerando abundancias
summary(mod4b)
drop1(mod4b, test="Chi")

# Tirando 'Zoanthus' chegamos ao mod4c:

f4c <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo  +
                 Crabs + Urchins + Fissu + (1|fSite))
mod4c <- glmmTMB(f4c,family = poisson, data=data) # Estou considerando abundancias
summary(mod4c)
drop1(mod4c, test="Chi")

# Tirando 'Fissu' chegamos ao mod4d:

f4d <- formula(Chiton_F ~ Area.lateral + Cov.Flu + Cov.Bryo  +
                 Crabs + Urchins  + (1|fSite))
mod4d <- glmmTMB(f4d,family = poisson, data=data) # Estou considerando abundancias
summary(mod4d)
drop1(mod4d, test="Chi")

# Vamos comparar os modelos:

library(MuMIn)
model.sel(mod4, mod4a, mod4b, mod4c, mod4d)

# vendo a diferenca entre os modelos pela ANOVA:
anova(mod4, mod4a, mod4b, mod4c, mod4d)

# Podemos considerar tanto o mod1c cquanto o mod1b, mas vamos coniderar o mod1c por ter o maior peso na explicação os modelos e ser o mais parcimonioso

r.squaredGLMM(mod4d)

#Vamos analisar o vif do modelo:
# Calculo da inflacao da variancia para fatores lineares vif - fator de inflacao - menor que 3

library(car)
vif(mod4d)


# Validacao:

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod4d)
plot(simulationOutput)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod4d,
                            n = 1000)
plot(simRes) # good

# Vou odar o mod4d co nbnomial:
mod4dNB1 <- glmmTMB(f4d,family = nbinom1, data=data) # Estou considerando abundancias
mod4dNB2 <- glmmTMB(f4d,family = nbinom2, data=data) # Estou considerando abundancias
summary(mod4dNB1)
summary(mod4dNB2)

simulationOutput <- simulateResiduals(fittedModel = mod4dNB2)
plot(simulationOutput)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod4dNB2,
                            n = 1000)
plot(simRes) # good

# one of the possible test, for other options see ?testResiduals / vignette
testDispersion(simRes)
testZeroInflation(simulationOutput)
testResiduals(simRes)
testUniformity(simulationOutput = simRes)
# b. Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,2))
plotResiduals(simulationOutput, data$Cov.Flu)
plotResiduals(simulationOutput, data$Cov.Bryo)

# Inflação de zero

mod4dZi <- glmmTMB(f4d,family = poisson, zi=~Samp.time, data=data)
summary(mod4dZi)

simulationOutput <- simulateResiduals(fittedModel = mod4dZi)
plot(simulationOutput)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod4dZi,
                            n = 1000)
plot(simRes) # good



+++++++++
# Modelo disconsiderando a inflacao de zero:
library(glmmTMB)

mod0 <- glmmTMB(Chiton_F ~ Temp + Wt.level + Weight + Area +
                  Cov.Flu + Cov.DeathCCAPey + Cov.Asc + Cov.Bryo + Cov.Spong +
                  F.sand + C.sand + Granule + Pebbles + f.Cobbles + Boulders + Rug.base +
                  Bryo + Asci + Bival + Crabs + Shrimp + Barna + Worms + Spong +
                  Urchins + Brittle + Coral + Zoanthus + Flatworm + Hermit + Gastro + (1|f.Site) + (1|Samp.time),
                family = nbinom1(link = "log"), zi=~0, data = data)

summary(mod0)

# Reduzindo modelo 0:

drop1(mod0, test="Chi")

mod0a <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.DeathCCAPey + Cov.Flu + Cov.Bryo + Wt.level + (1|fSite) + (1|Samp.time),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod0a)







# Modelo c/ zero inflado:

mod1 <- glmmTMB(Chiton_F ~ Temp + Wt.level + Weight + Area +
                  Cov.Flu + Cov.DeathCCAPey + Cov.Asc + Cov.Bryo + Cov.Spong +
                  F.sand + C.sand + Granule + Pebbles + Cobbles + Boulders + Rug.base +
                  Bryo + Asci + Bival + Crabs + Shrimp + Barna + Worms + Spong +
                  Urchins + Brittle + Coral + Zoanthus + Flatworm + Hermit + Gastro + (1|f.Site) + (1|Samp.time),
                family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1)

mod1 <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.DeathCCAPey + Cov.Flu + Cov.Asc + Cov.Bryo + Wt.level + (1|fSite) + (1|Samp.time),
                family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1)


# Modeloc/ zero inflado provocado pelas covariaveis:

mod2 <- glmmTMB(Chiton_F ~  Expo.area + Cov.Flu + Cov.Asc + Cov.Bryo + Wt.level + (1|fSite),
                family = nbinom1(link = "log"), zi=~Expo.area+Cov.Flu, data = data)
summary(mod2)


# Reduzindo Modelo 1:

drop1(mod1, test="Chi") # Retirando Chiton_NF

mod1a <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.DeathCCAPey + Cov.Flu + Cov.Bryo + Wt.level + (1|fSite) + (1|Samp.time),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1a)
drop1(mod1a, test="Chi") # Retirando Cov.Asc

mod1b <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.Flu + Cov.Bryo + Wt.level + (1|fSite) + (1|Samp.time),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1b)
drop1(mod1b, test="Chi") # Retirando Wt.level

mod1c <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.Flu + Cov.Bryo + Wt.level + (1|fSite) + (1|Samp.time),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1c)
drop1(mod1c, test="Chi") # Retirando Wt.level

plot(Chiton_F~Expo.area,data=mod1c)

mod1d <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.Flu + Cov.Bryo + (1|fSite) + (1|Samp.time),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1d)
drop1(mod1d, test="Chi")

mod1e <- glmmTMB(Chiton_F ~ Temp + Expo.area + Cov.Flu + Cov.Bryo + (1|fSite) + (1|Samp.time),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod1d)
drop1(mod1d, test="Chi")



# Comparando os modelos:
library(MuMIn)

model.sel(mod1, mod1a, mod1b, mod1c)

r.squaredGLMM(mod1c)

# vendo a diferenca entre os modelos pela ANOVA:
anova(mod1, mod1a, mod1b, mod1c)


install.packages("AICcmodavg")
library(AICcmodavg)

#setup a subset of models of Table 1
Cand.models <- list(mod1, mod1a, mod1b, mod1c)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)




#Vamos analisar o vif do modelo? Calculo da inflacao da variancia para fatores lineares
# vif - fator de inflacao - menor que 3

library(car)
vif(mod1c)

# todos deram menores que 3. Isso e otimo!

# Validacao:
install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod0)
plot(simulationOutput)


# Teste para inflacao de zero
testZeroInflation(simulationOutput)

# b. Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,2))
plotResiduals(simulationOutput, data$Cov.Flu)
plotResiduals(simulationOutput, data$Expo.area)


# Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod1c,
                            n = 1000)
plot(simRes) # good


# Grafico:

install.packages("ggplot2")
install.packages("rlang")
install.packages("ggeffects")
install.packages("effects")
library(ggplot2)
library(ggeffects)

pred_mod1c <- ggpredict(model=mod1c, terms=c("Expo.area","fSite"), type="random")
pred_mod1f <- ggpredict(model=mod1c, terms=c("Cov.Bryo"), type="fixed")
plot(pred_mod1f)

windows(12,6)

plot(pred_mod1c, add.data=T, show.title=T,facets = T, ci=T) +
  theme_classic() + ylab("N° Chitons") +  xlab("Expo.area") +
  labs(colour="") + guides(colour=guide_legend(ncol=4)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.text=element_text(size=12),axis.title=element_text(size=14),
        legend.position="bottom",legend.title=element_text(size=13),
        legend.text=element_text(size=12), strip.text=element_text(size=rel(1.3)))


#install.packages("glmtoolbox")
library(glmtoolbox)

adjR2(mod1c, digits = 4, verbose = TRUE)

install.packages("HH")
library(HH)
vif(data)

head(data[,c(2,8:13)])
vif(data[,c(2:13)])
vif(data[,c(2:13)],y.name="Cov.Others")


# MODELO 3:
mod3 <- glmmTMB(Chiton_F ~ Cov.DeathCCAPey + Expo.area + Cov.Others+ Cov.Flu + Cov.Asc + Cov.Bryo + Wt.level + (1|fSite),
                family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod3)
drop1(mod3, test="Chi") # Retirando Cov.Asc

mod3a <- glmmTMB(Chiton_F ~ Cov.DeathCCAPey + Expo.area + Cov.Others+ Cov.Flu + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod3a)
drop1(mod3a, test="Chi") # Retirando Cov.DeathCCAPEY

mod3b <- glmmTMB(Chiton_F ~  Expo.area + Cov.Others+ Cov.Flu + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod3b)
drop1(mod3b, test="Chi") # Retirando Cov. Others

mod3c <- glmmTMB(Chiton_F ~  Expo.area + Cov.Flu + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod3c)
drop1(mod3c, test="Chi") # Retirando W.Levl

mod3d <- glmmTMB(Chiton_F ~ Expo.area + Cov.Flu + Cov.Bryo + Expo.area*Cov.Bryo + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod3d)

drop1(mod3d, test="Chi") # Retirando W.Levl

model.sel(mod3, mod3a, mod3b, mod3c, mod3d)

# Modelo 4 (w/ Weight) ####

# MODELO 3:

mod4 <- glmmTMB(Chiton_F ~ Weight + Expo.area + Cov.Flu + Cov.DeathCCAPey + Cov.Others + Cov.Asc + Cov.Bryo + Wt.level + (1|fSite),
                family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4)
drop1(mod4, test="Chi") # Retirando Cov.Asc


mod4a <- glmmTMB(Chiton_F ~ Weight + Expo.area + Cov.Flu + Cov.DeathCCAPey + Cov.Others + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4a)
drop1(mod4a, test="Chi") # Retirando Cov.others

mod4b <- glmmTMB(Chiton_F ~ Weight + Expo.area + Cov.Flu + Cov.DeathCCAPey + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4b)
drop1(mod4b, test="Chi")

mod4c <- glmmTMB(Chiton_F ~ Weight + Expo.area + Cov.Flu + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4c)
drop1(mod4c, test="Chi")

mod4d <- glmmTMB(Chiton_F ~ Weight + Expo.area + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4d)
drop1(mod4d, test="Chi")

mod4e <- glmmTMB(Chiton_F ~ Weight +  + Cov.Bryo + Wt.level + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4e)
drop1(mod4e, test="Chi")

mod4f <- glmmTMB(Chiton_F ~ Weight +  + Cov.Bryo  + (1|fSite),
                 family = nbinom1(link = "log"), zi=~1, data = data)

summary(mod4f)
drop1(mod4f, test="Chi")


model.sel(mod4, mod4a, mod4b,mod4c, mod4d, mod4e,mod4f)
anova(mod4, mod4a, mod4b,mod4c, mod4d, mod4e,mod4f)




simulationOutput <- simulateResiduals(fittedModel = mod4b)
plot(simulationOutput)

simRes <- simulateResiduals(fittedModel = mod4b,
                            n = 1000)
plot(simRes) # good

pred_mod4b <- ggpredict(model=mod4b, terms=c("Weight","fSite"), type="random")
pred_mod4b <- ggpredict(model=mod4b, terms=c("Weight"), type="fixed")
plot(pred_mod4b)

windows(12,6)

plot(pred_mod4b, add.data=T, show.title=T,facets = F, ci=T)



testData = createData(sampleSize = 200, overdispersion = 1.5, family = poisson())
fittedModel <- glm(observedResponse ~  Environment1 , family = "poisson", data = testData)

simulationOutput <- simulateResiduals(fittedModel = fittedModel)
plot(simulationOutput)


# Exemplo dharma

testData = createData(sampleSize = 100, overdispersion = 0.5, randomEffectVariance = 0)
fittedModel <- glm(observedResponse ~ Environment1 , family = "poisson", data = testData)
simulationOutput <- simulateResiduals(fittedModel = fittedModel)

# the plot function runs 4 tests
# i) KS test i) Dispersion test iii) Outlier test iv) quantile test
plot(simulationOutput, quantreg = TRUE)

# testResiduals tests distribution, dispersion and outliers
# testResiduals(simulationOutput)

####### Individual tests #######

# KS test for correct distribution of residuals
testUniformity(simulationOutput)

# KS test for correct distribution within and between groups
testCategorical(simulationOutput, testData$group)

# Dispersion test - for details see ?testDispersion
testDispersion(simulationOutput) # tests under and overdispersion

# Outlier test (number of observations outside simulation envelope)
# Use type = "boostrap" for exact values, see ?testOutliers
testOutliers(simulationOutput, type = "binomial")

# testing zero inflation
testZeroInflation(simulationOutput)

# testing generic summaries
countOnes <- function(x) sum(x == 1)  # testing for number of 1s
testGeneric(simulationOutput, summary = countOnes) # 1-inflation
testGeneric(simulationOutput, summary = countOnes, alternative = "less") # 1-deficit

means <- function(x) mean(x) # testing if mean prediction fits
testGeneric(simulationOutput, summary = means)

spread <- function(x) sd(x) # testing if mean sd fits
testGeneric(simulationOutput, summary = spread)

library(R.utils)
install.packages(c("usethis", "renv", "tidyverse"))

# Se você quiser reproduzir os exemplos
install.packages("R.utils")

R.utils::getRelativePath(getwd())


# CHITON CHOICE EXPERIMENT ####

# database: chiton.experiment.csv
# "choice": experimental choice made (1) and not made (O)
# "subst": fluorescence substrate (1) and no fluorescence substrate (2)
# "treat": treatments, "Daylight" (1) and "Dawn light" (2)
# "orient": experimental arena orientation, North (1), South (2), East (3) and West (4)
# "time": time taken to make the choice (Unit: minutes)

data <- read.csv("chiton.experiment.csv",sep = ";",dec=".", header=T)

head(data)
dim(data)
class(data)
str(data)

summary(data)

# Successful choice:
sum(data$choice) # 29 successful

# Unsuccessful choice:
length(data$choice[data$choice!=1]) # 51 unsuccessful

# Percentage of successful trials:
trials <- length(data$choice)
successful <- length(data$choice[data$choice!=0])
percentage <- (successful*100)/trials
percentage # 36.25% successful trials

# Choice/ treatment:
tapply(data$choice,data$treat,sum, na.rm=T) # Daylight (8) & Dawn (21)

# Choice/ orientation:
tapply(data$choice,data$orient,sum,
       na.rm=T) # North (7), South (5), East (6) and West (11)

# Fluorescence substrate choice:
length(data$subst[data$subst=="FLUOR"& !is.na(data$subst)]) # 19

# Non-fluorescence substrate choice:
length(data$subst[data$subst=="UNFLU" & !is.na(data$subst)]) # 10

# Fluorescence choice/ treatment:
tapply(data$subst[data$subst=="FLUOR" & !is.na(data$subst)],
       data$treat[data$subst=="FLUOR" & !is.na(data$subst)],
       length) # Daylight (7) & Dawn (12)

# Fluorescence substrate / Orientation:
tapply(data$subst[data$subst=="FLUOR" & !is.na(data$subst)],
       data$orient[data$subst=="FLUOR" & !is.na(data$subst)],
       length)

# Non-fluorescence substrate / Orientation:
tapply(data$subst[data$subst=="UNFLU" & !is.na(data$subst)],
       data$orient[data$subst=="UNFLU" & !is.na(data$subst)],
       length)


# Spend time to choice/ treatment:
tapply(data$time,data$treat,mean,na.rm=T) # Daylight (6.67) & Dawn (7.92)
tapply(data$time,data$treat,sd,na.rm=T)
boxplot(time ~ treat,varwidth = F, data=data)

# Spend time / substrate choice:
tapply(data$time,data$subst,mean,na.rm=T) # Daylight (8.24) & Dawn (6.32)
tapply(data$time,data$subst,sd,na.rm=T)
boxplot(time ~ subst,varwidth = F, data=data)

library(lattice)
xyplot(time ~ factor(subst) | treat, data=data, type=c("p","r"))

# QQ-plot for 'time':
par(mfrow= c (1,1))
qqnorm(data$time)
qqline(data$time, col="blue", lwd=2)

# Normality test:
shapiro.test(data$time) # H0 rejected (non-normal)

# Homogeneity test:
variance <- tapply(data$time,data$subst,var)
variance

variance [1]/ variance [2] # 1.07

# Bartlett' test:
bartlett.test(data$time ~ data$subst) # H0 accepted (homoscedasticity)

# Fligner-Killeen' test:
fligner.test(data$time ~ data$subst) # H0 accepted (homoscedasticity)

# Wilcox.test:

wilcox.test(data$time[data$treat=="Daylight"],
            data$time[data$treat=="Dawn"],
            conf.int=TRUE,
            conf.level=0.95) # H0 accepted

wilcox.test(data$time[data$subst=="FLUOR"],
            data$time[data$subst=="UNFLU"],
            conf.int=TRUE,
            conf.level=0.95) # H0 accepted


## Treatment x Orientation ####

vector1 <- c(2,2,1,3,5,3,5,8)
orient <- matrix(vector1,byrow=T,nrow = 2)
colnames(orient)<-c("North", "South", "East", "West")
rownames (orient) <- c("Daylight","Dawn")
orient

chisq.test(orient)
chisq.test(orient)$expected

fisher.test(orient)

# Substrate choice x Orientation
vector2 <- c(3,5,3,8,3,2,2,3)
orient2 <- matrix(vector2,byrow=T,nrow = 2)
colnames(orient2)<-c("East","North","South","West")
rownames (orient2) <- c("Fluorescence","Non-Fluorescence")
orient2

chisq.test(orient2)
chisq.test(orient2)$expected

fisher.test(orient2)

## 2. Substrate choice x Treatment

vector3 <- c(7,12,1,9)
choices <- matrix(vector3,byrow=T,nrow = 2)
colnames(choices)<-c("Daylight","Dawn")
rownames (choices) <- c("Fluorecence","Non-Fluorescence")
choices

chisq.test(choices, correct=F)
chisq.test(choices)$expected

fisher.test(choices)

install.packages("Exact")
library(Exact)

exact.test(choices, method="z-pooled", model="Multinomial")

# "choice" FAZER

B1 <- glm(choice~subst*treat+time, data=data,
          family = binomial(link="cloglog"), na.action = na.omit)

summary(B1)
drop1(B1, tst="Chi")

B2 <- glm(choice~subst+treat+time, data=data,
          family = binomial(link="cloglog"), na.action = na.omit)

summary(B2)
drop1(B2, tst="Chi")

B3 <- glm(choice~subst+treat, data=data,
          family = binomial(link="cloglog"), na.action = na.omit)

summary(B3)

B4 <- glm(choice~subst, data=data,
          family = binomial, na.action = na.omit)

summary(B4)

B5 <- glm(choice~treat, data=data,
          family = binomial(link="cloglog"), na.action = na.omit)

summary(B5)

library(MuMIn)
model.sel(B1,B2,B3,B4,B5)
anova(G1,G2,G3,G3a)


# "time"

kruskal.test(time ~ subst,
             data = data)



G1 <- glm(time~fsubst*ftreat, data=data, family=gaussian)
summary(G1a)
drop1(G1,test="Chi")

G2 <- glm(time~fsubst+ftreat, data = data, family = gaussian)
summary(G2)
drop1(G2,test="Chi")

G3 <- glm(time~subst, data = data, family = gaussian)
summary(G3)

data$subst2 <- relevel(data$fsubst, ref="UNFLU")
levels(data$subst2)

G3a <- glm(time~subst2, data = data, family = gaussian)
summary(G3a)

library(MuMIn)
model.sel(G1,G1a,G2,G3,G3a)
anova(G1,G2,G3,G3a)

summary(G1)

GM1 <- glm(time~subst*treat, data=data, family=Gamma, na.action = na.omit)
summary(GM1)
drop1(GM1,test="Chi")

GM2 <- glm(time~fsubst+ftreat, data=data, family=Gamma)
summary(GM2)
drop1(GM2,test="Chi")

GM3 <- glm(time~subst, data=data, family=Gamma)
summary(GM3)
drop1(GM3,test="Chi")

(11.588/25) # 0.46

library(car)
vif(GM1)

library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = GM1)
plot(simulationOutput)

#Conclusions: there is no obvious patterns in the qq-plot or residuals, or at least there are no obvious trends remaining that would be indicative of overdispersion or non-linearity.

# Exploring the goodness of the fit of the model

testUniformity(GM1)

# Conclusions: neither Pearson residuals, Deviance or Kolmogorov-Smirnov test of uniformity indicate a lack of fit (p values greater than 0.05)

AIC(G1, GM1)

#Poisson deviance
G1$deviance
[1] 71.62288
#Gaussian deviance
GM1$deviance


testOverdispersion(simulateResiduals(GM1, refit=T))

# Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)

simRes <- simulateResiduals(fittedModel = GM1,
                            n = 1000)
plot(simRes) # bad!



par(mfrow=c(2,2))
plot(G3)

P1 <- glm(time~fsubst*ftreat, data=data, family=poisson)
summary(P1)

# deviance explained:
(((82.686-71.623)/82.686)*100) # 13%

# overdispersion?

(71.623/25) # 2.86492

library(glmmTMB)
mod1 <- glmmTMB(time~fsubst*ftreat, data=data, family=nbinom1)
summary(mod1)

mod2 <- glmmTMB(time~fsubst+ftreat, data=data, family=nbinom1)
summary(mod2)

mod2 <- glmmTMB(time~fsubst, data=data, family=nbinom1)
summary(mod2)

library(car)
vif(mod1)



Q1 <- glm(time~fsubst*ftreat, data=data, family=quasipoisson, na.action = na.omit)
summary(Q1)

Q2 <- glm(time~subst+ftreat, data=data, family=quasipoisson, na.action = na.omit)
summary(Q2)

drop1(Q2, test="Chi")

Q3 <- update(Q2,~.-treat)
summary(Q3)

library(MuMIn)
model.sel(Q1, Q2, Q3)
anova(Q1, Q2, Q3)

par(mfrow=c(2,2))
plot(mod1)

preditos <- predict(Q1, type="response")
RordQ <- data$time - preditos
RpeS <- RordQ / sqrt(7.630148 * preditos)
par(mfrow=c(1,1))
plot(x=preditos, y=RpeS, main="Pearson residuals scaled", cex=2)
abline(h=0, lty=2)


library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = G1)
plot(simulationOutput)

# b. Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,2))
plotResiduals(simulationOutput, data$time)
plotResiduals(simulationOutput, data$fsubst)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = P2,
                            n = 1000)
plot(simRes) # good




par(mfrow= c (2,2))
plot(P2q)




library(MASS)
NB1 <- glm.nb(time~fsubst*ftreat, data=data, link="log", na.action = na.omit)
summary(NB1)


NB1 <- glm.nb(time~fsubst+ftreat, data=data, link="log", na.action = na.omit)
summary(NB1)

data$subst2 <- relevel(designANCOVA$sex, ref="male")
levels(designANCOVA$sex2)

B1 <- glm(fsubst~time+ftreat+forient, family=binomial, data=data)
summary(B1)


