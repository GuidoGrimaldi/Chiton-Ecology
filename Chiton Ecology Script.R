#++++++++++++++++++++++++#
# Doutorado - UFSC
# Chapter 3: Ecology of chitons
#
# Script for: "chitons.eco.data"
#++++++++++++++++++++++++#

# Setting ####
R.utils::getRelativePath(getwd())
usethis::create_project("~/Documentos/R sessions/Chiton Ecology")

usethis::use_git() # 2 "Yeah" / 3 "Yes"


gitcreds::gitcreds_set()

usethis::use_github()


setwd("C:\\Users\\guido\\Documentos\\Git\\2023_Ecology_chitons")
getwd()
dir()

# Packages ####
install.packages("MuMIn")
library(DHARMa) # To validate model
library(glmmTMB) # To fit GLMM model zero-inflated
library(car) # For leveneTest and vif model
library(lattice)
library(AICcmodavg)
library(beeswarm) # For exploration graphs
library(MuMIn)

install.packages("ggplot2")
install.packages("rlang")
install.packages("ggeffects")
install.packages("effects")
library(ggplot2)
library(ggeffects)


#install.packages("glmtoolbox")
library(glmtoolbox)

adjR2(mod1c, digits = 4, verbose = TRUE)

install.packages("HH")
library(HH)
vif(data)

# Data ####

data <- read.csv("chitons.eco.data.csv", sep=";",dec=".", header=T)
head(data)
dim(data)
str(data)
data$fSite <- factor(data$Site)
dim(data)

# Variables ####

# a. Variable Y:
# + Chiton_F: Abundancia de quitons fluorescentes (Ischonplax pectinata)

# b. Variable X:
# + "Site" -> Locais de coleta (Fator)
# + "Samp.Time" -> Tempo de amostragem (min)
# + "Temp" -> Temperatura da agua (°C) - Var. discreta
# + "Wt.level" -> Altura da mare em relacao a base da rocha (cm) - Var. discreta
# + "Weight" -> Peso (Kg) - Var. continua
# + "Exposed.area -> Area superficial lateral exposta (cm2) - Var. continua
# + "Cov.Flu" -> % cobertura alga CCA + Peyssonnelia
# + "Cov.DeathCCAPey" -> % cobertura de alga CCA e Pey morta.
# + "Cov.Others" -> % cobertura de outros substratos (rocha nua, "mucos", other)
# + "Cov.Asci" -> % cobertura de Ascidias
# + "Cov.Bryo" -> % cobertura de Briozoarios
# + "Cov.Spon" -> % cobertura de Esponjas.

# Data Exploring ####

# 1. Amostragem balanceadas ?
tapply(data$Chiton_F,data$fSite,length) # YES

# 2. Missing data ?
summary(data) # NO

# 3. Zero trouble Y ?
table(data$Chiton_F)
plot(table(data$Chiton_F))
barplot(table(data$Chiton_F)) # Bastante zeros!

# Podemos saber a proporcao de zeros aplicando o seguinte calculo:
sum(data$Chiton_F==0)/length(data$Chiton_F)*100

# Proporcao dos dados que sao iguais a zero eh 63.33%

# Vamos olhar em um histograma com a distribuicao de probabilidade
hist(data$Chiton_F, prob = T, breaks=15)
rug(data$Chiton_F)
lines(density(data$Chiton_F), col="blue")

# Os dados apresentam uma distribuicao unimodal nao normal, ja esperada para dados discretos.
# Mas  vamos examinar...

# 4. Normality Y ?

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

# 5. Outliers Y & X ?

# A ferramenta grafica tipica utilizada para deteccao de outliers eh o boxplot.

# Vamos analisar por partes:
# Quitons:
boxplot(data[,c(2,3)])

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

# Modelo ####

library(glmmTMB)

# Modelo disconsiderando a inflacao de zero:

mod0 <- glmmTMB(Chiton_F ~ Samp.time + Temp + Expo.area + Cov.DeathCCAPey + Cov.Flu + Cov.Asc + Cov.Bryo + Wt.level + (1|fSite),
                family = nbinom1(link = "log"), zi=~0, data = data)

summary(mod0)


# Modelo c/ zero inflado:

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
