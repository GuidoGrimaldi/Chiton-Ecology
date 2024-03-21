#++++++++++++++++++++++++#
# Doutorado - UFSC
# Chapter 3: Ecology of chitons
#
# Script for: "chitons.eco.data"
#++++++++++++++++++++++++#

#install.packages("usethis")
usethis::create_project("~/Documents/R sessions/Chiton Ecology/")

# GIT & GITHUB ####

usethis::use_git() # 2 "Yeah" / 3 "Yes"

# gitcreds::gitcreds_set() # Cole o Token. Basta uma unica vez

usethis::use_github()

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

macro <- read.csv("rockside.csv", sep=";",dec=".", header=T)
str(macro)
macro$fSite <- factor(macro$Site)


# Data Sets Variables ####

## Data set: "data"
### A. Response Variable:
# + "Chiton_F": abundance of fluorescent chitons (Ischonoplax pectinata)

### B. Explanatory Variable:
# + "Site": study sites
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

## Data set:"macro"

# + Asci:
# + Barnacle:
# + Bival:
# + Brittle:
# + Bryo:
# + Cerithium:
# + Coral:
# + Crabs:
# + C.merca:
# + Fissu:
# + Flatworms:
# + Gastrop:
# + Hermit:
# + L.nassa:
# + Not.ident:
# + N.virgi
# + P.pusio:
# + Shrimp:
# + Spong:
# + T.viri:
# + Urchins:
# + Worms:
# + Zoanthus:

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

# De fato, Buzio eh diferente das demais praias


# ~Sampling Time ####
tapply(data$Samp.time,data$Site,sum)
tapply(data$Samp.time,data$Site,mean)
tapply(data$Samp.time,data$Site,sd)

# Normality test:
shapiro.test(data$Samp.time) # Nao  normal

# Homogeneity test:
variancia <- tapply(data$Samp.time,data$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [3]  # Buzios eh 2x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(data$Samp.time ~ data$fSite) # H0 aceita (Homocedastica)

# Teste de Fligner-Killeen:
fligner.test(data$Samp.time ~ data$Site) # H0 aceita (Homocedastica)

library(car)
leveneTest(Samp.time ~ fSite, data=data, center=mean) # H0 aceita (Homocedastica)

# Kruaskal-Wallis test:

kruskal.test(Samp.time ~ fSite,
             data = data) # H0 rejeitada, ha diferenca entre praias

# Dunn test:
library(FSA)
DT = dunnTest(Samp.time ~ fSite,
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

boxplot(Samp.time~fSite, data = data)

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

boxplot(Wt.level~fSite, data = data)

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


# ~ Boulder Area ####
summary(macro$Area.total)
mean(macro$Area.total)
sd(macro$Area.total)

tapply(macro$Area.total,macro$fSite,mean)
tapply(macro$Area.total,macro$fSite,sd)

# Normality test:
shapiro.test(macro$Area.total) # Nao  normal

# Homogeneity test:
variancia <- tapply(macro$Area.total,macro$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [1]  # Buzios eh 3x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(macro$Area.total ~ macro$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(macro$Area.total ~ macro$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Area.total ~ fSite, data=macro, center=mean) # H0 aceita (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Area.total ~ Site,
             data = macro) # H0 aceita, nao ha diferenca entre praias


# ~ Lateral Exposed area ####
summary(macro$Area.lateral)
mean(macro$Area.lateral)
sd(macro$Area.lateral)

tapply(macro$Area.lateral,macro$Site,mean)
tapply(macro$Area.lateral,macro$fSite,sd)

# Normality test:
shapiro.test(macro$Area.lateral) # Nao  normal

# Homogeneity test:
variancia <- tapply(macro$Area.lateral,macro$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [2]/ variancia [4]

# Teste de Bartlett:
bartlett.test(macro$Area.lateral ~ macro$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(macro$Area.lateral ~ macro$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Area.lateral ~ fSite, data=macro, center=mean) # H0 aceita (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Area.lateral ~ Site,
             data = macro) # H0 aceita, nao ha diferenca entre praias

# ~Cov.Flu ####
mean(macro$Cov.Flu)
sd(macro$Cov.Flu)
tapply(macro$Cov.Flu,macro$fSite,mean)
tapply(macro$Cov.Flu,macro$fSite,sd)

# Normality test:
shapiro.test(macro$Cov.Flu) # Nao  normal

# Homogeneity test:
variancia <- tapply(macro$Cov.Flu,macro$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [3]/ variancia [2]  # Pitangui eh 2.35x maior...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(macro$Cov.Flu ~ macro$fSite) # H0 aceita (Homocedastico)

# Teste de Fligner-Killeen:
fligner.test(macro$Cov.Flu ~ macro$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Cov.Flu ~ fSite, data=macro, center=mean) # H0 aceita (Homocedastico)

# Kruaskal-Wallis test:

kruskal.test(Cov.Flu ~ Site,
             data = macro) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(Cov.Flu ~ fSite,
              data=macro,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

boxplot(Cov.Flu ~ Site,
        data = macro)

# ~Cov.UNFlu ####

tapply(data$Cov.DeathCCAPey,data$Site,mean)
tapply(data$Cov.DeathCCAPey,data$Site,sd)

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

boxplot(Cov.Asc ~ Site,
        data = data)

# ~Cov.Bryo ####
tapply(macro$Cov.Bryo,macro$fSite,mean)
tapply(macro$Cov.Bryo,macro$fSite,sd)

# Normality test:
shapiro.test(macro$Cov.Bryo) # Nao  normal

# Homogeneity test:
variancia <- tapply(macro$Cov.Bryo,macro$fSite,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [4]/ variancia [2]  # Santa Rita eh 13x maior que Buzios

# Teste de Bartlett:
bartlett.test(macro$Cov.Bryo ~ macro$fSite) # H0 rejeitada (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(macro$Cov.Bryo ~ macro$fSite) # H0 aceita (Homocedastico)

library(car)
leveneTest(Cov.Bryo ~ fSite, data=macro, center=mean) # H0 Rejeitada (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(Cov.Bryo ~ Site,
             data = macro) # H0 aceita, nao ha diferenca entre praias


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

boxplot(Cov.Spong ~ Site,
           data = data)


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

# ~cov.cca ####
mean(macro$cov.cca)
sd(macro$cov.cca)
tapply(macro$cov.cca,macro$Site,mean)
tapply(macro$cov.cca,macro$Site,sd)

# Normality test:
shapiro.test(macro$cov.cca) # Normal

# Homogeneity test:
variancia <- tapply(macro$cov.cca,macro$Site,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [1]/ variancia [2]  # Baia Formosa excede 3.35x...provavelmente homocedastico

# Teste de Bartlett:
bartlett.test(macro$cov.cca ~ macro$Site) # H0 Rejected (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(macro$cov.cca ~ macro$Site) # H0 Accepted (Homocedastico)

library(car)
leveneTest(cov.cca ~ Site, data=macro, center=mean) # H0 Rejected (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(cov.cca ~ Site,
             data = macro) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(cov.cca ~ Site,
              data=macro,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

boxplot(cov.cca ~ Site,
        data = macro)


# ~cov.pey ####
mean(macro$cov.pey)
sd(macro$cov.pey)
tapply(macro$cov.pey,macro$Site,mean)
tapply(macro$cov.pey,macro$Site,sd)

# Normality test:
shapiro.test(macro$cov.pey) # Non Normal

# Homogeneity test:
variancia <- tapply(macro$cov.pey,macro$Site,var)
variancia

# Uma regra pratica eh que a maior variancia nao deve exceder de 3 a 4 vezes a
# menor variancia:

variancia [1]/ variancia [4]  # Baia Formosa excede 9x...provavelmente Heterocedastico

# Teste de Bartlett:
bartlett.test(macro$cov.pey ~ macro$Site) # H0 Rejected (Heterocedastico)

# Teste de Fligner-Killeen:
fligner.test(macro$cov.pey ~ macro$Site) # H0 Rejected (Heterocedastico)

library(car)
leveneTest(cov.pey ~ Site, data=macro, center=mean) # H0 Rejected (Heterocedastico)

# Kruaskal-Wallis test:

kruskal.test(cov.pey ~ Site,
             data = macro) # H0 rejeitada, ha diferenca entre praias

# Dunn test:

library(FSA)

DT = dunnTest(cov.pey ~ Site,
              data=macro,
              method="bh")      # Adjusts p-values for multiple comparisons;

DT

# Compact letter display:

PT = DT$res

PT

library(rcompanion)

cldList(P.adj ~ Comparison,
        data = PT,
        threshold = 0.05)

boxplot(cov.pey ~ Site,
        data = macro)


# 7. Relationship X & Y ? ####

# Vamos analisar por partes:

# Var.abioticas primeiro
head(data)
pairs(data[,c(2:8)])

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
pairs(data[,c(2:8)], diag.panel=panel.hist)

# Outra funcao para fornecer os coeficientes de correlacao ("KENDALL")
# com a fonte proporcional ao indice
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

# Agora plotamos com indices de correlacao e histogramas
pairs(data[,c(2:8)], diag.panel=panel.hist, lower.panel=panel.cor)

# Incluimos tambem a linha de suavizacao
pairs(data[,c(2:8)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Temos uma Correlacao Positiva entre Weight X Expo.area (0.54)
# Temos uma Correlacao Positiva entre Weight X Samp.time (0.50)
# Temos uma Correlacao Positiva fraca entre Samp.time X Expo.area (0.38)
# Temos uma Correlacao Negativa fraca entre Temp X Wt.level (-0.32)


cor.test(data$Temp,data$Wt.level, method = "k")

cor.test(data$Samp.time,data$Weight, method = "k")


cor.test(data$Samp.time,data$Chiton_F, method = "k")

# Vamos ver as coberturas laterais...
head(data)
pairs(data[,c(2,8:14)])

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
pairs(data[,c(2,8:14)], diag.panel=panel.hist)

# Outra funcao para fornecer os coeficientes de correlacao
# com a fonte proporcional ao indice
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

# Agora plotamos com indices de correlacao e histogramas
pairs(data[,c(2,8:14)], diag.panel=panel.hist, lower.panel=panel.cor)

# Incluimos tambem a linha de suavizacao
pairs(data[,c(2,8:14)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Correlacao entre: Cov.Flu x Cov.DeathCCAPey (-0.36)
# Correlacao entre: Cov.Flu x Cov.Others (-0.71)
# Correlacao entre: Cov.Flu x Cov.Spong (-0.39)

# Vamos em ver todas:
x11()
pairs(data[,c(3:14)])

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

pairs(data[,c(3:14)], diag.panel=panel.hist)

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

pairs(data[,c(3:14)], diag.panel=panel.hist, lower.panel=panel.cor)

pairs(data[,c(3:14)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Correlacao entre: Temp x Cov.Spong (-0.32)
# Correlacao entre: Wt.level x Cov.Others (-0.44)
# Correlacao entre: Weight x Cov.Spong (0.32)

# Vamos determinar a signifiancia das correlacoes...

# Para isso precisamos testar primero as normalidades

shapiro.test(data$Samp.time) # Nao normal
shapiro.test(data$Temp)   # Nao normal
shapiro.test(data$Wt.level) # Nao normal
shapiro.test(data$Weight)   # Nao normal
shapiro.test(data$Area.total) # Nao normal
shapiro.test(data$Area.lateral) # Nao normal
shapiro.test(data$Cov.Flu) # Nao normal
shapiro.test(data$Cov.DeathCCAPey) # Nao normal
shapiro.test(data$Cov.Others) # Nao normal
shapiro.test(data$Cov.Asc) # Nao normal
shapiro.test(data$Cov.Bryo) # Nao normal
shapiro.test(data$Cov.Spong) # Nao normal

# Kendall correlations:

# 1. Wt. level X Temp (-0.31)
cor.test(data$Temp, data$Wt.level, method= "k") # p<0.001
plot(data$Temp, data$Wt.level)

# 2. Samp.time X Weight (0.27)
cor.test(data$Weight, data$Samp.time, method= "k") # p<0.001
plot(data$Weight, data$Samp.time)

# 3. Area.lateral  X Weight (0.39)
cor.test(data$Weight, data$Area.lateral, method= "k") # p<0.001
plot(data$Weight, data$Area.lateral)

# 4. Area.lateral  X Samp. time (0.20)
cor.test(data$Samp.time, data$Area.lateral, method= "k") # p<0.01
plot(data$Sap.time, data$Area.lateral)

# 5. Area.lateral  X Area.total (0.47)
cor.test(data$Area.total, data$Area.lateral, method= "k") # p<0.001
plot(data$Area.total, data$Area.lateral)

# 6. Cov.Flu x Cov.DeathCCAPey (-0.27)
cor.test(data$Cov.Flu, data$Cov.DeathCCAPey, method= "kendall") # p<0.001
plot(data$Cov.Flu, data$Cov.DeathCCAPey)

# 7. Cov.Bryo x Cov.Asc (0.25)
cor.test(data$Cov.Bryo, data$Cov.Asc, method= "k") # p<0.001
plot(data$Cov.Bryo, data$Cov.Asc)

# 8. Temp x Cov.Spong (-0.20)
cor.test(data$Temp, data$Cov.Spong, method= "k") # p<0.001
plot(data$Temp, data$Cov.Spong)

# 9. Weight x Cov.Spong (0.29)
cor.test(data$Weight, data$Cov.Spong, method= "k") # p<0.001
plot(data$Weight, data$Cov.Spong)

# 10. Samp.time x Cov.Spong (0.22)
cor.test(data$Samp.time, data$Cov.Spong, method= "k")  #p<0.01
plot(data$Samp.time, data$Cov.Spong)

# 11. Cov.Spong X Cov.Asc (0.28)
cor.test(data$Cov.Spong, data$Cov.Asc, method= "k")  #p<0.001
plot(data$Cov.Spong, data$Cov.Asc)

# 12. Wt.level x Cov.Others (-0.35)
cor.test(data$Wt.level, data$Cov.Others, method= "k") #p<0.001
plot(data$Wt.level, data$Cov.Others)

# 13. Cov.Flu x Cov.Others (-0.50)
cor.test(data$Cov.Flu, data$Cov.Other, method= "k") #p<0.001
plot(data$Cov.Flu, data$Cov.Other)

# Vamos ver para a coocorrencia de invertebrados

dim(data)

x11()
pairs(data[,c(15:29)])

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

pairs(data[,c(3:14)], diag.panel=panel.hist)

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

pairs(data[,c(3:14)], diag.panel=panel.hist, lower.panel=panel.cor)

pairs(data[,c(3:14)], panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# ~ Macrofauna ####

# Tranformando a matriz de macrofauna em presenca/ausencia.

dim(macro)
str(macro)
library(vegan)
spe.pa <- decostand(macro[,15:34], method="pa")
str(spe.pa)
dim(spe.pa)
macro.pa <- c(macro [,1:14],spe.pa)
str(macro.pa)
dim(macro.pa)

# ~Bryozoans ####
sum(macro.pa$Bryo) # 37
tapply(macro.pa$Bryo,macro.pa$Site,sum)

# ~Ascidians ####
sum(macro.pa$Asci) # 98
tapply(macro.pa$Asci,macro.pa$Site,sum)

mod.logistico <- glm(Asci~Site, family=binomial(link="logit"), data=macro.pa)	# ajustar o modelo logistico (family=binomial)
summary(mod.logistico)
anova(mod.logistico,test="Chisq")

mod.logistico <- glm(Asci~Weight, family=binomial(link="logit"), data=macro.pa)	# ajustar o modelo logistico (family=binomial)
summary(mod.logistico)
anova(mod.logistico,test="Chisq")

plot(macro.pa$Asci~macro.pa$Weight, pch=16,
     xlab="Reefs", ylab="Presence of ascidians")

probabilidade <- exp(predict(mod.logistico)) / (1+exp(predict(mod.logistico))) # converter os valores preditos no modelo logit para probabilidades
lines(probabilidade~macro.pa$Weight, lwd=2, lty=2)


# ~Bivalves ####
sum(macro.pa$Bival) # 3
tapply(macro.pa$Bival,macro.pa$Site,sum)

# ~Crabs ####
sum(macro.pa$Crabs) # 14 occurence
tapply(macro.pa$Crabs,macro.pa$Site,sum)

sum(macro$Crabs) # 15 abundance
tapply(macro$Crabs,macro$Site,sum)

# ~Polychaetes ####
sum(macro.pa$Worm) # 17 occurence
tapply(macro.pa$Worm,macro.pa$Site,sum)

sum(macro$Worms) # 17 abundance
tapply(macro$Worms,macro$Site,sum)

# ~Sponges ####
sum(macro.pa$Spong) # 65
tapply(macro.pa$Spong,macro.pa$Site,sum)

# ~Zoanthus ####
sum(macro.pa$Zoanthus) # 10
tapply(macro.pa$Zoanthus,macro.pa$Site,sum)

# ~Sponges ####
sum(macro.pa$Spong) # 65
tapply(macro.pa$Spong,macro.pa$Site,sum)

# ~Hermit ####
sum(macro.pa$Hermit) # 9 occurence
tapply(macro.pa$Hermit,macro.pa$Site,sum)

sum(macro$Hermit) # 14 abundance
tapply(macro$Hermit,macro$Site,sum)

# ~Gastropods ####
sum(macro.pa$Gastrop) # 62 occurence
tapply(macro.pa$Gastrop,macro.pa$Site,sum)

sum(macro$Gastrop) # 178 abundance
tapply(macro$Gastrop,macro$Site,sum)

# ~Fissu ####
sum(macro.pa$Fissu) # 38 occurence
tapply(macro.pa$Fissu,macro.pa$Site,sum)

sum(macro$Fissu) # 114 abundance
tapply(macro$Fissu,macro$Site,sum)

# ~T.viri ####
sum(macro.pa$T.viri) # 17 occurence
tapply(macro.pa$T.viri,macro.pa$Site,sum)

sum(macro$T.viri) # 22 abundance
tapply(macro$T.viri,macro$Site,sum)


# ~ Diagnostico ####
# REFAZER

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




# MODELING ####

# ~at side scale ####

## 1° Model ####
# O objetivo eh investigar se a abundancia de quitons fluorescentes ('Chiton_F')
# eh influenciada pelo tipo de substrato que recobre a lateral dos boulders que
# os chitons habitam.

# Hipótese: maior cobertura de substratos fluorescentes (cca+peyssonnelia),
# maiores as são as abundancias dos quítons fluorescentes.

macro <- read.csv("rockside.csv", sep=";",dec=".", header=T)
str(macro)
macro$fSite <- factor(macro$Site)
head(macro)

# Variable Y: Chiton_F
# Variables X: Area.lateral,Cov.Flu,Cov.Asc,cov.unflu,Cov.Bryo,Cov.Spong,Cov.Others)
# Random factor: fSite

## Exploring variables ##

Z<-macro[,c("Area.lateral","Cov.Flu","Cov.Asc","cov.unflu","Cov.Bryo","Cov.Spong","Cov.Others")]

#install.packages("usdm")
library(usdm)
vif(Z)
vifcor(Z)
vifstep(Z, th=3,keep = "Cov.Flu") # cov.unflu excluded

boxplot(macro[,c("Area.lateral","Cov.Flu","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")])

pairs(macro[,c("Chiton_F","Area.lateral","Cov.Flu","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")],
      panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Transformando Area.lateral
macro$logArea.lateral <- log(macro$Area.lateral)
str(macro)

boxplot(macro[,c("logArea.lateral","Cov.Flu","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")])

pairs(macro[,c("Chiton_F","logArea.lateral","Cov.Flu","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")],
      panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)


# A variavel "area.lateral" apresenta dimensoes muito maiores que as demais.

dotplot(as.matrix(Z), groups = FALSE,
        strip = strip.custom(bg = 'white',
                             par.strip.text = list(cex = 0.8)),
        scales = list(x = list(relation = "free"),
                      y = list(relation = "free"),
                      draw = FALSE),
        col = 1, cex  = 0.5, pch = 16,
        xlab = "Value of the variable",
        ylab = "Order of the data")

# Vamos anasila-la em isolado

# Make the coplot
coplot(Chiton_F ~ Cov.Bryo | Site, data=macro, ylab = "Abudance",
       xlab = "Cover Fluoescence (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

library(beeswarm)
beeswarm::beeswarm(Area.lateral ~ Site,data=macro)
beeswarm::beeswarm(Area.lateral ~ Site,data=macro, col="red", pch=16, method="swarm")

dotchart(macro$Area.lateral, main = "AREA", group = macro$fSite)

stripchart(Area.lateral ~ Site,data=macro, method="stack")

library(lattice)
xyplot(Chiton_F ~ Area.lateral | Site, data=macro)
xyplot(Chiton_F ~ Cov.Bryo | Site, data=macro, type=c("p","r"))

xyplot(Area.lateral ~Cov.Flu| Site, data=macro, type=c("p","r"))


### ~GLMM Poisson ####

str(macro)
library(pscl)

f.full <- formula(Chiton_F ~ Area.lateral + Cov.Flu  + Cov.Asc + Cov.Bryo
              + Cov.Spong + (1|fSite))

library(glmmTMB)
mod1 <- glmmTMB(f.full,family = poisson, data=macro)
summary(mod1)

drop1(mod1, test="Chi") # Cov. Asc
f1 <- update(f.full,.~.-Cov.Asc)

mod1a <- glmmTMB(f1, family = poisson, data=macro)
summary(mod1a)

drop1(mod1a, test="Chi") # Cov. Spong
f2 <- update(f1,.~.-Cov.Spong)

mod1b <- glmmTMB(f2, family = poisson, data=macro)
summary(mod1b)

### ~AIC Selection ####

library(MuMIn)
model.sel(mod1, mod1a, mod1b) # mod1b (AIC: 305)
anova(mod1, mod1a, mod1b, test="Chi")

#install.packages("AICcmodavg")
library(AICcmodavg)

#setup a subset of models of Table 1
Cand.models <- list(mod1, mod1a, mod1b)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

r.squaredGLMM(mod1b)

### ~Model Validation ####

#install.packages("DHARMa")
library(DHARMa)

simulationOutput <- simulateResiduals(fittedModel = mod1b,
                                      plot=T)

plot(residuals(simulationOutput)) # Asymptotic is good!
plotResiduals(simulationOutput, form = macro$Cov.Flu)


# a. Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod1b,
                            n = 1000)
plot(simRes) # good

# b Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,3))
plotResiduals(simulationOutput, macro$Cov.Flu)
plotResiduals(simulationOutput, macro$Cov.Bryo)
plotResiduals(simulationOutput, macro$Area.lateral)

# c. Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)

# d. Dispersion
testDispersion(mod1b)
testDispersion(simRes)
testOutliers(simulationOutput = simulationOutput)

# e. spatial autocorrelation # DEU ERRO

testSpatialAutocorrelation(simulationOutput, x=macro$x, macro$y)
testSpatialAutocorrelation(simulationOutput, x=macro$Cov.Flu, macro$Chiton_F)
testSpatialAutocorrelation(mod1b, x=macro$Cov.Flu, macro$Chiton_F)
testSpatialAutocorrelation(mod1b, x=macro$x, macro$y)
testSpatialAutocorrelation(simRes, x=macro$x, macro$y)
testSpatialAutocorrelation(simRes, x=macro$Cov.Flu, macro$Chiton_F)
testSpatialAutocorrelation(simulationOutput, x=macro$Chiton_F, macro$Cov.Flu)

dM = as.matrix(dist(cbind(macro$x, macro$y)))
testSpatialAutocorrelation(simulationOutput, distMat = dM)


install.packages("sfsmisc")
library(sfsmisc)
testData = createData(sampleSize = 100, family = poisson(), spatialAutocorrelation = 5)
head(testData)
summary(testData)
fittedModel <- glmer(observedResponse ~ Environment1 + (1|group), data = testData, family = poisson() )
simulationOutput <- simulateResiduals(fittedModel = fittedModel)
testSpatialAutocorrelation(simulationOutput = simulationOutput, x = testData$x, y= testData$y)




## CONCLUIDO ###

## 2° MODEL ####

# O objetivo eh investigar a variavel "COV.FLU" como  especies distintas para
# saber se ha influencia especie-especifica.

# Hipótese: nao havera influencia especie especifica.

macro <- read.csv("rockside.csv", sep=";",dec=".", header=T)
str(macro)
macro$fSite <- factor(macro$Site)

# Variable Y: Chiton_F
# Variables X: Area.lateral, Cov.Asc,Cov.Bryo,Cov.Spong,Cov.Others, cov.cca, cov.pey)
# Random factor: fSite

## Exploring variables ##

Z<-macro[,c("Area.lateral","cov.cca","cov.pey","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")]

#install.packages("usdm")
library(usdm)
vif(Z)
vifcor(Z)
vifstep(Z, th=3,keep = "cov.cca") # Cov.Other excluded

boxplot(macro[,c("cov.cca","cov.pey","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")])

pairs(macro[,c("Chiton_F","Area.lateral","cov.cca","cov.pey","Cov.Asc","Cov.Bryo","Cov.Spong","Cov.Others")],
      panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)

# Transformando Area.lateral
macro$logArea.lateral <- log(macro$Area.lateral)
str(macro)

boxplot(macro[,c("logArea.lateral","cov.cca","cov.pey","Cov.Asc","Cov.Bryo","Cov.Spong")])

pairs(macro[,c("Chiton_F","logArea.lateral","cov.cca","cov.pey","Cov.Asc","Cov.Bryo","Cov.Spong")],
      panel=panel.smooth, diag.panel=panel.hist, lower.panel=panel.cor)


# Make the coplot
coplot(Chiton_F ~ cov.cca | Site, data=macro, ylab = "Abudance",
       xlab = "Cover Fluoescence (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

coplot(Chiton_F ~ cov.pey | Site, data=macro, ylab = "Abudance",
       xlab = "Cover Fluoescence (%)",
       panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
         abline(tmp)
         points(x, y) })

library(beeswarm)
beeswarm::beeswarm(cov.cca ~ Site,data=macro)
beeswarm::beeswarm(cov.pey ~ Site,data=macro, col="red", pch=16, method="swarm")

dotchart(macro$cov.cca, main = "CCA", group = macro$fSite)
dotchart(macro$cov.pey, main = "PEY", group = macro$fSite)

stripchart(cov.cca ~ Site,data=macro, method="stack")
stripchart(cov.pey ~ Site,data=macro, method="stack")

library(lattice)

xyplot(Chiton_F ~ cov.cca | Site, data=macro, type=c("p","r"))
xyplot(Area.lateral ~ cov.pey| Site, data=macro, type=c("p","r"))

### ~GLMM Poisson ####

str(macro)
library(pscl)

f.full <- formula(Chiton_F ~ Area.lateral + cov.cca + cov.pey + Cov.Asc + Cov.Bryo
                  + Cov.Spong + (1|fSite))

library(glmmTMB)
mod2 <- glmmTMB(f.full,family = poisson, data=macro)
summary(mod2)

drop1(mod2, test="Chi") # Cov. Asc
f1 <- update(f.full,.~.-Cov.Asc)

mod2a <- glmmTMB(f1, family = poisson, data=macro)
summary(mod2a)

drop1(mod2a, test="Chi") # Cov. Spong
f2 <- update(f1,.~.-Cov.Spong)

mod2b <- glmmTMB(f2, family = poisson, data=macro)
summary(mod2b)

drop1(mod2b, test="Chi") # cov.cca
f3 <- update(f2,.~.-cov.cca)

mod2c <- glmmTMB(f3, family = poisson, data=macro)
summary(mod2c)

drop1(mod2c, test="Chi") # cov.pey
f4 <- update(f3,.~.-cov.pey)

mod2d <- glmmTMB(f4, family = poisson, data=macro)
summary(mod2d)

### ~AIC Selection ####

library(MuMIn)
model.sel(mod2, mod2a, mod2b, mod2c, mod2d) # mod2b (AIC: 305.4)
anova(mod2, mod2a, mod2b, mod2c, mod2d, test="Chi")

#install.packages("AICcmodavg")
library(AICcmodavg)

#setup a subset of models of Table 1
Cand.models <- list(mod2, mod2a, mod2b, mod2c, mod2d)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

r.squaredGLMM(mod2b)
r.squaredGLMM(mod2c)

# O modelo de menor AIC foi o mod2b mas pela comparacao de deltaAIC eh possivel
# considerar mod2b = mod2c, portanto vamos testar a validade de ambos

### ~Model Validation ####

## mod2b: ~ cov.cca + cov.pey + Area.lateral + Cov.Bryo + (1|fSite)

#install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod2b)
plot(simulationOutput)

# a. Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod2b,
                            n = 1000)
plot(simRes) # good

# b Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,2))
plotResiduals(simulationOutput, macro$cov.cca)
plotResiduals(simulationOutput, macro$cov.pey)

# c. Teste para inflacao de zero
par(mfrow = c(1,1))
testZeroInflation(simulationOutput)

simulationOutput <- simulateResiduals(fittedModel = mod2b)
plot(simulationOutput)

# O mod2b esta otimo!

## mod2c: ~ cov.pey + Area.lateral + Cov.Bryo + (1|fSite)
simRes <- simulateResiduals(fittedModel = mod2c,
                            n = 1000)
plot(simRes) # bad

# Entao, ficamos com o mod2b
summary(mod2b)

# cobertura de Peyssonelia eh positivamente significativa

# Vamos ver o mod1 e mod2 juntos:
model.sel(mod1b,mod2b) # mod2b (AIC: 305.4)
anova(mod1b, mod2b, test="Chi")

# Os modelos nao sao diferentes ou seja por esse metodo nao eh possivel dessociar
# o efeito visual do espcie especifico. E a peysonelia mesmo em menor prporcao
# tem importancia na deerminacao do abundancia de quitons

model.sel(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod2c, mod2d)

evidence.ratio <- 0.255/0.208
evidence.ratio # 1.22 Ou seja o mod1b eh somente 1.22 vezes mais provavel de ser
# o melhor modelo do que o mod2b que considera o efeito das especies de algas calcarias.

# Mod.averaging
mod.sel <- model.sel(mod1, mod1a, mod1b, mod2, mod2a, mod2b, mod2c, mod2d)
summary(model.avg(mod.sel, subset = delta < 2))


## CONCLUIDO ###



## Model Plus ####
# NAO VAi ENTRAR NO ARTIGO

# O objetivo eh investigar se eh possivel determinar a abundância dos quítons com
# base na presença de outros grupos de taxons na lateral da rocha evidando assim
# o explorador a necessiade de turn over da rocha para possível contatação.

# Hipótese: a presença de fissurelideos é um forte inicativo de maiores abundância
# de quiton I.pectinata.

# Variable Y: Chiton_F
# Variables X: Area.lateral, Bryo, Asci, Bival, Crabs, Shrimp, Worms, Spong,
# Urchins, Brittle, Coral, Zoanthus,Flatworm, Hermit, Fissu.
# Random factor: fSite

macro <- read.csv("rockside.csv", sep=";",dec=".", header=T)
str(macro)


# Tranformando os dados de abundancias de coespecies em presenca-ausencia:

library(vegan) # decostand() function

macro.pa <- decostand(macro[,15:29], "pa")
macro.pa <- cbind(macro[,1:14],macro.pa)
str(macro.pa)
macro.pa$fSite <- factor(macro.pa$Site)

## Exploring variables ##

Z<-macro.pa[,c("Area.lateral","Bryo","Asci","Bival","Crabs","Shrimp","Worms","Spong",
            "Urchins","Brittle","Coral","Zoanthus","Flatworm","Hermit","Fissu")]

#install.packages("usdm")
library(usdm)
vif(Z)
vifcor(Z)

dotchart(Z$Fissu) # weight no eixo x
plot(Z$Fissu)  # Compare com a funcao generica 'plot'
plot(density(Z$Fissu))
plot(table(Z$Fissu))
plot()
hist(Z$Fissu)
stripchart(Fissu ~ Site,data=macro.pa, method="stack")
library(beeswarm)
beeswarm::beeswarm(Z$Fissu)
beeswarm::beeswarm(Fissu ~ Site,data=macro.pa)
beeswarm::beeswarm(Fissu ~ Site,data=macro.pa, col="red", pch=16, method="swarm")

head(macro.pa[,-c(2:14)])
class(macro.pa)

macro.mx <- matrix(macro.pa[,-c(2:14)])
class(macro.mx)
dim(macro.mx)
str(macro.mx)

specie <- c("Bryo","Asci","Bival","Crabs","Shrimp","Worms","Spong",
            "Urchins","Brittle","Coral","Zoanthus","Flatworm","Hermit","Fissu")


library(pscl)
f.full <- formula(Chiton_F ~ Bryo+Asci+Bival+Crabs+Shrimp+Worms+Spong+Urchins+
                    Brittle+Coral+Zoanthus+Flatworm+Hermit+Fissu + (1|fSite))

### ~GLMM Poisson ####

library(glmmTMB)
mod2 <- glmmTMB(f.full,family = poisson, data=macro.pa)
summary(mod2)

drop1(mod2, test="Chi") # Hermit
f1 <- update(f.full,.~.-Hermit)
mod2a <- glmmTMB(f1, family = poisson, data=macro.pa)
summary(mod2a)

drop1(mod2a, test="Chi") # Worms
f2 <- update(f1,.~.-Worms)
mod2b <- glmmTMB(f2, family = poisson, data=macro.pa)
summary(mod2b)

drop1(mod2b, test="Chi") # Coral
f3 <- update(f2,.~.-Coral)
mod2c <- glmmTMB(f3, family=poisson, data=macro.pa)
summary(mod2c)

drop1(mod2c, test="Chi") # Spong
f4 <- update(f3,.~.-Spong)
mod2d <- glmmTMB(f4, family=poisson, data=macro.pa)
summary(mod2d)

drop1(mod2d, test="Chi") # Asci
f5 <- update(f4,.~.-Asci)
mod2e <- glmmTMB(f5, family=poisson, data=macro.pa)
summary(mod2e)

drop1(mod2e, test="Chi") # Shrimp
f6 <- update(f5,.~.-Shrimp)
mod2f <- glmmTMB(f6, family=poisson, data=macro.pa)
summary(mod2f)

drop1(mod2f, test="Chi") # Bryo
f7 <- update(f6,.~.-Bryo)
mod2g <- glmmTMB(f7, family=poisson, data=macro.pa)
summary(mod2g)

drop1(mod2g, test="Chi") # Flatworm
f8 <- update(f7,.~.-Flatworm)
mod2h <- glmmTMB(f8, family=poisson, data=macro.pa)
summary(mod2h)

drop1(mod2h, test="Chi") # Bival
f9 <- update(f8,.~.-Bival)
mod2i <- glmmTMB(f9, family=poisson, data=macro.pa)
summary(mod2i)

drop1(mod2i, test="Chi") # Urchins
f10 <- update(f9,.~.-Urchins)
mod2j <- glmmTMB(f10, family=poisson, data=macro.pa)
summary(mod2j)

drop1(mod2j, test="Chi") # Zoanthus
f11 <- update(f10,.~.-Zoanthus)
mod2l <- glmmTMB(f11, family=poisson, data=macro.pa)
summary(mod2l)

#### + AIC Selection ####

library(MuMIn)
model.sel(mod2,mod2a,mod2b,mod2c,mod2d,mod2e,mod2f,mod2g,mod2h,mod2i,mod2j,mod2l)
anova(mod2,mod2a,mod2b,mod2c,mod2d,mod2e,mod2f,mod2g,mod2h,mod2i,mod2j,mod2l, test="Chi")

#install.packages("AICcmodavg")
library(AICcmodavg)

#setup a subset of models of Table 1
Cand.models <- list(mod2,mod2a,mod2b,mod2c,mod2d,mod2e,mod2f,mod2g,mod2h,mod2i,mod2j,mod2l)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

r.squaredGLMM(mod2i)

#### + Model Validation ####

#install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod2i)
plot(simulationOutput) # Not ajusted

# a. Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod2i,
                            n = 1000)
plot(simRes) # Not ajusted

### ~GLMM nbinom1 ####

library(glmmTMB)
nb1 <- glmmTMB(f.full,family = nbinom1, data=macro.pa)
summary(nb1)

drop1(nb1, test="Chi") # Hermit
f1 <- update(f.full,.~.-Hermit)
nb2 <- glmmTMB(f1, family = nbinom1, data=macro.pa)
summary(nb2)

drop1(nb2, test="Chi") # Spong
f2 <- update(f1,.~.-Spong)
nb3 <- glmmTMB(f2, family = nbinom1, data=macro.pa)
summary(nb3)

drop1(nb3, test="Chi") # Worms
f3 <- update(f2,.~.-Worms)
nb4 <- glmmTMB(f3, family=nbinom1, data=macro.pa)
summary(nb4)

drop1(nb4, test="Chi") # Coral
f4 <- update(f3,.~.-Coral)
nb5 <- glmmTMB(f4, family=nbinom1, data=macro.pa)
summary(nb5)

drop1(nb5, test="Chi") # Asci
f5 <- update(f4,.~.-Asci)
nb6 <- glmmTMB(f5, family=nbinom1, data=macro.pa)
summary(nb6)

drop1(nb6, test="Chi") # Bival
f6 <- update(f5,.~.-Bival)
nb7 <- glmmTMB(f6, family=nbinom1, data=macro.pa)
summary(nb7)

drop1(nb7, test="Chi") # Shrimp
f7 <- update(f6,.~.-Shrimp)
nb8 <- glmmTMB(f7, family=nbinom1, data=macro.pa)
summary(nb8)

drop1(nb8, test="Chi") # Flatworm
f8 <- update(f7,.~.-Flatworm)
nb9 <- glmmTMB(f8, family=nbinom1, data=macro.pa)
summary(nb9)

drop1(nb9, test="Chi") # Bryo
f9 <- update(f8,.~.-Bryo)
nb10 <- glmmTMB(f9, family=nbinom1, data=macro.pa)
summary(nb10)

drop1(nb10, test="Chi") # Zoanthus
f10 <- update(f9,.~.-Zoanthus)
nb11 <- glmmTMB(f10, family=nbinom1, data=macro.pa)
summary(nb11)

drop1(nb11, test="Chi") # Brittle
f11 <- update(f10,.~.-Brittle)
nb12 <- glmmTMB(f11, family=nbinom1, data=macro.pa)
summary(nb12)

drop1(nb12, test="Chi") # Urchins
f12 <- update(f11,.~.-Urchins)
nb13 <- glmmTMB(f12, family=nbinom1, data=macro.pa)
summary(nb13)

#### + AIC Selection ####

library(MuMIn)

model.sel(nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8,nb9,nb10,nb11,nb12,nb13)
anova(nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8,nb9,nb10,nb11,nb12,nb13, test="Chi")

#install.packages("AICcmodavg")
library(AICcmodavg)

#setup a subset of models of Table 1
Cand.models <- list(nb1,nb2,nb3,nb4,nb5,nb6,nb7,nb8,nb9,nb10,nb11,nb12,nb13)

##create a vector of names to trace back models in set
Modnames <- paste("mod", 1:length(Cand.models), sep = " ")

##generate AICc table
aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE)

##round to 4 digits after decimal point and give log-likelihood
print(aictab(cand.set = Cand.models, modnames = Modnames, sort = TRUE),
      digits = 4, LL = TRUE)

r.squaredGLMM(nb13)

#### + Model Validation ####

#install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = nb13)
plot(simulationOutput) # Ajusted

# a. Randomizacoes:
simRes <- simulateResiduals(fittedModel = nb13,
                            n = 1000)
plot(simRes) # Ajusted

model.sel(mod2h,nb13) # nb13 menor AIC

## CONCLUIDO ###

#========

# Corvif fuction

#Library files for courses provided by: Highland Statistics Ltd.
#To cite these functions, use:
#Mixed effects models and extensions in ecology with R. (2009).
#Zuur, AF, Ieno, EN, Walker, N, Saveliev, AA, and Smith, GM. Springer.

#Copyright Highland Statistics LTD.

#####################################################################
#VIF FUNCTION.
#To use:  corvif(YourDataFile)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  #cat("Correlations of the variables\n\n")
  #tmp_cor <- cor(dataz,use="complete.obs")
  #print(tmp_cor)

  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)

  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}


#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}



# ZEro inflated - NB

poiZi1 <- glmmTMB(f2,family = nbinom1, zi=~Chiton_F, data=macro)
summary(poiZi1)

simulationOutput <- simulateResiduals(fittedModel = poiZi1)
plot(simulationOutput)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = poiZi1,
                            n = 1000)
plot(simRes) # good


poiZi1 <- glmmTMB(f2,family = poisson, zi=~., data=macro)
summary(poiZi1)

simulationOutput <- simulateResiduals(fittedModel = poiZi1)
plot(simulationOutput)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = poiZi1,
                            n = 1000)
plot(simRes) # good


# Vamos comparar os modelos:

library(MuMIn)
model.sel(mod1, mod2, mod3, mod4) # mod4

# vendo a diferenca entre os modelos pela ANOVA:
anova(mod1, mod2, mod3, mod4)

# R2
r.squaredGLMM(mod4)

#install.packages("AICcmodavg")
detach(MuMIn)
library(AICcmodavg)

# setup a subset of models of Table 1
Cand.models <- list(mod1, mod2, mod3, mod4)

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
vif(mod4) # Nao deu certo

# Validacao:
#install.packages("DHARMa")
library(DHARMa)
simulationOutput <- simulateResiduals(fittedModel = mod4)
plot(simulationOutput)

# Teste para inflacao de zero
testZeroInflation(simulationOutput)

# b. Detecting missing predictors or wrong functional assumptions

testUniformity(simulationOutput = simulationOutput)
par(mfrow = c(1,2))
plotResiduals(simulationOutput, data$Cov.Bryo)
plotResiduals(simulationOutput, data$Area.total)

# Randomizacoes:
simRes <- simulateResiduals(fittedModel = mod4,
                            n = 1000)
plot(simRes) # good

# MODELO AJUSTADO COM POISSON!!!!
# Prinipais variáveis


# Load necessary libraries
library(MASS)

# Generate example data
set.seed(123)
n <- 100
covariate <- rbeta(n, shape1 = 2, shape2 = 2)  # Generating covariate following beta distribution
predictor <- rnorm(n)  # Generating predictor variable (just a random normal distribution for demonstration)
response <- 3 + 2 * covariate + rnorm(n)  # Generating response variable (linear relationship with covariate)

# Fit linear regression model with beta-distributed covariate
model <- lm(Chiton_F ~ Fissu + Crabs, data=macro.pa)

# Summarize model
summary(model)

# Plot data and fitted model
plot(macro.pa$Fissu, macro.pa$Chiton_F, main = "Linear Regression with Beta-Distributed Covariate", xlab = "Covariate", ylab = "Response")
abline(model, col = "red")







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
plotResiduals(simulationOutput, macro.pa$Fissu)
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


#CODIGOS SOLTOS ####

# Discritivo
summary(data)

aggregate(data$Samp.time, list(data$Site), FUN=mean)

anova(aov(Samp.time ~ Site, data=data))
TukeyHSD(aov(Samp.time ~ Site, data=data),ordered=T)


install.packages("AED")
library (AED)
data (Koalas)

install.packages("ncf")
library(ncf)
Correlog <- spline.correlog(x = macro$x,
                            y = macro$y,
                            z = macro$Chiton_F,  xmax = 10000)

plot(Correlog)
summary(Correlog)
