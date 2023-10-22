
library(ggplot2)
library(gridExtra)
library(fpc)
library(mclust)
library(PerformanceAnalytics)
library(MuMIn)
library(metasens)
library(tidyverse)
library(wildmeta)
library(dplyr)
library(tidybayes)
library(ggplot2)
library(ggridges)
library(glue)
library(stringr)
library(forcats)
library(brms)
library(pema)
library(mice)
library(loo)
library(ggthemes)

# VALENCIA	tomamos como referencia la escala -3 a +3
# AROUSAL	tomamos como referencia la escala 1 a 5
# FRECUENCIA	tomamos como referencia la métrica ftot/millón
# IMAGINABILIDAD	tomamos como referencia de 1 a 7
# EDAD ADQUISICIÓN	tomamos como referencia la escala de 1 a 7
# CONCRECIÓN	tomamos como referencia la escala de 1 a 7


#Importamos la matriz de datos
library(readxl)
VNN <- read_excel("H:/ALBERTO/Metaanalisis PALABRAS/RESULTADOS/Analisis.xlsx", sheet="OKPOSNEUrt")

# Calculamos el tamaño del efecto y su error de medida y lo incorporamos a la matriz
# Usamos standardized mean change con raw score standardization [SMCR], corrigiendo por la desviación típica de NEUTRO
# Usamos un dato de correlación conservador de 0.5 (cambiar el valor 27 por el tamaño de la matriz en cuestión)

VNN$r<-rep(0.5,nrow(VNN))
VNN$idx<-c(1:nrow(VNN))
VNN$mi<-VNN$N-1

library(meta)
library(dmetar)
library(metafor)

#dat <- escalc(measure="SMCR", m1i=Mneg, sd1i=SDneu, ni=N, m2i=Mneu, ri=r, slab=ID, data=VNN)

VNN$dif<-VNN$Mpos-VNN$Mneu
VNN$sdif<-sqrt((VNN$SDpos^2 + VNN$SDneu^2)-2*VNN$r*VNN$SDpos*VNN$SDneu)

VNN$yi_un<-VNN$dif/VNN$sdif
VNN$yi<-VNN$yi_un*(1-(3/((4*VNN$mi)-1)))
VNN$vi<-(1/VNN$N)+(VNN$yi^2/(2*VNN$N))

VNN$sei<-sqrt(VNN$vi)
VNN$lower <- VNN$yi - 1.96*VNN$sei
VNN$upper <- VNN$yi + 1.96*VNN$sei


dat<-VNN

#Realizmos el meta multinivel, incluyendo como estrato intermedio la variable ID, donde el mismo estudio puede aportar varios tamaños del efecto (para controlar variable predictora), por lo que los tamaños del efecto serían dependientes
full.model <- rma.mv(yi = yi, 
                     V = vi, 
                     slab = ID,
                     data = dat,
                     random = ~ 1 | ID/idx, 
                     test = "t", 
                     method = "REML")

summary(full.model)
forest(full.model, cex=.75 , header=c("Author","SMCR"))

# Componentes de varianza:
# sigma^2.1 indica la varianza inter-cluster-> heterogeneidad entre estudios
# sigma^2.2 indica la varianza intra-cluster


i2 <- var.comp(full.model)
summary(i2)
plot(i2)

## Influence analysis ###

inf <- cooks.distance(full.model)
plot(inf)
Number=sum(inf>1,na.rm=T)


#### Sesgo de publicación: funnel + Egger test

funnel(full.model, label="out", cex=0.5, level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), refline=0, legend=FALSE)


# Modelo de EFECTOS MIXTOS de TRES NIVELES

## Para implementar la meta-regresión, primero evaluamos el nivel de multicolinealidad entre los predictores.

#dat[,c("VALneg", "VALneu", "AROneg", "AROneg", "AROneu", "FREQneg", "FREQneu", "LENneg", "LENneu","CONneg", "CONneu")] %>% 
#chart.Correlation()

## Constatamos elevados niveles de multicolinealidad y pasamos a reformular los predictores como el promedio de cada dimensión/diferencia entre dimensiones

dat$VALm <-apply(dat[,8:9],1,mean)
dat$AROm <-apply(dat[,10:11],1,mean)
dat$FREQm <-apply(dat[,12:13],1,mean)
dat$LENm <-apply(dat[,14:15],1,mean)
dat$CONm <-apply(dat[,20:21],1,mean)

dat$VALd <-apply(dat[,8:9],1,diff)
dat$AROd <-apply(dat[,10:11],1,diff)
dat$FREQd <-apply(dat[,12:13],1,diff)
dat$LENd <-apply(dat[,14:15],1,diff)
dat$CONd <-apply(dat[,20:21],1,diff)


#dat[,c("VALm", "AROm", "FREQm", "LENm", "CONm")] %>% 
  #chart.Correlation()

#dat[,c("VALd", "AROd", "FREQd", "LENd", "CONd")] %>% 
  #chart.Correlation()



## BAYES ##

prior4 <- c(
  prior(normal(0, 1), coef = Intercept),
  prior(cauchy(0, 0.5), class = sd)
)


# Model
mod4 <- brm(
  yi | se(sqrt(vi) ) ~ 0 + Intercept + (1|EXPERIMENT) + (1|idx),
  data = dat,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)


## FOREST PLOT ##

ss_draws <- mod4 %>%
  spread_draws(b_Intercept, r_EXPERIMENT[EXPERIMENT,]) %>%
  # add the grand mean to the group-specific deviations
  mutate(mu = b_Intercept + r_EXPERIMENT) %>%
  dplyr::select(-b_Intercept) %>% 
  ungroup()


ss <- full_join(ss_draws,dat %>% dplyr::select(EXPERIMENT,ID,year), by = c("EXPERIMENT"))
avg_draws <- spread_draws(mod4, b_Intercept) %>% 
  mutate(EXPERIMENT = max(dat$EXPERIMENT) + 3,
         r_EXPERIMENT = NA,
         mu = b_Intercept,
         ID = "Overall Effect",
         year = 9999) %>% dplyr::select(-b_Intercept)

ss <- bind_rows(ss,avg_draws)
ss$year <- as.integer(ss$year)


col_pal<- c(
  "grey", ## summary effect
        "white", ## white space at bottom between summary effect size and study-specific effect sizes
        "white",
        "#FFCC33",
        "#FFBF33",
        "#FFB333",
        "#FFA633",
        "#FF9333",
        "#FF8C33",
        "#FF8033",
        "#FF7333",
        "#FF6633",
        "#FF3926",
        "#FF1A2B",
        "#FF0D46",
        "#FF0066",
        "#F80276",
        "#F00385",
        "#E90592",
        "#E2079F",
        "#DB08AA",
        "#D409B3",
        "#CD0BBC",
        "#B40DBF",
        "#9900CC",
        "#9C00F2",
        "#9A06FF",
        "#961AFF",
        "#9340FF",
        "#9966FF",
        "#8C66FF",
        "#7F66FF",
        "#6666FF",
        "#667FFF",
        "#668CFF",
        "#6699FF",
        "#55A6F6",
        "#47B4EB",
        "#3CC2DD",
        "#33CCCC",
        "#38CCC3",
        "#40CCBA",
        "#5ACCA1",
        "#66CC99")
        
overall_index <- max(dat$EXPERIMENT) + 3

test <- ss

forPlot1 <- test %>%
  # plot
  ggplot(aes(x = mu, y = reorder(ID,desc(EXPERIMENT)),fill = reorder(ID,desc(EXPERIMENT)))) +
  #ggplot(aes(x = mu, y = plot_label,fill = plot_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_point(aes(x = yi, y =  reorder(ID,desc(EXPERIMENT))), 
             
             ##add summary effect size
             data = dat %>%
               dplyr::add_row(EXPERIMENT = overall_index, ID = "Overall Effect", year = 9999) %>% 
               dplyr::add_row(EXPERIMENT = overall_index- 1, ID = "", year = 9998) %>% 
               dplyr::add_row(EXPERIMENT = overall_index- 2, ID = " ", year = 9997), 
             
             shape = 21, size = 1.5) +
  geom_density_ridges(scale = 2.2, rel_min_height = 0.01, alpha = .9) +
  #scale_fill_viridis(discrete = TRUE,option="magma") +
  scale_fill_manual(values = col_pal) +
  scale_x_continuous(limits = c(-2, 4.5), breaks = seq(-2,4,.5)) +
  scale_y_discrete(expand = expansion(add = c(1,1.5)), drop = F)+
  #scale_y_continuous(limits = c(-3,40)) +
  cowplot::theme_half_open(10) +
  theme(legend.position = "none", plot.background = element_rect(fill="white")) +
  cowplot::panel_border() +
  ylab(NULL) +
  xlab(expression(paste("Effect Size (Hedges' ", italic(g), ")"))) +
  coord_cartesian(xlim = c(-2,4.5))

forPlot1


#Tests for publication bias
dat$sei=sqrt(dat$vi)

mod4_egger <- brm(
  yi | se(sei) ~ 0 + Intercept + sei + (1|EXPERIMENT) + (1|idx),
  data = dat,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(mod4_egger, "sei > 0")


### Penalized Metanalysis ###
df <- dat[ , c("yi", "vi", "VALd",
               "AROd", "FREQd", "LENd")]
df <- mice(df)
fit1 <- brma(yi ~ ., data = df, vi = "vi", chains= 4, iter=2000, control=list(adapt_delta=0.99), seed = 1)


df <- dat[ , c("yi", "vi", "VALm",
                       "AROm", "FREQm", "LENm")]
df <- mice(df)
fit2 <- brma(yi ~ ., data = df, vi = "vi", chains= 4, iter=2000, control=list(adapt_delta=0.99), seed = 1)


#traceplot(as.stan(fit), pars = c("Intercept"))
#check_hmc_diagnostics(as.stan(fit))


## Bayesian Regresion ##
df <- dat[ , c("yi", "vi", "EXPERIMENT", "idx", "VALd",
               "AROd", "FREQd", "LENd")]
df <- mice(df)

imp1<-complete(df,1)

bm11 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + VALd + (1|EXPERIMENT) + (1|idx),
  data = imp1,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm11, "(Intercept + VALd) > Intercept")


df <- dat[ , c("yi", "vi", "EXPERIMENT", "idx", "VALm",
               "AROm", "FREQm", "LENm")]
df <- mice(df)

imp2<-complete(df,1)

bm12 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + VALm + (1|EXPERIMENT) + (1|idx),
  data = imp2,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm12, "(Intercept + VALm) > Intercept")


bm21 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + AROd + (1|EXPERIMENT) + (1|idx),
  data = imp1,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm21, "(Intercept + AROd) < Intercept")

bm22 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + AROm + (1|EXPERIMENT) + (1|idx),
  data = imp2,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm22, "(Intercept + AROm) > Intercept")

bm31 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + FREQd + (1|EXPERIMENT) + (1|idx),
  data = imp1,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm31, "(Intercept + FREQd) > Intercept")


bm32 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + FREQm + (1|EXPERIMENT) + (1|idx),
  data = imp2,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm32, "(Intercept + FREQm) > Intercept")

bm41 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + LENd + (1|EXPERIMENT) + (1|idx),
  data = imp1,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm41, "(Intercept + LENd) > Intercept")


bm42 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + LENm + (1|EXPERIMENT) + (1|idx),
  data = imp2,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99)
)

hypothesis(bm42, "(Intercept + LENm) > Intercept")


bm51 <- brm(
  yi | se(sqrt(vi)) ~ 0 + Intercept + VALd * AROd + (1|EXPERIMENT) + (1|idx),
  data = imp1,
  prior = prior4,
  save_all_pars = TRUE,
  warmup = 5000, iter = 1e4,
  cores = parallel::detectCores(),
  control = list(adapt_delta = .99, max_treedepth=12)
)

hypothesis(bm51, "(Intercept + AROd * VALd) > Intercept")



