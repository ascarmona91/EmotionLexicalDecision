Sys.setenv(JAGS_HOME="C:\\Program Files\\JAGS\\JAGS-4.3.1")

fit <- RoBMA(d = dat$yi_un, v = dat$vi, effect_direction="negative", study_ids = dat$EXPERIMENT,
             priors_effect_null        = NULL,
             priors_heterogeneity_null = NULL,
             priors_hierarchical_null= NULL,
             parallel = TRUE, seed = 1)

fit<-RoBMA(d=dat$yi, v=dat$vi, effect_direction="negative", study_ids = dat$EXPERIMENT,parallel = TRUE,seed = 1)