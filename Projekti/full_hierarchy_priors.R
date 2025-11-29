library(dplyr)
library(brms)


data <- read.csv("data.csv")

data$bacteremia <- factor(data$bacteremia)
data <- data %>%
  mutate(
    bact_num = ifelse(bacteremia == levels(bacteremia)[2], 1, 0),
    sex      = factor(sex)
  )


data <- data %>%
  mutate(
    age_group = cut(
      age,
      breaks = c(0, 39, 64, 79, Inf),
      labels = c("0-39", "40-64", "65-79", "80+"),
      right  = TRUE
    )
  )


data <- data %>%
  mutate(
    CRP_log = log(crp + 0.1)
  )

data_model <- data %>%
  mutate(
    CRP_log_s = as.numeric(scale(CRP_log)),
    WBC_s     = as.numeric(scale(wbc)),
    NEU_s     = as.numeric(scale(neu)),
    LYM_s     = as.numeric(scale(lym)),
    PLT_s     = as.numeric(scale(plt)),
    CREA_s    = as.numeric(scale(crea)),
    ALB_s     = as.numeric(scale(alb)),
    BUN_s     = as.numeric(scale(bun)),
    GBIL_s    = as.numeric(scale(gbil)),
    RBC_s     = as.numeric(scale(rbc)),
    SODIUM_s  = as.numeric(scale(sodium)),
    AGE_s     = as.numeric(scale(age))
  ) %>%
  select(
    bact_num, bacteremia, sex, age_group,
    CRP_log_s, WBC_s, NEU_s, LYM_s, PLT_s,
    CREA_s, ALB_s, BUN_s, GBIL_s, RBC_s,
    SODIUM_s, AGE_s
  )


data_model <- na.omit(data_model)

# Testasin full mallia
f_hier <- bf(
  bact_num ~ CRP_log_s + WBC_s + NEU_s + LYM_s + PLT_s +
    CREA_s + ALB_s + BUN_s + GBIL_s + RBC_s + SODIUM_s +
    AGE_s + sex +
    (1 + CREA_s + RBC_s | sex) +
    (1 + CREA_s + ALB_s + BUN_s + RBC_s | age_group)
)


prev <- mean(data_model$bact_num)
logit_prev <- qlogis(prev)
print(logit_prev)

priors <- c(

  prior(normal(-2.44, 0.5), class = "Intercept"), #Datasta napattu se baserate bakteremia, en tiiä onks laitonta (ei kai?)
  
  prior(normal(0.7, 0.3), coef = "CRP_log_s"),  # korkee CRP -> korkee riski
  prior(normal(0.3, 0.3), coef = "WBC_s"),      # korkee WBC -> korkee riski (mutta voi olla korrelaation takia tää pitäis muuttaa)
  prior(normal(1.2, 0.5), coef = "NEU_s"),
  prior(normal(-0.5, 0.4), coef = "LYM_s"),
  prior(normal(-0.3, 0.3), coef = "PLT_s"),
  prior(normal(0.3, 0.3), coef = "CREA_s"),
  prior(normal(0.2, 0.2), coef = "ALB_s"),
  prior(normal(0.4, 0.3), coef = "BUN_s"),
  prior(normal(0.3, 0.2), coef = "GBIL_s"),
  prior(normal(0, 0.3), coef = "RBC_s"),
  prior(normal(-0.5,0.2), coef = "SODIUM_s"),
  prior(normal(0.4, 0.2), coef = "AGE_s"),
  

  prior(exponential(2), class = "sd", group = "sex", coef = "Intercept"),
  prior(exponential(1.5), class = "sd", group = "sex", coef = "CREA_s"),
  prior(exponential(1.5), class = "sd", group = "sex", coef = "RBC_s"),
  
  prior(exponential(2), class = "sd", group = "age_group", coef = "Intercept"),
  prior(exponential(1.5), class = "sd", group = "age_group", coef = "CREA_s"),
  prior(exponential(1.5), class = "sd", group = "age_group", coef = "ALB_s"),
  prior(exponential(1.5), class = "sd", group = "age_group", coef = "BUN_s"),
  prior(exponential(1.5), class = "sd", group = "age_group", coef = "RBC_s")
)


fit_hier <- brm(
  formula = f_hier,
  data    = data_model,
  family  = bernoulli(),
  prior   = priors,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
)

summary(fit_hier)
plot(fit_hier)
pp_check(fit_hier)
