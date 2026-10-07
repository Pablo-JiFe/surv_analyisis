install.packages("survival")
install.packages("survminer")

library(ggplot2)
library(dplyr)
library(survival)
library(survminer)

# Para que sea lo mismo para todos

set.seed(37)

# El "dataset" simulado de pacientes

treatment_coh <- 
  data.frame(
  pacientes = paste0("p", seq(1, 200)),
  
  global = rep(x = 1, 200),
  
  grupo = sample(c("Medicamento", "Sin medicamento"), 200, replace = TRUE, prob = c(0.5, 0.5)),
  
  tiempo = rexp(200, rate = 0.5) * 20
  
) %>% 
  mutate(
    tiempo = ifelse(
      tiempo < 0,
                    abs(tiempo),
                    tiempo
                    ),
    evento = ifelse(grupo == "Medicamento",
                   sample(c(1, 0), 100, replace = TRUE, prob = c(0.25, 0.75)),
                   sample(c(1, 0), 100, replace = TRUE, prob = c(0.75, 0.25))
    ),
    enf_cardiaca = ifelse(evento == 1, 
                         sample(c(1, 0), 100, replace = TRUE, prob = c(0.45, 0.55)),
                         sample(c(1, 0), 100, replace = TRUE, prob = c(0.55, 0.45))
    ),
    grado = ifelse(evento == 1,
                   sample(c(1, 2, 3, 4), 100, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.4)),
                   sample(c(1, 2, 3, 4), 100, replace = TRUE, prob = c(0.1, 0.3, 0.4, 0.2))
    ),
    edad = ifelse(evento == 1,
                  rnorm(100, mean = 60, sd = 25),
                  rnorm(100, mean = 40, sd = 30)
    ),
    enf_renal =  ifelse(enf_cardiaca == 1,
                        sample(c(1, 0), 100, replace = TRUE, prob = c(0.95, 0.05)),
                        sample(c(1, 0), 100, replace = TRUE, prob = c(0.05, 0.95))
    ),
    marcador_sangre = ifelse(evento == 0 & edad > 50,
                      1,
                      0
                      ),
    marcador_clinico = ifelse(evento == 0 & tiempo > 80 | evento == 1 & tiempo < 80,
                      1,
                      0
    )
  )

# Para todos los subsequentes analysis se requiere de una formula
# La formula consiste en un objeto de supervivencia surv_obj
# y aquello que va a ser utilizado para comparaciones

# Primero para hacer el surv_obj, este sera util para todo lo que sigue

surv_obj <- Surv(time = treatment_coh$tiempo, event = treatment_coh$evento)

# Objeto de Kaplan meier independiente de grupo

fit_km <- survfit(formula = surv_obj ~ treatment_coh$global, data = treatment_coh)

# Graficar

ggsurvplot(fit_km,
           conf.int = FALSE,
           ggtheme = theme_classic(base_size = 30))

# Objeto Kaplan Meier para comparar entre grupos

fit_km_grupo <- survfit(formula = surv_obj ~ treatment_coh$grupo, data = treatment_coh)

# Graficar

ggsurvplot(
  fit_km_grupo,
  ggtheme = theme_classic(base_size = 25)
           )

# Log rank test comparando medicamento contra no medicamento

log_rank_grupo <- survdiff(surv_obj ~ grupo, data = treatment_coh)

print(log_rank_grupo)


# Log rank test comparando entre diuferentes grados de esta enfermedad

log_rank_grado <- survdiff(surv_obj ~ grado, data = treatment_coh)

print(log_rank_grado)

# Regresion tipo Cox simple

cox_univariado <- coxph(surv_obj ~ grupo, data = treatment_coh)

print(cox_univariado)

summary(cox_univariado)

# Regresion multivariada

cox_multivariado1 <- coxph(surv_obj ~ grupo + edad, data = treatment_coh)

summary(cox_multivariado1)

cox_multivariado2 <- coxph(surv_obj ~ grupo + enf_cardiaca + edad, data = treatment_coh)

summary(cox_multivariado2)

# Multicolinearidad

cox_multicolinearidad <- coxph(surv_obj ~ grupo + enf_cardiaca + enf_renal, data = treatment_coh)

summary(cox_multicolinearidad)

# Separación

coxph(surv_obj ~ grupo + marcador_sangre, data = treatment_coh)

# Supuesto de proporcionalidad

# Cumple

supuesto_ph <- cox.zph(cox_multivariado1)

ggcoxzph(supuesto_ph)

# No cumple

cox_supuesto <- coxph(surv_obj ~ grupo + marcador_clinico, data = treatment_coh)

supuesto_ph_no <- cox.zph(cox_supuesto)

print(supuesto_ph_no)

ggcoxzph(supuesto_ph_no)

# Arreglado

cox_supuesto_arreglado <- coxph(surv_obj ~ grupo + strata(marcador_clinico), data = treatment_coh)

supuesto_ph_arreglado <- cox.zph(cox_supuesto_arreglado)

print(supuesto_ph_arreglado)

ggcoxzph(supuesto_ph_arreglado)

# Transcriptomica

set.seed(37)

# Simulando datos de transcriptomica centrados y normalizados

transcriptomica_sim <- 
  data.frame(
  pacientes = treatment_coh$pacientes,
  CD20 = rnorm(200, mean = 0, sd = 1),
  TNF = rnorm(200, mean = 0, sd = 1),
  APOE = rnorm(200, mean = 0, sd = 1),
  EGFR = rnorm(200, mean = 0, sd = 1),
  ESR1 = rnorm(200, mean = 0, sd = 1)
  )

# Junto datos y agrego genes a proposito para que salgan positivos

treat_trans <- 
  treatment_coh %>% 
  select(pacientes,
         evento,
         tiempo) %>% 
  left_join(transcriptomica_sim, by = "pacientes") %>% 
  tibble::column_to_rownames("pacientes") %>% 
  mutate(
    NFKB1 = ifelse(evento == 1,
                   rnorm(100, mean = - 0.8, sd = 0.8),
                   rnorm(100, mean = 0.2, sd = 0.8)
    ),
    IL10 = ifelse(evento == 1,
                  rnorm(100, mean = 0.2, sd = 0.8),
                  rnorm(100, mean = - 0.8, sd = 0.8)
                  ),
    ESR1 = ifelse(evento == 1,
                  rnorm(100, mean = 0.1, sd = 0.9),
                  rnorm(100, mean = - 0.9, sd = 0.9)
    )
  )

# surv_obj

surv_obj_trans <- Surv(treat_trans$tiempo, treat_trans$evento)

# Quito los desenlaces porque si no sesgan el resultado

treat_trans <- 
  treat_trans %>% 
  select(-c(evento, 
            tiempo))
# Modelo cox

summary(coxph(surv_obj_trans ~ ., data = treat_trans))
