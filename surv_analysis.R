# 1. Cargar paquetes -------------------------------------------------------
pacman::p_load("survival",  "survminer", "ggplot2", "dplyr", "vroom" )


# 2.- Leer los datos para este ejercicio -----------------------------------

treatment_coh <- vroom( file = "https://data.biofreelancer.com/survivalcoh" )
treat_trans   <- vroom( file = "https://data.biofreelancer.com/survivaltrans" )

# Para todos los subsequentes analysis se requiere de una formula
# La formula consiste en un objeto de supervivencia surv_obj
# y aquello que va a ser utilizado para comparaciones

# Primero para hacer el surv_obj, este sera util para todo lo que sigue

surv_obj <- Surv(time = treatment_coh$tiempo, event = treatment_coh$evento)



# 3.- Curvas Kaplan Meier -------------------------------------------------


# Objeto de Kaplan meier independiente de grupo

fit_km <- survfit(formula = surv_obj ~ treatment_coh$global, data = treatment_coh)

# Graficar

ggsurvplot(fit = fit_km,
           conf.int = FALSE,
           ggtheme = theme_classic(base_size = 30))

# Objeto Kaplan Meier para comparar entre grupos

fit_km_grupo <- survfit(formula = surv_obj ~ treatment_coh$grupo, data = treatment_coh)

# Graficar

ggsurvplot(
  fit_km_grupo,
  ggtheme = theme_classic(base_size = 25),
  pval = TRUE,
  risk.table = TRUE,
  break.time.by = 50,
  break.y.by = 0.5,
  palette = c("red", "blue")
           )


# 4.- Log rank test -------------------------------------------------------


# Log rank test comparando medicamento contra no medicamento

log_rank_grupo <- survdiff(formula = surv_obj ~ grupo, data = treatment_coh)

print(log_rank_grupo)


# Log rank test comparando entre diuferentes grados de esta enfermedad

log_rank_grado <- survdiff(surv_obj ~ grado, data = treatment_coh)

print(log_rank_grado)


# 5.- Cox -----------------------------------------------------------------


# Regresion tipo Cox simple

cox_univariado <- coxph(formula = surv_obj ~ grupo, data = treatment_coh)

print(cox_univariado)

summary(cox_univariado)

# Regresion multivariada

cox_multivariado1 <- coxph(surv_obj ~ grupo + edad, data = treatment_coh)

summary(cox_multivariado1)

cox_multivariado2 <- coxph(surv_obj ~ grupo + enf_cardiaca + edad, data = treatment_coh)

summary(cox_multivariado2)


# surv_obj

surv_obj_trans <- Surv(treat_trans$tiempo, treat_trans$evento)

# Quito los desenlaces porque si no sesgan el resultado

treat_trans <- 
  treat_trans %>% 
  select(-c(evento, 
            tiempo))

# Modelo cox

summary(coxph(surv_obj_trans ~ ., data = treat_trans))

