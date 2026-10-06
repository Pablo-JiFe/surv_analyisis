library(ggplot2)
library(dplyr)

set.seed(37)

treatment <- 
  data.frame(
  pacientes = paste0("p", seq(1, 200)),
  
  grupo = sample(c("Medicamento", "Sin medicamento"), 200, replace = TRUE, prob = c(0.5, 0.5)),
  
  tiempo = rnorm(200, mean = 80, sd = 20)
  
) %>% 
  mutate(
    evento = ifelse(grupo == "Medicamento",
                   sample(c(1, 0), 100, replace = TRUE, prob = c(0.3, 0.7)),
                   sample(c(1, 0), 100, replace = TRUE, prob = c(0.7, 0.3))
    )
  )


treatment %>% 
  slice_head(n = 50) %>% 
  ggplot(aes(x =  tiempo, y = pacientes, colour = factor(evento))) +
  geom_segment(aes(x = 0, y = pacientes, xend = tiempo), 
               color = "gray") +
  geom_point() +
  facet_wrap(~ grupo, scales = "free_y", ncol = 2) + 
  theme_classic(base_size = 30) + 
  scale_color_manual(values = c("blue", "red"))

