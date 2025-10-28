TOV_dmts  # Gross motor development
TOV_pvts  # Visual perception
TOV_mfts  # Fine motor skills
TOV_cvts  # Verbal comprehension
TOV_evts  # Verbal expression
TOV_ssms   # MSEL standard score

# Make sure medgroup is a factor and SP is the reference
muller_study$medgroup <- relevel(factor(muller_study$medgroup), ref = "SP")
table(muller_study$medgroup)

model_std_score<- glm(TOV_ssms ~ medgroup, data = muller_study, family = "gaussian")
summary(model_std_score)
# Calculate 95% confidence intervals
confint(model_std_score, level = 0.95)
# Run each model separately with descriptive names
model_gross_motor <- glm(TOV_dmts ~ medgroup, data = muller_study, family = "gaussian")
summary(model_gross_motor)
confint(model_gross_motor, level = 0.95)
model_fine_motor      <- glm(TOV_mfts ~ medgroup, data = muller_study, family = "gaussian")
summary(model_fine_motor)
confint(model_fine_motor, level = 0.95)
model_visual_perc     <- glm(TOV_pvts ~ medgroup, data = muller_study, family = "gaussian")
summary(model_visual_perc)
confint(model_visual_perc, level = 0.95)

model_verbal_recept    <- glm(TOV_cvts ~ medgroup, data = muller_study, family = "gaussian")
summary(model_verbal_recept)
confint(model_verbal_recept, level = 0.95)
model_verbal_expr     <- glm(TOV_evts ~ medgroup, data = muller_study, family = "gaussian")
summary(model_verbal_expr)
confint(model_verbal_expr, level = 0.95)




