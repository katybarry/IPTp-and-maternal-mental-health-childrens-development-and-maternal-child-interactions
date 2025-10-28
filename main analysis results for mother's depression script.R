table(study_pop$medgroup)
class(study_pop$medgroup)
study_pop$medgroup=as.factor(study_pop$medgroup)
# Make sure medgroup is a factor and set SP as the reference
study_pop$medgroup <- relevel(factor(study_pop$medgroup), ref = "SP")

model_1 <- glm(EPDS_9 ~ medgroup, data = study_pop, family = "binomial")

# Get odds ratios
exp(coef(model_1))

# Get confidence intervals
confint_vals_1 <- confint(model_1)

# Exponentiate both coefficients and confidence intervals
ORs_1 <- exp(cbind(OR = coef(model_1), confint_vals_1))
ORs_1


model_2 <- glm(EPDS_12 ~ medgroup, data = study_pop, family = "binomial")

# Get odds ratios
exp(coef(model_2))

# Get confidence intervals
confint_vals_2 <- confint(model_2)

# Exponentiate both coefficients and confidence intervals
ORs_2 <- exp(cbind(OR = coef(model_2), confint_vals_2))
ORs_2
