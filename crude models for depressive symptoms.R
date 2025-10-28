
library(dplyr)
class(study_pop$EPDS_9)
# our question is: What is the effect of being assigned to a specific malaria medication group during pregnancy on maternal depression at 1 year postpartum?”
study_pop <- study_pop %>%
  mutate(
    EPDS_9 = case_when(
      EPDS_9 == 0 ~ "No",
      EPDS_9 == 1 ~ "Yes",
      TRUE ~ NA_character_  # handle any other/unexpected values
    ),
    EPDS_9 = factor(EPDS_9, levels = c("No", "Yes")),
    medgroup = factor(medgroup)
  )


study_pop <- study_pop %>%
  mutate(
    medgroup = relevel(factor(medgroup), ref = "SP")
  )

# Run logistic regression on depression with 9 as cut off
logit_model <- glm(EPDS_9 ~ medgroup, data = study_pop, family = binomial)

# View summary
summary(logit_model)

# Optional: exponentiate coefficients to get odds ratios
exp(coef(logit_model))

# ORs with confidence intervals
exp(confint(logit_model))  # might take a moment due to profiling



# Run logistic regression with 12 as cut off
logit_model <- glm(EPDS_12 ~ medgroup, data = study_pop, family = binomial)

# View summary
summary(logit_model)

# Optional: exponentiate coefficients to get odds ratios
exp(coef(logit_model))

# ORs with confidence intervals
exp(confint(logit_model))  # might take a moment due to profiling


# Run linear regression
lm_model <- lm(epds_total ~ medgroup, data = study_pop)

# View model summary
summary(lm_model)

# Calculate 95% confidence intervals
confint(lm_model, level = 0.95)




