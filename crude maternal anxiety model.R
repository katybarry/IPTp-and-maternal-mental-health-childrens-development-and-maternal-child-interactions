# Ensure EPDS_9 is binary factor with Yes as the event level
library(dplyr)
class(study_pop$maternal_anxiety)
study_pop$maternal_anxiety=as.factor(study_pop$maternal_anxiety)
# our question is: What is the effect of being assigned to a specific malaria medication group during pregnancy on maternal depression at 1 year postpartum?”
study_pop <- study_pop %>%
  mutate(
    maternal_anxiety = case_when(
      maternal_anxiety == 0 ~ "No",
      maternal_anxiety == 1 ~ "Yes",
      TRUE ~ NA_character_  # handle any other/unexpected values
    ),
    maternal_anxiety = factor(maternal_anxiety, levels = c("No", "Yes")),
    medgroup = factor(medgroup)
  )


study_pop <- study_pop %>%
  mutate(
    medgroup = relevel(factor(medgroup), ref = "SP")
  )

# Run logistic regression on depression with 9 as cut off
logit_model <- glm(maternal_anxiety ~ medgroup, data = study_pop, family = binomial)

# View summary
summary(logit_model)

# Optional: exponentiate coefficients to get odds ratios
exp(coef(logit_model))

# ORs with confidence intervals
exp(confint(logit_model))  # might take a moment due to profiling


# Run linear regression
lm_model <- lm(EPDS3a_total ~ medgroup, data = study_pop)

# View model summary
summary(lm_model)

# Calculate 95% confidence intervals
confint(lm_model, level = 0.95)




