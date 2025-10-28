# Set SP as reference level for medgroup
study_pop$medgroup <- relevel(factor(study_pop$medgroup), ref = "SP")

# Model 1: HOME I - Responsivity/sensitivity
model_resp <- glm(TOV_senstot ~ medgroup, data = study_pop, family = "gaussian")
summary(model_resp)
confint(model_resp, level = 0.95)
# Model 2: HOME II - Acceptance
model_accept <- glm(TOV_aceptot ~ medgroup, data = study_pop, family = "gaussian")
summary(model_accept)
confint(model_accept, level = 0.95)
# Model 3: HOME III - Organization
model_org <- glm(TOV_orgtot ~ medgroup, data = study_pop, family = "gaussian")
summary(model_org)
confint(model_org, level = 0.95)
# Model 4: HOME IV - Learning materials
model_mat <- glm(TOV_matot ~ medgroup, data = study_pop, family = "gaussian")
summary(model_mat)
confint(model_mat, level = 0.95)
# Model 5: HOME V - Involvement
model_eng <- glm(TOV_engtot ~ medgroup, data = study_pop, family = "gaussian")
summary(model_eng)
confint(model_eng, level = 0.95)
# Model 6: HOME Total Score
model_total <- glm(TOV_hometot ~ medgroup, data = study_pop, family = "gaussian")
summary(model_total)
# Calculate 95% confidence intervals
confint(model_total, level = 0.95)





