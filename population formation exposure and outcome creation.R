#load necessary libraries
library(dplyr)


# Calculate the number of unique mothers (midm)
num_mothers_MiPPAD <- tabmev5_kb %>%
  summarise(unique_mothers = n_distinct(midm))

# Print the result
print(num_mothers_MiPPAD)
# we have 1183 mothers enrolled in the MiPPAD study

table(tabmev5_kb$v05del_delivoutc_7)

table(tabmev5_kb$v05del_deloutc_7oth) # here is the other abortion that was hiding

# Create a new variable `birth_outcome` using only the first variable
tabmev5_kb <- tabmev5_kb %>%
  mutate(birth_outcome = case_when(
    v05del_delivoutc_7 == "Live birt" ~ "Live Birth",      # Group Live Birth
    v05del_delivoutc_7 == "Stillbirt" ~ "Still Birth",     # Group Still Birth
    v05del_delivoutc_7 %in% c("Spontaneo", "Other") ~ "Abortion", # Group Spontaneous and Other as Abortion
    v05del_delivoutc_7 %in% c(".i", ".a") ~ NA_character_  # Treat missing values
  ))


# Check the summary of the new variable
table(tabmev5_kb$birth_outcome, useNA = "ifany")


# Filter for distinct mothers based on `midm` and keep the first occurrence
distinct_mothers <- tabmev5_kb %>%
  distinct(midm, .keep_all = TRUE)

# Summarize the new `birth_outcome` variable for distinct mothers
birth_outcome_summary <- distinct_mothers %>%
  group_by(birth_outcome) %>%
  summarise(count = n()) %>%
  ungroup()

# Print the summary
print(birth_outcome_summary)
# we have 10 abortions
# we have 42 stillbirths 
# 1059 live births
# 72 NAs

# cross reference birth outcome with no study completion to understand


# Create the table with NA values considered
table_with_na <- table(
  distinct_mothers$vxxstu_notcompletersn_5, 
  distinct_mothers$v05del_delivoutc_7, 
  useNA = "ifany"
)

# Display the table
print(table_with_na)

table(tabmev5_kb$vxxstu_other_5)


# 1059 live births
#10 abortions
#42 stillbirths
#26 migrations
#14  loss to follow up
# 30 consent withdrawn
# 1 protocol not withheld
# 1 false pregnancy
# 28 twin births (there were 4 twin births which did not live)


mothers_with_twins <- tabmev5_kb %>%
  filter(
    birth_outcome == "Live Birth",             # Only live births
    v05del_multbirth_10 == "Yes"               # Only multiple births (twins)
  ) %>%
  distinct(midm)                               # Unique mothers




summary(mothers_with_twins) # here we have 32 mothers with twins

table(tabmev5_kb$v05del_multbirth_10, tabmev5_kb$birth_outcome)

TOVI_notwins <- tabmev5_kb %>%
  filter(
    v05del_multbirth_10 != "Yes",        # Exclude twins/multiples
    birth_outcome == "Live Birth"        # Keep only live births
  )
# 1027


# Filter for distinct mothers based on `midm` and keep the first occurrence
TOVI <- TOVI_notwins %>%
  filter(!is.na(i_TOVI)) %>%           # Keep only rows where the mother participated in TOVI
  distinct(midm, .keep_all = TRUE)    # Keep only one row per mother (based on `midm`)




#there is a question asking if the mother filled out the questionnaire with another person
#or if she was on her own... could see if this is an interaction
TOV_epds1
TOV_epds2
TOV_epds3
TOV_epds4
TOV_epds5
TOV_epds6
TOV_epds7
TOV_epds8
TOV_epds9
TOV_epds10

# exposure variable: type of malaria medication
#-------------------------------------------------------------------------------
table(TOVI$v01med_iptpgroup_1, useNA = "always")

TOVI$medgroup=TOVI$v01med_iptpgroup_1

# outcome: maternal depression
#-------------------------------------------------------------------------------
epds_items <- c("TOV_epds1", "TOV_epds2", "TOV_epds3", "TOV_epds4", "TOV_epds5",
                "TOV_epds6", "TOV_epds7", "TOV_epds8", "TOV_epds9", "TOV_epds10")


# Create total EPDS score
TOVI <- TOVI %>%
  mutate(
    epds_total = TOV_epds1 + TOV_epds2 + TOV_epds3 +
      TOV_epds4 + TOV_epds5 +
      TOV_epds6 + TOV_epds7 + TOV_epds8 + TOV_epds9 + TOV_epds10
  )

summary(TOVI$epds_total)  # range from 0 to 30 with 466 missing

# creating binary cut off based on african context:
#-------------------------------------------------------------------------------
# can explore two recommended cut offs 9 or higher
# if we want extreme cases, we can try 12 or more

# based on this paper: https://doi.org/10.1371/journal.pone.0082521

TOVI <- TOVI %>%
  mutate(
    EPDS_9 = ifelse(epds_total >= 9, 1, 0),
    EPDS_12 = ifelse(epds_total >= 12, 1, 0)
  )

# adding in anxiety outcome
# EPDS-3A total 
table(TOVI$TOV_epds3)
table(TOVI$TOV_epds4)
table(TOVI$TOV_epds5)

TOVI$EPDS3a_total <- with(TOVI,
                               TOV_epds3 + TOV_epds4 + TOV_epds5
)

# Maternal anxiety flag (≥6 = anxiety)
TOVI$maternal_anxiety <- as.integer(TOVI$EPDS3a_total >= 6)

# Quick check
table(TOVI$EPDS3a_total)
table(TOVI$maternal_anxiety)
#-------------------------------------------------------------------------------
# only those with exposure and outcome in population

TOVI_complete <- TOVI %>%
  filter(!is.na(medgroup), !is.na(epds_total))


# our population size is 753 mothers

write.csv(TOVI_complete, "study_pop.csv", row.names = FALSE)

# missing population

library(dplyr)

missing_pop <- TOVI_notwins %>%
  anti_join(TOVI_complete, by = "midm")

nrow(missing_pop)

write.csv(missing_pop, "missingpop_maternal.csv", row.names = FALSE)


summary(study_pop$TOV_dmts) # gross motor development

summary(study_pop$TOV_pvts) # visual perception

summary(study_pop$TOV_mfts) # fine motor skills

summary(study_pop$TOV_cvts) # verbal comprehension

summary(study_pop$TOV_evts) # verbal expression

summary(study_pop$TOV_ssms) # MSEL standard score

library(dplyr)

muller_study <- study_pop %>%
  filter(
    !is.na(TOV_dmts),  # Gross motor development
    !is.na(TOV_pvts),  # Visual perception
    !is.na(TOV_mfts),  # Fine motor skills
    !is.na(TOV_cvts),  # Verbal comprehension
    !is.na(TOV_evts),  # Verbal expression
    !is.na(TOV_ssms)   # MSEL standard score
  ) # 740 children

write.csv(muller_study, "muller_study.csv")

missing_pop_muller <- TOVI_notwins %>%
  anti_join(muller_study, by = "midm")

nrow(missing_pop_muller)

write.csv(missing_pop_muller, "missing_pop_muller.csv", row.names = FALSE)
# maybe i can do a sensitivity analysis where we look at only the children where there was not a problem
# but then again maybe the problems are important to have
table(study_pop$TOV_obser1)
table(study_pop$TOV_obser2)


# HOME score (all 753 are included)
summary(study_pop$TOV_senstot) #HOME I. (Responsivity/sensitivity) sub-total
summary(study_pop$TOV_aceptot) #HOME II. (Acceptance) sub-total
summary(study_pop$TOV_orgtot) #HOME III. (Organization) sub-total
summary(study_pop$TOV_matot) #HOME IV. (Learning materials) sub-total
summary(study_pop$TOV_engtot) #HOME V. (Involvement) sub-total
summary(study_pop$TOV_hometot) #HOME scale, total score

table(study_pop$TOV_obser9)
table(study_pop$TOV_obser10)


