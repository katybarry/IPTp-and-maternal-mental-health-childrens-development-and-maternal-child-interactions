# study population, exposure, outcome, and missing populations creation
#-------------------------------------------------------------------------------
library(dplyr)
# Filter for rows where the mother was reported dead
dead_mothers <- tabmev5_kb%>%
  filter(vxxstu_notcompletersn_5 == "Death") %>%
  select(midm, vxxstu_notcompletersn_5,v01med_iptpgroup_1,v01med_iptp1date_2m,v03ant_iptp2ddate_6m)

# Display the result
dead_mothers
# there is one mother coded as dead but she is not dead

# recoded this person who was said to be dead but they are not dead
tabmev5_kb <- tabmev5_kb %>%
  mutate(
    vxxstu_notcompletersn_5 = ifelse(
      midm == "" & vxxstu_notcompletersn_5 == "Death",
      "",  # or NA_character_ if you prefer it as missing
      vxxstu_notcompletersn_5
    )
  )

TOVI <- tabmev5_kb %>%
  filter(!is.na(TOV_cons)) %>%  # Keep only rows where the mother provided informed consent in TOVI
  distinct(midm, .keep_all = TRUE)  # keep only distinct mothers (this is ok because no twins are in TOVI anyway)
# 759 women provided consent for the TOVI follow up
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

# population for muller study
muller_study <- study_pop %>%
  filter(
    !is.na(TOV_dmts),  # Gross motor development
    !is.na(TOV_pvts),  # Visual perception
    !is.na(TOV_mfts),  # Fine motor skills
    !is.na(TOV_cvts),  # Verbal comprehension
    !is.na(TOV_evts),  # Verbal expression
    !is.na(TOV_ssms)   # MSEL standard score
  ) # 733 children

write.csv(muller_study, "muller_study.csv")

missing_pop_muller <- TOVI_notwins_relaxed %>%
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




# missing population will include those who were elligible for this study but went missing
#-------------------------------------------------------------------------------
missing_pop <- tabmev5_kb %>%
  anti_join(TOVI_complete, by = "midm")
nrow(missing_pop)
write.csv(missing_pop, "build_flowchart_missing.csv", row.names = FALSE)

# Filter for rows where the mother was reported dead
TOVI_test<- TOVI %>%
  filter(!is.na(vxxstu_notcompletersn_5)) %>%
  select(midm, vxxstu_notcompletersn_5,v01med_iptpgroup_1,v01med_iptp1date_2m,v03ant_iptp2ddate_6m,epds_total)
# 37 women who were marked as incomplete in MiPPAD were added back to TOVI


# here i can explore exactly why everyone is missing or not
tabmev5_kb <-  tabmev5_kb%>%
  mutate(birth_outcome = case_when(
    v05del_delivoutc_7 == "Live birt" ~ "Live Birth",      # Group Live Birth
    v05del_delivoutc_7 == "Stillbirt" ~ "Still Birth",     # Group Still Birth
    v05del_delivoutc_7 %in% c("Spontaneo", "Other") ~ "Abortion", # Group Spontaneous and Other as Abortion
    v05del_delivoutc_7 %in% c(".i", ".a") ~ NA_character_  # Treat missing values
  ))

tabmev5_kb <- tabmev5_kb %>%
  mutate(
    dose_one = ifelse(!is.na(v01med_iptp1date_2m), 1, 0),
    dose_two = ifelse(!is.na(v03ant_iptp2ddate_6m), 1, 0)
  )


epds_items <- c("TOV_epds1", "TOV_epds2", "TOV_epds3", "TOV_epds4", "TOV_epds5",
                "TOV_epds6", "TOV_epds7", "TOV_epds8", "TOV_epds9", "TOV_epds10")

# Create total EPDS score
tabmev5_kb <- tabmev5_kb %>%
  mutate(
    epds_total = TOV_epds1 + TOV_epds2 + TOV_epds3 +
      TOV_epds4 + TOV_epds5 +
      TOV_epds6 + TOV_epds7 + TOV_epds8 + TOV_epds9 + TOV_epds10
  )


tabmev5_kb <- tabmev5_kb %>% # Keep only rows where the mother provided informed consent in TOVI
  distinct(midm, .keep_all = TRUE) 

explore_it<- tabmev5_kb %>%
  select(midm,TOV_cons,birth_outcome, vxxstu_notcompletersn_5,v01med_iptpgroup_1,dose_one,dose_two,epds_total, vxxstu_other_5)

library(openxlsx)
write.xlsx(explore_it, "check_missing_all.xlsx", rowNames = FALSE)

# Now, i need to define the elligible people that were missing

### stopped here to be continued


# Find mothers with Live Birth in tabmev5_kb who are NOT in TOVI

library(dplyr)
table(explore_it$vxxstu_notcompletersn_5)

# 272 women who were elgibile but then had loss to follow up
livebirth_not_in_TOVI_blank <- explore_it %>%
  filter(
    birth_outcome == "Live Birth",
    is.na(TOV_cons),
    is.na(vxxstu_notcompletersn_5) | vxxstu_notcompletersn_5 == "Migration"| vxxstu_notcompletersn_5 == "Lost to f"| vxxstu_notcompletersn_5 == "Other"
  ) 


livebirth_not_in_TOVI_blank <- explore_it %>%
  mutate(
    meets_criteria = birth_outcome == "Live Birth" &
      is.na(TOV_cons) &
      (is.na(vxxstu_notcompletersn_5) |
         vxxstu_notcompletersn_5 %in% c("Migration", "Lost to f", "Other"))
  )
# View only those meeting criteria (without losing full dataset)
subset_livebirth_not_in_TOVI_blank <- livebirth_not_in_TOVI_blank %>%
  filter(meets_criteria)



# View relevant columns
livebirth_not_in_TOVI_blank %>%
  select(midm, birth_outcome, vxxstu_notcompletersn_5) %>%
  head()

# Step 1: Identify MiPPAD mothers who were NOT completers (had a non-NA reason)
mippad_noncompleters <- tabmev5_kb %>%
  filter(!is.na(vxxstu_notcompletersn_5))

# Step 2: Identify which of these reappeared in TOVI
rejoined_in_TOVI <- TOVI_complete %>%
  semi_join(mippad_noncompleters, by = "midm")

# Step 3: Remove those rejoined mothers from your previous dataset
livebirth_cleaned <- livebirth_not_in_TOVI_blank %>%
  anti_join(rejoined_in_TOVI, by = "midm")

# Step 4: Check result
nrow(livebirth_cleaned)

# those who went missing even though they could have been in the study (272)
library(openxlsx)
write.xlsx(subset_livebirth_not_in_TOVI_blank, "missingpop_maternal.xlsx", rowNames = FALSE)

# now make the subset for those in the Muller study
#-------------------------------------------------------------------------------

muller_study <- study_pop %>%
  filter(
    !is.na(TOV_dmts),  # Gross motor development
    !is.na(TOV_pvts),  # Visual perception
    !is.na(TOV_mfts),  # Fine motor skills
    !is.na(TOV_cvts),  # Verbal comprehension
    !is.na(TOV_evts),  # Verbal expression
    !is.na(TOV_ssms)   # MSEL standard score
  ) # 733 children

write.csv(muller_study, "muller_study.csv")

missing_pop_muller <- TOVI_notwins_relaxed %>%
  anti_join(muller_study, by = "midm")

nrow(missing_pop_muller)


# 1) Full data with a flag
muller_missing_pop <- tabmev5_kb %>%
  mutate(
    meets_criteria =
      birth_outcome == "Live Birth" &
      (is.na(TOV_cons) |                      # no TOV_cons
          (!is.na(TOV_cons) & is.na(TOV_ssms)|is.na(TOV_dmts)|is.na(TOV_pvts)|is.na(TOV_mfts)|is.na(TOV_cvts)|is.na(TOV_evts))   # has TOV_cons but missing TOV_ssms
      ) &(is.na(vxxstu_notcompletersn_5) |
          vxxstu_notcompletersn_5 %in% c("Migration", "Lost to f", "Other")
      )
  )

# 2) Subset with only those meeting the criteria
muller_missing_subset <- muller_missing_pop %>% 
  filter(meets_criteria)

write.csv(muller_missing_subset, "missing_muller_elligible.csv")



