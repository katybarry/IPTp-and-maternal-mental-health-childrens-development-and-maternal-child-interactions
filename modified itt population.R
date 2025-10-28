#load necessary libraries
library(dplyr)



# Calculate the number of unique mothers (midm)
num_mothers_MiPPAD <- tabmev5_kb %>%
  summarise(unique_mothers = n_distinct(midm))

# Print the result
print(num_mothers_MiPPAD)
# we have 1183 mothers enrolled in the MiPPAD study



# Create the table with NA values considered
table_with_na <- table(
  distinct_mothers$vxxstu_notcompletersn_5, 
  distinct_mothers$v05del_delivoutc_7, 
  useNA = "ifany"
)

# Display the table
print(table_with_na)




# 1059 live births
#10 abortions
#42 stillbirths
#26 migrations
#14  loss to follow up
# 30 consent withdrawn
# 1 protocol not withheld
# 1 false pregnancy
# 28 twin births (there were 4 twin births which did not live)


# participated in the TOVI follow up
TOVI <- tabmev5_kb %>%
  filter(!is.na(i_TOVI)) 

table(TOVI$v05del_multbirth_10) # only one twin pair included in TOVI but then why is the mother not logged twice if so?
TOVI_kb_duplicates <- TOVI %>%
  filter(mide %in% (TOVI %>%
                      count(midm) %>%
                      filter(n > 1) %>%
                      pull(midm))) %>%
  arrange(midm)

# There is no twin included here

distinct_mothers_TOVI <- TOVI %>%
  distinct(midm, .keep_all = TRUE)

# Subset rows where v05del_multbirth_10 is "Yes" and extract the midm column
midm_yes <- TOVI$midm[TOVI$v05del_multbirth_10 == "Yes"]

# View the result
print(midm_yes)


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


#-------------------------------------------------------------------------------
# only those with exposure and outcome in population

TOVI_complete <- TOVI %>%
  filter(!is.na(medgroup), !is.na(epds_total))


# our population size is 753 mothers

write.csv(TOVI_complete, "study_pop.csv", row.names = FALSE)




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



# there is also the raven scale

summary(study_pop$TOV_seriea)
summary(study_pop$TOV_serieb)
summary(study_pop$TOV_seriec)
summary(study_pop$TOV_seried)
summary(study_pop$TOV_seriee)
summary(study_pop$TOV_ratot)




















