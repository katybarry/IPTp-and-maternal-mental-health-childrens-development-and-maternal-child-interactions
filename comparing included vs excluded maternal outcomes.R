
library(dplyr)
# upload study pop and missingpop_maternal

#-------------------------------------------------------------------------------

# loss to follow up population

included_ids <- unique(study_pop$midm)


excluded_df <- missingpop_maternal[!missingpop_maternal$midm %in% included_ids, ]


excluded_unique_mothers <- excluded_df[!duplicated(excluded_df$midm), ]


# exposure variable: type of malaria medication
#-------------------------------------------------------------------------------

excluded_unique_mothers$medgroup=excluded_unique_mothers$v01med_iptpgroup_1
table(excluded_unique_mothers$v01med_iptpgroup_1, useNA = "always")

# set outcomes as missing
excluded_unique_mothers$EPDS_9 <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$EPDS_12 <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$TOV_orgtot  <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$TOV_senstot <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$TOV_aceptot <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$TOV_matot <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$TOV_engtot <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$TOV_hometot <- NA  # Add placeholder for missing outcome
excluded_unique_mothers$maternal_anxiety <- NA  # Add placeholder for missing outcome



# Compare baseline characteristics
# Step 1: Add the inclusion status to each dataset
study_pop$include <- "Yes"
excluded_unique_mothers$include <- "No"

library(tidyselect)

vars <- c("midm", "mide","include",
          "v01bas_birthday_1m",
          "v01bas_age_2",
          "TOV_sexe",
          "v01bas_canread_4",
          "v01bas_canwrite_5",
          "v01bas_prevgest_6",
          "v01bas_abort_7",
          "v01bas_livingchd_11",
          "v01bas_anydisease_13",
          "v01bas_diseasename_13yes",
          "v01bas_weight_14",
          "v01bas_gestage_17",
          "v01bas_height_15",
          "medgroup", "EPDS_9", "EPDS_12","maternal_anxiety", "TOV_orgtot", "TOV_senstot","TOV_aceptot","TOV_matot","TOV_engtot", "TOV_hometot")


study_pop_sub <- dplyr::select(study_pop, tidyselect::any_of(vars))
excluded_sub  <- dplyr::select(excluded_unique_mothers, tidyselect::any_of(vars))


compare <- rbind(study_pop_sub, excluded_sub)



compare$v01bas_canread_4=as.factor(compare$v01bas_canread_4)
compare$v01bas_canwrite_5=as.factor(compare$v01bas_canwrite_5)
compare$v01bas_abort_7=as.factor(compare$v01bas_abort_7)
compare$v01bas_livingchd_11 =as.numeric(compare$v01bas_livingchd_11)
compare$v01bas_anydisease_13=as.factor(compare$v01bas_anydisease_13)
compare$v01bas_diseasename_13yes=as.factor(compare$v01bas_diseasename_13yes)
compare$medgroup=as.factor(compare$medgroup)
table(compare$v01bas_canread_4, useNA = "always")
table(compare$v01bas_canwrite_5, useNA = "always")
table(compare$v01bas_abort_7, useNA = "always")
summary(compare$v01bas_livingchd_11, useNA = "always")
table(compare$v01bas_anydisease_13, useNA = "always")
table(compare$v01bas_diseasename_13yes, useNA = "always")


class(compare$v01bas_livingchd_11)
compare$v01bas_livingchd_11 <- as.numeric(compare$v01bas_livingchd_11)
table(compare$v01bas_livingchd_11)
class(compare$v01bas_prevgest_6)
table(compare$v01bas_prevgest_6)
compare <- compare %>%
  mutate(
    v01bas_prevgest_6_cat = case_when(
      v01bas_prevgest_6 == 0 ~ "0",
      v01bas_prevgest_6 == 1 ~ "1",
      v01bas_prevgest_6 == 2 ~ "2",
      v01bas_prevgest_6 == 3 ~ "3",
      v01bas_prevgest_6 >= 4 ~ "4 or more",
      TRUE ~ NA_character_  # catch any missing or odd values
    ),
    v01bas_prevgest_6_cat = factor(v01bas_prevgest_6_cat, levels = c("0", "1", "2", "3", "4 or more"))
  )


table(compare$v01bas_prevgest_6_cat, useNA = "always")


compare <- compare %>%
  mutate(
    v01bas_livingchd_cat = case_when(
      v01bas_livingchd_11 == 0 ~ "0",
      v01bas_livingchd_11 == 1 ~ "1",
      v01bas_livingchd_11 == 2 ~ "2",
      v01bas_livingchd_11 >= 3 ~ "3 or more",
      TRUE ~ NA_character_  # catch any missing or odd values
    ),
    v01bas_livingchd_cat = factor(v01bas_livingchd_cat, levels = c("0", "1", "2", "3 or more"))
  )


library(dplyr)
library(gtsummary)

compare <- compare %>%
  mutate(medgroup = recode(medgroup,
                           "MQ full d"= "Mefloquine full dose" ,
                           "MQ split"="Mefloquine split dose",
                           "SP"= "Sulfadoxine-pryimethamine dose" ))


table(compare$medgroup)




compare_subset <- compare %>%
  select(
    medgroup,
    include,
    v01bas_age_2,
    v01bas_canread_4,
    v01bas_canwrite_5,
    v01bas_prevgest_6_cat,
    v01bas_livingchd_cat,
    v01bas_abort_7,
    v01bas_anydisease_13,
    v01bas_weight_14,
    v01bas_gestage_17
  )

table1 <- compare_subset  %>%
  tbl_summary(
    label = list(
      medgroup~ "Treatment assignment",
      v01bas_age_2 ~ "Mother’s age at inclusion",
      v01bas_canread_4 ~ "Mother can read (yes)",
      v01bas_canwrite_5 ~ "Mother can write (yes)",  # corrected label here
      v01bas_prevgest_6_cat ~ "Number of previous pregnancies",
      v01bas_livingchd_cat ~ "Number of living children",
      v01bas_abort_7 ~ "Had previous abortion (yes)",
      v01bas_anydisease_13 ~ "Has a chronic disease (yes)",
      v01bas_gestage_17 ~ "Gestational age at inclusion (in weeks)",
      v01bas_weight_14 ~ "Weight at inclusion (kg)"
    ),
    by = include,
    missing = "ifany",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_p() %>%
  add_n() %>%
  bold_labels()

# Print the table
table1
# Load necessary libraries
library(gtsummary)
library(flextable)
library(officer)

table1 %>%
  as_flex_table() %>%
  save_as_docx(path = "Comparing included vs excluded population HOME and depression.docx")



write.csv(compare, "tipping_maternal_db.csv")


