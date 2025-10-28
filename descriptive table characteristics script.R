# load libraries
library(dplyr)

# List of your variables
vars <- c(
  "v01bas_birthday_1m",
  "v01bas_age_2",
  "v01bas_sex_3",
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
  "medgroup"
)

# Check the class of each variable
sapply(study_pop[vars], class)


study_pop$v01bas_canread_4=as.factor(study_pop$v01bas_canread_4)
study_pop$v01bas_canwrite_5=as.factor(study_pop$v01bas_canwrite_5)
study_pop$v01bas_abort_7=as.factor(study_pop$v01bas_abort_7)
study_pop$v01bas_livingchd_11 =as.numeric(study_pop$v01bas_livingchd_11)
study_pop$v01bas_anydisease_13=as.factor(study_pop$v01bas_anydisease_13)
study_pop$v01bas_diseasename_13yes=as.factor(study_pop$v01bas_diseasename_13yes)
study_pop$medgroup=as.factor(study_pop$medgroup)
study_pop$maternal_anxiety=as.factor(study_pop$maternal_anxiety)




study_pop <- study_pop %>%
  mutate(
    EPDS_9 = factor(EPDS_9, levels = c(0, 1), labels = c("No", "Yes")),
    EPDS_12 = factor(EPDS_12, levels = c(0, 1), labels = c("No", "Yes"))
  )

table(study_pop$v01bas_canread_4, useNA = "always")
table(study_pop$v01bas_canwrite_5, useNA = "always")
table(study_pop$v01bas_abort_7, useNA = "always")
summary(study_pop$v01bas_livingchd_11, useNA = "always")
table(study_pop$v01bas_anydisease_13, useNA = "always")
table(study_pop$v01bas_diseasename_13yes, useNA = "always")


class(study_pop$v01bas_livingchd_11)
study_pop$v01bas_livingchd_11 <- as.numeric(study_pop$v01bas_livingchd_11)
table(study_pop$v01bas_livingchd_11)
class(study_pop$v01bas_prevgest_6)
table(study_pop$v01bas_prevgest_6)

table(study_pop$EPDS_9, useNA = "always")

study_pop <- study_pop %>%
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


table(study_pop$v01bas_prevgest_6_cat, useNA = "always")


study_pop <- study_pop %>%
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

unique(study_pop$medgroup)                 # see exact raw values
dplyr::count(study_pop, medgroup, .drop=FALSE)  # keep zero-count factor levels
sum(is.na(study_pop$medgroup))             # any NAs introduced?


study_pop <- study_pop %>%
  mutate(medgroup = recode(medgroup,
                           "Mefloquine full dose" = "MQ full d",
                           "Mefloquine split dose" = "MQ split",
                           "Sulfadoxine-pryimethamine dose" = "SP"))

table(study_pop$medgroup)

study_pop_subset <- study_pop %>%
  select(
    medgroup,
    v01bas_age_2,
    v01bas_canread_4,
    v01bas_canwrite_5,
    v01bas_prevgest_6_cat,
    v01bas_livingchd_cat,
    v01bas_abort_7,
    v01bas_anydisease_13,
    v01bas_weight_14,
    v01bas_gestage_17,
    EPDS_9,
    EPDS_12,
    maternal_anxiety,
    TOV_ssms,
    TOV_hometot
  )
library(gtsummary)
# Create gtsummary table
table1 <- study_pop_subset %>%
  tbl_summary(
    label = list(
      v01bas_age_2 ~ "Mother’s age at inclusion",
      v01bas_canread_4 ~ "Mother can read (yes)",
      v01bas_canwrite_5 ~ "Mother can write (yes)",  # corrected label here
      v01bas_prevgest_6_cat ~ "Number of previous pregnancies",
      v01bas_livingchd_cat ~ "Number of living children",
      v01bas_abort_7 ~ "Had previous abortion (yes)",
      v01bas_anydisease_13 ~ "Has a chronic disease (yes)",
      v01bas_gestage_17 ~ "Gestational age at inclusion (in weeks)",
      v01bas_weight_14 ~ "Weight at inclusion (kg)",
      EPDS_9 ~ "Depressive symptoms one year post-birth (EPDS >=9)",
      EPDS_12 ~ "Depressive symptoms one year post-birth (EPDS3a >=12)",
      maternal_anxiety ~ "Anxiety symptoms at one year post-birth (EPDS3a >=6)",
      TOV_ssms ~ "Muller scale, mean (sd)",
      TOV_hometot ~ "HOME scale, mean (sd)"
    ),
    by = medgroup,
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
  save_as_docx(path = "Descriptive table.docx")

