
#' Preprocess Vaccine Survey Wave Data
#'
#' @param survey_data Data frame.
#' @param output_path String (optional). Path to save the resulting .RData file.
#'
#' @import readr dplyr here readxl forcats
#'
#' @return A cleaned data frame.
#' @export
process_survey_wave <- function(survey_data, output_path = NULL) {

  # 2. Standardize Column Names
  data <- survey_data |>
    # Outcomes ---

    # wave 1

    # MMR Doses
    rename_with(
      .fn = ~ rep("Oldest child received MMR doses", length(.)),
      .cols = matches("For the following question, by 'MMR vaccine' we mean the measles")
    ) |>

    # Flu Last Winter
    rename_with(
      .fn = ~ rep("Received flu vaccine last winter", length(.)),
      .cols = matches("Did you get the flu vaccine last winter")
    ) |>

    # Flu Next Year
    rename_with(
      .fn = ~ rep("Flu vaccine likelihood in next year", length(.)),
      .cols = matches("How likely, if at all, are you to get the flu vaccine this year?")
    ) |>

    # waves 1 & 2
    rename_with(~ "concerned_raw",
                matches("I am concerned about serious adverse effects")) |>
    # "Vaccines are effective"
    # "Childhood vaccines are important for a child's health"

    # Demographics (Wave 2 specific names)
    rename_with(~ "Gender", matches("^Are you...\\?"))

  if (("Ethnicity...3" %in% names(data)) && (!"Ethnicity" %in% names(data))) {
    data <- data |> rename_with(~ "Ethnicity", matches("^Ethnicity\\.\\.\\.3"))
  }

  # 3. Handle 'Concerned' Column Consistency
  # If Wave 1 name exists, rename it to 'concerned_raw' to match the Wave 2 logic above
  if ("I am concerned about serious adverse effects (i.e. negative side effects or unwanted results) of vaccines" %in% names(data)) {
    data <- data |>
      rename(concerned_raw = `I am concerned about serious adverse effects (i.e. negative side effects or unwanted results) of vaccines`)
  }

  likert_order <- c("Strongly disagree", "Tend to disagree",
                    "Don't know", "Neither agree nor disagree",
                    "Tend to agree", "Strongly agree")

  neg_hear_orig_colname <- "Not applicable - I've only heard negative comments about one or more of these campaigns"

  required_cols <- c(
    "Oldest child received MMR doses",
    "Received flu vaccine last winter",
    "Flu vaccine likelihood in next year",
    "concerned_raw",
    "Vaccines are effective",
    "Childhood vaccines are important for a child's health",
    neg_hear_orig_colname
  )

  for (col in required_cols) {
    if (!col %in% names(data)) {
      data[[col]] <- NA
    }
  }

  # 4. Main Cleaning Pipeline
  data <- data |>
    mutate(
      # Outcomes
      `Oldest child received MMR doses` = coalesce(`Oldest child received MMR doses`, "Not parent"),

      vacc_last_yr = as.factor(`Received flu vaccine last winter`),    # NLV_Q2
      vacc_next_yr = as.factor(`Flu vaccine likelihood in next year`), # NLV_Q3
      mmr_doses    = as.factor(`Oldest child received MMR doses`),     # NLV_Q5

      # Likert Scales
      vacc_eff       = factor(`Vaccines are effective`, levels = likert_order),  # NLV_Q6_1
      concerned      = factor(concerned_raw, levels = likert_order),             # NLV_Q6_2
      vacc_important = factor(`Childhood vaccines are important for a child's health`, levels = likert_order),  # NNLV_Q6_3

      negative_hear = .data[[neg_hear_orig_colname]]  # NLV_Q8a
    ) |>

    # --- Coalesce Matrix Questions ---
    # Checks for "y" (Wave 1), "selected" (Both), and "Yes" (Wave 2)

    # Working Status
    mutate(across(any_of(c("Working full time", "Working part time", "Full time student", "Retired", "Unemployed", "Not working/ Other")),
                  ~ if_else(.x %in% c("y", "Yes", "selected"), cur_column(), NA_character_))) %>%
    mutate(working_status = coalesce(!!!select(., any_of(c("Working full time", "Working part time", "Full time student", "Retired", "Unemployed", "Not working/ Other")))),
           working_status = factor(working_status, levels = c("Working full time", "Working part time", "Full time student", "Not working/ Other", "Retired", "Unemployed"))) |>

    # Children in Household
    mutate(across(any_of(c("0", "1", "2", "3+", "Refused")),
                  ~ if_else(.x %in% c("y", "Yes", "selected"), cur_column(), NA_character_))) %>%
    mutate(children_in_household = coalesce(!!!select(., any_of(c("0", "1", "2", "3+", "Refused")))),
           children_in_household = as.factor(children_in_household)) |>

    # Parent/Guardian
    mutate(across(any_of(c("Not parent/ guardian", "4 years and under", "5 to 11 years", "12 to 16 years", "17 to 18 years", "18 years and under", "Over 18 years")),
                  ~ if_else(.x %in% c("y", "Yes", "selected"), cur_column(), NA_character_))) %>%
    mutate(parent_guardian = coalesce(!!!select(., any_of(c("Not parent/ guardian", "4 years and under", "5 to 11 years", "12 to 16 years", "17 to 18 years", "18 years and under", "Over 18 years")))),
           parent_guardian = factor(parent_guardian, levels = c("4 years and under", "5 to 11 years", "12 to 16 years", "17 to 18 years", "Over 18 years", "Not parent/ guardian"))) |>

    # --- Demographics ---
    mutate(
      Age = factor(Age, levels = c("18-24", "25-34", "35-44", "45-54", "55+")),
      Gender = factor(Gender, levels = c("Female", "Male")),

      # ABC1: more affluent or middle-class groups (Grades A, B, and C1)
      # C2DE: less affluent or working-class groups (Grades C2, D, and E)
      social_grade = ifelse(`Social Grade` == "ABC1", "High", "Low"),

      marital_status = ifelse(is.na(`Marital Status`), "blank", `Marital Status`),
      marital_status = factor(marital_status, levels = c("Married/ Civil Partnership", "Living as married", "Never Married", "Separated/ Divorced", "Widowed", "blank"))
    )

  # 5. Housing Tenure (Check existence)
  if ("House Tenure" %in% names(data)) {
    data <- data |>
      mutate(
        house_tenure = ifelse(is.na(`House Tenure`), "blank", `House Tenure`),
        house_tenure = iconv(house_tenure, from = "", to = "UTF-8", sub = ""),
        house_tenure = gsub(x = house_tenure, pattern = "[^[:print:]]", replacement = ""),
        house_tenure = gsub(x = house_tenure, pattern = "\\?\\? ", replacement = ""),
        house_tenure = case_match(
          house_tenure,
          "Neither I live rent-free with my parents, family or friends" ~ "Rent-free",
          "Neither I live with my parents, family or friends but pay some rent to them" ~ "Family/ friends but pay rent",
          "Own (part-own) through shared ownership scheme (i.e. pay part mortgage, part rent)" ~ "Shared ownership scheme",
          .default = house_tenure),
        house_tenure = factor(house_tenure, levels = c("Own outright", "Own with a mortgage", "Rent-free", "Rent from a housing association", "Rent from a private landlord", "Rent from my local authority", "Family/ friends but pay rent", "Shared ownership scheme", "Other", "blank"))
      )
  } else {
    warning("Variable 'House Tenure' not found in dataset. Creating NA column.")
    data$house_tenure <- NA
    data$house_tenure <- as.factor(data$house_tenure)
  }

  # 6. Ethnicity Handling (Internal Column Only)
  if ("Ethnicity" %in% names(data)) {
    data <- data |>
      mutate(ethnicity_grp = case_when(
        Ethnicity %in% c("English / Welsh / Scottish / Northern Irish / British", "Irish", "Gypsy or Irish Traveller", "Any other White background") ~ "White",
        Ethnicity %in% c("Indian", "Pakistani", "Bangladeshi", "Chinese", "Any other Asian background") ~ "Asian",
        Ethnicity %in% c("African", "Caribbean", "Any other Black / African / Caribbean background") ~ "Black",
        Ethnicity %in% c("White and Black Caribbean", "White and Black African", "White and Asian", "Any other Mixed / Multiple ethnic background") ~ "Mixed",
        Ethnicity %in% c("Arab", "Any other ethnic group") ~ "Other",
        Ethnicity == "Prefer not to say" ~ "Prefer not to say",
        TRUE ~ "Unknown"),
        ethnicity_grp = as.factor(ethnicity_grp),
        ethnicity_grp = relevel(ref = "White", ethnicity_grp))
  } else {
    warning("'Ethnicity' column not found. 'ethnicity_grp' will not be created.")
  }

  # 7. Global Variable Creation
  data <- data |>
    mutate(
      # collapsed Marital Status
      marital_collapsed = fct_collapse(marital_status,
                                       "Separated/Divorced/Widowed" = c("Separated/ Divorced", "Widowed")),
      # collapsed Ethnicity
      ethnicity_collapsed = if("ethnicity_grp" %in% names(data)) {
        fct_collapse(ethnicity_grp, "Mixed" = c("Mixed", "Other"))
      } else { NA },

      Age_Group_Binary = factor(case_when(
        Age %in% c("18-24", "25-34", "35-44") ~ "under 45",
        Age %in% c("45-54", "55+") ~ "45+",
        TRUE ~ as.character(Age)), levels = c("under 45", "45+")),

      working_status_grouped = factor(case_when(
        working_status %in% c("Working full time", "Working part time") ~ "Working",
        TRUE ~ as.character(working_status)), levels = c("Working", "Full time student", "Retired", "Unemployed", "Not working/ Other")),

      parent_guardian18 = factor(case_when(
        parent_guardian %in% c("4 years and under", "5 to 11 years", "12 to 16 years", "17 to 18 years") ~ "18 and under",
        TRUE ~ as.character(parent_guardian)), levels = c("18 and under", "Over 18 years", "Not parent/ guardian")),

      house_tenure_grouped = if("house_tenure" %in% names(data)) {
        case_when(
          house_tenure %in% c("Own outright", "Own with a mortgage", "Shared ownership scheme") ~ "Owns",
          house_tenure %in% c("Rent from a private landlord", "Rent from my local authority", "Rent from a housing association") ~ "Rents",
          house_tenure %in% c("Family/ friends but pay rent", "Rent-free") ~ "Friends / family",
          TRUE ~ as.character(house_tenure))
      } else { NA }
    )

  if (!is.null(output_path)) {
    save(data, file = output_path)
    message(paste("Saved processed data to:", output_path))
  }

  return(data)
}
