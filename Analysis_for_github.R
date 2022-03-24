#####################################
# Analysis: outcome hospitalization #
#####################################

# Add libraries

  library(tidyverse)
  library(broom)
  library(glue)
  library(splines)
  library(cri.utils)

# model
  
  hosp_poiss_2_fit <- model_poisson(formula = formula_1_2_full, .data = data_hosp_poiss_2_prep)
  
# crude model
  
  hosp_poiss_2_fit_crude <- model_poisson(formula = formula_1_2_crude, .data = data_hosp_poiss_2_prep)
  
# present statistics adjusted
  
  hosp_poiss_2_summ <- summ_poisson(hosp_poiss_2_fit)
  
# present statistics crude 
  
  hosp_poiss_2_summ_crude <- summ_poisson_crude(hosp_poiss_2_fit_crude)
  

#####
# sensitivity analysis including only recently enhanced infection-induced immunity
#####
  
# model
  
  hosp_poiss_secondary_fit <- model_poisson(formula = formula_1_2_full, .data = data_hosp_poiss_secondary_prep)
  
# crude model
  
  hosp_poiss_secondary_fit_crude <- model_poisson(formula = formula_1_2_crude, .data = data_hosp_poiss_secondary_prep)
  
# present statistics adjusted
  
  hosp_poiss_secondary_summ <- summ_poisson(hosp_poiss_secondary_fit)
  
# present statistics crude 
  
  hosp_poiss_secondary_summ_crude <- summ_poisson_crude(hosp_poiss_secondary_fit_crude)
  
  
# create formulas for adjusted and crude models

    formula_1_2_full <- as.formula(
      "sw_outcome ~ offset(count_days) + immunity_status + splines::ns(week_number, df = 3) + splines::ns(pn_age, df = 2) +
      pc_sex + splines::ns(pc_ses_comb20, df = 2) + pc_sector_eng + residency_type_eng + splines::ns(covid_burden, df = 3) + pn_cdc_risk_factors_7_bins +
      splines::ns(pn_gp_visits_1yr_by_age, df = 2) + pn_flu_vacc_5yr_binary"
    )
    
    formula_1_2_crude <- as.formula(
      "sw_outcome ~ offset(count_days) + immunity_status"
    )

# model Poisson

    model_poisson <- function(formula, .data) {
      glm(formula, family = "poisson", data = .data)
    }

# summarise main analysis (adjusted) using 1 - IRR as point

    summ_poisson <- function(model) {
      summ <- broom::tidy(model, parametric = T) %>% 
        mutate(
          point = (1 - exp(estimate)) * 100,
          upper = (1 - exp(estimate - 1.96 * std.error)) * 100,
          lower = (1 - exp(estimate + 1.96 * std.error)) * 100
        )
      
      presentation <- summ %>% mutate_if(is.numeric, round, 0) %>%
        mutate(adjusted_estimate = glue::glue("{point} ({lower}-{upper})")) %>% 
        select(term, adjusted_estimate)
      
      return(list(summ, presentation))
    }

# summarise main analysis (crude) using 1- IRR as point

    summ_poisson_crude <- function(model) {
      summ <- broom::tidy(model, parametric = T) %>% 
        mutate(
          point = (1 - exp(estimate)) * 100,
          upper = (1 - exp(estimate - 1.96 * std.error)) * 100,
          lower = (1 - exp(estimate + 1.96 * std.error)) * 100
        )
      
      presentation <- summ %>% mutate_if(is.numeric, round, 0) %>%
        mutate(crude_estimate = glue::glue("{point} ({lower}-{upper})")) %>% 
        select(term, crude_estimate)
      
      return(list(summ, presentation))
    }


