## ---- warning=FALSE, results='hide',message=FALSE-----------------------------
## library packages
library(gompertztrunc)         ## calculate mortality differentials under double-truncation  
library(tidyverse)             ## data manipulation and visualization  
library(data.table)            ## fast data manipulation
library(cowplot)               ## publication-ready themes for ggplot
library(socviz)                ## helper functions for data visualization (Kieran Healy)  
library(broom)                 ## "tidy" model output 


## load data 
sim_data <- sim_data             ## simulated 
bunmd_demo <- bunmd_demo         ## real 
numident_demo <- gompertztrunc::numident_demo   ## real 

## -----------------------------------------------------------------------------
## Look at simulated data 
head(sim_data)

## What years do we have mortality coverage? 
sim_data %>% 
  summarize(min(dyear), max(dyear)) 

## -----------------------------------------------------------------------------
## run gompertz_mle function  
## returns a list 
simulated_example <- gompertz_mle(formula = aod ~ temp + as.factor(sex) + as.factor(isSouth),
                                                 left_trunc = 1888,
                                                 right_trunc = 1905,
                                                 data = sim_data)

## -----------------------------------------------------------------------------
## 1. starting value for coefficients (from linear regression)
simulated_example$starting_values

## 2. optim fit object
simulated_example$optim_fit

## 2. check model convergence (0 == convergence)
simulated_example$optim_fit$convergence

## 3. Look at model results 
simulated_example$results

## -----------------------------------------------------------------------------
## true coefficient values (we know because we simulated them)
mycoefs <- c("temp" = +.2, "sex" = -.5, "isSouth" = +.6)

## compare
simulated_example$results %>%
  filter(!stringr::str_detect(parameter, "gompertz")) %>%
  mutate(true_coef = mycoefs) %>%
  select(parameter, coef, coef_lower, coef_upper, true_coef)

## -----------------------------------------------------------------------------
## translate hazard rates to difference in e65
convert_hazards_to_ex(simulated_example$results, age = 65, use_model_estimates = T) %>% 
  select(parameter, hr, hr_lower, hr_upper, e65, e65_lower, e65_upper)

## -----------------------------------------------------------------------------
## look at data 
head(bunmd_demo)

## how many people per country? 
bunmd_demo %>%
  count(bpl_string)

## ---- fig.width = 6, fig.height = 4-------------------------------------------
## distribution of deaths?
ggplot(data = bunmd_demo) + 
  geom_histogram(aes(x = death_age),
                 fill = "grey",
                 color = "black",
                 binwidth = 1) + 
  cowplot::theme_cowplot() + 
  labs(x = "Age of Death",
       y = "N") + 
  facet_wrap(~bpl_string)

## -----------------------------------------------------------------------------
## run linear model 
lm_bpl <- lm(death_age ~ bpl_string + as.factor(byear), data = bunmd_demo)

## extract coefficients from model 
lm_bpl_tidy <- tidy(lm_bpl) %>%
  filter(str_detect(term, "bpl_string")) %>%
  mutate(term = prefix_strip(term, "bpl_string"))

## rename variables 
lm_bpl_tidy <- lm_bpl_tidy %>%
  mutate(
    e65 = estimate,
    e65_lower = estimate - 1.96 * std.error,
    e65_upper = estimate + 1.96 * std.error
  ) %>%
  rename(country = term) %>%
  mutate(method = "Regression on Age of Death")

## -----------------------------------------------------------------------------
## run gompertztrunc
## set truncation bounds to 1988-2005 because we are using BUNMD 
gompertz_mle_results <- gompertz_mle(formula = death_age ~ bpl_string, 
                                    left_trunc = 1988,
                                    right_trunc = 2005,
                                    data = bunmd_demo)

## convert to e65
## use model estimates — but can also set other defaults for Gompertz M and b. 
mle_results <- convert_hazards_to_ex(gompertz_mle_results$results, use_model_estimates = T)

## tidy up results 
mle_results <- mle_results %>% 
  rename(country = parameter) %>%
  filter(str_detect(country, "bpl_string")) %>%
  mutate(country = prefix_strip(country, "bpl_string")) %>%
  mutate(method = "Gompertz Parametric Estimate")

## look at results 
mle_results

## ---- fig.width = 7.2, fig.height = 5-----------------------------------------
## combine results from both models 
bpl_results <- lm_bpl_tidy %>%
  bind_rows(mle_results)

## calculate adjustment factor (i.e., how much bigger are Gompertz MLE results)
adjustment_factor <- bpl_results %>% 
  select(country, method, e65) %>%
  pivot_wider(names_from = method, values_from = e65) %>%
  mutate(adjustment_factor = `Gompertz Parametric Estimate` / `Regression on Age of Death`) %>%
  summarize(adjustment_factor_mean = round(mean(adjustment_factor), 3)) %>%
  as.vector()

## plot results
bpl_results %>%
  bind_rows(mle_results) %>%
  ggplot(aes(y = reorder(country, e65), x = e65, xmin = e65_lower, xmax = e65_upper, color = method)) +
  geom_pointrange(position = position_dodge(width = 0.2), shape = 1) +
  cowplot::theme_cowplot(font_size = 12) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(
    x = "Estimate",
    title = "Foreign-Born Male Ages at Death, BUNMD 1905-1914",
    y = "",
    subtitle = paste0("Gompertz MLE estimates are ~", adjustment_factor, " times larger than regression on age of death")
  ) +
  scale_color_brewer(palette = "Set1") +
  annotate("text", label = "Native Born Whites", x = 0.1, y = 3, angle = 90, size = 3, color = "black")

## -----------------------------------------------------------------------------
## create the dataset
bunmd_1915_cohort <- bunmd_demo %>% 
  filter(byear == 1915, death_age >= 65) %>% 
  filter(bpl_string %in% c("Native Born White", "Mexico"))

## run gompertz_mle()
bpl_results_1915_cohort <- gompertz_mle(formula = death_age ~ bpl_string, 
                                   data = bunmd_1915_cohort,
                                   left_trunc = 1988, 
                                   right_trunc = 2005)

## look at results
bpl_results_1915_cohort$results

## -----------------------------------------------------------------------------
diagnostic_plot(object = bpl_results_1915_cohort, data = bunmd_1915_cohort,
                covar = "bpl_string", death_var = "death_age")

## -----------------------------------------------------------------------------
diagnostic_plot_hazard(object = bpl_results_1915_cohort, data = bunmd_1915_cohort,
                covar = "bpl_string", death_var = "death_age", xlim=c(65,95))

## -----------------------------------------------------------------------------
## load in file 
numident_demo <- numident_demo

## recode categorical education variable to continuous "years of education" 
numident_demo <- numident_demo %>% 
  mutate(educ_yrs = case_when(
    educd == "No schooling completed" ~ 0,
    educd == "Grade 1" ~ 1,
    educd == "Grade 2" ~ 2,
    educd == "Grade 3" ~ 3,
    educd == "Grade 4" ~ 4,
    educd == "Grade 5" ~ 5,
    educd == "Grade 6" ~ 6,
    educd == "Grade 7" ~ 7,
    educd == "Grade 8" ~ 8,
    educd == "Grade 9" ~ 9,
    educd == "Grade 10" ~ 10,
    educd == "Grade 11" ~ 11,
    educd == "Grade 12" ~ 12,
    educd == "Grade 12" ~ 12,
    educd == "1 year of college" ~ 13,
    educd == "2 years of college" ~ 14,
    educd == "3 years of college" ~ 15,
    educd == "4 years of college" ~ 16,
    educd == "5+ years of college" ~ 17
  ))

## restrict to men 
data_numident_men <- numident_demo %>% 
  filter(sex == "Male") %>% 
  filter(byear %in% 1910:1920 & death_age > 65)


## -----------------------------------------------------------------------------
## look at person-level weights 
head(data_numident_men$weight)

## run gompertz model with person weights
education_gradient <- gompertz_mle(formula = death_age ~ educ_yrs, 
                                   data = data_numident_men,
                                   weights = weight, ## specify person-level weights 
                                   left_trunc = 1988, 
                                   right_trunc = 2005)

## look at results 
education_gradient$results 

## translate to e65
mle_results_educ <- convert_hazards_to_ex(education_gradient$results, use_model_estimates = T, age = 65) %>% 
  mutate(method = "Parametric Gompertz MLE")

## look at results
mle_results_educ

## ---- fig.width = 7.2, fig.height = 5-----------------------------------------
## run linear model 
lm_bpl <- lm(death_age ~ educ_yrs + as.factor(byear), data = data_numident_men, weights = weight)

## extract coefficients from model 
lm_bpl_tidy <- tidy(lm_bpl) %>%
  filter(str_detect(term, "educ_yrs"))

## rename variables 
ols_results <- lm_bpl_tidy %>%
  mutate(
    e65 = estimate,
    e65_lower = estimate - 1.96 * std.error,
    e65_upper = estimate + 1.96 * std.error
  ) %>%
  rename(parameter = term) %>%
  mutate(method = "Regression on Age of Death")

## Plot results
education_plot <- ols_results %>%
  bind_rows(mle_results_educ) %>%
  mutate(parameter = "Education (Years) Regression Coefficient") %>% 
  ggplot(aes(x = method, y = e65, ymin = e65_lower, ymax = e65_upper)) +
  geom_pointrange(position = position_dodge(width = 0.2), shape = 1) +
  cowplot::theme_cowplot(font_size = 12) +
  theme(legend.position = "bottom", legend.title = element_blank()) +
  labs(
    x = "",
    title = "Association Education (Years) and Longevity",
    subtitle = "Men, CenSoc-Numident 1910-1920",
    y = ""
  ) +
  scale_color_brewer(palette = "Set1") + 
  ylim(0, 0.5) 

education_plot

