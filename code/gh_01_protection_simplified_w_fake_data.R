# code for manuscript:
#    1st, 2nd, and 3rd+ SARS-CoV-2 infections: 
#    associations of prior infections with protection and severity

library(tidyverse) # data wrangling
library(lubridate) # work with dates
library(patchwork) # combine plots

source("./code/gh_00_fake_data_prep.R") # load the fake data
source("./code/gh_00_functions.R")      # load the functions 


# model examples----
glm(pcr ~ inf_til_now_cat,
    family = poisson, data = dat_long_ms #%>% filter(week < ymd("2024-01-01"))
)  %>%
  tidy_irr()


# *** nest data (to run multiple models at once)----
dat_nested <- dat_long_ms %>% 
  group_by(agecat_week) %>% 
  nest() %>% 
  mutate(variable = "age",
         period   = "all") %>% 
  relocate(variable, agecat_week, period) %>% 
  ungroup() %>% 
  mutate(across(c(variable, agecat_week, period), as.factor))
dat_nested


# MAIN EXPOSURE: # PRIOR INFECTIONS----

# *** plots of cumulative infection, vaccination----
p_cum_inf_adult <- dat_long_ms %>% plot_cum_inf() + facet_agecat_week()
p_cum_vax_adult <- dat_long_ms %>% plot_cum_vax() + facet_agecat_week()
p_cum_inf_adult / p_cum_vax_adult


# ~~~ INCIDENCE ~~~----
# *** long tables to plot----
inc <- dat_long_ms %>% 
  group_by(agecat_week, inf_til_now_cat) %>% # group by agecat (and groupvar / inf_til_now_cat)
  nest() %>% 
  mutate(variable = "age", 
         period = "all") %>% 
  get_ns() %>% 
  select(variable, agecat_week, period, inf_til_now_cat, py, starts_with("n."), starts_with("inc.")) %>% 
  mutate(across(starts_with("inc"), as.numeric)) %>% 
  
  pivot_longer(cols          = c(starts_with("n"), starts_with("inc")),
               names_to      = c("measure","outcome"),
               names_pattern = "(.*)\\.(.*)"
  ) %>% 
  pivot_wider(names_from  = measure,
              values_from = value) %>% 
  mutate(
    outcome = factor(outcome, levels = c("any","sympt","mod_sev"),
                     labels = c("pcr", "sympt", "mod/sev"))
  ) %>%
  ungroup()
inc

inc %>% count(outcome)



inc_adult        <- inc %>% filter(variable=="age" & agecat_week=="18+y")

# *** plot incidence (overall & age subsets)----
p_inc_adult <- plot_inc_inf(inc_adult) + facet_grid(rows = vars(agecat_week))

p_inc_adult


# ~~~ INCIDENCE RATIOS ~~~----
# ** run models----
irrs_inf_crude <- add_pcr_sympt_modsev(dat_nested, mod_prot_inf)     %>% unnest_mod() 
irrs_inf_crude

irrs_inf_adj   <- add_pcr_sympt_modsev(dat_nested, mod_prot_inf_vax) %>% unnest_mod_no_vax()
irrs_inf_adj



# *** plot protection (age subsets)----
irrs_inf_adj_adult <- irrs_inf_adj %>% filter(variable=="age" & agecat_week == "18+y")

p_prot_inf_adult <- plot_prot_inf_adj(irrs_inf_adj_adult) + facet_agecat_week()
p_prot_inf_adult


# ~~~ COMBINE FIGURES FOR MS ~~~----
# * figures - protection (by age)----
fig_prot_adult <- p_cum_inf_adult / p_cum_vax_adult / (p_inc_adult | p_prot_inf_adult) & 
  plot_annotation(
    tag_levels = "A",
    # added for example:
    title = "Example Figure 1 for manuscript: \n1st, 2nd, and 3rd+ SARS-CoV-2 infections: associations of prior infections with protection and severity",
    subtitle = "(using made up, fake data)"
  )

fig_prot_adult

# * tables - protection----
tab_prot <- make_table_n(
  nest_by_agecat_week_period_x(dat_long_ms, inf_til_now_cat), 
  irrs_inf_adj
)
tab_prot %>% count(variable, age, period)

tab_prot_1 <- tab_prot %>% filter(variable %in% c("summary", "age", "period"))
tab_prot_2 <- tab_prot %>% filter(variable %in% c("age and period"))

# format variable for printed display
tab_prot_1 <- tab_prot_1 %>% 
  mutate(
    subset = case_when( # this is just how I want labels formated in the table in Word
      `# prior infections`==0 & variable=="summary" ~ age,
      `# prior infections`==0 & variable=="age"     ~ age,
      `# prior infections`==0 & variable=="period"  ~ period
    )
  ) %>% relocate(subset, .after = period)

tab_prot_2 <- tab_prot_2 %>% 
  mutate(
    subset = case_when( # this is just how I want labels formated in the table in Word
      `# prior infections`==0 ~ paste0(age, ", "),
      `# prior infections`==1 ~ paste0(period)
    )
  ) %>% relocate(subset, .after = period) 


view(tab_prot_1)
view(tab_prot_2)


# save figures----

# # prior infections
ggsave(fig_prot_adult,        filename = "./results/gh_example_fig_1_prot_adult.pdf",         width = 10, height = 8)

# save tables (these are then formatted more in Word)----
write_csv(tab_prot_1, file = "./results/gh_example_protection_1.csv")
write_csv(tab_prot_2, file = "./results/gh_example_protection_2.csv")


