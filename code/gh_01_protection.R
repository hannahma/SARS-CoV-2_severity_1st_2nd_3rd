# code for manuscript:
#    1st, 2nd, and 3rd+ SARS-CoV-2 infections: 
#    associations of prior infections with protection and severity
#
# this is the full code used for analysis
# but without the data (which has protected health information), it can't be run
# - if you are interested in running the code, please see:
#   gh_01_protection_simplified_w_fake_data.R


library(tidyverse) # data wrangling
library(lubridate) # work with dates
library(patchwork) # combine plots

file_date_recent_data


source("./code/gh_00_functions.R")


# model examples----
glm(pcr_sympt ~ inf_til_now_cat,
    family = poisson, data = dat_long_ms #%>% filter(week < ymd("2024-01-01"))
)  %>%
  tidy_irr()

glm(pcr_sympt ~ inf_til_now_cat + vax_status ,
    family = poisson, 
    data = dat_long_ms 
)  %>%
  tidy_irr()

# *** nest data (to run multiple models at once)----
dat_nested <- nest_by_agecat_week_period_x(dat_long_ms) # nest, not including levels of prior infection
dat_nested

dat_nested_no     <- dat_nested %>% filter(period != "original") # no = no original
dat_nested_subset <- dat_nested %>% filter(! (variable=="age and period" & period=="original"))

dat_nested

# MAIN EXPOSURE: # PRIOR INFECTIONS----

# *** plots of cumulative infection, vaccination----
p_cum_inf_adult <- dat_long_ms %>% filter(agecat_week=="18+y") %>% plot_cum_inf() + facet_agecat_week()
p_cum_vax_adult <- dat_long_ms %>% filter(agecat_week=="18+y") %>% plot_cum_vax() + facet_agecat_week()
p_cum_inf_adult / p_cum_vax_adult

p_cum_inf_age <- dat_long_ms %>% filter(agecat_week != "18+y") %>% plot_cum_inf() + facet_agecat_week()
p_cum_vax_age <- dat_long_ms %>% filter(agecat_week != "18+y") %>% plot_cum_vax() + facet_agecat_week()
p_cum_inf_age / p_cum_vax_age


# ~~~ INCIDENCE ~~~----
# *** long tables to plot----
inc <- dat_long_ms %>% 
  nest_by_agecat_week_period_x(inf_til_now_cat) %>% # nest, including levels of prior infection
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


inc_sum          <- inc %>% filter(variable=="summary")
inc_age          <- inc %>% filter(variable=="age") %>% filter(agecat_week != "18+y")
inc_adult        <- inc %>% filter(variable=="age" & agecat_week=="18+y")
inc_period       <- inc %>% filter(variable=="period" & period != "original")
inc_period_age   <- inc %>% filter(variable=="age and period" & period != "original" & agecat_week != "18+y")
inc_period_adult <- inc %>% filter(variable=="age and period" & period != "original" & agecat_week == "18+y")

# *** plot incience (overall & age subsets)----
p_inc_adult <- plot_inc_inf(inc_adult) + facet_grid(rows = vars(agecat_week))
p_inc_age   <- plot_inc_inf(inc_age)  

p_inc_adult
p_inc_age


# *** plot incidence (by period & age)----
p_inc_period_adult <- plot_inc_inf(inc_period_adult) + facet_agecat_week_period()
p_inc_period_adult

# ~~~ INCIDENCE RATIOS ~~~----
# ** run models----
irrs_inf_adj   <- add_pcr_sympt_modsev(dat_nested_no, mod_prot_inf_vax) %>% unnest_mod_no_vax()
irrs_inf_adj

# *** plot protection (age subsets)----
irrs_inf_adj_adult <- irrs_inf_adj %>% filter(variable=="age" & agecat_week == "18+y")
irrs_inf_adj_age   <- irrs_inf_adj %>% filter(variable=="age" & agecat_week != "18+y")

p_prot_inf_adult <- plot_prot_inf_adj(irrs_inf_adj_adult) + facet_agecat_week()
p_prot_inf_age   <- plot_prot_inf_adj(irrs_inf_adj_age)   + facet_agecat_week()

p_prot_inf_adult
p_prot_inf_age

# *** plot protection by period (age subsets)----
irrs_inf_adj_period_adult <- irrs_inf_adj %>% filter(variable=="age and period" & agecat_week == "18+y")

p_prot_inf_period_adult <- plot_prot_inf_adj(irrs_inf_adj_period_adult) + facet_agecat_week_period()
p_prot_inf_period_adult


# ~~~ COMBINE FIGURES FOR MS ~~~----
# * figures - protection (by age)----
fig_prot_adult <- p_cum_inf_adult / p_cum_vax_adult / (p_inc_adult | p_prot_inf_adult) & plot_annotation(tag_levels = "A")

fig_prot_age <-  (
  (p_cum_inf_age + labs(fill = "prior\ninfections")) / 
    (p_cum_vax_age + labs(fill = "prior\nvaccinations")) |
    p_inc_age | p_prot_inf_age 
) &
  theme(legend.position = "bottom") & 
  plot_annotation(tag_levels = "A") 

fig_prot_adult
fig_prot_age

fig_prot_age[[1]] <- fig_prot_age[[1]] &
  scale_x_date(labels = scales::label_date_short(),
               date_breaks = "6 months",
               limits = c(ymd("2020-01-01"),
                          time_period_end
               )
  )
fig_prot_age

# * figures - protection (by period)----
fig_prot_period_adult <- (p_inc_period_adult / p_prot_inf_period_adult) & 
  theme(legend.position = "right") & plot_annotation(tag_levels = "A")
fig_prot_period_adult


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
    subset = case_when(
      `# prior infections`==0 & variable=="summary" ~ age,
      `# prior infections`==0 & variable=="age"     ~ age,
      `# prior infections`==0 & variable=="period"  ~ period
    )
  ) %>% relocate(subset, .after = period)

tab_prot_2 <- tab_prot_2 %>% 
  mutate(
    subset = case_when(
      `# prior infections`==0 ~ paste0(age, ", "),
      `# prior infections`==1 ~ paste0(period)
    )
  ) %>% relocate(subset, .after = period) 


view(tab_prot_1)
view(tab_prot_2)


# ALT EXPOSURE: PERIOD OF PRIOR INFECTION----
dat_long_ms %>% count(period_last_inf, period_li, period)
dat_long_ms %>% count(period_li)

glm(pcr_sympt ~  period_li +  vax_status ,
    family = poisson, 
    data = dat_long_ms %>% 
      filter(agecat_week=="18+y") %>% 
      filter(period=="omicron") #%>% count(period_li)
)  %>%
  tidy_irr()

irrs_per_pi <- add_pcr_sympt_modsev(dat_nested_subset, mod_prot_period_li) %>% unnest_mod_no_vax() %>% 
  mutate(exposure = factor(exposure, levels = c("original","pre-omicron variant", "pre-omicron", "omicron")))

#irrs_per_pi_sum   <- irrs_per_pi %>% filter(agecat_week == "all")
irrs_per_pi_age   <- irrs_per_pi %>% filter(variable=="age and period")
irrs_per_pi_adult <- irrs_per_pi %>% filter(variable=="age and period" & agecat_week == "18+y")

# noted in paper (convert to protection):
irrs_per_pi_adult %>% 
  #filter(str_detect(period, "pre")) %>% 
  filter(outcome=="sympt") %>% 
  select(agecat_week, period, outcome, exposure, irr, l, u) %>% 
  mutate(across(c(irr, l, u), ~ (1 - .x)*100)) 

# plot
p_prot_per_pi_adult  <- plot_prot_per_pi(irrs_per_pi_adult) + facet_agecat_week_period()
p_prot_per_pi_adult

# ALT EXPOSURE: TIME SINCE PRIOR INFECTION----
irrs_time_pi <- add_pcr_sympt_modsev(dat_nested_subset, mod_prot_time_li) %>% unnest_mod_no_vax() %>% 
  mutate(exposure = factor(exposure, levels = c("0-180 days", "180-365 days", "1-2 years",">= 2 years")))

irrs_time_pi_adult <- irrs_time_pi %>% filter(variable=="age and period" & agecat_week == "18+y")

# noted in paper (convert to protection):
irrs_time_pi_adult %>% 
  #filter(str_detect(period, "pre")) %>% 
  filter(outcome %in% c("sympt", "mod/sev")) %>% 
  select(agecat_week, period, outcome, exposure, irr, l, u) %>% 
  mutate(across(c(irr, l, u), ~ (1 - .x)*100)) 

p_prot_time_adult  <- plot_prot_time_pi(irrs_time_pi_adult) + facet_agecat_week_period()


# save figures----
file_suffix <- paste0(file_date_recent_data, ".pdf")


# # prior infections
ggsave(fig_prot_adult,        filename = paste0("./figures/fig_1_prot_adult_",         file_suffix), width = 10, height = 8)
ggsave(fig_prot_period_adult, filename = paste0("./figures/fig_s3_prot_period_adult_", file_suffix), width = 10, height = 8)
ggsave(fig_prot_age,          filename = paste0("./figures/fig_s6_prot_age_",          file_suffix), width = 17, height = 10)

# save tables----
write_csv(tab_prot_1, file = paste0("./tables/unformatted/protection_1_", file_date_recent_data,".csv"))
write_csv(tab_prot_2, file = paste0("./tables/unformatted/protection_2_", file_date_recent_data,".csv"))


