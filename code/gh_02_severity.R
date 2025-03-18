# code for manuscript:
#    1st, 2nd, and 3rd+ SARS-CoV-2 infections: 
#    associations of prior infections with protection and severity
#
# this is the full code used for analysis
# but without the data (which has protected health information), it can't be run
# - an example has been provided for the protection code, which is set up similarly:
#   gh_01_protection_simplified_w_fake_data.R

library(tidyverse)  # data wrangling
library(patchwork)  # combine plots
library(MMWRweek)   # epiweeks
library(broom)      # tidy model output
library(janitor)    # round_half_up function
 
source("./code/gh_01_protection.R")

# * subset data----
dat_sev_ms_adult <- dat_sev_ms %>% filter(agecat == "18+y")
dat_sev_ms_kid   <- dat_sev_ms %>% filter(agecat %in% c("0-9y","10-17y"))

# * nest data to run models----
dat_sev_nest <- dat_sev_ms %>% nest_by_agecat_period_x()
dat_sev_nest

# no original:
dat_sev_nest_no <- dat_sev_nest %>% filter(period != "original")



# MAIN EXPOSURE: NTH INFECTION----

# ~~ * RT-PCR-CONFIRMED INFECTIONS ~~~----
p_pcr_adult <- dat_sev_ms_adult %>% plot_rtpcr_inf(19) + labs(subtitle = "ages 18+y")
p_pcr_kid   <- dat_sev_ms_kid   %>% plot_rtpcr_inf(20) + labs(subtitle = "ages 0-17y")
p_pcr_adult 
p_pcr_kid


# ~~ * PERCENT SEVERITY ~~~----
sev_adult <- sev_prior_inf(dat_sev_ms_adult, "out", agecat)
sev_age   <- sev_prior_inf(dat_sev_ms_kid,    "in",  agecat) 

sev_age[[1]] <- sev_age[[1]] + facet_grid(cols = vars(agecat)) # switch to row facets
sev_adult
sev_age

# by period
sev_per_adult <- sev_prior_inf(dat_sev_ms_adult, "in", agecat, period)
sev_per_age   <- sev_prior_inf(dat_sev_ms_kid,   "in", agecat, period)
sev_per_adult
sev_per_age



# **** models examples ----
mod <- glm(severe ~  nth_infection_cat + vax_status, 
           family = poisson,
           data = dat_sev_ms 
)
tidy_rr(mod)


dat_sev_ms %>% count(vax_status, vax_bf_inf)


# ~~~ RELATIVE SEVERITIES ~~~----
rs_nth_inf_adj   <- dat_sev_nest_no %>% add_severity_mod(mod_sev_inf_vax) %>% format_model_results() 
rs_nth_inf_adj



# *** plot relative severities----
rs_nth_inf_adj_adult <- rs_nth_inf_adj %>% filter(variable=="age" & agecat == "18+y")
rs_nth_inf_adj_age   <- rs_nth_inf_adj %>% filter(variable=="age" & agecat != "18+y")

p_rs_nth_adult <- plot_rel_sev_nth_adj(rs_nth_inf_adj_adult) + facet_agecat()
p_rs_nth_age   <- plot_rel_sev_nth_adj(rs_nth_inf_adj_age)   + facet_grid(cols = vars(agecat))

p_rs_nth_adult
p_rs_nth_age



# *** plot relative severities, by period----
rs_nth_inf_adj_per_adult <- rs_nth_inf_adj %>% filter(variable=="age and period" & agecat == "18+y")
p_rs_inf_per_adult <- plot_rel_sev_nth_adj(rs_nth_inf_adj_per_adult) + facet_agecat_period()

p_rs_inf_per_adult


# ~~* COMBINE FIGURES FOR MS ~~~----
fig_sev_adult <- combine_sev_plots(p_pcr_adult, sev_adult[[1]], p_rs_nth_adult)
fig_sev_age   <- ((p_pcr_kid + theme(legend.position = "none")) | sev_age[[1]] / p_rs_nth_age) + plot_annotation(tag_levels = "A")
fig_sev_adult
fig_sev_age

fig_sev_per_adult <- (sev_per_adult[[1]] / p_rs_inf_per_adult) + plot_annotation(tag_levels = "A") & theme(legend.position = "right") 
fig_sev_per_adult

# create severity table ----
dat_sev_nest_nth_inf <- dat_sev_ms %>% nest_by_agecat_period_x(nth_infection_cat)

t_sev <- dat_sev_nest_nth_inf %>% get_ns_sev() %>% rename(infection_n = nth_infection_cat)


tab_sev <- t_sev %>% 
  left_join(rs_nth_inf_adj %>% widen_model_sev()) %>% 
  mutate(
    across(starts_with("rr_"), ~ case_when(
      infection_n=="1st" ~ "ref",
      is.na(.x)          ~ "not calculated",
      TRUE               ~ .x
    ))
  ) %>% 
  arrange(variable, agecat, period, infection_n)

tab_sev %>% print(n=50)

tab_sev_1 <- tab_sev %>% filter(variable %in% c("summary", "age", "period"))
tab_sev_2 <- tab_sev %>% filter(variable %in% c("age and period"))


# format variable for printed display
tab_sev_1 <- tab_sev_1 %>% 
  mutate(
    subset = case_when(
      infection_n=="1st" & variable=="summary" ~ agecat,
      infection_n=="1st" & variable=="age"     ~ agecat,
      infection_n=="1st" & variable=="period"  ~ period
    )
  ) %>% relocate(subset, .after = period)

tab_sev_2 <- tab_sev_2 %>% 
  mutate(
    subset = case_when(
      infection_n=="1st" ~ paste0(agecat, ", "),
      infection_n=="2nd" ~ paste0(period)
    )
  ) %>% relocate(subset, .after = period) 


view(tab_sev_1)
view(tab_sev_2)



# ALT EXPOSURE: PERIOD OF PRIOR INFECTION----
glm(subclinical ~  period_li , 
    family = poisson,
    data = dat_sev_ms %>% 
      filter(period=="pre-omicron variant") %>% filter(agecat=="18+y")
) %>% tidy_rr()

# ~~~ RELATIVE SEVERITIES ~~~----
rs_per_pi <- dat_sev_nest_no %>% add_severity_mod(mod_sev_per_li) %>% format_model_results() %>% 
  mutate(
    label_level = factor(label_level, 
                         levels = c("original",
                                    "pre-omicron variant",
                                    "pre-omicron",
                                    "omicron")),
  ) %>% arrange(variable, agecat, outcome, label_level)
rs_per_pi

rs_per_pi %>% print(n=50)

# *** plot relative severities----
rs_per_pi
rs_per_pi %>% count(variable, agecat, period)

rs_per_pi_adult <- rs_per_pi %>% filter(variable=="age and period" & agecat=="18+y")
p_rs_per_pi_adult <- plot_rel_sev_per_pi(rs_per_pi_adult)
p_rs_per_pi_adult


# ALT EXPOSURE: TIME SINCE PRIOR INFECTION----
glm(subclinical ~  days_since_last_inf_cat + vax_status , 
    family = poisson,
    data = dat_sev_ms 
) %>% tidy_rr()

# ~~~ RELATIVE SEVERITIES ~~~----
rs_time_pi <- dat_sev_nest_no %>% 
  add_severity_mod(mod_sev_time_li) %>% format_model_results() %>% 
  mutate(
    label_level = factor(label_level, 
                         levels = c("0-180 days",
                                    "180-365 days",
                                    "1-2 years",
                                    ">= 2 years"))
  ) %>% 
 arrange(variable, agecat, outcome, label_level) 
rs_time_pi

rs_time_pi %>% count(label_level)

rs_time_pi %>% print(n=50)

# *** plot relative severities----
rs_time_pi_adult <- rs_time_pi %>% filter(variable=="age and period" & agecat=="18+y")
p_rs_time_pi_adult <- plot_rel_sev_time_pi(rs_time_pi_adult)
p_rs_time_pi_adult


# combine protection & sev for supplement----
p_per <- (p_prot_per_pi_adult + guides(fill = "none", shape = "none")) / 
  p_rs_per_pi_adult +
  plot_layout(guides = "collect") &
  theme(legend.position = "top") &
  labs(title = NULL)

p_per 

p_prot_inf_adult

p_time <- (p_prot_time_adult + guides(fill = "none", shape = "none")) /
  p_rs_time_pi_adult + 
  plot_layout(guides = "collect") &
  theme(legend.position = "top") &
  labs(title = NULL)
p_time

p_per_time <- (p_per | p_time ) & plot_annotation(tag_levels = list(c("A","","B")))
p_per_time



# save figures----
as_tibble(ls(pattern = "fig")) 

file_suffix <- paste0(file_date_recent_data, ".pdf")
file_suffix

# # prior infections
ggsave(fig_sev_adult,     filename = paste0("./figures/fig_2_sev_adult_",           file_suffix), width = 10, height = 8)
ggsave(fig_sev_per_adult, filename = paste0("./figures/fig_s4_sev_period_adult_",   file_suffix), width = 10, height = 8)
# other exposures
ggsave(p_per_time,        filename = paste0("./figures/fig_s5_other_exposures_",    file_suffix), width = 18, height = 10)
# kids
ggsave(fig_sev_age,       filename = paste0("./figures/fig_s7_sev_age_",            file_suffix), width = 17, height = 10)


# save tables----
write_csv(tab_sev_1, file = paste0("./tables/unformatted/severity_1_", file_date_recent_data,".csv"))
write_csv(tab_sev_2, file = paste0("./tables/unformatted/severity_2_", file_date_recent_data,".csv"))







