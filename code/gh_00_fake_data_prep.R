# code for manuscript:
#    1st, 2nd, and 3rd+ SARS-CoV-2 infections: 
#    associations of prior infections with protection and severity

# this code creates a fake example dataset that has the same structure
# and variable names as used in the analysis


# 1) create a function to subset data by year (to show loss to follow-up) and make a dataset of participants----
sample_for_year <- function(data, group_to_lower, n) {
  data %>%
    filter(year %in% group_to_lower) %>%
    group_by(year) %>%
    slice_sample(n = n) %>%
    ungroup %>%
    bind_rows(data %>% filter(!year %in% group_to_lower))
}

weeks <- seq.Date(ymd("2020-01-05"),
                  ymd("2024-10-13"),
                  by = "week")
participants <- 2500

tot <- participants*length(weeks)



dat <- tibble(
  week = weeks,
  year = year(week)
) %>% 
  crossing(codigo = c(1:participants)) %>%  # number of participants
  # subset to show some loss to follow-up
  sample_for_year(2021, floor(tot/5*.95)) %>% 
  sample_for_year(2022, floor(tot/5*.95*.95)) %>% 
  sample_for_year(2023, floor(tot/5*.95*.95*.90)) %>% 
  sample_for_year(2024, floor(tot/5*.95*.95*.90*.80)) %>% 
  arrange(week, codigo) 

dat %>% 
  group_by(year) %>% 
  summarise(year_n = n())




# 2) define number of previous infections / vaccinations each year----
inf_til_now_cat <-   
c(
  sample(c("0","1","2"),      130000, prob = c(0.67, 0.32, 0.01),        replace = T), # 2020
  sample(c("0","1","2","3+"), 118750, prob = c(0.18, 0.75, 0.06, 0.004), replace = T), # 2021
  sample(c("0","1","2","3+"), 112812, prob = c(0.04, 0.65, 0.27, 0.04),  replace = T), # 2022
  sample(c("0","1","2","3+"), 101531, prob = c(0.01, 0.45, 0.40, 0.13),  replace = T), # 2023
  sample(c("0","1","2","3+"),  81225, prob = c(0.02, 0.39, 0.40, 0.18),  replace = T)  # 2024
)

vax_til_now_num_cat <-   
  c(
    sample(c("0"),                  130000, prob = c(1),                             replace = T), # 2020
    sample(c("0","1","2","3","4+"), 118750, prob = c(0.86, 0.08, 0.05, 0.01, 0.00),  replace = T), # 2021
    sample(c("0","1","2","3","4+"), 112812, prob = c(0.27, 0.17, 0.24, 0.22, 0.10),  replace = T), # 2022
    sample(c("0","1","2","3","4+"), 101531, prob = c(0.25, 0.12, 0.20, 0.23, 0.20),  replace = T), # 2023
    sample(c("0","1","2","3","4+"),  81225, prob = c(0.25, 0.11, 0.19, 0.23, 0.21),  replace = T)  # 2024
  )


# 3) sample some number of infections / symptomatic infections / moderate/severe infections per year----
# * pcr----
set.seed(1234)
pcr <- rbind(
  tibble(week = as.Date(rnorm(121, mean = 0, sd =  58), origin = "2020-05-24")),
  tibble(week = as.Date(rnorm(445, mean = 0, sd =  54), origin = "2021-08-22")),
  tibble(week = as.Date(rnorm(442, mean = 0, sd = 107), origin = "2022-05-08")),
  tibble(week = as.Date(rnorm(105, mean = 0, sd =  95), origin = "2023-06-25")),
  tibble(week = as.Date(rnorm( 31, mean = 0, sd =  19), origin = "2024-06-23"))
) %>% 
  mutate(week = floor_date(week, unit = "week"),
         year = year(week),
         pcr = 1) %>% 
  filter(year > 2019) %>% 
  arrange(week)

pcr %>% 
  group_by(year) %>% 
  count()

participants_w_pcr <- dat %>% filter(year==2020) %>% slice_sample(n = 121, replace = F) %>% 
  rbind(dat %>% filter(year==2021) %>% slice_sample(n = 492, replace = F)) %>% 
  rbind(dat %>% filter(year==2022) %>% slice_sample(n = 391, replace = F)) %>% 
  rbind(dat %>% filter(year==2023) %>% slice_sample(n = 108, replace = F)) %>% 
  rbind(dat %>% filter(year==2024) %>% slice_sample(n = 32, replace = F)) %>%
  arrange(week, codigo)
participants_w_pcr

dat_pcr <- pcr %>% 
  add_column(participants_w_pcr %>% select(codigo))

# * sympt pcr----
set.seed(5678)
pcr_sympt <- rbind(
  tibble(week = as.Date(rnorm(112, mean = 0, sd =  59), origin = "2020-05-17")),
  tibble(week = as.Date(rnorm(321, mean = 0, sd =  56), origin = "2021-08-22")),
  tibble(week = as.Date(rnorm(308, mean = 0, sd = 104), origin = "2022-05-22")),
  tibble(week = as.Date(rnorm( 76, mean = 0, sd =  97), origin = "2023-06-25")),
  tibble(week = as.Date(rnorm( 26, mean = 0, sd =  16), origin = "2024-06-23"))
) %>% 
  mutate(week = floor_date(week, unit = "week"),
         year = year(week),
         pcr_sympt = 1) %>% 
  filter(year > 2019) %>% 
  arrange(week)

pcr_sympt %>% 
  group_by(year) %>% 
  count()

participants_w_pcr_sympt <- participants_w_pcr %>% filter(year==2020) %>% slice_sample(n = 111, replace = F) %>% 
  rbind(participants_w_pcr %>% filter(year==2021) %>% slice_sample(n = 350, replace = F)) %>% 
  rbind(participants_w_pcr %>% filter(year==2022) %>% slice_sample(n = 279, replace = F)) %>% 
  rbind(participants_w_pcr %>% filter(year==2023) %>% slice_sample(n =  73, replace = F)) %>% 
  rbind(participants_w_pcr %>% filter(year==2024) %>% slice_sample(n =  29, replace = F)) %>%
  arrange(week, codigo)
participants_w_pcr_sympt

dat_pcr_sympt <- pcr_sympt %>% 
  add_column(participants_w_pcr_sympt %>% select(codigo))


# * pcr_mod_sev----
set.seed(9013)
pcr_mod_sev <- rbind(
  tibble(week = as.Date(rnorm( 29, mean = 0, sd =  67), origin = "2020-04-14")),
  tibble(week = as.Date(rnorm(137, mean = 0, sd =  58), origin = "2021-08-15")),
  tibble(week = as.Date(rnorm( 84, mean = 0, sd = 109), origin = "2022-03-20")),
  tibble(week = as.Date(rnorm( 26, mean = 0, sd =  89), origin = "2023-06-18")),
  tibble(week = as.Date(rnorm(  4, mean = 0, sd =   6), origin = "2024-06-23"))
) %>% 
  mutate(week = floor_date(week, unit = "week"),
         year = year(week),
         pcr_mod_sev = 1) %>% 
  filter(year > 2019) %>% 
  arrange(week)

pcr_mod_sev %>%  
  group_by(year) %>% 
  count()

participants_w_pcr_mod_sev <- participants_w_pcr_sympt %>% filter(year==2020) %>% slice_sample(n = 26, replace = F) %>% 
  rbind(participants_w_pcr_sympt %>% filter(year==2021) %>% slice_sample(n =  155, replace = F)) %>% 
  rbind(participants_w_pcr_sympt %>% filter(year==2022) %>% slice_sample(n =   67, replace = F)) %>% 
  rbind(participants_w_pcr_sympt %>% filter(year==2023) %>% slice_sample(n =   25, replace = F)) %>% 
  rbind(participants_w_pcr_sympt %>% filter(year==2024) %>% slice_sample(n =    4, replace = F)) %>%
  arrange(week, codigo)
participants_w_pcr_mod_sev

dat_pcr_mod_sev <- pcr_mod_sev %>% 
  add_column(participants_w_pcr_mod_sev %>% select(codigo))



dat_long_ms <- dat %>% 
  add_column(inf_til_now_cat) %>% 
  add_column(vax_til_now_num_cat) %>% 
  mutate(agecat_week = "18+y") %>% 
  left_join(dat_pcr) %>% 
  left_join(dat_pcr_sympt,   relationship = "many-to-many") %>% 
  left_join(dat_pcr_mod_sev, relationship = "many-to-many") %>% 
  mutate(
    across(c(pcr, pcr_sympt, pcr_mod_sev), ~ if_else(is.na(.x), 0, .x)),
    
    vax_status = case_when(
      vax_til_now_num_cat=="0"             ~ "unvax",
      vax_til_now_num_cat=="1"             ~ "partial", # keeping simple for example, but some vaccines were 1 dose = full
      vax_til_now_num_cat=="2"             ~ "full",
      vax_til_now_num_cat %in% c("3","4+") ~ "boosted"
    ),
    
    period = case_when(
      week < ymd("2021-01-01") ~ "original",
      week < ymd("2022-01-01") ~ "pre-omicron variant",
      TRUE                     ~ "omicron"
    ),
    period = factor(period, levels = c("original","pre-omicron variant","omicron"))
  )
dat_long_ms
dat_long_ms %>% count(inf_til_now_cat)

