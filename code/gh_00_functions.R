# code for manuscript:
#    1st, 2nd, and 3rd+ SARS-CoV-2 infections: 
#    associations of prior infections with protection and severity
#
# functions used in protection & severity analyses
# functions allow to make multiple versions of the same analysis / figure / table,
#  and only have to change the code in 1 place

library(broom)      # tidy model output
library(janitor)    # round_half_up function

basesize <- 12

time_period_end   <- "2024-10-14" # date through which data is subset and plots are cut off 
time_period_start <- "2020-01-01"

options(pillar.min_title_chars = Inf) # don't cut off variable names in console

# general (used in severity & protection analyses)----
# colors are from: https://colorbrewer2.org/
colors_severity <- c( # severe first
  "#49006a", # severe
  "#ae017e", # moderate
  "#f768a1", # mild
  "#fde0dd", # subclinical
  "grey60"   # missing
)

colors_num <- c(
  #"#ffffcc",
  "grey90",
  "#a1dab4",
  "#41b6c4",
  "#2c7fb8",
  "#253494"
)

# find out the hex code for default colors:
gg_color_hue <- function(n){
  
  hues = seq(15, 375, length = n + 1)
  
  hcl(h = hues, l = 65, c = 100)[1:n] # hue, chroma, luminance
}

gg_color_hue(5)

theme_hm <- function(basesize){
  theme_bw(base_size = {{ basesize }}) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(size = 13, face = "bold"), 
          
          
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          
          strip.text.y = element_text(angle = 0),
          
          axis.title = element_text(face = "bold"),
          title      = element_text(face = "bold")
          
    )
}

theme_hm_lines <- function(basesize){
  theme_bw(base_size = {{ basesize }}) +
    theme(strip.background = element_blank(),
          strip.text       = element_text(size = 13, face = "bold"), 
          
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank(),
          
          strip.text.y = element_text(angle = 0),
          
          axis.title = element_text(face = "bold"),
          title      = element_text(face = "bold")
          
    )
}


add_shading_blood <- function(color){
  alpha <- 0.1
  fill <- color
  
  list(
    # * add shading for blood sampling----
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf,# annual
             fill = color,
             xmin = ymd("2020-03-01"), 
             xmax = ymd("2020-04-17")),
    
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf,
             fill = color,
             xmin = ymd("2020-10-18") - 3.5, 
             xmax = ymd("2020-11-29") + 3.5),
    
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf, # annual
             fill = color,
             xmin = ymd("2021-02-15"), 
             xmax = ymd("2021-03-20")),
    
    # these are approximate timing:
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf,
             fill = color,
             xmin = ymd("2021-10-15"), 
             xmax = ymd("2021-11-30")),
    
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf, # annual
             fill = color,
             xmin = ymd("2022-02-15"), 
             xmax = ymd("2022-04-15")),
    
    
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf, # annual
             fill = color,
             xmin = ymd("2023-02-15"), 
             xmax = ymd("2023-04-05")),
    
    
    annotate("rect", alpha = alpha, ymin = -Inf, ymax = Inf, # annual
             fill = color,
             xmin = ymd("2024-02-15"), 
             xmax = ymd("2024-04-05")),
    
    geom_vline(xintercept = c(ymd("2020-01-01"),    # add lines for each year
                              ymd("2021-01-01"),
                              ymd("2022-01-01"),
                              ymd("2023-01-01"),
                              ymd("2024-01-01"))) 
  )
}

modify_date_plots <- function(plot){
  list(
    theme_hm(12),
    scale_x_date(labels = scales::label_date_short(),
                 date_minor_breaks = "1 months", 
                 date_breaks       = "3 months",
                 limits = c(ymd("2020-01-01"),
                            #today()
                            time_period_end
                            )
                 )
  )
}

# nest datat----
# for protection:
nest_by_agecat_week_period_x <- function(data, group_var){
  data %>% 
    group_by(agecat_week, {{ group_var }}) %>% # group by agecat (and groupvar / inf_til_now_cat)
    nest() %>% 
    mutate(period = "all") %>% 
    
    
    full_join(
      data %>% 
        group_by(period, {{ group_var }}) %>% # group by period (and groupvar / inf_til_now_cat)
        nest() %>% 
        mutate(agecat_week = "all")
    ) %>% 
    
    
    full_join(
      data %>% 
        group_by( {{ group_var }} ) %>% # summary (and grouped by groupvar / inf_til_now_cat)
        nest() %>% 
        mutate(agecat_week = "all",
               period      = "all")
    ) %>% 
    
    full_join(
      data %>% 
        # group by agecat AND period (and optionally groupvar / inf_til_now_cat)
        group_by(agecat_week, period, {{ group_var }}) %>% 
        nest() %>% 
        mutate(variable = "age and period")
    ) %>% 
    
    mutate(
      variable = case_when(agecat_week == "all" & period == "all" ~ "summary",
                           agecat_week != "all" & period == "all" ~ "age",
                           agecat_week == "all" & period != "all" ~ "period",
                           agecat_week != "all" & period != "all" ~ "age and period"),
      variable = factor(variable, levels = c("summary","age","period", "age and period"))
    ) %>% 
    
    # factor levels aren't preserved when grouping/nesting
    mutate(agecat_week = factor(agecat_week, levels = c("all", "18+y", "10-17y", "0-9y")),
           period = factor(period, levels = c("all", "original", "pre-omicron variant", "omicron"))) %>% 
    
    arrange(variable, agecat_week, period) %>% 
    relocate(variable, agecat_week, period)
  
}

# for severity
nest_by_agecat_period_x <- function(data, group_var){
  data %>% 
    group_by(agecat, {{ group_var }}) %>% # group by agecat (and groupvar / nth_infection_cat)
    nest() %>% 
    mutate(period = "all") %>% 
    
    
    full_join(
      data %>% 
        group_by(period, {{ group_var }}) %>% # group by period (and groupvar / nth_infection_cat)
        nest() %>% 
        mutate(agecat = "all")
    ) %>% 
    
    
    full_join(
      data %>% 
        group_by( {{ group_var }} ) %>% # summary (and grouped by groupvar / nth_infection_cat)
        nest() %>% 
        mutate(agecat = "all",
               period      = "all")
    ) %>% 
    
    full_join(
      data %>% 
        # group by agecat AND period (and optionally groupvar / nth_infection_cat)
        group_by(agecat, period, {{ group_var }}) %>% 
        nest() %>% 
        mutate(variable = "age and period")
    ) %>% 
    
    mutate(
      variable = case_when(agecat == "all" & period == "all" ~ "summary",
                           agecat != "all" & period == "all" ~ "age",
                           agecat == "all" & period != "all" ~ "period",
                           agecat != "all" & period != "all" ~ "age and period"),
      variable = factor(variable, levels = c("summary","age","period", "age and period"))
    ) %>% 
    
    mutate(agecat = factor(agecat, levels = c("all", "18+y", "10-17y", "0-9y")),
           period = factor(period, levels = c("all", "original", "pre-omicron variant", "omicron"))) %>% 
    
    arrange(variable, agecat, period) %>% 
    relocate(variable, agecat, period)
  
}


# get numbers----
get_ns     <- function(nested_data){
  nested_data %>% 
    mutate(
      py = map_dbl(data, ~ nrow(.x) / 52 %>% round_half_up(digits = 1)),
      'person-years' = sprintf("%.1f", py),
      n  = map_dbl(data, ~ nrow(distinct(.x, codigo))),
      
      n.any     = map_dbl(data, ~ nrow(.x %>% filter(pcr==1))),
      n.sympt   = map_dbl(data, ~ nrow(.x %>% filter(pcr_sympt==1))),
      n.mod_sev = map_dbl(data, ~ nrow(.x %>% filter(pcr_mod_sev==1))),
      
      inc.any     = n.any / py * 100,
      inc.sympt   = n.sympt / py * 100,
      inc.mod_sev = n.mod_sev / py * 100,
      
      across(contains("inc."), ~ sprintf("%.1f", round_half_up(.x, digits = 1)))
    ) %>% 
    select(-data) 
  
}
get_ns_sev <- function(nested_data){
  nested_data %>% 
    mutate(
      n  = map_dbl(data, ~ nrow(.x)),
      
      
      # n's
      subclin   = map_dbl(data, ~ nrow(.x %>% filter(subclinical==1))),
      mild      = map_dbl(data, ~ nrow(.x %>% filter(severity_clinic_imp=="mild"))),
      mod_sev   = map_dbl(data, ~ nrow(.x %>% filter(severe_moderate==1))),
      sev       = map_dbl(data, ~ nrow(.x %>% filter(severe==1))),
      
      
      
      across(c(subclin, mild, mod_sev, sev), ~ 
               paste0(.x, " (", 
                      #round_half_up(.x / n * 100, digits = 1),  
                      sprintf("%.1f", round_half_up(
                        .x / n * 100, 1
                      )),
                      "%)"))
      
    ) %>% 
    select(-data)
}


# add_bold_title_pw <- function(plot){
#   list(
#     plot_annotation(title = "All ages"),
#     theme()
#   )
# }



# plot cumulative infections & vaccinations----
plot_cum_inf <- function(data){
  data %>% 
    ggplot() +
    geom_bar(aes(x    = week,
                 fill = inf_til_now_cat)) + 
    labs(fill = "prior infections",
         title = "cumulative infections (serology & RT-PCR)") +
    modify_date_plots() +
    scale_fill_manual(values = colors_num, drop = F) +
    add_shading_blood("red") 
}
plot_cum_vax <- function(data){
  data %>% 
    ggplot() +
    geom_bar(aes(x    = week,
                 fill = vax_til_now_num_cat)) +
    labs(fill = "prior vaccinations",
         title = "cumulative vaccinations") +
    modify_date_plots() +
    scale_fill_manual(values = colors_num) +
    add_shading_blood("red") 
}


# facets----
# protection:
facet_agecat_week        <- function(plot){
  list(
    facet_grid(rows = vars(agecat_week), scales = "free_y")
  )
}
facet_agecat_week_period <- function(plot){
  list(
    facet_grid(rows = vars(agecat_week),
               cols = vars(period), scales = "free_y")
  )
}

# severity:
facet_agecat        <- function(plot){
  list(
    facet_grid(rows = vars(agecat))
  )
}
facet_agecat_period <- function(plot){
  list(
    facet_grid(rows = vars(agecat),
               cols = vars(period))
  )
}



# plot incidence----
plot_inc     <- function(data, fill){
  data %>% 
    mutate(
      l = (100 / py) * (n - (1.96*sqrt(n))), 
      u = (100 / py) * (n + (1.96*sqrt(n))),
      
      
      # agecat_week = fct_recode(agecat_week,  "all ages" = "all"),
      # period      = fct_recode(period,    "all periods" = "all")
    ) %>%
    
    ggplot(aes(x = outcome,
               y = inc,
               fill = {{ fill }},
               
               ymin = l,
               ymax = u)) +
    geom_col(position = position_dodge2(0.9),
             color = "black",
             linewidth = 0.1) +
    geom_errorbar(position = position_dodge(0.9),
                  width = 0.2) +
    scale_fill_manual(values = colors_num) +
    theme_hm_lines(basesize) +
    labs(y     = "incidence / 100 py",
         x     = "disease outcome",
         title = "incidence of RT-PCR") +
    
    facet_grid(rows = vars(agecat_week)
               #cols = vars(period)
    )
}
plot_inc_inf <- function(data){
  data %>% 
    plot_inc(inf_til_now_cat) +
    labs(fill = "prior \ninfections") 
}



# protection models----
mod_prot_inf        <- function(data, outcome){
  outcome <- data[[outcome]]
  
  glm(outcome ~ inf_til_now_cat,
      family = poisson,
      data = data)
} # crude model
# models adjusted for vaccination:
mod_prot_inf_vax    <- function(data, outcome){
  outcome <- data[[outcome]]
  
  glm(outcome ~ inf_til_now_cat + vax_status,
      family = poisson,
      data = data)
}
mod_prot_inf_yn_vax <- function(data, outcome){
  outcome <- data[[outcome]]
  
  glm(outcome ~ inf_til_now_yn + vax_status,
      family = poisson,
      data = data)
}
mod_prot_period_li  <- function(data, outcome){
  outcome <- data[[outcome]]
  
  glm(outcome ~ period_li + vax_status,
      family = poisson,
      data = data)
}
mod_prot_time_li    <- function(data, outcome){
  outcome <- data[[outcome]]
  
  glm(outcome ~ days_since_last_inf_cat + vax_status,
      family = poisson,
      data = data)
}



tidy_irr <- function(mod){# tidy & add CI's,
  tidy(mod) %>% 
  #tidy(mod, conf.int = T, exponentiate = T) %>% # 2024-11-06: added 
    # rename(
    #   irr = estimate,
    #   l = conf.low,
    #   u = conf.high
    # ) 
    # 
    mutate(l   = exp(estimate - 1.96*std.error),
            u   = exp(estimate + 1.96*std.error),
            
            irr = exp(estimate),

           across(where(is.numeric), \(x) round_half_up(x, digits = 3)),

           #across(c(l, u), \(x) if_else(abs(x) > 500, NA_real_, x)),
           across(c(irr, l, u, p.value),
                  ~ if_else(l==0 | u==Inf,
                            NA_real_,.x)
                  ),

           CI  = paste0("(",
                        sprintf("%.2f", round_half_up(l, 2)),
                        ", ",
                        sprintf("%.2f", round_half_up(u, 2)),
                        ")"),
           irr_CI = paste0(sprintf("%.2f", round_half_up(irr, 2)),
                           " ", CI),

           across(c(irr_CI, CI), ~ if_else(is.na(irr),
                                       NA_character_,
                                       .x)),
           
           
           prot   = (1-irr)*100,
           prot_l = (1-u)*100,
           prot_u = (1-l)*100,
           
           across(c(prot, prot_l, prot_u), ~
                  sprintf("%.1f", round_half_up(.x, 1))
                  ),
           
           protection = paste0(prot," (", prot_l, "-", prot_u,")")
           
    ) %>%
    filter(term != "(Intercept)") %>%
    select(term, irr, l, u, p.value, irr_CI, protection) %>% 
    mutate(aic = mod$aic)
  
}

add_model_n <- function(nested_data, outcome, model){
  outcome <- paste(outcome)
  
  nested_data %>% 
    mutate(
      model = map2(data, outcome, {{ model }}),
      model_output = map(model, tidy_irr),
      outcome = {{ outcome }}
    ) %>% 
    select(-c(model)) 
}

add_pcr_sympt_modsev <- function(nested_data, model){
  
  nested_data             %>% add_model_n("pcr",         {{ model }}) %>% 
    full_join(nested_data %>% add_model_n("pcr_sympt",   {{ model }})) %>%
    full_join(nested_data %>% add_model_n("pcr_mod_sev", {{ model }})) 
  
}

unnest_mod <- function(model_results){
  model_results %>% 
    unnest(model_output) %>% 
    relocate(variable, agecat_week, period,
             outcome, term) %>% 
    select(-data) %>% 
    mutate(
      exposure = case_when(
        str_detect(term, "inf_til_now")             ~ paste0(str_sub(term, 16,-1), " inf"),
        str_detect(term, "vax_til_now")             ~ paste0(str_sub(term, 16,-1), " vax"),
        str_detect(term, "period_li")               ~ paste0(str_sub(term, 10,-1)),
        str_detect(term, "days_since_last_inf_cat") ~ paste0(str_sub(term, 24,-1)),
        TRUE ~ term
      ),
      period = factor(period, levels = c("all", "original", "pre-omicron variant", "omicron")),
      
      agecat_week = factor(agecat_week, levels = c("all", "18+y", "10-17y", "0-9y")),
      
      outcome = factor(outcome, # format for plots
                       levels = c("pcr","pcr_sympt","pcr_mod_sev", "pcr_sev"),
                       labels = c("pcr","sympt","mod/sev","sev"))
      
    ) %>% 
    arrange(variable, agecat_week, period) 
}

# just keep main term model results
unnest_mod_no_vax <- function(model_results){
  model_results %>% 
    unnest_mod() %>% 
    #filter(str_detect(term, "inf")) 
    filter(! str_detect(term, "vax")) 
}



# plot protection----
plot_prot <- function(data){
  data %>% 
    ggplot(aes(x = outcome,
               y = (1 - irr)*100,
               ymin = (1 - l)*100,
               ymax = (1 - u)*100,
               fill  = exposure,
               shape = exposure)) +
    
    geom_rect(fill = "#fde0dd", xmin = 0.5, xmax = 1.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    geom_rect(fill = "#ae017e", xmin = 1.5, xmax = 2.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    geom_rect(fill = "#49006a", xmin = 2.5, xmax = 3.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    
    geom_pointrange(position = position_dodge(width = 0.5), 
                    size = 1
    ) +
    coord_cartesian(ylim = c(-5, 100)) +
    
    scale_shape_manual(# pick shapes that can have a fill
      values = c( 
        21, # circle
        24, # triangle
        22, # square 
        23  # diamond
      ) 
    ) +
   
    labs(title = "protection from RT-PCR", 
         x     = "disease outcome",
         y     = "% protection"
    ) +
    theme_hm_lines(basesize)
  
}

# * modify protection plot function based on exposure variables----
plot_prot_inf     <- function(data){
  data %>% 
    plot_prot() +
    scale_fill_manual(values = colors_num[2:5]) + # infection #
    
    labs(
      fill    = "prior \ninfections\n(vs uninfected)",
      shape   = "prior \ninfections\n(vs uninfected)",
    )
}    # sets colors for prior inf
plot_prot_inf_adj <- function(data){
  data %>% 
    plot_prot_inf() +
    labs(caption = "adjusted for vaccination")
}# adds caption "adjusted for vax"
plot_prot_per_pi  <- function(data){
  data %>% 
    plot_prot() +
    
    scale_fill_manual(
      values = c("#fdc086",
                 "#386cb0",
                 "#7fc97f",
                 "#beaed4")
    ) +
    
  
    labs(
         fill    = "period of last infection\n(vs uninfected)",
         shape   = "period of last infection\n(vs uninfected)",
         caption = "adjusted for vaccination"
    ) +
    
    facet_grid(rows = vars(agecat_week),
               cols = vars(period)) 
  
} # defualt ggplot colors & nests by age and period
plot_prot_time_pi <- function(data){
  data %>% 
    plot_prot() +
    
    scale_fill_manual(
      values = c("#fed98e",
                 "#fe9929",
                 "#d95f0e",
                 "#993404")
    ) +
    
    labs(
      fill    = "time since last infection\n(vs uninfected)",
      shape   = "time since last infection\n(vs uninfected)",
      caption = "adjusted for vaccination"
    ) +
    
    facet_grid(rows = vars(agecat_week),
               cols = vars(period)) 
  
} # defualt ggplot colors & nests by age and period


# make tables----
widen_model  <- function(irrs_long){
  irrs_long %>% 
    mutate(inf_til_now_cat = str_sub(term, 16, -1)) %>% 
    arrange(variable, agecat_week, period) %>% 
    select(variable, agecat_week, period,
           outcome, inf_til_now_cat, irr_CI, protection) %>% 
    
    pivot_wider(id_cols     = c(variable, agecat_week, period,
                                inf_til_now_cat), 
                names_from   = outcome, 
                #names_prefix = "irr_CI_",
                values_from = c(irr_CI, protection)) %>% 
    select(variable:inf_til_now_cat, 
           irr_CI_pcr,
           irr_CI_sympt,
           "irr_CI_mod/sev",
           
           protection_pcr,
           protection_sympt,
           "protection_mod/sev") 
}
make_table_n <- function(nested_data, irrs_long){
  
  tab_ns   <- get_ns( {{ nested_data }} ) 
  tab_irrs <- widen_model( {{ irrs_long }} )
  
  tab <- tab_ns %>% 
    left_join(tab_irrs) %>% 
    select(-c(py)) %>% 
    mutate(period = factor(period, levels = c("all", "original", "pre-omicron variant","omicron"))) %>% 
    arrange(variable, agecat_week, period, inf_til_now_cat) %>% 
    
    mutate(across(contains("irr_"), ~ case_when(is.na(.x) & inf_til_now_cat==0 ~ "ref", 
                                               is.na(.x) ~ "not calculated",
                                               TRUE ~ .x))) 
  
  tab %>% 
    rename(
      "age"                    = agecat_week,
      "# prior infections"     = inf_til_now_cat,
      
      "RT-PCR \n(all)"         = n.any,
      "RT-PCR \n(symptomatic)" = n.sympt,
      "RT-PCR \n(mod/sev)"     = n.mod_sev,

      "Incidence \n(all)"         = inc.any,
      "Incidence \n(symptomatic)" = inc.sympt,
      "Incidence \n(mod/sev)"     = inc.mod_sev,

      "IRR \n(all)"         = irr_CI_pcr,
      "IRR \n(symptomatic)" = irr_CI_sympt,
      "IRR \n(mod/sev)"     = "irr_CI_mod/sev",
      
    ) %>% 
    ungroup()
  
}





# severity analyses----
plot_rtpcr_inf <- function(data, death_ht){
  data %>% 
    ggplot() +
    add_shading_blood("grey20") +
    geom_bar(aes(x    = inf_week,
                 fill = severity_clinic_imp),
             color = "black",
             linewidth = 0.1
    ) +                                    # bar plot
    
    geom_point(data = data %>% filter(death==1), # add deaths
               aes(y = {{ death_ht }},
                   x = inf_week),  
               color = "#49006a",
               shape = 8, size = 2, show.legend = F) +
    
    facet_wrap(~ label_nth_inf_cat, ncol = 1) +     # facet by infection #
    
    theme_hm(basesize) +
    theme(legend.position = "bottom") +
    scale_fill_manual(values = colors_severity) +
    scale_x_date(labels = scales::label_date_short(),
                 date_minor_breaks = "1 months",
                 date_breaks = "3 months",
                 limits = c(ymd("2020-01-01"), max(dat_sev_ms$inf_week) + 30)
    ) +
    labs(x = "week", 
         y = "weekly counts", 
         fill = "severity",
         title = "infection timing") 
}


label_positions_sev <- tribble( # these are the positions to put totals text on plot
  ~label_position, ~text_pos, ~ylimit,
  "in", 104, 105,
  "out", 113, 100
)

# * count and plot severity percents & n's----
sev_percents <- function(data, inout, x_var, facet_var_1, facet_var_2){
  
  pos <- label_positions_sev %>% filter(label_position== {{ inout }} )
  
  counts_sev <- data %>%
    group_by( {{ x_var }}, {{ facet_var_1 }}, {{ facet_var_2 }} ) %>%
    count(severity_clinic_imp) %>%
    mutate(percent = n / sum(n)*100)  # percents of severity by groups
  
  counts <- data %>%
    group_by({{ facet_var_1 }}, {{ facet_var_2 }}) %>%
    count( {{ x_var }} ) %>%
    mutate(percent = n / sum(n)*100) # n's (totals) by group
  
  
  
  plot <- ggplot(counts_sev, aes(x    = {{ x_var }},
                                 y    = percent,
                                 fill = severity_clinic_imp)) +
    geom_col() +
    facet_grid(rows = vars( {{ facet_var_1 }} ),
               cols = vars( {{ facet_var_2 }} )
               ) +
    
    geom_text(data = counts,
              aes(x = {{ x_var }},
                  y = pos$text_pos, #104,#113,
                  label = paste0("n=", n),
                  fill = NULL),
              size = 4
    ) +
    coord_cartesian(
      ylim = c(0, pos$ylimit),
      clip = 'off'
    ) +
    
    theme_hm(basesize) +
    theme(
      plot.margin = unit(c(3,1,1,1), "lines"), # this extends the top margin
      plot.title = element_text(vjust = 5)     # raise title
    ) + 
    scale_fill_manual(values = colors_severity) +
    labs(#x     = "infection #",
      fill  = "severity",
      title = "severity")
  
  list(
    plot, counts_sev, counts
  )
  
  
} # enter inout input as "in" or "out", facet vars are optional

# modify function for exposure: # prior infection
sev_prior_inf <- function(data, inout, facet_var_1, facet_var_2){
  output <- sev_percents(data, inout,  nth_infection_cat,  {{ facet_var_1 }}, {{ facet_var_2 }})
  
  output[[1]] <- output[[1]] + 
    labs(x = "infection #")
  
  output
} # enter inout input as "in" or "out", facet vars are optional



# * severity models----
mod_sev_inf         <- function(data, severity_level){
  severity_level <- data[[severity_level]]
  
  glm(
    severity_level ~ nth_infection_cat,
    family = poisson,
    data = {{ data }}
  ) 
}
mod_sev_inf_vax     <- function(data, severity_level){
  severity_level <- data[[severity_level]]
  
  glm(
    severity_level ~ nth_infection_cat + vax_status,
    family = poisson,
    data = {{ data }}
  ) 
}
mod_sev_per_li      <- function(data, severity_level){
  severity_level <- data[[severity_level]]
  
  glm(
    severity_level ~ period_li + vax_status,
    family = poisson,
    data = {{ data }}
  ) 
}
mod_sev_time_li      <- function(data, severity_level){
  severity_level <- data[[severity_level]]
  
  glm(
    severity_level ~ days_since_last_inf_cat + vax_status,
    family = poisson,
    data = {{ data }}
  ) 
}


tidy_rr           <- function(mod){# tidy & add CI's
  tidy(mod) %>% 
    mutate(l   = exp(estimate - 1.96*std.error),
           u   = exp(estimate + 1.96*std.error),
           
           
           across(where(is.numeric), ~ round_half_up(.x, digits = 3)),
           ratio = ( exp(estimate)),
           
           across(c(ratio, l, u, p.value), ~ if_else(l==0 | u==Inf, 
                                            NA_real_,
                                            .x)),
           
           CI  = paste0("(", 
                        sprintf("%.2f", round_half_up(l,2)), 
                        ", ", 
                        sprintf("%.2f", round_half_up(u,2)), 
                        ")"),
           rr = paste0(
             sprintf("%.2f", round_half_up(ratio, 2)),
             " ", CI
           ),
           
           across(c(rr, CI), ~ if_else(is.na(ratio), 
                                            NA_character_,
                                            .x)),
    ) %>% 
    filter(term != "(Intercept)") %>% 
    select(term, rr, p.value, ratio, CI, l, u)
  
}


add_severity_mod <- function(data, model){
  data <- data
  
  add_mod <- function(data, model, outcome){
    data %>% 
      mutate(model = map(data, {{ model }}, {{ outcome }}),
             outcome = {{ outcome }},
             model_results = map(model, tidy_rr)
      ) 
  }
  
  data %>%             
    add_mod({{ model }},                    "severe") %>% 
    full_join(data %>% add_mod({{ model }}, "severe_moderate")) %>% 
    full_join(data %>% add_mod({{ model }}, "subclinical"))
  
}


format_model_results <- function(data){
  data %>% 
    unnest(model_results) %>% 
    mutate(
      exposure = case_when(
        str_detect(term, "nth_infection_cat")       ~ "inf_n",
        str_detect(term, "vax_status")              ~ "prior_vax",
        str_detect(term, "period_li")               ~ "period_li",
        str_detect(term, "days_since_last_inf_cat") ~ "days_since_last_inf_cat"
      ),
      
      outcome = factor(outcome, 
                       levels = c("subclinical","severe_moderate","severe"),
                       labels = c("subclin",
                                  "mod/sev",
                                  "severe")),
      
      term = case_when(
        str_detect(term, "nth_infection_cat")       ~ str_sub(term, 18, -1),
        str_detect(term, "vax_status")              ~ str_sub(term, 11, -1),
        str_detect(term, "period_li")               ~ str_sub(term, 10, -1),
        str_detect(term, "days_since_last_inf_cat") ~ str_sub(term, 24, -1),
        TRUE ~ term
      ),
    
    label_level = case_when(
      exposure=="inf_n"                   ~ paste0(term, " vs 1st"),
      exposure=="prior_vax"               ~ paste0(term, " vs 0"),
      exposure=="period_li"               ~ paste0(term),
      exposure=="days_since_last_inf_cat" ~ paste0(term)
      )
    ) %>% 
    filter(exposure != "prior_vax") %>% 
    rename(level = term)  %>%
    select(-c(data, model)) %>% 
    relocate(exposure, .after = outcome) %>% 
    arrange(variable, agecat, outcome) %>% 
    ungroup()
}


# plot severity----
plot_rel_sev <- function(data){
  
  data %>% 
    filter(outcome != "severe") %>% 
    
    ggplot(aes(
      x = outcome,
      y = ratio,
      ymin = l,
      ymax = u,
      fill  = label_level, 
      shape = label_level
    )) +
    
    geom_rect(fill = "#fde0dd", xmin = 0.5, xmax = 1.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    #geom_rect(fill = "#ae017e", xmin = 1.5, xmax = 2.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    #geom_rect(fill = "#49006a", xmin = 2.5, xmax = 3.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    geom_rect(fill = "#49006a", xmin = 1.5, xmax = 2.5, alpha = 0.05, ymin = -Inf, ymax = Inf) +
    
    geom_hline(yintercept = 1) +
    
    geom_pointrange(size = 0.9, position = position_dodge2(width = 0.4)) +
    
    scale_y_log10() +
    
    labs(y     = "relative severity",
         x     = "severity",
         title = "relative severity",
         shape = NULL,
         fill  = NULL) +
    theme_hm_lines(basesize) +
    theme(plot.title = element_text(vjust = 5))  # raise title, to match fig 2b
}

#  modify protection plot function based on exposure variables
plot_rel_sev_nth     <- function(data){
  data %>% 
    plot_rel_sev() +
    scale_fill_manual(values = colors_num[3:5]) +
    scale_shape_manual(
      values = c( 
        24, # triangle
        22  # square 
      ) 
    ) 
    
}
plot_rel_sev_nth_adj <- function(data){
  data %>% 
    plot_rel_sev_nth() +
    labs(caption = "adjusted for vaccination")
}
plot_rel_sev_per_pi  <- function(data){
  data %>% 
    plot_rel_sev() +
    scale_shape_manual(# pick shapes that can have a fill
      values = c( 
        21, # circle
        24, # triangle
        22, # square 
        23  # diamond
      ) 
    ) +
    
    scale_fill_manual(
      values = c("#fdc086",
                 "#386cb0",
                 "#7fc97f",
                 "#beaed4")
    ) +
    
    labs(caption = "adjusted for vaccination",
         fill    = "period of last infection\n(vs uninfected)",
         shape   = "period of last infection\n(vs uninfected)") +
    
    facet_grid(rows = vars(agecat),
               cols = vars(period)) 
}

plot_rel_sev_time_pi  <- function(data){
  data %>% 
    plot_rel_sev() +
    scale_shape_manual(# pick shapes that can have a fill
      values = c( 
        21, # circle
        24, # triangle
        22, # square 
        23  # diamond
      ) 
    ) +
    
    scale_fill_manual(
      values = c("#fed98e",
                 "#fe9929",
                 "#d95f0e",
                 "#993404")
    ) +
    
    labs(caption = "adjusted for vaccination",
         fill    = "time since last infection\n(vs uninfected)",
         shape   = "time since last infection\n(vs uninfected)") +
    
    facet_grid(rows = vars(agecat),
               cols = vars(period)) 
}



combine_sev_plots <- function(p1, p2_object, p3){
  ( {{ p1 }} + theme(legend.position = "none")) / 
    ( {{ p2_object }} | {{ p3 }} ) + 
    plot_layout(heights = c(1, 0.4)) & 
    plot_annotation(tag_levels = "A")
}


# severity tables----
widen_model_sev  <- function(rs_long){
  rs_long %>% 
    rename(infection_n = level) %>% 
    pivot_wider(id_cols = c(variable, agecat, period, infection_n),
                names_from = outcome,
                names_prefix = "rr_CI_",
                values_from = c(rr)) 
}



