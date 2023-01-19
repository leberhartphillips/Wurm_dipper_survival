#### European Dipper survival modelling
## January 18, 2022
## Data: Ulrich Knief
## Script: Luke Eberhart-Phillips

#### Libraries
if (!'RMark' %in% installed.packages()) install.packages('RMark') ; require(RMark)
if (!'tidyverse' %in% installed.packages()) install.packages('tidyverse') ; require(tidyverse)
if (!'dplyr' %in% installed.packages()) install.packages('dplyr') ; require(dplyr)
if (!'BaSTA' %in% installed.packages()) install.packages('BaSTA') ; require(BaSTA)
if (!'ggpubr' %in% installed.packages()) install.packages('ggpubr') ; require(ggpubr)
if (!'lubridate' %in% installed.packages()) install.packages('lubridate') ; require(lubridate)
if (!'viridis' %in% installed.packages()) install.packages('viridis') ; require(viridis)
if (!'patchwork' %in% installed.packages()) install.packages('patchwork') ; require(patchwork)
if (!'gt' %in% installed.packages()) install.packages('gt') ; require(gt)
if (!'RColorBrewer' %in% installed.packages()) install.packages('RColorBrewer') ; require(RColorBrewer)

# path <- "C:\\Users\\Knief\\Dropbox\\DBDipper\\SurvivalMARK\\"
path <- "data/"

#### Prepare data ---
dipper_df <- read.table(paste0(path,"data_MARK_dipper_2022.txt"), header = TRUE)
dipper_df$date <- as.Date(dipper_df$date, format = "%Y-%m-%d")

dipper_df_join <- 
  dipper_df %>% 
  dplyr::group_by(ring) %>% 
  dplyr::summarise(first_age = age[which.min(date)],
            first_site = site[which.min(date)]) %>% 
  dplyr::left_join(dipper_df, ., by = "ring")

# check if any individuals have more than one sex
dipper_df_join %>% 
  group_by(ring) %>% 
  summarise(n_sexes = n_distinct(sex)) %>% 
  filter(n_sexes > 1)

# check sampling distribution across sites
dipper_df_join %>% 
  group_by(age, first_site) %>% 
  summarise(n = n())

# extract year from date string
dipper_df_join$year <- format(dipper_df_join$date,"%Y")

# wrangle capture history for MARK
dipper_df_ch <- 
  dipper_df %>% 
  CensusToCaptHist(ID = dipper_df_join$ring, 
                   d = as.numeric(dipper_df_join$year),
                   timeInt = "Y") %>% 
  mutate(ring = row.names(.)) %>% 
  dplyr::left_join(., dplyr::select(dipper_df_join, ring, sex, first_age, first_site), by = "ring") %>%
  distinct() %>% 
  dplyr::rename(site = first_site,
         age = first_age)

# Check sampling distribution 
sum_dat <- 
  dipper_df %>% 
  dplyr::mutate(year = year(date)) %>% 
  dplyr::select(-date) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(site, ring, age, sex) %>%
  dplyr::summarise(n_encounters = n()) %>% 
  dplyr::group_by(site, age, sex) %>% 
  dplyr::summarise(total_encounters = sum(n_encounters),
            total_ind = n_distinct(ring))
sum_dat %>% print(n = Inf)

# make capture history for MARK
# 1. dataset to assess site and age-specific survival and encounter rates. Site
# in this case is "Wurm" vs. "not Wurm" due to the limited sampling at non-Wurm
# sites.
dipper_df_ch_age_site <- 
  data.frame(ch = apply(dipper_df_ch[, 2:7], 1, paste, collapse = "")) %>%
  bind_cols(., dplyr::select(dipper_df_ch, sex, site, age)) %>%
  mutate(across(everything(), ~str_trim(.x))) %>%
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>%
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  mutate(age = as.factor(age),
         site = as.factor(ifelse(site == "Wuerm", "Wuerm", "Not_Wuerm")),
         sex = as.factor(sex)) %>% 
  dplyr::select(-sex)

# 2. dataset to assess sex-specific survival and encounter rates. In this case,
# we have subset the data to individuals of known sex
dipper_df_ch_sex <- 
  data.frame(ch = apply(dipper_df_ch[, 2:7], 1, paste, collapse = "")) %>%
  bind_cols(., dplyr::select(dipper_df_ch, sex, site, age)) %>%
  mutate(across(everything(), ~str_trim(.x))) %>%
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>%
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  mutate(age = as.factor(age),
         sex = as.factor(sex)) %>% 
  dplyr::filter(sex != "U") %>% 
  dplyr::select(-site, -age)

# Create processed RMark data formatted as Cormack-Jolly_Seber with groups
# "Age" (age initally ringed), "sex", and "site", starting at year 2017, two age groups
# (first-years and adults) in which the first-year stage only lasts for 
# one year.
dipper_age_site.proc <- RMark::process.data(dipper_df_ch_age_site, model = "CJS",
                                        group = c("age", "site"),
                                        begin.time = 2017,
                                        initial.age = c(1, 0))

dipper_sex.proc <- RMark::process.data(dipper_df_ch_sex, model = "CJS",
                                           group = c("sex"),
                                           begin.time = 2017)

# Create the design matrix from the processed mark-recapture datasets
dipper_age_site.ddl <- RMark::make.design.data(dipper_age_site.proc)
dipper_sex.ddl <- RMark::make.design.data(dipper_sex.proc)

# adds juvenile / adult age field to design data in column "age"
dipper_age_site.ddl <- RMark::add.design.data(data = dipper_age_site.proc,
                                              ddl = dipper_age_site.ddl,
                                              parameter = "Phi",
                                              type = "age",
                                              bins = c(0, 1, dipper_age_site.proc$nocc),
                                              right = FALSE,
                                              name = "age", replace = TRUE)

dipper_age_site.ddl <- RMark::add.design.data(data = dipper_age_site.proc,
                                              ddl = dipper_age_site.ddl,
                                              parameter = "p",
                                              type = "age",
                                              bins = c(0, 1, dipper_age_site.proc$nocc),
                                              right = TRUE,
                                              name = "age", replace = TRUE)

# check the parameter index matricies to see if the dummy variable is
# is structured correctly...looks good
PIMS(mark(dipper_age_site.proc,
          dipper_age_site.ddl,
          model.parameters = list(Phi = list(formula = ~ age + time + site)),
          output = F, delete = TRUE, brief = TRUE,
          silent = TRUE, wrap = FALSE),
     "Phi")

PIMS(mark(dipper_age_site.proc,
          dipper_age_site.ddl,
          model.parameters = list(p = list(formula = ~ age + time + site)),
          output = F, delete = TRUE, brief = TRUE,
          silent = TRUE, wrap = FALSE),
     "p")

#### Model selection ####
# first assess variation in encounter probability while keeping 
# apparent survival constant
age_site_dipper_encounter_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula =  ~ 1)
    
    # age-dependent:
    p.age = list(formula =  ~ age)
    
    # site-dependent:
    p.site = list(formula =  ~ site)
    
    # age- and site-dependent:
    p.age_site = list(formula =  ~ age + site)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

sex_dipper_encounter_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula =  ~ 1)
    
    # age-dependent:
    p.sex = list(formula =  ~ sex)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

# run encounter models
age_site_dipper_encounter_analysis_out <-
  age_site_dipper_encounter_analysis_run(proc_data = dipper_age_site.proc, 
                                         design_data = dipper_age_site.ddl)

sex_dipper_encounter_analysis_out <-
  sex_dipper_encounter_analysis_run(proc_data = dipper_sex.proc, 
                                        design_data = dipper_sex.ddl)

# then assess variation in apparent survival with age-dependent encounter rate
age_site_dipper_survival_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # age-dependent:
    Phi.age = list(formula =  ~ age)
    
    # site-dependent:
    Phi.site = list(formula =  ~ site)
    
    # age- and site-dependent:
    Phi.age_site = list(formula =  ~ age + site)
    
    # Models exploring variation in encounter probability
    # age-dependent (based on encounter analysis, above, this is the best model):
    p.age = list(formula =  ~ age)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

sex_dipper_survival_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # sex-dependent:
    Phi.sex = list(formula =  ~ sex)
    
    # Models exploring variation in encounter probability
    # constant (based on encounter analysis, above):
    p.dot = list(formula =  ~ 1)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
                                       ddl = design_data, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

# run apparent survival models
age_site_dipper_survival_analysis_out <-
  age_site_dipper_survival_analysis_run(proc_data = dipper_age_site.proc, 
                                         design_data = dipper_age_site.ddl)

sex_dipper_survival_analysis_out <-
  sex_dipper_survival_analysis_run(proc_data = dipper_sex.proc, 
                                       design_data = dipper_sex.ddl)


#### Dredge: explore all model combinations ----
# (caution this is likely to produce false positives due to sample size issues)
# dredge_dipper_survival_analysis_run = 
#   function(proc_data, design_data){
#     # Models exploring variation in apparent survival
#     # constant:
#     Phi.dot = list(formula =  ~ 1)
#     
#     # age-dependent:
#     Phi.age = list(formula =  ~ age)
#     
#     # factorial variation across year:
#     Phi.year = list(formula =  ~ time)
#     
#     # additive effects of age and factorial year:
#     Phi.age_year = list(formula =  ~ age + time)
#     
#     # additive effects of age and linear year:
#     Phi.age_Time = list(formula =  ~ age + Time)
#     
#     # interactive effects of age and year:
#     Phi.age_x_Time = list(formula =  ~ age * Time)
#     
#     # Models exploring variation in encounter probability
#     # constant:
#     p.dot = list(formula =  ~ 1)
#     
#     # age-dependent:
#     p.age = list(formula =  ~ age)
#     
#     # factorial variation across year:
#     p.year = list(formula =  ~ time)
#     
#     # additive effects of sex and factorial year:
#     p.age_year = list(formula =  ~ age + time)
#     
#     # create a list of candidate models for all the a models above that begin with 
#     # either "Phi." or "p."
#     cml <-  RMark::create.model.list("CJS")
#     
#     # specify the data, design matrix, delete unneeded output files, and 
#     # run the models in Program MARK
#     model.list <-  RMark::mark.wrapper(cml, data = proc_data, 
#                                        ddl = design_data, delete = TRUE, 
#                                        wrap = FALSE, threads = 1, brief = TRUE,
#                                        silent = TRUE, output = FALSE)
#     
#     # output the model list and sotre the results
#     return(model.list)
#   }
# 
# dredge_dipper_survival_analysis_full_out <-
#   dredge_dipper_survival_analysis_run(proc_data = dipper_full.proc, 
#                                       design_data = dipper_full.ddl)
# dredge_dipper_survival_analysis_wurm_out <-
#   dredge_dipper_survival_analysis_run(proc_data = dipper_wurm.proc, 
#                                       design_data = dipper_wurm.ddl)
# dredge_dipper_survival_analysis_not_wurm_out <-
#   dredge_dipper_survival_analysis_run(proc_data = dipper_not_wurm.proc, 
#                                       design_data = dipper_not_wurm.ddl)

##### Results ####
# Extract the AIC model table from the model output
# not alot of support for a "significant-effect" of age, but it is slightly
# better than the null model and we have good a priori reasoning that there
# should be an effect
age_site_AIC_table_conservative <- 
  age_site_dipper_survival_analysis_out$model.table
sex_AIC_table_conservative <- 
  sex_dipper_survival_analysis_out$model.table

# Find the model number for the first ranked model of the AIC table
top_mod_num_sex <- 
  as.numeric(rownames(sex_dipper_survival_analysis_out$model.table[1,]))
top_mod_num_age_site <- 
  as.numeric(rownames(age_site_dipper_survival_analysis_out$model.table[1,]))

# extract and format the transformed parameter estimates
estimates_sex <- 
  sex_dipper_survival_analysis_out[[top_mod_num_sex]]$results$real %>% 
  bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                       n = 5)), .) %>% 
  dplyr::mutate(sex = as.factor(ifelse(unlist(str_extract_all(X2,"[FM]")) == "F", 
                                       "Female","Male")),
                parameter = as.factor(ifelse(unlist(str_extract_all(X1,"[pP]")) == "P", 
                                             "Phi","p")),
                analysis = "sex-specific") %>% 
  dplyr::select(analysis, parameter, sex, estimate, lcl, ucl) %>% 
  `rownames<-`( NULL )

estimates_age_site <- 
  age_site_dipper_survival_analysis_out[[top_mod_num_age_site]]$results$real %>% 
  bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                       n = 5)), .) %>% 
  dplyr::mutate(age = as.factor(ifelse(unlist(str_extract_all(X2,"[AJ]")) == "A", 
                                       "Adult","Juvenile")),
                site = as.factor(ifelse(str_detect(X2, "_"),
                                              "Other","Würm")),
                parameter = as.factor(ifelse(unlist(str_extract_all(X1,"[pP]")) == "P", 
                                             "Phi","p")),
                analysis = "age- and site-specific") %>% 
  dplyr::select(analysis, parameter, site, age, estimate, lcl, ucl) %>% 
  `rownames<-`( NULL )

estimates_combined <- 
  bind_rows(estimates_sex, estimates_age_site) %>% 
  mutate(coefString = ifelse(!is.na(lcl),
                             paste0("[", 
                                    round(lcl, 2), ", ", 
                                    round(ucl, 2), "]"),
                             NA))

betas_sex <- 
  sex_dipper_survival_analysis_out[[top_mod_num_sex]]$results$beta %>% 
  bind_cols(data.frame(X1 = rownames(.)), .) %>% 
  dplyr::mutate(component = as.factor(ifelse(str_detect(X1,"Intercept"), "Intercept", "Sex (Male)")), 
                parameter = as.factor(ifelse(str_detect(X1,"Phi"), "Phi", "p")),
                analysis = "sex-specific") %>% 
  dplyr::select(analysis, parameter, component, estimate, lcl, ucl) %>% 
  `rownames<-`( NULL )

betas_age_site <- 
  age_site_dipper_survival_analysis_out[[top_mod_num_age_site]]$results$beta %>% 
  bind_cols(data.frame(X1 = rownames(.)), .) %>% 
  dplyr::mutate(component = as.factor(ifelse(str_detect(X1,"Intercept"), "Intercept", 
                                             ifelse(str_detect(X1,"site"), "Site (Würm)","Age (Adult)"))), 
                parameter = as.factor(ifelse(str_detect(X1,"Phi"), "Phi", "p")),
                analysis = "age- and site-specific") %>% 
  dplyr::select(analysis, parameter, component, estimate, lcl, ucl) %>% 
  `rownames<-`( NULL )

betas_combined <- 
  bind_rows(betas_sex, betas_age_site) %>% 
  mutate(coefString = ifelse(!is.na(lcl),
                             paste0("[", 
                                    round(lcl, 2), ", ", 
                                    round(ucl, 2), "]"),
                             NA))

#### Uli plotting ####
# clr1 <- viridis(10)[6]
# clr2 <- viridis(11)[11]
# svg(filename=paste(path,"Survival1.svg",sep=""), height=90/25.4, width=70/25.4, family="Arial", pointsize = 10)
# par(mar=c(2.8, 2.8, 0.2, 0.2))
# par(mgp=c(1.2,0.2,0))
# plot(c(0,1),c(0,1), axes=FALSE, xlab="Altersklasse", ylab="Lokale Uberlebenswahrscheinlichkeit", tcl=-0.25, type="n")
# axis(2, at=seq(0,1,by=0.25), labels=c("0.00","0.25","0.50","0.75","1.00"), tcl=-0.25)
# axis(1, at=c(0.2,0.8), labels=c("Jungvogel","Altvogel"), tcl=-0.25)
# arrows(0.2,estimates$lcl[2],0.2,estimates$ucl[2],code=3,length=0.1,angle=90)
# arrows(0.8,estimates$lcl[1],0.8,estimates$ucl[1],code=3,length=0.1,angle=90)
# points(c(0.2), estimates$estimate[2], pch=23, bg=clr2, cex=sqrt(42)/2, type="b")
# points(c(0.8), estimates$estimate[1], pch=23, bg=clr1, cex=sqrt(33)/2, type="b")
# box()
# dev.off()

#  parameter      age  estimate       lcl       ucl
#1         p    Adult 0.4400000 0.2629157 0.6337978
#2         p Juvenile 0.2258065 0.1116412 0.4036688
#3       Phi    Adult 1.0000000 0.9999556 1.0000444
  
  
#  parameter      age  estimate       lcl       ucl
#1         p    Adult 0.4400000 0.2629157 0.6337978
#2         p Juvenile 0.2592593 0.1289173 0.4528708
#3       Phi    Adult 1.0000000 1.0000000 1.0000000

  
#### Plotting ####
# first setup the plotting theme
dipper_theme <- 
  theme(text = element_text(family = "Franklin Gothic Book"),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 10,margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.text.y  = element_text(size = 9),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.y = element_line(size = 0.2, colour = "grey40"),
        axis.ticks.length = unit(0.1, "cm"),
        axis.ticks.x = element_line(size = 0.2, colour = "grey40"),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.spacing = unit(0.3, "lines"),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        panel.border = element_rect(linetype = "solid", colour = "grey"),
        panel.background = element_rect(linetype = "solid", fill = "white"))


# get sample sizes of individuals used to estimate age-specific rates
dipper_df %>% 
  dplyr::select(ring, age) %>%
  distinct() %>% 
  group_by(age) %>% 
  summarise(N = n())

# age-specific apparent survival plot
dipper_df_ch_sex %>% 
  group_by(sex) %>% 
  summarise(N = n())

Sex_Phi_plot <- 
  ggplot(data = filter(estimates_sex, parameter == "Phi"),
         aes(x = sex, y = estimate, color = sex)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), 
               size = 0.3, linetype = "solid", width = 0.1) +
  geom_point(size = 3) +
  theme_bw() + 
  dipper_theme +
  ylab("Apparent survival (± 95% CI)") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Set1") +
  annotate(geom = "text", y = 0.6, x = 1, angle = 0,
           label = "N = 21", 
           color = "grey30", size = 2.5, fontface = 'italic') +
  annotate(geom = "text", y = 0.45, x = 2, angle = 0,
           label = "N = 25", 
           color = "grey30", size = 2.5, fontface = 'italic') 

dipper_df_ch_age_site %>% 
  mutate(recuit = ifelse(str_detect(ch, "11") & age == "J", "yes", "no")) %>% 
  group_by(site, age, recuit) %>% 
  summarise(N = n())

Age_Site_Phi_plot <- 
  ggplot(data = filter(estimates_age_site, parameter == "Phi"),
         aes(x = site, y = estimate, color = age)) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl), position = position_dodge(1), 
                size = 0.3, linetype = "solid", width = 0.1) +
  geom_point(size = 3, position = position_dodge(1)) +
  theme_bw() + 
  dipper_theme +
  theme(legend.position = "top",
        axis.title.y = element_blank(),
        axis.text.y  = element_blank()) +
  ylab("Apparent survival (± 95% CI)") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Dark2") +
  annotate(geom = "text", y = 0.73, x = 0.75, angle = 0,
           label = "N = 12", 
           color = "grey30", size = 2.5, fontface = 'italic') +
  annotate(geom = "text", y = 0.25, x = 1.25, angle = 0,
           label = "N = 66", 
           color = "grey30", size = 2.5, fontface = 'italic')  +
  annotate(geom = "text", y = 0.8, x = 1.75, angle = 0,
           label = "N = 33", 
           color = "grey30", size = 2.5, fontface = 'italic') +
  annotate(geom = "text", y = 0.33, x = 2.25, angle = 0,
           label = "N = 93", 
           color = "grey30", size = 2.5, fontface = 'italic')

# arrange multi-paneled plot
Wurm_plot <-
  Sex_Phi_plot + Age_Site_Phi_plot + 
  plot_annotation(
    title = 'CJS model of European Dippers (2017-2022)') + 
  plot_layout(widths = c(1, 1.5))

# export plot to working directory
ggsave(Wurm_plot,
       filename = "figs/age_site_sex_dipper_survival_2017-2022.jpeg",
       width = 5,
       height = 4, 
       units = "in")

#### Tabulation of results ####
dipper_CJS_effect_size_table <- 
  estimates_combined %>% 
  mutate(structure = ifelse(parameter != "Phi" & str_detect(analysis, "sex"), "-", 
                            ifelse(parameter == "p" & str_detect(analysis, "site"), as.character(age),
                                   ifelse(is.na(sex) & parameter == "Phi", 
                                          paste0(site, ", ", age), 
                                          as.character(sex))))) %>% 
  dplyr::select(parameter, structure, estimate, analysis, coefString) %>%
  gt(groupname_col = "analysis") %>% 
  cols_label(parameter = html("<i>Model</i>"),
             estimate = "Mean estimate",
             coefString = "95% confidence interval",
             structure = html("<i>structure</i>")) %>% 
  fmt_number(columns = vars(estimate),
             rows = 1:9,
             decimals = 2,
             use_seps = FALSE) %>% 
  cols_align(align = "right",
             columns = vars(parameter)) %>% 
  cols_align(align = "left",
             columns = vars(structure)) %>% 
  tab_options(row_group.font.weight = "bold",
              row_group.background.color = brewer.pal(9,"Greys")[3],
              table.font.size = 12,
              data_row.padding = 3,
              row_group.padding = 4,
              summary_row.padding = 2,
              column_labels.font.size = 14,
              row_group.font.size = 12,
              table.width = pct(50)) %>% 
  tab_header(
    title = "CJS model of European Dippers (2017 to 2022)",
    subtitle = "Effect sizes"
  )

dipper_CJS_effect_size_table %>%
  gtsave("dipper_CJS_effect_size_table.png", path = "tables/")

dipper_CJS_betas_table <- 
  betas_combined %>% 
  dplyr::select(parameter, component, estimate, analysis, coefString) %>%
  gt(groupname_col = "analysis") %>% 
  cols_label(parameter = html("<i>Model</i>"),
             estimate = "Mean estimate",
             coefString = "95% confidence interval",
             component = html("<i>structure</i>")) %>% 
  fmt_number(columns = vars(estimate),
             rows = 1:8,
             decimals = 2,
             use_seps = FALSE) %>% 
  cols_align(align = "right",
             columns = vars(parameter)) %>% 
  cols_align(align = "left",
             columns = vars(component)) %>% 
  tab_options(row_group.font.weight = "bold",
              row_group.background.color = brewer.pal(9,"Greys")[3],
              table.font.size = 12,
              data_row.padding = 3,
              row_group.padding = 4,
              summary_row.padding = 2,
              column_labels.font.size = 14,
              row_group.font.size = 12,
              table.width = pct(50)) %>% 
  tab_header(
    title = "CJS model of European Dippers (2017 to 2022)",
    subtitle = "Fixed effect coefficients"
  )

dipper_CJS_betas_table %>%
  gtsave("dipper_CJS_betas_table.png", path = "tables/")