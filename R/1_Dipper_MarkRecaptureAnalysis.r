#### European Dipper survival modelling
## September 28, 2020
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

path <- "C:\\Users\\Knief\\Dropbox\\DBDipper\\SurvivalMARK\\"

# (1) Prepare data ---
dipper_df <- read.table(paste0(path,"data_MARK_dipper_2022.txt"), header = TRUE)
dipper_df$date <- as.Date(dipper_df$date, format = "%Y-%m-%d")

dipper_df_join <- 
  dipper_df %>% 
  dplyr::group_by(ring) %>% 
  dplyr::summarise(first_age = age[which.min(date)],
            first_site = site[which.min(date)]) %>% 
  dplyr::left_join(dipper_df, ., by = "ring")

dipper_df_join$year <- format(dipper_df_join$date,"%Y")

dipper_df_ch <- 
  dipper_df %>% 
  CensusToCaptHist(ID = dipper_df_join$ring, 
                   d = as.numeric(dipper_df_join$year),
                   timeInt = "Y") %>% 
  dplyr::rename(ring = ID) %>%
  dplyr::left_join(., dplyr::select(dipper_df_join, ring, sex, first_age, first_site), by = "ring") %>%
  distinct() %>% 
  dplyr::rename(site = first_site,
         age = first_age)

# Check sampling distribution 
# conclusion: 
# best to model Wurm alone due to poor sampling at other sites
# drop sex due to limited high numbers of individuals of unknown sex
# focus on age (a priori most likely to have an effect)
sum_dat <- dipper_df %>% 
  dplyr::mutate(year = year(date)) %>% 
  dplyr::select(-date) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(site, ring, age, sex) %>%
  dplyr::summarise(n_encounters = n()) %>% 
  dplyr::group_by(site, age, sex) %>% 
  dplyr::summarise(total_counters = sum(n_encounters),
            total_ind = n_distinct(ring))
sum_dat %>% print(n = Inf)

# make capture history for MARK
dipper_df_ch_wurm <- 
  data.frame(ch = apply(dipper_df_ch[, 2:7], 1, paste, collapse = "")) %>%
  bind_cols(., dplyr::select(dipper_df_ch, sex, site, age)) %>%
  mutate(across(everything(), ~str_trim(.x))) %>%
  mutate(across(everything(), ~str_replace_all(.x, fixed(" "), ""))) %>%
  mutate(across(everything(), ~gsub("^$|^ $", NA, .x))) %>% 
  filter(site != "Wuerm") %>%
  mutate(age = as.factor(age)) %>% 
  dplyr::select(-sex, -site)

nrow(subset(dipper_df_ch_wurm, age=="J"))
nrow(subset(dipper_df_ch_wurm, age=="A"))

# Create processed RMark data formatted as Cormack-Jolly_Seber with 1 group
# (age initally ringed), starting at year 2017, two age groups
# (first-years and adults) in which the first-year stage only lasts for 
# one year.
dipper.proc <- RMark::process.data(dipper_df_ch_wurm, model = "CJS",
                                   group = "age",
                                   begin.time = 2017,
                                   initial.age = c(1, 0))

# Create the design matrix from the processed mark-recapture datasets
dipper.ddl <- RMark::make.design.data(dipper.proc)

# adds juvenile / adult age field to design data in column "age"
dipper.ddl <- RMark::add.design.data(data = dipper.proc,
                                     ddl = dipper.ddl,
                                     parameter = "Phi",
                                     type = "age",
                                     bins = c(0, 1, 4),
                                     right = FALSE,
                                     name = "age", replace = TRUE)

dipper.ddl$Phi$age[5] <- "[1,4]"

# create a dummy variable in the p design matrix called age2 which
# is "0" for all ages > 1 and "1" for all ages < 2 (i.e, the juvenile stage)
dipper.ddl$p$age2 = 0
dipper.ddl$p$age2[dipper.ddl$p$group == "J" & dipper.ddl$p$Age < 2] = 1

# check the parameter index matricies to see if the dummy variable is
# is structured correctly...looks good
PIMS(mark(dipper.proc,
          dipper.ddl,
          model.parameters = list(Phi = list(formula = ~ age + time)),
          output = F, delete = TRUE, brief = TRUE,
          silent = TRUE, wrap = FALSE),
     "Phi")

PIMS(mark(dipper.proc,
          dipper.ddl,
          model.parameters = list(p = list(formula = ~ age2 + time)),
          output = F, delete = TRUE, brief = TRUE,
          silent = TRUE, wrap = FALSE),
     "p")

#### Model selection ####
# first assess variation in encounter probability while keeping 
# apparent survival constant
dipper_encounter_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula =  ~ 1)
    
    # age-dependent:
    p.age2 = list(formula =  ~ age2)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = dipper.proc, 
                                       ddl = dipper.ddl, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

# run encounter models
# best model is the null model
dipper_encounter_analysis_out <-
  dipper_encounter_analysis_run()

# then assess variation in apparent survival with age-dependent encounter rate
dipper_survival_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # age-dependent:
    Phi.age = list(formula =  ~ age)
    
    # Models exploring variation in encounter probability
    # constant:
    p.age2 = list(formula =  ~ age2)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = dipper.proc, 
                                       ddl = dipper.ddl, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

# run apparent survival models
dipper_survival_analysis_out <-
  dipper_survival_analysis_run()

# explore all model combinations 
# (caution this is likely to produce false positives due to sample size issues)
dredge_dipper_survival_analysis_run = 
  function(proc_data, design_data){
    # Models exploring variation in apparent survival
    # constant:
    Phi.dot = list(formula =  ~ 1)
    
    # age-dependent:
    Phi.age = list(formula =  ~ age)
    
    # factorial variation across year:
    Phi.year = list(formula =  ~ time)
    
    # additive effects of age and factorial year:
    Phi.age_year = list(formula =  ~ age + time)
    
    # additive effects of age and linear year:
    Phi.age_Time = list(formula =  ~ age + Time)
    
    # interactive effects of age and year:
    Phi.age_x_Time = list(formula =  ~ age * Time)
    
    # Models exploring variation in encounter probability
    # constant:
    p.dot = list(formula =  ~ 1)
    
    # age-dependent:
    p.age2 = list(formula =  ~ age2)
    
    # factorial variation across year:
    p.year = list(formula =  ~ time)
    
    # additive effects of sex and factorial year:
    p.age2_year = list(formula =  ~ age2 + time)
    
    # create a list of candidate models for all the a models above that begin with 
    # either "Phi." or "p."
    cml <-  RMark::create.model.list("CJS")
    
    # specify the data, design matrix, delete unneeded output files, and 
    # run the models in Program MARK
    model.list <-  RMark::mark.wrapper(cml, data = dipper.proc, 
                                       ddl = dipper.ddl, delete = TRUE, 
                                       wrap = FALSE, threads = 1, brief = TRUE,
                                       silent = TRUE, output = FALSE)
    
    # output the model list and sotre the results
    return(model.list)
  }

full_dipper_survival_analysis_out <-
  dredge_dipper_survival_analysis_run()

##### Results ####
# Extract the AIC model table from the model output
# not alot of support for a "significant-effect" of age, but it is slightly
# better than the null model and we have good a priori reasoning that there
# should be an effect
AIC_table_conservative <- 
  dipper_survival_analysis_out$model.table

# Find the model number for the first ranked model of the AIC table
top_mod_num <- 
  as.numeric(rownames(dipper_survival_analysis_out$model.table[1,]))

# extract and format the transformed parameter estimates
estimates <- 
  dipper_survival_analysis_out[[top_mod_num]]$results$real %>% 
  bind_cols(data.frame(str_split_fixed(rownames(.), " ", 
                                       n = 5)), .) %>% 
  dplyr::mutate(age = as.factor(ifelse(unlist(str_extract_all(X2,"[AJ]")) == "A", 
                                "Adult","Juvenile")),
         parameter = as.factor(ifelse(unlist(str_extract_all(X1,"[pP]")) == "p", 
                                      "Phi","p"))) %>% 
  dplyr::select(parameter, age, estimate, lcl, ucl)

clr1 <- viridis(10)[6]
clr2 <- viridis(11)[11]
svg(filename=paste(path,"Survival1.svg",sep=""), height=90/25.4, width=70/25.4, family="Arial", pointsize = 10)
par(mar=c(2.8, 2.8, 0.2, 0.2))
par(mgp=c(1.2,0.2,0))
plot(c(0,1),c(0,1), axes=FALSE, xlab="Altersklasse", ylab="Lokale Uberlebenswahrscheinlichkeit", tcl=-0.25, type="n")
axis(2, at=seq(0,1,by=0.25), labels=c("0.00","0.25","0.50","0.75","1.00"), tcl=-0.25)
axis(1, at=c(0.2,0.8), labels=c("Jungvogel","Altvogel"), tcl=-0.25)
arrows(0.2,estimates$lcl[2],0.2,estimates$ucl[2],code=3,length=0.1,angle=90)
arrows(0.8,estimates$lcl[1],0.8,estimates$ucl[1],code=3,length=0.1,angle=90)
points(c(0.2), estimates$estimate[2], pch=23, bg=clr2, cex=sqrt(42)/2, type="b")
points(c(0.8), estimates$estimate[1], pch=23, bg=clr1, cex=sqrt(33)/2, type="b")
box()
dev.off()

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
  theme(text = element_text(family="Franklin Gothic Book"),
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

# age-specific apparent survival plot
Phi_plot <- 
  ggplot(data = filter(estimates, parameter == "Phi")) +
  geom_errorbar(aes(x = age, ymin = lcl, ymax = ucl), 
               color = "black", size = 0.3, linetype = "solid", width = 0.1) +
  geom_point(aes(x = age, y = estimate, color = age), size = 3) +
  theme_bw() + 
  dipper_theme +
  ylab("Apparent survival (± 95% CI)") +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = "Dark2") +
  annotate(geom = "text", y = 0.95, x = 1.5, angle = 0,
           label = "CJS model of European Dippers\n on the Wurm (2017-2020)", 
           color = "grey30", size = 3) +
  annotate(geom = "text", y = 0.75, x = 1, angle = 0,
           label = "N = 33 individuals", 
           color = "grey30", size = 2.5, fontface = 'italic') +
  annotate(geom = "text", y = 0.55, x = 2, angle = 0,
           label = "N = 46 individuals", 
           color = "grey30", size = 2.5, fontface = 'italic') 

# encounter probability plot
p_plot <- 
  ggplot(data = filter(estimates, parameter == "p")) +
  geom_errorbar(aes(x = year, ymin = lcl, ymax = ucl), 
                color = "black", size = 0.3, 
                linetype = "solid", width = 0.1) +
  geom_point(aes(x = year, y = estimate, color = age), color = "grey20", size = 3) +
  theme_bw() + 
  dipper_theme +
  theme(axis.title.x = element_blank(),
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()) +
  ylab("Encounter probability (± 95% CI)") +
  scale_y_continuous(limits = c(0, 1))

# arrange multi-panelled plot
Wurm_plot <-
  ggarrange(Phi_plot, 
            p_plot,
            ncol = 2, align = "h",
            widths = c(0.65, 0.35))

# export plot to working directory
ggsave(Wurm_plot,
       filename = "dipper_survival_plot.jpeg",
       width = 4,
       height = 4, 
       units = "in")