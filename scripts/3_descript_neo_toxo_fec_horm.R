###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############             3. Descriptive statsistics              #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 5 Oct 2020                 #############
#############               last updated: 6 Oct 2020              #############
###############################################################################


  ### PURPOSE: Do univariate and bivariate associations for associations 
             # of toxo. and neosp. with fecal horomones in spotted hyenas
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Univariate analysis 
    # 4: Bivariate analysis



###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 Global options
    ## a) clear global environment
      rm(list = ls())

    ## b) prevent R from automatically reading charater strins as factors
      options(stringsAsFactors = FALSE)
  

  ### 1.2 Install and load CRAN packages   
    ## a) Data Manipulation and Descriptive Stats Packages
      # Check for tidyverse and install if not already installed
        if (!'tidyverse' %in% installed.packages()[,1]){
          install.packages ('tidyyverse')
        }
      # load tidyverse packages
        library ('tidyverse')
      
      # Check for lubridate and install if not already installed
        if (!'lubridate' %in% installed.packages()[,1]){
          install.packages ('lubridate')
        }
      # load lubridate packages
        library ('lubridate') 
        
      # Check for here and install if not already installed
        if (!'here' %in% installed.packages()[,1]){
          install.packages ('here')
        }
      # load here packages
        library ('here')
     
    ## b) Graph Plotting and Visualization Packages
      # Check for ggplot2 and install if not already installed
      if (!'ggplot2' %in% installed.packages()[,1]){
        install.packages ('ggplot2')
      }
      # load ggplot2 packages
      library ('ggplot2')
      
      # Check for gridExtra and install if not already installed
      if (!'gridExtra' %in% installed.packages()[,1]){
        install.packages ('gridExtra')
      }
      # load gridExtra packages
      library ('gridExtra')
      
        
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.0.0 (2020-04-24)
    # Platform: x86_64-apple-darwin17.0 (64-bit)
    # Running under: macOS Mojave 10.14.6
    
  
  ### 1.4 Set working directory 
    setwd(here())
  
  
  ### 1.5 Set file paths for data importing and exporting
    ## a) The path to sample, normalized RNA expression, and transformed 
      # behavior data
      project_data_path <- paste0(here('data/'))
     
    ## b) The path for exporting to the output folder
      project_output_path <- paste0(here('output/'))
    
    ## c) Source scripts path
      source_path <- paste("~/Git/source_code/")
      
      
  ### 1.6 Source functions
    ## a) all_char_to_lower function
      source(file = paste0(source_path, "all_char_to_lower.R"))
      
    ## b) format_var_names function
      source(file = paste0(source_path, "format_var_names.R"))
      
    ## c) format_var_names_dash function
      source(file = paste0(source_path, "format_var_names_dash.R"))  


      
###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  
  
  ### 2.1 Load RData
    ## a) load RData: updated fec_horm_neosp_toxo_data joined to hyena 
      # data tables
      load(paste0(project_data_path, '2_tidy_data_neo_toxo_fec_horm.RData'))
     
      
      
###############################################################################
##############               3. Univariate analysis              ##############
###############################################################################
    
  ### 3.1 Univariate stats fecal testosterone (T)
    ## a) Descriptive stats fecal testosterone (T) data
      univar_T_stat <- fec_horm_neosp_toxo_data_12 %>%
        summarise (n.T = sum(!is.na(testosterone.ng.g)),
                   avg.T = round (mean(testosterone.ng.g, 
                                              na.rm = T),2),
                   stdev.T = round (sd(testosterone.ng.g, 
                                              na.rm = T), 2)) 
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_T_stat.pdf'), height = 4, width = 8)
      grid.table(univar_T_stat)
      dev.off()
      
    ## c) Descriptive stats fecal testosterone (T) data grouped by ID
      univar_T_id_stat <- fec_horm_neosp_toxo_data_12 %>%
        group_by(hy.id) %>%
        summarise (n.T = sum(!is.na(testosterone.ng.g)),
                   avg.T = round (mean(testosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.T = round (sd(testosterone.ng.g, 
                                       na.rm = T), 2)) 
      
    ## d) Histogram fecal testosterone
      hist_plot_T <- ggplot(data=fec_horm_neosp_toxo_data_12, 
                                aes(x=testosterone.ng.g)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 600, by = 5), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0,600)) +
        labs(title = 'Histogram of Fecal Testosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Fec. Testosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_T)
      
    ## e) Natural log transformation of T data 
      fec_horm_neosp_toxo_data_12$testosterone.ng.g.ln <- 
        log(fec_horm_neosp_toxo_data_12$testosterone.ng.g)
      
    ## f) Histogram of Nat. Log. of fecal testosterone
      hist_plot_ln_T <- ggplot(data=fec_horm_neosp_toxo_data_12, 
                            aes(x=testosterone.ng.g.ln)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 7.5, by = 0.5), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0,7.5)) +
        labs(title = 'Histogram of Nat. Log. Fecal Testosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Nat. Log. Fec. Testosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_ln_T)

    ## g) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_T.pdf', plot = hist_plot_ln_T, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
    
  ### 3.2 Univariate stats fecal corticosterone (cort)
    ## a) Descriptive stats fecal corticosterone (cort) data
      univar_cort_stat <- fec_horm_neosp_toxo_data_12 %>%
        summarise (n.cort = sum(!is.na(corticosterone.ng.g)),
                   avg.cort = round (mean(corticosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(corticosterone.ng.g, 
                                       na.rm = T), 2)) 
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_cort_stat.pdf'), height = 4, width = 8)
      grid.table(univar_cort_stat)
      dev.off()
      
    ## c) Descriptive stats fecal corticosterone (cort) data grouped by ID
      univar_cort_id_stat <- fec_horm_neosp_toxo_data_12 %>%
        group_by(hy.id) %>%
        summarise (n.cort = sum(!is.na(corticosterone.ng.g)),
                   avg.cort = round (mean(corticosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(corticosterone.ng.g, 
                                       na.rm = T), 2)) 
      
    ## d) Histogram fecal corticosterone
      hist_plot_T <- ggplot(data=fec_horm_neosp_toxo_data_12, 
                            aes(x=corticosterone.ng.g)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1000, by = 5), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0, 1000)) +
        labs(title = 'Histogram of Fecal Corticosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Fec. Corticosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_T)
      
    ## e) Natural log transformation of T data 
      fec_horm_neosp_toxo_data_12$corticosterone.ng.g.ln <- 
        log(fec_horm_neosp_toxo_data_12$corticosterone.ng.g)
      
    ## f) Histogram of Nat. Log. of fecal corticosterone
      hist_plot_ln_T <- ggplot(data=fec_horm_neosp_toxo_data_12, 
                               aes(x=corticosterone.ng.g.ln)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 7.5, by = 0.5), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0,7.5)) +
        labs(title = 'Histogram of Nat. Log. Fecal Corticosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Nat. Log. Fec. Corticosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_ln_T)
      
    ## g) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_cort.pdf', plot = hist_plot_ln_T, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
  

###############################################################################
##############               4. Bivariate analysis               ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate stats testosterone status by sex
    ## a) Sex ratio summary 
      T_sex_ratio_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, sex) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                       na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_sex_ratio_sum <- T_sex_ratio_sum %>%
        group_by (sex) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_sex_ratio_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_sex_ratio_sum)
      dev.off() 
      
      
  ### 4.2 Descriptive bivariate stats testosterone status by categorical age
      # when fecal sample was collected
    ## a) Age summary 
      T_age_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, fecal.age.cat) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_age_sum <- T_age_sum %>%
        group_by (fecal.age.cat) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_age_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_age_sum)
      dev.off() 
      
    ## c) Sex by Age summary 
      T_sex_age_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, sex, fecal.age.cat) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_sex_age_sum <- T_sex_age_sum %>%
        group_by (sex, fecal.age.cat) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
  
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_sex_age_sum.pdf'), 
          height = 5, width = 5)
      grid.table(T_sex_age_sum)
      dev.off() 
      
      
  ### 4.3 Descriptive bivariate stats testosterone status by reproductive state
      # when fecal sample was collected
    ## a) Reproductive state (for females) or males 
      T_repro_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, state) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_repro_sum <- T_repro_sum %>%
        group_by (state) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_repro_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_repro_sum)
      dev.off()     

      
  ### 4.4 Descriptive bivariate stats testosterone status by time of day
      # when fecal sample was collected
    ## a) AM vs PM summary 
      T_time_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, poop.am.pm) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_time_sum <- T_time_sum %>%
        group_by (poop.am.pm) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_time_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_time_sum)
      dev.off()     
      
    
  ### 4.5 Descriptive bivariate stats testosterone status by the season
      # when fecal sample was collected
    ## a) Migration present vs absent summary 
      T_migrtn_seas_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, migratn.seas.fec) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_migrtn_seas_sum <- T_migrtn_seas_sum %>%
        group_by (migratn.seas.fec) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_migrtn_seas_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_migrtn_seas_sum)
      dev.off()   
      

  ### 4.6 Descriptive bivariate stats testosterone status by human disturbance
      # when fecal sample was collected
    ## a) Human disturbance summary 
      T_hum_disturb_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, hum.pop.den) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   avg.Test = round (mean(testosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Test = round (sd(testosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      T_hum_disturb_sum <- T_hum_disturb_sum %>%
        group_by (hum.pop.den) %>%
        summarise (n.T = sum(!is.na(avg.Test)),
                   avg.T = round (mean(avg.Test, 
                                       na.rm = T),2),
                   stdev.T = round (sd(avg.Test, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_hum_disturb_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_hum_disturb_sum)
      dev.off()
      
      
  ### 4.7 Descriptive bivariate stats corticosterone status by sex
    ## a) Sex ratio summary 
      cort_sex_ratio_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, sex) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_sex_ratio_sum <- cort_sex_ratio_sum %>%
        group_by (sex) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_sex_ratio_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_sex_ratio_sum)
      dev.off() 
      
      
  ### 4.8 Descriptive bivariate stats corticosterone status by categorical age
      # when fecal sample was collected
    ## a) Age summary 
      cort_age_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, fecal.age.cat) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_age_sum <- cort_age_sum %>%
        group_by (fecal.age.cat) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_age_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_age_sum)
      dev.off() 
      
    ## c) Sex by Age summary 
      cort_sex_age_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, sex, fecal.age.cat) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_sex_age_sum <- cort_sex_age_sum %>%
        group_by (sex, fecal.age.cat) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_sex_age_sum.pdf'), 
          height = 5, width = 5)
      grid.table(cort_sex_age_sum)
      dev.off() 
      
      
  ### 4.9 Descriptive bivariate stats corticosterone status by reproductive state
      # when fecal sample was collected
    ## a) Reproductive state (for females) or males 
      cort_repro_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, state) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_repro_sum <- cort_repro_sum %>%
        group_by (state) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_repro_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_repro_sum)
      dev.off()     
      
      
  ### 4.10 Descriptive bivariate stats corticosterone status by time of day
      # when fecal sample was collected
    ## a) AM vs PM summary 
      cort_time_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, poop.am.pm) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_time_sum <- cort_time_sum %>%
        group_by (poop.am.pm) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_time_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_time_sum)
      dev.off()     
      
      
  ### 4.11 Descriptive bivariate stats corticosterone status by the season
      # when fecal sample was collected
    ## a) Migration present vs absent summary 
      cort_migrtn_seas_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, migratn.seas.fec) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_migrtn_seas_sum <- cort_migrtn_seas_sum %>%
        group_by (migratn.seas.fec) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_migrtn_seas_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_migrtn_seas_sum)
      dev.off()   
      
      
  ### 4.12 Descriptive bivariate stats corticosterone status by human disturbance
      # when fecal sample was collected
    ## a) Human disturbance summary 
      cort_hum_disturb_sum <- fec_horm_neosp_toxo_data_12 %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, hum.pop.den) %>%
        summarise (n.Cort = sum(!is.na(corticosterone.ng.g)),
                   avg.Cort = round (mean(corticosterone.ng.g, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(corticosterone.ng.g, 
                                          na.rm = T), 2)) 
      
      # summarize T averaged over individuals
      cort_hum_disturb_sum <- cort_hum_disturb_sum %>%
        group_by (hum.pop.den) %>%
        summarise (n.cort = sum(!is.na(avg.Cort)),
                   avg.cort = round (mean(avg.Cort, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(avg.Cort, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_hum_disturb_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_hum_disturb_sum)
      dev.off()
      

     
      
      