###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############       2. Descriptive statsistics plasma data        #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 5 Oct 2020                 #############
#############               last updated: 2 Nov 2020              #############
###############################################################################


  ### PURPOSE: Do univariate and bivariate associations for associations 
             # of toxo. and neosp. with plasma horomones in spotted hyenas
  
  
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
      # load tidyverse packages
        library ('tidyverse')

      # load lubridate packages
        library ('lubridate') 

      # load here packages
        library ('here')
     
    ## b) Graph Plotting and Visualization Packages
      # load ggplot2 packages
      library ('ggplot2')

      # load gridExtra packages
      library ('gridExtra')
      
        
  ### 1.3 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.0.2 (2020-06-22)
    # Platform: x86_64-apple-darwin17.0 (64-bit)
    # Running under: macOS Catalina 10.15.7
    
  
  ### 1.4 Set working directory 
    setwd(here())
  
  
  ### 1.5 Set file paths for data importing and exporting
    ## a) The path to sample, normalized RNA expression, and transformed 
      # behavior data
      project_data_path <- paste0(here('data/'))
     
    ## b) The path for exporting to the output folder
      project_output_path <- paste0(here('output/'))


      
###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  
  
  ### 2.1 Load RData
    ## a) load RData
      load(paste0(project_data_path, '2_tidy_data_neo_toxo_horm.RData'))
      
    ## b) Remove the fecal hormones data tables
      rm(list = c('fec_horm_neosp_toxo_data', 'fec_horm_neosp_toxo_data_4',
               'fec_horm_toxo_data_restrict', 
               'fec_horm_neosp_data_restrict'))
    
      
      
      
###############################################################################
##############               3. Univariate analysis              ##############
###############################################################################
    
  ### 3.1 Univariate stats plasma testosterone (T)
    ## a) Descriptive stats plasma testosterone (T) data
      univar_plasma_T_stat <- plasma_horm_neosp_toxo_data %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                              na.rm = T),2),
                   stdev.T = round (sd(t, 
                                              na.rm = T), 2),
                   med.T = round(median(t,
                                        na.rm = T), 2),
                   min.T = round(min(t,
                                         na.rm = T), 2),
                   max.T = round(max(t,
                                     na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_plasma_T_stat.pdf'), height = 4, width = 8)
      grid.table(univar_plasma_T_stat)
      dev.off()
      
      
    ## c) Histogram plasma testosterone
      hist_plot_plasma_T <- ggplot(data=plasma_horm_neosp_toxo_data , 
                                aes(x=t)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 6, by = 0.25), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0,6)) +
        labs(title = 'Histogram of Plasma Testosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Fec. Testosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_plasma_T)
      
    ## d) Natural log transformation of T data 
      plasma_horm_neosp_toxo_data $t.ln <- 
        log(plasma_horm_neosp_toxo_data $t)
      
    ## e) Histogram of Nat. Log. of plasma testosterone
      hist_plot_ln_plasma_T <- ggplot(data=plasma_horm_neosp_toxo_data , 
                            aes(x=t.ln)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 3, by = 0.125), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0, 3)) +
        labs(title = 'Histogram of Nat. Log. Plasma Testosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Nat. Log. Fec. Testosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_ln_plasma_T)

    ## g) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_T.pdf', plot = hist_plot_ln_T, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
    
  ### 3.2 Univariate stats plasma corticosterone (cort)
    ## a) Descriptive stats plasma corticosterone (cort) data
      univar_plasma_cort_stat <- plasma_horm_neosp_toxo_data  %>%
        summarise (n.cort = sum(!is.na(c )),
                   avg.cort = round (mean(c , 
                                       na.rm = T),2),
                   stdev.cort = round (sd(c , 
                                       na.rm = T), 2),
                   med.cort = round(median(c ,
                                        na.rm = T), 2),
                   min.cort = round(min(c ,
                                     na.rm = T), 2),
                   max.cort = round(max(c ,
                                     na.rm = T), 2)) 
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_plasma_cort_stat.pdf'), height = 4, width = 8)
      grid.table(univar_plasma_cort_stat)
      dev.off()
      
    ## c) Histogram plasma corticosterone
      hist_plot_plasma_c <- ggplot(data=plasma_horm_neosp_toxo_data , 
                            aes(x=c )) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 20, by = 1), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0, 20)) +
        labs(title = 'Histogram of Plasma Corticosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Fec. Corticosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_plasma_c)
      
    ## e) Natural log transformation of cort data 
      plasma_horm_neosp_toxo_data $c.ln <- 
        log(plasma_horm_neosp_toxo_data $c )
      
    ## f) Histogram of Nat. Log. of plasma corticosterone
      hist_plot_ln_plasma_c <- ggplot(data=plasma_horm_neosp_toxo_data , 
                               aes(x=c.ln)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 4, by = 0.125), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0, 4)) +
        labs(title = 'Histogram of Nat. Log. Plasma Corticosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Nat. Log. Fec. Corticosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_ln_plasma_c)
      
    ## g) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_cort.pdf', plot = hist_plot_ln_plasma_c, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
      
  ### 3.3 Univariate stats T. gondii 
    ## a) Descriptive stats T. gondii infection prevalence
      univar_toxo_plasma_stat <- plasma_horm_neosp_toxo_data  %>%
        group_by(toxo.status) %>%
        summarise(n.status = sum(!is.na(toxo.status))) %>%
        mutate(freq = n.status / sum(n.status, na.rm = T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_toxo_plasma_stat.pdf'), height = 4, width = 8)
      grid.table(univar_toxo_plasma_stat)
      dev.off()
      
      
  ### 3.4 Univariate stats N. caninum 
    ## a) Descriptive stats N. caninum  infection prevalence
      univar_neosp_plasma_stat <- plasma_horm_neosp_toxo_data  %>%
        group_by(neo.status) %>%
        summarise(n.status = sum(!is.na(neo.status))) %>%
        mutate(freq = n.status / sum(n.status, na.rm = T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_neosp_plasma_stat.pdf'), height = 4, width = 8)
      grid.table(univar_neosp_plasma_stat)
      dev.off()
      
      

###############################################################################
##############               4. Bivariate analysis               ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate stats testosterone status by sex
    ## a) Sex summary 
      plasma_T_sex_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(t)) %>%
        group_by (sex) %>%
        summarise (n.Test = sum(!is.na(t)),
                   avg.Test = round (mean(t, 
                                       na.rm = T),2),
                   stdev.Test = round (sd(t, 
                                       na.rm = T), 2)) 
      
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/ plasma_T_sex_sum.pdf'), 
          height = 4, width = 5)
      grid.table( plasma_T_sex_sum)
      dev.off() 
      
    ## c) Sex by toxo.status 
      plasma_T_sex_toxo_sum <- plasma_horm_neosp_toxo_data %>%
        group_by (sex, toxo.status) %>%
        summarise (n.id = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                       na.rm = T),2),
                   stdev.T = round (sd(t, 
                                       na.rm = T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
      ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_sex_toxo_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_T_sex_toxo_sum)
      dev.off() 
      
  ### 4.2 Descriptive bivariate stats testosterone status by categorical age
      # when plasma sample was collected
    ## a) Age summary 
      plasma_T_age_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(t)) %>%
        group_by (age.cat.dart) %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                          na.rm = T),2),
                   stdev.T = round (sd(t, 
                                          na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_age_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_T_age_sum)
      dev.off() 
      
    ## c) Age by toxo.status 
      plasma_T_age_toxo_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(t)) %>%
        group_by (age.cat.dart, toxo.status) %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                          na.rm = T),2),
                   stdev.T = round (sd(t, 
                                          na.rm = T), 2), 
                   toxo.status = first(toxo.status))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_age_toxo_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_T_age_toxo_sum)
      dev.off()
      
      
    ## e) Age by sex by toxo.status 
      plasma_T_age_sex_toxo_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(t)) %>%
        group_by (age.cat.dart, sex, toxo.status) %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                       na.rm = T),2),
                   stdev.T = round (sd(t, 
                                       na.rm = T), 2), 
                   toxo.status = first(toxo.status))%>%
        mutate(freq = n.T / sum(n.T))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_age_sex_toxo_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_T_age_sex_toxo_sum)
      dev.off()
      
      
  ### 4.3 Descriptive bivariate stats testosterone status by reproductive state
      # when plasma sample was collected
    ## a) Reproductive state (for females) or males 
      plasma_T_repro_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(t)) %>%
        group_by (dart.state) %>%
        summarise (n.T= sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                          na.rm = T),2),
                   stdev.T = round (sd(t, 
                                          na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_repro_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_T_repro_sum)
      dev.off()     

      
  ### 4.4 Descriptive bivariate stats testosterone status by the season
      # when plasma sample was collected
    ## a) Migration present vs absent summary 
      plasma_T_migrtn_seas_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(t)) %>%
        group_by (migratn.seas.dart) %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                          na.rm = T),2),
                   stdev.T = round (sd(t, 
                                          na.rm = T), 2))%>%
        mutate(freq = n.T / sum(n.T))
      
      ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_migrtn_seas_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_T_migrtn_seas_sum)
      dev.off()   
      

  ### 4.5 Descriptive bivariate stats corticosterone status by sex
    ## a) Sex summary 
      plasma_cort_sex_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (sex) %>%
        summarise (n.Cort = sum(!is.na(c)),
                   avg.Cort = round (mean(c, 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(c, 
                                          na.rm = T), 2))%>%
        mutate(freq = n.Cort / sum(n.Cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_sex_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_sex_sum)
      dev.off() 
      
    ## c) Sex by toxo.status 
      plasma_cort_sex_toxo_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (sex, toxo.status) %>%
        summarise (n.Cort = sum(!is.na(c )),
                   avg.Cort = round (mean(c , 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(c , 
                                          na.rm = T), 2))%>%
        mutate(freq = n.Cort / sum(n.Cort))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_sex_toxo_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_sex_toxo_sum)
      dev.off() 
      
      
  ### 4.6 Descriptive bivariate stats corticosterone status by categorical age
      # when plasma sample was collected
    ## a) Age summary 
      plasma_cort_age_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (age.cat.dart) %>%
        summarise (n.Cort = sum(!is.na(c )),
                   avg.Cort = round (mean(c , 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(c , 
                                          na.rm = T), 2))%>%
        mutate(freq = n.Cort / sum(n.Cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_age_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_age_sum)
      dev.off() 
      
    ## c) Age by toxo.status 
      plasma_cort_age_toxo_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (age.cat.dart, toxo.status) %>%
        summarise (n.cort = sum(!is.na(c )),
                   avg.cort = round (mean(c , 
                                          na.rm = T),2),
                   stdev.cort = round (sd(c , 
                                          na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_age_toxo_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_age_toxo_sum)
      dev.off() 
      
      
    ## e) Age by sex by toxo.status 
      plasma_cort_age_sex_toxo_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (age.cat.dart, sex, toxo.status) %>%
        summarise (n.cort = sum(!is.na(c )),
                   avg.cort = round (mean(c , 
                                          na.rm = T),2),
                   stdev.cort = round (sd(c , 
                                          na.rm = T), 2))%>%
        mutate(freq = n.cort / sum(n.cort))
      
      ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_age_sex_toxo_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_age_sex_toxo_sum)
      dev.off() 
      
  ### 4.7 Descriptive bivariate stats corticosterone status by reproductive 
      # state when plasma sample was collected
    ## a) Reproductive state (for females) or males 
      plasma_cort_repro_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (dart.state) %>%
        summarise (n.Cort = sum(!is.na(c )),
                   avg.Cort = round (mean(c , 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(c , 
                                          na.rm = T), 2))%>%
        mutate(freq = n.Cort / sum(n.Cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_repro_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_repro_sum)
      dev.off()     
      
      
  ### 4.8 Descriptive bivariate stats corticosterone status by the season
      # when plasma sample was collected
    ## a) Migration present vs absent summary 
      plasma_cort_migrtn_seas_sum <- plasma_horm_neosp_toxo_data  %>%
        filter(!is.na(c )) %>%
        group_by (migratn.seas.dart) %>%
        summarise (n.Cort = sum(!is.na(c )),
                   avg.Cort = round (mean(c , 
                                          na.rm = T),2),
                   stdev.Cort = round (sd(c , 
                                          na.rm = T), 2))%>%
        mutate(freq = n.Cort / sum(n.Cort))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_migrtn_seas_sum.pdf'), 
          height = 4, width = 5)
      grid.table(plasma_cort_migrtn_seas_sum)
      dev.off()   
      
      
      
###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
    ### 5.1 Export data to an RData file     
      ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '3_neo_toxo_plasma_horm.RData'), 
           list = c('plasma_horm_neosp_toxo_data '))
      

     
      
      