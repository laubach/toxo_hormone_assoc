###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############       2. Descriptive statsistics fecal data         #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 5 Oct 2020                 #############
#############               last updated: 2 Nov 2020              #############
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
      rm(list = c('plasma_horm_neosp_toxo_data'))
      
      
      
      
###############################################################################
##############               3. Univariate analysis              ##############
###############################################################################
    
  ### 3.1 Univariate stats fecal testosterone (T)
    ## a) Descriptive stats fecal testosterone (T) data
      univar_T_stat <- fec_horm_neosp_toxo_data_6 %>%
        summarise (n.T = sum(!is.na(testosterone.ng.g)),
                   avg.T = round (mean(testosterone.ng.g, 
                                              na.rm = T),2),
                   stdev.T = round (sd(testosterone.ng.g, 
                                              na.rm = T), 2),
                   med.T = round(median(testosterone.ng.g,
                                        na.rm = T), 2),
                   min.T = round(min(testosterone.ng.g,
                                         na.rm = T), 2),
                   max.T = round(max(testosterone.ng.g,
                                     na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_T_stat.pdf'), height = 4, width = 8)
      grid.table(univar_T_stat)
      dev.off()
      
    ## c) Descriptive stats fecal testosterone (T) data grouped by ID
      univar_T_id_stat <- fec_horm_neosp_toxo_data_6 %>%
      # univar_T_id_stat <- fec_horm_neosp_toxo_data_6 %>%
      # univar_T_id_stat <- fec_horm_neosp_data_restrict %>%
        group_by(hy.id) %>%
        summarise (n.T = sum(!is.na(testosterone.ng.g)),
                   avg.T = round (mean(testosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.T = round (sd(testosterone.ng.g, 
                                       na.rm = T), 2),
                   med.T = round(median(testosterone.ng.g,
                                        na.rm = T), 2),
                   min.T = round(min(testosterone.ng.g,
                                     na.rm = T), 2),
                   max.T = round(max(testosterone.ng.g,
                                     na.rm = T), 2)) 
      
    ## d) Histogram fecal testosterone
      hist_plot_T <- ggplot(data=fec_horm_neosp_toxo_data_6, 
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
      fec_horm_neosp_toxo_data_6$testosterone.ng.g.ln <- 
        log(fec_horm_neosp_toxo_data_6$testosterone.ng.g)
      
      fec_horm_toxo_data_restrict$testosterone.ng.g.ln <- 
        log(fec_horm_toxo_data_restrict$testosterone.ng.g)
      
      fec_horm_neosp_data_restrict$testosterone.ng.g.ln <- 
        log(fec_horm_neosp_data_restrict$testosterone.ng.g)
      
    ## f) Histogram of Nat. Log. of fecal testosterone
      hist_plot_ln_T <- ggplot(data=fec_horm_neosp_toxo_data_6, 
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
      univar_cort_stat <- fec_horm_neosp_toxo_data_6 %>%
        summarise (n.cort = sum(!is.na(corticosterone.ng.g)),
                   avg.cort = round (mean(corticosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(corticosterone.ng.g, 
                                       na.rm = T), 2),
                   med.cort = round(median(corticosterone.ng.g,
                                        na.rm = T), 2),
                   min.cort = round(min(corticosterone.ng.g,
                                     na.rm = T), 2),
                   max.cort = round(max(corticosterone.ng.g,
                                     na.rm = T), 2)) 
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_cort_stat.pdf'), height = 4, width = 8)
      grid.table(univar_cort_stat)
      dev.off()
      
    ## c) Descriptive stats fecal corticosterone (cort) data grouped by ID
      univar_cort_id_stat <- fec_horm_neosp_toxo_data_6 %>%
      # univar_cort_stat <- fec_horm_neosp_toxo_data_6 %>%
      # univar_cort_stat <- fec_horm_neosp_data_restrict %>%
        group_by(hy.id) %>%
        summarise (n.cort = sum(!is.na(corticosterone.ng.g)),
                   avg.cort = round (mean(corticosterone.ng.g, 
                                       na.rm = T),2),
                   stdev.cort = round (sd(corticosterone.ng.g, 
                                       na.rm = T), 2),
                   med.cort = round(median(corticosterone.ng.g,
                                           na.rm = T), 2),
                   min.cort = round(min(corticosterone.ng.g,
                                        na.rm = T), 2),
                   max.cort = round(max(corticosterone.ng.g,
                                        na.rm = T), 2)) 
      
    ## d) Histogram fecal corticosterone
      hist_plot_c <- ggplot(data=fec_horm_neosp_toxo_data_6, 
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
      
      print(hist_plot_c)
      
    ## e) Natural log transformation of cort data 
      fec_horm_neosp_toxo_data_6$corticosterone.ng.g.ln <- 
        log(fec_horm_neosp_toxo_data_6$corticosterone.ng.g)
      
      fec_horm_toxo_data_restrict$corticosterone.ng.g.ln <- 
        log(fec_horm_toxo_data_restrict$corticosterone.ng.g)
      
      fec_horm_neosp_data_restrict$corticosterone.ng.g.ln <- 
        log(fec_horm_neosp_data_restrict$corticosterone.ng.g)
      
    ## f) Histogram of Nat. Log. of fecal corticosterone
      hist_plot_ln_c <- ggplot(data=fec_horm_neosp_toxo_data_6, 
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
      
      print(hist_plot_ln_c)
      
    ## g) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_cort.pdf', plot = hist_plot_ln_c, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
      
  ### 3.3 Univariate stats T. gondii 
    ## a) Descriptive stats T. gondii infection prevalence
      univar_toxo_stat <- fec_horm_neosp_toxo_data_6 %>%
        group_by(hy.id, toxo.status) %>%
        summarise(n.status = sum(!is.na(toxo.status))) %>%
        mutate(freq = n.status / sum(n.status, na.rm = T))
      
      univar_toxo_stat <- univar_toxo_stat %>%
        group_by(toxo.status) %>%
        summarise(n.status = sum(!is.na(toxo.status))) %>%
        mutate(freq = n.status / sum(n.status, na.rm = T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_toxo_status.pdf'), height = 4, width = 8)
      grid.table(univar_toxo_stat)
      dev.off()
      
      
  # ### 3.4 Univariate stats N. caninum 
  #   ## a) Descriptive stats N. caninum  infection prevalence
  #     univar_neosp_stat <- fec_horm_neosp_toxo_data_6 %>%
  #       group_by(hy.id, neo.status) %>%
  #       summarise(n.status = sum(!is.na(neo.status))) %>%
  #       mutate(freq = n.status / sum(n.status, na.rm = T))
  #     
  #     univar_neosp_stat <- univar_neosp_stat %>%
  #       group_by(neo.status) %>%
  #       summarise(n.status = sum(!is.na(neo.status))) %>%
  #       mutate(freq = n.status / sum(n.status, na.rm = T))
  #     
  #   ## b) save the data frame of summary stats out as a pdf into output file
  #     pdf(here('output/univar_neosp_status.pdf'), height = 4, width = 8)
  #     grid.table(univar_neosp_stat)
  #     dev.off()
  #     
      

###############################################################################
##############               4. Bivariate analysis               ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate plots testosterone status by sex
    ## a) Scatter plot testosterone by age by sex  
      fec_T_scatter <- ggplot(fec_horm_neosp_toxo_data_6, aes(x = fecal.age.cat, 
                                              y = testosterone.ng.g, 
                                              shape = poop.state, color = sex)) +
        geom_point(position = position_jitter(h=0.0, w=0.3)) +
        geom_hline(yintercept = round(sd(
          fec_horm_neosp_toxo_data_6$testosterone.ng.g, 
                                          na.rm = T), 2) * 4, 
          linetype="solid", color = "red") +
        annotate(geom="text", label = '4 SD', x = 1, 
                 y = round(sd(fec_horm_neosp_toxo_data_6$testosterone.ng.g, 
                   na.rm = T), 2) * 4, vjust=-1) +
        
        geom_hline(yintercept = round(sd(
          fec_horm_neosp_toxo_data_6$testosterone.ng.g, 
          na.rm = T), 2) * 5,
          linetype="dashed", color = "red") +
        annotate(geom="text", label = '5 SD', x = 1, 
                 y = round(sd(fec_horm_neosp_toxo_data_6$testosterone.ng.g, 
                              na.rm = T), 2) * 5, vjust=-1) +
        labs(title = 'Scatter plot of fecal testoserone (repeated measures) 
by age by sex',
             x ='Hyena age class (fecal sampling date)', 
             y ='Fecal testosterone ng/g') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) 
        # add major axes

        print(fec_T_scatter)

    ## b) Save scatter plot
      ggsave('fec_T_scatter.pdf', plot = fec_T_scatter, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
    ## c) Scatter plot corticosterone by age by sex 
      fec_cort_scatter <- ggplot(fec_horm_neosp_toxo_data_6, 
                                 aes(x = fecal.age.cat, 
                                     y = corticosterone.ng.g,
                                     shape = poop.state, color = sex)) +
        geom_point(position = position_jitter(h=0.0, w=0.3)) +
        geom_hline(yintercept = round(sd(
          fec_horm_neosp_toxo_data_6$corticosterone.ng.g, 
          na.rm = T), 2) * 4, 
          linetype="solid", color = "red") +
        annotate(geom="text", label = '4 SD', x = 1, 
                 y = round(sd(fec_horm_neosp_toxo_data_6$corticosterone.ng.g, 
                              na.rm = T), 2) * 4, vjust=-1) +
        
        geom_hline(yintercept = round(sd(
          fec_horm_neosp_toxo_data_6$corticosterone.ng.g, 
          na.rm = T), 2) * 5,
          linetype="dashed", color = "red") +
        annotate(geom="text", label = '5 SD', x = 1, 
                 y = round(sd(fec_horm_neosp_toxo_data_6$corticosterone.ng.g, 
                              na.rm = T), 2) * 5, vjust=-1) +
        labs(title = 'Scatter plot of fecal corticosterone (repeated measures) 
by age by sex',
             x ='Hyena age class (fecal sampling date)', 
             y ='Fecal corticosterone ng/g') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) 
      # add major axes
      
      print(fec_cort_scatter)
      
    ## d) Save scatter plot
      ggsave('fec_cort_scatter.pdf', plot = fec_cort_scatter, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      

#********************** Data Inclusion/Exclusion Criteria **********************
  
  ### 4.2 Remove outliers
    ## a) Create a copy of full data set (including outliers) for sensitivity 
      # analyses
      fec_horm_neosp_toxo_data_6_full <- fec_horm_neosp_toxo_data_6
      
    ## b) Remove testosterone values > 5 SD 
      #set both the untransformed and log transformed values to NA)
      fec_horm_neosp_toxo_data_6 <- fec_horm_neosp_toxo_data_6_full  %>%
        mutate(testosterone.ng.g = ifelse
               ((fec_horm_neosp_toxo_data_6_full$testosterone.ng.g >
                 sd(fec_horm_neosp_toxo_data_6_full$testosterone.ng.g, 
                                        na.rm = T)* 5), 
                 NA, fec_horm_neosp_toxo_data_6_full$testosterone.ng.g)) %>%
        mutate(testosterone.ng.g.ln = ifelse
               ((fec_horm_neosp_toxo_data_6_full$testosterone.ng.g >
                   sd(fec_horm_neosp_toxo_data_6_full$testosterone.ng.g, 
                      na.rm = T)* 5), 
                 NA, fec_horm_neosp_toxo_data_6_full$testosterone.ng.g.ln))
      
    ## c) Remove corticosterone values > 5 SD 
      #set both the untransformed and log transformed values to NA)
      fec_horm_neosp_toxo_data_6 <- fec_horm_neosp_toxo_data_6  %>%
        mutate(corticosterone.ng.g = ifelse
               ((fec_horm_neosp_toxo_data_6$corticosterone.ng.g >
                   sd(fec_horm_neosp_toxo_data_6$corticosterone.ng.g, 
                      na.rm = T)* 5), 
                 NA, fec_horm_neosp_toxo_data_6$corticosterone.ng.g)) %>%
        mutate(corticosterone.ng.g.ln = ifelse
               ((fec_horm_neosp_toxo_data_6$corticosterone.ng.g >
                   sd(fec_horm_neosp_toxo_data_6$corticosterone.ng.g, 
                      na.rm = T)* 5), 
                 NA, fec_horm_toxo_data_restrict$corticosterone.ng.g.ln))
                 
#********************** Data Inclusion/Exclusion Criteria ********************** 

      
  ### 4.3 Descriptive bivariate stats testosterone
    ## a) Testosterone by toxo.status by sex by age
      T_toxo_sex_age_sum <- fec_horm_neosp_toxo_data_6 %>%
     #T_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   sex = first(sex),
                   fecal.age.cat = first(fecal.age.cat)) 
      
      # summarize T averaged over individuals
      T_toxo_sex_age_sum <- T_toxo_sex_age_sum %>%
        group_by (toxo.status, sex, fecal.age.cat) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_toxo_sex_age_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_toxo_sex_age_sum)
      dev.off()
      
    ## c) Testosterone by toxo.status by sex by poop.am.pm
      T_toxo_sex_am_pm_sum <- fec_horm_neosp_toxo_data_6 %>%
        #T_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   sex = first(sex),
                   poop.am.pm = first(poop.am.pm)) 
      
      # summarize T averaged over individuals
      T_toxo_sex_am_pm_sum <- T_toxo_sex_am_pm_sum %>%
        group_by (toxo.status, sex, poop.am.pm) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_toxo_sex_am_pm_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_toxo_sex_am_pm_sum)
      dev.off()
      
    ## e) Testosterone by toxo.status by sex by migratn.seas.fec
      T_toxo_migratn_sum <- fec_horm_neosp_toxo_data_6 %>%
        #T_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   sex = first(sex),
                   migratn.seas.fec = first(migratn.seas.fec)) 
      
      # summarize T averaged over individuals
      T_toxo_migratn_sum <- T_toxo_migratn_sum %>%
        group_by (toxo.status, sex, migratn.seas.fec) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## f) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_toxo_migratn_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_toxo_migratn_sum)
      dev.off()
      
    ## g) Testosterone by toxo.status by sex by hum.pop.poop
      T_toxo_hum_pop_sum <- fec_horm_neosp_toxo_data_6 %>%
        #T_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(testosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   sex = first(sex),
                   hum.pop.poop = first(hum.pop.poop)) 
      
      # summarize T averaged over individuals
      T_toxo_hum_pop_sum <- T_toxo_hum_pop_sum %>%
        group_by (toxo.status, sex, hum.pop.poop) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## h) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_toxo_hum_pop_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_toxo_hum_pop_sum)
      
    ## i) Testosterone by toxo.status by sex by poop.state (feamales)
      T_toxo_state_sum <- fec_horm_neosp_toxo_data_6 %>%
        #T_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(testosterone.ng.g) & sex == 'f') %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(testosterone.ng.g)),
                   sex = first(sex),
                   poop.state = first(poop.state)) 
      
      # summarize T averaged over individuals
      T_toxo_state_sum <- T_toxo_state_sum %>%
        group_by (toxo.status, sex, poop.state) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
      ## h) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/T_toxo_state_sum.pdf'), 
          height = 4, width = 5)
      grid.table(T_toxo_state_sum)
      dev.off()
      
      
  ### 4.3 Descriptive bivariate stats corticosterone
    ## a) Corticosterone by toxo.status by sex by age
      cort_toxo_sex_age_sum <- fec_horm_neosp_toxo_data_6 %>%
        #cort_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(corticosterone.ng.g)),
                   sex = first(sex),
                   fecal.age.cat = first(fecal.age.cat)) 
      
      # summarize T averaged over individuals
      cort_toxo_sex_age_sum <- cort_toxo_sex_age_sum %>%
        group_by (toxo.status, sex, fecal.age.cat) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_toxo_sex_age_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_toxo_sex_age_sum)
      dev.off()
      
    ## c) Corticosterone by toxo.status by sex by poop.am.pm
      cort_toxo_sex_am_pm_sum <- fec_horm_neosp_toxo_data_6 %>%
        #cort_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(corticosterone.ng.g)),
                   sex = first(sex),
                   poop.am.pm = first(poop.am.pm)) 
      
      # summarize T averaged over individuals
      cort_toxo_sex_am_pm_sum <- cort_toxo_sex_am_pm_sum %>%
        group_by (toxo.status, sex, poop.am.pm) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_toxo_sex_am_pm_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_toxo_sex_am_pm_sum)
      dev.off()
      
     ## e) Corticosterone by toxo.status by sex by migratn.seas.fec
      cort_toxo_migratn_sum <- fec_horm_neosp_toxo_data_6 %>%
        #cort_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(corticosterone.ng.g)),
                   sex = first(sex),
                   migratn.seas.fec = first(migratn.seas.fec)) 
      
      # summarize T averaged over individuals
      cort_toxo_migratn_sum <- cort_toxo_migratn_sum %>%
        group_by (toxo.status, sex, migratn.seas.fec) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## f) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_toxo_migratn_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_toxo_migratn_sum)
      dev.off()
      
    ## g) Corticosterone by toxo.status by sex by hum.pop.poop
      cort_toxo_hum_pop_sum <- fec_horm_neosp_toxo_data_6 %>%
        #cort_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(corticosterone.ng.g)) %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(corticosterone.ng.g)),
                   sex = first(sex),
                   hum.pop.poop = first(hum.pop.poop)) 
      
      # summarize T averaged over individuals
      cort_toxo_hum_pop_sum <- cort_toxo_hum_pop_sum %>%
        group_by (toxo.status, sex, hum.pop.poop) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## h) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_toxo_hum_pop_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_toxo_hum_pop_sum)
      
    ## i) Corticosterone by toxo.status by sex by poop.state (feamales)
      cort_toxo_state_sum <- fec_horm_neosp_toxo_data_6 %>%
        #cort_toxo_sum <- fec_horm_toxo_data_restrict %>%
        filter(!is.na(corticosterone.ng.g) & sex == 'f') %>%
        group_by (hy.id, toxo.status) %>%
        summarise (n.Test = sum(!is.na(corticosterone.ng.g)),
                   sex = first(sex),
                   poop.state = first(poop.state)) 
      
      # summarize T averaged over individuals
      cort_toxo_state_sum <- cort_toxo_state_sum %>%
        group_by (toxo.status, sex, poop.state) %>%
        summarise (n.id = sum(!is.na(n.Test)))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## h) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/cort_toxo_state_sum.pdf'), 
          height = 4, width = 5)
      grid.table(cort_toxo_state_sum)
      dev.off() 

      
      
###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
    ### 5.1 Export data to an RData file     
      ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '3_neo_toxo_fec_horm.RData'), 
           list = c('fec_horm_neosp_toxo_data', 'fec_horm_neosp_toxo_data_6',
                    'fec_horm_toxo_data_restrict',
                    'fec_horm_neosp_toxo_data_6_full'
                    #,'fec_horm_neosp_data_restrict'
                    ))
      

     
      
      