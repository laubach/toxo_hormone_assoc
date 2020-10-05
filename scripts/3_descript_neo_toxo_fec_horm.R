###############################################################################
#############            Spotted Hyena Neospora caninum:          #############
############# Determinants and behavior and fitness consequences  #############
#############                                                     #############
#############      3. Descriptive stats: Determinants Neosp.      #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 5 Oct 2020                 #############
#############               last updated: 5 Oct 2020              #############
###############################################################################

#**************************  Determinants of Neosp  **************************** 

  ### PURPOSE: Do univariate and bivariate associations for determinants of 
             # of Neospora caninum in spotted hyenas
  
  
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
    ## a) load RData: updated neosp_data joined to hyena data tables
      load(paste0(project_data_path, '2_tidy_data_neo_toxo_fec_horm.RData'))
     
      
      
###############################################################################
##############               3. Univariate analysis              ##############
###############################################################################
    
  ### 3.1 Univariate stats neosp 
    ## a) Descriptive stats neosp status
      univar_neosp_stat <- neosp_data %>%
        group_by(neo.status) %>%
        summarise(n.status = sum(!is.na(neo.status))) %>%
        mutate(freq = n.status / sum(n.status, na.rm = T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_neosp_status.pdf'), height = 4, width = 8)
      grid.table(univar_neosp_stat)
      dev.off()
      
    ## c) Histogram Neospora IFA
      hist_plot_neosp <- ggplot(data=neosp_data, aes(x=ifa.neospora)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1700, by = 100), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0,1700)) +
        labs(title = 'Histogram of Neospora caninum IFA') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='N. caninum IFA', y='Frequency') 
      
      print(hist_plot_neosp)
      
    ## d) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_neosp.pdf', plot = hist_plot_neosp, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
    ## e) Density plot Neospora IFA  
      density_plot_neosp <- ggplot(data=neosp_data, aes(x=ifa.neospora)) + 
        geom_density() +
        xlim(c(0,1700)) +
        labs(title = 'Density plot of Neospora caninum IFA') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='N. caninum IFA', y='Frequency') 
      
      print(density_plot_neosp)
      
    ## f) Save density plot
      # use ggsave to save the plot
      ggsave('density_plot_neosp.pdf', plot = density_plot_neosp, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE) 

    ## g) Histogram Neospora IFA by infection status
      hist_plot_neosp_by_status <- 
        ggplot(data=neosp_data, aes(x=ifa.neospora, 
                                          fill = neo.status)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 1700, by = 100), 
                       # col='black',
                       # fill = 'dark grey',
                       position='dodge') +
        xlim(c(0,1700)) +
        labs(title = 'Histogram of Neospora caninum IFA by infection status') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='N. caninum IFA', y='Frequency') 
      
      print(hist_plot_neosp_by_status)
      
    ## h) Save histogram plot
      # use ggsave to save the plot
      ggsave('hist_plot_neosp_by_status.pdf', plot = hist_plot_neosp_by_status, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      

###############################################################################
##############               4. Bivariate analysis               ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate stats neosp status by sex
    ## a) Sex ratio summary 
      neosp_sex_ratio_sum <- neosp_data %>%
        #group_by (id) %>%
        group_by (sex, neo.status) %>%
        summarise(n=n_distinct(hy.id)) %>%
        mutate(freq = n / sum(n))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_sex_ratio_sum.pdf'), 
          height = 3, width = 5)
      grid.table(neosp_sex_ratio_sum)
      dev.off() 
      
      
  ### 4.2 Descriptive bivariate stats neosp by age.mon.dart    
    ## a) Age summary  
      neosp_age_var_sum <- neosp_data %>%
        group_by (age.cat.dart) %>%
        summarise (n.age.dart = sum(!is.na(age.mon.dart)),
                   avg.age.dart = round (mean(age.mon.dart, na.rm = T),2),
                   stdev.age.dart = round (sd(age.mon.dart, na.rm = T), 2)) %>%
        mutate (freq_age =  (n.age.dart / sum(n.age.dart)))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_age_var_sum.pdf'), 
          height = 4, width = 8)
      grid.table(neosp_age_var_sum)
      dev.off() 
      
    ## c) Age categories summary 
      neosp_age_cat_sum <- neosp_data %>%
        #group_by (id) %>%
        group_by (age.cat.dart, neo.status) %>%
        summarise(n=n_distinct(hy.id)) %>%
        mutate(freq = n / sum(n))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_age_cat_sum.pdf'), 
          height = 3, width = 5)
      grid.table(neosp_age_cat_sum)
      dev.off() 
      
    ## e) Plot age at darting by neosp status
      age_dart_by_neosp_plot <- ggplot(data = neosp_data,
             aes (x = neo.status, y = age.mon.dart,
                  color = neo.status)) +
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs (title = 'Age at diagnosis by Neosp. 
infection status') +
        ylab ('Age at diagnosis (mon)') +
        xlab ('Neosp status (0 = uninfected)')
      
      print(age_dart_by_neosp_plot)
      
    ## f) Save Plot
      # use ggsave to save the boxplot
      ggsave('age_dart_by_neosp_plot.pdf', plot = age_dart_by_neosp_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 7, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
  ### 4.3 Descriptive bivariate stats neosp by rank
    ## a) Neosp. by rank summary, among females and males separately
      # separate sexes because rank calculate separately
      neosp_rank_sex_sum <- neosp_data %>%
        group_by  (sex, neo.status) %>%
        summarise (n.rank =  sum(!is.na(stan.rank.dart)),
                   #n_rank = n(), # n including na
                   avg.rank = round (mean(stan.rank.dart, na.rm = T),2),
                   stdev.rank = round (sd(stan.rank.dart, na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_rank_sex_sum.pdf'),
          height = 4, width = 8)
      grid.table(neosp_rank_sex_sum)
      dev.off()
      
    ## c) Neosp. by rank summary, among females and males separately by age
      neosp_rank_sex_age_sum <- neosp_data %>%
        group_by  (sex, age.cat.dart, neo.status) %>%
        summarise (n.rank =  sum(!is.na(stan.rank.dart)),
                   #n_rank = n(), # n including na
                   avg.rank = round (mean(stan.rank.dart, na.rm = T),2),
                   stdev.rank = round (sd(stan.rank.dart, na.rm = T), 2))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_rank_sex_age_sum.pdf'),
          height = 4, width = 8)
      grid.table(neosp_rank_sex_age_sum)
      dev.off()
      
    ## e) Plot standardize rank by neosp. by sex 
      neosp_rank_sex_plot <- ggplot(data = neosp_data,
                                        aes (x = sex, y = stan.rank.dart,
                                             color = neo.status)) +
        geom_boxplot() +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs (title = 'Rank by neosp. infection status by sex') +
        ylab ('Standardized rank') +
        xlab ('Neosp status (0 = uninfected)')
      
      print(neosp_rank_sex_plot)
      
    ## f) Save Plot
      # use ggsave to save the boxplot
      ggsave('neosp_rank_sex_plot.pdf', plot = neosp_rank_sex_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 7, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE)
      
    ## g) Plot standardize rank by neosp. by sex and by age
      neosp_rank_sex_age_plot <- ggplot(data = neosp_data,
                                       aes (x = sex, y = stan.rank.dart,
                                            color = neo.status)) +
        geom_boxplot() +
        facet_wrap(~age.cat.dart) +
        theme(text = element_text(size=20))+
        scale_colour_hue(l = 50) + # Use a slightly darker palette than normal
        labs (title = 'Rank by neosp. infection status
by sex by age') +
        ylab ('Standardized rank') +
        xlab ('Neosp status (0 = uninfected)')
      
      print(neosp_rank_sex_age_plot)
      
    ## h) Save Plot
      # use ggsave to save the boxplot
      ggsave('neosp_rank_sex_age_plot.pdf', plot = neosp_rank_sex_age_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 7, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      
  ### 4.4 Descriptive bivariate stats Neosp by human pop size
    ## a) Neosp. by human pop size (proxy) summary 
      neosp_hum_pop_sum <- neosp_data %>%
        #group_by (id) %>%
        group_by (hum.pop.den, neo.status) %>%
        summarise(n=n_distinct(hy.id)) %>%
        mutate(freq = n / sum(n))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_hum_pop_sum.pdf'), 
          height = 3, width = 5)
      grid.table(neosp_hum_pop_sum)
      dev.off() 
      
    ## c) Neosp. by human pop size (proxy) summary, among females and males 
      # separately by age
      neosp_hum_pop_sex_age_sum <- neosp_data %>%
        group_by  (sex, age.cat.dart, hum.pop.den, neo.status) %>%
        summarise(n=n_distinct(hy.id)) %>%
        mutate(freq = n / sum(n))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_hum_pop_sex_age_sum.pdf'),
          height = 8, width = 8)
      grid.table(neosp_hum_pop_sex_age_sum)
      dev.off()
      
    ## e) Neosp. by human pop size (proxy) summary, by sex in subadults 
        # and adult hyeans
      neosp_hum_pop_sex_subs_adlt_sum <- neosp_data %>%
        filter(age.cat.dart == 'subadult' | age.cat.dart == 'adult') %>%
        group_by (sex, hum.pop.den, neo.status) %>%
        summarise(n=n_distinct(hy.id)) %>%
        mutate(freq = n / sum(n))
      
    ## f) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/neosp_hum_pop_sex_subs_adlt_sum.pdf'),
          height = 7, width = 8)
      grid.table(neosp_hum_pop_sex_subs_adlt_sum)
      dev.off()
      
    
            
      
      

      
      