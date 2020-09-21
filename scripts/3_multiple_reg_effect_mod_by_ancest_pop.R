###############################################################################
##############    Guppy Predation Environment, RNA expression,   ##############
##############     and Behavior: Effect modification by Group    ##############
##############                 By: Zach Laubach                  ##############
##############               created: 28 July 2020               ##############
##############            last updated: 28 July 2020             ##############
###############################################################################


  ### PURPOSE: 
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Univariate analyses
    # 4: Effect modification models 
    # 5: Export data files
  



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
      
      # Check for forcats and install if not already installed
      if (!'forcats' %in% installed.packages()[,1]){
        install.packages ('forcats')
      }
      # load forcats packages
      library ('forcats')
        
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
      
      # Check for hrbrthemes and install if not already installed
      if (!'hrbrthemes' %in% installed.packages()[,1]){
        install.packages ('hrbrthemes')
      }
      # load hrbrthemes packages
      library ('hrbrthemes')
      
      # Check for viridis and install if not already installed
      if (!'viridis' %in% installed.packages()[,1]){
        install.packages ('viridis')
      }
      # load viridis packages
      library ('viridis')
      
      # Check for EnhancedVolcano and install if not already installed
      # used with EWAS data to make volcano plots
      if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')
      BiocManager::install('EnhancedVolcano')

      # load EnhancedVolcano package
      library ('EnhancedVolcano')
     
      # Check for dotwhisker and install if not already installed
      # used with broom to graph beta estimates
      if (!'dotwhisker' %in% installed.packages()[,1]){
        install.packages ('dotwhisker')
      }
      # load dotwhisker packages
      library ('dotwhisker')
      
      # Check for ggcorrplot and install if not already installed
      if (!'ggcorrplot' %in% installed.packages()[,1]){
        install.packages ('ggcorrplot')
      }
      # load ggcorrplot packages
      library ('ggcorrplot')
      
    ## c) Modeling Packages
      # Check for broom and install if not already installed
      if (!'broom' %in% installed.packages()[,1]){
        install.packages ('broom')
      }
      # load broom packages
      library ('broom')
      
      # install broom.mixed, a package under development by Ben Bolker
      # similar to broom but extends bey lm and glm
      install.packages("remotes")
      remotes::install_github("bbolker/broom.mixed")
      
      # # Check for nlme and install if not already installed
      # if (!'nlme' %in% installed.packages()[,1]){
      #   install.packages ('nlme')
      # }
      # # load nlme packages
      # library ('nlme')
      # 
      # # Check for lme4 and install if not already installed
      # if (!'lme4' %in% installed.packages()[,1]){
      #   install.packages ('lme4')
      # }
      # # load lme4 packages
      # library ('lme4')
      # 
      # # Check for bbmle and install if not already installed
      # # for AICtab
      # if (!'bbmle' %in% installed.packages()[,1]){
      #   install.packages ('bbmle')
      # }
      # # load bbmle packages
      # library ('bbmle')
      # 
      # # Check for car and install if not already installed
      # # includes vif function
      # if (!'car' %in% installed.packages()[,1]){
      #   install.packages ('car')
      # }
      # # load car packages
      # library ('car')
      
      # Check for qvalue and install if not already installed
      if (!'qvalue' %in% installed.packages()[,1]){
        install.packages ('qvalue')
      }
      # load qvalue package
      library ('qvalue')
      
        
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
    ## a) The path to sample, normalized RNA expression, and 
      # transformed behavior data
      project_data_path <- paste0(here('data/'))
     
    ## b) The path for exporting to the output folder
      project_output_path <- paste0(here('output/'))
    
    ## c) Source scripts path
      source_path <- paste("~/Git/source_code/")


      
###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  
  
  ### 2.1 Load RData
    ## a) Load RData (tables on guppy samples, RNA expression, and beahvior)
      load(paste0(project_data_path,'tidy_data_pred_rna_behav.RData'))
     
      
      
###############################################################################
##############               3. Univariate analyses              ##############
###############################################################################

  ### 3.1 Calculate basic descriptive statistics for each behavior 
      # outcome varialbe
      
    ## a) Subset aripo behavior variables and transform data from wide to long
      aripo_behav_long <- aripo_data %>%
        select(c(fish, sigtime:moving)) %>%
        gather(key = 'behav', value = 'transform', 'sigtime':'moving', 
               factor_key = T)      
      
    ## b) Summarize aripo_data: calculate the average, sd, etc. for each behavior
      summry_aripo_behav <- aripo_behav_long  %>%
        group_by(behav) %>%
        summarize(n.aripo = sum(!is.na(transform)),
                  avg.aripo = round(mean(transform, 
                                     na.rm = T), 3),
                  sd.aripo = round(sd(transform, 
                                  na.rm = T), 3),
                  med.aripo = round(median(transform, 
                                  na.rm = T), 3),
                  min.aripo = round(min(transform, 
                                     na.rm = T), 3),
                  max.aripo = round(max(transform, 
                                  na.rm = T), 3))
      
    ## c) Save the summary data table
      pdf(paste0(here(),"/output/summry_aripo_behav.pdf"),
          height = 5, width = 7)
      grid.table(summry_aripo_behav)
      dev.off()
      
    ## d) Subset quare behavior variables
      quare_behav_long <- quare_data %>%
        select(c(fish, sigtime:moving)) %>%
        gather(key = 'behav', value = 'transform', 'sigtime':'moving', 
               factor_key = T)       
      
    ## e) Summarize quare_data: calculate the average, sd, etc. for each behavior
      summry_quare_behav <- quare_behav_long  %>%
        group_by(behav) %>%
        summarize(n.quare = sum(!is.na(transform)),
                  avg.quare = round(mean(transform, 
                                       na.rm = T), 3),
                  sd.quare = round(sd(transform, 
                                    na.rm = T), 3),
                  med.quare = round(median(transform, 
                                         na.rm = T), 3),
                  min.quare = round(min(transform, 
                                      na.rm = T), 3),
                  max.quare = round(max(transform, 
                                      na.rm = T), 3))
      
    ## f) Save the summary data table
      pdf(paste0(here(),"/output/summry_quare_behav.pdf"),
          height = 5, width = 7)
      grid.table(summry_quare_behav)
      dev.off()
      
      
  ### 3.2 Plot Histogram Distributions
    ## a) Facet plot of Aripo behavior histograms
      aripo_hists <- aripo_behav_long %>%
        mutate(behav = fct_reorder(behav, transform)) %>%
        ggplot(aes(x=transform, color=behav, fill=behav)) +
        geom_histogram(alpha=0.8, binwidth = 0.5) +
        scale_fill_viridis(discrete=TRUE) +
        scale_color_viridis(discrete=TRUE) +
        theme_ipsum(base_family = "sans") +
        theme(
          legend.position='none',
          panel.spacing = unit(1, 'lines'),
          strip.text.x = element_text(size = 8)
        ) +
        xlab('') +
        ylab('Assigned Probability (%)') +
        facet_wrap(~behav)
      
    ## b) Print aripo_hists plots
      print(aripo_hists)
      
    ## c) Save aripo_hists plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_hists.pdf", plot = aripo_hists, device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)     
      
    ## d) Facet plot of Quare behavior histograms
      quare_hists <- quare_behav_long %>%
        mutate(behav = fct_reorder(behav, transform)) %>%
        ggplot(aes(x=transform, color=behav, fill=behav)) +
        geom_histogram(alpha=0.8, binwidth = 0.5) +
        scale_fill_viridis(discrete=TRUE) +
        scale_color_viridis(discrete=TRUE) +
        theme_ipsum(base_family = "sans") +
        theme(
          legend.position='none',
          panel.spacing = unit(1, 'lines'),
          strip.text.x = element_text(size = 8)
        ) +
        xlab('') +
        ylab('Assigned Probability (%)') +
        facet_wrap(~behav)
      
    ## e) Print quare_hists plots
      print(quare_hists)
      
    ## f) Save aripo_hists plots
      # use ggsave to save the plot as pdf
      ggsave("quare_hists.pdf", plot = quare_hists, device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)  
      
      
  ### 3.3 Correlation matrix of behaviors
    ## a) Make a correlation matrix of Aripo behaviors
      aripo_behav_corr<- aripo_data %>%
        select(c(sigtime:moving)) %>%
        do(as.data.frame(cor(., method = 'spearman', 
                             use="pairwise.complete.obs")))
      
    ## b) Graph the correlation matrix of Aripo behaviors
      aripo_behav_corr_plot <-ggcorrplot(aripo_behav_corr, method = 'circle',
                                         type = 'upper') + 
        labs(title = "Aripo behaviors correlation plot") +
        theme(plot.title = element_text(hjust = 0.5)) # center title
      
      print(aripo_behav_corr_plot)
    
    ## c) Save aripo_hists plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_behav_corr_plot.pdf", plot = aripo_behav_corr_plot, 
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
      
    ## d) Make a correlation matrix of Quare behaviors
      quare_behav_corr<- quare_data %>%
        select(c(sigtime:moving)) %>%
        do(as.data.frame(cor(., method = 'spearman', 
                             use="pairwise.complete.obs")))
      
    ## e) Graph the correlation matrix of Quare behaviors
      quare_behav_corr_plot <-ggcorrplot(quare_behav_corr, method = 'circle',
                                         type = 'upper') + 
        labs(title = "Quare behaviors correlation plot") +
        theme(plot.title = element_text(hjust = 0.5)) # center title
      
      print(quare_behav_corr_plot)
      
    ## f) Save quare_hists plots
      # use ggsave to save the plot as pdf
      ggsave("quare_behav_corr_plot.pdf", plot = quare_behav_corr_plot, 
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      

      
      
###############################################################################
##############           4. Effect modification models           ##############
###############################################################################      
      
  ### 4.1 Run Aripo Effect Modification Models  
    ## a) Identify the first and last columns containing RNA expression IDs 
      # in aripo_data
      colnames(aripo_data)
      tail(colnames(aripo_data))
     
    ## b) Transform aripo_data data from wide to long based on RNA expression
      aripo_rna_long <- aripo_data %>%
        gather(key = 'site', value = 'rna', 
               'ENSPREG00000006268':'ENSPREG00000007674', 
               factor_key = F)
      
    ## c) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'contact' - outcome and tidy lm objects using 'broom'
      aripo_contact_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(contact ~ rna*group, data = .)))
        
              
    ## d) Check if there are 8 estimates (intercept, rna, group & intx) per model
      summry_aripo_contact_coef <- aripo_contact_coef  %>%
        group_by(site) %>%
        summarize(n.terms = sum(!is.na(term)))
      
    ## e) filter to remove intercepts and retain only intx coefs and test stats
      aripo_contact_coef <- filter(aripo_contact_coef, 
                                    grepl('rna:', term)) 
      
    ## f) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'afaceoff' - outcome and tidy lm objects using 'broom'
      aripo_afaceoff_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(afaceoff ~ rna*group, data = .)))
      
    ## g) filter to remove intercepts and retain only intx coefs and test stats
      aripo_afaceoff_coef <- filter(aripo_afaceoff_coef, 
                                    grepl('rna:', term))
      
    ## h) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'acontact' - outcome and tidy lm objects using 'broom'
      aripo_acontact_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(acontact ~ rna*group, data = .)))
      
    ## i) filter to remove intercepts and retain only intx coefs and test stats
      aripo_acontact_coef <- filter(aripo_acontact_coef, 
                                    grepl('rna:', term)) 
      
    ## j) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'posturing' - outcome and tidy lm objects using 'broom'
      aripo_posturing_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(posturing ~ rna*group, data = .)))
      
    ## k) filter to remove intercepts and retain only intx coefs and test stats
      aripo_posturing_coef <- filter(aripo_posturing_coef, 
                                    grepl('rna:', term)) 
      
    ## l) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'thrust' - outcome and tidy lm objects using 'broom'
      aripo_thrust_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(thrust ~ rna*group, data = .)))
      
    ## m) filter to remove intercepts and retain only intx coefs and test stats
      aripo_thrust_coef <- filter(aripo_thrust_coef, 
                                    grepl('rna:', term)) 
      
    ## n) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'alunge' - outcome and tidy lm objects using 'broom'
      aripo_alunge_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(alunge ~ rna*group, data = .)))
      
    ## o) filter to remove intercepts and retain only intx coefs and test stats
      aripo_alunge_coef <- filter(aripo_alunge_coef, 
                                    grepl('rna:', term)) 
      
    ## p) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'athrust' - outcome and tidy lm objects using 'broom'
      aripo_athrust_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(athrust ~ rna*group, data = .)))
      
    ## q) filter to remove intercepts and retain only intx coefs and test stats
      aripo_athrust_coef <- filter(aripo_athrust_coef, 
                                    grepl('rna:', term)) 
      
    ## r) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'asigmoid' - outcome and tidy lm objects using 'broom'
      aripo_asigmoid_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(asigmoid ~ rna*group, data = .)))
      
    ## s) filter to remove intercepts and retain only intx coefs and test stats
      aripo_asigmoid_coef <- filter(aripo_asigmoid_coef, 
                                    grepl('rna:', term)) 
      
    ## t) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'surface' - outcome and tidy lm objects using 'broom'
      aripo_surface_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(surface ~ rna*group, data = .)))
      
    ## u) filter to remove intercepts and retain only intx coefs and test stats
      aripo_surface_coef <- filter(aripo_surface_coef, 
                                    grepl('rna:', term)) 
      
    ## v) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'sigmoid' - outcome and tidy lm objects using 'broom'
      aripo_sigmoid_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(sigmoid ~ rna*group, data = .)))
      
    ## w) filter to remove intercepts and retain only intx coefs and test stats
      aripo_sigmoid_coef <- filter(aripo_sigmoid_coef, 
                                    grepl('rna:', term)) 
      
    ## x) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'aswing' - outcome and tidy lm objects using 'broom'
      aripo_aswing_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(aswing ~ rna*group, data = .)))
      
    ## y) filter to remove intercepts and retain only intx coefs and test stats
      aripo_aswing_coef <- filter(aripo_aswing_coef, 
                                    grepl('rna:', term)) 
      
    ## z) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'swing' - outcome and tidy lm objects using 'broom'
      aripo_swing_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(swing ~ rna*group, data = .)))
      
    ## aa) filter to remove intercepts and retain only intx coefs and test stats
      aripo_swing_coef <- filter(aripo_swing_coef, 
                                    grepl('rna:', term)) 
      
    ## bb) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'center' - outcome and tidy lm objects using 'broom'
      aripo_center_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(center ~ rna*group, data = .)))
      
    ## cc) filter to remove intercepts and retain only intx coefs and test stats
      aripo_center_coef <- filter(aripo_center_coef, 
                                    grepl('rna:', term)) 
      
    ## dd) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'startle' - outcome and tidy lm objects using 'broom'
      aripo_startle_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(startle ~ rna*group, data = .)))
      
    ## ee) filter to remove intercepts and retain only intx coefs and test stats
      aripo_startle_coef <- filter(aripo_startle_coef, 
                                    grepl('rna:', term)) 
      
    ## ff) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'moving' - outcome and tidy lm objects using 'broom'
      aripo_moving_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(moving ~ rna*group, data = .)))
      
    ## gg) filter to remove intercepts and retain only intx coefs and test stats
      aripo_moving_coef <- filter(aripo_moving_coef, 
                                    grepl('rna:', term)) 
      
    ## hh) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'sigtime' - outcome and tidy lm objects using 'broom'
      aripo_sigtime_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(sigtime ~ rna*group, data = .)))
      
    ## ii) filter to remove intercepts and retain only intx coefs and test stats
      aripo_sigtime_coef <- filter(aripo_sigtime_coef, 
                                    grepl('rna:', term)) 
      
    ## jj) RNA expression site by site regression for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # 'asigtime' - outcome and tidy lm objects using 'broom'
      aripo_asigtime_coef<- aripo_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(asigtime ~ rna*group, data = .)))
      
    ## kk) filter to remove intercepts and retain only intx coefs and test stats
      aripo_asigtime_coef <- filter(aripo_asigtime_coef, 
                                    grepl('rna:', term)) 
      
      
  ### 4.2 Run Quare Effect Modification Models      
    ## a) Identify the first and last columns containing RNA expression IDs 
      # in quare_data
      colnames(quare_data)
      tail(colnames(quare_data))
      
    ## b) Transform quare_data data from wide to long based on RNA expression
      quare_rna_long <- quare_data %>%
        gather(key = 'site', value = 'rna', 
               'ENSPREG00000006268':'ENSPREG00000007674', 
               factor_key = F)
      
    ## c) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'contact' - outcome and tidy lm objects using 'broom'
      quare_contact_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(contact ~ rna*group, data = .)))
      
    ## d) Check if there are 8 estimates (intercept, rna, group & intx) per model
      summry_quare_contact_coef <- quare_contact_coef  %>%
        group_by(site) %>%
        summarize(n.terms = sum(!is.na(term)))
      
    ## e) filter to remove intercepts and retain only intx coefs and test stats
      quare_contact_coef <- filter(quare_contact_coef, 
                                   grepl('rna:', term)) 
      
    ## f) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'afaceoff' - outcome and tidy lm objects using 'broom'
      quare_afaceoff_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(afaceoff ~ rna*group, data = .)))
      
    ## g) filter to remove intercepts and retain only intx coefs and test stats
      quare_afaceoff_coef <- filter(quare_afaceoff_coef, 
                                    grepl('rna:', term))
      
    ## h) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'acontact' - outcome and tidy lm objects using 'broom'
      quare_acontact_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(acontact ~ rna*group, data = .)))
      
    ## i) filter to remove intercepts and retain only intx coefs and test stats
      quare_acontact_coef <- filter(quare_acontact_coef, 
                                    grepl('rna:', term)) 
      
    ## j) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'posturing' - outcome and tidy lm objects using 'broom'
      quare_posturing_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(posturing ~ rna*group, data = .)))
      
    ## k) filter to remove intercepts and retain only intx coefs and test stats
      quare_posturing_coef <- filter(quare_posturing_coef, 
                                     grepl('rna:', term)) 
      
    ## l) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'thrust' - outcome and tidy lm objects using 'broom'
      quare_thrust_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(thrust ~ rna*group, data = .)))
      
    ## m) filter to remove intercepts and retain only intx coefs and test stats
      quare_thrust_coef <- filter(quare_thrust_coef, 
                                  grepl('rna:', term)) 
      
    ## n) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'alunge' - outcome and tidy lm objects using 'broom'
      quare_alunge_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(alunge ~ rna*group, data = .)))
      
    ## o) filter to remove intercepts and retain only intx coefs and test stats
      quare_alunge_coef <- filter(quare_alunge_coef, 
                                  grepl('rna:', term)) 
      
    ## p) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'athrust' - outcome and tidy lm objects using 'broom'
      quare_athrust_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(athrust ~ rna*group, data = .)))
      
    ## q) filter to remove intercepts and retain only intx coefs and test stats
      quare_athrust_coef <- filter(quare_athrust_coef, 
                                   grepl('rna:', term)) 
      
    ## r) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'asigmoid' - outcome and tidy lm objects using 'broom'
      quare_asigmoid_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(asigmoid ~ rna*group, data = .)))
      
    ## s) filter to remove intercepts and retain only intx coefs and test stats
      quare_asigmoid_coef <- filter(quare_asigmoid_coef, 
                                    grepl('rna:', term)) 
      
    ## t) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'surface' - outcome and tidy lm objects using 'broom'
      quare_surface_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(surface ~ rna*group, data = .)))
      
    ## u) filter to remove intercepts and retain only intx coefs and test stats
      quare_surface_coef <- filter(quare_surface_coef, 
                                   grepl('rna:', term)) 
      
    ## v) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'sigmoid' - outcome and tidy lm objects using 'broom'
      quare_sigmoid_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(sigmoid ~ rna*group, data = .)))
      
    ## w) filter to remove intercepts and retain only intx coefs and test stats
      quare_sigmoid_coef <- filter(quare_sigmoid_coef, 
                                   grepl('rna:', term)) 
      
    ## x) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'aswing' - outcome and tidy lm objects using 'broom'
      quare_aswing_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(aswing ~ rna*group, data = .)))
      
    ## y) filter to remove intercepts and retain only intx coefs and test stats
      quare_aswing_coef <- filter(quare_aswing_coef, 
                                  grepl('rna:', term)) 
      
    ## z) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'swing' - outcome and tidy lm objects using 'broom'
      quare_swing_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(swing ~ rna*group, data = .)))
      
    ## aa) filter to remove intercepts and retain only intx coefs and test stats
      quare_swing_coef <- filter(quare_swing_coef, 
                                 grepl('rna:', term)) 
      
    ## bb) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'center' - outcome and tidy lm objects using 'broom'
      quare_center_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(center ~ rna*group, data = .)))
      
    ## cc) filter to remove intercepts and retain only intx coefs and test stats
      quare_center_coef <- filter(quare_center_coef, 
                                  grepl('rna:', term)) 
      
    ## dd) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'startle' - outcome and tidy lm objects using 'broom'
      quare_startle_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(startle ~ rna*group, data = .)))
      
    ## ee) filter to remove intercepts and retain only intx coefs and test stats
      quare_startle_coef <- filter(quare_startle_coef, 
                                   grepl('rna:', term)) 
      
    ## ff) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'moving' - outcome and tidy lm objects using 'broom'
      quare_moving_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(moving ~ rna*group, data = .)))
      
    ## gg) filter to remove intercepts and retain only intx coefs and test stats
      quare_moving_coef <- filter(quare_moving_coef, 
                                  grepl('rna:', term)) 
      
    ## hh) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'sigtime' - outcome and tidy lm objects using 'broom'
      quare_sigtime_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(sigtime ~ rna*group, data = .)))
      
    ## ii) filter to remove intercepts and retain only intx coefs and test stats
      quare_sigtime_coef <- filter(quare_sigtime_coef, 
                                   grepl('rna:', term)) 
      
    ## jj) RNA expression site by site regression for quare_rna_long
      # Estimate the RNA expression site by group interaction
      # 'asigtime' - outcome and tidy lm objects using 'broom'
      quare_asigtime_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(lm(asigtime ~ rna*group, data = .)))
      
    ## kk) filter to remove intercepts and retain only intx coefs and test stats
      quare_asigtime_coef <- filter(quare_asigtime_coef, 
                                    grepl('rna:', term)) 
      
      
  ## 4.3 Visualize Aripo distribution of interaction terms   
    ## a) Calculate the min and max estimates for x-axis scale
      aripo.contact.x.lim <- round(max(abs(aripo_contact_coef$estimate), 
                               na.rm = T), 1)
    ## b) Histogram Aripo rna x group intx terms
      aripo_contact_intx_hist <- ggplot(data = aripo_contact_coef, 
                                   aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.contact.x.lim, 
                                    aripo.contact.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.contact.x.lim, aripo.contact.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = contact)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## c) Print aripo_contact_intx_hist plots
      print(aripo_contact_intx_hist)
      
    ## d) Save aripo_contact_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_contact_intx_hist.pdf", plot = aripo_contact_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
    # ## e) Density plot Aripo rna x group intx terms
    #   ggplot(data = aripo_contact_coef, aes(x = estimate)) + 
    #     geom_density() +
    #     geom_vline(aes(xintercept = mean(estimate)),
    #                color = "blue", linetype = "dashed", size=1)
    #   
      
      
    ## f) Calculate the min and max estimates for x-axis scale
      aripo.afaceoff.x.lim <- round(max(abs(aripo_afaceoff_coef$estimate), 
                                        na.rm = T), 1)
    ## g) Histogram Aripo rna x group intx terms
      aripo_afaceoff_intx_hist <- ggplot(data = aripo_afaceoff_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.afaceoff.x.lim, 
                                    aripo.afaceoff.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.afaceoff.x.lim, aripo.afaceoff.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = afaceoff)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## h) Print aripo_afaceoff_intx_hist plots
      print(aripo_afaceoff_intx_hist)
      
    ## i) Save aripo_afaceoff_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_afaceoff_intx_hist.pdf", plot = aripo_afaceoff_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
    ## j) Calculate the min and max estimates for x-axis scale
      aripo.acontact.x.lim <- round(max(abs(aripo_acontact_coef$estimate), 
                                        na.rm = T), 1)
    ## k) Histogram Aripo rna x group intx terms
      aripo_acontact_intx_hist <- ggplot(data = aripo_acontact_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.acontact.x.lim, 
                                    aripo.acontact.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.acontact.x.lim, aripo.acontact.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = acontact)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## l) Print aripo_acontact_intx_hist plots
      print(aripo_acontact_intx_hist)
      
    ## m) Save aripo_acontact_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_acontact_intx_hist.pdf", plot = aripo_acontact_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## n) Calculate the min and max estimates for x-axis scale
      aripo.posturing.x.lim <- round(max(abs(aripo_posturing_coef$estimate), 
                                        na.rm = T), 1)
    ## o) Histogram Aripo rna x group intx terms
      aripo_posturing_intx_hist <- ggplot(data = aripo_posturing_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.posturing.x.lim, 
                                    aripo.posturing.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.posturing.x.lim, aripo.posturing.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = posturing)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## p) Print aripo_posturing_intx_hist plots
      print(aripo_posturing_intx_hist)
      
    ## q) Save aripo_posturing_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_posturing_intx_hist.pdf", plot = aripo_posturing_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## r) Calculate the min and max estimates for x-axis scale
      aripo.thrust.x.lim <- round(max(abs(aripo_thrust_coef$estimate), 
                                        na.rm = T), 1)
    ## s) Histogram Aripo rna x group intx terms
      aripo_thrust_intx_hist <- ggplot(data = aripo_thrust_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.thrust.x.lim, 
                                    aripo.thrust.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.thrust.x.lim, aripo.thrust.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = thrust)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## t) Print aripo_thrust_intx_hist plots
      print(aripo_thrust_intx_hist)
      
    ## u) Save aripo_thrust_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_thrust_intx_hist.pdf", plot = aripo_thrust_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## v) Calculate the min and max estimates for x-axis scale
      aripo.alunge.x.lim <- round(max(abs(aripo_alunge_coef$estimate), 
                                        na.rm = T), 1)
    ## w) Histogram Aripo rna x group intx terms
      aripo_alunge_intx_hist <- ggplot(data = aripo_alunge_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.alunge.x.lim, 
                                    aripo.alunge.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.alunge.x.lim, aripo.alunge.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = alunge)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## x) Print aripo_alunge_intx_hist plots
      print(aripo_alunge_intx_hist)
      
    ## y) Save aripo_alunge_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_alunge_intx_hist.pdf", plot = aripo_alunge_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## z) Calculate the min and max estimates for x-axis scale
      aripo.athrust.x.lim <- round(max(abs(aripo_athrust_coef$estimate), 
                                        na.rm = T), 1)
    ## aa) Histogram Aripo rna x group intx terms
      aripo_athrust_intx_hist <- ggplot(data = aripo_athrust_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.athrust.x.lim, 
                                    aripo.athrust.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.athrust.x.lim, aripo.athrust.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = athrust)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## bb) Print aripo_athrust_intx_hist plots
      print(aripo_athrust_intx_hist)
      
    ## cc) Save aripo_athrust_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_athrust_intx_hist.pdf", plot = aripo_athrust_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## dd) Calculate the min and max estimates for x-axis scale
      aripo.asigmoid.x.lim <- round(max(abs(aripo_asigmoid_coef$estimate), 
                                        na.rm = T), 1)
    ## ee) Histogram Aripo rna x group intx terms
      aripo_asigmoid_intx_hist <- ggplot(data = aripo_asigmoid_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.asigmoid.x.lim, 
                                    aripo.asigmoid.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.asigmoid.x.lim, aripo.asigmoid.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = asigmoid)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## ff) Print aripo_asigmoid_intx_hist plots
      print(aripo_asigmoid_intx_hist)
      
    ## gg) Save aripo_asigmoid_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_asigmoid_intx_hist.pdf", plot = aripo_asigmoid_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## hh) Calculate the min and max estimates for x-axis scale
      aripo.surface.x.lim <- round(max(abs(aripo_surface_coef$estimate), 
                                        na.rm = T), 1)
    ## ii) Histogram Aripo rna x group intx terms
      aripo_surface_intx_hist <- ggplot(data = aripo_surface_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.surface.x.lim, 
                                    aripo.surface.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.surface.x.lim, aripo.surface.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = surface)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## jj) Print aripo_surface_intx_hist plots
      print(aripo_surface_intx_hist)
      
    ## kk) Save aripo_surface_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_surface_intx_hist.pdf", plot = aripo_surface_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## ll) Calculate the min and max estimates for x-axis scale
      aripo.sigmoid.x.lim <- round(max(abs(aripo_sigmoid_coef$estimate), 
                                        na.rm = T), 1)
    ## mm) Histogram Aripo rna x group intx terms
      aripo_sigmoid_intx_hist <- ggplot(data = aripo_sigmoid_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.sigmoid.x.lim, 
                                    aripo.sigmoid.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.sigmoid.x.lim, aripo.sigmoid.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = sigmoid)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## nn) Print aripo_sigmoid_intx_hist plots
      print(aripo_sigmoid_intx_hist)
      
    ## oo) Save aripo_sigmoid_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_sigmoid_intx_hist.pdf", plot = aripo_sigmoid_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## pp) Calculate the min and max estimates for x-axis scale
      aripo.aswing.x.lim <- round(max(abs(aripo_aswing_coef$estimate), 
                                        na.rm = T), 1)
    ## qq) Histogram Aripo rna x group intx terms
      aripo_aswing_intx_hist <- ggplot(data = aripo_aswing_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.aswing.x.lim, 
                                    aripo.aswing.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.aswing.x.lim, aripo.aswing.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = aswing)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## rr) Print aripo_aswing_intx_hist plots
      print(aripo_aswing_intx_hist)
      
    ## ss) Save aripo_aswing_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_aswing_intx_hist.pdf", plot = aripo_aswing_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## tt) Calculate the min and max estimates for x-axis scale
      aripo.swing.x.lim <- round(max(abs(aripo_swing_coef$estimate), 
                                        na.rm = T), 1)
    ## uu) Histogram Aripo rna x group intx terms
      aripo_swing_intx_hist <- ggplot(data = aripo_swing_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.swing.x.lim, 
                                    aripo.swing.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.swing.x.lim, aripo.swing.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = swing)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## vv) Print aripo_swing_intx_hist plots
      print(aripo_swing_intx_hist)
      
    ## ww) Save aripo_swing_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_swing_intx_hist.pdf", plot = aripo_swing_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## xx) Calculate the min and max estimates for x-axis scale
      aripo.center.x.lim <- round(max(abs(aripo_center_coef$estimate), 
                                        na.rm = T), 1)
    ## yy) Histogram Aripo rna x group intx terms
      aripo_center_intx_hist <- ggplot(data = aripo_center_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.center.x.lim, 
                                    aripo.center.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.center.x.lim, aripo.center.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = center)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## zz) Print aripo_center_intx_hist plots
      print(aripo_center_intx_hist)
      
    ## aaa) Save aripo_center_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_center_intx_hist.pdf", plot = aripo_center_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## bbb) Calculate the min and max estimates for x-axis scale
      aripo.startle.x.lim <- round(max(abs(aripo_startle_coef$estimate), 
                                        na.rm = T), 1)
    ## ccc) Histogram Aripo rna x group intx terms
      aripo_startle_intx_hist <- ggplot(data = aripo_startle_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.startle.x.lim, 
                                    aripo.startle.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.startle.x.lim, aripo.startle.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = startle)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## ddd) Print aripo_startle_intx_hist plots
      print(aripo_startle_intx_hist)
      
    ## eee) Save aripo_startle_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_startle_intx_hist.pdf", plot = aripo_startle_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## fff) Calculate the min and max estimates for x-axis scale
      aripo.moving.x.lim <- round(max(abs(aripo_moving_coef$estimate), 
                                        na.rm = T), 1)
    ## ggg) Histogram Aripo rna x group intx terms
      aripo_moving_intx_hist <- ggplot(data = aripo_moving_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.moving.x.lim, 
                                    aripo.moving.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.moving.x.lim, aripo.moving.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = moving)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## hhh) Print aripo_moving_intx_hist plots
      print(aripo_moving_intx_hist)
      
    ## iii) Save aripo_moving_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_moving_intx_hist.pdf", plot = aripo_moving_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## jjj) Calculate the min and max estimates for x-axis scale
      aripo.sigtime.x.lim <- round(max(abs(aripo_sigtime_coef$estimate), 
                                        na.rm = T), 1)
    ## kkk) Histogram Aripo rna x group intx terms
      aripo_sigtime_intx_hist <- ggplot(data = aripo_sigtime_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.sigtime.x.lim, 
                                    aripo.sigtime.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.sigtime.x.lim, aripo.sigtime.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = sigtime)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## lll) Print aripo_sigtime_intx_hist plots
      print(aripo_sigtime_intx_hist)
      
    ## mmm) Save aripo_sigtime_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_sigtime_intx_hist.pdf", plot = aripo_sigtime_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## nnn) Calculate the min and max estimates for x-axis scale
      aripo.asigtime.x.lim <- round(max(abs(aripo_asigtime_coef$estimate), 
                                        na.rm = T), 1)
    ## ooo) Histogram Aripo rna x group intx terms
      aripo_asigtime_intx_hist <- ggplot(data = aripo_asigtime_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.asigtime.x.lim, 
                                    aripo.asigtime.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-aripo.asigtime.x.lim, aripo.asigtime.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = asigtime)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## ppp) Print aripo_asigtime_intx_hist plots
      print(aripo_asigtime_intx_hist)
      
    ## qqq) Save aripo_asigtime_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_asigtime_intx_hist.pdf", plot = aripo_asigtime_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
      
      
  ## 4.4 Visualize Quare distribution of interaction terms
    ## a) Calculate the min and max estimates for x-axis scale
      quare.contact.x.lim <- round(max(abs(quare_contact_coef$estimate), 
                                       na.rm = T), 1)
    ## b) Histogram Quare rna x group intx terms
      quare_contact_intx_hist <- ggplot(data = quare_contact_coef, 
                                        aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.contact.x.lim, 
                                    quare.contact.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.contact.x.lim, quare.contact.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = contact)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## c) Print quare_contact_intx_hist plots
      print(quare_contact_intx_hist)
      
    ## d) Save quare_contact_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_contact_intx_hist.pdf", plot = quare_contact_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
      # ## e) Density plot Quare rna x group intx terms
      #   ggplot(data = quare_contact_coef, aes(x = estimate)) + 
      #     geom_density() +
      #     geom_vline(aes(xintercept = mean(estimate)),
      #                color = "blue", linetype = "dashed", size=1)
      #   
      
      
    ## f) Calculate the min and max estimates for x-axis scale
      quare.afaceoff.x.lim <- round(max(abs(quare_afaceoff_coef$estimate), 
                                        na.rm = T), 1)
    ## g) Histogram Quare rna x group intx terms
      quare_afaceoff_intx_hist <- ggplot(data = quare_afaceoff_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.afaceoff.x.lim, 
                                    quare.afaceoff.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.afaceoff.x.lim, quare.afaceoff.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = afaceoff)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## h) Print quare_afaceoff_intx_hist plots
      print(quare_afaceoff_intx_hist)
      
    ## i) Save quare_afaceoff_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_afaceoff_intx_hist.pdf", plot = quare_afaceoff_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
    ## j) Calculate the min and max estimates for x-axis scale
      quare.acontact.x.lim <- round(max(abs(quare_acontact_coef$estimate), 
                                        na.rm = T), 1)
    ## k) Histogram Quare rna x group intx terms
      quare_acontact_intx_hist <- ggplot(data = quare_acontact_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.acontact.x.lim, 
                                    quare.acontact.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.acontact.x.lim, quare.acontact.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = acontact)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## l) Print quare_acontact_intx_hist plots
      print(quare_acontact_intx_hist)
      
    ## m) Save quare_acontact_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_acontact_intx_hist.pdf", plot = quare_acontact_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## n) Calculate the min and max estimates for x-axis scale
      quare.posturing.x.lim <- round(max(abs(quare_posturing_coef$estimate), 
                                         na.rm = T), 1)
    ## o) Histogram Quare rna x group intx terms
      quare_posturing_intx_hist <- ggplot(data = quare_posturing_coef, 
                                          aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.posturing.x.lim, 
                                    quare.posturing.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.posturing.x.lim, quare.posturing.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = posturing)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## p) Print quare_posturing_intx_hist plots
      print(quare_posturing_intx_hist)
      
    ## q) Save quare_posturing_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_posturing_intx_hist.pdf", plot = quare_posturing_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## r) Calculate the min and max estimates for x-axis scale
      quare.thrust.x.lim <- round(max(abs(quare_thrust_coef$estimate), 
                                      na.rm = T), 1)
    ## s) Histogram Quare rna x group intx terms
      quare_thrust_intx_hist <- ggplot(data = quare_thrust_coef, 
                                       aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.thrust.x.lim, 
                                    quare.thrust.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.thrust.x.lim, quare.thrust.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = thrust)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## t) Print quare_thrust_intx_hist plots
      print(quare_thrust_intx_hist)
      
    ## u) Save quare_thrust_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_thrust_intx_hist.pdf", plot = quare_thrust_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## v) Calculate the min and max estimates for x-axis scale
      quare.alunge.x.lim <- round(max(abs(quare_alunge_coef$estimate), 
                                      na.rm = T), 1)
    ## w) Histogram Quare rna x group intx terms
      quare_alunge_intx_hist <- ggplot(data = quare_alunge_coef, 
                                       aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.alunge.x.lim, 
                                    quare.alunge.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.alunge.x.lim, quare.alunge.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = alunge)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## x) Print quare_alunge_intx_hist plots
      print(quare_alunge_intx_hist)
      
    ## y) Save quare_alunge_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_alunge_intx_hist.pdf", plot = quare_alunge_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## z) Calculate the min and max estimates for x-axis scale
      quare.athrust.x.lim <- round(max(abs(quare_athrust_coef$estimate), 
                                       na.rm = T), 1)
    ## aa) Histogram Quare rna x group intx terms
      quare_athrust_intx_hist <- ggplot(data = quare_athrust_coef, 
                                        aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.athrust.x.lim, 
                                    quare.athrust.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.athrust.x.lim, quare.athrust.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = athrust)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## bb) Print quare_athrust_intx_hist plots
      print(quare_athrust_intx_hist)
      
    ## cc) Save quare_athrust_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_athrust_intx_hist.pdf", plot = quare_athrust_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## dd) Calculate the min and max estimates for x-axis scale
      quare.asigmoid.x.lim <- round(max(abs(quare_asigmoid_coef$estimate), 
                                        na.rm = T), 1)
    ## ee) Histogram Quare rna x group intx terms
      quare_asigmoid_intx_hist <- ggplot(data = quare_asigmoid_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.asigmoid.x.lim, 
                                    quare.asigmoid.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.asigmoid.x.lim, quare.asigmoid.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = asigmoid)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## ff) Print quare_asigmoid_intx_hist plots
      print(quare_asigmoid_intx_hist)
      
    ## gg) Save quare_asigmoid_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_asigmoid_intx_hist.pdf", plot = quare_asigmoid_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## hh) Calculate the min and max estimates for x-axis scale
      quare.surface.x.lim <- round(max(abs(quare_surface_coef$estimate), 
                                       na.rm = T), 1)
    ## ii) Histogram Quare rna x group intx terms
      quare_surface_intx_hist <- ggplot(data = quare_surface_coef, 
                                        aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.surface.x.lim, 
                                    quare.surface.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.surface.x.lim, quare.surface.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = surface)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## jj) Print quare_surface_intx_hist plots
      print(quare_surface_intx_hist)
      
    ## kk) Save quare_surface_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_surface_intx_hist.pdf", plot = quare_surface_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## ll) Calculate the min and max estimates for x-axis scale
      quare.sigmoid.x.lim <- round(max(abs(quare_sigmoid_coef$estimate), 
                                       na.rm = T), 1)
    ## mm) Histogram Quare rna x group intx terms
      quare_sigmoid_intx_hist <- ggplot(data = quare_sigmoid_coef, 
                                        aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.sigmoid.x.lim, 
                                    quare.sigmoid.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.sigmoid.x.lim, quare.sigmoid.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = sigmoid)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## nn) Print quare_sigmoid_intx_hist plots
      print(quare_sigmoid_intx_hist)
      
    ## oo) Save quare_sigmoid_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_sigmoid_intx_hist.pdf", plot = quare_sigmoid_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## pp) Calculate the min and max estimates for x-axis scale
      quare.aswing.x.lim <- round(max(abs(quare_aswing_coef$estimate), 
                                      na.rm = T), 1)
    ## qq) Histogram Quare rna x group intx terms
      quare_aswing_intx_hist <- ggplot(data = quare_aswing_coef, 
                                       aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.aswing.x.lim, 
                                    quare.aswing.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.aswing.x.lim, quare.aswing.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = aswing)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## rr) Print quare_aswing_intx_hist plots
      print(quare_aswing_intx_hist)
      
    ## ss) Save quare_aswing_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_aswing_intx_hist.pdf", plot = quare_aswing_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## tt) Calculate the min and max estimates for x-axis scale
      quare.swing.x.lim <- round(max(abs(quare_swing_coef$estimate), 
                                     na.rm = T), 1)
    ## uu) Histogram Quare rna x group intx terms
      quare_swing_intx_hist <- ggplot(data = quare_swing_coef, 
                                      aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.swing.x.lim, 
                                    quare.swing.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.swing.x.lim, quare.swing.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = swing)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## vv) Print quare_swing_intx_hist plots
      print(quare_swing_intx_hist)
      
    ## ww) Save quare_swing_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_swing_intx_hist.pdf", plot = quare_swing_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## xx) Calculate the min and max estimates for x-axis scale
      quare.center.x.lim <- round(max(abs(quare_center_coef$estimate), 
                                      na.rm = T), 1)
    ## yy) Histogram Quare rna x group intx terms
      quare_center_intx_hist <- ggplot(data = quare_center_coef, 
                                       aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.center.x.lim, 
                                    quare.center.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.center.x.lim, quare.center.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = center)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## zz) Print quare_center_intx_hist plots
      print(quare_center_intx_hist)
      
    ## aaa) Save quare_center_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_center_intx_hist.pdf", plot = quare_center_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## bbb) Calculate the min and max estimates for x-axis scale
      quare.startle.x.lim <- round(max(abs(quare_startle_coef$estimate), 
                                       na.rm = T), 1)
    ## ccc) Histogram Quare rna x group intx terms
      quare_startle_intx_hist <- ggplot(data = quare_startle_coef, 
                                        aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.startle.x.lim, 
                                    quare.startle.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.startle.x.lim, quare.startle.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = startle)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## ddd) Print quare_startle_intx_hist plots
      print(quare_startle_intx_hist)
      
    ## eee) Save quare_startle_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_startle_intx_hist.pdf", plot = quare_startle_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## fff) Calculate the min and max estimates for x-axis scale
      quare.moving.x.lim <- round(max(abs(quare_moving_coef$estimate), 
                                      na.rm = T), 1)
    ## ggg) Histogram Quare rna x group intx terms
      quare_moving_intx_hist <- ggplot(data = quare_moving_coef, 
                                       aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.moving.x.lim, 
                                    quare.moving.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.moving.x.lim, quare.moving.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = moving)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## hhh) Print quare_moving_intx_hist plots
      print(quare_moving_intx_hist)
      
    ## iii) Save quare_moving_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_moving_intx_hist.pdf", plot = quare_moving_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## jjj) Calculate the min and max estimates for x-axis scale
      quare.sigtime.x.lim <- round(max(abs(quare_sigtime_coef$estimate), 
                                       na.rm = T), 1)
    ## kkk) Histogram Quare rna x group intx terms
      quare_sigtime_intx_hist <- ggplot(data = quare_sigtime_coef, 
                                        aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.sigtime.x.lim, 
                                    quare.sigtime.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.sigtime.x.lim, quare.sigtime.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = sigtime)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## lll) Print quare_sigtime_intx_hist plots
      print(quare_sigtime_intx_hist)
      
    ## mmm) Save quare_sigtime_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_sigtime_intx_hist.pdf", plot = quare_sigtime_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## nnn) Calculate the min and max estimates for x-axis scale
      quare.asigtime.x.lim <- round(max(abs(quare_asigtime_coef$estimate), 
                                        na.rm = T), 1)
    ## ooo) Histogram Quare rna x group intx terms
      quare_asigtime_intx_hist <- ggplot(data = quare_asigtime_coef, 
                                         aes(x = estimate)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.asigtime.x.lim, 
                                    quare.asigtime.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.asigtime.x.lim, quare.asigtime.x.lim)) +
        geom_vline(aes(xintercept = mean(estimate, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Quare Histogram of RNA expression by Group Interaction 
Terms (Outcome = asigtime)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## ppp) Print quare_asigtime_intx_hist plots
      print(quare_asigtime_intx_hist)
      
    ## qqq) Save quare_asigtime_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("quare_asigtime_intx_hist.pdf", plot = quare_asigtime_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
 
  ### 4.5  Aripo Mutliple comparisons corrections
    ## a) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_contact_coef$fdr.p = p.adjust(aripo_contact_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_contact_coef$p.value))
      
    ## b) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.contact.0.3 <- 
        max(aripo_contact_coef$p.value[aripo_contact_coef$fdr.p <= 0.3],
            na.rm = T)
    
    ## c) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.contact <- 4 * round(sd(aripo_contact_coef$estimate, 
                                        na.rm = T), 3)
      
    ## d) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_afaceoff_coef$fdr.p = p.adjust(aripo_afaceoff_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_afaceoff_coef$p.value))
      
    ## e) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.afaceoff.0.3 <- 
        max(aripo_afaceoff_coef$p.value[aripo_afaceoff_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## f) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.afaceoff <- 4 * round(sd(aripo_afaceoff_coef$estimate, 
                                        na.rm = T), 3)
      
    ## g) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_acontact_coef$fdr.p = p.adjust(aripo_acontact_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_acontact_coef$p.value))
      
    ## h) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.acontact.0.3 <- 
        max(aripo_acontact_coef$p.value[aripo_acontact_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## i) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.acontact <- 4 * round(sd(aripo_acontact_coef$estimate, 
                                        na.rm = T), 3)
      
    ## j) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_posturing_coef$fdr.p = p.adjust(aripo_posturing_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_posturing_coef$p.value))
      
    ## k) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.posturing.0.3 <- 
        max(aripo_posturing_coef$p.value[aripo_posturing_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## l) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.posturing <- 4 * round(sd(aripo_posturing_coef$estimate, 
                                        na.rm = T), 3)
      
    ## m) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_thrust_coef$fdr.p = p.adjust(aripo_thrust_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_thrust_coef$p.value))
      
    ## n) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.thrust.0.3 <- 
        max(aripo_thrust_coef$p.value[aripo_thrust_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## o) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.thrust <- 4 * round(sd(aripo_thrust_coef$estimate, 
                                        na.rm = T), 3)
      
    ## p) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_alunge_coef$fdr.p = p.adjust(aripo_alunge_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_alunge_coef$p.value))
      
    ## q) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.alunge.0.3 <- 
        max(aripo_alunge_coef$p.value[aripo_alunge_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## r) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.alunge <- 4 * round(sd(aripo_alunge_coef$estimate, 
                                        na.rm = T), 3)
      
    ## s) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_athrust_coef$fdr.p = p.adjust(aripo_athrust_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_athrust_coef$p.value))
      
    ## t) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.athrust.0.3 <- 
        max(aripo_athrust_coef$p.value[aripo_athrust_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## u) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.athrust <- 4 * round(sd(aripo_athrust_coef$estimate, 
                                        na.rm = T), 3)
      
    ## v) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_asigmoid_coef$fdr.p = p.adjust(aripo_asigmoid_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_asigmoid_coef$p.value))
      
    ## x) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.asigmoid.0.3 <- 
        max(aripo_asigmoid_coef$p.value[aripo_asigmoid_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## y) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.asigmoid <- 4 * round(sd(aripo_asigmoid_coef$estimate, 
                                        na.rm = T), 3)
      
    ## z) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_surface_coef$fdr.p = p.adjust(aripo_surface_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_surface_coef$p.value))
      
    ## aa) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.surface.0.3 <- 
        max(aripo_surface_coef$p.value[aripo_surface_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## bb) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.surface <- 4 * round(sd(aripo_surface_coef$estimate, 
                                        na.rm = T), 3)
      
    ## cc) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_sigmoid_coef$fdr.p = p.adjust(aripo_sigmoid_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_sigmoid_coef$p.value))
      
    ## dd) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.sigmoid.0.3 <- 
        max(aripo_sigmoid_coef$p.value[aripo_sigmoid_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## ee) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.sigmoid <- 4 * round(sd(aripo_sigmoid_coef$estimate, 
                                        na.rm = T), 3)
      
    ## ff) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_aswing_coef$fdr.p = p.adjust(aripo_aswing_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_aswing_coef$p.value))
      
    ## gg) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.aswing.0.3 <- 
        max(aripo_aswing_coef$p.value[aripo_aswing_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## hh) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.aswing <- 4 * round(sd(aripo_aswing_coef$estimate, 
                                        na.rm = T), 3)
      
    ## ii) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_swing_coef$fdr.p = p.adjust(aripo_swing_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_swing_coef$p.value))
      
    ## jj) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.swing.0.3 <- 
        max(aripo_swing_coef$p.value[aripo_swing_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## kk) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.swing <- 4 * round(sd(aripo_swing_coef$estimate, 
                                        na.rm = T), 3)
      
    ## ll) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_center_coef$fdr.p = p.adjust(aripo_center_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_center_coef$p.value))
      
    ## mm) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.center.0.3 <- 
        max(aripo_center_coef$p.value[aripo_center_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## nn) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.center <- 4 * round(sd(aripo_center_coef$estimate, 
                                        na.rm = T), 3)
      
    ## oo) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_startle_coef$fdr.p = p.adjust(aripo_startle_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_startle_coef$p.value))
      
    ## pp) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.startle.0.3 <- 
        max(aripo_startle_coef$p.value[aripo_startle_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## qq) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.startle <- 4 * round(sd(aripo_startle_coef$estimate, 
                                        na.rm = T), 3)
      
    ## rr) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_moving_coef$fdr.p = p.adjust(aripo_moving_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_moving_coef$p.value))
      
    ## ss) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.moving.0.3 <- 
        max(aripo_moving_coef$p.value[aripo_moving_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## tt) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.moving <- 4 * round(sd(aripo_moving_coef$estimate, 
                                        na.rm = T), 3)
      
    ## uu) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_sigtime_coef$fdr.p = p.adjust(aripo_sigtime_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_sigtime_coef$p.value))
      
    ## vv) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.sigtime.0.3 <- 
        max(aripo_sigtime_coef$p.value[aripo_sigtime_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## ww) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.sigtime <- 4 * round(sd(aripo_sigtime_coef$estimate, 
                                        na.rm = T), 3)
      
    ## xx) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_asigtime_coef$fdr.p = p.adjust(aripo_asigtime_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_asigtime_coef$p.value))
      
    ## yy) extract the FDR p-value.threshold based on p.adjust <= 0.3
      aripo.p.asigtime.0.3 <- 
        max(aripo_asigtime_coef$p.value[aripo_asigtime_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## zz) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      aripo.fc.asigtime <- 4 * round(sd(aripo_asigtime_coef$estimate, 
                                        na.rm = T), 3)

      
  ### 4.6 Aripo Volcano Plots
    ## a) Volcano Plot 
      aripo.asigmoid_vol <- EnhancedVolcano(aripo_asigmoid_coef,
                                            lab = rownames(aripo_asigmoid_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.asigmoid.x.lim, aripo.asigmoid.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = asigmoid',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.asigmoid,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.asigmoid.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## b) Print aripo_asigmoid_vol plots    
      print(aripo.asigmoid_vol)
      
    ## c) Save aripo.asigmoid_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_asigmoid_vol.pdf", plot = aripo.asigmoid_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
    ## d) Volcano Plot 
      aripo.afaceoff_vol <- EnhancedVolcano(aripo_afaceoff_coef,
                                            lab = rownames(aripo_afaceoff_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.afaceoff.x.lim, aripo.afaceoff.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = afaceoff',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.afaceoff,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.afaceoff.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## e) Print aripo_afaceoff_vol plots    
      print(aripo.afaceoff_vol)
      
    ## f) Save aripo.afaceoff_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_afaceoff_vol.pdf", plot = aripo.afaceoff_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## g) Volcano Plot 
      aripo.posturing_vol <- EnhancedVolcano(aripo_posturing_coef,
                                            lab = rownames(aripo_posturing_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.posturing.x.lim, aripo.posturing.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = posturing',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.posturing,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.posturing.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## h) Print aripo_posturing_vol plots    
      print(aripo.posturing_vol)
      
    ## i) Save aripo.posturing_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_posturing_vol.pdf", plot = aripo.posturing_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## j) Volcano Plot 
      aripo.thrust_vol <- EnhancedVolcano(aripo_thrust_coef,
                                            lab = rownames(aripo_thrust_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.thrust.x.lim, aripo.thrust.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = thrust',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.thrust,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.thrust.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## k) Print aripo_thrust_vol plots    
      print(aripo.thrust_vol)
      
    ## l) Save aripo.thrust_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_thrust_vol.pdf", plot = aripo.thrust_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## m) Volcano Plot 
      aripo.contact_vol <- EnhancedVolcano(aripo_contact_coef,
                                            lab = rownames(aripo_contact_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.contact.x.lim, aripo.contact.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = contact',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.contact,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.contact.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## n) Print aripo_contact_vol plots    
      print(aripo.contact_vol)
      
    ## o) Save aripo.contact_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_contact_vol.pdf", plot = aripo.contact_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## p) Volcano Plot 
      aripo.sigmoid_vol <- EnhancedVolcano(aripo_sigmoid_coef,
                                            lab = rownames(aripo_sigmoid_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.sigmoid.x.lim, aripo.sigmoid.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = sigmoid',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.sigmoid,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.sigmoid.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## q) Print aripo_sigmoid_vol plots    
      print(aripo.sigmoid_vol)
      
    ## r) Save aripo.sigmoid_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_sigmoid_vol.pdf", plot = aripo.sigmoid_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## s) Volcano Plot 
      aripo.startle_vol <- EnhancedVolcano(aripo_startle_coef,
                                            lab = rownames(aripo_startle_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.startle.x.lim, aripo.startle.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = startle',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.startle,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.startle.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## t) Print aripo_startle_vol plots    
      print(aripo.startle_vol)
      
    ## u) Save aripo.startle_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_startle_vol.pdf", plot = aripo.startle_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## v) Volcano Plot 
      aripo.alunge_vol <- EnhancedVolcano(aripo_alunge_coef,
                                            lab = rownames(aripo_alunge_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.alunge.x.lim, aripo.alunge.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = alunge',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.alunge,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.alunge.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## w) Print aripo_alunge_vol plots    
      print(aripo.alunge_vol)
      
    ## x) Save aripo.alunge_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_alunge_vol.pdf", plot = aripo.alunge_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## y) Volcano Plot 
      aripo.acontact_vol <- EnhancedVolcano(aripo_acontact_coef,
                                            lab = rownames(aripo_acontact_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.acontact.x.lim, aripo.acontact.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = acontact',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.acontact,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.acontact.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## z) Print aripo_acontact_vol plots    
      print(aripo.acontact_vol)
      
    ## aa) Save aripo.acontact_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_acontact_vol.pdf", plot = aripo.acontact_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## bb) Volcano Plot 
      aripo.swing_vol <- EnhancedVolcano(aripo_swing_coef,
                                            lab = rownames(aripo_swing_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.swing.x.lim, aripo.swing.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = swing',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.swing,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.swing.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## cc) Print aripo_swing_vol plots    
      print(aripo.swing_vol)
      
    ## dd) Save aripo.swing_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_swing_vol.pdf", plot = aripo.swing_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## ee) Volcano Plot 
      aripo.athrust_vol <- EnhancedVolcano(aripo_athrust_coef,
                                            lab = rownames(aripo_athrust_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.athrust.x.lim, aripo.athrust.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = athrust',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.athrust,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.athrust.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## ff) Print aripo_athrust_vol plots    
      print(aripo.athrust_vol)
      
    ## gg) Save aripo.athrust_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_athrust_vol.pdf", plot = aripo.athrust_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## hh) Volcano Plot 
      aripo.asigtime_vol <- EnhancedVolcano(aripo_asigtime_coef,
                                            lab = rownames(aripo_asigtime_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.asigtime.x.lim, aripo.asigtime.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = asigtime',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.asigtime,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.asigtime.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## ii) Print aripo_asigtime_vol plots    
      print(aripo.asigtime_vol)
      
    ## jj) Save aripo.asigtime_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_asigtime_vol.pdf", plot = aripo.asigtime_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## kk) Volcano Plot 
      aripo.aswing_vol <- EnhancedVolcano(aripo_aswing_coef,
                                            lab = rownames(aripo_aswing_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.aswing.x.lim, aripo.aswing.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = aswing',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.aswing,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.aswing.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## ll) Print aripo_aswing_vol plots    
      print(aripo.aswing_vol)
      
    ## mm) Save aripo.aswing_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_aswing_vol.pdf", plot = aripo.aswing_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## nn) Volcano Plot 
      aripo.center_vol <- EnhancedVolcano(aripo_center_coef,
                                            lab = rownames(aripo_center_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.center.x.lim, aripo.center.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = center',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.center,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.center.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## oo) Print aripo_center_vol plots    
      print(aripo.center_vol)
      
    ## pp) Save aripo.center_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_center_vol.pdf", plot = aripo.center_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## qq) Volcano Plot 
      aripo.moving_vol <- EnhancedVolcano(aripo_moving_coef,
                                            lab = rownames(aripo_moving_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.moving.x.lim, aripo.moving.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = moving',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.moving,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.moving.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## rr) Print aripo_moving_vol plots    
      print(aripo.moving_vol)
      
    ## ss) Save aripo.moving_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_moving_vol.pdf", plot = aripo.moving_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## tt) Volcano Plot 
      aripo.sigtime_vol <- EnhancedVolcano(aripo_sigtime_coef,
                                            lab = rownames(aripo_sigtime_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.sigtime.x.lim, aripo.sigtime.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = sigtime',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.sigtime,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.sigtime.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## uu) Print aripo_sigtime_vol plots    
      print(aripo.sigtime_vol)
      
    ## vv) Save aripo.sigtime_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_sigtime_vol.pdf", plot = aripo.sigtime_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## ww) Volcano Plot 
      aripo.surface_vol <- EnhancedVolcano(aripo_surface_coef,
                                            lab = rownames(aripo_surface_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-aripo.surface.x.lim, aripo.surface.x.lim),
                                            title = 'Aripo Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = surface',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = aripo.fc.surface,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = aripo.p.surface.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## xx) Print aripo_surface_vol plots    
      print(aripo.surface_vol)
      
    ## yy) Save aripo.surface_vol plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_surface_vol.pdf", plot = aripo.surface_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    

  ### 4.7 Quare Mutliple comparisons corrections
    ## a) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_contact_coef$fdr.p = p.adjust(quare_contact_coef$p.value, 
                                          method = "BH", 
                                          length(quare_contact_coef$p.value))
      
    ## b) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.contact.0.3 <- 
        max(quare_contact_coef$p.value[quare_contact_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## c) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.contact <- 4 * round(sd(quare_contact_coef$estimate, 
                                       na.rm = T), 3)
      
    ## d) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_afaceoff_coef$fdr.p = p.adjust(quare_afaceoff_coef$p.value, 
                                           method = "BH", 
                                           length(quare_afaceoff_coef$p.value))
      
    ## e) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.afaceoff.0.3 <- 
        max(quare_afaceoff_coef$p.value[quare_afaceoff_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## f) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.afaceoff <- 4 * round(sd(quare_afaceoff_coef$estimate, 
                                        na.rm = T), 3)
      
    ## g) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_acontact_coef$fdr.p = p.adjust(quare_acontact_coef$p.value, 
                                           method = "BH", 
                                           length(quare_acontact_coef$p.value))
      
    ## h) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.acontact.0.3 <- 
        max(quare_acontact_coef$p.value[quare_acontact_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## i) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.acontact <- 4 * round(sd(quare_acontact_coef$estimate, 
                                        na.rm = T), 3)
      
    ## j) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_posturing_coef$fdr.p = p.adjust(quare_posturing_coef$p.value, 
                                            method = "BH", 
                                            length(quare_posturing_coef$p.value))
      
    ## k) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.posturing.0.3 <- 
        max(quare_posturing_coef$p.value[quare_posturing_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## l) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.posturing <- 4 * round(sd(quare_posturing_coef$estimate, 
                                         na.rm = T), 3)
      
    ## m) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_thrust_coef$fdr.p = p.adjust(quare_thrust_coef$p.value, 
                                         method = "BH", 
                                         length(quare_thrust_coef$p.value))
      
    ## n) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.thrust.0.3 <- 
        max(quare_thrust_coef$p.value[quare_thrust_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## o) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.thrust <- 4 * round(sd(quare_thrust_coef$estimate, 
                                      na.rm = T), 3)
      
    ## p) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_alunge_coef$fdr.p = p.adjust(quare_alunge_coef$p.value, 
                                         method = "BH", 
                                         length(quare_alunge_coef$p.value))
      
    ## q) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.alunge.0.3 <- 
        max(quare_alunge_coef$p.value[quare_alunge_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## r) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.alunge <- 4 * round(sd(quare_alunge_coef$estimate, 
                                      na.rm = T), 3)
      
    ## s) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_athrust_coef$fdr.p = p.adjust(quare_athrust_coef$p.value, 
                                          method = "BH", 
                                          length(quare_athrust_coef$p.value))
      
    ## t) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.athrust.0.3 <- 
        max(quare_athrust_coef$p.value[quare_athrust_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## u) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.athrust <- 4 * round(sd(quare_athrust_coef$estimate, 
                                       na.rm = T), 3)
      
    ## v) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_asigmoid_coef$fdr.p = p.adjust(quare_asigmoid_coef$p.value, 
                                           method = "BH", 
                                           length(quare_asigmoid_coef$p.value))
      
    ## x) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.asigmoid.0.3 <- 
        max(quare_asigmoid_coef$p.value[quare_asigmoid_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## y) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.asigmoid <- 4 * round(sd(quare_asigmoid_coef$estimate, 
                                        na.rm = T), 3)
      
    ## z) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_surface_coef$fdr.p = p.adjust(quare_surface_coef$p.value, 
                                          method = "BH", 
                                          length(quare_surface_coef$p.value))
      
    ## aa) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.surface.0.3 <- 
        max(quare_surface_coef$p.value[quare_surface_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## bb) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.surface <- 4 * round(sd(quare_surface_coef$estimate, 
                                       na.rm = T), 3)
      
    ## cc) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_sigmoid_coef$fdr.p = p.adjust(quare_sigmoid_coef$p.value, 
                                          method = "BH", 
                                          length(quare_sigmoid_coef$p.value))
      
    ## dd) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.sigmoid.0.3 <- 
        max(quare_sigmoid_coef$p.value[quare_sigmoid_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## ee) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.sigmoid <- 4 * round(sd(quare_sigmoid_coef$estimate, 
                                       na.rm = T), 3)
      
    ## ff) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_aswing_coef$fdr.p = p.adjust(quare_aswing_coef$p.value, 
                                         method = "BH", 
                                         length(quare_aswing_coef$p.value))
      
    ## gg) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.aswing.0.3 <- 
        max(quare_aswing_coef$p.value[quare_aswing_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## hh) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.aswing <- 4 * round(sd(quare_aswing_coef$estimate, 
                                      na.rm = T), 3)
      
    ## ii) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_swing_coef$fdr.p = p.adjust(quare_swing_coef$p.value, 
                                        method = "BH", 
                                        length(quare_swing_coef$p.value))
      
    ## jj) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.swing.0.3 <- 
        max(quare_swing_coef$p.value[quare_swing_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## kk) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.swing <- 4 * round(sd(quare_swing_coef$estimate, 
                                     na.rm = T), 3)
      
    ## ll) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_center_coef$fdr.p = p.adjust(quare_center_coef$p.value, 
                                         method = "BH", 
                                         length(quare_center_coef$p.value))
      
    ## mm) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.center.0.3 <- 
        max(quare_center_coef$p.value[quare_center_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## nn) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.center <- 4 * round(sd(quare_center_coef$estimate, 
                                      na.rm = T), 3)
      
    ## oo) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_startle_coef$fdr.p = p.adjust(quare_startle_coef$p.value, 
                                          method = "BH", 
                                          length(quare_startle_coef$p.value))
      
    ## pp) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.startle.0.3 <- 
        max(quare_startle_coef$p.value[quare_startle_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## qq) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.startle <- 4 * round(sd(quare_startle_coef$estimate, 
                                       na.rm = T), 3)
      
    ## rr) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_moving_coef$fdr.p = p.adjust(quare_moving_coef$p.value, 
                                         method = "BH", 
                                         length(quare_moving_coef$p.value))
      
    ## ss) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.moving.0.3 <- 
        max(quare_moving_coef$p.value[quare_moving_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## tt) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.moving <- 4 * round(sd(quare_moving_coef$estimate, 
                                      na.rm = T), 3)
      
    ## uu) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_sigtime_coef$fdr.p = p.adjust(quare_sigtime_coef$p.value, 
                                          method = "BH", 
                                          length(quare_sigtime_coef$p.value))
      
    ## vv) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.sigtime.0.3 <- 
        max(quare_sigtime_coef$p.value[quare_sigtime_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## ww) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.sigtime <- 4 * round(sd(quare_sigtime_coef$estimate, 
                                       na.rm = T), 3)
      
    ## xx) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_asigtime_coef$fdr.p = p.adjust(quare_asigtime_coef$p.value, 
                                           method = "BH", 
                                           length(quare_asigtime_coef$p.value))
      
    ## yy) extract the FDR p-value.threshold based on p.adjust <= 0.3
      quare.p.asigtime.0.3 <- 
        max(quare_asigtime_coef$p.value[quare_asigtime_coef$fdr.p <= 0.3],
            na.rm = T)
      
    ## zz) Calculate an FC cutoff based on 2 SD of the intx. term estimates
      quare.fc.asigtime <- 4 * round(sd(quare_asigtime_coef$estimate, 
                                        na.rm = T), 3)
      
      
      
  ### 4.8 Quare Volcano Plots
    ## a) Volcano Plot 
      quare.asigmoid_vol <- EnhancedVolcano(quare_asigmoid_coef,
                                            lab = rownames(quare_asigmoid_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-quare.asigmoid.x.lim, quare.asigmoid.x.lim),
                                            title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = asigmoid',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = quare.fc.asigmoid,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = quare.p.asigmoid.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## b) Print quare_asigmoid_vol plots    
      print(quare.asigmoid_vol)
      
    ## c) Save quare.asigmoid_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_asigmoid_vol.pdf", plot = quare.asigmoid_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
    ## d) Volcano Plot 
      quare.afaceoff_vol <- EnhancedVolcano(quare_afaceoff_coef,
                                            lab = rownames(quare_afaceoff_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-quare.afaceoff.x.lim, quare.afaceoff.x.lim),
                                            title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = afaceoff',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = quare.fc.afaceoff,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = quare.p.afaceoff.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## e) Print quare_afaceoff_vol plots    
      print(quare.afaceoff_vol)
      
    ## f) Save quare.afaceoff_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_afaceoff_vol.pdf", plot = quare.afaceoff_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## g) Volcano Plot 
      quare.posturing_vol <- EnhancedVolcano(quare_posturing_coef,
                                             lab = rownames(quare_posturing_coef),
                                             x = 'estimate',
                                             y = 'p.value',
                                             xlim = c(-quare.posturing.x.lim, quare.posturing.x.lim),
                                             title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                             subtitle = 'Behavior = posturing',
                                             caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                             xlab = 'Intx. β coefficients (group x RNA)',
                                             pCutoff = 0.05,
                                             FCcutoff = quare.fc.posturing,
                                             # transcriptPointSize = 1.5,
                                             # transcriptLabSize = 3.0,
                                             colAlpha = 1,
                                             cutoffLineType = 'longdash',
                                             cutoffLineCol = 'grey50',
                                             cutoffLineWidth = 0.8,
                                             hline = quare.p.posturing.0.3,
                                             hlineCol = 'red',
                                             hlineType = 'dotdash',
                                             hlineWidth = 0.8,
                                             gridlines.major = F,
                                             gridlines.minor = F)
      
    ## h) Print quare_posturing_vol plots    
      print(quare.posturing_vol)
      
    ## i) Save quare.posturing_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_posturing_vol.pdf", plot = quare.posturing_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## j) Volcano Plot 
      quare.thrust_vol <- EnhancedVolcano(quare_thrust_coef,
                                          lab = rownames(quare_thrust_coef),
                                          x = 'estimate',
                                          y = 'p.value',
                                          xlim = c(-quare.thrust.x.lim, quare.thrust.x.lim),
                                          title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                          subtitle = 'Behavior = thrust',
                                          caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                          xlab = 'Intx. β coefficients (group x RNA)',
                                          pCutoff = 0.05,
                                          FCcutoff = quare.fc.thrust,
                                          # transcriptPointSize = 1.5,
                                          # transcriptLabSize = 3.0,
                                          colAlpha = 1,
                                          cutoffLineType = 'longdash',
                                          cutoffLineCol = 'grey50',
                                          cutoffLineWidth = 0.8,
                                          hline = quare.p.thrust.0.3,
                                          hlineCol = 'red',
                                          hlineType = 'dotdash',
                                          hlineWidth = 0.8,
                                          gridlines.major = F,
                                          gridlines.minor = F)
      
    ## k) Print quare_thrust_vol plots    
      print(quare.thrust_vol)
      
    ## l) Save quare.thrust_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_thrust_vol.pdf", plot = quare.thrust_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## m) Volcano Plot 
      quare.contact_vol <- EnhancedVolcano(quare_contact_coef,
                                           lab = rownames(quare_contact_coef),
                                           x = 'estimate',
                                           y = 'p.value',
                                           xlim = c(-quare.contact.x.lim, quare.contact.x.lim),
                                           title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                           subtitle = 'Behavior = contact',
                                           caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                           xlab = 'Intx. β coefficients (group x RNA)',
                                           pCutoff = 0.05,
                                           FCcutoff = quare.fc.contact,
                                           # transcriptPointSize = 1.5,
                                           # transcriptLabSize = 3.0,
                                           colAlpha = 1,
                                           cutoffLineType = 'longdash',
                                           cutoffLineCol = 'grey50',
                                           cutoffLineWidth = 0.8,
                                           hline = quare.p.contact.0.3,
                                           hlineCol = 'red',
                                           hlineType = 'dotdash',
                                           hlineWidth = 0.8,
                                           gridlines.major = F,
                                           gridlines.minor = F)
      
    ## n) Print quare_contact_vol plots    
      print(quare.contact_vol)
      
    ## o) Save quare.contact_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_contact_vol.pdf", plot = quare.contact_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## p) Volcano Plot 
      quare.sigmoid_vol <- EnhancedVolcano(quare_sigmoid_coef,
                                           lab = rownames(quare_sigmoid_coef),
                                           x = 'estimate',
                                           y = 'p.value',
                                           xlim = c(-quare.sigmoid.x.lim, quare.sigmoid.x.lim),
                                           title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                           subtitle = 'Behavior = sigmoid',
                                           caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                           xlab = 'Intx. β coefficients (group x RNA)',
                                           pCutoff = 0.05,
                                           FCcutoff = quare.fc.sigmoid,
                                           # transcriptPointSize = 1.5,
                                           # transcriptLabSize = 3.0,
                                           colAlpha = 1,
                                           cutoffLineType = 'longdash',
                                           cutoffLineCol = 'grey50',
                                           cutoffLineWidth = 0.8,
                                           hline = quare.p.sigmoid.0.3,
                                           hlineCol = 'red',
                                           hlineType = 'dotdash',
                                           hlineWidth = 0.8,
                                           gridlines.major = F,
                                           gridlines.minor = F)
      
    ## q) Print quare_sigmoid_vol plots    
      print(quare.sigmoid_vol)
      
    ## r) Save quare.sigmoid_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_sigmoid_vol.pdf", plot = quare.sigmoid_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## s) Volcano Plot 
      quare.startle_vol <- EnhancedVolcano(quare_startle_coef,
                                           lab = rownames(quare_startle_coef),
                                           x = 'estimate',
                                           y = 'p.value',
                                           xlim = c(-quare.startle.x.lim, quare.startle.x.lim),
                                           title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                           subtitle = 'Behavior = startle',
                                           caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                           xlab = 'Intx. β coefficients (group x RNA)',
                                           pCutoff = 0.05,
                                           FCcutoff = quare.fc.startle,
                                           # transcriptPointSize = 1.5,
                                           # transcriptLabSize = 3.0,
                                           colAlpha = 1,
                                           cutoffLineType = 'longdash',
                                           cutoffLineCol = 'grey50',
                                           cutoffLineWidth = 0.8,
                                           hline = quare.p.startle.0.3,
                                           hlineCol = 'red',
                                           hlineType = 'dotdash',
                                           hlineWidth = 0.8,
                                           gridlines.major = F,
                                           gridlines.minor = F)
      
    ## t) Print quare_startle_vol plots    
      print(quare.startle_vol)
      
    ## u) Save quare.startle_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_startle_vol.pdf", plot = quare.startle_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## v) Volcano Plot 
      quare.alunge_vol <- EnhancedVolcano(quare_alunge_coef,
                                          lab = rownames(quare_alunge_coef),
                                          x = 'estimate',
                                          y = 'p.value',
                                          xlim = c(-quare.alunge.x.lim, quare.alunge.x.lim),
                                          title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                          subtitle = 'Behavior = alunge',
                                          caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                          xlab = 'Intx. β coefficients (group x RNA)',
                                          pCutoff = 0.05,
                                          FCcutoff = quare.fc.alunge,
                                          # transcriptPointSize = 1.5,
                                          # transcriptLabSize = 3.0,
                                          colAlpha = 1,
                                          cutoffLineType = 'longdash',
                                          cutoffLineCol = 'grey50',
                                          cutoffLineWidth = 0.8,
                                          hline = quare.p.alunge.0.3,
                                          hlineCol = 'red',
                                          hlineType = 'dotdash',
                                          hlineWidth = 0.8,
                                          gridlines.major = F,
                                          gridlines.minor = F)
      
    ## w) Print quare_alunge_vol plots    
      print(quare.alunge_vol)
      
    ## x) Save quare.alunge_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_alunge_vol.pdf", plot = quare.alunge_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## y) Volcano Plot 
      quare.acontact_vol <- EnhancedVolcano(quare_acontact_coef,
                                            lab = rownames(quare_acontact_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-quare.acontact.x.lim, quare.acontact.x.lim),
                                            title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = acontact',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = quare.fc.acontact,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = quare.p.acontact.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## z) Print quare_acontact_vol plots    
      print(quare.acontact_vol)
      
    ## aa) Save quare.acontact_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_acontact_vol.pdf", plot = quare.acontact_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## bb) Volcano Plot 
      quare.swing_vol <- EnhancedVolcano(quare_swing_coef,
                                         lab = rownames(quare_swing_coef),
                                         x = 'estimate',
                                         y = 'p.value',
                                         xlim = c(-quare.swing.x.lim, quare.swing.x.lim),
                                         title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                         subtitle = 'Behavior = swing',
                                         caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                         xlab = 'Intx. β coefficients (group x RNA)',
                                         pCutoff = 0.05,
                                         FCcutoff = quare.fc.swing,
                                         # transcriptPointSize = 1.5,
                                         # transcriptLabSize = 3.0,
                                         colAlpha = 1,
                                         cutoffLineType = 'longdash',
                                         cutoffLineCol = 'grey50',
                                         cutoffLineWidth = 0.8,
                                         hline = quare.p.swing.0.3,
                                         hlineCol = 'red',
                                         hlineType = 'dotdash',
                                         hlineWidth = 0.8,
                                         gridlines.major = F,
                                         gridlines.minor = F)
      
    ## cc) Print quare_swing_vol plots    
      print(quare.swing_vol)
      
    ## dd) Save quare.swing_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_swing_vol.pdf", plot = quare.swing_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## ee) Volcano Plot 
      quare.athrust_vol <- EnhancedVolcano(quare_athrust_coef,
                                           lab = rownames(quare_athrust_coef),
                                           x = 'estimate',
                                           y = 'p.value',
                                           xlim = c(-quare.athrust.x.lim, quare.athrust.x.lim),
                                           title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                           subtitle = 'Behavior = athrust',
                                           caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                           xlab = 'Intx. β coefficients (group x RNA)',
                                           pCutoff = 0.05,
                                           FCcutoff = quare.fc.athrust,
                                           # transcriptPointSize = 1.5,
                                           # transcriptLabSize = 3.0,
                                           colAlpha = 1,
                                           cutoffLineType = 'longdash',
                                           cutoffLineCol = 'grey50',
                                           cutoffLineWidth = 0.8,
                                           hline = quare.p.athrust.0.3,
                                           hlineCol = 'red',
                                           hlineType = 'dotdash',
                                           hlineWidth = 0.8,
                                           gridlines.major = F,
                                           gridlines.minor = F)
      
    ## ff) Print quare_athrust_vol plots    
      print(quare.athrust_vol)
      
    ## gg) Save quare.athrust_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_athrust_vol.pdf", plot = quare.athrust_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## hh) Volcano Plot 
      quare.asigtime_vol <- EnhancedVolcano(quare_asigtime_coef,
                                            lab = rownames(quare_asigtime_coef),
                                            x = 'estimate',
                                            y = 'p.value',
                                            xlim = c(-quare.asigtime.x.lim, quare.asigtime.x.lim),
                                            title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                            subtitle = 'Behavior = asigtime',
                                            caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                            xlab = 'Intx. β coefficients (group x RNA)',
                                            pCutoff = 0.05,
                                            FCcutoff = quare.fc.asigtime,
                                            # transcriptPointSize = 1.5,
                                            # transcriptLabSize = 3.0,
                                            colAlpha = 1,
                                            cutoffLineType = 'longdash',
                                            cutoffLineCol = 'grey50',
                                            cutoffLineWidth = 0.8,
                                            hline = quare.p.asigtime.0.3,
                                            hlineCol = 'red',
                                            hlineType = 'dotdash',
                                            hlineWidth = 0.8,
                                            gridlines.major = F,
                                            gridlines.minor = F)
      
    ## ii) Print quare_asigtime_vol plots    
      print(quare.asigtime_vol)
      
    ## jj) Save quare.asigtime_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_asigtime_vol.pdf", plot = quare.asigtime_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## kk) Volcano Plot 
      quare.aswing_vol <- EnhancedVolcano(quare_aswing_coef,
                                          lab = rownames(quare_aswing_coef),
                                          x = 'estimate',
                                          y = 'p.value',
                                          xlim = c(-quare.aswing.x.lim, quare.aswing.x.lim),
                                          title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                          subtitle = 'Behavior = aswing',
                                          caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                          xlab = 'Intx. β coefficients (group x RNA)',
                                          pCutoff = 0.05,
                                          FCcutoff = quare.fc.aswing,
                                          # transcriptPointSize = 1.5,
                                          # transcriptLabSize = 3.0,
                                          colAlpha = 1,
                                          cutoffLineType = 'longdash',
                                          cutoffLineCol = 'grey50',
                                          cutoffLineWidth = 0.8,
                                          hline = quare.p.aswing.0.3,
                                          hlineCol = 'red',
                                          hlineType = 'dotdash',
                                          hlineWidth = 0.8,
                                          gridlines.major = F,
                                          gridlines.minor = F)
      
    ## ll) Print quare_aswing_vol plots    
      print(quare.aswing_vol)
      
    ## mm) Save quare.aswing_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_aswing_vol.pdf", plot = quare.aswing_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## nn) Volcano Plot 
      quare.center_vol <- EnhancedVolcano(quare_center_coef,
                                          lab = rownames(quare_center_coef),
                                          x = 'estimate',
                                          y = 'p.value',
                                          xlim = c(-quare.center.x.lim, quare.center.x.lim),
                                          title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                          subtitle = 'Behavior = center',
                                          caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                          xlab = 'Intx. β coefficients (group x RNA)',
                                          pCutoff = 0.05,
                                          FCcutoff = quare.fc.center,
                                          # transcriptPointSize = 1.5,
                                          # transcriptLabSize = 3.0,
                                          colAlpha = 1,
                                          cutoffLineType = 'longdash',
                                          cutoffLineCol = 'grey50',
                                          cutoffLineWidth = 0.8,
                                          hline = quare.p.center.0.3,
                                          hlineCol = 'red',
                                          hlineType = 'dotdash',
                                          hlineWidth = 0.8,
                                          gridlines.major = F,
                                          gridlines.minor = F)
      
    ## oo) Print quare_center_vol plots    
      print(quare.center_vol)
      
    ## pp) Save quare.center_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_center_vol.pdf", plot = quare.center_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## qq) Volcano Plot 
      quare.moving_vol <- EnhancedVolcano(quare_moving_coef,
                                          lab = rownames(quare_moving_coef),
                                          x = 'estimate',
                                          y = 'p.value',
                                          xlim = c(-quare.moving.x.lim, quare.moving.x.lim),
                                          title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                          subtitle = 'Behavior = moving',
                                          caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                          xlab = 'Intx. β coefficients (group x RNA)',
                                          pCutoff = 0.05,
                                          FCcutoff = quare.fc.moving,
                                          # transcriptPointSize = 1.5,
                                          # transcriptLabSize = 3.0,
                                          colAlpha = 1,
                                          cutoffLineType = 'longdash',
                                          cutoffLineCol = 'grey50',
                                          cutoffLineWidth = 0.8,
                                          hline = quare.p.moving.0.3,
                                          hlineCol = 'red',
                                          hlineType = 'dotdash',
                                          hlineWidth = 0.8,
                                          gridlines.major = F,
                                          gridlines.minor = F)
      
    ## rr) Print quare_moving_vol plots    
      print(quare.moving_vol)
      
    ## ss) Save quare.moving_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_moving_vol.pdf", plot = quare.moving_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## tt) Volcano Plot 
      quare.sigtime_vol <- EnhancedVolcano(quare_sigtime_coef,
                                           lab = rownames(quare_sigtime_coef),
                                           x = 'estimate',
                                           y = 'p.value',
                                           xlim = c(-quare.sigtime.x.lim, quare.sigtime.x.lim),
                                           title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                           subtitle = 'Behavior = sigtime',
                                           caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                           xlab = 'Intx. β coefficients (group x RNA)',
                                           pCutoff = 0.05,
                                           FCcutoff = quare.fc.sigtime,
                                           # transcriptPointSize = 1.5,
                                           # transcriptLabSize = 3.0,
                                           colAlpha = 1,
                                           cutoffLineType = 'longdash',
                                           cutoffLineCol = 'grey50',
                                           cutoffLineWidth = 0.8,
                                           hline = quare.p.sigtime.0.3,
                                           hlineCol = 'red',
                                           hlineType = 'dotdash',
                                           hlineWidth = 0.8,
                                           gridlines.major = F,
                                           gridlines.minor = F)
      
    ## uu) Print quare_sigtime_vol plots    
      print(quare.sigtime_vol)
      
    ## vv) Save quare.sigtime_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_sigtime_vol.pdf", plot = quare.sigtime_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
    ## ww) Volcano Plot 
      quare.surface_vol <- EnhancedVolcano(quare_surface_coef,
                                           lab = rownames(quare_surface_coef),
                                           x = 'estimate',
                                           y = 'p.value',
                                           xlim = c(-quare.surface.x.lim, quare.surface.x.lim),
                                           title = 'Quare Volcano Plot of Group x RNA Interaction Terms',
                                           subtitle = 'Behavior = surface',
                                           caption = 'FC cutoff, 4-SD of β; p-value cutoff, 0.05; BH FDR cutoff, red dotdashline',
                                           xlab = 'Intx. β coefficients (group x RNA)',
                                           pCutoff = 0.05,
                                           FCcutoff = quare.fc.surface,
                                           # transcriptPointSize = 1.5,
                                           # transcriptLabSize = 3.0,
                                           colAlpha = 1,
                                           cutoffLineType = 'longdash',
                                           cutoffLineCol = 'grey50',
                                           cutoffLineWidth = 0.8,
                                           hline = quare.p.surface.0.3,
                                           hlineCol = 'red',
                                           hlineType = 'dotdash',
                                           hlineWidth = 0.8,
                                           gridlines.major = F,
                                           gridlines.minor = F)
      
    ## xx) Print quare_surface_vol plots    
      print(quare.surface_vol)
      
    ## yy) Save quare.surface_vol plots
      # use ggsave to save the plot as pdf
      ggsave("quare_surface_vol.pdf", plot = quare.surface_vol,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE)
      
      
      
###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
  ### 5.1 Export data to an RData file     
    ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      # save(file = paste0(project_data_path, 'tidy_data_pred_rna_behav.RData'), 
      #      list = c('aripo_data', 'quare_data'))
      # 
      