###############################################################################
##############    Guppy Predation Environment, RNA expression,   ##############
##############      and Behavior: MANOVA effect modification     ##############
##############                 By: Zach Laubach                  ##############
##############               created: 9 Sept 2020                ##############
##############             last updated: 9 Sept 2020             ##############
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
      
    ## WORKAROUND: https://github.com/rstudio/rstudio/issues/6692
    ## Revert to 'sequential' setup of PSOCK cluster in 
    ## RStudio Console on macOS and R 4.0.0
      if (Sys.getenv("RSTUDIO") == "1" && !nzchar(Sys.getenv("RSTUDIO_TERM")) && 
          Sys.info()["sysname"] == "Darwin" && getRversion() >= "4.0.0") {
        parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
      }
  

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
      # Check for car and install if not already installed
      # includes vif function
      if (!'car' %in% installed.packages()[,1]){
        install.packages ('car')
      }
      # load car packages
      library ('car')

      # Check for qvalue and install if not already installed
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      
      BiocManager::install("qvalue")
      # load qvalue package
      library ('qvalue')
      
      # Check for MVLM and install if not already installed
      if (!'MVLM' %in% installed.packages()[,1]){
        install.packages ('MVLM')
      }
      # load MVLM packages
      library ('MVLM')
      
      # Check for MANOVA.RM and install if not already installed
      if (!'MANOVA.RM' %in% installed.packages()[,1]){
        install.packages ('MANOVA.RM')
      }
      # load MANOVA.RM packages
      library ('MANOVA.RM')
      
        
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
      
    ## b) Remove missing data
      aripo_data <- aripo_data[complete.cases(aripo_data), ]
     
    ## b) Transform aripo_data data from wide to long based on RNA expression
      aripo_rna_long <- aripo_data %>%
        gather(key = 'site', value = 'rna', 
               'ENSPREG00000006268':'ENSPREG00000007674', 
               factor_key = F)
      
      
      test.x <- aripo_data %>%
        select(ENSPREG00000006268:ENSPREG00000006273) 
     
    ## c) RNA expression site by site MANOVA for aripo_rna_long
      # Estimate the RNA expression site by group interaction
      # all behaviors - outcome and tidy lm objects using 'broom'
      aripo_behaviors_coef <- aripo_rna_long%>% 
        nest(data = c(fish, family, week, pop, rear, group, basin, sigtime, sigmoid, 
                      thrust, swing, contact, posturing, asigtime, asigmoid, athrust, 
                      aswing, acontact, afaceoff, alunge, startle, surface, center, 
                      moving, rna)) %>% 
        mutate(model = map(data, 
                           ~manova(lm(cbind(sigtime, sigmoid, thrust, swing, 
                                            contact,posturing, asigtime, 
                                            asigmoid, athrust, aswing,
                                            acontact, afaceoff, alunge, 
                                            startle, surface, 
                                            center, moving) ~ rna*group, 
                                     .))), 
               tidy_data = map(model, broom::tidy)) %>% 
        select(site, tidy_data) %>% 
        unnest()
     
     
    ## d) Check if there are 4 estimates (intercept, rna, group & intx) per model
      summry_aaripo_behaviors_coef <- aripo_behaviors_coef  %>%
        group_by(site) %>%
        summarize(n.terms = sum(!is.na(term)))
      
    ## e) filter to remove intercepts and retain only intx coefs and test stats
      aripo_behaviors_coef <- filter(aripo_behaviors_coef, 
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
      # all behaviors - outcome and tidy lm objects using 'broom'
      quare_behaviors_coef<- quare_rna_long %>%
        group_by(site) %>%
        do(tidy(manova(cbind(sigtime, sigmoid, thrust, swing, contact,
                             posturing, asigtime, asigmoid, athrust, aswing,
                             acontact, afaceoff, alunge, startle, surface, 
                             center, moving) ~ rna*group, data = .)))
      
    ## d) Check if there are 8 estimates (intercept, rna, group & intx) per model
      summry_quare_behaviors_coef <- quare_behaviors_coef  %>%
        group_by(site) %>%
        summarize(n.terms = sum(!is.na(term)))
      
    ## e) filter to remove intercepts and retain only intx coefs and test stats
      quare_behaviors_coef <- filter(quare_behaviors_coef, 
                                   grepl('rna:', term)) 
      
      
      
  ## 4.3 Visualize Aripo distribution of interaction terms   
    ## a) Calculate the min and max estimates for x-axis scale
      aripo.behaviors.x.lim <- round(max(abs(aripo_behaviors_coef$pillai), 
                               na.rm = T), 1)
    ## b) Histogram Aripo rna x group intx terms
      aripo_behaviors_intx_hist <- ggplot(data = aripo_behaviors_coef, 
                                   aes(x = pillai)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-aripo.behaviors.x.lim, 
                                    aripo.behaviors.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(0, aripo.behaviors.x.lim)) +
        geom_vline(aes(xintercept = mean(pillai, na.rm = T)),
                   color = "blue", linetype = "dashed", size = 1) + 
        labs(title = "Aripo Histogram of RNA expression by Group Interaction 
Terms (Outcome = behaviors)") +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x = "Beta Intx. Terms", y = "Frequency")
      
    ## c) Print aripo_behaviors_intx_hist plots
      print(aripo_behaviors_intx_hist)
      
    ## d) Save aripo_behaviors_intx_hist plots
      # use ggsave to save the plot as pdf
      ggsave("aripo_behaviors_intx_hist.pdf", plot = aripo_behaviors_intx_hist,
             device = NULL, 
             path = paste0(here(),"/output"), 
             scale = 1, width = 8, height = 7,
             units = c("in"), dpi = 300, limitsize = TRUE) 
      
    # ## e) Density plot Aripo rna x group intx terms
    #   ggplot(data = aripo_behaviors_coef, aes(x = pillai)) + 
    #     geom_density() +
    #     geom_vline(aes(xintercept = mean(pillai)),
    #                color = "blue", linetype = "dashed", size=1)
    #   
      
  
      
      
  ## 4.4 Visualize Quare distribution of interaction terms
    ## a) Calculate the min and max estimates for x-axis scale
      quare.contact.x.lim <- round(max(abs(quare_contact_coef$pillai), 
                                       na.rm = T), 1)
    ## b) Histogram Quare rna x group intx terms
      quare_contact_intx_hist <- ggplot(data = quare_contact_coef, 
                                        aes(x = pillai)) + 
        geom_histogram(aes(y = ..count..),
                       breaks = seq(-quare.contact.x.lim, 
                                    quare.contact.x.lim, by = 0.025), 
                       col = "black",
                       fill = "dark grey") +
        xlim(c(-quare.contact.x.lim, quare.contact.x.lim)) +
        geom_vline(aes(xintercept = mean(pillai, na.rm = T)),
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
      #   ggplot(data = quare_contact_coef, aes(x = pillai)) + 
      #     geom_density() +
      #     geom_vline(aes(xintercept = mean(pillai)),
      #                color = "blue", linetype = "dashed", size=1)
      #   
      
      
  
  ### 4.5  Aripo Mutliple comparisons corrections
    ## a) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      aripo_contact_coef$fdr.p = p.adjust(aripo_contact_coef$p.value, 
                                           method = "BH", 
                                           length(aripo_contact_coef$p.value))
     
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
      
    
    

  ### 4.7 Quare Mutliple comparisons corrections
    ## a) Use p.adjust to calculate the Benjamini-Hocherberg FDR correction 
      quare_contact_coef$fdr.p = p.adjust(quare_contact_coef$p.value, 
                                          method = "BH", 
                                          length(quare_contact_coef$p.value))
      
    
      
      
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
      