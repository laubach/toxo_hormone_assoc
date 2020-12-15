###############################################################################
#############          Associations of Toxoplasma gondii          #############
#############             with steroid hormone levels             #############
#############                                                     #############
#############       3 Models: Toxo. and Neosp. associations       #############
#############                 w/plasma hormones                   #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 15 Oct 2020                 #############
#############              last updated: 15 Dec 2020              #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii infection with plasma 
             # testosterone and corticosterone levels in spotted hyenas

  ### NOTE: Cortisol data include only samples collected <= 13 minutes post
            # darting - measures baseline stress. Also only stress state 
            # categories 1 and 2 are included in analyses.
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Infection associations with testosterone
    # 4: Infection associations with corticosterone



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
 
      # load aod packages (used to for Wald test)
      library ('aod')
 
      # load car packages (used for type II and type III SS test and VIF)
      library ('car')
   
      # load dotwhisker packages; used with broom to graph beta estimates
      library ('dotwhisker')
      
      # load ggtext to use markdown syntax
      library ('ggtext')
      
    ## c) Modeling Packages
     
      # load broom packages
      library ('broom')
      
      # load broom.mixed, a package under development by Ben Bolker
      # similar to broom but extends bey lm and glm
      library ('broom.mixed')
   
      # load nlme packages
      library ('nlme')
  
      # load lme4 packages
      library ('lme4')
  
      # load boot packages; used to generate boot strap CI from LME4
      library ('boot')
   
      # load merTools packages
      library ('merTools')
      
      # load emmeans packages
      library ('emmeans')

        
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
    ## a) load RData: updated 3_neo_toxo_plasma_horm.RData joined to hyena 
      # data tables
      load(paste0(project_data_path, '3_neo_toxo_plasma_horm.RData'))
     
      
      
###############################################################################
##############   3: Infection associations with testosterone     ##############
############################################################################### 
      
  ### 3.1 Associations between T. gondii status and testosterone levels
    ## a) Unadjusted model: Male adults
      # testosterone by T. gondii infection
      T.toxo.unadj.mod.m.adult <- lm(t.ln ~ toxo.status, 
                              data = subset(plasma_horm_neosp_toxo_data,
                              sex == 'm' & 
                              age.cat.dart == 'adult' &
                                !is.na(x = toxo.status) &
                                !is.na(x = t.ln)))
      
      summary(T.toxo.unadj.mod.m.adult) # print model summary (ln scale)
      confint(T.toxo.unadj.mod.m.adult) # 95% CIs (ln scale)
      #plot(T.toxo.unadj.mod.m.adult) # view fitted vs residuals
      
      # Use emmeans to estimate marginal means
      T.toxo.unadj.mmean.m.adult <- emmeans(T.toxo.unadj.mod.m.adult, 
                                            'toxo.status')
      summary(T.toxo.unadj.mmean.m.adult)
      
      
      
  ### 3.2 Graph of testosterone levels by toxo      
      
    ## a) Unadjusted adult male testosterone model
      T.toxo.unadj.mmean.m.adult.tdy <- 
        broom::tidy(T.toxo.unadj.mmean.m.adult) %>%
        rename('term' = 'toxo.status') %>%
        relabel_predictors(negative = 'Negative',
                             positive = 'Positive')
      
    ## b) Add adult male model category variable
      T.toxo.unadj.mmean.m.adult.tdy$model <- 'Adult male model'
      
    ## c) # Note forcing dwplot to graph SE by specificying the upper and
      # lower bounds instead of automatically calculating actual
      # 95% CI from SE
      T.toxo.unadj.mmean.m.adult.tdy$conf.low <-
        T.toxo.unadj.mmean.m.adult.tdy$estimate - 
        T.toxo.unadj.mmean.m.adult.tdy$std.error
      
      T.toxo.unadj.mmean.m.adult.tdy$conf.high <-
        T.toxo.unadj.mmean.m.adult.tdy$estimate + 
        T.toxo.unadj.mmean.m.adult.tdy$std.error
      
    ## d) Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
      T.toxo.unadj.mmean.m.adult.plot <- dwplot(T.toxo.unadj.mmean.m.adult.tdy, 
             vline = geom_vline(xintercept = 0, colour = 'gray20', 
                                linetype = 2), # line at zero behind coefs
             dot_args = list(size = 3),
             whisker_args = list(size = 1),
             dodge_size = 1) + 
        coord_flip() + # flip x and y axes
        xlim (-3,1) +
        labs(title = expression(atop(paste('Associations between latent ', 
        italic('T. gondii '),'infection and'), 
        paste('plasma testosterone among adult male hyenas.')))) +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=18, face = 'bold')) +
        # theme(legend.justification = c(1,1), legend.position = c(1,0.25),
        #       legend.background = element_rect(fill = 'white'),
        #       legend.title = element_blank(), 
        #       legend.key = element_rect(fill = 'white')) +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        #theme(axis.line = element_line(colour = 'lightgrey', 
        #                               size = 1, linetype = 'solid')) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=18, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=18, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        #scale_color_grey (start = 0, end = 0) + # make color estimates black
        # color the dot and whiskers
        scale_color_manual(values=c('goldenrod3')) + 
        # NOTE: we flipped x and y axes above, so the 'xlab' is actually
        # 'ylab' and vice versa. 
        xlab(expression(atop(bold('Mean +/- SE'), 
                          paste(italic('Nat. Log. testosterone (ug/dL)'))))) +
        #scale_y_discrete(labels = c('Seropostive hyenas')) +
        ylab('')
    
      print(T.toxo.unadj.mmean.m.adult.plot)
      
    ## e) Save Plot
      # use ggsave to save the linearization plot
      ggsave('T_toxo_m_adult_plot.pdf', 
             plot = T.toxo.unadj.mmean.m.adult.plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 11,
             height = 6,
             units = c('in'), dpi = 300, limitsize = TRUE)
      

      
###############################################################################
##############   4: Infection associations with corticosterone   ##############
############################################################################### 
      
  ### 4.1 Associations between T. gondii status and corticosterone levels
    ## a) Unadjusted model: Female 
      # Corticosterone by T. gondii infection
      cort.toxo.unadj.mod.f <- lm(c.ln ~ toxo.status, 
                                    data = subset(plasma_horm_neosp_toxo_data,
                                                      sex == 'f' & 
                                                      dart.time.diff <= 13 & 
                                                      stressca <=2 &
                                                      !is.na(x = toxo.status) &
                                                      !is.na(x = c.ln)))
      
      summary(cort.toxo.unadj.mod.f) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.f) # 95% CIs (ln scale)
      #plot(cort.toxo.unadj.mod.f) # view fitted vs residuals
      
      # Use emmeans to estimate marginal means
      cort.toxo.unadj.mmean.f <- emmeans(cort.toxo.unadj.mod.f, 
                                               'toxo.status')
      
      summary(cort.toxo.unadj.mmean.f)
      
    ## b) Adjusted model: Female 
      # Corticosterone by T. gondii infection
      cort.toxo.adj.mod.f <- lm(c.ln ~ toxo.status + dart.am.pm, 
                                  data = subset(plasma_horm_neosp_toxo_data,
                                                    sex == 'f' & 
                                                    dart.time.diff <= 13 & 
                                                    stressca <=2 &
                                                    !is.na(x = toxo.status) &
                                                    !is.na(x = c.ln)))
      
      summary(cort.toxo.adj.mod.f) # print model summary (ln scale)
      confint(cort.toxo.adj.mod.f) # 95% CIs (ln scale)
      #plot(cort.toxo.adj.mod.f) # view fitted vs residuals
      
      # Use emmeans to estimate marginal means
      cort.toxo.adj.mmean.f <- emmeans(cort.toxo.adj.mod.f, 
                                             'toxo.status')
      
      summary(cort.toxo.adj.mmean.f)
      
  ### 4.1 Associations between T. gondii status and corticosterone levels
    ## a) Unadjusted model: Male 
      # Corticosterone by T. gondii infection
      cort.toxo.unadj.mod.m <- lm(c.ln ~ toxo.status, 
                                     data = subset(plasma_horm_neosp_toxo_data,
                                            sex == 'm' & 
                                            dart.time.diff <= 13 & 
                                            stressca <=2 &
                                            !is.na(x = toxo.status) &
                                            !is.na(x = c.ln)))
      
      summary(cort.toxo.unadj.mod.m) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.m) # 95% CIs (ln scale)
      #plot(cort.toxo.unadj.mod.m) # view fitted vs residuals
      
      # Use emmeans to estimate marginal means
      cort.toxo.unadj.mmean.m <- emmeans(cort.toxo.unadj.mod.m, 
                                            'toxo.status')
      
      summary(cort.toxo.unadj.mmean.m)
      
    ## b) Adjusted model: Male 
      # Corticosterone by T. gondii infection
      cort.toxo.adj.mod.m <- lm(c.ln ~ toxo.status + dart.am.pm, 
                                    data = subset(plasma_horm_neosp_toxo_data,
                                        sex == 'm' &
                                        dart.time.diff <= 13 &
                                        stressca <=2 &
                                        !is.na(x = toxo.status) &
                                        !is.na(x = c.ln)))
      
      summary(cort.toxo.adj.mod.m) # print model summary (ln scale)
      confint(cort.toxo.adj.mod.m) # 95% CIs (ln scale)
      #plot(cort.toxo.adj.mod.m) # view fitted vs residuals
      
      # Use emmeans to estimate marginal means
      cort.toxo.adj.mmean.m <- emmeans(cort.toxo.adj.mod.m, 
                                               'toxo.status')
      
      summary(cort.toxo.adj.mmean.m)
      
    
      
      
      
      
      
      
      
      
  