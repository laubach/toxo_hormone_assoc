###############################################################################
#############          Associations of Toxoplasma gondii          #############
#############             with steroid hormone levels             #############
#############                                                     #############
#############       3 Models: Toxo. and Neosp. associations       #############
#############                 w/plasma hormones                   #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 15 Oct 2020                 #############
#############             last updated: 24 Nov 2021              #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii infection with plasma 
             # testosterone and cortisol levels in spotted hyenas

  ### NOTE: Cortisol data include only samples collected <= 13 minutes post
            # darting - measures baseline stress. Also only stress state 
            # categories 1 and 2 are included in analyses.
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Infection associations with testosterone
    # 4: Infection associations with cortisol



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
    
      
  ### 1.6 Set seed for reproducible bootstrapping
      round(runif(1, min = 1, max = 999999),0) #generated 20 May 2021
      set.seed(512158)
      
      
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
      
      
  ### 3.1 Functions for extracting model coefficeints to be used in bootstrap
    ## a) function to extract lm beta estimates pass to car:Boot since data
      # are not normal even after ln transformation
      toxo_reg_coef <- function(mod) {
        reg_mod <- mod # run lm mod
        tdy_reg_mod <- coefficients(reg_mod) # extract estimates in data frame
        tdy_reg_mod # print estimates
      
      }
      
      
      # # test fucntion:
      # toxo_reg_coef(T.toxo.mod.m.cubsub)
      
      
    ## b) function to extract marginal means and pass to car:Boot since data
      # are not normal even after ln transformation
      toxo_mm_coef <- function(mod) {
        toxo.mod.mmean <- emmeans(mod, 'toxo.status') # estimate marg. means
        var.name <- c("Uninfected", "Infected")
        # extract marg.mean
        est.mmean <- summary(toxo.mod.mmean)$emmean
        names(est.mmean) <- var.name # rename variables
        est.mmean # print marg. mean
      }
      
      # # test fucntion:
      # toxo_mm_coef(T.toxo.mod.m.cubsub)
      
      
  ### 3.2 Associations between T. gondii status and testosterone levels
      # cub and subadult males
    ## a) Unadjusted model:testosterone by T. gondii infection
      T.toxo.mod.m.cubsub <- lm(t.ln
                                #t.imp.ln # sensitivity analyses imputed data
                                ~ toxo.status, 
                              data = subset(plasma_horm_neosp_toxo_data,
                              sex == 'm' & 
                              age.cat.dart != 'adult' &
                                !is.na(x = toxo.status) &
                                !is.na(x = t.ln)))
      
      summary(T.toxo.mod.m.cubsub) # print model summary (ln scale)
      confint(T.toxo.mod.m.cubsub) # 95% CIs (ln scale)
      # plot(T.toxo.mod.m.cubsub) # view fitted vs residuals
 
    ## b) Use emmeans to estimate marginal means
      T.toxo.mod.mmean.m.cubsub <- emmeans(T.toxo.mod.m.cubsub, 
                                            'toxo.status')
      summary(T.toxo.mod.mmean.m.cubsub)
      

    ## c) Bootstrap parameter lm estimates   
      # bootstrapping number of resampling simulations
      set.seed(512158)
      lm.boot.T.toxo.mod.m.cubsub <- car::Boot(obj = T.toxo.mod.m.cubsub, 
                                       f = toxo_reg_coef, R = 2000)
     
      # create a summary statistics data frame
      lm.summ.boot.T.toxo.mod.m.cubsub <- summary(lm.boot.T.toxo.mod.m.cubsub)
      
      # estimate permuted confidence intervals 
      lm.boot.ci.T.toxo.mod.m.cubsub <- confint(lm.boot.T.toxo.mod.m.cubsub, 
                                                type = "perc")
      
    ## d) Tidy bootstrapped regression estimates table 
      # create a data frame for graphing
      lm.boot.T.toxo.mod.m.cubsub <- cbind(lm.summ.boot.T.toxo.mod.m.cubsub, 
                                           lm.boot.ci.T.toxo.mod.m.cubsub)
      
      # label the estimates in data frame
      lm.boot.T.toxo.mod.m.cubsub$model <- c("Cub/Subadult Male")
      
      # back transform estimates of interest
      lm.boot.T.toxo.mod.m.cubsub$estimate <- 
        exp(lm.boot.T.toxo.mod.m.cubsub$bootMed)
      lm.boot.T.toxo.mod.m.cubsub$SE <- 
        exp(lm.boot.T.toxo.mod.m.cubsub$bootSE)
      lm.boot.T.toxo.mod.m.cubsub$low.conf <- 
        exp(lm.boot.T.toxo.mod.m.cubsub$`2.5 %`)
      lm.boot.T.toxo.mod.m.cubsub$hi.conf <- 
        exp(lm.boot.T.toxo.mod.m.cubsub$`97.5 %`)
      
    ## e) Bootstrap parameter marginal mean estimates   
      # bootstrapping number of resampling simulations
      # set.seed(512158)
      # mm.boot.T.toxo.mod.m.cubsub <- car::Boot(obj = T.toxo.mod.m.cubsub, 
      #                                       f = toxo_mm_coef, R = 2000)
      # 
      # # create a summary statistics data frame
      # mm.summ.boot.T.toxo.mod.m.cubsub <- summary(mm.boot.T.toxo.mod.m.cubsub)
      # 
      # # estimate permuted confidence intervals 
      # mm.boot.ci.T.toxo.mod.m.cubsub <- confint(mm.boot.T.toxo.mod.m.cubsub, 
      #                                        type = "perc")
      
    
      
    ## f) Logistic regression: cub and subadult males
      # t.bin (above vs below median) by toxo_status
      t.log.m.cubsub <- glm(t.bin ~ toxo.status, 
                         data = subset(plasma_horm_neosp_toxo_data,
                                  sex == 'm' & 
                                  age.cat.dart != 'adult' &
                                  !is.na(x = toxo.status) &
                                  !is.na(x = t.ln)),
                         family = binomial) 
      
      summary(t.log.m.cubsub)    # model summary (log odds scale)
      confint(t.log.m.cubsub)    # 95% CIs (log odds scale)
   
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(t.log.m.cubsub), confint (t.log.m.cubsub)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(t.log.m.cubsub), Sigma = vcov(t.log.m.cubsub), 
                estimates = 2) 
      
      
  ### 3.3 Associations between T. gondii status and testosterone levels
      # adult males
    ## a) Unadjusted model:testosterone by T. gondii infection
      T.toxo.mod.m.adult <- lm(t.ln
                               #t.imp.ln # sensitivity analyses imputed data
                               ~ toxo.status, 
                                     data = subset(plasma_horm_neosp_toxo_data,
                                                   sex == 'm' & 
                                                     age.cat.dart == 'adult' &
                                                     !is.na(x = toxo.status) &
                                                     !is.na(x = t.ln)))
      
      summary(T.toxo.mod.m.adult) # print model summary (ln scale)
      confint(T.toxo.mod.m.adult) # 95% CIs (ln scale)
      #plot(T.toxo.mod.m.adult) # view fitted vs residuals
      
    ## b) Use emmeans to estimate marginal means
      T.toxo.mmean.m.adult <- emmeans(T.toxo.mod.m.adult, 
                                            'toxo.status')
      summary(T.toxo.mmean.m.adult)
      
    ## c) Bootstrap parameter lm estimates   
      # bootstrapping number of resampling simulations
      set.seed(512158)
      lm.boot.T.toxo.mod.m.adult <- car::Boot(obj = T.toxo.mod.m.adult, 
                                               f = toxo_reg_coef, R = 2000)
      
      # create a summary statistics data frame
      lm.summ.boot.T.toxo.mod.m.adult <- summary(lm.boot.T.toxo.mod.m.adult)
      
      # estimate permuted confidence intervals 
      lm.boot.ci.T.toxo.mod.m.adult <- confint(lm.boot.T.toxo.mod.m.adult, 
                                               type = "perc")
      
    ## d) Tidy  bootstrapped regression estimates table 
      # create a data frame for graphing
      lm.boot.T.toxo.mod.m.adult <- cbind(lm.summ.boot.T.toxo.mod.m.adult, 
                                          lm.boot.ci.T.toxo.mod.m.adult)
      
      # label the estimates in data frame
      lm.boot.T.toxo.mod.m.adult$model <- c("Adult Male")
      
      # back transform estimates of interest
      lm.boot.T.toxo.mod.m.adult$estimate <- 
        exp(lm.boot.T.toxo.mod.m.adult$bootMed)
      lm.boot.T.toxo.mod.m.adult$SE <- 
        exp(lm.boot.T.toxo.mod.m.adult$bootSE)
      lm.boot.T.toxo.mod.m.adult$low.conf <- 
        exp(lm.boot.T.toxo.mod.m.adult$`2.5 %`)
      lm.boot.T.toxo.mod.m.adult$hi.conf <- 
        exp(lm.boot.T.toxo.mod.m.adult$`97.5 %`)
      
    # ## e) Bootstrap parameter marginal mean estimates   
    #   # bootstrapping number of resampling simulations
    #   set.seed(512158)
    #   mm.boot.T.toxo.mod.m.adult <- car::Boot(obj = T.toxo.mod.m.adult, 
    #                                            f = toxo_mm_coef, R = 2000)
    #   
    #   # create a summary statistics data frame
    #   mm.summ.boot.T.toxo.mod.m.adult <- summary(mm.boot.T.toxo.mod.m.adult)
    #   
    #   # estimate permuted confidence intervals 
    #   mm.boot.ci.T.toxo.mod.m.adult <- confint(mm.boot.T.toxo.mod.m.adult, 
    #                                             type = "perc")
      
    
      
    ## f) Logistic regression: adult male
      # t.bin (above vs below median) by toxo_status
      t.log.m.adult <- glm(t.bin ~ toxo.status, 
                   data = subset(plasma_horm_neosp_toxo_data,
                                 sex == 'm' & 
                                   age.cat.dart == 'adult' &
                                   !is.na(x = toxo.status) &
                                   !is.na(x = t.ln)),
                   family = binomial) 
      
      summary(t.log.m.adult)    # model summary (log odds scale)
      confint(t.log.m.adult)    # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(t.log.m.adult), confint (t.log.m.adult)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(t.log.m.adult), Sigma = vcov(t.log.m.adult), 
                estimates = 2) 

      
  ### 3.3 Associations between T. gondii status and testosterone levels
      # cub/subadult females
    ## a) Unadjusted model:testosterone by T. gondii infection   
      T.toxo.mod.f.cubsub <- lm(t.ln
                                #t.imp.ln # sensitivity analyses imputed data
                                ~ toxo.status, 
                                data = subset(plasma_horm_neosp_toxo_data,
                                              sex == 'f' & 
                                                age.cat.dart != 'adult' &
                                                !is.na(x = toxo.status) &
                                                !is.na(x = t.ln)))
      
      summary(T.toxo.mod.f.cubsub) # print model summary (ln scale)
      confint(T.toxo.mod.f.cubsub) # 95% CIs (ln scale)
      #plot(T.toxo.mod.f.cubsub) # view fitted vs residuals
  
      
    ## b) Use emmeans to estimate marginal means
      T.toxo.mod.mmean.f.cubsub <- emmeans(T.toxo.mod.f.cubsub, 
                                           'toxo.status')
      summary(T.toxo.mod.mmean.f.cubsub)
      
    ## c) Bootstrap parameter lm estimates   
      # bootstrapping number of resampling simulations
      set.seed(512158)
      lm.boot.T.toxo.mod.f.cubsub <- car::Boot(obj = T.toxo.mod.f.cubsub, 
                                               f = toxo_reg_coef, R = 2000)
      
      # create a summary statistics data frame
      lm.summ.boot.T.toxo.mod.f.cubsub <- summary(lm.boot.T.toxo.mod.f.cubsub)
      
      # estimate permuted confidence intervals 
      lm.boot.ci.T.toxo.mod.f.cubsub <- confint(lm.boot.T.toxo.mod.f.cubsub, 
                                                type = "perc")
      
    ## d) Tidy bootstrapped marginal means results table 
      # create a data frame for graphing
      lm.boot.T.toxo.mod.f.cubsub <- cbind(lm.summ.boot.T.toxo.mod.f.cubsub, 
                                           lm.boot.ci.T.toxo.mod.f.cubsub)
      
      # label the estimates in data frame
      lm.boot.T.toxo.mod.f.cubsub$model <- c("Cub/Subadult Female")
      
      # back transform estimates of interest
      lm.boot.T.toxo.mod.f.cubsub$estimate <- 
        exp(lm.boot.T.toxo.mod.f.cubsub$bootMed)
      lm.boot.T.toxo.mod.f.cubsub$SE <- 
        exp(lm.boot.T.toxo.mod.f.cubsub$bootSE)
      lm.boot.T.toxo.mod.f.cubsub$low.conf <- 
        exp(lm.boot.T.toxo.mod.f.cubsub$`2.5 %`)
      lm.boot.T.toxo.mod.f.cubsub$hi.conf <- 
        exp(lm.boot.T.toxo.mod.f.cubsub$`97.5 %`)
      
    ## e) Bootstrap parameter marginal mean estimates   
      # # bootstrapping number of resampling simulations
      # set.seed(512158)
      # mm.boot.T.toxo.mod.f.cubsub <- car::Boot(obj = T.toxo.mod.f.cubsub, 
      #                                          f = toxo_mm_coef, R = 2000)
      # 
      # # create a summary statistics data frame
      # mm.summ.boot.T.toxo.mod.f.cubsub <- summary(mm.boot.T.toxo.mod.f.cubsub)
      # 
      # # estimate permuted confidence intervals 
      # mm.boot.ci.T.toxo.mod.f.cubsub <- confint(mm.boot.T.toxo.mod.f.cubsub, 
      #                                           type = "perc")
      
    ## f) Logistic regression: cub and subadult females
      # t.bin (above vs below median) by toxo_status
      t.log.f.cubsub <- glm(t.bin ~ toxo.status, 
                            data = subset(plasma_horm_neosp_toxo_data,
                                          sex == 'f' & 
                                            age.cat.dart != 'adult' &
                                            !is.na(x = toxo.status) &
                                            !is.na(x = t.ln)),
                            family = binomial) 
      
      summary(t.log.f.cubsub)    # model summary (log odds scale)
      confint(t.log.f.cubsub)    # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(t.log.f.cubsub), confint (t.log.f.cubsub)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(t.log.f.cubsub), Sigma = vcov(t.log.f.cubsub), 
                estimates = 2) 
      
    
      
  ### 3.4 Graph of testosterone levels by toxo (regession coefficients)      
    ## a) Combine regression estimates into a tidy table
      test.toxo.lm.est <- bind_rows(lm.boot.T.toxo.mod.f.cubsub, 
                                       lm.boot.T.toxo.mod.m.cubsub, 
                                       lm.boot.T.toxo.mod.m.adult)
    
    ## b) Add a variable, 'term' to distinguish estimates 
      test.toxo.lm.est$term <- c('intercept', 'infected vs. uninfected', 
                                 'intercept', 'infected vs. uninfected',
                                 'intercept', 'infected vs. uninfected')
     
    ## c) Remove intercepts    
      test.toxo.lm.est <- test.toxo.lm.est %>%
        filter(term != 'intercept')
    

    ## d) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'model' variable 
      test.toxo.lm.est <- 
        transform(test.toxo.lm.est, 
                  model = factor(model,
                                 levels = c('Cub/Subadult Female', 
                                            'Cub/Subadult Male',
                                            'Adult Male')))  
    
    ## e) Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
      test.toxo.lm.est.plot <- 
        ggplot(test.toxo.lm.est, aes(x = model, y = estimate)) +
        geom_hline(yintercept = 1, colour = 'red',
                                linetype = 2) + # line at null behind coefs
        geom_point(size = 6, color = 'red', alpha = 0.8, 
                   position=position_dodge(width = 0.5)) +
        geom_errorbar(aes(ymin=low.conf, ymax=hi.conf), width=.1, color = 'red',
                      position=position_dodge(.5)) +
        #coord_flip() + # flip x and y axes
        labs(title = 'Sex and age stratified associations of T.gondii infection status 
             and plasma testosterone in spotted hyenas',
             subtitle = 'Geometric mean ratio and 95% CI are based on percentile bootstrap (2000 simulations)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=20, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        #theme(axis.line = element_line(colour = 'lightgrey', 
        #                               size = 1, linetype = 'solid')) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=20, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=20, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
    
        xlab(expression(italic("(Sex and age stratified models)"))) +
        ylab(expression(atop(bold("Testosterone geometric mean ratio and 95% CI"), 
        paste(italic("Infected vs. uninfected (reference)"))))) 
 
      
      print(test.toxo.lm.est.plot)
      
    ## f) Save Plot
      # use ggsave to save the linearization plot
      ggsave('test.toxo.lm.est.plot.pdf', plot = test.toxo.lm.est.plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 12,
             height = 9,
             units = c('in'), dpi = 300, limitsize = TRUE)
      
  #   ## a) Unadjusted adult male testosterone model
  #     T.toxo.mod.m.cubsub.tdy <- 
  #       broom::tidy(T.toxo.mod.mmean.m.cubsub) %>%
  #       rename('estimate' = 'toxo.status') %>%
  #       relabel_predictors(negative = 'Negative',
  #                            positive = 'Positive')
  #     
  #   ## b) Add adult male model category variable
  #     T.toxo.mod.m.cubsub.tdy$model <- 'cub and subadult male model'
  #     
  # 
  #     
  #   ## d) Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
  #     T.toxo.unadj.mmean.m.adult.plot <- dwplot(T.toxo.unadj.mmean.m.adult.tdy, 
  #            vline = geom_vline(xintercept = 0, colour = 'gray20', 
  #                               linetype = 2), # line at zero behind coefs
  #            dot_args = list(size = 3),
  #            whisker_args = list(size = 1),
  #            dodge_size = 1) + 
  #       coord_flip() + # flip x and y axes
  #       xlim (-3,1) +
  #       labs(title = expression(atop(paste('Associations between latent ', 
  #       italic('T. gondii '),'infection and'), 
  #       paste('plasma testosterone among adult male hyenas.')))) +
  #       theme(plot.title = element_text(hjust = 0.5)) + # center title
  #       theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
  #       # bold and size title and axes labels
  #       theme(text = element_text(size=18, face = 'bold')) +
  #       # theme(legend.justification = c(1,1), legend.position = c(1,0.25),
  #       #       legend.background = element_rect(fill = 'white'),
  #       #       legend.title = element_blank(), 
  #       #       legend.key = element_rect(fill = 'white')) +
  #       theme(axis.ticks = element_blank()) + # remove axis ticks
  #       # remove background color
  #       theme(panel.background = element_rect(fill = 'white')) +
  #       # add major axes
  #       #theme(axis.line = element_line(colour = 'lightgrey', 
  #       #                               size = 1, linetype = 'solid')) + 
  #       # change axes font style, color, size, angle, and margin
  #       theme(axis.text.x = element_text(face='bold', color='black', 
  #                                        size=18, angle=0,
  #                                        margin = margin(t = 10, r = 0, 
  #                                                        b = 10, l = 0)),
  #             axis.text.y = element_text(face='bold', color='black', 
  #                                        size=18, angle=0, 
  #                                        margin = margin(t = 0, r = 0, 
  #                                                        b = 0, l = 10))) +
  #       #scale_color_grey (start = 0, end = 0) + # make color estimates black
  #       # color the dot and whiskers
  #       scale_color_manual(values=c('goldenrod3')) + 
  #       # NOTE: we flipped x and y axes above, so the 'xlab' is actually
  #       # 'ylab' and vice versa. 
  #       xlab(expression(atop(bold('Mean +/- SE'), 
  #                         paste(italic('Nat. Log. testosterone'))))) +
  #       #scale_y_discrete(labels = c('Seropostive hyenas')) +
  #       ylab('')
  #   
  #     print(T.toxo.unadj.mmean.m.adult.plot)
  #     
  #   ## e) Save Plot
  #     # use ggsave to save the linearization plot
  #     ggsave('T_toxo_m_adult_plot.pdf', 
  #            plot = T.toxo.unadj.mmean.m.adult.plot, 
  #            device = NULL,
  #            path = paste0(here(),'/output'), 
  #            scale = 1, width = 11,
  #            height = 6,
  #            units = c('in'), dpi = 300, limitsize = TRUE)
  #     

      
###############################################################################
##############       4: Infection associations with cortisol     ##############
############################################################################### 
      
  ### 4.1 Associations between T. gondii status and cortisol levels
    ## a) Unadjusted model: Female 
      # cortisol by T. gondii infection
      cort.toxo.unadj.mod.f <- lm(c.ln
                                  #c.imp.ln # sensitivity analyses imputed data
                                  ~ toxo.status,
                                    data = subset(plasma_horm_neosp_toxo_data,
                                                      sex == 'f' & 
                                                      dart.time.diff <= 13 & 
                                                      stressca <=2 &
                                                      !is.na(x = spratio) &
                                                      !is.na(x = c.ln)))
      
      summary(cort.toxo.unadj.mod.f) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.f) # 95% CIs (ln scale)
      #plot(cort.toxo.unadj.mod.f) # view fitted vs residuals
      
    ## b) Use emmeans to estimate marginal means
      cort.toxo.unadj.mmean.f <- emmeans(cort.toxo.unadj.mod.f, 
                                               'toxo.status')
      
      summary(cort.toxo.unadj.mmean.f)

    ## c) Bootstrap parameter lm estimates   
      # bootstrapping number of resampling simulations
      set.seed(512158)
      lm.boot.cort.toxo.unadj.mod.f <- car::Boot(obj = cort.toxo.unadj.mod.f, 
                                               f = toxo_reg_coef, R = 2000)
      
      # create a summary statistics data frame
      lm.summ.boot.cort.toxo.unadj.mod.f <- summary(lm.boot.cort.toxo.unadj.mod.f)
      
      # estimate permuted confidence intervals 
      lm.boot.ci.cort.toxo.unadj.mod.f <- confint(lm.boot.cort.toxo.unadj.mod.f, 
                                                type = "perc")
      
    ## d) Tidy bootstrapped regression estimates table 
      # create a data frame for graphing
      lm.boot.cort.toxo.unadj.mod.f <- cbind(lm.summ.boot.cort.toxo.unadj.mod.f, 
                                           lm.boot.ci.cort.toxo.unadj.mod.f)
      
      # label the estimates in data frame
      lm.boot.cort.toxo.unadj.mod.f$model <- c("Female")
      
      # back transform estimates of interest
      lm.boot.cort.toxo.unadj.mod.f$estimate <- 
        exp(lm.boot.cort.toxo.unadj.mod.f$bootMed)
      lm.boot.cort.toxo.unadj.mod.f$SE <- 
        exp(lm.boot.cort.toxo.unadj.mod.f$bootSE)
      lm.boot.cort.toxo.unadj.mod.f$low.conf <- 
        exp(lm.boot.cort.toxo.unadj.mod.f$`2.5 %`)
      lm.boot.cort.toxo.unadj.mod.f$hi.conf <- 
        exp(lm.boot.cort.toxo.unadj.mod.f$`97.5 %`)
      
    ## e) Bootstrap parameter marginal mean estimates   
      # bootstrapping number of resampling simulations
      # set.seed(512158)
      # mm.boot.cort.toxo.unadj.mod.fb <- car::Boot(obj = cort.toxo.unadj.mod.f, 
      #                                       f = toxo_mm_coef, R = 2000)
      # 
      # # create a summary statistics data frame
      # mm.summ.boot.cort.toxo.unadj.mod.f <- summary(mm.boot.cort.toxo.unadj.mod.f)
      # 
      # # estimate permuted confidence intervals 
      # mm.boot.ci.cort.toxo.unadj.mod.f <- confint(mm.boot.cort.toxo.unadj.mod.f, 
      #                                        type = "perc")
      
    ## f) Adjusted model: Female 
      # cortisol by T. gondii infection
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

    ## g) Use emmeans to estimate marginal means
      cort.toxo.adj.mmean.f <- emmeans(cort.toxo.adj.mod.f, 
                                             'toxo.status')
 
### NOTE: Not enough PM meausrements to bootstrap estimates           
    # ## h) Bootstrap parameter lm estimates   
    #   # bootstrapping number of resampling simulations
    #   set.seed(512158)
    #   lm.boot.cort.toxo.adj.mod.f <- car::Boot(obj = cort.toxo.adj.mod.f,
    #                                              f = toxo_reg_coef, R = 2000)
    # 
    #   # create a summary statistics data frame
    #   lm.summ.boot.cort.toxo.adj.mod.f <- summary(lm.boot.cort.toxo.adj.mod.f)
    # 
    #   # estimate permuted confidence intervals
    #   lm.boot.ci.cort.toxo.adj.mod.f <- confint(lm.boot.cort.toxo.adj.mod.f,
    #                                               type = "perc")
    # 
    # ## i) Tidy bootstrapped regression estimates table
    #   # create a data frame for graphing
    #   lm.boot.cort.toxo.adj.mod.f <- cbind(lm.summ.boot.cort.toxo.adj.mod.f,
    #                                          lm.boot.ci.cort.toxo.adj.mod.f)

    #   # label the estimates in data frame
    #   lm.boot.cort.toxo.adj.mod.f$model <- c("Female")
    #   
    #   # back transform estimates of interest
    #   lm.boot.cort.toxo.adj.mod.f$estimate <- 
    #     exp(lm.boot.cort.toxo.adj.mod.f$bootMed)
    #   lm.boot.cort.toxo.adj.mod.f$SE <- 
    #     exp(lm.boot.cort.toxo.adj.mod.f$bootSE)
    #   lm.boot.cort.toxo.adj.mod.f$low.conf <- 
    #     exp(lm.boot.cort.toxo.adj.mod.f$`2.5 %`)
    #   lm.boot.cort.toxo.adj.mod.f$hi.conf <- 
    #     exp(lm.boot.cort.toxo.adj.mod.f$`97.5 %`)
      
    ## j) Bootstrap parameter marginal mean estimates   
      # bootstrapping number of resampling simulations
      # set.seed(512158)
      # mm.boot.cort.toxo.adj.mod.fb <- car::Boot(obj = cort.toxo.adj.mod.f, 
      #                                       f = toxo_mm_coef, R = 2000)
      # 
      # # create a summary statistics data frame
      # mm.summ.boot.cort.toxo.adj.mod.f <- summary(mm.boot.cort.toxo.adj.mod.f)
      # 
      # # estimate permuted confidence intervals 
      # mm.boot.ci.cort.toxo.adj.mod.f <- confint(mm.boot.cort.toxo.adj.mod.f, 
      #                                        type = "perc")
      
    ## k) Undjusted logistic regression: Females
      # c.bin (above vs below median) by toxo_status
      cort.unadj.log.f <- glm(c.bin ~ toxo.status #+ dart.am.pm
                            ,
                        data = subset(plasma_horm_neosp_toxo_data,
                                      sex == 'f' & 
                                        dart.time.diff <= 13 & 
                                        stressca <=2 &
                                        !is.na(x = toxo.status) &
                                        !is.na(x = c.ln)), 
                        family = binomial)
      
      summary(cort.unadj.log.f)    # model summary (log odds scale)
      confint(cort.unadj.log.f)    # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(cort.unadj.log.f), confint (cort.unadj.log.f)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(cort.unadj.log.f), Sigma = vcov(cort.unadj.log.f), 
                estimates = 2) 
      
      
  ### 4.2 Associations between T. gondii status and cortisol levels
    ## a) Unadjusted model: Male 
      # cortisol by T. gondii infection
      cort.toxo.unadj.mod.m <- lm(c.ln
                                  #c.imp.ln # sensitivity analyses imputed data
                                  ~ toxo.status, 
                                     data = subset(plasma_horm_neosp_toxo_data,
                                            sex == 'm' & 
                                            dart.time.diff <= 13 & 
                                            stressca <=2 &
                                            !is.na(x = toxo.status) &
                                            !is.na(x = c.ln)))
      
      summary(cort.toxo.unadj.mod.m) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.m) # 95% CIs (ln scale)
      #plot(cort.toxo.unadj.mod.m) # view fitted vs residuals
      
    ## b) Use emmeans to estimate marginal means
      cort.toxo.unadj.mmean.m <- emmeans(cort.toxo.unadj.mod.m, 
                                            'toxo.status')
      
      summary(cort.toxo.unadj.mmean.m)
      
    ## c) Bootstrap parameter lm estimates   
      # bootstrapping number of resampling simulations
      set.seed(512158)
      lm.boot.cort.toxo.unadj.mod.m <- car::Boot(obj = cort.toxo.unadj.mod.m, 
                                                 f = toxo_reg_coef, R = 2000)
      
      # create a summary statistics data frame
      lm.summ.boot.cort.toxo.unadj.mod.m <- summary(lm.boot.cort.toxo.unadj.mod.m)
      
      # estimate permuted confidence intervals 
      lm.boot.ci.cort.toxo.unadj.mod.m <- confint(lm.boot.cort.toxo.unadj.mod.m, 
                                                  type = "perc")
      
    ## d) Tidy bootstrapped regression estimates table 
      # create a data frame for graphing
      lm.boot.cort.toxo.unadj.mod.m <- cbind(lm.summ.boot.cort.toxo.unadj.mod.m, 
                                             lm.boot.ci.cort.toxo.unadj.mod.m)
      
      # label the estimates in data frame
      lm.boot.cort.toxo.unadj.mod.m$model <- c("Male")
      
      # back transform estimates of interest
      lm.boot.cort.toxo.unadj.mod.m$estimate <- 
        exp(lm.boot.cort.toxo.unadj.mod.m$bootMed)
      lm.boot.cort.toxo.unadj.mod.m$SE <- 
        exp(lm.boot.cort.toxo.unadj.mod.m$bootSE)
      lm.boot.cort.toxo.unadj.mod.m$low.conf <- 
        exp(lm.boot.cort.toxo.unadj.mod.m$`2.5 %`)
      lm.boot.cort.toxo.unadj.mod.m$hi.conf <- 
        exp(lm.boot.cort.toxo.unadj.mod.m$`97.5 %`)
      
    ## e) Bootstrap parameter marginal mean estimates   
      # bootstrapping number of resampling simulations
      # set.seed(512158)
      # mm.boot.cort.toxo.unadj.mod.m <- car::Boot(obj = cort.toxo.unadj.mod.m, 
      #                                       f = toxo_mm_coef, R = 2000)
      # 
      # # create a summary statistics data frame
      # mm.summ.boot.cort.toxo.unadj.mod.m <- summary(mm.boot.cort.toxo.unadj.mod.m)
      # 
      # # estimate permuted confidence intervals 
      # mm.boot.ci.cort.toxo.unadj.mod.m <- confint(mm.boot.cort.toxo.unadj.mod.m, 
      #                                        type = "perc")    
      
    ## f) Adjusted model: Male 
      # cortisol by T. gondii infection
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
      
    ## g) Use emmeans to estimate marginal means
      cort.toxo.adj.mmean.m <- emmeans(cort.toxo.adj.mod.m, 
                                               'toxo.status')
      
      summary(cort.toxo.adj.mmean.m)
      
    # ## h) Bootstrap parameter lm estimates   
    #   # bootstrapping number of resampling simulations
    #   set.seed(512158)
    #   lm.boot.cort.toxo.adj.mod.m <- car::Boot(obj = cort.toxo.adj.mod.m, 
    #                                            f = toxo_reg_coef, R = 2000)
    #   
    #   # create a summary statistics data frame
    #   lm.summ.boot.cort.toxo.adj.mod.m <- summary(lm.boot.cort.toxo.adj.mod.m)
    #   
    #   # estimate permuted confidence intervals 
    #   lm.boot.ci.cort.toxo.adj.mod.m <- confint(lm.boot.cort.toxo.adj.mod.m, 
    #                                             type = "perc")
    #   
    # ## i) Tidy bootstrapped regression estimates table 
    #   # create a data frame for graphing
    #   lm.boot.cort.toxo.adj.mod.m <- cbind(lm.summ.boot.cort.toxo.adj.mod.m, 
    #                                        lm.boot.ci.cort.toxo.adj.mod.m)
    #   
    #   # label the estimates in data frame
    #   lm.boot.cort.toxo.adj.mod.m$model <- c("Cub/Subadult Male")
    #   
    #   # back transform estimates of interest
    #   lm.boot.cort.toxo.adj.mod.m$estimate <- 
    #     exp(lm.boot.cort.toxo.adj.mod.m$bootMed)
    #   lm.boot.cort.toxo.adj.mod.m$SE <- 
    #     exp(lm.boot.cort.toxo.adj.mod.m$bootSE)
    #   lm.boot.cort.toxo.adj.mod.m$low.conf <- 
    #     exp(lm.boot.cort.toxo.adj.mod.m$`2.5 %`)
    #   lm.boot.T.toxo.mod.m.cubsub$hi.conf <- 
    #     exp(lm.boot.cort.toxo.adj.mod.m$`97.5 %`)
    #   
      ## j) Bootstrap parameter marginal mean estimates   
      # bootstrapping number of resampling simulations
      # set.seed(512158)
      # mm.boot.cort.toxo.adj.mod.m <- car::Boot(obj = cort.toxo.adj.mod.m, 
      #                                       f = toxo_mm_coef, R = 2000)
      # 
      # # create a summary statistics data frame
      # mm.summ.boot.cort.toxo.adj.mod.m <- summary(mm.boot.cort.toxo.adj.mod.m)
      # 
      # # estimate permuted confidence intervals 
      # mm.boot.ci.cort.toxo.adj.mod.m <- confint(mm.boot.cort.toxo.adj.mod.m, 
      #                                        type = "perc")
      
      
    ## k) Unadjusted logistic regression: Males
      # c.bin (above vs below median) by toxo_status
      cort.unadj.log.m <- glm(c.bin ~ toxo.status #+ dart.am.pm
                              ,
                            data = subset(plasma_horm_neosp_toxo_data,
                                          sex == 'm' & 
                                            dart.time.diff <= 13 & 
                                            stressca <=2 &
                                            !is.na(x = toxo.status) &
                                            !is.na(x = c.ln)))
      
      summary(cort.unadj.log.m)    # model summary (log odds scale)
      confint(cort.unadj.log.m)    # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(cort.unadj.log.m), confint (cort.unadj.log.m)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(cort.unadj.log.m), Sigma = vcov(cort.unadj.log.m), 
                estimates = 2) 
      
    
      
      
  ### 4.3 Graph of testosterone levels by toxo (regession coefficients)      
    ## a) Combine regression estimates into a tidy table
      cort.toxo.lm.est <- bind_rows(lm.boot.cort.toxo.unadj.mod.f, 
                                    lm.boot.cort.toxo.unadj.mod.m)
      
    ## b) Add a variable, 'term' to distinguish estimates 
      cort.toxo.lm.est$term <- c('intercept', 'infected vs. uninfected', 
                                 'intercept', 'infected vs. uninfected')
      
    ## c) Remove intercepts    
      cort.toxo.lm.est <- cort.toxo.lm.est %>%
        filter(term != 'intercept')
      
      
    ## d) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of 'model' variable 
      cort.toxo.lm.est <- 
        transform(cort.toxo.lm.est, 
                  model = factor(model,
                                 levels = c('Female', 
                                            'Male')))  
      
    ## e) Graph results using dotwhisker, broom, dplyr, and ggplot2 packages
      cort.toxo.lm.est.plot <- 
        ggplot(cort.toxo.lm.est, aes(x = model, y = estimate)) +
        geom_hline(yintercept = 1, color = 'red',
                   linetype = 2) + # line at null behind coefs
        geom_point(size = 6, color = 'red', 
                   position=position_dodge(width = 0.5), alpha =0.8) +
        geom_errorbar(aes(ymin=low.conf, ymax=hi.conf), width=.1,
                      position=position_dodge(.5), color = 'red') +
        #coord_flip() + # flip x and y axes
        labs(title = 'Sex and age stratified associations of T.gondii infection status 
             and plasma cortisol in spotted hyenas',
             subtitle = 'Geometric mean ratio and 95% CI are based on percentile bootstrap (2000 simulations)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(plot.subtitle = element_text(hjust = 0.5, size = 14)) + 
        # bold and size title and axes labels
        theme(text = element_text(size=20, face = 'bold')) +
        theme(legend.position = 'none') +
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) +
        # add major axes
        #theme(axis.line = element_line(colour = 'lightgrey', 
        #                               size = 1, linetype = 'solid')) + 
        # change axes font style, color, size, angle, and margin
        theme(axis.text.x = element_text(face='bold', color='black', 
                                         size=20, angle=0,
                                         margin = margin(t = 10, r = 0, 
                                                         b = 10, l = 0)),
              axis.text.y = element_text(face='bold', color='black', 
                                         size=20, angle=0, 
                                         margin = margin(t = 0, r = 0, 
                                                         b = 0, l = 10))) +
        
        xlab(expression(italic("(Sex stratified models)"))) +
        ylab(expression(atop(bold("Cortisol geometric mean ratio and 95% CI"), 
                             paste(italic("Infected vs. uninfected (reference)"))))) 
      
      
      print(cort.toxo.lm.est.plot)
      
    ## f) Save Plot
      # use ggsave to save the linearization plot
      ggsave('cort.toxo.lm.est.plot.pdf', plot = cort.toxo.lm.est.plot, 
             device = NULL,
             path = paste0(here(),'/output'), 
             scale = 1, width = 12,
             height = 9,
             units = c('in'), dpi = 300, limitsize = TRUE)
      
      
      
      
      
  