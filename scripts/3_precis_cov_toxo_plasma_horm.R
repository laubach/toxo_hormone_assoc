###############################################################################
#############          Associations of Toxoplasma gondii          #############
#############             with steroid hormone levels             #############
#############                                                     #############
#############     3 Models: Precision covariate associations      #############
#############                   plasma hormones                   #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 6 Oct 2020                 #############
#############              last updated: 15 Dec 2020              #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii infection with fecal 
             # testosterone and corticosterone levels in spotted hyenas

  ### NOTE: Cortisol data include only samples collected <= 13 minutes post
          # darting - measures baseline stress. Also only stress state 
          # categories 1 and 2 are included in analyses.

  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Model precision covariates



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
    ## a) load RData: updated 3_neo_toxo_plasma_horm joined to hyena 
      # data tables
      load(paste0(project_data_path, '3_neo_toxo_plasma_horm.RData'))
     
      
      
###############################################################################
##############           3. Model precision covariates           ##############
###############################################################################
    
  ### 3.1 Associations between precision covariates and testosterone levels
    ## a) Unadjusted model: Females
      # testosterone by age 
      T.age.mod.f <- lm(t.ln ~ age.cat.dart, 
                                     data = subset(plasma_horm_neosp_toxo_data,
                                                   sex == 'f' & 
                                                     !is.na(x = t.ln)))
      
      summary(T.age.mod.f) # print model summary (ln scale)
      confint(T.age.mod.f) # 95% CIs (ln scale)
      #plot(T.age.mod.f) # view fitted vs residuals
      Anova(T.age.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.age.mmean.f <- emmeans(T.age.mod.f, 
                                            'age.cat.dart')
      summary(T.age.mmean.f)
  
    ## b) Unadjusted model: Males
      # testosterone by age 
      T.age.mod.m <- lm(t.ln ~ age.cat.dart, 
                        data = subset(plasma_horm_neosp_toxo_data,
                                      sex == 'm' & 
                                        !is.na(x = t.ln)))
      
      summary(T.age.mod.m) # print model summary (ln scale)
      confint(T.age.mod.m) # 95% CIs (ln scale)
      #plot(T.age.mod.m) # view fitted vs residuals
      Anova(T.age.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.age.mmean.m <- emmeans(T.age.mod.m, 
                               'age.cat.dart')
      summary(T.age.mmean.m)
      
    ## c) Unadjusted model: Female adults
      # testosterone by reproductive state 
      T.state.mod.f.adult <- lm(t.ln ~ dart.state, 
                        data = subset(plasma_horm_neosp_toxo_data,
                                      sex == 'f' & age.cat.dart == 'adult'
                                      & !is.na(x = t.ln)))
      
      summary(T.state.mod.f.adult) # print model summary (ln scale)
      confint(T.state.mod.f.adult) # 95% CIs (ln scale)
      #plot(T.state.mod.f.adult) # view fitted vs residuals
      Anova(T.state.mod.f.adult, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.state.mmean.f.adult <- emmeans(T.state.mod.f.adult, 
                               'dart.state')
      summary(T.state.mmean.f.adult)
      
    ## d) Unadjusted model: Male adults
      # testosterone by residency status
      T.status.mod.m.adult <- lm(t.ln ~ status, 
                                data = subset(plasma_horm_neosp_toxo_data,
                                          sex == 'm' & age.cat.dart == 'adult'
                                              & !is.na(x = t.ln)))
      
      summary(T.status.mod.m.adult) # print model summary (ln scale)
      confint(T.status.mod.m.adult) # 95% CIs (ln scale)
      #plot(T.status.mod.m.adult) # view fitted vs residuals
      Anova(T.status.mod.m.adult, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.status.mmean.m.adult <- emmeans(T.status.mod.m.adult, 
                                       'status')
      summary(T.status.mmean.m.adult)
      
    ## e) Unadjusted model: Male adults
      # testosterone by time of day
      #*** NOTE: no female model becaus all samples from am
      T.am.pm.mod.m.adult <- lm(t.ln ~ dart.am.pm, 
                                 data = subset(plasma_horm_neosp_toxo_data,
                                          sex == 'm' & age.cat.dart == 'adult'
                                               & !is.na(x = t.ln)))
      
      summary(T.am.pm.mod.m.adult) # print model summary (ln scale)
      confint(T.am.pm.mod.m.adult) # 95% CIs (ln scale)
      #plot(T.am.pm.mod.m.adult) # view fitted vs residuals
      Anova(T.am.pm.mod.m.adult, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.am.pm.mmean.m.adult <- emmeans(T.am.pm.mod.m.adult, 
                                        'dart.am.pm')
      summary(T.am.pm.mmean.m.adult)
      
    ## f) Unadjusted model: Female adults
      # testosterone by migration season
      #*** NOTE: no female model becaus all samples from am
      T.migrtn.mod.f.adult <- lm(t.ln ~ migratn.seas.dart, 
                                data = subset(plasma_horm_neosp_toxo_data,
                                          sex == 'f' & age.cat.dart == 'adult'
                                              & !is.na(x = t.ln)))
      
      summary(T.migrtn.mod.f.adult) # print model summary (ln scale)
      confint(T.migrtn.mod.f.adult) # 95% CIs (ln scale)
      #plot(T.migrtn.mod.f.adult) # view fitted vs residuals
      Anova(T.migrtn.mod.f.adult, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.migrtn.mmean.f.adult <- emmeans(T.migrtn.mod.f.adult, 
                                       'migratn.seas.dart')
      summary(T.migrtn.mmean.f.adult)
      
    ## g) Unadjusted model: Male adults
      # testosterone by migration season
      #*** NOTE: no female model becaus all samples from am
      T.migrtn.mod.m.adult <- lm(t.ln ~ migratn.seas.dart, 
                                 data = subset(plasma_horm_neosp_toxo_data,
                                        sex == 'm' & age.cat.dart == 'adult'
                                               & !is.na(x = t.ln)))
      
      summary(T.migrtn.mod.m.adult) # print model summary (ln scale)
      confint(T.migrtn.mod.m.adult) # 95% CIs (ln scale)
      #plot(T.migrtn.mod.m.adult) # view fitted vs residuals
      Anova(T.migrtn.mod.m.adult, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.migrtn.mmean.m.adult <- emmeans(T.migrtn.mod.m.adult, 
                                        'migratn.seas.dart')
      summary(T.migrtn.mmean.m.adult)
      
    
  ### 3.2 Associations between precision covariates and corticosterone levels
    ## a) Unadjusted model: Females
      # corticosterone by age 
      cort.age.mod.f <- lm(c.ln ~ age.cat.dart, 
                        data = subset(plasma_horm_neosp_toxo_data,
                                      sex == 'f' & 
                                      dart.time.diff <= 13 & stressca <=2 &
                                        !is.na(x = c.ln)))
      
      summary(cort.age.mod.f) # print model summary (ln scale)
      confint(cort.age.mod.f) # 95% CIs (ln scale)
      #plot(cort.age.mod.f) # view fitted vs residuals
      Anova(cort.age.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.age.mmean.f <- emmeans(cort.age.mod.f, 
                               'age.cat.dart')
      summary(cort.age.mmean.f)
      
    ## b) Unadjusted model: Males
      # corticosterone by age 
      cort.age.mod.m <- lm(c.ln ~ age.cat.dart, 
                        data = subset(plasma_horm_neosp_toxo_data,
                                      sex == 'm' & 
                                        dart.time.diff <= 13 & stressca <=2 &
                                        !is.na(x = c.ln)))
      
      summary(cort.age.mod.m) # print model summary (ln scale)
      confint(cort.age.mod.m) # 95% CIs (ln scale)
      #plot(cort.age.mod.m) # view fitted vs residuals
      Anova(cort.age.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.age.mmean.m <- emmeans(cort.age.mod.m, 
                               'age.cat.dart')
      summary(cort.age.mmean.m)
      
    ## c) Unadjusted model: Female 
      # corticosterone by reproductive state 
      cort.state.mod.f <- lm(c.ln ~ dart.state, 
                                data = subset(plasma_horm_neosp_toxo_data,
                                        sex == 'f' &
                                        dart.time.diff <= 13 & stressca <=2 &
                                        !is.na(x = c.ln)))
      
      summary(cort.state.mod.f) # print model summary (ln scale)
      confint(cort.state.mod.f) # 95% CIs (ln scale)
      #plot(cort.state.mod.f) # view fitted vs residuals
      Anova(cort.state.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.state.mmean.f <- emmeans(cort.state.mod.f, 
                                       'dart.state')
      summary(cort.state.mmean.f)
      
    ## d) Unadjusted model: Male 
      # corticosterone by residency status
      cort.status.mod.m <- lm(c.ln ~ status, 
                                 data = subset(plasma_horm_neosp_toxo_data,
                                        sex == 'm' & 
                                        dart.time.diff <= 13 & stressca <=2 &
                                        !is.na(x = c.ln)))
      
      summary(cort.status.mod.m) # print model summary (ln scale)
      confint(cort.status.mod.m) # 95% CIs (ln scale)
      #plot(cort.status.mod.m) # view fitted vs residuals
      Anova(cort.status.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.status.mmean.m <- emmeans(cort.status.mod.m, 
                                        'status')
      summary(cort.status.mmean.m)
      
    ## e) Unadjusted model: Feale 
      # corticosterone by time of day
      #*** NOTE: no female model becaus all samples from am
      cort.am.pm.mod.f <- lm(c.ln ~ dart.am.pm, 
                             data = subset(plasma_horm_neosp_toxo_data,
                                           sex == 'f' &
                                             dart.time.diff <= 13 & 
                                             stressca <=2 &
                                             !is.na(x = c.ln)))
      
      summary(cort.am.pm.mod.f) # print model summary (ln scale)
      confint(cort.am.pm.mod.f) # 95% CIs (ln scale)
      #plot(cort.am.pm.mod.f) # view fitted vs residuals
      Anova(cort.am.pm.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.am.pm.mmean.f <- emmeans(cort.am.pm.mod.f, 
                                    'dart.am.pm')
      summary(cort.am.pm.mmean.f)
      
    ## f) Unadjusted model: Male 
      # corticosterone by time of day
      #*** NOTE: no female model becaus all samples from am
      cort.am.pm.mod.m <- lm(c.ln ~ dart.am.pm, 
                                data = subset(plasma_horm_neosp_toxo_data,
                                       sex == 'm' &
                                       dart.time.diff <= 13 & stressca <=2 &
                                               !is.na(x = c.ln)))
      
      summary(cort.am.pm.mod.m) # print model summary (ln scale)
      confint(cort.am.pm.mod.m) # 95% CIs (ln scale)
      #plot(cort.am.pm.mod.m) # view fitted vs residuals
      Anova(cort.am.pm.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.am.pm.mmean.m <- emmeans(cort.am.pm.mod.m, 
                                       'dart.am.pm')
      summary(cort.am.pm.mmean.m)
      
    ## g) Unadjusted model: Female 
      # corticosterone by migration season
      #*** NOTE: no female model becaus all samples from am
      cort.migrtn.mod.f <- lm(c.ln ~ migratn.seas.dart, 
                                 data = subset(plasma_horm_neosp_toxo_data,
                                        sex == 'f' & 
                                        dart.time.diff <= 13 & stressca <=2 &
                                                !is.na(x = c.ln)))
      
      summary(cort.migrtn.mod.f) # print model summary (ln scale)
      confint(cort.migrtn.mod.f) # 95% CIs (ln scale)
      #plot(cort.migrtn.mod.f) # view fitted vs residuals
      Anova(cort.migrtn.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.migrtn.mmean.f <- emmeans(cort.migrtn.mod.f, 
                                        'migratn.seas.dart')
      summary(cort.migrtn.mmean.f)
      
    ## h) Unadjusted model: Male 
      # corticosterone by migration season
      #*** NOTE: no female model becaus all samples from am
      cort.migrtn.mod.m <- lm(c.ln ~ migratn.seas.dart, 
                                 data = subset(plasma_horm_neosp_toxo_data,
                                        sex == 'm' & 
                                        dart.time.diff <= 13 & stressca <=2 &
                                                !is.na(x = c.ln)))
      
      summary(cort.migrtn.mod.m) # print model summary (ln scale)
      confint(cort.migrtn.mod.m) # 95% CIs (ln scale)
      #plot(cort.migrtn.mod.m) # view fitted vs residuals
      Anova(cort.migrtn.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.migrtn.mmean.m <- emmeans(cort.migrtn.mod.m, 
                                        'migratn.seas.dart')
      summary(cort.migrtn.mmean.m)
