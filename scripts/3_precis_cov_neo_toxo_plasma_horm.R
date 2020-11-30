###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############     3 Models: Precision covariate associations      #############
#############                  w/fecal hormones                   #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 6 Oct 2020                 #############
#############              last updated: 12 Nov 2020              #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii infection with fecal 
             # testosterone and corticosterone levels in spotted hyenas
  
  
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
    
  ### 3.1 Testosterone precision covariates     
    ## a)  Unadjusted mixed-model: female testosterone by fecal.age.cat
      T.age.mod.f <- nlme::lme(testosterone.ng.g.ln ~ fecal.age.cat, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_6,
                                           sex == 'f' &
                                        !is.na(x = fecal.age.cat) & 
                                        !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.age.mod.f) # print model summary (ln scale)
      intervals(T.age.mod.f) # 95% CIs (ln scale)
      plot(T.age.mod.f) # view fitted vs residuals
      Anova(T.age.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.age.mmean.f <- emmeans(T.age.mod.f, 'fecal.age.cat')
      summary(T.age.mmean.f)
  
    ## b) Unadjusted mixed-model: male testosterone by fecal.age.cat 
      T.age.mod.m <- nlme::lme(testosterone.ng.g.ln ~ fecal.age.cat, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_6,
                                           sex == 'm' &
                                        !is.na(x = fecal.age.cat) & 
                                        !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.age.mod.m) # print model summary (ln scale)
      intervals(T.age.mod.m) # 95% CIs (ln scale)
      plot(T.age.mod.m) # view fitted vs residuals
      Anova(T.age.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.age.mmean.m <- emmeans(T.age.mod.m, 'fecal.age.cat')
      summary(T.age.mmean.m)
      
    ## c) Unadjusted mixed-model: female testosterone by poop.am.pm
      T.poop.am.pm.mod.f <- nlme::lme(testosterone.ng.g.ln ~ poop.am.pm, 
                               random = ~ 1|hy.id, 
                               data = subset(fec_horm_neosp_toxo_data_6,
                                             sex == 'f' &
                                          !is.na(x = poop.am.pm) & 
                                          !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.poop.am.pm.mod.f) # print model summary (ln scale)
      intervals(T.poop.am.pm.mod.f) # 95% CIs (ln scale)
      plot(T.poop.am.pm.mod.f) # view fitted vs residuals
      Anova(T.poop.am.pm.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.poop.am.pm.mmean.f <- emmeans(T.poop.am.pm.mod.f, 'poop.am.pm')
      summary(T.poop.am.pm.mmean.f)
      
    ## d) Unadjusted mixed-model: male testosterone by poop.am.pm 
      T.poop.am.pm.mod.m <- nlme::lme(testosterone.ng.g.ln ~ poop.am.pm, 
                               random = ~ 1|hy.id, 
                               data = subset(fec_horm_neosp_toxo_data_6,
                                             sex == 'm' &
                                        !is.na(x = poop.am.pm) & 
                                        !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.poop.am.pm.mod.m) # print model summary (ln scale)
      intervals(T.poop.am.pm.mod.m) # 95% CIs (ln scale)
      plot(T.poop.am.pm.mod.m) # view fitted vs residuals
      Anova(T.poop.am.pm.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.poop.am.pm.mmean.m <- emmeans(T.poop.am.pm.mod.m, 'poop.am.pm')
      summary(T.poop.am.pm.mmean.m)
      
    ## e) Unadjusted mixed-model: female testosterone by migratn.seas.fec
      T.migratn.seas.fec.mod.f <- nlme::lme(testosterone.ng.g.ln ~ 
                                              migratn.seas.fec, 
                                      random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                    sex == 'f' &
                                              !is.na(x = migratn.seas.fec) & 
                                            !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.migratn.seas.fec.mod.f) # print model summary (ln scale)
      intervals(T.migratn.seas.fec.mod.f) # 95% CIs (ln scale)
      plot(T.migratn.seas.fec.mod.f) # view fitted vs residuals
      Anova(T.migratn.seas.fec.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.migratn.seas.fec.mmean.f <- emmeans(T.migratn.seas.fec.mod.f, 
                                        'migratn.seas.fec')
      summary(T.migratn.seas.fec.mmean.f)
      
    ## f) Unadjusted mixed-model: male testosterone by migratn.seas.fec 
      T.migratn.seas.fec.mod.m <- nlme::lme(testosterone.ng.g.ln ~ 
                                              migratn.seas.fec, 
                                      random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                    sex == 'm' &
                                            !is.na(x = migratn.seas.fec) & 
                                            !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.migratn.seas.fec.mod.m) # print model summary (ln scale)
      intervals(T.migratn.seas.fec.mod.m) # 95% CIs (ln scale)
      plot(T.migratn.seas.fec.mod.m) # view fitted vs residuals
      Anova(T.migratn.seas.fec.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.migratn.seas.fec.mmean.m <- emmeans(T.migratn.seas.fec.mod.m, 
                                            'migratn.seas.fec')
      summary(T.migratn.seas.fec.mmean.m)
     
    ## g) Unadjusted mixed-model: female testosterone by hum.pop.poop
      T.hum.pop.poop.mod.f <- nlme::lme(testosterone.ng.g.ln ~ 
                                              hum.pop.poop, 
                                            random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                          sex == 'f' &
                                      !is.na(x = hum.pop.poop) & 
                                      !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.hum.pop.poop.mod.f) # print model summary (ln scale)
      intervals(T.hum.pop.poop.mod.f) # 95% CIs (ln scale)
      plot(T.hum.pop.poop.mod.f) # view fitted vs residuals
      Anova(T.hum.pop.poop.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.hum.pop.poop.mmean.f <- emmeans(T.hum.pop.poop.mod.f, 
                                            'hum.pop.poop')
      summary(T.hum.pop.poop.mmean.f)
      
    ## h) Unadjusted mixed-model: male testosterone by hum.pop.poop 
      T.hum.pop.poop.mod.m <- nlme::lme(testosterone.ng.g.ln ~ 
                                              hum.pop.poop, 
                                            random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                          sex == 'm' &
                                      !is.na(x = hum.pop.poop) & 
                                      !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.hum.pop.poop.mod.m) # print model summary (ln scale)
      intervals(T.hum.pop.poop.mod.m) # 95% CIs (ln scale)
      plot(T.hum.pop.poop.mod.m) # view fitted vs residuals
      Anova(T.hum.pop.poop.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.hum.pop.poop.mmean.m <- emmeans(T.hum.pop.poop.mod.m, 
                                            'hum.pop.poop')
      summary(T.hum.pop.poop.mmean.m)
      
      
    ## i) Unadjusted mixed-model: female testosterone by poop.state
      T.poop.state.mod.f <- nlme::lme(testosterone.ng.g.ln ~ 
                                          poop.state, 
                                        random = ~ 1|hy.id, 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                                      sex == 'f' &
                                  !is.na(x = poop.state) & 
                                  !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.poop.state.mod.f) # print model summary (ln scale)
      intervals(T.poop.state.mod.f) # 95% CIs (ln scale)
      plot(T.poop.state.mod.f) # view fitted vs residuals
      Anova(T.poop.state.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.poop.state.mmean.f <- emmeans(T.poop.state.mod.f, 
                                        'poop.state')
      summary(T.poop.state.mmean.f)
      
      
  ### 3.2 Corticosterone precision covariates     
    ## a)  Unadjusted mixed-model: female corticosterone by fecal.age.cat
      cort.age.mod.f <- nlme::lme(corticosterone.ng.g.ln ~ fecal.age.cat, 
                               random = ~ 1|hy.id, 
                               data = subset(fec_horm_neosp_toxo_data_6,
                                             sex == 'f' &
                                        !is.na(x = fecal.age.cat) & 
                                        !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.age.mod.f) # print model summary (ln scale)
      intervals(cort.age.mod.f) # 95% CIs (ln scale)
      plot(cort.age.mod.f) # view fitted vs residuals
      Anova(cort.age.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.age.mmean.f <- emmeans(cort.age.mod.f, 'fecal.age.cat')
      summary(cort.age.mmean.f)
      
    ## b) Unadjusted mixed-model: male corticosterone by fecal.age.cat 
      cort.age.mod.m <- nlme::lme(corticosterone.ng.g.ln ~ fecal.age.cat, 
                               random = ~ 1|hy.id, 
                               data = subset(fec_horm_neosp_toxo_data_6,
                                             sex == 'm' &
                                          !is.na(x = fecal.age.cat) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.age.mod.m) # print model summary (ln scale)
      intervals(cort.age.mod.m) # 95% CIs (ln scale)
      plot(cort.age.mod.m) # view fitted vs residuals
      Anova(cort.age.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.age.mmean.m <- emmeans(cort.age.mod.m, 'fecal.age.cat')
      summary(cort.age.mmean.m)
      
    ## c) Unadjusted mixed-model: female corticosterone by poop.am.pm
      cort.poop.am.pm.mod.f <- nlme::lme(corticosterone.ng.g.ln ~ poop.am.pm, 
                                      random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                    sex == 'f' &
                                          !is.na(x = poop.am.pm) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.poop.am.pm.mod.f) # print model summary (ln scale)
      intervals(cort.poop.am.pm.mod.f) # 95% CIs (ln scale)
      plot(cort.poop.am.pm.mod.f) # view fitted vs residuals
      Anova(cort.poop.am.pm.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.poop.am.pm.mmean.f <- emmeans(cort.poop.am.pm.mod.f, 'poop.am.pm')
      summary(cort.poop.am.pm.mmean.f)
      
    ## d) Unadjusted mixed-model: male corticosterone by poop.am.pm 
      cort.poop.am.pm.mod.m <- nlme::lme(corticosterone.ng.g.ln ~ poop.am.pm, 
                                      random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                    sex == 'm' &
                                          !is.na(x = poop.am.pm) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.poop.am.pm.mod.m) # print model summary (ln scale)
      intervals(cort.poop.am.pm.mod.m) # 95% CIs (ln scale)
      plot(cort.poop.am.pm.mod.m) # view fitted vs residuals
      Anova(cort.poop.am.pm.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.poop.am.pm.mmean.m <- emmeans(cort.poop.am.pm.mod.m, 'poop.am.pm')
      summary(cort.poop.am.pm.mmean.m)
      
    ## e) Unadjusted mixed-model: female corticosterone by migratn.seas.fec
      cort.migratn.seas.fec.mod.f <- nlme::lme(corticosterone.ng.g.ln ~ 
                                              migratn.seas.fec, 
                                            random = ~ 1|hy.id, 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                                          sex == 'f' &
                                  !is.na(x = migratn.seas.fec) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.migratn.seas.fec.mod.f) # print model summary (ln scale)
      intervals(cort.migratn.seas.fec.mod.f) # 95% CIs (ln scale)
      plot(cort.migratn.seas.fec.mod.f) # view fitted vs residuals
      Anova(cort.migratn.seas.fec.mod.f, type = 'II') # type II SS from Car 
      
      # Use emmeans to estimate marginal means
      cort.migratn.seas.fec.mmean.f <- emmeans(cort.migratn.seas.fec.mod.f, 
                                            'migratn.seas.fec')
      summary(cort.migratn.seas.fec.mmean.f)
      
    ## f) Unadjusted mixed-model: male corticosterone by migratn.seas.fec 
      cort.migratn.seas.fec.mod.m <- nlme::lme(corticosterone.ng.g.ln ~ 
                                              migratn.seas.fec, 
                                            random = ~ 1|hy.id, 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                                          sex == 'm' &
                                  !is.na(x = migratn.seas.fec) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.migratn.seas.fec.mod.m) # print model summary (ln scale)
      intervals(cort.migratn.seas.fec.mod.m) # 95% CIs (ln scale)
      plot(cort.migratn.seas.fec.mod.m) # view fitted vs residuals
      Anova(cort.migratn.seas.fec.mod.m, type = 'II') # type II SS from Car 
      
      # Use emmeans to estimate marginal means
      cort.migratn.seas.fec.mmean.m <- emmeans(cort.migratn.seas.fec.mod.m, 
                                            'migratn.seas.fec')
      summary(cort.migratn.seas.fec.mmean.m)
      
    ## g) Unadjusted mixed-model: female corticosterone by hum.pop.poop
      cort.hum.pop.poop.mod.f <- nlme::lme(corticosterone.ng.g.ln ~ 
                                          hum.pop.poop, 
                                        random = ~ 1|hy.id, 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                                      sex == 'f' &
                                  !is.na(x = hum.pop.poop) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.hum.pop.poop.mod.f) # print model summary (ln scale)
      intervals(cort.hum.pop.poop.mod.f) # 95% CIs (ln scale)
      plot(cort.hum.pop.poop.mod.f) # view fitted vs residuals
      Anova(cort.hum.pop.poop.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.hum.pop.poop.mmean.f <- emmeans(cort.hum.pop.poop.mod.f, 
                                        'hum.pop.poop')
      summary(cort.hum.pop.poop.mmean.f)
      
    ## h) Unadjusted mixed-model: male corticosterone by hum.pop.poop 
      cort.hum.pop.poop.mod.m <- nlme::lme(corticosterone.ng.g.ln ~ 
                                          hum.pop.poop, 
                                        random = ~ 1|hy.id, 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                                      sex == 'm' &
                                  !is.na(x = hum.pop.poop) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.hum.pop.poop.mod.m) # print model summary (ln scale)
      intervals(cort.hum.pop.poop.mod.m) # 95% CIs (ln scale)
      plot(cort.hum.pop.poop.mod.m) # view fitted vs residuals
      Anova(cort.hum.pop.poop.mod.m, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.hum.pop.poop.mmean.m <- emmeans(cort.hum.pop.poop.mod.m, 
                                        'hum.pop.poop')
      summary(cort.hum.pop.poop.mmean.m)
      
      
    ## i) Unadjusted mixed-model: female corticosterone by poop.state
      cort.poop.state.mod.f <- nlme::lme(corticosterone.ng.g.ln ~ 
                                        poop.state, 
                                      random = ~ 1|hy.id, 
                                      data = subset(fec_horm_neosp_toxo_data_6,
                                                    sex == 'f' &
                                                      !is.na(x = poop.state) & 
                                                      !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.poop.state.mod.f) # print model summary (ln scale)
      intervals(cort.poop.state.mod.f) # 95% CIs (ln scale)
      plot(cort.poop.state.mod.f) # view fitted vs residuals
      Anova(cort.poop.state.mod.f, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.poop.state.mmean.f <- emmeans(cort.poop.state.mod.f, 
                                      'poop.state')
      summary(cort.poop.state.mmean.f)
      
      
      