###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
############# 3 Models: Toxo. and Neosp. associations w/ hormones #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 6 Oct 2020                 #############
#############              last updated: 27 Oct 2020              #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii infection with fecal 
             # testosterone and corticosterone levels in spotted hyenas
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Model precision covariates
    # 4: Infection associations with testosterone
    # 5: Infection associations with corticosterone



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
    ## a) load RData: updated fec_horm_neosp_toxo_data joined to hyena 
      # data tables
      load(paste0(project_data_path, '3_neo_toxo_fec_horm.RData'))
     
      
      
###############################################################################
##############           3. Model precision covariates           ##############
###############################################################################
    
### 3.1 Testosterone precision covariates     
  ## a)  Unadjusted mixed-model: testosterone by sex
      T.sex.mod <- lme4::lmer(testosterone.ng.g.ln ~ sex 
                              + (1|hy.id) , 
                     data = subset(fec_horm_neosp_toxo_data_12,
                            !is.na(x = sex))) 
      
      summary(T.sex.mod) # print model summary (ln scale)
      confint(T.sex.mod) # 95% CIs (ln scale)
      plot(T.sex.mod) # view fitted vs residuals
 
    # same mode using nlme instead of lme4;the latter produces p-values  
      T.sex.mod <- nlme::lme(testosterone.ng.g.ln ~ sex, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                           !is.na(x = sex) & 
                                             !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.sex.mod) # print model summary (ln scale)
      intervals(T.sex.mod) # 95% CIs (ln scale)
      plot(T.sex.mod) # view fitted vs residuals
      Anova(T.sex.mod, type = 'II') # type II SS from Car package
      
    # Use emmeans to estimate marginal means
      T.sex.mmean <- emmeans(T.sex.mod, 'sex')
      summary(T.sex.mmean)
      
      
    ## b) Unadjusted mixed-model: testosterone by age 
      T.age.mod <- nlme::lme(testosterone.ng.g.ln ~ fecal.age.cat, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                           !is.na(x = fecal.age.mon) & 
                                             !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.age.mod) # print model summary (ln scale)
      intervals(T.age.mod) # 95% CIs (ln scale)
      plot(T.age.mod) # view fitted vs residuals
      Anova(T.age.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.age.mmean <- emmeans(T.age.mod, 'fecal.age.cat')
      summary(T.age.mmean)
      
    ## c) Unadjusted mixed-model: testosterone by reproductive state 
      T.repro.state.mod <- nlme::lme(testosterone.ng.g.ln ~ state, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                           !is.na(x = state) & 
                                             !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.repro.state.mod) # print model summary (ln scale)
      intervals(T.repro.state.mod) # 95% CIs (ln scale)
      plot(T.repro.state.mod) # view fitted vs residuals
      Anova(T.repro.state.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.repro.state.mmean <- emmeans(T.repro.state.mod, 'state')
      summary(T.repro.state.mmean)
      
      
    ## d) Unadjusted mixed-model: testosterone by fecal sample time of day
      T.poop.am.pm.mod <- nlme::lme(testosterone.ng.g.ln ~ poop.am.pm, 
                                     random = ~ 1|hy.id, 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                            !is.na(x = poop.am.pm) & 
                                            !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.poop.am.pm.mod) # print model summary (ln scale)
      intervals(T.poop.am.pm.mod) # 95% CIs (ln scale)
      plot(T.poop.am.pm.mod) # view fitted vs residuals
      Anova(T.poop.am.pm.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.poop.am.pm.mmean <- emmeans(T.poop.am.pm.mod, 'poop.am.pm')
      summary(T.poop.am.pm.mmean)
      
      
    ## e) Unadjusted mixed-model: testosterone by fecal sample collection season
      T.migratn.seas.fec.mod <- nlme::lme(testosterone.ng.g.ln ~ 
                                            migratn.seas.fec,
                                    random = ~ 1|hy.id, 
                                    data = subset(fec_horm_neosp_toxo_data_12,
                                           !is.na(x = migratn.seas.fec) & 
                                           !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.migratn.seas.fec.mod) # print model summary (ln scale)
      intervals(T.migratn.seas.fec.mod) # 95% CIs (ln scale)
      plot(T.migratn.seas.fec.mod) # view fitted vs residuals
      Anova(T.migratn.seas.fec.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.migratn.seas.mmean <- emmeans(T.migratn.seas.fec.mod, 'migratn.seas.fec')
      summary(T.migratn.seas.mmean)
      
      
    ## f) Unadjusted mixed-model: testosterone by human disturance when fecal
      # sample was collected
      T.hum.pop.poop.mod <- nlme::lme(testosterone.ng.g.ln ~ hum.pop.poop,
                                          random = ~ 1|hy.id, 
                                    data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = hum.pop.poop) & 
                                          !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.hum.pop.poop.mod) # print model summary (ln scale)
      intervals(T.hum.pop.poop.mod) # 95% CIs (ln scale)
      plot(T.hum.pop.poop.mod) # view fitted vs residuals
      Anova(T.hum.pop.poop.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      T.hum.pop.poop.mmean <- emmeans(T.hum.pop.poop.mod, 'hum.pop.poop')
      summary(T.hum.pop.poop.mmean)
      
      
  ### 3.2 Corticosterone precision covariates     
    ## a)  Unadjusted mixed-model: corticosterone by sex
      cort.sex.mod <- nlme::lme(corticosterone.ng.g.ln ~ sex, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = sex) & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.sex.mod) # print model summary (ln scale)
      intervals(cort.sex.mod) # 95% CIs (ln scale)
      plot(cort.sex.mod) # view fitted vs residuals
      Anova(cort.sex.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.sex.mmean <- emmeans(cort.sex.mod, 'sex')
      summary(cort.sex.mmean)
      
    ## b) Unadjusted mixed-model: corticosterone by age 
      cort.age.mod <- nlme::lme(corticosterone.ng.g.ln ~ fecal.age.cat, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = fecal.age.mon) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.age.mod) # print model summary (ln scale)
      intervals(cort.age.mod) # 95% CIs (ln scale)
      plot(cort.age.mod) # view fitted vs residuals
      Anova(cort.age.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.age.mmean <- emmeans(cort.age.mod, 'fecal.age.cat')
      summary(cort.age.mmean)
      
    ## c) Unadjusted mixed-model: corticosterone by reproductive state 
      cort.repro.state.mod <- nlme::lme(corticosterone.ng.g.ln ~ state, 
                                     random = ~ 1|hy.id, 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = state) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.repro.state.mod) # print model summary (ln scale)
      intervals(cort.repro.state.mod) # 95% CIs (ln scale)
      plot(cort.repro.state.mod) # view fitted vs residuals
      Anova(cort.repro.state.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.repro.state.mmean <- emmeans(cort.repro.state.mod, 'state')
      summary(cort.repro.state.mmean)
      
 
    ## d) Unadjusted mixed-model: corticosterone by fecal sample time of day
      cort.poop.am.pm.mod <- nlme::lme(corticosterone.ng.g.ln ~ poop.am.pm, 
                                    random = ~ 1|hy.id, 
                                    data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = poop.am.pm) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.poop.am.pm.mod) # print model summary (ln scale)
      intervals(cort.poop.am.pm.mod) # 95% CIs (ln scale)
      plot(cort.poop.am.pm.mod) # view fitted vs residuals
      Anova(cort.poop.am.pm.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.poop.am.pm.mmean <- emmeans(cort.poop.am.pm.mod, 'poop.am.pm')
      summary(cort.poop.am.pm.mmean)
      
      
    ## e) Unadjusted mixed-model: corticosterone by fecal sample collection season
      cort.migratn.seas.fec.mod <- nlme::lme(corticosterone.ng.g.ln ~ 
                                            migratn.seas.fec,
                                    random = ~ 1|hy.id, 
                                    data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = migratn.seas.fec) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.migratn.seas.fec.mod) # print model summary (ln scale)
      intervals(cort.migratn.seas.fec.mod) # 95% CIs (ln scale)
      plot(cort.migratn.seas.fec.mod) # view fitted vs residuals
      Anova(cort.migratn.seas.fec.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.migratn.seas.fec.mmean <- emmeans(cort.migratn.seas.fec.mod, 
                                       'migratn.seas.fec')
      summary(cort.migratn.seas.fec.mmean)
      
      
    ## f) Unadjusted mixed-model: corticosterone by human disturance when fecal
      # sample was collected
      cort.hum.pop.poop.mod <- nlme::lme(corticosterone.ng.g.ln ~ hum.pop.poop,
                                     random = ~ 1|hy.id, 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = hum.pop.poop) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.hum.pop.poop.mod) # print model summary (ln scale)
      intervals(cort.hum.pop.poop.mod) # 95% CIs (ln scale)
      plot(cort.hum.pop.poop.mod) # view fitted vs residuals
      Anova(cort.hum.pop.poop.mod, type = 'II') # type II SS from Car package
      
      # Use emmeans to estimate marginal means
      cort.hum.pop.poop.mmean <- emmeans(cort.hum.pop.poop.mod, 
                                             'hum.pop.poop')
      summary(cort.hum.pop.poop.mmean)
      
   

###############################################################################
##############   4: Infection associations with testosterone     ##############
############################################################################### 
      
  ### 4.1 Raw associations of hormones by Toxo status
    ## a) Unadjusted mixed-model: testosterone by T. gondii infection
      T.toxo.unadj.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status 
                              + (1|hy.id), 
                              data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = toxo.status) & 
                                          !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.unadj.mod) # print model summary (ln scale)
      confint(T.toxo.unadj.mod) # 95% CIs (ln scale)
      plot(T.toxo.unadj.mod) # view fitted vs residuals
      
      
    ## b) Unadjusted mixed-model: corticosterone by T. gondii infection
      cort.toxo.unadj.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status 
                                     + (1|hy.id), 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = toxo.status) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.unadj.mod) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod) # 95% CIs (ln scale)
      plot(cort.toxo.unadj.mod) # view fitted vs residuals
      
    # ## b) Effect modification mixed-model: 
    #   # testosterone by T. gondii infection * sex
    #   T.toxo.sex.mod <- nlme::lme(corticosterone.ng.g.ln ~ toxo.status*sex,
    #                              random = ~ 1|hy.id, 
    #                              data = subset(fec_horm_neosp_toxo_data_12,
    #                                       !is.na(x = sex) & 
    #                                       !is.na(x = corticosterone.ng.g.ln)))
    #   
    #   summary(T.toxo.sex.mod) # print model summary (ln scale)
    #   intervals(T.toxo.sex.mod) # 95% CIs (ln scale)
    #   plot(T.toxo.sex.mod) # view fitted vs residuals
    #   Anova(T.toxo.sex.mod, type = 'II') # type II SS from Car package
    #   
    #   
    # ## b) Effect modification mixed-model: 
    #   # testosterone by T. gondii infection * fecal.age.cat
    #   T.toxo.age.mod <- nlme::lme(corticosterone.ng.g.ln ~ toxo.status*
    #                                 fecal.age.cat,
    #                                 random = ~ 1|hy.id, 
    #                                 data = subset(fec_horm_neosp_toxo_data_12,
    #                                       !is.na(x = fecal.age.cat) & 
    #                                       !is.na(x = corticosterone.ng.g.ln)))
    #   
    #   summary(T.toxo.age.mod) # print model summary (ln scale)
    #   intervals(T.toxo.age.mod) # 95% CIs (ln scale)
    #   plot(T.toxo.age.mod) # view fitted vs residuals
    #   Anova(T.toxo.age.mod, type = 'II') # type II SS from Car package
      
      
  ### 4.2 Associations between T. gondi and testosterone: Sex stratified
    ## a) Female, sex stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.f.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                   + (1|hy.id) , 
                                   data = subset(fec_horm_neosp_toxo_data_12,
                                            sex == 'f' & 
                                            !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.f.mod) # print model summary (ln scale)
      confint(T.toxo.f.mod) # 95% CIs (ln scale)
      plot(T.toxo.f.mod) # view fitted vs residuals
      
  ## b) Female, sex stratified, adjusted mixed-model: 
      # testosterone by T. gondii infection   
      T.toxo.f.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status + state 
                                 + fecal.age.cat
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'f' & 
                                          !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.f.mod) # print model summary (ln scale)
      confint(T.toxo.f.mod) # 95% CIs (ln scale)
      plot(T.toxo.f.mod) # view fitted vs residuals
      
    ## c) Male, sex stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.m.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'm' & 
                                          !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.m.mod) # print model summary (ln scale)
      confint(T.toxo.m.mod) # 95% CIs (ln scale)
      plot(T.toxo.m.mod) # view fitted vs residuals
      
    ## d) Male, sex stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.m.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                 + fecal.age.cat
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'm' & 
                                        !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.m.mod) # print model summary (ln scale)
      confint(T.toxo.m.mod) # 95% CIs (ln scale)
      plot(T.toxo.m.mod) # view fitted vs residuals
    
      
      
  ### 4.3 Associations between T. gondi and corticosterone: Sex stratified
    ## a) Female, sex stratified, unadjusted mixed-model: 
      # corticosterone by T. gondii infection 
      cort.toxo.f.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'f' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.f.mod) # print model summary (ln scale)
      confint(cort.toxo.f.mod) # 95% CIs (ln scale)
      plot(cort.toxo.f.mod) # view fitted vs residuals
      
    ## b) Female, sex stratified, adjusted mixed-model: 
      # corticosterone by T. gondii infection   
      cort.toxo.f.adj.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status 
                                 + state 
                                 + fecal.age.cat + poop.am.pm 
                                 + migratn.seas.fec + hum.pop.poop
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'f' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.f.adj.mod) # print model summary (ln scale)
      confint(cort.toxo.f.adj.mod) # 95% CIs (ln scale)
      plot(cort.toxo.f.adj.mod) # view fitted vs residuals
      
    ## c) Male, sex stratified, unadjusted mixed-model: 
      # corticosterone by T. gondii infection 
      cort.toxo.m.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'm' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.m.mod) # print model summary (ln scale)
      confint(cort.toxo.m.mod) # 95% CIs (ln scale)
      plot(cort.toxo.m.mod) # view fitted vs residuals
      
    ## d) Male, sex stratified, unadjusted mixed-model: 
      # corticosterone by T. gondii infection 
      cort.toxo.m.adj.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status
                                 + fecal.age.cat + poop.am.pm 
                                 + migratn.seas.fec + hum.pop.poop
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'm' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.m.adj.mod) # print model summary (ln scale)
      confint(cort.toxo.m.adj.mod) # 95% CIs (ln scale)
      plot(cort.toxo.m.adj.mod) # view fitted vs residuals
      
      
  ### 4.4 Associations between T. gondi and testosterone: Age stratified
    ## a) Cub, age stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.c.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               fecal.age.cat == 'cub' & 
                                            !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.c.mod) # print model summary (ln scale)
      confint(T.toxo.c.mod) # 95% CIs (ln scale)
      plot(T.toxo.c.mod) # view fitted vs residuals
      
    ## b) Cub, age stratified,, adjusted mixed-model: 
      # testosterone by T. gondii infection   
      T.toxo.c.adj.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status + sex
                                 # note sex is part of repro state 
                                 # all female cubs nulliporous
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               fecal.age.cat == 'cub' & 
                                            !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.c.adj.mod) # print model summary (ln scale)
      confint(T.toxo.c.adj.mod) # 95% CIs (ln scale)
      plot(T.toxo.c.adj.mod) # view fitted vs residuals
      
    ## c) Subadult, age stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.sa.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               fecal.age.cat == 'subadult' & 
                                            !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.sa.mod) # print model summary (ln scale)
      confint(T.toxo.sa.mod) # 95% CIs (ln scale)
      plot(T.toxo.sa.mod) # view fitted vs residuals
      
    ## d) Subadult, age stratified,, adjusted mixed-model: 
      # testosterone by T. gondii infection   
      T.toxo.sa.adj.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status + sex
                                     # note sex is part of repro state 
                                     # all female subadults nulliporous
                                     + (1|hy.id) , 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                                fecal.age.cat == 'subadult' & 
                                          !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.sa.adj.mod) # print model summary (ln scale)
      confint(T.toxo.sa.adj.mod) # 95% CIs (ln scale)
      plot(T.toxo.sa.adj.mod) # view fitted vs residuals
      
    ## e) Adult, age stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.a.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                  + (1|hy.id) , 
                                  data = subset(fec_horm_neosp_toxo_data_12,
                                                fecal.age.cat == 'adult' & 
                                            !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.a.mod) # print model summary (ln scale)
      confint(T.toxo.a.mod) # 95% CIs (ln scale)
      plot(T.toxo.a.mod) # view fitted vs residuals
      
    ## d) Adult, age stratified,, adjusted mixed-model: 
      # testosterone by T. gondii infection   
      T.toxo.a.adj.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status 
                                     + state
                                      # note sex is part of repro states
                                      + (1|hy.id) , 
                                      data = subset(fec_horm_neosp_toxo_data_12,
                                                fecal.age.cat == 'adult' & 
                                              !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.a.adj.mod) # print model summary (ln scale)
      confint(T.toxo.a.adj.mod) # 95% CIs (ln scale)
      plot(T.toxo.a.adj.mod) # view fitted vs residuals
      
      
      
      
      
  ### 4.5 Associations between T. gondi and corticosterone: Age stratified
    ## a) Cub, age stratified, unadjusted mixed-model: 
      # corticosterone by T. gondii infection 
      cort.toxo.c.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               fecal.age.cat == 'cub' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.c.mod) # print model summary (ln scale)
      confint(cort.toxo.c.mod) # 95% CIs (ln scale)
      plot(cort.toxo.c.mod) # view fitted vs residuals
      
    ## b) Cub, age stratified,, adjusted mixed-model: 
      # corticosterone by T. gondii infection   
      cort.toxo.c.adj.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status 
                                  + sex + poop.am.pm 
                                  + migratn.seas.fec + hum.pop.poop
                                     # note sex is part of repro state 
                                     # all female cubs nulliporous
                                     + (1|hy.id) , 
                                  data = subset(fec_horm_neosp_toxo_data_12,
                                              fecal.age.cat == 'cub' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.c.adj.mod) # print model summary (ln scale)
      confint(cort.toxo.c.adj.mod) # 95% CIs (ln scale)
      plot(cort.toxo.c.adj.mod) # view fitted vs residuals
      
    ## c) Subadult, age stratified, unadjusted mixed-model: 
      # corticosterone by T. gondii infection 
      cort.toxo.sa.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status
                                  + (1|hy.id) , 
                                  data = subset(fec_horm_neosp_toxo_data_12,
                                                fecal.age.cat == 'subadult' & 
                                                  !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.sa.mod) # print model summary (ln scale)
      confint(cort.toxo.sa.mod) # 95% CIs (ln scale)
      plot(cort.toxo.sa.mod) # view fitted vs residuals
      
    ## d) Subadult, age stratified,, adjusted mixed-model: 
      # corticosterone by T. gondii infection   
      cort.toxo.sa.adj.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status 
                                      + sex + poop.am.pm 
                                      + migratn.seas.fec + hum.pop.poop
                                      # note sex is part of repro state 
                                      # all female subadults nulliporous
                                      + (1|hy.id) , 
                                      data = subset(fec_horm_neosp_toxo_data_12,
                                          fecal.age.cat == 'subadult' & 
                                          !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.sa.adj.mod) # print model summary (ln scale)
      confint(cort.toxo.sa.adj.mod) # 95% CIs (ln scale)
      plot(cort.toxo.sa.adj.mod) # view fitted vs residuals
      
    ## e) Adult, age stratified, unadjusted mixed-model: 
      # corticosterone by T. gondii infection 
      cort.toxo.a.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               fecal.age.cat == 'adult' & 
                                                 !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.a.mod) # print model summary (ln scale)
      confint(cort.toxo.a.mod) # 95% CIs (ln scale)
      plot(cort.toxo.a.mod) # view fitted vs residuals
      
    ## d) Adult, age stratified,, adjusted mixed-model: 
      # corticosterone by T. gondii infection   
      cort.toxo.a.adj.mod <- lme4::lmer(corticosterone.ng.g.ln ~ toxo.status 
                                     + state + poop.am.pm 
                                     + migratn.seas.fec + hum.pop.poop
                                     # note sex is part of repro states
                                     + (1|hy.id) , 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                                   fecal.age.cat == 'adult' & 
                                                     !is.na(x = corticosterone.ng.g.ln))) 
      
      summary(cort.toxo.a.adj.mod) # print model summary (ln scale)
      confint(cort.toxo.a.adj.mod) # 95% CIs (ln scale)
      plot(cort.toxo.a.adj.mod) # view fitted vs residuals
      

     