###############################################################################
#############            Spotted Hyena Neospora caninum:          #############
############# Determinants and behavior and fitness consequences  #############
#############                                                     #############
############# 4 Models: Toxo. and Neosp. associations w/ hormones #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 6 Oct 2020                #############
#############              last updated: 9 Oct 2020             #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii and Neospora caninum 
             # infection with fecal testosterone and corticosterone levels 
             # in spotted hyenas
  
  
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
      
      # Check for aod and install if not already installed
      if (!'aod' %in% installed.packages()[,1]){
        install.packages ('aod')
      }
      # load aod packages (used to for Wald test)
      library ('aod')
      
      # Check for car and install if not already installed
      if (!'car' %in% installed.packages()[,1]){
        install.packages ('car')
      }
      # load car packages (used for type II and type III SS test)
      library ('car')
      #options(contrasts = c('contr.sum', 'contr.poly'))
      
      # Check for dotwhisker and install if not already installed
      # used with broom to graph beta estimates
      if (!'dotwhisker' %in% installed.packages()[,1]){
        install.packages ('dotwhisker')
      }
      # load dotwhisker packages
      library ('dotwhisker')
      
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
      
      # Check for nlme and install if not already installed
      if (!'nlme' %in% installed.packages()[,1]){
        install.packages ('nlme')
      }
      # load nlme packages
      library ('nlme')
      
      # Check for lme4 and install if not already installed
      if (!'lme4' %in% installed.packages()[,1]){
        install.packages ('lme4')
      }
      # load lme4 packages
      library ('lme4')
      
      # Check for boot and install if not already installed
      # used to generate boot strap CI from LME4
      if (!'boot' %in% installed.packages()[,1]){
        install.packages ('boot')
      }
      # load boot packages
      library ('boot')
      
      # Check for car and install if not already installed
      # includes vif function
      if (!'car' %in% installed.packages()[,1]){
        install.packages ('car')
      }
      # load car packages
      library ('car')
      
      # Check for merTools and install if not already installed
      if (!'merTools' %in% installed.packages()[,1]){
        install.packages ('merTools')
      }
      # load merTools packages
      library ('merTools')

        
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
      
      
      
      
      
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (beta_est = fixef(T.sex.mod), 
                 exp(intervals(T.sex.mod, which = c('fixed'))) )
          
          T.sex.mod.tdy <- tidy(T.sex.mod) %>%
            filter(term != '(Intercept)') %>%
            relabel_predictors(c(toxo.status1 = 'Seropositive'))
          
          
          
          
          
          
      
    ## b) Unadjusted mixed-model: testosterone by age 
      T.age.mod <- nlme::lme(testosterone.ng.g.ln ~ fecal.age.cat, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                           !is.na(x = fecal.age.mon) & 
                                             !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.age.mod) # print model summary (ln scale)
      intervals(T.age.mod) # 95% CIs (ln scale)
      plot(T.age.mod) # view fitted vs residuals
      
    ## c) Unadjusted mixed-model: testosterone by reproductive state 
      T.repro.state.mod <- nlme::lme(testosterone.ng.g.ln ~ state, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                           !is.na(x = state) & 
                                             !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.repro.state.mod) # print model summary (ln scale)
      intervals(T.repro.state.mod) # 95% CIs (ln scale)
      plot(T.repro.state.mod) # view fitted vs residuals
      
      
    ## d) Unadjusted mixed-model: testosterone by fecal sample time of day
      T.poop.am.pm.mod <- nlme::lme(testosterone.ng.g.ln ~ poop.am.pm, 
                                     random = ~ 1|hy.id, 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                            !is.na(x = poop.am.pm) & 
                                            !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.poop.am.pm.mod) # print model summary (ln scale)
      intervals(T.poop.am.pm.mod) # 95% CIs (ln scale)
      plot(T.poop.am.pm.mod) # view fitted vs residuals
      
      
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
      
    ## b) Unadjusted mixed-model: corticosterone by age 
      cort.age.mod <- nlme::lme(corticosterone.ng.g.ln ~ fecal.age.cat, 
                             random = ~ 1|hy.id, 
                             data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = fecal.age.mon) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.age.mod) # print model summary (ln scale)
      intervals(cort.age.mod) # 95% CIs (ln scale)
      plot(cort.age.mod) # view fitted vs residuals
      
    ## c) Unadjusted mixed-model: corticosterone by reproductive state 
      cort.repro.state.mod <- nlme::lme(corticosterone.ng.g.ln ~ state, 
                                     random = ~ 1|hy.id, 
                                     data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = state) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.repro.state.mod) # print model summary (ln scale)
      intervals(cort.repro.state.mod) # 95% CIs (ln scale)
      plot(cort.repro.state.mod) # view fitted vs residuals
      
 
    ## d) Unadjusted mixed-model: corticosterone by fecal sample time of day
      cort.poop.am.pm.mod <- nlme::lme(corticosterone.ng.g.ln ~ poop.am.pm, 
                                    random = ~ 1|hy.id, 
                                    data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = poop.am.pm) & 
                                          !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.poop.am.pm.mod) # print model summary (ln scale)
      intervals(cort.poop.am.pm.mod) # 95% CIs (ln scale)
      plot(cort.poop.am.pm.mod) # view fitted vs residuals
      
      
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
      
   

###############################################################################
##############   4: Infection associations with testosterone     ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate stats neosp status by sex
    ## a) Unadjusted mixed-model: testosterone by T. gondii infection
      T.toxo.unadj.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status 
                              + (1|hy.id), 
                              data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = toxo.status) & 
                                          !is.na(x = testosterone.ng.g.ln)))
      
      
      summary(T.toxo.unadj.mod) # print model summary (ln scale)
      confint(T.toxo.unadj.mod) # 95% CIs (ln scale)
      plot(T.toxo.unadj.mod) # view fitted vs residuals
      
    ## b) Effect modification mixed-model: 
      # testosterone by T. gondii infection * sex
      T.toxo.sex.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                     * sex
                                     + (1|hy.id) , 
                                   data = subset(fec_horm_neosp_toxo_data_12,
                                          !is.na(x = toxo.status) & 
                                          !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.sex.mod) # print model summary (ln scale)
      confint(T.toxo.sex.mod) # 95% CIs (ln scale)
      plot(T.toxo.sex.mod) # view fitted vs residuals
            
      
    ## c) Female, sex stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.f.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                   + (1|hy.id) , 
                                   data = subset(fec_horm_neosp_toxo_data_12,
                                            sex == 'f' & 
                                            !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.f.mod) # print model summary (ln scale)
      confint(T.toxo.f.mod) # 95% CIs (ln scale)
      plot(T.toxo.f.mod) # view fitted vs residuals
      
    ## d) Male, sex stratified, unadjusted mixed-model: 
      # testosterone by T. gondii infection 
      T.toxo.m.mod <- lme4::lmer(testosterone.ng.g.ln ~ toxo.status
                                 + (1|hy.id) , 
                                 data = subset(fec_horm_neosp_toxo_data_12,
                                               sex == 'm' & 
                                          !is.na(x = testosterone.ng.g.ln))) 
      
      summary(T.toxo.m.mod) # print model summary (ln scale)
      confint(T.toxo.m.mod) # 95% CIs (ln scale)
      plot(T.toxo.m.mod) # view fitted vs residuals
    
      