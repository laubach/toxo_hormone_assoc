###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############       4 Models: Toxo. and Neosp. associations       #############
#############                  w/fecal hormones                   #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 6 Oct 2020                 #############
#############               last updated: 2 Nov 2020              #############
###############################################################################



  ### PURPOSE: Model associations between T. gondii infection with fecal 
             # testosterone and corticosterone levels in spotted hyenas
  
  
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
##############   3: Infection associations with testosterone     ##############
############################################################################### 
      
      
  ### 3.1 Female associations between T. gondii status and testosterone levels
    ## a) Unadjusted mixed-model: Female cub and subadult 
      # testosterone by T. gondii infection
      T.toxo.unadj.mod.f.cub.sub <- lme4::lmer(testosterone.ng.g.ln ~ 
                                                 toxo.status 
                              + (1|hy.id), 
                              data = subset(fec_horm_neosp_toxo_data_6,
                              #data = subset(fec_horm_toxo_data_restrict,
                              sex == 'f' & 
                              (fecal.age.cat == 'cub' | 
                               fecal.age.cat == 'subadult') &
                                          !is.na(x = toxo.status) & 
                                          !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.unadj.mod.f.cub.sub) # print model summary (ln scale)
      confint(T.toxo.unadj.mod.f.cub.sub) # 95% CIs (ln scale)
      plot(T.toxo.unadj.mod.f.cub.sub) # view fitted vs residuals
      
  ## b) Adjusted mixed-model: Female cub and subadult 
      # testosterone by T. gondii infection
      T.toxo.adj.mod.f.cub.sub <- lme4::lmer(testosterone.ng.g.ln ~ 
                                                 toxo.status + fecal.age.cat
                                               + (1|hy.id), 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                  #data = subset(fec_horm_toxo_data_restrict,
                                  sex == 'f' & 
                                  (fecal.age.cat == 'cub' | 
                                  fecal.age.cat == 'subadult') &
                                  !is.na(x = toxo.status) & 
                                  !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.adj.mod.f.cub.sub) # print model summary (ln scale)
      confint(T.toxo.adj.mod.f.cub.sub) # 95% CIs (ln scale)
      plot(T.toxo.adj.mod.f.cub.sub) # view fitted vs residuals
      
    ## c) Unadjusted mixed-model: Female adult
      # testosterone by T. gondii infection
      T.toxo.unadj.mod.f.adult <- lme4::lmer(testosterone.ng.g.ln ~ 
                                               toxo.status 
                                       + (1|hy.id), 
                            data = subset(fec_horm_neosp_toxo_data_6,
                            #data = subset(fec_horm_toxo_data_restrict,
                            sex == 'f' & 
                            (fecal.age.cat == 'adult') &
                            !is.na(x = toxo.status) & 
                             !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.unadj.mod.f.adult) # print model summary (ln scale)
      confint(T.toxo.unadj.mod.f.adult) # 95% CIs (ln scale)
      plot(T.toxo.unadj.mod.f.adult) # view fitted vs residuals
      
    ## d) Adjusted mixed-model: Female adult 
      # testosterone by T. gondii infection
      T.toxo.adj.mod.f.adult <- lme4::lmer(testosterone.ng.g.ln ~ 
                                               toxo.status + poop.state
                                             + (1|hy.id), 
                                data = subset(fec_horm_neosp_toxo_data_6,
                                #data = subset(fec_horm_toxo_data_restrict,
                                             sex == 'f' & 
                                    (fecal.age.cat == 'adult') &
                                    !is.na(x = toxo.status) & 
                                    !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.adj.mod.f.adult) # print model summary (ln scale)
      confint(T.toxo.adj.mod.f.adult) # 95% CIs (ln scale)
      plot(T.toxo.adj.mod.f.adult) # view fitted vs residuals
      
      
      
      
  ### 3.2 Male associations between T. gondii status and testosterone levels
    ## a) Unadjusted mixed-model: Male cub and subadult 
      # testosterone by T. gondii infection
      T.toxo.unadj.mod.m.cub.sub <- lme4::lmer(testosterone.ng.g.ln ~ 
                                                 toxo.status 
                                               + (1|hy.id), 
                                    data = subset(fec_horm_neosp_toxo_data_6,
                                    #data = subset(fec_horm_toxo_data_restrict,
                                    sex == 'm' & 
                                    (fecal.age.cat == 'cub' | 
                                    fecal.age.cat == 'subadult') &
                                    !is.na(x = toxo.status) & 
                                    !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.unadj.mod.m.cub.sub) # print model summary (ln scale)
      confint(T.toxo.unadj.mod.m.cub.sub) # 95% CIs (ln scale)
      plot(T.toxo.unadj.mod.m.cub.sub) # view fitted vs residuals
      
    ## b) Adjusted mixed-model: Male cub and subadult 
      # testosterone by T. gondii infection
      T.toxo.adj.mod.m.cub.sub <- lme4::lmer(testosterone.ng.g.ln ~ 
                                               toxo.status + migratn.seas.fec
                                             + (1|hy.id), 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                  #data = subset(fec_horm_toxo_data_restrict,
                                  sex == 'm' & 
                                  (fecal.age.cat == 'cub' | 
                                  fecal.age.cat == 'subadult') &
                                  !is.na(x = toxo.status) & 
                                  !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.adj.mod.m.cub.sub) # print model summary (ln scale)
      confint(T.toxo.adj.mod.m.cub.sub) # 95% CIs (ln scale)
      plot(T.toxo.adj.mod.m.cub.sub) # view fitted vs residuals
      
    ## c) Unadjusted mixed-model: Male adult
      # testosterone by T. gondii infection
      T.toxo.unadj.mod.m.adult <- lme4::lmer(testosterone.ng.g.ln ~ 
                                               toxo.status 
                                             + (1|hy.id), 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                  #data = subset(fec_horm_toxo_data_restrict,
                                  sex == 'm' & 
                                  (fecal.age.cat == 'adult') &
                                  !is.na(x = toxo.status) & 
                                  !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.unadj.mod.m.adult) # print model summary (ln scale)
      confint(T.toxo.unadj.mod.m.adult) # 95% CIs (ln scale)
      plot(T.toxo.unadj.mod.m.adult) # view fitted vs residuals
      
    ## d) Adjusted mixed-model: Male adult 
      # testosterone by T. gondii infection
      T.toxo.adj.mod.m.adult <- lme4::lmer(testosterone.ng.g.ln ~ 
                                             toxo.status + migratn.seas.fec
                                           + (1|hy.id), 
                                data = subset(fec_horm_neosp_toxo_data_6,
                                #data = subset(fec_horm_toxo_data_restrict,
                                sex == 'm' & 
                                (fecal.age.cat == 'adult') &
                                !is.na(x = toxo.status) & 
                                !is.na(x = testosterone.ng.g.ln)))
      
      summary(T.toxo.adj.mod.m.adult) # print model summary (ln scale)
      confint(T.toxo.adj.mod.m.adult) # 95% CIs (ln scale)
      plot(T.toxo.adj.mod.m.adult) # view fitted vs residuals
      
      
  ### 3.3 Female associations between T. gondii status and corticosterone levels
    ## a) Unadjusted mixed-model: Female cub and subadult 
      # corticosterone by T. gondii infection
      cort.toxo.unadj.mod.f.cub.sub <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                                 toxo.status 
                                               + (1|hy.id), 
                                    data = subset(fec_horm_neosp_toxo_data_6,
                                    #data = subset(fec_horm_toxo_data_restrict,
                                    sex == 'f' & 
                                    (fecal.age.cat == 'cub' | 
                                    fecal.age.cat == 'subadult') &
                                    !is.na(x = toxo.status) & 
                                    !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.unadj.mod.f.cub.sub) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.f.cub.sub) # 95% CIs (ln scale)
      plot(cort.toxo.unadj.mod.f.cub.sub) # view fitted vs residuals
      
    ## b) Adjusted mixed-model: Female cub and subadult 
      # corticosterone by T. gondii infection
      cort.toxo.adj.mod.f.cub.sub <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                               toxo.status + fecal.age.cat
                                             + (1|hy.id), 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                  #data = subset(fec_horm_toxo_data_restrict,
                                  sex == 'f' & 
                                  (fecal.age.cat == 'cub' | 
                                   fecal.age.cat == 'subadult') &
                                  !is.na(x = toxo.status) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.adj.mod.f.cub.sub) # print model summary (ln scale)
      confint(cort.toxo.adj.mod.f.cub.sub) # 95% CIs (ln scale)
      plot(cort.toxo.adj.mod.f.cub.sub) # view fitted vs residuals
      
    ## c) Unadjusted mixed-model: Female adult
      # corticosterone by T. gondii infection
      cort.toxo.unadj.mod.f.adult <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                               toxo.status 
                                             + (1|hy.id), 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                  #data = subset(fec_horm_toxo_data_restrict,
                                  sex == 'f' & 
                                  (fecal.age.cat == 'adult') &
                                  !is.na(x = toxo.status) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.unadj.mod.f.adult) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.f.adult) # 95% CIs (ln scale)
      plot(cort.toxo.unadj.mod.f.adult) # view fitted vs residuals
      
    ## d) Adjusted mixed-model: Female adult 
      # corticosterone by T. gondii infection
      cort.toxo.adj.mod.f.adult <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                             toxo.status + poop.state
                                           + (1|hy.id), 
                                data = subset(fec_horm_neosp_toxo_data_6,
                                #data = subset(fec_horm_toxo_data_restrict,
                                sex == 'f' & 
                                (fecal.age.cat == 'adult') &
                                !is.na(x = toxo.status) & 
                                !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.adj.mod.f.adult) # print model summary (ln scale)
      confint(cort.toxo.adj.mod.f.adult) # 95% CIs (ln scale)
      plot(cort.toxo.adj.mod.f.adult) # view fitted vs residuals
      
      
  ### 3.4 Male associations between T. gondii status and corticosterone levels
    ## a) Unadjusted mixed-model: Male cub and subadult 
      # corticosterone by T. gondii infection
      cort.toxo.unadj.mod.m.cub.sub <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                                 toxo.status 
                                               + (1|hy.id), 
                                    data = subset(fec_horm_neosp_toxo_data_6,
                                    #data = subset(fec_horm_toxo_data_restrict,
                                    sex == 'm' & 
                                    (fecal.age.cat == 'cub' | 
                                    fecal.age.cat == 'subadult') &
                                    !is.na(x = toxo.status) & 
                                    !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.unadj.mod.m.cub.sub) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.m.cub.sub) # 95% CIs (ln scale)
      plot(cort.toxo.unadj.mod.m.cub.sub) # view fitted vs residuals
      
    ## b) Adjusted mixed-model: Male cub and subadult 
      # corticosterone by T. gondii infection
      cort.toxo.adj.mod.m.cub.sub <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                               toxo.status + migratn.seas.fec
                                               
                                             + (1|hy.id), 
                                  data = subset(fec_horm_neosp_toxo_data_6,
                                  #data = subset(fec_horm_toxo_data_restrict,
                                  sex == 'm' & 
                                  (fecal.age.cat == 'cub' | 
                                  fecal.age.cat == 'subadult') &
                                  !is.na(x = toxo.status) & 
                                  !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.adj.mod.m.cub.sub) # print model summary (ln scale)
      confint(cort.toxo.adj.mod.m.cub.sub) # 95% CIs (ln scale)
      plot(cort.toxo.adj.mod.m.cub.sub) # view fitted vs residuals
      
    ## c) Unadjusted mixed-model: Male adult
      # corticosterone by T. gondii infection
      cort.toxo.unadj.mod.m.adult <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                               toxo.status 
                                             + (1|hy.id), 
                                 data = subset(fec_horm_neosp_toxo_data_6,
                                 #data = subset(fec_horm_toxo_data_restrict,
                                 sex == 'm' & 
                                 (fecal.age.cat == 'adult') &
                                 !is.na(x = toxo.status) & 
                                 !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.unadj.mod.m.adult) # print model summary (ln scale)
      confint(cort.toxo.unadj.mod.m.adult) # 95% CIs (ln scale)
      plot(cort.toxo.unadj.mod.m.adult) # view fitted vs residuals
      
    ## d) Adjusted mixed-model: Male adult 
      # corticosterone by T. gondii infection
      cort.toxo.adj.mod.m.adult <- lme4::lmer(corticosterone.ng.g.ln ~ 
                                             toxo.status + migratn.seas.fec
                                           + (1|hy.id), 
                                data = subset(fec_horm_neosp_toxo_data_6,
                                #data = subset(fec_horm_toxo_data_restrict,
                                sex == 'm' & 
                                (fecal.age.cat == 'adult') &
                                !is.na(x = toxo.status) & 
                                !is.na(x = corticosterone.ng.g.ln)))
      
      summary(cort.toxo.adj.mod.m.adult) # print model summary (ln scale)
      confint(cort.toxo.adj.mod.m.adult) # 95% CIs (ln scale)
      plot(cort.toxo.adj.mod.m.adult) # view fitted vs residuals


     