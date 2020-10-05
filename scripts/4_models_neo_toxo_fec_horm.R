###############################################################################
#############            Spotted Hyena Neospora caninum:          #############
############# Determinants and behavior and fitness consequences  #############
#############                                                     #############
#############           5_1 Models: Determinants Neosp.           #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 24 Sept 2020                #############
#############              last updated: 24 Sept 2020             #############
###############################################################################

#**************************  Determinants of Neosp. **************************** 

  ### PURPOSE: Model associations of determinants of Neospora caninum 
             # in spotted hyenas
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Determinants of N. caninum infection models



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
      load(paste0(project_data_path, '3_1_determnts_neosp.RData'))
     
      
      
###############################################################################
##############  3. Determinants of N. caninum infection models   ##############
###############################################################################
    
### 3.1 Neosp. by sex     
  ## a) Unadjusted logistic regression neosp_status by Sex
      sex.log <- glm(neo.status ~ sex , 
                     subset(neosp_data,
                            !is.na(x = sex)),family = binomial) 
      
      summary(sex.log) # print model summary (log odds scale)
      confint(sex.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(sex.log), confint (sex.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(sex.log), Sigma = vcov(sex.log), Terms = 2) 
      
    ## b) Adjusted logistic regression neosp_status by sex
      sex.log.adj <- glm(neo.status ~ sex + age.mon.dart 
                         , 
                         data = neosp_data,
                         #data = neosp_data_no_gil_baj, # sensitivity
                         family = binomial) 
      #****NOTE...control for continuous age to save power and because ****
      summary(sex.log.adj) # print model summary (log odds scale)
      confint(sex.log.adj) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(sex.log.adj), 
                 confint (sex.log.adj)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(sex.log.adj), 
                Sigma = vcov(sex.log.adj), 
                Terms = 2) 
      
      
  ### 3.2 Neosp. by age 
    ## a) Unadjusted logistic regression neosp_status by age
      age.log <- glm(neo.status ~ age.cat.dart , 
                     subset(neosp_data,
                            !is.na(x = age.cat.dart)),family = binomial) 
      
      summary(age.log) # print model summary (log odds scale)
      confint(age.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(age.log), confint (age.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(age.log), Sigma = vcov(age.log), Terms = 2:3) 
      
    ## b) Adjusted logistic regression neosp_status by age.cat.dart
      age.cat.dart.log.adj <- glm(neo.status ~ age.cat.dart +  sex  
                                  , 
                                  data = neosp_data,
                                  #data = neosp_data_no_gil_baj, # sensitivity
                                  family = binomial)  
      
      summary(age.cat.dart.log.adj) # print model summary (log odds scale)
      confint(age.cat.dart.log.adj) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(age.cat.dart.log.adj), 
                 confint (age.cat.dart.log.adj)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(age.cat.dart.log.adj), 
                Sigma = vcov(age.cat.dart.log.adj), 
                Terms = 2) 
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(age.cat.dart.log.adj), 
                Sigma = vcov(age.cat.dart.log.adj), 
                Terms = 3) 
      
      
  ### 3.3 Neosp. by standardized rank 
    ## a) Adult female stratified logistic regression neosp_status by 
      # stanrank.dart
      stanrank.ad.f.log <- glm(neo.status ~ stan.rank.dart,
                          data = neosp_data,
                          subset = (sex == 'f' & age.cat.dart == 'adult'),
                          family = binomial)
      
      
      summary(stanrank.ad.f.log) # print model summary (log odds scale)
      confint(stanrank.ad.f.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(stanrank.ad.f.log),
                 confint (stanrank.ad.f.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(stanrank.ad.f.log),
                Sigma = vcov(stanrank.ad.f.log), Terms = 2)
      
    ## b) Adult male stratified logistic regression neosp_status by 
      # stanrank.dart
      stanrank.ad.m.log <- glm(neo.status ~ stan.rank.dart,
                               data = neosp_data,
                               subset = (sex == 'm' & age.cat.dart == 'adult'),
                               family = binomial)
      
      
      summary(stanrank.ad.m.log) # print model summary (log odds scale)
      confint(stanrank.ad.m.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(stanrank.ad.m.log),
                 confint (stanrank.ad.m.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(stanrank.ad.m.log),
                Sigma = vcov(stanrank.ad.m.log), Terms = 2)
      
      
    ## b) Cub and subadult stratified logistic regression neosp_status by 
      # stanrank.dart (both sexes)
      stanrank.cub.sub.log <- glm(neo.status ~ stan.rank.dart,
                               data = neosp_data,
                               subset = (age.cat.dart == 'cub' | 
                                           age.cat.dart == 'subadult'),
                               family = binomial)
      
      
      summary(stanrank.cub.sub.log) # print model summary (log odds scale)
      confint(stanrank.cub.sub.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(stanrank.cub.sub.log),
                 confint (stanrank.cub.sub.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(stanrank.cub.sub.log),
                Sigma = vcov(stanrank.cub.sub.log), Terms = 2)
      
      
  ### 3.4 Neosp. by human disturbance 
      
#***NOTE: All cubs are from low disturbance, so cubs removed due to their 
      # strong influence which leads to an apparent association btwn 
      # hum.pop.den low and lower infection prevalence; cannot disentagle
      # effect of age from hum.pop.den
      
    ## a) Interaction logistic regression neosp_status by hum_pop_den * sex
      hum.pop.sex.log <- glm(neo.status ~ hum.pop.den * sex, 
                         data = neosp_data,
                         subset = (age.cat.dart == 'adult' | 
                                     age.cat.dart == 'subadult'),
                         family = binomial) 
      
      summary(hum.pop.sex.log) # print model summary (log odds scale)
      confint(hum.pop.sex.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(hum.pop.sex.log), confint (hum.pop.sex.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(hum.pop.sex.log), 
                Sigma = vcov(hum.pop.sex.log), Terms = 4) 
      
      
    ## b) Unadjusted logistic regression neosp_status by hum_pop_den * sex
      hum.pop.log <- glm(neo.status ~ hum.pop.den, 
                             data = neosp_data,
                             subset = (age.cat.dart == 'adult' | 
                                         age.cat.dart == 'subadult'),
                             family = binomial) 
      
      summary(hum.pop.log) # print model summary (log odds scale)
      confint(hum.pop.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(hum.pop.log), confint (hum.pop.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(hum.pop.log), 
                Sigma = vcov(hum.pop.log), Terms = 2) 
      
      
    ## c) Adjusted logistic regression neosp_status by hum_pop_den * sex
      hum.pop.log <- glm(neo.status ~ hum.pop.den + sex + age.mon.dart, 
                         data = neosp_data,
                         subset = (age.cat.dart == 'adult' | 
                                     age.cat.dart == 'subadult'),
                         family = binomial) 
      
      summary(hum.pop.log) # print model summary (log odds scale)
      confint(hum.pop.log) # 95% CIs (log odds scale)
      
      # exponentiate estimates to get onto odds scale
      exp(cbind (O.R. = coef(hum.pop.log), confint (hum.pop.log)))
      
      # Wald Chi-square test of significance using 'aod'
      wald.test(b = coef(hum.pop.log), 
                Sigma = vcov(hum.pop.log), Terms = 2) 
      

###############################################################################
##############               4. Bivariate analysis               ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate stats neosp status by sex

            
      
      

      
      