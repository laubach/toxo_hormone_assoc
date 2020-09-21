###############################################################################
##############       Associations of Toxoplasma gondii and       ##############
##############       Neospora caninum with hormone levels        ##############
##############              2. Tidy and Join Data                ##############
##############                 By: Zach Laubach                  ##############
##############               created: 19 Sept 2020               ##############
##############             last updated: 19 Sept 2020            ##############
###############################################################################


  ### PURPOSE: 
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Tidy data tables
    # 4: Join tables and re-tidy data 
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
        
      # Check for here and install if not already installed
        if (!'here' %in% installed.packages()[,1]){
          install.packages ('here')
        }
      # load here packages
        library ('here')

        
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


      
###############################################################################
##############                    2. Load RData                  ##############
###############################################################################  
  
  ### 2.1 Load RData
    ## a) Load RData (tables on guppy samples, RNA expression, and beahvior)
      load(paste0(project_data_path,'raw_data_pred_rna_behav.RData'))
     
      
      
###############################################################################
##############              3. Tidy individual tables            ##############
###############################################################################

      
      
      setwd("~/Desktop/NEOSPORA_ANALYSIS/Neospora_Analyses")
      coccid_data <- read.csv("neospora_and_toxo_diagnostics.csv")
      rbind(coccid_data[coccid_data$diagnosis_neo=="high",],coccid_data[coccid_data$diagnosis_neo=="low",])->coccid_data
      as.factor(coccid_data$diagnosis_neo)->coccid_data$diagnosis_neo
      #coccid_data$diagnosis_toxo[coccid_data$diagnosis_toxo=="doubtful",]<-"negative"
      as.factor(coccid_data$diagnosis_toxo)->coccid_data$diagnosis_toxo
      as.factor(coccid_data$sex)->coccid_data$diagnosis_sex
      as.factor(coccid_data$age.cat.dart)->coccid_data$age
      coccid_data$age<-factor(coccid_data$age, levels=c("cub","subadult","adult"))
      coccid_data$neo_status<-factor(coccid_data$neo_status, levels=c("negative","positive"))
      coccid_data$toxo_status<-factor(coccid_data$toxo_status, levels=c("negative","positive"))
      
      as.numeric(coccid_data$age.mon.dart)->coccid_data$agem
      rename(coccid_data, clan = dob.event.data)->coccid_data
      #NOTE: this is clan of birth; some males may have dispersed?
      
      
      
      ##add columns for neospora diagnosis at two cutoffs
      cutofflow<-rep(NA,nrow(coccid_data))
      cutoffhigh<-rep(NA,nrow(coccid_data))
      cbind(coccid_data,cutoffhigh,cutofflow)
      coccid_data$cutofflow <- ifelse(coccid_data$IFA_Neospora<=160, 0, NA)
      coccid_data$cutofflow <- ifelse(coccid_data$IFA_Neospora>160, 1,coccid_data$cutofflow)
      coccid_data$cutoffhigh <- ifelse(coccid_data$IFA_Neospora<=640, 0, NA)
      coccid_data$cutoffhigh <- ifelse(coccid_data$IFA_Neospora>640, 1,coccid_data$cutoffhigh)
      
      
      
      
      
      
      
      
      
      
      
      
  ### 3.1 Tidy _samp_data
    ## a) Add Aripo basin identifer variable to aripo_samp_data
      aripo_samp_data$basin <- paste('aripo')
      
    ## b) Add Quare basin identifer variable to quare_samp_data
      quare_samp_data$basin <- paste('quare') 
      
    ## c) Make a variable group by collapsing 'pop' and 'rear in order to
      # create a 4-level ordinal factor for Aripo data
      aripo_samp_data <- aripo_samp_data  %>%
        mutate(group = case_when(aripo_samp_data$pop %in% c('LP') & 
                                   aripo_samp_data$rear %in% c ('NP')
                                    ~ c('LP_NP'),
                                 aripo_samp_data$pop %in% c('LP') & 
                                   aripo_samp_data$rear %in% c ('P')
                                 ~ c('LP_P'),
                                 aripo_samp_data$pop %in% c('HP') & 
                                   aripo_samp_data$rear %in% c ('NP')
                                 ~ c('HP_NP'),
                                 aripo_samp_data$pop %in% c('HP') & 
                                   aripo_samp_data$rear %in% c ('P')
                                 ~ c('HP_P'))) %>%
        transform(group = factor(group,
                                 levels = c('LP_NP', 'LP_P', 'HP_NP',
                                                     'HP_P')))

    ## d) Make a variable group by collapsing 'pop' and 'rear in order to
      # create a 4-level ordinal factor for Quare data 
      quare_samp_data <- quare_samp_data  %>%
        mutate(group = case_when(quare_samp_data$pop %in% c('LP') & 
                                   quare_samp_data$rear %in% c ('NP')
                                 ~ c('LP_NP'),
                                 quare_samp_data$pop %in% c('LP') & 
                                   quare_samp_data$rear %in% c ('P')
                                 ~ c('LP_P'),
                                 quare_samp_data$pop %in% c('HP') & 
                                   quare_samp_data$rear %in% c ('NP')
                                 ~ c('HP_NP'),
                                 quare_samp_data$pop %in% c('HP') & 
                                   quare_samp_data$rear %in% c ('P')
                                 ~ c('HP_P'))) %>%
        transform(group = factor(group,
                                 levels = c('LP_NP', 'LP_P', 'HP_NP',
                                            'HP_P')))
      
          
  ### 3.2 Tidy _behav_data
    ## a) Rename fish ids by adding 'X' to aripo data
      aripo_behav_data$fish <- sub('^', 'X', aripo_behav_data$fish)
         
    ## b) Select repeat variables to drop from aripo_behave
      aripo_behav_data <- aripo_behav_data %>%
        select(-c(pop, rear, group))

    ## c) Select repeat variables to drop from quare_behave
      quare_behav_data <- quare_behav_data %>%
        select(-c(pop, rear, group))
      
    ## d) Reorder column names to match aripo vs quare
      quare_behav_data <- quare_behav_data[, colnames(aripo_behav_data)]
      
    ## e) Remove rows from aripo_behav_data where 'fish' contains NA
      aripo_behav_data <- aripo_behav_data %>% 
        drop_na('fish')
      
    ## f) Remove rows from quare_behav_data where 'fish' contains NA
      quare_behav_data <- quare_behav_data %>% 
        drop_na('fish')
      
    ## g) Z-score standardize Aripo behavior counts
      aripo_behav_data <- aripo_behav_data %>%
        mutate_at(vars(sigtime:moving), as.numeric) %>%
        mutate_at(vars(sigtime:moving), scale)
      
    ## h) Z-score standardize Quare behavior counts
      quare_behav_data <- quare_behav_data %>%
        mutate_at(vars(sigtime:moving), as.numeric) %>%
        mutate_at(vars(sigtime:moving), scale)

      
  ### 3.3 Tidy _rna_data   
    ## a) View column names of aripo_rna_data
      colnames(aripo_rna_data)
      
    ## b) use gather and spread function to transform data into wide format 
      aripo_rna_data <- aripo_rna_data %>%
        gather(key = 'fish', value = 'rna', 'X001':'X058') %>% 
        pivot_wider(names_from = X1, values_from = rna) 
      
    ## c) View column names of aripo_rna_data
      colnames(quare_rna_data)
      
    ## d) use gather and spread function to transformdata into wide format 
      quare_rna_data <- quare_rna_data %>%
        gather(key = 'fish', value = 'rna', 'M01':'M62') %>% 
        pivot_wider(names_from = X1, values_from = rna) 
      
      
      
###############################################################################
##############          4. Join tables and re-tidy data          ##############
###############################################################################        

  ### 4.1 Make the aripo_data dataframe
    ## a) Left join aripo_behav_data to aripo_samp_data and rename aripo_data
      aripo_data <- aripo_samp_data %>%
        left_join(aripo_behav_data, by = "fish")
      
    ## b) Left join aripo_rna_data to aripo_data
      aripo_data <- aripo_data %>%
        left_join(aripo_rna_data, by = "fish")
      
  
  ### 4.2 Make the quare_data dataframe
      ## a) Left join quare_behav_data to quare_samp_data and rename quare_data
      quare_data <- quare_samp_data %>%
        left_join(quare_behav_data, by = "fish")
      
      ## b) Left join quare_rna_data to aripo_data
      quare_data <- quare_data %>%
        left_join(quare_rna_data, by = "fish")
      
  

###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
  ### 5.1 Export data to an RData file     
    ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 'tidy_data_pred_rna_behav.RData'), 
           list = c('aripo_data', 'quare_data'))
      
      

      
      