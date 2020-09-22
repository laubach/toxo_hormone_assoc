###############################################################################
##############       Associations of Toxoplasma gondii and       ##############
##############       Neospora caninum with hormone levels        ##############
##############              2. Tidy and Join Data                ##############
##############                 By: Zach Laubach                  ##############
##############               created: 19 Sept 2020               ##############
##############             last updated: 21 Sept 2020            ##############
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
    ## a) Load RData (diognotistic, hyena lion interactions, 
      # and hyena data base)
      load(paste0(project_data_path,'raw_data_neosp_toxo_hormone.RData'))
     
      
      
###############################################################################
##############              3. Tidy individual tables            ##############
###############################################################################
  
  ### 3.1 Tidy tblFecalHormones
    ## a) Convert all text to lower case
      fecal_horm <- AllCharactersToLower(fecal_horm)
      
    ## b) Format variable names (lowercase and separated by '.')
      fecal_horm <- FormatVarNames(fecal_horm)
      
    ## c) Convert hormone concentrations to numeric
      # make list of variable namges that contain 'ng.g'
      horm_columns <- fecal_horm %>%
        select(contains('ng.g')) %>%
        colnames()
      # convert variables (from horm_columns list) to character first 
      fecal_horm <- fecal_horm %>%
        mutate_at(horm_columns, list(as.character)) 
      
    ## d) Check
      class(fecal_horm$androgens.ng.g)
      
    ## e) Convert variables (from horm_columns list) to numeric 
      fecal_horm <- fecal_horm %>%
        mutate_at(horm_columns, list(as.numeric))
    
    ## f) Check again
      class(fecal_horm$androgens.ng.g)
      
      
  ### 3.2 Tidy fecal_repos
    ## a) Convert all text to lower case
      fecal_repos <- AllCharactersToLower(fecal_repos)
      
    ## b) Format variable names (lowercase and separated by '.')
      fecal_repos <- FormatVarNames(fecal_repos)
      
    ## c) Rename hyena.id as hy.id
      fecal_repos <- fecal_repos %>%
        rename('hy.id' = 'hyena.id')
      
      
  ### 3.3 Tidy tblReproStates
    ## a) Convert all text to lower case
      repro_state <- AllCharactersToLower(repro_state)
      
    ## b) Format variable names (lowercase and separated by '.')
      repro_state <- FormatVarNames(repro_state)
      
    ## c) Rename mom as hy.id 
      repro_state <- repro_state %>%
        rename('hy.id' = 'mom')
      
    
  ### 3.4 Tidy neosp_toxo_data
    ## a) Convert all text to lower case
      neosp_toxo_data <- AllCharactersToLower(neosp_toxo_data)
      
    ## b) Update 'toxo_status'
      neosp_toxo_data <- neosp_toxo_data  %>%
        transform(toxo_status = case_when(!is.na(neosp_toxo_data$toxo_status) & 
                                            toxo_status == 'positive'
                                         ~ 1,
                                         !is.na(neosp_toxo_data$toxo_status) & 
                                           toxo_status == 'negative'
                                         ~ 0))
     
    ## c) Drop redundant variable, 'toxo.status'   
      neosp_toxo_data <- neosp_toxo_data  %>%
        select(-c(toxo.status))
      
    ## d) Update 'neo_status'
      neosp_toxo_data <- neosp_toxo_data  %>%
        transform(neo_status = case_when(!is.na(neosp_toxo_data$neo_status) & 
                                           neo_status == 'positive'
                                          ~ 1,
                                          !is.na(neosp_toxo_data$neo_status) & 
                                           neo_status == 'negative'
                                          ~ 0))
      
    ## e) Format variable names (lowercase and separated by '.')
      neosp_toxo_data <- FormatVarNames(neosp_toxo_data)  
      

      
###############################################################################
##############     4. Join tables and re-tidy fecal luma data    ##############
###############################################################################         
      
  ### 4.1 Make the fecal_luma_data (fecal cort as outcome) dataframe
    ## a) Left join fecal_repos to fecal_horm, retains all columns from both
      fecal_data <- fecal_horm %>%
        left_join(select(fecal_repos, c(fecal.sample.id, hy.id, kaycode,
                                        poop.date, poop.time)),
                  by = "fecal.sample.id")
      
      
  ### 4.2 Combine repro_state w fecal_data
    ## a) Make an empty data frame
      # This is an empty data frame that can store the overlaping 
      # fecal_data and repro_states data
      fecal_repro_data  <- c()   
      
    ## b) For loop to find overlap
      # Iterate over fecal data (hy.id and poop.date), to find interseciton
      # with repro_state data
      for (i in 1:nrow(fecal_data)) { 
        
        # loop through 1:n IDs in fecal_data
        id = paste (fecal_data$hy.id[i])  
        
        # loop through 1:n dates in fecal_data
        poop.date <-(fecal_data$poop.date[i])
        
        # loop through 1:n poop.time in fecal_data
        poop.time <-(fecal_data$poop.time[i])
        
        # create a dataframe to store the rows from repro_states where the
        # id matches hy.id, and poop.date is in between cycle start and stop 
        overlap_poop_repro <- filter(repro_state, id == hy.id & 
                                       poop.date >= cycle.start & 
                                       poop.date <= cycle.stop)
        
        # Control flow
        # if there is no id match and date overlap, 
        # then go to next loop iteration in fecal_data
        if (nrow(overlap_poop_repro) < 1) {
          next
        }
        
        # add the poop.date onto the overlap_poop_repro data
        overlap_poop_repro <- cbind(overlap_poop_repro, poop.date, poop.time)
        
        # add the filtered overlap_poop_repro data to a new dataframe
        # over each iteration of the loop
        fecal_repro_data <- rbind(fecal_repro_data, 
                                  overlap_poop_repro)
      }
      
    ## c) Join repro state data to fecal data 
      fecal_data <- fecal_data %>%
        left_join(select(fecal_repro_data, c(hy.id, poop.date, poop.time,  
                                             state, cycle.start, cycle.stop, 
                                             trimester, parity)),
                  by = c('hy.id' = 'hy.id',
                         'poop.date' = 'poop.date',
                         'poop.time' = 'poop.time')) 
      
    # ## d) Tidy fecal_data after join to repro states
    #   # NOTE: Boom fec samp 5208 ; lyle fec samp id 3617 both
    #   # have overlapping duplicate parities in the same year.
    #   # HACK: Keep the first entry and drop the second
    #   fecal_data <- fecal_data %>%
    #     distinct(fecal.sample.id, .keep_all = T)
      
      
      
    ### 4.3 STOPPED 4.4 Tidy Fecal Data of 3_hy_luma_stress.R...need to go through 4.6
      
      
      
      
      
      
      
      
      ## a) Left join life_hist to fecal_data,
      fecal_data <- fecal_data %>% 
        left_join(life_hist, by = c('hy.id' = 'hy.id'))     
      
      
      
      
      
      
      
      
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
      
      

      
      