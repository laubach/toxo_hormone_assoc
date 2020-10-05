###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############               2. Tidy and Join Data                 #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 19 Sept 2020                #############
#############              last updated: 5 Oct 2020               #############
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
    ## a) Load RData (diognotistics and hyena data base)
      load(paste0(project_data_path,'1_raw_data_neo_toxo_fec_horm.RData'))
     
      
      
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
      
    ## f) Format the darting dart.date in neosp_toxo_data    
      neosp_toxo_data <- neosp_toxo_data %>%
        mutate(dart.date = as.Date(neosp_toxo_data$dart.date,
                                   format = '%m/%d/%y'))
      
    ## g) Format the darting dob.date in neosp_toxo_data    
      neosp_toxo_data <- neosp_toxo_data %>%
        mutate(dob.date = as.Date(neosp_toxo_data$dob.date,
                                   format = '%m/%d/%y'))
      
    ## h) Format the darting disappeared.date in neosp_toxo_data    
      neosp_toxo_data <- neosp_toxo_data %>%
        mutate(disappeared.date = as.Date(neosp_toxo_data$disappeared.date,
                                   format = '%m/%d/%y'))
      

      
###############################################################################
##############          4. Join tables and re-tidy data          ##############
###############################################################################         
      
  ### 4.1 Make the fecal_data dataframe
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
      

  ### 4.3 Combine fecal_data with neosp_toxo_data into a long dataframe
    ## a) Left join neosp_toxo_data to fecal_data, retains all columns from both
      fec_horm_neosp_toxo_data <- fecal_data %>%
        left_join(neosp_toxo_data, by = 'hy.id')
      
    ## b) Remove fecal hormone data that does not have corresponding neosp
        # and/or toxo data
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data %>%
        filter(!is.na(toxo.status) | !is.na(neo.status))
      
      
  ### 4.4 Tidy fec_horm_neosp_toxo_data including preciion covariates 
    ## a) Create a varialbes, 'fecal.age.days', and 'fecal.age.mon,' 
      # which indicates how old hyena was on the poop.date using lubridate
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(fecal.age.days = round(interval(dob.date,
                                               poop.date) %/% days(1), 1)) %>%
        mutate(fecal.age.mon = round((interval(dob.date,
                                               poop.date) %/% days(1) 
                                      / 30.44), 1))
      
    ## b) Create a categorical age variable using fecal.age.mon 
      # and based on Holekamp and Smale 1998
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(fecal.age.cat = case_when(sex == "m" & fecal.age.mon <= 12 
                                         ~ c("cub"),
                                         sex == "m" & fecal.age.mon > 12 & 
                                           fecal.age.mon <=24 
                                         ~ c("subadult"),
                                         sex == "m" & fecal.age.mon > 24 
                                         ~ c("adult"),
                                         sex == "f" & fecal.age.mon <= 12 
                                         ~ c("cub"),
                                         sex == "f" & fecal.age.mon > 12 & 
                                           fecal.age.mon <=24 
                                         ~ c("subadult"),
                                         sex == "f" & fecal.age.mon > 24 
                                         ~ c("adult")))
      
    ## c) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) fecal.age.cat variable and sets the 
      # reference level to cub 
      # NOTE: model output is difference in means btwn reference and each 
      # level; factorial contrasts are differences in group level means
      # stored internally as 1, 2, 3 (equal spacing)
      fec_horm_neosp_toxo_data <- transform(fec_horm_neosp_toxo_data, 
                              fecal.age.cat = factor(fecal.age.cat,
                                                     levels = c("cub", 
                                                                "subadult", 
                                                                "adult")))
      
    ## d)  Extract month and year from poop.date
      # Use lubridate to extract the month/yr during which a poop sample was 
      # collected and make a new variable  
      fec_horm_neosp_toxo_data$poop.mon <- 
        month(fec_horm_neosp_toxo_data$poop.date)
      fec_horm_neosp_toxo_data$poop.yr <- 
        year(fec_horm_neosp_toxo_data$poop.date)
      
    ## e) Create a varialbe, 'migratn.seas.fec,' which indicates if a poop 
      # sample was collected in migration (June 1 - Oct 31)
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(migratn.seas.fec = ifelse(poop.mon >= 6 & poop.mon<= 10, 
                                         'migration', 'none'))
      
    ## f) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of migratn.seas.fec variable and sets the reference  
      # level to 'none'
      fec_horm_neosp_toxo_data <- transform(fec_horm_neosp_toxo_data, 
                              migratn.seas.fec = 
                                factor(migratn.seas.fec,
                                       levels = c("none", "migration")))  
      
    ## g) Convert poop.time to a datetime class
      fec_horm_neosp_toxo_data$poop.time <- 
        as.POSIXct(paste(fec_horm_neosp_toxo_data$poop.date,
                         fec_horm_neosp_toxo_data$poop.time), 
                   format = '%Y-%m-%d %H:%M:%S')
      
    ## h) Extract am vs. pm from poop.time
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(poop.am.pm = ifelse(lubridate::am(poop.time), 'am', 
                                   'pm'))
      
    ## i) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of poop.am.pm variable and sets the reference  
      # level to 'am' 
      fec_horm_neosp_toxo_data <- transform(fec_horm_neosp_toxo_data, 
                              poop.am.pm = factor(poop.am.pm,
                                                  levels = c("am", 
                                                             "pm"))) 
      
    ## j) Change state to character
      fec_horm_neosp_toxo_data$state <- 
        as.character(fec_horm_neosp_toxo_data$state)
      
    ## k) Replaces NA with repro state
      # *** NOTE *** Hyena's less than ~750 days old are not in tblReprostates,  
      # but are by default n = nulliparous. Animals older than ~750 days
      # sometimes have missing data on repro state, possibly because
      # cub goes missing - HERE WE MADE DECISION to classify these
      # animals' rerpro state as o = other
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(state = case_when(!is.na(fec_horm_neosp_toxo_data$state)
                                 ~ state,
                                 sex == 'f' & 
                                   is.na(fec_horm_neosp_toxo_data$state) &
                                   fecal.age.days < 750
                                 ~ c('n'),
                                 sex == 'f' & 
                                   is.na(fec_horm_neosp_toxo_data$state) &
                                   fecal.age.days > 750
                                 ~ c('o'),
                                 sex == 'm'
                                 ~ c('m')))
      
    ## l) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of state variable and sets the reference level 
      # to 'n' makes this
      fec_horm_neosp_toxo_data <- transform( fec_horm_neosp_toxo_data, 
                               state = factor(state,
                                              levels = c("n", "p", 
                                                         "l", "o", "m")))
     
      

###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
  ### 5.1 Export data to an RData file     
    ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '2_tidy_data_neo_toxo_fec_horm.RData'), 
           list = c('fec_horm_neosp_toxo_data'))
      
      

      
      