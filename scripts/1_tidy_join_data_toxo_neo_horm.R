###############################################################################
#############        Associations of Toxoplasma gondii and        #############
#############        Neospora caninum with hormone levels         #############
#############                                                     #############
#############               1. Tidy and Join Data                 #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 19 Sept 2020                #############
#############              last updated: 30 Nov 2020              #############
###############################################################################


  ### PURPOSE: Tidy and join data in preparation for downstream analyses
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Tidy data tables
    # 4: Join tables and re-tidy: Fecal data 
    # 5: Join tables and re-tidy: Plasma data 
    # 6: Export data files
  


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
     
      # load tidyverse package
        library('tidyverse')
      
      # load lubridate package
        library('lubridate') 
     
      # load here package
        library('here')
      
      # load lubridate package
        library ('lubridate')
      

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
        select(contains("ng.g")) %>%
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
      
      
  ### 3.3 Tidy tblDarting
    ## a) Convert all text to lower case
      darting <- AllCharactersToLower(darting)
      
    ## b) Format variable names (lowercase and separated by '.')
      darting <- FormatVarNames(darting)
      
    
  ### 3.5 Tidy neosp_toxo_data
    ## a) Convert all text to lower case
      neosp_toxo_data <- AllCharactersToLower(neosp_toxo_data)
      
    ## b) Format as toxo_status as a factor
      neosp_toxo_data <- transform(neosp_toxo_data,
                                   toxo_status = 
                                     factor(toxo_status,
                                            levels = c("negative", 
                                                       "positive")))
     
    ## c) Drop redundant variable, 'toxo.status'   
      neosp_toxo_data <- neosp_toxo_data  %>%
        select(-c(toxo.status))
      
    ## d) Format as neo_status as a factor
      neosp_toxo_data <- transform(neosp_toxo_data,
                                   neo_status = 
                                     factor(neo_status,
                                            levels = c("negative", 
                                                       "positive")))
      
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
      
    ## i) Format sex as a factor
      neosp_toxo_data <- transform(neosp_toxo_data,
                                   sex = 
                                     factor(sex,
                                            levels = c("f", 
                                                       "m")))
    
  ### 3.6 Tidy plasma_horm data
    ## a) Convert all text to lower case
      plasma_horm <- AllCharactersToLower(plasma_horm)
      
    ## b) Format variable names (lowercase and separated by '.')
      plasma_horm <- FormatVarNames(plasma_horm) 
      
    ## c) Rename darting.date as dart.date 
      plasma_horm <- plasma_horm %>%
        rename('dart.date' = 'darting.date')
      
    ## d) Rename id as hy.id
      plasma_horm <- plasma_horm %>%
        rename('hy.id' = 'id')
      
    ## e) Format the darting dart.date in plasma_horm    
      plasma_horm <- plasma_horm %>%
        mutate(dart.date = as.Date(plasma_horm$dart.date,
                                   format = '%m/%d/%y'))
   
 

###############################################################################
##############      4. Join tables and re-tidy: Fecal data       ##############
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
      
    ## d) Rename state as poop.state (repro state when poop sample collected)
      fecal_data <- fecal_data %>%
        rename('poop.state' = 'state')
      

  ### 4.3 Combine fecal_data with neosp_toxo_data into a long dataframe
    ## a) Left join neosp_toxo_data to fecal_data, retains all columns from both
      fec_horm_neosp_toxo_data <- fecal_data %>%
        left_join(neosp_toxo_data, by = 'hy.id')
      
    ## b) Remove fecal hormone data that does not have corresponding neosp
        # and/or toxo data
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data %>%
        filter(!is.na(toxo.status) | !is.na(neo.status))
      
    ## c) Get the clan status for each hyena on their poop date
      clan_status <- hyenadata::get_clan_status(fec_horm_neosp_toxo_data$hy.id,
                                      fec_horm_neosp_toxo_data$poop.date)
    
    ## d) Rename clan as poop.clan
      clan_status <- clan_status %>%
        rename('poop.clan' = 'clan') %>%
        rename('poop.status' = 'status')
      
    ## e) Left join clan_status to fec_horm_neosp_toxo_data
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data %>%
        left_join(select(clan_status, c(ids, dates, poop.clan, poop.status)),
                  by = c('hy.id' = 'ids',
                         'poop.date' = 'dates')) %>%
        distinct(fecal.sample.id, .keep_all = TRUE) # add this because join 
                                                    # causes a duplicate
      
      
  ### 4.4 Tidy fec_horm_neosp_toxo_data including precision covariates 
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
      
    ## c) Fill in missing ages for (Serena animals and immigrant males, 
      # whose dob is unknown)
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(fecal.age.cat = case_when
                 (!is.na(fec_horm_neosp_toxo_data$fecal.age.cat)
                   ~ fecal.age.cat,
                   is.na(fec_horm_neosp_toxo_data$fecal.age.cat) 
                  & sex == "m" & status == 'i'
                  ~ c("adult"),
                 is.na(fec_horm_neosp_toxo_data$fecal.age.cat) 
                  & sex == "f" & (grepl('serena', dob.event.data) |
                                    grepl('happy', dob.event.data))
                  ~ c("adult")))  
      
      
    ## d) Re-code *nominal* factor (with ordered levels)  
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
      
    ## e)  Extract month and year from poop.date
      # Use lubridate to extract the month/yr during which a poop sample was 
      # collected and make a new variable  
      fec_horm_neosp_toxo_data$poop.mon <- 
        month(fec_horm_neosp_toxo_data$poop.date)
      fec_horm_neosp_toxo_data$poop.yr <- 
        year(fec_horm_neosp_toxo_data$poop.date)
      
    ## f) Create a varialbe, 'migratn.seas.fec,' which indicates if a poop 
      # sample was collected in migration (June 1 - Oct 31)
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(migratn.seas.fec = ifelse(poop.mon >= 6 & poop.mon<= 10, 
                                         'migration', 'none'))
      
    ## g) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of migratn.seas.fec variable and sets the reference  
      # level to 'none'
      fec_horm_neosp_toxo_data <- transform(fec_horm_neosp_toxo_data, 
                              migratn.seas.fec = 
                                factor(migratn.seas.fec,
                                       levels = c("none", "migration")))  
      
    ## h) Convert poop.time to a datetime class
      fec_horm_neosp_toxo_data$poop.time <- 
        as.POSIXct(paste(fec_horm_neosp_toxo_data$poop.date,
                         fec_horm_neosp_toxo_data$poop.time), 
                   format = '%Y-%m-%d %H:%M:%S')
      
    ## i) Extract am vs. pm from poop.time
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(poop.am.pm = ifelse(lubridate::am(poop.time), 'am', 
                                   'pm'))
      class(fec_horm_neosp_toxo_data$poop.time)
      
    ## j) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of poop.am.pm variable and sets the reference  
      # level to 'am' 
      fec_horm_neosp_toxo_data <- transform(fec_horm_neosp_toxo_data, 
                              poop.am.pm = factor(poop.am.pm,
                                                  levels = c("am", 
                                                             "pm"))) 
      
    ## k) Change poop.state to character
      fec_horm_neosp_toxo_data$poop.state <- 
        as.character(fec_horm_neosp_toxo_data$poop.state)
      
    ## l) Replaces NA with repro state
      # *** NOTE *** Hyena's less than ~750 days old are not in tblReprostates,  
      # but are by default n = nulliparous. Animals older than ~750 days
      # sometimes have missing data on repro state, possibly because
      # cub goes missing - HERE WE MADE DECISION to classify these
      # animals' rerpro state as o = other
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(poop.state = 
                 case_when(!is.na(fec_horm_neosp_toxo_data$poop.state)
                           ~ poop.state,
                           sex == 'f' & 
                             is.na(fec_horm_neosp_toxo_data$poop.state) &
                             fecal.age.days < 750
                           ~ c('n'),
                           sex == 'f' & 
                             is.na(fec_horm_neosp_toxo_data$poop.state) &
                             fecal.age.days > 750
                           ~ c('o'),
                           sex == 'm'
                           ~ c('m')))
      
    ## m) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of poop.state variable and sets the reference level 
      # to 'n' makes this
      fec_horm_neosp_toxo_data <- transform( fec_horm_neosp_toxo_data, 
                               poop.state = factor(poop.state,
                                              levels = c("n", "p", 
                                                         "l", "o", "m")))
      
    ## n) Rename existing human disturbance variable (which is based on dob)
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data %>%
        rename('hum.pop.dob' = 'hum.pop.den')
      
    ## o) Create a 2-level ordinal factor indicating human pastoralist presence
      # /disturbance based Green et. al 2018 when fecal sample was collected
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data %>%
        mutate(hum.pop.poop = case_when(fec_horm_neosp_toxo_data$poop.clan
                                        %in% c('talek', 'talek.e', 'kcm', 
                                               'fig tree') &
                                    fec_horm_neosp_toxo_data$poop.yr >= 2000
                                        ~ c('hi'),
                                    fec_horm_neosp_toxo_data$poop.clan
                                        %in% c('talek', 'talek.e', 'kcm',
                                               'fig tree') &
                                    fec_horm_neosp_toxo_data$poop.yr < 2000 
                                        ~ c('low'),
                                    fec_horm_neosp_toxo_data$poop.clan 
                                        %in% c('serena.n', 'serena.s',
                                               'happy.zebra')
                                        ~ c('low')))
      
    ## p) Re-code hum.dist as nominal factor and set level (order)
      fec_horm_neosp_toxo_data <- transform(fec_horm_neosp_toxo_data,
                                   hum.pop.poop = factor(hum.pop.poop,
                                                   levels = c('hi','low')))

        
  ### 4.5 Consider data inclusion cut-offs for fec_horm_neosp_toxo_data
    ## a) Calcuate time between diagnosis and fecal samlple colleciton
      fec_horm_neosp_toxo_data <- fec_horm_neosp_toxo_data  %>%
        mutate(mon.btwn.diag.fec = round((interval(dart.date,
                                                   poop.date) %/% days(1) 
                                      / 30.44), 1))
      
      
    ## b) Generate a flag to indicate if the fecal sample collection comes 
      # before or after diagnosis (darting date) +/- 4 months (122 days)
      
      fec_horm_neosp_toxo_data_restrict <- fec_horm_neosp_toxo_data %>%
        mutate(poop.before.after = 
                 case_when(fec_horm_neosp_toxo_data$poop.date >=
                             (fec_horm_neosp_toxo_data$dart.date)
                           ~ c('after'),
                  # fuzzy window setting dart.date 3 months earlier to 
                  # include more fecal samples for more stable average
                           fec_horm_neosp_toxo_data$poop.date <=
                             (fec_horm_neosp_toxo_data$dart.date)
                           ~ c('before')))
      
      
#********************** Data Inclusion/Exclusion Criteria **********************
      
    ## c) Subset the data to include only samples collected within 4 months 
      # of the diagnosis
      fec_horm_neosp_toxo_data_6 <- fec_horm_neosp_toxo_data  %>%
        filter(abs(mon.btwn.diag.fec) <= 6)
      
    ## d) Fecal hormone measues based on toxo diagnosis
      # Subset data to include only fecal hormone data before negative 
      # infection and after postive infection
      fec_horm_toxo_data_restrict <- fec_horm_neosp_toxo_data_restrict %>%
        filter ((grepl('positive', toxo.status) &
                   grepl('after', poop.before.after)) | 
                  (grepl('negative', toxo.status) & 
                     grepl('before', poop.before.after)))
      
    # ## e) Fecal hormone measues based on neosp diagnosis
    #   fec_horm_neosp_data_restrict <- fec_horm_neosp_toxo_data %>%
    #     filter ((grepl('positive', neo.status) &
    #                grepl('after', poop.before.after)) | 
    #               (grepl('negative', neo.status) & 
    #                  grepl('before', poop.before.after)))
      
#********************** Data Inclusion/Exclusion Criteria **********************  
      
 
      
###############################################################################
##############      5. Join tables and re-tidy: Plasma data      ##############
###############################################################################     
      
  ### 5.1  Combine repro_state w plasma_horm
    ## a) Make an empty data frame
      # This is an empty data frame that can store the overlaping
      # plasma_repro_data and repro_states data
      plasma_repro_data  <- c()

      ## b) For loop to find overlap
      # Iterate over plasma data (hy.id and dart.date), to find interseciton
      # with repro_state data
      for (i in 1:nrow(plasma_horm)) {

        # loop through 1:n IDs in plasma_horm
        id <- paste (plasma_horm$hy.id[i])

        # loop through 1:n dates in plasma_horm
        dart.date <- (plasma_horm$dart.date[i])

        # create a dataframe to store the rows from repro_states where the
        # id matches hy.id, and dart.date is in between cycle start and stop
        overlap_plasma_repro <- filter(repro_state, id == hy.id &
                                         dart.date >= cycle.start &
                                         dart.date <= cycle.stop)

        # Control flow
        # if there is no id match and date overlap,
        # then go to next loop iteration in fecal_data
        if (nrow(overlap_plasma_repro) < 1) {
          next
        }

        # add the poop.date onto the overlap_plasma_repro data
        overlap_plasma_repro <- cbind(overlap_plasma_repro, dart.date)

        # add the filtered overlap_plasma_repro data to a new dataframe
        # over each iteration of the loop
        plasma_repro_data <- rbind(plasma_repro_data,
                                   overlap_plasma_repro)
      }

    ## c) Join repro state data to plasma_horm
      plasma_horm <- plasma_horm %>%
        left_join(select(plasma_repro_data, c(hy.id, dart.date,
                                             state, cycle.start, cycle.stop,
                                             trimester, parity)),
                  by = c('hy.id' = 'hy.id',
                         'dart.date' = 'dart.date'))

    ## d) Rename state as dart.state
      plasma_horm <- plasma_horm %>%
        rename('dart.state' = 'state')


  ### 5.2 Combine plasma_horm with neosp_toxo_data into a dataframe
    ## a) Left join neosp_toxo_data to plasma_horm, retains all columns from both
      plasma_horm_neosp_toxo_data <- plasma_horm %>%
        select(-c(sex, kay.code)) %>%
        left_join(neosp_toxo_data, by = 'hy.id')

    ## b) Remove plasma hormone data that does not have corresponding neosp
      # and/or toxo data
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(toxo.status) | !is.na(neo.status)) %>%
        filter(dart.date.x == dart.date.y)

    ## c) Remove extra dart.date.y and rename
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        select(-c(dart.date.y)) %>%
        rename('dart.date' = 'dart.date.x')

    ## d) Get the clan status for each hyena on their poop date
      clan_status <- hyenadata::get_clan_status(
        plasma_horm_neosp_toxo_data$hy.id,
        plasma_horm_neosp_toxo_data$dart.date)

    ## e) Rename clan as poop.clan
      clan_status <- clan_status %>%
        rename('dart.clan' = 'clan') %>%
        rename('dart.status' = 'status')

    ## f) Left join clan_status to plasma_horm_neosp_toxo_data
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        left_join(select(clan_status, c(ids, dates, dart.clan, dart.status)),
                  by = c('hy.id' = 'ids',
                         'dart.date' = 'dates'))
      
    ## g) Left join darting info to plasma_horm_neosp_toxo_data
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        left_join(select(darting, c(id, collection.date, found.dead, 
                                    darting.time, time.down, 
                                    blood.sampling.time, gnrh.challenge)),
                  by = c('hy.id' = 'id',
                         'dart.date' = 'collection.date'))  


  ### 5.3 Tidy plasma_horm_neosp_toxo_data including precision covariates
    ## a) Create a varialbes, 'dart.age.days', and 'dart.age.mon,'
      # which indicates how old hyena was on the dart.date using lubridate
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(dart.age.days = round(interval(dob.date,
                                              dart.date) %/% days(1), 1))

    ## b) Extract month and year from dart.date
      # Use lubridate to extract the month/yr during which a plasma sample was
      # collected and make a new variable
      plasma_horm_neosp_toxo_data$dart.mon <-
        month(plasma_horm_neosp_toxo_data$dart.date)
      plasma_horm_neosp_toxo_data$dart.yr <-
        year(plasma_horm_neosp_toxo_data$dart.date)

    ## c) Create a varialbe, 'migratn.seas.dart,' which indicates if a plasma
      # sample was collected in migration (June 1 - Oct 31)
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(migratn.seas.dart = ifelse(dart.mon >= 6 & dart.mon<= 10,
                                         'migration', 'none'))

    ## d) Re-code *nominal* factor (with ordered levels)
      # Set levels (odering) of migratn.seas.dart variable and sets the 
      # reference level to 'none'
      plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data,
                                               migratn.seas.dart =
                                              factor(migratn.seas.dart,
                                                     levels = c("none",
                                                                "migration")))

    ## e) Change dart.state (repro state when dartede) to character
      plasma_horm_neosp_toxo_data$dart.state <-
        as.character(plasma_horm_neosp_toxo_data$dart.state)

    ## f) Replaces NA with repro state
      # *** NOTE *** Hyena's less than ~750 days old are not in tblReprostates,
      # but are by default n = nulliparous. Animals older than ~750 days
      # sometimes have missing data on repro state, possibly because
      # cub goes missing - HERE WE MADE DECISION to classify these
      # animals' rerpro state as o = other
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(dart.state =
                 case_when(!is.na(plasma_horm_neosp_toxo_data$dart.state)
                           ~ dart.state,
                           sex == 'f' &
                             is.na(plasma_horm_neosp_toxo_data$dart.state) &
                             dart.age.days < 750
                           ~ c('n'),
                           sex == 'f' &
                             is.na(plasma_horm_neosp_toxo_data$dart.state) &
                            dart.age.days > 750
                           ~ c('o'),
                           sex == 'm'
                           ~ c('m')))

    ## g) Re-code *nominal* factor (with ordered levels)
      # Set levels (odering) of state variable and sets the reference level
      # to 'n' makes this
      plasma_horm_neosp_toxo_data <- transform( plasma_horm_neosp_toxo_data,
                                             state = factor(dart.state,
                                                        levels = c("n", "p",
                                                                   "l", "o",
                                                                   "m")))

    ## h) Rename existing human disturbance variable (which is based on dob)
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        rename('hum.pop.dob' = 'hum.pop.den')

    ## i) Create a 2-level ordinal factor indicating human pastoralist presence
      # /disturbance based Green et. al 2018 when darting sample was collected
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        mutate(hum.pop.dart = case_when(plasma_horm_neosp_toxo_data$dart.clan
                                        %in% c('talek', 'talek.e', 'kcm',
                                               'fig tree') &
                                  plasma_horm_neosp_toxo_data$dart.yr >= 2000
                                        ~ c('hi'),
                                  plasma_horm_neosp_toxo_data$dart.clan
                                        %in% c('talek', 'talek.e', 'kcm',
                                               'fig tree') &
                                  plasma_horm_neosp_toxo_data$dart.yr < 2000
                                        ~ c('low'),
                                  plasma_horm_neosp_toxo_data$dart.clan
                                        %in% c('serena.n', 'serena.s',
                                               'happy.zebra')
                                        ~ c('low')))

    ## j) Re-code hum.dist as nominal factor and set level (order)
      plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data,
                                            hum.pop.dart = factor(hum.pop.dart,
                                                      levels = c('hi','low')))
      
    ## k) Re-code *nominal* factor (with ordered levels)
      # Set levels (odering) of state variable and sets the reference level
      # to 'n' makes this
      plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data,
                                                age.cat.dart = 
                                                 factor(age.cat.dart,
                                                        levels = c('cub',
                                                                   'subadult',
                                                                   'adult')))
      
    ## l) Convert darting.time to a datetime class
      plasma_horm_neosp_toxo_data$darting.time <- 
        as.POSIXct(paste(plasma_horm_neosp_toxo_data$dart.date,
                         plasma_horm_neosp_toxo_data$darting.time), 
                   format = '%Y-%m-%d %H:%M:%S')
      
    ## m) Convert time.down to a datetime class
      plasma_horm_neosp_toxo_data$time.down <- 
        as.POSIXct(paste(plasma_horm_neosp_toxo_data$dart.date,
                         plasma_horm_neosp_toxo_data$time.down), 
                   format = '%Y-%m-%d %H:%M:%S')
      
    ## n) Convert blood.sampling.time to a datetime class
      plasma_horm_neosp_toxo_data$blood.sampling.time <- 
        as.POSIXct(paste(plasma_horm_neosp_toxo_data$dart.date,
                         plasma_horm_neosp_toxo_data$blood.sampling.time), 
                   format = '%Y-%m-%d %H:%M:%S')
      
    ## o) calculate the number of mintues between darting and blood draw
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(dart.time.diff = difftime(blood.sampling.time, darting.time,
                                         units = 'mins'))
      
    ## p) Extract am vs. pm from darting.time
      # Use lubridate to extract the month during which a poop sample was 
      # collected and make a new variable
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(dart.am.pm = ifelse(lubridate::am(darting.time), 'am', 
                                   'pm'))
      
    ## q) Re-code *nominal* factor (with ordered levels)  
      # Set levels (odering) of dart.am.pm variable and sets the reference  
      # level to 'am' 
      plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data, 
                                            dart.am.pm = factor(dart.am.pm,
                                                  levels = c('am', 'pm'))) 
    
    
  ### 5.3 Tidy plasma_horm_neosp_toxo_data by removing extra variables
    ## a) Make list of variables to keep
      var_list <- c('hy.id', 'dart.date', 'tucb', 'a4ucb', 't', 'p', 'c', 'e',
                    'a' , 'lh', 'stressca', 'testes', 'plate', 'ifa.neospora', 
                    'diagnosis.neo', 'neo.status','toxo.status', 'spratio', 
                    'diagnosis.toxo', 'hum.pop.dob', 'kay.code',
                    'sample.origin', 'notes.appearance', 'dart.year', 
                    'dob.date', 'dob.event.data', 'sex', 'status', 'mom', 
                    'dad', 'dob.yr', 'age.cat.dart', 'hum.dist.dob', 
                    'rank.dart', 'stan.rank.dart', 'dart.clan', 'dart.status', 
                    'found.dead', 'darting.time', 'time.down', 
                    'blood.sampling.time', 'gnrh.challenge', 'dart.age.days',
                    'dart.mon', 'dart.yr', 'migratn.seas.dart', 'dart.state', 
                    'hum.pop.dart', 'dart.time.diff','dart.am.pm')
      
    ## b) Select variables according to var_list
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        select(all_of(var_list))
   
      
###############################################################################
##############               6. Export data files                ##############
###############################################################################
      
  ### 5.1 Export data to an RData file     
    ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '2_tidy_data_neo_toxo_horm.RData'), 
           list = c('fec_horm_neosp_toxo_data', 'fec_horm_neosp_toxo_data_6',
                    'fec_horm_toxo_data_restrict' 
                    #,'fec_horm_neosp_data_restrict', 
                    ,'plasma_horm_neosp_toxo_data'
                    ))

    