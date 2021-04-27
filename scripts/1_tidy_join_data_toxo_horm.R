###############################################################################
#############          Associations of Toxoplasma gondii          #############
#############             with steroid hormone levels             #############
#############                                                     #############
#############               1. Tidy and Join Data                 #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 19 Sept 2020                #############
#############              last updated: 20 April 2021            #############
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
     
      # load tidyverse package
        library('tidyverse')
      
      # load lubridate package
        library('lubridate') 
     
      # load here package
        library('here')
      
      # load hyenadata package
        library('hyenadata')
      
      # load mice package (for imputation)
        library('mice')
      # prevent mice from masking base cbind and rbind
        cbind <- base::cbind
        rbind <- base::rbind
      
      # load naniar package (graph missing data)
        library(naniar)
      
      

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
  
  ### 3.1  Tidy tblDarting
    ## a) Convert all text to lower case
      darting <- AllCharactersToLower(darting)
      
    ## b) Format variable names (lowercase and separated by '.')
      darting <- FormatVarNames(darting)
      
      
  ### 3.2 Tidy tblReproStates
    ## a) Convert all text to lower case
      repro_state <- AllCharactersToLower(repro_state)
      
    ## b) Format variable names (lowercase and separated by '.')
      repro_state <- FormatVarNames(repro_state)
      
    ## c) Rename mom as hy.id 
      repro_state <- repro_state %>%
        rename('hy.id' = 'mom')
      
      
  ### 3.3 Tidy neosp_toxo_data
    ## a) Convert all text to lower case
      neosp_toxo_data <- AllCharactersToLower(neosp_toxo_data)
      
    ## b) Format as toxo_status as a factor
      neosp_toxo_data <- transform(neosp_toxo_data,
                                   toxo_status = 
                                     factor(toxo_status,
                                            levels = c("negative", 
                                                       "positive")))
      
    ## c) Drop redundant variable, 'toxo.status'   
      neosp_toxo_data <- neosp_toxo_data %>%
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
    
  ### 3.4 Tidy plasma_horm data
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
      
#********************** Data Inclusion/Exclusion Criteria ********************** 
      
    ## f) Remove hy 'oak' from plasma_horm becuase it has two kaycodes listed
      # for a darting that occurs on a single date
      plasma_horm <- plasma_horm %>%
        dplyr::filter(hy.id != 'oak')
      
#********************** Data Inclusion/Exclusion Criteria **********************    

      
      
###############################################################################
##############          4. Impute plasma hormone data            ##############
###############################################################################
  
  ### 4.1 Select and tidy variables to use in imputation
    ## a) Select variables  
      plasma_horm_imp <- plasma_horm %>%
        select(c(hy.id, dart.date, t, p, c, e, a, lh, igf, igf06, nips, 
                 testes, body, shoulder, wt, rank, agemo))
    
    ## b) Replace 0 with NA for T and Cort in order to impute 
      plasma_horm_imp$t <- ifelse(plasma_horm_imp$t == 0, NA, 
                                  plasma_horm_imp$t)
      
      plasma_horm_imp$c <- ifelse(plasma_horm_imp$c == 0, NA, 
                                  plasma_horm_imp$c)
      
    ## c) Extract id for missing T and Cort
      t_id <- plasma_horm_imp %>%
        dplyr::filter(is.na(t)) %>%
        select(hy.id, dart.date)
      
      c_id <- plasma_horm_imp %>%
        dplyr::filter(is.na(c)) %>%
        select(hy.id, dart.date)
      
    ## d) Strip hy.id from plasma_horm_imp
      plasma_horm_imp <- plasma_horm_imp %>%
        select(-c(hy.id, dart.date))
      
    
  ### 4.2 Imputation
    ## a) Determine the amount of missing data
      perc_miss_plot <- gg_miss_var(plasma_horm_imp, show_pct = TRUE)
      print(perc_miss_plot)
      
    ## b) Save histogram plot
      # use ggsave to save the plot
      ggsave('perc_miss_plot.pdf', plot = perc_miss_plot, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE)  
      
      # Graham et al. 2007 recommend 20 imputation for 10-30% missing,
      # and 40 imputations for 50% missing
    
    ## b) Use 'mice' to impute data   
      imputed_hormone_data <- mice(plasma_horm_imp, m=40, 
                           maxit = 50, method = 'pmm', seed = 500)
      
      # m = nos imputed data sets
      # maxit = nos iterations to impute missing data
      # method = 1 of 4 methods for imputation (PMM - numeric vars, 
                # logreg - binar vars, polyreg - factor vars, 
                # proportional odds model - ordered factor vars)
     
    ## c) Extract the imputed t values
      t_imp <- imputed_hormone_data$imp$t
      
    ## d) Combine imputed values with hy.id
      t_imp <- cbind(t_id, t_imp)
      
    ## e) Calculate row averages for each imputed data
      t_imp$t.imp <- rowMeans(t_imp[ , c(3,42)], na.rm=TRUE)
      
    ## f) Extract the imputed c values
      c_imp <- imputed_hormone_data$imp$c
      
    ## g) Combine imputed values with hy.id
      c_imp <- cbind(c_id, c_imp)
      
    ## h) Calculate row averages for each imputed data
      c_imp$c.imp <- rowMeans(c_imp[ , c(3,42)], na.rm=TRUE)
      
  ### 4.3 Join imputed values to back to plasma_horm and tidy data
    ## a) join the average t imputed values to plasma hormone  
      plasma_horm <- plasma_horm %>%
      left_join(select(t_imp, c(hy.id, dart.date, t.imp)),
                by = c('hy.id' = 'hy.id',
                       'dart.date' = 'dart.date'))
      
    ## b) join the average t imputed values to plasma hormone  
      plasma_horm <- plasma_horm %>%
        left_join(select(c_imp, c(hy.id, dart.date, c.imp)),
                  by = c('hy.id' = 'hy.id',
                         'dart.date' = 'dart.date'))
      
    ## c) Fill in NA in t.imp with actual t values 
      plasma_horm$t.imp <- ifelse(is.na(plasma_horm$t.imp), plasma_horm$t,
                                  plasma_horm$t.imp)
      
    ## d) Fill in NA in c.imp with actual c values 
      plasma_horm$c.imp <- ifelse(is.na(plasma_horm$c.imp), plasma_horm$c,
                                  plasma_horm$c.imp)
      
      
      
###############################################################################
##############      4. Join tables and re-tidy: Plasma data      ##############
###############################################################################     
      
  ### 4.1  Combine repro_state w plasma_horm
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


  ### 4.2 Combine plasma_horm with neosp_toxo_data into a dataframe
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

    # ## d) Get the clan status for each hyena on their poop date
    #   clan_status <- hyenadata::get_clan_status(
    #     plasma_horm_neosp_toxo_data$hy.id,
    #     plasma_horm_neosp_toxo_data$dart.date)
    # 
    # ## e) Rename clan as poop.clan
    #   clan_status <- clan_status %>%
    #     rename('dart.clan' = 'clan') %>%
    #     rename('dart.status' = 'status')
    # 
    # ## f) Left join clan_status to plasma_horm_neosp_toxo_data
    #   plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
    #     left_join(select(clan_status, c(ids, dates, dart.clan, dart.status)),
    #               by = c('hy.id' = 'ids',
    #                      'dart.date' = 'dates'))
      
    ## g) Left join darting info to plasma_horm_neosp_toxo_data
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        left_join(select(darting, c(id, collection.date, found.dead, 
                                    darting.time, time.down, 
                                    blood.sampling.time, gnrh.challenge)),
                  by = c('hy.id' = 'id',
                         'dart.date' = 'collection.date'))  


  ### 4.3 Tidy plasma_horm_neosp_toxo_data including precision covariates
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

    # ## i) Create a 2-level ordinal factor indicating human pastoralist presence
    #   # /disturbance based Green et. al 2018 when darting sample was collected
    #   plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
    #     mutate(hum.pop.dart = case_when(plasma_horm_neosp_toxo_data$dart.clan
    #                                     %in% c('talek', 'talek.e', 'kcm',
    #                                            'fig tree') &
    #                               plasma_horm_neosp_toxo_data$dart.yr >= 2000
    #                                     ~ c('hi'),
    #                               plasma_horm_neosp_toxo_data$dart.clan
    #                                     %in% c('talek', 'talek.e', 'kcm',
    #                                            'fig tree') &
    #                               plasma_horm_neosp_toxo_data$dart.yr < 2000
    #                                     ~ c('low'),
    #                               plasma_horm_neosp_toxo_data$dart.clan
    #                                     %in% c('serena.n', 'serena.s',
    #                                            'happy.zebra')
    #                                     ~ c('low')))
    # 
    # ## j) Re-code hum.dist as nominal factor and set level (order)
    #   plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data,
    #                                         hum.pop.dart = factor(hum.pop.dart,
    #                                                   levels = c('hi','low')))
      
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
    
    
  ### 4.4 Tidy plasma_horm_neosp_toxo_data
    ## a) Make list of variables to keep
      var_list <- c('hy.id', 'dart.date', 'tucb', 'a4ucb', 't', 'p', 'c', 'e',
                    'a', 'lh', 'stressca', 'testes', 'plate', 'ifa.neospora', 
                    'diagnosis.neo', 'neo.status','toxo.status', 'spratio', 
                    'diagnosis.toxo', 'hum.pop.dob', 'kay.code',
                    'sample.origin', 'notes.appearance', 'dart.year', 
                    'dob.date', 'dob.event.data', 'sex', 'status', 'mom', 
                    'dad', 'dob.yr', 'age.cat.dart', 'hum.dist.dob', 
                    'rank.dart', 'stan.rank.dart', #'dart.clan', 'dart.status', 
                    'found.dead', 'darting.time', 'time.down', 
                    'blood.sampling.time', 'gnrh.challenge', 'dart.age.days',
                    'dart.mon', 'dart.yr', 'migratn.seas.dart', 'dart.state', 
                    #'hum.pop.dart', 
                    'dart.time.diff','dart.am.pm',
                    't.imp', 'c.imp')
      
    ## b) Select variables according to var_list.
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        select(all_of(var_list))
      
    ## c) Remove two hyenas, baj and gil which have no toxo.status
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(toxo.status))
      

      
###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
  ### 5.1 Export data to an RData file     
    ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '2_tidy_data_neo_toxo_horm.RData'), 
           list = c('plasma_horm_neosp_toxo_data'))

    