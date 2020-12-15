###############################################################################
#############          Associations of Toxoplasma gondii          #############
#############                  with hormone levels                #############
#############                                                     #############
#############                   0. Data Import                    #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                created: 19 Sept 2020                #############
#############              last updated: 15 Dec 2020              #############
###############################################################################


  ### PURPOSE: Load Access Fisi backend from MHP R package
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Import data 
    # 3: Export data files
  



###############################################################################
##############             1.  Configure workspace               ##############
###############################################################################

  ### 1.1 Global options
    ## a) clear global environment
      rm(list = ls())

    ## b) prevent R from automatically reading charater strins as factors
      options(stringsAsFactors = FALSE)
      
    # ## c) Setup packrat for reproducibility
    #   library('packrat')
    #   packrat::init('.') #initiate packrat in the current working directory

      
  ### 1.2 Install and load CRAN packages   
    ## a) Data Manipulation and Descriptive Stats Packages
      # load tidyverse packages
        library ('tidyverse')
 
      # load here packages
        library ('here')
      
      
  ### 1.3 Install and load Mara Hyena Project packages 
    ## a) Load the Mara Hyena Project data files from github
    
      # load hyenadata package
      library('hyenadata')

        
  ### 1.4 Get Version and Session Info
    R.Version()
    sessionInfo()
    
    # Developed in:   
    # R version 4.0.2 (2020-06-22)
    # Platform: x86_64-apple-darwin17.0 (64-bit)
    # Running under: macOS Catalina 10.15.7
    
  
  ### 1.5 Set working directory 
    setwd(here())
  
  
  ### 1.6 Set file paths for data importing and exporting
    ## a) The path to sample, normalized RNA expression, and transformed behavior data
      project_data_path <- paste0(here('data/'))
     
  
      
###############################################################################
##############                  2. Import data                   ##############
###############################################################################    
      
  ### 2.1 Import sample data files
    ## a) Import nesop_toxo_data.
      neosp_toxo_data <- read_csv(paste0(project_data_path,
                                  'neosp_toxo_data.csv'))
      
  ### 2.2 Import Access Fisi data files
    ## a) Import Access data backend from Mara Hyena Project data package
      # Use hyenadata package 'load_all_tables' function to load all tables
      hyenadata::load_all_tables()
      # Versiion 1.2.82
      
    ## b) Rename data tables
      # fecal_horm <- tblFecalHormones
      # fecal_repos <- tblFecalRepository
      repro_state <- tblReproStates
      darting <- tblDarting
  
    ## c) Manually load tblHormones2005 csv file
      plasma_horm <- read_csv(paste0(project_data_path,
                                         'tblHormones2005.csv'))
  
  

###############################################################################
##############                3. Export data files               ##############
###############################################################################
      
  ### 3.1 Export data to an RData file     
    ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '1_raw_data_neo_toxo_fec_horm.RData'), 
           list = c('neosp_toxo_data', 
                    #'fecal_horm', 'fecal_repos', 
                    'repro_state', 'plasma_horm', 'darting',
                    'tblClanMembership'))
    