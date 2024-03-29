###############################################################################
#############          Associations of Toxoplasma gondii          #############
#############             with steroid hormone levels             #############
#############                                                     #############
#############       2. Descriptive statsistics plasma data        #############
#############                                                     #############
#############                  By: Zach Laubach                   #############
#############                 created: 5 Oct 2020                 #############
#############             last updated: 22 April 2021             #############
###############################################################################


  ### PURPOSE: Do univariate and bivariate associations for associations 
             # of toxo. and neosp. with plasma horomones in spotted hyenas

  ### NOTE: Cortisol data include only samples collected <= 13 minutes post
            # darting - measures baseline stress. Also only stress state
            #  categories 1 and 2 are included in analyses.
  
  
  # Code Blocks
    # 1: Configure workspace
    # 2: Load RData
    # 3: Univariate analysis 
    # 4: Bivariate analysis



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
    ## a) load RData
      load(paste0(project_data_path, '2_tidy_data_neo_toxo_horm.RData'))
      
      
      
###############################################################################
##############               3. Univariate analysis              ##############
###############################################################################
    
  ### 3.1 Univariate stats plasma testosterone (T)
    ## a) Descriptive stats plasma testosterone (T) data
      univar_plasma_T_stat <- plasma_horm_neosp_toxo_data %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                              na.rm = T),2),
                   stdev.T = round (sd(t, 
                                              na.rm = T), 2),
                   med.T = round(median(t,
                                        na.rm = T), 2),
                   min.T = round(min(t,
                                         na.rm = T), 2),
                   max.T = round(max(t,
                                     na.rm = T), 2))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_plasma_T_stat.pdf'), height = 4, width = 8)
      grid.table(univar_plasma_T_stat)
      dev.off()
      
      
    ## c) Sample size (includes sample below assay detection limit)
      sex_age_plasma_T_tot_samp_size <- plasma_horm_neosp_toxo_data %>%
        group_by(sex, age.cat.dart) %>%
        summarise (n.T = sum(!is.na(t))) %>%
        ungroup()
      
    ## d) Descriptive stats plasma testosterone (T) data within sex and age
      sex_age_sum_plasma_T_stat <- plasma_horm_neosp_toxo_data %>%
        filter(t > 0) %>%
        group_by(sex, age.cat.dart) %>%
        summarise (n.T = sum(!is.na(t)),
                   avg.T = round (mean(t, 
                                       na.rm = T),2),
                   stdev.T = round (sd(t, 
                                       na.rm = T), 2),
                   med.T = round(median(t,
                                        na.rm = T), 2),
                   min.T = round(min(t,
                                     na.rm = T), 2),
                   max.T = round(max(t,
                                     na.rm = T), 2)) %>%
        ungroup()
      
      
    ## e) Histogram plasma testosterone
      hist_plot_plasma_T <- ggplot(data=plasma_horm_neosp_toxo_data , 
                                aes(x=t)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 6, by = 0.25), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0,6)) +
        labs(title = 'Histogram of Plasma Testosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Fec. Testosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_plasma_T)
      
      
#********************** Data Inclusion/Exclusion Criteria ********************** 
      
    ## f) Set Zero's to 1/2 the minimum detected value
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(t.adj = ifelse(t == 0, 
           (0.5*min(plasma_horm_neosp_toxo_data[,'t']
                    [which(plasma_horm_neosp_toxo_data[,'t']> 0)])),
            plasma_horm_neosp_toxo_data$t))
      
#********************** Data Inclusion/Exclusion Criteria **********************      

      
    ## g) Natural log transformation of T data and T.imp
      plasma_horm_neosp_toxo_data $t.ln <- 
        log(plasma_horm_neosp_toxo_data$t.adj)
      
      plasma_horm_neosp_toxo_data $t.imp.ln <- 
        log(plasma_horm_neosp_toxo_data$t.imp)
      
    ## h) Histogram of Nat. Log. of plasma testosterone
      hist_plot_ln_plasma_T <- ggplot(data=plasma_horm_neosp_toxo_data , 
                            aes(x=t.ln)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(-5, 3, by = 0.5), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(-5, 3)) +
        labs(title = 'Histogram of Nat. Log. Plasma Testosterone 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Nat. Log. Fec. Testosterone (ng/g)', y='Frequency') 
      
      print(hist_plot_ln_plasma_T)

    ## i) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_T.pdf', plot = hist_plot_ln_plasma_T, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE)  

      
#************************** Data Manipulation Criteria ************************* 
      
    ## j) Descriptive stats plasma testosterone (T) to determine cutoffs
      sex_age_plasma_T_cutoffs <- plasma_horm_neosp_toxo_data %>%
        group_by(sex, age.cat.dart) %>%
                  summarise (n.T = sum(!is.na(t)),
                  avg.T = round (mean(t, 
                                       na.rm = T),2),
                  stdev.T = round (sd(t, 
                                       na.rm = T), 2),
                  med.T = round(median(t,
                                        na.rm = T), 2),
                  min.T = round(min(t,
                                     na.rm = T), 2),
                  max.T = round(max(t,
                                     na.rm = T), 2)) %>%
        ungroup()
      
    ## k) Within sex and age strata, categorize t-levels as low vs hi
      #*** BASED on within age/sex median values***
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%     
        mutate(t.bin = case_when(plasma_horm_neosp_toxo_data$sex == 'f' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart == 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t <= 0.03
                                 ~ 0,
                                 plasma_horm_neosp_toxo_data$sex == 'f' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart == 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t > 0.03
                                 ~ 1,
                                 plasma_horm_neosp_toxo_data$sex == 'f' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart != 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t <= 0
                                 ~ 0,
                                 plasma_horm_neosp_toxo_data$sex == 'f' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart != 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t > 0
                                 ~ 1,
                                 plasma_horm_neosp_toxo_data$sex == 'm' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart == 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t <= 0.24
                                 ~ 0,
                                 plasma_horm_neosp_toxo_data$sex == 'm' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart == 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t > 0.24
                                 ~ 1,
                                 plasma_horm_neosp_toxo_data$sex == 'm' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart != 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t <= 0
                                 ~ 0,
                                 plasma_horm_neosp_toxo_data$sex == 'm' &
                                   plasma_horm_neosp_toxo_data$age.cat.dart != 
                                   'adult' & 
                                   plasma_horm_neosp_toxo_data$t > 0
                                 ~ 1)) 
      
    # ## l) Re-code t.bin as nominal factor and set level (order)
    #   plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data,
    #                                            t.bin = factor(t.bin,
    #                                            levels = c('hi','low')))
      
#************************** Data Manipulation Criteria *************************       
      
    
  ### 3.2 Univariate stats plasma Cortisol (cort)
    ## a) Descriptive stats plasma Cortisol (cort) data
      univar_plasma_cort_stat <- plasma_horm_neosp_toxo_data  %>%
        filter(stressca <=2 & dart.time.diff <= 13 & 
                 !is.na(x = toxo.status)) %>%
        summarise (n.cort = sum(!is.na(c )),
                   avg.cort = round (mean(c , 
                                       na.rm = T),2),
                   stdev.cort = round (sd(c , 
                                       na.rm = T), 2),
                   med.cort = round(median(c ,
                                        na.rm = T), 2),
                   min.cort = round(min(c ,
                                     na.rm = T), 2),
                   max.cort = round(max(c ,
                                     na.rm = T), 2)) 
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_plasma_cort_stat.pdf'), height = 4, width = 8)
      grid.table(univar_plasma_cort_stat)
      dev.off()
      
    ## c) Sample size (includes sample below assay detection limit)
      sex_age_plasma_cort_tot_samp_size <- plasma_horm_neosp_toxo_data %>%
        filter(stressca <=2 & dart.time.diff <= 13 & 
                 !is.na(x = toxo.status)) %>%
        group_by(sex, age.cat.dart) %>%
        summarise (n.C = sum(!is.na(c))) %>%
        ungroup()
      
    ## d) Descriptive stats plasma corticosterone (C) data within sex and age
      sex_age_sum_plasma_C_stat <- plasma_horm_neosp_toxo_data %>%
        filter(c > 0 & stressca <=2 & dart.time.diff <= 13 & 
                 !is.na(x = toxo.status)) %>%
        group_by(sex, age.cat.dart) %>%
        summarise (n.C = sum(!is.na(c)),
                   avg.C = round (mean(c, 
                                       na.rm = T),2),
                   stdev.C = round (sd(c, 
                                       na.rm = T), 2),
                   med.C = round(median(c,
                                        na.rm = T), 2),
                   min.C = round(min(c,
                                     na.rm = T), 2),
                   max.C = round(max(c,
                                     na.rm = T), 2)) %>%
        ungroup()
      
    ## e) Histogram plasma Cortisol
      hist_plot_plasma_c <- ggplot(data=subset(plasma_horm_neosp_toxo_data,
                                                 dart.time.diff <= 13 & 
                                                 stressca <=2 &
                                                 !is.na(x = spratio) &
                                                 !is.na(x = c.ln)),
                                   aes(x=c )) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 20, by = 1), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0, 20)) +
        labs(title = 'Histogram of Plasma Cortisol 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Fec. Cortisol (ng/g)', y='Frequency') 
      
      print(hist_plot_plasma_c)
      
#********************** Data Inclusion/Exclusion Criteria ********************** 
      
    ## f) Set Zero's to 1/2 the minimum detected value
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(c.adj = ifelse(c == 0, 
                              (0.5*min(plasma_horm_neosp_toxo_data[,'c']
                                [which(plasma_horm_neosp_toxo_data[,'c']> 0)])),
                              plasma_horm_neosp_toxo_data$c))
      
#********************** Data Inclusion/Exclusion Criteria **********************      
      
    ## g) Natural log transformation of cort data 
      plasma_horm_neosp_toxo_data $c.ln <- 
        log(plasma_horm_neosp_toxo_data$c.adj )
      
      plasma_horm_neosp_toxo_data $c.imp.ln <- 
        log(plasma_horm_neosp_toxo_data$c.imp )
      
    ## h) Histogram of Nat. Log. of plasma Cortisol
      hist_plot_ln_plasma_c <- ggplot(data=subset(plasma_horm_neosp_toxo_data,
                                                  dart.time.diff <= 13 & 
                                                    stressca <=2 &
                                                    !is.na(x = spratio) &
                                                    !is.na(x = c.ln)), 
                               aes(x=c.ln)) + 
        geom_histogram(aes(y = ..count..),
                       breaks=seq(0, 4, by = 0.125), 
                       col='black',
                       fill = 'dark grey') +
        xlim(c(0, 4)) +
        labs(title = 'Histogram of Nat. Log. Plasma Cortisol 
(repeated measures)') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        labs(x='Nat. Log. Fec. Cortisol (ng/g)', y='Frequency') 
      
      print(hist_plot_ln_plasma_c)
      
    ## i) Save histogram plot
      # use ggsave to save the plot
      ggsave('histogram_ln_cort.pdf', plot = hist_plot_ln_plasma_c, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 5, 
             height = 3, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
#************************** Data Manipulation Criteria ************************* 
    ## j) Descriptive stats plasma corticosterone (C) to determine cutoffs
      sex_age_plasma_C_cutoffs <- plasma_horm_neosp_toxo_data %>%
        filter(c > 0 & stressca <=2 & dart.time.diff <= 13 & 
                 !is.na(x = toxo.status)) %>%
        group_by(sex) %>%
        summarise (n.C = sum(!is.na(c)),
                   avg.C = round (mean(c, 
                                       na.rm = T),2),
                   stdev.C = round (sd(c, 
                                       na.rm = T), 2),
                   med.C = round(median(c,
                                        na.rm = T), 2),
                   min.C = round(min(c,
                                     na.rm = T), 2),
                   max.C = round(max(c,
                                     na.rm = T), 2)) %>%
        ungroup()
      
    ## k) Categorize cort-levels as low vs hi
      #*** BASED on within sex median values***
      plasma_horm_neosp_toxo_data <- plasma_horm_neosp_toxo_data  %>%
        mutate(c.bin = case_when(plasma_horm_neosp_toxo_data$sex == 'f' &
                                   plasma_horm_neosp_toxo_data$c <= 1.92
                                 ~ 0,
                                 plasma_horm_neosp_toxo_data$sex == 'f' &
                                   plasma_horm_neosp_toxo_data$c > 1.92
                                 ~ 1,
                                 plasma_horm_neosp_toxo_data$sex == 'm' &
                                   plasma_horm_neosp_toxo_data$c <= 0.90
                                 ~ 0,
                                 plasma_horm_neosp_toxo_data$sex == 'm' &
                                   plasma_horm_neosp_toxo_data$c > 0.90
                                 ~ 1))
    # ## l) Re-code c.bin as nominal factor and set level (order)
    #   plasma_horm_neosp_toxo_data <- transform(plasma_horm_neosp_toxo_data,
    #                                            c.bin = factor(c.bin,
    #                                                   levels = c('hi','low')))
      
#************************** Data Manipulation Criteria *************************
      
  ### 3.3 Univariate stats T. gondii 
    ## a) Descriptive stats T. gondii infection prevalence
      univar_toxo_plasma_stat <- plasma_horm_neosp_toxo_data  %>%
        group_by(toxo.status) %>%
        summarise(n.status = sum(!is.na(toxo.status))) %>%
        mutate(freq = n.status / sum(n.status, na.rm = T))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/univar_toxo_plasma_stat.pdf'), height = 4, width = 8)
      grid.table(univar_toxo_plasma_stat)
      dev.off()
      
      
  # ### 3.4 Univariate stats N. caninum 
  #   ## a) Descriptive stats N. caninum  infection prevalence
  #     univar_neosp_plasma_stat <- plasma_horm_neosp_toxo_data  %>%
  #       group_by(neo.status) %>%
  #       summarise(n.status = sum(!is.na(neo.status))) %>%
  #       mutate(freq = n.status / sum(n.status, na.rm = T))
  #     
  #   ## b) save the data frame of summary stats out as a pdf into output file
  #     pdf(here('output/univar_neosp_plasma_stat.pdf'), height = 4, width = 8)
  #     grid.table(univar_neosp_plasma_stat)
  #     dev.off()
      
      

###############################################################################
##############               4. Bivariate analysis               ##############
############################################################################### 
      
  ### 4.1 Descriptive bivariate plots testosterone status by sex
    ## a) Scatter plot plasma testosterone by age by sex  
      plasma_T_scatter <- ggplot(plasma_horm_neosp_toxo_data, 
                              aes(x = age.cat.dart, y = t, 
                              shape = dart.state, color = sex)) +
        geom_point(position = position_jitter(h=0.0, w=0.3), cex = 2) +
        # geom_hline(yintercept = round(sd(
        #   plasma_horm_neosp_toxo_data$t, 
        #   na.rm = T), 2) * 4, 
        #   linetype="solid", color = "red") +
        # annotate(geom="text", label = '4 SD', x = 1, 
        #          y = round(sd(plasma_horm_neosp_toxo_data$t, 
        #                       na.rm = T), 2) * 4, vjust=-1) +
        # 
        # geom_hline(yintercept = round(sd(
        #   plasma_horm_neosp_toxo_data$t, 
        #   na.rm = T), 2) * 5,
        #   linetype="dashed", color = "red") +
        # annotate(geom="text", label = '5 SD', x = 1, 
        #          y = round(sd(plasma_horm_neosp_toxo_data$t, 
        #                       na.rm = T), 2) * 5, vjust=-1) +
        labs(title = 'Scatter plot of plasma testosterone by age and sex',
             x ='Hyena age class (darting date)', 
             y ='Plasma testosterone ug/dL') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) 
      # add major axes
      
      print(plasma_T_scatter)
      
      ## b) Save scatter plot
      ggsave('plasma_T_scatter.pdf', plot = plasma_T_scatter, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      
      
      
  ### 4.2 Descriptive bivariate plots Cortisol status by sex
    ## a) Scatter plot plasma Cortisol by age by sex  
      plasma_cort_scatter <- ggplot(data = filter(plasma_horm_neosp_toxo_data, 
                                                  dart.time.diff <= 13 & 
                                                    stressca <=2), 
                                 aes(x = age.cat.dart, y = c, 
                                     shape = dart.state, color = sex)) +
        geom_point(position = position_jitter(h=0.0, w=0.3), cex = 2) +
        # geom_hline(yintercept = round(sd(
        #   plasma_horm_neosp_toxo_data$c, 
        #   na.rm = T), 2) * 4, 
        #   linetype="solid", color = "red") +
        # annotate(geom="text", label = '4 SD', x = 1, 
        #          y = round(sd(plasma_horm_neosp_toxo_data$c, 
        #                       na.rm = T), 2) * 4, vjust=-1) +
        # 
        # geom_hline(yintercept = round(sd(
        #   plasma_horm_neosp_toxo_data$c, 
        #   na.rm = T), 2) * 5,
      #   linetype="dashed", color = "red") +
      # annotate(geom="text", label = '5 SD', x = 1, 
      #          y = round(sd(plasma_horm_neosp_toxo_data$c, 
      #                       na.rm = T), 2) * 5, vjust=-1) +
      labs(title = 'Scatter plot of plasma cortisol by age and sex',
           x ='Hyena age class (darting date)', 
           y ='Plasma cortisol ug/dL') +
        theme(plot.title = element_text(hjust = 0.5)) + # center title
        theme(axis.ticks = element_blank()) + # remove axis ticks
        # remove background color
        theme(panel.background = element_rect(fill = 'white')) 
      # add major axes
      
      print(plasma_cort_scatter)
      
    ## b) Save scatter plot
      ggsave('plasma_cort_scatter.pdf', plot = plasma_cort_scatter, 
             device = NULL, 
             path = here('output/'), scale = 1, width = 8, 
             height = 5, 
             units = c('in'), dpi = 300, limitsize = TRUE) 
      

  ### 4.3 Descriptive bivariate stats testosterone
    ## a) Testosterone by toxo.status by sex by age
      plasma_T_toxo_sex_age_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(t)) %>%
        group_by (sex, age.cat.dart, toxo.status) %>%
        summarise (n.id = sum(!is.na(t)),
                   avg.T = round(mean(t, na.rm = T),2),
                   stdev.T = round(sd(t, na.rm = T), 2),
                   min.T = round(min(t, na.rm = T), 2),
                   med.T = round(median(t, na.rm = T), 2),
                   max.T = round(max(t, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_toxo_sex_age_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_toxo_sex_age_sum)
      dev.off()
      
    ## c) Testosterone by sex by age where t > 0
      plasma_T_non_zero_sex_age_sum <- plasma_horm_neosp_toxo_data %>%
        filter(t > 0 & !is.na(x = toxo.status)) %>%
        #filter(!is.na(t) & sex == 'm' & age.cat.dart == 'subadult') %>%
        group_by (sex, age.cat.dart) %>%
        summarise (n.id = sum(!is.na(t)),
                   avg.T = round(mean(t, na.rm = T),2),
                   stdev.T = round(sd(t, na.rm = T), 2),
                   min.T = round(min(t, na.rm = T), 2),
                   med.T = round(median(t, na.rm = T), 2),
                   max.T = round(max(t, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_non_zero_sex_age_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_non_zero_sex_age_sum)
      dev.off()
      
    ## e) Nat. log testosterone by sex by age
      plasma_T_sex_age_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(t.ln)) %>%
       #filter(!is.na(t.ln) & sex == 'm' & age.cat.dart == 'subadult') %>%
        group_by (sex, age.cat.dart) %>%
        summarise (n.id = sum(!is.na(t)),
                   avg.T = round(mean(t.ln, na.rm = T),2),
                   stdev.T = round(sd(t.ln, na.rm = T), 2),
                   min.T = round(min(t.ln, na.rm = T), 2),
                   med.T = round(median(t.ln, na.rm = T), 2),
                   max.T = round(max(t.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## f) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_sex_age_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_sex_age_sum)
      dev.off()
      
    ## g) Nat. log Testosterone by reproductive state (adult females only)
      plasma_T_state_sum_adult_f <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(t.ln) & sex == 'f' & 
                 age.cat.dart == 'adult') %>%
        group_by (dart.state) %>%
        summarise (n.id = sum(!is.na(t.ln)),
                   avg.T = round(mean(t.ln, na.rm = T),2),
                   stdev.T = round(sd(t.ln, na.rm = T), 2),
                   min.T = round(min(t.ln, na.rm = T), 2),
                   med.T = round(median(t.ln, na.rm = T), 2),
                   max.T = round(max(t.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## h) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_state_sum_adult_f.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_state_sum_adult_f)
      dev.off()
      
    ## i) Nat. log Testosterone by residency (adult males only)
      plasma_T_status_sum_adult_m <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(t.ln) & sex == 'm' & 
                 age.cat.dart == 'adult') %>%
        group_by (status) %>%
        summarise (n.id = sum(!is.na(t.ln)),
                   avg.T = round(mean(t.ln, na.rm = T),2),
                   stdev.T = round(sd(t.ln, na.rm = T), 2),
                   min.T = round(min(t.ln, na.rm = T), 2),
                   med.T = round(median(t.ln, na.rm = T), 2),
                   max.T = round(max(t.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## j) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_status_sum_adult_m.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_status_sum_adult_m)
      dev.off()
      
    ## k) Nat. log Testosterone by time of day
      plasma_T_time_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(t.ln)) %>%
        group_by (sex, dart.am.pm) %>%
        summarise (n.id = sum(!is.na(t.ln)),
                   avg.T = round(mean(t.ln, na.rm = T),2),
                   stdev.T = round(sd(t.ln, na.rm = T), 2),
                   min.T = round(min(t.ln, na.rm = T), 2),
                   med.T = round(median(t.ln, na.rm = T), 2),
                   max.T = round(max(t.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## l) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_time_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_time_sum)
      dev.off()
      
    ## m) Nat. log Testosterone by migration season (adult females and males only)
      plasma_T_migrtn_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(t.ln)) %>%
        group_by (sex, migratn.seas.dart) %>%
        summarise (n.id = sum(!is.na(t.ln)),
                   avg.T = round(mean(t.ln, na.rm = T),2),
                   stdev.T = round(sd(t.ln, na.rm = T), 2),
                   min.T = round(min(t.ln, na.rm = T), 2),
                   med.T = round(median(t.ln, na.rm = T), 2),
                   max.T = round(max(t.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## n) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_T_migrtn_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_T_migrtn_sum)
      dev.off()  
      
      
  ### 4.4 Descriptive bivariate stats Cortisol    
      ## a) Cortisol by toxo.status by sex by age
      plasma_cort_toxo_sex_age_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(c) & dart.time.diff <= 13 & stressca <=2) %>%
        group_by (sex, age.cat.dart, toxo.status) %>%
        summarise (n.id = sum(!is.na(c)),
                   avg.cort = round(mean(c, na.rm = T),2),
                   stdev.cort = round(sd(c, na.rm = T), 2),
                   min.cort = round(min(c, na.rm = T), 2),
                   med.cort = round(median(c, na.rm = T), 2),
                   max.cort = round(max(c, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## b) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_toxo_sex_age_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_toxo_sex_age_sum)
      dev.off() 
      
    ## c) Cortisol by sex by age where c > 0
      plasma_cort_non_zero_sex_age_sum <- plasma_horm_neosp_toxo_data %>%
        filter(c > 0 & stressca <=2 & dart.time.diff <= 13 & 
                 !is.na(x = toxo.status)) %>%
        #filter(!is.na(c) & sex == 'm' & age.cat.dart == 'subadult') %>%
        group_by (sex, age.cat.dart) %>%
        summarise (n.id = sum(!is.na(c)),
                   avg.cort = round(mean(c, na.rm = T),2),
                   stdev.cort = round(sd(c, na.rm = T), 2),
                   min.cort = round(min(c, na.rm = T), 2),
                   med.cort = round(median(c, na.rm = T), 2),
                   max.cortT = round(max(c, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_non_zero_sex_age_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_non_zero_sex_age_sum)
      dev.off()
      
    ## c) Cortisol by sex by age
      plasma_cort_sex_age_sum <- plasma_horm_neosp_toxo_data %>%
        filter(stressca <=2 & dart.time.diff <= 13 & 
                 !is.na(x = toxo.status)) %>%
        group_by (sex, age.cat.dart) %>%
        summarise (n.id = sum(!is.na(c.ln)),
                   avg.cort = round(mean(c.ln, na.rm = T),2),
                   stdev.cort = round(sd(c.ln, na.rm = T), 2),
                   min.cort = round(min(c.ln, na.rm = T), 2),
                   med.cort = round(median(c.ln, na.rm = T), 2),
                   max.cort = round(max(c.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## d) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_sex_age_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_sex_age_sum)
      dev.off()
      
     ## e) Cortisol by reproductive state (adult females only)
      plasma_cort_state_sum_adult_f <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(c.ln) & sex == 'f' & age.cat.dart == 'adult'
               & dart.time.diff <= 13 & stressca <=2) %>%
        group_by (dart.state) %>%
        summarise (n.id = sum(!is.na(c.ln)),
                   avg.cort = round(mean(c.ln, na.rm = T),2),
                   stdev.cort = round(sd(c.ln, na.rm = T), 2),
                   min.cort = round(min(c.ln, na.rm = T), 2),
                   med.cort = round(median(c.ln, na.rm = T), 2),
                   max.cort = round(max(c.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## f) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_state_sum_adult_f.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_state_sum_adult_f)
      dev.off()
      
    ## g) Cortisol by residency (adult males only)
      plasma_cort_status_sum_adult_m <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(c.ln) & sex == 'm' & age.cat.dart == 'adult'
               & dart.time.diff <= 13 & stressca <=2) %>%
        group_by (status) %>%
        summarise (n.id = sum(!is.na(c.ln)),
                   avg.cort = round(mean(c.ln, na.rm = T),2),
                   stdev.cort = round(sd(c.ln, na.rm = T), 2),
                   min.cort = round(min(c.ln, na.rm = T), 2),
                   med.cort = round(median(c.ln, na.rm = T), 2),
                   max.cort = round(max(c.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## h) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_status_sum_adult_m.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_status_sum_adult_m)
      dev.off()
      
      
    ## i) Cortisol by time of day
      plasma_cort_time_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(c.ln)
               & dart.time.diff <= 13 & stressca <=2) %>%
        group_by (sex, dart.am.pm) %>%
        summarise (n.id = sum(!is.na(c.ln)),
                   avg.cort = round(mean(c.ln, na.rm = T),2),
                   stdev.cort = round(sd(c.ln, na.rm = T), 2),
                   min.cort = round(min(c.ln, na.rm = T), 2),
                   med.cort = round(median(c.ln, na.rm = T), 2),
                   max.cort = round(max(c.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## j) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_time_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_time_sum)
      dev.off()
      
    ## k) Testosterone by migration season 
      plasma_cort_migrtn_sum <- plasma_horm_neosp_toxo_data %>%
        filter(!is.na(c.ln) 
               & dart.time.diff <= 13 & stressca <=2) %>%
        group_by (sex, migratn.seas.dart) %>%
        summarise (n.id = sum(!is.na(c.ln)),
                   avg.cort = round(mean(c.ln, na.rm = T),2),
                   stdev.cort = round(sd(c.ln, na.rm = T), 2),
                   min.cort = round(min(c.ln, na.rm = T), 2),
                   med.cort = round(median(c.ln, na.rm = T), 2),
                   max.cort = round(max(c.ln, na.rm =T), 2))%>%
        mutate(freq = n.id / sum(n.id))
      
    ## l) save the data frame of summary stats out as a pdf into output file
      pdf(here('output/plasma_cort_migrtn_sum.pdf'), 
          height = 4, width = 8)
      grid.table(plasma_cort_migrtn_sum)
      dev.off()  
   
      
      
###############################################################################
##############               5. Export data files                ##############
###############################################################################
      
    ### 5.1 Export data to an RData file     
      ## a) Save and export raw data tables 
      # Files are saved in the 'data' folder in the working directory as an
      # RData file.
      save(file = paste0(project_data_path, 
                         '3_neo_toxo_plasma_horm.RData'), 
           list = c('plasma_horm_neosp_toxo_data'))
      

     
      
      