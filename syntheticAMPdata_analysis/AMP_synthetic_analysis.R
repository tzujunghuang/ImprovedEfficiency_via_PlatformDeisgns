##### This script illustrates the code of a single Monte Carlo run 
##### under the setting in our submitted manuscript 
##### or a selected setting for a quick coding check.

##### Should be in the same working directory	
wk_dir = 'D:/Covid19_vaccine efficacy/public_R_codes/' # where files are saved
setwd(wk_dir)
to_dir = wk_dir
source('AMPdata_analysis_functions.R')

num_digits = 7; options(digits = num_digits)


# Methods for efficacy comparison
meth = 'intersection_testing';

# Configuration
t = 560; k = 2; m = 4; # k: Number of interventions; m: Number of windows
pre_trts = c(1) # Pre-selected treatments to be compared
pre_trt_type = 'low'
tolerances = c(1) # Selected tolerance criteria
prob_vec = seq(0.1, 0.6, 0.1)

setting_index = 'selected_setting'
setting_choice = switch(setting_index,
  ### The full setting in the manuscript
  full_setting = { nr = 500 # Number of replication
  },   
  ### The selected simplified setting 
  selected_setting = { nr = 50 # Number of replication.
  },
  stop("setting_choice: choose full_setting or selected_setting"))

# set-up of windows
window_setup = expand.grid(1:m, 1:m, 1:m, 1:m); 
colnames(window_setup) = c('absence_window_1', 'fix_window_1', 'absence_window_2', 'fix_window_2') 
window_setup = ( window_setup %>% filter(absence_window_1 != fix_window_1 & absence_window_1 != absence_window_2
                                         & absence_window_2 != fix_window_2 & fix_window_1 == fix_window_2)
                 %>% arrange(absence_window_1, absence_window_2) )
window_setup$adj_window = apply(window_setup, 1, function(row){setdiff(1:m, row)})

settings = expand.grid(t, pre_trts, tolerances, prob_vec); 
colnames(settings) = c('t', 'pre_trt', 'tolerance', 'prob')
settings = settings[with(settings, order(t, pre_trt, tolerance, prob)),]

input_alpha = 0.025 # Significance level


# Load in the simulated data
setwd('D:/Covid19_vaccine efficacy/public_R_codes/syntheticAMPdata_analysis/')
filename = 'syntheticAMPdata.Rdata'
sim_data = as.data.frame( get(load(filename)) )
title_text = 'simdata_'

# Run the analysis based on the built functions
measure_type = 'stratified'


ret = do.call(rbind, lapply(unique(sim_data$Z), function(z0){ 
  data_z = sim_data %>% filter(Z==z0)
  ret_z = Test_under_random_configuration(data_z, z0, settings, window_setup, num_rep=nr)
  ret_z$Z = z0
  return( ret_z ) }))

save(ret, file=paste0(title_text, 'ciwidths_', measure_type, '_', pre_trt_type, '.Rdata'))

