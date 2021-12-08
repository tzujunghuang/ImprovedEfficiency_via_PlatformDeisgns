##### Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'dplyr', 'PWEALL', 'tidyr')
for (package in package_list){
  require(package, character.only = TRUE); library(package, character.only = TRUE) }
num_digits = 7; options(digits = num_digits)


##### Should be in the same working directory
wk_dir = 'D:/Covid19_vaccine efficacy/public R codes/' # where files are saved
setwd(wk_dir)
to_dir = wk_dir
source('preliminary_functions.R')
source('generate_survivaldata.R')

##### Data generation settings
num_rep = 10 # Number of replications
trial_types = c('platform', 'separate') # Trial types
k = 10 # Number of interventions
m = 5 # Number of windows
p = 1 # Number of predictors
input_n_m=c(1000, 900, 1200, 1400, 1000) # Basic unit sizes of windows  
# Treatment labels in windows
input_treatment_labels_in_windows = list('window1' = c(2, 3, 4, 7, 8, 9, 10),
                                         'window2' = c(1, 4, 5, 6, 7, 9, 10),
                                         'window3' = c(4, 5, 6, 7, 9, 10),
                                         'window4' = c(1, 3, 6, 8, 9, 10),
                                         'window5' = c(1, 2, 3, 7, 8, 10))
# Incidence rates over windows for placebo and treatments
input_incidence_rates = matrix((1/100)*c(  2.1,   2.1,    1,    0.7,    0.7,
                                           2.1,   2.1,    1,    0.7,    0.7,
                                           1.89,  1.89,  0.9,   0.63,   0.63,
                                           1.68,  1.68,  0.8,   0.56,   0.56,
                                           1.47,  1.47,  0.7,   0.49,   0.49,
                                           1.26,  1.26,  0.6,   0.42,   0.42,
                                           1.05,  1.05,  0.5,   0.35,   0.35,
                                           0.84,  0.84,  0.4,   0.28,   0.28,
                                           0.84,  0.84,  0.4,   0.28,   0.28,
                                           0.63,  0.63,  0.3,   0.21,   0.21,
                                           0.63,  0.63,  0.3,   0.21,   0.21), byrow=TRUE, ncol=m)
rownames(input_incidence_rates) = 0:k
# Bands and length of windows
input_window_bands = data.frame(lower_end = c(0, 1, 1.5, 2, 2.5), upper_end = c(1, 1.5, 2, 2.5, 3))
rownames(input_window_bands) = 1:m
input_window_bands$length = input_window_bands$upper_end - input_window_bands$lower_end

##### Grids at which we generate true parameters
preferred_trt = 10
t = c(3, 6); # Selected t's
z = 1; # Selected z
tz_combns = expand.grid(t,z); colnames(tz_combns) = c('t','z')
pre_trts = c(3, 5, 7, 9); # Pre-selected treatments to be compared
trts = 1:k # All the treatments 
trt_combns = do.call(rbind, lapply(pre_trts, function(pre_trt){ data.frame('a1'=pre_trt, 'a2'=trts[trts!=pre_trt]) } ))

settings = do.call(rbind, lapply(1:dim(tz_combns)[1], function(i){ 
                     data.frame('t'=tz_combns[i,'t'], 'z'=tz_combns[i,'z'], trt_combns) }))



true_pars0 = NULL; true_pars1 = NULL
for (trial_type in trial_types) {
  
  if (trial_type=='platform') { NN=250 } else if (trial_type=='separate'){ NN=150 }
  
  ret0 = NULL; ret1 = NULL
  for (j in 1:num_rep) {
    
    true_dat_obj = Generate_SurvivalData$new(treatment_labels_in_windows=input_treatment_labels_in_windows,
                                             incidence_rates=input_incidence_rates, window_bands=input_window_bands, 
                                             n_m=NN*input_n_m, p=1)
    true_dat = true_dat_obj$Generate_data(trial_type, loss_to_follow_rate=NA, adm_censoring_calendartime=NA, true_data=TRUE) 
    print( c(trial_type, j, dim(true_dat)[1]) )
    
##### Generate parameters for Section 5, ignoring the role of Z. 
    if (trial_type=='separate') { 
    } else if (trial_type=='platform') {
      settings0 = subset(settings, select=-c(z))
      true_RRs = do.call(rbind, lapply(1:dim(settings0)[1], function(i){
      t = settings0[i,'t']; a1 = settings0[i,'a1']; ws1 = Get_Wset_treatment(a1, true_dat)
      a2 = settings0[i,'a2']; ws2 = Get_Wset_treatment(a2, true_dat)
        
      rr1 = True_RelativeRisk(t, ws1, a1, trial_type, data=true_dat)  
      rr2 = True_RelativeRisk(t, ws2, a2, trial_type, data=true_dat)
      #rr_diff = True_contrast_RRs(t, ws1, ws2, a1, a2, contrast_func='difference', trial_type, data=true_dat)
        
      return( data.frame('t'=t, 'a1'=a1, 'a2'=a2, 'rr1'=rr1, 'rr2'=rr2) ) }))
      
      ret0 = rbind(ret0, data.frame('rep'=j, true_RRs))  }  
    
##### Generate parameters for Section 4.
    settings1 = settings[settings[,'a1']==7,]; a1 = 7; ws1 = Get_Wset_treatment(a1, true_dat)
    true_RRs = do.call(rbind, lapply(1:dim(settings1)[1], function(i){
      t = settings1[i,'t']; z = settings1[i,'z']; a2 = settings1[i,'a2']
      ws2 = Get_Wset_treatment(a2, true_dat)
    
      srr1 = True_StratRelativeRisk(t, ws1, a1, z, trial_type, data=true_dat)
      srr2 = True_StratRelativeRisk(t, ws2, a2, z, trial_type, data=true_dat)
      #srr_ratio = True_contrast_StratRRs(t, ws1, ws2, a1, a2, z, contrast_func='ratio', trial_type, data=true_dat)
                          
      rr1 = True_RelativeRisk(t, ws1, a1, trial_type, data=true_dat)  
      rr2 = True_RelativeRisk(t, ws2, a2, trial_type, data=true_dat)
      #rr_ratio = True_contrast_RRs(t, ws1, ws2, a1, a2, contrast_func='ratio', trial_type, data=true_dat)
                          
      return( data.frame('t'=t, 'a1'=a1, 'a2'=a2, 'srr1'=srr1, 'srr2'=srr2, 'rr1'=rr1, 'rr2'=rr2) ) }))
    
    ret1 = rbind(ret1, data.frame('rep'=j, true_RRs)) 
  }
  ret1$trial_type = trial_type
  
  true_pars0 = rbind(ret0, true_pars0); true_pars1 = rbind(ret1, true_pars1)
}

save(true_pars0, file='pars_efficacycompetition.Rdata')
save(true_pars1, file='pars_efficiencygain.Rdata')


##### Collect and summarize the results
load('pars_efficiencygain.Rdata')
true_pars_dat = true_pars1
settings1 = settings[settings[,'a1']==7,]

true_pars = do.call(rbind, lapply(trial_types, function(trial_type){
  do.call(rbind, lapply(1:dim(settings1)[1], function(i){  
    t = settings1[i,'t']; a1 = settings1[i,'a1']; a2 = settings1[i,'a2']
    condition = (true_pars_dat$trial_type==trial_type & true_pars_dat$t==t & true_pars_dat$a1==a1 & 
                 true_pars_dat$a2==a2)
    true_pars_vals = as.data.frame( t(apply(true_pars_dat[condition, c('srr1','srr2','rr1','rr2')], 2, mean)) )
    return( data.frame( 'trial_type'=trial_type, 't'=t, 'a1'=a1, 'a2'=a2, 
                        'strat_RRratio'=true_pars_vals[,'srr1']/true_pars_vals[,'srr2'], 
                        'adj_RRratio'=true_pars_vals[,'rr1']/true_pars_vals[,'rr2'],
                        'strat_RRdiff'=true_pars_vals[,'srr1']-true_pars_vals[,'srr2'],
                        'adj_RRdiff'=true_pars_vals[,'rr1']-true_pars_vals[,'rr2'] ) ) }))
 }))

save(true_pars, file='true_pars_efficiencygain.Rdata')


load('pars_efficacycompetition.Rdata')
true_pars_dat = true_pars0
true_pars = do.call(rbind, lapply(1:dim(settings)[1], function(i){  
    t = settings[i,'t']; a1 = settings[i,'a1']; a2 = settings[i,'a2']
    condition = (true_pars_dat$t==t & true_pars_dat$a1==a1 & true_pars_dat$a2==a2)
    
    true_pars_vals = as.data.frame( t(apply(true_pars_dat[condition, c('rr1','rr2')], 2, mean)) )
    return( data.frame( 't'=t, 'a1'=a1, 'a2'=a2, 
                        'RRratio'=true_pars_vals[,'rr1']/true_pars_vals[,'rr2'],
                        'RRdiff'=true_pars_vals[,'rr1']-true_pars_vals[,'rr2'] ) ) }))

save(true_pars, file='true_pars_efficacycompetition.Rdata')

