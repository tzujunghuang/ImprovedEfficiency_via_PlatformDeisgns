##### This script provides the code of a single Monte Carlo run 
##### under the full setting in our submitted manuscript 
##### or a selected setting for a quick coding check.

set.seed(r)

##### Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'dplyr', 'PWEALL', 'tidyr')
for (package in package_list){
  require(package, character.only=TRUE); library(package, character.only=TRUE) }
num_digits = 7; options(digits = num_digits)


##### Should be in the same working directory	
wk_dir = 'D:/Covid19_vaccine efficacy/public R codes/' # where files are saved
setwd(wk_dir)
to_dir = wk_dir
source('preliminary_functions.R')
source('generate_survivaldata.R')
source('code_KMbased_Coxbased_InfluFuncs.R')

##### Methods to evaluate efficacy
meth_vecs = c('cov_stratRR', 'cov_adjRR', 'cov_stratRR_nullZ') 

##### Simulation settings
trial_types = c('platform', 'separate') # Trial types
loss_to_follow_rate = '10%' # Loss_to_follow_rates = c('2%', '10%')
loss_to_follow_pars = c(600, 120) # Loss_to_follow_parameters corresponding to '2%' and '10%'
z = 1 # Selected z
trt0 = 7;  rest_trts = c(1:6, 8:10) # Example treatment and the rest

setting_index = 'selected_setting'
setting_choice = switch(setting_index,
### The full setting in the manuscript
  full_setting = { adm_censoring_calendartimes = c(6, 9, 12, 18) # Administrative censoring times 
                   t = c(3, 6) # Selected t's
                 },   
### The selected simplified setting 
  selected_setting = { adm_censoring_calendartimes = c(6) # Administrative censoring times 
                         t = c(6) # Selected t's
                     },
  stop("setting_choice: choose full_setting or selected_setting"))

settings = expand.grid(t, z, trt0, rest_trts)
colnames(settings) = c('t', 'z', 'a1', 'a2')
settings = settings[with(settings, order(t, a2)),]
alpha_one_side = 0.025 # Significance level for one-side

##### Data generation settings
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


sim = data.frame( method=NA, trial_type=NA, adm_censoring_calendartime=NA, t=NA, z=NA, a1=NA, a2=NA,
                  contrast=NA, est=NA, se=NA, ci_lb=NA, ci_ub=NA )

for (trial_type in trial_types){
  for (adm_censoring_calendartime in adm_censoring_calendartimes){
    
    dat_obj = Generate_SurvivalData$new(treatment_labels_in_windows=input_treatment_labels_in_windows,
                                        incidence_rates=input_incidence_rates, window_bands=input_window_bands, 
                                        n_m=input_n_m, p=1)
    dat = dat_obj$Generate_data(trial_type, loss_to_follow_rate, adm_censoring_calendartime, true_data=FALSE, 
                                loss_to_follow_pars=loss_to_follow_pars)
    obj <- Estimation$new(input_data=dat)
    print(c(trial_type, loss_to_follow_rate, adm_censoring_calendartime))
    
    out = list()
    for (meth in meth_vecs) {
      if (meth=='cov_stratRR') {
        
        print(meth)
        
        CSR_ret = do.call(rbind, lapply(1:dim(settings)[1], function(i){
          t = settings[i,'t']; z = settings[i,'z']; a1 = settings[i,'a1']; a2 = settings[i,'a2']
          ws1 = Get_Wset_treatment(a1, dat); ws2 = Get_Wset_treatment(a2, dat)
          CSR_est = obj$Estimate_contrast_stratRRs(t, ws1, ws2, a1, a2, z, trial_type, contrast_func='ratio', 
                                                   alpha=alpha_one_side, Delta_val=1)
          t(c(t, a1, a2, 'ratio', CSR_est$est, CSR_est$se, CSR_est$ci, z)) }))
      
        out$cov_stratRR = CSR_ret
        
      } else if (meth=='cov_stratRR_nullZ'){
        
        print(meth)
        
        dat1 = dat; dat1$Z = 1
        print(paste0('Z-value = ', unique(dat1$Z)))
        obj1 <- Estimation$new(input_data=dat1)
        
        CSR_ret_nullZ = do.call(rbind, lapply(1:dim(settings)[1], function(i){
          t = settings[i,'t']; z = settings[i,'z']; a1 = settings[i,'a1']; a2 = settings[i,'a2']
          ws1 = Get_Wset_treatment(a1, dat1); ws2 = Get_Wset_treatment(a2, dat1)
          CSR_est1 = obj1$Estimate_contrast_stratRRs(t, ws1, ws2, a1, a2, z, trial_type, contrast_func='ratio', 
                                                     alpha=alpha_one_side, Delta_val=1)
          t(c(t, a1, a2, 'ratio', CSR_est1$est, CSR_est1$se, CSR_est1$ci, z)) }))
        
        out$cov_stratRR_nullZ = CSR_ret_nullZ
        
      } else if (meth=='cov_adjRR') {
        
        print(meth)
        
        settings0 = subset(settings, select=-c(z))
        CAR_ret = do.call(rbind, lapply(1:dim(settings0)[1], function(i){
          t = settings0[i,'t']; a1 = settings0[i,'a1']; a2 = settings0[i,'a2']
          ws1 = Get_Wset_treatment(a1, dat); ws2 = Get_Wset_treatment(a2, dat)
          CAR_est = obj$Estimate_contrast_adjustRRs(t, ws1, ws2, a1, a2, trial_type, contrast_func='ratio', 
                                                    alpha=alpha_one_side, Delta_val=1)
          t(c(t, a1, a2, 'ratio', CAR_est$est, CAR_est$se, CAR_est$ci)) }))
        
        out$cov_adjRR = CAR_ret
        
      } else stop ('Invalid method.') }
    gc()
    
    sim = rbind(sim, do.call(rbind,lapply(meth_vecs, function(meth) {
      if ( meth=='cov_stratRR' ) {
        
        result = data.frame( method=meth, trial_type=trial_type, adm_censoring_calendartime=adm_censoring_calendartime, 
                             t=out[[meth]][,1], z=out[[meth]][,9], a1=out[[meth]][,2], a2=out[[meth]][,3],
                             contrast=out[[meth]][,4], est=out[[meth]][,5], se=out[[meth]][,6],
                             ci_lb=out[[meth]][,7], ci_ub=out[[meth]][,8] )
      
      } else if ( meth=='cov_stratRR_nullZ' ) {
        
        result = data.frame( method=meth, trial_type=trial_type, adm_censoring_calendartime=adm_censoring_calendartime, 
                             t=out[[meth]][,1], z=out[[meth]][,9], a1=out[[meth]][,2], a2=out[[meth]][,3],
                             contrast=out[[meth]][,4], est=out[[meth]][,5], se=out[[meth]][,6],
                             ci_lb=out[[meth]][,7], ci_ub=out[[meth]][,8] )
      
      } else if (meth=='cov_adjRR') {
        
        result = data.frame( method=meth, trial_type=trial_type, adm_censoring_calendartime=adm_censoring_calendartime, 
                             t=out[[meth]][,1], z=NA, a1=out[[meth]][,2], a2=out[[meth]][,3],
                             contrast=out[[meth]][,4], est=out[[meth]][,5], se=out[[meth]][,6],
                             ci_lb=out[[meth]][,7], ci_ub=out[[meth]][,8] ) } 
      return( result )} )))
  }
}

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
