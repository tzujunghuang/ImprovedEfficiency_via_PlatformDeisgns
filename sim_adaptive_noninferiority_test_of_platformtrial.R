##### This script provides the code of a single Monte Carlo run 
##### under the setting in our submitted manuscript 
##### or a selected setting for a quick coding check.

set.seed(r)

##### Make sure we have every package needed
package_list = c('survival', 'MASS', 'Matrix', 'dplyr', 'PWEALL', 'tidyr', 'mvtnorm', 'nloptr')
for (package in package_list){
  require(package, character.only=TRUE); library(package, character.only=TRUE) }
options(nloptr.show.inequality.warning=FALSE)
num_digits = 7; options(digits = num_digits)


##### Should be in the same working directory	
wk_dir = 'D:/Covid19_vaccine efficacy/public R codes/' # where files are saved
setwd(wk_dir)
to_dir = wk_dir
source('preliminary_functions.R')
source('generate_survivaldata.R')
source('code_KMbased_Coxbased_InfluFuncs.R')

##### Methods for efficacy comparison
meth_vecs = c('intersection_test', 'likelihood_ratio_test')

##### Simulation settings
trial_types = c('platform', 'separate')
loss_to_follow_rate = '10%' # Loss_to_follow_rates = c('2%', '10%')
loss_to_follow_pars = c(600, 120) # Loss_to_follow_parameters corresponding to '2%' and '10%'
tolerances = c(0.7) # Selected tolerance criteria delta = 0.7

setting_index = 'selected_setting'
setting_choice = switch(setting_index,
### The full setting in the manuscript
  full_setting = { adm_censoring_calendartimes = c(6, 18) # Administrative censoring times 
                   t = c(3, 6) # Selected t's
                   pre_trts = c(7, 9) # Pre-selected treatments to be compared
                   margins = seq(0, 0.4, 0.1) # Selected non-inferiority margins
                 },   
### The selected simplified setting 
  selected_setting = { adm_censoring_calendartimes = c(6) # Administrative censoring times 
                       t = c(6) # Selected t's
                       pre_trts = c(7) # Pre-selected treatments to be compared
                       margins = seq(0, 0.2, 0.1) # Selected non-inferiority margins
                     },
  stop("setting_choice: choose full_setting or selected_setting"))

settings = expand.grid(t, pre_trts, margins, tolerances); 
colnames(settings) = c('t', 'pre_trt', 'margin', 'tolerance')
settings = settings[with(settings, order(t, pre_trt, margin, tolerance)),]
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
input_incidence_rates = matrix((1/100)*c( 2.1,   2.1,    1,    0.7,    0.7,
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


col_names = c('trial_type', 'adm_censoring_calendartime', 't', 'pre_trt', 'margin',
              'tolerance', 'method', 'init_rej', 'rej', paste0('rest_rej', 1:(k-1)))
sim = data.frame( matrix(nrow=1, ncol=length(col_names)) ); colnames(sim) = col_names

for (trial_type in trial_types){
  for (adm_censoring_calendartime in adm_censoring_calendartimes){
  
    dat_obj = Generate_SurvivalData$new(treatment_labels_in_windows = input_treatment_labels_in_windows,
                                        incidence_rates = input_incidence_rates, 
                                        window_bands = input_window_bands, n_m = input_n_m, p=1)
    dat = dat_obj$Generate_data(trial_type, loss_to_follow_rate, adm_censoring_calendartime, 
                                loss_to_follow_pars = loss_to_follow_pars, true_data=FALSE)
    dat$Z = 1

    obj <- Estimation$new(input_data=dat)
    N_val = obj$N;
    print(c(trial_type, loss_to_follow_rate, adm_censoring_calendartime))
   
    if (trial_type=='platform'){
      
      est_cov_mats_list = lapply(unique(settings[,'t']), function(time){ 
        est_cov = obj$Estimate_cov_mat(t=time, z=1) 
        est_cov })
      names(est_cov_mats_list) = unique(settings[,'t'])
      
      trt_arms = sort(unique(dat$A))[-1]; 
      r_n_list = lapply(unique(settings[,'t']), function(time){ 
        r_n = sapply(trt_arms, function(trt){ w_trt = Get_Wset_treatment(trt, dat)
          SR_est = obj$Estimate_stratRRs(time, w_trt, trt, z=1, trial_type, alpha=alpha_one_side, Delta_val=1)
          est_val = SR_est$est; 
          return( est_val ) })
        r_n })
      names(r_n_list) = unique(settings[,'t'])
      
      ret = do.call(rbind, lapply(1:dim(settings)[1], function(i){
        t = settings[i,'t']; pre_trt = settings[i,'pre_trt']; margin = settings[i,'margin'];
        tolerance = settings[i,'tolerance']; w_pre_trt = Get_Wset_treatment(pre_trt, dat); 
        trt_arms = sort(unique(dat$A))[-1]; rest_trts = sort(trt_arms[trt_arms != pre_trt]) 
        
        ret0 = do.call(rbind, lapply(meth_vecs, function(meth){
          if (meth=='intersection_test'){
            SR_est = obj$Estimate_stratRRs(t, w_pre_trt, pre_trt, z=1, trial_type, alpha=alpha_one_side, 
                                           Delta_val=1)
            init_test = 1*( SR_est$ci[2] < tolerance ); init_rej = init_test 
            
            CSRests = do.call(rbind, lapply(rest_trts, function(rest_trt){
              w_rest_trt = Get_Wset_treatment(rest_trt, dat)
              CSR_est = obj$Estimate_contrast_stratRRs(t, w_pre_trt, w_rest_trt, pre_trt, rest_trt, 
                                                       z=1, trial_type, contrast_func='difference', 
                                                       alpha=alpha_one_side, Delta_val=1)
              ci = CSR_est$ci; rej = 1*( ci[2] < margin )
              return( rej ) }))
            rest_tests = CSRests
            
            rej = 1*( all(rest_tests==1) & init_test==1 )
            
            temp = data.frame( t(c(t, pre_trt, margin, tolerance, meth, init_rej, rej, rest_tests)) )
            colnames(temp) = col_names[3:length(col_names)]
            return( temp )
            
          } else if (meth=='likelihood_ratio_test'){
            Dmat = solve( est_cov_mats_list[[as.character(t)]]/N_val )
            r_n = r_n_list[[as.character(t)]]
            
            obj_fn = function(g, r0=r_n, Dmat0=Dmat){ t(r0-g) %*% Dmat0 %*% (r0-g)/2 }
            constraint_hn0 = function(g, tol=tolerance, margn=margin){ 
              g[pre_trt] - min(tol, min(g[-pre_trt])+margn) }          
            quad_sol0 = cobyla(x0=r_n, fn=obj_fn, hin=constraint_hn0, 
                               control=list(xtol_rel=1e-4, maxeval=1000)) 
            g_0 = quad_sol0$par
            Tn =  t(r_n-g_0) %*% Dmat %*% (r_n-g_0) 
            rej = 1*( Tn >= qchisq(1-alpha_one_side, df=k, ncp=0, lower.tail=TRUE, log.p=FALSE) )
            
            temp = data.frame( t(c(t, pre_trt, margin, tolerance, meth, init_rej=NA, rej, rep(NA,k-1))) )
            colnames(temp) = col_names[3:length(col_names)]
            return( temp )
          }
        }))
        return( ret0 )
      }))
    
    } else if (trial_type=='separate'){
      
      trt_arms = sort(unique(dat$A))[-1]; 
      rets_list = lapply(unique(settings[,'t']), function(time){ 
        ret = do.call(rbind, lapply(trt_arms, function(trt){ w_trt = Get_Wset_treatment(trt, dat)
          SR_est = obj$Estimate_stratRRs(time, w_trt, trt, z=1, trial_type, alpha=alpha_one_side, Delta_val=1)
          est_val = SR_est$est; var_val = (SR_est$sd)^2;
          return( data.frame('est_val'=est_val, 'var_val'=var_val) ) }))
        return( ret ) })
      names(rets_list) = unique(settings[,'t'])
    
      ret = do.call(rbind, lapply(1:dim(settings)[1], function(i){
        t = settings[i,'t']; pre_trt = settings[i,'pre_trt']; margin = settings[i,'margin'];
        tolerance = settings[i,'tolerance']; w_pre_trt = Get_Wset_treatment(pre_trt, dat);
        trt_arms = sort(unique(dat$A))[-1]; rest_trts = sort(trt_arms[trt_arms != pre_trt])  
        
        ret0 = do.call(rbind, lapply(meth_vecs, function(meth){
          if (meth=='intersection_test'){
            SR_est = obj$Estimate_stratRRs(t, w_pre_trt, pre_trt, z=1, trial_type, alpha=alpha_one_side, 
                                           Delta_val=1)
            init_test = 1*( SR_est$ci[2] < tolerance ); init_rej = init_test 
            
            CSRests = do.call(rbind, lapply(rest_trts, function(rest_trt){
              w_rest_trt = Get_Wset_treatment(rest_trt, dat)
              CSR_est = obj$Estimate_contrast_stratRRs(t, w_pre_trt, w_rest_trt, pre_trt, rest_trt, 
                                                       z=1, trial_type, contrast_func='difference', 
                                                       alpha=alpha_one_side, Delta_val=1)
              ci = CSR_est$ci; rej = 1*( ci[2] < margin )
              return( rej ) }))
            rest_tests = CSRests
            rej = 1*( all(rest_tests==1) & init_test==1 )
            
            temp = data.frame( t(c(t, pre_trt, margin, tolerance, meth, init_rej, rej, rest_tests)) )
            colnames(temp) = col_names[3:length(col_names)]  
            return( temp )
            
          } else if (meth=='likelihood_ratio_test'){
            Dmat = diag( N_val/unlist(rets_list[[as.character(t)]]['var_val']) )
            r_n = unlist(rets_list[[as.character(t)]][['est_val']])
            
            obj_fn = function(g, r0=r_n, Dmat0=Dmat){ t(r0-g) %*% Dmat0 %*% (r0-g)/2 }
            constraint_hn0 = function(g, tol=tolerance, margn=margin){ 
              g[pre_trt] - min(tol, min(g[-pre_trt])+margn) }
            quad_sol0 = cobyla(x0=r_n, fn=obj_fn, hin=constraint_hn0, 
                               control=list(xtol_rel=1e-4, maxeval=1000)); 
            g_0 = quad_sol0$par
            Tn =  t(r_n-g_0) %*% Dmat %*% (r_n-g_0) 
            rej = 1*( Tn >= qchisq(1-alpha_one_side, df=k, ncp=0, lower.tail=TRUE, log.p=FALSE) )
            
            temp = data.frame( t(c(t, pre_trt, margin, tolerance, meth, init_rej=NA, rej, rep(NA,k-1))) ) 
            colnames(temp) = col_names[3:length(col_names)]  
            return( temp )
          }  
        })) 
        return( ret0 )
      })) 
    }   
    gc()
    
    sim0 = ret; sim0$trial_type = trial_type; sim0$adm_censoring_calendartime = adm_censoring_calendartime
    sim = rbind(sim, sim0)
 }
}

sim = sim[-1,]
rownames(sim) = 1:nrow(sim)
