#setwd('D:/Covid19_vaccine efficacy/public_R_codes/')
source('preliminary_functions.R')
source('code_KMbased_Coxbased_InfluFuncs.R')

package_list = c('survival', 'MASS', 'Matrix', 'dplyr', 'lubridate', 'tidyr', 'purrr')
for (package in package_list){
  require(package, character.only = TRUE); library(package, character.only = TRUE) }

num_digits = 7; options(digits = num_digits)

'%ni%' <- Negate('%in%')


Compute_final_pval = function(est_val, se_val, mu_val){
  pval = pnorm(abs((est_val-mu_val)/se_val), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
  return( pval )
}

Do_random_intersection_test = function(dat, t_val, z_val, pre_trt, fix_window, adj_window, prob, num_rep, 
                                       alpha_val=input_alpha){
  
  trt_arms = sort(unique(dat$A))[-1]; rest_trts = trt_arms[trt_arms != pre_trt]; alt_trt = rest_trts[1]
  w_pre_trt = Get_Wset_treatment(pre_trt, dat); adj_windows = setdiff(w_pre_trt, fix_window) 
  w_alt_trt = Get_Wset_treatment(alt_trt, dat); alt_window = setdiff(w_alt_trt, c(fix_window, adj_window))
  pre_window = setdiff(adj_windows, adj_window)
  
  sample_df_pre = as.data.frame( dat %>% filter(W %in% adj_windows) )
  n = dim(sample_df_pre %>% filter(A==pre_trt))[1]; 
  fixed_df = as.data.frame( dat %>% filter(W %in% fix_window) )
  
  ret1 = do.call(rbind, lapply(1:num_rep, function(i){ set.seed(i)
    resample_prep = sample_df_pre %>% filter(W %in% adj_window) %>% mutate(A_label=A) %>%
      group_by(A) %>% nest() %>% arrange(A) %>%
      mutate(size = map_int(data, nrow)) %>% mutate(keep_size = min(size,round(n*prob,0))) %>%
      mutate(kick_size = size-keep_size) %>% mutate(keep_sample = map2(data, keep_size, sample_n)) %>%
      mutate(kick_sample = map2(data, keep_sample, function(df, alt_df){
                                                     df[!row.names(df) %in% row.names(alt_df),]}))
    
    temp1 = as.data.frame(do.call(rbind, resample_prep$keep_sample)) %>% rename(A = A_label)
    kick_sample = as.data.frame(do.call(rbind, resample_prep$kick_sample)) %>% rename(A = A_label)
    temp2 = as.data.frame( bind_rows( kick_sample %>% filter(A %in% c(0,pre_trt)) %>% mutate(W=pre_window),
                                      sample_df_pre %>% filter(W %in% pre_window) ) )
    temp3 = as.data.frame( bind_rows( kick_sample %>% filter(A %in% c(0,alt_trt)) %>% mutate(W=alt_window),
                                      dat %>% filter(W %in% alt_window) ) )
    
    resample_df = rbind(temp1, temp2, temp3, fixed_df)
    
    obj <- Estimation$new(input_data=resample_df)
    
    rest_tests = do.call(rbind, lapply(rest_trts, function(rest_trt){
      w_rest_trt = Get_Wset_treatment(rest_trt, dat=resample_df)
      
      CSR_est = obj$Estimate_contrast_stratRRs(t=t_val, w_pre_trt, w_rest_trt, pre_trt, rest_trt,
                                               z=z_val, trial_type='platform', contrast_func='difference', 
                                               alpha=alpha_val, Delta_val=1)
      est = CSR_est$est; se = CSR_est$se; ci = CSR_est$ci;
      
      w_intersect= intersect(w_pre_trt, w_rest_trt); w_union = union(w_pre_trt, w_rest_trt)
      shared_placebo_n = dim(resample_df %>% filter(W %in% w_intersect & A==0))[1]
      total_placebo_n = dim(resample_df %>% filter(W %in% w_union & A==0))[1]
      
      sharing_rate = shared_placebo_n/total_placebo_n
      
      return( data.frame('est'=est, 'se'=se, 'lb'=ci[1], 'ub'=ci[2], 'sharing_rate'=sharing_rate) ) }))
    
    mean_est = mean(rest_tests[,'est']); mean_se = mean(rest_tests[,'se']);
    mean_sharing_rate = mean(rest_tests[,'sharing_rate'])  
    min_lb = min(rest_tests[,'lb']) ; max_ub = max(rest_tests[,'ub']); sample_used = dim(resample_df)[1]
    
    return( data.frame('sharing_rate'=mean_sharing_rate, 'sample_used'=sample_used, 'N'=obj$N, 'est'=mean_est, 
                       'ci_width'=mean_se, 'lb'=min_lb, 'ub'=max_ub, 'idx'=i) ) }))

  return( as.data.frame(ret1) )
}

Test_under_random_configuration = function(dat, z_val, settings, window_setup, num_rep){ 
  ret0 = do.call(rbind, lapply(1:nrow(settings), function(i){ #print(i)
    # Extract the needed variables from one input data frame
    for(name in names(settings)){ assign(paste0(name), settings[i,name]) }
    
    ret1 = do.call(rbind, lapply(1:nrow(window_setup), function(idx){
      # if(idx %in% c(6,18)) { print(idx) }
      # Extract the needed variables from the other input data frame
      for(name in names(window_setup)){ assign(paste0(name), window_setup[idx,name]) }
      
      keep_condition = ((dat$W!=absence_window_1 & dat$A==1) | (dat$W!=absence_window_2 & dat$A==2) | (dat$A==0) )
      dat0 = dat[keep_condition,]
      ret2 = Do_random_intersection_test(dat0, t, z_val, pre_trt, fix_window_1, adj_window, prob, num_rep) 
      return( ret2 ) }))
    
    ret3 = cbind( data.frame('t'=t, 'pre_trt'=pre_trt, 'prob'=prob),
                  as.data.frame( t(colMeans(ret1)) ) )
    return( ret3 ) }))
  return( ret0 )
}

