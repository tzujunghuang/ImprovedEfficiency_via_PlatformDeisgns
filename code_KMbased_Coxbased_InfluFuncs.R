source('preliminary_functions.R')
num_digits = 7
library('survival')

Estimation <- setRefClass( "Estimation",
  fields = list(
    input_data = "data.frame",
    N = "numeric"
  ),
                               
  methods = list(
  initialize = function(input_data = dat){ input_data <<- input_data; N <<- dim(input_data)[1] },

  Strat_KMSurF = function(t, ws, a, z, trial_type, label_a, Delta_val=1){ 
    # ws: a column vector with >=1 element
    if (trial_type=='platform') { data_used = subset(input_data, W %in% ws & A==a & Z==z, select=c(X, Delta))
    } else if (trial_type=='separate') { 
      data_used = subset(input_data, W %in% ws & A==a & Z==z & Label_A==label_a, select=c(X, Delta)) }
    
    data_km = data.frame(X=data_used$X, Delta=1*(data_used$Delta==Delta_val))
    km = survfit(Surv(X, Delta)~1, data=data_km); rm(data_km);
    survest = cbind(km$time, km$surv)
      
    if ( length(which(survest[, 1] <= t)) > 0 ) {
      return( survest[max(which(survest[, 1] <= t)), 2] )  
    } else { stop('Survival probability too close to zero') }
  },
  
  Strat_CHF = function(t, ws, a, z, trial_type, label_a, Delta_val=1){ # ws: a column vector with >=1 element
    if (trial_type=='platform') { data_used = subset(input_data, W %in% ws & A==a & Z==z) 
                                  data_used = data_used[order(data_used$X),]; N0 = N
    } else if (trial_type=='separate') { data_used0 = subset(input_data, Label_A==label_a)
                                         data_used = subset(data_used0, W %in% ws & A==a & Z==z) 
                                         data_used = data_used[order(data_used$X),]
                                         N0 = dim(data_used0)[1] }
    
    X = as.vector(data_used$X); Delta = as.vector(1*(data_used$Delta==Delta_val))
    A = as.vector(data_used$A); W = as.vector(data_used$W); Z = as.vector(data_used$Z)
    
    time_comparison1 = diag(length(X))
    time_comparison2 = outer(X, X, '>=')
    
    event_counting = colSums(time_comparison1*Delta*(W %in% ws & A==a & Z==z)) / N0; 
    risk_set = colSums(time_comparison2*(W %in% ws & A==a & Z==z)) / N0
    return( sum((event_counting / risk_set)*(X <= t & risk_set > 0)) )
  },
  
  Strat_RelativeRisk = function(t, ws, a, z, trial_type, Delta_val=1){ # ws: a column vector with >=1 element
    SurF0_z = Strat_KMSurF(t, ws, a=0, z, trial_type, label_a=a, Delta_val); 
    SurF1_z = Strat_KMSurF(t, ws, a, z, trial_type, label_a=a, Delta_val)
    return( (1-SurF1_z)/(1-SurF0_z) )
  },  
  
  Adjust_RelativeRisk = function(t, ws, a, trial_type, Delta_val=1){ # ws: a column vector with >=1 element
    if (trial_type=='platform') { data_used = input_data[order(input_data$X),]
    } else if (trial_type=='separate') { 
      data_used = subset(input_data, Label_A==a); data_used = data_used[order(data_used$X),] }

    W = as.vector(data_used$W); Z = as.vector(data_used$Z); z_vals = sort(as.numeric(unique(Z)))
    
    temp = do.call(rbind, lapply(1:length(z_vals), function(i) { 
      z = z_vals[i]; p_zw = mean(Z==z & W %in% ws); p_w = mean(W %in% ws); 
      p_z_given_w = p_zw/p_w 
      SurF0_z = Strat_KMSurF(t, ws, a=0, z, trial_type, label_a=a, Delta_val); 
      SurF1_z = Strat_KMSurF(t, ws, a, z, trial_type, label_a=a, Delta_val)
      data.frame('top'=(1-SurF1_z)*p_z_given_w, 'bottom'=(1-SurF0_z)*p_z_given_w) }))
    return( sum(temp['top'])/sum(temp['bottom']) )
  },
  
  IF_stratCHF = function(t, ws, a, z, trial_type, label_a, result_type, Delta_val=1){ 
    # ws: a column vector with >=1 element
    if (trial_type=='platform') { data_used = subset(input_data, W %in% ws & A==a & Z==z) 
                                  data_used = data_used[order(data_used$X),]; N0 = N
    } else if (trial_type=='separate') { data_used0 = subset(input_data, Label_A==label_a)
                                         data_used = subset(data_used0, W %in% ws & A==a & Z==z) 
                                         data_used = data_used[order(data_used$X),] 
                                         N0 = dim(data_used0)[1] }
    
    X = as.vector(data_used$X); Delta = as.vector(1*(data_used$Delta==Delta_val))
    A = as.vector(data_used$A); W = as.vector(data_used$W); Z = as.vector(data_used$Z)
    
    time_comparison1 = diag(length(X))
    time_comparison2 = outer(X, X, '>=')
      
    event_counting = colSums(time_comparison1*Delta*(W %in% ws & A==a & Z==z))/N0; 
    risk_set = colSums(time_comparison2*(W %in% ws & A==a & Z==z))/N0
    if_1 = 1*(X<=t & Delta==Delta_val & W %in% ws & A==a & Z==z)/risk_set
    t0 = pmin(t, X)
    if_2 = cumsum((event_counting / (risk_set^2))*(X <= t0 & risk_set > 0)) #elements in risk_set >=1
      
    ret <- switch(result_type,
                  iid = { data.frame('a'=a, 'time'=X, 'iid_vals'=if_1-if_2) },
                  sum = { sum(if_1-if_2) },
                  square_sum = { sum((if_1-if_2)^2) },
                  stop("IF_stratCHF: choose iid, sum or square_sum"))
    return( ret ) 
  },
  
  IF_stratKM = function(t, ws, a, z, trial_type, label_a, result_type, Delta_val=1){ 
    # ws: a column vector with >=1 element
    S1_z = Strat_KMSurF(t, ws, a, z, trial_type, label_a, Delta_val);
    ret <- switch(result_type,
                  iid = { iid_data = IF_stratCHF(t, ws, a, z, trial_type, label_a, 'iid', Delta_val) 
                          data.frame('a'=iid_data[,'a'], 'time'=iid_data[,'time'], 
                                     'iid_vals'=S1_z*iid_data[,'iid_vals']) },
                  sum = { S1_z*IF_stratCHF(t, ws, a, z, trial_type, label_a, 'sum', Delta_val) },
                  square_sum = { ((S1_z)^2)*IF_stratCHF(t, ws, a, z, trial_type, label_a, 'square_sum', 
                                                        Delta_val) },
                  stop("IF_stratKM: choose iid, sum or square_sum"))
    return( ret )
  },
    
  IF_zw = function(ws, a, z, trial_type, result_type){ # ws: a column vector with >=1 element
    if (trial_type=='platform') { data_used = input_data[order(input_data$X),]
    } else if (trial_type=='separate') { data_used = subset(input_data, Label_A==a) 
    data_used = data_used[order(data_used$X),] }
    
    X = as.vector(data_used$X); A = as.vector(data_used$A); W = as.vector(data_used$W); 
    Z = as.vector(data_used$Z)
    
    p_zw = mean(Z==z & W %in% ws); p_w = mean(W %in% ws); 
    p_z_given_w = p_zw/p_w
    if_zw = (1*(Z==z & W %in% ws) - p_zw)/p_w -(1*(W %in% ws) - p_w)*p_z_given_w/p_w
    
    condition = (W %in% ws & A %in% c(a,0) & Z==z)
    A0 = A[condition]; time = X[condition]; iid_vals = if_zw[condition]
    
    condition1 = (W %in% ws & A %in% c(a,0))
    A1 = A[condition1]; time1 = X[condition1]; iid_vals1 = if_zw[condition1]
    
    ret <- switch(result_type,
                  iid = { temp = data.frame('a'=A0, 'time'=time, 'iid_vals'=iid_vals) 
                          temp[with(temp, order(-a, time)),] },
                  iid1 = { temp = data.frame('a'=A1, 'time'=time1, 'iid_vals1'=iid_vals1) 
                           temp[with(temp, order(-a, time)),] },
                  sum = { sum(iid_vals) },
                  square_sum = { sum((iid_vals)^2) },
                  stop("IF_zw: choose iid, sum or square_sum"))
    return( ret ) 
  },
    
  IF_stratRR = function(t, ws, a, z, trial_type, result_type, Delta_val=1){ 
    # ws: a column vector with >=1 element
    SurF0_z = Strat_KMSurF(t, ws, a=0, z, trial_type, label_a=a, Delta_val); 
    SurF1_z = Strat_KMSurF(t, ws, a, z, trial_type, label_a=a, Delta_val)
    R1 = SurF1_z/(1-SurF0_z); R2 = (1-SurF1_z)*SurF0_z/(1-SurF0_z)^2
    
    ret <- switch(result_type,
                  iid = { iid_data_a = IF_stratCHF(t, ws, a, z, trial_type, label_a=a, 'iid', Delta_val) 
                          iid_data_0 = IF_stratCHF(t, ws, a=0, z, trial_type, label_a=a, 'iid', Delta_val) 
                          data.frame( 'a' = c(iid_data_a[,'a'],iid_data_0[,'a']),
                                      'time' = c(iid_data_a[,'time'], iid_data_0[,'time']),
                                      'iid_vals' = c(-R1*iid_data_a[,'iid_vals'], 
                                                     +R2*iid_data_0[,'iid_vals']) ) },
                  sum = { ( -R1*IF_stratCHF(t, ws, a, z, trial_type, label_a=a, 'sum', Delta_val) 
                            +R2*IF_stratCHF(t, ws, a=0, z, trial_type, label_a=a, 'sum', Delta_val) ) },
                  square_sum = { ( (R1^2)*IF_stratCHF(t, ws, a, z, trial_type, label_a=a, 
                                                      'square_sum', Delta_val) 
                                   +(R2^2)*IF_stratCHF(t, ws, a=0, z, trial_type, label_a=a, 'square_sum', 
                                                       Delta_val) ) },
                  stop("IF_stratRR: choose iid, sum or square_sum"))
    return( ret )
  },
    
  IF_stratLogRR = function(t, ws, a, z, trial_type, result_type, Delta_val=1){ 
    # ws: a column vector with >=1 element
    srr = Strat_RelativeRisk(t, ws, a, z, Delta_val); 
    ret <- switch(result_type,
                  iid = { iid_data = IF_stratRR(t, ws, a, z, trial_type, 'iid', Delta_val) 
                          data.frame( 'a' = iid_data[,'a'],
                                      'time' = iid_data[,'time'],
                                      'iid_vals' = iid_data[,'iid_vals']/srr ) },
                  sum = { IF_stratRR(t, ws, a, z, trial_type, 'sum', Delta_val)/srr },
                  square_sum = { IF_stratRR(t, ws, a, z, trial_type, 'square_sum', Delta_val)/(srr)^2 },
                  stop("IF_stratLogRR: choose iid, sum or square_sum"))
    return( ret )
  },
  
  IF_adjustRR = function(t, ws, a, trial_type, result_type, Delta_val=1){ 
    # ws: a column vector with >=1 element
    if (trial_type=='platform') { data_used = input_data[order(input_data$X),]
    } else if (trial_type=='separate') { data_used = subset(input_data, Label_A==a) 
    data_used = data_used[order(data_used$X),] }
    
    W = as.vector(data_used$W); Z = as.vector(data_used$Z); z_vals = sort(as.numeric(unique(Z)))
    arr = Adjust_RelativeRisk(t, ws, a, trial_type, Delta_val)
    
    temp0 = do.call(rbind, lapply(1:length(z_vals), function(i) {
               z = z_vals[i]; p_zw = mean(Z==z & W %in% ws); p_w = mean(W %in% ws); p_z_given_w = p_zw/p_w 
               SurF0_z = Strat_KMSurF(t, ws, a=0, z, trial_type, label_a=a, Delta_val); 
               SurF1_z = Strat_KMSurF(t, ws, a, z, trial_type, label_a=a, Delta_val)
               
               iid_data_0 = IF_stratCHF(t, ws, a=0, z, trial_type, label_a=a, 'iid', Delta_val)
               iid_data_a = IF_stratCHF(t, ws, a, z, trial_type, label_a=a, 'iid', Delta_val)
               iid_data_zw = IF_zw(ws, a, z, trial_type, 'iid')
               
               condition1 = ( iid_data_zw[,'a'] == c(iid_data_a[,'a'], iid_data_0[,'a']) )
               condition2 = ( iid_data_zw[,'time'] == c(iid_data_a[,'time'], iid_data_0[,'time']) )
                 
               if (all(condition1 & condition2)) {
                 iid_vals_z = ( c( -p_z_given_w*SurF1_z*iid_data_a[,'iid_vals'],
                                   +arr*p_z_given_w*SurF0_z*iid_data_0[,'iid_vals'] )
                                + ((1-SurF1_z)-arr*(1-SurF0_z))*iid_data_zw[,'iid_vals'] )
               } else { stop ('Treatment labels and observed times unmatch.') }   
               
               return( data.frame('a'=iid_data_zw[,'a'], 'time'=iid_data_zw[,'time'], 'z'=z,
                                  'iid_vals_z'=iid_vals_z) ) }))
    
    temp1 = do.call(rbind, lapply(1:length(z_vals), function(i) { 
               z = z_vals[i]; p_zw = mean(Z==z & W %in% ws); p_w = mean(W %in% ws); p_z_given_w = p_zw/p_w 
               SurF0_z = Strat_KMSurF(t, ws, a=0, z, trial_type, label_a=a, Delta_val);
               return( data.frame('z'=z, 'p_z_given_w'=p_z_given_w, 'SurF0_z'=SurF0_z) ) }))
    Q0 = 1/(1-sum((temp1$SurF0_z)*(temp1$p_z_given_w)))
      
    ret <- switch(result_type,
                  iid = { data.frame( 'a'=temp0[,'a'], 'time'=temp0[,'time'], 
                                      'iid_vals'=Q0*temp0[,'iid_vals_z'] ) },
                  sum = { sum(Q0*temp0[,'iid_vals_z']) },
                  square_sum = { sum((Q0*temp0[,'iid_vals_z'])^2) },
                  stop("IF_adjustRR: choose iid, sum or square_sum"))
    return( ret ) 
  },
  
  Estimate_cov_mat = function(t, z) { 
    trt_arms = sort(unique(input_data$A))[-1]
    est_vars = do.call(rbind, lapply(trt_arms, function(trt){
                  w_trt = Get_Wset_treatment(trt, input_data)
                  est_var_trt = ( IF_stratRR(t, w_trt, trt, z, trial_type='platform', 'square_sum', 
                                             Delta_val=1)/N 
                                  - (IF_stratRR(t, w_trt, trt, z, trial_type='platform', 'sum', 
                                                Delta_val=1)/N)^2 )
                  return( est_var_trt ) }))
    var_mat = diag(as.vector(est_vars))
  
    trt_combns = t(combn(trt_arms,2)) 
    est_covs = do.call(rbind, lapply(1:dim(trt_combns)[1], function(i){
                  a1 = trt_combns[i,1]; a2 = trt_combns[i,2]
                  ws1 = Get_Wset_treatment(a1, input_data); ws2 = Get_Wset_treatment(a2, input_data)
                  iid_data_a1 = IF_stratRR(t, ws1, a1, z, trial_type='platform', 'iid', Delta_val=1)
                  iid_data_a2 = IF_stratRR(t, ws2, a2, z, trial_type='platform', 'iid', Delta_val=1)
                  unique_common_times = sort(unique(subset(input_data, W %in% intersect(ws1, ws2) & A==0)$X))
      
                  if (length(unique_common_times)==0) { est_cov_a1a2 = 0
                  } else { 
                    temp = do.call(rbind, lapply(1:length(unique_common_times), function(i){
                      time = unique_common_times[i]
                      condition1 = (iid_data_a1$a==0 & iid_data_a1$time==time)
                      condition2 = (iid_data_a2$a==0 & iid_data_a2$time==time)
                      
                      iid_vals_a1 = iid_data_a1[condition1, 'iid_vals']
                      iid_vals_a2 = iid_data_a2[condition2, 'iid_vals']
                      min_length = min(length(iid_vals_a1), length(iid_vals_a2))
                      
                      return( data.frame( 'time' = time,
                                          'iid_a1' = iid_vals_a1[1:min_length],
                                          'iid_a2' = iid_vals_a2[1:min_length] ) ) }))
                    est_cov_a1a2 = sum(temp[,'iid_a1']*temp[,'iid_a2'])/N - 
                                   (sum(temp[,'iid_a1'])/N)*(sum(temp[,'iid_a2'])/N) }
                  return( data.frame( 'loc1'=a1, 'loc2'=a2, 'val'=est_cov_a1a2 ) ) }))
  
    cov_mat = do.call(rbind, lapply(1:length(trt_arms), function(i){
                 trt = trt_arms[i]; 
                 if ((i+1) <= dim(var_mat)[2]) { row = c(rep(0,i), est_covs[est_covs['loc1']==trt,'val'])
                 } else { row = c(rep(0, dim(var_mat)[2])) }
                 return( row ) }))
    return( var_mat + cov_mat + t(cov_mat) )
  },
  
  Estimate_stratRRs = function(t, ws, a, z, trial_type, alpha, Delta_val=1){
    # ws1, ws2: column vectors with >=1 element
    if (trial_type=='platform') { N1=N;
    } else if (trial_type=='separate') {
      data_used1 = subset(input_data, Label_A==a); N1=dim(data_used1)[1] }
    
    est = Strat_RelativeRisk(t, ws, a, z, trial_type, Delta_val)
    est_var0 = ( IF_stratRR(t, ws, a, z, trial_type, 'square_sum', Delta_val)/N1 
                 - (IF_stratRR(t, ws, a, z, trial_type, 'sum', Delta_val)/N1)^2 )

    if (trial_type=='platform') { est_var = est_var0
    } else if (trial_type=='separate'){
      probs_separate_trials = Get_probs_separate_trials(input_data)
      sum_probs_separate_trials = sum(probs_separate_trials[,'prob'])
      prob_separate_trial = probs_separate_trials[probs_separate_trials$a==a,'prob']
      est_var = (sum_probs_separate_trials/prob_separate_trial)*est_var0 }
    
    ci = c( est-qnorm(1-alpha)*sqrt(est_var)/sqrt(N), est+qnorm(1-alpha)*sqrt(est_var)/sqrt(N) )
    pval = pnorm(abs(sqrt(N)*est/sqrt(est_var)), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
    
    return( list(est=est, sd=sqrt(est_var), se=sqrt(est_var)/sqrt(N), ci=ci, pval=pval) )
  },
  
  Estimate_contrast_stratRRs = function(t, ws1, ws2, a1, a2, z, trial_type, contrast_func, alpha, 
                                        Delta_val=1){
    # ws1, ws2: column vectors with >=1 element
    if (trial_type=='platform') { N1=N2=N;
    } else if (trial_type=='separate') {
        data_used1 = subset(input_data, Label_A==a1); N1=dim(data_used1)[1]
        data_used2 = subset(input_data, Label_A==a2); N2=dim(data_used2)[1] }
    
    est_stratRR1 = Strat_RelativeRisk(t, ws1, a1, z, trial_type, Delta_val)
    est_stratRR2 = Strat_RelativeRisk(t, ws2, a2, z, trial_type, Delta_val)
    est = if (contrast_func=='difference') { est_stratRR1-est_stratRR2 
    } else if (contrast_func=='ratio') { est_stratRR1/est_stratRR2 }
    
    if (contrast_func=='difference') { Theta1=1; Theta2=-1;
    } else if (contrast_func=='ratio') { Theta1=1/est_stratRR2; Theta2=-est_stratRR1/(est_stratRR2)^2 }
    
    est_var_a1 = ( IF_stratRR(t, ws1, a1, z, trial_type, 'square_sum', Delta_val)/N1 
                   - (IF_stratRR(t, ws1, a1, z, trial_type, 'sum', Delta_val)/N1)^2 )
    est_var_a2 = ( IF_stratRR(t, ws2, a2, z, trial_type, 'square_sum', Delta_val)/N2 
                   - (IF_stratRR(t, ws2, a2, z, trial_type, 'sum', Delta_val)/N2)^2 )
    
    if (trial_type=='platform') {
      iid_data_a1 = IF_stratRR(t, ws1, a1, z, trial_type, 'iid', Delta_val)
      iid_data_a2 = IF_stratRR(t, ws2, a2, z, trial_type, 'iid', Delta_val)
      unique_common_times = sort(unique(subset(input_data, W %in% intersect(ws1, ws2) & A==0 & Z==z)$X))
      if (length(unique_common_times)==0) { est_cov_a1a2 = 0
      } else{ temp = do.call(rbind, lapply(1:length(unique_common_times), function(i){
                        time = unique_common_times[i]
                        condition1 = (iid_data_a1$a==0 & iid_data_a1$time==time)
                        condition2 = (iid_data_a2$a==0 & iid_data_a2$time==time)
                        
                        iid_vals_a1 = iid_data_a1[condition1, 'iid_vals']
                        iid_vals_a2 = iid_data_a2[condition2, 'iid_vals']
                        min_length = min(length(iid_vals_a1), length(iid_vals_a2))
                      
                        return( data.frame( 'time' = time,
                                            'iid_a1' = iid_vals_a1[1:min_length],
                                            'iid_a2' = iid_vals_a2[1:min_length] ) ) }))
      est_cov_a1a2 = sum(temp[,'iid_a1']*temp[,'iid_a2'])/N - (sum(temp[,'iid_a1'])/N)*(sum(temp[,'iid_a2'])/N) }
      est_var = ((Theta1)^2)*est_var_a1 + ((Theta2)^2)*est_var_a2 + 2*Theta1*Theta2*est_cov_a1a2
    } else if (trial_type=='separate'){
      probs_separate_trials = Get_probs_separate_trials(input_data)
      sum_probs_separate_trials = sum(probs_separate_trials[,'prob'])
      prob_a1_separate_trial = probs_separate_trials[probs_separate_trials$a==a1,'prob'] 
      prob_a2_separate_trial = probs_separate_trials[probs_separate_trials$a==a2,'prob'] 
      
      est_var = ( (sum_probs_separate_trials/prob_a1_separate_trial)*((Theta1)^2)*est_var_a1 
                  + (sum_probs_separate_trials/prob_a2_separate_trial)*((Theta2)^2)*est_var_a2 ) }

    ci = c( est-qnorm(1-alpha)*sqrt(est_var)/sqrt(N), est+qnorm(1-alpha)*sqrt(est_var)/sqrt(N) )
    pval = pnorm(abs(sqrt(N)*est/sqrt(est_var)), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
    
    return( list(est=est, sd=sqrt(est_var), se=sqrt(est_var)/sqrt(N), ci=ci, pval=pval) )
  },
  
  Estimate_contrast_adjustRRs = function(t, ws1, ws2, a1, a2, trial_type, contrast_func, alpha, Delta_val=1){
    # ws1, ws2: column vectors with >=1 element
    if (trial_type=='platform') { N1=N2=N
    } else if (trial_type=='separate') {
    data_used1 = subset(input_data, Label_A==a1); N1=dim(data_used1)[1]
    data_used2 = subset(input_data, Label_A==a2); N2=dim(data_used2)[1] }
    
    est_adjRR1 = Adjust_RelativeRisk(t, ws1, a1, trial_type, Delta_val)
    est_adjRR2 = Adjust_RelativeRisk(t, ws2, a2, trial_type, Delta_val)
    est = if (contrast_func=='difference') { est_adjRR1-est_adjRR2
      } else if (contrast_func=='ratio') { est_adjRR1/est_adjRR2 }
    
    if (contrast_func=='difference') { Theta1=1; Theta2=-1;
      } else if (contrast_func=='ratio') { Theta1=1/est_adjRR2; Theta2=-est_adjRR1/(est_adjRR2)^2 }
    
    est_var_a1 = ( IF_adjustRR(t, ws1, a1, trial_type, 'square_sum', Delta_val)/N1 
                   - (IF_adjustRR(t, ws1, a1, trial_type, 'sum', Delta_val)/N1)^2 )
    est_var_a2 = ( IF_adjustRR(t, ws2, a2, trial_type, 'square_sum', Delta_val)/N2 
                   - (IF_adjustRR(t, ws2, a2, trial_type, 'sum', Delta_val)/N2)^2 )
    
    if (trial_type=='platform') {
      iid_data_a1 = IF_adjustRR(t, ws1, a1, trial_type, 'iid', Delta_val)
      iid_data_a2 = IF_adjustRR(t, ws2, a2, trial_type, 'iid', Delta_val)
      unique_common_times = sort(unique(subset(input_data, W %in% intersect(ws1, ws2) & A==0)$X))
      
      temp = do.call(rbind, lapply(1:length(unique_common_times), function(i){
                time = unique_common_times[i]
                condition1 = (iid_data_a1$a==0 & iid_data_a1$time==time)
                condition2 = (iid_data_a2$a==0 & iid_data_a2$time==time)
        
                iid_vals_a1 = iid_data_a1[condition1, 'iid_vals']
                iid_vals_a2 = iid_data_a2[condition2, 'iid_vals']
                min_length = min(length(iid_vals_a1), length(iid_vals_a2))
        
                return( data.frame( 'time' = time,
                                    'iid_a1' = iid_vals_a1[1:min_length],
                                    'iid_a2' = iid_vals_a2[1:min_length] ) ) }))
      est_cov_a1a2 = sum(temp[,'iid_a1']*temp[,'iid_a2'])/N - (sum(temp[,'iid_a1'])/N)*(sum(temp[,'iid_a2'])/N)
      est_var = ((Theta1)^2)*est_var_a1 + ((Theta2)^2)*est_var_a2 + 2*Theta1*Theta2*est_cov_a1a2
    } else if (trial_type=='separate'){
      probs_separate_trials = Get_probs_separate_trials(input_data)
      sum_probs_separate_trials = sum(probs_separate_trials[,'prob'])
      prob_a1_separate_trial = probs_separate_trials[probs_separate_trials$a==a1,'prob'] 
      prob_a2_separate_trial = probs_separate_trials[probs_separate_trials$a==a2,'prob'] 
      
      est_var = ( (sum_probs_separate_trials/prob_a1_separate_trial)*((Theta1)^2)*est_var_a1 
                   + (sum_probs_separate_trials/prob_a2_separate_trial)*((Theta2)^2)*est_var_a2 ) }
    
    ci = c( est-qnorm(1-alpha)*sqrt(est_var)/sqrt(N), est+qnorm(1-alpha)*sqrt(est_var)/sqrt(N) )
    pval = pnorm(abs(sqrt(N)*est/sqrt(est_var)), mean=0, sd=1, lower.tail=FALSE, log.p=FALSE)
    
    return( list(est=est, sd=sqrt(est_var), se=sqrt(est_var)/sqrt(N), ci=ci, pval=pval) )
  }
 )
)  
    