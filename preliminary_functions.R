num_digits=7

Create_window_intervals_bydays = function(start_date, widths_in_days){
  cum_widths = c(0, cumsum(widths_in_days))
  
  d_lb = do.call(rbind, lapply(cum_widths, function(day_num){
    d = ymd(as.Date(start_date)) + days(day_num)
    return( format(d, format = "%Y-%m-%d") ) }))
  
  d_ub = do.call(rbind, lapply(d_lb[-1], function(d0){
    d = ymd(as.Date(d0)) - days(1)
    return( format(d, format = "%Y-%m-%d") ) }))
  
  intervals = cbind(d_lb[-length(d_lb)], d_ub)
  colnames(intervals) = c('start', 'end')
  return( intervals )
}

Create_window_intervals_bymonth = function(start_date, widths_in_months){
  cum_widths = c(0, cumsum(widths_in_months))
  
  d_lb = do.call(rbind, lapply(cum_widths, function(month_num){
    d = ymd(as.Date(start_date)) %m+% months(month_num)
    return( format(d, format = "%Y-%m") ) }))
  
  d_ub = do.call(rbind, lapply(d_lb[-1], function(d0){
    d = ymd(as.Date(paste0(d0,"-01"))) %m-% months(1)
    return( format(d, format = "%Y-%m") ) }))
  
  intervals = cbind(d_lb[-length(d_lb)], d_ub)
  colnames(intervals) = c('start', 'end')
  return( intervals )
}

Assign_window = function(target_date, window_intervals){
  window_indx = sum( (1:nrow(window_intervals))*(window_intervals[,'start']<=target_date & 
                                                 target_date<=window_intervals[,'end']) )
  return( window_indx )
}

round_df <- function(df, digits) {
  numeric_columns <- sapply(df, mode) == 'numeric'
  df[numeric_columns] <- round(df[numeric_columns], digits)
  return( df )
}

True_StratCumIncidence = function(t, ws, a, z, data){ 
  # ws: a column vector with >=1 element
  data1 = subset(data, W %in% ws & A==a & Z==z, select=c(T0))
  return( mean(data1$T0 <= t) )
}

True_StratRelativeRisk = function(t, ws, a, z, trial_type, data){ 
  # ws: a column vector with >=1 element
  if (trial_type=='platform'){
    data0 = subset(data, W %in% ws & A==0 & Z==z, select=c(T0))
    data1 = subset(data, W %in% ws & A==a & Z==z, select=c(T0))
  } else if (trial_type=='separate'){
    data0 = subset(data, W %in% ws & A==0 & Z==z & Label_A==a, select=c(T0))
    data1 = subset(data, W %in% ws & A==a & Z==z & Label_A==a, select=c(T0)) }
  
  True_F0_z = mean(data0$T0 <= t); True_F1_z = mean(data1$T0 <= t)
  return( True_F1_z/True_F0_z )
}

True_contrast_StratRRs = function(t, ws1, ws2, a1, a2, z, contrast_func, trial_type, data){ 
  # ws1, ws2: column vectors with >=1 element
  srr1 = True_StratRelativeRisk(t, ws1, a1, z, trial_type, data)
  srr2 = True_StratRelativeRisk(t, ws2, a2, z, trial_type, data)
  val = if (contrast_func=='difference') { srr1-srr2 } else if (contrast_func=='ratio') { srr1/srr2 }    
  return( val )
}

True_CumIncidence = function(t, ws, a, data){ # ws: a column vector with >=1 element
  data1 = subset(data, W %in% ws & A==a, select=c(T0))
  return( mean(data1$T0 <= t) )
}

True_RelativeRisk = function(t, ws, a, trial_type, data){ 
  # ws: a column vector with >=1 element
  if (trial_type=='platform'){
    data0 = subset(data, W %in% ws & A==0, select=c(T0))
    data1 = subset(data, W %in% ws & A==a, select=c(T0))
  } else if (trial_type=='separate'){
    data0 = subset(data, W %in% ws & A==0 & Label_A==a, select=c(T0))
    data1 = subset(data, W %in% ws & A==a & Label_A==a, select=c(T0)) }
    
  True_F0 = mean(data0$T0 <= t); True_F1 = mean(data1$T0 <= t)
  return( True_F1/True_F0 )
}

True_contrast_RRs = function(t, ws1, ws2, a1, a2, contrast_func, trial_type, data){ 
  # ws1, ws2: column vectors with >=1 element
  rr1 = True_RelativeRisk(t, ws1, a1, trial_type, data)
  rr2 = True_RelativeRisk(t, ws2, a2, trial_type, data)
  val = if (contrast_func=='difference') { rr1-rr2 } else if (contrast_func=='ratio') { rr1/rr2 }    
  return( val )
}

Get_Wset_treatment = function(a, data){
  return( sort(unique(subset(data, A==a)$W)) )
}

Get_probs_separate_trials = function(data){
  treatments = sort(unique(subset(data, A!=0)$A))
  probs_separate_trials = do.call(rbind, lapply(1:length(treatments), function(i){
    a = treatments[i]; Wset = Get_Wset_treatment(a, data)
    return( data.frame('a'=a, 'prob'=dim(subset(data, W %in% Wset & A %in% c(a,0)))[1]/dim(data)[1]) )}))
  return( probs_separate_trials )
}


##### To compute the event numbers for treatment and placebo groups at different time points
Compute_EventsOfGroups_ByTime = function(t, incidence_rates, data){
  
  temp0 = do.call(rbind, lapply(1:k, function(i, t0=t){ 
    
    windows_taken = which(incidence_rates[i+1,]>0)
    temp_df0 = data[data$A %in% c(0,i) & data$W %in% windows_taken,] 
    temp_df1 = temp_df0[!is.na(temp_df0$Event_CalendarTime),]
    temp = temp_df1[temp_df1$Event_CalendarTime<=t0,]
    
    result_df = data.frame('group'=paste0(0, ' v.s ', as.character(i)), 
                           'size'=dim(temp_df0)[1], 
                           'size_p'=dim(temp_df0[temp_df0$A==0,])[1],
                           'size_i'=dim(temp_df0[temp_df0$A==i,])[1],
                           'uncensored_num_p'=dim(temp_df1[temp_df1$A==0,])[1],
                           'uncensored_num_i'=dim(temp_df1[temp_df1$A==i,])[1],
                           'events_p'=dim(temp[temp$A==0,])[1], 
                           'events_i'=dim(temp[temp$A==i,])[1],
                           'events'=dim(temp)[1])
    
    colnames(result_df) = c('group','size','size_p','size_i','uncensored_num_p',
                            'uncensored_num_i', paste0('events_p_',t0), 
                            paste0('events_i_',t0), paste0('events_',t0))
    
    return (result_df) } )) 
}


##### To find the time when the event number start >=150 per group
Get_Time_ToEventNums = function(num, incidence_rates, data, ts){ 
  
  temp_events = do.call(cbind, lapply(1:length(ts), function(i){ 
    temp_df = Compute_EventsOfGroups_ByTime(t=ts[i], incidence_rates, data=data)
    temp_df[c('group', paste0('events_',ts[i]))] }))
  
  temp_events = temp_events[, !duplicated(colnames(temp_events))]
  
  result_df = do.call(rbind, lapply(1:dim(temp_events)[1], function(row){
    month_string = colnames(temp_events[min(which(temp_events[row,-1]>=num))]) 
    data.frame('group'=temp_events[row, 'group'],
               'month_to_num'=as.numeric(sub(".*_", "", month_string))) }))
  
  return (result_df) 
}
