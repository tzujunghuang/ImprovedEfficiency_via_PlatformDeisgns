source('preliminary_functions.R')

##### Function to generate the rate parameters of the piecewise exponential distribution 
Get_piecewise_rates = function(complete_incidence_rates, num_interventions, num_windows, window_bands){
  k = num_interventions; m = num_windows
  rates0 = do.call(rbind, lapply(0:k, function(i){ data.frame('A'=rep(i, times=m), 'W'=1:m, 
                                                              'rate'=complete_incidence_rates[i+1,],
                                                              'tchange'= window_bands$lower_end) }))
  
  rates = do.call(rbind, lapply(0:k, function(i){ 
    rate_vec0 = do.call(rbind, lapply(1:m, function(j){ list(rates0[rates0$A==i,'rate']) }))
    tchange_vec0 = do.call(rbind, lapply(1:m, function(j){ list(rates0[rates0$A==i,'tchange']) }))
    data.frame('A'=rep(i, times=m), 'W'=1:m, 'rate_vec'=rate_vec0, 'tchange_vec'=tchange_vec0) }))
  
  return (rates)
}

##### Class to format the required parameters and to simulate survival data
Generate_SurvivalData <- setRefClass( "Generate_SurvivalData",
  fields = list(
    treatment_labels_in_windows = "list",
    incidence_rates = "matrix",
    window_bands = "data.frame",
    k = "numeric",
    m = "numeric",
    entry_info_by_windows = "data.frame",
    rates = "data.frame",
    n_m = "numeric",
    data_p = "data.frame",
    n_p = "numeric",
    data_s = "data.frame",
    n_s = "numeric",
    p = "numeric"
  ),
                               
  methods = list(
    initialize = function(treatment_labels_in_windows = input_treatment_labels_in_windows,
                          incidence_rates = input_incidence_rates, window_bands = input_window_bands, 
                          n_m = input_n_m, p=1){
      treatment_labels_in_windows <<- treatment_labels_in_windows
      incidence_rates <<- incidence_rates
      window_bands <<- window_bands
      k <<- dim(incidence_rates)[1]-1
      m <<- dim(incidence_rates)[2]
      entry_info_by_windows <<- do.call(rbind, lapply(1:m, function(i){ 
        data.frame('W'=i, 'entry_pt'=window_bands[i,'lower_end'], 
                   'enroll_length'=window_bands[i,'length']) }))
      rates <<- Get_piecewise_rates(input_incidence_rates, k, m, input_window_bands)
      n_m <<- n_m
      data_p <<- do.call(rbind, lapply(1:m, function(i) { 
        data.frame('A'=rep(c(0, treatment_labels_in_windows[[i]]), each=n_m[i]),
                   'W'=rep(i, times=(length(treatment_labels_in_windows[[i]])+1)*n_m[i])) }))
      n_p <<- dim(data_p)[1]
      data_s <<- do.call(rbind, lapply(1:m, function(i) { 
        data.frame('A'=c(rep(0, length(treatment_labels_in_windows[[i]])*n_m[i]), 
                         rep(treatment_labels_in_windows[[i]], each=n_m[i])),
                   'W'=rep(i, times=2*length(treatment_labels_in_windows[[i]])*n_m[i]),
                   'Label_A'=c(rep(treatment_labels_in_windows[[i]], each=n_m[i]),
                               rep(treatment_labels_in_windows[[i]], each=n_m[i])))  }))
      n_s <<- dim(data_s)[1]
      p <<- p
  },

  
##### Simulate data
  Generate_data = function(trial_type, loss_to_follow_rate, adm_censoring_calendartime, true_data=FALSE, 
                           loss_to_follow_pars=c(600, 120)){
    # if true_data=TRUE, loss_to_follow_rate=NA; adm_censoring_calendartime=NA; loss_to_follow_pars=NA
    
    if (trial_type=='platform'){ data1 = data_p } else if (trial_type=='separate'){ data1 = data_s }  
    n = dim(data1)[1]
    data1 = merge(data1, entry_info_by_windows, by=c('W'))
    data1 = data1[!duplicated(as.list(data1))]; colnames(data1) = gsub(".x", "", names(data1), fixed = TRUE)
    data1$Entry_CalendarTime = data1$entry_pt + runif(n, 0, data1$enroll_length)
      
    data1$Z = sapply(data1$W, function(w){ if(w %in% c(1,2)){ U=runif(1,0,0.4); 
                                                              z=(U<=0.04)+(U<=0.12)+(U<=0.24)+(U<=0.4)
                                } else{ U=runif(1,0.4,1); z=(U<=0.46)+(U<=0.58)+(U<=0.76)+(U<=1) } })
    
    AW_pairs = unique(data1[,c('A','W')])
    data1$T00 = do.call(rbind, lapply(1:dim(AW_pairs)[1], function(i) { 
                         a = AW_pairs[i,'A']; w = AW_pairs[i,'W']; nr = dim(data1[data1$A==a & data1$W==w,])[1] 
                         t0 = rpwe(nr=nr, rate=unlist(rates[rates$A==a & rates$W==w, 'rate_vec']), 
                                   tchange=unlist(rates[rates$A==a & rates$W==w, 'tchange_vec']))[[1]]
                         matrix(round(t0, num_digits), ncol=1) } ))
    scale_num = ifelse(data1$Z %in% c(3,4), 0.5, 1)
    data1$T0 = data1$T00*scale_num
    
    if (!(true_data)){
      for (i in 1:length(loss_to_follow_pars)) { assign( paste('par', i, sep=''), loss_to_follow_pars[i] ) }
      ltf_par = par1*(loss_to_follow_rate=='2%') + par2*(loss_to_follow_rate=='10%')
      
      data1$C0 = runif(n, 0, ltf_par)
      data1$C = pmin(adm_censoring_calendartime-data1$Entry_CalendarTime, data1$C0)
      data1$X = round(pmin(data1$T0, data1$C), num_digits);
      data1$Delta = 1*(data1$T0 <= data1$C);
      
      data1$Event_CalendarTime = data1$Entry_CalendarTime + ifelse(data1$Delta==1, data1$T0, NA)
      data1$Observed_CalendarTime = data1$Entry_CalendarTime + data1$X
    } else{ data1$Event_CalendarTime0 = data1$Entry_CalendarTime + data1$T0 }
      
    return( data1[complete.cases(data1[c('A', 'W', 'Z', 'T0')]),] )
  }
 )
)