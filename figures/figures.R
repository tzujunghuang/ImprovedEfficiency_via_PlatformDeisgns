
package_list = c('ggplot2', 'latex2exp', 'grid', 'tidyr', 'dplyr')
for (package in package_list){
  require(package, character.only=TRUE); library(package, character.only=TRUE) }

wk_dir = 'D:/Covid19_vaccine efficacy/public R codes/figures/' #where files are saved
setwd(wk_dir)
to_dir = wk_dir

gg_color_hue = function(n) { hues=seq(15, 375, length=n + 1); hcl(h=hues, l=65, c=100)[1:n] }


##### Main Text
### Figure 1: Relative risk (RR) constancy over windows
# Generate various paris of event rates for arms, given a value of RR
data_example = data.frame(intervention_label=rep(c('Intervention 1', 'Intervention 2'), each=10),
                          ws=rep(c('1', '2', '3', '4', '5'), 4),
                          arm_label=rep(rep(c('Control', 'Intervention'), each=5), 2),
                          ars=c(10, 15, 5, NA, 5, 4, 6, 2, NA, 2,
                                NA, 15, 7.5, 10, NA, NA, 9, 4.5, 6, NA))
data_example$ws = factor(data_example$ws, levels=c('1', '2', '3', '4', '5'))
data_RR = data.frame(intervention_label=c('Intervention 1', 'Intervention 2'), RRs=c('RR=0.4', 'RR=0.6'))

plt = ggplot(data_example, aes(x=ws, y=ars, group=arm_label, label=as.character(ars))) +
      geom_point(aes(color=arm_label, shape=arm_label), size=5) + 
      geom_label(data=data_RR, aes(group=intervention_label, label=RRs), x=Inf, y=Inf, hjust=0.9, 
                 vjust=0.9, label.padding=unit(0.5, 'lines'), size=7) +         
      scale_x_discrete('Window of Enrollment', labels = parse(text=levels(data_example$ws))) + 
      scale_y_continuous('Event Rate (%)', limits=c(0, 16.5)) + 
      geom_text(vjust=0.2, nudge_x=0.4, nudge_y=0.01, check_overlap=TRUE, size=7) + 
      scale_shape_manual('Six-month\nAttack Rates', values=c(17, 15)) + 
      scale_color_manual( values=c(gg_color_hue( length(unique(data_example$arm_label)) )) ) + 
      facet_grid(intervention_label~ .) + 
      labs(y='Event Rate (%)', x='Window of Enrollment', color='Six-month\nAttack Rates') + 
      theme_bw() + 
      theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
            axis.ticks.length = unit(2, 'pt'), 
            axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), 
            legend.text = element_text(size=22), legend.title = element_blank(),
            strip.text.y = element_text(size=20))


### Figure 2: Ratio of confidence interval widths
meth_vecs = c('cov_stratRR', 'cov_adjRR')
select_vars_titles = list('mse_ratios'=c('Ratio of Mean Squared Errors', 'MSEratios'), 
                          'ci_width_ratios'=c('Ratio of Confidence Interval Widths', 'CIWidthratios'))
varname = 'ci_width_ratios'
title_text = select_vars_titles[[varname]][1]
filename_text = select_vars_titles[[varname]][2]

d1 = get(load('ratios_efficiencygain_ltfr10.Rdata')) %>%
     mutate(meth_name = factor(meth, levels=c('cov_stratRR', 'cov_stratRR_nullZ', 'cov_adjRR') , 
                               labels=c('Conditional Relative Risk Ratio', 'Unadjusted Relative Risk Ratio', 
                                        'Marginal Relative Risk Ratio'))) %>%
     filter(as.character(meth) %in% meth_vecs)

plt = ggplot(d1, aes(x=d1[,varname])) + 
      geom_vline(xintercept=1, linetype=2, size=0.8) + 
      geom_histogram(binwidth=0.05, alpha=.8, position='identity', color='black', fill='cornflower blue') +
      labs(x=title_text, y='Count') + theme_bw() + facet_grid(~meth_name, scales='free') +
      coord_cartesian(xlim=c(0.68,1.02), ylim=c(0, 30), expand=TRUE, clip='off') +
      theme(axis.text.x = element_text(size=20), axis.text.y = element_text(size=20),
            axis.ticks.length = unit(2, 'pt'),
            axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), 
            legend.position = '', legend.text = element_text(), legend.title = element_blank(),
            panel.spacing.x = unit(0.2, 'cm'),
            strip.text.x = element_text(size=22), strip.text.y = element_text(size=22))


### Figure 3: Empirical rejeciton rates of the intersection test and the likelihood ratio test
# Get true parameters
true_pars = get(load('true_pars_efficacycompetition.Rdata'))

# Obtain the boundry in margin between the null and complementary alternatives
diff_RRs = aggregate(true_pars[,'RRdiff'], list(true_pars$t, true_pars$a1), FUN=max)
names(diff_RRs) = c('t', 'pre_trt', 'diff_RRs')
margins = seq(0, 0.4, 0.1)
for (m in margins) { col_name = paste0('diff_vs_', m)
                     diff_RRs[ ,col_name] = m*(diff_RRs[,'diff_RRs'] >= m) }
diff_RRs$max_margin_null = as.character(apply(diff_RRs[,grepl('diff_vs_', names(diff_RRs))], 1, max))
diff_RRs = diff_RRs[diff_RRs$pre_trt!=3,-(c(1:dim(diff_RRs)[2])*grepl('diff_vs_', names(diff_RRs)))]

# Prepare data for plot 
acc_time_val = 6; delta_val = 0.7; pre_trts = c(7,9); 
meth_vecs = c('intersection_test','likelihood_ratio_test')

d1 = get(load('stats_efficacycompetition_ltfr10_alpha025_comparison.Rdata')) %>%
     filter(acc_time==acc_time_val & tolerance==delta_val & pre_trt%in%pre_trts) %>%
     left_join(., diff_RRs, by=c('t', 'pre_trt')) %>%
     mutate(type_name = factor(trial_type, levels=c('platform', 'separate'), labels=c('platform', 'separate')),
            marginf = factor(margin, levels=c(0, 0.1, 0.2, 0.3, 0.4), labels=c('0','0.1','0.2','0.3','0.4')),
            acc_time_f = factor(acc_time, levels=c(6,18), labels=c('censored at t=6', 'censored at t=18')))

d1 = d1[with(d1, order(t, pre_trt, trial_type, acc_time, margin, tolerance)),]

d1_l = tidyr::gather(d1, test, rej, rej_intersection, rej_oracle, rej_LRT, factor_key=TRUE) %>%
       mutate(test = gsub('rej_', '', test),
              trial_test = as.factor(paste0('(',type_name,', ',test,')'))) %>%
       filter(test!='oracle')
  
# Start plot   
t.labs = c("t = 3", "t = 6")
names(t.labs) = c(3, 6)

pre_trt.labs = c("Intervention 7", "Intervention 9")
names(pre_trt.labs) = c(7, 9)

ann_arrows = data.frame(xmid=2, xmin=1, xmax=3, y=0.6, ymin=0.57, ymax=0.63, t=factor(3, levels=c(3,6)), 
                        pre_trt=factor(7, levels=c(3,5,7,9)))
ann_texts = data.frame(x=c(1.5,2.9), y=c(0.52,0.7), t=factor(3, levels=c(3,6)), 
                       pre_trt=factor(7, levels=c(3,5,7,9)), label=c("null", "alternative"))

plt = ggplot(d1_l, aes(x=marginf, y=rej, colour=type_name, shape=test, group=trial_test)) +
      geom_line(size=0.6) + geom_point(size=1.5) + scale_shape_manual(values=c(15, 17)) + 
      facet_grid(pre_trt~t, scales='free', labeller = labeller(pre_trt=pre_trt.labs, t=t.labs)) + 
      geom_hline(yintercept=0.025, linetype=2, size=0.4) + 
      geom_vline(aes(xintercept=max_margin_null), linetype=6, size=0.4) + 
      geom_segment(data=ann_arrows, aes(x=xmid,xend=xmin,y=ymin,yend=ymin), size=0.6, 
                   arrow=arrow(length=unit(0.2,"cm")), show.legend=FALSE, inherit.aes=FALSE) + 
      geom_point(data=ann_arrows, mapping=aes(x=xmid, y=ymin), size=1.5, shape=19, inherit.aes=FALSE) +
      geom_segment(data=ann_arrows, aes(x=xmid,xend=xmax,y=ymax,yend=ymax), size=0.6, 
                   arrow=arrow(length=unit(0.2,"cm")), show.legend=FALSE, inherit.aes=FALSE) +
      geom_point(data=ann_arrows, aes(x=xmid, y=ymax), size=1.5, shape=1, inherit.aes=FALSE) +
      geom_text(data=ann_texts,aes(x=x, y=y, label=label), size=4, show.legend=FALSE, inherit.aes=FALSE) +
      scale_x_discrete('Margin', labels = parse(text=levels(d1_l$marginf)) ) +
      scale_y_continuous(limits=c(0, 1)) +
      labs(y='Empirical Rejection Rate', x='Margin') + theme_bw() +
      guides(shape = guide_legend(override.aes = list(size=1.7))) +
      theme(axis.text.x = element_text(size=14), axis.text.y = element_text(size=14),
            axis.ticks.length = unit(2, 'pt'), axis.title.x = element_text(size=18), 
            axis.title.y = element_text(size=18), legend.position = 'bottom', 
            legend.text = element_text(size=16), legend.title = element_blank(),
            panel.spacing.x = unit(0.2, 'lines'), panel.spacing.y = unit(0.2, 'lines'),  
            strip.text.x = element_text(size=16), strip.text.y = element_text(size=16))


### Figure 4: Confidence width vs. Proportion of controls shared in HVTN 703 and 704
result = get(load('efficacycompetition_cis_AMPdata_stratified_low_fixed500.Rdata')) %>%
         filter(prob!=0.5)
Z.labs = c('HVTN 703', 'HVTN 704'); names(Z.labs) = c('HVTN 703', 'HVTN 704')

plt = ggplot( data=result, aes(x=sharing_rate, y=ci_width, group=Z)) + 
      geom_line(color='blue') + geom_point(size=2, colour='blue') + 
      facet_grid(~Z, scales='free', labeller = labeller(Z=Z.labs)) + 
      scale_x_continuous('Proportion of Controls Shared', breaks=c(0.2,0.3,0.4,0.5), limits=c(0.2,0.51)) + 
      scale_y_continuous('Confidence Interval Width', breaks=c(0.25,0.3,0.35), limits=c(0.25, 0.35)) + 
      theme_bw() + 
      theme(axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=22), axis.title.y = element_text(size=22), 
            legend.position = 'none', legend.text = element_text(), panel.spacing.x = unit(0.6, 'cm'),
            strip.text.x = element_text(size=20))



##### Supplementary Materials 
meth_vecs = c('cov_stratRR', 'cov_adjRR')
select_vars_titles = list('mse_ratios'=c('Ratio of Mean Squared Errors', 'MSEratios'), 
                          'ci_width_ratios'=c('Ratio of Confidence Interval Widths', 'CIWidthratios'))
varname = 'ci_width_ratios'
title_text = select_vars_titles[[varname]][1]
filename_text = select_vars_titles[[varname]][2]

d1 = get(load('ratios_efficiencygain_ltfr10.Rdata')) %>%  filter(as.character(meth) %in% meth_vecs) %>%
     mutate(meth_name = factor(meth, levels=c('cov_stratRR', 'cov_stratRR_nullZ', 'cov_adjRR') , 
                               labels=c('Conditional Relative Risk Ratio', 'Unadjusted Relative Risk Ratio', 
                                        'Marginal Relative Risk Ratio')),
            acc_time_f = factor(acc_time, levels=c(6,9,12,18), labels=c('censored at 6', 'censored at 9', 
                                                                        'censored at 12', 'censored at 18')))

acc_time.labs = c("censored at 6", "censored at 9", "censored at 12", "censored at 18")
names(acc_time.labs) = c(6, 9, 12, 18)
t.labs = c("t = 3", "t = 6")
names(t.labs) = c(3,6)

### Figure 1: Ratio of confidence interval widths for conditional relative risk ratio - Full plot
d2 = d1 %>% filter(meth=='cov_stratRR')
plt = ggplot(d2, aes(x=d2[,varname])) + geom_vline(xintercept=1, linetype=2, size=0.8) +
      geom_histogram(binwidth=0.05, alpha=.8, position="identity", color='black', fill='cornflower blue') +
      facet_grid(acc_time~t, scales='free', labeller = labeller(acc_time=acc_time.labs, t=t.labs)) +
      scale_y_continuous(limits=c(0, 8)) + labs(x=title_text, y='Count') + theme_bw() +
      theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), 
            axis.ticks.length = unit(2, 'pt'),
            axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
            legend.position = 'none', 
            panel.spacing.x = unit(0.2, 'lines'), panel.spacing.y = unit(0.2, 'lines'),
            strip.text.x = element_text(size=14), strip.text.y = element_text(size=12))

 
### Figure 2: Ratio of confidence interval widths for marginal relative risk ratio - Full plot
d2 = d1 %>% filter(meth=='cov_adjRR')
plt = ggplot(d2, aes(x=d2[,varname])) + geom_vline(xintercept=1, linetype=2, size=0.8) +
      geom_histogram(binwidth=0.05, alpha=.8, position='identity', color='black', fill='cornflower blue') +
      facet_grid(acc_time~t, scales='free', labeller = labeller(acc_time=acc_time.labs, t=t.labs)) +
      scale_y_continuous(limits=c(0, 8)) + labs(x=title_text, y='Count') + theme_bw() +
      theme(axis.text.x = element_text(size=13), axis.text.y = element_text(size=13), 
            axis.ticks.length = unit(2, 'pt'),
            axis.title.x = element_text(size=20), axis.title.y = element_text(size=20),
            legend.position = 'none', 
            panel.spacing.x = unit(0.2, 'lines'), panel.spacing.y = unit(0.2, 'lines'),
            strip.text.x = element_text(size=14), strip.text.y = element_text(size=12))


### Figure 3: The numbers of participants enrolled in various groups and windows across trials
p.labs = c('HVTN 703', 'HVTN 704'); names(p.labs) = c('HVTN 703', 'HVTN 704')
gp_freq = get(load('frequency_across_windows.Rdata'))
plt = ggplot(gp_freq, aes(x=W, y=freq, fill=txf, label=freq, group=txf)) + 
      geom_bar(stat='identity') + geom_text(size=5, position=position_stack(vjust=0.5)) + 
      guides(alpha='none') + labs(x='Window', y='Frequency', title='', fill='Group') + 
      theme_bw() + facet_grid(.~protocol, scales='free', labeller = labeller(protocol=p.labs)) +
      theme(plot.title = element_text(hjust=0.5), 
            axis.text.x = element_text(size=16), axis.text.y = element_text(size=16),
            axis.ticks.length = unit(2, 'pt'), axis.title.x = element_text(size=22),
            axis.title.y = element_text(size=22), legend.text = element_text(size=16), 
            legend.title = element_text(size=18), panel.spacing.y = unit(0.2, "lines"), 
            strip.text.x = element_text(size=18))
