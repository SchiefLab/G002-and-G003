library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure S29 comparisons
DIR = file.path(SOURCEDIR, 'Sup-Metrics/S25')
STATS = './Stats.codebase'

# fn = list.files(DIR)
fn = c(
  "figA_prime_v_heavy_nt_mutation.csv",
  "figB_prime_v_light_nt_mutation.csv"
)

dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n), check.names = FALSE)

  if( 'cottrell_focused_v_common_score' %in% names(df) ) {
    stopifnot(all(df$cottrell_focused_v_common_score == df$residues))
    df$residues = NULL
  }
  stopifnot(ncol(df)==5)
  stopifnot(all(names(df)[1:4] == c('pubID', 'trial', 'pseudogroup', 'weeks')))

  oc = names(df)[5]
  df$y_label = oc
  df$y = df[[oc]]
  df[[oc]] = NULL

  df$panel = substr(n, 4, 4)

  return(df)
})
dat = do.call(rbind, dat)


# Comparison table,

stat = c()

for( p in unique(dat$panel) ) {
  for( wk in c(4, 8, 16) ) {
    sss = subset(dat, panel==p & weeks==wk & trial %in% c('G001', 'G002'))

    tmp = pairwise_test_cont(
      x = sss$y, group = sss$trial,
      method = 'wilcox', paired = FALSE,
      alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
      verbose = FALSE)
    tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
    names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

    tmp$panel = p
    tmp$weeks = wk
    tmp$comp = 'G001 vs. G002'
    stat = rbind(stat, tmp)
  }
}

stat = stat %>% arrange(panel, weeks)

# add FDR adjustment
stat = stat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()


# across time comparisons
# Panels A and B
# G001 wk8 vs wk16
# G002 wk8 vs wk16
# G002 wk16 vs wk24
# G003 wk8 vs wk16
# G003 wk16 vs wk21
comp = data.frame(
  trial = c('G001', 'G002', 'G002', 'G003', 'G003'),
  wk1 = c(8, 8, 16, 8, 16),
  wk2 = c(16, 16, 24, 16, 21)
)
tstat = c()

for( p in c('A','B') ) {
  for( i in 1:nrow(comp) ) {
    sss = subset(dat, panel==p & weeks %in% c(comp$wk1[i], comp$wk2[i]) & trial == comp$trial[i])
    sss$wk = factor(sss$weeks,
                    levels = c(comp$wk1[i], comp$wk2[i]),
                    labels = paste0('wk', c(comp$wk1[i], comp$wk2[i])))

    tmp = pairwise_test_cont(
      x = sss$y, group = sss$wk, id = sss$pubID,
      method = 'wilcox', paired = TRUE,
      alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
      verbose = FALSE)
    tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
    names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

    tmp$panel = p
    tmp$trial = comp$trial[i]
    tstat = rbind(tstat, tmp)
  }
}

# add FDR adjustment
tstat = tstat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()


write.csv(stat, file.path(STATS, 'fig3_alt_pvals.csv'), row.names = FALSE)
write.csv(tstat, file.path(STATS, 'fig3_alt_time_pvals.csv'), row.names = FALSE)

q(save='no')
