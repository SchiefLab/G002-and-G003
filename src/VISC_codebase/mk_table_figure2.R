library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

DIR = file.path(SOURCEDIR, 'Main-Metrics/fig2')
STATS = './Stats.codebase'

### Get the figure 2 data and combine

# The data for each panel is in a separate file.  The list `map` maps
# the panel label (A, B, C, etc.) to the file name
map = list(
  a = "figA_percent_gt8$^{++}$_amongigg$^+$_b_cells.csv",
  b = "figB_percent_cd4bs-specific_amongigg$^{+}$_b_cells.csv" ,
  c = "figC_percent_of_gt8$^{++}$igg$^{+}$b_cells_that_are_ko$^{-}$.csv",
  d = "figD_number_of_vrc01-classigg$^{+}$_b_cells.csv",
  e = "figE_percent_of_igg$^{+}$_memory_b_cellsdetected_as_vrc01-class.csv",
  f = "figF_percent_vrc01-class_amongcd4bs-specificigg$^{+}$_memory_sequences.csv",
  g = "figG_percent_vrc01-class_amonggt8-specificigg$^{+}$_memory_sequences.csv",
  h = "figH_percent_of_vrc01-class_responders.csv"
)

# read in the data files
dat = lapply(names(map), function(panel) {
  read.csv(file.path(DIR, map[[panel]]), check.names = FALSE)
})
names(dat) = names(map)

# process each file so that we can rbind them
# the data are trial, pubID, weeks and y-axis label name
# convert to a format where I can rbind the files and use
# outcome to process the comparisons.

# NOTE, panel H has a different column order
dat = lapply( names(dat), function(panel) {
  dat_panel = dat[[panel]]
  if(panel=='h') {
    stopifnot(all(names(dat_panel) == c("pubID", "trial", "percent_of_vrc01-class_responders", "weeks")))
    dat_panel = dat_panel[,c("trial", "pubID", "weeks", "percent_of_vrc01-class_responders")]
  }
  stopifnot(all(names(dat_panel)[1:3] == c('trial', 'pubID', 'weeks')))
  oc = names(dat_panel)[4]
  dat_panel$panel = toupper(panel)
  dat_panel$y_label = oc
  dat_panel$y = dat_panel[[oc]]
  dat_panel[[oc]] = NULL

  return(dat_panel)
})

# rbind the data files
dat = do.call(rbind, dat)

# response call CIs panels A and B
thold = subset(dat, weeks == -5 & panel %in% c('A', 'B')) %>%
  group_by(trial, panel) %>%
  summarize(thold = mean(y) + 3*sd(y), .groups = 'drop')

out = dat %>%
  filter( panel %in% c('A', 'B') ) %>%
  merge(thold) %>%
  mutate(Response = as.numeric(y > thold)) %>%
  group_by(trial, panel, weeks) %>%
  summarize(
    binom_ci(Response),
    x = sum(Response==1),
    n = n(),
    sum_info = paste0(x, ' of ', n, ' = ',
                      round_away_0(100*mean, 1, trailing_zeros = TRUE), '\\% (',
                      round_away_0(100*lower, 1, trailing_zeros = TRUE),'\\%, ',
                      round_away_0(100*upper, 1, trailing_zeros = TRUE),'\\%)'),
    `.groups` = "drop")
out = merge(out, thold) # add back the thresolds applied

write.csv(out, file.path(STATS,'G00X_response_call_panelsAB.csv'), row.names = FALSE)

# remove out so not inadvertently used below
out = NULL


# magnitude comparisons of baseline vs. post-vaccination
out = c()
for( t in c('G001','G002','G003') ) {
  for( p in c('A', 'B', 'C', 'E') ) {
    ss = subset(dat, trial==t & panel==p)
    wks = sort(unique(ss$weeks))
    wks = wks[wks > 0]

    for( w in wks ) {
      sss = subset(ss, weeks %in% c(-5, w))
      sss$timepoint = paste0('Wk', sss$weeks)
      tmp = pairwise_test_cont(
        x = sss$y, group = sss$timepoint, id = sss$pubID,
        method = 'wilcox', paired = TRUE,
        alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
        verbose = FALSE)
      tmp$trial = t
      tmp$panel = p
      out = rbind(out, tmp)

    }
  }
}

# save results
write.csv(out, file.path(STATS,'G00X_mag_baseline_v_post.csv'), row.names = FALSE)

# remove out, tmp and indexing variables so they aren't inadvertently used below
out = tmp = t = p = w = NULL



# magnitude comparisons of post-vaccination time points
# requested comparisons
# wk8 vs wk16 for G001
# wk8 vs wk16 for G002
# wk8 vs wk24 for G002
# wk8 vs wk10 for G003
# wk8 vs wk16 for G003
# wk8 vs wk21 for G003
comp = data.frame(t1=rep(8,6),
                  t2 = c(16, 16, 24, 10, 16, 21),
                  trial=c('G001','G002', 'G002', 'G003', 'G003', 'G003'))
comp$order = 1:nrow(comp)

out = c()

for( p in c('A', 'B', 'C', 'E') ) {
  for( i in 1:nrow(comp) ) {
    sss = subset(dat, trial==comp$trial[i] & weeks %in% c(comp$t1[i], comp$t2[i]) & panel==p)

    sss$timept = factor(sss$weeks,
                        levels = c(comp$t1[i], comp$t2[i]),
                        labels = paste0('Wk', c(comp$t1[i], comp$t2[i])))
    tmp = pairwise_test_cont(
      x = sss$y, group = sss$timept, id = sss$pubID,
      method = 'wilcox', paired = TRUE,
      alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
      verbose = FALSE)
    tmp$trial = comp$trial[i]
    tmp$panel = p
    tmp$order = comp$order[i]
    out = rbind(out, tmp)
  }
}

# save results
write.csv(out, file.path(STATS,'G00X_mag_postbaseline_v_postbaseline.csv'), row.names = FALSE)

# remove out, tmp and indexing variables so they aren't inadvertently used below
out = tmp = comp = p = i = NULL

# Fig 2F magnitude testing postbaseline vs. postbaseline
# wk8 vs wk16 for G001
# wk4 vs wk24 for G002 --- special case
# wk8 vs wk16 for G002
# wk8 vs wk24 for G002
# wk8 vs wk10 for G003
# wk8 vs wk16 for G003
# wk8 vs wk21 for G003
comp = data.frame(t1=c(8, 4, rep(8,5)),
                  t2 = c(16, 24, 16, 24, 10, 16, 21),
                  trial=c('G001','G002', 'G002', 'G002', 'G003', 'G003', 'G003'))
comp$order = 1:nrow(comp)

out = c()

for( p in c('F', 'G') ) {
  for( i in 1:nrow(comp) ) {
    sss = subset(dat, trial==comp$trial[i] & weeks %in% c(comp$t1[i], comp$t2[i]) & panel==p)

    sss$timept = factor(sss$weeks,
                        levels = c(comp$t1[i], comp$t2[i]),
                        labels = paste0('Wk', c(comp$t1[i], comp$t2[i])))
    tmp = pairwise_test_cont(
      x = sss$y, group = sss$timept, id = sss$pubID,
      method = 'wilcox', paired = TRUE,
      alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
      verbose = FALSE)
    tmp$trial = comp$trial[i]
    tmp$panel = p
    tmp$order = comp$order[i]
    out = rbind(out, tmp)
  }
}

# save results
write.csv(out, file.path(STATS,'G00X_mag_postbaseline_v_postbaseline_panelFG.csv'), row.names = FALSE)

# remove out, tmp and indexing variables so they aren't inadvertently used below
out = tmp = comp = p = i = NULL


# detectable response call CIs based on panel D
out = dat %>%
  filter( panel == 'D' ) %>%
  mutate(Response = as.numeric( y > 0)) %>%
  group_by(trial, weeks) %>%
  summarize(
    binom_ci(Response),
    x = sum(Response==1),
    n = n(),
    sum_info = paste0(x, ' of ', n, ' = ',
                      round_away_0(100*mean, 1, trailing_zeros = TRUE), '\\% (',
                      round_away_0(100*lower, 1, trailing_zeros = TRUE),'\\%, ',
                      round_away_0(100*upper, 1, trailing_zeros = TRUE),'\\%)'),
    `.groups` = "drop")

# save results
write.csv(out, file.path(STATS,'G00X_detectable_response_panelD.csv'), row.names = FALSE)

# remove out so not inadvertently used below
out = NULL

# comparisons of post-vaccination time points

out = c()

for( t in c('G001', 'G002', 'G003') ) {
  ss = subset(dat, trial==t & panel=='D')
  ss$Response = as.numeric( ss$y> 0 )
  wks = sort(unique(ss$weeks))
  wks = wks[wks > 0]

  for( w in wks ) {
    sss = subset(ss, weeks %in% c(-5, w))
    sss$timepoint = paste0('Wk', sss$weeks)
    tmp = pairwise_test_bin(
      x = sss$Response, group = sss$timepoint, id = sss$pubID,
      method = 'mcnemar',
      alternative = 'two.sided', num_needed_for_test = 3, digits = 1,
      latex_output = TRUE, verbose = FALSE
    )
    tmp$trial = t
    out = rbind(out, tmp)
  }
}

# save results
write.csv(out, file.path(STATS,'G00X_detectable_response_baseline_v_post.csv'), row.names = FALSE)

# remove out, tmp, wks and indexing variables so they aren't inadvertently used below
out = tmp = w = wks = NULL



q(save='no')
