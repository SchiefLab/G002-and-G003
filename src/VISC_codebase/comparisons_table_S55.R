library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

DIR3 = file.path(SOURCEDIR, 'Main-Metrics/fig3')
STATS = './Stats.codebase'


# The data for each panel is in a separate file.  The list `map` maps
# the panel label (A, B, etc.) to the file name
# the different methods for the 90th percentile computation are in
# different directories (e.g. fig3-nearest)
map3 = list(
  a = "figA_prime_v_heavy_aa_mutation.csv",
  b = "figB_prime_v_light_aa_mutation.csv",
  e = "figE_90th_percentile_number_of_key_VRC01-class_HC_residues.csv"
)


# Figure 3 panels A, B and E
# G001 wk8 vs wk16
# G002 wk8 vs wk16
# G002 wk8 vs wk24
# G002 wk16 vs wk24
# G003 wk8 vs wk16
# G003 wk8 vs wk21
# G003 wk16 vs wk21
comp = data.frame(t1 = c(8, 8, 8, 16, 8, 8, 16),
                  t2 = c(16, 16, 24, 24, 16, 21, 21),
                  trial=c('G001','G002', 'G002', 'G002', 'G003', 'G003', 'G003'))
comp$order = 1:nrow(comp)

out = c()


for( SUFF in c('') ) {
  METHOD = ifelse(SUFF=='', 'nearest', sub('^-', '', SUFF))

  # read in the data files
  dat3 = lapply(names(map3), function(panel) {
    D = paste0(DIR3, SUFF)
    read.csv(file.path(D, map3[[panel]]), check.names = FALSE)
  })
  names(dat3) = names(map3)

  # process each file so that we can rbind them
  # the data are trial, pubID, weeks and y-axis label name
  # convert to a format where I can rbind the files and use
  # outcome to process the comparisons.
  dat3 = lapply( names(dat3), function(panel) {
    dat_panel = dat3[[panel]]
    stopifnot(all(names(dat_panel)[1:4] == c('pubID', 'trial', 'pseudogroup', 'weeks')))
    oc = names(dat_panel)[5]
    if( oc == 'cottrell_focused_v_common_score' ) {
      # remove duplicated column
      stopifnot(all(dat_panel$cottrell_focused_v_common_score == dat_panel$residues))
      dat_panel$residues = NULL
    }
    dat_panel$panel = toupper(panel)
    dat_panel$y_label = oc
    dat_panel$y = dat_panel[[oc]]
    dat_panel[[oc]] = NULL

    return(dat_panel)
  })

  # rbind the data files
  dat3 = do.call(rbind, dat3)

  for( p in c('A', 'B', 'E') ) {
    for( i in 1:nrow(comp) ) {
      sss = subset(dat3, trial==comp$trial[i] & weeks %in% c(comp$t1[i], comp$t2[i]) & panel==p)

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
      tmp$method = METHOD
      tmp$order = comp$order[i]
      out = rbind(out, tmp)
    }
  }
}

# save results
write.csv(out, file.path(STATS,'G00X_mag_postbaseline_v_postbaseline_fig3ABE_by_method.csv'), row.names = FALSE)

q(save='no')
