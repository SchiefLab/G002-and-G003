library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure 5 comparisons by method
# for panels A, B, E, and F
DIR = file.path(SOURCEDIR, 'Main-Metrics/fig5')
STATS = './Stats.codebase'

stat_comb = c()

for( SUFF in c('') ) {
  METHOD = ifelse(SUFF=='', 'nearest', sub('^-', '', SUFF))
  D = paste0(DIR, SUFF)

  fn = list.files(D)

  # File names
  fn = fn[grepl('^figA_|figB_|figE_|figF_', fn)]

  # pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
  # the first two are week 8; sencond two are week 16; and last three are week 24
  pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')


  dat = lapply(1:4, function(i) {
    n = fn[i]
    df = read.csv(file.path(D, n), check.names = FALSE)
    if( i %in% 1:2 ) { # panels A and B
      stopifnot(all(names(df)[1:4] == c('pubID', 'trial', 'pseudogroup', 'weeks')))

      df$panel = substr(n, 4, 4)
      df$y_label = names(df)[5] # 5th column is the outcome variable
      df$y = df[,5]
      df = df[,c('pubID','pseudogroup','weeks','panel','y_label','y')]
    } else if( i==3 ) {
      stopifnot(ncol(df)==6)
      stopifnot(all(names(df)[1:6] == c('pubID', 'psname', 'weeks', 'pseudogroup','cottrell_focused_v_common_score', 'residues')))
      stopifnot(all(df$cottrell_focused_v_common_score==df$residues))

      df$panel = substr(n, 4, 4)
      df$y_label = names(df)[5] # 5th column is the outcome variable
      df$y = df[,5]
      df = df[,c('pubID','pseudogroup','weeks','panel','y_label','y')]

    } else if( i==4 ) {
      stopifnot(ncol(df)==6)
      stopifnot(all(names(df)[1:6] == c('pubID', 'psname', 'weeks', 'pseudogroup','num_hcdr2_mutations', 'residues')))
      stopifnot(all(df$num_hcdr2_mutations==df$residues))

      df$panel = substr(n, 4, 4)
      df$y_label = names(df)[5] # 5th column is the outcome variable
      df$y = df[,5]
      df = df[,c('pubID','pseudogroup','weeks','panel','y_label','y')]
    }

    return(df)
  })
  dat = do.call(rbind, dat)
  dat$grpname = factor(dat$pseudogroup,
                       levels=1:7,
                       labels=pmap)

  # Only Weeks 16 and 24 are comparable.
  dat = dat %>%
    filter(weeks %in% c(16, 24))

  stopifnot(!any(is.na(dat$pseudogroup)))
  stopifnot(all(dat$pseudogroup %in% 3:7)) # because we are only looking at weeks 16 and 24

  # check for NAs
  table(dat$panel, is.na(dat$y))

  # Comparison table
  compt = data.frame( pgrp1=c(3,5,5,6),
                      pgrp2=c(4,6,7,7),
                      weeks = c(16, 24, 24, 24),
                      compn=1:4)

  stat = c()

  for( p in unique(dat$panel) ) {
    for( i in 1:nrow(compt) ) {
      sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))
      stopifnot(all(sss$weeks==compt$weeks[i]))

      tmp = pairwise_test_cont(
        x = sss$y, group = sss$grpname,
        method = 'wilcox', paired = FALSE,
        alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
        verbose = FALSE)
      tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
      names(tmp) = c('comp', 'SampleSizes', 'Median_Min_Max', 'p.value')

      tmp$panel = p
      tmp$weeks = compt$weeks[i]
      stat = rbind(stat, tmp)
    }
  }

  stat$comp = factor(stat$comp,
                     levels=c("eOD->eOD vs. eOD->core", "eOD->eOD vs. eOD->eOD->core", "eOD->core vs. eOD->eOD->core"))
  stat = stat %>% arrange(panel, weeks, comp)


  # additional inter time point analysis
  # Comparison table
  compt = data.frame( pgrp1=c(3, 4, 3, 4),
                      pgrp2=c(5, 6, 7, 7),
                      wk1 = c(16, 16, 16, 16),
                      wk2 = c(24, 24, 24, 24),
                      paired = c(TRUE, TRUE, TRUE, FALSE))

  stat_wk = c()

  for( p in unique(dat$panel) ) {
    for( i in 1:nrow(compt) ) {
      ss1 = subset(dat, panel==p & pseudogroup == compt$pgrp1[i] & weeks == compt$wk1[i])
      ss1$grpname = sprintf('%s at wk%d', pmap[compt$pgrp1[i]], compt$wk1[i])
      ss2 = subset(dat, panel==p & pseudogroup == compt$pgrp2[i] & weeks == compt$wk2[i])
      ss2$grpname = sprintf('%s at wk%d', pmap[compt$pgrp2[i]], compt$wk2[i])
      sss = rbind(ss1, ss2)

      if( compt$paired[i] ) {
        tmp = pairwise_test_cont(
          x = sss$y, group = sss$grpname, id = sss$pubID,
          method = 'wilcox', paired = TRUE,
          alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
          verbose = FALSE)
      } else {
        tmp = pairwise_test_cont(
          x = sss$y, group = sss$grpname,
          method = 'wilcox', paired = FALSE,
          alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
          verbose = FALSE)
      }
      tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
      names(tmp) = c('comp', 'SampleSizes', 'Median_Min_Max', 'p.value')

      tmp$panel = p
      tmp$weeks = compt$weeks[i]
      stat_wk = rbind(stat_wk, tmp)
    }
  }

  stat_wk$comp = factor(stat_wk$comp,
                     levels=c("eOD->eOD at wk16 vs. eOD->eOD at wk24",
                              "eOD->core at wk16 vs. eOD->core at wk24",
                              "eOD->eOD at wk16 vs. eOD->eOD->core at wk24",
                              "eOD->core at wk16 vs. eOD->eOD->core at wk24"),
                     labels = c("eOD->eOD",
                                "eOD->core",
                                "eOD->eOD vs. eOD->eOD->core",
                                "eOD->core vs. eOD->eOD->core"))
  stat_wk$weeks = '16 vs. 24'

  # add FDR adjustment
  stat = rbind(stat, stat_wk) %>%
    arrange(panel) %>%
    group_by(panel) %>%
    mutate(fdr = p.adjust(p.value, method='fdr')) %>%
    ungroup()
  stat$method = METHOD

  stat_comb = rbind(stat_comb, stat)
}


write.csv(stat_comb, file.path(STATS, 'fig5_by_method_pvals.csv'), row.names = FALSE)

q(save='no')
