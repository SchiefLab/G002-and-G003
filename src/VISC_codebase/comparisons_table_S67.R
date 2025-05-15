library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure 7 comparisons
DIR = file.path(SOURCEDIR, 'Main-Metrics/fig7')
STATS = './Stats.codebase'

fn = list.files(DIR)
fn = fn[!grepl('^figX',fn)]

# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')

dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n), check.names = FALSE)
  df$psname = NULL # not needied
  stopifnot(ncol(df)==4)
  stopifnot( names(df)[1] == c('pubID'))
  stopifnot(all(names(df)[3:4] == c('pseudogroup', 'weeks')))

  df$panel = substr(n, 4, 4)
  df$y_label = names(df)[2] # 2th column is the outcome variable
  df$y = df[,2]
  df[,2] = NULL
  df = df[,c('pubID','pseudogroup','weeks','panel','y_label','y')]

  return(df)
})
dat = do.call(rbind, dat)
dat$grpname = factor(dat$pseudogroup,
                     levels=1:7,
                     labels=pmap)

# remove week -5 and 4 data (not in figure)
# Week 8, 16, and 24 should remain.
dat = dat %>%
  filter( !(weeks %in% c(-5, 4)) ) # not used

# various checks
stopifnot(!any(is.na(dat$pseudogroup)))
stopifnot(all(dat$pseudogroup %in% 1:7))
stopifnot(all(dat$weeks %in% c(8, 16, 24)))

# Comparison table
compt = data.frame( pgrp1=c(1,3,5,5,6),
                    pgrp2=c(2,4,6,7,7),
                    weeks = c(8, 16, 24, 24, 24),
                    compn=1:5)

stat = c()

for( p in unique(dat$panel) ) {
  for( i in 1:nrow(compt) ) {
    sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))

    if( p %in% c('B', 'F') ) {
      samplesize = sss %>% group_by(grpname) %>% summarize(n=sum(!is.na(y)), pos=sum(y, na.rm=TRUE), rate=sprintf("%d/%d", pos, n))
      tmp = pairwise_test_bin(
        x = sss$y, group = sss$grpname,
        method = 'barnard', paired = FALSE,
        alternative = 'two.sided', num_needed_for_test = 3, digits = 1,
        latex_output = TRUE, verbose = FALSE)

      # check sample size calc against ResponseStats
      stopifnot(grepl(sprintf("^%s", samplesize$rate[1]), tmp$ResponseStats))
      stopifnot(grepl(sprintf("vs. %s", samplesize$rate[2]), tmp$ResponseStats))
      tmp$SampleSizes = sprintf('%d vs. %d', samplesize$n[1], samplesize$n[2])

      tmp = tmp %>% dplyr::select(Comparison, SampleSizes, ResponseStats, ResponseTest)

      names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')
    } else {
      tmp = pairwise_test_cont(
        x = sss$y, group = sss$grpname,
        method = 'wilcox', paired = FALSE,
        alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
        verbose = FALSE)
      if( is.null(tmp) ) next # note enough observations to run test
      tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
      names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')
    }

    tmp$panel = p
    tmp$weeks = compt$weeks[i]
    stat = rbind(stat, tmp)
  }
}

# Comparison table for across time comparisons
compt = data.frame( pgrp1=c(3,4,3,4),
                    pgrp2=c(5,6,7,7),
                    wk1 = c(16, 16, 16, 16),
                    wk2 = c(24, 24, 24, 24))

# tests are paired if the pseudogroup name is the same across group 1 and 2 and
# week 1 is 16 and week 2 is 24 or eOD->eOD at wk16 vs. eOD->eOD->core at wk24
compt$paired = (pmap[compt$pgrp1] == pmap[compt$pgrp2] & compt$wk1==16 & compt$wk2==24) |
  (pmap[compt$pgrp1] == 'eOD->eOD' & pmap[compt$pgrp2] == 'eOD->eOD->core' & compt$wk1==16 & compt$wk2==24)

stat_wk = c()

for( p in unique(dat$panel) ) {
  for( i in 1:nrow(compt) ) {
    sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))

    if( compt$paired[i] ) {
      if( p %in% c('B', 'F') ) {
        sss.paired = sss %>%
          pivot_wider(id_cols = pubID, values_from = y, names_from = weeks, names_prefix = 'wk') %>%
          filter( !is.na(wk16) & !is.na(wk24))

        tmp = pairwise_test_bin(
          x = sss$y, group = sss$weeks, id = sss$pubID,
          method = 'mcnemar', paired = TRUE,
          alternative = 'two.sided', num_needed_for_test = 3, digits = 1,
          latex_output = TRUE, verbose = FALSE)

        # check sample size calc against ResponseStats
        stopifnot(grepl(sprintf("^%d/%d", sum(sss.paired$wk16), nrow(sss.paired)), tmp$ResponseStats))
        stopifnot(grepl(sprintf("%d/%d", sum(sss.paired$wk24), nrow(sss.paired)), tmp$ResponseStats))
        tmp$SampleSizes = sprintf('%d', nrow(sss.paired))

        tmp = tmp %>% dplyr::select(Comparison, SampleSizes, ResponseStats, ResponseTest)

        names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')
      } else {
        tmp = pairwise_test_cont(
          x = sss$y, group = sss$weeks, id = sss$pubID,
          method = 'wilcox', paired = TRUE,
          alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
          verbose = FALSE)
        if( is.null(tmp) ) next # note enough observations to run test
        tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
        names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')
      }
    } else {
      if( p %in% c('B', 'F') ) {
        samplesize = sss %>% group_by(grpname) %>% summarize(n=sum(!is.na(y)), pos=sum(y, na.rm=TRUE), rate=sprintf("%d/%d", pos, n))
        tmp = pairwise_test_bin(
          x = sss$y, group = sss$grpname,
          method = 'barnard', paired = FALSE,
          alternative = 'two.sided', num_needed_for_test = 3, digits = 1,
          latex_output = TRUE, verbose = FALSE)

        # check sample size calc against ResponseStats
        stopifnot(grepl(sprintf("^%s", samplesize$rate[1]), tmp$ResponseStats))
        stopifnot(grepl(sprintf("vs. %s", samplesize$rate[2]), tmp$ResponseStats))
        tmp$SampleSizes = sprintf('%d vs. %d', samplesize$n[1], samplesize$n[2])

        tmp = tmp %>% dplyr::select(Comparison, SampleSizes, ResponseStats, ResponseTest)

        names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')
      } else {
        tmp = pairwise_test_cont(
          x = sss$y, group = sss$grpname,
          method = 'wilcox', paired = FALSE,
          alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
          verbose = FALSE)
        if( is.null(tmp) ) next # note enough observations to run test
        tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
        names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')
      }
    }

    stopifnot( (tmp$comp == sprintf('%s vs. %s', pmap[compt$pgrp1[i]], pmap[compt$pgrp2[i]])) |
                 (tmp$comp == sprintf('%d vs. %d', compt$wk1[i], compt$wk2[i])) )

    tmp$panel = p
    tmp$comp = ifelse(pmap[compt$pgrp1[i]] == pmap[compt$pgrp2[i]],
                      pmap[compt$pgrp1[i]],
                      sprintf('%s vs. %s', pmap[compt$pgrp1[i]], pmap[compt$pgrp2[i]]))
    tmp$weeks = ifelse(compt$wk1[i] == compt$wk2[i],
                       compt$wk1[i],
                       sprintf('%d vs. %d', compt$wk1[i], compt$wk2[i]))

    stat_wk = rbind(stat_wk, tmp)
  }
}


# add FDR adjustment
stat = rbind(stat, stat_wk) %>%
  arrange(panel) %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'fig7_pvals.csv'), row.names = FALSE)

q(save='no')
