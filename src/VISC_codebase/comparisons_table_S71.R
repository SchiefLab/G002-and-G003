library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure 4 comparisons
DIR = file.path(SOURCEDIR, 'Main-Metrics/fig4')
STATS = './Stats.codebase'

fn = list.files(DIR)
fn = fn[grepl('^fig',fn)]

# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')

dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n), check.names = FALSE)
  stopifnot(ncol(df)==5)
  stopifnot(all(names(df)[1:4] == c('pubID', 'trial', 'pseudogroup', 'weeks')))

  df$panel = substr(n, 4, 4)
  df$y_label = names(df)[5] # 5th column is the outcome variable
  df$y = df[,5]
  df[,5] = NULL
  df = df[,c('pubID','trial','pseudogroup','weeks','panel','y_label','y')]

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
stopifnot(all(dat$trial=='G002'))

# Account for 3 NAs in panel C
# NAs in panel C apply to all VRC01-panels (C, D, E, G, and H)
# remaining NAs, panel D only, converted to zero
ss = subset(dat, panel=='C' & is.na(y), select=c('pubID', 'weeks'))
ss$set.na = 1
tmp = dat
tmp$order = 1:nrow(tmp)
ss = merge(tmp, ss, all.x=TRUE)
ss = ss[order(ss$order),]
set.na = which(ss$set.na==1 & ss$panel %in% c('C', 'D', 'E', 'G', 'H'))
# stopifnot(length(set.na)==15)
set.zero = setdiff(which(is.na(ss$y) & ss$panel=='D'), set.na)

# all set.na values should already be NA except panel E where the response should be 0
stopifnot(all(is.na(dat$y)[set.na] | (dat$panel[set.na]=='E' & dat$y[set.na] == 0)))
dat$y[set.na] = NA

# all set.zero values are expected to be NA
stopifnot(all(is.na(dat$y[set.zero])))
dat$y[set.zero] = 0

# data should have three each in panels C, D, G, and H and
# have three unique samples (i.e. pubID and weeks)
# stopifnot(all(table(dat$panel, is.na(dat$y))[c('C','D','E','G','H'),2] == 3))
stopifnot(nrow(unique(subset(dat, is.na(y), select=c('pubID','weeks'))))==3)


# Comparison table
compt = data.frame( pgrp1=c(1,3),
                    pgrp2=c(4,7),
                    wk1=c(8,16),
                    wk2=c(16,24))

# tests are paired will only use ppts common across pseudogroup
compt$paired = TRUE

# create wide data set to compare fold change of group 2 at week 16 versus week 8
# to fold change for group 3 at week 16 versus 24
dat2 = dat %>%
  filter(pseudogroup %in% c(compt$pgrp1[1], compt$pgrp2[1])) %>%
  pivot_wider(id_cols = c(pubID, panel), names_from = weeks, names_prefix = 'y', values_from = y) %>%
  mutate( y = y16/y8, grpname='eOD->core') %>%
  select(-y16, -y8) %>%
  filter(!is.na(y))

dat3 = dat %>%
  filter(pseudogroup %in% c(compt$pgrp1[2], compt$pgrp2[2])) %>%
  pivot_wider(id_cols = c(pubID, panel), names_from = weeks, names_prefix = 'y', values_from = y) %>%
  mutate( y = y24/y16, grpname='eOD->eOD->core') %>%
  select(-y24, -y16) %>%
  filter(!is.na(y))



stat = c()

for( p in c('A', 'B', 'D') ) {
  for( i in 1:nrow(compt) ) {
    sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))

    tmp = pairwise_test_cont(
      x = sss$y, group = sss$weeks, id = sss$pubID,
      method = 'wilcox', paired = TRUE,
      alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
      verbose = FALSE)
    tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
    names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

    tmp$panel = p
    tmp$comp = sprintf('%s at wk%d vs. %s at wk%d', pmap[compt$pgrp1[i]], compt$wk1[i], pmap[compt$pgrp2[i]], compt$wk2[i])
    stat = rbind(stat, tmp)
  }

  sss = rbind(dat2, dat3) %>% filter(panel==p)

  tmp = pairwise_test_cont(
    x = sss$y, group = sss$grpname,
    method = 'wilcox', paired = FALSE,
    alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
    verbose = FALSE)
  tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
  names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

  tmp$panel = p
  tmp$comp = sprintf('%s (fold-change)', tmp$comp)
  stat = rbind(stat, tmp)
}

# add FDR adjustment
stat = stat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'fig4masking_pvals.csv'), row.names = FALSE)

q(save='no')
