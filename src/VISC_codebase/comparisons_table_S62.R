library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure S31 comparisons
# DIR = file.path(SOURCEDIR, 'Supp-Metrics/S31')
DIR = file.path(SOURCEDIR, 'Sup-Metrics/S39')
STATS = './Stats.codebase'

fn = list.files(DIR)
fn = fn[!grepl('figA', fn)]
fn = fn[!grepl('figC', fn)]
fn = fn[!grepl('figB_lambda', fn)]
fn = fn[!grepl('_seqs.csv', fn)]


# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')

dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n))
  oc = names(df)[5]
  df$oc = df[,5]
  df$oc_label = oc
  stopifnot(names(df)[6] %in% c('freq', 'frequency'))
  names(df)[6] = 'y'
  df[,5] = NULL

  df$panel = substr(n, 4, 4)
  if(grepl('kappa', n)) df$panel = paste0(df$panel, '_kappa')
  if(grepl('lambda', n)) df$panel = paste0(df$panel, '_lambda')

  return(df)
})
dat = do.call(rbind, dat)
dat$grpname = factor(dat$pseudogroup,
                     levels=1:7,
                     labels=pmap)


# Comparison table
compt = data.frame( pgrp1=c(3,5,5,6),
                    pgrp2=c(4,6,7,7),
                    weeks = c(16, 24, 24, 24) )


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


# add FDR adjustment
stat = stat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'figS31_pvals.csv'), row.names = FALSE)

q(save='no')
