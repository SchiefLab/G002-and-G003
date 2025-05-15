library(tidyverse)
library(VISCfunctions)

print('Comparisons table for figure S35')

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure S35 comparisons
# DIR = file.path(SOURCEDIR, 'Supp-Metrics/S29') 29 was the old name
DIR = file.path(SOURCEDIR, 'Sup-Metrics/S35')
if (!dir.exists(DIR)) {
  stop("Directory does not exist: ", DIR)
}

STATS = './Stats.codebase'
if (!dir.exists(STATS)) {
  dir.create(STATS, recursive = TRUE, showWarnings = FALSE)
  message("Created directory: ", STATS)
}

fn = list.files(DIR)
fn = fn[grepl('boost_clonality', fn)]

# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')

dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n), check.names = FALSE)
  df$psname = NULL # not needed
  stopifnot(ncol(df)==4)
  stopifnot(all(names(df)[1:3] %in% c('pubID', 'pseudogroup', 'weeks')))

  oc = names(df)[4]
  df$y_label = oc
  df$y = df[[oc]]
  df[[oc]] = NULL

  df$panel = substr(n, 4, 4)

  return(df)
})
dat = do.call(rbind, dat)

stopifnot(all(dat$pseudogroup %in% 1:7))
dat$grpname = factor(dat$pseudogroup,
                     levels=1:7,
                     labels=pmap)

# Comparison table
compt = data.frame( pgrp1=c(1,3,5,5,6),
                    pgrp2=c(2,4,6,7,7),
                    weeks = c(8, 16, 24, 24, 24),
                    compn=1:5)

stat = c()

for( p in unique(dat$panel) ) {
  for( i in 1:nrow(compt) ) {
    sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))

    tmp = pairwise_test_cont(
      x = sss$y, group = sss$grpname,
      method = 'wilcox', paired = FALSE,
      alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
      verbose = FALSE)
    tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
    names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

    tmp$panel = p
    tmp$weeks = compt$weeks[i]
    stat = rbind(stat, tmp)
  }
}

stat$comp = factor(stat$comp,
                   levels=c("eOD vs. core", "eOD->eOD vs. eOD->core", "eOD->eOD vs. eOD->eOD->core", "eOD->core vs. eOD->eOD->core"))
stat = stat %>% arrange(panel, weeks, comp)


# add FDR adjustment
stat = stat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'figS29_pvals.csv'), row.names = FALSE)

q(save='no')
