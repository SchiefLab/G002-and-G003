library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure S29 comparisons
# DIR = file.path(SOURCEDIR, 'Supp-Metrics/S29')
DIR = file.path(SOURCEDIR, 'Sup-Metrics/S35')
STATS = './Stats.codebase'

fn = list.files(DIR)
fn = fn[grepl('boost_clonality', fn)]

# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')

dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n), check.names = FALSE)
  df$psname = NULL
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

# pseudogroup corresponds to a grpname and week
label = dat %>%
  dplyr::select(pseudogroup, grpname, weeks) %>%
  unique() %>%
  mutate(grpwk = sprintf('%s wk%d', grpname, weeks)) %>%
  select(pseudogroup, grpwk) %>%
  arrange(pseudogroup)


# Comparison table,
# a.	eOD wk8 vs eOD->eOD wk16
# b.	eOD wk8 vs eOD->core wk16
# c.	eOD wk8 vs eOD->eOD wk24
# d.	eOD wk8 vs eOD->core wk24
# e.	eOD wk8 vs eOD->eOD->core wk24
# f.	eOD->eOD wk16 vs eOD->core wk16
# g.	eOD->eOD wk 24 vs eOD->core wk24
compt = data.frame( pgrp1 = c(1,1,1,1,1,3,5),
                    pgrp2 = c(3,4,5,6,7,4,6),
                    is.paired = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE) )
compt$lbl1 = label[compt$pgrp1, 'grpwk']
compt$lbl2 = label[compt$pgrp2, 'grpwk']
compt$comp = sprintf('%s vs. %s', compt$lbl1, compt$lbl2)

# check which comparisons are paired and which are unpaired
# tmp = c()
# for( p in unique(dat$panel) ) {
#   for( i in 1:nrow(compt) ) {
#     sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))
#     tmp = rbind(tmp, data.frame( panel=p, i=i, paired = as.numeric(max(as.numeric(names(table(table(sss$pubID))))) == 2)))
#   }
# }

stat = c()

for( p in unique(dat$panel) ) {
  for( i in 1:nrow(compt) ) {
    sss = subset(dat, panel==p & pseudogroup %in% c(compt$pgrp1[i], compt$pgrp2[i]))

    tmp = c()
    if( compt$is.paired[i] ) {
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
    names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

    tmp$is.paired = compt$is.paired[i]
    tmp$panel = p
    tmp$comp = compt$comp[i]
    stat = rbind(stat, tmp)
  }
}

stat$comp = factor(stat$comp, levels=compt$comp)
stat = stat %>% arrange(panel, comp)


# add FDR adjustment
stat = stat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'figS29_v2_pvals.csv'), row.names = FALSE)

q(save='no')
