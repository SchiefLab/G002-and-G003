args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

library(tidyverse)
library(VISCfunctions)

# Figure 2F vs. 4G (different probesets)
DIR2 = file.path(SOURCEDIR, 'Main-Metrics/fig2')
DIR4 = file.path(SOURCEDIR, 'Main-Metrics/fig4')
STATS = 'Stats.codebase'

fn2 = c(
  "figF_percent_vrc01-class_amongcd4bs-specificigg$^{+}$_memory_sequences.csv",
  "figG_percent_vrc01-class_amonggt8-specificigg$^{+}$_memory_sequences.csv"
)
fn4 = c(
  "figG_percent VRC01-class among CD4bs-specific IgG$^{+}$ memory sequences.csv",
  "figH_percent VRC01-class among core-specific IgG$^{+}$ memory sequences.csv"
)


# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')

# Figure 2 data; restricting to panel F, trial G002, and weeks 16 and 24
# pseudogroup is 3 at week 16 and 5 at week 24
dat2f = read.csv(file.path(DIR2, fn2[1]), check.names = FALSE)
dat2f$panel = '2F'
dat2f$y = dat2f$`Percent of IGHG sequences that are VRC01-class`
dat2f$`Percent of IGHG sequences that are VRC01-class` = NULL
dat2g = read.csv(file.path(DIR2, fn2[2]), check.names = FALSE)
dat2g$panel = '2G'
dat2g$y = dat2g$`Percent VRC01-class among eOD-specific IgG+ memory BCR sequences`
dat2g$`Percent VRC01-class among eOD-specific IgG+ memory BCR sequences` = NULL
dat2 = rbind(dat2f, dat2g)
dat2 = subset(dat2, trial=='G002' & weeks %in% c(16, 24))
dat2$pseudogroup = case_when(
  dat2$weeks==16 ~ 3,
  dat2$weeks==24 ~ 5
)



# Figure 4 data
dat4g = read.csv(file.path(DIR4, fn4[1]), check.names = FALSE)
dat4g$panel = '4G'
dat4g$y = dat4g$`Percent of IGHG sequences that are VRC01-class`
dat4g$`Percent of IGHG sequences that are VRC01-class` = NULL
dat4h = read.csv(file.path(DIR4, fn4[2]), check.names = FALSE)
dat4h$panel = '4H'
dat4h$y = dat4h$`Percent VRC01-class among Core-specific IgG+ memory BCR sequences`
dat4h$`Percent VRC01-class among Core-specific IgG+ memory BCR sequences` = NULL
dat4 = rbind(dat4g, dat4h)

# Restrict to week 8, 16 and 24 data
dat4 = dat4 %>%
  filter( weeks %in% c(8, 16, 24) )

# various checks
stopifnot(!any(is.na(dat4$pseudogroup)))
stopifnot(all(dat4$pseudogroup %in% 1:7))
stopifnot(all(dat4$weeks %in% c(8, 16, 24)))
stopifnot(all(dat4$trial=='G002'))

# Figure 4, account for 3 NAs in panel C
# NAs in panel C apply to all VRC01-panels (C, D, E, G, and H)
# remaining NAs, panel D only, converted to zero
# we only have panel G and H for this analysis so there should be 3 NAs
set.na = which(is.na(dat4$y))
# stopifnot(length(set.na)==6)  # no longer true for figure 4


# combine data and set pgroup name
dat2$probe = 'eOD'
dat4$probe = 'core'
dat = rbind(dat2, dat4)
dat$grpname = factor(dat$pseudogroup,
                      levels=1:7,
                      labels=pmap)

# set outcome name
dat$outcome_name = factor(case_when(
  dat$panel %in% c('4G','2F') ~ 'VRC01-class among CD4bs-specific IgG B cells',
  dat$panel %in% c('4H','2G') ~ 'VRC01-class among antigen-specific IgG B cells'
),
  levels = c('VRC01-class among CD4bs-specific IgG B cells',
             'VRC01-class among antigen-specific IgG B cells')
)



# Comparison table
compt = data.frame( pgrp1=c(3, 4, 5, 6, 7), # core sort data
                    pgrp2=c(3, 3, 5, 5, 5), # eOD sort data
                    weeks = c(16, 16, 24, 24, 24) )
compt$is.paired = compt$pgrp1 == compt$pgrp2
compt$grp1 = sprintf('%s at wk%d', pmap[compt$pgrp1], compt$weeks)
compt$grp2 = sprintf('%s at wk%d', pmap[compt$pgrp2], compt$weeks)


stat = c()

for( oc in levels(dat$outcome_name) ) {
  for( i in 1:nrow(compt) ) {
    ss.core = subset(dat, pseudogroup == compt$pgrp1[i] & probe=='core' & outcome_name==oc)
    ss.eod  = subset(dat, pseudogroup == compt$pgrp2[i] & probe=='eOD'  & outcome_name==oc)
    sss = rbind(ss.core, ss.eod)

    # make sure we have the right order for groups
    sss$probe = factor(sss$probe, levels=c('core', 'eOD'))

    if( compt$is.paired[i] ) {
      tmp = pairwise_test_cont(
        x = sss$y, group = sss$probe, sss$pubID,
        method = 'wilcox', paired = TRUE,
        alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
        verbose = FALSE)
    } else {
      tmp = pairwise_test_cont(
        x = sss$y, group = sss$probe,
        method = 'wilcox', paired = FALSE,
        alternative = 'two.sided', num_needed_for_test = 3, digits = 4,
        verbose = FALSE)
    }

    tmp = tmp %>% dplyr::select(Comparison, SampleSizes, Median_Min_Max, MagnitudeTest)
    names(tmp) = c('comp', 'SampleSizes', 'summary', 'p.value')

    # use the compt grp names instead of probe names for 'comp'
    tmp$comp = sprintf('%s vs. %s', compt$grp1[i], compt$grp2[i])

    tmp$is.paired = compt$is.paired[i]
    tmp$weeks = compt$weeks[i]
    tmp$outcome_name = oc
    stat = rbind(stat, tmp)
  }
}

# add FDR adjustment
stat = stat %>%
  group_by(outcome_name) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  ungroup()


write.csv(stat, file.path(STATS, 'probe_comparison_pvals.csv'), row.names = FALSE)

q(save='no')
