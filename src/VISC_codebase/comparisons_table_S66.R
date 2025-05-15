library(tidyverse)
library(gee)

# DIR = '~/Repo/Schief856_G002/data/dropbox/NeutData'
# STATS = '~/Repo/Schief856_G002/data/dropbox/Stats.codebase'

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

DIR = file.path(SOURCEDIR, 'NeutData')
# DIR = file.path(SOURCEDIR, 'Table-Metrics')
STATS = './Stats.codebase'

# FN="Table1.xlsx"
FN = 'G002_neuts_v3_two-mAbs-from-v2-deleted_06April2024--add non-nAb controls 09June2024--add gp3 mAbs 22Jan2025_CCfixedcolors add 16 more gp1+3 mAbs 05Feb2025.xlsx'

dat.main = readxl::read_excel(file.path(DIR, FN),
                         sheet = 'main', range = 'C2:S59' )
dat.supp = readxl::read_excel(file.path(DIR, FN),
                              sheet = 'supp', range = 'C3:S49' )

dat = rbind(dat.main, dat.supp) %>%
  fill(Group, Week, .direction = 'down')

# count data and binary outcomes greater than 1 'g1' and greater than or equal to 1 'ge1'
dat = dat %>%
  mutate(nneg = rowSums(across(5:17, ~ str_detect(.x, ">")), na.rm = TRUE),
         npos = length(5:17) - nneg,
         g1 = as.numeric(npos > 1),
         ge1 = as.numeric(npos > 0)) %>%
  arrange(Group, Week, Participant)

dat = dat %>%
  mutate(across(10:17, ~ ifelse(grepl("^>", .), 2*as.numeric(sub('>', '', .)), as.numeric(.)), .names = "converted_{.col}"),
         across(10:17, ~ ifelse(grepl("^>", .), as.numeric(sub('>', '', .)), as.numeric(.)), .names = "converted2_{.col}"))


dat$Group = factor(
  case_when(
    grepl('Group 1\\+3', dat$Group) ~ 'eOD->eOD',
    grepl('Group 2', dat$Group) ~ 'eOD->core',
    grepl('Group 3', dat$Group) ~ 'eOD->eOD->core'),
  levels = c('eOD->eOD', 'eOD->core', 'eOD->eOD->core'))


comp = data.frame(pgrp1 = c('eOD->core', 'eOD->core', 'eOD->core'),
                  pgrp2 = c('eOD->eOD', 'eOD->eOD->core', 'eOD->eOD->core'),
                  wk1 = c(16, 24, 16),
                  wk2 = c(16, 24, 24))

stat = c()

for(oc in c('g1', 'ge1', 'ic50', 'ic50_v2') ) {
  for( i in 1:nrow(comp) ) {

    if( oc %in% c('g1', 'ge1') ) {

      ss = subset(dat, (Group == comp$pgrp1[i] & Week==comp$wk1[i]) | (Group == comp$pgrp2[i] & Week==comp$wk2[i]))
      ss = ss %>%
        mutate(
          Participant = factor(Participant),
          pgrp = factor(Group, levels=c(comp$pgrp1[i], comp$pgrp2[i]))) %>%
        arrange(Participant)

      ss$y = ss[[oc]]

      fit = try(gee(y ~ pgrp, id=Participant, data=ss, family=binomial, corstr = 'independence'))
      if(!inherits(fit, 'try-error')) {
        est = summary(fit)$coefficients[2,1]
        se = summary(fit)$coefficients[2,"Robust S.E."]
        lci = est - qnorm(0.975)*se
        uci = est + qnorm(0.975)*se
        z = summary(fit)$coefficients[2,"Robust z"]
        pval = 2 * pnorm(-abs(z), lower.tail = TRUE)
      } else {
        est = lci = uci = pval = NA
      }
      rm(fit)

      est = exp(est)
      lci = exp(lci)
      uci = exp(uci)
    } else {

      prefix = ifelse(oc=='ic50', 'converted', 'converted2' )

      ss = subset(dat, (Group == comp$pgrp1[i] & Week==comp$wk1[i]) | (Group == comp$pgrp2[i] & Week==comp$wk2[i]))
      ss = ss %>%
        mutate(
          Participant = factor(Participant),
          pgrp = factor(Group, levels=c(comp$pgrp1[i], comp$pgrp2[i]))) %>%
          select(Group, pgrp, Week, Participant, Antibody, starts_with('converted')) %>%
          pivot_longer(cols=starts_with(prefix), values_to = 'y') %>%
          arrange(Participant, Antibody, name)

      ss$y = log10(ss$y)

      fit = gee(y ~ pgrp, id=Participant, data=ss, family=gaussian, corstr = 'independence')
      est = summary(fit)$coefficients[2,1]
      se = summary(fit)$coefficients[2,"Robust S.E."]
      lci = est - qnorm(0.975)*se
      uci = est + qnorm(0.975)*se
      z = summary(fit)$coefficients[2,"Robust z"]
      pval = 2 * pnorm(-abs(z), lower.tail = TRUE)
      rm(fit)

      est = 10^est
      lci = 10^lci
      uci = 10^uci
    }

    tmp = data.frame( outcome = oc,
                      comp = sprintf('%s vs. %s', comp$pgrp1[i], comp$pgrp2[i]),
                      week = ifelse(comp$wk1[i] == comp$wk2[i], comp$wk1[i], sprintf('%d vs. %d', comp$wk1[i], comp$wk2[i])),
                      est = est,
                      lci = lci,
                      uci = uci,
                      p.value = pval)

    stat = rbind(stat, tmp)
  }
}


stat = stat %>%
  group_by(outcome) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'Neut_pvalues.csv'))

q(save='no')
