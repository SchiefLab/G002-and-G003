args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

library(tidyverse)
library(gee)

# comparisons for eOD at wk8 verus eOD->eOD at wk16
# Fig. 6A comparisons:
# random mAbs: eOD at week 8 vs eOD-eOD at week 16
# selected mAbs: eOD at week 8 vs eOD-eOD at week 16

# DIR = file.path(SOURCEDIR, 'SPR-Metrics')
DIR = file.path(SOURCEDIR, 'Main-Metrics/fig6')
STATS = 'Stats.codebase'

# pseudogroup names
pmap = c('eOD', 'eOD->eOD', 'eOD->core','eOD->eOD', 'eOD->core','eOD->eOD->core')

# SPR for core-g28v2 and core-N276
spr.g28 = read.csv(file.path(DIR, 'figA-spr-g28v2.csv'))
spr.g28$panel = 'A'
if ("ptid" %in% names(spr.g28)) {
  spr.g28 <- spr.g28 %>% rename(pubID = ptid)
}
spr.g28 = spr.g28 %>% dplyr::select(panel, pubID, pseudogroup, is_cp, weeks, KD_fix)

spr.n276 = read.csv(file.path(DIR, 'figB-spr-N276.csv'))
spr.n276$panel = 'B'
if ("ptid" %in% names(spr.n276)) {
  spr.n276 <- spr.n276 %>% rename(pubID = ptid)
}
spr.n276 = spr.n276 %>% dplyr::select(panel, pubID, pseudogroup, is_cp, weeks, KD_fix)

spr = rbind(spr.g28, spr.n276) %>%
  mutate(pgrp = factor(pseudogroup, levels=c(1:6), labels=pmap),
         cp = factor(is_cp, levels=c('False', 'True'), labels=c('Random','Selected')))
spr$y = log10(spr$KD_fix)

# run analyses conditional on KD < 10^-4
spr = subset(spr, KD_fix < 10^-4)


comp = data.frame(pgrp1 = c(1),
                  pgrp2 = c(2),
                  wk1 = c(8), # these have to match pseudogroup 1
                  wk2 = c(16)) # these have to match pseudogroup 2

stat = c()

for(p in c('A') ) {
  for( set in levels(spr$cp) ) {
    for( i in 1:nrow(comp) ) {
      ss = subset(spr, pseudogroup %in% c(comp$pgrp1[i],comp$pgrp2[i]) &
                       cp==set &
                       panel==p )
      # sanity check for pseudogroup and week
      stopifnot(all((ss$pseudogroup == comp$pgrp1[i] & ss$weeks==comp$wk1[i]) |
                    (ss$pseudogroup == comp$pgrp2[i] & ss$weeks==comp$wk2[i]) ))

      ss = ss %>%
        mutate(pubID = factor(pubID),
               pgrp = factor(weeks,
                             levels=c(comp$wk1[i], comp$wk2[i]),
                             labels = paste0('wk', c(comp$wk1[i], comp$wk2[i])))) %>%
        arrange(pubID)

      fit = gee(y ~ pgrp, id=pubID, data=ss, family=gaussian, corstr = 'independence')
      est = summary(fit)$coefficients[2,1]
      se = summary(fit)$coefficients[2,"Robust S.E."]
      lci = est - qnorm(0.975)*se
      uci = est + qnorm(0.975)*se
      z = summary(fit)$coefficients[2,"Robust z"]
      pval = 2 * pnorm(-abs(z), lower.tail = TRUE)


      tmp = data.frame( panel = p,
                        comp = ifelse(pmap[comp$pgrp1[i]]==pmap[comp$pgrp2[i]],
                                      pmap[comp$pgrp1[i]],
                                      sprintf('%s vs. %s', pmap[comp$pgrp1[i]], pmap[comp$pgrp2[i]])),
                        week = ifelse(comp$wk1[i]==comp$wk2[i], comp$wk1[i], sprintf('%d vs. %d', comp$wk1[i], comp$wk2[i])),
                        set = set,
                        est = 10^est,
                        lci = 10^lci,
                        uci = 10^uci,
                        p.value = pval)

      stat = rbind(stat, tmp)
    }
  }
}

stat = stat %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>%
  ungroup()

write.csv(stat, file.path(STATS, 'SPR_eOD_pvalues.csv'))

q(save='no')
