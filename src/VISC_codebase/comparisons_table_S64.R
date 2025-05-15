library(tidyverse)
library(gee)

# DIR = c('~/Repo/Schief856_G002/data/dropbox/Supp-Metrics/S46',
#         '~/Repo/Schief856_G002/data/dropbox/Supp-Metrics/S47')
# STATS = '~/Repo/Schief856_G002/data/dropbox/Stats.codebase'

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

DIR = c(file.path(SOURCEDIR, 'Sup-Metrics/S46'),
        file.path(SOURCEDIR, 'Sup-Metrics/S47'))
STATS = './Stats.codebase'

# pseudogroup names
pmap = c('eOD', 'eOD->eOD', 'eOD->core','eOD->eOD', 'eOD->core','eOD->eOD->core')


for( D in DIR ) {
  SUFF = ifelse(grepl('S46$', D), 'S46', 'S47')
  fn = list.files(D)
  print(fn)
  spr = c()
  for( n in fn ) {
    if (grepl("table-spr-", n)) next
    tmp = read.csv(file.path(D, n), check.names = FALSE)
    tmp$panel = substr(n, 4, 4)
    spr = rbind(spr, tmp)
  }


  spr = spr %>%
    select( panel, pubID, pseudogroup, is_cp, weeks, KD_fix) %>%
    mutate(pgrp = factor(pseudogroup, levels=c(1:6), labels=pmap),
           cp = factor(is_cp, levels=c('False', 'True'), labels=c('Random','Selected')))
  spr$y = log10(spr$KD_fix)

  # run analyses conditional on KD < 50uM
  spr = subset(spr, KD_fix < 5*(10^-5))


  comp = data.frame(pgrp1 = c(2,4,4,5,2,3,2,3),
                    pgrp2 = c(3,5,6,6,4,5,6,6),
                    wk1 = c(16, 24, 24, 24, 16, 16, 16, 16), # these have to match pseudogroup 1
                    wk2 = c(16, 24, 24, 24, 24, 24, 24, 24)) # these have to match pseudogroup 2

  stat = c()

  for(p in c('A', 'B', 'C', 'D') ) {
    for( i in 1:nrow(comp) ) {

      ss = subset(spr, pseudogroup %in% c(comp$pgrp1[i],comp$pgrp2[i]) &
                       panel==p )
      # sanity check for pseudogroup and week
      stopifnot(all((ss$pseudogroup == comp$pgrp1[i] & ss$weeks==comp$wk1[i]) |
                    (ss$pseudogroup == comp$pgrp2[i] & ss$weeks==comp$wk2[i]) ))

      if( comp$wk1[i] == comp$wk2[i] ) {
        ss = ss %>%
          mutate(pubID = factor(pubID),
                 pgrp = factor(pseudogroup,
                               levels=c(comp$pgrp1[i], comp$pgrp2[i]),
                               labels = pmap[c(comp$pgrp1[i], comp$pgrp2[i])]))%>%
          arrange(pubID)
      } else {
        ss = ss %>%
          mutate(pubID = factor(pubID),
                 pgrp = factor(weeks,
                               levels=c(comp$wk1[i], comp$wk2[i]),
                               labels = paste0('wk', c(comp$wk1[i], comp$wk2[i])))) %>%
          arrange(pubID)
      }

      fit = try(gee(y ~ pgrp, id=pubID, data=ss, family=gaussian, corstr = 'independence'))
      if( !inherits(fit, 'try-error') ) {
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

      tmp = data.frame( panel = p,
                        comp = ifelse(pmap[comp$pgrp1[i]]==pmap[comp$pgrp2[i]],
                                      pmap[comp$pgrp1[i]],
                                      sprintf('%s vs. %s', pmap[comp$pgrp1[i]], pmap[comp$pgrp2[i]])),
                        week = ifelse(comp$wk1[i]==comp$wk2[i], comp$wk1[i], sprintf('%d vs. %d', comp$wk1[i], comp$wk2[i])),
                        est = 10^est,
                        lci = 10^lci,
                        uci = 10^uci,
                        p.value = pval)

      stat = rbind(stat, tmp)
    }
  }

  stat = stat %>%
    group_by(panel) %>%
    mutate(fdr = p.adjust(p.value, 'fdr')) %>%
    ungroup()

  write.csv(stat, file.path(STATS, sprintf('SPR_supp_%s_pvalues.csv', SUFF)))
}

q(save='no')
