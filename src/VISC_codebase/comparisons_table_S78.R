library(tidyverse)
library(gee)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure 8A comparisons
# DIR = file.path(SOURCEDIR, 'External')
DIR = file.path(SOURCEDIR, 'Main-Metrics/fig8')
STATS = './Stats.codebase'

# fn = "Fig8 data in Data S10 simplified and study week added 16Feb2025.xlsx"
fn = "fig8A.csv"

# pseudogroup map; pseudogroups are labelled 1 to 7 and define by the following 7 x-axis labels
# the first two are week 8; sencond two are week 16; and last three are week 24
pmap = c('eOD', 'core', 'eOD->eOD', 'eOD->core', 'eOD->eOD', 'eOD->core','eOD->eOD->core')
# dat = readxl::read_excel(file.path(DIR, fn), range = 'A4:R66')
dat = read_csv(file.path(DIR, fn))

map = c(
  "core-g28v2-N276",
  "001428-core",
  "BG505-core",
  "BG505-cd4bsHxB2",
  "BG505-cd4bsHxB2-M278",
  "BG505-T278M",
  "1HD2-N276Q",
  "001428-T278M",
  "V703-0537-T278M",
  "V703-0739-T278M",
  "CAP260-T278M",
  "CNE40-T278M",
  "235-T278M")


dat$pubID = substr(dat$cellid, 1, 8)
dat$week = dat$`study week`
table(dat$pubID, dat$`study week`)


stat = c()

for( oc in map ) {
  ss = dat
  ss = ss %>%
    mutate(
      pubID = factor(pubID),
      pgrp = factor(week, levels=c(16, 24))) %>%
    arrange(pubID)

  ss$y = ss[[oc]]
  ss = subset(ss, y < 5*10^-5) # analysis is conditional on KD less than 50 uM
  ss$y = log10(ss$y)

  fit = gee(y ~ pgrp, id=pubID, data=ss, family=gaussian, corstr = 'independence')
  est = summary(fit)$coefficients[2,1]
  wk16 = summary(fit)$coefficients[1,1]
  wk24 = sum(wk16, est)
  se = summary(fit)$coefficients[2,"Robust S.E."]
  lci = est - qnorm(0.975)*se
  uci = est + qnorm(0.975)*se
  z = summary(fit)$coefficients[2,"Robust z"]
  pval = 2 * pnorm(-abs(z), lower.tail = TRUE)
  rm(fit)

  wk16 = 10^wk16
  wk24 = 10^wk24
  est = 10^est
  lci = 10^lci
  uci = 10^uci

  tmp = data.frame( outcome = oc,
                    comp = '16 vs. 24',
                    wk16,
                    wk24,
                    est = est,
                    lci = lci,
                    uci = uci,
                    p.value = pval)

  stat = rbind(stat, tmp)
}

stat = stat %>%
  mutate(fdr = p.adjust(p.value, 'fdr'))

write.csv(stat, file.path(STATS, 'Fig8A_pvalues.csv'))

q(save='no')
