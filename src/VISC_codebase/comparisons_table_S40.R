library(tidyverse)
library(VISCfunctions)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# Figure S19 comparisons
# DIR = '~/Repo/Schief856_G002/data/dropbox/Supp-Metrics/S19-alt'
# DIR = '/Users/tmsincomb/Dropbox/repos/G002-and-G003/G00X-plots/Sup-Metrics/S19'
# TODO: update to ./stats
# STATS = '~/Repo/Schief856_G002/data/dropbox/Stats.codebase'

DIR = file.path(SOURCEDIR, 'Sup-Metrics/S19')
STATS = './Stats.codebase'

fn = list.files(DIR)


dat = lapply(fn, function(n) {
  df = read.csv(file.path(DIR, n), check.names = FALSE)
  # If PTID column exists, rename it to pubID
  if ("PTID" %in% names(df)) {
    names(df)[names(df) == "PTID"] = "pubID"
  }
  stopifnot(ncol(df)==4)
  stopifnot(all(sort(c('trial', 'method', 'pubID')) %in% sort(names(df)[1:3])))
  stopifnot(all(df$trial == df$method))

  oc = names(df)[4]
  df$y_label = oc
  df$y = df[[oc]]
  df[[oc]] = NULL
  df$panel = substr(n, 4, 4)

  return(df)
})
dat = do.call(rbind, dat)
stopifnot(length(unique(dat$pubID))==10)
stopifnot(all(dat$trial == dat$method))
dat$trial = NULL

# check for NAs
stopifnot(!any(is.na(dat$y)))

# check for zeros in data that will be log transformed panels A, B, D, G, H, and J
# w.log = which(dat$panel %in% c('A','B','D','G','H','J'))
w.log = which(dat$panel %in% c('A','B','C','F','G','H'))
stopifnot(all(dat$y[w.log] > 0))


### run the comparisons ties create warnings which are suppressed
# G001 vs G002 in A, B, C, D, E, F
stat_col1 = c()
for( p in c('A','B','C','D','E','F') ) {
# for( p in c('A','B','C','D','E') ) {
  sss = subset(dat, panel==p)
  # if( p %in% c('A', 'B', 'D') ) sss$y = log10(sss$y)
  if( p %in% c('A', 'B', 'C') ) sss$y = log10(sss$y)
  sss = sss %>% pivot_wider(id_cols = pubID, names_from = method, values_from = y)

  t = t.test(sss$G002, sss$G001, paired = TRUE)
  tmp = data.frame(n = t$parameter + 1,
                   estimate=t$estimate,
                   lci = t$conf.int[1],
                   uci = t$conf.int[2],
                   p.value = t$p.value,
                   comp = 'G001 vs. G002',
                   panel = p)

  # if( p %in% c('A', 'B', 'D') ) {
  if( p %in% c('A', 'B', 'C') ) {
    tmp$estimate = 10^tmp$estimate
    tmp$lci = 10^tmp$lci
    tmp$uci = 10^tmp$uci
  }

  stat_col1 = rbind(stat_col1, tmp)
}

# G001 vs G003 in G, H, I, J, K, L
stat_col2_g001 = c()
# for( p in c('G', 'H', 'I', 'J', 'K', 'L') ) {
for( p in c('G', 'H', 'I', 'J') ) {
  sss = subset(dat, panel==p & method %in% c('G001', 'G003'))
  # if( p %in% c('G', 'H', 'J') ) sss$y = log10(sss$y)
  if( p %in% c('F', 'G', 'H') ) sss$y = log10(sss$y)
  sss = sss %>% pivot_wider(id_cols = pubID, names_from = method, values_from = y)

  t = t.test(sss$G003, sss$G001, paired = TRUE)
  tmp = data.frame(n = t$parameter + 1,
                   estimate=t$estimate,
                   lci = t$conf.int[1],
                   uci = t$conf.int[2],
                   p.value = t$p.value,
                   comp = 'G001 vs. G003',
                   panel = p)

  # if( p %in% c('G', 'H', 'J') ) {
  if( p %in% c('F', 'G', 'H') ) {
    tmp$estimate = 10^tmp$estimate
    tmp$lci = 10^tmp$lci
    tmp$uci = 10^tmp$uci
  }

  stat_col2_g001 = rbind(stat_col2_g001, tmp)
}

# G002 vs G003 in G, H, I, J, K, L
stat_col2_g002 = c()
# for( p in c('G', 'H', 'I', 'J', 'K', 'L') ) {
for( p in c('F', 'G', 'H', 'I', 'J') ) {
  sss = subset(dat, panel==p & method %in% c('G002', 'G003'))
  # if( p %in% c('G', 'H', 'J') ) sss$y = log10(sss$y)
  if( p %in% c('F', 'G', 'H') ) sss$y = log10(sss$y)
  sss = sss %>% pivot_wider(id_cols = pubID, names_from = method, values_from = y)

  t = t.test(sss$G003, sss$G002, paired = TRUE)
  tmp = data.frame(n = t$parameter + 1,
                   estimate=t$estimate,
                   lci = t$conf.int[1],
                   uci = t$conf.int[2],
                   p.value = t$p.value,
                   comp = 'G002 vs. G003',
                   panel = p)

  # if( p %in% c('G', 'H', 'J') ) {
  if( p %in% c('F', 'G', 'H') ) {
    tmp$estimate = 10^tmp$estimate
    tmp$lci = 10^tmp$lci
    tmp$uci = 10^tmp$uci
  }

  stat_col2_g001 = rbind(stat_col2_g001, tmp)
}



# add FDR adjustment
stat = rbind(stat_col1, stat_col2_g001, stat_col2_g002) %>%
  dplyr::select(comp, panel, n, estimate, lci, uci, p.value) %>%
  group_by(comp) %>%
  mutate(fdr = p.adjust(p.value, method='fdr')) %>%
  arrange(comp, panel) %>%
  ungroup()


write.csv(stat, file.path(STATS, 'figS19_pvals.csv'), row.names = FALSE)
# write.csv(stat, file.path(STATS, 'table_S40.csv'), row.names = FALSE)

q(save='no')
