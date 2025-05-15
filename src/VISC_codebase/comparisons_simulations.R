# Run parametric bootstrap method

library(tidyverse)
library(VISCfunctions)

# Accept SOURCEDIR as a command line argument, default to '../../G00x-plots'
args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')
DIR = file.path(SOURCEDIR, '/Main-Metrics/fig2')
DIR3 = file.path(SOURCEDIR, '/Main-Metrics/fig3')
CDIR = file.path(SOURCEDIR, '/Sup-Metrics/S19')
STATS = './Stats.codebase'

# Validate if the folders exist
if (!dir.exists(SOURCEDIR)) {
  stop(paste("Source directory does not exist:", SOURCEDIR))
}

if (!dir.exists(DIR)) {
  stop(paste("Main metrics fig2 directory does not exist:", DIR))
}

if (!dir.exists(DIR3)) {
  stop(paste("Main metrics fig3 directory does not exist:", DIR3))
}

if (!dir.exists(CDIR)) {
  stop(paste("Supplementary metrics S19 directory does not exist:", CDIR))
}

if (!dir.exists(STATS)) {
  stop(paste("Stats codebase directory does not exist:", STATS))
}

# DIR = '~/Repo/Schief856_G002/data/dropbox/Main-Metrics/fig2'
# DIR3 = '~/Repo/Schief856_G002/data/dropbox/Main-Metrics/fig3'
# CDIR = '~/Repo/Schief856_G002/data/dropbox/Supp-Metrics/S19-alt'
# STATS = '~/Repo/Schief856_G002/data/dropbox/Stats.codebase'


# create panel map between figures 2 and 3 and S19
# list names are figure 2 or 3 panels and value is S19 panels
# panel_map = list(
#   'Fig. 2A' = c('A', 'G'),
#   'Fig. 2B' = c('B', 'H'),
#   'Fig. 2C' = c('C', 'I'),
#   'Fig. 2E' = c('D', 'J'),
#   'Fig. 3A' = c('E', 'K'),
#   'Fig. 3B' = c('F', 'L')
# )

panel_map = list(
  # 'Fig. 2A' = c('A', 'F'),
  # 'Fig. 2B' = c('B', 'G'),
  # 'Fig. 2C' = c('', ''),
  'Fig. 2E' = c('C', 'H'),
  'Fig. 3A' = c('D', 'I'),
  'Fig. 3B' = c('E', 'J')
)

# get method data
fn = list.files(CDIR)

# read individual files
dat = lapply(fn, function(n) {
  df = read.csv(file.path(CDIR, n), check.names = FALSE)

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
# stopifnot(length(unique(dat$pubID))==10)
stopifnot(all(dat$trial == dat$method))
dat$trial = NULL

# check for NAs
stopifnot(!any(is.na(dat$y)))

# check for zeros in data that will be log transformed panels A, B, D, G, H, and J
w.log = which(dat$panel %in% c('A','B','C','F','G','H'))
stopifnot(all(dat$y[w.log] > 0))

# Remove duplicated data.  In figure S19 panels G-L duplicate data in panels A-H
# for method G001 and G002 so only the G003 data is needed from those panels
dat = rbind(subset(dat, panel %in% c('A', 'B', 'C', 'D', 'E')),
            subset(dat, panel %in% c('F', 'G', 'H', 'I', 'J') & method=='G003'))

# convert panel names to main figure panel names
for( i in 1:length(panel_map) ) {
  dat$panel[ dat$panel %in% panel_map[[i]] ] = names(panel_map)[i]
}

# rename dat to cdat
cdat = dat
dat = NULL

### Get the figure 2 data and combine

# The data for each panel is in a separate file.  The list `map` maps
# the panel label (A, B, C, etc.) to the file name
map = list(
  a = "figA_percent_gt8$^{++}$_amongigg$^+$_b_cells.csv",
  b = "figB_percent_cd4bs-specific_amongigg$^{+}$_b_cells.csv" ,
  c = "figC_percent_of_gt8$^{++}$igg$^{+}$b_cells_that_are_ko$^{-}$.csv",
  d = "figD_number_of_vrc01-classigg$^{+}$_b_cells.csv",
  e = "figE_percent_of_igg$^{+}$_memory_b_cellsdetected_as_vrc01-class.csv",
  f = "figF_percent_vrc01-class_amongcd4bs-specificigg$^{+}$_memory_sequences.csv",
  g = "figG_percent_vrc01-class_amonggt8-specificigg$^{+}$_memory_sequences.csv",
  h = "figH_percent_of_vrc01-class_responders.csv"
)

# read in the data files
dat = lapply(names(map), function(panel) {
  read.csv(file.path(DIR, map[[panel]]), check.names = FALSE)
})
names(dat) = names(map)

# process each file so that we can rbind them
# the data are trial, pubID, weeks and y-axis label name
# convert to a format where I can rbind the files and use
# outcome to process the comparisons.

# NOTE, panel H has a different column order
dat = lapply( names(dat), function(panel) {
  dat_panel = dat[[panel]]
  if(panel=='h') {
    stopifnot(all(names(dat_panel) == c("pubID", "trial", "percent_of_vrc01-class_responders", "weeks")))
    dat_panel = dat_panel[,c("trial", "pubID", "weeks", "percent_of_vrc01-class_responders")]
  }
  stopifnot(all(names(dat_panel)[1:3] == c('trial', 'pubID', 'weeks')))
  oc = names(dat_panel)[4]
  dat_panel$panel = toupper(panel)
  dat_panel$y_label = oc
  dat_panel$y = dat_panel[[oc]]
  dat_panel[[oc]] = NULL

  return(dat_panel)
})

# rbind the data files
dat = do.call(rbind, dat)


### Get the figure 3 data and combine

# The data for each panel is in a separate file.  The list `map` maps
# the panel label (A, B, etc.) to the file name
map3 = list(
  a = "figA_prime_v_heavy_aa_mutation.csv",
  b = "figB_prime_v_light_aa_mutation.csv",
  e = "figE_90th_percentile_number_of_key_VRC01-class_HC_residues.csv"
)

# read in the data files
dat3 = lapply(names(map3), function(panel) {
  read.csv(file.path(DIR3, map3[[panel]]), check.names = FALSE)
})
names(dat3) = names(map3)

# process each file so that we can rbind them
# the data are trial, pubID, weeks and y-axis label name
# convert to a format where I can rbind the files and use
# outcome to process the comparisons.
dat3 = lapply( names(dat3), function(panel) {
  dat_panel = dat3[[panel]]
  stopifnot(all(names(dat_panel)[1:4] == c('pubID', 'trial', 'pseudogroup', 'weeks')))
  oc = names(dat_panel)[5]
  if( oc == 'cottrell_focused_v_common_score' ) {
    # remove duplicated column
    stopifnot(all(dat_panel$cottrell_focused_v_common_score == dat_panel$residues))
    dat_panel$residues = NULL
  }
  dat_panel$panel = toupper(panel)
  dat_panel$y_label = oc
  dat_panel$y = dat_panel[[oc]]
  dat_panel[[oc]] = NULL

  return(dat_panel)
})

# rbind the data files
dat3 = do.call(rbind, dat3)

# combine figure 2 A, B, C, and E data with figure 3 A and B
# dat = subset(dat, panel %in% c('A', 'B', 'C', 'E'))
dat = subset(dat, panel %in% c('A', 'B', 'E'))
dat$panel = paste0('Fig. 2', dat$panel)

dat3 = subset(dat3, panel %in% c('A', 'B'))
dat3$pseudogroup = NULL
dat3$panel = paste0('Fig. 3', dat3$panel)

dat = rbind(dat, dat3)
dat$panel = factor(dat$panel)


### run analysis

# comparisons
comp = data.frame( g1 = c('G001', 'G001', 'G002'),
                   g2 = c('G002', 'G003', 'G003'))

# store the results in res
res = c()

for( i in 1:nrow(comp) ) {
  comp_name = sprintf('%s vs. %s', comp$g1[i], comp$g2[i])

  for( p in levels(dat$panel) ) {
    df = dat %>%
      filter( panel==p ) %>% # filter for panel
      filter( trial %in% c(comp$g1[i], comp$g2[i]) )
    dfc = cdat %>%
      filter( panel == p ) %>%
      filter( method %in% c(comp$g1[i], comp$g2[i]) )

    if(p %in% c('Fig. 2A', 'Fig. 2B', 'Fig. 2E')) {

      # stopifnot(all(df$y > 0))
      df = subset(df, y>0)

      stopifnot(all(dfc$y > 0))

      df$y = log10(df$y)
      dfc$y = log10(dfc$y)
    }


    # get paired concordance data for these two trials and this panel
    df.conc = pivot_wider(dfc, id_cols = pubID, names_from = method, values_from = y)

    v = df.conc[[comp$g2[i]]] - df.conc[[comp$g1[i]]]
    v = v[!is.na(v)] # NAs come up because the G003 method was run on a subset of samples
    n.d = length(v)
    md = mean(v)
    sdd = sd(v)

    out = c()

    # fit repeatedly for each time point
    for( wk in c(-5, 4, 8, 10, 16) ) {
      if( grepl('G002', comp_name) && wk==10 ) next
      if( grepl('G003', comp_name) && wk==4 ) next
      if( wk==-5 && p %in% c('Fig. 2E', 'Fig. 3A','Fig. 3B') ) next

      ss1 = subset(df, weeks==wk & trial==comp$g1[i] & !is.na(y))
      ss2 = subset(df, weeks==wk & trial==comp$g2[i] & !is.na(y))

      N=10000

      n1 = nrow(ss1)
      mu1 = mean(ss1$y)
      sd1 = sd(ss1$y)
      n2 = nrow(ss2)
      mu2 = mean(ss2$y)
      sd2 = sd(ss2$y)

      d_obs = mu2 - mu1 - md
      lci = uci = pval = NA

      if( n1 >=2 & n2 >=2 ) {
        set.seed(1)
        d_samp = replicate(N, {
          d1 = rnorm(n1, mu1, sd1)
          d2 = rnorm(n2, mu2, sd2)
          m = md + sdd*rt(n.d, n.d-1)

          mean(d2) - mean(d1) - mean(m)
        })

        pval = (2*min(sum(d_samp <0 ), sum(d_samp > 0)) + 0.5)/(N + 1)

        lci = quantile(d_samp, 0.025)
        uci = quantile(d_samp, 0.975)
      }

      out = rbind(out, data.frame(comp=comp_name,
                                  panel=p,
                                  week = wk,
                                  n1 = n1,
                                  n2 = n2,
                                  d = d_obs,
                                  lci=lci,
                                  uci=uci,
                                  p.value=pval,
                                  d.m = md))
    }

    # backtransform estimates for log-transformed outcomes
    if(p %in% c('Fig. 2A', 'Fig. 2B', 'Fig. 2E')) {
      out$d = 10^out$d
      out$lci = 10^out$lci
      out$uci = 10^out$uci
    }


    res = rbind(res, out)
  }
}

# add fdr adjusted p-values by panel
res = res %>%
  group_by(panel) %>%
  mutate(fdr = p.adjust(p.value, 'fdr')) %>%
  ungroup() %>%
  arrange(panel, comp)


# write to file only creating file with all g00X comparisons
write.csv(res, file.path(STATS, 'simulations_pvals_g00X.csv'), row.names = FALSE)

q(save='no')
