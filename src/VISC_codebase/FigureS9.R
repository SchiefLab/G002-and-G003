#Figure S9
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggplot2)
library(DescTools)
library(latex2exp)

args <- commandArgs(trailingOnly = TRUE)
SOURCEDIR <- ifelse(length(args) >= 1, args[1], '../../G00x-plots')

# DIR = file.path(SOURCEDIR, 'Sup/S9')
# STATS = file.path(SOURCEDIR, 'Stats.codebase')

# path= "path-to-data"
path = SOURCEDIR
file_name= "2025-02-06G00X_manuscript_FigS9.csv"
full_path = file.path(path, file_name)
df= read.csv(full_path)

###################################################
time2= c(2,8,10,16)

varname1= c('auc_scharp','delta')
var1_name= c('AUTC','net MFI')

# antigen_select= c( 'eOD-GT8 60mer',
#                    'eOD-GT8.1_His-AvimC (CR479) 293F',
#                    'eOD-GT8.1_KO11(CR176) 293F',
#                    'Lumazine Synthase',
#                    'CD4bs')

antigen_select=c('eOD-GT8 60mer',
                 'eOD-GT8 monomer',
                 'eOD-GT8 KO',
                 'Lumazine Synthase',
                 'CD4bs')
antigen_label = antigen_select
names(antigen_label)= antigen_select

##############################
CCC_AUTC_plot= function(df, i2, time1, space1=0, wt=0){

  #function for AUTC truncated at 100
  antigen1= antigen_select[i2]
  df_sub = df %>% filter(outcome %in% antigen1, week %in% time1 )

  dat1 = pivot_wider( df_sub, names_from = labid, values_from = `response.magnitude` )
  dat1 = as.data.frame(dat1)

  yname= c('GT','KAVI')
  y1= log10(dat1[,yname[1]] ) #y1= duke
  y2= log10(dat1[,yname[2]] ) #y2= kavi

  ccc = DescTools::CCC(y1, y2, ci = "z-transform", conf.level = 0.95, na.rm =TRUE)

  tmp.mean <- mean( ccc$blalt$delta)
  tmp.sd <- sqrt(var( ccc$blalt$delta))

  #SpearmanRho
  #sc= SpearmanRho(y1, y2, use='complete.obs', conf.level= .95)

  ccc_res = ccc$rho.c
  ccc_res = ccc_res %>% mutate_if(is.numeric, round, 4)

  a= ccc_res$est
  ccc_lower= ccc_res$lwr.ci
  ccc_upper= ccc_res$upr.ci
  #b=sc['rho'] #precision

  #peaerson's rho
  pc= confintr::ci_cor(y1, y2, method = c("pearson"),type = c("normal"))

  min1=  min(c(y2)) #- .1
  max1=  ( max(c(y2))) #+.1
  min2=  min(c(y1)) # -.1
  max2=  ( max(c(y1))) #+.1
  diff= max1-min1
  diff2= max2-min2
  min0= min(min1,min2)
  min00 = trunc(min0)
  max0= max(y1,y2)

  dat2= dat1
  dat2$y1= y1
  dat2$y2= y2

  p1= ggplot( dat2, aes(x=y1, y=y2) )+
    geom_abline(  intercept = 0, slope = 1, color="gray60", linetype="dashed" , linewidth=.5) +
    stat_smooth(  method = "lm", se = FALSE , color='dodgerblue2', linewidth=.5 ) +
    geom_point( color=adjustcolor('gray30', alpha=.7), size=2 ) +
    #scale_colour_manual(name = "Visit (Week)", values = col1) +
    theme_bw() +
    coord_fixed(ratio = 1)+
    annotation_logticks( base = 10)+
    scale_y_continuous(#breaks = seq(1,3,by = 1),
      limits = c( min(min0,1.7), max(max1,max2,5.1) ) ,# c(  min0-.1-wt , max(max1, max2)+.1+wt ),
      #labels = scales::math_format(),
      breaks = c(2,3,4,5),
      labels = function(y) ifelse(y == 2, parse(text = TeX('$\\leq 10^2$')),  scales::math_format()(y) )
    ) +
    scale_x_continuous(#breaks = seq(1,3,by = 1),
      limits = c(min(min0,1.7), max(max1,max2,5.1)), #c(   min0-.1-wt , max(max1, max2)+.1+wt ),
      #labels = scales::math_format(),
      breaks = c(2,3,4,5),
      labels = function(x) ifelse(x == 2, parse(text = TeX('$\\leq 10^2$')),  scales::math_format()(x) ) ) +
    theme(plot.title = element_text(size = 14),
          plot.subtitle = element_text(size = 12),
          strip.text.x = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.title=element_text(size= 12),
          legend.text=element_text(size= 12) )

  if (i2==5){
    p1= p1 +labs(title=  antigen_label[i2] ,
                 subtitle = paste(  "CCC = ",format( round(a,4), nsmall=4), ', ',
                                    "Accuracy = ",format( round(ccc$C.b,4), nsmall=4) ,', ',
                                    "Precision = ",format( round(pc$estimate,4), nsmall=4), sep=''),
                 y =parse(text = TeX('$KAVI (\\Delta AUTC)$')), x =  TeX('$Duke (\\Delta AUTC)$'))
  } else{
    p1 = p1 +labs(title=  antigen_label[i2] ,
                  subtitle = paste(  "CCC = ",format( round(a,4), nsmall=4), ', ',
                                     "Accuracy = ",format( round(ccc$C.b,4), nsmall=4) ,', ',
                                     "Precision = ",format( round(pc$estimate,4), nsmall=4), sep=''),
                  y ="KAVI (AUTC)", x = "Duke (AUTC)")

  }

  return(p1)
}


plot_list2=list()
plot_list2[[1]]= CCC_AUTC_plot(df, i2=1, time1=time2, space1=.03, wt=.21)
plot_list2[[2]]= CCC_AUTC_plot(df, i2=2, time1=time2, space1=.03, wt=.21)
plot_list2[[3]]= CCC_AUTC_plot(df, i2=3, time1=time2, space1=.03, wt=.21)
plot_list2[[4]]= CCC_AUTC_plot(df, i2=4, time1=time2, space1=.03, wt=.21)
plot_list2[[5]]= CCC_AUTC_plot(df, i2=5, time1=time2, space1=.03, wt=.21)


layout= rbind( c(1,1,2,2), c(3,3,4,4), c(NA,5,5,NA))
pdf("pdfs/Reproducibility_postbaseline.pdf", width = 9, height = 13 )
grid.arrange(grobs = plot_list2, ncol = 2, nrow = 3,
             layout_matrix = layout)
#dev.off()
