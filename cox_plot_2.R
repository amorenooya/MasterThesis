library(dplyr)
library(glue)
#library(stringr)
library(ggplot2)
#library(cowplot)
#library(showtext)
#library(ggrepel)

setwd("/home/amorenoo/GWASP/NMIBC/EPICURO/output_models_all")
#numchr<-1:22


# file_recurrence<-glue('cox_recurrence_{numchr}_2.txt')
# recurrence<-lapply(file_recurrence, function(x) read.table(x, header=TRUE))
# recurrence<-do.call("rbind", recurrence) %>% as_tibble()
# 
# 
# 
# recurrence<-as_tibble(recurrence) %>%
#   rename(SNP=var) %>%
#   mutate(chr=word(SNP_new, 1, sep='_')) %>%
#   mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
#   relocate(chr, .after=SNP_new) %>%
#   mutate(chr=as.numeric(chr)) %>%
#   arrange(chr) %>%
#   mutate(highlight=case_when(p_value<1e-5~1,
#                              TRUE~0)) %>%
#   mutate(highlight=as.factor(highlight)) %>%
#   mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
#                       TRUE~'')) %>%
#   mutate(axis=1:nrow(recurrence))
# 
# breaks <- recurrence %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()
# 
# 
# 
theme_set(theme_minimal())
theme_update(plot.title=element_text(color='black', face='bold', size=30),#,
             #family='Poppins'),
             plot.title.position='plot',
             plot.subtitle=element_text(color='black', size=25),
             axis.title=element_text(color='black', size=15),
             #family='Poppins'),
             axis.text=element_text(size=10, color='black'),#,
             #family='Poppins'),
             plot.caption=element_text(size=20, color='black'),#,
             #family='Poppins'),
             legend.position='none',
             legend.justification='right',
             legend.title=element_text(size=50),#,
             # family='Poppins'),
             legend.text=element_text(size=35),
             #family='Poppins'),
             panel.grid = element_line(color="gray93", size=0.3))

# 
# g1<-ggplot(recurrence,
#            aes(x=axis, y=-log10(p_value),
#                color=as.factor(chr), label=rs)) +
#   geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
#   scale_color_manual(values=rep(c('#7d3c98',
#                                   '#ffcc00'),11)) +
#   scale_x_continuous(label=unique(recurrence$chr),
#                      breaks=breaks, expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
#   labs(x='Chromosome',
#        y='-log10(p-value)')+#,
#       # title='Recurrence free survival') +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x= element_blank())
# 
# 
# 
# tiff('recurrence_pvalues_notitle.tif',units='in', res=500,
#      compression='lzw',
#      width=15,
#      height=8)
# 
# g1
# 
# dev.off()
# 
# rm(recurrence)
# 
# 
# file_progression<-glue('cox_progression_{numchr}.txt')
# progression<-lapply(file_progression, function(x) read.table(x, header=TRUE))
# progression<-do.call("rbind", progression) %>% as_tibble()
# progression<-as_tibble(progression) %>%
#   rename(SNP=var) %>%
#   mutate(chr=word(SNP_new, 1, sep='_')) %>%
#   relocate(chr, .after=SNP_new)  %>%
#   mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
#   relocate(chr, .after=SNP_new) %>%
#   mutate(chr=as.numeric(chr)) %>%
#   arrange(chr) %>%
#   mutate(highlight=case_when(p_value<1e-5~1,
#                              TRUE~0)) %>%
#   mutate(highlight=as.factor(highlight)) %>%
#   mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
#                       TRUE~'')) %>%
#   mutate(axis=1:nrow(progression))
# 
# breaks <- progression %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()
# 
# 
# g2<-ggplot(progression,
#            aes(x=axis, y=-log10(p_value),
#                color=as.factor(chr), label=rs)) +
#   geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
#   scale_color_manual(values=rep(c('#7d3c98',
#                                   '#ffcc00'),11)) +
#   scale_x_continuous(label=unique(progression$chr),
#                      breaks=breaks, expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
#   labs(x='Chromosome',
#        y='-log10(p-value)') +#,
#        #title='Progression free survival') +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x= element_blank())
# 
# tiff('progression_pvalues_notile.tif',units='in', res=500,
#      compression='lzw',
#      width=15,
#      height=8)
# 
# g2
# 
# dev.off()

# 
# file_disease<-glue('cox_disease_{numchr}.txt') 
# disease<-lapply(file_disease, function(x) read.table(x, header=TRUE))
# disease<-do.call("rbind", disease) %>% as_tibble()
# 
# 

disease <- read.table('cox_disease_all.txt', header=TRUE, fill=TRUE)

disease<-as_tibble(disease) %>%
  #rename(SNP=var) %>%
  mutate(chr=word(SNP_new, 1, sep='_')) %>%
  mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
  relocate(chr, .after=SNP_new) %>%
  mutate(chr=as.numeric(chr)) %>%
  arrange(chr) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
                      TRUE~'')) %>%
  mutate(axis=1:nrow(disease))

breaks <- disease %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()


g1<-ggplot(disease,
           aes(x=axis, y=-log10(p_value),
               color=as.factor(chr), label=rs)) +
  geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
  scale_color_manual(values=rep(c('#7d3c98',
                                  '#ffcc00'),11)) +
  scale_x_continuous(label=unique(disease$chr),
                     breaks=breaks, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
  labs(x='Chromosome',
       y='-log10(p-value)') +#,
     #  title='Disease free survival') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())



tiff('relapse_all.tif',units='in', res=500,
     compression='lzw',
     width=15,
     height=8)

g1

dev.off()

rm(disease)