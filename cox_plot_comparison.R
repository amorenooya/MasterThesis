library(dplyr)
library(glue)
library(stringr)
library(ggplot2)
library(cowplot)
#library(showtext)
library(ggrepel)

setwd("/home/amorenoo/GWASP/NMIBC/EPICURO/output_models_all")
#numchr<-1:22

#file_recurrence<-glue('cox_recurrence_all.txt')  
#file_lowrisk<-glue('cox_disease_lowrisk_{numchr}.txt')
#recurrence<-lapply(file_recurrence, function(x) read.table(x, header=TRUE))
#recurrence<-do.call("rbind", recurrence) 

recurrence<-read.table('cox_recurrence_all.txt', header=TRUE, fill=TRUE)%>% as_tibble()

# 
#file_highrisk<-glue('cox_disease_highrisk_{numchr}.txt')
#highrisk<-lapply(file_highrisk, function(x) read.table(x, header=TRUE))
#highrisk<-do.call("rbind", highrisk) %>% as_tibble()

# 
# file_progression<-glue('cox_progression_{numchr}.txt')
# progression<-lapply(file_progression, function(x) read.table(x, header=TRUE))
# progression<-do.call("rbind", progression) %>% as_tibble()



progression<-read.table('cox_progression_all.txt', header=TRUE, fill=TRUE)%>% as_tibble()

#relapse_lowrisk <- read.table('cox_disease_lowrisk_all.txt', header=TRUE, fill=TRUE)%>% as_tibble()
#relapse_highrisk <- read.table('cox_disease_highrisk_all.txt', header=TRUE, fill=TRUE)%>% as_tibble()

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


# relapse_lowrisk<-as_tibble(relapse_lowrisk) %>%
#   #rename(SNP=var) %>%
#   mutate(chr=word(SNP_new, 1, sep='_')) %>%
#   mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
#   relocate(chr, .after=SNP_new) %>%
#   mutate(highlight=case_when(p_value<1e-5~1,
#                              TRUE~0)) %>%
#   mutate(highlight=as.factor(highlight)) %>%
#   mutate(p_value=log10(p_value)) %>%
#   mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
#                       TRUE~'')) %>%
#   mutate(chr=as.numeric(chr)) %>%
#   arrange(chr) %>%
#   mutate(axis=1:nrow(relapse_lowrisk)) %>%
#   mutate(type='Low_risk')
# 
# 
# relapse_highrisk<-as_tibble(relapse_highrisk) %>%
#   #rename(SNP=var) %>%
#   mutate(chr=word(SNP_new, 1, sep='_')) %>%
#   mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
#   relocate(chr, .after=SNP_new) %>%
#   mutate(highlight=case_when(p_value<1e-5~1,
#                              TRUE~0)) %>%
#   mutate(highlight=as.factor(highlight)) %>%
#   mutate(p_value=-log10(p_value)) %>%
#   mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
#                       TRUE~'')) %>%
#   mutate(chr=as.numeric(chr)) %>%
#   arrange(chr) %>%
#   mutate(axis=1:nrow(relapse_highrisk)) %>%
#   mutate(type='High_risk')
# 



progression<-as_tibble(progression) %>%
  #rename(SNP=var) %>%
  mutate(chr=word(SNP_new, 1, sep='_')) %>%
  mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
  relocate(chr, .after=SNP_new) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(p_value=log10(p_value)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
                      TRUE~'')) %>%
  mutate(chr=as.numeric(chr)) %>%
  arrange(chr, position) %>%
  mutate(axis=1:nrow(progression)) %>%
  mutate(type='Progression')


recurrence<-as_tibble(recurrence) %>%
  #rename(SNP=var) %>%
  mutate(chr=word(SNP_new, 1, sep='_')) %>%
  mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
  relocate(chr, .after=SNP_new) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(p_value=-log10(p_value)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
                      TRUE~'')) %>%
  mutate(chr=as.numeric(chr)) %>%
  arrange(chr, position) %>%
  mutate(axis=1:nrow(recurrence)) %>%
  mutate(type='recurrence')

#relapse_highrisk<-rbind(relapse_highrisk, relapse_lowrisk) %>% mutate(chr=as.numeric(chr))

recurrence <- rbind(recurrence, progression) %>% mutate(chr=as.numeric(chr))

recurrence_snps <- read.csv2(file='recurrence_FUMA_job237685_table.csv', sep=',')
progression_snps <- read.csv2(file='progression_FUMA_job238384_table.csv', sep=',')

recurrence <- recurrence %>% mutate(nombre=case_when(rs %in% recurrence_snps$rsID~1,
                                                     rs %in% progression_snps$rsID~1,
                                                     TRUE~0))


# progression<- rbind(progression, lowrisk) %>%
#   mutate(chr=as.numeric(chr))

# progression<- rbind(progression, highrisk) %>%
#   mutate(chr=as.numeric(chr))
# # 
# 
# lowrisk<-rbind(lowrisk, highrisk) %>%
#   mutate(chr=as.numeric(chr))

#breaks <- lowrisk %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()

#breaks <- relapse_highrisk %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()
breaks <- recurrence %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()

# 
# g1<-ggplot(progression,
#            aes(x=axis, y=p_value, 
#                color=as.factor(chr):as.factor(type), label=rs)) +
#   geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
#   scale_color_manual(values=rep(c('#7d3c98', 
#                                   '#ffcc00', 
#                                   '#693480',
#                                   '#d8ae08'),11)) +
#   scale_x_continuous(label=unique(progression$chr),
#                      breaks=breaks) +
#   labs(x='Chromosome',
#        y='-log10(p-value)',
#        title='Progression free survival',
#        subtitle='All vs high risk') +
#   geom_hline(yintercept=5, colour='#2471a3')+
#   geom_hline(yintercept=-5, colour='#2471a3')+
#   scale_y_continuous(expand=c(0,0), 
#                      limits = c(-9,9), breaks=c(-8,-6,-4,-2,0,2,4,6,8),
#                      labels = c('8','6','4','2','0','2','4','6','8'))+
#   # geom_label_repel(data=lowrisk %>% filter(highlight==1),
#   #                  size=1.5, max.overlaps = 55, color='black',
#   #                  fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
#   #                  segment.size=0.2, segment.linetype=2,
#   #                  min.segment.length = 0.5) +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.x= element_blank())




g1<-ggplot(recurrence,
           aes(x=axis, y=p_value, 
               color=as.factor(chr):as.factor(type), label=rs)) +
  geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
  scale_color_manual(values=rep(c('#7d3c98', 
                                  '#ffcc00', 
                                  '#693480',
                                  '#d8ae08'),11)) +
  scale_x_continuous(label=unique(recurrence$chr),
                     breaks=breaks, expand=c(0,0)) +
  geom_label_repel(data=recurrence %>% filter(nombre==1 & highlight==1),
                                   size=2, max.overlaps = 55, color='black',
                                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
                                   segment.size=0.2, segment.linetype=2,
                                   min.segment.length = 0.5) +
  labs(x='Chromosome',
       y='-log10(p-value)')+#,
       #title='Progression free survival',
       #subtitle='All vs low risk') +
  scale_y_continuous(expand=c(0,0), 
                     limits = c(-9,9), breaks=c(-8:8),
                     labels = c('8','7','6','5','4','3','2', '1', '0', '1','2', '3','4', '5','6', '7','8'))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())




# 
# 
# tiff('progression_pvalues_allvslow.tif',units='in', res=500,
#      compression='lzw',
#      width=15,
#      height=8)
# 
# g1
# 
# dev.off() 
# 
# 

tiff('recurrence_vs_progression_label.tif',units='in', res=500,
     compression='lzw',
     width=15,
     height=8)

g1

dev.off()


# tiff('relapse_high_vs_low.tif',units='in', res=500,
#      compression='lzw',
#      width=15,
#      height=8)
# 
# g1
# 
# dev.off() 
