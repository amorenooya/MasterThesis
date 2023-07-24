library(dplyr)
library(glue)
library(ggplot2)


setwd("/home/amorenoo/GWASP/NMIBC/EPICURO/output_models_all")
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

disease <- read.table('cox_disease_all.txt', header=TRUE, fill=TRUE)

disease<-as_tibble(disease) %>%
  #rename(SNP=var) %>%
  #mutate(chr=word(SNP_new, 1, sep='_')) %>%
  #mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
  #relocate(chr, .after=SNP_new) %>%
  #mutate(chr=as.numeric(chr)) %>%
  #arrange(chr) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(axis=1:nrow(disease))

disease <- disease %>%
  mutate(HR=exp(coef))

#breaks <- disease %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()


g1<-ggplot(disease,
           aes(y=-log10(p_value), x=coef,
               color=highlight)) +
  geom_point(alpha=0.3, width = 0.4, size=2) +
  scale_color_manual(values=c('#b2b2b2','#85b678')) +
  scale_x_continuous(limits=c(-3,3), breaks=c(-2,-1,0,1,2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
  labs(x='Beta coefficient',
       y='-log10(p-value)') +#,
  #  title='Disease free survival') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())



tiff('relapse_volcano.tif',units='in', res=500,
     compression='lzw',
     width=8,
     height=8)

g1

dev.off()



recurrence <- read.table('cox_recurrence_all.txt', header=TRUE, fill=TRUE)

recurrence<-as_tibble(recurrence) %>%
  #rename(SNP=var) %>%
  #mutate(chr=word(SNP_new, 1, sep='_')) %>%
  #mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
  #relocate(chr, .after=SNP_new) %>%
  #mutate(chr=as.numeric(chr)) %>%
  #arrange(chr) %>%
  mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
                      TRUE~'')) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(axis=1:nrow(recurrence)) 

recurrence_snps <- read.csv2(file='recurrence_FUMA_job237685_table.csv', sep=',')
progression_snps <- read.csv2(file='progression_FUMA_job238384_table.csv', sep=',')

recurrence <- recurrence %>% mutate(nombre=case_when(rs %in% recurrence_snps$rsID~1,
                                                     TRUE~0))


#breaks <- disease %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()


g1<-ggplot(recurrence,
           aes(y=-log10(p_value), x=coef, label=rs,
               color=highlight)) +
  geom_point(alpha=0.3, width = 0.4, size=2) +
  scale_color_manual(values=c('#b2b2b2','#ffcc00')) +
  scale_x_continuous(limits=c(-3,3), breaks=c(-2,-1,0,1,2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
  geom_label_repel(data=recurrence %>% filter(nombre==1 & highlight==1),
                   size=2, max.overlaps = 55, color='black',
                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
                   segment.size=0.2, segment.linetype=2,
                   min.segment.length = 0.5) +
  labs(x='ln(Hazard Ratio)',
       y='-log10(p-value)') +#,
  #  title='Disease free survival') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x=element_blank())



tiff('recurrence_volcano.tif',units='in', res=500,
     compression='lzw',
     width=8,
     height=7)

g1

dev.off()






progression <- read.table('cox_progression_all.txt', header=TRUE, fill=TRUE)

progression<-as_tibble(progression) %>%
  #rename(SNP=var) %>%
  #mutate(chr=word(SNP_new, 1, sep='_')) %>%
  #mutate(chr=readr::parse_number(chr)) %>% #extract the numerical value
  #relocate(chr, .after=SNP_new) %>%
  #mutate(chr=as.numeric(chr)) %>%
  #arrange(chr) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(axis=1:nrow(progression)) %>%
  mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
                                                       TRUE~''))

progression <- progression %>% mutate(nombre=case_when(rs %in% progression_snps$rsID~1,
                                                     TRUE~0))


#breaks <- disease %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()


g1<-ggplot(progression,
           aes(y=-log10(p_value), x=coef,label=rs,
               color=highlight)) +
  geom_point(alpha=0.3, width = 0.4, size=2) +
  scale_color_manual(values=c('#b2b2b2','#7d3c98')) +
  scale_x_continuous(limits=c(-3,3), breaks=c(-2,-1,0,1,2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
  geom_label_repel(data=progression %>% filter(nombre==1  & highlight==1),
                   size=2, max.overlaps = 55, color='black',
                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
                   segment.size=0.2, segment.linetype=2,
                   min.segment.length = 0.5) +
  labs(x='ln(Hazard Ratio)',
       y='-log10(p-value)') +#,
  #  title='Disease free survival') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x=element_blank())



tiff('progression_volcano.tif',units='in', res=500,
     compression='lzw',
     width=8,
     height=7)

g1

dev.off()
