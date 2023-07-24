library(dplyr)
library(glue)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggrepel)

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


miami_plot <- function(outcome1, outcome2){
 result_1<-read.table(glue('cox_{outcome1}_all.txt'), #output from Cox_PH.R
                      header=TRUE, fill=TRUE)%>% as_tibble() 
 result_2<-read.table(glue('cox_{outcome2}_all.txt'), #output from Cox_PH.R
                      header=TRUE, fill=TRUE)%>% as_tibble() 
  
 result_1<-as_tibble(result_1) %>%
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
  mutate(axis=1:nrow(result_1)) %>%
  mutate(type='result_1')

  result_2<-as_tibble(result_2) %>%
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
  mutate(axis=1:nrow(result_2)) %>%
  mutate(type='result_2')
  
 result_2 <- rbind(result_2, result_1) %>% mutate(chr=as.numeric(chr))
 breaks <- result_2 %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()

 

g1<-ggplot(result_2,
           aes(x=axis, y=p_value, 
               color=as.factor(chr):as.factor(type), label=rs)) +
  geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
  scale_color_manual(values=rep(c('#7d3c98', 
                                  '#ffcc00', 
                                  '#693480',
                                  '#d8ae08'),11)) +
  scale_x_continuous(label=unique(result_2$chr),
                     breaks=breaks, expand=c(0,0)) +
  geom_label_repel(data=result_2 %>% filter(nombre==1 & highlight==1),
                                   size=2, max.overlaps = 55, color='black',
                                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
                                   segment.size=0.2, segment.linetype=2,
                                   min.segment.length = 0.5) +
  labs(x='Chromosome',
       y='-log10(p-value)')+
  scale_y_continuous(expand=c(0,0), 
                     limits = c(-9,9), breaks=c(-8:8),
                     labels = c('8','7','6','5','4','3','2', '1', '0', '1','2', '3','4', '5','6', '7','8'))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())

  return(g1)
  
}


miami_plot_withlabel <- function(outcome1, outcome2){
 result_1<-read.table(glue('cox_{outcome1}_all.txt'), #output from Cox_PH.R
                      header=TRUE, fill=TRUE)%>% as_tibble() 
 result_2<-read.table(glue('cox_{outcome2}_all.txt'), #output from Cox_PH.R
                      header=TRUE, fill=TRUE)%>% as_tibble() 
  
 result_1<-as_tibble(result_1) %>%
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
  mutate(axis=1:nrow(result_1)) %>%
  mutate(type='result_1')

  result_2<-as_tibble(result_2) %>%
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
  mutate(axis=1:nrow(result_2)) %>%
  mutate(type='result_2')
  
 result_2 <- rbind(result_2, result_1) %>% mutate(chr=as.numeric(chr))
 breaks <- result_2 %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()

  
 result_1_snps <- read.csv2(file=glue('{result_1}_FUMA.csv'), sep=',') #output from FUMA GWAS 
 result_2_snps <- read.csv2(file=glue('{result_2}_FUMA.csv'), sep=',') #output from FUMA GWAS

 result_2 <- result_2 %>% mutate(nombre=case_when(rs %in% result_2_snps$rsID~1,
                                                     rs %in% result_1_snps$rsID~1,
                                                     TRUE~0))
 

  g1<-ggplot(result_2,
           aes(x=axis, y=p_value, 
               color=as.factor(chr):as.factor(type), label=rs)) +
  geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
  scale_color_manual(values=rep(c('#7d3c98', 
                                  '#ffcc00', 
                                  '#693480',
                                  '#d8ae08'),11)) +
  scale_x_continuous(label=unique(result_2$chr),
                     breaks=breaks, expand=c(0,0)) +
  geom_label_repel(data=result_2 %>% filter(nombre==1 & highlight==1),
                                   size=2, max.overlaps = 55, color='black',
                                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
                                   segment.size=0.2, segment.linetype=2,
                                   min.segment.length = 0.5) +
  labs(x='Chromosome',
       y='-log10(p-value)')+
  scale_y_continuous(expand=c(0,0), 
                     limits = c(-9,9), breaks=c(-8:8),
                     labels = c('8','7','6','5','4','3','2', '1', '0', '1','2', '3','4', '5','6', '7','8'))+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())

  return(g1)
  
}
