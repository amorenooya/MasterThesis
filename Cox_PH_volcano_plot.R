library(dplyr)
library(glue)
library(ggplot2)

theme_set(theme_minimal())
theme_update(plot.title=element_text(color='black', face='bold', size=30),
             plot.title.position='plot',
             plot.subtitle=element_text(color='black', size=25),
             axis.title=element_text(color='black', size=15),
             axis.text=element_text(size=10, color='black'),
             plot.caption=element_text(size=20, color='black'),
             legend.position='none',
             legend.justification='right',
             legend.title=element_text(size=50),
             legend.text=element_text(size=35),
             panel.grid = element_line(color="gray93", size=0.3))


volcano_plot<- function(outcome){

result <- read.table(glue('cox_{outcome}_all.txt'), 
                     header=TRUE, fill=TRUE) #result from Cox PH

result<-as_tibble(result) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(axis=1:nrow(result)) 

g1<-ggplot(result,
           aes(y=-log10(p_value), x=coef,
               color=highlight)) +
  geom_point(alpha=0.3, width = 0.4, size=2) +
  scale_color_manual(values=c('#b2b2b2','#85b678')) +
  scale_x_continuous(limits=c(-3,3), breaks=c(-2,-1,0,1,2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
  labs(x='Beta coefficient',
       y='-log10(p-value)') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())

  return(g1)
}



volcano_plot_withlabel<- function(outcome){

result <- read.table(glue('cox_{outcome}_all.txt'), 
                     header=TRUE, fill=TRUE) #result from Cox PH
  
result_snps <- read.csv2(file=glue('{result_1}_FUMA.csv'), sep=',') #output from FUMA GWAS 

result <- result %>% mutate(nombre=case_when(rs %in% result_2_snps$rsID~1,
                                             rs %in% result_1_snps$rsID~1,
                                             TRUE~0))

result<-as_tibble(result) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(axis=1:nrow(result)) 

g1<-ggplot(result,
           aes(y=-log10(p_value), x=coef,
               color=highlight)) +
  geom_point(alpha=0.3, width = 0.4, size=2) +
  scale_color_manual(values=c('#b2b2b2','#85b678')) +
  scale_x_continuous(limits=c(-3,3), breaks=c(-2,-1,0,1,2)) +
  scale_y_continuous(expand = c(0, 0),limits = c(0,9), breaks=c(0:8))+
  geom_label_repel(data=recurrence %>% filter(nombre==1 & highlight==1),
                   size=2, max.overlaps = 55, color='black',
                   fill=rgb(red = 1, green = 1, blue = 1, alpha = 0.6),
                   segment.size=0.2, segment.linetype=2,
                   min.segment.length = 0.5) +
  labs(x='Beta coefficient',
       y='-log10(p-value)') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())

  return(g1)
}

