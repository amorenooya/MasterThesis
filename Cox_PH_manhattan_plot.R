library(dplyr)
library(glue)
library(ggplot2)

theme_set(theme_minimal())
theme_update(plot.title=element_text(color='black', face='bold', size=30),
             plot.title.position='plot',
             plot.subtitle=element_text(color='black', size=25),
             axis.title=element_text(color='black', size=15)
             axis.text=element_text(size=10, color='black'),
             plot.caption=element_text(size=20, color='black'),
             legend.position='none',
             legend.justification='right',
             legend.title=element_text(size=50),
             legend.text=element_text(size=35),
             panel.grid = element_line(color="gray93", size=0.3))

#outcome: recurrence, progression or relapse

manhattan_plot<- function(outcome){
  
result <- read.table(glue('cox_{outcome}_all.txt'), header=TRUE, fill=TRUE) #output from Cox_PH.R
result<-as_tibble(result) %>%
  rename(SNP=var) %>%
  mutate(chr=word(SNP_new, 1, sep='_'),
         chr=readr::parse_number(chr)) %>% #extract the numerical value
  mutate(chr=as.numeric(chr)) %>%
  arrange(chr) %>%
  mutate(highlight=case_when(p_value<1e-5~1,
                             TRUE~0)) %>%
  mutate(highlight=as.factor(highlight)) %>%
  mutate(rs=case_when(str_detect(SNP, 'rs')~word(SNP, 1, sep='_'),
                      TRUE~'')) %>%
  mutate(axis=1:nrow(result))

breaks <- result %>% group_by(chr) %>% mutate(breaks=median(axis)) %>% ungroup() %>% pull(breaks) %>% unique()

g1<-ggplot(result,
           aes(x=axis, y=-log10(p_value),
               color=as.factor(chr), label=rs)) +
  geom_jitter(alpha=0.3, width = 0.4, size=1, stroke = 0.5) +
  scale_color_manual(values=rep(c('#7d3c98',
                                  '#ffcc00'),11)) +
  scale_x_continuous(label=unique(result$chr),
                     breaks=breaks, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0,9), 
                     breaks=c(0:8))+
  labs(x='Chromosome',
       y='-log10(p-value)') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x= element_blank())

  return(g1)
  
  }



  

