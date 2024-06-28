###
###   Benjamin Harrison
###


# variable importance and intersection in the A.Type and C.Type clocks of Hubert, Greenspan, et al. Phillip's lab
# install.packages('ggvenn')
library(ggvenn)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)

rm(list=ls())



b <- read.csv("06_19_24_raw_features.csv")
head(b)
b <- b[-c(1), -c(1)]
head(b)
colnames(b)[1] <- 'mz'

table(b$A.Type==0)
table(b$C.Type==0)
table(b$C.Type!=0, b$A.Type!=0) # 9 in common

fisher.test(table(b$C.Type!=0, b$A.Type!=0) ) # the intersecting metabolites are not an over/under-representation of the total number of features. Fisher test odds ratio 2.69 P=0.081
 
###################
# venn diagram
x <- list('A clock'=b$mz[b$A.Type!=0], 'C clock'=b$mz[b$C.Type!=0])
ggvenn(x, show_percentage = F, text_size = 10, set_name_size = 10, auto_scale=T)

intersecting.features <- intersect(x[[1]], x[[2]])
write.table(x[[1]][!x[[1]]%in% intersecting.features], quote=F, file='A features', row.names = F, col.names = F)
write.table(x[[2]][!x[[2]]%in% intersecting.features], quote=F, file='C features', row.names = F, col.names = F)
write.table(intersecting.features, quote=F, file='intersecting features', row.names = F, col.names = F)


####################

## plots of betas:
b$intersecting <- ifelse(rowSums(b[ ,2:3]==0)==0, 'intersecting', 'unique')

ggplot(b, aes(y=A.Type, x=C.Type, color=intersecting))+
  geom_point()+
  theme_bw(base_size = 16)+
  scale_color_manual(values=c('black', 'grey'))

b[rowSums(b[ ,2:3]==0)==0, ]

p1 <- ggplot(b[rowSums(b[ ,2:3]==0)==0, ], aes(y=A.Type, x=C.Type, label=mz))+
  geom_point()+
  theme_bw(base_size = 16)+
  ylab(expression(beta * ' A clock'))+
  xlab(expression(beta * ' C clock'))+
  ggrepel::geom_text_repel()

p1  

ggplot(b, aes(y=A.Type, x=C.Type, color=intersecting, label=mz))+
  geom_point()+
  theme_bw(base_size = 18)+
  scale_color_manual(values=c('black', 'grey40'))+
  ylab(expression(beta * ' A-clock'))+
  xlab(expression(beta * ' C-clock'))+
  ggrepel::geom_text_repel(size=3, max.overlaps = 8)+
  facet_wrap(~intersecting)+
  theme(legend.position = 'none')


b[b==0] = NA
cor.test(b$A.Type, b$C.Type)

# are the intersecting metabolites any more or less imporant (magnitude of beta) than the unique metabolties in either clock?
head(b)

colnames(b)[2:3] <- c('A clock', 'C clock')
l <- b %>% pivot_longer(-c(mz, intersecting), values_to = 'coef', names_to = 'clock')
head(l)

wilcox.test(abs(coef) ~ intersecting, l[l$clock=='A clock', ]) # nope
wilcox.test(abs(coef) ~ intersecting, l[l$clock=='C clock', ]) # nope

p2 <- ggplot(l, aes(y=coef, x=intersecting, color=intersecting))+
  geom_jitter(width=0.05)+
  theme_bw(base_size = 16)+
  ylab(expression(beta))+
  facet_wrap(~clock, scales='free')+
  theme(axis.text.x=element_blank())+
  xlab('feature type')+ 
  labs(color='feature type') 

p2



ggarrange(p1, p2)



