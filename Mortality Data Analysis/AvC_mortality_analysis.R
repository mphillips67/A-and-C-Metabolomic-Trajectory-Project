dat<-read.csv("Winter 2023 Mortality.csv",header=TRUE)

#change for different sex
sex<-"M"

#change for convergence
conv<-"A"

#add intervals to full data
dat.int<-cbind(dat,ceiling((dat$Age-1)/3))
colnames(dat.int)[10]<-"Interval"
dat.int<-dat.int[!dat.int$Interval==0 & !dat.int$Alive<=20,]
dat.int$AgeSpecficMort<-log(dat.int$AgeSpecficMort)
dat.int$AgeSpecficMort<-as.numeric(gsub(-Inf,-7.5,as.numeric(dat.int$AgeSpecficMort)))

#Convergence testing

#### USE FOR CONVERGENCE WITHIN TREATMENT ####
#dat.con<-dat.int[dat.int$Sex==sex & dat.int$SimpleSel==conv,]
####                                      ####

#### USE FOR DIVERGENCE BETWEEN TREATMENT ####
dat.con<-dat.int[dat.int$Sex==sex,]
####                                      ####


mini.set<-NULL
for (i in unique(dat.con$Population)){
  mini.set<-c(mini.set,max(dat.con[dat.con$Population==i,]$Interval))
}
mini<-min(mini.set)
fin.dat.con<-dat.con[dat.con$Interval<=mini,]

#analysis
library(nlme)
conv.dat<-groupedData(AgeSpecficMort~Age|Population,data=fin.dat.con)
conv.dat$Interval<-as.factor(conv.dat$Interval)

#### USE FOR CONVERGENCE WITHIN TREATMENT ####
#conv.lme<-lme(AgeSpecficMort~Age*Interval + Interval*Selection,data=conv.dat,random=~1)
####                                      ####

#### USE FOR DIVERGENCE BETWEEN TREATMENT ####
conv.lme<-lme(AgeSpecficMort~Age*Interval + Interval*SimpleSel,data=conv.dat,random=~1)
####                                      ####


summary(conv.lme)

#standard error
stde<-summary(conv.lme)$tTable[(mini+2),2]
for (i in 0:(mini-2)){
  s<-conv.lme$varFix[(mini+2),(mini+2)]
  a<-conv.lme$varFix[(i+((2*mini)+2)),(i+((2*mini)+2))]
  b<-conv.lme$varFix[(mini+2),(i+((2*mini)+2))]
  stde<-c(stde,sqrt(s+a-2*b))
}

#intercepts
type1.inters<-fixed.effects(conv.lme)[1]
for (i in 1:(mini-1)){
  type1.inters<-c(type1.inters,(fixed.effects(conv.lme)[1]+fixed.effects(conv.lme)[i+2]))
}
type2.inters<-fixed.effects(conv.lme)[1]+fixed.effects(conv.lme)[mini+2]
for (i in 1:(mini-1)){
  type2.inters<-c(type2.inters,type1.inters[i+1]+fixed.effects(conv.lme)[mini+2]+fixed.effects(conv.lme)[i+(2*mini)+1])
}

#p-values
p.done<-data.frame(1:mini,rep(0,mini))
colnames(p.done)<-c("Interval","p-value")
for (i in 1:mini){
  p.done[i,2]<-(1-pt(abs((type1.inters[i]-type2.inters[i]))/stde[i],8))*2
}
write.csv(p.done,"table.csv",row.names=FALSE)

