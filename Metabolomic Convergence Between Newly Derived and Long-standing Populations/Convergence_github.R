
library(ggplot2)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(pheatmap)
library(svglite)


rm(list=ls())

d <- read.csv('A-C Normalized Metabolomics Time Series Data.csv')
#change col head SimpleSel to regime
colnames(d)[colnames(d) == "SimpleSel"] <- "regime"
d$regime <- gsub("-Type", "-type", d$regime)
head(d)
str(d)

mzs <- colnames(d)[8:ncol(d)]

d$Pop <- as.factor(d$Pop)
d$Day
d$Sel <- as.factor(d$Sel)
d$regime <- as.factor(d$regime)
d$Batch <- as.factor(d$Batch)

boxplot(d[ ,mzs])
boxplot(t(d[ ,mzs]), col=d$regime)
boxplot(t(d[ ,mzs]), col=d$Day)

boxplot(t(d[order(d$Day), mzs]), col=d$Day[order(d$Day)]) # variance by biological variable/confounded with sampling session (a technical batching). 


pheatmap(cor(d[ ,mzs]), fontsize=5)
pheatmap(cor(t(d[ ,mzs])), fontsize=5)


# AO and nCO having been derived more recently than ACO and CO
d$history <- factor(ifelse(d$Sel %in% c('ACO', 'CO'), 'longstanding', 'recent'), levels=c('longstanding', 'recent'))
d$history

colnames(d) %in%mzs
l <- pivot_longer(d, all_of(mzs), names_to = 'mz')
head(l)

table(l$Day)

# in these plots, each point is a replicate, so, 5 reps per 4 pops per age
ggplot(subset(l, Day%in% c(21, 28, 35)), aes(y=value, x=as.factor(Day), color=mz, group=mz))+
  geom_jitter(width=0.05)+
  geom_line(alpha=0.2)+
  theme_classic(base_size = 14)+
  facet_wrap(~Sel)+
  theme(legend.position = 'none')

table(l$mz, l$Day, l$Sel)

# in these plots, each line connects the mean of the 5 reps per pops per age, points removed to clarify.
ggplot(subset(l, Day%in% c(21, 28, 35)), aes(y=value, x=as.factor(Day), color=mz, group=mz))+
  stat_summary(fun = "mean", geom = "line", alpha=0.5)+
  facet_grid(~regime ~history)+
  ylab('mean metabolite (n=5)')+
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')

# variance in each metabolite by selection regime & history & day
o <- list()
for(i in 1:length(mzs)){
o[[i]] <- aggregate(d[ ,mzs[i]] ~ Day * regime * history, d, var) 
names(o[[i]])[4] <- 'var'
o[[i]]$mz <- mzs[i] }
v <- do.call(rbind, o)

head(v)

# one line per mz, lines connect the variance among 5 reps by age
ggplot(subset(v, Day%in% c(21, 28, 35)), aes(y=var, x=as.factor(Day), color=regime, group=mz))+
  geom_line(alpha=0.5)+
  facet_grid(~regime ~history)+
  ylab('metabolite variance (n=5)')+
  theme_bw(base_size = 14)+
  theme(legend.position = 'none')
# I would not say that there is strong evidence of inc. var within an age among the A and C, nor a diff btw longstanding and recent.  Although it does appear that the A type may have higher var over all ages.


## use PCA to reveal relationships that might exist for the metabolome with age, regime (A and C), history of selection, and the consistency/inconsistency of replicates:
pca <- prcomp(d[ ,mzs], scale=T)

plot(pca)
summary(pca)

loadings <- as.data.frame(pca$rotation)

eigs <- pca$sdev^2
o <- AssocTests::tw(eigs, eigenL=length(eigs), criticalpoint=0.9793) # alpha at 5%
sigEigs <- print(o$SigntEigenL) # there are 12 PCs that explain/detect variation above that expected from random data at an alpha of <5%
PCs <- as.data.frame(pca$x[ ,1:sigEigs])
pcs <- colnames(PCs)

head(PCs)

pca <- data.frame(Pop=d$Pop, Day=d$Day, Sel=d$Sel, regime=d$regime, history=d$history, PCs)

l <- pivot_longer(pca, all_of(pcs), names_to = 'PC')
str(l)

l$Pop
l$Sel
l$regime
l$Sel <- factor(l$Sel, levels=c('ACO', 'AO', 'nCO', 'CO'))
l$PC <- factor(l$PC, levels=pcs)

ggplot(l, aes(y=value, x=Sel, fill=history, color=regime))+
  geom_boxplot(alpha=0.5, outliers = F)+
  geom_jitter(width=0.1, size=0.8)+
  theme_classic()+
  scale_fill_manual(values=c('white', 'grey'))+
  facet_wrap(~PC, scales='free')+
  ylab('')+
  ggtitle('these plots do not consider age')

head(l)
PCApanels <- ggplot(subset(l, Day %in% c(21, 28, 35) & PC %in% pcs[1:sigEigs] ), aes(y=value, x=as.factor(Day), color=regime, group=Pop))+
  geom_line(alpha=0.5, aes(linetype=history))+
  geom_point(aes(shape=history), size = 3)+
  scale_shape_manual(values=c(1, 19))+
  scale_color_manual(values=c("A-type" = "red", "C-type" = "blue")) +  # Color coding
  theme_bw(base_size = 18)+
  facet_wrap(~PC)+
  ylab('')+
  xlab('age from egg (days)')


PCApanels

#ggsave('plots/PCApannels.svg', PCApanels)

# this version uses separate Y-axis scales for each panel, this can allow you to see into the data that has less variance (which by definition are the higher PCS), but can be a little misleading:
ggplot(subset(l, Day %in% c(21, 28, 35) & PC %in% pcs[1:sigEigs] ), aes(y=value, x=as.factor(Day), color=regime, group=Pop))+
  geom_line(alpha=0.5, aes(linetype=history))+
  geom_point(aes(shape=history))+
  scale_shape_manual(values=c(1, 19))+
  theme_bw(base_size = 14)+
  facet_wrap(~PC, scales='free')+
  ylab('')+
  xlab('age from egg (days)')



# test of convergence: does the extend (main effect of PC) or the trajectory (PC ~ age) of the more recently derived populations (history) differ from those of the long standing populations: 
library(lme4)
library(lmerTest)
d$Pop

aList <- list()
sub <- subset(pca, Day %in% c(21, 28, 35)) # subset to the sampling times (age) that are directly comparable across C and A regimes

for(i in 1:length(pcs)){
 aList[[i]] <- anova(lmer(sub[ ,pcs[i]] ~ regime * as.factor(Day) * history + (1|Pop), sub)) } # ANOVA with random effect of rep, to increase power to detect effects, given the non-independence of rep over age

pees <- sapply(aList, function(x) x[ ,6])
rownames(pees) <- rownames(aList[[1]])
colnames(pees) <- pcs
pees <- print(t(pees)) # [could correct for multicomp, ie, bonferroni] ....although a more liberal test is appropriate to ask if there is any difference (ie, to more strongly test the hypothesis of convergence)

# summary of test results: 
# this test was made using PCs that revealed highly significant variation in PC ~ age (PCs 3,4,8), and by regime (PCs 1,2,4), even after bonferroni correction (not shown in summary).
# :None of the first 12 PC differed in main effect (detectable by: regime:history), nor in PC ~ age (detectable by regime:as.factor(Day):history) between the more-recently derived populations and the longs-standing populations.
# there is no evidence that recent populations are any different than long-standing populations in their association with the metabolome nor the agr-trajectory of the metabolome.


pca$Day
daynames <- c(`9`="day 9", `21`="day 21", `28` = "day 28", `35`= "day 35", `70`= "day 70")

# in case people want to see a descriptive plot of PCs that includes the early and late un-matched ages:
PCA_AGE<-ggplot(pca, aes(y=PC2, x=PC1, color=Sel, fill=Sel))+
  stat_ellipse(geom='polygon', alpha=0.2, col=NA)+
  geom_point()+
  theme_bw()+
  facet_wrap(~Day, nrow=1, labeller = as_labeller(daynames))

PCA_AGE
## measure convergence with:

## as recommended by an MBE reviewer: "using PCoA, defining an axis of maximum differentiation between A-type and C-type populations and assessing how the old versus more recent populations score on that axis. "

# I think the reviewer meant a supervised method, when they said "defining an axis of maximum differentiation between A-type and C-type populations", they meant something like supervised PCA.  I don't know of a way to do supervised PCA with a categorical outcome, which is what i think the reviewer was suggesting.  The goal being to see how well the more-recently selected populations resemble the long-standing populations.  The idea to do it in PCoA "space" is a means to overcome the co-linearity problem with the metabolites, which is indeed a good idea.  I suggest two ways to ask this question: 

# PCoA/MDS is NOT a supervised method, and so is not explicitly deriving a metric to maximize the differences between A and C, upon which to evaluate the position of more-recent populations along that axis, which the I think is what the reviewer had in mind.  We do however have the PCA result and the ANOVA which failed to find an axis (PC) that was either directly modified by recent selection, or whose trajectory over age (PCS ~ age) was affected by recent selection ..which leads to the second method:

# a simple learning model to discriminate A and C in ling standing and ask how well it discriminates A and C among the more recently selected pops. This could include age and a covar, but may not need to in order to achieve the resolution/accuracy to ask the question.



##############################################################
### regime discrimination model
##############################################################
library(caret)
library(mlbench)
library(pROC)
library(plotROC)

sub <- subset(d, Day %in% c(21, 28, 35)) # subset to the sampling times (age) that are directly comparable across C and A regimes
old <- sub[sub$history=='longstanding', ]
new <- sub[sub$history=='recent', ]

m <- plsda(y=old$regime, x=old[ ,mzs], ncomp=3)
plot(m)

cm <- print(confusionMatrix(predict(m, new[ ,mzs]), new$regime))
cm$table
cm$byClass
cm$overall

sub$plsda <- predict(m, sub[ ,mzs])
table(sub$plsda, sub$regime, sub$history)


# have to make names compatible with train function in caret
sub$regime <- as.factor(make.names(sub$regime))
old <- sub[sub$history=='longstanding', ]
new <- sub[sub$history=='recent', ]


?trainControl

set.seed(1)
ctrl <- trainControl(method="cv", number=5, summaryFunction=twoClassSummary, classProbs=T, savePredictions = T)
m <- train(y=old$regime, x=old[ ,mzs], method="pls", trControl=ctrl)

m$results

head(m$pred)
bestTune <- as.numeric(m$bestTune) 

posplot<- ggplot(m$pred[m$pred$ncomp == bestTune, ], aes(m = C.type, d = factor(obs, levels = c("A.type", "C.type")), color= 'train')) + 
  geom_roc(labels=F, pointsize = 0)+
  theme_bw(base_size = 16)+
  xlab('false positive rate')+
  ylab('true positite rate')+
  coord_equal()+
  theme(legend.position='none')

posplot

confusionMatrix(predict(m, newdata=old[ ,mzs]), old$regime)
confusionMatrix(predict(m, newdata=new[ ,mzs]), new$regime)

probs <- predict(m, newdata=new[ ,mzs], type='prob') # if you want probabilities, prob A, prob C, you'll need to match the rows with the true obs in the data
head(probs)


probs <- predict(m, newdata=sub, type='prob')
head(probs)
probs$set <- as.factor(ifelse(rownames(probs) %in% rownames(new), 'testing', 'training'))

tmp <- data.frame(sub, probs)
tmp$set <- factor(tmp$set, levels=c('training', 'testing'))
levels(tmp$set) <- c('training\n(long-standing)', 'testing\n(recently derived)')

head(tmp)
levels(tmp$regime) <- c('A-type', 'C-type') # rename for plotting
probplot<- ggplot(tmp, aes(x=regime, y=`A.type`, color=set))+
  geom_point()+
  facet_wrap(~set)+
  theme_bw(base_size = 16)+
  ylab('probability of A-type')+
  xlab('regime')+
  geom_abline(intercept = 0.5, slope=0, linetype='dashed', color='grey40')+
  theme(legend.position = 'none')

probplot

confusionMatrix(predict(m, newdata=new[ ,mzs]), new$regime)

##############################################################
blank <- 0
ggarrange(posplot, blank, probplot, ncol = 3, nrow = 1, widths = c(1, 0.1, 0.75))
ggarrange(blank, PCA_AGE, blank, PCApanels, ncol = 1, nrow = 4, heights = c(0.1, 1.25, 0.15, 2.25))



