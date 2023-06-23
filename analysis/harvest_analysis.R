library('ggplot2')
library('gridExtra')
library('reshape2')
library('ggsignif')
library('car')
library('agricolae')
library('ggpubr')
library('ggfortify')
library('multcomp')
library('emmeans')
library('igraph')
library('ggrepel')
library('tidyverse')
library('FactoMineR')
library('factoextra')
library('corrplot')
library('ggbeeswarm')
library('ggbreak')

# Necessary functions
se <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))
central.median <- function(x) median(x, na.rm=TRUE)
central.mean <- function(x) mean(x, na.rm=TRUE)

# Aggregate, calculate mean, and se

# Generic function to automate plot
response.plot <- function(dataset, y, x, order=NULL, fill=NULL){
  if (is.null(order) != T & is.null(fill) != T){
    temp.dataset <-dataset[, c(y, x, order, fill)]
    colnames(temp.dataset) <- c("y", "x1", "o", 'f')
  } else {
    temp.dataset <- dataset[, c(y, x)]
    colnames(temp.dataset) <- c("y", "x1")
  }
  
  agg.mean <- aggregate(y ~ x1, data = temp.dataset, mean, 
                        na.rm = T)
  agg.se <- aggregate(y ~ x1, data = temp.dataset, se)
  colnames(agg.se)[2] <- "se"
  
  if (is.null(order) != T){
    agg.o <- aggregate(o ~ x1+f, data=temp.dataset, mean)
    colnames(agg.o)[2] <- 'f'
    colnames(agg.o)[3] <- 'o'
  }
  
  agg <- merge(agg.mean, agg.se, by = "x1")
  
  # Reorder dataset
  if (is.null(order) != T){
    agg <- merge(agg, agg.o, by='x1')
    agg <- transform(agg, x1=reorder(x1, -o))
  }
  
  
  # Base-plot
  if (is.null(order) != T){
    ggplot() +
      geom_bar(data = agg, aes(x=x1, y=y, fill=f), stat = "identity", position = "dodge", na.rm = T) + 
      geom_errorbar(data = agg, aes(x=x1, y=y, ymin = y - se, ymax = y + se), position = position_dodge(width = 0.9), width = 0.2) +
      xlab(x) + ylab(y) + theme(legend.title = element_blank(), axis.text.x = element_text(angle=90))
  } else {
    ggplot() +
      geom_bar(data = agg, aes(x=x1, y=y, fill=f), stat = "identity", position = "dodge", na.rm = T) + 
      geom_errorbar(data = agg, aes(x=x1, y=y, ymin = y - se, ymax = y + se), position = position_dodge(width = 0.9), width = 0.2) +
      xlab(x) + ylab(y) + theme(legend.title = element_blank(), axis.text.x = element_text(angle=90))
  }
  
}

# Generic function to automate plot
response.plot2 <- function(dataset, y, x, order=NULL, fill=NULL){
  if (is.null(order) != T & is.null(fill) != T){
    temp.dataset <-dataset[, c(y, x, order, fill)]
    colnames(temp.dataset) <- c("y", "x1", "o", 'f')
  } else {
    temp.dataset <- dataset[, c(y, x)]
    colnames(temp.dataset) <- c("y", "x1")
  }
  
  agg.mean <- aggregate(y ~ x1, data = temp.dataset, mean, 
                        na.rm = T)
  agg.se <- aggregate(y ~ x1, data = temp.dataset, se)
  colnames(agg.se)[2] <- "se"
  
  if (is.null(order) != T){
    agg.o <- aggregate(o ~ x1+f, data=temp.dataset, mean)
    colnames(agg.o)[2] <- 'f'
    colnames(agg.o)[3] <- 'o'
  }
  
  agg <- merge(agg.mean, agg.se, by = "x1")
  
  # Reorder dataset
  if (is.null(order) != T){
    agg <- merge(agg, agg.o, by='x1')
    agg <- transform(agg, x1=reorder(x1, -o))
  }
  
  
  # Base-plot
  if (is.null(order) != T){
    ggplot(data = agg, aes(x=x1, y=y, fill=f)) +
      geom_bar(stat = "identity", position = "dodge", na.rm = T) + 
      geom_errorbar(data = agg, aes(x=x1, y=y, ymin = y - se, ymax = y + se), position = position_dodge(width = 0.9), width = 0.2) +
      xlab(x) + ylab(y) + theme(legend.title = element_blank(), axis.text.x = element_text(angle=90))
  } else {
    ggplot(data = agg, aes(x=x1, y=y, fill=f)) +
      geom_bar(stat = "identity", position = "dodge", na.rm = T) + 
      geom_errorbar(data = agg, aes(x=x1, y=y, ymin = y - se, ymax = y + se), position = position_dodge(width = 0.9), width = 0.2) +
      xlab(x) + ylab(y) + theme(legend.title = element_blank(), axis.text.x = element_text(angle=90))
  }
  
}

# Get legend for multiple plots in single panel
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# read sampling harvest csv
data <- read.csv('../data/Sampling_harvest_v0.8.csv', header = TRUE)

# Subset only harvest relevant variables
harvest <- data[, c('plant_pos', 'plant_id', 'root_mass', 'shoot_mass', 'date', 'total_nodules', 'mean_weight')]

##############################
## Calculating new variables #
##############################

# Decouple treatment and replicate id
harvest$treatment <- as.character(gsub("_\\d+", '', harvest$plant_id))
harvest$replicate <- as.integer(gsub(".*_", '', harvest$plant_id))
length(unique(harvest$treatment))

# Seperating batches into variable
harvest$batch <- '1'
harvest[which(harvest$replicate>=6),]$batch <- '2'
#batch2 <- data[which(data$replicate>=6),]

# Calculate elapsed DPI
harvest$elapsed_dpi <- as.Date(harvest$date, format='%Y-%m-%d') - as.Date('2019-04-5', format='%Y-%m-%d')
harvest[which(harvest$batch==2),]$elapsed_dpi <- as.Date(harvest[which(harvest$batch==2),]$date, format='%Y-%m-%d') - as.Date('2019-04-11', format='%Y-%m-%d')

# converting weight to mg
harvest$root_mass <- harvest$root_mass*1000
harvest$shoot_mass <- harvest$shoot_mass*1000
harvest$mean_weight <- harvest$mean_weight*1000

# Relative Growth
DM_I_minus <- mean(harvest[which(harvest$treatment=='H2O'),]$shoot_mass + 
                     harvest[which(harvest$treatment=='H2O'),]$root_mass, na.rm = TRUE)

DM_I_plus <- harvest$shoot_mass + harvest$root_mass
harvest$RG <- DM_I_plus / DM_I_minus

# Shoot Relative Growth
harvest$ShootRG <- harvest$shoot_mass / mean(harvest[which(harvest$treatment=='H2O'),]$shoot_mass, na.rm=T)

# Symbiotic Efficiency
harvest$SEff <- DM_I_plus - DM_I_minus

# Investment to nodulation
harvest$investment <- (harvest$mean_weight*harvest$total_nodules) / (harvest$shoot_mass + harvest$root_mass)

head(harvest)


# Seperate mulitple plants from the data, ie. the plants coming from same pot
numbers_only <- function(x) !grepl("\\D", x)
harvest <- harvest[ which(numbers_only(harvest$plant_pos) == TRUE), ]

##################################################
## data breakdown by Fix+/-, T3SS, and abundance #
##################################################

# Loading trait file
traits.single <- as.data.frame(read.csv('../data/treatment_traits_single.csv', header=T))
traits.coinoc <- as.data.frame(read.csv('../data/treatment_traits_coinoculation.csv', header=T))
traits.all <- as.data.frame(read.csv('../data/treatment_traits_all.csv', header=T))

# merge harvest data and trait
harvest.single <- merge(harvest, traits.single, by.x='treatment', by.y='Treatment')
harvest.coinoc <- merge(harvest, traits.coinoc, by.x='treatment', by.y='Treatment')
harvest <- merge(harvest, traits.all, by.x='treatment', by.y='Treatment')

################################################
# Having two fixing strain is not always great #
################################################
three.single <- harvest.single[which(harvest.single$treatment %in% c(131, 184, 4)),]

psingle <- ggbarplot(three.single, x='treatment', y='RG', add = "mean_se", 
          position = position_dodge(0.8), fill='#069A8E') + xlab('')

three.coinoc <- harvest.coinoc[which(harvest.coinoc$treatment %in% c('131+184', '4+131', '4+184')),]

pcoinoc <- ggbarplot(three.coinoc, x='treatment', y='RG', add = "mean_se", 
          position = position_dodge(0.8), fill=c('#E8630A', '#069A8E', '#069A8E'))  + xlab('') +
          scale_y_continuous(limits = c(0, 10))


ggarrange(psingle, pcoinoc)
############################
# Working with nodule area #
############################
nod <- read.csv('../data/area_est_vc.tsv', header = TRUE, sep='\t')

nod$id <- as.character(gsub(".xls:*.\\d", '', nod$Nodules))
nod$number <- as.integer(gsub(".*xls:", '', nod$Nodules))
nod$treatment <- as.character(gsub("_\\d*", '', nod$id))
nod$AreaMM <- nod$Area*100

nod.clean <- (nod[nod$treatment %in% harvest$treatment, -1])

ggplot(nod.clean, aes(treatment, AreaMM)) + 
  geom_bar(stat='summary', fun.y=mean) + 
  theme(axis.text.x = element_text(angle=90))

## Pairwise t-test for independent group
#pwc <- nod.clean %>%
#  pairwise_t_test(AreaMM ~ treatment, pool.sd=F,
#                  p.adjust.method = 'bonferroni')#
#
#pwc <- pwc %>% add_xy_position(x='treatment', step.increase = 1)

#ggboxplot(nod.clean, x='treatment', y='AreaMM')



# Dry nodule biomass based on Wendlandt's formula 2019
nod$dryMass <- (nod$AreaMM - 0.9097853)/5.525844
nodule_biomass <- aggregate(nod[, c(9, 10)], by=list(nod$id), FUN=sum)

harvest <- merge(harvest, nodule_biomass, by.y='Group.1', by.x='plant_id', all.x=T, all.y=F)
harvest.coinoc <- merge(harvest.coinoc, nodule_biomass, by.y='Group.1', by.x='plant_id', all.x=T, all.y=F)
harvest.single <- merge(harvest.single, nodule_biomass, by.y='Group.1', by.x='plant_id', all.x=T, all.y=F)

dim(harvest)
#write.csv(harvest, '../data/harvest_master.csv', row.names = F, quote = F)
#######################################################
# Pearson correlation of fresh nodule biomass and HGR #
#######################################################
cor.test(harvest.single$RG, harvest.single$mean_weight, method='pearson')
cor.test(harvest.single$RG, harvest.single$total_nodules, method='pearson')


cor.test(harvest.coinoc$RG, harvest.coinoc$total_nodules, method='pearson')
cor.test(harvest.coinoc$RG, harvest.coinoc$mean_weight, method='pearson')
cor.test(harvest.coinoc$RG, harvest.coinoc$investment, method='pearson')



#################
# Order dataset #
#################
harvest.single.ordered <- harvest.single[order(harvest.single$Fix), ]
harvest.single.ordered$treatment <- factor(harvest.single.ordered$treatment, levels=unique(harvest.single.ordered$treatment))

harvest.coinoc.ordered <- harvest.coinoc[order(harvest.coinoc$Fix), ]
harvest.coinoc.ordered$treatment <- factor(harvest.coinoc.ordered$treatment, levels=unique(harvest.coinoc.ordered$treatment))



###########################
# Analysis for manuscript #
###########################




##############################
# Are there effect of batch? #
# Testing Mixed Models #######
########################
library(car)
library(MASS)
library(lme4)

######################
# Single inoculation #
######################

# RG
qqp(log10(harvest.single.ordered$RG), "norm")
mod1 <- lmer(log(RG)~treatment + (1|batch),  data=harvest.single.ordered, REML=F)
mod2 <- lm(log(RG)~treatment,  data=harvest.single.ordered)
anova(mod1, mod2)




# Total nodules
qqp(sqrt(harvest.single.ordered$total_nodules), "norm")
mod1 <- lmer(sqrt(total_nodules)~treatment + (1|batch),  data=harvest.single.ordered, REML=F)
mod2 <- lm(sqrt(total_nodules)~treatment,  data=harvest.single.ordered)
anova(mod1, mod2)

# Nodule biomass
qqp(log(harvest.single.ordered$mean_weight+0.1), "norm")
mod1 <- lmer(log(harvest.single.ordered$mean_weight+0.1)~treatment + (1|batch),  data=harvest.single.ordered, REML=F)
mod2 <- lm(log(harvest.single.ordered$mean_weight+0.1)~treatment,  data=harvest.single.ordered)
anova(mod1, mod2)

# Investment
qqp(log10(harvest.single.ordered$investment+0.2), "norm")
mod1 <- lmer(log10(harvest.single.ordered$investment+0.2)~treatment + (1|batch),  data=harvest.single.ordered, REML=F)
mod2 <- lm(log10(harvest.single.ordered$investment+0.2)~treatment,  data=harvest.single.ordered)
anova(mod1, mod2)

#################
# Coinoculation #
#################

# RG
qqp(log10(harvest.coinoc.ordered$RG), "norm")
model1 <- lmer(log10(RG)~treatment + (1|batch), data=harvest.coinoc.ordered)
model2 <- lm(log10(RG)~treatment, data=harvest.coinoc.ordered)
anova(model1, model2)

# Total nodules
qqp(sqrt(harvest.coinoc.ordered$total_nodules+2), "norm")
mod1 <- lmer(sqrt(total_nodules+2)~treatment + (1|batch),  data=harvest.coinoc.ordered, REML=F)
mod2 <- lm(sqrt(total_nodules+2)~treatment,  data=harvest.coinoc.ordered)
anova(mod1, mod2)
summary(mod1)

# Nodule biomass
qqp(log(harvest.coinoc.ordered$mean_weight+0.1), "norm")
mod1 <- lmer(log(harvest.coinoc.ordered$mean_weight+0.1)~treatment + (1|batch),  data=harvest.coinoc.ordered, REML=F)
mod2 <- lm(log(harvest.coinoc.ordered$mean_weight+0.1)~treatment,  data=harvest.coinoc.ordered)
anova(mod1, mod2)

# Investment
qqp(log10(harvest.coinoc.ordered$investment+0.2), "norm")
mod1 <- lmer(log10(harvest.coinoc.ordered$investment+0.2)~treatment + (1|batch),  data=harvest.coinoc.ordered, REML=F)
mod2 <- lm(log10(harvest.coinoc.ordered$investment+0.2)~treatment,  data=harvest.coinoc.ordered)
anova(mod1, mod2)

#############################################################
#############################################################

# Single Innoculation
# By Fix Trait

## Unpaired 2-sample t-test for log10(RG), total_nodules, mean_weight

# For RG
t.test(log10(harvest.single.ordered[which(harvest.single.ordered$Fix == '-'), ]$ShootRG), 
       log10(harvest.single.ordered[which(harvest.single.ordered$Fix == '+'), ]$ShootRG))
# For total_nodules
t.test(harvest.single.ordered[which(harvest.single.ordered$Fix == '-'), ]$total_nodules, 
       harvest.single.ordered[which(harvest.single.ordered$Fix == '+'), ]$total_nodules)

# For total_nodules
t.test(harvest.single.ordered[which(harvest.single.ordered$Fix == '-'), ]$mean_weight, 
       harvest.single.ordered[which(harvest.single.ordered$Fix == '+'), ]$mean_weight)


## Visualization by Fix
p1 <- ggplot(harvest.single.ordered[which(harvest.single.ordered$Fix != 'Control'), ], aes(x=Single_effectivity, y=log10(ShootRG), fill=Single_effectivity)) +
  geom_boxplot(notch=T) + xlab('') + 
  geom_signif(comparison=list(c('-', '+')), map_signif_level=T)

p2 <- ggplot(harvest.single.ordered[which(harvest.single.ordered$Fix != 'Control'), ], aes(x=Single_effectivity, y=total_nodules, fill=Single_effectivity)) +
  geom_boxplot(notch=T) + ylab('Total Nodules') + xlab('') +
  geom_signif(comparison=list(c('-', '+')), map_signif_level=T)

p3 <- ggplot(harvest.single.ordered[which(harvest.single.ordered$Fix != 'Control'), ], aes(x=Single_effectivity, y=mean_weight, fill=Single_effectivity)) +
  geom_boxplot(notch=T) + xlab ('') + ylab('Mean fresh nodule biomass') +
  geom_signif(comparison=list(c('-', '+')), map_signif_level=T) 

ggarrange(p1, p2, p3, nrow=1, common.legend = T, legend='bottom')



## Plants with nodules, reject all plants w/o nodules
single.nod <- harvest.single.ordered[which(harvest.single.ordered$total_nodules > 0),]

## Anova and lm for RG
hist(log10(single.nod$ShootRG+1))
mod.rg <- lm(log10(ShootRG+0.5) ~ elapsed_dpi + treatment, data=single.nod)
summary(mod.rg)

autoplot(mod.rg)
shapiro.test(residuals(mod.rg))

res <- aov(mod.rg)
summary(res)

Anova(mod.rg, type='III')

(HSD.test(aov(mod.rg), 'treatment', group=T))
emm <- emmeans(mod.rg, pairwise~treatment:elapsed_dpi, type='response')
emm$emmeans


## anova and lm for nodules
hist(sqrt(single.nod$total_nodules))
mod.nod <- lm(sqrt(total_nodules) ~ elapsed_dpi + treatment, data=single.nod)
summary(mod.nod)
autoplot(mod.nod)

Anova(mod.nod, type='III')
shapiro.test(residuals(mod.nod))

res <- aov(mod.nod)
summary(res)


(HSD.test(aov(mod.nod), 'treatment', group=T))




emm <- emmeans(mod.nod, pairwise~treatment:elapsed_dpi, type='response')
emm$emmeans

#p.nod <-symbio.plot(harvest.single.ordered[which(harvest.single.ordered$Fix != 'Control'), ], 
#                    'total_nodules', 'treatment', 'Fix') +
#  geom_text(aes(label=c('a', 'ab', 'ab', 'abc', 'bc', 'c', 'cd', 'bc')), vjust=-5) +
#  ggtitle('Total Nodule')

## anoval and lm for mean fresh nodule biomass
hist(log(single.nod$mean_weight+0.1))
mod.nod.biomass <- lm(log(mean_weight+0.1) ~ elapsed_dpi + treatment, data=single.nod)
summary(mod.nod.biomass)
autoplot(mod.nod.biomass)
shapiro.test(residuals(mod.nod.biomass))

Anova(mod.nod.biomass, type='III')

res <- aov(mod.nod.biomass)
summary(res)
(HSD.test(aov(mod.nod.biomass), 'treatment', group=T))


# Comparing nodule dry mass and area by t-test
t.test(single.nod[which(single.nod$Single_effectivity == 'Effective'), ]$dryMass, 
       single.nod[which(single.nod$Single_effectivity == 'Ineffective'), ]$dryMass)

t.test(single.nod[which(single.nod$Single_effectivity == 'Effective'), ]$AreaMM, 
       single.nod[which(single.nod$Single_effectivity == 'Ineffective'), ]$AreaMM)

t.test(single.nod[which(single.nod$Single_effectivity == 'Effective'), ]$total_nodules, 
       single.nod[which(single.nod$Single_effectivity == 'Ineffective'), ]$total_nodules)

emm <- emmeans(mod.nod.biomass, pairwise~treatment:elapsed_dpi, type='response')
emm$emmeans

#p.biomass <- symbio.plot(harvest.single.ordered[which(harvest.single.ordered$Fix != 'Control'), ],
#                         'mean_weight', 'treatment', 'Fix') +
#  geom_text(aes(label=c('a', 'ab', 'ab', 'ab', 'cd', 'ab', 'cd', 'bc')), vjust=-5) +
#  ggtitle('Mean fresh nodule biomass')

#grid.arrange(p1, p2, p3, ncol=1)

## anova and lm for investment
hist(log10(single.nod$investment + 0.2))
mod.invest <- lm(log10(investment + 0.2) ~ elapsed_dpi + treatment, data=single.nod)
summary(mod.invest)
autoplot(mod.invest)
shapiro.test(residuals(mod.invest))

res <- aov(mod.invest)
summary(res)

Anova(mod.invest, type='III')

(HSD.test(aov(mod.invest), 'treatment', group=T))


###########################
# Single Inoculation Data #
###########################
# Main result figure
single.nod[which(single.nod$Fix == '+'),]$Fix <- 'Fix+'
single.nod[which(single.nod$Fix == '-'),]$Fix <- 'Fix-'
single.nod$Fix <- factor(single.nod$Fix, levels = c('Fix+', 'Fix-'))

q1 <- response.plot2(single.nod, 'ShootRG', 'treatment', 'ShootRG', 'Fix')
#q1$data
#q1$data[which(q1$data$f == '-'),]$f <- 'Fix-'
#q1$data[which(q1$data$f == '+'),]$f <- 'Fix+'


rg_sig <- data.frame(treatment=c(131,184,4,156,2,186,200,187),
                     sig= c('*', '*', '*', 'NS', 'NS', 'NS', 'NS', '·'))



q1 <- q1 + geom_hline(yintercept = 1, color='blue1') + ggtitle('A. Host Growth Benefit')+  ylab('Shoot RG') + xlab('Treatment') + 
  geom_text(aes(y=y+se, label=sig), col='#8b0000', vjust= -0.05) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

q2 <- response.plot2(single.nod, 'total_nodules', 'treatment', 'ShootRG', 'Fix') +
  ggtitle('B. Nodules Per Plant') + ylab('Average Nodules') + xlab('Treatment') + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

q3 <- response.plot2(single.nod, 'mean_weight', 'treatment', 'ShootRG', 'Fix') + 
  ggtitle('C. Mean Fresh Nodule Biomass') + ylab('Mean Weight (mg)') + xlab('Treatment') + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

q4 <- response.plot2(single.nod, 'investment', 'treatment', 'ShootRG', 'Fix') + 
  ggtitle('D. Investment to Nodulation') + ylab('Nodule Biomass / Shoot Biomass') + xlab('Treatment') + theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

ggarrange(q1, q2, q3, q4, nrow=2, ncol=2, common.legend=T, legend='top')


# Main result figures with violin plot/jittered dataset
q0 <- response.plot(single.nod, 'ShootRG', 'treatment', 'ShootRG', 'Fix')

rg_sig <- data.frame(treatment=c('131','184','4','156','2','186','200','187'),
                     rg =c(10.5, 9, 6, 4, 4, 3.5, 3, 1),
                     sig= c('*', '*', '•', '*', 'NS', 'NS', 'NS', 'NS'))


q1 <- q0 + 
  geom_beeswarm(data=single.nod, aes(x=treatment, y=RG),size=1, alpha=0.3) +
  geom_hline(yintercept = 1, color='blue1') + 
  ggtitle('A')+  ylab('Shoot RG') + xlab('Treatment') + 
  geom_text(data=rg_sig, aes(x=treatment, y=rg, label=sig), col='#8b0000', vjust= -0.05) +
  theme_bw() + theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

q2 <- response.plot(single.nod, 'total_nodules', 'treatment', 'ShootRG', 'Fix') + 
  geom_beeswarm(data=single.nod, aes(x=treatment, y=total_nodules),size=1, alpha=0.3) +
  ggtitle('B') + ylab('Average Nodules') + xlab('Treatment') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

q3 <- response.plot(single.nod, 'mean_weight', 'treatment', 'ShootRG', 'Fix') + 
  geom_beeswarm(data=single.nod, aes(x=treatment, y=mean_weight),size=1, alpha=0.3) +
  ggtitle('C') + ylab('Mean Weight (mg)') + xlab('Treatment') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

q4 <- response.plot(single.nod, 'investment', 'treatment', 'ShootRG', 'Fix') + 
  geom_beeswarm(data=single.nod, aes(x=treatment, y=investment),size=1, alpha=0.3) +
  ggtitle('D') + ylab('Nodule Biomass / Shoot Biomass') + xlab('Treatment') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c('#2E8B57', '#F8766D'))

ggarrange(q1, q2, q3, q4, nrow=2, ncol=2, common.legend=T, legend='top')

# Estimated marginal means
emm <- emmeans(mod.invest, pairwise~treatment:elapsed_dpi, type='response')
emm$emmeans



## Testing N-fixing effectiveness
t.test(single.nod[which(single.nod$treatment == '131'), ]$ShootRG, mu=1, alternative = 'greater')


#################
# Coinoculation #
#################

coinoc.nod <- harvest.coinoc.ordered[which(harvest.coinoc.ordered$total_nodules > 0),]


coinoc.nod[which(coinoc.nod$Fix == '+/+'),]$Fix <- 'Fix+/Fix+'
coinoc.nod[which(coinoc.nod$Fix == '+/-'),]$Fix <- 'Fix+/Fix-'
coinoc.nod[which(coinoc.nod$Fix == '-/-'),]$Fix <- 'Fix-/Fix-'

coinoc.nod$Fix <- factor(coinoc.nod$Fix, levels = c('Fix+/Fix+', 'Fix+/Fix-', 'Fix-/Fix-'))

# Main result figure
p1 <- response.plot(coinoc.nod, 'ShootRG', 'treatment', 'ShootRG', 'Fix')
p1$data$sig <- c('*','NS','*','.','NS','NS','.','*','NS','.','NS','NS','*','NS', '.','NS','NS','NS','.','NS','NS', 'NS', 'NS','*','NS','*','*','*')
p1 <- p1 + 
  geom_hline(yintercept = 1, color='blue1') + ggtitle('a. Host Growth Benefit') + xlab('Treatment') + ylab('Shoot RG') +
  geom_text(aes(y=y+se, label=sig), col='#8b0000', vjust= -0.25) + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=c('#228B22', '#364ADD', '#F8766D'))
p2 <- response.plot(coinoc.nod, 'total_nodules', 'treatment', 'ShootRG', 'Fix') +
  ggtitle('b. Nodules Per Plant')  + xlab('Treatment') +
  theme(legend.position='none') + ylab('Average Nodules') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=c('#228B22', '#364ADD', '#F8766D'))
p3 <- response.plot(coinoc.nod, 'mean_weight', 'treatment', 'ShootRG', 'Fix') +
  ggtitle('c. Mean Fresh Nodule Biomass')  + xlab('Treatment') +
  theme(legend.position='none') + ylab('Mean Weight (mg)') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=c('#228B22', '#364ADD', '#F8766D'))
p4 <- response.plot(coinoc.nod, 'investment', 'treatment', 'ShootRG', 'Fix') +
  ggtitle('d. Investment to Nodulation')  +
  theme(legend.position='none') + ylab('Nodule Biomass / Shoot Biomass') + xlab('Treatment') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_manual(values=c('#228B22', '#0096FF', '#F8766D'))


ggarrange(p1, p2, p3, p4, nrow=2, ncol=2, common.legend = T, legend='top')

p5 <- response.plot(coinoc.nod, 'AreaMM', 'treatment', 'ShootRG', 'Fix')
p6 <- response.plot(coinoc.nod, 'dryMass', 'treatment', 'ShootRG', 'Fix')

# Main result figure with datapoints
p0 <- response.plot(coinoc.nod, 'ShootRG', 'treatment', 'ShootRG', 'Fix')
coinoc_rg <- data.frame(signf <- c('•', '*', '•', '*', '*', '*', '*', '*', '•', '•', '*', '*', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS', 'NS' ),
                        treatment <- c('184+200', '4+200', '187+184', '186+4', '131+187','4+131',
                                       '4+184', '4+187', '131+200', '186+184', '131+156', '186+131', '2+184', '2+131',
                                       '156+184', '4+156', '2+4', '156+200', '131+184', '186+187',
                                       '187+156', '186+200', '2+186', '186+156', '2+156', '2+200',
                                       '187+200', '2+187'),
                        y <- c(7.5, 6, 6.5, 6.5, 5.5, 5.5, 5, 5, 5, 5, 
                               4.5, 4.4, 4, 4, 3.5, 3, 2.8, 2.5, 2.3, 2.2, 
                               2, 2, 1, 1, 1, 1, 1, 1))

p1 <- p0 + 
  geom_beeswarm(data=coinoc.nod, aes(x=treatment, y=ShootRG), size=1, alpha=0.3) +
  geom_hline(yintercept = 1, color='blue1') + ggtitle('A') + xlab('Treatment') + ylab('Shoot RG') +
  geom_text(data=coinoc_rg, aes(y=y, x=treatment, label=signf), col='#8b0000', vjust= -0.25) +
  theme_bw() + theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position='none') + 
  scale_fill_manual(values=c('#228B22', '#0096FF', '#F8766D'))
p2 <- response.plot(coinoc.nod, 'total_nodules', 'treatment', 'ShootRG', 'Fix') +
  geom_beeswarm(data=coinoc.nod, aes(x=treatment, y=total_nodules), size=1, alpha=0.3) +
  ggtitle('B')  + xlab('Treatment') +
  ylab('Average Nodules') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position='none') + 
  scale_fill_manual(values=c('#228B22', '#0096FF', '#F8766D'))

p3 <- response.plot(coinoc.nod, 'mean_weight', 'treatment', 'ShootRG', 'Fix') +
  geom_beeswarm(data=coinoc.nod, aes(x=treatment, y=mean_weight), size=1, alpha=0.3) +
  ggtitle('C')  + xlab('Treatment') +
  scale_y_break(breaks = c(8, 18)) +
  ylab('Mean Weight (mg)') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position='none') + 
  scale_fill_manual(values=c('#228B22', '#0096FF', '#F8766D'))


p4 <- response.plot(coinoc.nod, 'investment', 'treatment', 'ShootRG', 'Fix') +
  geom_beeswarm(data=coinoc.nod, aes(x=treatment, y=investment), size=1, alpha=0.3) +
  ggtitle('D')  +
  ylab('Nodule Biomass / Shoot Biomass') + xlab('Treatment') + theme_bw() +
  theme(plot.title = element_text(face = "bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(legend.position='none') + 
  scale_fill_manual(values=c('#228B22', '#0096FF', '#F8766D'))

ggarrange(print(p1), print(p2), print(p3), print(p4), nrow=2, ncol=2, common.legend = T, legend='top')

# Heatmap 
# Main result figure
p1 <- response.plot2(coinoc.nod, 'ShootRG', 'treatment', 'ShootRG', 'Fix')
p2 <- response.plot2(coinoc.nod, 'total_nodules', 'treatment', 'ShootRG', 'Fix')
p3 <- response.plot2(coinoc.nod, 'mean_weight', 'treatment', 'ShootRG', 'Fix')
p4 <- response.plot2(coinoc.nod, 'investment', 'treatment', 'ShootRG', 'Fix')
p5 <- response.plot2(coinoc.nod, 'AreaMM', 'treatment', 'ShootRG', 'Fix')
p6 <- response.plot2(coinoc.nod, 'dryMass', 'treatment', 'ShootRG', 'Fix')


colnames(p1$data) <- c('treatment', 'ShootRG_mean', 'ShootRG_se', 'group', 'ShootRG')
colnames(p2$data) <- c('treatment', 'Nodules_mean', 'Nodules_se', 'group', 'Nodules')
colnames(p3$data) <- c('treatment', 'NoduleWeight_mean', 'NoduleWeight_se', 'group', 'NoduleWeight')
colnames(p4$data) <- c('treatment', 'Investment_mean', 'Investment_se', 'group', 'Investment')
colnames(p5$data) <-  c('treatment', 'NoduleArea_mean', 'NoduleArea_se', 'group', 'NoduleArea')
colnames(p6$data) <-  c('treatment', 'DryMass_mean', 'DryMass_se', 'group', 'DryMass')

group <- merge(p1$data, p2$data, by=c('treatment', 'group'))
group <- merge(group, p3$data, by=c('treatment', 'group'))
group <- merge(group, p4$data, by=c('treatment', 'group'))
group <- merge(group, p5$data, by=c('treatment', 'group'))
group <- merge(group, p6$data, by=c('treatment', 'group'))


group.heatmap <- group[,c('treatment', 'group', 'ShootRG_mean', 'Nodules_mean', 'NoduleWeight_mean', 'Investment_mean', 'NoduleArea_mean', 'DryMass_mean')]
colnames(group.heatmap) <- c('Treatment', 'Group', 'ShootRG', 'Total Nodules', 'Fresh Nodule Weight', 'Investment', 'Nodule Area', 'Dry Nodule Weight')

group.heatmap[,c('ShootRG', 'Total Nodules', 'Fresh Nodule Weight', 'Investment', 'Nodule Area', 'Dry Nodule Weight')] <- 
  scale(group.heatmap[,c('ShootRG', 'Total Nodules', 'Fresh Nodule Weight', 'Investment', 'Nodule Area', 'Dry Nodule Weight')])

group.heatmap <- melt(group.heatmap, id.vars = c('Treatment', 'Group'))
variable_order <- c('ShootRG', 'Total Nodules', 'Fresh Nodule Weight', 'Dry Nodule Weight', 'Nodule Area', 'Investment')



ggplot(group.heatmap, aes(y=factor(variable, level=variable_order), x=Treatment, fill=value)) + 
  geom_tile(color='black') + 
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  ggtitle('E') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
        strip.text.x = element_text(size = 14)) +
  theme(plot.title = element_text(face = "bold")) +
  facet_grid(~Group, scales="free_x", space='free_x') +
  labs(fill='Z-Score') + ylab('') 

# Grouping the dataset by Single-effectivity combination
p1 <- response.plot2(coinoc.nod, 'ShootRG', 'treatment', 'ShootRG', 'Fix')
p2 <- response.plot2(coinoc.nod, 'total_nodules', 'treatment', 'ShootRG', 'Fix')
p3 <- response.plot2(coinoc.nod, 'mean_weight', 'treatment', 'ShootRG', 'Fix')
p4 <- response.plot2(coinoc.nod, 'investment', 'treatment', 'ShootRG', 'Fix')


p1.summary <- as.data.frame(as.matrix(aggregate(y~f, data= p1$data, FUN=function(x) c(mean=mean(x), se= se(x)))))
p2.summary <- as.data.frame(as.matrix(aggregate(y~f, data= p2$data, FUN=function(x) c(mean=mean(x), se= se(x)))))
p3.summary <- as.data.frame(as.matrix(aggregate(y~f, data= p3$data, FUN=function(x) c(mean=mean(x), se= se(x)))))
p4.summary <- as.data.frame(as.matrix(aggregate(y~f, data= p4$data, FUN=function(x) c(mean=mean(x), se= se(x)))))


p1s <- ggplot(p1.summary, aes(x=f, y=as.numeric(y.mean), fill=f)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=as.numeric(y.mean)-as.numeric(y.se), ymax=as.numeric(y.mean)+as.numeric(y.se)), width=0.2, colour='blue') +
  ggtitle('A. RG') +  ylab('') + theme(legend.position='none')
  

p2s <- ggplot(p2.summary, aes(x=f, y=as.numeric(y.mean), fill=f)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=as.numeric(y.mean)-as.numeric(y.se), ymax=as.numeric(y.mean)+as.numeric(y.se)), width=0.2, colour='blue') +
  ggtitle('B. Total Nodules') +  ylab('') + theme(legend.position='none')
  

p3s <- ggplot(p3.summary, aes(x=f, y=as.numeric(y.mean), fill=f)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=as.numeric(y.mean)-as.numeric(y.se), ymax=as.numeric(y.mean)+as.numeric(y.se)), width=0.2, colour='blue') +
  ggtitle('C. Mean Fresh Nodule Biomass') +  ylab('') + theme(legend.position='none')
  

p4s <- ggplot(p4.summary, aes(x=f, y=as.numeric(y.mean), fill=f)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=as.numeric(y.mean)-as.numeric(y.se), ymax=as.numeric(y.mean)+as.numeric(y.se)), width=0.2, colour='blue') +
  ggtitle('D. Investment to Nodulation') +  ylab('') + theme(legend.position='none')
  

ggarrange(p1s, p2s, p3s, p4s, nrow=2, ncol=2, common.legend = T, legend='top')

# t-test
# Main result figure
t_test_comparison.q1 <- compare_means(y~f, data=p1$data)
t_test_comparison.q1 <- t_test_comparison.q1 %>% mutate(y.position=c(4.5, 5, 5.5))

q1 <- ggbarplot(p1$data, x='f', y='y', fill='f', add='mean_se') +
  stat_pvalue_manual(t_test_comparison.q1, label="p.adj") +
  #stat_compare_means(method='t.test', comparisons= list(c('Fix+/Fix+', 'Fix+/Fix-'),
  #                                                      c('Fix+/Fix+', 'Fix-/Fix-'),
  #                                                      c('Fix+/Fix-', 'Fix-/Fix-'))) +
  theme(legend.position='none', legend.title=element_blank(), plot.title=element_text(face="bold")) +
  xlab('Strain Fix Trait Combinations') + ylab('Shoot RG') + ggtitle('A') +
  scale_fill_manual(values=c('#00BA38', '#619CFF', '#F8766D'))

t_test_comparison.q2 <- compare_means(y~f, method='t.test', data=p2$data)
t_test_comparison.q2 <- t_test_comparison.q2 %>% mutate(y.position=c(25, 28, 31))

q2 <- ggbarplot(p2$data, x='f', y='y', fill='f', add='mean_se') + 
  stat_pvalue_manual(t_test_comparison.q2, label="p.adj") +
  #stat_compare_means(method='t.test', comparisons= list(c('Fix+/Fix+', 'Fix+/Fix-'),
  #                                                      c('Fix+/Fix+', 'Fix-/Fix-'),
  #                                                      c('Fix+/Fix-', 'Fix-/Fix-'))) +
  theme(legend.position='none', legend.title=element_blank(), plot.title=element_text(face="bold")) +
  xlab('Strain Fix Trait Combinations') + ylab('Total Nodules') + ggtitle('B') +
  scale_fill_manual(values=c('#00BA38', '#619CFF', '#F8766D'))

t_test_comparison.q3 <- compare_means(y~f, data=p3$data)
t_test_comparison.q3 <- t_test_comparison.q3 %>% mutate(y.position=c(2.3, 2.6, 3.1))

q3 <- ggbarplot(p3$data, x='f', y='y', fill='f', add='mean_se') + 
  stat_pvalue_manual(t_test_comparison.q3, label="p.adj") +
  #stat_compare_means(method='t.test', comparisons= list(c('Fix+/Fix+', 'Fix+/Fix-'),
  #                                                      c('Fix+/Fix+', 'Fix-/Fix-'),
  #                                                      c('Fix+/Fix-', 'Fix-/Fix-'))) +
  theme(legend.position='none', legend.title=element_blank(), plot.title=element_text(face="bold")) +
  xlab('Strain Fix Trait Combinations') + ylab('. Mean Fresh Nodule Biomass (mg)') + ggtitle('C') +
  scale_fill_manual(values=c('#00BA38', '#619CFF', '#F8766D'))

t_test_comparison.q4 <- compare_means(y~f, data=p4$data)
t_test_comparison.q4 <- t_test_comparison.q4 %>% mutate(y.position=c(1.4, 1.7, 2))

q4 <- ggbarplot(p4$data, x='f', y='y', fill='f', add='mean_se') +  
  stat_pvalue_manual(t_test_comparison.q4, label="p.adj") +
  #stat_compare_means(method='t.test', comparisons= list(c('Fix+/Fix+', 'Fix+/Fix-'),
  #                                                      c('Fix+/Fix+', 'Fix-/Fix-'),
  #                                                      c('Fix+/Fix-', 'Fix-/Fix-'))) +
  theme(legend.position='none', legend.title=element_blank(), plot.title=element_text(face="bold")) +
  xlab('Strain Fix Trait Combinations') + ylab('Investment to Nodulation') + ggtitle('D') +
  scale_fill_manual(values=c('#00BA38', '#619CFF', '#F8766D'))

  
  
ggarrange(q1, q2, q3,q4, nrow=2, ncol=2, common.legend = T, legend='top')



## Modeling RG
hist(log(coinoc.nod$ShootRG+0.5))
mod.rg.2 <- lm(log(ShootRG+0.5) ~ elapsed_dpi+treatment, data=coinoc.nod)
mod.rg.2 <- lm(log(ShootRG+0.5) ~ treatment + elapsed_dpi, data=coinoc.nod)
summary(mod.rg.2)

autoplot(mod.rg.2)
shapiro.test(residuals(mod.rg.2))

Anova(mod.rg.2, type='III')

## Modeling total_nodules
hist(sqrt(coinoc.nod$total_nodules))
mod.nod.2 <- lm(sqrt(total_nodules) ~ elapsed_dpi + treatment, data=coinoc.nod)
# Exclude outlier
#mod.nod.2 <- lm(sqrt(total_nodules) ~ treatment*elapsed_dpi, data=coinoc.nod[rownames(coinoc.nod) != '58',])
summary(mod.nod.2)

autoplot(mod.nod.2)
shapiro.test(residuals(mod.nod.2))

Anova(mod.nod.2, type='III')

## Modeling mean_weight
hist(log(coinoc.nod$mean_weight+0.1))
mod.weight.2 <- lm(log(mean_weight+0.1) ~ elapsed_dpi + treatment, data=coinoc.nod)
# Exclude outlier
#mod.weight.2 <- lm(log(mean_weight+0.1) ~ treatment * elapsed_dpi, data=coinoc.nod[!rownames(coinoc.nod) %in% c('135', '240', '14'),])
summary(mod.weight.2)

autoplot(mod.weight.2)
shapiro.test(residuals(mod.weight.2))

Anova(mod.weight.2, type='III')

contrasts.mean_weight <- emmeans(mod.weight.2, pairwise~treatment:elapsed_dpi, type='response')
contrasts.mean_weight <- as.data.frame(contrasts.mean_weight$contrasts)
contrasts.mean_weight[which(contrasts.mean_weight$p.value <= 0.05), ]
# Modeling investment
hist(log(coinoc.nod$investment+0.2))
mod.invest.2 <- lm(log(investment+0.2) ~ elapsed_dpi + treatment, data=coinoc.nod)
summary(mod.invest.2)

autoplot(mod.invest.2)
shapiro.test(residuals(mod.invest.2))

Anova(mod.invest.2, type='III')
summary(aov(mod.invest.2))
# 
p1 <- ggplot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
       aes(x=Fix, y=log10(RG), fill=Fix)) + geom_boxplot(notch=T) + 
  stat_compare_means(comparison=list(c('-/-', '+/-'), c('+/-', '+/+'), c('-/-', '+/+')), label.y=c(1.3, 1.5, 1.7))

p2 <- ggplot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
             aes(x=Fix, y=total_nodules, fill=Fix)) + geom_boxplot(notch=T) + 
  stat_compare_means(comparison=list(c('-/-', '+/-'), c('+/-', '+/+'), c('-/-', '+/+')), label.y=c(40, 42, 45))

p3 <- ggplot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
             aes(x=Fix, y=mean_weight, fill=Fix)) + geom_boxplot(notch=T) + 
  stat_compare_means(comparison=list(c('-/-', '+/-'), c('+/-', '+/+'), c('-/-', '+/+')), label.y=c(4, 6, 8))

# Contrast on traits as participant isolate's trait in pariwise combination
head(harvest.coinoc)

mod.rg.pair <- lm(log10(ShootRG+0.5) ~ elapsed_dpi + Single_effectivity, data=coinoc.nod)
summary(mod.rg.pair)
autoplot(mod.rg.pair)
Anova(mod.rg.pair, type='III')

emm <- emmeans(mod.rg.pair, pairwise~Single_effectivity:elapsed_dpi, type='response')
emm$emmeans

# Contrast on traits as a whole in pariwise combination
head(harvest.coinoc)

mod.rg.whole <- lm(log10(ShootRG+0.5) ~ elapsed_dpi + Co_effectivity, data=coinoc.nod)
summary(mod.rg.whole)
autoplot(mod.rg.whole)
Anova(mod.rg.whole, type='III')

emm <- emmeans(mod.rg.whole, pairwise~Co_effectivity:elapsed_dpi, type='response')
emm$emmeans

#################
# All together #
###############

harvest.nod <- harvest[which(harvest$total_nodules > 0), ]

harvest.nod$harvest_batch  <- 0
harvest.nod[which(harvest.nod$elapsed_dpi > 40), ]$harvest_batch <- 1
harvest.nod[which(harvest.nod$elapsed_dpi > 55), ]$harvest_batch <- 2
#barplot(harvest.nod$harvest_batch)


# RG modeling
hist(log10(harvest.nod$ShootRG+0.5))
mod.rg.3 <- lm(log10(ShootRG+0.5) ~ elapsed_dpi + treatment , data=harvest.nod)
summary(mod.rg.3)
shapiro.test(residuals(mod.rg.3))
autoplot(mod.rg.3)

summary(aov(mod.rg.3))
Anova(mod.rg.3, type='III')


#Anova(mod.rg.3, type='III', singular.ok=T)

#pwc <- emmeans_test(data=harvest.coinoc,
#  ShootRG ~ treatment, covariate = elapsed_dpi,
#  p.adjust.method='bonferroni'
#)
#cld(get_emmeans(pwc))

# Posthoc test
cld(emmeans(mod.rg.3, ~treatment:elapsed_dpi, type='response'), decreasing=T)
cld(emmeans(mod.rg.3, ~treatment, type='response'), decreasing=T)
emm <- emmeans(mod.rg.3, pairwise~treatment:elapsed_dpi, type='response')
plot(emm, comparison=T)
contrasts <- summary(emm$contrasts)
contrasts[which(contrasts$p.value <= 0.05), ]

p1 <- response.plot(harvest.nod[order(harvest.nod$ShootRG),], 'ShootRG', 'treatment', 'ShootRG', 'Co_effectivity') +
  #geom_text(aes(label=c(9,789,23456789,89,456789,12345678,456789,23456789,6789,6789,1234,6789,12345678,6789,1234567,
  #                      12345678,56789,1,12345678,23456789,123456,123,23456789,12345678,3456789,123456,
  #                      12345,1234567,23456789,12,23456789,789,23456789,89,6789,9), vjust=-2)) + 
  ggtitle('A. Symbiotic Benefit') + xlab('')



# Total nodules
hist(sqrt(harvest.nod$total_nodules))
mod.nod.3 <- lm(sqrt(total_nodules) ~ elapsed_dpi + treatment, data=harvest.nod)
summary(mod.nod.3)
shapiro.test(residuals(mod.nod.3))
autoplot(mod.nod.3)

summary(aov(mod.nod.3))
Anova(mod.nod.3, type='III', singular.ok=T)

# Posthoc
emm <- emmeans(mod.nod.3, pairwise~treatment:elapsed_dpi, type='response')
plot(emm, comparison=T)
contrasts <- summary(emm$contrasts)
contrasts[which(contrasts$p.value <= 0.05), ]

p2 <- response.plot(harvest.nod[order(harvest.nod$ShootRG),], 'total_nodules', 'treatment', 'ShootRG', 'Co_effectivity') +
  ggtitle('B. Nodulation') + xlab('')

grid.arrange(p1, p2, nrow=2)

# Mean weight
hist(log(harvest.nod$mean_weight+0.5))
mod.weight.3 <- lm(log(mean_weight+0.1) ~ elapsed_dpi + treatment, data=harvest.nod)
summary(mod.weight.3)
shapiro.test(residuals(mod.weight.3))
hist(residuals(mod.weight.3))
autoplot(mod.weight.3)

summary(aov(mod.weight.3))
Anova(aov(mod.weight.3), type='III')

# Posthoc
emm <- emmeans(mod.weight.3, pairwise~treatment:elapsed_dpi, type='response')
plot(emm, comparison=T)
contrasts <- summary(emm$contrasts)
dim(contrasts[which(contrasts$p.value <= 0.05), ])

p3 <- response.plot(harvest.nod[order(harvest.nod$ShootRG),], 'mean_weight', 'treatment', 'ShootRG', 'Co_effectivity') +
  ggtitle('C. Mean Fresh Nodule Weight') + xlab('')

# Investment
hist(harvest.nod$investment)
hist(log(harvest.nod$investment+0.2))
mod.invest.3 <- lm(log(investment+0.2) ~ elapsed_dpi + treatment, data=harvest.nod)
summary(mod.invest.3)
shapiro.test(residuals(mod.invest.3))
hist(residuals(mod.invest.3))
autoplot(mod.invest.3)

summary(aov(mod.invest.3))
Anova(mod.invest.3, type='III')

# Posthoc
emm <- emmeans(mod.invest.3, pairwise~treatment:elapsed_dpi, type='response')
plot(emm, comparison=T)
contrasts <- summary(emm$contrasts)
contrasts[which(contrasts$p.value <= 0.05), ]

p4 <- response.plot(harvest.nod[order(harvest.nod$ShootRG),], 'investment', 'treatment', 'ShootRG', 'Co_effectivity') +
  ggtitle('D. Investment to Nodulation') + xlab('Treatment')


grid.arrange(p3, p4, nrow=2)

## Testing N-fixing effectiveness
library(reshape2)

t.test(harvest.nod[which(harvest.nod$treatment == '186+187'), ]$ShootRG, mu=1, alternative = 'greater')

harvest.RG <- harvest.nod[,c('treatment','ShootRG')]
harvest.RG$exp <- 'Co'
harvest.RG[which(harvest.RG$treatment %in% c(2, 4, 131, 156, 184, 187, 187, 200)),]$exp <- 'Single'

t_test_summary <- data.frame(treatment=NA, p_val=NA)

for (t in unique(harvest.RG$treatment)){
  p_value <- t.test(harvest.RG[which(harvest.RG$treatment == t), ]$ShootRG, mu=1, alternative = 'greater')$p.value
  t_test_summary <- rbind(t_test_summary, data.frame(treatment=t, p_val=p_value))
}

t_test_summary$exp <- 'Co'
t_test_summary[which(t_test_summary$treatment %in% c(2, 4, 131, 156, 184, 187, 187, 200)),]$exp <- 'Single'
t_test_summary <- t_test_summary[-1,]

# subset
t_test_summary.single <- t_test_summary[which(t_test_summary$exp == 'Single'),]
t_test_summary.co <- t_test_summary[which(t_test_summary$exp != 'Single'),]

t_test_summary.single$p_val
p.adjust(t_test_summary.single$p_val, method="bonferroni")

t_test_summary.co$p_val
p.adjust(t_test_summary.co$p_val, method="bonferroni")
#
glm.invest <- glm(investment ~ treatment + harvest_batch, data=harvest.nod, family=poisson(link='log'))
summary(glm.invest)
autoplot(glm.invest)
shapiro.test(residuals(glm.invest))
Anova(glm.invest, type='III', singular.ok = T)
  
###########
# ANCOVA
###########
  
  # Testing linearity
  ggscatter(harvest.nod, x=elapsed_dpi, y=ShootRG, 
            color='treatment', add='reg.line') 
  
  # Testing interaction between covariate and categoriacal variable
  harvest.nod %>% anova_test(ShootRG ~ treatment*elapsed_dpi)
  
  # Testing normality of residuals
  lm.rg.4 <- lm(log(ShootRG) ~ elapsed_dpi + treatment, data=harvest.nod)
  model.metrics <- augment(lm.rg.4) %>%
    select(-.hat, -.sigma, -.fitted, -.se.fit)
  head(model.metrics, 3)
  
  # Assess normality of residuals
  shapiro_test(model.metrics$.resid)
  
  # Assess homogeneity of variences
  model.metrics %>% levene_test(.resid ~ treatment)
  
  # Analyzing ANCOVA
  res.aov <- harvest.nod %>% anova_test(log(ShootRG) ~ elapsed_dpi + treatment)
  get_anova_table(res.aov)
  
  # Post-hoc test
  library(emmeans)
  pwc <- harvest.nod %>%
    emmeans_test(
      ShootRG ~ treatment, covariate = elapsed_dpi,
      p.adjust.method='bonferroni'
    )
  pwc %>% arrange(p.adj.signif)
  
  #Adjusted means
  get_emmeans(pwc)
  
  #Reporting
  pwc <- pwc %>% add_xy_position(x='treatment', fun='mean_se')

    ggline(get_emmeans(pwc), x= 'treatment', y='emmean') +
  geom_errorbar(aes(ymin=conf.low, ymax=conf.high), width=0.2)
    +
  stat_pvalue_manual(pwc, hide.ns=T, tip.length=F) +
  labs(
    subtitle=get_test_label(res.aov, detailed=T), 
    caption=get_pwc_label(pwc)
  )

################################
# Coinoculation data analysis #
###############################
library(symbioUtil)
    
symbio.plot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
            'ShootRG', 'treatment', 'Fix')

symbio.plot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
            'ShootRG', 'treatment', 'Fix')

symbio.plot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
            'total_nodules', 'treatment', 'Fix')

symbio.plot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
            'mean_weight', 'treatment', 'Fix')

symbio.plot(harvest.coinoc.ordered[which(harvest.coinoc.ordered$Fix != 'Control'), ], 
            'investment', 'treatment', 'Fix')
# LM for 
mod.rg <- lm(log10(RG) ~ treatment + elapsed_dpi, data=harvest.single.ordered)
autoplot(mod.rg)
shapiro.test(residuals(mod.rg))


#######################################################
# Comparing single predicted with actual RG in coinoc #
#######################################################
single.response <- response.plot2(single.nod, 'ShootRG', 'treatment', 'ShootRG', 'Single_effectivity') + 
  geom_hline(yintercept = 1, color='coral')
single.response <- single.response$data

#factor in abundance
abundance = read.csv( '../data/competition_abundance.csv', header=T)
abundance[is.na(abundance)] <- 0
abundance$sum <- rowSums(abundance[, c('A', 'A.B', 'B')])
head(abundance)
abundance$percentA <- (abundance$A / abundance$sum) + ((abundance$A.B / abundance$sum)/2) 
abundance$percentB <- (abundance$B / abundance$sum) + ((abundance$A.B / abundance$sum)/2)


co.response.rg <- response.plot2(coinoc.nod, 'ShootRG', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
co.response.rg <- co.response.rg$data

co.response.rg$decouple <- strsplit(as.character(co.response.rg$x1), '[+]' )
co.response.rg$predict <- NA
co.response.rg$percent_predict <- NA
co.response.rg$percent_predictA <- NA
co.response.rg$percent_predictB <- NA



for (i in co.response.rg$decouple){
  a <- single.response[which(single.response$x1 == i[1]),]$y
  b <- single.response[which(single.response$x1 == i[2]),]$y
  t <- paste(i[1], i[2], sep='+' )
  co.response.rg[which(co.response.rg$x1 == t), ]$predict <- (a+b)/2
  # abundance weighted
  a_percent <- abundance[which(abundance$treatment == t), ]$percentA
  b_percent <- abundance[which(abundance$treatment == t), ]$percentB
  co.response.rg[which(co.response.rg$x1 == t), ]$percent_predict <- (a*a_percent+b*b_percent)
  co.response.rg[which(co.response.rg$x1 == t), ]$percent_predictA <- a*a_percent
  co.response.rg[which(co.response.rg$x1 == t), ]$percent_predictB <- b*b_percent
  
}

head(co.response.rg)  
co.response.rg$Fix <- c('Effective', 'Ineffective', 'Effective', 'Effective', 'Ineffective', 
                        'Ineffective', 'Effective', 'Effective', 'Ineffective', 'Effective', 
                        'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 'Effective',
                        'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                        'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                        'Effective', 'Effective', 'Effective')
co.response.rg$p <- c(0.05532, 1.18E-07, 0.7838, 0.6344, 0.02495, 0.948, 0.4636, 0.04459, 0.000892, 0.7219, 0.5681, 0.8965, 0.2928, 0.4187,
                   0.6981, 0.4803, 0.07239, 7.58E-05, 0.059, 5.18E-05, 0.06739, 0.0005016, 0.002652, 0.1788, 0.06682, 0.2787, 0.1822,
                   0.1033)

co.response.rg$p_abundance <- c(0.05532, 1.18E-07, 0.7838, 0.6344, 0.02495, 0.948, 0.4636, 0.04459, 0.000892, 0.7219, 0.5681, 0.8965, 0.2928, 0.4187,
                      0.6981, 0.4803, 0.07239, 7.58E-05, 0.059, 5.18E-05, 0.06739, 0.0005016, 0.002652, 0.1788, 0.06682, 0.2787, 0.1822,
                      0.1033)
                   
co.response.rg$optimize <- c('Optimum', 'Suboptimum','Optimum', 'Optimum','Suboptimum', 'Optimum', 'Optimum', 'Suboptimum',
                          'Suboptimum','Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum',
                          'Suboptimum', 'Optimum', 'Suboptimum', 'Optimum', 'Suboptimum', 'Suboptimum', 'Optimum', 'Optimum',
                          'Optimum', 'Optimum', 'Optimum')

# T-test for observed co-inoculation values significantly different from %predicted RG or not
#t = '156+200'
#test.group <- coinoc.nod[which(coinoc.nod$treatment == t), 'RG']
#test = co.response.rg[which(co.response.rg$x1 == t),  'percent_predict']
#t.out <- t.test(test.group, mu=test, alternative="two.sided")

co.response.rg$percent_tstat <- NA
co.response.rg$percent_pvalue <- NA
co.response.rg$percent_confint <- NA


for (t in co.response.rg$x1){
  test.group <- coinoc.nod[which(coinoc.nod$treatment == t), 'RG']
  test = co.response.rg[which(co.response.rg$x1 == t),  'percent_predict']
  t.out <- t.test(test.group, mu=test, alternative="less")
  t.stat <- t.out$statistic
  t.pval <- t.out$p.value
  co.response.rg[which(co.response.rg$x1 == t),]$percent_tstat  <- t.stat
  co.response.rg[which(co.response.rg$x1 == t),]$percent_pvalue <- t.pval
  co.response.rg[which(co.response.rg$x1 == t), ]$percent_predictB <- b*b_percent
  #co.response.rg$percent_confint <- t.out$conf.int
  co.response.rg[which(co.response.rg$x1 == t), ]$percent_confint <- t.out$stderr
  }

p0 <- ggplot(co.response.rg, aes(x=predict, y=y)) + 
  #geom_point(aes(shape=Fix)) + 
  geom_point(size=3, aes( color=ifelse(p <= 0.05, 'green', 'red')), show.legend = T) + 
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  scale_color_manual(labels = c('Suboptimal', 'Optimal'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Effective', 'Ineffective'), values=c(3, 16, 17)) + labs(colour='Interaction') +
    xlab('Expected Relative Growth') + ylab('Observed Relative Growth') + 
  ggtitle('Relative growth/response comparison') +
  geom_abline(slope=1, intercept=0)

p1 <- ggplot(co.response.rg, aes(x=percent_predict, y=y)) + 
  #geom_point(aes(shape=Fix)) + 
  geom_point(size=3, aes( color=ifelse(percent_pvalue <= 0.05, 'green', 'red')), show.legend = T) + 
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  scale_color_manual(labels = c('Suboptimal', 'Optimal'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Effective', 'Ineffective'), values=c(3, 16, 17)) + labs(colour='Interaction') +
  xlab('Abundance Weighted Relative Growth') + ylab('Observed Relative Growth') + 
  ggtitle('Abundance Weighted Relative growth/response comparison') +
  geom_abline(slope=1, intercept=0)

grid.arrange(p0, p1, nrow=1)

  
###########################################
# Fitting model for abundance weighted RG #
###########################################
head(co.response.rg)
head(coinoc.nod)
head(abundance)

# Merge co.response.rg and abundance
abundance_rg <- merge(co.response.rg, abundance, by.x = 'x1', by.y='treatment')
mod = lm(o ~ percent_predict, data=abundance_rg)
summary(mod)
mod$coefficients
autoplot(mod)
shapiro.test(residuals(mod))

p1 <- p1 + geom_abline(intercept=0.8129, slope=0.37177, color="#E41A1C")

grid.arrange(p0, p1, nrow=2)

# Misc
t.test(coinoc.nod[which(coinoc.nod$treatment=='4+200'), ]$ShootRG, mu=2.0080448)
#pw.comb <- combn(single.response$data$x1, 2)
#paste(as.vector(pw.comb[,1]), sep='+')
dim(coinoc.nod)
mod.co.seperate <- lm(log10(ShootRG) ~ elapsed_dpi +  Single_effectivity, data=coinoc.nod)
summary(mod.co.seperate)

autoplot(mod.co.seperate)

shapiro.test(residuals(mod.co.seperate))
Anova(mod.co.seperate, type='III')

posthocs <- glht(mod.co.seperate, linfct=mcp(Single_effectivity='Tukey'))
summary(posthocs)

#
mod.co.single <- lm(log10(ShootRG) ~ elapsed_dpi +  Co_effectivity, data=coinoc.nod)
summary(mod.co.single)
autoplot(mod.co.single)

shapiro.test(residuals(mod.co.single))
Anova(mod.co.single, type='III')

posthocs <- glht(mod.co.single, linfct=mcp(Co_effectivity='Tukey'))
summary(posthocs)



############################################################
# Comparing single predicted with actual nodules in coinoc #
############################################################

single.response <- response.plot2(single.nod, 'total_nodules', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
single.response <- single.response$data

co.response.nod <- response.plot2(coinoc.nod, 'total_nodules', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
co.response.nod <- co.response.nod$data

co.response.nod$decouple <- strsplit(as.character(co.response.nod$x1), '[+]' )
co.response.nod$predict <- NA
co.response.nod$percent_predict <- NA
#co.response[co.response$decouple[[1]][2], 'y']



for (i in co.response.nod$decouple){
  a <- single.response[which(single.response$x1 == i[1]),]$y
  b <- single.response[which(single.response$x1 == i[2]),]$y
  t <- paste(i[1], i[2], sep='+' )
  a_percent <- abundance[which(abundance$treatment == t), ]$percentA
  b_percent <- abundance[which(abundance$treatment == t), ]$percentB
  
  co.response.nod[which(co.response.nod$x1 == t), ]$predict <- (a+b)/2
  co.response.nod[which(co.response.nod$x1 == t), ]$percent_predict <- (a*a_percent)+(b*b_percent)
}

head(co.response.nod)  
co.response.nod$Fix <- c('Effective', 'Ineffective', 'Effective', 'Effective', 'Ineffective', 
                         'Ineffective', 'Effective', 'Effective', 'Ineffective', 'Effective', 
                         'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 'Effective',
                         'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                         'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                         'Effective', 'Effective', 'Effective')
# t-test
co.response.nod$p <- NA
co.response.nod$t_p_percent <- NA
co.response.nod$t_stat_percent <- NA
co.response.nod$t_confint_percent <- NA

for (i in co.response.nod$x1){
  predicted <- co.response.nod[which(co.response.nod$x1 == i), 'predict']
  percent_predicted <- co.response.nod[which(co.response.nod$x1 == i), 'percent_predict']
  p <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$total_nodules, mu=predicted)$p.value
  t_test <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$total_nodules, mu=percent_predicted)
  t_p_percent <- t_test$p.value
  t_stat_percent <- t_test$stat
  #t_confint_percent <- t_test$conf.int[1]
  t_confint_percent <- t_test$stderr
  co.response.nod[which(co.response.nod$x1 == i), 'p'] <- p
  co.response.nod[which(co.response.nod$x1 == i), 't_p_percent'] <- t_p_percent
  co.response.nod[which(co.response.nod$x1 == i), 't_stat_percent'] <- t_stat_percent
  co.response.nod[which(co.response.nod$x1 == i), 't_confint_percent'] <- t_confint_percent
}


#co.response.nod$optimize <- c('Optimum', 'Suboptimum','Optimum', 'Optimum','Suboptimum', 'Optimum', 'Optimum', 'Suboptimum',
#                          'Suboptimum','Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum',
#                          'Suboptimum', 'Optimum', 'Suboptimum', 'Optimum', 'Suboptimum', 'Suboptimum', 'Optimum', 'Optimum',
#                         'Optimum', 'Optimum', 'Optimum')



p1 <- ggplot(co.response.nod, aes(x=predict, y=y)) + 
  #geom_point(aes(shape=optimize)) + 
  geom_point(size=3, aes(shape=optimize, color=ifelse(p <= 0.05, 'green', 'red')), show.legend = T) + 
  scale_color_manual(labels = c('Different', 'No Difference'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Optimal', 'Suboptimal'), values=c(3, 16)) + labs(colour='Significance', shape='Interaction') +
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  xlab('Expected Nodule Number') + ylab('Observed Nodule Number') +
  ggtitle('A. Nodulation comparison') +
  geom_abline(slope=1, intercept=0)



nod.tradeoff <- co.response.rg

# Checking trade-off
tradeoff <- merge(RG.tradeoff, nod.tradeoff,  by=c('x1', 'decouple', 'Fix', 'optimize', 'Single_effectivity'))
colnames(tradeoff) <- c('treatment', 'pairs', 'Fix', 'optimize', 'RG', 'RGse', 'oRG', 'predictedRG', 'pRG','Nod',
                                'Nodse', 'oNod', 'predictedNod', 'pNod')

head(tradeoff)


library(ggpmisc)

ggplot(tradeoff, aes(x=Nod, y=RG, color=optimize)) + 
  geom_point() + geom_smooth(method='lm', formula=y~x) +
  stat_smooth_func(geom='text', method='lm', hjust=0, parse=T)

mod.optimum <- lm(RG~Nod, data=tradeoff[which(tradeoff$optimize == 'Optimum'), ])
mod.suboptimum <- lm(RG~Nod, data=tradeoff[which(tradeoff$optimize == 'Suboptimum'), ])

Anova(mod.optimum, mod.suboptimum, type='III')

##################################################################
# Comparing single predicted with observed mean_weight in coinoc #
##################################################################

single.response <- response.plot(single.nod, 'mean_weight', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
single.response <- single.response$data

co.response.mw <- response.plot(coinoc.nod, 'mean_weight', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
co.response.mw <- co.response.mw$data

co.response.mw$decouple <- strsplit(as.character(co.response.mw$x1), '[+]' )
co.response.mw$predict <- NA
co.response.mw$percent_predict <- NA
#co.response[co.response$decouple[[1]][2], 'y']



for (i in co.response.mw$decouple){
  a <- single.response[which(single.response$x1 == i[1]),]$y
  b <- single.response[which(single.response$x1 == i[2]),]$y
  a_percent <- abundance[which(abundance$treatment == t), ]$percentA
  b_percent <- abundance[which(abundance$treatment == t), ]$percentB
  t <- paste(i[1], i[2], sep='+' )
  co.response.mw[which(co.response.mw$x1 == t), ]$predict <- (a+b)/2
  co.response.mw[which(co.response.nod$x1 == t), ]$percent_predict <- (a*a_percent)+(b*b_percent)
}

head(co.response.mw)  
co.response.mw$Fix <-c('Effective', 'Ineffective', 'Effective', 'Effective', 'Ineffective', 
                       'Ineffective', 'Effective', 'Effective', 'Ineffective', 'Effective', 
                       'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 'Effective',
                       'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                       'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                       'Effective', 'Effective', 'Effective')

# t-test
co.response.mw$p <- NA
co.response.mw$t_p_percent <- NA
co.response.mw$t_stat_percent <- NA
co.response.mw$t_confint_percent <- NA

for (i in co.response.mw$x1){
  predicted <- co.response.mw[which(co.response.mw$x1 == i), 'predict']
  p <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$mean_weight, mu=predicted)$p.value
  t_test <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$mean_weight, mu=percent_predicted)
  t_p_percent <- t_test$p.value
  t_stat_percent <- t_test$stat
  #t_confint_percent <- t_test$conf.int[1]
  t_confint_percent <- t_test$stderr
  co.response.mw[which(co.response.mw$x1 == i), 'p'] <- p
  co.response.mw[which(co.response.mw$x1 == i), 't_p_percent'] <- t_p_percent
  co.response.mw[which(co.response.mw$x1 == i), 't_stat_percent'] <- t_stat_percent
  co.response.mw[which(co.response.mw$x1 == i), 't_confint_percent'] <- t_confint_percent
}

co.response.mw$optimize <- c('Optimum', 'Suboptimum','Optimum', 'Optimum','Suboptimum', 'Optimum', 'Optimum', 'Suboptimum',
                          'Suboptimum','Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum',
                          'Suboptimum', 'Optimum', 'Suboptimum', 'Optimum', 'Suboptimum', 'Suboptimum', 'Optimum', 'Optimum',
                          'Optimum', 'Optimum', 'Optimum')


ggplot(co.response.mw, aes(x=predict, y=y)) + 
  #geom_point(aes(shape=optimize)) + 
  geom_point(size=3, aes(shape=optimize, color=ifelse(p <= 0.05, 'green', 'red')), show.legend = T) + 
  scale_color_manual(labels = c('Different', 'No Difference'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Optimal', 'Suboptimal'), values=c(3, 16)) + labs(colour='Significance', shape='Interaction') +
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  xlab('Expected Nodule Number') + ylab('Observed Nodule Number') +
  ggtitle('Mean fresh nodule weight compare') +
  geom_abline(slope=1, intercept=0)



#################################################################
# Comparing single predicted with observed investment in coinoc #
#################################################################

single.response <- response.plot(single.nod, 'investment', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
single.response <- single.response$data

co.response.inv <- response.plot(coinoc.nod, 'investment', 'treatment', 'ShootRG', 'Single_effectivity') + geom_hline(yintercept = 1, color='coral')
co.response.inv <- co.response.inv$data

co.response.inv$decouple <- strsplit(as.character(co.response.inv$x1), '[+]' )
co.response.inv$predict <- NA
#co.response[co.response$decouple[[1]][2], 'y']



for (i in co.response.inv$decouple){
  a <- single.response[which(single.response$x1 == i[1]),]$y
  b <- single.response[which(single.response$x1 == i[2]),]$y
  t <- paste(i[1], i[2], sep='+' )
  co.response.inv[which(co.response.inv$x1 == t), ]$predict <- (a+b)/2
}

head(co.response.inv)  
co.response.inv$Fix <- c('Effective', 'Ineffective', 'Effective', 'Effective', 'Ineffective', 
                         'Ineffective', 'Effective', 'Effective', 'Ineffective', 'Effective', 
                         'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 'Effective',
                         'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                         'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                         'Effective', 'Effective', 'Effective')


# t-test
co.response.inv$p <- NA

for (i in co.response.inv$x1){
  predicted <- co.response.inv[which(co.response.inv$x1 == i), 'predict']
  p <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$investment, mu=predicted)$p.value
  co.response.inv[which(co.response.inv$x1 == i), 'p'] <- p
}



co.response.inv$optimize <- c('Optimum', 'Suboptimum','Optimum', 'Optimum','Suboptimum', 'Optimum', 'Optimum', 'Suboptimum',
                          'Suboptimum','Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum',
                          'Suboptimum', 'Optimum', 'Suboptimum', 'Optimum', 'Suboptimum', 'Suboptimum', 'Optimum', 'Optimum',
                          'Optimum', 'Optimum', 'Optimum')

p3 <- ggplot(co.response.inv, aes(x=predict, y=y)) + 
  #geom_point(aes(shape=optimize)) + 
  geom_point(size=3, aes(shape=optimize, color=ifelse(p <= 0.05, 'green', 'red')), show.legend = T) + 
  scale_color_manual(labels = c('Different', 'No Difference'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Optimal', 'Suboptimal'), values=c(3, 16)) + labs(colour='Significance', shape='Interaction') +
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  xlab('Expected Investment') + ylab('Observed Investment') +
  ggtitle('C. Investment comparison') +
  geom_abline(slope=1, intercept=0)


##################################################
# Comparing nodule area: prediction vs. observed #
##################################################
library(dplyr)
library(stringr)

#nod.clean.single <- nod.clean %>% filter(str_detect(treatment, "\\+", negate=T))
#nod.clean.coinoc <- nod.clean %>% filter(str_detect(treatment, "\\+"))

single.response <- response.plot2(single.nod, 'AreaMM', 'treatment', 'ShootRG', 'Single_effectivity')
single.response <- single.response$data

co.response.area <- response.plot2(coinoc.nod, 'AreaMM', 'treatment', 'ShootRG', 'Single_effectivity')
co.response.area <- co.response.area$data

co.response.area$decouple <- strsplit(as.character(co.response.area$x1), '[+]' )
co.response.area$predict <- NA
co.response.area$percent_predict <- NA
#co.response[co.response$decouple[[1]][2], 'y']



for (i in co.response.area$decouple){
  a <- single.response[which(single.response$x1 == i[1]),]$y
  b <- single.response[which(single.response$x1 == i[2]),]$y
  t <- paste(i[1], i[2], sep='+' )
  a_percent <- abundance[which(abundance$treatment == t), ]$percentA
  b_percent <- abundance[which(abundance$treatment == t), ]$percentB
  
  co.response.area[which(co.response.area$x1 == t), ]$predict <- (a+b)/2
  co.response.area[which(co.response.area$x1 == t), ]$percent_predict <- (a*a_percent)+(b*b_percent)
}

co.response.area$Fix <- c('Effective', 'Ineffective', 'Effective', 'Effective', 'Ineffective', 
                          'Ineffective', 'Effective', 'Effective', 'Ineffective', 'Effective', 
                          'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 'Effective',
                          'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                          'Ineffective', 'Ineffective', 'Ineffective', 'Effective', 'Ineffective', 
                          'Effective', 'Effective', 'Effective')

# t-test
co.response.area$p <- NA
co.response.area$t_p_percent <- NA
co.response.area$t_stat_percent <- NA
co.response.area$t_confint_percent <- NA


for (i in co.response.area$x1){
  predicted <- co.response.area[which(co.response.area$x1 == i), 'predict']
  p <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$AreaMM, mu=predicted)$p.value
  t_test <- t.test(coinoc.nod[which(coinoc.nod$treatment==i), ]$AreaMM, mu=predicted)
  t_p_percent <- t_test$p.value
  t_stat_percent <- t_test$stat
  #t_confint_percent <- t_test$conf.int[1]
  t_confint_percent <- t_test$stderr
  co.response.area[which(co.response.area$x1 == i), 'p'] <- p
  co.response.area[which(co.response.area$x1 == i), 't_p_percent'] <- t_p_percent
  co.response.area[which(co.response.area$x1 == i), 't_stat_percent'] <- t_stat_percent
  co.response.area[which(co.response.area$x1 == i), 't_confint_percent'] <- t_confint_percent
}

co.response.area$optimize <- c('Optimum', 'Suboptimum','Optimum', 'Optimum','Suboptimum', 'Optimum', 'Optimum', 'Suboptimum',
                          'Suboptimum','Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum', 'Optimum',
                          'Suboptimum', 'Optimum', 'Suboptimum', 'Optimum', 'Suboptimum', 'Suboptimum', 'Optimum', 'Optimum',
                          'Optimum', 'Optimum', 'Optimum')



head(co.response.area)  

p4 <- ggplot(co.response.area, aes(x=predict, y=y)) + 
  geom_point(size=3, aes(shape=optimize, color=ifelse(p <= 0.05, 'green', 'red')), show.legend = T) + 
  scale_color_manual(labels = c('Different', 'No Difference'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Optimal', 'Suboptimal'), values=c(3, 16)) + labs(colour='Significance', shape='Interaction') +
  #geom_point(size=3, aes(shape=optimize)) + 
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  xlab('Expected Mean Nodule Area') + ylab('Observed Mean Nodule Area') +
  ggtitle('D. Nodule Area comparison') +
  geom_abline(slope=1, intercept=0)



p2 <- ggplot(co.response.mw, aes(x=predict, y=y)) + 
  #geom_point(aes(shape=optimize)) + 
  geom_point(size=3, aes(shape=optimize, color=ifelse(p <= 0.05, 'green', 'red')), show.legend = T) + 
  scale_color_manual(labels = c('Different', 'No difference'), values=c('tomato2', 'turquoise3')) +
  scale_shape_manual(labels = c('Optimal', 'Suboptimal'), values=c(3, 16)) + labs(colour='Significance', shape='Interaction') +
  geom_label_repel(aes(label=x1), box.padding=0.35, point.padding = 0.5) +
  xlab('Expected Mean Fresh Nodule Weight') + ylab('Observed Mean Fresh Nodule Weight') +
  ggtitle('B. Fresh nodule weight comparison') +
  geom_abline(slope=1, intercept=0)

grid.arrange(p1, p2, p3, p4, nrow=2)

###########################################
# Compare observed vs expected difference #
###########################################
# Main result figure
co.response.rg$diff <- co.response.rg$y - co.response.rg$percent_predict
co.response.rg$x1 <- factor(co.response.rg$x1, levels=c('4+200', '186+4', '187+200', '4+187', 
                                                     '2+187', '186+187', '187+184', '186+200', 
                                                     '156+200', '2+186', '187+156', '4+156',  
                                                     '156+184', '2+200', '4+184', '186+156', 
                                                     '2+156', '184+200', '131+156', '186+184', 
                                                     '2+4', '131+187', '4+131', '2+184', '131+200', 
                                                     '186+131', '2+131',  '131+184'))
  
a0 <- ggplot(co.response.rg[which(co.response.rg$x1 != '131+156'),], aes(x=x1, y=diff, fill=f)) + 
  #geom_bar(stat='identity', aes( fill =ifelse(percent_pvalue <= 0.05, 'green', 'red'))) +
  #geom_bar(stat='identity', fill='#228B22') +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=diff - percent_confint, ymax=diff+percent_confint), color='grey', width=0.2) +
  geom_text(aes(label=ifelse(percent_pvalue <= 0.05, '*', '')),color='red', position = position_dodge(width = .9), vjust = 1, size = 20 / .pt) +
  scale_fill_manual(labels = c("Fix+/Fix+", "Fix+/Fix-", "Fix-/Fix-"), values=c('#228B22', '#0096FF', '#F8766D')) +
  #scale_fill_manual(labels = c('Suboptimal', 'Optimal'), values=c('tomato2', 'turquoise3')) +
  labs(fill='Performance')+
  ylab('Shoot RG Difference') + xlab('Treatments') +
  ggtitle('A') + theme_bw() + theme(legend.position='none') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 

co.response.nod$diff <- co.response.nod$y - co.response.nod$percent_predict
co.response.nod$x1 <- factor(co.response.nod$x1, levels=c('4+200', '186+4', '187+200', '4+187', 
                                                        '2+187', '186+187', '187+184', '186+200', 
                                                        '156+200', '2+186', '187+156', '4+156',  
                                                        '156+184', '2+200', '4+184', '186+156', 
                                                        '2+156', '184+200', '131+156', '186+184', 
                                                        '2+4', '131+187', '4+131', '2+184', '131+200', 
                                                        '186+131', '2+131',  '131+184'))


a1 <- ggplot(co.response.nod[which(co.response.nod$x1 != '131+156'),], aes(x=x1, y=diff, fill=f)) + 
  #geom_bar(stat='identity', aes( fill =ifelse(t_p_percent <= 0.05, 'green', 'red'))) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=diff - t_confint_percent, ymax= diff + t_confint_percent), width=.2, colour='grey') +
  geom_text(aes(label=ifelse(t_p_percent <= 0.05, '*', '')),color='red', position = position_dodge(width = .9), vjust = 1, size = 20 / .pt) +
  scale_fill_manual(labels = c("Fix+/Fix+", "Fix+/Fix-", "Fix-/Fix-"), values=c('#228B22', '#0096FF', '#F8766D')) +
  labs(fill='Difference')+
  ylab('Nodule Number Difference') + xlab('Treatments') +
  ggtitle('B') + theme_bw() + theme(legend.position='none') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 



co.response.mw$diff <- co.response.mw$y - co.response.mw$percent_predict
co.response.mw$x1 <- factor(co.response.mw$x1, levels=c('4+200', '186+4', '187+200', '4+187', 
                                                        '2+187', '186+187', '187+184', '186+200', 
                                                        '156+200', '2+186', '187+156', '4+156',  
                                                        '156+184', '2+200', '4+184', '186+156', 
                                                        '2+156', '184+200', '131+156', '186+184', 
                                                        '2+4', '131+187', '4+131', '2+184', '131+200', 
                                                        '186+131', '2+131',  '131+184'))

a2 <- ggplot(co.response.mw[which(co.response.mw$x1 != '131+156'),], aes(x=x1, y=diff, fill=f)) + 
  #geom_bar(stat='identity', aes( fill =ifelse(t_p_percent <= 0.05, 'green', 'red'))) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=diff - t_confint_percent, ymax= diff + t_confint_percent), width=.2, colour='grey') +
  geom_text(aes(label=ifelse(t_p_percent <= 0.05, '*', '')),color='red', position = position_dodge(width = .9), vjust = 1, size = 20 / .pt) +
  scale_fill_manual(labels = c("Fix+/Fix+", "Fix+/Fix-", "Fix-/Fix-"), values=c('#228B22', '#0096FF', '#F8766D')) +
  labs(fill='Difference')+
  ylab('Difference') + xlab('Treatments') +
  ggtitle('D. Observed mean weight - Expected mean weight')  + theme_bw() + theme(legend.position='none') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 

co.response.area$diff <- co.response.area$y - co.response.area$percent_predict
co.response.area$x1 <- factor(co.response.area$x1, levels=c('4+200', '186+4', '187+200', '4+187', 
                                                        '2+187', '186+187', '187+184', '186+200', 
                                                        '156+200', '2+186', '187+156', '4+156',  
                                                        '156+184', '2+200', '4+184', '186+156', 
                                                        '2+156', '184+200', '131+156', '186+184', 
                                                        '2+4', '131+187', '4+131', '2+184', '131+200', 
                                                        '186+131', '2+131',  '131+184'))

a3 <- ggplot(co.response.area[which(co.response.area$x1 != '131+156'),], aes(x=x1, y=diff, fill=f)) + 
  #geom_bar(stat='identity', aes( fill =ifelse(t_p_percent <= 0.05, 'green', 'red'))) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=diff - t_confint_percent, ymax= diff + t_confint_percent), width=.2, colour='grey') +
  geom_text(aes(label=ifelse(t_p_percent <= 0.05, '*', '')),color='red', position = position_dodge(width = .9), vjust = 1, size = 20 / .pt) +
  scale_fill_manual(labels = c("Fix+/Fix+", "Fix+/Fix-", "Fix-/Fix-"), values=c('#228B22', '#0096FF', '#F8766D')) +
  labs(fill='Difference')+ ylab('Nodule Area Difference') + xlab('Treatments') +
  ggtitle('C')  + theme_bw() + theme(legend.position='none') +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 

#grid.arrange(a0, a1, a3, a2, nrow=4)
ggarrange(a0, a1, a3, nrow=3, common.legend = T, legend='top')

# Mean fresh weight
co.response.mw$diff <- co.response.mw$y - co.response.mw$percent_predict
co.response.mw$x1 <- factor(co.response.mw$x1, levels=c('4+200', '186+4', '187+200', '4+187', 
                                                            '2+187', '186+187', '187+184', '186+200', 
                                                            '156+200', '2+186', '187+156', '4+156',  
                                                            '156+184', '2+200', '4+184', '186+156', 
                                                            '2+156', '184+200', '131+156', '186+184', 
                                                            '2+4', '131+187', '4+131', '2+184', '131+200', 
                                                            '186+131', '2+131',  '131+184'))

a4 <- ggplot(co.response.mw, aes(x=x1, y=diff)) + 
  geom_bar(stat='identity', aes( fill =ifelse(t_p_percent <= 0.05, 'green', 'red'))) +
  geom_errorbar(aes(ymin=diff - t_confint_percent, ymax= diff + t_confint_percent), width=.2, colour='grey') +
  scale_fill_manual(labels = c('Significant', 'Non-significant'), values=c('tomato2', 'turquoise3')) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + labs(fill='Difference')+
  ylab('Difference') + xlab('') +
  ggtitle('D. Observed nodule mean weight - Expected nodule mean weight') +
  theme(plot.title = element_text(face = "bold")) 


grid.arrange(a0, a1, a3, a4,  nrow=4)

########################
# Peer-review analysis #
########################

# total biomass
harvest$total_biomass <- harvest$root_mass + harvest$root_mass

# total nodule biomass
harvest$total_nodmass <- harvest$mean_weight*harvest$total_nodules

# per nodule biomass
harvest$per_nod_benefit <- harvest$ShootRG/harvest$total_nodules

treated <- harvest %>%
  filter(treatment != "H2O") %>%
  droplevels(.)

## rename Fix treatment
treated$type <- ifelse(treated$Fix == "+", "Single_Effective", 
                       ifelse(treated$Fix == "-", "Single_Ineffective",
                              ifelse(treated$Fix == "+/+", "Co_Effective",
                                     ifelse(treated$Fix == "+/-", "Co_Mixed",
                                            "Co_Ineffective"))))

treated$type <- factor(treated$type, levels = c("Single_Effective", "Co_Effective",
                                                "Co_Mixed", "Single_Ineffective", "Co_Ineffective"))

### correlation between shoot_mass and total nodule weight
ggplot(treated, aes(x=total_nodules, y = ShootRG)) +
  geom_point(aes(color = type)) +
  geom_smooth(method="lm", formula = "y~x")

### get residuals
lm <- lm(log(ShootRG) ~ total_nodules*type, data = treated)
summary(lm) ## sig
aov <- aov(lm)
summary(aov)



res_investment <- as.data.frame(lm$residuals)
res_investment$ID <- row.names(res_investment)
colnames(res_investment) <- c("res","ID")

### add residuals to master
treated$ID <- row.names(treated)
treated <- left_join(treated, res_investment, by = "ID")

## correlations between shoot_mass and all other traits
traits <- c("root_mass","shoot_mass", "total_nodules", "mean_weight",
            "RG","ShootRG","investment","AreaMM","dryMass","total_biomass",
            "total_nodmass","res","per_nod_benefit")

## gather traits to be correlated
corrs <- treated %>%
  select(treatment, plant_id, all_of(traits)) %>%
  gather(., "trait", "value", -treatment, -ShootRG, -plant_id)

ggplot(corrs, aes(x=value, y = ShootRG)) +
  geom_point() +
  geom_smooth(method = "loess", formula = "y~x") +
  facet_wrap(~trait, scales = "free_x") +
  theme_bw()

## filter to just single inoculation treatments
single <- treated %>%
  filter(treatment %in% c("131", "156", "4", "184", 
                          "186", "187", "2", "200")) %>%
  droplevels(.)

### reorder treatments
single$treatment <- factor(single$treatment, 
                           levels = c("131", "156", "4", "184", 
                                      "186", "187", "2", "200")
)

### corrplot between relative shoot growth and total nodule mass
ggplot(single, aes(x=total_nodules, y = log(ShootRG), color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", size = 2) +
  theme_bw()

### corrplot between relative shoot growth and total nodule mass
ggplot(single, aes(x=total_nodmass, y = log(ShootRG), color = treatment)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", size = 2) +
  theme_bw()

### ANCOVA
lm <- lm(log(ShootRG) ~ total_nodules*treatment, data = single)
lm_sum <- summary(lm)
(aov <- Anova(lm, type = 3)) ## sig

## specify traits
traits <- c("root_mass","shoot_mass", "total_nodules", "mean_weight",
            "RG","ShootRG","investment","AreaMM","dryMass","total_biomass",
            "total_nodmass")

### get means
single_sum <- single %>%
  group_by(treatment) %>%
  summarize_at(traits, mean, na.rm = TRUE)

## calculate benefit per nod
single_sum$per_nod_benefit <- single_sum$ShootRG/single_sum$total_nodules

## quick graph
ggplot(single_sum, aes(x=reorder(treatment, per_nod_benefit, median),
                       y = per_nod_benefit)) +
  geom_bar(stat = "identity")

## add in slopes
all_coef <- as.data.frame(lm_sum$coefficients)
treatment <- single_sum$treatment
###
treat131 <- all_coef["total_nodules","Estimate"]
treat156 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment156","Estimate"]
treat4 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment4","Estimate"]
treat184 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment184","Estimate"]
treat186 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment186","Estimate"]
treat187 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment187","Estimate"]
treat2 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment2","Estimate"]
treat200 <- all_coef["total_nodules","Estimate"] + 
  all_coef["total_nodules:treatment200","Estimate"]
slopes <- c(treat131, treat156, treat4, treat184,
            treat186, treat187, treat2, treat200)

single_sum$slopes <- slopes

## specify
traits <- c("root_mass","shoot_mass", "total_nodules", "mean_weight",
            "RG","ShootRG","investment","AreaMM","dryMass","total_biomass",
            "total_nodmass")

## gather traits
treated.l <- treated %>%
  select(type, plant_id, treatment, all_of(traits)) %>%
  gather(., "trait", "value", -type, -plant_id, -treatment)

## plot of all traits
ggplot(treated.l, aes(x=type, y = value)) +
  geom_boxplot() +
  geom_jitter(aes(color = treatment)) +
  facet_wrap(~trait, scales = "free_y") +
  theme_bw() +
  guides(color = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)
  )

## test whether shootRG is highest in + vs +/-
lm <- lm(log(ShootRG) ~ type, data = treated)
aov <- aov(lm)
summary(aov)
TukeyHSD(aov)

## test whether total_nodules is highest in + vs +/-
lm <- lm(sqrt(total_nodules) ~ type, data = treated)
aov <- aov(lm)
summary(aov)
TukeyHSD(aov)

## plot of shootRG and total_nodles
p <- ggplot(treated.l %>% filter(trait %in% c("ShootRG","total_nodules")), 
            aes(x=type, y = value)) +
  geom_boxplot() +
  geom_jitter(aes(color = treatment)) +
  facet_wrap(~trait, scales = "free_y") +
  theme_bw() +
  guides(color = "none") +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 0.5, vjust = 0.5)
  )

#save(single_sum, file = "./single_sum.Rdata")

### corrplot between relative shoot growth and total nodule mass
ggplot(single_sum, aes(x=total_nodmass, y = log(ShootRG))) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x", size = 2) +
  theme_bw()

## Nodule occupancy freq
nod_occ <- read_csv('../data/competition_abundance.csv')

## total nodules
nod_occ$total_nods <- rowSums(nod_occ[,c("A","B","A+B")], na.rm = TRUE)

## make dup treatment
nod_occ$treatment2 <- nod_occ$treatment

## split strains
nod_occ <- nod_occ %>%
  separate(treatment2, c("strain_A", "strain_B"))

## replaces NAs with zeros
nod_occ <- replace(nod_occ, is.na(nod_occ), 0)

## frequency
nod_occ <- nod_occ %>%
  rowwise() %>%
  mutate(prop_A_single = (A/total_nods),
         prop_A_mixed = `A+B`/(2*total_nods),
         fA = prop_A_single + prop_A_mixed,
         fB = 1 - fA)

head(nod_occ)

## filter to co-inoculation treatments
mixed <- treated %>%
  filter(!treatment %in% c("131", "156", "4", "184", 
                           "186", "187", "2", "200")) %>%
  droplevels(.)

### add in fA and fB given treatments
mixed <- left_join(mixed, nod_occ[,c("strain_A", "strain_B", "fA","fB","treatment")], 
                   by = "treatment")

### expect values from single inoculation
mixed$exp_strain_A <- single_sum$ShootRG[match(mixed$strain_A, 
                                               single_sum$treatment)]
mixed$exp_strain_B <- single_sum$ShootRG[match(mixed$strain_B, 
                                               single_sum$treatment)]

### calculate according to their paper
mixed$exp <- ((mixed$exp_strain_A)*(mixed$fA)) + ((mixed$exp_strain_B)*(mixed$fB))

### t.test for each treatment separately
treat.list <- unique(mixed$treatment)

t.test_fun <- function(treat, df){
  
  ## filter for each treatment
  mixed_treat <- df %>%
    filter(treatment %in% treat) %>%
    droplevels(.)
  
  predicted <- mean(mixed_treat$exp, na.rm = TRUE)
  
  ttest <- t.test(mixed_treat$ShootRG, mu = predicted)
  
  return(ttest)
}

out <- sapply(treat.list, FUN = t.test_fun, df = mixed,
              simplify = FALSE, USE.NAMES = TRUE)

### expect values from single inoculation
mixed$exp2_strain_A <- single_sum$per_nod_benefit[match(mixed$strain_A, 
                                                        single_sum$treatment)]
mixed$exp2_strain_B <- single_sum$per_nod_benefit[match(mixed$strain_B, 
                                                        single_sum$treatment)]

### total nodules of each strain
mixed$total_nods_A <- mixed$fA * mixed$total_nodules
mixed$total_nods_B <- mixed$fB * mixed$total_nodules

### calculate according to per benefit 
# Proposed by reviewr 3
mixed$exp2 <- ((mixed$exp2_strain_A)*(mixed$total_nods_A)) + ((mixed$exp2_strain_B)*(mixed$total_nods_B))

## compare graphically
mixed.l <- mixed %>%
  select(ShootRG, exp, exp2, treatment) %>%
  gather(., method, value, -treatment, -ShootRG)

p0 <- ggplot(mixed.l[which(mixed.l$method=='exp'),], aes(x=value, y = ShootRG, color = method)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_discrete(breaks = c("exp","exp2"),
                       labels = c("Model I",
                                  "Model II")) +
  ylab("Observed relative shoot growth") +
  xlab("Expected relative shoot growth") +
  theme_bw() + theme(legend.title=element_blank()) 

p <- ggplot(mixed.l, aes(x=value, y = ShootRG, color = method)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  scale_color_discrete(breaks = c("exp","exp2"),
                       labels = c("Model I",
                                  "Model II")) +
  ylab("Observed relative shoot growth") +
  xlab("Expected relative shoot growth") +
  theme_bw() + theme(legend.title=element_blank()) 

## compare fit
lm <- lm(ShootRG ~ value*method, data = mixed.l)
summary(lm)

lma <- lm(ShootRG ~ exp, data = mixed)
summary(lma) ## Adjusted R-squared:  0.06295 

lmb <- lm(ShootRG ~ exp2, data = mixed)
summary(lmb) ## Adjusted R-squared:  0.6824 

## compare among coinoculation treatments for exp2 (by the reviewer)
ggplot(mixed, aes(x=exp2, y = ShootRG, color = type)) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y~x") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  ylab("Predicted relative shoot growth") +
  xlab("Observed relative shoot growth")


lmc <- lm(ShootRG ~ exp2*type, data = mixed)
summary(lmc) ## sig less than 1, but treatments are not sig different

## compare among approaches

## compare graphically
mixed.l <- mixed %>%
  select(ShootRG, exp, exp2, treatment) %>%
  gather(., method, value, -treatment)

p <- ggplot(mixed.l, aes(x=method, y = value)) +
  geom_boxplot() +
  geom_jitter() +
  scale_x_discrete(breaks = c("ShootRG","exp","exp2"),
                   labels = c("Obsesrved values",
                              "Current expected values",
                              "Alternative expected values")) +
  xlab("Comparisons") +
  ylab("Relative shoot growth") +
  theme_bw()


lm <- lm(sqrt(value) ~ method, data = mixed.l)
aov <- aov(lm)
summary(aov)
TukeyHSD(aov)


##########################################
# Summarize observed-expected difference #
##########################################
rgdiff <- a0$data[,c('x1', 'diff', 'percent_pvalue')]
noddiff <- a1$data[,c('x1', 'diff', 't_p_percent')]
nodareadiff <- a3$data[,c('x1', 'diff', 't_p_percent')]
nodweightdiff <- a4$data[,c('x1', 'diff', 't_p_percent')]

colnames(rgdiff) <- c('treatments', 'rg_diff', 'rg_pval')
colnames(noddiff) <- c('treatments', 'nod_diff', 'nod_pval')
colnames(nodareadiff) <- c('treatments', 'nodarea_diff', 'nodarea_pval')
colnames(nodweightdiff) <- c('treatments', 'nodweight_diff', 'nodweight_pval')

diff_data <- merge(rgdiff, noddiff, by = 'treatments')
diff_data <- merge(diff_data, nodareadiff, by = 'treatments')
diff_data <- merge(diff_data, nodweightdiff, by = 'treatments')
View(diff_data)


###################################
# Correlation between differences #
###################################
library(ggpmisc)

rg <- co.response.rg[,c('x1', 'diff')]
area <- co.response.area[,c('x1', 'diff')]
nod <- co.response.nod[,c('x1', 'diff')]

colnames(rg)[2] <- 'rg.diff'
colnames(area)[2] <- 'area.diff'
colnames(nod)[2] <- 'nod.diff'

diff <- merge(rg, area, by='x1')
diff <- merge(diff, nod, by='x1')

head(diff)

p1 <- ggplot(diff, aes(y=rg.diff, x=nod.diff)) + 
  geom_point() + geom_smooth(method='lm', formula=y~x) + 
  stat_poly_eq(formula=y~x, aes(label=paste(..eq.label.., ..rr.label.., sep='~~')), parse=T) +
  xlab('ΔNodule Number') + ylab('ΔRG') + ggtitle('b.')
  
p2 <- ggplot(diff, aes(x=rg.diff, y=area.diff)) + geom_point() +
  geom_smooth(method='lm', formula=y~x) + 
  stat_poly_eq(formula=y~x, aes(label=paste(..eq.label.., ..rr.label.., sep='~~')), parse=T) +
  xlab('ΔRG') + ylab('ΔNodule Area')

p3 <- ggplot(diff, aes(x=nod.diff, y=area.diff)) + geom_point() +
  geom_smooth(method='lm', formula=y~x) + 
  stat_poly_eq(formula=y~x, aes(label=paste(..eq.label.., ..rr.label.., sep='~~')), parse=T) +
  xlab('ΔNodule Number') + ylab('ΔNodule Area') + ggtitle('a.')


grid.arrange(p3, p1, nrow=1)

cor.test(diff$rg.diff, diff$nod.diff)
cor.test(diff$rg.diff, diff$area.diff)

# Test for autocorrelation
model <- lm(rg.diff~nod.diff, data=diff)
durbinWatsonTest(model)

model <- lm(rg.diff~area.diff, data=diff)
durbinWatsonTest(model)

# mod.area <- lm(rg.diff ~ area.diff, data=diff)
# summary(mod.nod)
# mod.nod <- lm(rg.diff ~ nod.diff, data=diff)
# summary(aov(mod.nod))

 ####
legend <- get_legend(p1)
p1 <- p1 + theme(legend.position = 'none')
p2 <- p2 + theme(legend.position = 'none')
p3 <- p3 + theme(legend.position = 'none')
p4 <- p4 + theme(legend.position = 'none')
blankplot <- ggplot()+geom_blank(aes(1,1)) +
  cowplot::theme_nothing()

grid.arrange(p1, p2, legend, p3, p4, blankplot, ncol=3, widths=c(2.3, 2.3, 0.8))


#########################
# Analyzing interaction #
#########################
co.response.rg$ratio <- co.response.rg$y/co.response.rg$predict
q <- ggplot(co.response.rg, aes(x=optimize, color=Fix, y=y)) + geom_boxplot() + ylab('Relative Growth')  + xlab('')
r <- ggplot(co.response.rg, aes(x=ratio, y=y, color=Fix)) + geom_point() + ylab('Relative Growth') + xlab('') + 
  geom_smooth(method='lm') + theme(legend.position = c(.1, .90))

co.response.nod$ratio <- co.response.rg$ratio
q1 <- ggplot(co.response.nod, aes(x=optimize, y=y, color=Fix)) + geom_boxplot() + ylab('Total Nodules')  + xlab('')
r1 <- ggplot(co.response.nod, aes(x=ratio, y=y, color=Fix)) + geom_point() +  ylab('Total Nodules')  + xlab('') + 
  geom_smooth(method='lm') + theme(legend.position = c(.1, .90)) + stat_regline_equation(label.x=1, label.y=10)

r1 + facet_grid()
co.response.mw$ratio <- co.response.rg$ratio
q2 <- ggplot(co.response.mw, aes(x=optimize, y=y, color=Fix)) + geom_boxplot() + ylab('Mean Nodule Weight') + xlab('')
r2 <- ggplot(co.response.mw, aes(x=ratio, y=y, color=Fix)) + geom_point() + ylab('Mean Nodule Weight') + xlab('') + 
  geom_smooth(method='lm') + theme(legend.position = c(.1, .90))

co.response.inv$ratio <- co.response.rg$ratio

q3 <- ggplot(co.response.inv, aes(x=optimize, y=y, color=Fix)) + geom_boxplot() + ylab('Investment') + xlab('')
r3 <- ggplot(co.response.inv, aes(x=ratio, y=y, color=Fix)) + geom_point() + ylab('Investment') + xlab('Intensity of Interaction (Observed RG / Expected RG) ') + 
  geom_smooth(method='lm') + theme(legend.position = c(.1, .90))

co.response.area$ratio <- co.response.rg$ratio

q4 <- ggplot(co.response.area, aes(x=optimize, y=y, color=Fix)) + geom_boxplot() +
  xlab('Interaction Type') + ylab('Nodule Area')
r4 <- ggplot(co.response.area, aes(x=ratio, y=y, color=Fix)) + 
  geom_point() +xlab('Intensity of Interaction (Observed RG / Expected RG) ') + ylab('Nodule Area') + 
  geom_smooth(method='lm') + theme(legend.position = c(.1, .90))

r
grid.arrange( r1, r2, r3, r4)
# 





#########################
# Visualizing Residuals #
#########################
head(co.response.rg)

co.response.rg$diff <- co.response.rg$predict - co.response.rg$y
hist(co.response.rg$diff)

p1 <- ggplot(co.response.rg, aes(x=x1, y=diff, fill=optimize)) + 
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=90)) +
  ylab('Expected RG - Observed RG') + xlab('Treatments')

co.response.rg$residuals <- sqrt((co.response.rg$predict - co.response.rg$y)^2/2)
hist(co.response.rg$residuals)

p2 <- ggplot(co.response.rg, aes(x=x1, y=residuals, fill=optimize)) + 
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle=90)) +
  ylab('Residuals') + xlab('Treatments')
  

grid.arrange(p1, p2, nrow=2)

###################################################
# Initializing co.response.rg dataset for network #
###################################################
co.response.rg$A <- NA
co.response.rg$B <- NA

for (i in co.response.rg$decouple) {
  A <- i[1]
  B <- i[2]
  co.response.rg[which(co.response.rg$x1 == paste(i[1], i[2], sep='+')), 'A'] <- A
  co.response.rg[which(co.response.rg$x1 == paste(i[1], i[2], sep='+')), 'B'] <- B
}

##########################################
# Calculating Relative Competition Index #
##########################################
single.nod.melt <- response.plot(single.nod, 'RG', 'treatment')$data
coinoc.nod.melt <- response.plot(coinoc.nod, 'RG', 'treatment')$data

single.nod.melt$RCI <- NA

for (i in single.nod.melt$x1){
  RG.coinoc.sum <- sum(co.response.rg[which(co.response.rg$A == i | co.response.rg$B == i), 'y'])
  n <- length(co.response.rg[which(co.response.rg$A == i | co.response.rg$B == i), 'y'])
  RG.single <- single.nod.melt[which(single.nod.melt$x1 == i), 'y']
  single.nod.melt[which(single.nod.melt$x1 == i), 'RCI'] <- (RG.single - (RG.coinoc.sum/n) ) / RG.single
}

rci <- ggplot(single.nod.melt, aes(x=factor(x1, levels=c(131, 184, 2, 156, 4, 186, 200, 187)), y= RCI)) + geom_bar(stat='identity') +
  xlab('Treatments')

rci
#######################
# Visualizing Network #
#######################

# Create nodes list
sources <- co.response.rg %>%
  distinct(A) %>%
  rename(label = A)

destination <- co.response.rg %>%
  distinct(B) %>%
  rename(label = B)

nodes <- full_join(sources, destination, by='label')
nodes <- nodes %>%
  mutate(id=1:nrow(nodes)) %>%
  select(id, everything())
nodes$size <- c(0.6390008, 0.3490073, 0.5136868, -0.2462751, -4.4573908, 0.4894246, 0.1188352, -2.7133815) + 5

# Create edges list
co.response.net <- co.response.rg %>%
  rename(weight=y, optimum_p=p, fix=Fix)

edges <- co.response.net %>%
  left_join(nodes, by=c('A'='label')) %>%
  rename(from=id)

edges <- edges %>%
  left_join(nodes, by=c('B'='label')) %>%
  rename(to=id)

edges <- select(edges, from, to, weight, optimum_p, fix)

# Create network/graph

g <- graph_from_data_frame(d=edges, vertices=nodes, directed=F)

is_weighted(g)

E(g)$color <- ifelse(E(g)$fix == 'Effective', 'blue', 'red2')

plot(g, layout=layout_with_graphopt, edge.width=E(g)$weight, vertex.size=5*nodes$size+2)





#####################
# Frequency in pair #
#####################
#p <- response.plot(coinoc.nod, 'ShootRG', 'treatment')

pairs <- data.frame(Treatment=NA, A=NA, B=NA)

for (i in co.response.rg$decouple) {
  A <- i[1]
  B <- i[2]
  Treatment <- paste(as.character(i[1]), as.character(i[2]), sep='+')
  if (A>B) {
    pairs[nrow(pairs)+1, ] <- c(Treatment, A, B)
  } else {
    pairs[nrow(pairs)+1, ] <- c(Treatment, B, A)
    }
  
  if (dim(pairs[which(pairs$A == B & pairs$B == A), ])[1] != 0) {
    pairs[nrow(pairs)+1, ] <- c(Treatment, B, A)
  } 
}

pairs <- merge(pairs[-1,], p0$data[, c('x1', 'p', 'Fix')], by.x='Treatment', by.y='x1')
pairs$opt <- 'Optimum'
pairs[which(pairs$p < 0.05), ]$opt <- 'Suboptimum'

mat <- ggplot(pairs, aes(as.character(A), as.character(B), fill=opt)) + geom_tile(color='gray') +
  scale_fill_manual(values=c('#4cb3b8', '#ef6358')) + xlab('')  + ylab('') + theme(legend.title = element_blank())

mat

ggplot(pairs, aes(as.character(A), as.character(B), fill=Fix)) + geom_tile(color='gray') +
  scale_fill_manual(values=c('#4cb3b8', '#ef6358'))

#unexp <- pairs[which(pairs$p <= 0.05), ]
#exp <- pairs[which(pairs$p > 0.05), ]

#table(c(unexp$A, unexp$B))
#table(c(exp$A, exp$B))

df1 <- data.frame(treatment=c(131, 156, 184, 186, 187, 2, 200, 4),
                  Suboptimum=c(2, 3, 2, 3, 0, 4, 1, 1),
                  Optimum=c(5, 4, 5, 4, 7, 3, 6, 6))

df1 <- melt(df1, id.vars='treatment', measure.vars=c('Suboptimum', 'Optimum'))

ct <- ggplot(df1, aes(x=as.factor(treatment), y=value, fill=variable)) + 
  geom_bar(stat='identity') + xlab('Treatment') + ylab('#') + ggtitle('Times Appeared in Co-inoculation') +
  theme(legend.position='none')

mat + annotation_custom(ggplotGrob(ct), xmin=as.character(156), xmax=as.character(187), ymin=as.character(187), ymax=as.character(200))

#grid.arrange(rci, mat)


#pairs <- merge(pairs, tradeoff[, c('treatment', 'optimize')], by.x='Treatment', by.y='treatment')

#sorted.pairs <- pairs %>% 
#  rowwise() %>% 
#  mutate(pair=sort(c(A, B)) %>% paste(collapse=",")) %>%
#  group_by(pair) %>%
#  distinct(pair, .keep_all=T)
  
ggplot(pairs, aes(A, B, fill=optimize)) + geom_tile() + 
  scale_fill_manual(values=c('#4cb3b8', '#ef6358'))





##############
## Doing PCA #
##############

# Harvest variables dataset
coinoc.pca.data <- coinoc.nod[,c('treatment', 'Single_effectivity', 'Co_effectivity', 'total_nodules', 'ShootRG', 'investment', 'AreaMM', 'dryMass')] %>%
  drop_na() 

coinoc.active_var <- coinoc.pca.data[, c('total_nodules', 'ShootRG', 'investment', 'AreaMM', 'dryMass', 'Single_effectivity', 'Co_effectivity')]

res.pca <- PCA(coinoc.active_var[, 1:5], graph=F)

# Scree plot
fviz_eig(res.pca, addlabels = T)

var <- get_pca_var(res.pca)
# plot variables
fviz_pca_var(res.pca, col.var='black')
# correlation plot
corrplot(var$cos2, is.corr=FALSE)

# Biplot fill by Single_effectivity
pca1 <- fviz_pca_biplot(res.pca,
                pointshape=21,
                pointsize=2.5,
                fill.ind = coinoc.active_var$Single_effectivity,
                palette =  c("#00AFBB", "#E7B800", "#FC4E07"),
                addEllipses = T, label='var',
                repel=T, 
                title="Harvest variables (pair-wise Fix trait)")

# Biplot fill by Single_effectivity
pca2 <- fviz_pca_biplot(res.pca,
                pointshape=21,
                pointsize=2.5,
                fill.ind = coinoc.active_var$Co_effectivity,
                palette =  c("#00AFBB", "#E7B800"),
                addEllipses = T, label='var',
                repel=T,
                title='A. Harvest variables (community Fix trait)')



#######################################
# PCA on observed dataset by optimity #
#######################################

# Observed dataset
merged <- merge(p0$data[, c('x1' ,'y', 'optimize')], p1$data[, c('x1', 'y')], by='x1')
colnames(merged) <- c('treatment', 'RG', 'Optimum', 'nodules')
merged <- merge(merged, p2$data[, c('x1', 'y')], by.x='treatment', by.y='x1')
colnames(merged) <- c('treatment', 'RG', 'Optimum', 'nodules', 'weight')
merged <- merge(merged, p3$data[, c('x1', 'y')], by.x='treatment', by.y='x1')
colnames(merged) <- c('treatment', 'RG', 'Optimum', 'nodules', 'weight', 'investment')
merged <- merge(merged, p4$data[, c('x1', 'y')], by.x='treatment', by.y='x1')
colnames(merged) <- c('treatment', 'RG', 'Optimum', 'nodules', 'weight', 'investment', 'area')

rownames(merged) <- merged$treatment

head(merged)

#PCA
optimity.pca <- PCA(merged[, c(2,4,5,6,7)], graph=F)

# Scree plot
fviz_eig(optimity.pca, addlabels = T)

var <- get_pca_var(optimity.pca)
# plot variables
fviz_pca_var(optimity.pca, col.var='black')
# correlation plot
corrplot(var$cos2, is.corr=FALSE)

# Biplot fill by optimum or not
pca3 <- fviz_pca_biplot(optimity.pca,
                pointshape=21,
                pointsize=2.5,
                fill.ind = merged$Optimum,
                palette =  c("#00AFBB", "#E7B800"),
                addEllipses = T, label='var',
                repel=T,
                title='B. Harvest variables (optimicity)')

grid.arrange(pca2, pca3, nrow=1)
#######################################
# PCA on residual dataset by optimity #
#######################################

# Diff = overved - expected
p0$data$Diff <- p0$data$y - p0$data$predict # RG
p1$data$Diff <- p1$data$y - p1$data$predict # Nodules
p2$data$Diff <- p2$data$y - p2$data$predict # Nod weight
p3$data$Diff <- p3$data$y - p3$data$predict # Investment
p4$data$Diff <- p4$data$y - p4$data$predict # nod Area

# Residuals = |obs - exp| / 2
p0$data$residuals <- sqrt((p0$data$predict - p0$data$y)^2/2)
p1$data$residuals <- sqrt((p1$data$predict - p1$data$y)^2/2)
p3$data$residuals <- sqrt((p3$data$predict - p3$data$y)^2/2)
p2$data$residuals <- sqrt((p2$data$predict - p2$data$y)^2/2)
p4$data$residuals <- sqrt((p4$data$predict - p4$data$y)^2/2)

# Residuals dataset
residuals <- merge(p0$data[, c('x1' ,'Diff', 'residuals', 'optimize', 'Fix')], p1$data[, c('x1', 'Diff', 'residuals')], by='x1')
colnames(residuals) <- c('treatment', 'RG.diff', 'RG.res', 'Optimum', 'Fix', 'nod.diff', 'nod.res')
residuals <- merge(residuals, p2$data[, c('x1', 'Diff', 'residuals')], by.x='treatment', by.y='x1')
colnames(residuals) <- c('treatment', 'RG.diff', 'RG.res', 'Optimum', 'Fix', 'nod.diff', 'nod.res', 'weight.diff', 'weight.res')
residuals <- merge(residuals, p3$data[, c('x1', 'Diff', 'residuals')], by.x='treatment', by.y='x1')
colnames(residuals) <- c('treatment', 'RG.diff', 'RG.res', 'Optimum', 'Fix', 'nod.diff', 'nod.res', 'weight.diff', 'weight.res',
                         'inv.diff', 'inv.res')
residuals <- merge(residuals, p4$data[, c('x1', 'Diff', 'residuals')], by.x='treatment', by.y='x1')
colnames(residuals) <- c('treatment', 'RG.diff', 'RG.res', 'Optimum', 'Fix', 'nod.diff', 'nod.res', 'weight.diff', 'weight.res',
                         'inv.diff', 'inv.res', 'area.diff', 'area.res')

head(residuals)

#PCA
residuals.pca <- PCA(residuals[, c(2,3,6,7, 8, 9, 10, 11, 12, 13)], graph=F)

# Scree plot
fviz_eig(residuals.pca, addlabels = T)

var <- get_pca_var(residuals.pca)
# plot variables
fviz_pca_var(residuals.pca, col.var='black')
# correlation plot
corrplot(var$cos2, is.corr=FALSE)

# Biplot fill by optimum or not
pca4 <- fviz_pca_biplot(residuals.pca,
                pointshape=21,
                pointsize=2.5,
                fill.ind = residuals$Optimum,
                palette =  c("#00AFBB", "#E7B800"),
                addEllipses = T, 
                label='var',
                repel=T, 
                title= 'A. Residuals variables (optimicity) ')

# Biplot fill by Fix
pca5 <- fviz_pca_biplot(residuals.pca,
                pointshape=21,
                pointsize=2.5,
                fill.ind = residuals$Fix,
                palette =  c("#00AFBB", "#E7B800"),
                addEllipses = T, label='var',
                repel=T,
                title='B. Residuals variables (community Fix trait)'
              )


grid.arrange(pca4, pca5, nrow=1)

#################################
# MDS on obsserved vs expeceted #
#################################
merged <-  merge(p1$data[, c('x1', 'p')], p2$data[, c('x1', 'p')], by='x1' )
merged <- merge(merged, p3$data[, c('x1', 'p')], by='x1')
merged <- merge(merged, p4$data[, c('x1', 'p')], by='x1')
colnames(merged) <- c('Treatment', 'RG', 'Nodulation', 'Nodule_Weight', 'Investment')
rownames(merged) <- merged$Treatment

mds <- merged[, c(2,3,5)] %>%
  dist() %>%
  cmdscale() %>%
  as_tibble()
colnames(mds) <- c('Dim.1', 'Dim.2')

#k-mean clustering
clust <- kmeans(mds, 4)$cluster %>%
  as.factor()

mds <- mds %>%
  mutate(group=clust)

ggscatter(mds, x='Dim.1', y='Dim.2', label=rownames(merged), 
          size=1, repel=T, color='group', pallette='jco', ellipse=T)



#########
## GLM ##
#########


glm.rg <- glm(RG~ treatment + elapsed_dpi, data=harvest.single.ordered, family=poisson(link='log'))
plot(residuals(glm.rg))

res <- aov(mod.rg)
summary(res)

#############################
# Data visualization by Fix #
#############################
nod.mean <- aggregate(harvest[,c("total_nodules","mean_weight", "RG")], by=list(Fix=harvest$Fix), FUN=central.mean)
colnames(nod.mean) <-  c('Fix', 'Total Nodules', 'Nodule Mean Weight (mg)', 'Relative Growth')
nod.mean <- melt(nod.mean)

nod.se <- aggregate(harvest[,c("total_nodules","mean_weight", "RG")], by=list(Fix=harvest$Fix), FUN=se)
colnames(nod.se) <- c('Fix', 'Total Nodules', 'Nodule Mean Weight (mg)', 'Relative Growth')
nod.se <- melt(nod.se)

nod.stat <- merge(nod.mean, nod.se, by = c("Fix", "variable")) # may be start using cbind(nod.mean, nod.se)
colnames(nod.stat) <- c('Fix', 'variable','mean','se')

ggplot(nod.stat[which(nod.stat$Fix != 'Control'),], aes(fill=Fix, x=Fix, y=mean)) + 
  geom_bar(stat='identity', position=position_dodge(width=0.9)) + 
  geom_errorbar(aes( ymin=mean-se, ymax=mean+se), position = position_dodge(width = 0.90), width=0.2) +
  facet_wrap(.~variable, scales='free')
  

###############################
# ANOVA, TukeyHSD test by Fix #
###############################

# New FixTuk variable with letter assignment

# ANOVA models

# Assessing ANOVA models

# TukeyHSD tests

# Function for extracting Tukey letters

# Merging Tukey labels

# Boxplot with Tukey labels




####################
# Comparison graph #
####################



comparison.graph <- function(combination, ab_line){
  if (missing(ab_line)){
    ab_line<- 0
  }
  
  coinoc_comparison <- aggregate( cbind(total_nodules, RG)~treatment, data=combination, FUN=mean, na.rm=TRUE)
  nodule_se <- aggregate( combination$total_nodules, by=list(treatment=combination$treatment), FUN=se)
  colnames(nodule_se)[2] <- 'se' 
  RG_se <- aggregate( combination$RG, by=list(treatment=combination$treatment), FUN=se)
  colnames(RG_se)[2] <-'se'
  coinoc_comparison$nodule_se <- nodule_se$se
  coinoc_comparison$RG_se <- RG_se$se
  #head(coinoc_comparison)
  
  # Distribution
  nodules_dist <- ggplot(data=combination, aes(y=total_nodules, x=treatment)) + 
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2), aes(color=batch))
  
  biomass_dist <- ggplot(data=combination, aes(y=RG, x=treatment)) + 
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2), aes(color=batch))
  
  # Calculated means
  #nodules_mean <- ggplot(data=coinoc_comparison, aes(y=total_nodules, x=treatment)) + 
  #  geom_bar(stat = 'summary', fun.y=central.mean) +
  #  geom_errorbar(aes(ymax=total_nodules+nodule_se, ymin=total_nodules-nodule_se), width=.1, position=position_dodge(0.9))
  ggplot(data=coinoc_comparison, aes(y=RG, x=treatment)) + 
    geom_bar(stat = 'summary', fun.y=central.median) +
    geom_hline(yintercept=ab_line) +
    geom_errorbar(aes(ymax=RG+RG_se, ymin=RG-RG_se), width=.1, position=position_dodge(0.9))
  
  #plot_grid(nodules_dist, biomass_dist, nodules_mean, biomass_mean, labels='AUTO')
  #grid.arrange(nodules_dist, biomass_dist, nodules_mean, biomass_mean, ncol=2, nrow=2)
  #grid.arrange( biomass_dist,  biomass_mean, ncol=2, nrow=1)
  #biomass_mean
}
#################
# Finalized Viz #
#################

# I need photos of combination plots


coinoc <- c('131+156', '131+184', '131+187', '131+200', '156+200', '156+184', '184+200', 
            '186+131', '186+156', '186+184', '186+200', '186+4', '187+156', '187+184', 
            '187+200', '2+131', '2+156', '2+184', '2+187', '2+200', '2+4', '4+131', '4+156', 
            '4+184', '4+187', '4+200', '186+187', '2+186')

png(filename='comparison_graph.png')

for (treatment in coinoc){
  a = unlist(strsplit(treatment, '\\+'))
  combination <- harvest[which(harvest$treatment==a[1] | harvest$treatment==a[2] | harvest$treatment==treatment),]
  comparison.graph(combination)
}

dev.off()
##################################
# Initial Harvest Visualizations #
##################################

# Calculate SE for error bar
se <- function(x) sqrt(sd(x, na.rm=TRUE)/sqrt(length(x)))
central.median <- function(x) median(x, na.rm=TRUE)
central.mean <- function(x) mean(x, na.rm=TRUE)

# Single Inoc
##############

# Box plot
harvest.single.ordered <- harvest.single[order(harvest.single$Fix), ]
harvest.single.ordered$treatment <- factor(harvest.single.ordered$treatment, levels=unique(harvest.single.ordered$treatment))

ggplot(data=harvest.single.ordered, aes(x=treatment, y=shoot_mass, fill=Fix)) + 
  #geom_boxplot() +
  geom_bar(stat='summary', fun.y=central.mean)

ggplot(data=harvest.single.ordered, aes(x=treatment, y=root_mass, fill=Fix)) + 
  geom_boxplot() 

ggplot(data=harvest.single.ordered, aes(x=treatment, y=total_nodules, fill=Fix)) + 
  geom_boxplot() +
  geom_jitter()

ggplot(data=harvest.single.ordered, aes(x=treatment, y=AreaMM, fill=Fix)) + 
  geom_boxplot() +
  geom_jitter()

ggplot(data=harvest.single.ordered, aes(x=treatment, y=mean_weight, fill=Fix)) + 
  geom_boxplot() +
  geom_jitter()

ggplot(data=harvest.single.ordered, aes(x=treatment, y=RG, fill=Fix)) + 
  geom_boxplot() +
  geom_jitter()

ggplot(data=harvest.single.ordered, aes(x=treatment, y=SEff, fill=Fix)) + 
  geom_boxplot() +
  geom_jitter()


# Shoot/Root mass bar plot
shoot_se <- aggregate(harvest.single$shoot_mass, by=list(treatment=harvest.single$treatment), FUN=se)
colnames(shoot_se)[2] <- 'shoot_se'
shoot_mass <- aggregate(harvest.single$shoot_mass, by=list(treatment=harvest.single$treatment), FUN=central.mean)
colnames(shoot_mass)[2] <- 'shoot_mass'
root_se <- aggregate(harvest.single$root_mass, by=list(treatment=harvest.single$treatment), FUN=se)
colnames(root_se)[2] <- 'root_se'
root_mass <- aggregate(harvest.single$root_mass, by=list(treatment=harvest.single$treatment), FUN=central.mean)
colnames(root_mass)[2] <- 'root_mass'

growth_bar <- Reduce(function(x, y) merge(x, y, all=T), list(shoot_mass, shoot_se, root_mass, root_se))
growth_bar <- merge(growth_bar, traits.single, by.x='treatment', by.y='Treatment')

ggplot(growth_bar, aes(x=treatment, y=shoot_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=shoot_mass-shoot_se, ymax=shoot_mass+shoot_se)) +
  ggtitle('Shoot mass statistics (mean)') +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(growth_bar, aes(x=treatment, y=root_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=root_mass-root_se, ymax=root_mass+root_se)) +
  ggtitle('Root mass statistics (mean)') +
  theme(axis.text.x = element_text(angle=80, hjust=1))

# Nodules bar plot
nodule_count <-  aggregate(harvest.single$total_nodules, by=list(treatment=harvest.single$treatment), FUN=central.mean)
colnames(nodule_count)[2] <- 'nodule_count'
nodule_se <- aggregate(harvest.single$total_nodules, by=list(treatment=harvest.single$treatment), FUN=se)
colnames(nodule_se)[2] <- 'nodule_se'
nodule_fresh_mass <- aggregate(harvest.single$mean_weight, by=list(treatment=harvest.single$treatment), FUN=central.mean)
colnames(nodule_fresh_mass)[2] <- 'nodule_fresh_mass'
nodule_mass_se <- aggregate(harvest.single$mean_weight, by=list(treatment=harvest.single$treatment), FUN=se)
colnames(nodule_mass_se)[2] <- 'nodule_mass_se'

nodule_bar <- Reduce(function(x, y) merge(x, y, all=T), list(nodule_count, nodule_se, nodule_fresh_mass, nodule_mass_se))
nodule_bar <- merge(nodule_bar, traits.single, by.x='treatment', by.y='Treatment')

ggplot(nodule_bar, aes(x=treatment, y=nodule_count, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=nodule_count-nodule_se, ymax=nodule_count+nodule_se)) +
  ggtitle('Total nodule count statistics (mean)') +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(nodule_bar, aes(x=treatment, y=nodule_fresh_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=nodule_fresh_mass-nodule_mass_se, ymax=nodule_fresh_mass+nodule_mass_se)) +
  ggtitle('Fresh nodule biomass statistics (mean)') + 
  theme(axis.text.x = element_text(angle=80, hjust=1))

# Co-inoc
##############

# Box plot
harvest.coinoc.ordered <- harvest.coinoc[order(harvest.coinoc$Fix), ]
harvest.coinoc.ordered$treatment <- factor(harvest.coinoc.ordered$treatment, levels=unique(harvest.coinoc.ordered$treatment))

ggplot(data=harvest.coinoc.ordered, aes(x=treatment, y=shoot_mass, fill=Fix)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(data=harvest.coinoc.ordered, aes(x=treatment, y=root_mass, fill=Fix)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(data=harvest.coinoc.ordered, aes(x=treatment, y=total_nodules, fill=Fix)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(data=harvest.coinoc.ordered, aes(x=treatment, y=AreaMM, fill=Fix)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(data=harvest.coinoc.ordered, aes(x=treatment, y=RG, fill=Fix)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(data=harvest.coinoc.ordered, aes(x=treatment, y=SEff, fill=Fix)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle=80, hjust=1))


# Shoot mass
shoot_se <- aggregate(harvest.coinoc$shoot_mass, by=list(treatment=harvest.coinoc$treatment), FUN=se)
colnames(shoot_se)[2] <- 'shoot_se'
shoot_mass <- aggregate(harvest.coinoc$shoot_mass, by=list(treatment=harvest.coinoc$treatment), FUN=central.median)
colnames(shoot_mass)[2] <- 'shoot_mass'
root_se <- aggregate(harvest.coinoc$root_mass, by=list(treatment=harvest.coinoc$treatment), FUN=se)
colnames(root_se)[2] <- 'root_se'
root_mass <- aggregate(harvest.coinoc$root_mass, by=list(treatment=harvest.coinoc$treatment), FUN=central.median)
colnames(root_mass)[2] <- 'root_mass'

growth_bar <- Reduce(function(x, y) merge(x, y, all=T), list(shoot_mass, shoot_se, root_mass, root_se))
growth_bar <- merge(growth_bar, traits.coinoc, by.x='treatment', by.y='Treatment')

ggplot(growth_bar, aes(x=treatment, y=shoot_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=shoot_mass-shoot_se, ymax=shoot_mass+shoot_se)) +
  ggtitle('Shoot mass statistics (median)') +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(growth_bar, aes(x=treatment, y=root_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=root_mass-root_se, ymax=root_mass+root_se)) +
  ggtitle('Root mass statistics (median)') +
  theme(axis.text.x = element_text(angle=80, hjust=1))

# Nodules
nodule_count <-  aggregate(harvest.coinoc$total_nodules, by=list(treatment=harvest.coinoc$treatment), FUN=central.median)
colnames(nodule_count)[2] <- 'nodule_count'
nodule_se <- aggregate(harvest.coinoc$total_nodules, by=list(treatment=harvest.coinoc$treatment), FUN=se)
colnames(nodule_se)[2] <- 'nodule_se'
nodule_fresh_mass <- aggregate(harvest.coinoc$mean_weight, by=list(treatment=harvest.coinoc$treatment), FUN=central.median)
colnames(nodule_fresh_mass)[2] <- 'nodule_fresh_mass'
nodule_mass_se <- aggregate(harvest.coinoc$mean_weight, by=list(treatment=harvest.coinoc$treatment), FUN=se)
colnames(nodule_mass_se)[2] <- 'nodule_mass_se'

nodule_bar <- Reduce(function(x, y) merge(x, y, all=T), list(nodule_count, nodule_se, nodule_fresh_mass, nodule_mass_se))
nodule_bar <- merge(nodule_bar, traits.coinoc, by.x='treatment', by.y='Treatment')


ggplot(nodule_bar, aes(x=treatment, y=nodule_count, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=nodule_count-nodule_se, ymax=nodule_count+nodule_se)) +
  ggtitle('Total nodule count statistics (median)') +
  theme(axis.text.x = element_text(angle=80, hjust=1))

ggplot(nodule_bar, aes(x=treatment, y=nodule_fresh_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=nodule_fresh_mass-nodule_mass_se, ymax=nodule_fresh_mass+nodule_mass_se)) +
  ggtitle('Fresh nodule biomass statistics (median)') + 
  theme(axis.text.x = element_text(angle=80, hjust=1))





###############################
# Visualizing Relative Growth #
###############################

# Single Inoc
#############
RG_se <- aggregate(harvest.single$RG, by=list(treatment=harvest.single$treatment), FUN=se)
colnames(RG_se)[2] <- 'RG_se'
RG_mass <- aggregate(harvest.single$RG, by=list(treatment=harvest.single$treatment), FUN=central.median)
colnames(RG_mass)[2] <- 'RG_mass'

growth_bar <- Reduce(function(x, y) merge(x, y, all=T), list(RG_mass, RG_se, root_mass, root_se))
growth_bar <- merge(growth_bar, traits.single, by.x='treatment', by.y='Treatment')

ggplot(growth_bar, aes(x=treatment, y=RG_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RG_mass-RG_se, ymax=RG_mass+RG_se)) +
  ggtitle('RG statistics (median)')

##########
# Coinoc #
##########
DM_I_minus <- mean(harvest.coinoc[which(harvest.coinoc$treatment=='H2O'),]$shoot_mass + 
                     harvest.coinoc[which(harvest.coinoc$treatment=='H2O'),]$root_mass, na.rm = TRUE)
DM_I_plus <- harvest.coinoc$shoot_mass + harvest.coinoc$root_mass
harvest.coinoc$RG <- DM_I_plus / DM_I_minus

# Visualize relative growth
RG_se <- aggregate(harvest.coinoc$RG, by=list(treatment=harvest.coinoc$treatment), FUN=se)
colnames(RG_se)[2] <- 'RG_se'
RG_mass <- aggregate(harvest.coinoc$RG, by=list(treatment=harvest.coinoc$treatment), FUN=central.median)
colnames(RG_mass)[2] <- 'RG_mass'

growth_bar <- Reduce(function(x, y) merge(x, y, all=T), list(RG_mass, RG_se, root_mass, root_se))
growth_bar <- merge(growth_bar, traits.coinoc, by.x='treatment', by.y='Treatment')

ggplot(growth_bar, aes(x=treatment, y=RG_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RG_mass-RG_se, ymax=RG_mass+RG_se)) +
  ggtitle('RG statistics (median)') +
  theme(axis.text.x = element_text(angle=80, hjust=1)) +
  geom_hline(yintercept=0.4415794, color='red')


####################################
# Calculating symbiotic efficiency #
####################################


# Single Inoc
#############
DM_I_minus <- mean(harvest.single[which(harvest.single$treatment=='H2O'),]$shoot_mass + 
                     harvest.single[which(harvest.single$treatment=='H2O'),]$root_mass, na.rm = TRUE)
DM_I_plus <- harvest.single$shoot_mass + harvest.single$root_mass + harvest.single$total_nodules*harvest.single$mean_weight
harvest.single$SEff <- (DM_I_plus - DM_I_minus) 

# Visualize relative growth
SE_se <- aggregate(harvest.single$SEff, by=list(treatment=harvest.single$treatment), FUN=se)
colnames(SE_se)[2] <- 'SE_se'
SE_mass <- aggregate(harvest.single$SEff, by=list(treatment=harvest.single$treatment), FUN=central.mean)
colnames(SE_mass)[2] <- 'SE_mass'

growth_bar <- Reduce(function(x, y) merge(x, y, all=T), list(SE_mass, SE_se, root_mass, root_se))
growth_bar <- merge(growth_bar, traits.single, by.x='treatment', by.y='Treatment')

ggplot(growth_bar, aes(x=treatment, y=SE_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=SE_mass-SE_se, ymax=SE_mass+SE_se)) +
  ggtitle('Symbiotic efficiency statistics (mean)')

# Coinoc
##############
DM_I_minus <- mean(harvest.coinoc[which(harvest.coinoc$treatment=='H2O'),]$shoot_mass + 
                     harvest.coinoc[which(harvest.coinoc$treatment=='H2O'),]$root_mass, na.rm = TRUE)
DM_I_plus <- harvest.coinoc$shoot_mass + harvest.coinoc$root_mass + harvest.coinoc$total_nodules*harvest.coinoc$mean_weight
harvest.coinoc$SEff <- (DM_I_plus - DM_I_minus)

# Visualize relative growth
SE_se <- aggregate(harvest.coinoc$SEff, by=list(treatment=harvest.coinoc$treatment), FUN=se)
colnames(SE_se)[2] <- 'SE_se'
SE_mass <- aggregate(harvest.coinoc$SEff, by=list(treatment=harvest.coinoc$treatment), FUN=central.mean)
colnames(SE_mass)[2] <- 'SE_mass'

growth_bar <- Reduce(function(x, y) merge(x, y, all=T), list(SE_mass, SE_se, root_mass, root_se))
growth_bar <- merge(growth_bar, traits.coinoc, by.x='treatment', by.y='Treatment')

ggplot(growth_bar, aes(x=treatment, y=SE_mass, fill=Fix)) + 
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=SE_mass-SE_se, ymax=SE_mass+SE_se)) +
  ggtitle('Symbiotic efficiency statistics (mean)') + 
  theme(axis.text.x = element_text(angle=80, hjust=1))


ggplot(harvest.coinoc, aes(x=treatment, y=SEff, fill=Fix.x)) + 
  geom_boxplot() + theme(axis.text.x = element_text(angle=80, hjust=1))
###############################
# Calculated projected growth #
###############################


###############################
# Post-Hoc test vizualization #
###############################
library(plyr)
library(multcompView)


generate_label_df <- function(HSD, flev){
  # Extract labels and factor levels from Tukey post-hoc 
  Tukey.levels <- HSD[[flev]][,4]
  Tukey.labels <- multcompLetters(Tukey.levels)['Letters']
  plot.labels <- names(Tukey.labels[['Letters']])
  
  # Get highest quantile for Tukey's 5 number summary and add a bit of space to buffer between    
  # upper quantile and label placement
  boxplot.df <- ddply(d, flev, function (x) max(fivenum(x$y)) + 0.2)
  
  # Create a data frame out of the factor levels and Tukey's homogenous group letters
  plot.levels <- data.frame(plot.labels, labels = Tukey.labels[['Letters']],
                            stringsAsFactors = FALSE)
  
  # Merge it with the labels
  labels.df <- merge(plot.levels, boxplot.df, by.x = 'plot.labels', by.y = flev, sort = FALSE)
  
  return(labels.df)
}


# tHSD for RG
mod.coinc <-  aov(log(RG)~treatment, data=harvest.coinoc.ordered)
tHSD<- TukeyHSD(mod.coinc)

d <- data.frame(lev=harvest.coinoc.ordered$treatment, y=log(harvest.coinoc.ordered$RG))
treatment <- unique(harvest.coinoc.ordered$treatment)

ggplot(harvest.coinoc.ordered, aes(x=treatment, y=RG)) + 
  geom_boxplot(data=harvest.coinoc.ordered, aes(fill=Fix)) +
  geom_text(data = generate_label_df(tHSD, 'treatment'), aes(x = plot.labels, y = 10, label = labels)) + 
  theme(axis.text.x = element_text(angle=80, hjust=1))


# tHSD for SEff
mod.coincSEff <-  aov(SEff~treatment, data=harvest.coinoc.ordered)
tHSD.SEff <- TukeyHSD(mod.coincSEff)

d <- data.frame(lev=harvest.coinoc.ordered$treatment, y=harvest.coinoc.ordered$SEff)
treatment <- unique(harvest.coinoc.ordered$treatment)

ggplot(harvest.coinoc.ordered, aes(x=treatment, y=SEff)) + 
  geom_boxplot(data=harvest.coinoc.ordered, aes(fill=Fix)) +
  geom_text(data = generate_label_df(tHSD.SEff, 'treatment'), aes(x = plot.labels, y = 200, label = labels)) + 
  theme(axis.text.x = element_text(angle=80, hjust=1))

##############
# Using glht #
##############

library(multcomp)

mod <-  aov(log(RGperDayCorrected)~treatment, data=harvest.coinoc)
tuk <- glht(mod, linfct=mcp(treatment='Tukey'))
tuk.cld <- cld(tuk)
plot(tuk.cld)

mod.coinc <-  aov(log(RGperDay)~treatment, data=harvest.coinoc)
tuk.coinc <- glht(mod.coinc, linfct=mcp(treatment='Tukey'))
tuk.coinc.cld <- cld(tuk.coinc)
plot(tuk.coinc.cld)


################
# Innoculation #
################
innoc.count <- read.csv('../data/inoculation_count.csv', header = TRUE)


#########################
# Comparison test Graph #
#########################
se <- function(x) sqrt(sd(x, na.rm=TRUE)/sqrt(length(x)))

comparison.graph <- function(combination, ab_line){
  if (missing(ab_line)){
    ab_line<- 0
  }
  
  coinoc_comparison <- aggregate( cbind(total_nodules, RG)~treatment, data=combination, FUN=mean, na.rm=TRUE)
  nodule_se <- aggregate( combination$total_nodules, by=list(treatment=combination$treatment), FUN=se)
  colnames(nodule_se)[2] <- 'se' 
  RG_se <- aggregate( combination$RG, by=list(treatment=combination$treatment), FUN=se)
  colnames(RG_se)[2] <-'se'
  coinoc_comparison$nodule_se <- nodule_se$se
  coinoc_comparison$RG_se <- RG_se$se
  head(coinoc_comparison)
  
  # Distribution
  nodules_dist <- ggplot(data=combination, aes(y=total_nodules, x=treatment)) + 
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2), aes(color=batch))
  
  biomass_dist <- ggplot(data=combination, aes(y=RG, x=treatment)) + 
    geom_boxplot() + 
    geom_jitter(shape=16, position=position_jitter(0.2), aes(color=batch))
  
  # Calculated means
  #nodules_mean <- ggplot(data=coinoc_comparison, aes(y=total_nodules, x=treatment)) + 
  #  geom_bar(stat = 'summary', fun.y=central.mean) +
  #  geom_errorbar(aes(ymax=total_nodules+nodule_se, ymin=total_nodules-nodule_se), width=.1, position=position_dodge(0.9))
  biomass_mean <- ggplot(data=coinoc_comparison, aes(y=RG, x=treatment)) + 
    geom_bar(stat = 'summary', fun.y=central.median) +
    geom_hline(yintercept=ab_line) +
    geom_errorbar(aes(ymax=RG+RG_se, ymin=RG-RG_se), width=.1, position=position_dodge(0.9))
  
  #plot_grid(nodules_dist, biomass_dist, nodules_mean, biomass_mean, labels='AUTO')
  #grid.arrange(nodules_dist, biomass_dist, nodules_mean, biomass_mean, ncol=2, nrow=2)
  grid.arrange( biomass_dist,  biomass_mean, ncol=2, nrow=1)
}

#segmenting harvest dataset by batchs
harvest.a <- harvest[which(harvest$batch==1), ] 
harvest.b <- harvest[which(harvest$batch==2), ] 

#combination <- harvest[which(harvest$treatment==2 | 
#                            harvest$treatment==187 | 
#                            harvest$treatment=='2+187' | 
#                            harvest$treatment=='H2O'),]#

#comparison.graph(combination)
#mu = (central.mean(combination[which(combination$treatment==2),]$shoot_mass) +
#central.mean(combination[which(combination$treatment==187),]$shoot_mass))/2
#test <-na.omit(combination[which(combination$treatment=='2+187'),]$shoot_mass)
#wilcox.test(test, mu=mu)

coinoc <- c('131+156', '131+184', '131+187', '131+200', '156+200', '156+184', '184+200', 
            '186+131', '186+156', '186+184', '186+200', '186+4', '187+156', '187+184', 
            '187+200', '2+131', '2+156', '2+184', '2+187', '2+200', '2+4', '4+131', '4+156', 
            '4+184', '4+187', '4+200', '186+187', '2+186')

table <- NULL

for (treatment in coinoc){
  a = unlist(strsplit(treatment, '\\+'))
  #print(treatment)
  combination <- harvest[which(harvest$treatment==a[1] | harvest$treatment==a[2] | harvest$treatment==treatment | harvest$treatment=='H2O'),]
  combination.a <- harvest.a[which(harvest.a$treatment==a[1] | harvest.a$treatment==a[2] | harvest.a$treatment==treatment | harvest.a$treatment=='H2O'),]
  combination.b <- harvest.b[which(harvest.b$treatment==a[1] | harvest.b$treatment==a[2] | harvest.b$treatment==treatment | harvest.b$treatment=='H2O'),]
  #comparison.graph(combination)
  
  
  ## test 50%
  mu <- (central.median(combination[which(combination$treatment==a[1]),]$RG) +
          central.median(combination[which(combination$treatment==a[2]),]$RG))/2
  test <-na.omit(combination[which(combination$treatment==treatment),]$RG)
  stat <- wilcox.test(test, mu=mu)
  
  # Proportion by batch
  total <- innoc.count[innoc.count$Inocula == a[1],] + innoc.count[innoc.count$Inocula == a[2],]
  proportion.strA <- innoc.count[innoc.count$Inocula == a[1],] / total
  proportion.strB <- innoc.count[innoc.count$Inocula == a[2],] / total
  
  # Proportion by treatment
  strA.total <- sum(innoc.count[innoc.count$Inocula == a[1],][2:3])
  strB.total <- sum(innoc.count[innoc.count$Inocula == a[2],][2:3])
  
  proportion.strA.total <- strA.total / (strA.total + strB.total)
  proportion.strB.total <- strB.total / (strA.total + strB.total)
  
  
  ## test by batch
  mu.a <- median(combination.a[which(combination.a$treatment==a[1]),]$RG, na.rm=T) * proportion.strA[2] +
    + median(combination.a[which(combination.a$treatment==a[2]),]$RG, na.rm=T) * proportion.strB[2]
  #mu.a <- (central.median(combination.a[which(combination.a$treatment==a[1]),]$RG) +
  #         central.median(combination.a[which(combination.a$treatment==a[2]),]$RG))/2
  test.a <-na.omit(combination.a[which(combination.a$treatment==treatment),]$RG)
  stat.a <- wilcox.test(test.a, mu=mu.a[[1]])
  
  
  mu.b <- median(combination.b[which(combination.b$treatment==a[1]),]$RG, na.rm=T) * proportion.strA[3] +
    + median(combination.b[which(combination.b$treatment==a[2]),]$RG, na.rm=T) * proportion.strB[3]
  #mu.b <- (central.median(combination.b[which(combination.b$treatment==a[1]),]$RG) +
  #           central.median(combination.b[which(combination.b$treatment==a[2]),]$RG))/2
  test.b <-na.omit(combination.b[which(combination.b$treatment==treatment),]$RG)
  stat.b <- wilcox.test(test.b, mu=mu.b[[1]])
  
  
  ## test all
  mu.ab <- median(combination[which(combination$treatment==a[1]),]$RG, na.rm=T) * proportion.strA.total +
    + median(combination[which(combination$treatment==a[2]),]$RG, na.rm=T) * proportion.strB.total
  test.ab <- na.omit(combination[which(combination$treatment==treatment),]$RG)
  stat.ab <- wilcox.test(test.ab, mu=mu.ab[[1]])
  
  #print(c(treatment, stat$p.value, stat.a$p.value, stat.b$p.value))
  table <- rbind(table, data.frame(treatment, mu.ab[[1]], stat.ab$p.value, mu, stat$p.value, mu.a[[1]], stat.a$p.value, mu.b[[1]], stat.b$p.value))
  #table[nrow(table)+1,] <- list(treatment, stat$p.value, stat.a$p.value, stat.b$p.value)
}




#######################################################
# Testing if the batches has significant differences  #
#######################################################
all_treatments <- unique(data$treatment)
report_t <- data.frame(matrix(ncol=4, nrow=0), stringsAsFactors = FALSE)
colnames(report_t) <- c('treatment','root_t','shoot_t','nodule_t')

for (i in all_treatments) {
  shoot_t <- t.test(data[which(data$batch=='1' & data$treatment==i),]$shoot_mass, 
                    data[which(data$batch=='2' & data$treatment==i),]$shoot_mass, 
                    var.equal = T, paired=F)
  root_t <- t.test(data[which(data$batch=='1' & data$treatment==i),]$root_mass, 
                   data[which(data$batch=='2' & data$treatment==i),]$root_mass, 
                   var.equal = T, paired=F)
  nodule_t <- t.test(data[which(data$batch=='1' & data$treatment==i),]$total_nodule, 
                     data[which(data$batch=='2' & data$treatment==i),]$total_nodule, 
                     var.equal = T, paired=F)
  report_t <- rbind(report_t, data.frame('treatment'=i, 'root_t'=shoot_t$p.value, 'shoot_t'=root_t$p.value, 'nodule_t'=nodule_t$p.value))
  #report_t[nrow(report_t)+1,] <= list(i, shoot_t$p.value, root_t$p.value, nodule_t$p.value)
  #print(i, shoot_t$p.value, root_t$p.value, nodule_t$p.value)
}

# Printing treatments that have significant differences between batch
report_t[which(report_t$root_t<=.05 | report_t$shoot_t <=0.05 | report_t$nodule_t <= 0.05), ]

p1 <- ggplot(data=data[which(data$treatment=='131+184'),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('131+184') 
p2 <- ggplot(data=data[which(data$treatment=='131+187'),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('131+187') 
p3 <- ggplot(data=data[which(data$treatment=='131+200'),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('131+200') 
p4 <- ggplot(data=data[which(data$treatment=='184+200'),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('184+200') 
p5 <- ggplot(data=data[which(data$treatment=='186+187'),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('186+187')
p6 <- ggplot(data=data[which(data$treatment=='2+200'),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('2+200') 
p7 <- ggplot(data=data[which(data$treatment==200),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('200') 
p8 <- ggplot(data=data[which(data$treatment==156),], aes(y=shoot_mass, x=elapsed_dpi, color=batch)) + geom_point() + ggtitle('156') 

grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8, nrow=4, ncol=2)

###################################################
# Testing if there are significant effect of time #
###################################################


#####################
# Correcting counts #
#####################

ggplot(data=data, aes(x=elapsed_dpi, y=shoot_mass)) + geom_point() + geom_smooth(method='lm')
fit <- lm(shoot_mass~elapsed_dpi+total_nodules, data=data)
summary(fit)


fit <- glm(shoot_mass~elapsed_dpi, data=data[which(data$treatment=='H2O'),], family=Gamma())
summary(fit)
ggplot(data=data[which(data$treatment=='H2O'),], aes(x=elapsed_dpi, y=shoot_mass, color=batch)) + 
  geom_point() + 
  geom_smooth(method='lm')


########################
# Regression with Stat #
########################

# Function that produce ggplot with regression statistics
ggplotRegression <- function (fit, i) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       "Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)
                       ),
         x= paste(i, 'elapsed_dpi'),
         y='RG')
}

# Append plots 
plot_list = list()

for (i in coinoc){
  p <- ggplotRegression(lm(RG~elapsed_dpi, data=harvest.coinoc[which(harvest.coinoc$treatment==i), ]), i)
  plot_list[[i]] = p
}

# Save plots in pdf
pdf("RG with time regression plots.pdf")
for (i in coinoc){
  print(plot_list[[i]])
}
dev.off()


# data frame to extract statistics from each regression model
dpi_regression.df <- data.frame(treatment=as.character(), RSquared=as.numeric(), 
                                Intercept=as.numeric(), Slope=as.numeric(), P=as.numeric())


for (i in coinoc){
       fit <- lm(RG~elapsed_dpi, data=harvest.coinoc[which(harvest.coinoc$treatment==i), ])
       newrow <- data.frame(treatment=i, 
                   RSquared=signif(summary(fit)$adj.r.squared, 5), 
                   Intercept=signif(fit$coef[[1]],5 ), 
                   Slope=signif(fit$coef[[2]],5 ), 
                   P=summary(fit)$coef[2,4])
       dpi_regression.df <- rbind(dpi_regression.df, newrow)
       }
# merge trait data with regression stat
dpi_regression.df <- merge(dpi_regression.df, traits.coinoc, by.x='treatment', by.y='Treatment')
dpi_regression.df[which(dpi_regression.df$P<= 0.05), ]

# Estimate corrected RGperDay
calculateRGperDay <- function(i) {
  #print(i)
  if (harvest.coinoc[i,]$treatment=='H2O'){
    return(harvest.coinoc[i,]$RG/45)
  }
  if (dpi_regression.df[harvest.coinoc[i,]$treatment, ]$P <= 0.05){
    return(harvest.coinoc[i,]$RG/(1.5*as.numeric(harvest.coinoc[i,]$elapsed_dpi)))
    #return(harvest.coinoc[i,]$RG/(1.2*as.numeric(harvest.coinoc[i,]$elapsed_dpi)))
  }else{
    return(harvest.coinoc[i,]$RG/45)
  }
}

harvest.coinoc$RGperDayCorrected <- sapply(1:nrow(harvest.coinoc), calculateRGperDay)


