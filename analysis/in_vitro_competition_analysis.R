# In vitro competition analysis

library('ggplot2')
library('reshape2')

#Exp 1
data <- read.csv('~/GDrive/Sachs/Chapter1_Competition/in_vitro_test/in_vitro_competition_klett.csv')


colnames(data) <- c('Strain', 'Replicate', 'Treatment', 19, 34, 42, 58.75, 69.76, 82.5, 93.33, 107, 119.5, 134.32,	157, 190.5)
data.series <- melt(data, id=c('Strain', 'Replicate', 'Treatment'))


# Growth Curve by individual replicats
ggplot(data.series[which(data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=value, group=Replicate)) + 
  geom_point() + geom_line() +
  facet_wrap(~Strain) + ggtitle('A. Clonal')


ggplot(data.series[which(!data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=value, group=Replicate)) + 
  geom_point() + geom_line() +
  facet_wrap(~Strain) + ggtitle('B. Competition')


# Growth Curve by individual replicats 
data.series$variable <- as.numeric(as.character(data.series$variable))


ggplot(data.series[which(data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=log10(value), group=Replicate)) + 
  geom_point() + geom_line() +
  facet_wrap(~Strain) + ggtitle('A. Clonal')


ggplot(data.series[which(!data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=log10(value), group=Replicate)) + 
  geom_point() + geom_line() +
  facet_wrap(~Strain) + ggtitle('B. Competition')



# Convert Klett numbers to CFU counts
# But it gives negative CFU counts for lower Klett values, which does not make any sense
#a = 4576468.3
#b = 46322867
#(data[,4:15]*a)-b
#(data[,4:15]*5.8*10^6)+2.96*10^7

# Growth curve measurement
library('growthcurver')

gc.data <- data[, 3:dim(data)[2]]
gc.data.t <- t(gc.data)
colnames(gc.data.t) <- gc.data.t[1,]
gc.data.t <- cbind(time=rownames(gc.data.t), gc.data.t)
gc.data.t <- gc.data.t[2:dim(gc.data.t)[1],]
gc.data.t <- as.data.frame(apply(gc.data.t, 2, as.numeric))

head(gc.data.t)

gc_fit <- SummarizeGrowth(gc.data.t$time, gc.data.t$`156+184.5`, bg_correct = 'none' )
plot(gc_fit)
gc_fit

gc_out <- SummarizeGrowthByPlate(gc.data.t, bg_correct='none')
gc_out




# Plotting growth stat
library(ggpubr)

gw_stat <- read.csv('~/GDrive/Sachs/Chapter1_Competition/in_vitro_test/growth_stat.csv')

head(gw_stat)

# Plot Carrying capacity
ggplot(gw_stat, aes(x=Strain, y=K)) + geom_boxplot() + ggtitle('Carrying Capacity') +
  facet_wrap(~Experiment, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggbarplot(gw_stat, x='Strain', y='K', label=F, add='mean_se', fill = 'Experiment', scales='free')

# Plot DT
ggplot(gw_stat, aes(x=Strain, y=DT)) + geom_boxplot() + ggtitle('Doubling Time') +
  facet_wrap(~Experiment, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggbarplot(gw_stat, x='Strain', y='DT', label=F, add='mean_se', fill = 'Experiment', scales='free')

# Plot r
ggplot(gw_stat, aes(x=Strain, y=r)) + geom_boxplot() + ggtitle('Growth Rate') +
  facet_wrap(~Experiment, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggbarplot(gw_stat, x='Strain', y='r', label=F, add='mean_se', fill = 'Experiment', scales='free')

# Median Statistics
ggplot(data.series[which(data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=value, group=Strain)) + 
  stat_summary(fun=mean, geom="line") +
  facet_wrap(~Strain) + ggtitle('A. Clonal')

ggplot(data.series[which(!data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=value, group=Strain)) + 
  stat_summary(fun=mean, geom="line") +
  facet_wrap(~Strain) + ggtitle('B. Competition')


# Exp2

data <- read.csv('~/GDrive/Sachs/Chapter1_Competition/in_vitro_test/in_vitro_competition_klett_2nd_exp.csv')

colnames(data) <- c('Strain', 'Replicate', 'Treatment', 1,	5.5,	8,	11,	14,	17.5,	19.5,	23,	28,	31,	33,	36,	41,	46,	54,	63.5,	68.5,	78.5,	85,	90,	101.5,	109,	114.5,	125.5,	150)
data.series <- melt(data, id=c('Strain', 'Replicate', 'Treatment'))

data.series$variable <- as.numeric(as.character(data.series$variable))
# Growth Curve by individual replicats
ggplot(data.series[which(data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=value, group=Replicate)) + 
  geom_point() + geom_line() + 
  scale_y_continuous(limits = c(0, 250)) +
  facet_wrap(~Strain) + ggtitle('A. Clonal')


ggplot(data.series[which(!data.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=value, group=Replicate)) + 
  geom_point() + geom_line() +
  scale_y_continuous(limits = c(0, 250)) +
  facet_wrap(~Strain) + ggtitle('B. Competition')

# Convert Klett numbers to CFU counts
# Use two equation for lower and higher Klett
# But it gives negative CFU counts for lower Klett values, which does not make any sense

cfu <- data[4:dim(data)[2]]


# If klett value is more than 50:
a = 4576468.3
b = 46322867
cfu <- cfu*a - b

# If CFU is negative, put 0
cfu[cfu < 0] <- 6.66E6


cfu <- cbind(data[1:3], cfu)



cfu.series <- melt(cfu, id=c('Strain', 'Replicate', 'Treatment'))
cfu.series$variable <- as.numeric(as.character(cfu.series$variable))

# Growth Curve by individual replicats for CFU counts
ggplot(cfu.series[which(cfu.series$Strain %in% c(4, 131, 156, 184)),], aes(x=variable, y=log10(value), group=Replicate)) + 
  geom_point() + geom_smooth() + 
  facet_wrap(~Strain) + ggtitle('A. Clonal') + xlab('Hours Post Inoculation') + ylab('log10(cells/mL)')


ggplot(cfu.series[which(!cfu.series$Strain %in% c(4, 131, 156, 184, 'blank')),], aes(x=variable, y=log10(value), group=Replicate)) + 
  geom_point() + geom_smooth() + 
  facet_wrap(~Strain) + ggtitle('B. Competition') + xlab('Hours Post Inoculation') + ylab('log10(cells/mL)')


# Growth curve measurement
library('growthcurver')

gc.data <- cfu[, 3:dim(cfu)[2]]
gc.data.t <- t(gc.data)
colnames(gc.data.t) <- gc.data.t[1,]
gc.data.t <- cbind(time=rownames(gc.data.t), gc.data.t)
gc.data.t <- gc.data.t[2:dim(gc.data.t)[1],]
gc.data.t <- as.data.frame(apply(gc.data.t, 2, as.numeric))

head(gc.data.t)

gc_fit <- SummarizeGrowth(gc.data.t$time, gc.data.t$`4.2`, bg_correct = 'none')
plot(gc_fit)
gc_fit

gc_out <- SummarizeGrowthByPlate(gc.data.t, bg_correct='none')
gc_out



# Growth curve measurement for klett values
gc.data <- data[, 3:dim(data)[2]]
gc.data.t <- t(gc.data)
colnames(gc.data.t) <- gc.data.t[1,]
gc.data.t <- cbind(time=rownames(gc.data.t), gc.data.t)
gc.data.t <- gc.data.t[2:dim(gc.data.t)[1],]
gc.data.t <- as.data.frame(apply(gc.data.t, 2, as.numeric))

head(gc.data.t)

gc_fit <- SummarizeGrowth(gc.data.t$time, gc.data.t$`4.5` )
plot(gc_fit)
gc_fit

# Plotting 24hr data
cfu.24 <- cfu[5:15]
# Exponential model fitting
#test = 
y = log2(as.numeric(cfu.24[4,4:dim(cfu.24)[2]]))
t = as.numeric(colnames(cfu.24)[4: dim(cfu.24)[2]])

test = data.frame(y, t)
plot(t, y)
mod <- nls(y~ exp(a + b*t), data=test, start=list(a=0, b=0))
mod <- lm(y~t, test)

lines(test$t, predict(mod, list(x=test$t)))

1/coef(mod)[[2]]

########################
# Plotting growth stat #
########################

library(ggpubr)
library(dplyr)

gw_stat <- read.csv('~/GDrive/Sachs/Chapter1_Competition/in_vitro_test/growth_stat.csv')

# Summary stat
# Mean and SE of K, r, DT by each strain(s)
gw_summary <- gw_stat
gw_summary$K <- log10(gw_summary$K)

gw_summary.mean <- gw_summary %>% group_by(Strain) %>%
  summarise(across(c(K,r,DT), mean))

gw_summary.se <- gw_summary %>% group_by(Strain) %>%
  summarise(across(c(K,r,DT), se))

colnames(gw_summary.se) <- c('Strain', 'K.se', 'r.se', 'DT.se')

gw_summary <- merge(gw_summary.mean, gw_summary.se, by='Strain')

View(gw_summary)

# Plot Carrying capacity
ggplot(gw_stat, aes(x=Strain, y=K)) + geom_boxplot() + ggtitle('Carrying Capacity') +
  facet_wrap(~Experiment, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggbarplot(gw_stat, x='Strain', y='K', label=F, add='mean_se', fill = 'Experiment', scales='free')

# Plot DT
ggplot(gw_stat, aes(x=Strain, y=DT)) + geom_boxplot() + ggtitle('Doubling Time') +
  facet_wrap(~Experiment, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggbarplot(gw_stat, x='Strain', y='DT', label=F, add='mean_se', fill = 'Experiment', scales='free')

# Plot r
ggplot(gw_stat, aes(x=Strain, y=r)) + geom_boxplot() + ggtitle('Growth Rate') +
  facet_wrap(~Experiment, scales = "free") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggbarplot(gw_stat, x='Strain', y='r', label=F, add='mean_se', fill = 'Experiment', scales='free')


################
# ANOVA Clonal #
################

# Carrying capacity
res.aov <- aov(log10(K)~Strain, data=gw_stat[which(gw_stat$Experiment == 'Clonal'),] )

summary(res.aov)
View(TukeyHSD(res.aov)$Strain)

# Growth Rate
res.aov <- aov(r~Strain, data=gw_stat[which(gw_stat$Experiment == 'Clonal'),] )

summary(res.aov)
View(TukeyHSD(res.aov)$Strain)


# Doubling Time
res.aov <- aov(DT~Strain, data=gw_stat[which(gw_stat$Experiment == 'Clonal'),])
summary(res.aov)
View(TukeyHSD(res.aov)$Strain)

#####################
# ANOVA Competition #
#####################

# Carrying capacity
res.aov <- aov(log10(K)~Strain, data=gw_stat[which(gw_stat$Experiment != 'Clonal'),] )

summary(res.aov)
View(TukeyHSD(res.aov)$Strain)

# Growth Rate
res.aov <- aov(r~Strain, data=gw_stat[which(gw_stat$Experiment != 'Clonal'),] )

summary(res.aov)
TukeyHSD(res.aov)


# Doubling Time
res.aov <- aov(DT~Strain, data=gw_stat[which(gw_stat$Experiment != 'Clonal'),])
summary(res.aov)
TukeyHSD(res.aov)

##################################
# Growth Stat Hypothesis Testing #
##################################
gw_stat <- read.csv('~/GDrive/Sachs/Chapter1_Competition/in_vitro_test/growth_stat.csv')

gw_stat$K <- log10(gw_stat$K)
# Carrying capacity, K
# Comparing expectation based on clonal vs observation in competition
t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 131)), ]$K),
       gw_stat[which(gw_stat$Strain == '4+131'),]$K)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 156)), ]$K),
       gw_stat[which(gw_stat$Strain == '4+156'),]$K)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 184)), ]$K),
       gw_stat[which(gw_stat$Strain == '4+184'),]$K)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(131, 156)), ]$K),
       gw_stat[which(gw_stat$Strain == '131+156'),]$K)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(131, 184)), ]$K),
       gw_stat[which(gw_stat$Strain == '131+184'),]$K)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(156, 184)), ]$K),
       gw_stat[which(gw_stat$Strain == '156+184'),]$K)

# Doubling Time, DT
# Comparing expectation based on clonal vs observation in competition
t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 131)), ]$DT),
       gw_stat[which(gw_stat$Strain == '4+131'),]$DT)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 156)), ]$DT),
       gw_stat[which(gw_stat$Strain == '4+156'),]$DT)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 184)), ]$DT),
       gw_stat[which(gw_stat$Strain == '4+184'),]$DT)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(131, 156)), ]$DT),
       gw_stat[which(gw_stat$Strain == '131+156'),]$DT)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(131, 184)), ]$DT),
       gw_stat[which(gw_stat$Strain == '131+184'),]$DT)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(156, 184)), ]$DT),
       gw_stat[which(gw_stat$Strain == '156+184'),]$DT)

# Growth Rate, r
# Comparing expectation based on clonal vs observation in competition
t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 131)), ]$r),
       gw_stat[which(gw_stat$Strain == '4+131'),]$r)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 156)), ]$r),
       gw_stat[which(gw_stat$Strain == '4+156'),]$r)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(4, 184)), ]$r),
       gw_stat[which(gw_stat$Strain == '4+184'),]$r)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(131, 156)), ]$r),
       gw_stat[which(gw_stat$Strain == '131+156'),]$r)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(131, 184)), ]$r),
       gw_stat[which(gw_stat$Strain == '131+184'),]$r)

t.test(mu=mean(gw_stat[which(gw_stat$Strain %in% c(156, 184)), ]$r),
       gw_stat[which(gw_stat$Strain == '156+184'),]$r)

##############
# Solid RDM  #
##############
solid_exp <- data.frame(isolate=c(4, 131, 156, 184, '156+184', '4+156', '4+184', '131+184', '4+131', '131+156'),
           population = c(2.25E8, 2.83E9, 1.77E9, 3.22E9, 2.88E9, 2.95E9, 3.42E9, 3.05E9, 3.36E9, 2.63E9),
           test=c('clonal','clonal','clonal','clonal','competition','competition','competition','competition','competition','competition'))

ggbarplot(solid_exp, x='isolate', y='population', fill='test', yscale='log10')

ggplot(solid_exp, aes(x=isolate, y=log10(population))) + 
  geom_bar(stat='identity') + 
  #scale_y_continuous(limits = c(5, 10)) +
  facet_wrap(~test, scales = 'free_x') 
  
##################
# Max's Equation #
##################
#

# If klett value is between 1-50:
cfu <- data[4:dim(data)[2]]

a = 5850000
b = 29600000
cfu[cfu < 50] <- cfu[cfu < 50]*a + b


cfu <- cbind(data[1:3], cfu)

head(cfu)
