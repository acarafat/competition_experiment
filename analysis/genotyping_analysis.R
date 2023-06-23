library(ggplot2)
library(ggrepel)
library(reshape2)
library(tidyr)

genotype = '../data/updated_count_Nov16.csv'
genotype = read.csv(genotype)

# Some nodules were labeled incorrectly. Lets update them (i.e. 184+2 > 2+184, 131+186 > 186+131)
genotype[which(genotype$id_a == '184' & genotype$id_b == '2'),] <- c(833, '351_184_2___6_1_S339', 2, 184, 0, 4062)

genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),]
new_id_a <- genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),]$id_b
new_id_b <- genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),]$id_a
new_count_a <- genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),]$count_b
new_count_b <- genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),]$count_a

genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),c('id_a', 'id_b', 'count_a', 'count_b')] <- data.frame(new_id_a, new_id_b, new_count_a, new_count_b)

# check update, should be empty
genotype[which(genotype$id_a == '184' & genotype$id_b == '2'),]
genotype[which(genotype$id_a == '131' & genotype$id_b == '186'),]


genotype[which(genotype$id_a == '186' & genotype$id_b =='131'),]

makeSwap <- function(genotype, a, b){
  new_id_a <- genotype[which(genotype$id_a == a & genotype$id_b == b),]$id_b
  new_id_b <- genotype[which(genotype$id_a == a & genotype$id_b == b),]$id_a
  new_count_a <- genotype[which(genotype$id_a == a & genotype$id_b == b),]$count_b
  new_count_b <- genotype[which(genotype$id_a == a & genotype$id_b == b),]$count_a
  upadteGenotypeSubset = data.frame(new_id_a, new_id_b, new_count_a, new_count_b)
  #subset <- genotype[which(genotype$id_a == a & genotype$id_b == b),c('id_a', 'id_b', 'count_a', 'count_b')]
  #genotype[which(genotype$id_a == a & genotype$id_b == b),c('id_a', 'id_b', 'count_a', 'count_b')] <- data.frame(new_id_a, new_id_b, new_count_a, new_count_b)\
  return(upadteGenotypeSubset)
}

# Flip all the treatments
genotype[which(genotype$id_a == 186 & genotype$id_b == 131),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '186', '131')
genotype[which(genotype$id_a == 186 & genotype$id_b == 156),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '186', '156')
genotype[which(genotype$id_a == 186 & genotype$id_b == 184),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '186', '184')
genotype[which(genotype$id_a == 186 & genotype$id_b == 131),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '186', '131')
genotype[which(genotype$id_a == 186 & genotype$id_b == 4),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '186', '4')
genotype[which(genotype$id_a == 187 & genotype$id_b == 156),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '187', '156')
genotype[which(genotype$id_a == 187 & genotype$id_b == 184),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '187', '184')
genotype[which(genotype$id_a == 2 & genotype$id_b == 131),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '2', '131')
genotype[which(genotype$id_a == 2 & genotype$id_b == 156),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '2', '156')
genotype[which(genotype$id_a == 2 & genotype$id_b == 184),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '2', '184')
genotype[which(genotype$id_a == 2 & genotype$id_b == 4),c('id_a', 'id_b', 'count_a', 'count_b')] <- makeSwap(genotype, '2', '4')


# Create treatment column
genotype$treatment = paste(genotype$id_a, genotype$id_b, sep='+')

#genotype = genotype[which(genotype$count_a != 0 | genotype$count_b != 0), ]
dim(genotype)
head(genotype)

genotype$class <- 'NA'

cut = 10

genotype[which(genotype$count_a >= cut & genotype$count_b >= cut),]$class = 'A+B'
genotype[which(genotype$count_a >= cut & genotype$count_b < cut),]$class = 'A'
genotype[which(genotype$count_a < cut & genotype$count_b >= cut),]$class = 'B'
genotype$count = 1

dim(genotype[which(genotype$class != 'NA'), ])



genoSummary <-  aggregate(count~class+treatment, data=genotype, FUN=sum)

genoSummary <- genoSummary %>% separate(treatment, c('A', 'B'))
genoSummary$treatment <- paste(genoSummary$A, genoSummary$B, sep='+')


both_fix <- c('131+184', '156+184', '4+131', '4+156', '4+184')
no_fix <- c('186+187', '186+200', '187+200', '2+186', '2+187', '2+200')
  
genoSummary$fix_combination <- '+/-'
genoSummary[which(genoSummary$treatment %in% both_fix),]$fix_combination <- '+/+'
genoSummary[which(genoSummary$treatment %in% no_fix),]$fix_combination <- '-/-'
  
# Flip some of the treatments 
genoSummary[which(genoSummary$A == '186' & genoSummary$B == '131'),]

# Final plot with facet_wrap
# Main result figure
ggplot(genoSummary[which(genoSummary$class != 'NA'),], aes(fill=class, x=treatment, y=count)) + 
  geom_bar(position="fill", stat="identity", show.legend = F) +
  scale_fill_manual(values=c('#B55A30', '#F7E0D4', '#0072B5'), name = "Nodule Infected By", labels = c("Strain A", "Coinfected", "Strain B")) +
  ylab('Nodule Occupancy') + xlab('Treatments') +
  geom_text(aes(label=A, y=1.04), color='#B55A30', angle=45) + 
  geom_text(aes(label=B, y=-0.04), color='#0072B5', angle=45) + 
  #scale_y_continuous(labels = scales::percent) +
  scale_y_continuous(labels=c('0%/100%', '25%/75%', '50%/50%', '75%/25%', '100%/0%'), breaks=c(0,0.25,0.50,0.75,1)) +
  facet_grid(~factor(fix_combination, levels=c('+/+', '+/-', '-/-')), scales="free_x", space='free_x', ) +
  #scale_fill_discrete(name = "", labels = c("Strain A", "Coinfected", "Strain B"), values=c('#DD4444', 'grey', '#555599')) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        strip.text.x = element_text(size = 14)) +
  ggtitle('A') + 
  theme(plot.title = element_text(face = "bold")) 
  

#theme(axis.text.x = element_blank()) 


abundance = genoSummary[which(genoSummary$class != 'NA'),]
abundance = dcast(abundance, treatment~class)

#write.csv(abundance, file= '/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition Experiment/AmpliconSeq/Amplicon_Fragments/competition_abundanc.csv', row.names=F)


##################################
# Analysis of coinfected nodules #
##################################
genoSummary <- genoSummary[which(genoSummary$class != 'NA'),]
genoSummary <- rbind(genoSummary, c('A+B', 131, 200, 0, '131+200', '+/-'))
genoSummary <- rbind(genoSummary, c('A+B', 4, 200, 0, '4+200', '+/-'))

total_nodules <- aggregate(as.numeric(genoSummary$count),  list(genoSummary$treatment), FUN=sum)
colnames(total_nodules) <- c('treatment', 'total_nodules')


#coinfected <- merge(coinfected, total_nodules, by='treatment')
coinfected <- merge(genoSummary, total_nodules, by='treatment')
coinfected$percent_coinfected <- (as.numeric(coinfected$count)/as.numeric(coinfected$total_nodules))*100
coinfected <- coinfected[which(coinfected$class == 'A+B'),]



ggplot(coinfected, aes(x=reorder(treatment, -percent_coinfected), y=percent_coinfected, fill=fix_combination)) + 
  geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 45, vjust = 0.1, hjust=0.1)) +
  xlab("Treatment") + ylab("% Coinfected nodules")

# Categorizing and comparing by ++, +-, and --
t_test_comparison <- compare_means(percent_coinfected~fix_combination, data=coinfected)
t_test_comparison <- t_test_comparison %>% mutate(y.position=c(90, 100, 110))

aggregate(percent_coinfected~fix_combination, data=coinfected, FUN=mean)
aggregate(percent_coinfected~fix_combination, data=coinfected, FUN=se)

ggbarplot(coinfected, x='fix_combination', y='percent_coinfected', fill='fix_combination', add='mean_se') +
  stat_pvalue_manual(t_test_comparison, label="p.adj") +
  scale_x_discrete(labels=c('Fix+/Fix+', 'Fix+/Fix-', 'Fix-/Fix-')) +
  theme(plot.title = element_text(face = "bold"), legend.position='bottom', legend.title=element_blank()) +
  xlab('') + ylab('% Coinfected nodules') + ggtitle('B') +
  scale_fill_manual(values=c('#2E8B57', '#619CFF', '#F8766D')) 



# Focusing only on strain 156

# 200 can not coinfect against 4 ...

####################################  
# Compare observed RG and expected #
####################################
abundance = read.csv( '../data/competition_abundanc.csv', header=T)
mean_RG_coinoc = read.csv( '../data/mean_RG_calculation.csv', header=T)
mean_RG_single = read.csv( '../data/mean_RG_single.csv', header=T)



# Estimate competence
library('tidyverse')
library('gridExtra')
library("cowplot")
library('ggpubr')
require(grid)

genoSummary.AB <- genoSummary[which(genoSummary$class != 'NA'),]
genoSummary.AB <- genoSummary.AB %>% separate(treatment, c('A', 'B'))

# Function to subset all combination of isolates with a key one.

genotype_subset <- function(base){
  subset <- genoSummary.AB[which(genoSummary.AB$A == base | genoSummary.AB$B == base),]
  subset$strain_A <- base
  subset$strain_B <- NA
  subset[which(subset$A != base),]$strain_B <- subset[which(subset$A != base),]$A
  subset[which(subset$B != base),]$strain_B <- subset[which(subset$B != base),]$B
  
  # adjust A/B order label for A+B in original treatment 
  if (dim(subset[which(subset$A != base & subset$class == 'B'),])[1] != 0) {
    subset[which(subset$A != base & subset$class == 'B'),]$class <- 'x'
  }
  if (dim(subset[which(subset$A != base & subset$class == 'A'),])[1] != 0){
    subset[which(subset$A != base & subset$class == 'A'),]$class <- 'y'
  }
  if (dim(subset[which(subset$A != base & subset$class == 'x'),])[1] != 0){
    subset[which(subset$A != base & subset$class == 'x'),]$class <- 'A'
  }
  if (dim(subset[which(subset$A != base & subset$class == 'y'),])[1] != 0){
    subset[which(subset$A != base & subset$class == 'y'),]$class <- 'B'
  }
  
  # get %occupance
  #agg.set <- aggregate(count~strain_B, data=subset, FUN = sum)
  
  return(subset)
}


plot_genotype_subset <- function(subset, base='', color_vector=c('#DD4444', 'grey', '#555599')){
  ggplot(subset, aes(fill=class, y=count, x=strain_B)) + 
    geom_bar(position="fill", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(values=color_vector) + 
    xlab('') + ylab('') + ggtitle(paste('Strain', base, sep=' ')) +
    theme(legend.position='none', legend.title = element_blank(), plot.title=element_text(color='#DD4444'), axis.text.x =element_text(color="#555599"))
}


# Whole panel
subset.2  <- genotype_subset('2')
subset.4 <- genotype_subset('4')
subset.131 <- genotype_subset('131')
subset.156 <- genotype_subset('156')
subset.184 <- genotype_subset('184')
subset.186 <- genotype_subset('186')
subset.187 <- genotype_subset('187')
subset.200 <- genotype_subset('200')

p2 <- plot_genotype_subset(subset.2, '2')
p4 <- plot_genotype_subset(subset.4, '4')
p131 <- plot_genotype_subset(subset.131, '131')
p156 <- plot_genotype_subset(subset.156, '156')
p184 <- plot_genotype_subset(subset.184, '184')
p186 <- plot_genotype_subset(subset.186, '186')
p187 <- plot_genotype_subset(subset.187, '187')
p200 <- plot_genotype_subset(subset.200, '200', c('grey', '#555599'))



genotype.summary <- ggarrange(p4, p131, p184, p156,  p2,p186, p187, p200, nrow=2, ncol=4, 
          common.legend = T, legend='top')

annotate_figure(genotype.summary, left=textGrob('Occupancy', rot=90), 
                fig.lab='A. Genotype based nodule occupancy')

relative_occupancy <- function(subset){
  group1 <- aggregate(count~strain_B, data=subset, FUN = sum)
  group2 <- aggregate(count~strain_B+class, data=subset, FUN = sum)
  group <- merge(group1, group2, by='strain_B', all.y = T)
  group <- group[which(group$class != 'B'), ]
  # A+B occupancy should be half
  group[which(group$class == 'A+B'), ]$count.y <- group[which(group$class == 'A+B'), ]$count.y/2
  #aggregate(count.y~strain_B, data=test, FUN=sum)
  group$percent <- group$count.y / group$count.x
  group = aggregate(percent ~ strain_B, data=group, FUN=sum)
  
  return(list(mean(group$percent), se(group$percent)))
}

# Mean Relative Growth
se <- function(x) sd(x, na.rm=TRUE)/sqrt(length(x))

get_mean_var <- function(base, var){
  mean_var = mean(rg[which(rg$strain_A == base | rg$strain_B == base), var], na.rm=T)
  se_var = se(rg[which(rg$strain_A == base | rg$strain_B == base), var])
  return(list(mean_var, se_var))
}

RO.value <- c(relative_occupancy(subset.131)[[1]],
              relative_occupancy(subset.156)[[1]],
        relative_occupancy(subset.4)[[1]],
        relative_occupancy(subset.184)[[1]],
        relative_occupancy(subset.186)[[1]],
        relative_occupancy(subset.187)[[1]],
        relative_occupancy(subset.2)[[1]],
        relative_occupancy(subset.200)[[1]])

RO.se <- c(relative_occupancy(subset.131)[[2]],
           relative_occupancy(subset.156)[[2]],
              relative_occupancy(subset.4)[[2]],
              relative_occupancy(subset.184)[[2]],
              relative_occupancy(subset.186)[[2]],
              relative_occupancy(subset.187)[[2]],
              relative_occupancy(subset.2)[[2]],
              relative_occupancy(subset.200)[[2]])

rg <- read.csv('../data/mean_RG_calculation.csv')

RG.value <- c(get_mean_var('131', 'RG')[[1]],
              get_mean_var('156', 'RG')[[1]],
              get_mean_var('4', 'RG')[[1]],
              get_mean_var('184', 'RG')[[1]],
              get_mean_var('186', 'RG')[[1]],
              get_mean_var('187', 'RG')[[1]],
              get_mean_var('2', 'RG')[[1]],
              get_mean_var('200', 'RG')[[1]])

RG.se <- c(get_mean_var('131', 'RG')[[2]],
           get_mean_var('156', 'RG')[[2]],
              get_mean_var('4', 'RG')[[2]],
              get_mean_var('184', 'RG')[[2]],
              get_mean_var('186', 'RG')[[2]],
              get_mean_var('187', 'RG')[[2]],
              get_mean_var('2', 'RG')[[2]],
              get_mean_var('200', 'RG')[[2]])

nodule.value <- c(get_mean_var('131', 'total_nodules')[[1]],
                  get_mean_var('156', 'total_nodules')[[1]],
              get_mean_var('4', 'total_nodules')[[1]],
              get_mean_var('184', 'total_nodules')[[1]],
              get_mean_var('186', 'total_nodules')[[1]],
              get_mean_var('187', 'total_nodules')[[1]],
              get_mean_var('2', 'total_nodules')[[1]],
              get_mean_var('200', 'total_nodules')[[1]])

nodule.se <- c(get_mean_var('131', 'total_nodules')[[2]],
               get_mean_var('156', 'total_nodules')[[2]],
                  get_mean_var('4', 'total_nodules')[[2]],
                  get_mean_var('184', 'total_nodules')[[2]],
                  get_mean_var('186', 'total_nodules')[[2]],
                  get_mean_var('187', 'total_nodules')[[2]],
                  get_mean_var('2', 'total_nodules')[[2]],
                  get_mean_var('200', 'total_nodules')[[2]])

RO.names <- c('131', '156', '4', '184',  '186', '187', '2', '200')
fix <- c('+', '+', '+', '+', '-', '-', '-', '-')
RO <- data.frame(RO.value, RO.se, RO.names, fix)
RO$RO.names <- factor(RO$RO.names, levels=RO$RO.names)
RO$RG.value <- RG.value
RO$RG.se <- RG.se
RO$nodule.value <- nodule.value
RO$nodule.se <- nodule.se
RO$dominance <- c(6, 6, 5, 4,  3, 2, 1, 0)

RO[which(RO$fix == '+'),]$fix <- "Fix+"
RO[which(RO$fix == '-'),]$fix <- "Fix-"
RO$fix <- factor(RO$fix, levels = c('Fix+', 'Fix-'))


Anova(lm(RO.value~fix, data=RO), type='II')
compare_means(RO.value~fix, data=RO)

# Main result figure


b = ggplot(RO, aes(x=RO.names, y=dominance, fill=fix)) + geom_bar(stat='identity') +
  xlab('Strain') + ylab('#Times dominant in coinoculation') + ggtitle('C') + 
  theme(legend.position='none') +
  theme_bw() + scale_fill_manual(values=c('#2E8B57', '#F8766D')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
theme(plot.title = element_text(face = "bold")) 


c = ggplot(RO, aes(x=RO.names, y=RO.value, fill=fix)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RO.value - RO.se, ymax=RO.value + RO.se), width=0.2) +
  xlab('Strain') + ylab('Nodule Occupancy Mean Count Value') + ggtitle('D') + 
  theme(legend.position='none') +
  #annotate('text', x='186', y=1, label='Wilcoxon p-value: 0.0286', col='blue') +
  theme_bw() + scale_fill_manual(values=c('#2E8B57', '#F8766D')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 

d = ggplot(RO, aes(x=RO.names, y=nodule.value, fill=fix)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=nodule.value - nodule.se, ymax=nodule.value + nodule.se), width=0.2) +
  xlab('Strain') + ylab('Nodule Mean Count Value') + ggtitle('E') + 
  theme(legend.position='none') +
  theme_bw() + scale_fill_manual(values=c('#2E8B57', '#F8766D')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 

e = ggplot(RO, aes(x=RO.names, y=RG.value, fill=fix)) + geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=RG.value - RG.se, ymax=RG.value + RG.se), width=0.2) +
  xlab('Strain') + ylab('Shoot RG Mean Count Value') + ggtitle('F') + theme(legend.position='none') +
  theme_bw() + scale_fill_manual(values=c('#2E8B57', '#F8766D')) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), legend.title = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(plot.title = element_text(face = "bold")) 



ggarrange(b, c, d, e, ncol=4, common.legend = T, legend='bottom')

t.test(e$data[which(e$data$RO.names %in% c('131', '145', '4', '184')), 'RG.value'], 
       e$data[which(e$data$RO.names %in% c('2', '186', '187', '200')), 'RG.value'])

t.test(e$data[which(e$data$RO.names %in% c('131', '145', '4', '184')), 'nodule.value'], 
       e$data[which(e$data$RO.names %in% c('2', '186', '187', '200')), 'nodule.value'])

t.test(e$data[which(e$data$RO.names %in% c('131', '145', '4', '184')), 'RO.value'], 
       e$data[which(e$data$RO.names %in% c('2', '186', '187', '200')), 'RO.value'])

m <- lm(RG.value ~ RO.value, RO)


eqn = paste('y', '=', format(unname(coef(m)[1]), digits=2), '+', format(unname(coef(m)[2]), digits=2), '*', 'x', sep='')

ggplot(RO, aes(x=RO.value, y=RG.value)) + stat_smooth(method='lm', se=F) + 
  geom_point(aes(color=fix)) + 
  geom_label_repel(aes(label=RO.names, size=NULL), show.legend = F) +
  xlab('Mean Occupancy in co-inoculation') + ylab('Mean RG in co-inoculation') +
  geom_text(x=0.25, y=2.9, label=eqn)



##########################
# Fix- with Fix- or Fix+ #
##########################
relative_occupancy(subset.186[which(subset.186$fix_combination == '+-'),])
relative_occupancy(subset.186[which(subset.186$fix_combination == '--'),])

relative_occupancy(subset.187[which(subset.187$fix_combination == '+-'),])
relative_occupancy(subset.187[which(subset.187$fix_combination == '--'),])

relative_occupancy(subset.200[which(subset.200$fix_combination == '+-'),])
relative_occupancy(subset.200[which(subset.200$fix_combination == '--'),])

relative_occupancy(subset.2[which(subset.2$fix_combination == '+-'),])
relative_occupancy(subset.2[which(subset.2$fix_combination == '--'),])

########################
# Goodness of fit test #
########################

class.subset.2  <- aggregate(count~strain_B+class, data=subset.2, FUN=sum)
class.subset.4 <- aggregate(count~strain_B+class, data=subset.4, FUN=sum)
class.subset.131 <- aggregate(count~strain_B+class, data=subset.131, FUN=sum)
class.subset.156 <- aggregate(count~strain_B+class, data=subset.156, FUN=sum)
class.subset.184 <- aggregate(count~strain_B+class, data=subset.184, FUN=sum)
class.subset.186 <- aggregate(count~strain_B+class, data=subset.186, FUN=sum)
class.subset.187 <- aggregate(count~strain_B+class, data=subset.187, FUN=sum)
class.subset.200 <- aggregate(count~strain_B+class, data=subset.200, FUN=sum)

# dcast all subset by the class/category
class.summary.4 <- dcast(class.subset.4, strain_B~class, value.var = 'count')
class.summary.4$strain_A <- '4'
class.summary.2 <- dcast(class.subset.2, strain_B~class, value.var = 'count')
class.summary.2$strain_A <- '2'
class.summary.131 <- dcast(class.subset.131, strain_B~class, value.var = 'count')
class.summary.131$B <- NA
class.summary.131$strain_A <- '131'

class.summary.156 <- dcast(class.subset.156, strain_B~class, value.var = 'count')
class.summary.156$strain_A <- '156'
class.summary.184 <- dcast(class.subset.184, strain_B~class, value.var = 'count')
class.summary.184$strain_A <- '184'
class.summary.186 <- dcast(class.subset.186, strain_B~class, value.var = 'count')
class.summary.186$strain_A <- '186'
class.summary.187 <- dcast(class.subset.187, strain_B~class, value.var = 'count')
class.summary.187$strain_A <- '187'
class.summary.200 <- dcast(class.subset.200, strain_B~class, value.var = 'count')
class.summary.200$A <- NA
class.summary.200$strain_A <- '200'

class.summary <- do.call("rbind", list(class.summary.4, class.summary.2, class.summary.131, class.summary.156,
                                       class.summary.184, class.summary.186, class.summary.187, class.summary.200))

#

class.summary[is.na(class.summary)] <- 0


class.summary$chi.p.value <- apply(class.summary, 1, function(x) chisq.test(c(as.numeric(x[2]), as.numeric(x[3]), as.numeric(x[4])), p=c(1/3, 1/3, 1/3))$p.value)

class.summary <- class.summary %>% filter(!duplicated(paste0(pmax(strain_A, strain_B), pmin(strain_A, strain_B))))
class.summary$dominant <- class.summary$strain_B
class.summary[which(class.summary$B < class.summary$A),]$dominant <- class.summary[which(class.summary$B < class.summary$A),]$strain_A

dominant.graph <- class.summary[, c('strain_A', 'strain_B', 'A', 'A+B', 'B', 'chi.p.value', 'dominant')]

# Use inoculum ratio from quantitative culture in the class.summary to test X^2 values
inoculum = data.frame(strain=c(2, 186, 4, 131, 187, 156, 184, 200),
           meanLog=c(8.438, 8.088, 8.5, 8.544, 8.287, 8.151, 8.389, 8.088))

# Update class.summary for dividing coinfected nodules into only two groups
class.summary$calcA <- class.summary$A + class.summary$`A+B`/2
class.summary$calcB <- class.summary$B + class.summary$`A+B`/2

class.summary <- merge(class.summary, inoculum, by.x='strain_B', by.y='strain')
colnames(class.summary)[10] <- 'meanLog_B'
class.summary <- merge(class.summary, inoculum, by.x='strain_A', by.y='strain')
colnames(class.summary)[11] <- 'meanLog_A'

class.summary$inoculum.chi.p.value <- apply(class.summary, 1, function(x) chisq.test(c(as.numeric(x[8]), as.numeric(x[9])), 
                                               p=c(as.numeric(x[10])/(as.numeric(x[10])+as.numeric(x[11])), as.numeric(x[11])/(as.numeric(x[10])+as.numeric(x[11]))))$p.value)

class.summary[which(class.summary$inoculum.chi.p.value >= 0.05),]
########################################
# Prokka gene-presence absence summary #
########################################
pangenome <- read.csv('../data/gene_presence_absence.csv')
colnames(pangenome) <- gsub('X', '', colnames(pangenome))
panmat <- pangenome[, c('Gene', 'Annotation', '131', '156', '184', '186', '187', '2', '200', '4')]

# Genes present in 131 only
unique.131 <- panmat[which(panmat$`131` != '' &
               panmat$`156` == '' &
               panmat$`184` == '' &
               panmat$`186` == '' &
               panmat$`187` == '' &
               panmat$`200` == '' &
               panmat$`2` == '' &
               panmat$`4` == ''),]
dim(unique.131)
unique(unique.131$Annotation)
View(unique.131[, c('Gene', 'Annotation')])
# Genes present in 156 only
unique.156 <- panmat[which(panmat$`131` == '' &
                             panmat$`156` != '' &
                             panmat$`184` == '' &
                             panmat$`186` == '' &
                             panmat$`187` == '' &
                             panmat$`200` == '' &
                             panmat$`2` == '' &
                             panmat$`4` == ''),]

unique(unique.156$Annotation)
View(unique.156[, c('Gene', 'Annotation')])
dim(unique.156)
# Genes present in 131+156 only
unique.131.156 <- panmat[which(panmat$`131` != '' &
                                            panmat$`156` != '' &
                                            panmat$`184` == '' &
                                            panmat$`186` == '' &
                                            panmat$`187` == '' &
                                            panmat$`200` == '' &
                                            panmat$`2` == '' &
                                            panmat$`4` == ''),]

unique(unique.131.156$Annotation)
dim(unique.131.156)
View(unique.131.156[, c('Gene', 'Annotation')])
# Make a list of genes for each category to use for GO enrichment for functions


# GO Enrichment of biological process terms of unique genes
library('rbioapi')

# What datasets are available
annots <- rba_panther_info(what = "datasets")
orgs <- rba_panther_info(what = "organisms")

genes131 <- c('AbaF','aes','caiD', 'ChaA','cnrB', 'CusA', 'DnaK', 'egtE','etfA','etfB',
           'EthA','LigD','mpa','nagX','phdht','PolX','queD','scoB','XerC','YflN')

genes131 <- c("A0A1U9YWS8", "A0A1U9YWS8", "A0A1U9Z131", "A0A1U9Z7T7", "A0A2A5KZJ9", "A0A2A6M516", "A0A6P2N956", "A0A7Z8ZAR0", "A0A7Z8ZAR0","O34409",
              "P9WQN5", "Q10725", "Q7D515", "A0A192TJ09", "A0A192TJ09", "A0A1B9TKE5", "A0A1U9YWT3", "A0A1U9Z1M4", "A0A1U9Z3S0", "A0A2A6I5L7", "A0A2A6I9W7",
              "A0A410VAL4", "A0A7S7Q177", "A0A7W6UKT1", "A0A7Z2SPG3")

enriched <- rba_panther_enrich(genes = genes131,
                               organism = 224911,
                               annot_dataset = "ANNOT_TYPE_ID_PANTHER_PATHWAY",
                               cutoff = 0.05)

enriched$result

