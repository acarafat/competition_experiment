################################################################################
# Visualize competition gene matrix based on Wheatley 2020 from BLAST analysis #
################################################################################

# Load competition gene matrix
matBacteroidSpecific <- read.csv(file='../data/competition_matrix_BacteroidSpecific.csv', sep=',', header=T)
matNoduleBacteria <- read.csv(file='../data/competition_matrix_NoduleBacteria.csv', sep=',', header=T)
matNoduleGeneral <- read.csv(file='../data/competition_matrix_NoduleGeneral.csv', sep=',', header=T)
matRhizosphereRoot <- read.csv(file='../data/competition_matrix_RhizosphereRoot.csv', sep=',', header=T)
matRhizosphereProgressive <- read.csv(file='../data/competition_matrix_RhizosphereProgressive.csv', sep=',', header=T)
matRootSpecific <- read.csv(file='../data/competition_matrix_RootSpecific.csv', sep=',', header=T)
matSuarez <- read.csv(file='../data/competition_matrix_suarez.csv', sep=',', header=T)

colnames(matBacteroidSpecific) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
colnames(matNoduleBacteria) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
colnames(matNoduleGeneral) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
colnames(matRhizosphereRoot) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
colnames(matRootSpecific) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
colnames(matRhizosphereProgressive) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
colnames(matSuarez) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')

matWheatley <- rbind(matBacteroidSpecific, matNoduleBacteria)
matWheatley <- rbind(matWheatley, matNoduleGeneral)
matWheatley <- rbind(matWheatley, matRhizosphereRoot)
matWheatley <- rbind(matWheatley, matRootSpecific)
matWheatley <- rbind(matWheatley, matRhizosphereProgressive)

  
plotMatHeatmap <- function(matrix, title) {
  colnames(matrix)[1] <- 'geneName'
  colnames(matrix) <- gsub('X', '', colnames(matrix))
  colnames(matrix) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
  #rownames(matrix) <- matrix$geneName
  matrix = matrix[which(rowSums(matrix[2:9]) != 8 & rowSums(matrix[2:9]) != 0), ]
  heatmap(data.matrix(matrix[2:9]), scale='none', main=title, 
          labRow = matrix$geneDesc)
}


# plot heatmap
#pdf('/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/wheatley.pdf',
#    width=25, height=20)
plotMatHeatmap(matWheatley, '')
plotMatHeatmap(matSuarez, '')
#dev.off()

plotMatHeatmap(matBacteroidSpecific, 'Bacteroid Specific')
plotMatHeatmap(matNoduleBacteria, 'Nodule Bacteria')
plotMatHeatmap(matNoduleGeneral, ' (Nodule General)')
plotMatHeatmap(matRhizosphereRoot, 'Rhizoshere and Root')
plotMatHeatmap(matRhizosphereProgressive, 'Rhizosphere Progressive')
plotMatHeatmap(matRootSpecific, 'Root Specific')
plotMatHeatmap(matSuarez, '')
matSuarez[which(matSuarez$geneDesc == 'M'),'geneDesc'] = 'MOSC domain-containing protein'
plotMatHeatmap(matSuarez, '')


# Combine Suarez and Wheatley
matCombined <- rbind(matSuarez, matWheatley)
dim(matCombined)
# Investigate particular contrasts
# 131 vs 156
contrast131_156 <- matCombined[which(rowSums(matCombined[,c('131','156')]) == 1) ,c('131','156', 'geneDesc')]
contrast131 <- matCombined[which(matCombined$`131`== 1 & rowSums(matCombined[,2:9]) == 1),c('131','156', 'geneDesc')]
contrast156 <- matCombined[which(matCombined$`156`== 1 & rowSums(matCombined[,2:9]) == 1),c('131','156', 'geneDesc')]
heatmap(data.matrix(contrast131_156[,c('131','156')]), labRow=contrast131_156$geneDesc)
plotMatHeatmap(matCombined, 'Combined')

# 156 vs. rest excluding 131
head(matCombined)
matCombined[which(matCombined$`156`== 1 & rowSums(matCombined[,4:9]) == 0),]
matCombined[which(matCombined$`131`== 1 & rowSums(matCombined[,4:9]) == 0),]
matCombined[which(matCombined$`131`== 1 & rowSums(matCombined[,4:9]) == 0),]


################
# PCA Analysis #
################
library(ggfortify)

colnames(matBacteroidSpecific) <- c('geneName', '131', '156', '184', '186', '187', '2', '4', '200', 'geneDesc')
colnames(matNoduleBacteria) <- c('geneName', '131', '156', '184', '186', '187', '2', '4', '200', 'geneDesc')
colnames(matNoduleGeneral) <- c('geneName', '131', '156', '184', '186', '187', '2', '4', '200', 'geneDesc')
colnames(matRhizosphereRoot) <- c('geneName', '131', '156', '184', '186', '187', '2', '4', '200', 'geneDesc')
colnames(matRhizosphereProgressive) <- c('geneName', '131', '156', '184', '186', '187', '2', '4', '200', 'geneDesc')
colnames(matRootSpecific) <- c('geneName', '131', '156', '184', '186', '187', '2', '4', '200', 'geneDesc')

matBacteroidSpecific$cat <- 'Bacteroid Specific'
matNoduleBacteria$cat <- 'Nodule Bacteria Specific'
matNoduleGeneral$cat <- 'Nodule General'
matRhizosphereRoot$cat <- 'Rhizosphere Root'
matRhizosphereProgressive$cat <- 'Rhizosphere Progressive'
matRootSpecific$cat <- 'Rhizosphere Specific'

finalMat <- rbind(matBacteroidSpecific, matNoduleBacteria)
finalMat <- rbind(finalMat, matNoduleGeneral)
finalMat <- rbind(finalMat, matRhizosphereRoot)
finalMat <- rbind(finalMat, matRhizosphereProgressive)
finalMat <- rbind(finalMat, matRootSpecific)
finalMat <- finalMat[which(rowSums(finalMat[2:9]) != 8),]

finalMatPCA <- prcomp(finalMat[2:9], center = TRUE)

autoplot(finalMatPCA, data=finalMat, loadings=T, loadings.label = TRUE, colour='cat')




