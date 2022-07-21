################################################################################
# Visualize competition gene matrix based on Wheatley 2020 from BLAST analysis #
################################################################################

# Load competition gene matrix
matBacteroidSpecific <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/Wheatley2022/competition_matrix_BacteroidSpecific.csv', sep=',', header=T)
matNoduleBacteria <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/Wheatley2022/competition_matrix_NoduleBacteria.csv', sep=',', header=T)
matNoduleGeneral <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/Wheatley2022/competition_matrix_NoduleGeneral.csv', sep=',', header=T)
matRhizosphereRoot <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/Wheatley2022/competition_matrix_RhizosphereRoot.csv', sep=',', header=T)
matRhizosphereProgressive <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/Wheatley2022/competition_matrix_RhizosphereProgressive.csv', sep=',', header=T)
matRootSpecific <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/Wheatley2022/competition_matrix_RootSpecific.csv', sep=',', header=T)
matSuarez <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/competition_matrix_suarez.csv', sep=',', header=T)

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

colnames(finalMat)


############################### 
# Create pan-gene matrix file #
###############################

# Mar 11, 2022

# Open gene presence-absence matrix from pan-genome
matrix <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter3_Epidemic_Genotype/GWAS_analysis/pan_gene_presence.csv', sep=',', header=T)
genes <- matrix$Gene

# Remove columns not important
matrix = matrix[,16:280]

# Update row and column names
colnames(matrix) <- gsub('X', '', colnames(matrix))
rownames(matrix) <- genes  

# Subset so that it contain only 8 genomes from competition experiment
genomes <- c('13LoS28_1', '11LoS34_4', '11LoS34_10', '11LoS6_2', '11LoS7_1', '13LoS78_1', 'inoc2', 'inoc4.2')

matrix <- matrix[, genomes]

matrix <- matrix[rowSums(matrix[])>0,]

# update rownames 
genes <- gsub('.fna', '', rownames(matrix))
genes <- stringr::str_split_fixed(genes, "_", 2)[,2]
matrix$gene <- genes

matrix <- matrix[which(matrix$gene != 'hypothetical_protein'),]
matrix$gene <- stringr::str_split_fixed(matrix$gene, "_", 2)[,1]

matrix[1:20,]

############################################
## Searching gene as string to find match  #
############################################


# Genetic components for Rhizosphere colonization
rcol_genes <- c('mot', 'flg', 'fli', 'aap', 'liv', 'acpXL', 'mosA', 'mosB', 'phaC', 'bdhA', 'rhaR', 'rhaS',
                'iolE', 'iolA', 'iolC', 'glpD', 'glpQ', 'glpK')


mat1 <- matrix[grep('cobF', matrix$gene), ]

for (i in rcol_genes){
  mat1 <- rbind(mat1, matrix[grep(i, matrix$gene),])
}

colnames(mat1) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
col<- colorRampPalette(c("white", "blue"))(1000)
heatmap(data.matrix(mat1[,1:8]), scale="none", col = col)

# Genetic components to prevent the growth of other bacterial cells
pgb_genes <- c('prsD', 'prsE', 'nopP', 'nopX', 'virB')

mat2 <- matrix[grep('cobF', matrix$gene), ]

for (i in pgb_genes){
  mat2 <- rbind(mat2, matrix[grep(i, matrix$gene),])
}


colnames(mat2) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
col<- colorRampPalette(c("white", "blue"))(500)
heatmap(data.matrix(mat2[,1:8]), scale="none", col = col)

# Genetic component to establish an efficient symmbiosis
mat3 <- matrix[grep('nodD', matrix$gene), ] # nodulation
colnames(mat3) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat3[,1:8]), scale='none', col=col)

# Genetic components to promote plant growth
pgp_genes <- c('trpE', 'trpF', 'trpC', 'trpB', 'phoR', 'phoU', 'phoB', 'phoC', 'phoA', 'rhtA', 'rhtB')

mat4 <- matrix[grep('cobF', matrix$gene), ]

for (i in pgp_genes){
  mat4 <- rbind(mat4, matrix[grep(i, matrix$gene),])
}

colnames(mat4) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat4[,1:8]), scale='none', col=col)

#####################################
# Based on Wheatly et al 2020 paper #
#####################################

# Rhizosphere specific
rs_genes <- c('sodB', 'hfq', 'cspA', 'pth')

mat5 <- matrix[grep('cobF', matrix$gene), ]

for (i in rs_genes){
  mat5 <- rbind(mat5, matrix[grep(i, matrix$gene),])
}


colnames(mat5) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat5[,1:8]), scale='none', col=col, main='Rhizosphere Specific Genes')


# Rhizosphere progressive
rp_genes <-c('cobF', 'fixN', 'trpB', 'trpA', 'aroA', 'gpmA', 'dapB', 'lepA', 'rpoN', 'fmt', 'metZ', 'edd',
          'gmd', 'bioY', 'ctaC', 'ctaD', 'rnhA', 'rosR', 'feuQ', 'relA', 'pepA', 'ppx', 'glyA', 'ribG', 
          'ribC', 'pyrC', 'pyrB', 'ilvD', 'icdB', 'eno', 'dus', 'ntrB', 'ntrC', 'ccdA', 'ppiD', 'trpD',
          'trpC', 'rLuC', 'deaD', 'rpe', 'purB', 'gor', 'ropA', 'cobU', 'covV', 'cobO', 'ilvC', 'ilvI', 
          'greA', 'petA', 'leuA', 'trpE', 'bacA', 'pssN', 'noeK', 'rLuD', 'purA', 'rubA', 'ruvB', 'cbbT',
          'ptsP', 'argJ', 'ftsE', 'ftsX', 'tyrC', 'hisC', 'cobT', 'cobS', 'typA', 'ccmA', 'fbpB', 'metA',
          'ctpA', 'leuB', 'leuS')


mat6 <- matrix[grep('cobF', matrix$gene), ]

for (i in rp_genes){
  mat6 <- rbind(mat6, matrix[grep(i, matrix$gene),])
}

colnames(mat6) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat6[,1:8]), scale='none', col=col, main='Rhizosphere Progressive')

# Rhizosphere and root specific
matrix[grep('eda', matrix$gene), ]

# Root specific genes

r_genes <- c('repA', 'repB', 'recC', 'pssO', 'tlpA', 'hslR')

mat7 <- matrix[grep('cobF', matrix$gene), ]

for (i in r_genes){
  mat7 <- rbind(mat7, matrix[grep(i, matrix$gene),])
}

colnames(mat6) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat6[,1:8]), scale='none', col=col, main='Root Specific')


# Root progressive
rp_genes <- c('hemN', 'rph', 'truA', 'hisH', 'guaB', 'amn', 'nodM', 'lspL', 'gap', 'ccmB',
              'cycZ', 'cycY')

mat7 <- matrix[grep('cobF', matrix$gene), ]

for (i in rp_genes){
  mat7 <- rbind(mat7, matrix[grep(i, matrix$gene),])
}

colnames(mat7) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat7[,1:8]), scale='none', col=col, main='Root Progressive')

# Nodule General
ng_genes <- c('cspA', 'gabD', 'nifN', 'nifK', 'nifD', 'nodL', 'nodE', 'nodF', 'nodD', 'nodA', 'nodB', 'nodC', 
              'nodI', 'nodJ', 'fixB', 'fixA', 'fixN', 'rhaD', 'eryB', 'eryC', 'truB', 'rbfA', 'recF', 'fur',
              'ptsN', 'pyrC', 'frk', 'phoU', 'cheX', 'cheY', 'cheA', 'cheW', 'cheR', 'cheB', 'cheY', 'cheD',
              'fliF', 'flhB', 'motA', 'flgF', 'flgI', 'flgH', 'motB', 'motC', 'flgL', 'flaF', 'flbT', 'flDE',
              'flhA', 'iolA', 'fcl', 'gph', 'moaB', 'purD', 'cycH', 'cycJ', 'cycK', 'cycL', 'degP', 'ropA',
              'purF', 'rplI', 'purN', 'purM', 'ihfA', 'rpmG', 'scpB', 'clpS', 'cysG', 'cysI', 'queA', 'tgt', 
              'tyrS', 'purC', 'purL', 'recA', 'ilvH', 'lpcB', 'gspA', 'csaA', 'proC', 'petB', 'glnII', 'pssI',
              'pssH', 'thiE', 'purE', 'purK', 'hss', 'glgA', 'xerD', 'dacC', 'cheW', 'gpsA', 'argG', 'ispZ', 'nagA')

mat8 <- matrix[grep('cobF', matrix$gene), ]

for (i in ng_genes){
  mat8 <- rbind(mat8, matrix[grep(i, matrix$gene),])
}

colnames(mat8) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat8[,1:8]), scale='none', col=col, main='Nodule General')

# Nodule Specific
ns_genes <- c('nodM', 'fixG', 'mccc1', 'accC', 'accB', 'traI', 'hisD', 'secB', 'glcB', 'fliQ', 'dgoA', 'dgoK',
              'kup', 'aroA', 'ppdK', 'tam', 'ubiG', 'hupA', 'aat', 'cspA', 'clpA', 'ecfE', 'mefG', 'kptA', 
              'purQ', 'dksA', 'hflK', 'recN', 'pfp', 'rnsA', 'ilvD', 'rpoH', 'ruvC', 'pgk', 'rpmE', 'talB', 
              'cysZ', 'fdsB')

mat9 <- matrix[grep('cobF', matrix$gene), ]

for (i in ns_genes){
  mat9 <- rbind(mat9, matrix[grep(i, matrix$gene),])
}

colnames(mat9) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat9[,1:8]), scale='none', col=col, main='Nodule Specific')


# Bacteroid specific
bs_genes <- c('nifE', 'nifH', 'nifB', 'nifA', 'fixX', 'fixC', 'pckA', 'smf', 'dctB', 'dctD', 'petC', 'guaD')

mat10 <- matrix[grep('cobF', matrix$gene), ]

for (i in bs_genes){
  mat10 <- rbind(mat10, matrix[grep(i, matrix$gene),])
}

colnames(mat10) <- c('131', '156', '184', '186', '187', '200', '2', '4', 'gene')
heatmap(data.matrix(mat10[,1:8]), scale='none', col=col, main='Bacteroid Specific')


#############################
# Marge all matrix together #
#############################
#mat <- Reduce(function(x, y) merge(x, y, 1, all=T), list(mat1, mat2, mat3, mat4, mat5, mat6, mat7, mat8, mat9, mat10))


# Heatmap
# Table
mat <- rbind(mat5, mat6)
mat <- rbind(mat, mat7)
mat <- rbind(mat, mat8)
mat <- rbind(mat, mat9)
mat <- rbind(mat, mat10)
dim(mat)

matX <- rbind(mat1, mat2)
matX <- rbind(matX, mat3)
matX <- rbind(matX, mat4)
  
head(mat)
length(unique(mat$gene))

# Venn-diagram
library("gplots")
mat[mat == 0] <- NA

s131 <- unique(mat[!is.na(mat$`131`), 'gene'])
s156 <- unique(mat[!is.na(mat$`156`), 'gene'])
s184 <- unique(mat[!is.na(mat$`184`), 'gene'])
s186 <- unique(mat[!is.na(mat$`186`), 'gene'])
s187 <- unique(mat[!is.na(mat$`187`), 'gene'])
s200 <- unique(mat[!is.na(mat$`200`), 'gene'])
s2 <- unique(mat[!is.na(mat$`2`), 'gene'])
s4 <- unique(mat[!is.na(mat$`4`), 'gene'])

input <- list('156'=s156, '131'=s131, '184'=s184, '4'=s4)
venn(input)

setdiff(s131, s156)
setdiff(s131, s184)
setdiff(s184, s156)
setdiff(s131, s4)
setdiff(s4, s156)



# For matX
matX[matX == 0] <- NA

s131 <- unique(matX[!is.na(matX$`131`), 'gene'])
s156 <- unique(matX[!is.na(matX$`156`), 'gene'])
s184 <- unique(matX[!is.na(matX$`184`), 'gene'])
s186 <- unique(matX[!is.na(matX$`186`), 'gene'])
s187 <- unique(matX[!is.na(matX$`187`), 'gene'])
s200 <- unique(matX[!is.na(matX$`200`), 'gene'])
s2 <- unique(matX[!is.na(matX$`2`), 'gene'])
s4 <- unique(matX[!is.na(matX$`4`), 'gene'])

input <- list('156'=s156, '131'=s131, '184'=s184, '4'=s4)
venn(input)

setdiff(s131, s156)
setdiff(s131, s184)
setdiff(s184, s156)
setdiff(s131, s4)
setdiff(s4, s156)

# All (Wheatly+Suarez)
matAll <- rbind(mat, matX)

s131 <- unique(matAll[!is.na(matAll$`131`), 'gene'])
s156 <- unique(matAll[!is.na(matAll$`156`), 'gene'])
s184 <- unique(matAll[!is.na(matAll$`184`), 'gene'])
s186 <- unique(matAll[!is.na(matAll$`186`), 'gene'])
s187 <- unique(matAll[!is.na(matAll$`187`), 'gene'])
s200 <- unique(matAll[!is.na(matAll$`200`), 'gene'])
s2 <- unique(matAll[!is.na(matAll$`2`), 'gene'])
s4 <- unique(matAll[!is.na(matAll$`4`), 'gene'])

input <- list('Strain 156'=s156, 'Strain 131'=s131, 'Strin 184'=s184, 'Strain 4'=s4)
#venn(input)

setdiff(s131, s156)
setdiff(s131, s184)
setdiff(s184, s156)
setdiff(s131, s4)
setdiff(s4, s156)
setdiff(s184, s156)


library(UpSetR)

input <- list('Strain 156'=s156, 'Strain 131'=s131, 'Strain 184'=s184, 'Strain 4'=s4, 'Strain 186'=s186, 'Strain 187'=s187, 'Strain 200'=s200, 'Strain 2'=s2)

upset(fromList(input), order.by='degree', nsets=8)



matAll = rbind(matX, mat)
# Reductionist heatplot   
matAll[ is.na(matAll)] = 0

# Convert all numbers ggreater than 1 into 1 
# so that end matrix will be only presence-absence, not count
matAll.pamat <- matAll[,1:8]
matAll.pamat[matAll.pamat > 1] = 1
matAll.pamat$genes <- matAll$gene

# Remove rows with all 1
matAll.pamat <- matAll.pamat[which(rowSums(matAll.pamat[,1:8]) != 8), ]

# Need to check repeated genes have same presence/absence pattern or not
table(matAll.pamat$genes)

matAll.pamat[which(matAll.pamat$genes == 'xerD'),]

length((matAll.pamat$genes))
heatmap(as.matrix(matAll.pamat[,1:8]), scale='none',
        col=col, main='Unique Competition Genes', labRow = matAll.pamat$genes)


########################################################################
########################################################################
# Pretty Venn Diagram
# Does not work with more then 4 sets
library(VennDiagram)

input <- list('Strain 156'=s156, 'Strain 131'=s131, 'Strin 184'=s184, 'Strain 4'=s4)


display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
display_venn(input)

display_venn(input, 
             lwd = 2,
             lty = 'blank',
             fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))


##############################################################################
# Visualize competition gene matrix based on Suarez 2021 from BLAST analysis #
##############################################################################

# Load competition gene matrix
matrix <- read.csv(file='/Users/arafat/GDrive/Sachs/Chapter1_Competition/Competition_study/Competition_genes/competition_matrix_suarez.csv', sep=',', header=T)
colnames(matrix)[1] <- 'geneName'
colnames(matrix) <- gsub('X', '', colnames(matrix))
colnames(matrix) <- c('geneName', '131', '156', '184', '186', '187', '2', '200', '4', 'geneDesc')
rownames(matrix) <- matrix$geneName

matrix = matrix[which(rowSums(matrix[2:9]) != 8 & rowSums(matrix[2:9]) != 0), ]

# plot heatmap
col<- colorRampPalette(c("white", "blue"))(10)
heatmap(data.matrix(matrix[2:9]), scale='none', main='Suarez 2021', 
        labRow = matrix$geneDesc, labCol = )


