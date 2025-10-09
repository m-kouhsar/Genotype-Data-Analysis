## script to plot PCA output and colour by 1000G ethnicity
library(ggplot2)
library(ggpubr)
library(dplyr)


calcPopDist<-function(dat.pca, ref){
# dat.pca is test individual's data across pcs
# ref is matrix with 1 row per population and 1 column per pca
	popDist<-rep(NA, nrow(ref))
	names(popDist)<-rownames(ref)
	for(i in 1:nrow(ref)){
		diffs<-dat.pca - ref[i,]
		sqdiffs<-diffs^2
		popDist[i]<-sqrt(sum(sqdiffs))
	}
	return(popDist)
}


args<-commandArgs(TRUE)

#args <- c("C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/Genotyping/QC/GRCh37",
#          "NIH37","C:/Users/mk693/OneDrive - University of Exeter/Desktop/2021/NIH/Data/Genotyping/Ref/GRCh37_HG19_Referance")

setwd(args[1])
prefix<-args[2]
refSamples<-args[3]
#refPed <- args[4]
#refInfo <- args[5]

pcas<-read.table(paste0(prefix,"_mergedw1000G.pca.eigenvec"), stringsAsFactors = FALSE)
pop <- read.csv(refSamples, header = T, stringsAsFactors = F)
names(pop) = c("Sample" , "Population" , "SuperPopulation")
pcas = pcas[,-1]
colnames(pcas) = c("Sample" , paste0("PC" , c(1:(ncol(pcas)-1))))

pcas <- left_join(pcas , pop , by="Sample")
pcas$Population[is.na(pcas$Population)] = "OurSamples"
pcas$SuperPopulation[is.na(pcas$SuperPopulation)]="OurSamples"




# pcs <- pcas
# pcs <- pcs[,c(-1)]
# colnames(pcs)[2:dim(pcs)[2]] <- paste0("PC",c(1:(dim(pcs)[2]-1)))
# colnames(pcs)[1] <- "sample"
# pcs$population <- pop$SuperPopulation[match(pcs$sample, pop$Sample)]
# pcs[is.na(pcs)] <- " "
# pcs$population[which(pcs$population== " ")]<-"OurSamples"
# pcs$population <- as.factor(pcs$population)
pcs_our <- pcas[pcas$Population=="OurSamples",]
p1 = ggplot(pcs_our, aes(PC1, PC2)) + geom_point() 
p2 = ggplot(pcs_our, aes(PC1, PC3)) + geom_point() 
p3 = ggplot(pcs_our, aes(PC2, PC3)) + geom_point() 
pdf(paste0(prefix,"_PCAplot.pdf"), width = 10, height = 10)
ggarrange(p1,p2,p3 ,nrow = 2 , ncol = 2)
dev.off()

p1 = ggplot(pcas, aes(PC1, PC2,colour =Population)) + geom_point() 
p2 = ggplot(pcas, aes(PC1, PC3,colour =Population)) + geom_point() 
p3 = ggplot(pcas, aes(PC1, PC4,colour =Population)) + geom_point() 
p4 = ggplot(pcas, aes(PC2, PC3,colour =Population)) + geom_point() 
p5 = ggplot(pcas, aes(PC2, PC4,colour =Population)) + geom_point() 
p6 = ggplot(pcas, aes(PC3, PC4,colour =Population)) + geom_point()

p7 = ggplot(pcas, aes(PC1, PC2,colour =SuperPopulation)) + geom_point() 
p8 = ggplot(pcas, aes(PC1, PC3,colour =SuperPopulation)) + geom_point() 
p9 = ggplot(pcas, aes(PC1, PC4,colour =SuperPopulation)) + geom_point() 
p10 = ggplot(pcas, aes(PC2, PC3,colour =SuperPopulation)) + geom_point() 
p11 = ggplot(pcas, aes(PC2, PC4,colour =SuperPopulation)) + geom_point() 
p12 = ggplot(pcas, aes(PC3, PC4,colour =SuperPopulation)) + geom_point()

pdf(paste0(prefix,"_PCAplotwith1KG.pdf"), width = 10, height = 10)
ggarrange(p1,p2,p3,p4,p5,p6 ,nrow = 3 , ncol = 2 , common.legend = T , legend = "right")
ggarrange(p7,p8,p9,p10,p11,p12 ,nrow = 3 , ncol = 2 , common.legend = T , legend = "right")
dev.off()


# KGped<-read.table(refPed, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
# popInfo<-read.table(refInfo, stringsAsFactors = FALSE, header = TRUE, sep = "\t") ## table made from 1000G website


# 
# # KGped<-KGped[match(pcas[,2], KGped[,2]),]
# nPops<-length(table(pcas$Population))
# # popInfo<-popInfo[match(popInfo$Population.Code,levels(as.factor(KGped$Population))),]
# # 
# # KGped<-cbind(KGped,popInfo$Super.Population.Code[match(KGped$Population, popInfo$Population.Code)])
# # colnames(KGped)[ncol(KGped)]<-"SuperPopulation"
# nSuperPops<-length(table(pcas$SuperPopulation))
# 
# 
# 
# ptType<-ifelse(pcas$Population == "OurSamples", 3, 20)
# ptCol<-rainbow(nPops)[as.factor(pcas$Population)]
# ptCol[pcas$Population == "OurSamples"]<-"black"
# 
# pdf(paste0(prefix,"_PCAplotwith1KG.pdf"), width = 10, height = 10)
# layout(matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3), widths = c(3, 3, 1))
# par(mar = c(4, 4, 2, 1))
# plot(pcas[,2], pcas[,3], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
# plot(pcas[,2], pcas[,4], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
# plot(pcas[,2], pcas[,5], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
# plot(pcas[,3], pcas[,4], xlab = "PC2", ylab = "PC3", pch = ptType, col = ptCol)
# plot(pcas[,4], pcas[,5], xlab = "PC3", ylab = "PC4", pch = ptType, col = ptCol)
# par(mar = c(0, 0, 0, 0))
# plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
# legend("center", pch = 16, col = rainbow(nPops), levels(as.factor(pcas$Population)), cex = 1.5, ncol = 2)
# 
# 
# ## alternatively plot "super populations"
# 
# 
# ptCol<-rainbow(nSuperPops)[as.factor(pcas$SuperPopulation)]
# ptCol[pcas$SuperPopulation == "OurSamples"]<-"black"
# layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
# plot(pcas[,2], pcas[,3], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
# plot(pcas[,2], pcas[,4], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
# plot(pcas[,2], pcas[,5], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
# plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
# legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(pcas$SuperPopulation)), cex = 1.5)
# 
# dev.off()



## for each super population calculate cluster medians
message("Calculating population means")
nMatches<-rep(NA,20) 
for(nPCs in 2:20){
	pop.medians<-apply(pcas[,-1][,1:nPCs], 2,aggregate, by = list(pcas$SuperPopulation), median)
	pop.medians<-cbind.data.frame(pop.medians)
	rownames(pop.medians)<-pop.medians[,1]
	pop.medians<-pop.medians[,seq(2,nPCs*2,2)]

	## for each individual compare to each super population and find most similar
	popDistsAll<-apply(pcas[,-1][,1:nPCs], 1, calcPopDist, pop.medians)
	popDistsAll<-t(popDistsAll)
	predPop<-colnames(popDistsAll)[apply(popDistsAll, 1, which.min)]

	compTrue<-table(predPop, pcas$SuperPopulation)
	nMatches[nPCs]<-sum(diag(compTrue))
  print(nPCs)
}


pdf(paste0(prefix,"_SelectOptimalnPCsForPopulationPrediction.pdf"))
plot(1:20,nMatches/sum(!is.na(pcas$SuperPopulation))*100, xlab = "nPCs", ylab = "Percentage Correct")
dev.off()

nPCs<-which.max(nMatches)
pop.medians<-apply(pcas[,-1][,1:nPCs], 2,aggregate, by = list(pcas$SuperPopulation), median)
pop.medians<-cbind.data.frame(pop.medians)
rownames(pop.medians)<-pop.medians[,1]
pop.medians<-pop.medians[,seq(2,nPCs*2,2)]

## for each individual compare to each super population and find most similar
popDistsAll<-apply(pcas[,-1][,1:nPCs], 1, calcPopDist, pop.medians)
popDistsAll<-t(popDistsAll)
predPop<-colnames(popDistsAll)[apply(popDistsAll, 1, which.min)]

## calculate a quality score for prediction
## ideally want one population much closer than the others
rangeDist<-t(diff(apply(popDistsAll,1,range)))
qsPred<-(apply(popDistsAll,1,quantile, 0.25)-apply(popDistsAll,1,min))/rangeDist
pdf(paste0(prefix,"_BoxplotPrePopQCScores.pdf"), width = 12, height = 6)
par(mfrow = c(1,2))
boxplot(qsPred ~ pcas$SuperPopulation, col = rainbow(5), xlab = "Known populations")
boxplot(qsPred ~ predPop, col = rainbow(5), xlab = "Predicted populations")
dev.off()

## can define thresholds for each population based on 1000 genomes samples
#pop99<-aggregate(popDistsAll, by = list(KGped$SuperPopulation), quantile, 0.95)
#pop99Thres<-diag(as.matrix(pop99[,-1]))
#names(pop99Thres)<-colnames(popDistsAll)
#table(apply(popDistsAll, 1, min) < pop99Thres[predPop],as.factor(predPop))

print("Plotting predicted populations")
nSuperPops<-length(table(pcas$SuperPopulation))
ptType<-ifelse(pcas$Population == "OurSamples", 3, 20)
ptCol<-rainbow(nSuperPops)[as.factor(predPop)]
pdf(paste0(prefix,"_PCAplotwith1KGpredictedPopulations.pdf"), width = 10, height = 10)
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,2], pcas[,3], xlab = "PC1", ylab = "PC2", pch = ptType, col = ptCol)
plot(pcas[,2], pcas[,4], xlab = "PC1", ylab = "PC3", pch = ptType, col = ptCol)
plot(pcas[,2], pcas[,5], xlab = "PC1", ylab = "PC4", pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(pcas$SuperPopulation)), cex = 1.5)
plot(pcas[,2], pcas[,3], xlab = "PC1", ylab = "PC2", type = "n")
points(pcas[which(ptType == 20),2], pcas[which(ptType == 20),3], pch = 20, col = ptCol[which(ptType == 20)])
plot(pcas[,2], pcas[,3], xlab = "PC1", ylab = "PC2", type = "n")
points(pcas[which(ptType == 3),2], pcas[which(ptType == 3),3], pch = 3, col = ptCol[which(ptType == 3)])

dev.off()


## look at predications of our sample
outPred<-cbind(pcas[,1:2], predPop, qsPred, popDistsAll)
outPred<-outPred[which(ptType == 3),]
write.csv(table(outPred$predPop), paste0(prefix,"_TablePredictedPopulations.csv"))
write.csv(outPred, paste0(prefix,"_PredictedPopulations.csv"), quote = FALSE, row.names = FALSE)

#fam <- read.table(paste0(prefix,".fam"), stringsAsFactors = FALSE)

#outliers <- pcs_our$sample[pcs_our$PC2< 0.015]
#outliers1 <- fam[fam$V2 %in% outliers,]
#write.table(outliers1[,c(1,2)],file = paste0(prefix,"_EthnicityOutliers.txt"),quote = F,row.names = F,col.names = F)		
