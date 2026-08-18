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

setwd(trimws(args[1]))
prefix<-trimws(args[2])
refSamples<-trimws(args[3])

pcas<-read.table(paste0(prefix,".ld.prune.pca.eigenvec"), stringsAsFactors = FALSE)
pop <- read.csv(refSamples, header = T, stringsAsFactors = F)
names(pop) = c("Sample" , "Population" , "SuperPopulation")
pcas = pcas[,-1]
colnames(pcas) = c("Sample" , paste0("PC" , c(1:(ncol(pcas)-1))))

pcas <- left_join(pcas , pop , by="Sample")
pcas$Population[is.na(pcas$Population)] = "InputSamples"
pcas$SuperPopulation[is.na(pcas$SuperPopulation)]="InputSamples"

eigenvals <- read.table(paste0(prefix,".ld.prune.pca.eigenval"))$V1
n_samples <- nrow(pcas)
pct_variance <- (eigenvals / n_samples) * 100

p1 = ggplot(pcas, aes(x = PC1, y = PC2,colour =Population)) + geom_point() +
  xlab(paste0("PC1 (",round(pct_variance[1],2),"%)"))+ylab(paste0("PC2 (",round(pct_variance[2],2),"%)"))

p2 = ggplot(pcas, aes(x = PC1, y = PC3,colour =Population)) + geom_point() +
  xlab(paste0("PC1 (",round(pct_variance[1],2),"%)"))+ylab(paste0("PC3 (",round(pct_variance[3],2),"%)"))

p3 = ggplot(pcas, aes(x = PC1, y = PC4,colour =Population)) + geom_point() +
  xlab(paste0("PC1 (",round(pct_variance[1],2),"%)"))+ylab(paste0("PC4 (",round(pct_variance[4],2),"%)"))

p4 = ggplot(pcas, aes(x = PC2, y = PC3,colour =Population)) + geom_point() +
  xlab(paste0("PC2 (",round(pct_variance[2],2),"%)"))+ylab(paste0("PC3 (",round(pct_variance[3],2),"%)"))

p5 = ggplot(pcas, aes(x = PC2, y = PC4,colour =Population)) + geom_point() +
  xlab(paste0("PC2 (",round(pct_variance[2],2),"%)"))+ylab(paste0("PC4 (",round(pct_variance[4],2),"%)"))

p6 = ggplot(pcas, aes(x = PC3, y = PC4,colour =Population)) + geom_point()+
  xlab(paste0("PC3 (",round(pct_variance[3],2),"%)"))+ylab(paste0("PC4 (",round(pct_variance[4],2),"%)"))

p7 = ggplot(pcas, aes(x = PC1, y = PC2,colour =SuperPopulation)) + geom_point() +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC1 (",round(pct_variance[1],2),"%)"))+ylab(paste0("PC2 (",round(pct_variance[2],2),"%)"))

p8 = ggplot(pcas, aes(x = PC1, y = PC3,colour =SuperPopulation)) + geom_point() +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC1 (",round(pct_variance[1],2),"%)"))+ylab(paste0("PC3 (",round(pct_variance[3],2),"%)"))

p9 = ggplot(pcas, aes(x = PC1, y =  PC4,colour =SuperPopulation)) + geom_point() +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC1 (",round(pct_variance[1],2),"%)"))+ylab(paste0("PC4 (",round(pct_variance[4],2),"%)"))

p10 = ggplot(pcas, aes(x = PC2, y = PC3,colour =SuperPopulation)) + geom_point() +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC2 (",round(pct_variance[2],2),"%)"))+ylab(paste0("PC3 (",round(pct_variance[3],2),"%)"))

p11 = ggplot(pcas, aes(x = PC2, y = PC4,colour =SuperPopulation)) + geom_point() +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC2 (",round(pct_variance[2],2),"%)"))+ylab(paste0("PC4 (",round(pct_variance[4],2),"%)"))

p12 = ggplot(pcas, aes(x = PC3, y = PC4,colour =SuperPopulation)) + geom_point() +
  scale_color_brewer(palette = "Set1") +
  xlab(paste0("PC3 (",round(pct_variance[3],2),"%)"))+ylab(paste0("PC4 (",round(pct_variance[4],2),"%)"))


pdf(paste0(prefix,"_PCAplotwith1KG.pdf"), width = 10, height = 10)
ggarrange(p1,p2,p3,p4,p5,p6 ,nrow = 3 , ncol = 2 , common.legend = T , legend = "right")
ggarrange(p7,p8,p9,p10,p11,p12 ,nrow = 3 , ncol = 2 , common.legend = T , legend = "right")
dev.off()

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

print("Plotting predicted populations")
nSuperPops<-length(table(pcas$SuperPopulation))
ptType<-ifelse(pcas$Population == "InputSamples", 3, 20)
ptCol<-rainbow(nSuperPops)[as.factor(predPop)]
pdf(paste0(prefix,"_PCAplotwith1KGpredictedPopulations.pdf"), width = 10, height = 10)
layout(matrix(c(1,2,3,4), ncol = 2), widths = c(3,3,3,0.75))
plot(pcas[,2], pcas[,3], xlab = paste0("PC1 (",round(pct_variance[1],2),"%)"), ylab = paste0("PC2 (",round(pct_variance[2],2),"%)"), pch = ptType, col = ptCol)
plot(pcas[,2], pcas[,4], xlab = paste0("PC1 (",round(pct_variance[1],2),"%)"), ylab = paste0("PC3 (",round(pct_variance[3],2),"%)"), pch = ptType, col = ptCol)
plot(pcas[,2], pcas[,5], xlab = paste0("PC1 (",round(pct_variance[1],2),"%)"), ylab = paste0("PC4 (",round(pct_variance[4],2),"%)"), pch = ptType, col = ptCol)
plot(0,1,type = "n", axes = FALSE, xlab = "", ylab = "")
legend("center", pch = 16, col = rainbow(nSuperPops), levels(as.factor(pcas$SuperPopulation)), cex = 1.5)
plot(pcas[,2], pcas[,3], xlab = paste0("PC1 (",round(pct_variance[1],2),"%)"), ylab = paste0("PC2 (",round(pct_variance[2],2),"%)"), type = "n")
points(pcas[which(ptType == 20),2], pcas[which(ptType == 20),3], pch = 20, col = ptCol[which(ptType == 20)])
plot(pcas[,2], pcas[,3], xlab = paste0("PC1 (",round(pct_variance[1],2),"%)"), ylab = paste0("PC2 (",round(pct_variance[2],2),"%)"), type = "n")
points(pcas[which(ptType == 3),2], pcas[which(ptType == 3),3], pch = 3, col = ptCol[which(ptType == 3)])

dev.off()


## look at predications of our sample
outPred<-cbind(pcas[,1:2], predPop, qsPred, popDistsAll)
outPred<-outPred[which(ptType == 3),]
write.csv(table(outPred$predPop), paste0(prefix,"_TablePredictedPopulations.csv"))
write.csv(outPred, paste0(prefix,"_PredictedPopulations.csv"), quote = FALSE, row.names = FALSE)
