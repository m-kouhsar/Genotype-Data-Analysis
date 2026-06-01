## Written by Eilis
## takes .info files (1 per chromosome) and creates plots to summarise imputation


args<-commandArgs(trailingOnly = TRUE)

#args <- c("/lustre/projects/Research_Project-191391/Genotyping/Imputation_Results/All","/lustre/projects/Research_Project-191391/Genotyping/Ref/GRCh37_HG19_Referance/1000GP_Phase3_combined.legend","ALL")

library(data.table)
InDir<-args[1]
OutDir<-args[2]
refFile<-args[3]

print("Reading the reference data...")
refPanel<-fread(refFile, data.table = FALSE,nThread=8)
if(!"id" %in% colnames(refPanel)){
	refPanel$id<-paste(refPanel$"#CHROM", refPanel$POS, refPanel$REF, refPanel$ALT, sep = ":")
}


compRefPanel<-matrix(data = NA, ncol = 2, nrow = 22)
colnames(compRefPanel)<-c("cor", "MAD")
nVar<-matrix(data = 0, nrow = 10, ncol = 10)
for(chr in 1:22){
  print(paste("Plotting Chr",chr,"..."))
	png(paste(OutDir,"/ImputationQualityPlots_chr", chr, ".png", sep = ""), width = 12, height = 8, units = "in", res = 200)
	par(mfrow = c(2,2))
	
	imputScores<-read.table(gzfile(paste(InDir,"/chr", chr, ".info.gz", sep = "")), header = TRUE, na.strings = "-")
	index<-match(imputScores$SNP,refPanel$id)

	plot(density(imputScores$Rsq, na.rm = TRUE), main = paste("chr", chr, sep = ""), xlab = "Rsq")

	## compare freq against refPanel
	plot(refPanel[index,"ALL"],imputScores$ALT_Frq,xlab = "RefPanel MAF", ylab = "Imputed sample MAF", pch = 16)
	## calculate median absolute deviation
	mtext(side = 3, line = 0.5, adj = 1, paste("MAD =", signif(median(abs(refPanel[index,"ALL"] -imputScores$ALT_Frq), na.rm = TRUE), 3), sep = ""))
	
	compRefPanel[chr,1]<-cor(refPanel[index,"ALL"],imputScores$ALT_Frq, use = "p")
	compRefPanel[chr,2]<-median(abs(refPanel[index,"ALL"] -imputScores$ALT_Frq), na.rm = TRUE)	

	## summarise INFO score distribution by MAF	
	tabRsqMAF<-table(cut(imputScores$Rsq, seq(0,1,0.1)), cut(imputScores$MAF, seq(0,0.5,0.05)))
	barplot(tabRsqMAF,ylab = "N variants", legend = rownames(tabRsqMAF), xlab = "MAF", col = colorRampPalette(c("white", "navy"))(10), main = paste("chr", chr, sep = ""))
	tabRsqMAFPer<-t(t(tabRsqMAF)/colSums(tabRsqMAF)*100)
	barplot(tabRsqMAFPer,ylab = "%", xlab = "MAF", col = colorRampPalette(c("white", "navy"))(10))

	nVar<-nVar+tabRsqMAF
	

	dev.off()
}
print("Writing summary reports...")
write.csv(nVar, paste(OutDir,"/nVariantsSummary.csv", sep = ""))
write.csv(compRefPanel, paste(OutDir,"/CompRefPanelSummary.csv", sep = "")) 


