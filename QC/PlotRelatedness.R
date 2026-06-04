library(ggplot2)

args<-commandArgs(TRUE)

wd <- trimws(args[1])
Kinship <- as.numeric(trimws(args[2]))
prefix <- trimws(args[3])

setwd(wd)

ibd <- read.table(file =paste0(prefix,"_ibd.genome"),stringsAsFactors = F,header = T )
ibd$Relatedness[ibd$RT=="UN"] <- "Unrelated"
ibd$Relatedness[ibd$RT=="HS"] <- "Half Sibs"
ibd$Relatedness[ibd$RT=="OT"] <- "Other related"
ibd$Relatedness[ibd$RT=="PO"] <- "Parent offspring"

if(file.exists(paste0(prefix,"_king.kin0"))){
  kin <- read.table(file =paste0(prefix,"_king.kin0") ,stringsAsFactors = F,header = T )
}else {
  if(file.exists(paste0(prefix,"_king.kin"))){
    kin <- read.table(file =paste0(prefix,"_king.kin") ,stringsAsFactors = F,header = T )
  }
}
  

outliers <- kin[kin$Kinship >= Kinship,] 

pdf(paste0(prefix,".relatedness.pdf"))

hist(kin$Kinship,main =paste0(prefix, " KinShip Coefficient"),xlab = "Kinship Coefficient")
abline(v = Kinship, col = 'red', lwd = 2, lty = 'dashed')

hist(ibd$PI_HAT,main =paste0(prefix, " IBD Pi_hat"),xlab = "Pi_hat")

ggplot(data = ibd,aes(Z0,Z1,color=Relatedness))+geom_point()+ggtitle(paste0(prefix," IBD with outliers"))
if(nrow(outliers)>0){
  ggplot(data = ibd[!((ibd$IID1 %in% outliers$ID1)&(ibd$IID2 %in% outliers$ID2)), ],aes(Z0,Z1,color=Relatedness))+geom_point()+
    ggtitle(paste0(prefix," IBD without outliers"))
}

dev.off()

if(nrow(outliers) > 0){
  message("[INFO] Finding the samples with lowest call rate (highest missingness) to remove as outliers...")
  setwd("..")
  system(paste0("plink --bfile ",prefix , " --missing --out " , wd , "/",prefix))
  setwd(wd)
  
  sample_missing <- read.table(file = paste0(prefix , ".imiss") , stringsAsFactors = F , header = T)
  
  index1 <- match(paste(outliers$FID1 , outliers$ID1),paste(sample_missing$FID , sample_missing$IID))
  outliers$N_MISS_1 <- sample_missing[index1 , "N_MISS"]
  outliers$F_MISS_1 <- sample_missing[index1 , "F_MISS"]
  
  index2 <- match(paste(outliers$FID2 , outliers$ID2),paste(sample_missing$FID , sample_missing$IID))
  outliers$N_MISS_2 <- sample_missing[index2 , "N_MISS"]
  outliers$F_MISS_2 <- sample_missing[index2 , "F_MISS"]
  
  outliers1 <- data.frame(FID = character(nrow(outliers)),
                          IID = character(nrow(outliers)))
  index <- apply(outliers[, c("F_MISS_1", "F_MISS_2")], 1, which.max)
  for(i in seq_along(index)) {
    outliers1[i, ] <- outliers[i,c(paste0("FID", index[i]),
                                 paste0("ID", index[i]))]
  }
  
  write.table(outliers1,file = paste0(prefix,".RelatednessOutliers.txt"),quote = F,row.names = F,col.names = F)
  write.table(outliers,file = paste0(prefix,".RelatedSamples.txt"),quote = F,row.names = F,col.names = T)
  
}else{
  warning("There is no related samples based on Kinship threshold = ",Kinship)
}

outliers1 <- cbind.data.frame(outliers$FID2,outliers$ID2)
outliers1 <- outliers1[!duplicated(outliers1),]
if(nrow(outliers1)==0){
  warning("There is no related samples based on Kinship threshold = ",Kinship)
}
		
