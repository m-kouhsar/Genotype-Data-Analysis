library(ggplot2)

args<-commandArgs(TRUE)

setwd(trimws(args[1]))
Kinship <- as.numeric(trimws(args[2]))
prefix <- trimws(args[3])

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

hist(ibd$PI_HAT,main =paste0(prefix, " IBD Pi_hat"),xlab = "Pi_hat")

ggplot(data = ibd,aes(Z0,Z1,color=Relatedness))+geom_point()+ggtitle(paste0(prefix," IBD with outliers"))
if(nrow(outliers)>0){
  ggplot(data = ibd[!((ibd$IID1 %in% outliers$ID1)&(ibd$IID2 %in% outliers$ID2)), ],aes(Z0,Z1,color=Relatedness))+geom_point()+
    ggtitle(paste0(prefix," IBD without outliers"))
}

dev.off()

outliers1 <- cbind.data.frame(outliers$FID,outliers$ID2)
outliers1 <- outliers1[!duplicated(outliers1),]
if(nrow(outliers1)==0){
  warning("There is no related samples based on Kinship threshold = ",Kinship)
}
write.table(outliers1,file = paste0(prefix,".RelatednessOutliers.txt"),quote = F,row.names = F,col.names = F)		
