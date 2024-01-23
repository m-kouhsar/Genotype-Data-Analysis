library(ggplot2)

args<-commandArgs(TRUE)
#args <- c("/lustre/projects/Research_Project-191391/Genotyping/Imputation_Results/PostImpQCoutput_NIH37","NIH37_PostImp_1-22_BiAllelic")
setwd(args[1])
prefix <- args[2]

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
  

outliers <- kin[kin$Kinship >= 0.25,]

pdf(paste0(prefix,"relatedness.pdf"))

hist(kin$Kinship,main =paste0(prefix, " KinShip Coefficient"),xlab = "Kinship Coefficient")

hist(ibd$PI_HAT,main =paste0(prefix, " IBD Pi_hat"),xlab = "Pi_hat")

ggplot(data = ibd,aes(Z0,Z1,color=Relatedness))+geom_point()+ggtitle(paste0(prefix," IBD with outliers"))
ggplot(data = ibd[!((ibd$IID1 %in% outliers$ID1)&(ibd$IID2 %in% outliers$ID2)), ],aes(Z0,Z1,color=Relatedness))+geom_point()+
  ggtitle(paste0(prefix," IBD without outliers"))
dev.off()

outliers1 <- cbind(outliers$FID1,outliers$ID1)
write.table(outliers1,file = paste0(prefix,"_RelatednessOutliers.txt"),quote = F,row.names = F,col.names = F)		
