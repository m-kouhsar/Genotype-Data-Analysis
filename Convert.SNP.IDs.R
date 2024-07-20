library(biomaRt)
library(data.table)
arguments <- commandArgs(T)

Input.File = arguments[1]
Out.Prefix <- arguments[2]
ID.col <- as.numeric(arguments[3])
Genome.Ver <- arguments[4]
Sourc.ID <- arguments[5]

Genome.Ver <- match.arg(Genome.Ver , choices = c("GRCh37","GRCh38"))
Sourc.ID <- match.arg(Sourc.ID , choices = c("rsid" , "chr:pos"))

if(Genome.Ver=="GRCh37"){
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp",host = "https://grch37.ensembl.org", verbose = T)
}else if(Genome.Ver=="GRCh38"){
  snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", verbose = T)
}
snp_attributes = c("refsnp_id", "chr_name", "chrom_start","chrom_end","allele","allele_1")
snp_ids = fread(file = Input.File , header = F , stringsAsFactors = F)

if(Sourc.ID == "rsid" ){
  
  snp_filters = c("snp_filter")
  results = getBM(attributes=snp_attributes, filters=snp_filters, values=snp_ids[,ID.col], mart=snp_mart)
  
} else if(Sourc.ID == "chr:pos"){
  
  snp_filters = c("chr_name","start","end")
  results = getBM(attributes=snp_attributes, filters=snp_filters, values=snp_ids[,ID.col], mart=snp_mart)
  
}

fwrite(results, file = paste0(Out.Prefix , "SNP.ID.txt"), row.names = F , sep = "\t",quote = F)





