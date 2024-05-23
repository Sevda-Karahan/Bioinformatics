library(GEOquery)

if(!file.exists("geo_downloads")) 
  dir.create("geo_downloads")

my.gse <- "GSE19188"  

if(!file.exists(paste0("./geo_downloads/",my.gse)))
  getGEOSuppFiles(my.gse, makeDirectory=T, baseDir="geo_downloads") 

my.geo.gse <- getGEO(GEO=my.gse, filename=NULL, destdir="./geo_downloads", GSElimits=NULL, GSEMatrix=TRUE, AnnotGPL=FALSE, getGPL=FALSE)
my.geo.gse <- my.geo.gse[[1]]

untar(paste0("geo_downloads/",my.gse,"/",my.gse,"_RAW.tar"), exdir=paste0("geo_downloads/",my.gse,"/CEL"))
my.cels <- list.files(paste0("geo_downloads/",my.gse,"/CEL"))

#Preparing the Phenodata
my.pdata <- as.data.frame(pData(my.geo.gse), stringsAsFactors=F)
table(rownames(my.pdata) == my.cels) #kontrol yap1yoruz
temp.rownames <- paste(rownames(my.pdata), ".CEL.gz", sep="")
table(temp.rownames == my.cels)
rownames(my.pdata) <- temp.rownames
rm(temp.rownames)
table(rownames(my.pdata) == my.cels)
write.table(my.pdata, file=paste0("geo_downloads/",my.gse,"/CEL/",my.gse,"_PhenoData.txt"), sep="\t", quote=F)

#Reading the CEL Files
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(affy)

BiocManager::install("affy")

library(affy)
cel.path <- paste0("geo_downloads/",my.gse,"/CEL")
my.affy <- ReadAffy(celfile.path=cel.path, phenoData=paste(cel.path, paste0(my.gse,"_PhenoData.txt"), sep="/"))
show(my.affy)
exprs(my.affy)
head(exprs(my.affy))

colnames(pData(my.affy))
pData(my.affy)
pData(my.affy)$title
pData(my.affy)$description

#Affy nesnesini rma nesnesine d??n????t??r?? normalize eden fonksiyon:
my.rma <- rma(my.affy, normalize=T, background=T)  #quantile normalizasyonu.
head(exprs(my.rma))

#Annotation#
my.rma@annotation
library(hgu133plus2.db) #BiocManager::install("hgu133plus2.db")
library(annotate)
library(R2HTML)

ID <- featureNames(my.rma)
Symbol<-getSYMBOL(ID,"hgu133plus2.db")
sym <- as.data.frame(Symbol)

data <- as.data.frame(exprs(my.rma)) 
data <- cbind(sym,data)  #sC<tun ekleme

i <- which(is.na(data$Symbol) == TRUE)
data<-data[-c(i),]

rownames(data) <- data[,1] #bu satD1rda amac satir isimlerini sembol sutunundaki degerlere atamak ancak hata verecek, cunku tekrarli sembol isimleri var.

X <- data.table::as.data.table(data)  #library(data.table), data.table, data.frame'e gore daha fazla islem yapabilmeye olanak saglar.
final_data <- X[,lapply(.SD,mean),"Symbol"] #Symbol sC<tununda tekrar eden genler var, bu satD1rlarD1n ortalamasD1 alD1nD1p tek bir satD1r olarak yazD1lacak.
final_data <- as.data.frame(final_data)
rownames(final_data) <- final_data[,1] 
final_data <- final_data[,-c(1)]

saveRDS(final_data,"geo_downloads/GSE19188/GSE19188_raw.RDS") #okumak iC'in final_data = readRDS("geo_downloads/GSE17536/GSE17536_raw.RDS")
final_data = readRDS("geo_downloads/GSE1988/GSE1988_raw.RDS")

final_data  = t(final_data) 
metadata = pData(my.affy)
table (rownames(final_data) == rownames(metadata)) #her iki datasetteki samplelarin ayni sirada oldugunu kontrol ettik.

final_data = as.data.frame(final_data)
final_data$state = metadata$tissue.type.ch1

# Karakter de??erlere kar????l??k gelen say??sal de??erleri tan??mlay??n
state_numeric <- ifelse(final_data$state == "tumor", 1, ifelse(final_data$state == "healthy", 0, NA))
final_data$state_numeric <- state_numeric

# final_data veri ??er??evesinin s??tun isimlerini al??n
column_names <- colnames(final_data)

# Son 5 s??tun ismini yazd??r??n
last_5_columns <- tail(column_names, 5)
print(last_5_columns)
head(final_data$state)
head(final_data$state_numeric)

write.csv(final_data,file="geo_downloads/GSE19188/GSE19188.csv")



