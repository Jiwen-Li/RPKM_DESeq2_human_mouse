# RPKM_DESeq2_human_mouse
#####RPKM
##### Set working directory and load the packages

setwd("~/Desktop/RNAseqworkshop/data") # for the workshop where you have your data
library(biomaRt)
library(limma)
library(edgeR)
library(RColorBrewer)
library(ggplot2)
library(Glimma)
library(gplots)

###### load the data
directory <- "./" # "/path to your/files/"
sampleFiles=grep('count',list.files(directory),value=TRUE)
sampleFiles
# my files are not ordered 	
sampleFiles=sampleFiles[order(nchar(sampleFiles),sampleFiles)]

# Read the meta data
metadata=read.table("./samples.txt",header=T,sep="\t")
head(metadata) # 
dim(metadata)
# create a group variable
group=metadata$group

#create a DGElist for downstream analyses
x <- readDGE(sampleFiles,header=F,group=group) # be sure that group and sampleFiles are in the same order!
# Meta tags detected: __no_feature, __ambiguous, __too_low_aQual, __not_aligned, __alignment_not_unique
MetaTags <- grep("^__", rownames(x))
MetaTags
x <- x[-MetaTags, ]# we remove the MetaTags from the dataframe


## we want to have the same name for count and metadata
sampleNames<-metadata$sample
colnames(x) <- sampleNames
x
## Get gene names from our Dataframe
gData = data.frame(gencode_id = rownames(x))
gData$ensembl_gene_id = sub("\\.\\d+", "", gData$gencode_id)# we transform our gencode_id to ensembl id
rownames(x)<-gData$ensembl_gene_id #conversion of the row name into ensembl id

# Get gene name and biotype with BioMart
ensembl = useMart("ENSEMBL_MART_ENSEMBL")
ensembl = useDataset("mmusculus_gene_ensembl", ensembl) 
#For human data, do not use the line above. Use this line instead: ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
View(listAttributes(ensembl))
View(listFilters(ensembl))

martData = getBM(attributes = c("ensembl_gene_id", 	"external_gene_name", "gene_biotype", "start_position", "end_position"),					# attributes
                 filters = "ensembl_gene_id",		# identifier to filter
                 values = gData$ensembl_gene_id,	# values to filter
                 mart  = ensembl)

dim(martData)
head(martData)
length(unique(martData$ensembl_gene_id))
length(unique(martData$external_gene_name))

mart<-martData[order(martData$ensembl_gene_id),] # order by ensembl id 
head(mart)
mart <- mart[!duplicated(mart$external_gene_name),]#remove duplicate gene name
dim(mart)
head(mart)

#Make gene length file
length<-data.frame(mart$end_position - mart$start_position)
write.csv(length,file="length.csv")

# merging your annotation with your DGE list count (x$counts)
NewDGE<-merge(mart,x$counts,by.x="ensembl_gene_id",by.y="row.names") ##

head(NewDGE)# count file with name and counts
count<-NewDGE[,-c(1,2,3,4,5)] 
rownames(count)<-NewDGE[,1] # We add row names
head(count)
symbol<-NewDGE[,-c(1,4:9)] # we remove the counts ! the number in parenthesis can change based on your experiment! 4 is always constant. The number after 4: should be the number of your samples+5
rownames(symbol)<-NewDGE[,1]


#### recreate a new DGE list 
x<-DGEList(count,genes=NewDGE[,1,drop=FALSE],group=group)
x # check your data
## now add gene annotation
x$genes$symbol<-symbol
x

### Normalization with TMM
x <- calcNormFactors(x, method = "TMM")
x$samples$norm.factors


########## LIMMA DE analyses
#### 1. creating design matrix and contrast

design <- model.matrix(~ 0 + x$samples$group)
design
#change the name
colnames(design) <- levels(x$samples$group)
design

#### 2. VOOM
v_norm=voom(x,design,plot=T)
v_norm
## save the norm count for downstream visualization analyses
v_norm1=v_norm$E
write.table(v_norm1,"v_norm1.txt",sep="\t",quote=F)
# add the gene name 
logCPM<-merge(v_norm1,v_norm$genes, by.x="row.names",by.y="ensembl_gene_id" )
head(logCPM) # it is a dataframe and we can not save it as it is
colnames(logCPM)[1]<-c("ensembl_ID")
head(logCPM)
logCPM$gene_name<-logCPM$symbol$external_gene_name
logCPM$biotype<-logCPM$symbol$gene_biotype
head(logCPM)
# remove the symbol. columns
dim(logCPM)
logCPM<-logCPM[,-c(19)] #The number in c( ) is the number of your samples +2 
head(logCPM)

write.table(logCPM,file="logCPM.txt",sep="\t",quote=F,col.names=NA)

###Calculate rpkm
length<-read.csv("length.csv",header=T) 
#add row names to x to save it as csv
tags=rownames(x)
tagsfile=cpm(x)[tags,]
write.csv(tagsfile,file="cpm.csv")

####Automated (for loop) to calculate RPKM
#works only AFTER you generate logCPM file (see sections LIMMA DE Analysis and VOOM)
headerNames<-"gene name" #creates variable to add sample names to later
logCPM<-read.table("logCPM.txt",header=T) 	# reads logCPM file into workspace
rpkm_data <- data.frame(logCPM[,7])	#number used after comma is number of samples + 2, creates list with gene names as first column
# for loop to generate list of sample names into variable headerNames; calculates rpkm for all samples into variable rpkm_data
for(i in 1:length(sampleFiles)){
  headerNames <- c(headerNames,tools::file_path_sans_ext(sampleFiles[i]))  
  countsData<-data.frame(tagsfile[,i]) #Read counts data for each data into a variable
  rpkm_data<-cbind(rpkm_data,rpkm(countsData,length)) #Calculate rpkm
}
colnames(rpkm_data)<-headerNames 
write.csv(rpkm_data,file="rpkm_data.csv") #writes file “rpkm_data.csv” which contains all RPKM data of samples with corresponding gene and sample names

####DESeq2
library(DESeq2)
countData <- read.csv('countdata.csv', header = TRUE, sep = ",")
metaData <- read.csv('metadata.csv', header = TRUE, sep = ",")
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~rep + treatment, tidy = TRUE)
ddsColl <- collapseReplicates(dds, dds$condition, dds$run) #for technical replicate, sum rpkm of all technical rep
colData(ddsColl)
colnames(ddsColl)
factor(ddsColl$rep)
factor(ddsColl$treatment)
ddsColl$Group <- relevel( ddsColl$treatment, "control" )
dds <- DESeq(ddsColl)
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
bm <- getBM(
  attributes=
    c("ensembl_gene_id", "external_gene_name", "gene_biotype"),
  filter="ensembl_gene_id",
  values= rownames(res),
  mart=ensembl)

library(plyr)
library(LSD)
library(stringr)
library(genefilter)
library(dplyr)
bm <- arrange(bm, ensembl_gene_id)
head(bm)
DESeq2Features <- data.frame(ensembl_gene_id = rownames(res))
DESeq2Features$ensembl_gene_id <- as.character(DESeq2Features$ensembl_gene_id)
rowData <- dplyr::left_join(DESeq2Features, bm, by = "ensembl_gene_id")
rowData <- as(rowData, "DataFrame")
write.csv(res, file="res_paired.csv")

