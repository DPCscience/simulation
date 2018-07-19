#################################################
###KEEP MARKERS SELECTED LINES FROM POPULATIONS 
#################################################

#F3MASres<-F3MAS(
#F3hapfile="F3soy2.hap"
#mapfileF3="F3soy2.map"
#outF3self=2
#F3MAS="F3MAS2.txt"
#index="../index0.txt"
#valsel=10
#length(indexMAS)*2

F3MAS<-function() {
args <- commandArgs(trailingOnly = TRUE)
  hapfile<-args[1]#"soy0.hap"
  mapfile<-args[2]#"soy0.map"
  outF3self<-args[3]
  candidates<-args[4]#"random_base0.cad"
  nsels<-args[5]
  phenoF3<-args[6]  

#hapfile<-"F3soy2.hap" 
#mapfile<-"F3soy2.map" 
#outF3self<-2
#candidates<-"F3MAS2.txt" 
#nsels<-6000 
#phenoF3<-"soyF3pheno2.pht" #only in cycle 0 otherwise .add

####
library(data.table)
traingen1<-fread(file=hapfile,stringsAsFactors=FALSE, nrows=-1,header=FALSE,data.table=FALSE)
trainnames<-read.table(file=candidates,stringsAsFactors=FALSE,header=FALSE)
vec<-as.character(trainnames$V1)
namevec<- t(traingen1[1,])
rownames(namevec)<-NULL
common<-match(vec,namevec[,1])
indsub<-sort(c(common,common+1))
traingen2<-traingen1[,indsub]
traingen3<-traingen2[-1,]
traingen3 <- apply(traingen3, 2, as.numeric)

####INDEX FOR MAS##REPLACE BY READING INDEX####
qtlMAS<-read.table(file=mapfile,header=FALSE,stringsAsFactors = FALSE)
colnames(qtlMAS)<-c("mapsoy.V1","mapsoy.V2","mapsoy.V3")

####can be replaced by index value of MAS QTLs from step 1
index<-read.table(file="index0.txt",header=FALSE)
index<-c(index$V1)

###two marker selection
####select position of MAS markers i-1,i+1 if they are on the same chromos
MASpos<-data.frame(row.names=1:length(index))
seq<-1:length(index)
for (i in seq){
  if(qtlMAS[(index[i]-1),"mapsoy.V1"] == qtlMAS[(index[i]+1),"mapsoy.V1"]){ #row number= number of MAS markers and column total number of haplotypes=geno*2
    a<-qtlMAS[c(index[i]-1,index[i]+1),]
    MASpos<-rbind(MASpos,a)
  }else if (qtlMAS[(index[i]-1),"mapsoy.V1"] == qtlMAS[(index[i]),"mapsoy.V1"]){
    b<-qtlMAS[c(index[i]-1),]
    MASpos<-rbind(MASpos,b)
  }else if (qtlMAS[(index[i]+1),"mapsoy.V1"] == qtlMAS[(index[i]),"mapsoy.V1"]){
    d<-qtlMAS[c(index[i]+1),]
    MASpos<-rbind(MASpos,d)
  }else{NULL}
}

MASidx<-unique(MASpos$mapsoy.V2)
MAS<-qtlMAS[qtlMAS$mapsoy.V2%in%MASidx,]
indexMAS<-as.numeric(rownames(MAS))# row MARKER position
valsel<-length(indexMAS)*2

######## select two markers if values are 1 on each column 1 1 haplotype check!!
genoF3MAS<-traingen3
F3MASmarkers<-genoF3MAS[indexMAS,]
genoF3MASres <- apply(F3MASmarkers, 2, as.numeric)
seq<-seq(1,(ncol(genoF3MASres)-1),2)
res<-sapply(seq,function(i) sum(genoF3MASres[,i:(i+1)]))
ID<-which(res==valsel)

#how many lines were selected after MAS in each cycle
write.table(ID,file=paste0("MAS",outF3self,".txt"),quote = F,row.names = F,sep=" ",col.names = FALSE)
####subsetting
indres<-seq[ID]
arr<-sort(c(indres,(indres+1)))
genoMAS<-traingen2[,arr]

#####AFTER MAS prepare FOR SELFING hap and map
namesF3MASsel<-as.character(unlist(genoMAS[1,], use.names=FALSE))

######PREPARE F3 self file
F3MASsel<-data.frame(unique(namesF3MASsel),unique(namesF3MASsel))
F3MASself<-sample(as.character(F3MASsel[,1]),nsels)
F3MASselfed<-data.frame(F3MASself,F3MASself)
vec1<-as.character(F3MASselfed$F3MASself)
common2<-match(vec1,namesF3MASsel)
indsubpop2<-sort(c(common2,common2+1))
genoF3MASsel1<-genoMAS[,indsubpop2]
namesF3MASsel1<-as.character(unlist(genoF3MASsel1[1,], use.names=FALSE))
genoF3MASsel1<-genoF3MASsel1[-1,]

#####HAP
write(namesF3MASsel1,file=paste0("F3soysub",outF3self,".hap"),sep=" ",ncolumns = length(namesF3MASsel1))
write.table(genoF3MASsel1,file=paste0("F3soysub",outF3self,".hap"),quote = F,row.names = F,sep=" ",col.names = FALSE,append=TRUE)

#####MAP
write.table(qtlMAS,file=paste0("F3soysub",outF3self,".map"),quote = F,row.names = FALSE,sep=" ",col.names = FALSE)

#####SELFING
write.table(F3MASselfed,file=paste0("F3MASself",outF3self,".txt"),quote = F,row.names = F,col.names = F,sep=" ")

#####PHENOS modified to take many traits
F3MASpheno<-read.table(file=phenoF3,header=TRUE,stringsAsFactors = FALSE)
colnames(F3MASselfed)<-c("Sample","SampleID")
F3soysub_add<-merge(F3MASselfed,F3MASpheno,by="SampleID")
F3soysub_add1<-F3soysub_add[,2:ncol(F3soysub_add)]
colnames(F3soysub_add1)<-colnames(F3MASpheno)
write.table(F3soysub_add1,file=paste0("F3soysub",outF3self,".add"),quote = F,row.names = F,col.names = T,sep=" ")

}

F3MAS()

