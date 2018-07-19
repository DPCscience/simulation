#################################################
###KEEP MARKERS SELECTED LINES FROM POPULATIONS 
#################################################


#F3hapfile="F3soy2.hap" in two column format
#mapfileF3="F3soy2.map"
#candidates="F3MAS2.txt" list of selected individuals to go through MAS filtered from original hapfile
#index="../index0.txt" row index of MAS position
#nsels how many individuals should be selected
#phenoF3 phenotype file for further output generation
#outF3self = descriptor can be any word

F3MAS<-function() {
args <- commandArgs(trailingOnly = TRUE)
  hapfile<-args[1]#"soy0.hap"
  mapfile<-args[2]#"soy0.map"
  outF3self<-args[3]
  candidates<-args[4]
  nsels<-args[5]
  phenoF3<-args[6]  


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

####remove markers with marker effects
qtlMAS<-read.table(file=mapfile,header=FALSE,stringsAsFactors = FALSE)
colnames(qtlMAS)<-c("mapsoy.V1","mapsoy.V2","mapsoy.V3")

####markers selected for MAS
index<-read.table(file="index0.txt",header=FALSE)
index<-c(index$V1)

###two marker selection based on position of MAS markers use the adjacent markers for scoring
###check adjacent markers
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

##valsel is the value of homozygous scoring to filter from hap file  
##i.e if 2 marker position = 4 markers (2 adjacent markers/position) = scoring = 8

MASidx<-unique(MASpos$mapsoy.V2)
MAS<-qtlMAS[qtlMAS$mapsoy.V2%in%MASidx,]
indexMAS<-as.numeric(rownames(MAS))# row MARKER position
valsel<-length(indexMAS)*2 #the *2 indicates the haplotype column = 1 1 = homozygous for the locus

########  haplotype check
genoF3MAS<-traingen3
F3MASmarkers<-genoF3MAS[indexMAS,]
genoF3MASres <- apply(F3MASmarkers, 2, as.numeric)
seq<-seq(1,(ncol(genoF3MASres)-1),2)
res<-sapply(seq,function(i) sum(genoF3MASres[,i:(i+1)])) ##what is the sum value of the haplotype loci
ID<-which(res==valsel) ##if the sum value of the haplotype loci = valsel then pass MAS

#how many lines were selected after MAS in each cycle
write.table(ID,file=paste0("MAS",outF3self,".txt"),quote = F,row.names = F,sep=" ",col.names = FALSE)
