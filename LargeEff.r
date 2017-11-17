
  LGeff<-function(NoQTL,numLG,ratio,NoSNP,traitLG,traits,mapfile,output,mas){
#
#numLG=c(5,5) #vector of number of large effect QTL for each trait with large effects
#traitLG<-c("DiseaseA","DiseaseB") #i.e vector of trait names with large effects
#traits<-c("maturity","grainyield") #i.e vector of trait names with other effects
#NoSNP=6909
#NoQTL=200
#ratio=0.5 #value of large effect vs small effect
#mapfile="input_pheno/soy.map"
#output rep number
#mas rep number

args <- commandArgs(trailingOnly = TRUE)
numLG<-args[1]
traitLG<-args[2]
traits<-args[3]
NoSNP<-as.numeric(args[4])
NoQTL<-as.numeric(args[5])
ratio<-as.numeric(args[6])
mapfile<-args[7]
output<-args[8]
mas<-args[9]


numLG<-as.numeric(strsplit(numLG,",")[[1]])
traitLG<-as.character(strsplit(traitLG,",")[[1]])
traits<-as.character(strsplit(traits,",")[[1]])

NoPheno=length(c(traitLG,traits))

res<-data.frame()
for (i in 1:length(traitLG))
{
ScaledVarLargeEffs = (ratio/numLG[i])/((1-ratio)/(NoQTL-numLG[i]))
res<-rbind(res,ScaledVarLargeEffs)
} 

NoSmallEffs<-data.frame()
for (i in 1:length(traitLG))
{
small = NoQTL-numLG[i]
NoSmallEffs<-rbind(NoSmallEffs,small)
} 

smallEffs<-apply(NoSmallEffs,1,function(x) rnorm(x,0,1))

largeEffs<-list()
for (i in 1:length(traitLG))
{
largeE<-rnorm(numLG[i],0,res[i,])
largeEffs[[i]]<-largeE
} 

largeEtraits<-data.frame(row.names=1:NoQTL)
for (i in 1:length(traitLG)) {
traitsLG=data.frame(c(smallEffs[,1],largeEffs[[i]]))
largeEtraits<-cbind(largeEtraits,traitsLG)
}
largeEtraits<-as.matrix(largeEtraits)

RandomEffects = matrix(rnorm(NoQTL*(NoPheno-length(traitLG))), nrow=NoQTL, ncol = NoPheno-length(traitLG))
RandomEffects = cbind(largeEtraits,RandomEffects)
RescaledEffects = RandomEffects %*% solve(chol(var(RandomEffects)))
colnames(RescaledEffects)=c(traitLG,traits)
Positions=matrix(sample(NoSNP,(NoPheno*NoQTL)),ncol=NoPheno)
colnames(Positions)=c(traitLG,traits)

EffectSNPs=matrix(0,nrow = NoSNP,ncol=NoPheno)


for (i in seq(1:NoPheno)) {
  EffectSNPs[Positions[,i],i]=RescaledEffects[,i]
}
colnames(EffectSNPs)=c(traitLG,traits)

#######qtleffect format for c++
mapsoy<-read.table(file=mapfile,header=FALSE,stringsAsFactors = FALSE)
qtleff<-data.frame(mapsoy$V1,mapsoy$V2,EffectSNPs)


qtlfile<-data.frame()
for (j in 1:NoPheno) {
 for (i in 3:ncol(qtleff))
 {
qtl<-qtleff[qtleff[i]!= 0, ]
}
qtlset<-data.frame(mapsoy.V1=paste0("Trait_",j),mapsoy.V2="additive",EffectSNPs...1.=as.character(nrow(qtl)), stringsAsFactors = FALSE)
nam<-qtl[,c(1,2,i)]
colnames(nam)<-c("mapsoy.V1","mapsoy.V2","EffectSNPs...1.")
mrk<-rbind(qtlset,nam)
qtlfile<-rbind(qtlfile,mrk)
  }

write(c("[#set]",NoPheno),file=paste0("qtlfilerep",output,".txt"),sep=" ",ncolumns = 2)
write.table(qtlfile,paste0("qtlfilerep",output,".txt"),quote = F,row.names = F,sep=" ",col.names = FALSE,append=TRUE)

#######index values for all LargeEff trait///sorting to select few for index is missing, for each trait select 5 markers
index<-as.vector(Positions[,traitLG])

largemrk<-data.frame()
for (i in 1:length(traitLG)) {
mar<-qtleff[,colnames(qtleff)%in%traitLG]
sort<-mar[order(-mar[i]),]
names<-as.numeric(rownames(sort[1:5,]))
largemrk<-rbind(largemrk,names)
}

index<-unlist(largemrk)

nom<-c(traitLG,traits)
weight<-data.frame(nom,1)

write.table(weight,file=paste0("weights",mas,".txt"),sep = ' ',row.names= FALSE,quote = F,col.names=FALSE)

write.table(index,file=paste0("index",mas,".txt"),row.names= FALSE,quote = F,col.names=FALSE)

}

LGeff()
