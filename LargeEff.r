LGeff<-function(NoQTL,numLG,ratio,NoSNP,traitLG,traits){
#numLG vector of number of large effect QTL for each trait with large effects
#traitLG<-c("DiseaseA","DiseaseB") i.e vector of trait names with large effects
#traits<-c("maturity","grainyield") i.e vector of trait names with other effects
#NoSNP
#NoQTL
#ratio value of large effect vs small effect

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
Positions=matrix(sample(NoSNP,NoPheno*NoQTL),ncol=NoPheno)
EffectSNPs=matrix(0,nrow = NoSNP,ncol=NoPheno)

for (i in seq(1:NoPheno)) {
  EffectSNPs[Positions[,1],i]=RescaledEffects[,i]
}

return(EffectSNPs)
}

