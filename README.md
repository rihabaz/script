# Auteur : Rihab AZMANI 
# script Zhang data (R)
rm(list=ls())

library(tidyverse)
library(dplyr)
library(stats)
library(ggplot2)
library(RColorBrewer)
library(S4Vectors)

setwd("/home/orlab/Desktop/Rihab/DatasetZhangEtal/") 

Mydata <- read.table("listOligo.csv" , header = T , sep <- "," , dec = ".")
str(Mydata)
## Max all 

max_all=apply(Mydata[,2:6], 1, max)
maxAll=round(max_all,3)
Mydata= cbind(Mydata , maxAll)

max_oligo=apply(Mydata[,4:6], 1, max)  
maxOligo= round(max_oligo,3)
Mydata= cbind(Mydata , maxOligo)

max_non_oligo=apply(Mydata[,2:3], 1, max)
maxNonOligo= round(max_non_oligo,3)
Mydata= cbind(Mydata , maxNonOligo)

sum_oligo=apply(Mydata[,4:6], 1, sum)
sumOligo= round(sum_oligo,3)
Mydata= cbind(Mydata , sumOligo)

sum_non_oligo=apply(Mydata[,2:3], 1, sum)
sumNonOligo= round(sum_non_oligo,3)
Mydata= cbind(Mydata , sumNonOligo)

mean_oligo=apply(Mydata[,4:6], 1, mean)
meanOligo=round(mean_oligo,3)
Mydata= cbind(Mydata , meanOligo)

mean_non_oligo=apply(Mydata[,2:3], 1, mean)
meanNonOligo=round(mean_non_oligo,3)
Mydata= cbind(Mydata , meanNonOligo)

Fc= maxOligo/maxNonOligo
FcMax=round(Fc,3)
Mydata<- cbind(Mydata , FcMax)

logfc= log2(FcMax)
logFC=round(logfc,3)
Mydata<- cbind(Mydata , logFC)


Mydata$FC_Max = NULL
FC_Max= ifelse((Mydata$logFC <0) ,(-(0.5)^logFC) , (2^logFC))
FCRealMax=round(FC_Max,3)
Mydata<- cbind(Mydata , FCRealMax)


Fc_sum= sumOligo/sumNonOligo
FcSum=round(Fc_sum,3)
Mydata= cbind(Mydata ,FcSum )

FC_mean= meanOligo/meanNonOligo 
FCMean=round(FC_mean,3)
Mydata= cbind(Mydata , FCMean)

#### Calcul du FC des neurones et Astrocytes 
Fc_neuron = maxNonOligo/maxOligo
FcMax_neuron=round(Fc_neuron,3)
Mydata <- cbind(Mydata , FcMax_neuron)

logfc_neuron= log2(FcMax_neuron)
logFC_neuron=round(logfc_neuron,3)
Mydata<- cbind(Mydata , logFC_neuron)


# Mydata$FcMax_neuron = NULL
FcMax_neuron= ifelse((Mydata$logFC_neuron <0) ,(-(0.5)^logFC_neuron) , (2^logFC_neuron))
FCRealMax_neuron=round(FcMax_neuron,3)
Mydata<- cbind(Mydata , FCRealMax_neuron)

write.csv(Mydata , file="list_oligoNonOligo.csv" , row.names = F , quote = F)
dim(Mydata) # 22458 18

##Filtre 1:
dataFiltre1 <- filter(Mydata, Mydata$maxAll >= 0.5)
dim(dataFiltre1) ## 14156 
write.csv(dataFiltre1 , file = "list_oligo_NonOligo_0.5.csv" , row.names = F , quote = F)

## suppression des génes des pericytes 

exclu=read.table("toExcludeFromZhang.csv" , header=T , sep=",")
dim(exclu) ## 170 gènes 

## exclusion 
setdiff(exclu$Name, dataFiltre1$Gene.symbol) ## Slc13a4

NewdataFiltre1 <-dataFiltre1 [!(dataFiltre1$Gene.symbol %in% exclu$Name), ] ## 13987 ==> 169 gènes éliminés 
write.table(NewdataFiltre1, file= "Oligo_NonOligo_Zhang.all.csv" ,row.names = F , quote = F)

## Fitre 2 FcMaxReal <=-1.8 | FcMaxReal >=1.8 
dataFiltre2Up <- filter(NewdataFiltre1 , FCRealMax >= 1.8)
dim(dataFiltre2Up) ## 2577 
write.csv(dataFiltre2Up , file= "Up_Zhang.csv" , row.names = F , quote = F )
dataFiltre2Down <- filter(NewdataFiltre1 , FCRealMax <= -1.8)
dim(dataFiltre2Down) #3410 
write.csv(dataFiltre2Down , file= "down_Zhang.csv" ,  row.names = F , quote = F)
dataFiltre2UpDown <- filter(NewdataFiltre1 , NewdataFiltre1$FCRealMax>= 1.8 | NewdataFiltre1$FCRealMax <= -1.8)
dim(dataFiltre2UpDown) # 5987 
write.csv(dataFiltre2UpDown , file= "Up_down_Zhang.csv" ,  row.names = F , quote = F)

## heatmap 

d= arrange(dataFiltre2Up , desc(dataFiltre2Up$FCRealMax))
donnes <- (d[,1:6])
coul = colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(donnes[,-1]) , col = coul ,  Colv = NA, Rowv = NA  ,  main ="heatmap NA FC:1.8")

## OPCs 
selec1 <- subset(dataFiltre2UpDown, select= c(Gene.symbol , Astrocytes, Neuron, OPCs, imOLs, maOLs))
dim(selec1)

FCimLOS= ifelse((selec1$OPCs < selec1$imOLs) , (-(0.5)^(log2(selec1$OPCs/selec1$imOLs))) , (2^(log2(selec1$OPCs/selec1$imOLs))))
FCVSimLOs =round(FCimLOS,3)
selec1<- cbind(selec1 , FCVSimLOs )

FCmaOLS= ifelse((selec1$OPCs < selec1$maOLs) , (-(0.5)^(log2(selec1$OPCs/selec1$maOLs))) , (2^(log2(selec1$OPCs/selec1$maOLs))))
FCVSmaOLs =round(FCmaOLS,3)
selec1<- cbind(selec1 , FCVSmaOLs )

FCimLOSMeanFCmaOLS=apply(selec1[,5:6], 1, mean)
FCmeanimLOsmaOLS =round(FCimLOSMeanFCmaOLS,3)
selec1<- cbind(selec1 , FCmeanimLOsmaOLS )

FcoPCsvsaLL= selec1$OPCs/FCimLOSMeanFCmaOLS
FcoPCSvsALL=round(FcoPCsvsaLL,3)
selec1<- cbind(selec1 , FcoPCSvsALL )

write.csv(selec1 , file ="OPCs.csv" ,  row.names = F , quote = F)

## Top markers 
TopOPCS <- arrange(selec1, desc(FcoPCSvsALL))
TopOPCs= slice( TopOPCS , 1:500)
TopOPCS10= slice(TopOPCS , 1:10)
write.csv(TopOPCs , file="topOPCS00.csv", row.names = F , quote = F)
write.csv(TopOPCS10 , file="topOPCs10.csv" , row.names = F , quote = F)
barplot(TopOPCS10$FcoPCSvsALL,horiz=T,names = TopOPCS10$Gene.symbol,las=1,col="blue")

## enrichement top 10 markers 

par(mfrow=c(1,2))
n = length(TopOPCS10$Gene.symbol)
for(i in 1:2)
{
  SecLOPCS<-subset(TopOPCS10, Gene.symbol ==TopOPCS10$Gene.symbol[i])
  OPCS = c( SecLOPCS$Astrocytes, SecLOPCS$Neuron, SecLOPCS$OPCs ,SecLOPCS$imOLs , SecLOPCS$maOLs)
  barplot( OPCS, horiz = T , col=c("darkred" , "red" , "palegreen2","palegreen3","darkgreen"), beside= TRUE , main = TopOPCS10$Gene.symbol[i] )#, names.arg = c("Astrocytes" ,"Neuron", "OPCs","imOLs" ,"maOLs"))
  i <- i+1
}


## imOLs

selec2 <- subset(dataFiltre2UpDown, select =c(Gene.symbol , Astrocytes, Neuron, OPCs, imOLs, maOLs))

FCOPCs= ifelse(( selec2$imOLs < selec2$OPCs) , (-(0.5)^(log2(selec2$imOLs/selec2$OPCs))) , (2^(log2(selec2$imOLs/selec2$OPCs))))
FCimLOsvsOPCs =round(FCOPCs,3)
selec2<- cbind(selec2 , FCimLOsvsOPCs )

FCmaOLS= ifelse((selec2$imOLs < selec2$maOLs) , (-(0.5)^(log2(selec2$imOLs/selec2$maOLs))) , (2^(log2(selec2$imOLs/selec2$maOLs))))
FCimLOsVSmaOLs =round(FCmaOLS,3)
selec2<- cbind(selec2 ,  FCimLOsVSmaOLs )

FCiOPCsMeanFCmaOLS=apply(selec2[,c(4,6)], 1, mean)
FCmeanOPCsMaOLS =round(FCiOPCsMeanFCmaOLS,3)
selec2<- cbind(selec2 , FCmeanOPCsMaOLS )

FcimOLsvsaLL=  selec2$imOLs/FCmeanOPCsMaOLS
FcimOLsvsAll=round(FcimOLsvsaLL,3)
selec2<- cbind(selec2 , FcimOLsvsAll )

write.table(selec2 , file="imOLS.csv" , row.names = F , quote = F)

### Top markers 
TopimaOLS= arrange(selec2, desc(FcimOLsvsAll))
TopimaOLS= slice( TopimaOLS , 1:500)
write.csv(TopimaOLS , file="topimaOLS500.csv" , row.names = F , quote = F)

TopimaOLS10= slice( TopimaOLS , 1:10)
write.csv(TopimaOLS10 , file="topimaOLS10.csv" , row.names = F , quote = F)
barplot(TopimaOLS10$FcimOLsvsAll,horiz=T,names = TopimaOLS10$Gene.symbol,las=1,col="orange")
## enrichement top 10 markers 
par(mfrow=c(1,2))
n = length(TopimaOLS10$Gene.symbol)
for(i in 1:n)
{
  SecLimaOLS<-subset(TopimaOLS10, Gene.symbol ==TopimaOLS10$Gene.symbol[i])
  imaOLS = c( SecLimaOLS$Astrocytes, SecLimaOLS$Neuron,  SecLimaOLS$OPCs , SecLimaOLS$imOLs , SecLimaOLS$maOLs)
  barplot( imaOLS, horiz = T , col= c("darkred" , "red" , "palegreen2","palegreen3","darkgreen"), beside= TRUE , main = TopimaOLS10$Gene.symbol[i] )#, names.arg = c("Astrocytes" ,"Neuron", "OPCs","imOLs" ,"maOLs"))
  i <- i+1
}


### mOLs

selec3 <- subset(dataFiltre2UpDown, select=c(Gene.symbol , Astrocytes, Neuron,  OPCs, imOLs, maOLs))

FCopcs= ifelse(( selec3$maOLs < selec3$OPCs) , (-(0.5)^(log2(selec3$maOLs/selec3$OPCs))) , (2^(log2(selec3$maOLs/selec3$OPCs))))
FCmLOsvsOPCs =round(FCopcs,3)
selec3<- cbind(selec3 , FCmLOsvsOPCs )

FcimaOLS= ifelse((selec3$maOLs < selec3$imOLs) , (-(0.5)^(log2(selec3$maOLs/selec3$imOLs))) , (2^(log2(selec3$maOLs/selec3$imOLs))))
FCmaOLsVSimLOs =round(FcimaOLS,3)
selec3<- cbind(selec3 , FCmaOLsVSimLOs )

FCiOPCsMeanFCimaOLS=apply(selec3[,4:5], 1, mean)
FCmeanOPCsIMaOLS =round(FCiOPCsMeanFCimaOLS,3)
selec3<- cbind(selec3 , FCmeanOPCsIMaOLS )

FcmaOLsvsaLL= (selec3$maOLs/FCmeanOPCsIMaOLS)
FcmaOLsvsAll=round(FcmaOLsvsaLL,3)
selec3<- cbind(selec3 , FcmaOLsvsAll )

write.csv(selec3 , file="maOLS.csv" , row.names = F , quote = F)


### Top markers 
TopmaOLS = arrange(selec3, desc(FcmaOLsvsAll))
TopmaOLS= slice( TopmaOLS , 1:500)
write.csv(TopmaOLS , file="topmaOLS500.csv" , row.names = F , quote = F)
TopmaOLS10= slice( TopmaOLS , 1:10)
write.csv(TopmaOLS10 , file="topmaOLS10.csv" , row.names = F , quote = F)
barplot(TopmaOLS10$FcmaOLsvsAll,horiz=T,names = TopmaOLS10$Gene.symbol,las=1,col="grey")

## enrichement top 10 markers 
par(mfrow=c(1,2))
n = length(TopmaOLS10$Gene.symbol)
for(i in 1:2)
{
  SecLmaOLS<-subset(TopmaOLS10, Gene.symbol ==TopmaOLS10$Gene.symbol[i])
  maOLS = c( SecLmaOLS$Astrocytes, SecLmaOLS$Neuron,  SecLmaOLS$OPCs , SecLmaOLS$imOLs , SecLmaOLS$maOLs)
  barplot( maOLS, horiz = T , col=c("darkred" , "red" , "palegreen2","palegreen3","darkgreen"), beside= TRUE , main = TopmaOLS10$Gene.symbol[i] )#, names.arg = c("Astrocytes" ,"Neuron", "OPCs","imOLs" ,"maOLs"))
  i <- i+1
}


## FC 4 
dataFiltre3Up <- filter(NewdataFiltre1 , FCRealMax >= 4)
dim(dataFiltre3Up) ## 1040 
write.csv(dataFiltre3Up , file= "Up_Fc4_Zhang.csv" , row.names = F , quote = F)
dataFiltre3Down <- filter(NewdataFiltre1 , FCRealMax <= -4)
dim(dataFiltre3Down) #1379
write.csv(dataFiltre3Down , file= "down_Fc4_Zhang.csv" , row.names = F , quote = F)
dataFiltre3UpDown <- filter(NewdataFiltre1 , NewdataFiltre1$FCRealMax>= 4 | NewdataFiltre1$FCRealMax <= -4)
dim(dataFiltre3UpDown) # 2419 
write.csv(dataFiltre3UpDown , file= "Up_down_Fc4_Zhang.csv" , row.names = F , quote = F)

## heatmap 

library(stats)
library(ggplot2)
library(RColorBrewer)
donnes <- (dataFiltre3Up[,1:6])
coul = colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(donnes[,-1]) , col = coul ,  Colv = NA, Rowv = NA  ,  main ="heatmap NA FC:4")

## Fc8 
dataFiltre4Up <- filter(NewdataFiltre1 , FCRealMax >= 8)
dim(dataFiltre4Up) ## 587
write.csv(dataFiltre4Up , file= "Up_Fc8_Zhang.csv" , row.names = F , quote = F)
dataFiltre4Down <- filter(NewdataFiltre1 , FCRealMax <= -8)
dim(dataFiltre4Down) #648 
write.csv(dataFiltre4Down , file= "down_Fc8_Zhang.csv" , row.names = F , quote = F)
dataFiltre4UpDown <- filter(NewdataFiltre1 , NewdataFiltre1$FCRealMax>= 8 | NewdataFiltre1$FCRealMax <= -8)
dim(dataFiltre4UpDown) # 1235 
write.csv(dataFiltre4UpDown , file= "Up_down_FC8_Zhang.csv" , row.names = F , quote = F)


### heatmap 

library(stats)
library(ggplot2)
library(RColorBrewer)
donnes <- (dataFiltre4Up[,1:6])
coul = colorRampPalette(brewer.pal(8, "PiYG"))(25)
heatmap(as.matrix(donnes[,-1]) , col = coul ,  Colv = NA, Rowv = NA  ,  main ="heatmap NA FC:8")


### Neurone - Astrocytes
#### FC1.8 
UpNeuronAstrocytesFc1.8<- filter(NewdataFiltre1 , FCRealMax_neuron >= 1.8)
dim(UpNeuronAstrocytesFc1.8) # 3412 18
write.csv(UpNeuronAstrocytesFc1.8  , file= "up_NA_ZhangFC1.8.csv" , 
            quote=F, col.names= T, row.names=F)

DownNeuronAstrocytesFC1.8<- filter(NewdataFiltre1 , FCRealMax_neuron <=-1.8)
dim(DownNeuronAstrocytesFC1.8) #2576  18
write.csv(DownNeuronAstrocytesFC1.8 , file= "down_NA_Zhang.csv" , 
            quote=F, col.names=T, row.names=F)


#### FC4 
UpNeuronAstrocytesFc4<- filter(NewdataFiltre1 , FCRealMax_neuron >= 4)
dim(UpNeuronAstrocytesFc4) # 1377 18
write.csv(UpNeuronAstrocytesFc4  , file= "up_NA_ZhangFC4.csv" , 
            quote=F, col.names= T, row.names=F)

DownNeuronAstrocytesFC4<- filter(NewdataFiltre1 , FCRealMax_neuron <=-4)
dim(DownNeuronAstrocytesFC4) #1042  18
write.csv(DownNeuronAstrocytesFC4, file= "down_NA_ZhangFC4.csv" ,
            quote=F, col.names=T, row.names=F)

### FC8 
UpNeuronAstrocytesFc8<- filter(NewdataFiltre1 , FCRealMax_neuron >= 8)
dim(UpNeuronAstrocytesFc8) # 646 18
write.csv(UpNeuronAstrocytesFc8  , file= "up_NA_ZhangFC8.csv" ,
            quote=F, col.names= T, row.names=F)

DownNeuronAstrocytesFC8<- filter(NewdataFiltre1 , FCRealMax_neuron <= -8)
dim(DownNeuronAstrocytesFC8) #589  18
write.csv(DownNeuronAstrocytesFC8 , file= "down_NA_ZhangFC8.csv" ,
            quote=F, col.names=T, row.names=F)
