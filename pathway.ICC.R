library("stringr")

profile.all <- read.csv("L3.profile",header=T,row.names=1,check.names=F,quote="",sep="\t") 
#rownames(profile.all)=str_extract(rownames(profile.all),"s__[\\w\\W]+")
group       <- as.data.frame(read.table("group.txt",header=F,row.names=1,sep="\t",check.names=F))
colnames(group) <- 'group'
profile.all <- profile.all[,rownames(group)]

profile.all.filt = profile.all
profile.all.filt = as.data.frame(profile.all[which (rowMeans(profile.all) > 1e-4),])  #### 524
#mean(colSums(as.data.frame(profile.all[which (rowMeans(profile.all) > 1e-4),]))/colSums(profile.all)) #### 97.86%

profile.offspring = profile.all.filt[,grep("Offspring",group[,1])]
profile.Longevous = profile.all.filt[,grep("Longevous",group[,1])]

mean.offspring = rowMeans(profile.offspring)
mean.Longevous = rowMeans(profile.Longevous)

library("multilevel")
ICC1=matrix(0,nrow=nrow(profile.all.filt),ncol=1)
rownames(ICC1)=rownames(profile.all.filt)
colnames(ICC1)="ICC1"
ICC2=matrix(0,nrow=nrow(profile.all.filt),ncol=1)
rownames(ICC2)=rownames(profile.all.filt)
colnames(ICC2)="ICC2"
for (i in 1:nrow(profile.all.filt)){
  hrs.mod <- aov(unlist(profile.all.filt[i,]) ~ as.factor(group[,1]))
  ICC1[i,1]=ICC1(hrs.mod)
  ICC2[i,1]=ICC2(hrs.mod)
}

mean.diff=abs(mean.Longevous - mean.offspring)/(mean.offspring+mean.Longevous)
log2=abs(log2(mean.Longevous/mean.offspring))

###########Longevous ??offspring֮??????????
pvalue.1=matrix(0,nrow(profile.offspring),1)
for (i in 1:nrow(profile.offspring)){
  pvalue.1[i,1]=wilcox.test(unlist(profile.offspring[i,]),unlist(profile.Longevous[i,]))$p.value
}

fdr.1 = p.adjust(pvalue.1, n = length(pvalue.1), method = "fdr")

######
out = cbind(mean.offspring, mean.Longevous, pvalue.1,fdr.1,mean.diff,log2,ICC1,ICC2)
colnames(out) = c("mean abundance of offspring", "mean abundance of Longevous", 
                  "Pvalue of comparison between offspring and Longevous",
                  "FDR of comparison between offspring and Longevous",
                  "mean.diff","log2(mean.Longevous/mean.offspring","ICC1","ICC2")
write.table(out, file = "pathway.compare.1e-4.ICC.xls", sep = "\t")


########################################################################################
# diff ggplot boxplot
# https://cmdlinetips.com/2019/02/how-to-make-grouped-boxplots-with-ggplot2/
# http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually ????????????
library(tidyverse)
library(gapminder)
#cols = c("#006d2c", "#fd8d3c", "#a50f15", "#33A02C","#E31A1C")
cols <- c("#33A02C","#E31A1C","#1F78B4","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99","#DA1C5C","#F8D5E1","#35978F","#C7EAE5","#8C510A","#F6E8C3","#8E0152","#BC6898","#737373","#A6CEE3")

index = Reduce(intersect,list(A=which (fdr.1>0.1),B=which (ICC1<0)))
write.table(profile.all[index,], file = "pathway.nodiff.xls", sep = "\t",quote = F)
#ko = read.table("map_title.tab",sep="\t",head=F,check.names = F, row.names=2, quote="")
profile.all.filt = profile.all.filt[index,]
sort.index = sort(rowMeans(profile.all.filt[,grep("Longevous",group[,1])]),decreasing = T, index.return =T)$ix
sort.index = sort.index[1:20]
profile.all.filt.sort = profile.all.filt[sort.index,]

write.table(profile.all.filt.sort, file = "pathway.nodiff.forplot.xls", sep = "\t",quote = F)

pathway.sort = rownames(profile.all.filt.sort)
profile.all.filt.nodiff = cbind(pathway = pathway.sort, profile.all.filt.sort)


library("reshape2")
data = melt(profile.all.filt.nodiff, id = c("pathway"))
data.group = c(replicate(length(sort.index)*116,"Longevous"),replicate(length(sort.index)*232,"Offspring"))
data = as.data.frame(cbind(data,group=data.group),stringAsFactor=F)
data$group = factor(data$group, levels = c("Longevous","Offspring"))
data$pathway = factor(data$pathway, levels = pathway.sort)

library(ggplot2)

pdf(paste("pathway.nodiff.ICC1_less_0.fdr_more_0.1.pdf",sep=''),8,6.5)
cols <- c("#E31A1C","#1F78B4","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99","#DA1C5C","#F8D5E1","#35978F","#C7EAE5","#8C510A","#F6E8C3","#8E0152","#BC6898","#737373","#A6CEE3")

ggplot(data,aes(x=pathway, y=value, fill=factor(group))) +
  geom_boxplot() + scale_fill_manual(values = cols) +
  labs(fill = "Pathway") + 
  ylab("Abundance") + 
  xlab("") +
  #  geom_point(position=position_jitterdodge(),alpha=0.3) +
  theme_bw(base_size = 16)+
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

dev.off()



