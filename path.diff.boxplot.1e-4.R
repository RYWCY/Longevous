
profile.all <- read.csv("L3.profile",header=T,row.names=1,check.names=F,quote="",sep="\t") 

group       <- as.data.frame(read.table("group.txt",header=F,row.names=1,sep="\t",check.names=F))
colnames(group) <- 'group'
profile.all <- profile.all[,rownames(group)]

profile.offspring = profile.all[,grep("Offspring",group[,1])]
profile.Longevous = profile.all[,grep("Longevous",group[,1])]

mean.offspring = rowMeans(profile.offspring)
mean.Longevous = rowMeans(profile.Longevous)

###########Longevous ??offspring֮??????????
pvalue.1=matrix(0,nrow(profile.offspring),1)
for (i in 1:nrow(profile.offspring)){
  pvalue.1[i,1]=wilcox.test(unlist(profile.offspring[i,]),unlist(profile.Longevous[i,]))$p.value
}

fdr.1 = p.adjust(pvalue.1, n = length(pvalue.1), method = "fdr")


##########
out = cbind(mean.offspring, mean.Longevous,
            pvalue.1,fdr.1)
colnames(out) = c("mean abundance of offspring", "mean abundance of Longevous",
                  "Pvalue of comparison between offspring and Longevous",
                  "FDR of comparison between offspring and Longevous")
write.table(out, file = "pathway.compare.xls", sep = "\t")

###### increase, Longevous ????offspring???????ӵ?????fdr<0.05??
increase.index = Reduce(intersect,list(A=as.data.frame(which (mean.Longevous > mean.offspring))[,1], 
                                       #                                       B=as.data.frame(which (mean.Longevous > mean.inLongevous))[,1], 
                                       C=which (fdr.1<0.05))) 
length(increase.index)

increase.out = cbind(mean.offspring[increase.index], mean.Longevous[increase.index],
                     pvalue.1[increase.index],fdr.1[increase.index])
colnames(increase.out) = c("mean abundance of offspring", "mean abundance of Longevous",
                           "Pvalue of comparison between offspring and Longevous",
                           "FDR of comparison between offspring and Longevous")
write.table(increase.out, file = "pathway.increase_in_Longevous.xls", sep = "\t")

###### decrease, Longevous????offspring???????ٵ?????fdr<0.05??
decrease.index = Reduce(intersect,list(A=as.data.frame(which (mean.Longevous < mean.offspring))[,1], 
                                       #                                       B=as.data.frame(which (mean.Longevous < mean.inLongevous))[,1], 
                                       C=which (fdr.1<0.05)))
length(decrease.index)

decrease.out = cbind(mean.offspring[decrease.index], mean.Longevous[decrease.index], 
                     pvalue.1[decrease.index],fdr.1[decrease.index])
colnames(decrease.out) = c("mean abundance of offspring", "mean abundance of Longevous", 
                           "Pvalue of comparison between offspring and Longevous",
                           "FDR of comparison between offspring and Longevous")
write.table(decrease.out, file = "pathway.decrease_in_Longevous.xls", sep = "\t")


########################################################################################
# diff ggplot boxplot
# https://cmdlinetips.com/2019/02/how-to-make-grouped-boxplots-with-ggplot2/
# http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually ????????????
library(tidyverse)
library(gapminder)
#cols = c("#006d2c", "#fd8d3c", "#a50f15", "#33A02C","#E31A1C")
cols <- c("#33A02C","#E31A1C","#1F78B4","#FB9A99","#FF7F00","#FDBF6F","#6A3D9A","#CAB2D6","#B15928","#FFFF99","#DA1C5C","#F8D5E1","#35978F","#C7EAE5","#8C510A","#F6E8C3","#8E0152","#BC6898","#737373","#A6CEE3")

#diff.index = Reduce(intersect,list(A=which (fdr.1<0.05)))
diff.index = Reduce(intersect,list(A=which (fdr.1<0.05),B=which(mean.Longevous>1e-4),C=which(mean.offspring>1e-4)))
write.table(profile.all[diff.index,], file = "pathway.diff.xls", sep = "\t",quote = F)

#ko = read.table("map_title.tab",sep="\t",head=F,check.names = F, row.names=2, quote="")
increase.index = Reduce(intersect,list(A=increase.index,B=which(mean.Longevous>1e-4),C=which(mean.offspring>1e-4)))
profile.all.filt.increase = profile.all[increase.index,]
sort.index.increase = sort(rowMeans(profile.all.filt.increase[,grep("Longevous",group[,1])]),decreasing = T, index.return =T)$ix
profile.all.filt.increase.sort = profile.all.filt.increase[sort.index.increase,]

decrease.index = Reduce(intersect,list(A=decrease.index,B=which(mean.Longevous>1e-4),C=which(mean.offspring>1e-4)))
profile.all.filt.decrease = profile.all[decrease.index,]
sort.index.decrease = sort(rowMeans(profile.all.filt.decrease[,grep("Offspring",group[,1])]),decreasing = T, index.return =T)$ix
profile.all.filt.decrease.sort = profile.all.filt.decrease[sort.index.decrease,]

profile.all.filt.diff = rbind(profile.all.filt.increase.sort,profile.all.filt.decrease.sort)
pathway.sort = rownames(profile.all.filt.diff)
profile.all.filt.diff = cbind(pathway = pathway.sort, profile.all.filt.diff)

library("reshape2")
data = melt(profile.all.filt.diff, id = c("pathway"))
data.group = c(replicate(length(diff.index)*116,"Longevous"),replicate(length(diff.index)*232,"Offspring"))
data = as.data.frame(cbind(data,group=data.group),stringAsFactor=F)
data$group = factor(data$group, levels = c("Longevous","Offspring"))
data$pathway = factor(data$pathway, levels = pathway.sort)

library(ggplot2)

pdf(paste("pathway.diff.1e-4.pdf",sep=''),8,7)
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

