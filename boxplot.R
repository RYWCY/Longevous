profile.all <- read.table("Metabolome_blood_groupfilt.xls",header=T,row.names=1,check.names=F,quote="",sep="\t", na.strings = "", stringsAsFactors=F) #### 821
group       <- as.data.frame(read.table("group_metabolome_blood_clinical.txt",header=T,row.names=1,sep="\t",check.names=F)) ## 382
profile.all <- profile.all[,rownames(group)]

description = read.table("Metabolome_blood.description.xls", header = T, row.names = 1,sep="\t", check.names = F, quote = "")
description = description[rownames(profile.all),]

####################################
cols <- c("#F34800","#00447E")

profile.all.filt= profile.all[grep("Androgenic Steroids",description[,3]),]
write.table(profile.all.filt, file = "Androgenic Steroids.xls", sep = "\t",quote = F)

sort.index = sort(rowMeans(profile.all.filt[,grep("Longevous",group[,2])]),decreasing = F, index.return =T)$ix
profile.all.filt.sort = profile.all.filt[sort.index,]

rownames(profile.all.filt.sort) = gsub("s__","",rownames(profile.all.filt.sort))
profile.ql = profile.all.filt.sort[,grep("Longevous",group[,2])]
profile.qy = profile.all.filt.sort[,grep("Offspring",group[,2])]

pdf(paste("Androgenic_Steroids.boxplot.pdf",sep=''),5,4)

layout(c(1))
par(mar=c(5,20,1,1))
#par(mar=c(3,12,2,2),mgp=c(1.5,0.8,0))

nshow = nrow(profile.ql)
boxplot(log10(t(profile.ql+1e-4)),at=0:(nshow-1)*3+2,las=1,xlab="log10 (Relative abundance +1e-4)",
        col="#F34800",border="#fa6a59",cex.lab=0.9,cex.axis=0.7,horizontal = T,yaxt="n",outcex=0.1)
boxplot(log10(t(profile.qy+1e-4)),at=0:(nshow-1)*3+1,xaxt='n',las=1,
        col="#00447E",border="#5a9ac6",cex.lab=0.9,cex.axis=0.7,horizontal = T,yaxt="n",outcex=0.1,add=T)#,notch=F,)
axis(2,at=0:(nshow-1)*3+1.5,labels=FALSE)
spe=0.3
width=1
text(y=seq(spe+width*1.5, length = nshow, by = 0*spe+3*width), par("usr")[3]-1.5,labels=rownames(profile.ql),cex=0.9,xpd=TRUE,srt=0,adj=1)
legend("bottomright",legend=c("Longevous","Offspring"),col=cols,pch=15,bty="n",cex=0.9)
dev.off()

