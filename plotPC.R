#!/usr/bin/env Rscript
# Date: Mar 4 2019
# This script will plot the first 2 PCs from the plink pca output 

setwd("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gmmat/pca")

options(bitmapType='cairo')


args = commandArgs(trailingOnly=TRUE)

# test if there are at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one arguments must be supplied (plink.eigenvec pca file).\n", call.=FALSE)
} 

pheno4<-read.delim(args[1], sep=" ", header=T, stringsAsFactors=F)
pheno1<-read.delim(args[2], sep=" ", header=T, stringsAsFactors=F)
pheno2<-read.delim(args[3], sep=" ", header=T, stringsAsFactors=F)
pheno3<-read.delim(args[4], sep=" ", header=T, stringsAsFactors=F)

pheno4.fam<-read.delim("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/data_files/merge1.fam", sep="\t")
pheno1.fam<-read.delim("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/data_files/pheno1.fam", sep=" ")
pheno2.fam<-read.delim("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/data_files/pheno2.fam", sep=" ")
pheno3.fam<-read.delim("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/data_files/pheno3.fam", sep=" ")

colnames(pheno4.fam)<-c("FID","IID","PID","MID","SEX","bTB_status")
colnames(pheno1.fam)<-c("FID","IID","PID","MID","SEX","bTB_status")
colnames(pheno2.fam)<-c("FID","IID","PID","MID","SEX","bTB_status")
colnames(pheno3.fam)<-c("FID","IID","PID","MID","SEX","bTB_status")

df1<-merge(pheno4,pheno4.fam,by="IID")
df2<-merge(pheno1,pheno1.fam,by="IID")
df3<-merge(pheno2,pheno2.fam,by="IID")
df4<-merge(pheno3,pheno3.fam,by="IID")

df1<-df1[,c(1,3,4,27)]
df2<-df2[,c(1,3,4,27)]
df3<-df3[,c(1,3,4,27)]
df4<-df4[,c(1,3,4,27)]

df1$status<-ifelse(df1$bTB_status == 0, "Control", "Case")
df2$status<-ifelse(df2$bTB_status == 0, "Control", "Case")
df3$status<-ifelse(df3$bTB_status == 0, "Control", "Case")
df4$status<-ifelse(df4$bTB_status == 0, "Control", "Case")

library(ggplot2)
library(gridExtra)

p1<-ggplot(df1, aes(x=PC1, y=PC2, col=as.factor(status))) + geom_point() + ggtitle("PhenoGrp4") + xlab("PC1 (22.8%)") + ylab("PC2 (17.3%)") + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                                       panel.grid.minor = element_blank()) + scale_color_manual("bTB status",values = c("Control" = "blue", "Case" = "red")) + theme(legend.justification=c(1,0), 
                                                                                                                                                                                            legend.position=c(0.95, 0.75),
                                                                                                                                                                                            legend.background = element_blank(),
                                                                                                                                                                                            legend.key = element_blank())


p2<-ggplot(df2, aes(x=PC1, y=PC2, col=as.factor(status))) + geom_point() + ggtitle("PhenoGrp1") + xlab("PC1 (14.7%)") + ylab("PC2 (12.4%)") + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                                   panel.grid.minor = element_blank()) + scale_color_manual("bTB status",values = c("Control" = "blue", "Case" = "red")) + theme(legend.justification=c(1,0), 
                                                                                                                                                                                                                                                 legend.position=c(0.95, 0.75),
                                                                                                                                                                                                                                                 legend.background = element_blank(),
                                                                                                                                                                                                                                                 legend.key = element_blank())


p3<-ggplot(df3, aes(x=PC1, y=PC2, col=as.factor(status))) + geom_point() + ggtitle("PhenoGrp2") + xlab("PC1 (20.5%)") + ylab("PC2 (15.1%)")  + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                                     panel.grid.minor = element_blank()) + scale_color_manual("bTB status",values = c("Control" = "blue", "Case" = "red")) + theme(legend.justification=c(1,0), 
                                                                                                                                                                                                                                                   legend.position=c(0.95, 0.75),
                                                                                                                                                                                                                                                   legend.background = element_blank(),
                                                                                                                                                                                                                                                   legend.key = element_blank())


p4<-ggplot(df4, aes(x=PC1, y=PC2, col=as.factor(status))) + geom_point() + ggtitle("PhenoGrp3") + xlab("PC1 (20.8%)") + ylab("PC2 (15.3%)") + theme_bw() + theme(panel.grid.major = element_blank(), 
                                                                                                                     panel.grid.minor = element_blank()) + scale_color_manual("bTB status",values = c("Control" = "blue", "Case" = "red")) + theme(legend.justification=c(1,0), 
                                                                                                                                                                                                                                                   legend.position=c(0.95, 0.75),
                                                                                                                                                                                                                                                   legend.background = element_blank(),
                                                                                                                                                                                                                                                   legend.key = element_blank())


g<-grid.arrange(p2, p3, p4, p1, ncol=2)
ggsave("allPhenoGrps.pdf", plot=g, device = pdf(), scale = 1, width = 11, height = 8.5, units = "cm", dpi = 600)  