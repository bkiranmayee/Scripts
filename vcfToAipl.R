#!/usr/bin/env Rscript
#Date:Sep 21 2018
#Author:Kiranmayee Bakshy

# A program to convert the plink genotype recodings (A-transpose) output to AIPL format ped and map files

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Usage: vcfToAipl.R plink.traw output filename.\n.At least one argument must be supplied (plink.traw).\n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.aipl"
}


dat <- read.delim(args[1], header=T, stringsAsFactors=FALSE)
dat.map<-dat[,c(1:6)]
names(dat.map)<-c("CHR","SNP","Dist","POS","REF","ALT")
write.table(x = dat.map[,c(1:4)], file = paste(args[2],"map", sep="."), sep=" ", row.names = F, col.names = F, quote=F)
samples<-names(dat[,c(7:ncol(dat))])
samples.list<-gsub("X0_","",samples)
names(dat[,c(7:ncol(dat))])<-samples.list
dat.ped<-t(dat[,c(7:ncol(dat))])
dat1.ped<-as.data.frame(cbind(samples.list,dat.ped),stringsAsFactors = F)
row.names(dat1.ped)<-samples.list
dat.ped<-dat1.ped
dat.ped[2:ncol(dat.ped)]<-lapply(dat1.ped[2:ncol(dat1.ped)], function(x) ifelse(x == 0, 2, ifelse(x == 2, 0, x)))
#dat.ped[2:ncol(dat.ped)]<-lapply(dat.ped[2:ncol(dat.ped)], function(x) ifelse(x == "NA", 5, x))
write.table(dat.ped, file = paste(args[2],"ped", sep="."), sep=" ",col.names=F, quote=F, NA = "5")

