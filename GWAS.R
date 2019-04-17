#!/bin/Rscript
# Date: Nov 30 2018
# This script runs GWAS using GMMAT package on the list of snps that the user inputs
# GDS file can be generated from the binary plink files, for example: SeqArray::seqBED2GDS("file.bed","file.fam","file.bim", "file.gds") 
# GRM can be generated from the binary plink files using the GEMMA package but it should have the column and row names same as the ID column of the phenotype file
# Phenotype file should have the ID column in the current GMMAT version  1.0.3
# for more information please visit this page "https://github.com/hanchenphd/GMMAT/blob/master/inst/doc/GMMAT.pdf"

library(GMMAT)

setwd("/beegfs/project/rumen_longread_metagenome_assembly/kiranmayee/IGC/prelim_gwas/round2/gmmat")

args = commandArgs(trailingOnly=TRUE)

# test if there are at least 5 argument: if not, return an error
if (length(args)<5) {
  stop("At least 5 arguments must be supplied gds file, phenotype file, GRM, snp.ids and output file name.\n", call.=FALSE)
} 


infile<-args[1]

pheno<-read.table(args[2], header=T, stringsAsFactors=F, sep="\t")

GRM<-as.matrix(read.table(args[3]))

colnames(GRM)<-pheno$id

row.names(GRM)<-pheno$id

snps <- scan(args[4], what="character")

results<-glmm.wald(fixed = bTB_status ~ breed + age + yr + season + reason + prevalence, data = pheno, kins = GRM, id = "id", family = binomial(link = "logit"), infile = infile, snps = snps)

write.table(results, file=args[5], append=T, row.names=F, col.names=F, quote=F, sep = "\t")
