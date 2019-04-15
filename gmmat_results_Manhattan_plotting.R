#!/usr/bin/env Rscript
#Date:Mar 18 2019
#Author:Kiranmayee Bakshy

# A program to make Manhattan plot of GWAS results obtained using GMMAT package
# This program also outputs the significant SNPs and also IGC marker selections

options(bitmapType='cairo')

library(ggplot2)
library(dplyr)
library(qqman)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

# test if there are at least 4 argument: if not, return an error
if (length(args)<4) {
  stop("At least 4 arguments must be supplied input gwas results.\n", call.=FALSE)
} 


gwasResults1<-read.delim(args[1], sep="\t", header=F, stringsAsFactors=F)
gwasResults2<-read.delim(args[2], sep="\t", header=F, stringsAsFactors=F)
gwasResults3<-read.delim(args[3], sep="\t", header=F, stringsAsFactors=F)
gwasResults4<-read.delim(args[4], sep="\t", header=F, stringsAsFactors=F)

colnames(gwasResults1)<-c("SNP", "CHR", "BP", "ref", "alt", "n", "af", "beta", "se", "pval", "converged")
colnames(gwasResults2)<-c("SNP", "CHR", "BP", "ref", "alt", "n", "af", "beta", "se", "pval", "converged")
colnames(gwasResults3)<-c("SNP", "CHR", "BP", "ref", "alt", "n", "af", "beta", "se", "pval", "converged")
colnames(gwasResults4)<-c("SNP", "CHR", "BP", "ref", "alt", "n", "af", "beta", "se", "pval", "converged")


#SnpsOfInterest
soi1 <- c("ARS_PIRBRIGHT_5_98985645","ARS_PIRBRIGHT_5_99036250","ARS_PIRBRIGHT_5_99333595","ARS_PIRBRIGHT_5_99445508","ARS_PIRBRIGHT_5_99780983","ARS_PIRBRIGHT_18_62460291","ARS_PIRBRIGHT_18_62559417","ARS_PIRBRIGHT_18_62644240","ARS_PIRBRIGHT_18_62670367","ARS_PIRBRIGHT_18_62774779","ARS_PIRBRIGHT_18_62812297","ARS_PIRBRIGHT_18_63089154","ARS_PIRBRIGHT_23_28535354","ARS_PIRBRIGHT_23_28651088","ARS_PIRBRIGHT_CH240_391K10_KIR_41480","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_61352","ARS_PIRBRIGHT_LIB14413_LRC_65014","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_72198","ARS_PIRBRIGHT_LIB14413_LRC_81028","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_98785","ARS_PIRBRIGHT_LIB14413_LRC_106729","ARS_PIRBRIGHT_CH240_370M3_LILR_LRC_174904","ARS_PIRBRIGHT_LIB14427_MHC_9213","ARS_PIRBRIGHT_LIB14427_MHC_43656","ARS_PIRBRIGHT_LIB14427_MHC_59013","ARS_PIRBRIGHT_LIB14427_MHC_86084","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_115082","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_143922","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_154399","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_208321","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_260913","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_286137","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_317666","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_324231","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_380177","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_400205")
soi <- c("ARS_PIRBRIGHT_5_99094858","ARS_PIRBRIGHT_5_99190989","ARS_PIRBRIGHT_5_99194167","ARS_PIRBRIGHT_5_99203402","ARS_PIRBRIGHT_5_99412564","ARS_PIRBRIGHT_5_99437819","ARS_PIRBRIGHT_5_99749976","ARS_PIRBRIGHT_18_62572950","ARS_PIRBRIGHT_18_62665698","ARS_PIRBRIGHT_18_62716825","ARS_PIRBRIGHT_18_62766196","ARS_PIRBRIGHT_18_62994372","ARS_PIRBRIGHT_18_63020431","ARS_PIRBRIGHT_18_63036451","ARS_PIRBRIGHT_18_63067935","ARS_PIRBRIGHT_18_63082203","ARS_PIRBRIGHT_18_63084493","ARS_PIRBRIGHT_18_63114542","ARS_PIRBRIGHT_18_63122963","ARS_PIRBRIGHT_18_63131111","ARS_PIRBRIGHT_18_63141688","ARS_PIRBRIGHT_18_63156196","ARS_PIRBRIGHT_18_63185710","ARS_PIRBRIGHT_18_63186702","ARS_PIRBRIGHT_18_63284810","ARS_PIRBRIGHT_18_63294645","ARS_PIRBRIGHT_18_63340289","ARS_PIRBRIGHT_18_63385249","ARS_PIRBRIGHT_18_63399311","ARS_PIRBRIGHT_18_63417698","ARS_PIRBRIGHT_23_28648192","ARS_PIRBRIGHT_LIB14427_MHC_2500","ARS_PIRBRIGHT_LIB14427_MHC_3538","ARS_PIRBRIGHT_LIB14427_MHC_36271","ARS_PIRBRIGHT_LIB14427_MHC_45690","ARS_PIRBRIGHT_LIB14427_MHC_57505","ARS_PIRBRIGHT_LIB14427_MHC_73766","ARS_PIRBRIGHT_LIB14427_MHC_116810","ARS_PIRBRIGHT_LIB14427_MHC_122678","ARS_PIRBRIGHT_LIB14427_MHC_126886","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_120784","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_181938","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_302664","ARS_PIRBRIGHT_TPI4222_A14_MHCclassI_MHC_395846")


# Prepare the datasets
gwasResults<-list(gwasResults1,gwasResults2,gwasResults3,gwasResults4)

don1 <- gwasResults1 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwasResults1, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%

  # Add highlight and annotation information
  mutate( is_highlight1=ifelse(SNP %in% soi1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(SNP %in% soi, "yes", "no")) 

# Prepare X axis
axisdf1 <- don1 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


don2 <- gwasResults2 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(gwasResults2, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  mutate( is_highlight1=ifelse(SNP %in% soi1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(SNP %in% soi, "yes", "no")) 
  axisdf2 <- don2 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don3 <- gwasResults3 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(gwasResults3, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  mutate( is_highlight1=ifelse(SNP %in% soi1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(SNP %in% soi, "yes", "no")) 
  axisdf3 <- don3 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

don4 <- gwasResults4 %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(gwasResults4, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  mutate( is_highlight1=ifelse(SNP %in% soi1, "yes", "no")) %>%
  mutate( is_highlight2=ifelse(SNP %in% soi, "yes", "no")) 
  axisdf4 <- don4 %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )



# Make the plots
p1<-ggplot(don1, aes(x=BPcum, y=-log10(pval))) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +

    # custom X axis:
    scale_x_continuous( label = axisdf1$CHR, breaks= axisdf1$center ) +
    scale_y_continuous(limits = c(0,15), expand = c(0, 0) ) +     # remove space between plot area and x axis
   
    # Add highlighted points
    geom_point(data=subset(don1, is_highlight1=="yes"), color="red", size=2) +

    # Add second set of highlighted points
    geom_point(data=subset(don1, is_highlight2=="yes"), color="green", size=2) +

    # Add suggestive and genome wide significance threshold lines
    geom_hline(yintercept=-log10(5.3e-06), linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(2.7e-07), linetype="dashed", color = "blue") +

   # Add title and axis labels
   ggtitle("PhenoGrp1") + xlab("Chromosome") + ylab("-log10(p)") +

   # Change angle of X-axis text 
   theme(axis.text.x = element_text(angle=45)) +
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

p2<-ggplot(don2, aes(x=BPcum, y=-log10(pval))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous( label = axisdf2$CHR, breaks= axisdf2$center ) +
    scale_y_continuous(limits = c(0,15), expand = c(0, 0) ) +     # remove space between plot area and x axis
    geom_point(data=subset(don2, is_highlight1=="yes"), color="red", size=2) +
    geom_point(data=subset(don2, is_highlight2=="yes"), color="green", size=2) +
    geom_hline(yintercept=-log10(5.3e-06), linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(2.7e-07), linetype="dashed", color = "blue") +
    ggtitle("PhenoGrp2") + xlab("Chromosome") + ylab("-log10(p)") +
    theme(axis.text.x = element_text(angle=45)) +
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

p3<-ggplot(don3, aes(x=BPcum, y=-log10(pval))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous( label = axisdf3$CHR, breaks= axisdf3$center ) +
    scale_y_continuous(limits = c(0,15), expand = c(0, 0) ) +     # remove space between plot area and x axis
    geom_point(data=subset(don3, is_highlight1=="yes"), color="red", size=2) +
    geom_point(data=subset(don3, is_highlight2=="yes"), color="green", size=2) +
    geom_hline(yintercept=-log10(5.3e-06), linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(2.7e-07), linetype="dashed", color = "blue") +
    ggtitle("PhenoGrp3") + xlab("Chromosome") + ylab("-log10(p)") +
    theme(axis.text.x = element_text(angle=45)) +
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

p4<-ggplot(don4, aes(x=BPcum, y=-log10(pval))) +
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
    scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
    scale_x_continuous( label = axisdf4$CHR, breaks= axisdf4$center ) +
    scale_y_continuous(limits = c(0,15), expand = c(0, 0) ) +     # remove space between plot area and x axis
    geom_point(data=subset(don4, is_highlight1=="yes"), color="red", size=2) +
    geom_point(data=subset(don4, is_highlight2=="yes"), color="green", size=2) +
    geom_hline(yintercept=-log10(5.3e-06), linetype="dashed", color = "black") +
    geom_hline(yintercept=-log10(2.7e-07), linetype="dashed", color = "blue") +
    ggtitle("PhenoGrp4") + xlab("Chromosome") + ylab("-log10(p)") +
    theme(axis.text.x = element_text(angle=45)) +
    theme_bw() +
    theme( 
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )

g<-grid.arrange(p1,p2,p3,p4, ncol = 2)
ggsave("Manhattan.pdf", g, height=7, width=14)

pdf("QQplot.pdf", w=14, h=7)
par(mfrow=c(2,2))
qq(gwasResults1$pval, main = "PhenoGrp1", xlim = c(0, 7), ylim = c(0, 21), pch = 1, col = "blue4")
text(2,10,expression(paste(lambda, "=1.001")))
qq(gwasResults2$pval, main = "PhenoGrp2", xlim = c(0, 7), ylim = c(0, 21), pch = 1, col = "blue4")
text(2,10,expression(paste(lambda, "=1.013")))
qq(gwasResults3$pval, main = "PhenoGrp3", xlim = c(0, 7), ylim = c(0, 21), pch = 1, col = "blue4")
text(2,10,expression(paste(lambda, "=1.010")))
qq(gwasResults4$pval, main = "PhenoGrp4", xlim = c(0, 7), ylim = c(0, 21), pch = 1, col = "blue4")
text(2,10,expression(paste(lambda, "=1.014")))
dev.off()


