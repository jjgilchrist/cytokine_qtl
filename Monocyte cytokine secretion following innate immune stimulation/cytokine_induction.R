rm(list = ls(all=TRUE))
library(ppcor)
library(patchwork)
library(ggplot2)
library(stringr)
library(RColorBrewer)

#cytokine induction
###1st normalize raw data###
###raw data is background (negative control) subtracted###

#read in raw cytokine quantification in supernatants
raw.data <- read.table("cytokine_data.raw.txt", row.names = 1, header = T)

#read in covariate data
covs <- read.table("covar.txt", row.names = 3,header = T)

#read in genotype principal components
geno.pcs <- read.table("geno_pcs.txt", row.names = 1, header = T)

#read in cytokine exclusions
cyto.rm <- read.table("cyto_rm.txt", header = F)

#reformat raw data
cytokines <- NA
for (i in c(1:200)){
  cytokines[i] <- str_split(colnames(raw.data)[i], "_")[[1]][2]
}

condition.index <- c(rep("IFN", 345), rep("LPS2", 345), rep("LPS24",345), rep("UT", 345))
colnames(raw.data) <- cytokines

raw.data2 <- rbind(raw.data[,c(1:50)], raw.data[,c(51:100)], raw.data[,c(101:150)], raw.data[,c(151:200)])

row.names(raw.data2) <- c(paste0(rownames(raw.data), "IFN"),paste0(rownames(raw.data), "LPS2"),paste0(rownames(raw.data), "LPS24"),paste0(rownames(raw.data), "UT"))

#remove cytokines not passing QC
raw.data.qc <- raw.data2[,-which(colnames(raw.data2) %in% c(cyto.rm[,1]))]

#normalise across all samples (for cytokine induction analysis)
raw.data.norm <- matrix(NA,ncol=28,nrow=1380)

for (i in c(1:28)){
  keep <- !is.na(raw.data.qc[,i])
  raw.data.norm[keep,i] <- qnorm((0.5 + rank(raw.data.qc[keep,i]))/(1 + sum(!is.na(raw.data.qc[keep,i]))))
}

#then regress out covariates
raw.data.norm2 <- raw.data.norm*NA

all(row.names(raw.data2)==row.names(covs))

for (i in 1:28){
  keep <- !is.na(raw.data.norm[,i]) & !is.na(covs$age) & !is.na(covs$sex) & !is.na(covs$batch) & !is.na(covs$monocyte)
  temp <- predict(lm(raw.data.norm[keep,i] ~  covs$age[keep] + covs$sex[keep] + covs$batch[keep] + covs$monocyte[keep]))
  raw.data.norm2[keep,i] <- raw.data.norm[keep,i]  - temp
}

raw.data.norm2 <- data.frame(raw.data.norm2)
colnames(raw.data.norm2) <- colnames(raw.data.qc)
rownames(raw.data.norm2) <- rownames(raw.data2)

#cytokine data normalised across samples

cyt.ind.data <- raw.data.norm2

#calculate cytokine induction for each of three stimulations - collect p-values

p.ifn <- c()
p.lps2 <- c()
p.lps24 <- c()
for (i in c(1:28)){
  p.ifn[i] <- summary(lm(cyt.ind.data[c(1:345, 1036:1380),i]~c(rep(1,345), rep(0,345))))$coef[2,4]
  p.lps2[i] <- summary(lm(cyt.ind.data[c(346:690, 1036:1380),i]~c(rep(1,345), rep(0,345))))$coef[2,4]
  p.lps24[i] <- summary(lm(cyt.ind.data[c(691:1035, 1036:1380),i]~c(rep(1,345), rep(0,345))))$coef[2,4]
  
}

#calculate fold-change cytokine induction for each of three stimulations - using raw data
raw.data <- raw.data.qc

#set negative values to zero
raw.data[raw.data<0] <- 0

fc.ifn <- c()
fc.lps2 <- c()
fc.lps24 <- c()
for (i in c(1:28)){
  fc.ifn[i] <- (mean(na.omit(raw.data[c(1:345),i]))+1)/(mean(na.omit(raw.data[c(1036:1380),i]))+1)
  fc.lps2[i] <- (mean(na.omit(raw.data[c(346:690),i]))+1)/(mean(na.omit(raw.data[c(1036:1380),i]))+1)
  fc.lps24[i] <- (mean(na.omit(raw.data[c(691:1035),i]))+1)/(mean(na.omit(raw.data[c(1036:1380),i]))+1)
}

#compare RNA induction with cytokine secretion
#read in RNA expression data (log2 transformed) for 28 post-qc cytokines

naive.rna.28 <- read.table("naive_rna_cyto.txt", header = T, row.names = 1)
lps2.rna.28 <- read.table("lps2_rna_cyto.txt", header = T, row.names = 1)
lps24.rna.28 <- read.table("lps24_rna_cyto.txt", header = T, row.names = 1)
ifn.rna.28 <- read.table("ifn_rna_cyto.txt", header = T, row.names = 1)

#calculate fold-change for 3 stimulation conditions c.f. baseline (naive)

fc.ifn.rna <- c()
fc.lps2.rna <- c()
fc.lps24.rna <- c()

for (i in c(1:28)){
  fc.ifn.rna[i]  <- 2^(mean(na.omit(ifn.rna.28[,i]))-mean(na.omit(naive.rna.28[,i])))
  fc.lps2.rna[i]  <- 2^(mean(na.omit(lps2.rna.28[,i]))-mean(na.omit(naive.rna.28[,i])))
  fc.lps24.rna[i]  <- 2^(mean(na.omit(lps24.rna.28[,i]))-mean(na.omit(naive.rna.28[,i])))
}

#calculate log2 fold-change for 3 stimulation conditions c.f. baseline (naive)

log.fc.ifn.rna <- c()
log.fc.lps2.rna <- c()
log.fc.lps24.rna <- c()

for (i in c(1:28)){
  log.fc.ifn.rna[i]  <- mean(na.omit(ifn.rna.28[,i]))-mean(na.omit(naive.rna.28[,i]))
  log.fc.lps2.rna[i]  <- mean(na.omit(lps2.rna.28[,i]))-mean(na.omit(naive.rna.28[,i]))
  log.fc.lps24.rna[i]  <- mean(na.omit(lps24.rna.28[,i]))-mean(na.omit(naive.rna.28[,i]))
}

#calculate p-values for 3 stimulation conditions c.f. baseline (naive)

p.ifn.rna <- c()
p.lps2.rna <- c()
p.lps24.rna <- c()

for (i in c(1:28)){
  p.ifn.rna[i]  <- summary(lm(c(ifn.rna.28[,i], naive.rna.28[,i])~c(rep(1,367), rep(0,414))))$coef[2,4]
  p.lps2.rna[i]  <- summary(lm(c(lps2.rna.28[,i], naive.rna.28[,i])~c(rep(1,261), rep(0,414))))$coef[2,4]
  p.lps24.rna[i]  <- summary(lm(c(lps24.rna.28[,i], naive.rna.28[,i])~c(rep(1,322), rep(0,414))))$coef[2,4]
}

#plot out cytokine induction (Fig S1)
#calculate induction FDRs - baseline vs 3 stimulations (for RNA induction)

fdr.lps2.rna <- p.adjust(p.lps2.rna, method = "fdr")
fdr.lps24.rna <- p.adjust(p.lps24.rna, method = "fdr")
fdr.ifn.rna <- p.adjust(p.ifn.rna, method = "fdr")

#calculate induction FDRs - baseline vs 3 stimulations (for cytokine protein secretion)

fdr.lps2 <- p.adjust(p.lps2, method = "fdr")
fdr.lps24 <- p.adjust(p.lps24, method = "fdr")
fdr.ifn <- p.adjust(p.ifn, method = "fdr")

#time and vector of ordered fold-changes
time <- rep(rep(c(0,12, 12,24),4), 28)
fc <- c(1,fc.lps2[1], fc.lps2[1], fc.lps24[1],1,NA, NA, fc.ifn[1], 1,fc.lps2.rna[1], fc.lps2.rna[1], fc.lps24.rna[1],1,NA, NA, fc.ifn.rna[1])

for (i in c(2:28)){
  fc <- c(fc, c(1,fc.lps2[i], fc.lps2[i], fc.lps24[i],1,NA, NA, fc.ifn[i], 1,fc.lps2.rna[i], fc.lps2.rna[i], fc.lps24.rna[i],1,NA, NA, fc.ifn.rna[i]))
}


stim <- rep(c("lps","lps","lps","lps","ifn","ifn","ifn", "ifn","lps","lps","lps","lps","ifn","ifn", "ifn","ifn"), 28)
assay <- rep(c(rep("cyt",8), rep("rna",8)), 28)
#highlight significant RNA/cytokine induction
sig <- c(sum(fdr.lps2[1]<0.05 & abs(log2(fc.lps2[1]))>1), sum(fdr.lps2[1]<0.05 & abs(log2(fc.lps2[1]))>1), 
         sum(fdr.lps24[1]<0.05 & abs(log2(fc.lps24[1]))>1), sum(fdr.lps24[1]<0.05 & abs(log2(fc.lps24[1]))>1), 
         sum(fdr.ifn[1]<0.05 & abs(log2(fc.ifn[1]))>1), NA, NA, sum(fdr.ifn[1]<0.05 & abs(log2(fc.ifn[1]))>1), 
         sum(fdr.lps2.rna[1]<0.05& abs(log2(fc.lps2.rna[1]))>1), sum(fdr.lps2.rna[1]<0.05& abs(log2(fc.lps2.rna[1]))>1), 
         sum(fdr.lps24.rna[1]<0.05& abs(log2(fc.lps24.rna[1]))>1), sum(fdr.lps24.rna[1]<0.05& abs(log2(fc.lps24.rna[1]))>1), 
         sum(fdr.ifn.rna[1]<0.05& abs(log2(fc.ifn.rna[1]))>1), NA, NA, sum(fdr.ifn.rna[1]<0.05& abs(log2(fc.ifn.rna[1]))>1))

for (i in c(2:28)){
  sig <- c(sig, c(sum(fdr.lps2[i]<0.05 & abs(log2(fc.lps2[i]))>1), sum(fdr.lps2[i]<0.05 & abs(log2(fc.lps2[i]))>1), 
                  sum(fdr.lps24[i]<0.05 & abs(log2(fc.lps24[i]))>1), sum(fdr.lps24[i]<0.05 & abs(log2(fc.lps24[i]))>1), 
                  sum(fdr.ifn[i]<0.05 & abs(log2(fc.ifn[i]))>1), NA, NA, sum(fdr.ifn[i]<0.05 & abs(log2(fc.ifn[i]))>1), 
                  sum(fdr.lps2.rna[i]<0.05& abs(log2(fc.lps2.rna[i]))>1), sum(fdr.lps2.rna[i]<0.05& abs(log2(fc.lps2.rna[i]))>1), 
                  sum(fdr.lps24.rna[i]<0.05& abs(log2(fc.lps24.rna[i]))>1), sum(fdr.lps24.rna[i]<0.05& abs(log2(fc.lps24.rna[i]))>1), 
                  sum(fdr.ifn.rna[i]<0.05& abs(log2(fc.ifn.rna[i]))>1)), NA, NA, sum(fdr.ifn.rna[i]<0.05& abs(log2(fc.ifn.rna[i]))>1))
}


#read in Illumina probe IDs
ilmn.ids <- read.table("probe_lookup.txt", header = T)

#line up c.f. fold-changes
cyto <- rep(ilmn.ids$cyto[1],16)
for (i in c(2:28)){
  cyto <- c(cyto, rep(ilmn.ids$cyto[i],16))
}

#create data.frame for plotting
for.plot <- data.frame(cbind(time, fc, stim, assay, sig, cyto))
for.plot$time <- as.numeric(for.plot$time)
for.plot$fc <- as.numeric(for.plot$fc)
for.plot$sig <- as.numeric(for.plot$sig)

for.plot <- na.omit(for.plot)


cols <- brewer.pal(9, "Set1")

p<-ggplot(for.plot, aes(x=time, y=log2(fc))) +
  geom_line(data = subset(for.plot, assay=="cyt" & stim=="lps" & sig==1), color=cols[1], linetype=1, linewidth=2)+
  geom_line(data = subset(for.plot, assay=="rna" & stim=="lps" & sig==1), color=cols[1], linetype=2, linewidth=2)+
  geom_line(data = subset(for.plot, assay=="cyt" & stim=="lps" & sig==0), color=cols[1], linetype=1, linewidth=2, alpha=0.3)+
  geom_line(data = subset(for.plot, assay=="rna" & stim=="lps" & sig==0), color=cols[1], linetype=2, linewidth=2, alpha=0.3)+
  geom_line(data = subset(for.plot, assay=="cyt" & stim=="ifn" & sig==0), color=cols[2], linetype=1, linewidth=2, alpha=0.3)+
  geom_line(data = subset(for.plot, assay=="rna" & stim=="ifn" & sig==0), color=cols[2], linetype=2, linewidth=2, alpha=0.3)+
  geom_line(data = subset(for.plot, assay=="cyt" & stim=="ifn" & sig==1), color=cols[2], linetype=1, linewidth=2)+
  geom_line(data = subset(for.plot, assay=="rna" & stim=="ifn" & sig==1), color=cols[2], linetype=2, linewidth=2)+
  geom_point(data = subset(for.plot, assay=="cyt" & stim=="lps" & sig==1), color=cols[1], fill=cols[1], shape=16, size=5)+
  geom_point(data = subset(for.plot, assay=="cyt" & stim=="ifn" & sig==1), color=cols[2], fill=cols[2], shape=16, size=5)+
  geom_point(data = subset(for.plot, assay=="rna" & stim=="lps" & sig==1), color=cols[1], fill=cols[1], shape=16, size=5)+
  geom_point(data = subset(for.plot, assay=="rna" & stim=="ifn" & sig==1), color=cols[2], fill=cols[2], shape=16, size=5)+
  geom_point(data = subset(for.plot, assay=="cyt" & stim=="lps" & sig==0), color=cols[1], fill=cols[1], shape=16, size=5, alpha = 0.3)+
  geom_point(data = subset(for.plot, assay=="cyt" & stim=="ifn" & sig==0), color=cols[2], fill=cols[2], shape=16, size=5, alpha = 0.3)+
  geom_point(data = subset(for.plot, assay=="rna" & stim=="lps" & sig==0), color=cols[1], fill=cols[1], shape=16, size=5, alpha = 0.3)+
  geom_point(data = subset(for.plot, assay=="rna" & stim=="ifn" & sig==0), color=cols[2], fill=cols[2], shape=16, size=5, alpha = 0.3)+
  annotate("point", x=0, y=0, color=cols[9], size = 3)+
  xlab("hours")+
  ylab("logFC")+
  facet_wrap(~cyto, ncol=7)+
  scale_x_continuous(breaks=c(12,24),
                   labels=c("2", "24"))+
  theme_bw() +
  theme(axis.text.y=element_text(size=25), 
        axis.text.x=element_text(size=25), 
        axis.title=element_text(size=40), 
        strip.background = element_rect(color="black", fill=cols[9], size=1, linetype="solid"), strip.text.x = element_text(
          size = 30, color = "white", face = "bold"))
p


ggplot2::ggsave(
  "cytokine_induction.jpg",
  width = 20,
  height = 20,
  dpi = 300
)


#plot out bar plots of cytokine induction - Figure 1A

for.bar_plot <-  data.frame(rbind(cbind(ilmn.ids$cyto, as.integer(p.adjust(p.lps2, method = "fdr")<0.05 & abs(log2(fc.lps2))>1), log2(fc.lps2), "LPS - 2 hours"),
cbind(ilmn.ids$cyto, as.integer(p.adjust(p.lps24, method = "fdr")<0.05 & abs(log2(fc.lps24))>1), log2(fc.lps24), "LPS - 24 hours"),
cbind(ilmn.ids$cyto, as.integer(p.adjust(p.ifn, method = "fdr")<0.05 & abs(log2(fc.ifn))>1), log2(fc.ifn), "IFN - 24 hours")))

colnames(for.bar_plot) <- c("Cytokine", "sig", "log2FC", "cond")
for.bar_plot$sig <- as.numeric(for.bar_plot$sig)
for.bar_plot$log2FC <- as.numeric(for.bar_plot$log2FC)

for.bar_plot$cond <- factor(for.bar_plot$cond, levels = c("LPS - 2 hours", "LPS - 24 hours", "IFN - 24 hours"))

for.bar_plot$for_colour <- 1
for.bar_plot$for_colour[which(for.bar_plot$cond=="LPS - 2 hours" & for.bar_plot$sig==1)] <- 2
for.bar_plot$for_colour[which(for.bar_plot$cond=="LPS - 24 hours" & for.bar_plot$sig==0)] <- 3
for.bar_plot$for_colour[which(for.bar_plot$cond=="LPS - 24 hours" & for.bar_plot$sig==1)] <- 4
for.bar_plot$for_colour[which(for.bar_plot$cond=="IFN - 24 hours" & for.bar_plot$sig==0)] <- 5
for.bar_plot$for_colour[which(for.bar_plot$cond=="IFN - 24 hours" & for.bar_plot$sig==1)] <- 6

bp <- ggplot(for.bar_plot, aes(x=Cytokine, y=log2FC, fill=factor(for_colour))) +
  geom_bar(stat="identity") +
  facet_wrap(~cond, ncol=1) +
  scale_fill_manual(values = c(brewer.pal(8, "Set2")[8], brewer.pal(12, "Paired")[5],
                               brewer.pal(8, "Set2")[8], brewer.pal(12, "Paired")[6],
                               brewer.pal(8, "Set2")[8], brewer.pal(12, "Paired")[2]))+
theme_bw() +
  theme(axis.text.y=element_text(size=40), 
        axis.text.x=element_text(size=25, angle = 45, vjust = 1, hjust=1), 
        axis.title=element_text(size=40), 
        strip.background = element_rect(color="black", fill=cols[9], size=1, linetype="solid"), strip.text.x = element_text(
          size = 50, color = "white", face = "bold"),
        legend.position="none")
bp

ggplot2::ggsave(
  "cytokine_induction_barplot.jpg",
  width = 17,
  height = 18,
  dpi = 300
) 

#tabulate cytokine  and RNA fold-changes and significance across 3 stimulation conditions

table.out <- cbind(ilmn.ids$cyto, 
      fc.lps2.rna, p.lps2.rna, p.adjust(p.lps2.rna, method = "fdr"),
      fc.lps2, p.lps2, p.adjust(p.lps2, method = "fdr"),
      fc.lps24.rna, p.lps24.rna, p.adjust(p.lps24.rna, method = "fdr"),
      fc.lps24, p.lps24, p.adjust(p.lps24, method = "fdr"),
      fc.ifn.rna, p.ifn.rna, p.adjust(p.ifn.rna, method = "fdr"),
      fc.ifn, p.ifn, p.adjust(p.ifn, method = "fdr"))

colnames(table.out) <- c("cytokine", "lps2.rna.fc", "lps2.rna.p", "lps2.rna.fdr", "lps2.pro.fc", "lps2.pro.p", "lps2.pro.fdr",
                         "lps24.rna.fc", "lps24.rna.p", "lps24.rna.fdr", "lps24.pro.fc", "lps24.pro.p", "lps24.pro.fdr",
                         "ifn.rna.fc", "ifn.rna.p", "ifn.rna.fdr", "ifn.pro.fc", "ifn.pro.p", "ifn.pro.fdr")

#next explore the relationship between cytokine responses and transcription
#1st normalise data within stimulation condition
#first split into 4 conditions

ifn.raw.qc <- raw.data.qc[c(1:345),]
lps2.raw.qc <- raw.data.qc[c(346:690),]
lps24.raw.qc <- raw.data.qc[c(691:1035),]
ut.raw.qc <- raw.data.qc[c(1036:1380),]

ifn.covs <- covs[c(1:345),]
lps2.covs <- covs[c(346:690),]
lps24.covs <- covs[c(691:1035),]
ut.covs <- covs[c(1036:1380),]


#quantile normalize the supernatant (within batch and condition)
#ifn
ifn.norm <- matrix(NA,ncol=28,nrow=345)

for(j in c(4:19)){
  for (i in c(1:28)){
    plate <- which(ifn.covs$plate==j & !is.na(ifn.raw.qc[,i]))
    ifn.norm[plate,i] <- qnorm((0.5 + rank(ifn.raw.qc[plate,i]))/(1 + sum(!is.na(ifn.raw.qc[plate,i]))))
  }
}

#lps2
lps2.norm <- matrix(NA,ncol=28,nrow=345)

for(j in c(4:19)){
  for (i in c(1:28)){
    plate <- which(lps2.covs$plate==j & !is.na(lps2.raw.qc[,i]))
    lps2.norm[plate,i] <- qnorm((0.5 + rank(lps2.raw.qc[plate,i]))/(1 + sum(!is.na(lps2.raw.qc[plate,i]))))
  }
}

#lps24
lps24.norm <- matrix(NA,ncol=28,nrow=345)

for(j in c(4:19)){
  for (i in c(1:28)){
    plate <- which(lps24.covs$plate==j & !is.na(lps24.raw.qc[,i]))
    lps24.norm[plate,i] <- qnorm((0.5 + rank(lps24.raw.qc[plate,i]))/(1 + sum(!is.na(lps24.raw.qc[plate,i]))))
  }
}

#ut
ut.norm <- matrix(NA,ncol=28,nrow=345)

for(j in c(4:19)){
  for (i in c(1:28)){
    plate <- which(ut.covs$plate==j & !is.na(ut.raw.qc[,i]))
    ut.norm[plate,i] <- qnorm((0.5 + rank(ut.raw.qc[plate,i]))/(1 + sum(!is.na(ut.raw.qc[plate,i]))))
  }
}

#then regress out other confounders
#ifn
ifn.norm2 <- ifn.norm*NA

for (i in 1:28){
  keep <- !is.na(ifn.norm[,i]) & !is.na(ifn.covs$age) & !is.na(ifn.covs$sex) & !is.na(ifn.covs$batch) & !is.na(ifn.covs$monocyte)
  temp <- predict(lm(ifn.norm[keep,i] ~  ifn.covs$age[keep] + ifn.covs$sex[keep] + ifn.covs$batch[keep] + ifn.covs$monocyte[keep]))
  ifn.norm2[keep,i] <- ifn.norm[keep,i]  - temp
}

ifn.norm2 <- data.frame(ifn.norm2)
colnames(ifn.norm2) <- colnames(ifn.raw.qc)
rownames(ifn.norm2) <- rownames(ifn.raw.qc)


#lps2
lps2.norm2 <- lps2.norm*NA

for (i in 1:28){
  keep <- !is.na(lps2.norm[,i]) & !is.na(lps2.covs$age) & !is.na(lps2.covs$sex) & !is.na(lps2.covs$batch) & !is.na(lps2.covs$monocyte)
  temp <- predict(lm(lps2.norm[keep,i] ~  lps2.covs$age[keep] + lps2.covs$sex[keep] + lps2.covs$batch[keep] + lps2.covs$monocyte[keep]))
  lps2.norm2[keep,i] <- lps2.norm[keep,i]  - temp
}

lps2.norm2 <- data.frame(lps2.norm2)
colnames(lps2.norm2) <- colnames(lps2.raw.qc)
rownames(lps2.norm2) <- rownames(lps2.raw.qc)

#lps24
lps24.norm2 <- lps24.norm*NA

for (i in 1:28){
  keep <- !is.na(lps24.norm[,i]) & !is.na(lps24.covs$age) & !is.na(lps24.covs$sex) & !is.na(lps24.covs$batch) & !is.na(lps24.covs$monocyte)
  temp <- predict(lm(lps24.norm[keep,i] ~  lps24.covs$age[keep] + lps24.covs$sex[keep] + lps24.covs$batch[keep] + lps24.covs$monocyte[keep]))
  lps24.norm2[keep,i] <- lps24.norm[keep,i]  - temp
}

lps24.norm2 <- data.frame(lps24.norm2)
colnames(lps24.norm2) <- colnames(lps24.raw.qc)
rownames(lps24.norm2) <- rownames(lps24.raw.qc)

#ut
ut.norm2 <- ut.norm*NA

for (i in 1:28){
  keep <- !is.na(ut.norm[,i]) & !is.na(ut.covs$age) & !is.na(ut.covs$sex) & !is.na(ut.covs$batch) & !is.na(ut.covs$monocyte)
  temp <- predict(lm(ut.norm[keep,i] ~  ut.covs$age[keep] + ut.covs$sex[keep] + ut.covs$batch[keep] + ut.covs$monocyte[keep]))
  ut.norm2[keep,i] <- ut.norm[keep,i]  - temp
}

ut.norm2 <- data.frame(ut.norm2)
colnames(ut.norm2) <- colnames(ut.raw.qc)
rownames(ut.norm2) <- rownames(ut.raw.qc)

cyt.ind.lps2 <- lps2.norm2
cyt.ind.lps24 <- lps24.norm2
cyt.ind.ifn <- ifn.norm2

#read in gene expression (RNA) principal components
ifn.pc <- read.table("ifn_pc.txt", row.names = 1, header = T)
lps2.pc <- read.table("lps2_pc.txt", row.names = 1, header = T)
lps24.pc <- read.table("lps24_pc.txt", row.names = 1, header = T)
cd14.pc <- read.table("cd14_pc.txt", row.names = 1, header = T)

#start with IFN stimulation - line up datasets
ifn.rna.28.qc_ifn <- ifn.rna.28[rownames(cyt.ind.ifn),]
cd14.rna.28.qc_ifn <- naive.rna.28[rownames(cyt.ind.ifn),]
ifn.pc.qc_ifn <- data.frame(ifn.pc)[rownames(cyt.ind.ifn),]
cd14.pc.qc_ifn <- data.frame(cd14.pc)[rownames(cyt.ind.ifn),]


#collect r2, beta and p-value for association between stimulated cytokine secretion and stimulated RNA transcription

ifn.ifn.r2 <- c()
ifn.ifn.beta <- c()
ifn.ifn.p <- c()

for (i in c(1:28)){
  ifn.ifn.r2[i] <- summary(lm(cyt.ind.ifn[,i]~ifn.rna.28.qc_ifn[,i]))$r.squared
  ifn.ifn.beta[i] <- summary(lm(cyt.ind.ifn[,i]~ifn.rna.28.qc_ifn[,i]))$coef[2,1]
  ifn.ifn.p[i] <- summary(lm(cyt.ind.ifn[,i]~ifn.rna.28.qc_ifn[,i]))$coef[2,4]
}

#collect r2, beta and p-value for association between stimulated cytokine secretion and baseline RNA transcription

ifn.cd14.r2 <- c()
ifn.cd14.beta <- c()
ifn.cd14.p <- c()


for (i in c(1:28)){
  ifn.cd14.r2[i] <- summary(lm(cyt.ind.ifn[,i]~cd14.rna.28.qc_ifn[,i]))$r.squared
  ifn.cd14.beta[i] <- summary(lm(cyt.ind.ifn[,i]~cd14.rna.28.qc_ifn[,i]))$coef[2,1]
  ifn.cd14.p[i] <- summary(lm(cyt.ind.ifn[,i]~cd14.rna.28.qc_ifn[,i]))$coef[2,4]
}



#LPS 2 hours - line up datasets

lps2.rna.28.qc_lps2 <- lps2.rna.28[rownames(cyt.ind.lps2),]
cd14.rna.28.qc_lps2 <- naive.rna.28[rownames(cyt.ind.lps2),]
lps2.pc.qc_lps2 <- data.frame(lps2.pc)[rownames(cyt.ind.lps2),]
cd14.pc.qc_lps2 <- data.frame(cd14.pc)[rownames(cyt.ind.lps2),]


#collect r2, beta and p-value for association between stimulated cytokine secretion and stimulated RNA transcription

lps2.lps2.r2 <- c()
lps2.lps2.beta <- c()
lps2.lps2.p <- c()

for (i in c(1:28)){
  lps2.lps2.r2[i] <- summary(lm(cyt.ind.lps2[,i]~lps2.rna.28.qc_lps2[,i]))$r.squared
  lps2.lps2.beta[i] <- summary(lm(cyt.ind.lps2[,i]~lps2.rna.28.qc_lps2[,i]))$coef[2,1]
  lps2.lps2.p[i] <- summary(lm(cyt.ind.lps2[,i]~lps2.rna.28.qc_lps2[,i]))$coef[2,4]
}


#collect r2, beta and p-value for association between stimulated cytokine secretion and baseline RNA transcription

lps2.cd14.r2 <- c()
lps2.cd14.beta <- c()
lps2.cd14.p <- c()

for (i in c(1:28)){
  lps2.cd14.r2[i] <- summary(lm(cyt.ind.lps2[,i]~cd14.rna.28.qc_lps2[,i]))$r.squared
  lps2.cd14.beta[i] <- summary(lm(cyt.ind.lps2[,i]~cd14.rna.28.qc_lps2[,i]))$coef[2,1]
  lps2.cd14.p[i] <- summary(lm(cyt.ind.lps2[,i]~cd14.rna.28.qc_lps2[,i]))$coef[2,4]
}


#LPS 24 hours - line up datasets
lps24.rna.28.qc_lps24 <- lps24.rna.28[rownames(cyt.ind.lps24),]
lps2.rna.28.qc_lps24 <- lps2.rna.28[rownames(cyt.ind.lps24),]
cd14.rna.28.qc_lps24 <- naive.rna.28[rownames(cyt.ind.lps24),]
lps24.pc.qc_lps24 <- data.frame(lps24.pc)[rownames(cyt.ind.lps24),]
lps2.pc.qc_lps24 <- data.frame(lps2.pc)[rownames(cyt.ind.lps24),]
cd14.pc.qc_lps24 <- data.frame(cd14.pc)[rownames(cyt.ind.lps24),]

#collect r2, beta and p-value for association between stimulated cytokine secretion and stimulated RNA transcription

lps24.lps24.r2 <- c()
lps24.lps24.beta <- c()
lps24.lps24.p <- c()


for (i in c(1:28)){
  lps24.lps24.r2[i] <- summary(lm(cyt.ind.lps24[,i]~lps24.rna.28.qc_lps24[,i]))$r.squared
  lps24.lps24.beta[i] <- summary(lm(cyt.ind.lps24[,i]~lps24.rna.28.qc_lps24[,i]))$coef[2,1]
  lps24.lps24.p[i] <- summary(lm(cyt.ind.lps24[,i]~lps24.rna.28.qc_lps24[,i]))$coef[2,4]
}

#collect r2, beta and p-value for association between stimulated cytokine secretion and stimulated RNA transcription (LPS 2 hours)

lps24.lps2.r2 <- c()
lps24.lps2.beta <- c()
lps24.lps2.p <- c()

for (i in c(1:28)){
  lps24.lps2.r2[i] <- summary(lm(cyt.ind.lps24[,i]~lps2.rna.28.qc_lps24[,i]))$r.squared
  lps24.lps2.beta[i] <- summary(lm(cyt.ind.lps24[,i]~lps2.rna.28.qc_lps24[,i]))$coef[2,1]
  lps24.lps2.p[i] <- summary(lm(cyt.ind.lps24[,i]~lps2.rna.28.qc_lps24[,i]))$coef[2,4]
}

#collect r2, beta and p-value for association between stimulated cytokine secretion and basline RNA transcription

lps24.cd14.r2 <- c()
lps24.cd14.beta <- c()
lps24.cd14.p <- c()

for (i in c(1:28)){
  lps24.cd14.r2[i] <- summary(lm(cyt.ind.lps24[,i]~cd14.rna.28.qc_lps24[,i]))$r.squared
  lps24.cd14.beta[i] <- summary(lm(cyt.ind.lps24[,i]~cd14.rna.28.qc_lps24[,i]))$coef[2,1]
  lps24.cd14.p[i] <- summary(lm(cyt.ind.lps24[,i]~cd14.rna.28.qc_lps24[,i]))$coef[2,4]
}


#plot out correlation between RNA transcription and cytokine secretion - LPS 2 hours (Figure 2)

lps2.plot <- data.frame(rbind(cbind(ilmn.ids$cyto, lps2.lps2.r2, lps2.lps2.p, as.integer(p.adjust(lps2.lps2.p, method = "fdr")<0.05), "LPS2", "LPS - 2 hours"), 
      cbind(ilmn.ids$cyto, lps2.cd14.r2, lps2.cd14.p, as.integer(p.adjust(lps2.cd14.p, method = "fdr")<0.05), "Naive", "LPS - 2 hours")))

colnames(lps2.plot) <- c("Cytokine", "R2", "p.value", "sig", "RNA", "facet")

lps2.plot$R2 <- as.numeric(lps2.plot$R2)
lps2.plot$p.value <- as.numeric(lps2.plot$p.value)
lps2.plot[,4][which(lps2.plot[,4]==0)] <- "no"
lps2.plot[,4][which(lps2.plot[,4]==1)] <- "yes"

cols2 <- brewer.pal(12, "Paired")

lps2Plot<- ggplot(lps2.plot, aes(x = Cytokine, y = R2, colour = sig)) +
  geom_point(aes(size = -log10(p.value), pch = RNA))+
scale_size_continuous(range = c(3, 10)) +
  guides(size = "none") +
  scale_colour_manual(values = c(cols[9], cols2[5]), name = "FDR<0.05")+  expand_limits(y=c(0,0.6)) +
  #scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  ylab("R2") +
  theme_bw() +
  theme(axis.text.y=element_text(size=15), 
        axis.text.x=element_text(size=15), 
        axis.title=element_text(size=20), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = c(0.97, .28),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        strip.background = element_rect(color="black", fill=cols[9], size=1, linetype="solid"), strip.text.x = element_text(
          size = 20, color = "white", face = "bold")) +
  facet_wrap(~facet, ncol=1) +
  coord_flip()
lps2Plot

ggplot2::ggsave(
  "cytokine_RNA_relationship_lps2.jpg",
  width = 5,
  height = 9,
  dpi = 300
) 


#plot out correlation between RNA transcription and cytokine secretion - LPS 24 hours (Figure 2)

lps24.plot <- data.frame(rbind(cbind(ilmn.ids$cyto, lps24.lps24.r2, lps24.lps24.p, as.integer(p.adjust(lps24.lps24.p, method = "fdr")<0.05), "LPS24", "LPS - 24 hours"),
                              cbind(ilmn.ids$cyto, lps24.lps2.r2, lps24.lps2.p, as.integer(p.adjust(lps24.lps2.p, method = "fdr")<0.05), "LPS2", "LPS - 24 hours"), 
                              cbind(ilmn.ids$cyto, lps24.cd14.r2, lps24.cd14.p, as.integer(p.adjust(lps24.cd14.p, method = "fdr")<0.05), "Naive", "LPS - 24 hours")))

colnames(lps24.plot) <- c("Cytokine", "R2", "p.value", "sig", "RNA", "facet")

lps24.plot$R2 <- as.numeric(lps24.plot$R2)
lps24.plot$p.value <- as.numeric(lps24.plot$p.value)
lps24.plot[,4][which(lps24.plot[,4]==0)] <- "no"
lps24.plot[,4][which(lps24.plot[,4]==1)] <- "yes"

lps24.plot$RNA <- factor(lps24.plot$RNA, levels = c("LPS24", "Naive", "LPS2"))


lps24Plot<- ggplot(lps24.plot, aes(x = Cytokine, y = R2, colour = sig)) +
  geom_point(aes(size = -log10(p.value), pch = RNA))+
  scale_size_continuous(range = c(3, 10)) +
  guides(size = "none") +
  scale_colour_manual(values = c(cols[9], cols2[6]), name = "FDR<0.05")+  expand_limits(y=c(0,0.6)) +
  #scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  ylab("R2") +
  theme_bw() +
  theme(axis.text.y=element_text(size=15), 
        axis.text.x=element_text(size=15), 
        axis.title=element_text(size=20), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = c(0.97, .28),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        strip.background = element_rect(color="black", fill=cols[9], size=1, linetype="solid"), strip.text.x = element_text(
          size = 20, color = "white", face = "bold")) +
  facet_wrap(~facet, ncol=1) +
  coord_flip()
lps24Plot

ggplot2::ggsave(
  "cytokine_RNA_relationship_lps24.jpg",
  width = 5,
  height = 9,
  dpi = 300
) 


#plot out correlation between RNA transcription and cytokine secretion - IFN 24 hours (Figure 2)


ifn.plot <- data.frame(rbind(cbind(ilmn.ids$cyto, ifn.ifn.r2, ifn.ifn.p, as.integer(p.adjust(ifn.ifn.p, method = "fdr")<0.05), "IFN", "IFN - 24 hours"), 
                              cbind(ilmn.ids$cyto, ifn.cd14.r2, ifn.cd14.p, as.integer(p.adjust(ifn.cd14.p, method = "fdr")<0.05), "Naive", "IFN - 24 hours")))

colnames(ifn.plot) <- c("Cytokine", "R2", "p.value", "sig", "RNA", "facet")

ifn.plot$R2 <- as.numeric(ifn.plot$R2)
ifn.plot$p.value <- as.numeric(ifn.plot$p.value)
ifn.plot[,4][which(ifn.plot[,4]==0)] <- "no"
ifn.plot[,4][which(ifn.plot[,4]==1)] <- "yes"


ifnPlot<- ggplot(ifn.plot, aes(x = Cytokine, y = R2, colour = sig)) +
  geom_point(aes(size = -log10(p.value), pch = RNA))+
  scale_size_continuous(range = c(3, 10)) +
  guides(size = "none") +
  scale_colour_manual(values = c(cols[9], cols2[2]), name = "FDR<0.05")+  expand_limits(y=c(0,0.6)) +
  #scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  ylab("R2") +
  theme_bw() +
  theme(axis.text.y=element_text(size=15), 
        axis.text.x=element_text(size=15), 
        axis.title=element_text(size=20), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.key.size = unit(0.5, 'cm'),
        legend.position = c(0.97, .28),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        strip.background = element_rect(color="black", fill=cols[9], size=1, linetype="solid"), strip.text.x = element_text(
          size = 20, color = "white", face = "bold")) +
  facet_wrap(~facet, ncol=1) +
  coord_flip()
ifnPlot

ggplot2::ggsave(
  "cytokine_RNA_relationship_ifn.jpg",
  width = 5,
  height = 9,
  dpi = 300
) 


#Pairwise correlation of cytokines (Figure 1E)
library(ClassDiscovery)

#Start with IFN stimulations
cyt.ind.ifn.no_missing <- na.omit(cyt.ind.ifn)
#calculate pairwise distance matrix (absolute Pearson correlation coeff as metric)
dist.mat.ifn <- distanceMatrix(as.matrix(cyt.ind.ifn.no_missing), metric = "absolute pearson")

cor.ifn <- matrix(, nrow = 28, ncol = 28)

for(j in c(1:28)){
  for(k in c(1:28)){
    cor.ifn[j,k]  <- abs(cor(cyt.ind.ifn.no_missing[,j], cyt.ind.ifn.no_missing[,k], method = 'pearson'))
  }
}

#plot out heatmap
jpeg("ifn_heatmap.jpg")
heatmap(cor.ifn[hclust(dist.mat.ifn)$order,hclust(dist.mat.ifn)$order], Rowv=NA, Colv=NA, 
        scale="none", 
        labRow = colnames(cyt.ind.ifn.no_missing)[hclust(dist.mat.ifn)$order],
        labCol = colnames(cyt.ind.ifn.no_missing)[hclust(dist.mat.ifn)$order],
        col = rev(heat.colors(6)))
dev.off()

#LPS 2 hours
cyt.ind.lps2.no_missing <- na.omit(cyt.ind.lps2)
#calculate pairwise distance matrix (absolute Pearson correlation coeff as metric)

dist.mat.lps2 <- distanceMatrix(as.matrix(cyt.ind.lps2.no_missing), metric = "absolute pearson")

cor.lps2 <- matrix(, nrow = 28, ncol = 28)

for(j in c(1:28)){
  for(k in c(1:28)){
    cor.lps2[j,k]  <- abs(cor(cyt.ind.lps2.no_missing[,j], cyt.ind.lps2.no_missing[,k], method = 'pearson'))
  }
}

#plot out heatmap
jpeg("lps2_heatmap.jpg")
heatmap(cor.lps2[hclust(dist.mat.lps2)$order,hclust(dist.mat.lps2)$order], Rowv=NA, Colv=NA, 
        scale="none", 
        labRow = colnames(cyt.ind.lps2.no_missing)[hclust(dist.mat.lps2)$order],
        labCol = colnames(cyt.ind.lps2.no_missing)[hclust(dist.mat.lps2)$order],
        col = rev(heat.colors(6)))
dev.off()

#LPS 24 hours
cyt.ind.lps24.no_missing <- na.omit(cyt.ind.lps24)
#calculate pairwise distance matrix (absolute Pearson correlation coeff as metric)

dist.mat.lps24 <- distanceMatrix(as.matrix(cyt.ind.lps24.no_missing), metric = "absolute pearson")

cor.lps24 <- matrix(, nrow = 28, ncol = 28)

for(j in c(1:28)){
  for(k in c(1:28)){
    cor.lps24[j,k]  <- abs(cor(cyt.ind.lps24.no_missing[,j], cyt.ind.lps24.no_missing[,k], method = 'pearson'))
  }
}
#plot out heatmap

jpeg("lps24_heatmap.jpg")
heatmap(cor.lps24[hclust(dist.mat.lps24)$order,hclust(dist.mat.lps24)$order], Rowv=NA, Colv=NA, 
        scale="none", 
        labRow = colnames(cyt.ind.lps24.no_missing)[hclust(dist.mat.lps24)$order],
        labCol = colnames(cyt.ind.lps24.no_missing)[hclust(dist.mat.lps24)$order],
        col = rev(heat.colors(6)))
dev.off()

#next step is then to look at expression predictors genome-wide
#here we provide data to reproduce the the JUN/IL1b, CAVIN2/BDNF & PDGF-BB, and CTS6 & MMP9 associations
#complete gene expression data is available here: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2232

#JUN and IL-1beta association
for.plot.il1b_jun <- read.table("for.plot.il1b_jun.txt", row.names = 1, header = T)
#resid = baseline RNA expression corrected for gene expression PCs (20)

summary(lm(for.plot.il1b_jun$IL1B~for.plot.il1b_jun$resid))


p_labels6.2 = data.frame(facet = c("IFN - 24 hours"), 
                         label = c("italic(P)==7.8%*%10^-9"))

cyto_ifn.il1b_jun <- ggplot(for.plot.il1b_jun, aes(y=IL1B, x=resid)) +
  geom_point(col=cols[2], size = 2, alpha=0.7) +
  geom_smooth(method=lm, col=cols[2]) +
  xlab(expression(paste(italic("JUN"), " RNA expression"))) + ylab("IL-1b secretion") +
  facet_wrap(~facet, ncol = 1) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=20),
                     axis.title=element_text(size=20), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 20, color = "white", face = "bold")) + # ylim(NA, 11.3) + #xlim(-2.5,NA) +
  geom_text(x=1, y=2.5, aes(label=label), data=p_labels6.2, parse=TRUE, inherit.aes=F, size = 8, hjust = 0)

cyto_ifn.il1b_jun

ggplot2::ggsave(
  "jun_il1b.jpg",
  width = 6,
  height = 6,
  dpi = 300
) 

#CAVIN2 expression and PDGF-BB secretion

for.plot.pdgfb_cavin2 <- read.table("for.plot.pdgfb_cavin2.txt", row.names = 1, header =T)

p_labels.cb = data.frame(facet = c("LPS - 2 hours", "LPS - 24 hours", "IFN - 24 hours"), 
                         label = c("italic(P)==3.5%*%10^-8", "italic(P)==6.3%*%10^-13", "italic(P)==1.8%*%10^-9"))

p_labels.cb$facet <- factor(p_labels.cb$facet, levels = c("LPS - 2 hours", "LPS - 24 hours", "IFN - 24 hours"))

cyto.cavin2_pdgfb <- ggplot(for.plot.pdgfb_cavin2, aes(y=PDGFBB, x=resid, colour = facet)) +
  geom_point(size = 2, alpha=0.7) +
  geom_smooth(method=lm) +
  scale_colour_manual(values = c(brewer.pal(12, "Paired")[5]
                               , brewer.pal(12, "Paired")[6]
                               , brewer.pal(12, "Paired")[2]))+
  xlab(expression(paste(italic("CAVIN2"), " RNA expression"))) + ylab("PDGF-BB secretion") +
  facet_wrap(~facet, ncol = 3) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=20),
                     axis.title=element_text(size=20), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 20, color = "white", face = "bold")) + ylim(NA, 3) + #xlim(-2.5,NA) +
  geom_text(x=0.1, y=2.8, aes(label=label), data=p_labels.cb, parse=TRUE, inherit.aes=F, size = 6, hjust = 0)

cyto.cavin2_pdgfb

ggplot2::ggsave(
  "cavin2_pdgfbb.jpg",
  width = 9,
  height = 3,
  dpi = 300
) 

#CAVIN2 expression and BDNF secretion
for.plot.bdnf_cavin2 <- read.table("for.plot.bdnf_cavin2.txt", row.names = 1, header =T)

p_labels.cbd = data.frame(facet = c("LPS - 2 hours", "LPS - 24 hours", "IFN - 24 hours"), 
                         label = c("italic(P)==3.7%*%10^-10", "italic(P)==5.4%*%10^-18", "italic(P)==18.5%*%10^-20"))

p_labels.cbd$facet <- factor(p_labels.cb$facet, levels = c("LPS - 2 hours", "LPS - 24 hours", "IFN - 24 hours"))

cyto.cavin2_bdnf <- ggplot(for.plot.bdnf_cavin2, aes(y=BDNF, x=resid, colour = facet)) +
  geom_point(size = 2, alpha=0.7) +
  geom_smooth(method=lm) +
  scale_colour_manual(values = c(brewer.pal(12, "Paired")[5]
                                 , brewer.pal(12, "Paired")[6]
                                 , brewer.pal(12, "Paired")[2]))+
  xlab(expression(paste(italic("CAVIN2"), " RNA expression"))) + ylab("BDNF secretion") +
  facet_wrap(~facet, ncol = 3) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=20),
                     axis.title=element_text(size=20), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 20, color = "white", face = "bold")) + ylim(NA, 3) + #xlim(-2.5,NA) +
  geom_text(x=0, y=2.8, aes(label=label), data=p_labels.cbd, parse=TRUE, inherit.aes=F, size = 6, hjust = 0)


ggplot2::ggsave(
  "cavin2_bdnf.jpg",
  width = 9,
  height = 3,
  dpi = 300
) 

#Stimulated MMP9 expression and cytokine secretion
plot.mmp9 <- read.table("plot.mmp9.txt", row.names = 1, header = T)
plot.mmp9$cytokine <- factor(plot.mmp9$cytokine, levels = c("GRO.alpha", "HGF", "IL.21", "IL.6", "MCP.1", "MIP.1.alpha", "PIGF.1", "SDF.1.alpha"))


mmp9.6.pvals <- data.frame(cytokine = levels(plot.mmp9$cytokine), 
                           label = c("italic(P)==3.4%*%10^-8", 
                                     "italic(P)==9.4%*%10^-11", 
                                     "italic(P)==2.7%*%10^-7",
                                     "italic(P)==1.6%*%10^-10",
                                     "italic(P)==1.5%*%10^-8",
                                     "italic(P)==4.7%*%10^-10",
                                     "italic(P)==3.0%*%10^-8",
                                     "italic(P)==3.1%*%10^-8"))

mmp9.6.pvals$cytokine <- factor(mmp9.6.pvals$cytokine, levels = c("GRO.alpha", "HGF", "IL.21", "IL.6", "MCP.1", "MIP.1.alpha", "PIGF.1", "SDF.1.alpha"))

cols1 <- brewer.pal(9, "Set2")

mmp9 <- ggplot(plot.mmp9, aes(y=cyto, x=RNA)) +
  geom_point(size = 1, alpha=0.7, color = cols1[8]) +
  geom_smooth(method=lm, color = "black") +
  xlab(expression(paste(italic("MMP9"), " RNA expression"))) + ylab("Cytokine secretion") +
  facet_wrap(~cytokine, ncol = 8) +
  scale_x_continuous(breaks = c(-1,0,1)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=10),
                     axis.title=element_text(size=10), strip.background = element_rect(color="black", fill=cols1[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 10, color = "white", face = "bold")) + #ylim(NA, 3) + #xlim(-2.5,NA) +
  geom_text(x=-0, y=2.8, aes(label=label), data=mmp9.6.pvals, parse=TRUE, inherit.aes=F, size = 3, hjust = 0.5)
mmp9

ggplot2::ggsave(
  "MMP9.jpg",
  width = 9,
  height = 3,
  dpi = 300
) 


#Stimulated CTS6 expression and cytokine secretion

plot.cts6 <- read.table("plot.cts6.txt", row.names = 1, header = T)

plot.cts6$cyto <- as.numeric(plot.cts6$cyto)
plot.cts6$RNA <- as.numeric(plot.cts6$RNA)
plot.cts6$cytokine <- factor(plot.cts6$cytokine, levels = c("GRO.alpha", "IL.1.beta", "IL.21", "MCP.1", "MIP.1.alpha", "PIGF.1", "SDF.1.alpha", "MCP.1.lps2"))


cts.6.pvals <- data.frame(cytokine = levels(plot.cts6$cytokine), 
                          label = c("italic(P)==2.1%*%10^-8", 
                                    "italic(P)==2.8%*%10^-10", 
                                    "italic(P)==2.2%*%10^-9",
                                    "italic(P)==3.9%*%10^-10",
                                    "italic(P)==6.4%*%10^-11",
                                    "italic(P)==1.3%*%10^-9",
                                    "italic(P)==4.7%*%10^-9",
                                    "italic(P)==2.0%*%10^-8"))


cts.6.pvals$cytokine <- factor(cts.6.pvals$cytokine, levels = c("GRO.alpha", "IL.1.beta", "IL.21", "MCP.1", "MIP.1.alpha", "PIGF.1", "SDF.1.alpha", "MCP.1.lps2"))

cts6 <- ggplot(plot.cts6, aes(y=cyto, x=RNA)) +
  geom_point(size = 1, alpha=0.7, color = cols1[8]) +
  geom_smooth(method=lm, color = "black") +
  xlab(expression(paste(italic("CTS6"), " RNA expression"))) + ylab("Cytokine secretion") +
  facet_wrap(~cytokine, ncol = 8) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=10),
                     axis.title=element_text(size=10), strip.background = element_rect(color="black", fill=cols1[8], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 10, color = "white", face = "bold")) + #ylim(NA, 3) + #xlim(-2.5,NA) +
  geom_text(x=-0, y=2.8, aes(label=label), data=cts.6.pvals, parse=TRUE, inherit.aes=F, size = 3, hjust = 0.5)
cts6

ggplot2::ggsave(
  "CTS6.jpg",
  width = 9,
  height = 3,
  dpi = 300
) 
