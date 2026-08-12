rm(list = ls(all=TRUE))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(patchwork)
library(ggbio)
library(data.table)
library(GenomicRanges)
library(biomaRt)
library(stringr)


#First perform multivariate GWAS (MANOVA) and univariate GWAS for each cytokine passing quality control

#args 1-22 for each chromosome
args = commandArgs(trailingOnly=TRUE)

#read in genotype (dosages) for each chromosome - individual-level genotypes are available at the European Genome-Phenome Archive (EGA) with accession ID EGAS00000000109. 
geno <- read.table(paste0("chr", args[1], ".dosage.gz"), header = T)

#collect SNP information from dosage file
chrom <- geno$chromosome
pos <- geno$position
pos <- geno$position
rsid <- geno$rsid
alleleA <- geno$alleleA
alleleB <- geno$alleleB

#remove SNP information from genotype matrix
geno$chromosome <- NULL
geno$SNPID <- NULL
geno$rsid <- NULL
geno$position <- NULL
geno$alleleA <- NULL
geno$alleleB <- NULL

#read in covariate-corrected, normalised cytokine secretion data (here following 2 hours of LPS stimulation - data for IFN 24 hour and LPS 24 hour stimulated samples are also provided)  
lps2.cyto <- read.table("./data/cyto_data.norm_lps2.geno_corr.txt", row.names = 1, header = T)

#transpose genotype matrix
geno.t <- t(geno)

#restrict genotype matrix to samples with cytokine secrection data in that condition
geno.t.incl <- data.frame(geno.t[rownames(lps2.cyto),])

#perform multivariate GWAS with MANOVA

#empty vectors to collect p-values and F-statistics
fstats <- c()
pvals <- c()

#recursively perform MANOVA for each SNP
for (i in c(1:dim(geno.t.incl)[2])){
  fstats[i] <- summary(manova(as.matrix(lps2.cyto)~geno.t.incl[,i]))$stats[1,3]
  pvals[i] <- summary(manova(as.matrix(lps2.cyto)~geno.t.incl[,i]))$stats[1,6]
}

#output MANOVA results
write.table(cbind(chrom, pos, rsid, alleleA, alleleB,  fstats, pvals), paste0("chr",args[1],"manova.lps2.out"), row.names = F, col.names = F, quote = F)

#then perform univariate GWAS for each of 28 cytokines passing QC

#empty vectors to collect betas, standard errors and p-values
beta <- c()
se <- c()
p <- c()

#recursively perform linear regression for each SNP for the first cytokine in the phenotype matrix
for (i in c(1:dim(geno.t.incl)[2])){
  beta[i] <- summary(lm(lps2.cyto[,1]~geno.t.incl[,i]))$coef[2,1]
  se[i] <- summary(lm(lps2.cyto[,1]~geno.t.incl[,i]))$coef[2,2]
  p[i] <- summary(lm(lps2.cyto[,1]~geno.t.incl[,i]))$coef[2,4]
}

#collect output for first cytokine into a dataframe
out <- data.frame(cbind(beta, se, p))

#recursively perform the same GWAS analysis for each of the remaining 27 cytokines - appending association statistics to the output dataframe
for (j in c(2:28)){
  beta <- c()
  se <- c()
  p <- c()
  for (i in c(1:dim(geno.t.incl)[2])){
    beta[i] <- summary(lm(lps2.cyto[,j]~geno.t.incl[,i]))$coef[2,1]
    se[i] <- summary(lm(lps2.cyto[,j]~geno.t.incl[,i]))$coef[2,2]
    p[i] <- summary(lm(lps2.cyto[,j]~geno.t.incl[,i]))$coef[2,4]
  }
  out <- cbind(out, beta, se, p)
  print(j)
}

#add in SNP information to output matrix
out2 <- cbind(chrom, pos, rsid, alleleA, alleleB,  out)

#add column names including cytokine information
colnames(out2) <- c("Chr", "POS", "SNP", "alleleA", "alleleB", paste0(c("beta.", "se.", "p."), rep(colnames(lps2.cyto), each=3)))

#output univariate GWAS results
write.table(out2, paste0("chr",args[1],"univariate.lps2.out"), row.names = F, col.names = T, quote = F)

#Next visualise results of interest

#read in summary statistics: Mulvariate GWAS (MANOVA) of cytokine secretion in monocyte (IFN-gamma 24 hours)
#summary statistics are randomly thinned (1 in 10) with the exception of chr3 and chr22 regions
#full summary statistics will be deposited with the GWAS Catalog (https://www.ebi.ac.uk/gwas/)

#Manhattan plot for multivariate GWAS (MANOVA) of monocyte cytokine secretion following 24 hours of IFN-gamma stimulation
disc.add <- read.table("./summ_stats/manova.IFN24h_stim.thinned_stats.txt.gz", header = T)

disc.add <- na.omit(disc.add)

head(disc.add)

str(disc.add)
disc.add$chr <- factor(disc.add$chr, levels = c(1:22))
levels(disc.add$chr) <- c(1:22)


#cumulative base pair position
disc.add2 <- disc.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(disc.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = disc.add2 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#need to highlight associated loci
disc.add2$label <- NA
disc.add2$label[which(disc.add2$BPcum==538676599)] <- "rs113010081"
disc.add2$label[which(disc.add2$BPcum==2867186538)] <- "rs5757584"

table(disc.add$label)

disc.add2$anno <- NA
disc.add2$anno[which(disc.add2$BPcum==538676599)] <- 1
disc.add2$anno[which(disc.add2$BPcum==2867186538)] <- 1


#set colour palette
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")


#plot manhattan plot
ifn24_manova.manh <- ggplot(disc.add2, aes(x=BPcum, y=-log10(p))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  
  # Add test using ggrepel to avoid overlapping
  geom_text_repel( data=subset(disc.add2, anno==1), aes(label=label), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 2) +
  
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = disc.add2, linetype="dashed", color = "red") +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","4","8","12","16","20"), breaks = c(0,4,8,12,16,20), expand = c(0, 0), limits= c(0,26)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20)
  )

ggsave(
  "~/Documents/cytokine_project/sept24/for_github/ifn24_manova.manh.jpg",
  width = 14,
  height = 3.5,
  dpi = 300
)


#Manhattan plot for multivariate GWAS (MANOVA) of monocyte cytokine secretion following 24 hours of LPS stimulation
disc.add <- read.table("./summ_stats/manova.LPS24h_stim.thinned_stats.txt.gz", header = T)

disc.add <- na.omit(disc.add)

head(disc.add)

str(disc.add)
disc.add$chr <- factor(disc.add$chr, levels = c(1:22))
levels(disc.add$chr) <- c(1:22)


#cumulative base pair position
disc.add2 <- disc.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(disc.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = disc.add2 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#need to highlight associated loci
disc.add2$label <- NA
disc.add2$label[which(disc.add2$BPcum==538676599)] <- "rs113010081"
disc.add2$label[which(disc.add2$BPcum==2867186538)] <- "rs5757584"
disc.add2$label[which(disc.add2$BPcum==363029430)] <- "rs11123160"

disc.add2$anno <- NA
disc.add2$anno[which(disc.add2$BPcum==538676599)] <- 1
disc.add2$anno[which(disc.add2$BPcum==2867186538)] <- 1
disc.add2$anno[which(disc.add2$BPcum==363029430)] <- 1


#set colour palette
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")


#plot manhattan plot
lps24_manova.manh <- ggplot(disc.add2, aes(x=BPcum, y=-log10(p))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  
  # Add test using ggrepel to avoid overlapping
  geom_text_repel( data=subset(disc.add2, anno==1), aes(label=label), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 2) +
  
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = disc.add2, linetype="dashed", color = "red") +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","4","8","12","16","20"), breaks = c(0,4,8,12,16,20), expand = c(0, 0), limits= c(0,26)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20)
  )

ggsave(
  "~/Documents/cytokine_project/sept24/for_github/lps24_manova.manh.jpg",
  width = 14,
  height = 3.5,
  dpi = 300
)


#Manhattan plot for univariate GWAS of monocyte IP-10 secretion following 24 hours of LPS stimulation
disc.add <- read.table("./summ_stats/IP10.LPS24h_stim.thinned_stats.txt.gz", header = T)

disc.add <- na.omit(disc.add)

head(disc.add)

str(disc.add)
disc.add$chr <- factor(disc.add$chr, levels = c(1:22))
levels(disc.add$chr) <- c(1:22)


#cumulative base pair position
disc.add2 <- disc.add %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(disc.add, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = disc.add2 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


#need to highlight associated loci
disc.add2$label <- NA
disc.add2$label[which(disc.add2$BPcum==1558848251)] <- "rs13296842"

disc.add2$anno <- NA
disc.add2$anno[which(disc.add2$BPcum==1558848251)] <- 1



#set colour palette
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")


#plot manhattan plot
ip10_univar.manh <- ggplot(disc.add2, aes(x=BPcum, y=-log10(p))) +
  
  # Show all points
  geom_point( aes(color=as.factor(chr)), size=2) +
  scale_color_manual(values = rep(c(cols3[2], cols3[1]), 11 )) +
  
  
  # Add test using ggrepel to avoid overlapping
  geom_text_repel( data=subset(disc.add2, anno==1), aes(label=label), size=5, min.segment.length = unit(0, 'lines'),
                   nudge_y = 2) +
  
  #gwas sig line
  geom_segment(aes(x = 10177, y = 7.30103, xend = 2879943885, yend = 7.30103), data = disc.add2, linetype="dashed", color = "red") +
  
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:18,20,22)], breaks= axisdf$center[c(1:18,20,22)] ) +
  scale_y_continuous( labels = c("0","4","8","12","16","20"), breaks = c(0,4,8,12,16,20), expand = c(0, 0), limits= c(0,26)) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20)
  )

ggsave(
  "~/Documents/cytokine_project/sept24/for_github/ip10_univar.manh.jpg",
  width = 14,
  height = 3.5,
  dpi = 300
)

#plot out regions of interest, illustrate effects on individual cytokines and consider evidence for colocalisation with human disease/traits
#we have used the CCR5del locus as an example here but the code and approach is identical for other loci
#1. Chromosome 3 region - IFN-gamma 24 hour stimulated monocytes - MANOVA summary statistics
#read in summary statistics
chr3.locus <- read.table("./summ_stats/manova.IFN24h_stim.chr3_locus.txt", header = T)

#download LD matrix for region SNPs
ld.out1 <- LDmatrix(chr3.locus$rsid2, pop = "CEU", r2d = "r2", token = "[token]", file = FALSE)

#extract pairwise LD for all SNP to peak associated SNP
chr3.locus$r2 <- ld.out1$rs113010081[match(chr3.locus$rsid2, ld.out1$RS_number)]

#categorise pairwise LD into bins 
chr3.locus$bin_r2 <- 1
chr3.locus$bin_r2[which(chr3.locus$r2>0.2 & chr3.locus$r2 <= 0.5)] <- 2
chr3.locus$bin_r2[which(chr3.locus$r2>0.5 & chr3.locus$r2 <= 0.8)] <- 3
chr3.locus$bin_r2[which(chr3.locus$r2>0.8)] <- 4

#set colour palettes
cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(11,"Spectral")
cols3 <- brewer.pal(8,"Paired")

#A: plot regional association
chr3.locus$annotate <- 0
chr3.locus$annotate[which(chr3.locus$rsid2=="rs113010081")] <- 1

ra.plot <- ggplot(chr3.locus, aes(x=bp, y=-log10(p))) + 
  xlim(min(chr3.locus$bp), max(chr3.locus$bp)) +
  geom_point(data=subset(chr3.locus, bin_r2==1), color=cols[8], size=3) + 
  geom_point(data=subset(chr3.locus, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(chr3.locus, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(chr3.locus, bin_r2==4), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,5,10,15,20,25)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(chr3.locus, annotate==1), aes(label=rsid2), size=5, col = c("black"), nudge_x = 50000) +
  annotate("text", x = (min(chr3.locus$bp)+20000+400000), y = 20, label = quote(r^2), hjust = 0.5,size=3) +
  annotate("point", x = (min(chr3.locus$bp)+60000+400000), y = 18, size = 3, colour = cols2[1]) +
  annotate("text", x = (min(chr3.locus$bp)+20000+400000), y = 18, label = c("0.8-1.0"), hjust = 0.5,size=3) +
  annotate("point", x = (min(chr3.locus$bp)+60000+400000), y =16, size = 3, colour = cols2[3]) +
  annotate("text", x = (min(chr3.locus$bp)+20000+400000), y = 16, label = c("0.5-0.8"), hjust = 0.5,size=3) +
  annotate("point", x = (min(chr3.locus$bp)+60000+400000), y = 14, size = 3, colour = cols2[5]) +
  annotate("text", x = (min(chr3.locus$bp)+20000+400000), y = 14, label = c("0.2-0.5"), hjust = 0.5,size=3) +
  annotate("point", x = (min(chr3.locus$bp)+60000+400000), y = 12, size = 3, colour = cols[8]) +
  annotate("text", x = (min(chr3.locus$bp)+20000+400000), y = 12, label = c("<0.2"), hjust = 0.5,size=3)
ra.plot


#B: plot out genes
#download coordinates and strand of local genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=3
sel.pos=min(chr3.locus$bp)+(max(chr3.locus$bp)-min(chr3.locus$bp))/2
range=1000000

out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'strand'), 
  filters = c('chromosome_name','start','end'), 
  values = list(sel.chr, sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)

#define midpoint of each gene
out.bm.genes.region$mid <- out.bm.genes.region$start_position+(out.bm.genes.region$end_position-out.bm.genes.region$start_position)/2

#restrict to protein coding genes
genes <- subset(out.bm.genes.region, gene_biotype=="protein_coding")

genes$start <- genes$start_position
genes$end <- genes$end_position

genes$start[which(genes$strand==-1)] <- genes$end_position[which(genes$strand==-1)]
genes$end[which(genes$strand==-1)] <- genes$start_position[which(genes$strand==-1)]

#limit to genes in region of peak
genes <- genes[c(10:19),]

#plot out genes
plot.range <- c(min(chr3.locus$bp), max(chr3.locus$bp))
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(min(chr3.locus$bp), max(chr3.locus$bp)) +
  ylim(c(1.9,2.2)) +
  geom_segment(data = genes[seq(1,10,3),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(2,10,3),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(3,10,3),],
               aes(x=start, xend=end, y=order+1.2, yend=order+1.2), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(1,10,3),], aes(x=mid,y=order+1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(2,10,3),], aes(x=mid,y=order+1.1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(3,10,3),], aes(x=mid,y=order+1.2, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
genes.plot


#C: plot out recombination rate
#read in data
recomb <- read.table("./genetic_maps/genetic_map_chr3_combined_b37.txt", header = T)
recomb <- subset(recomb, position>min(chr3.locus$bp) & position<max(chr3.locus$bp))

recomb_rate <- ggplot(recomb, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 3") +
  scale_x_continuous(breaks=c(46300000,46400000,46500000,46600000,46700000),
                     labels=c("46.3Mb", "46.4Mb", "46.5Mb", "46.6Mb", "46.7Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

p1 <- (ra.plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 1.5, 0.5))

ggplot2::ggsave(
  "ccr5del_ra.jpg",
  width = 5,
  height = 5,
  dpi = 300
) 

#Forest plots to depict the effect size of CCR5del on RANTES and MIP-1beta secretion
label <- c("LPS - 2h", "LPS - 24h", "IFN - 24h", "LPS - 2h", "LPS - 24h", "IFN - 24h")
mean  <- c(-0.05, 0.10, 0.73, 0.25, 0.36, 0.06)
se  <- c(0.14, 0.13, 0.12, 0.14, 0.13, 0.12)
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- c(rep("MIP1b/CCL4",3), rep("RANTES/CCL5", 3))

df <- data.frame(label, mean, lower, upper, facet)

#plot forest plot of effects at rs113010081

ccr5.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2, size=2) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols3[c(2,6,5)])) + scale_colour_manual(values = c(cols3[c(2,6,5)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("beta") + scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_y_continuous(breaks=c(0.0, 0.5),
  labels=c("0.0", "0.5")) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=15),
                     axis.title.x=element_text(size=15), axis.title.y=element_blank()) +
  facet_wrap(~facet, ncol = 1) +
  theme(strip.background = element_rect(color="black", fill=cols[8], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
ccr5.fp


ggplot2::ggsave(
  "ccr5del_effects.jpg",
  width = 3,
  height = 5,
  dpi = 300
)

#consider evidence of colocalisation between CCR5del locus and human disease/traits - we have used the CCR5del locus as an example here but the code and approach is identical for other loci
library(coloc)

#start with UK-Biobank summary statistics: https://www.nealelab.is/uk-biobank
#index of 47 UK Biobank traits
ukb.index <- read.table("./gwas_studies/ukbb.cont", header = T)

#vector to catch colocalisation statistics
peak.snp <- c()
peak.p <- c()
peak.beta <- c()
peak.se <- c()
ccr5d.p <- c()
ccr5d.beta <- c()
ccr5d.se <- c()
abf <- c()

#run colocalisation analysis across 47 traits against summary statistics describing association of CCR5del region with MIP-1b secretion (following IFN-gamma stimulation)

for (i in c(1:47)){
  #read in summary stats for each UKB GWAS
  ukb <- read.table(paste0("./ukb_stats/", ukb.index$phenotype[i], ".gwas.chr3_region"), header = T)
  
  #construct shared regional associations
  chr3.locus.coloc <- chr3.locus
  chr3.locus.coloc$ukb.beta <- ukb$beta[match(chr3.locus.coloc$rsid2, ukb$rsid)]
  chr3.locus.coloc$ukb.se <- ukb$se[match(chr3.locus.coloc$rsid2, ukb$rsid)]
  chr3.locus.coloc$ukb.p <- ukb$pval[match(chr3.locus.coloc$rsid2, ukb$rsid)]
  chr3.locus.coloc$minor_AF.1 <- ukb$minor_AF.1[match(chr3.locus.coloc$rsid2, ukb$rsid)]
  
  chr3.locus.coloc <- na.omit(chr3.locus.coloc)
  #collect peak UKB GWAS association in region
  peak.snp[i] <- chr3.locus.coloc$rsid2[which(chr3.locus.coloc$ukb.p==min(chr3.locus.coloc$ukb.p))]
  peak.p[i] <- min(chr3.locus.coloc$ukb.p)
  peak.beta[i] <- chr3.locus.coloc$ukb.beta[which(chr3.locus.coloc$ukb.p==min(chr3.locus.coloc$ukb.p))]
  peak.se[i] <- chr3.locus.coloc$ukb.se[which(chr3.locus.coloc$ukb.p==min(chr3.locus.coloc$ukb.p))]
  #collect CCR5del association for the UKB GWAS
  ccr5d.p[i] <- chr3.locus.coloc$ukb.p[which(chr3.locus.coloc$rsid2=="rs113010081")]
  ccr5d.beta[i] <- chr3.locus.coloc$ukb.beta[which(chr3.locus.coloc$rsid2=="rs113010081")]
  ccr5d.se[i] <- chr3.locus.coloc$ukb.beta[which(chr3.locus.coloc$rsid2=="rs113010081")]
  #set up lists for coloc
  cyt.list <- list(beta = chr3.locus.coloc$mip1b.beta, varbeta = (chr3.locus.coloc$mip1b.se)^2, type = "quant", N = 345, MAF=chr3.locus.coloc$minor_AF.1, snp=chr3.locus.coloc$rsid2)
  gwas.list <- list(beta = chr3.locus.coloc$ukb.beta, varbeta = (chr3.locus.coloc$ukb.se)^2, type = "quant", N = ukb.index$nos[i], MAF=chr3.locus.coloc$minor_AF.1, snp=chr3.locus.coloc$rsid2)
  #run colocalisation
  abf.res <- coloc.abf(cyt.list, gwas.list)
  #collect ABF
  abf[i] <- abf.res$summary[6]
}

#output UKB colocalisation results
chr3.ukb.out <- cbind(ukb.index, peak.snp, peak.p, peak.beta, peak.se, ccr5d.p, ccr5d.beta, ccr5d.se,  abf)
write.table(chr3.ukb.out, "coloc.ukb_chr3.txt")

#run colocalisation analysis across 100 GWAS catalog traits against summary statistics describing association of CCR5del region with MIP-1b secretion (following IFN-gamma stimulation)
#read in index describing GWAS catalogue traits
gwas_cat.index <- read.table("./gwas_studies/gwas.cat3", header = F)

#empty vectors to catch colocalisation statistics
peak.snp <- rep(NA,100)
peak.p <- rep(NA,100)
peak.beta <- rep(NA,100)
peak.se <- rep(NA,100)
ccr5d.p <- rep(NA,100)
ccr5d.beta <- rep(NA,100)
ccr5d.se <- rep(NA,100)
abf <- rep(NA,100)

#run colocalisation analysis across 100 GWAS catalog traits against summary statistics describing association of CCR5del region with MIP-1b secretion (following IFN-gamma stimulation)
for (i in c(1:100)){
  #read in summary stats for each GWAS
  gwas.cat <- read.table(paste0("./gwas_cat_stats/", gwas_cat.index$V1[i], ".gwas.chr3_region"), header = T)
  #if no variants in CCR5 region - skip
  if(dim(gwas.cat)[1]==0) next
  #construct shared regional associations
  gwas.cat$rsid2 <- paste0(gwas.cat$chr, ":", gwas.cat$bp)
  
  chr3.locus.coloc <- chr3.locus
  chr3.locus.coloc$gwas.beta <- gwas.cat$beta[match(chr3.locus.coloc$rsid, gwas.cat$rsid2)]
  chr3.locus.coloc$gwas.se <- gwas.cat$se[match(chr3.locus.coloc$rsid, gwas.cat$rsid2)]
  chr3.locus.coloc$gwas.p <- gwas.cat$p[match(chr3.locus.coloc$rsid, gwas.cat$rsid2)]
  chr3.locus.coloc$minor_AF.1 <- ukb$minor_AF.1[match(chr3.locus.coloc$rsid2, ukb$rsid)]

  chr3.locus.coloc <- na.omit(chr3.locus.coloc)
  #skip if no overlapping variants between QTL and GWAS
  if (dim(chr3.locus.coloc)[1]==0) next
  if (sum(chr3.locus.coloc$gwas.beta==0)>0){chr3.locus.coloc <- chr3.locus.coloc[-which(chr3.locus.coloc$gwas.beta==0),]}
  if (sum(chr3.locus.coloc$gwas.se==0)>0){chr3.locus.coloc <- chr3.locus.coloc[-which(chr3.locus.coloc$gwas.se==0),]}
  
  #collect peak GWAS association in region
  peak.snp[i] <- chr3.locus.coloc$rsid2[which(chr3.locus.coloc$gwas.p==min(chr3.locus.coloc$gwas.p))]
  peak.p[i] <- min(chr3.locus.coloc$gwas.p)
  peak.beta[i] <- chr3.locus.coloc$gwas.beta[which(chr3.locus.coloc$gwas.p==min(chr3.locus.coloc$gwas.p))]
  peak.se[i] <- chr3.locus.coloc$gwas.se[which(chr3.locus.coloc$gwas.p==min(chr3.locus.coloc$gwas.p))]
  #collect CCR5del association for the GWAS
  if (sum(which(chr3.locus.coloc$rsid2=="rs113010081"))==0) next
  ccr5d.p[i] <- chr3.locus.coloc$gwas.p[which(chr3.locus.coloc$rsid2=="rs113010081")]
  ccr5d.beta[i] <- chr3.locus.coloc$gwas.beta[which(chr3.locus.coloc$rsid2=="rs113010081")]
  ccr5d.se[i] <- chr3.locus.coloc$gwas.se[which(chr3.locus.coloc$rsid2=="rs113010081")]
  #set up GWAS and cytokine QTL lists for coloc
  cyt.list <- list(beta = chr3.locus.coloc$mip1b.beta, varbeta = (chr3.locus.coloc$mip1b.se)^2, type = "quant", N = 345, MAF=chr3.locus.coloc$minor_AF.1, snp=chr3.locus.coloc$rsid)
  gwas.list <- list(beta = chr3.locus.coloc$gwas.beta, varbeta = (chr3.locus.coloc$gwas.se)^2, type = "cc", N = gwas_cat.index$V2[i], s = gwas_cat.index$V2[i]/(gwas_cat.index$V2[i]+gwas_cat.index$V3[i]), snp=chr3.locus.coloc$rsid)
  #run coloc
  abf.res <- coloc.abf(cyt.list, gwas.list)
  abf[i] <- abf.res$summary[6]
}
#output GWAS catalog colocalisation results
chr3.gwas.out <- na.omit(cbind(gwas_cat.index, peak.snp, peak.p, peak.beta, peak.se, ccr5d.p, ccr5d.beta, ccr5d.se, abf))
write.table(chr3.gwas.out, "coloc.gwas_cat_chr3.txt")

