rm(list = ls(all=TRUE))

#code to illustrate co-expression QTL analysis and visualisation for IFN-gamma (24 hours) stimulated monocytes

#read in genotype data. n.b. genotypes have been randomly shuffled - and so results here will not recapiulate study results - code is provided to illustrated the analysis approach/pipeline.
#Individual-level genotypes are available at the European Genome-Phenome Archive (EGA) with accession ID EGAS00000000109. 
geno.ifn.qc <- read.table("./data/IFN_genotypes_shuffled.txt", header = T, row.names = 1)
#read in covariate-corrected cytokine secretion from IFN-gamma (24 hours) stimulated monocytes
cyto.ifn.t.work <- read.table("./data/IFN_cytokine_secretion_covar_corrected.txt", header = T, row.names = 1)
#read in covariate-corrected RNA expression from IFN-gamma (24 hours) stimulated monocytes
resid_rna.ifn.t.cyt <- read.table("./data/IFN_RNA_expn_covar_corrected.txt", header = T, row.names = 1)
#read in paired peak eSNP and eQTL gene
ifn.pairs <- read.table("./data/ifn.fdr_sig.1e5.txt", header = T)

#check cytokine and RNA sample IDs in same order
all(colnames(cyto.ifn.t.work)==colnames(resid_rna.ifn.t.cyt))

#check genotype sample IDs in same order
all(colnames(cyto.ifn.t.work)==colnames(geno.ifn.qc))

#add in RNA expression data for the first eQTL gene
cyto.rna.ifn <- rbind(cyto.ifn.t.work, resid_rna.ifn.t.cyt[paste0("X",ifn.pairs$probe[1]),])
rownames(cyto.rna.ifn)[29] <- "RNA"

#transpose data frame
cyto.rna.ifn.t <- data.frame(t(cyto.rna.ifn))

#extract genotypes for first eSNP:eGene pair
rsid <- data.frame(t(geno.ifn.qc[ifn.pairs$rsid[1],]))

#run co-expression analysis for first eSNP:eGene pair across all 28 cytokines - extracting interaction p-value
p <- c()
for (j in (c(1:28))){
  p[j] <- summary(lm(cyto.rna.ifn.t$RNA~round(rsid[,1],0)*cyto.rna.ifn.t[,j]))$coef[4,4]
}

#collect output
out <- cbind(colnames(cyto.rna.ifn.t)[c(1:28)], ifn.pairs$probe[1], p)
total <- out

#recapitulate analysis across all eSNP:eGene pairs - excluding instances where the peak eSNP has fewer that 5 minor alleles
for (i in c(2:3319)){
  #add ith eGene RNA expression to cytokine data
  cyto.rna.ifn <- rbind(cyto.ifn.t.work, resid_rna.ifn.t.cyt[paste0("X",ifn.pairs$probe[i]),])
  rownames(cyto.rna.ifn)[29] <- "RNA"
  #transpose data frame
  cyto.rna.ifn.t <- data.frame(t(cyto.rna.ifn))
  #extract genotypes for ith eSNP:eGene pair
  rsid <- data.frame(t(geno.ifn.qc[ifn.pairs$rsid[i],]))
  
  #skip where SNP minor allele count <5
  if(length(which(round(rsid[,1],0)==0))<5) next
  if(length(which(round(rsid[,1],0)==2))<5) next
  
  #run co-expression analysis for ith eSNP:eGene pair across all 28 cytokines - extracting interaction p-value
  p <- c()
  for (j in (c(1:28))){
    p[j] <- summary(lm(cyto.rna.ifn.t$RNA~round(rsid[,1],0)*cyto.rna.ifn.t[,j]))$coef[4,4]
  }
  
  colnames(cyto.rna.ifn.t)[c(1:28)]
  #collect output
  out <- cbind(colnames(cyto.rna.ifn.t)[c(1:28)], ifn.pairs$probe[i], p)
  #add output to total
  total <- rbind(total, out)
  print(i)
}

#apply FDR
total.ifn <- data.frame(total)
total.ifn$p <- as.numeric(total.ifn$p)
total.ifn$fdr <- p.adjust(total.ifn$p, method = "fdr")

#output significant (FDR<0.05) results
write.table(total.ifn[which(total.ifn$fdr<0.05),], "ifn_coexpn_sig.txt", col.names = T, row.names= F, quote = F)

#then plot out interaction plots, e.g. for PLTP interaction (again, genotypes are shuffled and so code illustrates pipeline but will not recapitulate results)

#first add PLTP expression (probe = 7040035) to cytokine secretion data
cyto.rna.ifn <- rbind(cyto.ifn.t.work, resid_rna.ifn.t.cyt[paste0("X",ifn.pairs$probe[which(ifn.pairs$probe=="7040035")]),])
rownames(cyto.rna.ifn)[29] <- "RNA"

#transpose dataframe
cyto.rna.ifn.t <- data.frame(t(cyto.rna.ifn))

#extract eSNP for PLTP expression (rs7270354 - 12:113357193)
rsid <- data.frame(t(geno.ifn.qc[which(rownames(geno.ifn.qc)=="12:113357193"),]))

#read in colour palette
cols <- brewer.pal(9,"Greens") 
cols2 <- brewer.pal(8,"Set2") 
cols1 <- brewer.pal(8,"Accent") 
cols3 <- brewer.pal(12, "Paired")

#constract dataframe of RNA expression, cytokine secretion and genotype for each significant co-expression QTL
for.int.plot <- data.frame(rbind(cbind(cyto.rna.ifn.t$RNA, round(rsid[,1],0), cyto.rna.ifn.t[,3], "GRO-a"),
                                 cbind(cyto.rna.ifn.t$RNA, round(rsid[,1],0), cyto.rna.ifn.t[,14], "IL-6"),
                                 cbind(cyto.rna.ifn.t$RNA, round(rsid[,1],0), cyto.rna.ifn.t[,5], "IL-10"),
                                 cbind(cyto.rna.ifn.t$RNA, round(rsid[,1],0), cyto.rna.ifn.t[,18], "IL-21"),
                                 cbind(cyto.rna.ifn.t$RNA, round(rsid[,1],0), cyto.rna.ifn.t[,23], "PIGF1"),
                                 cbind(cyto.rna.ifn.t$RNA, round(rsid[,1],0), cyto.rna.ifn.t[,25], "VEGF-A")))
colnames(for.int.plot) <- c("PLTP", "rs7270354", "cytokine", "cyto")

#extract p-values for labels
p_labels.int = data.frame(cyto = c("GRO-a", "IL-6", "IL-10", "IL-21", "PIGF1", "VEGF-A"), 
                          label = c("italic(P)==2%*%10^-7", "italic(P)==7.4%*%10^-12", "italic(P)==1.9%*%10^-6", "italic(P)==5.8%*%10^-9", "italic(P)==9.0%*%10^-9", "italic(P)==5.7%*%10^-7"))
p_labels.int$cyto <- factor(p_labels.int$cyto, levels = c("GRO-a", "IL-6", "IL-10", "IL-21", "PIGF1", "VEGF-A"))

for.int.plot$PLTP <- as.numeric(for.int.plot$PLTP)
for.int.plot$cytokine <- as.numeric(for.int.plot$cytokine)

#re-label genotypes
for.int.plot$rs7270354[which(for.int.plot$rs7270354==0)] <- "GG"
for.int.plot$rs7270354[which(for.int.plot$rs7270354==1)] <- "GA"
for.int.plot$rs7270354[which(for.int.plot$rs7270354==2)] <- "AA"
for.int.plot$rs7270354 <- factor(for.int.plot$rs7270354, levels = c("GG", "GA", "AA"))

#plot facetted interaction plots across 6 cytokines
int.pltp = ggplot(for.int.plot, aes(x=cytokine, y=PLTP)) +
  geom_point(size = 0.5) +
  geom_smooth(method = 'lm', se = FALSE) +
  aes(fill = rs7270354, col=rs7270354) + scale_fill_manual(values = cols1[c(1,2,3)]) + scale_colour_manual(values = cols1[c(1,2,3)]) +
  theme_bw() + theme(axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) +
  ylab(expression(italic(PLTP)~ RNA~ expression)) + xlab("cytokine secretion") +
facet_wrap( ~ cyto, ncol=3) + #ylim(NA, 1) +
  geom_text(x=0.5, y=1.1, aes(label=label), data=p_labels.int, parse=TRUE, inherit.aes=F, size = 4) +
  theme(strip.background = element_rect(color="black", fill=cols3[2], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
  
int.pltp

ggplot2::ggsave("pltp_int.jpg",
                width = 5,
                height = 5,
                dpi = 300)



