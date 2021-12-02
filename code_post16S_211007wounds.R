#!/usr/bin/env Rscript
#####################################################################
########## To run this code in scrum ################################
#####################################################################
# module load R
# Rscript code.R ~/inputDir/

#####################################################################
########## Modifications made to Taxonomy ###########################
### Escherichia/Shigella -> Escherichia
### Cutibacterium -> Propionibacterium
#####################################################################

#####################################################################
########## Processing Set −up #######################################
#####################################################################
#Rpath <- "~/R/x86_64-redhat-linux-gnu-library/3.6/"
########## Call libraries for use ##########
# Packages from CRAN
#library("Cairo") # allow use of Times font for graphics # version 1.5-10
library("ggplot2") # version 3.2.1
#library("ggpubr") # for stat_compare_means command = diversity calculations of plot_richness # version 0.2.3
library("reshape2") # version 1.4.3
library("data.table") # version 1.12.6
library("RColorBrewer") # design new color palette for data # version 1.1-2
library("scales") # for scientific notation in plots # version 1.0.0
#library("cowplot") # used to prepare final figuresMd # version 1.0.0
library("tidyr") # version 1.0.0
library("dplyr") # version 0.8.3
#library("pheatmap") # version 1.0.12
library("vegan") # version 2.5-6

# Packages from Bioconductor
library("phyloseq") # version 1.28.0
library("decontam") # identify contaminant ASVs # version 1.4.0
#library("pairwiseAdonis") # for PERMANOVA calculations # version 0.0.1
library("DESeq2") # The DESeq2 model internally corrects for library size, so transformed or normalized values such as counts scaled by library size should not be used as input. # version 1.24.0
library("dada2") # version 1.12.1
library("Biostrings") # 2.52.0
library("microbiome") # version 1.6.0

########## Call functions for use ##########
source("~/Masters/taxa_summary.R",local=TRUE) # load fast_melt function

########## Define color palettes ##########
pal.pairBiome<-c("#1c541c","#2f8f2f","#49c349","#bfeabf",
                 "#B80000","#F00000","#FF7777","#ffcccc",
                 "#000080","#0000cd","#8282ff","#cfcfff",
                 "#623800","#c47000","#ff9914","#ffddb1",
                 "#430059","#7d00a7","#cc32ff","#eebbff",
                 "#626200","#c4c400","#ffff14","#ffffb1",
                 "#005f6c","#00a4bb","#1ee3ff","#bbf7ff",
                 "#750063","#c400a5","#ff13da","#ffb0f3",
                 "#1a1a1a","#808080","#d9d9d9","#ffffff","#808080")
pal.pairMini<-c("#e6e6e6","#c0c0c0","#404040",
                "#A6CEE3","#1F78B4",
                "#B2DF8A","#33A02C",
                "#FB9A99","#E31A1C",
                "#FDBF6F","#FF7F00",
                "#CAB2D6","#6A3D9A",
                "#FFFF99","#eded09",
                "#9aebff","#01cdff"
                )
pal.CB<-c("#e69f00","#56b4e9","#009e73","#f0e442","#0072b2","#D55e00","#cc79a7","#999999","#000000","#ffffff") # colorblind-friendly color palette

### Starting Data 
setwd("~/Desktop/waniganR_emily211007/")
load("output/physeq_start.RData")
phyN.pruS
taxa_names(phyN.pruS)<-tax_table(phyN.pruS)[,"ASV"]
phyN.wound <- subset_samples(phyN.pruS,Project%in%c("Controls","DoDwound"))
phyN.wound <- subset_taxa(phyN.wound,taxa_sums(phyN.wound)>1)

### Edit with modified table
MAP = sample_data(data.frame(read.csv("table_map_xtend.csv", header = TRUE, row.names = 1, check.names=FALSE))) # reads csv file into data.frame with row names in column 1
phyR.wound <- phyloseq(otu_table(phyN.wound),tax_table(phyN.wound),MAP)

### Modify taxonomy for easier reading / medical terminology
df.tax <- as.data.frame(tax_table(phyR.wound))
df.tax$Genus <- with(df.tax,
                   ifelse(Genus=="Corynebacterium_1","Corynebacterium",
                          ifelse(Genus=="Escherichia/Shigella","Escherichia",
                                 ifelse(Genus=="Cutibacterium","Propionibacterium",
                                        Genus))))
phyR.wound <- phyloseq(otu_table(phyR.wound),tax_table(as.matrix(df.tax)),sample_data(data.frame(read.csv("table_map_xtend.csv", header = TRUE, row.names = 1, check.names=FALSE))))
  
#phyN.pruS <- subset_samples(phyN.pruS,!SpecimenType=="Fasia")
### Remove Animalia reads ###
phyN.animal <- subset_taxa(phyR.wound,Kingdom=="Animalia")
phyN.bacteria <- subset_taxa(phyR.wound,!Kingdom=="Animalia")

### Evaluate animal contamination ###
phyG.wound <- tax_glom(phyR.wound,taxrank="Genus")
taxa_names(phyG.wound) <- tax_table(phyG.wound)[,"Genus"]

df.nReads <- as.data.frame(cbind(sample_sums(phyN.bacteria),
                                 sample_sums(phyR.wound),
                                 otu_table(phyG.wound)
                                 ))
names(df.nReads)[1:2]=c("Bacteria","Total")
df.nReads$MousePer<- with(df.nReads, ifelse(Total=="0","0", Mus / Total))
df.nReads$HumanPer<- with(df.nReads, ifelse(Total=="0","0", Homo / Total))
df.nReads$BactPer<- with(df.nReads, ifelse(Total=="0","0", Bacteria / Total))

df.nReads<-cbind(df.nReads,
                 sample_data(phyR.wound)$DNAconc,
                 sample_data(phyR.wound)$ExtractDate,
                 sample_data(phyR.wound)$MiSeqPlate,
                 sample_data(phyR.wound)$RossID
                 )
names(df.nReads)[eval(length(names(df.nReads))-3):length(names(df.nReads))]=c("MiSeqPlate","ExtractDate","DNAconc","Pt")

write.csv(df.nReads,"nReads.csv")

ggplot(df.nReads, aes(x=MousePer, y=Bacteroides)) + geom_point()

### Divide bacterial counts by human counts
phy.homo <- subset_taxa(subset_samples(phyG.wound,otu_table(phyG.wound)>0),taxa_sums(phyG.wound)>0)
df.homo <- as.data.frame(otu_table(phy.homo))
phy.HOMO <- phyloseq(otu_table(df.homo/df.homo[,"Homo"],taxa_are_rows = FALSE),tax_table(phy.homo),MAP)
write.csv(df.homo/df.homo[,"Homo"],"perHomo.csv")

phy.homo11 <- subset_samples(phy.HOMO,RossID=="MB-11")
phy.homo11 <- subset_taxa(phy.homo11,taxa_sums(phy.homo11)>1e-3)
mdf.otherMB11 = psmelt(subset_samples(phy.homo11,!SpecimenType=="RemoteSwab"))
pBW.otherMB11 <- ggplot(mdf.otherMB11, aes(as.factor(CollectionDate), Abundance, color = Genus)) + 
  theme_bw() +
  geom_boxplot(outlier.shape=16,outlier.size=2, notch=FALSE) +
  facet_grid(StoredSampleType~SpecimenType, scales="free_x", space="free_x") +
  xlab("") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=pal.pairBiome)

phyA <- subset_taxa(phy.homo,Genus=="Mus")
plot_bar(subset_samples(phyA,Plate=="Ross210823A"),fill="Project")+
  theme_bw() +
  facet_grid(Column~Row,scales="free_x") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_blank(),panel.spacing = unit(0, "lines"))




### Keep samples with >100 bacterial reads ###
phyN.bact100 <- subset_samples(phyN.bacteria,sample_sums(phyN.bacteria)>2) 
#phyN.bact100 <-  phyN.bacteria
### Remove unique taxa ###
dt.bact100 = fast_melt(phyN.bact100) # make data table from phyloseq object (data table. physeq Rarefied)
prev.bact100 = dt.bact100[, list(Prevalence = sum(count > 1),TotalPer = sum(count),
                             MinCount = min(count), MaxCount = max(count)),
                      by = ASV] # make simple table listing 'ASV, Prevalence, TotalPer, MinCount, and MaxCount' (prevalence . physeq Rarefied)
ls.bact100 = prev.bact100[(Prevalence > 1 & TotalPer > 10), ASV] # Make list of ASVs presMdent in dataset at least once (list . PresMdent)
phyN.pruBact <- subset_taxa(phyN.bact100,ASV%in%c(ls.bact100)) # remove taxa not in ls.PresMd from phyloseq object (physeq Raw . pruned 2)

### Transform reads to percentage ###
phyT.bact <- transform_sample_counts(phyN.pruBact, function(x) x/sum(x)) # transform raw counts to fraction (physeq Transform . leech)

### Separate into "culturable" and "other" groups
sample_data(phyT.bact)$RossID <- factor(sample_data(phyT.bact)$RossID, levels = c("R-02","R-03","R-04","R-05","MB-06","MB-07","MB-08","MB-09","MB-10","MB-11","MB-12"))
sample_data(phyT.bact)$xxx <- paste(sample_data(phyT.bact)$DHMCspecimen,sample_data(phyT.bact)$CollectionDate,sample_data(phyT.bact)$StoredSampleType)
phyT.cultTax <- subset_taxa(phyT.bact,Genus%in%c("Acinetobacter","Aeromonas","Arcanobacterium",
                                                 "Bacillus","Bacteroides",
                                                 "Chryseobacterium","Citrobacter","Clostridium","Corynebacterium","Cutibacterium",
                                                 "Enterobacter","Enterococcus","Escherichia/Shigella","Escherichia",
                                                 "Fusobacterium","Graunilicitella","Hafnia","Klebsiella",
                                                 "Lactobacillus","Listeria","Micrococcus","Morganella",
                                                 "Peptostreptococcus","Prevotella","Propionibacterium","Proteus","Providencia","Pseudomonas",
                                                 "Serratia","Sphingomonas","Stenotrophomonas","Staphylococcus","Streptococcus"))
phyT.otherTax <- subset_taxa(phyT.bact,!Genus%in%c("Acinetobacter","Aeromonas","Arcanobacterium",
                                                   "Bacillus","Bacteroides",
                                                   "Chryseobacterium","Citrobacter","Clostridium","Corynebacterium","Cutibacterium",
                                                   "Enterobacter","Enterococcus","Escherichia/Shigella","Escherichia",
                                                   "Fusobacterium","Graunilicitella","Hafnia","Klebsiella",
                                                   "Lactobacillus","Listeria","Micrococcus","Morganella",
                                                   "Peptostreptococcus","Prevotella","Propionibacterium","Proteus","Providencia","Pseudomonas",
                                                   "Serratia","Sphingomonas","Stenotrophomonas","Staphylococcus","Streptococcus"))

### Evaluation of MB-11
phyT.cultMB11 <- subset_samples(phyT.cultTax,RossID == "MB-11")
phyT.cultMB11 <- subset_taxa(phyT.cultMB11,taxa_sums(phyT.cultMB11)>0)
phyG.cultMB11 <- tax_glom(phyT.cultMB11,taxrank="Genus")
boxplot_abundance(phyG.cultMB11,x="CollectionDate",y="ASV31")

mdf.cultMB11 = psmelt(subset_samples(phyG.cultMB11,!SpecimenType=="RemoteSwab"))
pBW.cultMB11 <- ggplot(mdf.cultMB11, aes(as.factor(CollectionDate), Abundance, color = Genus)) + 
  theme_bw() +
  geom_boxplot(outlier.shape=16,outlier.size=2, notch=FALSE) +
  facet_grid(StoredSampleType~SpecimenType, scales="free_x", space="free_x") +
  xlab("") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=pal.CB)
ggsave(pBW.cultMB11,filename="wounds_plotBW_MB11cultTax.pdf", dpi="retina",width=15,height=8,units="in") # Save figure to .pdf file

phyT.otherMB11 <- subset_samples(phyT.otherTax,RossID == "MB-11")
phyT.otherMB11 <- subset_taxa(phyT.otherMB11,taxa_sums(phyT.otherMB11)>0)
phyG.otherMB11 <- tax_glom(phyT.otherMB11,taxrank="Genus")
boxplot_abundance(phyG.otherMB11,x="CollectionDate",y="ASV31")

mdf.otherMB11 = psmelt(subset_samples(phyG.otherMB11,!SpecimenType=="RemoteSwab"))
pBW.otherMB11 <- ggplot(mdf.otherMB11, aes(as.factor(CollectionDate), Abundance, color = Genus)) + 
  theme_bw() +
  geom_boxplot(outlier.shape=16,outlier.size=2, notch=FALSE) +
  facet_grid(StoredSampleType~SpecimenType, scales="free_x", space="free_x") +
  xlab("") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=pal.pairMini)
ggsave(pBW.otherMB11,filename="wounds_plotBW_MB11otherTax.pdf", dpi="retina",width=15,height=8,units="in") # Save figure to .pdf file

### Evaluation of MB-10
phyT.cultMB10 <- subset_samples(phyT.cultTax,RossID == "MB-10")
phyT.cultMB10 <- subset_taxa(phyT.cultMB10,taxa_sums(phyT.cultMB10)>0)
phyG.cultMB10 <- tax_glom(phyT.cultMB10,taxrank="Genus")
mdf.cultMB10 = psmelt(subset_samples(phyG.cultMB10,!SpecimenType=="RemoteSwab"))

pBW.cultMB10 <- ggplot(mdf.cultMB10, aes(as.factor(CollectionDate), Abundance, color = Genus)) + 
  theme_bw() +
  geom_boxplot(outlier.shape=16,outlier.size=2, notch=FALSE) +
  facet_grid(StoredSampleType~SpecimenType, scales="free_x", space="free_x") +
  xlab("") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=pal.pairMini)
ggsave(pBW.cultMB10,filename="wounds_plotBW_MB10cultTax.pdf", dpi="retina",width=15,height=8,units="in") # Save figure to .pdf file

phyT.otherMB10 <- subset_samples(phyT.otherTax,RossID == "MB-10")
phyT.otherMB10 <- subset_taxa(phyT.otherMB10,taxa_sums(phyT.otherMB10)>0)
phyG.otherMB10 <- tax_glom(phyT.otherMB10,taxrank="Genus")
mdf.otherMB10 = psmelt(subset_samples(phyG.otherMB10,!SpecimenType=="RemoteSwab"))

pBW.otherMB10 <- ggplot(mdf.otherMB10, aes(as.factor(CollectionDate), Abundance, color = Genus)) + 
  theme_bw() +
  geom_boxplot(outlier.shape=16,outlier.size=2, notch=FALSE) +
  facet_grid(StoredSampleType~SpecimenType, scales="free_x", space="free_x") +
  xlab("") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_color_manual(values=pal.pairBiome)
ggsave(pBW.otherMB10,filename="wounds_plotBW_MB10otherTax.pdf", dpi="retina",width=15,height=8,units="in") # Save figure to .pdf file





#mdf.cultMB11[which(mdf.cultMB11$Abundance>0),]




phyG.bact2 <- tax_glom(phyN.pruBact,"Genus") 
phyG.bact2 <- subset_samples(phyG.bact2,sample_sums(phyG.bact2)>10)
phyG.bact <- subset_taxa(phyG.bact2,taxa_names(phyG.bact2)%in%names(sort(taxa_sums(phyG.bact2),decreasing=TRUE)[1:50]))


phyT.bact <- subset_samples(phyT.bact,!RossID%in%c("none","unk"))


#sample_data(phyT.cultTax)$RossID <- factor(sample_data(phyT.cultTax)$RossID, levels = c("R-02","R-03","R-04","R-05","MB-06","MB-07","MB-08","MB-09","MB-10","MB-11","MB-12"))
pBar.cultTax <- plot_bar(phyT.cultTax, x="xxx",fill="Genus") + 
  theme_bw() +
  #scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(0.1),expression(1),expression(10),""),limits=c(.001,1)) +
  xlab("") +
  facet_grid(SpecimenType~RossID, scales="free_x", space="free_x") +
  ylab("Abundance (%)") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.pairBiome)
pBar.cultTax
ggsave(pBar.cultTax,filename="wounds_plotBar_cultTax.pdf", dpi="retina",width=15,height=8,units="in") # Save figure to .pdf file

pBar.otherTax <- plot_bar(phyT.otherTax, x="xxx",fill="Genus") + 
  theme_bw() +
  #scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(0.1),expression(1),expression(10),""),limits=c(.001,1)) +
  xlab("") +
  facet_grid(SpecimenType~RossID, scales="free_x", space="free_x") +
  ylab("Abundance (%)") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.pairMini)
#pBar.otherTax
ggsave(pBar.otherTax,filename="wounds_plotBar_otherTax.pdf", dpi="retina",width=15,height=8,units="in") # Save figure to .pdf file


df.otu <- t(as.data.frame(otu_table(phyT.cultTax)))
df.otuXtend <- cbind(tax_table(phyT.cultTax),df.otu)
rownames(df.otuXtend) <- df.otuXtend[,"ASV"]
write.csv(df.otuXtend,"output/table_otuXtend_cultTax_count.csv")

aaa<-subset_taxa(phyG.bact,taxa_names(phyG.bact)%in%taxa_names(phyT.cultTax))




### Evaluation of MB-10
phyT.cultMB10 <- subset_samples(phyT.cultTax,RossID == "MB-10")
phyT.cultMB10 <- subset_taxa(phyT.cultMB10,taxa_sums(phyT.cultMB10)>1)
phyT.cultMB10 <- subset_samples(phyT.cultMB10,SpecimenType)
pBar.cultMB10 <- plot_bar(phyT.cultMB10, x="OriginalSample",fill="Genus") + 
  theme_bw() +
  #scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(0.1),expression(1),expression(10),""),limits=c(.001,1)) +
  xlab("") +
  facet_grid(SpecimenType+StoredSampleType~CollectionDate, scales="free_x", space="free_x") +
  ylab("Abundance (%)") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.CB)
pBar.cultMB10

phyT.otherMB10 <- subset_samples(phyT.otherTax,RossID == "MB-10")
#phyT.otherMB10 <- subset_samples(phyT.otherMB10,SpecimenType!)
phyT.otherMB10 <- subset_taxa(phyT.otherMB10,taxa_names(phyT.otherMB10)%in%names(sort(taxa_sums(phyT.otherMB10),decreasing=TRUE)[1:36]))

pBar.otherMB10 <- plot_bar(phyT.otherMB10, x="OriginalSample",fill="Genus") + 
  theme_bw() +
  #scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(0.1),expression(1),expression(10),""),limits=c(.001,1)) +
  xlab("") +
  facet_grid(SpecimenType+StoredSampleType~CollectionDate, scales="free_x", space="free_x") +
  ylab("Abundance (%)") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.pairBiome)
pBar.otherMB10

### Evaluation of MB-11
phyT.cultMB11 <- subset_samples(phyT.cultTax,RossID == "MB-11")
phyT.cultMB11 <- subset_taxa(phyT.cultMB11,taxa_sums(phyT.cultMB11)>1)
phyT.cultMB11 <- subset_samples(phyT.cultMB11,SpecimenType)
pBar.cultMB11 <- plot_bar(phyT.cultMB11, x="OriginalSample",fill="Genus") + 
  theme_bw() +
  #scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(0.1),expression(1),expression(10),""),limits=c(.001,1)) +
  xlab("") +
  facet_grid(SpecimenType+StoredSampleType~CollectionDate, scales="free_x", space="free_x") +
  ylab("Abundance (%)") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.CB)
pBar.cultMB11

phyT.otherMB11 <- subset_samples(phyT.otherTax,RossID == "MB-11")
#phyT.otherMB11 <- subset_samples(phyT.otherMB11,SpecimenType!)
phyT.otherMB11 <- subset_taxa(phyT.otherMB11,taxa_names(phyT.otherMB11)%in%names(sort(taxa_sums(phyT.otherMB11),decreasing=TRUE)[1:36]))

pBar.otherMB11 <- plot_bar(phyT.otherMB11, x="OriginalSample",fill="Genus") + 
  theme_bw() +
  #scale_y_log10(breaks=c(.001,.01,.1,1),labels=c(expression(0.1),expression(1),expression(10),""),limits=c(.001,1)) +
  xlab("") +
  facet_grid(SpecimenType+StoredSampleType~CollectionDate, scales="free_x", space="free_x") +
  ylab("Abundance (%)") +
  theme(text=element_text(size=12),strip.text.y = element_text(size=8),axis.text.x = element_text(size=8,angle = 45, hjust=1),panel.spacing = unit(0, "lines")) +
  scale_fill_manual(values=pal.pairBiome)
pBar.otherMB11

sdt = data.table(as(sample_data(phyT.cultMB11), "data.frame"),
                 TotalReads = sample_sums(mp), keep.rownames = TRUE)






mdf.dplyr = psmelt(phyT.cultMB11)


ggplot(mdf.dplyr, aes(CollectionDate, Abundance, color = Genus)) + 
  geom_point(size = 5)  +
  geom_boxplot(outlier.colour="black", outlier.shape=16,
               outlier.size=2, notch=FALSE)
  #geom_smooth(method = lm) +
  #ggtitle("Sequencing Depth vs. Time") +
  #scale_y_log10() +
  #facet_grid(StoredSampleType ~ .)









