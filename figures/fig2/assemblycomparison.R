library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(hrbrthemes)
library(GGally)
library(viridis)
library('stringr')
library(ggtext)
library(cowplot)
install.packages("magick")
library(magick)
install.packages('patchwork')
library('patchwork')
install.packages('ggpubr')
library(ggpubr)
theme_set(theme_pubr())


setwd("~/R_assemblycomparison/")
getwd()

# Need to locate input files (both .txt and .png in google drive) to working directory

##################### Panel a; Kmer duplication #####################

KDDat <- read.table("Kmer_Dupl.txt", sep = '\t')
KDDat$V1_reorder <- factor(KDDat$V1, levels=c("CLR", "HiFi-only", "HiFi-Hic", "HiFi-trio"))
KDDat_Bar <- ggplot(KDDat, aes(x=V1_reorder, y=V2)) + 
  geom_bar(stat="identity", fill="#666666", width = 0.8, alpha=1, color="black", size=0.4) + #, alpha=.8) + 
  scale_y_continuous(limits = c(0,2.2), expand = c(0, 0)) +
  geom_segment(aes(x = 0.6, xend = 1.4, y = 2.0, yend = 2.0), size = 0.4) +
  geom_segment(aes(x = 1.6, xend = 4.4, y = 2.0, yend = 2.0), size = 0.4) +
  annotate("text", x= 1, y= 2.13, label="CLR", size = 3) +
  annotate("text", x= 3, y= 2.13, label="HiFi", size = 3) +
  theme_classic() + 
  ylab("Proportion of k-mer (%)") +
  #ylab("Proportion of *k*-mer (%)") +
  theme(axis.title.y = element_text(size = 8)) +
  theme(axis.title.y = ggtext::element_markdown()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 8, angle=45, hjust= 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 0.4)) +
  theme(axis.ticks = element_line(colour = 'black', size = 0.4)) +
  theme(axis.title.x=element_blank()) 


KDDat_Bar
ggsave(file="KD_Bar.svg", plot=KDDat_Bar, width=1.8, height=2.1)


##################### Panel b; Kmer exp/collap #####################

KECDat <- read.table("Kmer_ExpCollap.txt", sep = '\t')
KECDat$V1_reorder <- factor(KECDat$V1, levels=c("CLR", "HiFi-only", "HiFi-Hic", "HiFi-trio"))
KECDat$V3_reorder <- factor(KECDat$V3, levels=c("Expansion", "Collapse"))
KECDat_Bar <- ggplot(KECDat, aes(x=V1_reorder, y = V2, fill=V3_reorder)) + 
  geom_bar(stat="identity", position = "dodge", width = 0.8, alpha=1, color="black", size=0.4) + #, alpha=.8) + 
  scale_y_continuous(limits = c(0,13), expand = c(0, 0)) +
  geom_segment(aes(x = 0.6, xend = 1.4, y = 11.7, yend = 11.7), size = 0.4) +
  geom_segment(aes(x = 1.6, xend = 4.4, y = 11.7, yend = 11.7), size = 0.4) +
  annotate("text", x= 1, y= 12.5, label="CLR", size = 3) +
  annotate("text", x= 3, y= 12.5, label="HiFi", size = 3) +
  theme_classic() + 
  scale_fill_grey(start=0.55, end=0.85, name=NULL) + 
  ylab("Proportion of k-mer (%)") +
  #ylab("Proportion of *k*-mer (%)") +
  theme(legend.position = "top",) +
  theme(legend.key.size = unit(0.4, 'cm')) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) +
  theme(axis.title.y = element_text(size = 8)) +
  theme(axis.title.y = ggtext::element_markdown()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 8, angle=45, hjust= 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 0.4)) +
  theme(axis.ticks = element_line(colour = 'black', size = 0.4)) +
  theme(axis.title.x=element_blank()) 


KECDat_Bar
ggsave(file="KEC_Bar.svg", plot=KECDat_Bar, width=2.5, height=2.57)



##################### Panel c; Rebinned merged #####################
Paternal <- read.table("Paternal.txt", sep = '\t',header=TRUE)
Paternal$X_reorder <- factor(Paternal$X, levels=c("Default", "Rebinned"))
Paternal

Pat <- ggplot(Paternal, aes(x=X_reorder, group = 1)) +
  geom_line(aes(y=KmerDuplication), color="#CC0000", size = 1) +
  geom_point(aes(y=KmerDuplication), color="#CC0000", size = 2) + 
  geom_line(aes(y=KmerCollapse^10 / 2.1), color="#003366", size = 1) +
  geom_point(aes(y=KmerCollapse^10/ 2.1), color="#003366", size = 2) + 
  scale_y_continuous(limits = c(0.825,0.873), "K-mer duplication (red)", sec.axis = sec_axis(~(.*2.1)^(1/10),  name = "K-mer collapse (blue)",  labels = number_format(accuracy = 0.0001), )) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(axis.title.y=element_text(size = 8)) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(panel.grid.minor = element_blank())

Symbol <- file.path("./Zev_male.png")
Symbol
Pat_Sym <- ggdraw() +
  draw_plot(Pat) +
  draw_image(Symbol, x = 1, y = 1, hjust = 1.8, vjust = 1.25, width = 0.28, height = 0.28) 

#Pat_Sym
#ggsave(file="PaternalKmer.svg", plot=Pat_Sym, width=2.1, height=1.4)



Maternal <- read.table("Maternal.txt", sep = '\t',header=TRUE)
Maternal$X_reorder <- factor(Maternal$X, levels=c("Default", "Rebinned"))
Maternal

Mat <- ggplot(Maternal, aes(x=X_reorder, group = 1)) +
  geom_line(aes(y=KmerDuplication), color="#CC0000", size = 1) +
  geom_point(aes(y=KmerDuplication), color="#CC0000", size = 2) + 
  geom_line(aes(y=KmerCollapse^100/130 ), color="#003366", size = 1) +
  geom_point(aes(y=KmerCollapse^100/130), color="#003366", size = 2) + 
  scale_y_continuous( labels = number_format(accuracy = 0.01), limits = c(3.65,4.45), "K-mer duplication (red)", sec.axis = sec_axis(~(.*130)^(1/100),  name = "K-mer collapse (blue)")) +
  theme_bw() +
  theme(panel.grid.major = element_blank()) +
  theme(axis.title.y=element_text(size = 8)) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.text.x = element_text(size = 8)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(panel.grid.minor = element_blank())

Mat
Symbol2 <- file.path("./Zev_female.png")
Symbol2
Mat_Sym <- ggdraw() +
  draw_plot(Mat) +
  draw_image(Symbol2, x = 1, y = 1, hjust = 1.8, vjust = 1.25, width = 0.28, height = 0.28) 

#Mat_Sym
#ggsave(file="MaternalKmer.svg", plot=Mat_Sym, width=2.1, height=1.4)

PatMat <- Pat_Sym / Mat_Sym

PatMat
ggsave(file="BinningIssueKmer.svg", plot=PatMat, width=2.5, height=3.2)


##################### Panel d; False duplication ##################### 

FDDat <- read.table("False_Duplication.txt", sep = '\t')
FDDat$V1_reorder <- factor(FDDat$V1, levels=c("CLR", "HiFi-only", "HiFi-Hic", "HiFi-trio"))
FDDat_Bar <- ggplot(FDDat, aes(x=V1_reorder, y=V2)) + 
  geom_bar(stat="identity", fill="#CC0000", width = 0.6, alpha=1, color="white", size=0) + #, alpha=.8) +
  geom_text(aes(label=V3, vjust = -0.8), size = 3) +
  scale_y_continuous(limits = c(0,1.7), expand = c(0, 0)) +
  geom_segment(aes(x = 0.7, xend = 1.3, y = 1.55, yend = 1.55), size = 0.4) +
  geom_segment(aes(x = 1.7, xend = 4.3, y = 1.55, yend = 1.55), size = 0.4) +
  annotate("text", x= 1, y= 1.65, label="CLR", size = 3) +
  annotate("text", x= 3, y= 1.65, label="HiFi", size = 3) +
  theme_classic() + 
  ylab("Proportion of false duplications (%)") +
  theme(axis.title.y = element_text(size = 8)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(size = 8, angle=45, hjust= 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 0.4)) +
  theme(axis.ticks = element_line(colour = 'black', size = 0.4)) +
  theme(axis.title.x=element_blank()) 


FDDat_Bar
ggsave(file="FD_Bar.svg", plot=FDDat_Bar, width=1.7, height=2.6)



#####################  Panel e; HEAT MAP #####################
### False Loss Heat map

FLDat <- read.table("False_Loss_HM.txt", sep = '\t')
#FLDat$V2_reorder <- factor(FLDat$V2, levels=c("HiFi-trio", "HiFi-Hic", "HiFi-only", "CLR"))
#FLDat$V1_reorder <- factor(FLDat$V1, levels=c("CLR", "HiFi-only", "HiFi-Hic", "HiFi-trio"))
FLDat$V2_reorder <- factor(FLDat$V2, levels=c("CLR", "HiFi-only", "HiFi-Hic", "HiFi-trio"))
FLDat$V1_reorder <- factor(FLDat$V1, levels=c("HiFi-trio", "HiFi-Hic", "HiFi-only", "CLR"))
colnames(FLDat)[7] <- "Mbp"
FLDat_HM <- ggplot(FLDat, aes(V2_reorder, V1_reorder, fill=Mbp)) + 
  scale_fill_distiller(direction = +1) +
  theme_bw() +
  theme(legend.key.size = unit(0.3, 'cm')) +
  theme(legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8)) +
  theme(legend.position = "right") +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title.y=element_blank()) +
  theme(axis.text.y = element_text(size = 8)) +
  theme(axis.text.x = element_text(size = 8, angle=45, hjust= 1)) +
  theme(legend.title = element_text(size= 8)) + 
  geom_tile()

#FLDat_HM
#ggsave(file="FL_HM.svg", plot=FLDat_HM, width=2.515, height=1.70)

### False Loss Bar

FLDat2 <- read.table("False_Loss.txt", sep = '\t')
FLDat2$V1_reorder <- factor(FLDat2$V1, levels=c("CLR", "HiFi-only", "HiFi-Hic", "HiFi-trio"))
FLDat_Bar <- ggplot(FLDat2, aes(x=V1_reorder, y=V3)) + 
  geom_bar(stat="identity", fill="dodgerblue4", width = 0.95, alpha=1) + #, alpha=.8) + 
  theme_classic() + 
  geom_text(aes(label=V4, vjust = -0.1), size = 3) +
  scale_y_continuous(limits = c(0,9.5), labels = scales::number_format(accuracy = 0.1,decimal.mark = '.')) +
  geom_segment(aes(x = 0.6, xend = 1.4, y = 8.5, yend = 8.5), size = 0.4) +
  geom_segment(aes(x = 1.6, xend = 4.4, y = 8.5, yend = 8.5), size = 0.4) +
  annotate("text", x= 1, y= 9.3, label="CLR", size = 3) +
  annotate("text", x= 3, y= 9.3, label="HiFi", size = 3) +
  ylab("Proportion of\nfalse losses (%)") +
  theme(axis.title.y = element_text(size = 8)) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank()) +
  #theme(axis.text.x = element_text(size = 8, angle=45, hjust= 1)) +
  theme(axis.text.y = element_text(size = 8)) +
  #theme(panel.border = element_blank()) +
  theme(axis.line = element_line(colour = 'black', size = 0.4)) +
  theme(axis.ticks = element_line(colour = 'black', size = 0.4)) +
  theme(axis.title.x=element_blank()) 

#FLDat_Bar
#ggsave(file="FL_Bar.svg", plot=FLDat_Bar, width=1.9, height=1.0)

FL <- FLDat_Bar / FLDat_HM
ggsave(file="FL_BarHM.svg", plot=FL, width=2.9, height=3.0)

