fancy_scientific_b = function(l) {
  l = as.character(l)
  l = unlist(strsplit(l, "e"))
  if(nchar(l[1]) > 4) l[1] = round(as.numeric(l[1]), 2)
  l = paste0(l[1], "e", l[2])
  l = format(l, scientific = TRUE)
  l = gsub("\\+", "", l)
  l = gsub("^(.*)e", "'\\1'e", l)
  l = gsub("e", "%*%10^", l)
  if(substr(l, 1,3) == "'1'") l = substr(l, 7, nchar(l))
  l
}
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(Hmisc)
library(magick)
library(readxl)
library(ggsci)
library(ggplot2)
mypal = pal_npg("nrc", alpha = 1)(9)
mypal[1] = "#FC4100"
mypal[4] = "#FFC55A"
mypal[7] = "#2C4E80"

#Panels A to C

A = image_read_pdf("BioRender_plots.pdf", pages = 1)
A = image_ggplot(A)
A = A + labs(tag = "A", face = "bold") + theme(plot.tag = element_text(face = "bold"))
B = image_read_pdf("BioRender_plots.pdf", pages = 2)
B = image_ggplot(B)
B = B + labs(tag = "B", face = "bold") + theme(plot.tag = element_text(face = "bold"))
C = image_read_pdf("BioRender_plots.pdf", pages = 3)
C = image_ggplot(C)
C = C + labs(tag = "C", face = "bold") + theme(plot.tag = element_text(face = "bold"))

#Panel D

promiscuity = readRDS("promiscuity.RDS")
promiscuity$Allele.name = sapply(promiscuity$Allele.name, FUN = function(x) {
  temp = unlist(strsplit(as.character(x), "/"))
  temp = sapply(temp, FUN = function(y) {
    paste0(substr(y, 1, 4), "\\*", substr(y, 5, 6), ":", substr(y, 7, nchar(y)))
  })
  if(length(temp) == 1) temp else paste0(temp[1], "/", temp[2])
})

common = readLines("hla_ref_set.class_ii.txt") #downloaded from https://help.iedb.org/hc/en-us/article_attachments/360047583591
common = gsub("HLA-", "", common)
common = gsub("\\*", "\\\\*", common)
common = c(common, "DRB1\\*04:03", "DRB1\\*15:02")

promiscuity$Allele.name[promiscuity$Allele.name %in% common] = paste0("**", promiscuity$Allele.name[promiscuity$Allele.name %in% common], "**")
promiscuity$Allele.name = factor(promiscuity$Allele.name, levels = promiscuity$Allele.name[order(promiscuity$Allele.promiscuity)])

D = ggplot(promiscuity, aes(x = Allele.name, y = Allele.promiscuity, fill = HLA.gene)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  scale_fill_manual(values = substr(mypal[c(1,4,7)], 1, 7)) + 
  theme(legend.position = c(0.05, 0.8), legend.title = element_blank()) + 
  ylab("Allele promiscuity") + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  theme(plot.tag = element_text(face = "bold")) +
  labs(tag = "D", face = "bold") +
  theme(axis.text.x = ggtext::element_markdown())

#Panel E

E = ggplot(promiscuity, aes(x=HLA.gene, y=Allele.promiscuity)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = HLA.gene), position=position_jitter(0.2), size = 2.5) + 
  scale_color_manual(values = substr(mypal[c(1,4,7)], 1, 7)) +
  theme_classic() +
  annotate(geom = "text", x = 1, y = 4, label = "italic(p) == 0.12", parse = T) +
  ylim(1.5, 4.3) +
  ylab("Allele promiscuity") +
  xlab("HLA type") +
  theme(legend.position = "n") +
  labs(tag = "E", face = "bold") +
  theme(plot.tag = element_text(face = "bold"))

#Panel F

load("aav_diversity_lowerthan0.RData")
aav = read_xlsx("Table_1_The HLA class-II immunopeptidomes of AAV capsids proteins.xlsx", sheet = 1)
alleles = table(aav$`Predicted HLA binding`)
cores = unique(aav[,c(9:10)])
seqnum = table(cores$`Predicted HLA binding`)
diversity_median = diversity_median[seqnum[names(diversity_median)] > 25 & grepl("DQ", names(diversity_median))]
names(diversity_median) = gsub("HLA-", "", names(diversity_median))
names(diversity_median) = gsub("_", "", names(diversity_median))
names(diversity_median) = gsub("-", "/", names(diversity_median))
df = data.frame(diversity = diversity_median, promiscuity = promiscuity$Allele.promiscuity[match(names(diversity_median), rownames(promiscuity))])
correl = rcorr(df$promiscuity, df$diversity, type = "spearman")
rh = round(correl$r[1,2], 2)
p = round(correl$P[1,2], 3)
label_spearman = paste0("rho: '", rh, "; italic(p) = '", as.character(p))
df = df[!is.na(df$promiscuity),]
f = ggplot(df, aes(x = promiscuity, y = diversity)) + 
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", color = "black", linetype = "dashed", fill = "lightgrey", linewidth = 0.75) +
  geom_point(colour =  "#FFC55A", size = 2) +
  ylab("Fraction of highly divergent sequence pairs") +
  xlab("HLA-DQ promiscuity") +
  theme_classic() +
  theme(legend.position = "n") +
  labs(tag = "F", face = "bold") +
  annotate(geom = "text", x = 2.4, y = 0.83, label = "paste(rho, ' = 0.83; ', italic(p), ' = 0.042')", parse = T)+
  theme(plot.tag = element_text(face = "bold"))

#Panel G

load("aav_diversity_lowerthan0.RData")
aav = read_xlsx("Table_1_The HLA class-II immunopeptidomes of AAV capsids proteins.xlsx", sheet = 1)
alleles = table(aav$`Predicted HLA binding`)
cores = unique(aav[,c(9:10)])
seqnum = table(cores$`Predicted HLA binding`)
diversity_median = diversity_median[seqnum[names(diversity_median)] > 25 & grepl("DRB1", names(diversity_median))]
names(diversity_median) = gsub("HLA-", "", names(diversity_median))
names(diversity_median) = gsub("_", "", names(diversity_median))
names(diversity_median) = gsub("-", "/", names(diversity_median))
df = data.frame(diversity = diversity_median, promiscuity = promiscuity$Allele.promiscuity[match(names(diversity_median), rownames(promiscuity))])
library(Hmisc)
correl = rcorr(df$promiscuity, df$diversity, type = "spearman")
rh = round(correl$r[1,2], 2)
p = round(correl$P[1,2], 3)
label_spearman = paste0("rho: '", rh, "; italic(p) = '", as.character(p))
G = ggplot(df, aes(x = promiscuity, y = diversity)) + 
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", color = "black", linetype = "dashed", fill = "lightgrey", linewidth = 0.75) +
  geom_point(colour = "#00215E", size = 2) +
  ylab("Fraction of highly divergent sequence pairs") +
  xlab("HLA-DRB1 allele promiscuity") +
  theme_classic() +
  theme(legend.position = "n") +
  labs(tag = "G", face = "bold") +
  annotate(geom = "text", x = 2.57, y = 0.8, label = "paste(rho, ': 0.9; ', italic(p), ' = 0.037')", parse = T)+
  theme(plot.tag = element_text(face = "bold"))

#Panels H and I

load("df_bound.Rdata")
correl = rcorr(df$Allele.promiscuity, df$score_iedb, type = "spearman")
rh = round(correl$r[1,2], 2)
p = fancy_scientific_b(correl$P[1,2])
label_spearman = paste0("rho: '", rh, "; p = ", substr(p, 2, nchar(p)))
H = ggplot(df, aes(x = Allele.promiscuity, y = score_iedb)) + 
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", color = "black", linetype = "dashed", fill = "lightgrey", linewidth = 0.75) +
  geom_point(aes(colour = HLA.gene)) +
  ylab("Median binding score for IEDB peps.") +
  xlab("Allele promiscuity") +
  theme_classic() +
  scale_colour_manual(values = substr(mypal[c(1,4,7)], 1, 7)) +
  theme(legend.position = c(0.85,0.17), legend.title = element_blank()) +
  labs(tag = "H", face = "bold") +
  scale_y_log10(labels = scales::number_format(accuracy = 0.001)) +
  annotate(geom = "text", x = 2.1, y = 0.019, label = label_spearman, parse = TRUE)+
  theme(plot.tag = element_text(face = "bold"))

correl = rcorr(df$Allele.promiscuity, df$score_vfdb, type = "spearman")
rh = round(correl$r[1,2], 2)
p = fancy_scientific_b(correl$P[1,2])
label_spearman = paste0("rho: '", rh, "; p = ", substr(p, 2, nchar(p)))

I = ggplot(df, aes(x = Allele.promiscuity, y = score_vfdb)) + 
  stat_smooth(method = "lm", 
              formula = y ~ x, 
              geom = "smooth", color = "black", linetype = "dashed", fill = "lightgrey", linewidth = 0.75) +
  geom_point(aes(colour = HLA.gene)) +
  ylab("Median binding score for VFDB peps.") +
  xlab("Allele promiscuity") +
  theme_classic() +
  theme(legend.position = c(0.85,0.17), legend.title = element_blank()) +
  scale_colour_manual(values = substr(mypal[c(1,4,7)], 1, 7)) +
  labs(tag = "I", face = "bold") +
  scale_y_log10(labels = scales::number_format(accuracy = 0.001)) +
  annotate(geom = "text", x = 2.1, y = 0.031, label = label_spearman, parse = TRUE) +
  theme(plot.tag = element_text(face = "bold"))


layout = rbind(c(1,1,1,1,1,2,2,2,2,2,2,3,3,3),
               c(1,1,1,1,1,2,2,2,2,2,2,3,3,3),
               c(1,1,1,1,1,2,2,2,2,2,2,3,3,3),
               c(4,4,4,4,4,4,4,4,4,4,4,4,4,4),
               c(4,4,4,4,4,4,4,4,4,4,4,4,4,4),
               c(5,5,6,6,6,7,7,7,8,8,8,9,9,9),
               c(5,5,6,6,6,7,7,7,8,8,8,9,9,9))

figure_1 =  arrangeGrob(A,B,C,D,E,f,G,H,I, layout_matrix = layout)
ggsave(plot = figure_1, filename = "Figure_1.png", width = 38, height = 32, units = "cm", dpi = 300)

