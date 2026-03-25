fancy_scientific_b = function(l) {
  # turn in to character string in scientific notation
  l = as.character(l)
  l = unlist(strsplit(l, "e"))
  if(nchar(l[1]) > 4) l[1] = round(as.numeric(l[1]), 2)
  l = paste0(l[1], "e", l[2])
  l = format(l, scientific = TRUE)
  l = gsub("\\+", "", l)
  # quote the part before the exponent to keep all the digits
  l = gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l = gsub("e", "%*%10^", l)
  if(substr(l, 1,3) == "'1'") l = substr(l, 7, nchar(l))
  # return this as an expression
  l
}
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(Hmisc)
library(magick)
library(readxl)

A = image_read_pdf("../../results/plots/FA_BioRender_plots.pdf", pages = 1)
A = image_ggplot(A)
A = A + labs(tag = "A", face = "bold") + theme(plot.tag = element_text(face = "bold"))
B = image_read_pdf("../../results/plots/FA_BioRender_plots.pdf", pages = 2)
B = image_ggplot(B)
B = B + labs(tag = "B", face = "bold") + theme(plot.tag = element_text(face = "bold"))

load("../../data/FA_promiscuity_table_all12_new.RData")
promiscuities[promiscuities == "-"] = NA
promiscuity = promiscuities[,"gfeller_100"]
names(promiscuity) = rownames(promiscuities)
promiscuity = promiscuity[!is.na(promiscuity)]

promiscuity = data.frame(Allele.name = names(promiscuity), Allele.promiscuity = promiscuity, HLA.gene = NA)
promiscuity$HLA.gene[grepl("DRB1", promiscuity$Allele.name)] = "DRB1"
promiscuity$HLA.gene[grepl("DPB1", promiscuity$Allele.name)] = "DP"
promiscuity$HLA.gene[grepl("DQB1", promiscuity$Allele.name)] = "DQ"
promiscuity = promiscuity[!is.na(promiscuity$HLA.gene),]

library(ggsci)
library(ggplot2)
mypal = pal_npg("nrc", alpha = 1)(9)
mypal[1] = "#FC4100"
mypal[4] = "#FFC55A"
mypal[7] = "#2C4E80"


Bb = image_read_pdf("../../results/plots/FA_BioRender_plots.pdf", pages = 3)
Bb = image_ggplot(Bb)
Bb = Bb + labs(tag = "C", face = "bold") + theme(plot.tag = element_text(face = "bold"))

promiscuity$Allele.name = sapply(promiscuity$Allele.name, FUN = function(x) {
  temp = unlist(strsplit(as.character(x), "/"))
  temp = sapply(temp, FUN = function(y) {
    paste0(substr(y, 1, 4), "\\*", substr(y, 5, 6), ":", substr(y, 7, nchar(y)))
  })
  if(length(temp) == 1) temp else paste0(temp[1], "/", temp[2])
})

common = readLines("../../data/hla_ref_set.class_ii.txt symlink")
common = gsub("HLA-", "", common)
common = gsub("\\*", "\\\\*", common)
common = c(common, "DRB1\\*04:03", "DRB1\\*15:02")

promiscuity$Allele.name[promiscuity$Allele.name %in% common] = paste0("**", promiscuity$Allele.name[promiscuity$Allele.name %in% common], "**")

promiscuity$Allele.name = factor(promiscuity$Allele.name, levels = promiscuity$Allele.name[order(promiscuity$Allele.promiscuity)])



C = ggplot(promiscuity, aes(x = Allele.name, y = Allele.promiscuity, fill = HLA.gene)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  scale_fill_manual(values = substr(mypal[c(1,4,7)], 1, 7)) + 
  theme(legend.position = c(0.05, 0.8), legend.title = element_blank()) + 
  ylab("Allele promiscuity") + theme(axis.title.x = element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) + 
  theme(plot.tag = element_text(face = "bold")) +
  labs(tag = "D", face = "bold") +
  theme(axis.text.x = ggtext::element_markdown())

library(ggpubr)
D = ggplot(promiscuity, aes(x=HLA.gene, y=Allele.promiscuity)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(color = HLA.gene), position=position_jitter(0.2), size = 2.5) + 
  scale_color_manual(values = substr(mypal[c(1,4,7)], 1, 7)) +
  theme_classic() +
  annotate(geom = "text", x = 1, y = 4, label = parse(text = expression("~italic(p)~': 0.12'"))) +
  ylim(1.5, 4.3) +
  ylab("Allele promiscuity") +
  xlab("HLA type") +
  theme(legend.position = "n") +
  labs(tag = "E", face = "bold") +
  theme(plot.tag = element_text(face = "bold"))

load("../../data/MM_df_bound.Rdata")
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
  annotate(geom = "text", x = 2.1, y = 0.019, label = parse(text = label_spearman))+
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
  annotate(geom = "text", x = 2.1, y = 0.031, label = parse(text = label_spearman)) +
  theme(plot.tag = element_text(face = "bold"))

load("../../data/MM_aav_diversity_lowerthan0.RData")
aav = read_xlsx("../../data/Table_1_The HLA class-II immunopeptidomes of AAV capsids proteins.xlsx", sheet = 1)
alleles = table(aav$`Predicted HLA binding`)
cores = unique(aav[,c(9:10)])
seqnum = table(cores$`Predicted HLA binding`)
diversity_median = diversity_median[seqnum[names(diversity_median)] > 25 & grepl("DQ", names(diversity_median))]
names(diversity_median) = gsub("HLA-", "", names(diversity_median))
names(diversity_median) = gsub("_", "", names(diversity_median))
names(diversity_median) = gsub("-", "/", names(diversity_median))
df = data.frame(diversity = diversity_median, promiscuity = promiscuity$Allele.promiscuity[match(names(diversity_median), 
                                                                                                 rownames(promiscuity))])
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
  annotate(geom = "text", x = 2.4, y = 0.83, label = parse(text = expression("rho: '0.83; '~italic(p)~': 0.042'")))+
  theme(plot.tag = element_text(face = "bold"))


load("../../data/MM_aav_diversity_lowerthan0.RData")
aav = read_xlsx("../../data/Table_1_The HLA class-II immunopeptidomes of AAV capsids proteins.xlsx", sheet = 1)
alleles = table(aav$`Predicted HLA binding`)
cores = unique(aav[,c(9:10)])
seqnum = table(cores$`Predicted HLA binding`)
diversity_median = diversity_median[seqnum[names(diversity_median)] > 25 & grepl("DRB1", names(diversity_median))]
names(diversity_median) = gsub("HLA-", "", names(diversity_median))
names(diversity_median) = gsub("_", "", names(diversity_median))
names(diversity_median) = gsub("-", "/", names(diversity_median))
df = data.frame(diversity = diversity_median, promiscuity = promiscuity$Allele.promiscuity[match(names(diversity_median), 
                rownames(promiscuity))])
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
  annotate(geom = "text", x = 2.57, y = 0.8, label = parse(text = expression("rho: '0.9; '~italic(p)~': 0.037'")))+
  theme(plot.tag = element_text(face = "bold"))

load("../../data/FA_HLA_motif_heatmap_data_0925.RData")
load("../../data/FA_HLA_motif_heatmap_annotation.RData")
library(pheatmap)
breaks = seq(min(positions_rowordered), max(positions_rowordered), length.out = 100)
positions_rowordered = positions_rowordered[grepl("DRB1", rownames(positions_rowordered)) | grepl("DP", rownames(positions_rowordered)) |
                                        grepl("DQ", rownames(positions_rowordered)),]

hist(log2(unlist(positions_rowordered)))

rownames(positions_rowordered) = gsub("-", "/", rownames(positions_rowordered))
positions_rowordered = positions_rowordered[,order(apply(positions_rowordered, 2, max))]

J = pheatmap(t(positions_rowordered[,apply(positions_rowordered, 2, max)>0.4]), cluster_cols = F, 
         cluster_rows = F, color = colorRampPalette(c("#2C4E80", "red", "yellow"), bias = 2)(length(breaks) - 1), 
         breaks = breaks)$gtable

# J = J + labs(tag = "A", face = "bold")

test_this = t(positions_rowordered[,apply(positions_rowordered, 2, max)>0.4])
apply(test_this, 1, FUN = function(x) {
  
  fisher.test(rbind(c(sum(x[1:32] > 2), sum(x[1:32] <= 2)),
                    c(sum(x[33:64] > 2), sum(x[33:64] <= 2))))$p.value
})

library(dplyr)

gfeller = read.csv("../../data/TF_Gfeller_full.csv")
gfeller = gfeller[,-1]
gfeller <- gfeller %>%
  dplyr::select(Allele, Core) %>%
  distinct()

alleles = readLines("../../data/TF_HLA_II_alleles")
load("../../data/TF_aa_prevalence_in_human_proteins.RData")

library(universalmotif)

motifs_full <- setNames(
  lapply(alleles, function(x) {
    seq_cor = gfeller$Core[gfeller$Allele == x]
    mot1 = create_motif(seq_cor, alphabet = "AA", type = "ICM", pseudocount = 1,  bkg = stat_all[[2]][,5])
    view_motifs(mot1, method = "KL", use.type = "ICM", min.height = 0, y.spacer = 0.001) + ylim(c(0,1.75))
  }),
  alleles
)

seq_cor = gfeller$Core[gfeller$Allele == "DRB11602"]
mot1 = create_motif(seq_cor, alphabet = "AA", type = "ICM", pseudocount = 1,  bkg = stat_all[[2]][,5])
# pheatmap(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif)
# hist(as.vector(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif))
klvals_drb11602 = as.vector(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif)
drb11602 = view_motifs(mot1, use.type = "ICM", min.height = 0, y.spacer = 0.001) + ylim(c(0,2.75))
drb11602 = drb11602 + 
  annotate(geom = "text", x = 5, y = 2.75, label = "HLA-DRB1*16:02", size = 4) +
  theme(axis.title=element_text(size=14)) + 
  labs(tag = "B", face = "bold") +
  theme(plot.tag = element_text(face = "bold", size = 12))



seq_cor = gfeller$Core[gfeller$Allele == "DRB10804"]
mot1 = create_motif(seq_cor, alphabet = "AA", type = "ICM", pseudocount = 1,  bkg = stat_all[[2]][,5])
# pheatmap(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif)
klvals_drb10804 = as.vector(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif)
drb10804 = view_motifs(mot1, use.type = "ICM", min.height = 0, y.spacer = 0.001) + ylim(c(0,2.75))
drb10804 = drb10804 + 
  annotate(geom = "text", x = 5, y = 2.75, label = "HLA-DRB1*08:04", size = 4) +
  theme(axis.title=element_text(size=14)) + 
  labs(tag = " ", face = "bold") +
  theme(plot.tag = element_text(face = "bold", size = 12))
df_dens = data.frame(Allele = c(rep("DRB1*16:02", 180), rep("DRB1*08:02", 180)), Values = c(klvals_drb11602, klvals_drb10804))
Kb = ggplot(df_dens, aes(x = Values, colour = Allele)) + geom_density(bw = 0.4) + theme_bw() + 
  scale_color_manual(values = substr(mypal[c(7,1)], 1, 7)) + xlab("log2(observed/background)") + ylab("Density") +
  theme(legend.position="none") + annotate(geom = "text", x = -6, y = 0.15, label = "HLA-DRB1*16:02", size = 3, color = "#FC4100") +
  annotate(geom = "text", x = -4, y = 0.3, label = "HLA-DRB1*08:04", size = 3, color = "#2C4E80") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(tag = "C", face = "bold") +
  theme(plot.tag = element_text(face = "bold", size = 12))

K = grid.arrange(drb11602, drb10804, nrow = 2)


load("../../data/FA_HLA_motif_max.RData")
plottable = data.frame(Max_motif = apply(positions_rowordered, 1, max), 
                       Promiscuity = promiscuity$Allele.promiscuity[match(rownames(positions_rowordered), rownames(promiscuity))])

plottable = plottable[grepl("DRB1", rownames(plottable)) | grepl("DP", rownames(plottable)) |
                                        grepl("DQ", rownames(plottable)),]
plottable = as.data.frame(plottable)

correl = rcorr(plottable$Promiscuity, plottable$Max_motif, type = "spearman")
rh = round(correl$r[1,2], 2)
p = fancy_scientific_b(correl$P[1,2])
label_spearman = paste0("rho: '", rh, "; p = ", substr(p, 2, nchar(p)))
L = ggplot(plottable, aes(x = Promiscuity, y = Max_motif)) + 
  stat_smooth(method = "loess", 
              formula = y ~ x, 
              geom = "smooth", color = "black", linetype = "dashed", fill = "lightgrey", linewidth = 0.75, span = 1.5) +
  geom_point(colour = "#00215E", size = 2) +
  ylab("Maximum amino acid specificity") +
  xlab("HLA-II allele promiscuity") +
  theme_classic() +
  theme(legend.position = "n") +
  labs(tag = "D", face = "bold") +
  annotate(geom = "text", x = 2.9, y = 3.5, label = label_spearman, parse = TRUE)+
  theme(plot.tag = element_text(face = "bold"))


layout = rbind(c(1,1,1,1,1,2,2,2,2,2,2,3,3,3),
               c(1,1,1,1,1,2,2,2,2,2,2,3,3,3),
               c(1,1,1,1,1,2,2,2,2,2,2,3,3,3),
               c(4,4,4,4,4,4,4,4,4,4,4,4,4,4),
               c(4,4,4,4,4,4,4,4,4,4,4,4,4,4),
               c(5,5,6,6,6,7,7,7,8,8,8,9,9,9),
               c(5,5,6,6,6,7,7,7,8,8,8,9,9,9))#,
               # c(10,10,10,10,10,10,10,10,10,10,10,11,11,11),
               # c(10,10,10,10,10,10,10,10,10,10,10,11,11,11),
               # c(10,10,10,10,10,10,10,10,10,10,10,12,12,12),
               # c(10,10,10,10,10,10,10,10,10,10,10,12,12,12))
               # 

figure_1 =  arrangeGrob(A,B,Bb,C,D,f,G,H,I, layout_matrix = layout)
ggsave(plot = figure_1, filename = "../../results/plots/MM_Figure_1.png", width = 38, height = 32, units = "cm", dpi = 300)

layout = rbind(c(10,10,10,10,10,10,10,10,10,10,10,11,11,11),
c(10,10,10,10,10,10,10,10,10,10,10,11,11,11),
c(10,10,10,10,10,10,10,10,10,10,10,12,12,12),
c(10,10,10,10,10,10,10,10,10,10,10,13,13,13),
c(10,10,10,10,10,10,10,10,10,10,10,13,13,13))


figure_2 =  arrangeGrob(J,K,Kb, L, layout_matrix = layout)
ggsave(plot = figure_2, filename = "../../results/plots/MM_Figure_2_revision_ppm.pdf", width = 38, height = 22, units = "cm", dpi = 300)
