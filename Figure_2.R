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
library(pheatmap)
library(ggsci)
library(ggplot2)
library(universalmotif)

mypal = pal_npg("nrc", alpha = 1)(9)
mypal[1] = "#FC4100"
mypal[4] = "#FFC55A"
mypal[7] = "#2C4E80"

#Panel A

positions_rowordered = readRDS("heatmap_mtx.RDS")
breaks = seq(min(positions_rowordered), max(positions_rowordered), length.out = 100)
A = pheatmap(t(positions_rowordered[,apply(positions_rowordered, 2, max)>0.4]), cluster_cols = F, 
         cluster_rows = F, color = colorRampPalette(c("#2C4E80", "red", "yellow"), bias = 2)(length(breaks) - 1), 
         breaks = breaks)$gtable

#Panel B

mot1 = readRDS("motif_DRB11602.RDS")
klvals_drb11602 = as.vector(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif)
drb11602 = view_motifs(mot1, use.type = "ICM", min.height = 0, y.spacer = 0.001) + ylim(c(0,2.75))
drb11602 = drb11602 + 
  annotate(geom = "text", x = 5, y = 2.75, label = "HLA-DRB1*16:02", size = 4) +
  theme(axis.title=element_text(size=14)) + 
  labs(tag = "B", face = "bold") +
  theme(plot.tag = element_text(face = "bold", size = 12))


mot1 = readRDS("motif_DRB10804.RDS")
klvals_drb10804 = as.vector(view_motifs(mot1, use.type = "PWM", min.height = 0, y.spacer = 0.001, return.raw = T)$motif)
drb10804 = view_motifs(mot1, use.type = "ICM", min.height = 0, y.spacer = 0.001) + ylim(c(0,2.75))
drb10804 = drb10804 + 
  annotate(geom = "text", x = 5, y = 2.75, label = "HLA-DRB1*08:04", size = 4) +
  theme(axis.title=element_text(size=14)) + 
  labs(tag = " ", face = "bold") +
  theme(plot.tag = element_text(face = "bold", size = 12))

df_dens = data.frame(Allele = c(rep("DRB1*16:02", 180), rep("DRB1*08:02", 180)), Values = c(klvals_drb11602, klvals_drb10804))

B = grid.arrange(drb11602, drb10804, nrow = 2)

#Panel C

C = ggplot(df_dens, aes(x = Values, colour = Allele)) + geom_density(bw = 0.4) + theme_bw() + 
  scale_color_manual(values = substr(mypal[c(7,1)], 1, 7)) + xlab("log2(observed/background)") + ylab("Density") +
  theme(legend.position="none") + annotate(geom = "text", x = -6, y = 0.15, label = "HLA-DRB1*16:02", size = 3, color = "#FC4100") +
  annotate(geom = "text", x = -4, y = 0.3, label = "HLA-DRB1*08:04", size = 3, color = "#2C4E80") + 
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(tag = "C", face = "bold") +
  theme(plot.tag = element_text(face = "bold", size = 12))

#Panel D

load("HLA_motif_max.RData")
promiscuity = readRDS("promiscuity.RDS")
plottable = data.frame(Max_motif = apply(positions_rowordered, 1, max), 
                       Promiscuity = promiscuity$Allele.promiscuity[match(rownames(positions_rowordered), rownames(promiscuity))])

plottable = plottable[grepl("DRB1", rownames(plottable)) | grepl("DP", rownames(plottable)) |
                                        grepl("DQ", rownames(plottable)),]
plottable = as.data.frame(plottable)

correl = rcorr(plottable$Promiscuity, plottable$Max_motif, type = "spearman")
rh = round(correl$r[1,2], 2)
p = fancy_scientific_b(correl$P[1,2])
label_spearman = paste0("rho: '", rh, "; p = ", substr(p, 2, nchar(p)))
D = ggplot(plottable, aes(x = Promiscuity, y = Max_motif)) + 
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

layout = rbind(c(10,10,10,10,10,10,10,10,10,10,10,11,11,11),
c(10,10,10,10,10,10,10,10,10,10,10,11,11,11),
c(10,10,10,10,10,10,10,10,10,10,10,12,12,12),
c(10,10,10,10,10,10,10,10,10,10,10,13,13,13),
c(10,10,10,10,10,10,10,10,10,10,10,13,13,13))


figure_2 =  arrangeGrob(A,B,C,D, layout_matrix = layout)
ggsave(plot = figure_2, filename = "Figure_2.png", width = 38, height = 22, units = "cm", dpi = 300)
