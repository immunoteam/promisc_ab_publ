trace(ggpubr:::.cor_test, edit=TRUE) #edit cor.coef.name to rho

library(ggplot2)
library(ggpubr)
library(ggimage)

subjects = read.delim("SDY314-DR54_Tab/SDY314-DR54_Tab/Tab/subject.txt") #Downloaded from ImmPort
gender = data.frame(table(subjects$GENDER))
ggplot(gender, aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0) +
  geom_text(aes(label = paste(Var1, "\n", paste0("n = ", Freq))), 
            position = position_stack(vjust = 0.5)) + 
  scale_fill_manual("legend", values = c("Female" = "#C255BF", "Male" = "#7287D7")) +
  labs(title = "Sex distribution") +
  theme_void() + 
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 


arm_subjects = read.delim("SDY314-DR54_Tab/SDY314-DR54_Tab/Tab/arm_2_subject.txt") #Downloaded from ImmPort
age = data.frame(table(arm_subjects$ARM_ACCESSION))
ggplot(age[1:3,], aes(x = "", y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = c(paste(c("18-30", "60-79", "80-90"), "\n", paste0("n = ", Freq)))), 
            position = position_stack(vjust = 0.5)) + 
  scale_fill_manual("legend", 
                    values = c("ARM2134" = "#27B7EE", 
                               "ARM2135" = "#FFC55A", 
                               "ARM2136" = "#FC4100")) +
  labs(title = "Age distribution") +
  theme_void() + 
  theme(legend.position="none", plot.title = element_text(hjust = 0.5)) 


promiscuity_df = readRDS("promiscuity.RDS")
promiscuity = promiscuity_df$Allele.promiscuity
names(promiscuity) = promiscuity_df$Allele.name
hla_data = read.delim("SDY314-DR54_Tab/SDY314-DR54_Tab/Tab/hla_typing_result.txt") #Downloaded from ImmPort
subjects = unique(hla_data$SUBJECT_ACCESSION)
df = data.frame(subject_id = subjects)
df$drb1_promiscuity = sapply(df$subject_id, FUN = function(x) {
  hla = c(hla_data$ALLELE_1[hla_data$SUBJECT_ACCESSION == x & hla_data$LOCUS_NAME == "HLA-DRB1"], 
          hla_data$ALLELE_2[hla_data$SUBJECT_ACCESSION == x & hla_data$LOCUS_NAME == "HLA-DRB1"])
  pr = promiscuity[substr(gsub("\\*|:", "", hla), 1, 8)]
  mean(pr)
})

arms = c("ARM2134", "ARM2135", "ARM2136")
years = c("18-30 year olds", "60-79 year olds", "80-90 year olds")
plots = list(NULL)
n = 1
for(i in 1:3){
  hai = read.delim("SDY314-DR54_Tab/SDY314-DR54_Tab/Tab/hai_result.txt") #Downloaded from ImmPort
  hai = hai[hai$ARM_ACCESSION == arms[i],]
  df$hai_0 = sapply(df$subject_id, FUN = function(x) {
    temp = hai[hai$SUBJECT_ACCESSION == x,]
    mean(temp$VALUE_REPORTED[temp$STUDY_TIME_COLLECTED == 0])
  })
  df$hai_28 = sapply(df$subject_id, FUN = function(x) {
    temp = hai[hai$SUBJECT_ACCESSION == x,]
    mean(temp$VALUE_REPORTED[temp$STUDY_TIME_COLLECTED == 28])
  })
  plots[[n]] = ggplot(df, aes(drb1_promiscuity, hai_28/hai_0)) + 
    geom_smooth(method = lm, color = "black", fill = "lightgrey") + 
    geom_point(color = "#00215E") + 
    labs(title = years[i], x = "HLA-DRB1 allele promiscuity", 
         y = "Fold change in HAI-titer \n (day 28/day 0)") +
    stat_cor(method = "spearman", 
             aes(label = paste("rho==", ..r.., "*paste(\";\")~~~italic(p) ==",
                               ifelse(..p..<0.01, paste(10^(log10(..p..) %% 1), "%*% 10^", 
                                                        floor(log10(..p..))), ..p..)))) +
    scale_y_continuous(trans = 'log2') +
    theme_classic()  
  n = n + 1
}

for(i in 1:3){
  neut_ab = read.delim("SDY314-DR54_Tab/SDY314-DR54_Tab/Tab/neut_ab_titer_result.txt") #Downloaded from ImmPort
  neut_ab = neut_ab[neut_ab$ARM_ACCESSION == arms[i],]
  df$neut_ab_0 = sapply(df$subject_id, FUN = function(x) {
    temp = neut_ab[neut_ab$SUBJECT_ACCESSION == x ,]
    mean(temp$VALUE_REPORTED[temp$STUDY_TIME_COLLECTED == 0])
  })
  df$neut_ab_28 = sapply(df$subject_id, FUN = function(x) {
    temp = neut_ab[neut_ab$SUBJECT_ACCESSION == x ,]
    mean(temp$VALUE_REPORTED[temp$STUDY_TIME_COLLECTED == 28])
  })
  plots[[n]] = ggplot(df, aes(drb1_promiscuity, neut_ab_28/neut_ab_0)) + 
    geom_smooth(method = lm, color = "black", fill = "lightgrey") + 
    geom_point(color = "#00215E") + 
    labs(x = "HLA-DRB1 allele promiscuity", 
         y = "Fold change in neutr. ab. titer \n (day 28/day 0)") +
    stat_cor(method = "spearman", 
             aes(label = paste("rho==", ..r.., "*paste(\";\")~~~italic(p) ==",
                               ifelse(..p..<0.01, paste(10^(log10(..p..) %% 1), "%*% 10^", 
                                                        floor(log10(..p..))), ..p..)))) +
    scale_y_continuous(trans = 'log2') +
    theme_classic()  
  n = n + 1
}


panel2 <- ggarrange(plots[[1]], plots[[2]], plots[[3]], 
                    plots[[4]], plots[[5]], plots[[6]], 
                    labels = c("B", "", "", "C", "", "", ""))


