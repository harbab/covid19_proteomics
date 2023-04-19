# Load the libraries and set the working directory where you saved the folder
library(ggplot2)
library(scales)
library(ggsci)
setwd("code_resubmission/data/")

s1 = as.data.frame(read_csv("table_s1.csv"))
det_prot = read.csv("detected_proteins.csv", row.names = 1)

# Plot the number of identified proteins in the dataset obtained from public resources or the authors
setwd("../results/")
fit = lm(s1$Proteins_identified ~ s1$N_participants)
sfit = summary(fit)
p = ggplot(s1, aes(N_participants, Proteins_identified)) +
  geom_point() + theme_classic(base_size=12) +
  stat_smooth(method = "loess", col = "red", se=F) +
  geom_text_repel(label=s1$Author) +
  scale_x_continuous("N participants", breaks=seq(0,300,30)) +
  scale_y_continuous("Proteins identified", breaks=seq(0,2200,200))
p

ggsave("participants_vs_proteins.pdf", p, width=6.25, height=5)

# Plot in how many cohorts a protein has been identified
t = table(det_prot$nstud)
tf = as.data.frame(cbind(names(t), t))
tf$t = as.numeric(tf$t)
tf$V1 = factor(tf$V1, levels=1:nrow(det_prot))
tf$perc = round(tf$t/sum(tf$t)*100, 1)

bp <- ggplot(tf, aes(V1, perc))+
  geom_bar(stat="identity", fill="black") + theme_classic() +
  scale_y_continuous("% of identified proteins", breaks=seq(0,100,5)) +
  xlab("N of studies") +
  geom_text(aes(V1, perc+5), label=paste(tf$t, " (", tf$perc, "%)", sep=""), 
            angle=45, hjust=0.2)
bp

ggsave("pie_studies.pdf", bp, dpi=300, width=6, height=4.8)

# Plot the number of identified proteins and number of proteins included in the meta-analysis
table(s1$Year)
s1[,c(7:11)] = apply(s1[,c(7:11)], 2, as.numeric)
sdf = as.data.frame(cbind(c(s1$Proteins_identified, s1$`Proteins_in_Meta-analysis`),
                          c(rep("Identified", nrow(s1)), rep("Included in MA", nrow(s1))),
                          c(rep(s1$Author, 2))))
colnames(sdf) = c("prot", "cat", "author")
sdf$prot = as.numeric(sdf$prot)
sdf = sdf[order(sdf$prot, decreasing = F),]
sdf$author = as.factor(sdf$author)
sdf$author = ordered(sdf$author, levels=s1$Author[order(s1$`Proteins_in_Meta-analysis`)])

p = ggplot(sdf, aes(prot, author, fill=cat)) + geom_bar(stat="identity", width=0.75, position="dodge") +
  theme_classic() +
  scale_x_continuous("N of protein", breaks=seq(0,2000, 200)) +
  theme(legend.position = c(0.8, 0.2)) +
  scale_fill_manual("Proteins", values=c("#009E73", "#CC79A7"))
p

ggsave("studies_proteins.pdf", p, width=6.25, height=5)

aut = list()
for (i in colnames(det_prot)[2:12]) {
  aut[[i]] = det_prot$gene_name[which(det_prot[,i]==1)]
}

det_prot$included = 0
det_prot$included[which(det_prot$nstud>1)] = 1
det_prot = det_prot[order(det_prot$included, det_prot$nstud, det_prot$Babacic, decreasing = T),]

library(RColorBrewer)
display.brewer.all()
cols = brewer.pal(n = 12, name = "Paired")
cols = cols[-11]
cols = c(cols, "#D55E00")

ha <- rowAnnotation(
  "Included in MA" = det_prot$included,
  "Babacic" = det_prot$Babacic,
  "Shen" = det_prot$Shen,
  "Di" = det_prot$Di,
  "Suvarna" = det_prot$Suvarna,
  "Overmyer" = det_prot$Overmyer,
  "Shu" = det_prot$Shu,
  "Geyer" = det_prot$Geyer,
  "Sullivan" = det_prot$Sullivan,
  "Zhang" = det_prot$Zhang,
  "Messner - Discovery" = det_prot$Messner.Discovery,
  "Messner - Validation" = det_prot$Messner.Validation,
  col=list("Included in MA"=c("1"=cols[1],"0"="white"), "Babacic"=c("1"=cols[2],"0"="white"),
           "Shen"=c("1"=cols[3],"0"="white"), "Di"=c("1"=cols[4],"0"="white"), 
           "Suvarna"=c("1"=cols[5],"0"="white"), "Overmyer"=c("1"=cols[6],"0"="white"),
           "Shu"=c("1"=cols[7],"0"="white"), "Geyer"=c("1"=cols[8],"0"="white"), 
           "Sullivan"=c("1"=cols[9],"0"="white"), "Zhang"=c("1"=cols[10],"0"="white"),
           "Messner - Discovery"=c("1"=cols[11],"0"="white"), "Messner - Validation"=c("1"=cols[12],"0"="white"))
)

dev.off()
draw(ha)

pdf("protein_overlap_studies.pdf", width=6.25, height=3)
draw(ha)
dev.off()

```