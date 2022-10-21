# Set the working directory where you saved the folder
setwd("code_submission/data/study_summary")
s1 = as.data.frame(read_csv("table_s1.csv"))
det_prot = read.csv("detected_proteins.csv", row.names = 1)

setwd("../../results/")
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

t = table(det_prot$nstud)
tf = as.data.frame(cbind(names(t), t))
tf$t = as.numeric(tf$t)
tf$V1 = factor(tf$V1, levels=1:11)

library(scales)
library(ggsci)

bp <- ggplot(tf, aes(x=1, y=t, fill=V1))+
  geom_bar(width = 1, stat = "identity")
pie = bp + coord_polar("y", start=0) +
  scale_fill_brewer("N of studies", palette="Paired") +
  theme_void() +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  geom_text_repel(aes(x=1+seq(0,1,0.1), y = t, 
                      label = percent(tf$t/sum(tf$t))), size=5)
pie

ggsave("pie_studies.pdf", pie, dpi=300, width=8, height=6.4)

table(s1$Year)

s1[,c(7:11)] = apply(s1[,c(7:11)], 2, as.numeric)
s1 = s1[-c(1,2,6,14),]
sdf = as.data.frame(cbind(c(s1$Proteins_identified, s1$`Proteins_in_Meta-analysis`),
                          c(rep("Identified", 11), rep("Included in MA", 11)),
                          c(rep(s1$Author, 2))))
colnames(sdf) = c("prot", "cat", "author")
sdf$prot = as.numeric(sdf$prot)
sdf = sdf[order(sdf$prot, decreasing = F),]
sdf$author = as.factor(sdf$author)
sdf$author = ordered(sdf$author, levels=s1$Author[order(s1$`Proteins_in_Meta-analysis`)])

p = ggplot(sdf, aes(prot, author, fill=cat)) + geom_bar(stat="identity", width=0.75, position="dodge") +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,2000, 200)) +
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