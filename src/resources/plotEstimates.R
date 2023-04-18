library(stringr)
library(reshape2)
library(dplyr)
library(ggplot2)
library(cowplot)

t_mapFastq=read.delim("grandslam_t15_noresc/noresc.binomEstimated.tsv")
t_mapFastq$Condition=factor(as.character(t_mapFastq$Condition),levels=as.character(t_mapFastq$Condition))
t_mapFastq$Status <- paste(t_mapFastq$Status, "unrescued", sep="_")

t_resc=read.delim("grandslam_t15_resc/resc.binomEstimated.tsv")
t_resc$Condition=factor(as.character(t_resc$Condition),levels=as.character(t_resc$Condition))
t_resc$Status <- paste(t_resc$Status, "rescued", sep="_")
t_resc$Condition <- str_replace(t_resc$Condition, "_rescued", "")

t_merge <- rbind(t_mapFastq,t_resc)
t_merge=t_merge[!grepl("no4sU",t_merge$Condition),]
t_merge$conv_old <- NULL
r_merge <- melt(t_merge[,c(1:2,6)])

ggplot(r_merge,aes(Condition,value,col=variable,shape=Status))+
  geom_point(size=2, color = "#D95F02")+scale_shape_manual(values=c(16,1,17,2), labels = c("rescued","no rescue"))+scale_color_brewer(labels = c("new RNA"))+
  theme_bw()+theme(text=element_text(size=24),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ylab("T-to-C Incorporation frequency\ninto new RNA")+xlab(NULL) +
  theme_cowplot()+theme(legend.title = element_blank()) +
  xlab("T-to-C conversion rate") + ylim(c(0,max(r_merge$value)))
ggsave("plotEstimates.png", bg="white")

