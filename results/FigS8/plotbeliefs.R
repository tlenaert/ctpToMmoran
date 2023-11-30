library(magrittr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggtext)
library(tidyverse)

flist<-c("beliefanalytical.csv", "avgbelief00001.csv","avgbelief0001.csv","avgbelief001.csv")

dt_list <- lapply(flist, read.delim, header = TRUE)

#prepare data for plotting
data<-data.frame(t(do.call(cbind,dt_list)))
names(data)<-c("mutation","eps","1","2","3","4","5")
data <-subset(data,select = -c(eps))

data[1,1]<-"analytic"
data[2,1]<-"*\u03BC=10^-5*"
data[3,1]<-"*\u03BC=10^-4*"
data[4,1]<-"*\u03BC=10^-3*"


df <- melt(data, id.vars="mutation")
df[,3] <- round (df[,3],3)

plot<-ggplot(data=df, aes(x=fct_inorder(variable), y=value, fill=fct_inorder(factor(mutation)))) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  labs(x="belief(*t*)", y = "Frequency", fill="mutation")+
  #geom_text(aes(label=value), vjust=-0.8, color="black", position = position_dodge(0.9), size=6)+
  scale_fill_brewer( palette="Blues", direction=-1)+
  theme_minimal()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position = c(0.2, 0.8),
        legend.text = element_markdown(size=27),
        legend.title = element_blank(),
        legend.spacing.y = unit(0.5,"cm"),
        axis.title.x = element_markdown(size=35),
        axis.title.y=element_text(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_text(size=27),
        axis.text.y=element_text(size=27))
plot

