library(magrittr)
library(tidyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(ggtext)
library(tibble)


#read all the files
flist<-c("misbeliefanalytical.csv", "avgmisbelief0.csv","avgmisbelief0001.csv","avgmisbelief001.csv")

dt_list <- lapply(flist, read.delim, sep=",",header = TRUE)

#prepare data for plotting
data<-data.frame(do.call(rbind,dt_list))
names(data)<-c("tmp","correct","negative","optimistic","mb")
data <- cbind(data,level=data$tmp)
data <-data[, c("tmp","level","correct","negative","optimistic","mb")]

data[1:5,1]<-"*A*"
data[6:10,1]<-"0%"
data[11:15,1]<-"0.01%"
data[16:20,1]<-"0.1%"

datasub<-select(data,"tmp", "level", "mb")
df <- melt(datasub, id.vars=c("level", "tmp"))
df[,4] <- round (df[,4],3)


plot<-ggplot(data=df, aes(x=fct_inorder(factor(tmp)), y=value, fill=fct_inorder(factor(tmp)))) +
  geom_bar(stat="identity", color="black")+
  labs(title="*k*-levels", x="Analytic (*A*) versus simulation (*\u03BC*)", y = "*t-T* Frequency")+
  scale_fill_grey(start=0.8, end=0.2)+
  facet_grid(~level)+
  theme_minimal()+
  guides(fill = guide_legend(byrow = TRUE))+
  theme(legend.position ="none",
        legend.text = element_markdown(size=27),
        legend.title = element_markdown(size=30),
        legend.spacing.y = unit(0.1,"cm"),
        plot.title = element_markdown(size=35, hjust = 0.5),
        axis.title.x=element_markdown(size=35),
        axis.title.y=element_markdown(size=35),
        panel.border = element_blank(),
        axis.line = element_line(colour="black", size=1),
        axis.text.x = element_markdown(size=20),
        axis.text.y=element_markdown(size=27),
        strip.text.x=element_text(size=27))
plot


