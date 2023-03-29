
##########################
### Cells Characteristics
##########################

library(ggplot2)
library(tcltk2)

select.cellsChar <- tk_choose.files(caption = "Select Cells characteristics file")
cells <- read.csv(select.cellsChar,header=T,sep="\t")
dir <- dirname(select.cellsChar)


pdf(file = paste(dir,'\\CellsIsolationCharacteristics.pdf',sep =''),
    width = 8,
    height = 6)

ggplot(cells, aes(x = factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted")),
                  y = Diameter)) +
  geom_violin(position = position_dodge(),aes(color=factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted")),
              fill=factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted"))),alpha=0.3,size=1)+
  geom_point(position = position_jitter(),pch = 21,aes(fill=factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted"))),
  alpha=0.5,width=0.07,size=2, color='darkgray')+
  geom_boxplot(position = position_dodge(),color="black",fill="white",alpha=0.8,size=0.5,width=0.07,outlier.shape = NA)+
  theme_bw() + xlab("")+ ylab("Diameter") + labs(title = "Diameter by Cell Type")+
  scale_color_manual(values = c("#FFE082","#EF9A9A", "#C5E1A5", "#90CAF9", "#AF7AC5","#BDBDBD" ))+
  scale_fill_manual(values = c("#FFE082","#EF9A9A", "#C5E1A5", "#90CAF9", "#AF7AC5","#BDBDBD" ))+
  theme(legend.position = "none",legend.title = element_text(size=10),
        legend.text = element_text(size=8), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=13), axis.ticks = element_blank()) + 
  guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE))


ggplot(cells, aes(x = factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted")),
                  y = Elongation)) +
  geom_violin(position = position_dodge(),aes(color=factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted")),
                                              fill=factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted"))),alpha=0.3,size=1)+
  geom_point(position = position_jitter(),pch = 21,aes(fill=factor(CellType, level=c('mESCs', 'Endoderm', 'Mesoderm', 'Ectoderm', "Unsorted"))),
             alpha=0.5,width=0.07,size=2, color='darkgray')+
  geom_boxplot(position = position_dodge(),color="black",fill="white",alpha=0.8,size=0.5,width=0.07,outlier.shape = NA)+
  theme_bw() + xlab("")+ ylab("Elongation") + labs(title = "Elongation by Cell Type")+
  scale_color_manual(values = c("#FFE082","#EF9A9A", "#C5E1A5", "#90CAF9", "#AF7AC5","#BDBDBD" ))+
  scale_fill_manual(values = c("#FFE082","#EF9A9A", "#C5E1A5", "#90CAF9", "#AF7AC5","#BDBDBD" ))+
  theme(legend.position = "none",legend.title = element_text(size=10),
        legend.text = element_text(size=8), 
        axis.text = element_text(size=12), 
        axis.title = element_text(size=13), axis.ticks = element_blank()) + 
  guides(fill = "none", color = guide_legend(nrow=3,byrow=TRUE))


dev.off()
