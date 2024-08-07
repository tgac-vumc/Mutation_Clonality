library(ggpubr)

plot_data_LUAD <-
    data.frame(
        Panel = factor(c(rep('MayoComplete',6),rep('InhousePanel',6),rep('TSO500',6),rep('F1CDX',6)), levels = c('MayoComplete','InhousePanel','F1CDX','TSO500')),
        Percentages = c(57,3,26,71,5,38,   #64,4,32, # Mayo
                        84,4,6,90,6,10,    #87,5,8, # Inhouse
                        92,4,4,95,2,3,     #93.5,3,3.5,# TSO500
                        91,4,5,95,2,3),     # 93,3,4), # F1CDX  
        Label =
            factor(rep(rep(c('Correct','Incorrect','Inconclusive'),2),4), levels = c('Incorrect','Inconclusive','Correct')),
        Subtype = 'LUAD')


plot_data_LUSC <-
    data.frame(
        Panel = factor(c(rep('MayoComplete',6),rep('InhousePanel',6),rep('TSO500',6),rep('F1CDX',6)), levels = c('MayoComplete','InhousePanel','F1CDX','TSO500')),
        Percentages = c(78,0,22,41,0,59,   #64,4,32, # Mayo
                        95,1,4,97,0,3,    #87,5,8, # Inhouse
                        95,1,4,97,0,3,     #93.5,3,3.5,# TSO500
                        95,1,4,97,0,3),     # 93,3,4), # F1CDX  
        Label = factor(rep(rep(c('Correct','Incorrect','Inconclusive'),2),4), levels = c('Incorrect','Inconclusive','Correct')),
        Subtype = 'LUSC')

plot_data_LUAD
plot_data <- rbind(plot_data_LUAD,plot_data_LUSC)

pdf('Barchart_VisualAbstract_LUAD.pdf',height = 4,width = 5.5)
ggbarplot(
  plot_data_LUAD, x = "Panel", y = "Percentages", add = 'mean_se',
   fill = "Label", palette = c('#ff4d55ff','#ccccccff','#51be60ff')
) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.title= element_blank()) +
     scale_y_continuous(limits = c(0,120),expand = c(0, 0), breaks = c(25,50,75,100)) 
dev.off()

pdf('Barchart_VisualAbstract_LUSC.pdf',height = 4,width = 5.5)
ggbarplot(
  plot_data_LUSC, x = "Panel", y = "Percentages", add = 'mean_se',
   fill = "Label", palette = c('#ff4d55ff','#ccccccff','#51be60ff')
) +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.title= element_blank()) +
    scale_y_continuous(limits = c(0,120),expand = c(0, 0), breaks = c(25,50,75,100)) 
dev.off()


plot_data <-
    data.frame(
64        panel =  factor(c('MayoComplete','InhousePanel','TSO500','F1CDX'), levels = c('MayoComplete','InhousePanel','F1CDX','TSO500')),
        size = c(34338,15944,1284224,3273215),
        Ngenes = c(12,22,485,259)
    )


pdf('LinePlot_size.pdf', height = 1.25 , width = 5.5)
ggline(plot_data,x='panel',y='Ngenes') + ylab('Number of genes') + scale_y_log10() + annotation_logticks() + theme(axis.title.x = element_blank()) +
    font("ylab", size = 8)
dev.off()

