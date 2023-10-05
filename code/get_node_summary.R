library(ggplot2)
library(reshape2)
library(ggalluvial)
library(ggExtra)

##功能统计
##获得边的统计结果,以桑葚图表示
node_function <- read.csv("../data/RMT-ASV_0.00001_0.01_node-all.csv", header=T, row.names=1)  ##4992个ASV
nrow(node_function[node_function$v_phylum %in% c("Proteobacteria","Bacteroidota","Myxococcota","Chloroflexi"),])/nrow(node_function)  ##0.6752804
##以phylum和function分类
top_phylum_list <- c("Proteobacteria","Bacteroidota","Myxococcota","Chloroflexi","Planctomycetota","Acidobacteriota",
                     "Patescibacteria","Verrucomicrobiota","Actinobacteriota","Bdellovibrionota","Desulfobacterota")
phylum_function <- node_function[,c(5,24:39)]
phylum_function[!phylum_function$v_phylum %in% top_phylum_list,]$v_phylum <- "Others"
phylum_function_long <- melt(phylum_function,id.vars = c('v_phylum'), variable.name='Function', value.name='weight')
phylum_function_long <- phylum_function_long[phylum_function_long$weight !=0,]
phylum_function_summary <- as.data.frame(table(phylum_function_long[c('v_phylum','Function')]))
write.csv(phylum_function_summary, './igraph/pearson/RMT-ASV_0.00001_0.01_phylum_to_function-summary.csv',quote = FALSE,row.names = FALSE)
sum(phylum_function_summary[phylum_function_summary$Function=="Nitrifiers",]$Freq) ##55
sum(phylum_function_summary[phylum_function_summary$Function=="Denitrifiers",]$Freq) ##430
sum(phylum_function_summary[phylum_function_summary$Function=="Anammox",]$Freq) ##1
sum(phylum_function_summary[phylum_function_summary$Function=="DNRA",]$Freq) ##85
sum(phylum_function_summary[phylum_function_summary$Function=="PAOs",]$Freq) ##46
sum(phylum_function_summary[phylum_function_summary$Function=="GAOs",]$Freq) ##60
sum(phylum_function_summary[phylum_function_summary$Function=="SRB",]$Freq) ##61
sum(phylum_function_summary[phylum_function_summary$Function=="SOB",]$Freq) ##34
sum(phylum_function_summary[phylum_function_summary$Function=="Filamentous",]$Freq) ##246
sum(phylum_function_summary[phylum_function_summary$Function=="Acetogen",]$Freq) ##7
sum(phylum_function_summary[phylum_function_summary$Function=="Methanogen",]$Freq) ##3
sum(phylum_function_summary[phylum_function_summary$Function=="Fermenter",]$Freq) ##345
sum(phylum_function_summary[phylum_function_summary$Function=="Sugars",]$Freq) ##964
sum(phylum_function_summary[phylum_function_summary$Function=="Lipids",]$Freq) ##584
sum(phylum_function_summary[phylum_function_summary$Function=="Proteins",]$Freq) ##835
sum(phylum_function_summary[phylum_function_summary$Function=="Unassigned",]$Freq) ##3648

phylum_function_form <- to_lodes_form(data.frame(phylum_function_long),
                              key = "Phylum_func",
                              axes = 1:2)
phylum_function_form$Phylum_func <- factor(phylum_function_form$Phylum_func, levels =c("v_phylum","Function"))
phylum_function_form$stratum <- factor(phylum_function_form$stratum, levels = c("Proteobacteria","Bacteroidota","Myxococcota","Chloroflexi","Planctomycetota","Acidobacteriota","Patescibacteria","Verrucomicrobiota","Actinobacteriota","Bdellovibrionota","Desulfobacterota",
                                                                                "Others","Nitrifiers","Denitrifiers","Anammox","DNRA","PAOs","GAOs","SRB","SOB","Filamentous","Acetogen","Methanogen","Fermenter","Sugars","Lipids","Proteins","Unassigned"))
mycol <- c("#006C71", "#FF8003","#FF82FF","#3CDB00","#814A18","#00DCFF","#C8AF00","#00D19A","#72466B","#4F8D32","#FF5E71","#C0C0C0",'#008B8B','#556B2F','#228B22','#20B2AA','#483D8B','#6495ED','#DDA0DD','#D87093','#FF1493','#A0522D','#BC8F8F','#0066ff','#FFD700','#006699','#666699','#999999')
ggplot(phylum_function_form,
       aes(x = Phylum_func, y = weight, stratum = stratum, alluvium = alluvium,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  geom_flow(width = 1/6) + #线跟方块间空隙的宽窄
  geom_stratum(alpha = 0.5,width = 1/6) + #方块的透明度、宽度
  geom_text(stat = "stratum", size = 2,color="black") + #文字大小、颜色
  #不喜欢默认的配色方案，用前面自己写的配色方案
  scale_fill_manual(values = mycol) +
  xlab("") + ylab("") +
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(panel.border = element_blank()) + #去除外层边框
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank()) + #去掉坐标轴
  ggtitle("")+
  guides(colour = "colorbar") 
dev.off()
ggsave("./igraph/pearson/RMT-ASV_0.00001_0.01_phylum_to_function-sankey.pdf",width = 9, height = 6)



