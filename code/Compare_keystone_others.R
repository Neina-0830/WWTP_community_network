library(ggplot2)
library(igraph)
library(dplyr)
library(reshape2)

Node_all <- read.csv('./Gephi/pearson/RMT-ASV_0.00001_0.01_node_all.csv',header=T,row.names = 1)
Node_select <- Node_all[,c("v_degree","v_betweenness","v_frequency","MRA")]
Node_select$Group <- "Others"
Node_select[Node_select$v_degree>100&Node_select$v_betweenness<5000,]$Group <- "Keystone"

par(mar=c(3,5,2,2))
p<-ggplot(Node_select, aes(x=v_degree, y=v_betweenness, color=Group)) +
  geom_point(size=1,alpha=0.8) +
  scale_color_manual(values = c("darkred","gray70"),limits = c('Keystone','Others')) +
  ylim(0,250000) +
  xlim(0,300) +
  geom_vline(aes(xintercept=100), colour = "black", linetype=2, size =0.3) +
  geom_hline(aes(yintercept=5000), colour = "black", linetype=2, size =0.3)
p + labs(x="Node degree", y = "Node betweenness") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
dev.off()
ggsave("./Gephi/pearson/Keystone_degree_betweenness_plot.pdf", width = 8, height = 6)

##keystone的丰度和频率分布情况
par(mar=c(3,5,2,2))
p<-ggplot(Node_select, aes(x=v_frequency, y=MRA, color=Group)) +
  geom_point(size=0.5,alpha=0.6) +
  scale_color_manual(values = c("darkgreen","gray70"),limits = c('Keystone','Others')) +
  ylim(0,0.015) +
  xlim(0,1) +
  geom_vline(aes(xintercept=0.1), colour = "red", linetype=2, size =0.3)
p + labs(x="Occurrence frequency", y = "Mean relative abundance") + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 14, color = "black"),
        axis.title = element_text(size = 14, color = "black"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
dev.off()
ggsave("./Gephi/pearson/Keystone_frequency_MRA_plot.pdf", width = 8, height = 6)

Keystone_select <- Node_select[Node_select$Group=="Keystone",]
library(ggalluvial)
library(ggExtra)
pdf(file="./Gephi/pearson/Node_keystone_MRA_frequency_plot_bar.pdf", width = 4, height = 3.4)
p <- ggplot(Keystone_select, aes(x=v_frequency, y=MRA)) +
  geom_point(size=1) +
  labs(x="Occurrence frequency", y = "Average relative abundance") + 
  theme_bw()+
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 12, color = "black"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm')) +
  theme(legend.position="none")
# Set relative size of marginal plots (main plot 10x bigger than marginals)
p1 <- ggMarginal(p, type="histogram", xparams = list(bins=20))
print(p1)
dev.off()


##比较不同keystone和others类群的信息（物种组成、功能组成、出现频率及相对丰度）
Node_all <- read.csv('./Gephi/pearson/RMT-ASV_0.00001_0.01_node_all_MRA.csv',header=T,row.names = 1)
keystone <- read.csv("./keystone/pearson/RMT-ASV_0.00001_0.01_node_keystone.csv",header=T,row.names = 1)
Node_all$Key_type <- "Others"
Node_all[rownames(Node_all)%in%keystone$Label,]$Key_type <- "Keystone"
write.csv(Node_all,"./keystone/pearson/RMT-ASV_0.00001_0.01_node_all_keystone_others.csv", quote = FALSE)
##phylum
all_phylum <- as.data.frame(table(Node_all$Key_type, Node_all$v_phylum))
colnames(all_phylum) <- c("Key_type", "Phylum", "ASV_num")
all_phylum_wide <- dcast(all_phylum, Key_type ~ Phylum, value.var="ASV_num")
write.csv(all_phylum_wide, "./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_Phylum_ASVs_num.csv",quote = FALSE, row.names = FALSE)
all_phylum_wide <- read.csv("./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_Phylum_ASVs_num.csv",header=T,row.names=1)
all_phylum_prop <- all_phylum_wide/rowSums(all_phylum_wide)
write.csv(all_phylum_prop, "./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_Phylum_proportion.csv",quote = FALSE)
all_phylum_prop['Key_type'] <- rownames(all_phylum_prop)
all_phylum_prop_long <- melt(all_phylum_prop, id = 'Key_type')
all_phylum_prop_long$Key_type <- factor(all_phylum_prop_long$Key_type,levels = c("Keystone","Others"))
##部分门涂色
all_phylum_prop_long$variable<-factor(all_phylum_prop_long$variable,levels = c("Caldatribacteriota","Caldisericota","Cloacimonadota","Coprothermobacterota","Dadabacteria","Deinococcota","Halobacterota","Margulisbacteria","MBNT15","midas_p_12349","midas_p_49974","Euryarchaeota","FCPU426","midas_p_31827",
                                                                               "midas_p_7082","Synergistota","WPS-2","Sumerlaeota","midas_p_25076","Fusobacteriota","Hydrogenedentes","Calditrichota","Campylobacterota","Nitrospirota","Fibrobacterota","Dependentiae","SAR324_cladeMarine_group_B","Elusimicrobiota","Spirochaetota","Armatimonadota","Cyanobacteria","Latescibacterota",
                                                                               "Gemmatimonadota","Firmicutes","Unknown","Desulfobacterota","Bdellovibrionota","Actinobacteriota","Verrucomicrobiota","Patescibacteria","Acidobacteriota","Planctomycetota","Chloroflexi","Myxococcota","Bacteroidota","Proteobacteria"))
par(mar=c(3,5,2,2))
ggplot(all_phylum_prop_long, aes(Key_type, value*100, fill= variable)) +
  geom_bar(stat = "identity")+
  labs(y = 'Proportion of Phylum(%)') +
  scale_fill_manual(values = c("#C0C0C0","#C0C0C0", "#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0", "#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0", "#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0"
                               ,"#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#C0C0C0","#FF5E71","#4F8D32","#72466B","#00D19A","#C8AF00","#00DCFF","#814A18","#3CDB00","#FF82FF","#FF8003","#006C71"),
                    limits = c("Caldatribacteriota","Caldisericota","Cloacimonadota","Coprothermobacterota","Dadabacteria","Deinococcota","Halobacterota","Margulisbacteria","MBNT15","midas_p_12349","midas_p_49974","Euryarchaeota","FCPU426","midas_p_31827",
                               "midas_p_7082","Synergistota","WPS-2","Sumerlaeota","midas_p_25076","Fusobacteriota","Hydrogenedentes","Calditrichota","Campylobacterota","Nitrospirota","Fibrobacterota","Dependentiae","SAR324_cladeMarine_group_B","Elusimicrobiota","Spirochaetota","Armatimonadota","Cyanobacteria","Latescibacterota",
                               "Gemmatimonadota","Firmicutes","Unknown","Desulfobacterota","Bdellovibrionota","Actinobacteriota","Verrucomicrobiota","Patescibacteria","Acidobacteriota","Planctomycetota","Chloroflexi","Myxococcota","Bacteroidota","Proteobacteria"))+
  ylim(0,100) +
  theme(axis.text = element_text(size = 10,face = "plain",color ='black'), axis.title = element_text(size = 10,color ='black')) +
  theme(legend.text = element_text(size = 8,color ='black'))+
  geom_col(position = 'stack', width = 0.6)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
dev.off()
ggsave("./keystone/pearson/Key_type_Phylum_proportion-ASVs.pdf",height=6,width=11)

##其他分类水平组成
Node_all <- read.csv("./keystone/pearson/RMT-ASV_0.00001_0.01_node_all_keystone_others.csv",header = T, row.names = 1)
##family
all_family <- as.data.frame(table(Node_all$Key_type, Node_all$v_family))
colnames(all_family) <- c("Key_type", "Phylum", "ASV_num")
all_family_wide <- dcast(all_family, Key_type ~ Phylum, value.var="ASV_num")
write.csv(all_family_wide, "./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_family_ASVs_num.csv",quote = FALSE, row.names = FALSE)
all_family_wide <- read.csv("./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_family_ASVs_num.csv",header=T,row.names=1)
all_family_prop <- all_family_wide/rowSums(all_family_wide)
write.csv(all_family_prop, "./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_family_proportion.csv",quote = FALSE)
all_family_prop <- read.csv('./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_Family_proportion-new.csv',header=T,row.names = 1) ##**-new数据是将前面数据表中的midas..属信息都替换成Others之后的结果
all_family_prop['Key_type'] <- rownames(all_family_prop)
all_family_prop_long <- melt(all_family_prop, id = 'Key_type')
all_family_prop_long$Key_type <- factor(all_family_prop_long$Key_type,levels = c("Keystone","Others"))



##其他分布特征
Node_all$Key_type <- factor(Node_all$Key_type, levels = c("Keystone","Others"))
##frequency
par(mar=c(3,5,2,2))
p<- ggplot(Node_all, aes(x=Key_type, y=v_frequency, fill=Key_type)) + 
  geom_boxplot(show.legend=FALSE) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_fill_manual(values = c("#F0E442", "#0072B2"),
                     limits = c("Keystone","Others")) +
  ylim(0,1.0)
p + labs(y = "Occurrence frequency") +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Key_type_frequency_compare-boxplot.pdf", width = 4, height = 3)
t.test(Node_all[Node_all$Key_type=="Keystone",]$v_frequency, Node_all[Node_all$Key_type=="Others",]$v_frequency) ##p-value < 2.2e-16(0.02737933 0.07387103)

##MRA
par(mar=c(3,5,2,2))
p<- ggplot(Node_all, aes(x=Key_type, y=MRA, fill=Key_type)) + 
  geom_boxplot(show.legend=FALSE) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_fill_manual(values = c("#F0E442", "#0072B2"),
                    limits = c("Keystone","Others")) +
  ylim(0,0.015)
p + labs(y = "Average relative abundance") +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Key_type_MRA_compare-boxplot.pdf", width = 4, height = 3)
t.test(Node_all[Node_all$Key_type=="Keystone",]$MRA, Node_all[Node_all$Key_type=="Others",]$MRA) ##p-value = 3.98e-16(4.809765e-05 1.456905e-04)

##功能比较
Node_func_all <- read.csv("./Function/pearson/RMT-ASV_0.00001_0.01_node-all-function.csv",header=T,row.names = 1)
Node_func_all$Key_type <- "Others"
Node_func_all[rownames(Node_func_all)%in%keystone$Label,]$Key_type <- "Keystone"
write.csv(Node_func_all, "./keystone/pearson/RMT-ASV_0.00001_0.01_node_all_Key_type-function.csv",quote = FALSE)
group_func_sum <- aggregate(Node_func_all[,c(24:39)], by=list(group=Node_func_all$Key_type),sum)
write.csv(group_func_sum, "./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_function_num.csv",quote = FALSE, row.names = FALSE)
group_func_sum <- read.csv("./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_function_num.csv", header=T, row.names=1)
# rowSums(group_func_sum)
# Keystone Others
# 115         7289
group_func_prop <- group_func_sum/rowSums(group_func_sum)
write.csv(group_func_prop, "./keystone/pearson/RMT-ASV_0.00001_0.01_Key_type_function_proportion.csv",quote = FALSE)
group_func_prop['group'] <- rownames(group_func_prop)
group_func_prop_long <- melt(group_func_prop, id = 'group')
group_func_prop_long$group<-factor(group_func_prop_long$group,levels = c("Keystone","Others"))
group_func_prop_long$variable<-factor(group_func_prop_long$variable,levels = c('Nitrifiers','Denitrifiers','Anammox','DNRA','PAOs','GAOs','SRB','SOB','Filamentous','Acetogen','Methanogen','Fermenter','Sugars','Lipids','Proteins','Unassigned'))

par(mar=c(3,5,2,2))
ggplot(group_func_prop_long, aes(group, value*100,fill= variable)) +
  geom_bar(stat = "identity")+
  labs(y = 'Proportion of functioner(%)') +
  scale_fill_manual(values = c('#008B8B','#556B2F','#228B22','#20B2AA','#483D8B','#6495ED','#DDA0DD','#D87093','#FF1493','#A0522D','#BC8F8F','#0066ff','#FFD700','#006699','#666699','#999999'),
                    limits = c('Nitrifiers','Denitrifiers','Anammox','DNRA','PAOs','GAOs','SRB','SOB','Filamentous','Acetogen','Methanogen','Fermenter','Sugars','Lipids','Proteins','Unassigned'))+
  ylim(-0.1,100.1) +
  theme(axis.text = element_text(size = 14,face = "plain",color ='black'), axis.title = element_text(size = 14,color ='black')) +
  theme(legend.text = element_text(size = 8,color ='black'))+
  #theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
  geom_col(position = 'stack', width = 0.6)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent')) +
  theme(legend.title = element_blank())
dev.off()
ggsave("./keystone/pearson/Key_type_function_proportion.pdf",height=6,width=8)

##node centrality
node_centrality <- read.csv("./igraph/pearson/RMT-ASV_0.00001_0.01_node_all_centrality.csv",header=T) ##共4992个ASVs
Node_all$Label <- rownames(Node_all)
all_node_centrality <- merge(Node_all, node_centrality, by="Label")
all_node_centrality$Key_type <- factor(all_node_centrality$Key_type, levels = c("Keystone","Others"))
##Degree_centrality
par(mar=c(3,5,2,2))
p<- ggplot(all_node_centrality, aes(x=Key_type, y=Degree_centrality, fill=Key_type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#F0E442", "#0072B2"),
                    limits = c("Keystone","Others")) +
  ylim(0,0.06)
p + labs(y = "Degree centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Key_Degree_centrality_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_node_centrality[all_node_centrality$Key_type=="Keystone",]$Degree_centrality, all_node_centrality[all_node_centrality$Key_type=="Others",]$Degree_centrality) ##p-value < 2.2e-16(0.029466523 0.004814908)

##Betweenness_centrality
par(mar=c(3,5,2,2))
p<- ggplot(all_node_centrality, aes(x=Key_type, y=Betweenness_centrality, fill=Key_type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#F0E442", "#0072B2"),
                    limits = c("Keystone","Others")) +
  ylim(0,0.04)
p + labs(y = "Betweenness centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Key_Betweenness_centrality_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_node_centrality[all_node_centrality$Key_type=="Keystone",]$Betweenness_centrality, all_node_centrality[all_node_centrality$Key_type=="Others",]$Betweenness_centrality) ##p-value = 0.6986(0.0008313496 0.0007819795)

##Closeness_centrality
par(mar=c(3,5,2,2))
p<- ggplot(all_node_centrality, aes(x=Key_type, y=Closeness_centrality, fill=Key_type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#F0E442", "#0072B2"),
                    limits = c("Keystone","Others")) +
  ylim(0,0.4)
p + labs(y = "Closeness centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Key_Closeness_centrality_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_node_centrality[all_node_centrality$Key_type=="Keystone",]$Closeness_centrality, all_node_centrality[all_node_centrality$Key_type=="Others",]$Closeness_centrality) ##p-value < 2.2e-16(0.2753382 0.1987975)

##Eigenvector_centrality
par(mar=c(3,5,2,2))
p<- ggplot(all_node_centrality, aes(x=Key_type, y=Eigenvector_centrality, fill=Key_type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#F0E442", "#0072B2"),
                    limits = c("Keystone","Others")) +
  ylim(0,1)
p + labs(y = "Eigenvector centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Key_Eigenvector_centrality_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_node_centrality[all_node_centrality$Key_type=="Keystone",]$Eigenvector_centrality, all_node_centrality[all_node_centrality$Key_type=="Others",]$Eigenvector_centrality) ##p-value < 2.2e-16(0.47949529 0.02590317)

##相关性
library(ggcor)
##env_influx
env <- read.csv("E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all-influx.csv", header = T, row.names = 1)
mantel_res <- read.csv("./keystone/pearson/ASV-data-all-environment-influx-all-keystone_others-mantel-test.csv", header = T)
mantel<-mutate(mantel_res,r = cut(abs(r), breaks = c(-Inf, 0.1, 0.2, Inf),
                                  labels = c("<0.1","0.1-0.2", ">=0.2"),
                                  right = FALSE),
               p.value = cut(p.value, breaks = c(-Inf,0.01, 0.05, Inf),
                             labels = c("<0.01", "0.01-0.05", ">=0.05"),
                             right = FALSE))
#绘制组合图
par(mar=c(3,5,2,2))
p1<-quickcor(env, type = "upper") + 
  geom_square(inherit.aes = TRUE)+ 
  anno_link(mantel, mapping = aes(colour = p.value, size = r),alpha = 0.6)+
  scale_size_manual(values = c(0.75, 1.5, 2),limits=c("<0.1","0.1-0.2",">=0.2"))+   #线段size
  scale_color_manual(values=c("darkgreen", "darkorange","darkgray"),
                     limits=c("<0.01", "0.01-0.05",">=0.05"))+  #线段颜色
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1), # 数据上下限
                       breaks = c(-1.0, -0.5, 0, 0.5, 1.0), # 分段点
                       low = "darkblue", mid = "white",high = "darkred", space = "Lab" ) #热图颜色
#geom_diag_label(geom = "text", remove.axis = TRUE) ##这个句决定env是否以横纵坐标展示
#定义图例标题
p1+theme(legend.position="right", legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.width = unit(1,"cm"))+
  guides(size=guide_legend(title="Mantel's r",override.aes=list(colour="grey35"),order=1),
         colour=guide_legend(title="Mantel's p",override.aes = list(size=0.75),order=2),
         fill=guide_colorbar(title="Pairwise correlations",barheight = 4,  order=2))
dev.off()
ggsave("./keystone/pearson/ASV-data-all-environment-influx-all-keystone_others-mantel-test.pdf", width = 9, height = 6)

##performance
env <- read.csv("E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_performance.csv", header = T, row.names = 1)
mantel_res <- read.csv("./keystone/pearson/ASV-data-all-performance-keystone_others-mantel-test.csv", header = T)
mantel<-mutate(mantel_res,r = cut(abs(r), breaks = c(-Inf, 0.1, 0.2,Inf),
                                  labels = c("<0.1","0.1-0.2",">=0.2"),
                                  right = FALSE),
               p.value = cut(p.value, breaks = c(-Inf, 0.01, 0.05, Inf),
                             labels = c("<0.01","0.01-0.05", ">=0.05"),
                             right = FALSE))
#绘制组合图
par(mar=c(3,5,2,2))
p1<-quickcor(env, type = "upper") + 
  geom_square(inherit.aes = TRUE)+ 
  anno_link(mantel, mapping = aes(colour = p.value, size = r),alpha = 0.6)+
  scale_size_manual(values = c(0.75, 1.5, 2),limits=c("<0.1","0.1-0.2",">=0.2"))+   #线段size
  scale_color_manual(values=c("darkgreen", "darkorange","darkgray"),
                     limits=c("<0.01", "0.01-0.05",">=0.05"))+  #线段颜色
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1), # 数据上下限
                       breaks = c(-1.0, -0.5, 0, 0.5, 1.0), # 分段点
                       low = "darkblue", mid = "white",high = "darkred", space = "Lab" ) #热图颜色
#geom_diag_label(geom = "text", remove.axis = TRUE) ##这个句决定env是否以横纵坐标展示
#定义图例标题
p1+theme(legend.position="right", legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.width = unit(1,"cm"))+
  guides(size=guide_legend(title="Mantel's r",override.aes=list(colour="grey35"),order=1),
         colour=guide_legend(title="Mantel's p",override.aes = list(size=0.75),order=2),
         fill=guide_colorbar(title="Pairwise correlations",barheight = 4,  order=3))
dev.off()
ggsave("./keystone/pearson/ASV-data-all-performance-keystone_others-mantel-test.pdf", width = 9, height = 6)
