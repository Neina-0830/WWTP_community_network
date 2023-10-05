##获得边的统计结果
library(ggplot2)
library(statnet)
library(circlize)

my_edge <- read.csv("../data/RMT_ASV_0.00001_0.01_edge.csv",header=T) ##65457个边
my_node <- read.csv("../data/RMT-ASV_0.00001_0.01_node_all.csv",header=T) ##4992个ASV
edge_all <- my_edge[,c("Source","Target")]
node_source <- my_node[my_node$Id %in% edge_all$Source,] 
node_source_temp <- merge(edge_all, node_source, by.x = 'Source',by.y = 'Id', all.x = TRUE)
edge_all_new <- node_source_temp[,c("Source","Target","v_phylum","v_genus","v_partition","v_group","v_taxa")]
colnames(edge_all_new) <- c("Source","Target","Source_phylum","Source_genus","Source_partition","Source_group","Source_taxa")
node_target <- my_node[my_node$Id %in% edge_all$Target,] 
node_target_temp <- merge(edge_all_new, node_target, by.x = 'Target',by.y = 'Id', all.x = TRUE)
edge_all_new <- node_target_temp[,c("Source","Target","Source_phylum","Source_genus","Source_partition","Source_group","Source_taxa","v_phylum","v_genus","v_partition","v_group","v_taxa")]
colnames(edge_all_new) <- c("Source","Target","Source_phylum","Source_genus","Source_partition","Source_group","Source_taxa","Target_phylum","Target_genus","Target_partition","Target_group","Target_taxa")
write.csv(edge_all_new, './igraph/pearson/RMT-ASV_0.00001_0.01_edge_taxonomy.csv',quote = FALSE,row.names = FALSE)

edge_all_new <- read.csv("./igraph/pearson/RMT-ASV_0.00001_0.01_edge_taxonomy.csv",header=T)
##Phylum统计
my_node_phylum <- as.data.frame(table(my_node$v_phylum))
write.csv(my_node_phylum, './igraph/pearson/RMT-ASV_0.00001_0.01_node_Phylum-summary.csv',quote = FALSE,row.names = FALSE)
top_phylum_list <- c("Proteobacteria","Bacteroidota","Myxococcota","Chloroflexi","Planctomycetota","Acidobacteriota",
                     "Patescibacteria","Verrucomicrobiota","Actinobacteriota","Bdellovibrionota","Desulfobacterota")
edge_all_new[!edge_all_new$Source_phylum %in% top_phylum_list,]$Source_phylum <- "Others"
edge_all_new[!edge_all_new$Target_phylum %in% top_phylum_list,]$Target_phylum <- "Others"
all_phylum <- c("Proteobacteria","Bacteroidota","Myxococcota","Chloroflexi","Planctomycetota","Acidobacteriota","Patescibacteria","Verrucomicrobiota","Actinobacteriota","Bdellovibrionota","Desulfobacterota","Others")
pair <- as.data.frame(t(combn(all_phylum, 2)))
colnames(pair) <- c("Node1","Node2")
same <- as.data.frame(cbind(all_phylum,all_phylum))
colnames(same) <- c("Node1","Node2")
edge_phylum <- as.data.frame(table(edge_all_new[c("Source_phylum", "Target_phylum")]))
##不区分边的方向重新整理边
pair$Freq <- 0
for(i in 1:nrow(pair)){
  pair[i,]$Freq <- edge_phylum[edge_phylum$Source_phylum==pair[i,1] & edge_phylum$Target_phylum==pair[i,2],]$Freq + edge_phylum[edge_phylum$Source_phylum==pair[i,2] & edge_phylum$Target_phylum==pair[i,1],]$Freq
}
same$Freq <- 0
for(i in 1:nrow(same)){
  same[i,]$Freq <- edge_phylum[edge_phylum$Source_phylum==same[i,1] & edge_phylum$Target_phylum==same[i,2],]$Freq
}
all_pair <- rbind(pair,same)
write.csv(all_pair, './igraph/pearson/RMT-ASV_0.00001_0.01_edge_Phylum-summary.csv',quote = FALSE,row.names = FALSE)
sum(all_pair[all_pair$Node1=="Proteobacteria",]$Freq,all_pair[all_pair$Node2=="Proteobacteria",]$Freq) #34462
sum(all_pair[all_pair$Node1=="Bacteroidota",]$Freq,all_pair[all_pair$Node2=="Bacteroidota",]$Freq) #33853
sum(all_pair[all_pair$Node1=="Myxococcota",]$Freq,all_pair[all_pair$Node2=="Myxococcota",]$Freq) #8531
sum(all_pair[all_pair$Node1=="Chloroflexi",]$Freq,all_pair[all_pair$Node2=="Chloroflexi",]$Freq) #7932
sum(all_pair[all_pair$Node1=="Planctomycetota",]$Freq,all_pair[all_pair$Node2=="Planctomycetota",]$Freq) #7340
sum(all_pair[all_pair$Node1=="Acidobacteriota",]$Freq,all_pair[all_pair$Node2=="Acidobacteriota",]$Freq) #6485
sum(all_pair[all_pair$Node1=="Patescibacteria",]$Freq,all_pair[all_pair$Node2=="Patescibacteria",]$Freq) #5150
sum(all_pair[all_pair$Node1=="Verrucomicrobiota",]$Freq,all_pair[all_pair$Node2=="Verrucomicrobiota",]$Freq) #4370
sum(all_pair[all_pair$Node1=="Actinobacteriota",]$Freq,all_pair[all_pair$Node2=="Actinobacteriota",]$Freq) #3045
sum(all_pair[all_pair$Node1=="Bdellovibrionota",]$Freq,all_pair[all_pair$Node2=="Bdellovibrionota",]$Freq) #3777
sum(all_pair[all_pair$Node1=="Desulfobacterota",]$Freq,all_pair[all_pair$Node2=="Desulfobacterota",]$Freq) #2922
sum(all_pair[all_pair$Node1=="Others",]$Freq,all_pair[all_pair$Node2=="Others",]$Freq) #13047


# 设置颜色、图片文件名、长宽和字体大小
grid.col = NULL
grid.col[all_phylum] = c("#006C71", "#FF8003","#FF82FF","#3CDB00","#814A18","#00DCFF","#C8AF00","#00D19A","#72466B","#4F8D32","#FF5E71","#C0C0C0")
##画图
pdf(file="./igraph/pearson/RMT-ASV_0.00001_0.01_edge_Phylum-circlize-new.pdf", width=8, height=5, pointsize=8)
chordDiagram(all_pair, order = c("Proteobacteria","Bacteroidota","Myxococcota","Chloroflexi","Planctomycetota","Acidobacteriota",
                                    "Patescibacteria","Verrucomicrobiota","Actinobacteriota","Bdellovibrionota","Desulfobacterota","Others"),
             grid.col = grid.col, 
             transparency = 0.5)
legend("right",pch=20,legend=all_phylum,
       col=grid.col[all_phylum],bty="n",
       cex=1,pt.cex=3,border="black")
circos.clear()
dev.off()
