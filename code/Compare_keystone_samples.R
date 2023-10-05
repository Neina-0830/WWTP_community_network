library(ggplot2)
library(igraph)
library(dplyr)
library(reshape2)
data_summary <-  function(data=NULL, varname, groupnames=NULL, na.rm=FALSE,conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  #group_by + summarise
  datac <- ddply(data, groupnames, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 varname
  )
  datac <- plyr::rename(datac, c("mean" = varname))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}

##keystone类群在样本中的分布情况
keystone <- read.csv('./keystone/pearson/RMT-ASV_0.00001_0.01_node_keystone.csv',header = T, row.names = 1)
ASV_data <- read.csv('../data/ASV_data_0.00001_0.01.csv', row.names = 1,header=TRUE)
keystone_data <- ASV_data[keystone$Label,]
#包含keystone的样本
keystone_sample <- keystone_data[,colSums(keystone_data)>0]
write.csv(keystone_sample, './keystone/pearson/RMT-ASV_0.00001_0.01_keystone_samples.csv',quote=FALSE)
ncol(keystone_sample)/ncol(keystone_data)
##[1] 0.4418212
no_keystone_sample <- keystone_data[,colSums(keystone_data)==0]
write.csv(no_keystone_sample, './keystone/pearson/RMT-ASV_0.00001_0.01_no_keystone_samples.csv',quote=FALSE)
keystone_sample_data <- ASV_data[,colnames(ASV_data)%in%colnames(keystone_sample)]
keystone_sample_data <- keystone_sample_data[rowSums(keystone_sample_data)>0,] ##6600个ASVs
no_keystone_sample_data <- ASV_data[,colnames(ASV_data)%in%colnames(no_keystone_sample)]
no_keystone_sample_data <- no_keystone_sample_data[rowSums(no_keystone_sample_data)>0,] ##6535个ASVs
write.csv(keystone_sample_data, './keystone/pearson/RMT-ASV_0.00001_0.01_keystone_sample_data.csv',quote=FALSE)
write.csv(no_keystone_sample_data, './keystone/pearson/RMT-ASV_0.00001_0.01_no_keystone_sample_data.csv',quote=FALSE)


##从总网络中筛选包含/不包含keystone的样本的子网络
#输入数据示例，邻接矩阵
#这是一个微生物互作网络，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adjacency_unweight <- read.csv('./igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_adj.csv', row.names = 1, check.names = FALSE)
#邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
mynet <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
keystone_sample_node <- rownames(keystone_sample_data[rownames(keystone_sample_data) %in% V(mynet)$name,]) ##V(mynet)$name获得节点名
keystone_sample_net <- induced_subgraph(mynet, keystone_sample_node)
taxonomy = read.csv("./rawdata/ASV-MiDAS-MRA-function-taxonomy-all-new.csv",header=T,row.names=1)  ##Taxonomy与data的顺序问题
##按照V(keystone_sample_net)$name筛选且排序rownames(taxonomy)
loc = match(V(keystone_sample_net)$name, rownames(taxonomy))
taxonomy <- taxonomy[loc,]
##set vertices size and color
V(keystone_sample_net)$Kindom <- taxonomy$Kindom
V(keystone_sample_net)$Phylum <- taxonomy$Phylum
V(keystone_sample_net)$Class <- taxonomy$Class
V(keystone_sample_net)$Order <- taxonomy$Order
V(keystone_sample_net)$Family <- taxonomy$Family
V(keystone_sample_net)$Genus <- taxonomy$Genus
V(keystone_sample_net)$Species <- taxonomy$Species
V(keystone_sample_net)$Partition <- taxonomy$Partition
V(keystone_sample_net)$Group <- taxonomy$Group
V(keystone_sample_net)$Taxa <- taxonomy$Taxa
V(keystone_sample_net)$Function <- taxonomy$Function
V(keystone_sample_net)$Process <- taxonomy$Process
V(keystone_sample_net)$Frequency <- taxonomy$frequency
V(keystone_sample_net)$Niche <- taxonomy$niche_width
V(keystone_sample_net)$MRA <- taxonomy$MRA
n <- length(V(keystone_sample_net))
mycolors <- colorRampPalette(colors())(n)
V(keystone_sample_net)$color <- mycolors
##remove isolated nodes
bad.vs = V(keystone_sample_net)[degree(keystone_sample_net) == 0]
keystone_sample_net = delete.vertices(keystone_sample_net, bad.vs)
matrix_rep <- as.matrix(get.adjacency(keystone_sample_net)) #this gives you the adjacency
node <- length(V(keystone_sample_net))###Number of nodes
edge <- length(E(keystone_sample_net))###Number of edges
set.seed(111)
V(keystone_sample_net)$modularity <- membership(cluster_louvain(keystone_sample_net))
V(keystone_sample_net)$degree <- degree(keystone_sample_net)####Degree
V(keystone_sample_net)$betweenness <- betweenness(keystone_sample_net)####Betweenness centrality
V(keystone_sample_net)$closeness <- closeness(keystone_sample_net)####Closeness centrality
write.csv(matrix_rep, "./keystone/pearson/RMT-ASV_level_0.00001_0.01_network_keystone_samples_pearson_adj.csv",quote = FALSE)
write_graph(keystone_sample_net, "./keystone/pearson/RMT_ASV_0.00001_0.01-Network-keystone_samples_pearson.graphml","graphml")
###Caculation of node-level topological features
modularity.class <- V(keystone_sample_net)$modularity
node.degree <- degree(keystone_sample_net)####Degree
betweenness.centrality <- betweenness(keystone_sample_net)####Betweenness centrality
closeness.centrality <- closeness(keystone_sample_net)####Closeness centrality
node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality,modularity.class)
write.csv(node.topology, file="./keystone/pearson/RMT-Net_0.00001_0.01_node_topology_keystone_samples_pearson.csv",quote = FALSE)
###Caculation of network-level topological features
node <- length(V(keystone_sample_net))###Number of nodes  5503
edge <- length(E(keystone_sample_net))###Number of edges 86753
ave_degree <- mean(degree(keystone_sample_net))###Average degree 
clust_coeff <- transitivity(keystone_sample_net)###Clustering coefficient
ave_path_len <- average.path.length(keystone_sample_net)###Average path length
diameter <- diameter(keystone_sample_net)###Network diameter
density <- graph.density(keystone_sample_net)###Graph density
modular <- modularity(cluster_louvain(keystone_sample_net))  ##Graph modularity
betweenness_centr <- centralization.betweenness(keystone_sample_net)$centralization
global.topology <- data.frame(node,edge,ave_degree,clust_coeff,ave_path_len,diameter,density,modular,betweenness_centr)
write.csv(global.topology, file="./keystone/pearson/RMT-Net_0.00001_0.01_global_topology_keystone_samples_pearson.csv",quote = FALSE)

no_keystone_sample_node <- rownames(no_keystone_sample_data[rownames(no_keystone_sample_data) %in% V(mynet)$name,]) ##V(mynet)$name获得节点名
no_keystone_sample_net <- induced_subgraph(mynet, no_keystone_sample_node)
taxonomy = read.csv("../data/ASV-MiDAS-MRA-function-taxonomy-all-new.csv",header=T,row.names=1)  ##Taxonomy与data的顺序问题
##按照V(no_keystone_sample_net)$name筛选且排序rownames(taxonomy)
loc = match(V(no_keystone_sample_net)$name, rownames(taxonomy))
taxonomy <- taxonomy[loc,]
##set vertices size and color
V(no_keystone_sample_net)$Kindom <- taxonomy$Kindom
V(no_keystone_sample_net)$Phylum <- taxonomy$Phylum
V(no_keystone_sample_net)$Class <- taxonomy$Class
V(no_keystone_sample_net)$Order <- taxonomy$Order
V(no_keystone_sample_net)$Family <- taxonomy$Family
V(no_keystone_sample_net)$Genus <- taxonomy$Genus
V(no_keystone_sample_net)$Species <- taxonomy$Species
V(no_keystone_sample_net)$Partition <- taxonomy$Partition
V(no_keystone_sample_net)$Group <- taxonomy$Group
V(no_keystone_sample_net)$Taxa <- taxonomy$Taxa
V(no_keystone_sample_net)$Function <- taxonomy$Function
V(no_keystone_sample_net)$Process <- taxonomy$Process
V(no_keystone_sample_net)$Frequency <- taxonomy$frequency
V(no_keystone_sample_net)$Niche <- taxonomy$niche_width
V(no_keystone_sample_net)$MRA <- taxonomy$MRA
n <- length(V(no_keystone_sample_net))
mycolors <- colorRampPalette(colors())(n)
V(no_keystone_sample_net)$color <- mycolors
##remove isolated nodes
bad.vs = V(no_keystone_sample_net)[degree(no_keystone_sample_net) == 0]
no_keystone_sample_net = delete.vertices(no_keystone_sample_net, bad.vs)
matrix_rep <- as.matrix(get.adjacency(no_keystone_sample_net)) #this gives you the adjacency
node <- length(V(no_keystone_sample_net))###Number of nodes
edge <- length(E(no_keystone_sample_net))###Number of edges
set.seed(111)
V(no_keystone_sample_net)$modularity <- membership(cluster_louvain(no_keystone_sample_net))
V(no_keystone_sample_net)$degree <- degree(no_keystone_sample_net)####Degree
V(no_keystone_sample_net)$betweenness <- betweenness(no_keystone_sample_net)####Betweenness centrality
V(no_keystone_sample_net)$closeness <- closeness(no_keystone_sample_net)####Closeness centrality
write.csv(matrix_rep, "./keystone/pearson/RMT-ASV_level_0.00001_0.01_network_no_keystone_samples_pearson_adj.csv",quote = FALSE)
write_graph(no_keystone_sample_net, "./keystone/pearson/RMT_ASV_0.00001_0.01-Network-no_keystone_samples_pearson.graphml","graphml")
###Caculation of node-level topological features
modularity.class <- V(no_keystone_sample_net)$modularity
node.degree <- degree(no_keystone_sample_net)####Degree
betweenness.centrality <- betweenness(no_keystone_sample_net)####Betweenness centrality
closeness.centrality <- closeness(no_keystone_sample_net)####Closeness centrality
node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality,modularity.class)
write.csv(node.topology, file="./keystone/pearson/RMT-Net_0.00001_0.01_node_topology_no_keystone_samples_pearson.csv",quote = FALSE)
###Caculation of network-level topological features
node <- length(V(no_keystone_sample_net))###Number of nodes  5503
edge <- length(E(no_keystone_sample_net))###Number of edges 86753
ave_degree <- mean(degree(no_keystone_sample_net))###Average degree 
clust_coeff <- transitivity(no_keystone_sample_net)###Clustering coefficient
ave_path_len <- average.path.length(no_keystone_sample_net)###Average path length
diameter <- diameter(no_keystone_sample_net)###Network diameter
density <- graph.density(no_keystone_sample_net)###Graph density
modular <- modularity(cluster_louvain(no_keystone_sample_net))  ##Graph modularity
betweenness_centr <- centralization.betweenness(no_keystone_sample_net)$centralization
global.topology <- data.frame(node,edge,ave_degree,clust_coeff,ave_path_len,diameter,density,modular,betweenness_centr)
write.csv(global.topology, file="./keystone/pearson/RMT-Net_0.00001_0.01_global_topology_no_keystone_samples_pearson.csv",quote = FALSE)

##计算模块内连通度（Zi）和模块间连通度（Pi）
zi.pi<-function(adj.matrix,node.topology){ ##节点属性node.topology，邻接矩阵adj.matrix
  module<-which(colnames(node.topology)=='modularity.class')
  module.max<-max(node.topology[,module])
  degree<-which(colnames(node.topology)=='node.degree')
  #按照模块将相关矩阵分割
  bulk.module<-list(NA)
  length(bulk.module)<-module.max
  for(i in 1:max(node.topology[,module])){
    bulk.module[[i]]<-adj.matrix[which(node.topology[,module]==i),which(node.topology[,module]==i)]
    bulk.module[[i]]<-as.data.frame(bulk.module[[i]])
    rownames(bulk.module[[i]])<-rownames(adj.matrix)[which(node.topology[,module]==i)]
    colnames(bulk.module[[i]])<-colnames(adj.matrix)[which(node.topology[,module]==i)]
  }
  # within-module degree z
  z_bulk<-list(NA)
  length(z_bulk)<-module.max
  for(i in 1:length(z_bulk)){
    z_bulk[[i]]<-bulk.module[[i]][,1]
    z_bulk[[i]]<-as.data.frame(z_bulk[[i]])
    colnames(z_bulk[[i]])<-"z"
    rownames(z_bulk[[i]])<-rownames(bulk.module[[i]])
  }
  #计算z值
  for(i in 1:max(node.topology[,module])){
    if(length(bulk.module[[i]])==1){
      z_bulk[[i]][,1]<-0
    }else if(sum(bulk.module[[i]])==0){
      z_bulk[[i]][,1]<-0
    }else{
      k<-rowSums(bulk.module[[i]])
      mean<-mean(k)
      sd<-sd(k)
      if (sd==0){
        z_bulk[[i]][,1]<-0
      }else{
        z_bulk[[i]][,1]<-(k-mean)/sd
      }
    }
  }
  #z值合并
  for(i in 2:max(node.topology[,module])) {
    z_bulk[[i]]<-rbind(z_bulk[[i-1]],z_bulk[[i]])
  }
  z_bulk<-z_bulk[[module.max]]
  #按照模块将相关矩阵列分割
  bulk.module1<-list(NA)
  length(bulk.module1)<-module.max
  for(i in 1:max(node.topology[,module])){
    bulk.module1[[i]]<-adj.matrix[,which(node.topology[,module]==i)]
    bulk.module1[[i]]<-as.data.frame(bulk.module1[[i]])
    rownames(bulk.module1[[i]])<-rownames(adj.matrix)
    colnames(bulk.module1[[i]])<-colnames(adj.matrix)[which(node.topology[,module]==i)]
  }
  #among-module connectivity c
  c_bulk<-list(NA)
  length(c_bulk)<-module.max
  for(i in 1:length(c_bulk)){
    c_bulk[[i]]<-adj.matrix[,1]
    c_bulk[[i]]<-as.matrix(c_bulk[[i]])
    colnames(c_bulk[[i]])<-"c"
    rownames(c_bulk[[i]])<-rownames(adj.matrix)
    c_bulk[[i]][,1]<-NA
  }
  #每个节点各模块连接数平方
  for(i in 1:max(node.topology[,module])){
    c_bulk[[i]]<-rowSums(bulk.module1[[i]])
    c_bulk[[i]]<-as.matrix(c_bulk[[i]])
    c_bulk[[i]]<-c_bulk[[i]]*c_bulk[[i]]
    colnames(c_bulk[[i]])<-"c"
    rownames(c_bulk[[i]])<-rownames(adj.matrix)
  }
  #平方和
  for(i in 2:max(node.topology[,module])){
    c_bulk[[i]]<-c_bulk[[i]]+c_bulk[[i-1]]
  }
  c_bulk<-c_bulk[[module.max]]
  
  c_bulk1<-1-(c_bulk/(node.topology[,degree]*node.topology[,degree]))
  colnames(c_bulk1)<-"c"
  #z,c整合
  z_c_bulk<-c_bulk1
  z_c_bulk<-as.data.frame(z_c_bulk)
  z_c_bulk$z<-z_bulk[match(rownames(c_bulk1),rownames(z_bulk)),]
  z_c_bulk<-z_c_bulk[,c(2,1)]
  names(z_c_bulk)[1:2]<-c('within_module_connectivities','among_module_connectivities')
  z_c_bulk$nodes_id<-rownames(z_c_bulk)
  node.topology$nodes_id<-rownames(node.topology)
  z_c_bulk<-merge(z_c_bulk,node.topology,by='nodes_id')
  z_c_bulk
}
adj.matrix <- read.csv('./keystone/pearson/RMT-ASV_level_0.00001_0.01_network_keystone_samples_pearson_adj.csv', row.names = 1, check.names = FALSE)
node.topology <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_node_topology_keystone_samples_pearson.csv', row.names = 1, check.names = FALSE)
#两个文件的节点顺序要一致
node.topology <- node.topology[rownames(adj.matrix), ]
#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(adj.matrix, node.topology)
head(zi_pi)
write.csv(zi_pi, './keystone/pearson/RMT-Net_0.00001_0.01_node_topology_keystone_samples_pearson_zi_pi.csv',quote=FALSE)
##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)
zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.csv(zi_pi, './keystone/pearson/RMT-Net_0.00001_0.01_node_topology_keystone_samples_pearson_zi_pi-type.csv',quote=FALSE, row.names = FALSE)
par(mar=c(3,5,2,2))
ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  xlim(0,0.8) +
  ylim(-2.5,7.5) +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-ASV_level_0.00001_0.01_network_keystone_samples-Hubs.pdf", width = 5, height = 3.5)

adj.matrix <- read.csv('./keystone/pearson/RMT-ASV_level_0.00001_0.01_network_no_keystone_samples_pearson_adj.csv', row.names = 1, check.names = FALSE)
node.topology <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_node_topology_no_keystone_samples_pearson.csv', row.names = 1, check.names = FALSE)
#两个文件的节点顺序要一致
node.topology <- node.topology[rownames(adj.matrix), ]
#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(adj.matrix, node.topology)
head(zi_pi)
write.csv(zi_pi, './keystone/pearson/RMT-Net_0.00001_0.01_node_topology_no_keystone_samples_pearson_zi_pi.csv',quote=FALSE)
##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)
zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.csv(zi_pi, './keystone/pearson/RMT-Net_0.00001_0.01_node_topology_no_keystone_samples_pearson_zi_pi-type.csv',quote=FALSE, row.names = FALSE)
par(mar=c(3,5,2,2))
ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 1) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
                     limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  xlim(0,0.8) +
  ylim(-2.5,7.5) +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-ASV_level_0.00001_0.01_network_no_keystone_samples-Hubs.pdf", width = 5, height = 3.5)



##包含/不包含keystone类群的样本的差异分析
keystone_sample <- read.csv('./keystone/pearson/RMT-ASV_0.00001_0.01_keystone_samples.csv',header = T, row.names = 1)
no_keystone_sample <- read.csv('./keystone/pearson/RMT-ASV_0.00001_0.01_no_keystone_samples.csv',header = T, row.names = 1)
keystone_sample <- as.data.frame(t(keystone_sample))
keystone_sample$group <- "Keystone_sample"
no_keystone_sample <- as.data.frame(t(no_keystone_sample))
no_keystone_sample$group <- "No_keystone_sample"
all_keystone_sample <- rbind(keystone_sample, no_keystone_sample)
group=as.data.frame(all_keystone_sample$group)
rownames(group)=rownames(all_keystone_sample)
group$sample=rownames(group)
colnames(group)=c("group","sample")
write.csv(group, "./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group.csv", quote=FALSE)

##环境和结构变异
##包含Keystone类群的样本
keystone_sample <- read.csv('./keystone/pearson/RMT-ASV_0.00001_0.01_keystone_samples.csv', row.names = 1,header=TRUE)
performance <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_performance.csv', row.names = 1,header=TRUE)
performance_keystone_sample <- performance[colnames(keystone_sample),]
write.csv(performance_keystone_sample, "./keystone/pearson/WWTP_meta_predata_all_performance-keystone_sample.csv", quote = FALSE)
SubNet_all <- read.csv('./igraph/pearson/RMT-Subnet_0.00001_0.01_pearson_global_topology.csv', row.names = 1,header=TRUE)
SubNet_keystone_sample <- SubNet_all[colnames(keystone_sample),]
write.csv(SubNet_keystone_sample, "./keystone/pearson/RMT-Subnet_0.00001_0.01_pearson_global_topology-keystone_sample.csv", quote = FALSE)
ENV_influx <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env-influx.csv', row.names = 1,header=TRUE)
ENV_keystone_sample <- ENV_influx[colnames(keystone_sample),]
write.csv(ENV_keystone_sample, "./keystone/pearson/WWTP_meta_predata_all_env-influx-keystone_sample.csv", quote = FALSE)
ENV_keystone_sample_standard <- decostand(ENV_keystone_sample, 'standardize')   #可选 z-score 标准化，消除量纲差异
ENV_keystone_sample_dis <- as.matrix(vegdist(ENV_keystone_sample_standard, method = 'euclidean'))
write.csv(ENV_keystone_sample_dis, "./keystone/pearson/WWTP_meta_predata_all_env-influx-keystone_sample-standardize-euclidean-matrix.csv", quote = FALSE)
ENV_keystone_sample_dis[!upper.tri(ENV_keystone_sample_dis, diag = FALSE)] <- 100000
ENV_keystone_sample_dis_long <- melt(ENV_keystone_sample_dis, id = 'sample')
colnames(ENV_keystone_sample_dis_long) = c('sample1','sample2','Env')
ENV_keystone_sample_dis_long <- ENV_keystone_sample_dis_long[ENV_keystone_sample_dis_long$Env!=100000,]
write.table(ENV_keystone_sample_dis_long, "./keystone/pearson/WWTP_meta_predata_all_env-influx-keystone_sample-standardize-euclidean-tran.txt", quote = FALSE, row.names = FALSE, sep = '\t')
##读取AS系统微生物群落的weighted_unifrac距离
Weight_Unifrac <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_weighted_uniFrac.csv', header=TRUE, row.names=1)
Weight_Unifrac_keystone_sample <- Weight_Unifrac[colnames(keystone_sample),colnames(keystone_sample)]
write.csv(Weight_Unifrac_keystone_sample, "./keystone/pearson/WWTP_ASV_weighted_uniFrac-keystone_sample.csv", quote = FALSE)
Weight_Unifrac_keystone_sample[!upper.tri(Weight_Unifrac_keystone_sample, diag = FALSE)] <- 100000
Weight_Unifrac_keystone_sample_long <- melt(as.matrix(Weight_Unifrac_keystone_sample), id = 'sample')
colnames(Weight_Unifrac_keystone_sample_long) = c('sample1','sample2','weighted_Unifrac')
Weight_Unifrac_keystone_sample_long <- Weight_Unifrac_keystone_sample_long[Weight_Unifrac_keystone_sample_long$weighted_Unifrac!=100000,]
##合并环境距离和群落w_unifrac距离
Env_weight_Unifrac <- merge(ENV_keystone_sample_dis_long,Weight_Unifrac_keystone_sample_long, by=c('sample1','sample2'))
write.table(Env_weight_Unifrac, "./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-weighted_Unifrac-keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
##读取AS系统微生物群落的bray_curtis距离
bray_curtis <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_Bray.csv', header=TRUE, row.names=1)
bray_curtis_keystone_sample <- bray_curtis[colnames(keystone_sample),colnames(keystone_sample)]
write.csv(bray_curtis_keystone_sample, "./keystone/pearson/WWTP_ASV_bray_curtis-keystone_sample.csv", quote = FALSE)
bray_curtis_keystone_sample[!upper.tri(bray_curtis_keystone_sample, diag = FALSE)] <- 100000
bray_curtis_keystone_sample_long <- melt(as.matrix(bray_curtis_keystone_sample), id = 'sample')
colnames(bray_curtis_keystone_sample_long) = c('sample1','sample2','bray_curtis')
bray_curtis_keystone_sample_long <- bray_curtis_keystone_sample_long[bray_curtis_keystone_sample_long$bray_curtis!=100000,]
write.table(bray_curtis_keystone_sample_long, "./keystone/pearson/ASV-bray_curtis-keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
##读取AS系统微生物群落的jaccard距离
jaccard <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_Jaccard.csv', header=TRUE, row.names=1)
jaccard_keystone_sample <- jaccard[colnames(keystone_sample),colnames(keystone_sample)]
write.csv(jaccard_keystone_sample, "./keystone/pearson/WWTP_ASV_jaccard-keystone_sample.csv", quote = FALSE)
jaccard_keystone_sample[!upper.tri(jaccard_keystone_sample, diag = FALSE)] <- 100000
jaccard_keystone_sample_long <- melt(as.matrix(jaccard_keystone_sample), id = 'sample')
colnames(jaccard_keystone_sample_long) = c('sample1','sample2','jaccard')
jaccard_keystone_sample_long <- jaccard_keystone_sample_long[jaccard_keystone_sample_long$jaccard!=100000,]
write.table(jaccard_keystone_sample_long, "./keystone/pearson/ASV-jaccard-keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
##读取AS系统微生物群落的unweighted_Unifrac距离
unweighted_Unifrac <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_unweight_uniFrac.csv', header=TRUE, row.names=1)
unweighted_Unifrac_keystone_sample <- unweighted_Unifrac[colnames(keystone_sample),colnames(keystone_sample)]
write.csv(unweighted_Unifrac_keystone_sample, "./keystone/pearson/WWTP_ASV_unweighted_Unifrac-keystone_sample.csv", quote = FALSE)
unweighted_Unifrac_keystone_sample[!upper.tri(unweighted_Unifrac_keystone_sample, diag = FALSE)] <- 100000
unweighted_Unifrac_keystone_sample_long <- melt(as.matrix(unweighted_Unifrac_keystone_sample), id = 'sample')
colnames(unweighted_Unifrac_keystone_sample_long) = c('sample1','sample2','unweighted_Unifrac')
unweighted_Unifrac_keystone_sample_long <- unweighted_Unifrac_keystone_sample_long[unweighted_Unifrac_keystone_sample_long$unweighted_Unifrac!=100000,]
write.table(unweighted_Unifrac_keystone_sample_long, "./keystone/pearson/ASV-unweighted_Unifrac-keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)


##不包含Keystone类群的样本
no_keystone_sample <- read.csv('./keystone/pearson/RMT-ASV_0.00001_0.01_no_keystone_samples.csv', row.names = 1,header=TRUE)
performance <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_performance.csv', row.names = 1,header=TRUE)
performance_no_keystone_sample <- performance[colnames(no_keystone_sample),]
write.csv(performance_no_keystone_sample, "./keystone/pearson/WWTP_meta_predata_all_performance-no_keystone_sample.csv", quote = FALSE)
SubNet_all <- read.csv('./igraph/pearson/RMT-Subnet_0.00001_0.01_pearson_global_topology.csv', row.names = 1,header=TRUE)
SubNet_no_keystone_sample <- SubNet_all[colnames(no_keystone_sample),]
write.csv(SubNet_no_keystone_sample, "./keystone/pearson/RMT-Subnet_0.00001_0.01_pearson_global_topology-no_keystone_sample.csv", quote = FALSE)
ENV_influx <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env-influx.csv', row.names = 1,header=TRUE)
ENV_no_keystone_sample <- ENV_influx[colnames(no_keystone_sample),]
write.csv(ENV_no_keystone_sample, "./keystone/pearson/WWTP_meta_predata_all_env-influx-no_keystone_sample.csv", quote = FALSE)
ENV_no_keystone_sample_standard <- decostand(ENV_no_keystone_sample, 'standardize')   #可选 z-score 标准化，消除量纲差异
ENV_no_keystone_sample_dis <- as.matrix(vegdist(ENV_no_keystone_sample_standard, method = 'euclidean'))
write.csv(ENV_no_keystone_sample_dis, "./keystone/pearson/WWTP_meta_predata_all_env-influx-no_keystone_sample-standardize-euclidean-matrix.csv", quote = FALSE)
ENV_no_keystone_sample_dis[!upper.tri(ENV_no_keystone_sample_dis, diag = FALSE)] <- 100000
ENV_no_keystone_sample_dis_long <- melt(ENV_no_keystone_sample_dis, id = 'sample')
colnames(ENV_no_keystone_sample_dis_long) = c('sample1','sample2','Env')
ENV_no_keystone_sample_dis_long <- ENV_no_keystone_sample_dis_long[ENV_no_keystone_sample_dis_long$Env!=100000,]
write.table(ENV_no_keystone_sample_dis_long, "./keystone/pearson/WWTP_meta_predata_all_env-influx-no_keystone_sample-standardize-euclidean-tran.txt", quote = FALSE, row.names = FALSE, sep = '\t')
##读取AS系统微生物群落的weighted_unifrac距离
Weight_Unifrac <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_weighted_uniFrac.csv', header=TRUE, row.names=1)
Weight_Unifrac_no_keystone_sample <- Weight_Unifrac[colnames(no_keystone_sample),colnames(no_keystone_sample)]
write.csv(Weight_Unifrac_no_keystone_sample, "./keystone/pearson/WWTP_ASV_weighted_uniFrac-no_keystone_sample.csv", quote = FALSE)
Weight_Unifrac_no_keystone_sample[!upper.tri(Weight_Unifrac_no_keystone_sample, diag = FALSE)] <- 100000
Weight_Unifrac_no_keystone_sample_long <- melt(as.matrix(Weight_Unifrac_no_keystone_sample), id = 'sample')
colnames(Weight_Unifrac_no_keystone_sample_long) = c('sample1','sample2','weighted_Unifrac')
Weight_Unifrac_no_keystone_sample_long <- Weight_Unifrac_no_keystone_sample_long[Weight_Unifrac_no_keystone_sample_long$weighted_Unifrac!=100000,]
##合并环境距离和群落w_unifrac距离
Env_weight_Unifrac <- merge(ENV_no_keystone_sample_dis_long,Weight_Unifrac_no_keystone_sample_long, by=c('sample1','sample2'))
write.table(Env_weight_Unifrac, "./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-weighted_Unifrac-no_keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
##读取AS系统微生物群落的bray_curtis距离
bray_curtis <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_Bray.csv', header=TRUE, row.names=1)
bray_curtis_no_keystone_sample <- bray_curtis[colnames(no_keystone_sample),colnames(no_keystone_sample)]
write.csv(bray_curtis_no_keystone_sample, "./keystone/pearson/WWTP_ASV_bray_curtis-no_keystone_sample.csv", quote = FALSE)
bray_curtis_no_keystone_sample[!upper.tri(bray_curtis_no_keystone_sample, diag = FALSE)] <- 100000
bray_curtis_no_keystone_sample_long <- melt(as.matrix(bray_curtis_no_keystone_sample), id = 'sample')
colnames(bray_curtis_no_keystone_sample_long) = c('sample1','sample2','bray_curtis')
bray_curtis_no_keystone_sample_long <- bray_curtis_no_keystone_sample_long[bray_curtis_no_keystone_sample_long$bray_curtis!=100000,]
write.table(bray_curtis_no_keystone_sample_long, "./keystone/pearson/ASV-bray_curtis-no_keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
##读取AS系统微生物群落的jaccard距离
jaccard <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_Jaccard.csv', header=TRUE, row.names=1)
jaccard_no_keystone_sample <- jaccard[colnames(no_keystone_sample),colnames(no_keystone_sample)]
write.csv(jaccard_no_keystone_sample, "./keystone/pearson/WWTP_ASV_jaccard-no_keystone_sample.csv", quote = FALSE)
jaccard_no_keystone_sample[!upper.tri(jaccard_no_keystone_sample, diag = FALSE)] <- 100000
jaccard_no_keystone_sample_long <- melt(as.matrix(jaccard_no_keystone_sample), id = 'sample')
colnames(jaccard_no_keystone_sample_long) = c('sample1','sample2','jaccard')
jaccard_no_keystone_sample_long <- jaccard_no_keystone_sample_long[jaccard_no_keystone_sample_long$jaccard!=100000,]
write.table(jaccard_no_keystone_sample_long, "./keystone/pearson/ASV-jaccard-no_keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
##读取AS系统微生物群落的unweighted_Unifrac距离
unweighted_Unifrac <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_unweight_uniFrac.csv', header=TRUE, row.names=1)
unweighted_Unifrac_no_keystone_sample <- unweighted_Unifrac[colnames(no_keystone_sample),colnames(no_keystone_sample)]
write.csv(unweighted_Unifrac_no_keystone_sample, "./keystone/pearson/WWTP_ASV_unweighted_Unifrac-no_keystone_sample.csv", quote = FALSE)
unweighted_Unifrac_no_keystone_sample[!upper.tri(unweighted_Unifrac_no_keystone_sample, diag = FALSE)] <- 100000
unweighted_Unifrac_no_keystone_sample_long <- melt(as.matrix(unweighted_Unifrac_no_keystone_sample), id = 'sample')
colnames(unweighted_Unifrac_no_keystone_sample_long) = c('sample1','sample2','unweighted_Unifrac')
unweighted_Unifrac_no_keystone_sample_long <- unweighted_Unifrac_no_keystone_sample_long[unweighted_Unifrac_no_keystone_sample_long$unweighted_Unifrac!=100000,]
write.table(unweighted_Unifrac_no_keystone_sample_long, "./keystone/pearson/ASV-unweighted_Unifrac-no_keystone_sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)


library(ggplot2)
library(vegan)
library(plyr)
group <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group.csv',header = T, row.names = 1)
##ASVs表差异
ASV_data <- read.csv('./rawdata/ASV_data_0.00001_0.01.csv', row.names = 1,header=TRUE)
ASV_data <- t(ASV_data)
##PCoA分析
dis1 <- vegdist(ASV_data, method = 'bray')
pcoa <- cmdscale(dis1, k = (nrow(ASV_data) - 1), eig = TRUE)
point <- data.frame(pcoa$point)
#坐标轴解释量（前两轴）
pcoa_eig <- (pcoa$eig)[1:2] / sum(pcoa$eig)
#提取样本点坐标（前两轴）
sample_site <- data.frame({pcoa$point})[1:2]
sample_site$sample <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1', 'PCoA2')
#为样本点坐标添加分组信息
sample_site <- merge(sample_site, group,by = 'sample')
# #将各分组类型转化为因子数据，方便作图识别，顺便进行因子排序。
sample_site$group <- factor(sample_site$group)

par(mar=c(3,5,2,2))
ggplot(sample_site, aes(PCoA1, PCoA2)) +
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        axis.text= element_text(size=12,color="black"),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        legend.text = element_text(size = 8, colour = "black")) + #去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  #geom_polygon(data = group_border, aes(fill = group), alpha = 0.5)+  #绘制多边形区域
  geom_point(aes(color = group), size = 1.5, alpha = 0.8)+ #可在这里修改点的透明度、大小
  stat_ellipse(geom="polygon", aes(fill=group), alpha=0.2, level = 0.95, show.legend = FALSE) +
  scale_shape_manual(values = c(15, 17, 16)) +  #可在这里修改点的形状
  scale_color_manual(values = c('#556B2F','#800000','#D87093','#0066ff')) + #可在这里修改点的颜色
  scale_fill_manual(values = c('#556B2F','#800000','#D87093','#0066ff'))  +#可在这里修改区块的颜色
  labs(x=paste("PCoA1 (", round(pcoa_eig[1]*100,2), "%)", sep = ''), y = paste("PCoA2 (", round(pcoa_eig[2]*100,2), "%)", sep = ''))
dev.off() 
ggsave("./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group-ASVs-PCoA.pdf", height=6,width=8)
#置换多元方差分析（PERMANOVA），vegan 包 adonis()，详情使用 ?adonis 查看帮助
dis1 <- vegdist(ASV_data, method = 'bray')
adonis2(dis1~group, group, permutations = 999) 
'''
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dis1 ~ group, data = group, permutations = 999)
           Df SumOfSqs      R2      F Pr(>F)    
group       1    13.80 0.03007 36.701  0.001 ***
Residual 1184   445.32 0.96993                  
Total    1185   459.12 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'''
##PCA分析
pca1 <- prcomp(ASV_data,center = TRUE,scale. = TRUE)
df1 <- pca1$x # 提取PC score
df1 <- as.data.frame(df1) # 注意：如果不转成数据框形式后续绘图时会报错
# 提取主成分的方差贡献率,生成坐标轴标题
summ1 <- summary(pca1)
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
df1$sample <- rownames(df1)
df2 <- merge(df1, group, by="sample")
# 绘制PCA得分图
p.pca1 <- ggplot(data = df2,aes(x = PC1,y = PC2))+
  stat_ellipse(geom="polygon", aes(fill=group), alpha=0.2, level = 0.95, show.legend = FALSE) +
  scale_shape_manual(values = c(15, 17, 16)) +  #可在这里修改点的形状
  scale_color_manual(values = c('#556B2F','#800000','#D87093','#0066ff')) + #可在这里修改点的颜色
  scale_fill_manual(values = c('#556B2F','#800000','#D87093','#0066ff'))  +#可在这里修改区块的颜色
  geom_point(aes(color = group), size = 1.5, alpha = 0.8)+ #可在这里修改点的透明度、大小
  labs(x = xlab1,y = ylab1,color = "Condition")+
  theme(panel.grid = element_line(color = 'gray', linetype = 2, size = 0.1), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        axis.text= element_text(size=12,color="black"),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        legend.text = element_text(size = 8, colour = "black")) + #去掉背景框
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.4)
ggsave(p.pca1,filename ="./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group-ASVs-PCA.pdf", height=6,width=8)

#多样性差异(PERManove)
alpha <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/rarefied_ASV_diversity-new.csv', row.names = 1,header=TRUE)
alpha$sample <- rownames(alpha)
alpha_group <- merge(alpha, group, by="sample")
##shannon
par(mar=c(3,5,2,2))
p<- ggplot(alpha_group, aes(x=group, y=shannon, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,8)
p + labs(y = "Shannon-Wiener index") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_shannon_compare-boxplot.pdf", width = 4, height = 3)
##evenness
par(mar=c(3,5,2,2))
p<- ggplot(alpha_group, aes(x=group, y=evenness, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "Pielou's evenness index") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_evenness_compare-boxplot.pdf", width = 4, height = 3)
##richness
par(mar=c(3,5,2,2))
p<- ggplot(alpha_group, aes(x=group, y=richness, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,2500)
p + labs(y = "Species richness") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_richness_compare-boxplot.pdf", width = 4, height = 3)
##PD
par(mar=c(3,5,2,2))
p<- ggplot(alpha_group, aes(x=group, y=PD, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,200)
p + labs(y = "Phylogenetic diversity") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_PD_compare-boxplot.pdf", width = 4, height = 3)
t.test(alpha_group[alpha_group$group=="Keystone_sample",]$shannon, alpha_group[alpha_group$group=="No_keystone_sample",]$shannon) ##p-value < 2.2e-16(5.606224  5.174574)
t.test(alpha_group[alpha_group$group=="Keystone_sample",]$evenness, alpha_group[alpha_group$group=="No_keystone_sample",]$evenness) ##p-value < 2.2e-16(0.8327265 0.8006012)
t.test(alpha_group[alpha_group$group=="Keystone_sample",]$richness, alpha_group[alpha_group$group=="No_keystone_sample",]$richness) ##p-value < 2.2e-16(858.9027  666.8353)
t.test(alpha_group[alpha_group$group=="Keystone_sample",]$PD, alpha_group[alpha_group$group=="No_keystone_sample",]$PD) ##p-value < 2.2e-16(79.08164  65.64136)

##群落结构变异性差异（Beta多样性）
##bray_curtis
keytsone_bray_curtis <- read.table("./keystone/pearson/ASV-bray_curtis-keystone_sample.txt", sep = '\t', header=TRUE)
keytsone_bray_curtis$group <-"Keystone_samples"
no_keytsone_bray_curtis <- read.table("./keystone/pearson/ASV-bray_curtis-no_keystone_sample.txt", sep = '\t', header=TRUE)
no_keytsone_bray_curtis$group <- "No_keystone_samples"
all_bray_curtis <- rbind(keytsone_bray_curtis, no_keytsone_bray_curtis)
par(mar=c(3,5,2,2))
p<- ggplot(all_bray_curtis, aes(x=group, y=bray_curtis, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "Bray curtis") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_bray_curtis_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_bray_curtis[all_bray_curtis$group=="Keystone_samples",]$bray_curtis, all_bray_curtis[all_bray_curtis$group=="No_keystone_samples",]$bray_curtis) ##p-value < 2.2e-16(0.8583213 0.8884534)

##jaccard
keytsone_jaccard <- read.table("./keystone/pearson/ASV-jaccard-keystone_sample.txt", sep = '\t', header=TRUE)
keytsone_jaccard$group <-"Keystone_samples"
no_keytsone_jaccard <- read.table("./keystone/pearson/ASV-jaccard-no_keystone_sample.txt", sep = '\t', header=TRUE)
no_keytsone_jaccard$group <- "No_keystone_samples"
all_jaccard <- rbind(keytsone_jaccard, no_keytsone_jaccard)
par(mar=c(3,5,2,2))
p<- ggplot(all_jaccard, aes(x=group, y=jaccard, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "Jaccard") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_jaccard_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_jaccard[all_jaccard$group=="Keystone_samples",]$jaccard, all_jaccard[all_jaccard$group=="No_keystone_samples",]$jaccard) ##p-value < 2.2e-16(0.9202201 0.9379143)

##unweighted_Unifrac
keytsone_unweighted_Unifrac <- read.table("./keystone/pearson/ASV-unweighted_Unifrac-keystone_sample.txt", sep = '\t', header=TRUE)
keytsone_unweighted_Unifrac$group <-"Keystone_samples"
no_keytsone_unweighted_Unifrac <- read.table("./keystone/pearson/ASV-unweighted_Unifrac-no_keystone_sample.txt", sep = '\t', header=TRUE)
no_keytsone_unweighted_Unifrac$group <- "No_keystone_samples"
all_unweighted_Unifrac <- rbind(keytsone_unweighted_Unifrac, no_keytsone_unweighted_Unifrac)
par(mar=c(3,5,2,2))
p<- ggplot(all_unweighted_Unifrac, aes(x=group, y=unweighted_Unifrac, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "Unweighted Unifrac") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_unweighted_Unifrac_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_unweighted_Unifrac[all_unweighted_Unifrac$group=="Keystone_samples",]$unweighted_Unifrac, all_unweighted_Unifrac[all_unweighted_Unifrac$group=="No_keystone_samples",]$unweighted_Unifrac) ##p-value < 2.2e-16(0.6653606 0.7107541)

##群落变异性差异（以weighted_unifrac距离为例）
keytsone_unifrac <- read.table("./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-weighted_Unifrac-keystone_sample.txt", sep = '\t', header=TRUE)
keytsone_unifrac$group <-"Keystone_samples"
no_keytsone_unifrac <- read.table("./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-weighted_Unifrac-no_keystone_sample.txt", sep = '\t', header=TRUE)
no_keytsone_unifrac$group <- "No_keystone_samples"
all_unifrac <- rbind(keytsone_unifrac, no_keytsone_unifrac)
samples_compare_group <-all_unifrac[,c("sample1","sample2","group")]
write.table(samples_compare_group, "./keystone/pearson/WWTP_samples_compare-group.txt", sep = '\t', quote = FALSE, row.names = FALSE)
Env_compare_all <-all_unifrac[,c("sample1","sample2","Env")]
write.table(Env_compare_all, "./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean_all-samples.txt", sep = '\t', quote = FALSE, row.names = FALSE)

par(mar=c(3,5,2,2))
p<- ggplot(all_unifrac, aes(x=group, y=weighted_Unifrac, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "Weighted Unifrac") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_weighted_unifrac_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_unifrac[all_unifrac$group=="Keystone_samples",]$weighted_Unifrac, all_unifrac[all_unifrac$group=="No_keystone_samples",]$weighted_Unifrac) ##p-value < 2.2e-16(0.3401960 0.3804729)

par(mar=c(3,5,2,2))
p<- ggplot(all_unifrac, aes(x=group, y=Env, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,30)
p + labs(y = "Environmental distance") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_env_distance_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_unifrac[all_unifrac$group=="Keystone_samples",]$Env, all_unifrac[all_unifrac$group=="No_keystone_samples",]$Env) ##p-value = 0.000596(8.979899  9.006808)

##环境因子校正的weighted_unifrac距离
all_unifrac$weighted_Unifrac_EnvNormed <- all_unifrac$weighted_Unifrac/all_unifrac$Env
all_unifrac_new <- all_unifrac[all_unifrac$weighted_Unifrac_EnvNormed != "Inf",]
par(mar=c(3,5,2,2))
p<- ggplot(all_unifrac_new, aes(x=group, y=weighted_Unifrac_EnvNormed, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,800)
p + labs(y = "Corrected weighted unifrac") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_weighted_unifrac_envNormed_compare-boxplot.pdf", width = 4, height = 3)
t.test(all_unifrac_new[all_unifrac_new$group=="Keystone_samples",]$weighted_Unifrac_EnvNormed, all_unifrac_new[all_unifrac_new$group=="No_keystone_samples",]$weighted_Unifrac_EnvNormed) ##p-value = 0.0003056(0.04408236 0.07764000)

##计算并比较校正的bray_curtis距离
Env_compare_all <- read.table("./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean_all-samples.txt", header = T, sep="\t")
Env_bray_curtis<- merge(Env_compare_all,all_bray_curtis, by=c('sample1','sample2'))
write.table(Env_bray_curtis, "./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-bray_curtis-all-sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
Env_bray_curtis$bray_curtis_EnvNormed <- Env_bray_curtis$bray_curtis/Env_bray_curtis$Env
Env_bray_curtis_new <- Env_bray_curtis[Env_bray_curtis$bray_curtis_EnvNormed != "Inf",]
par(mar=c(3,5,2,2))
p<- ggplot(Env_bray_curtis_new, aes(x=group, y=bray_curtis_EnvNormed, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1500)
p + labs(y = "Corrected bray curtis") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_bray_curtis_envNormed_compare-boxplot.pdf", width = 4, height = 3)
t.test(Env_bray_curtis_new[Env_bray_curtis_new$group=="Keystone_samples",]$bray_curtis_EnvNormed, Env_bray_curtis_new[Env_bray_curtis_new$group=="No_keystone_samples",]$bray_curtis_EnvNormed) ##p-value = 0.0008815(0.1123922 0.1767616)

##计算并比较校正的jaccard距离
Env_jaccard<- merge(Env_compare_all,all_jaccard, by=c('sample1','sample2'))
write.table(Env_jaccard, "./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-jaccard-all-sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
Env_jaccard$jaccard_EnvNormed <- Env_jaccard$jaccard/Env_jaccard$Env
Env_jaccard_new <- Env_jaccard[Env_jaccard$jaccard_EnvNormed != "Inf",]
par(mar=c(3,5,2,2))
p<- ggplot(Env_jaccard_new, aes(x=group, y=jaccard_EnvNormed, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1600)
p + labs(y = "Corrected jaccard") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_jaccard_envNormed_compare-boxplot.pdf", width = 4, height = 3)
t.test(Env_jaccard_new[Env_jaccard_new$group=="Keystone_samples",]$jaccard_EnvNormed, Env_jaccard_new[Env_jaccard_new$group=="No_keystone_samples",]$jaccard_EnvNormed) ##p-value = 0.001417(0.1240167 0.1940706)

##计算并比较校正的unweighted_unifrac距离
Env_compare_all <- read.table("./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean_all-samples.txt", header = T, sep="\t")
Env_unweighted_unifrac<- merge(Env_compare_all,all_unweighted_Unifrac, by=c('sample1','sample2'))
write.table(Env_unweighted_unifrac, "./keystone/pearson/WWTP_meta_predata_all_env-influx-standardize-euclidean-unweighted_unifrac-all-sample.txt", sep = '\t', quote = FALSE, row.names = FALSE)
Env_unweighted_unifrac$unweighted_unifrac_EnvNormed <- Env_unweighted_unifrac$unweighted_Unifrac/Env_unweighted_unifrac$Env
Env_unweighted_unifrac_new <- Env_unweighted_unifrac[Env_unweighted_unifrac$unweighted_unifrac_EnvNormed != "Inf",]
par(mar=c(3,5,2,2))
p<- ggplot(Env_unweighted_unifrac_new, aes(x=group, y=unweighted_unifrac_EnvNormed, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1500)
p + labs(y = "Corrected unweighted_unifrac") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_unweighted_unifrac_envNormed_compare-boxplot.pdf", width = 4, height = 3)
t.test(Env_unweighted_unifrac_new[Env_unweighted_unifrac_new$group=="Keystone_samples",]$unweighted_unifrac_EnvNormed, Env_unweighted_unifrac_new[Env_unweighted_unifrac_new$group=="No_keystone_samples",]$unweighted_unifrac_EnvNormed) ##p-value = 0.0007268(0.09179023 0.14518546)


#功能差异(PERManove)
functioner <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/qiime分析/ASV-new/OTU_analysis/ASV-MiDAS-function-all-new.csv', row.names = 1,header=TRUE)
functioner$sample <- rownames(functioner)
functioner_group <- merge(functioner, group, by="sample")
##Nitrifiers
par(mar=c(3,5,2,2))
p<- ggplot(functioner_group, aes(x=group, y=Nitrifiers, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,0.15)
p + labs(y = "MRA of Nitrifiers") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_Nitrifiers_compare-boxplot.pdf", width = 4, height = 3)
##Denitrifiers
par(mar=c(3,5,2,2))
p<- ggplot(functioner_group, aes(x=group, y=Denitrifiers, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "MRA of Denitrifiers") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_Denitrifiers_compare-boxplot.pdf", width = 4, height = 3)
##PAOs
par(mar=c(3,5,2,2))
p<- ggplot(functioner_group, aes(x=group, y=PAOs, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,0.25)
p + labs(y = "MRA of PAOs") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_PAOs_compare-boxplot.pdf", width = 4, height = 3)
##GAOs
par(mar=c(3,5,2,2))
p<- ggplot(functioner_group, aes(x=group, y=GAOs, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,0.2)
p + labs(y = "MRA of GAOs") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_GAOs_compare-boxplot.pdf", width = 4, height = 3)
##Filamentous
par(mar=c(3,5,2,2))
p<- ggplot(functioner_group, aes(x=group, y=Filamentous, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,0.8)
p + labs(y = "MRA of Filamentous") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_Filamentous_compare-boxplot.pdf", width = 4, height = 3)
t.test(functioner_group[functioner_group$group=="Keystone_sample",]$Nitrifiers, functioner_group[functioner_group$group=="No_keystone_sample",]$Nitrifiers) ##p-value = 3.841e-13(0.02934436 0.02062913)
t.test(functioner_group[functioner_group$group=="Keystone_sample",]$Denitrifiers, functioner_group[functioner_group$group=="No_keystone_sample",]$Denitrifiers) ##p-value = 3.368e-16(0.1694936 0.2039769)
t.test(functioner_group[functioner_group$group=="Keystone_sample",]$PAOs, functioner_group[functioner_group$group=="No_keystone_sample",]$PAOs) ##p-value = 0.0003765(0.03129674 0.03735164)
t.test(functioner_group[functioner_group$group=="Keystone_sample",]$GAOs, functioner_group[functioner_group$group=="No_keystone_sample",]$GAOs) ##p-value = 0.0001626(0.02223992 0.02764428)
t.test(functioner_group[functioner_group$group=="Keystone_sample",]$Filamentous, functioner_group[functioner_group$group=="No_keystone_sample",]$Filamentous) ##p-value = 1.582e-05(0.04212089 0.05296374)

#性能差异(PERManove)
performance <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_performance.csv', row.names = 1,header=TRUE)
performance$sample <- rownames(performance)
performance_group <- merge(performance, group, by="sample")
write.csv(performance_group, './keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group-meta_predata_all_performance.csv',quote = FALSE, row.names = FALSE)
performance_group <- read.csv("./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group-meta_predata_all_performance.csv", header = T, row.names = 1)
##BOD
par(mar=c(3,5,2,2))
p<- ggplot(performance_group, aes(x=group, y=BOD_remove_rate, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,6)
p + labs(y = "BOD removal rate") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"))
        #panel.background = element_rect(fill = "transparent"),
        #panel.border = element_rect(fill = "transparent", color = "transparent"),
        #axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_BOD_removal_rate_compare-boxplot-new.pdf", width = 3, height = 3.8)
##COD
par(mar=c(3,5,2,2))
p<- ggplot(performance_group, aes(x=group, y=COD_remove_rate, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(-0.1,20)
p + labs(y = "COD removal rate") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"))
#panel.background = element_rect(fill = "transparent"),
#panel.border = element_rect(fill = "transparent", color = "transparent"),
#axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_COD_removal_rate_compare-boxplot-new.pdf", width = 3, height = 3.8)
##NH4N
par(mar=c(3,5,2,2))
p<- ggplot(performance_group, aes(x=group, y=NH4N_remove_rate, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(-0.2,1)
p + labs(y = "NH4N removal rate") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"))
#panel.background = element_rect(fill = "transparent"),
#panel.border = element_rect(fill = "transparent", color = "transparent"),
#axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_NH4N_removal_rate_compare-boxplot-new.pdf", width = 3, height = 3.8)
##TN
par(mar=c(3,5,2,2))
p<- ggplot(performance_group, aes(x=group, y=TN_remove_rate, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(-0.1,1.5)
p + labs(y = "TN removal rate") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"))
#panel.background = element_rect(fill = "transparent"),
#panel.border = element_rect(fill = "transparent", color = "transparent"),
#axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_TN_removal_rate_compare-boxplot-new.pdf", width = 3, height = 3.8)
##TP
par(mar=c(3,5,2,2))
p<- ggplot(performance_group, aes(x=group, y=TP_remove_rate, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(-0.01,0.4)
p + labs(y = "TP removal rate") +
  theme_bw() +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"))
#panel.background = element_rect(fill = "transparent"),
#panel.border = element_rect(fill = "transparent", color = "transparent"),
#axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_TP_removal_rate_compare-boxplot-new.pdf", width = 3, height = 3.8)
t.test(performance_group[performance_group$group=="Keystone_sample",]$BOD_remove_rate, performance_group[performance_group$group=="No_keystone_sample",]$BOD_remove_rate) ##p-value = 8.525e-12(0.2043972 0.4701664)
t.test(performance_group[performance_group$group=="Keystone_sample",]$COD_remove_rate, performance_group[performance_group$group=="No_keystone_sample",]$COD_remove_rate) ##p-value = 2.885e-10(0.4358145 0.9618106)
t.test(performance_group[performance_group$group=="Keystone_sample",]$NH4N_remove_rate, performance_group[performance_group$group=="No_keystone_sample",]$NH4N_remove_rate) ##p-value = 5.289e-07(0.02649756 0.04411104)
t.test(performance_group[performance_group$group=="Keystone_sample",]$TN_remove_rate, performance_group[performance_group$group=="No_keystone_sample",]$TN_remove_rate) ##p-value = 0.001195(0.04367684 0.06106122)
t.test(performance_group[performance_group$group=="Keystone_sample",]$TP_remove_rate, performance_group[performance_group$group=="No_keystone_sample",]$TP_remove_rate) ##p-value = 0.003475(0.005860259 0.009207885)

##合并不同matter的性能一起画图
library(reshape2)
performance_group_long <- melt(performance_group,id.vars = c('group'), variable.name='Matter', value.name='Performance')
par(mar=c(3,5,2,2))
p<- ggplot(performance_group_long, aes(x=Matter, y=Performance, fill=group)) + 
  geom_boxplot(alpha=0.8) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(-0.2,16)
p + labs(y = "Performance(removal rate)") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_all-matters_removal_rate_compare-boxplot.pdf", width = 6, height = 3)


performance <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_NA_performance.csv', row.names = 1,header=TRUE)
performance$sample <- rownames(performance)
performance_group <- merge(performance, group, by="sample")
t.test(as.numeric(performance_group[performance_group$group=="Keystone_sample"&performance_group$BOD_removal!="#VALUE!",]$BOD_removal), as.numeric(performance_group[performance_group$group=="No_keystone_sample"&performance_group$BOD_removal!="#VALUE!",]$BOD_removal)) ##p-value = 2.517e-11(0.1617599 0.3192662)
t.test(as.numeric(performance_group[performance_group$group=="Keystone_sample"&performance_group$COD_removal!="#VALUE!",]$COD_removal), as.numeric(performance_group[performance_group$group=="No_keystone_sample"&performance_group$COD_removal!="#VALUE!",]$COD_removal)) ##p-value = 2.539e-06(0.3678888 0.7145210)
t.test(as.numeric(performance_group[performance_group$group=="Keystone_sample"&performance_group$NH4N_removal!="#VALUE!",]$NH4N_removal), as.numeric(performance_group[performance_group$group=="No_keystone_sample"&performance_group$NH4N_removal!="#VALUE!",]$NH4N_removal)) ##p-value = 0.0003487(0.02200588 0.03418323)
t.test(as.numeric(performance_group[performance_group$group=="Keystone_sample"&performance_group$TN_removal!="#VALUE!",]$TN_removal), as.numeric(performance_group[performance_group$group=="No_keystone_sample"&performance_group$TN_removal!="#VALUE!",]$TN_removal)) ##p-value = 5.291e-05(0.02499996 0.04638631)
t.test(as.numeric(performance_group[performance_group$group=="Keystone_sample"&performance_group$TP_removal!="#VALUE!",]$TP_removal), as.numeric(performance_group[performance_group$group=="No_keystone_sample"&performance_group$TP_removal!="#VALUE!",]$TP_removal)) ##p-value = 0.5145(0.005761176 0.006232365)

performance <- decostand(performance, 'standardize')   #可选 z-score 标准化，消除量纲差异
dis <- vegdist(performance, method = 'euclidean')      #获取欧式距离矩阵（这里以欧式距离为例，实际分析中视情况选择合适的距离测度）
adonis2(dis~group, group, permutations = 999)  #999 次置换获得 p 值（其实就是通过置换检验来实现）
'''
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dis ~ group, data = group, permutations = 999)
           Df SumOfSqs      R2    F Pr(>F)  
group       1     25.6 0.00432 5.14  0.016 *
Residual 1184   5899.4 0.99568              
Total    1185   5925.0 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'''


#环境差异(PERManove)
env <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env.csv', row.names = 1,header=TRUE)
env <- decostand(env, 'standardize')   #可选 z-score 标准化，消除量纲差异
dis <- vegdist(env, method = 'euclidean')      #获取欧式距离矩阵（这里以欧式距离为例，实际分析中视情况选择合适的距离测度）
adonis2(dis~group, group, permutations = 999)  #999 次置换获得 p 值（其实就是通过置换检验来实现）
'''
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dis ~ group, data = group, permutations = 999)
           Df SumOfSqs      R2      F Pr(>F)    
group       1     2472 0.04347 53.803  0.001 ***
Residual 1184    54408 0.95653                  
Total    1185    56880 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'''
env <- env[group$sample,]
env_group <- cbind(env, group)
para_list <- colnames(env)
t_test_result <- data.frame(para = character(),Keystone_mean=double(),No_keystone_mean=double(),p_value=double(),Sig=character())
for(i in 1:length(para_list)){
  para <- para_list[i]
  t_test <- t.test(env_group[env_group$group=="Keystone_sample",para], env_group[env_group$group=="No_keystone_sample",para])
  Keystone_mean <- t_test[["estimate"]][[1]]
  No_keystone_mean <- t_test[["estimate"]][[2]]
  p_value <- t_test[["p.value"]]
  if (p_value <= 0.001) Sig <- '***' else if (p_value <= 0.01) Sig <- '**' else if (p_value <= 0.05) Sig <- '*' else if (p_value > 0.05) Sig <- 'n.s'
  para_result <- data.frame(para, Keystone_mean, No_keystone_mean, p_value, Sig)
  t_test_result <- rbind(t_test_result, para_result)
}
write.csv(t_test_result, "./keystone/pearson/Keystone_samples_group-diffenv-t_test-sig.csv", quote = FALSE, row.names = FALSE)

t_test_result <- read.csv("./keystone/pearson/Keystone_samples_group-diffenv-t_test-sig.csv", header = T, row.names = 1)
t_test_result_sig <- t_test_result[t_test_result$Sig != "n.s",]
env_group_performance <- env_group[,c("InfBOD","EffBOD","InfCOD","EffCOD","InfNH4N","EffNH4N","InfTN","EffTN","InfTP","EffTP","group")]
env_group_long <- melt(env_group_performance,id.vars = c('group'), variable.name='Env', value.name='para')
#加载包
library(ggplot2)#绘图包
library(ggpubr)#基于ggplot2的可视化包，主要用于绘制符合出版要求的图形
library(ggsignif)#用于P值计算和显著性标记
library(tidyverse)#数据预处理
par(mar=c(3,5,2,2))
p <- ggplot(env_group_long,aes(x=Env,y=para))+#指定数据
  geom_boxplot(aes(fill=group), outlier.colour="white", width=0.5,size=0.8) +##异常点去除
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1000)
p + labs(y = "Concentration of pollutants(mg/l)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_Concentration_pollutants_compare-boxplot.pdf", width = 9, height = 6)

par(mar=c(3,5,2,2))
p <- ggplot(env_group,aes(x=group,y=InfR))+#指定数据
  geom_boxplot(aes(fill=group), outlier.colour="white", width=0.5,size=0.8) +##异常点去除
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,500000)
p + labs(y = "Actual influx rate(m3/d)") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"),
        legend.position="none")
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_InfR_compare-boxplot.pdf", width = 4, height = 3)

par(mar=c(3,5,2,2))
p <- ggplot(env_group,aes(x=group,y=Vol))+#指定数据
  geom_boxplot(aes(fill=group), outlier.colour="white", width=0.5,size=0.8) +##异常点去除
  scale_fill_manual(values = c('#556B2F','#800000')) + 
  ylim(0,100000)
p + labs(y = "Volume of aeration tanks(m3)") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"),
        legend.position="none")
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_Vol_compare-boxplot.pdf", width = 4, height = 3)

par(mar=c(3,5,2,2))
p <- ggplot(env_group,aes(x=group,y=MLSS))+#指定数据
  geom_boxplot(aes(fill=group), outlier.colour="white", width=0.5,size=0.8) +##异常点去除
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,6000)
p + labs(y = "MLSS(mg/l)") +
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"),
        legend.position="none")
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_MLSS_compare-boxplot.pdf", width = 4, height = 3)

##网络结构差异
SubNet_all <- read.csv('./igraph/pearson/RMT-Subnet_0.00001_0.01_pearson_global_topology.csv', row.names = 1,header=TRUE)
SubNet_stat <- colnames(SubNet_all)
SubNet_all$sample <- rownames(SubNet_all)
SubNet_all_group <- merge(SubNet_all, group,by="sample")
##网络参数比较
t_test_result <- data.frame(para = character(),Keystone_mean=double(),No_keystone_mean=double(),p_value=double(),Sig=character())
for(i in 1:length(SubNet_stat)){
  para <- SubNet_stat[i]
  t_test <- t.test(SubNet_all_group[SubNet_all_group$group=="Keystone_sample",para], SubNet_all_group[SubNet_all_group$group=="No_keystone_sample",para])
  Keystone_mean <- t_test[["estimate"]][[1]]
  No_keystone_mean <- t_test[["estimate"]][[2]]
  p_value <- t_test[["p.value"]]
  if (p_value <= 0.001) Sig <- '***' else if (p_value <= 0.01) Sig <- '**' else if (p_value <= 0.05) Sig <- '*' else if (p_value > 0.05) Sig <- 'n.s'
  para_result <- data.frame(para, Keystone_mean, No_keystone_mean, p_value, Sig)
  t_test_result <- rbind(t_test_result, para_result)
}
write.csv(t_test_result, "./keystone/pearson/Keystone_samples_group-diffSubNet_stat-t_test-sig.csv", quote = FALSE, row.names = FALSE)

##比较包含/不包含keystone类群，环境因子与群落结构的相关性
keystone_mantel <- read.csv("./keystone/pearson/ASV-data-all-environment-influx-all-ASVs-keystone_sample-mantel-test.csv",header = T)
no_keystone_mantel <- read.csv("./keystone/pearson/ASV-data-all-environment-influx-all-ASVs-no_keystone_sample-mantel-test.csv",header = T)
#matter <- c("F.M","MLSS","DO","Ph","MIT","Con","SVI","ReInfBOD","AtInfBOD","ReInfCOD","AtInfCOD","AtInfNH4N","AtInfTN","AtInfTP","NH4N","NO3N")
keystone_all_env <- keystone_mantel[!keystone_mantel$env%in%c("Lat","Lon"),]
no_keystone_all_env <- no_keystone_mantel[!no_keystone_mantel$env%in%c("Lat","Lon"),]
keystone_all_env[keystone_all_env$p.value>0.05,]$r <- 0
no_keystone_all_env[no_keystone_all_env$p.value>0.05,]$r <- 0
keystone_all_env <- keystone_all_env[,c("env","r")]
colnames(keystone_all_env) <- c("env","Keystone_samples")
no_keystone_all_env <- no_keystone_all_env[,c("env","r")]
colnames(no_keystone_all_env) <- c("env","No_keystone_samples")
all_group_env <- merge(keystone_all_env, no_keystone_all_env, by="env")
write.csv(all_group_env, "./keystone/pearson/ASV-data-all-environment-influx-all-ASVs-keystone_sample-groups-mantel-test.csv",quote = FALSE, row.names = FALSE)
##全部环境因子热图
all_group_env <- read.csv("./keystone/pearson/ASV-data-all-environment-influx-all-ASVs-keystone_sample-groups-mantel-test.csv", header = T, row.names = 1)
##热图
library(pheatmap)
bk <- c(seq(0,0.4,by=0.001))
p <- pheatmap(as.matrix(all_group_env),cluster_row = TRUE,cluster_col = FALSE, show_colnames= TRUE, angle_col = 0,
              fontsize=5,clustering_method = "ward.D2",color = c(colorRampPalette(colors = c("white","darkred"))(length(bk))),
              legend_breaks=seq(0,0.4,by=0.1),breaks=bk,filename="./keystone/pearson/ASV-data-all-environment-influx-all-ASVs-keystone_sample-groups-mantel-test-heatmap.pdf",width = 8,height = 6)


##全部Influx环境因子相关图
all_group_env <- read.csv("./keystone/pearson/ASV-data-all-environment-influx-all-ASVs-keystone_sample-groups-mantel-test.csv", header = T, row.names = 1)
Influx = c("IndConInf", "IndPer","InfBOD","InfCOD", "InfNH4N", "InfTN", "InfTP")##influx
all_group_influx <- all_group_env[Influx,]
library(pheatmap)
bk <- c(seq(0,0.3,by=0.001))
p <- pheatmap(as.matrix(all_group_influx),cluster_row = TRUE,cluster_col = FALSE, show_colnames= TRUE, angle_col = 0,
              fontsize=5,clustering_method = "ward.D2",color = c(colorRampPalette(colors = c("white","darkred"))(length(bk))),
              legend_breaks=seq(0,0.3,by=0.1),breaks=bk,filename="./keystone/pearson/Influx-all-ASVs-keystone_sample-groups-mantel-test-heatmap.pdf",width = 4,height = 3)
##mantel test 组合图
env <- read.csv("E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all-influx.csv", header = T, row.names = 1)
env_influx <- env[,Influx]
all_group_influx$para <- rownames(all_group_influx)
all_group_influx_long <- melt(all_group_influx,id.vars = c('para'), variable.name='group', value.name='r')
mantel<-mutate(all_group_influx_long,r = cut(abs(r), breaks = c(-Inf, 0.1, 0.2, Inf),
                                  labels = c("<0.1","0.1-0.2", ">=0.2"),
                                  right = FALSE))
mantel_adj <- mantel[,c("group","para","r")]
#绘制组合图
par(mar=c(3,5,2,2))
p1<-quickcor(env_influx, type = "upper") + 
  geom_square(inherit.aes = TRUE)+ 
  anno_link(mantel_adj, mapping = aes(size = r), colour ="#FF9664" ,alpha = 0.6)+
  scale_size_manual(values = c(0.5, 1.5, 2.5),limits=c("<0.1","0.1-0.2",">=0.2"))+   #线段size
  scale_fill_gradient2(midpoint = 0, limits = c(-1, 1), # 数据上下限
                       breaks = c(-1.0, -0.5, 0.0, 0.5, 1.0), # 分段点
                       low = "darkblue", mid = "white",high = "darkred", space = "Lab" ) #热图颜色
#geom_diag_label(geom = "text", remove.axis = TRUE) ##这个句决定env是否以横纵坐标展示
#定义图例标题
p1+theme(legend.position="right", legend.title = element_text(size = 10), legend.text = element_text(size = 10),legend.key.width = unit(1,"cm"))+
  guides(size=guide_legend(title="Mantel's r",override.aes=list(colour="#FF9664"),order=1),
         fill=guide_colorbar(title="Pairwise correlations",barheight = 5,  order=2))
dev.off()
ggsave("./keystone/pearson/Influx-all-ASVs-keystone_sample-groups-mantel-test.pdf", width = 6, height = 4)

