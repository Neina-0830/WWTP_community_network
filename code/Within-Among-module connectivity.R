library(igraph)

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

adj.matrix <- read.csv('./igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_adj.csv', row.names = 1, check.names = FALSE)
node.topology <- read.csv('./igraph/pearson/RMT-Net_0.00001_0.01_node_topology_pearson.csv', row.names = 1, check.names = FALSE)
#两个文件的节点顺序要一致
node.topology <- node.topology[rownames(adj.matrix), ]
#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(adj.matrix, node.topology)
head(zi_pi)
write.csv(zi_pi, './igraph/pearson/RMT-Net_0.00001_0.01_node_topology_pearson_zi_pi.csv',quote=FALSE)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)
zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'
write.csv(zi_pi, './igraph/pearson/RMT-Net_0.00001_0.01_node_topology_pearson_zi_pi-type.csv',quote=FALSE, row.names = FALSE)

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
ggsave("./igraph/pearson/RMT-ASV_level_0.00001_0.01_network-Hubs.pdf", width = 5, height = 3.5)

##从zi_pi中定义不同类型的数据
Peripherals=zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities <= 0.62),]
Connectors=zi_pi[which(zi_pi$within_module_connectivities <= 2.5 & zi_pi$among_module_connectivities > 0.62),] 
Module_hubs=zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities <= 0.62),]
Network_hubs=zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),]
#导出数据
write.csv(Peripherals, './igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_Peripherals.csv',quote = FALSE,row.names = FALSE)
write.csv(Connectors, './igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_Connectors.csv',quote = FALSE,row.names = FALSE)
write.csv(Module_hubs, './igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_Module_hubs.csv',quote = FALSE,row.names = FALSE)
write.csv(Network_hubs, './igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_Network_hubs.csv',quote = FALSE,row.names = FALSE)
