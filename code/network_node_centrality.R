##“中心性(Centrality)"来判断网络中节点重要性或影响力。
library(igraph)

mynet.adj<- read.csv('./igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_adj.csv', row.names = 1, check.names = FALSE)
mynet <- graph_from_adjacency_matrix(as.matrix(mynet.adj), mode = 'undirected', weighted = NULL, diag = FALSE)
node_all <- read.csv("./Gephi/pearson/RMT-ASV_0.00001_0.01_node_all.csv",header=T)
##1、Degree centrality
node_degree_centrality <- node_all[,c("Label", "v_degree")]
node <- nrow(node_all)
node_degree_centrality$v_degree <- node_degree_centrality$v_degree/(node-1) ##标准化度中心性Dnorm=D/(n-1)
colnames(node_degree_centrality) <- c("Label","Degree_centrality")
write.csv(node_degree_centrality,"./igraph/pearson/RMT-ASV_0.00001_0.01_node_degree_centrality.csv",quote = FALSE, row.names = FALSE)

##2、Betweeness centrality
node_betweenness_centrality <- as.data.frame(betweenness(mynet, normalized=TRUE)) ##标准化介数中心性Bnorm=2*B/((n-1)*(n-2)).
node_betweenness_centrality$Label <- rownames(node_betweenness_centrality)
colnames(node_betweenness_centrality) <- c("Betweenness_centrality","Label")
write.csv(node_betweenness_centrality,"./igraph/pearson/RMT-ASV_0.00001_0.01_node_betweenness_centrality.csv",quote = FALSE, row.names = FALSE)

##3、Closeness centrality
node_closeness_centrality <- node_all[,c("Label", "norm_closeness_w")] ##标准化紧密中心性Cnorm=C/(n-1)
colnames(node_closeness_centrality) <- c("Label","Closeness_centrality")
write.csv(node_closeness_centrality,"./igraph/pearson/RMT-ASV_0.00001_0.01_node_closeness_centrality.csv",quote = FALSE, row.names = FALSE)

##4、Eigenvector centrality
node_eigenvector_centrality <- as.data.frame(eigen_centrality(mynet)$vector) ##标准化特征向量中心性，Enorm=E/max(E)
node_eigenvector_centrality$Label <- rownames(node_eigenvector_centrality)
colnames(node_eigenvector_centrality) <- c("Eigenvector_centrality","Label")
write.csv(node_eigenvector_centrality,"./igraph/pearson/RMT-ASV_0.00001_0.01_node_eigenvector_centrality.csv",quote = FALSE, row.names = FALSE)

##Combine centrality
node_centrality <- merge(node_degree_centrality, node_betweenness_centrality, by="Label")
node_centrality <- merge(node_centrality, node_closeness_centrality, by="Label")
node_centrality <- merge(node_centrality, node_eigenvector_centrality, by="Label")
write.csv(node_centrality,"./igraph/pearson/RMT-ASV_0.00001_0.01_node_all_centrality.csv",quote = FALSE, row.names = FALSE)
