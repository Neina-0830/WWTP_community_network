library(RMThreshold)
library(igraph)
library(psych)
#install.packages("tnet")
library(tnet)  ##计算不连通图中closeness的方法
##comun: A community table with samples as rows and taxa as columns. 
comun = read.csv("../data/ASV_data_0.00001_0.01.csv",header=T,row.names=1) ##6752个ASV
comun = t(comun)
comun[comun==0] <- 0.01  ##0丰度的转化（这里矩阵反映绝对序列数，最小非0值是1）
log_comun <- log(comun)
set.seed(111)
# 取相关性矩阵R值
# comun行为sample，列为OTUs
occor.r = cor(log_comun,method = "pearson")  ##考虑正相关和负相关
occor.r <- as.matrix(occor.r)
diag(occor.r) <- 0
res <- rm.get.threshold(occor.r, nr.thresholds = 26, interval = c(0.3, 0.8), discard.zeros = TRUE, unfold.method = "spline")
thre <- 0.44 
cleaned.matrix <- rm.denoise.mat(occor.r, threshold = thre) #denoise adjacency matrix by threshold
cleaned.matrix <- rm.discard.zeros(cleaned.matrix) #delete all-zero rows
mynet = graph_from_adjacency_matrix(cleaned.matrix ,mode="undirected",weighted=TRUE,diag=FALSE)
node <- length(V(mynet))###Number of nodes
edge <- length(E(mynet))###Number of edges
taxonomy = read.csv("../data/ASV-MiDAS-MRA-function-taxonomy-all-new.csv",header=T,row.names=1)  ##Taxonomy与data的顺序问题
##按照V(mynet)$name筛选且排序rownames(taxonomy)
loc = match(V(mynet)$name, rownames(taxonomy))
taxonomy <- taxonomy[loc,]
##set vertices size and color
V(mynet)$Kindom <- taxonomy$Kindom
V(mynet)$Phylum <- taxonomy$Phylum
V(mynet)$Class <- taxonomy$Class
V(mynet)$Order <- taxonomy$Order
V(mynet)$Family <- taxonomy$Family
V(mynet)$Genus <- taxonomy$Genus
V(mynet)$Species <- taxonomy$Species
V(mynet)$Partition <- taxonomy$Partition
V(mynet)$Group <- taxonomy$Group
V(mynet)$Taxa <- taxonomy$Taxa
V(mynet)$Function <- taxonomy$Function
V(mynet)$Process <- taxonomy$Process
V(mynet)$Frequency <- taxonomy$frequency
V(mynet)$Niche <- taxonomy$niche_width
V(mynet)$MRA <- taxonomy$MRA
n <- length(V(mynet))
mycolors <- colorRampPalette(colors())(n)
V(mynet)$color <- mycolors
E(mynet)$correlation<-E(mynet)$weight
E(mynet)$weight<-abs(E(mynet)$weight)
E(mynet)$Relation <- "Pos"
E(mynet)[E(mynet)$correlation<0]$Relation <- "Neg" ##90
set.seed(111)
V(mynet)$modularity <- membership(cluster_louvain(mynet))
V(mynet)$degree <- degree(mynet)####Degree
V(mynet)$betweenness <- betweenness(mynet)####Betweenness centrality
V(mynet)$closeness <- closeness(mynet)####Closeness centrality

##remove isolated nodes
bad.vs = V(mynet)[degree(mynet) == 0]
mynet = delete.vertices(mynet, bad.vs)
matrix_rep <- as.matrix(get.adjacency(mynet)) #this gives you the adjacency
node <- length(V(mynet))###Number of nodes
edge <- length(E(mynet))###Number of edges
write.csv(matrix_rep, "./igraph/RMT-ASV_level_0.00001_0.01_network_pearson_adj.csv",quote = FALSE)
write_graph(mynet, "./Gephi/RMT_ASV_0.00001_0.01-Network-pearson.graphml","graphml")

###Caculation of node-level topological features
modularity.class <- V(mynet)$modularity
node.degree <- degree(mynet)####Degree
betweenness.centrality <- betweenness(mynet)####Betweenness centrality
closeness.centrality <- closeness(mynet)####Closeness centrality
node.topology <- data.frame(node.degree, betweenness.centrality, closeness.centrality,modularity.class)
write.csv(node.topology, file="./igraph/RMT-Net_0.00001_0.01_node_topology_pearson.csv",quote = FALSE)

###Caculation of network-level topological features
node <- length(V(mynet))###Number of nodes  5503
edge <- length(E(mynet))###Number of edges 86753
ave_degree <- mean(degree(mynet))###Average degree 
clust_coeff <- transitivity(mynet)###Clustering coefficient
ave_path_len <- average.path.length(mynet)###Average path length
diameter <- diameter(mynet)###Network diameter
density <- graph.density(mynet)###Graph density
modular <- modularity(cluster_louvain(mynet))  ##Graph modularity
betweenness_centr <- centralization.betweenness(mynet)$centralization
global.topology <- data.frame(node,edge,ave_degree,clust_coeff,ave_path_len,diameter,density,modular,betweenness_centr)
write.csv(global.topology, file="./igraph/RMT-Net_0.00001_0.01_global_topology_pearson.csv",quote = FALSE)

##Random network
#Erdös−Réyni random network
library(igraph)
node <- length(V(mynet))###Number of nodes
edge <- length(E(mynet))###Number of edges 
ave_degree <- 0
clust_coeff <- 0
ave_path_len <- 0
diameter <- 0
density <- 0
modular <- 0
betweenness_centr <- 0
for (i in 1:100) {
  rand_mynet <- erdos.renyi.game(node,edge,'gnm') ##点，边
  rand_mynet <- simplify(rand_mynet)
  ave_degree[i] <- mean(degree(rand_mynet))###Average degree 
  clust_coeff[i] <- transitivity(rand_mynet)###Clustering coefficient
  ave_path_len[i] <- average.path.length(rand_mynet)###Average path length
  diameter[i] <- diameter(rand_mynet)###mynetwork diameter
  density[i] <- graph.density(rand_mynet)###Graph density
  modular[i] <- modularity(cluster_louvain(rand_mynet))  ##Graph modularity
  betweenness_centr[i] <- centralization.betweenness(rand_mynet)$centralization
} 
mean <- c(mean(node),mean(edge),mean(ave_degree),mean(clust_coeff),mean(ave_path_len),mean(diameter),mean(density),mean(modular),mean(betweenness_centr))
sd <- c(sd(node),sd(edge),sd(ave_degree),sd(clust_coeff),sd(ave_path_len),sd(diameter),sd(density),sd(modular),sd(betweenness_centr))#计算标准偏差并转换为向量
para<-rbind(mean,sd)
colnames(para)<-c("node","edge","ave_degree","clust_coeff","ave_path_len","diameter","density","modular","betweenness_centr")
write.csv(para,"./igraph/RMT-ASV_0.00001_0.01_random_global_topology_pearson.csv",quote = FALSE)
