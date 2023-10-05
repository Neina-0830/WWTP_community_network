setwd("/media/wxl/Run1/lxn/Network/ASV-new/")
library(igraph)
library(ggplot2)
#定义函数，计算网络中节点的重要性，并根据重要性依次移出节点后，计算网络的自然连通度
nc <- function(g) {
  natcon <- function(g) {
    N   <- vcount(g)
    adj <- get.adjacency(g)
    evals <- eigen(adj)$value
    nc <- log(mean(exp(evals)))
    nc / (N - log(N))
  }
  nc.attack <- function(g) {
    hubord <- order(rank(betweenness(g)), rank(degree(g)), decreasing=TRUE)
    list <- seq(from=1, to=vcount(g), by=100)
    sapply(list, function(i) {
      ind <- hubord[1:i]
      tmp <- delete_vertices(g, V(g)$name[ind])
      natcon(tmp)
    })
  }
  nc <- nc.attack(g)
  nc
}

##AS network vs. random network
g <- read.csv('./igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_adj.csv', row.names = 1, check.names = FALSE)
#转换网络数据（邻接矩阵 -> R 的网络数据结构）
g <- graph_from_adjacency_matrix(as.matrix(g), mode = 'undirected', weighted = NULL, diag = FALSE)
g_nc <- nc(g)
#构建随机网络，并计算自然连通度
g_rand <- erdos.renyi.game(n = vcount(g), p = ecount(g), type = 'gnm', directed = FALSE)
g_rand_nc <- nc(g_rand)

dat <- data.frame(
  network = c(rep('microbial network', length(g_nc)), rep('random network', length(g_rand_nc))), 
  'Proportion of removes nodes' = c((1:length(g_nc))/length(g_nc), (1:length(g_rand_nc))/length(g_rand_nc)),
  'Natural Connectivity' = c(g_nc, g_rand_nc), check.names = FALSE
)
dat <- subset(dat, `Proportion of removes nodes` <= 0.8)
write.csv(dat,"./igraph/RMT-order-removal-compare-nc.csv",quote = FALSE)

dat <- read.csv("./igraph/pearson/RMT-order-removal-compare-nc.csv",header=T, row.names = 1)

ggplot(dat, aes(Proportion.of.removes.nodes, Natural_Connectivity_adj, color = network)) +
  geom_point() +
  theme_bw() +
  labs(x = 'Proportion of removes nodes', y = 'Natural connectivity') +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
      axis.title = element_text(face = "bold", size = 10, color = "black"),
      panel.background = element_rect(fill = "transparent"),
      panel.border = element_rect(fill = "transparent", color = "transparent"),
      axis.line = element_line(color = "black"))
dev.off()
ggsave("./igraph/pearson/RMT-order-removal-compare-nc.pdf", width = 5.5, height = 3)

##keystone_sample vs. no_keystone_sample
dat <- read.csv("./keystone/pearson/RMT-keystone-no-keystone_sample-order-removal-compare-nc.csv",header=T, row.names = 1)

ggplot(dat, aes(Proportion.of.removes.nodes, Natural_Connectivity_adj, color = network)) +
  geom_point() +
  ylim(0,150)+
  theme_bw() +
  scale_color_manual(values=c('#556B2F','#800000'))+
  labs(x = 'Proportion of removes nodes', y = 'Natural connectivity') +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-keystone-no-keystone_sample-order-removal-compare-nc.pdf", width = 5.5, height = 3)
##两子网络自然连通度曲线的F检验
##方差是否有显著差异
var.test(dat[dat$network=="no_keystone_network",]$Natural_Connectivity_adj, dat[dat$network=="keystone_network",]$Natural_Connectivity_adj)
'''
F = 0.24603, num df = 47, denom df = 48, p-value = 3.862e-06
alternative hypothesis: true ratio of variances is not equal to 1
95 percent confidence interval:
 0.1384708 0.4379111
sample estimates:
ratio of variances 
         0.2460308 
'''
##F检验显示，两个样本的方差有显著性差异，因此应该使用非等方差的t检验：
t.test(dat[dat$network=="no_keystone_network",]$Natural_Connectivity_adj, dat[dat$network=="keystone_network",]$Natural_Connectivity_adj, var.equal = F)## p-value = 0.006909（25.28948  47.94620）

