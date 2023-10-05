library(igraph)

mynet.adj<- read.csv('./igraph/pearson/RMT-ASV_level_0.00001_0.01_network_pearson_adj.csv', row.names = 1, check.names = FALSE)
mynet <- graph_from_adjacency_matrix(as.matrix(mynet.adj), mode = 'undirected', weighted = NULL, diag = FALSE)
# 做幂律分布图Ploting node degreee distribution in a log-log plot
#计算节点度（Degree）以及获得其分布
#计算节点度
V(mynet)$degree <-degree(mynet) 
#度分布统计
degree_dist <- table(V(mynet)$degree)
degree_num <- as.numeric(names(degree_dist))
degree_count <- as.numeric(degree_dist)
dat <- data.frame(degree = degree_num, count = degree_count)
head(dat)
#查看度分布
#可观察到微生物相关网络通常服从幂律分布
par(mfrow = c(1, 3))
hist(V(mynet)$degree, xlab = 'Degree', ylab = 'Frequency', 
     main = 'Degree distribution')
plot(degree_count, degree_num, xlab = 'Degree', ylab = 'Count', 
     main = 'Degree distribution')
plot(degree_count, degree_num, log = 'xy', xlab = 'Log-degree', 
     ylab = 'Log-count', main = 'Log-log degree distribution')

####幂律分布
#拟合，a 和 b 的初始值手动指定，数据不同需要多加尝试
mod <- nls(count ~ a*degree^b, data = dat, start = list(a = 2, b = 0.5))
summary(mod)
#提取关键值
a <- round(coef(mod)[1], 3)
b <- round(coef(mod)[2], 3)
a; b
#使用构建好的方程，通过预测变量拟合响应变量，再根据和观测值的差异，获得 R2
#SSre：the sum of the squares of the distances of the points from the fit
#SStot：the sum of the squares of the distances of the points from a horizontal line through the mean of all Y values
fit <- fitted(mod)
SSre <- sum((dat$count-fit)^2)
SStot <- sum((dat$count-mean(dat$count))^2)
R2 <- round(1 - SSre/SStot, 3)
R2
#p 值可以根据置换检验的原理获得
#将 count 的值分别随机置换 N 次（例如 999 次），通过随机置换数据后数据获取 R2（R2'）
#比较随机置换后值的 R2' 大于观测值的 R2 的频率，即为 p 值
p_num <- 1
dat_rand <- dat
for (i in 1:999) {
  dat_rand$count <- sample(dat_rand$count)
  SSre_rand <- sum((dat_rand$count-fit)^2)
  SStot_rand <- sum((dat_rand$count-mean(dat_rand$count))^2)
  R2_rand <- 1 - SSre_rand/SStot_rand
  if (R2_rand > R2) p_num <- p_num + 1
}
p_value <- p_num / (999+1)
p_value

#作图展示
library(ggplot2)
p <- ggplot(dat, aes(x = degree, y = count)) +
  geom_point(color = 'red') +
  stat_smooth(method = 'nls', formula = y ~ a*x^b, method.args = list(start = list(a = 2, b = 0.5)), se = FALSE) +
  labs(x = 'Degree', y = 'Count')
#添加公式拟合的注释
label <- data.frame(
  formula = sprintf('italic(Y) == %.3f*italic(X)^%.3f', a, b),
  R2 = sprintf('italic(R^2) == %.3f', R2),
  p_value = sprintf('italic(P) < %.3f', p_value)
)
p1=p + geom_text(x = 25, y = 650, aes(label = formula), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 25, y = 600, aes(label = R2), data = label, parse = TRUE, hjust = 0) +
  geom_text(x = 25, y = 550, aes(label = p_value), data = label, parse = TRUE, hjust = 0)
p1
#去除背景和网络
p2=p1+theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2


##随机网络
node <- length(V(mynet))###Number of nodes
edge <- length(E(mynet))###Number of edges 
all_dat <- data.frame(degree=integer(),count=integer(),group=character(),stringsAsFactors=FALSE)
for(i in 1:100){
  randnet <- erdos.renyi.game(node,edge,'gnm') ##点，边
  #计算节点度
  V(randnet)$degree <-degree(randnet) 
  #度分布统计
  degree_dist <- table(V(randnet)$degree)
  degree_num <- as.numeric(names(degree_dist))
  degree_count <- as.numeric(degree_dist)
  dat <- data.frame(degree = degree_num, count = degree_count)
  dat$group <- paste("Network",i,sep = '')
  all_dat <- rbind(all_dat, dat)
}

head(all_dat)

#作图展示
library(ggplot2)
ggplot(all_dat, aes(x = degree, y = count),groups=group) +
  geom_point(aes(color = group),size=0.5,alpha=0.6) +
  labs(x = 'Degree', y = 'Number of nodes') +
  ylim(0,500) +
  xlim(0,60) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +
  theme(legend.position="none")
dev.off()
ggsave("./igraph/pearson/Random-poisson-fit.pdf",height=3,width=4)
