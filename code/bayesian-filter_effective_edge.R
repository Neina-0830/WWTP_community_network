library(igraph)

# 设置工作目录（请修改为您的目录路径）
setwd("E:/Desktop/WWTP-network_analysis/MiDAS/ASV/")

# 读取数据
Env <- read.csv("E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env.csv", header=TRUE, row.names=1)
Env_names <- colnames(Env)

# 定义分类
class1 <- c("Lat", "Lon")  # climate
class2 <- c("MAT", "AMMaxT", "AMMinT", "MMT", "SMT", "AP", "SMP")  # climate
class3 <- c("GDP", "Pop")  # climate
class4 <- c("BY", "DC", "Vol", "InfR", "HRT", "AtHRT", "SRT", "Nitri", "Denitri")  # control
class5 <- c("IndConInf", "IndPer", "InfBOD", "InfCOD", "InfNH4N", "InfTN", "InfTP")  # influx
class6 <- c("F.M", "MLSS", "DO", "Ph", "MIT", "Con", "SVI", "ReInfBOD", "AtInfBOD", "ReInfCOD", "AtInfCOD", "AtInfNH4N", "AtInfTN", "AtInfTP", "NH4N", "NO3N", "group")  # matter
class7 <- c("EffBOD", "EffCOD", "EffNH4N", "EffTN", "EffTP")

# 读取边数据
arcs_data <- read.csv("./keystone/pearson/Keystone_samples_group_meta_normed-all-h2pc-arcs.csv", header=TRUE)

# 过滤 Lat/Lon 只能为根节点
arcs_data <- arcs_data[!(arcs_data$to %in% class1), ]

# 进行层级依赖删除
arcs_data <- arcs_data[!(arcs_data$to %in% class2 & arcs_data$from %in% c(class3, class4, class5, class6, class7)), ]
arcs_data <- arcs_data[!(arcs_data$to %in% class3 & arcs_data$from %in% c(class4, class5, class6, class7)), ]
arcs_data <- arcs_data[!(arcs_data$to %in% class4 & arcs_data$from %in% c(class5, class6, class7)), ]
arcs_data <- arcs_data[!(arcs_data$to %in% class5 & arcs_data$from %in% c(class6, class7)), ]
arcs_data <- arcs_data[!(arcs_data$from %in% class7), ]

# 生成输出文件名
if("group" %in% arcs_data$to){
  newname <- "./keystone/pearson/Keystone_samples_group_meta_predata_all-h2pc-arcs-effective.csv"
} else {
  newname <- "./keystone/pearson/Keystone_samples_group_meta_predata_all-h2pc-arcs-error.csv"
}

# 保存过滤后的数据
write.csv(arcs_data, newname, row.names=FALSE)
arcs_data <- read.csv(newname, header=TRUE)

# 生成邻接矩阵
edge_list <- as.matrix(arcs_data[, c("from", "to")])
Graph <- graph_from_edgelist(edge_list, directed=TRUE)

# 确保节点集合完整
node_set <- c(Env_names, "group")
Graph <- add_vertices(Graph, nv=length(setdiff(node_set, V(Graph)$name)), name=setdiff(node_set, V(Graph)$name))

# 计算邻接矩阵
adj_matrix <- as.matrix(as_adj(Graph, sparse=FALSE))
colnames(adj_matrix) <- node_set
rownames(adj_matrix) <- node_set

# 保存邻接矩阵
netname <- "./keystone/pearson/Keystone_samples_group_meta_predata_all-h2pc-arcs-adj_matrix.csv"
write.csv(adj_matrix, netname)
