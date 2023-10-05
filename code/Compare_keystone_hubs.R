library(ggplot2)
library(igraph)
library(dplyr)
library(reshape2)

Node_key <- read.csv('./Gephi/pearson/RMT-ASV_0.00001_0.01_Node_key_MRA.csv',header=T,row.names = 1)
zi_pi_type <- read.csv("./hubs/pearson/RMT-Net_0.00001_0.01_node_topology_pearson_zi_pi-type.csv",header=T,row.names=1)
Node_key$type <- zi_pi_type[rownames(Node_key),]$type
Node_hubs <- Node_key[Node_key$type !="Keystone",]
keystone <- read.csv("./keystone/pearson/RMT-ASV_0.00001_0.01_node_keystone.csv",header=T,row.names = 1)
Node_keystone <- Node_key[rownames(Node_key)%in%keystone$Label,]
Node_keystone$type <- "Keystone"
Node_key <- rbind(Node_hubs, Node_keystone)
write.csv(Node_key, "./key/pearson/RMT-ASV_0.00001_0.01_node_key_type-all.csv",quote = FALSE)


##module hubs与connectors均为hubs
Node_key <- read.csv('./key/pearson/RMT-ASV_0.00001_0.01_node_key_type-all.csv',header=T,row.names = 1)
Node_key[Node_key$type != "Keystone",]$type <- "Hubs" 
write.csv(Node_key, "./key/pearson/RMT-ASV_0.00001_0.01_node_key_type-all-new.csv",quote = FALSE)

##其他分布特征
Node_key$type <- factor(Node_key$type, levels = c("Hubs","Keystone"))
##frequency
par(mar=c(3,5,2,2))
p<- ggplot(Node_key, aes(x=type, y=v_frequency, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"),
                    limits = c("Hubs","Keystone")) +
  ylim(0,1.0)
p + labs(y = "Occurrence frequency") +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/type_frequency_compare-boxplot-new.pdf", width = 4, height = 3)
t.test(Node_key[Node_key$type=="Hubs",]$v_frequency, Node_key[Node_key$type=="Keystone",]$v_frequency) ##p-value = 3.137e-12(0.07897570 0.02737933)

##MRA
par(mar=c(3,5,2,2))
p<- ggplot(Node_key, aes(x=type, y=MRA, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"),
                    limits = c("Hubs","Keystone")) +
  ylim(0,0.005)
p + labs(y = "Average relative abundance") +
  #theme(axis.text.x = element_text(angle=30, hjust=0.5, vjust=0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/type_MRA_compare-boxplot-new.pdf", width = 4, height = 3)
t.test(Node_key[Node_key$type=="Hubs",]$MRA, Node_key[Node_key$type=="Keystone",]$MRA) ##p-value = 1.889e-06(1.992021e-04 4.809765e-05)

##node centrality
node_centrality <- read.csv("./igraph/pearson/RMT-ASV_0.00001_0.01_Node_all_centrality.csv",header=T, row.names = 1) ##共4992个ASVs
node_centrality_key <- node_centrality[rownames(Node_key),]
node_centrality_key$v_name <- rownames(node_centrality_key)
key_type_centrality <- merge(node_centrality_key, Node_key, by="v_name")
key_type_centrality$type <- factor(key_type_centrality$type, levels = c("Hubs","Keystone"))
##Degree_centrality
par(mar=c(3,5,2,2))
p<- ggplot(key_type_centrality, aes(x=type, y=Degree_centrality, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"),
                    limits = c("Hubs","Keystone")) +
  ylim(0,0.06)
p + labs(y = "Degree centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/Difftype_Degree_centrality_compare-boxplot-new.pdf", width = 4, height = 3)
t.test(key_type_centrality[key_type_centrality$type=="Hubs",]$Degree_centrality, key_type_centrality[key_type_centrality$type=="Keystone",]$Degree_centrality) ##p-value < 2.2e-16(0.008695116 0.029466523)

##Betweenness_centrality
par(mar=c(3,5,2,2))
p<- ggplot(key_type_centrality, aes(x=type, y=Betweenness_centrality, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"),
                    limits = c("Hubs","Keystone")) +
  ylim(0,0.02)
p + labs(y = "Betweenness centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/Difftype_Betweenness_centrality_compare-boxplot-new.pdf", width = 4, height = 3)
t.test(key_type_centrality[key_type_centrality$type=="Hubs",]$Betweenness_centrality, key_type_centrality[key_type_centrality$type=="Keystone",]$Betweenness_centrality) ##p-value = 5.146e-16(0.0027597618 0.0008313472)

##Closeness_centrality
par(mar=c(3,5,2,2))
p<- ggplot(key_type_centrality, aes(x=type, y=Closeness_centrality, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"),
                    limits = c("Hubs","Keystone")) +
  ylim(0,0.4)
p + labs(y = "Closeness centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/Difftype_Closeness_centrality_compare-boxplot-new.pdf", width = 4, height = 3)
t.test(key_type_centrality[key_type_centrality$type=="Hubs",]$Closeness_centrality, key_type_centrality[key_type_centrality$type=="Keystone",]$Closeness_centrality) ##p-value < 2.2e-16(0.2282138 0.2753382)

##Eigenvector_centrality
par(mar=c(3,5,2,2))
p<- ggplot(key_type_centrality, aes(x=type, y=Eigenvector_centrality, fill=type)) + 
  geom_boxplot(show.legend=FALSE) +
  scale_fill_manual(values = c("#CC79A7", "#F0E442"),
                    limits = c("Hubs","Keystone")) +
  ylim(0,1)
p + labs(y = "Eigenvector centrality") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/Difftype_Eigenvector_centrality_compare-boxplot-new.pdf", width = 4, height = 3)
t.test(key_type_centrality[key_type_centrality$type=="Hubs",]$Eigenvector_centrality, key_type_centrality[key_type_centrality$type=="Keystone",]$Eigenvector_centrality) ##p-value < 2.2e-16(0.003470863 0.479495288)

