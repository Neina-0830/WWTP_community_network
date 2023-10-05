##比较两条回归线斜率和截距是否显著性差异
#install.packages("D:/Program Files/R/package/simba_0.3-5.tar.gz", repos = NULL, type = "source")
library(simba)
library(ggplot2)
##keystone
dat1 <- read.csv("./keystone/pearson/RMT-keystone-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
keystone_data <- dat1[dat1$Network=="keystone",]
randnode_data <- dat1[dat1$Network=="randnode",]
keystone_down_rate <- (keystone_data[keystone_data$Number.hub.removed==nrow(keystone_data)-1,]$remain.mean-keystone_data[keystone_data$Number.hub.removed==0,]$remain.mean)/(nrow(keystone_data)-1)
#-0.5077114
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==nrow(randnode_data)-1,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/(nrow(randnode_data)-1)
#-0.03300759
ratio <- keystone_down_rate/rand_down_rate
#15.38166
#斜率的差异比较，差异显著
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.4766 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0388 0.0517 0.0606 0.0711 
'''
#截距的差异比较，差异显著
diffic(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.1731 
Significance: 0.462 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
 3.12  3.90  4.55  5.17
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-Robustness-keystone-random_deletion-compare-nc-new.pdf", width = 5, height = 3)

##non-keystone
dat1 <- read.csv("./keystone/pearson/RMT-0.00001-0.01-no_keystone-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
non_keystone_data <- dat1[dat1$Network=="Others",]
randnode_data <- dat1[dat1$Network=="randnode",]
non_keystone_down_rate <- (non_keystone_data[non_keystone_data$Number.hub.removed==4900,]$remain.mean-non_keystone_data[non_keystone_data$Number.hub.removed==0,]$remain.mean)/4900
#-0.02058521
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==4900,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/4900
#-0.02906287
ratio <- non_keystone_down_rate/rand_down_rate
##0.7082993
#斜率的差异比较，差异显著
diffslope(non_keystone_data$Number.hub.removed, non_keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: 0.008902 
Significance: 0.001 

Empirical upper confidence limits of r:
     90%      95%    97.5%      99% 
0.000721 0.000899 0.001100 0.001350    
'''
#截距的差异比较，差异显著
diffic(non_keystone_data$Number.hub.removed, non_keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: -0.07986 
Significance: 0.487 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
 2.89  3.70  4.22  4.75
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,5000) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-Robustness-non-keystone-random_deletion-compare-nc-new.pdf", width = 5, height = 3)


##module hubs
dat1 <- read.csv("./hubs/pearson/RMT-hubs-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
hubs_data <- dat1[dat1$Network=="hubs",]
randnode_data <- dat1[dat1$Network=="randnode",]
hubs_down_rate <- (hubs_data[hubs_data$Number.hub.removed==nrow(hubs_data)-1,]$remain.mean-hubs_data[hubs_data$Number.hub.removed==0,]$remain.mean)/(nrow(hubs_data)-1)
#8.136965e-05
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==nrow(randnode_data)-1,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/(nrow(randnode_data)-1)
#-0.02877757
ratio <- hubs_down_rate/rand_down_rate
#-0.002827538
#斜率的差异比较，差异显著
diffslope(hubs_data$Number.hub.removed, hubs_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: 0.02994 
Significance: 0.001 

Empirical upper confidence limits of r:
    90%     95%   97.5%     99% 
0.00209 0.00272 0.00336 0.00378
'''
#截距的差异比较，差异显著
diffic(hubs_data$Number.hub.removed, hubs_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: 0.009214 
Significance: 0.485 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.223 0.290 0.358 0.409
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,150) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./hubs/pearson/RMT-Robustness-hubs-random_deletion-compare-nc-new.pdf", width = 5, height = 3)

##connectors
dat1 <- read.csv("./hubs/pearson/RMT-connectors-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
connectors_data <- dat1[dat1$Network=="connectors",]
randnode_data <- dat1[dat1$Network=="randnode",]
connectors_down_rate <- (connectors_data[connectors_data$Number.hub.removed==nrow(connectors_data)-1,]$remain.mean-connectors_data[connectors_data$Number.hub.removed==0,]$remain.mean)/(nrow(connectors_data)-1)
#9.82883e-05
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==nrow(randnode_data)-1,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/(nrow(randnode_data)-1)
# -0.03110075
ratio <- connectors_down_rate/rand_down_rate
#-0.003160319
#斜率的差异比较，差异显著
diffslope(connectors_data$Number.hub.removed, connectors_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: 0.03078 
Significance: 0.001 

Empirical upper confidence limits of r:
    90%     95%   97.5%     99% 
0.00262 0.00328 0.00380 0.00467  
'''
#截距的差异比较，差异显著
diffic(connectors_data$Number.hub.removed, connectors_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: -0.005083 
Significance: 0.483 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.190 0.245 0.289 0.359
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./hubs/pearson/RMT-Robustness-connectors-random_deletion-compare-nc-new.pdf", width = 5, height = 3)

##peripherals
dat1 <- read.csv("./hubs/pearson/RMT-peripherals-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
peripherals_data <- dat1[dat1$Network=="peripherals",]
randnode_data <- dat1[dat1$Network=="randnode",]
peripherals_down_rate <- (peripherals_data[peripherals_data$Number.hub.removed==4750,]$remain.mean-peripherals_data[peripherals_data$Number.hub.removed==0,]$remain.mean)/4750
#-0.02817848
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==4750,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/4750
#-0.02954263
ratio <- peripherals_down_rate/rand_down_rate
##0.9538244
#斜率的差异比较，差异显著
diffslope(peripherals_data$Number.hub.removed, peripherals_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: 5.448e-05 
Significance: 0.407 

Empirical upper confidence limits of r:
     90%      95%    97.5%      99% 
0.000379 0.000501 0.000583 0.000719
'''
#截距的差异比较，差异显著
diffic(peripherals_data$Number.hub.removed, peripherals_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: -2.138 
Significance: 0.001 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.633 0.817 0.930 1.152 
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,5000) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./hubs/pearson/RMT-Robustness-peripherals-random_deletion-compare-nc.pdf", width = 5, height = 3)

##core
dat1 <- read.csv("./igraph/pearson/RMT-core-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
core_data <- dat1[dat1$Network=="core",]
randnode_data <- dat1[dat1$Network=="randnode",]
core_down_rate <- (core_data[core_data$Number.node.removed==nrow(core_data)-1,]$remain.mean-core_data[core_data$Number.node.removed==0,]$remain.mean)/(nrow(core_data)-1)
#0.0002002215
rand_down_rate <- (randnode_data[randnode_data$Number.node.removed==nrow(randnode_data)-1,]$remain.mean-randnode_data[randnode_data$Number.node.removed==0,]$remain.mean)/(nrow(randnode_data)-1)
# -0.02964143
ratio <- core_down_rate/rand_down_rate
#-0.006754785
#斜率的差异比较，差异显著
diffslope(core_data$Number.node.removed, core_data$remain.mean, randnode_data$Number.node.removed, randnode_data$remain.mean)
'''
Difference in Slope: 0.03043 
Significance: 0.001 

Empirical upper confidence limits of r:
    90%     95%   97.5%     99% 
0.00161 0.00189 0.00231 0.00277 
'''
#截距的差异比较，差异显著
diffic(core_data$Number.node.removed, core_data$remain.mean, randnode_data$Number.node.removed, randnode_data$remain.mean)
'''
Difference in Slope: -0.01032 
Significance: 0.506 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.301 0.385 0.458 0.551
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.node.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,250) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./igraph/pearson/RMT-Robustness-core-random_deletion-compare-nc-new.pdf", width = 5, height = 3)

##non-core
dat1 <- read.csv("./igraph/pearson/RMT-non-core-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
non_core_data <- dat1[dat1$Network=="non-core",]
randnode_data <- dat1[dat1$Network=="randnode",]
non_core_down_rate <- (non_core_data[non_core_data$Number.hub.removed==4750,]$remain.mean-non_core_data[non_core_data$Number.hub.removed==0,]$remain.mean)/4750
#-0.02867832
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==4750,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/4750
#-0.02957319
ratio <- non_core_down_rate/rand_down_rate
##0.9697403
#斜率的差异比较，差异显著
diffslope(non_core_data$Number.hub.removed, non_core_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: -0.0007835 
Significance: 0.001 

Empirical upper confidence limits of r:
     90%      95%    97.5%      99% 
0.000330 0.000398 0.000465 0.000548   
'''
#截距的差异比较，差异显著
diffic(non_core_data$Number.hub.removed, non_core_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Difference in Slope: -1.146 
Significance: 0.01 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.620 0.798 0.924 1.026 
'''
##plot
par(mar=c(3,5,2,2))
ggplot(dat1, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red","blue"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,5000) +
  ylim(0,150) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./igraph/pearson/RMT-Robustness-non-core-random_deletion-compare-nc-new.pdf", width = 5, height = 3)


##将non-keystone抽平到keystone数量一样探究对网络稳定性的影响
dat1 <- read.csv("./keystone/pearson/RMT-keystone-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
keystone_data <- dat1[dat1$Network=="keystone",]
randnode_data <- dat1[dat1$Network=="randnode",]
dat2 <- read.csv("./keystone/pearson/RMT-no_keystone-rand89-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
Others_data <- dat2[dat2$Network=="Others",]
all_data <- rbind(keystone_data, Others_data, randnode_data)
all_data$Network <- factor(all_data$Network, levels = c("keystone","Others","randnode"))
write.csv(all_data, "./keystone/pearson/RMT-keystone-others-random_deletion_result-compare-nc.csv",quote = FALSE, row.names = FALSE)
keystone_down_rate <- (keystone_data[keystone_data$Number.hub.removed==nrow(keystone_data)-1,]$remain.mean-keystone_data[keystone_data$Number.hub.removed==0,]$remain.mean)/(nrow(keystone_data)-1)
#-0.5077114
Others_down_rate <- (Others_data[Others_data$Number.hub.removed==nrow(Others_data)-1,]$remain.mean-Others_data[Others_data$Number.hub.removed==0,]$remain.mean)/(nrow(Others_data)-1)
#--0.02172106
rand_down_rate <- (randnode_data[randnode_data$Number.hub.removed==nrow(randnode_data)-1,]$remain.mean-randnode_data[randnode_data$Number.hub.removed==0,]$remain.mean)/(nrow(randnode_data)-1)
#-0.03300759
ratio1 <- keystone_down_rate/rand_down_rate
#15.38166
ratio2 <- Others_down_rate/rand_down_rate
#0.6580627
#斜率的差异比较，差异显著
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.4766 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0388 0.0517 0.0606 0.0711 
'''
diffslope(Others_data$Number.hub.removed, Others_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = Others_data$Number.hub.removed, y1 = Others_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: 0.009873 
Significance: 0.001 

Empirical upper confidence limits of r:
    90%     95%   97.5%     99% 
0.00112 0.00150 0.00172 0.00199 
'''
#截距的差异比较，差异显著
diffic(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.1731 
Significance: 0.451 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
 2.99  3.84  4.41  5.28
'''
diffic(Others_data$Number.hub.removed, Others_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = Others_data$Number.hub.removed, y1 = Others_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.03307 
Significance: 0.264 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0626 0.0766 0.0891 0.1070
'''
##plot
par(mar=c(3,5,2,2))
ggplot(all_data, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line()+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("red", "blue", "gray"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(80,160) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-Robustness-keystone-others-random_deletion-compare-nc-new.pdf", width = 5, height = 3)


##hubs抽平到与keystone数量一样(module hubs与connectors均为hubs)
dat1 <- read.csv("./keystone/pearson/RMT-keystone-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
keystone_data <- dat1[dat1$Network=="keystone",]
randnode_data <- dat1[dat1$Network=="randnode",]
dat2 <- read.csv("./hubs/pearson/RMT-hubs-all-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
hubs_data <- dat2[dat2$Network=="hubs",][1:90,]
all_data <- rbind(hubs_data, keystone_data, randnode_data)
write.csv(all_data, "./keystone/pearson/RMT-keystone-hubs_deletion_result-compare-nc.csv",quote = FALSE, row.names = FALSE)
all_data$Network <- factor(all_data$Network, levels = c("hubs","keystone","randnode"))

hubs_down_rate <- (hubs_data[hubs_data$Number.hub.removed==nrow(hubs_data)-1,]$remain.mean-hubs_data[hubs_data$Number.hub.removed==0,]$remain.mean)/(nrow(hubs_data)-1)
#8.730899e-05
ratio <- hubs_down_rate/rand_down_rate
#-0.002645119
#斜率的差异比较，差异显著
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.4766 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0388 0.0517 0.0606 0.0711 
'''
diffslope(hubs_data$Number.hub.removed, hubs_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = hubs_data$Number.hub.removed, y1 = hubs_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: 0.03167 
Significance: 0.001 

Empirical upper confidence limits of r:
    90%     95%   97.5%     99% 
0.00268 0.00356 0.00423 0.00538 
'''
#截距的差异比较，差异显著
diffic(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.1731 
Significance: 0.462 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
 3.12  3.90  4.55  5.17
'''
diffic(hubs_data$Number.hub.removed, hubs_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = hubs_data$Number.hub.removed, y1 = hubs_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.03971 
Significance: 0.415 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.182 0.227 0.265 0.330 
'''
par(mar=c(3,5,2,2))
ggplot(all_data, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line(size=0.2,alpha=0.1)+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue", "red", "gray"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(80,160) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-Robustness-keystone-hubs-random_deletion-compare-nc-new.pdf", width = 5, height = 3)

##只比较keystone和hubs
all_data <- rbind(hubs_data, keystone_data)
write.csv(all_data, "./keystone/pearson/RMT-keystone-hubs_deletion_result-compare-nc-new.csv",quote = FALSE, row.names = FALSE)
all_data$Network <- factor(all_data$Network, levels = c("hubs","keystone"))

#斜率的差异比较，差异显著
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, hubs_data$Number.hub.removed, hubs_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = hubs_data$Number.hub.removed, y2 = hubs_data$remain.mean) 

Difference in Slope: -0.5083 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0454 0.0581 0.0654 0.0858 
'''
#截距的差异比较，差异显著
diffic(keystone_data$Number.hub.removed, keystone_data$remain.mean, hubs_data$Number.hub.removed, hubs_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.1334 
Significance: 0.485 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
 3.04  3.80  4.28  5.22
'''
par(mar=c(3,5,2,2))
ggplot(all_data, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line(size=0.2,alpha=0.1)+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue", "red"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(80,160) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-Robustness-keystone-hubs_deletion-compare-nc-new.pdf", width = 5, height = 3)



##core抽平到与keystone数量一样
dat1 <- read.csv("./keystone/pearson/RMT-keystone-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
keystone_data <- dat1[dat1$Network=="keystone",]
randnode_data <- dat1[dat1$Network=="randnode",]
dat2 <- read.csv("./igraph/pearson/RMT-core-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
core_data <- dat2[dat2$Network=="core",][1:90,]
all_data <- rbind(core_data, keystone_data, randnode_data)
write.csv(all_data, "./keystone/pearson/RMT-keystone-core_deletion_result-compare-nc.csv",quote = FALSE, row.names = FALSE)
all_data$Network <- factor(all_data$Network, levels = c("core","keystone","randnode"))

core_down_rate <- (core_data[core_data$Number.hub.removed==nrow(core_data)-1,]$remain.mean-core_data[core_data$Number.hub.removed==0,]$remain.mean)/(nrow(core_data)-1)
#0.0001890865
ratio <- core_down_rate/rand_down_rate
#-0.005728578
#斜率的差异比较，差异显著
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.4766 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0388 0.0517 0.0606 0.0711 
'''
diffslope(core_data$Number.hub.removed, core_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = core_data$Number.hub.removed, y1 = core_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: 0.03178 
Significance: 0.001 

Empirical upper confidence limits of r:
    90%     95%   97.5%     99% 
0.00290 0.00382 0.00453 0.00515
'''
#截距的差异比较，差异显著
diffic(keystone_data$Number.hub.removed, keystone_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.1731 
Significance: 0.462 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
 3.12  3.90  4.55  5.17
'''
diffic(core_data$Number.hub.removed, core_data$remain.mean, randnode_data$Number.hub.removed, randnode_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffic(x1 = core_data$Number.hub.removed, y1 = core_data$remain.mean,      x2 = randnode_data$Number.hub.removed, y2 = randnode_data$remain.mean) 

Difference in Slope: -0.0397 
Significance: 0.403 

Empirical upper confidence limits of r:
  90%   95% 97.5%   99% 
0.195 0.254 0.286 0.343
'''
par(mar=c(3,5,2,2))
ggplot(all_data, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line(size=0.2,alpha=0.1)+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue", "red", "gray"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(80,160) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/RMT-Robustness-keystone-core-random_deletion-compare-nc-new.pdf", width = 5, height = 3)


##module_hubs/connectors/keystone抽平后的比较
dat1 <- read.csv("./hubs/pearson/RMT-hubs-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
hubs_data <- dat1[dat1$Network=="hubs",][1:90,]
dat1 <- read.csv("./hubs/pearson/RMT-connectors-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
connectors_data <- dat1[dat1$Network=="connectors",][1:90,]
dat1 <- read.csv("./keystone/pearson/RMT-keystone-random_deletion_result-compare-nc.csv", header=TRUE, row.names = 1)
keystone_data <- dat1[dat1$Network=="keystone",]
all_data <- rbind(hubs_data, connectors_data, keystone_data)
write.csv(all_data, "./key/pearson/RMT-keystone-hubs-connectors_deletion_result-compare-nc.csv",quote = FALSE, row.names = FALSE)
all_data$Network <- factor(all_data$Network, levels = c("hubs","connectors","keystone"))
hubs_down_rate <- (hubs_data[hubs_data$Number.hub.removed==nrow(hubs_data)-1,]$remain.mean-hubs_data[hubs_data$Number.hub.removed==0,]$remain.mean)/(nrow(hubs_data)-1)
#8.218627e-05
connectors_down_rate <- (connectors_data[connectors_data$Number.hub.removed==nrow(connectors_data)-1,]$remain.mean-connectors_data[connectors_data$Number.hub.removed==0,]$remain.mean)/(nrow(connectors_data)-1)
#9.782012e-05
keystone_down_rate <- (keystone_data[keystone_data$Number.hub.removed==nrow(keystone_data)-1,]$remain.mean-keystone_data[keystone_data$Number.hub.removed==0,]$remain.mean)/(nrow(keystone_data)-1)
#-0.5077114
#斜率的差异比较，差异显著
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, hubs_data$Number.hub.removed, hubs_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = hubs_data$Number.hub.removed, y2 = hubs_data$remain.mean) 

Difference in Slope: -0.5083 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0443 0.0554 0.0681 0.0781  
'''
diffslope(keystone_data$Number.hub.removed, keystone_data$remain.mean, connectors_data$Number.hub.removed, connectors_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = keystone_data$Number.hub.removed, y1 = keystone_data$remain.mean,      x2 = connectors_data$Number.hub.removed, y2 = connectors_data$remain.mean) 

Difference in Slope: -0.5083 
Significance: 0.001 

Empirical upper confidence limits of r:
   90%    95%  97.5%    99% 
0.0454 0.0588 0.0677 0.0836 
'''
diffslope(hubs_data$Number.hub.removed, hubs_data$remain.mean, connectors_data$Number.hub.removed, connectors_data$remain.mean)
'''
Is difference in slope significant? 
Significance is based on 1000 permutations 

Call:
diffslope(x1 = hubs_data$Number.hub.removed, y1 = hubs_data$remain.mean,      x2 = connectors_data$Number.hub.removed, y2 = connectors_data$remain.mean) 

Difference in Slope: -1.977e-05 
Significance: 0.001 

Empirical upper confidence limits of r:
     90%      95%    97.5%      99% 
2.19e-06 2.74e-06 3.36e-06 4.05e-06 
'''
par(mar=c(3,5,2,2))
ggplot(all_data, aes(x=Number.hub.removed, y=remain.mean, group=Network, color=Network)) + 
  geom_line(size=0.2,alpha=0.1)+
  geom_pointrange(aes(ymin=remain.mean-remain.sd, ymax=remain.mean+remain.sd),size=0.2)+
  scale_color_manual(values=c("blue", "red", "purple"))+
  xlab("Number of node removed")+
  #ylab("Proportion of species remained")+
  ylab("Natural connectivity")+
  xlim(0,100) +
  ylim(80,160) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./key/pearson/RMT-Robustness-keystone-hubs-connectors-random_deletion-compare-nc-new.pdf", width = 5, height = 3)