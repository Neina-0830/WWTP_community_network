##consumer-resource model analysis
library(reshape2)
library(ggplot2)
##底物浓度100
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_100_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_100_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_100<- read.table("./result/resource_100_all_composition.txt", header=TRUE, sep=',')

test_otu_100 <- all_100[all_100$Type=="consumer" & all_100$ID==4,]
par(mar=c(3,5,2,2))
p<-ggplot(test_otu_100, aes(Transfer, Abundance),groups=Well) +
  geom_line(aes(color = Well)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black")) +
  theme(legend.position = 'none')
p
dev.off()
ggsave("./resource_100-consumer4_change_line.pdf", width = 4, height = 3)
##群落稳定性计算
otu_100 <- all_100[all_100$Type=="consumer"&all_100$Transfer==30,]
otu_100$ID <- paste("Spe", otu_100$ID, sep = "")
otu_100_table <- dcast(otu_100, Well ~ ID, value.var="Abundance")
otu_100_table[is.na(otu_100_table)] <- 0
write.csv(otu_100_table, "resource_100_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_100_table <- read.csv("resource_100_otu_table.csv", header=TRUE, row.names = 1)
otu_100_rarefied_table <- t(apply(otu_100_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_100_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_100_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.2435762
##总物质浓度的利用比例计算
resource_100 <- all_100[all_100$Type=="resource",]
well_list <- unique(resource_100$Well)
Utilization_100 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_100[resource_100$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_100 <- rbind(Utilization_100, result)
}
colnames(Utilization_100) <- c("Well","Utilization_proportion")
write.csv(Utilization_100, "resource_100_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_100$Utilization_proportion) ##0.4526373
Utilization_100 <- read.csv("resource_100_all_utilization.csv", header=TRUE, row.names = 1)
result_100 <- cbind(AVD, Utilization_100)
result_100$Resource <- 100
write.csv(result_100, "resource_100_all_result.csv", quote = FALSE) 


##底物浓度250
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_250_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_250_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_250<- read.table("./result/resource_250_all_composition.txt", header=TRUE, sep=',')
##群落稳定性计算
otu_250 <- all_250[all_250$Type=="consumer"&all_250$Transfer==30,]
otu_250$ID <- paste("Spe", otu_250$ID, sep = "")
otu_250_table <- dcast(otu_250, Well ~ ID, value.var="Abundance")
otu_250_table[is.na(otu_250_table)] <- 0
write.csv(otu_250_table, "resource_250_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_250_table <- read.csv("resource_250_otu_table.csv", header=TRUE, row.names = 1)
otu_250_rarefied_table <- t(apply(otu_250_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_250_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_250_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.3567483
##总物质浓度的利用比例计算
resource_250 <- all_250[all_250$Type=="resource",]
well_list <- unique(resource_250$Well)
Utilization_250 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_250[resource_250$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_250 <- rbind(Utilization_250, result)
}
colnames(Utilization_250) <- c("Well","Utilization_proportion")
write.csv(Utilization_250, "resource_250_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_250$Utilization_proportion) ##0.8990928
Utilization_250 <- read.csv("resource_250_all_utilization.csv", header=TRUE, row.names = 1)
result_250 <- cbind(AVD, Utilization_250)
result_250$Resource <- 250
write.csv(result_250, "resource_250_all_result.csv", quote = FALSE) 


##底物浓度500
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_500_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_500_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_500<- read.table("./result/resource_500_all_composition.txt", header=TRUE, sep=',')
##群落稳定性计算
otu_500 <- all_500[all_500$Type=="consumer"&all_500$Transfer==30,]
otu_500$ID <- paste("Spe", otu_500$ID, sep = "")
otu_500_table <- dcast(otu_500, Well ~ ID, value.var="Abundance")
otu_500_table[is.na(otu_500_table)] <- 0
write.csv(otu_500_table, "resource_500_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_500_table <- read.csv("resource_500_otu_table.csv", header=TRUE, row.names = 1)
otu_500_rarefied_table <- t(apply(otu_500_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_500_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_500_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.3915243
##总物质浓度的利用比例计算
resource_500 <- all_500[all_500$Type=="resource",]
well_list <- unique(resource_500$Well)
Utilization_500 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_500[resource_500$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_500 <- rbind(Utilization_500, result)
}
colnames(Utilization_500) <- c("Well","Utilization_proportion")
write.csv(Utilization_500, "resource_500_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_500$Utilization_proportion) ##0.9617039
Utilization_500 <- read.csv("resource_500_all_utilization.csv", header=TRUE, row.names = 1)
result_500 <- cbind(AVD, Utilization_500)
result_500$Resource <- 500
write.csv(result_500, "resource_500_all_result.csv", quote = FALSE) 


##底物浓度750
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_750_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_750_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_750<- read.table("./result/resource_750_all_composition.txt", header=TRUE, sep=',')
##群落稳定性计算
otu_750 <- all_750[all_750$Type=="consumer"&all_750$Transfer==30,]
otu_750$ID <- paste("Spe", otu_750$ID, sep = "")
otu_750_table <- dcast(otu_750, Well ~ ID, value.var="Abundance")
otu_750_table[is.na(otu_750_table)] <- 0
write.csv(otu_750_table, "resource_750_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_750_table <- read.csv("resource_750_otu_table.csv", header=TRUE, row.names = 1)
otu_750_rarefied_table <- t(apply(otu_750_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_750_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_750_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.3998138
##总物质浓度的利用比例计算
resource_750 <- all_750[all_750$Type=="resource",]
well_list <- unique(resource_750$Well)
Utilization_750 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_750[resource_750$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_750 <- rbind(Utilization_750, result)
}
colnames(Utilization_750) <- c("Well","Utilization_proportion")
write.csv(Utilization_750, "resource_750_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_750$Utilization_proportion) ##0.973132
Utilization_750 <- read.csv("resource_750_all_utilization.csv", header=TRUE, row.names = 1)
result_750 <- cbind(AVD, Utilization_750)
result_750$Resource <- 750
write.csv(result_750, "resource_750_all_result.csv", quote = FALSE) 


##底物浓度1000
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_1000_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_1000_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_1000<- read.table("./result/resource_1000_all_composition.txt", header=TRUE, sep=',')
##群落稳定性计算
otu_1000 <- all_1000[all_1000$Type=="consumer"&all_1000$Transfer==30,]
otu_1000$ID <- paste("Spe", otu_1000$ID, sep = "")
otu_1000_table <- dcast(otu_1000, Well ~ ID, value.var="Abundance")
otu_1000_table[is.na(otu_1000_table)] <- 0
write.csv(otu_1000_table, "resource_1000_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_1000_table <- read.csv("resource_1000_otu_table.csv", header=TRUE, row.names = 1)
otu_1000_rarefied_table <- t(apply(otu_1000_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_1000_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_1000_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.404984
##总物质浓度的利用比例计算
resource_1000 <- all_1000[all_1000$Type=="resource",]
well_list <- unique(resource_1000$Well)
Utilization_1000 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_1000[resource_1000$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_1000 <- rbind(Utilization_1000, result)
}
colnames(Utilization_1000) <- c("Well","Utilization_proportion")
write.csv(Utilization_1000, "resource_1000_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_1000$Utilization_proportion) ##0.9778099
Utilization_1000 <- read.csv("resource_1000_all_utilization.csv", header=TRUE, row.names = 1)
result_1000 <- cbind(AVD, Utilization_1000)
result_1000$Resource <- 1000
write.csv(result_1000, "resource_1000_all_result.csv", quote = FALSE) 


##底物浓度1250
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_1250_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_1250_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_1250<- read.table("./result/resource_1250_all_composition.txt", header=TRUE, sep=',')

test_otu_1250 <- all_1250[all_1250$Type=="consumer" & all_1250$ID==101,]
par(mar=c(3,5,2,2))
p<-ggplot(test_otu_1250, aes(Transfer, Abundance),groups=Well) +
  geom_line(aes(color = Well)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black")) +
  theme(legend.position = 'none')
p
dev.off()
ggsave("./resource_1250-consumer101_change_line.pdf", width = 4, height = 3)

##群落稳定性计算
otu_1250 <- all_1250[all_1250$Type=="consumer"&all_1250$Transfer==30,]
otu_1250$ID <- paste("Spe", otu_1250$ID, sep = "")
otu_1250_table <- dcast(otu_1250, Well ~ ID, value.var="Abundance")
otu_1250_table[is.na(otu_1250_table)] <- 0
write.csv(otu_1250_table, "resource_1250_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_1250_table <- read.csv("resource_1250_otu_table.csv", header=TRUE, row.names = 1)
otu_1250_rarefied_table <- t(apply(otu_1250_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_1250_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_1250_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.4009965
##总物质浓度的利用比例计算
resource_1250 <- all_1250[all_1250$Type=="resource",]
well_list <- unique(resource_1250$Well)
Utilization_1250 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_1250[resource_1250$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_1250 <- rbind(Utilization_1250, result)
}
colnames(Utilization_1250) <- c("Well","Utilization_proportion")
write.csv(Utilization_1250, "resource_1250_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_1250$Utilization_proportion) ## 0.9819779
Utilization_1250 <- read.csv("resource_1250_all_utilization.csv", header=TRUE, row.names = 1)
result_1250 <- cbind(AVD, Utilization_1250)
result_1250$Resource <- 1250
write.csv(result_1250, "resource_1250_all_result.csv", quote = FALSE) 



##底物浓度1500
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_1500_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_1500_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_1500<- read.table("./result/resource_1500_all_composition.txt", header=TRUE, sep=',')
##群落稳定性计算
otu_1500 <- all_1500[all_1500$Type=="consumer"&all_1500$Transfer==30,]
otu_1500$ID <- paste("Spe", otu_1500$ID, sep = "")
otu_1500_table <- dcast(otu_1500, Well ~ ID, value.var="Abundance")
otu_1500_table[is.na(otu_1500_table)] <- 0
write.csv(otu_1500_table, "resource_1500_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_1500_table <- read.csv("resource_1500_otu_table.csv", header=TRUE, row.names = 1)
otu_1500_rarefied_table <- t(apply(otu_1500_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_1500_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_1500_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.4227345
##总物质浓度的利用比例计算
resource_1500 <- all_1500[all_1500$Type=="resource",]
well_list <- unique(resource_1500$Well)
Utilization_1500 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_1500[resource_1500$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_1500 <- rbind(Utilization_1500, result)
}
colnames(Utilization_1500) <- c("Well","Utilization_proportion")
write.csv(Utilization_1500, "resource_1500_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_1500$Utilization_proportion) ##0.980254
Utilization_1500 <- read.csv("resource_1500_all_utilization.csv", header=TRUE, row.names = 1)
result_1500 <- cbind(AVD, Utilization_1500)
result_1500$Resource <- 1500
write.csv(result_1500, "resource_1500_all_result.csv", quote = FALSE) 


##底物浓度2000
alldata <-  data.frame(exp_id=character(),Transfer=numeric(),Type=character(),ID=numeric(),Well=character(),Abundance=double(),stringsAsFactors=FALSE)
for(i in 1:100){
  name <- paste("./result/resource_2000_",as.character(i),"_composition.txt", sep = '')
  data <- read.csv(name,head=T)
  data$Well <- paste("W",as.character(i),sep = '')
  alldata <- rbind(alldata,data) 
}
write.table(alldata, "./result/resource_2000_all_composition.txt", quote=FALSE, row.names=FALSE, sep=',')
##预处理OTU表
all_2000<- read.table("./result/resource_2000_all_composition.txt", header=TRUE, sep=',')
##群落稳定性计算
otu_2000 <- all_2000[all_2000$Type=="consumer"&all_2000$Transfer==30,]
otu_2000$ID <- paste("Spe", otu_2000$ID, sep = "")
otu_2000_table <- dcast(otu_2000, Well ~ ID, value.var="Abundance")
otu_2000_table[is.na(otu_2000_table)] <- 0
write.csv(otu_2000_table, "resource_2000_otu_table.csv", quote = FALSE, row.names = FALSE)  
otu_2000_table <- read.csv("resource_2000_otu_table.csv", header=TRUE, row.names = 1)
otu_2000_rarefied_table <- t(apply(otu_2000_table, 1, function(x)  x/sum(x))) ##抽平操作（丰度归一化)
deviation_otu <- apply(otu_2000_rarefied_table, 2, function(x)  abs(x-mean(x))/sd(x))
AVD <- as.data.frame(rowSums(deviation_otu)/ncol(deviation_otu))
colnames(AVD) <- "AVD"
write.csv(AVD, "resource_2000_AVD.csv", quote = FALSE)  
mean(AVD$AVD) ##0.4351667
##总物质浓度的利用比例计算
resource_2000 <- all_2000[all_2000$Type=="resource",]
well_list <- unique(resource_2000$Well)
Utilization_2000 <- data.frame(Well = character(), Utilization_proportion=double())
for(i in 1:length(well_list)){
  well <- well_list[i]
  resource_select <- resource_2000[resource_2000$Well==well,]
  utilization <- (sum(resource_select[resource_select$Transfer==0,]$Abundance) - sum(resource_select[resource_select$Transfer==30,]$Abundance))/sum(resource_select[resource_select$Transfer==0,]$Abundance)
  result <- data.frame(well, utilization)
  Utilization_2000 <- rbind(Utilization_2000, result)
}
colnames(Utilization_2000) <- c("Well","Utilization_proportion")
write.csv(Utilization_2000, "resource_2000_all_utilization.csv", quote = FALSE, row.names = FALSE) 
mean(Utilization_2000$Utilization_proportion) ##0.9835639
Utilization_2000 <- read.csv("resource_2000_all_utilization.csv", header=TRUE, row.names = 1)
result_2000 <- cbind(AVD, Utilization_2000)
result_2000$Resource <- 2000
write.csv(result_2000, "resource_2000_all_result.csv", quote = FALSE) 


all_result <- rbind(result_100, result_250, result_500, result_750, result_1000, result_1500, result_2000)
##spearman相关性
cor.test(all_result$Resource, all_result$AVD,method="spearman")
'''
Spearmans rank correlation rho

data:  all_result$Resource and all_result$AVD
S = 9411995, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.8353584
'''
cor.test(all_result$Resource, all_result$Utilization_proportion,method="spearman")
'''
Spearmans rank correlation rho

data:  all_result$Resource and all_result$Utilization_proportion
S = 13218635, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
  rho 
0.7687698 
'''

all_result$Resource <- factor(all_result$Resource)
##AVD
par(mar=c(3,5,2,2))
p<- ggplot(all_result, aes(x=Resource, y=AVD)) + 
  geom_boxplot(alpha=0.8,fill="grey", show.legend=FALSE) +
  ylim(0,1)
p + labs(y = "AVD") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("All_resource_AVD-boxplot.pdf", width = 4, height = 3)
t.test(all_result[all_result$Resource==100,]$AVD, all_result[all_result$Resource==250,]$AVD) ##p-value < 2.2e-16(0.2435762 0.3567483)
t.test(all_result[all_result$Resource==250,]$AVD, all_result[all_result$Resource==500,]$AVD) ##p-value < 2.2e-16(0.3567483 0.3915243)
t.test(all_result[all_result$Resource==500,]$AVD, all_result[all_result$Resource==750,]$AVD) ##p-value = 0.006573(0.3915243 0.3998138)
t.test(all_result[all_result$Resource==750,]$AVD, all_result[all_result$Resource==1000,]$AVD) ##p-value = 0.06855(0.3998138 0.4049840)
t.test(all_result[all_result$Resource==1000,]$AVD, all_result[all_result$Resource==1500,]$AVD) ##p-value = 1.823e-08(0.4049840 0.4227345)
t.test(all_result[all_result$Resource==1500,]$AVD, all_result[all_result$Resource==2000,]$AVD) ##p-value = 0.0001011(0.4227345 0.4351667)


##Utilization_proportion
par(mar=c(3,5,2,2))
p<- ggplot(all_result, aes(x=Resource, y=Utilization_proportion)) + 
  geom_boxplot(alpha=0.8,fill="grey", show.legend=FALSE) +
  ylim(0,1)
p + labs(y = "Utilization_proportion") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("All_resource_Utilization_proportion-boxplot.pdf", width = 4, height = 3)
t.test(all_result[all_result$Resource==100,]$Utilization_proportion, all_result[all_result$Resource==250,]$Utilization_proportion) ##p-value < 2.2e-16(0.4526373 0.8990928)
t.test(all_result[all_result$Resource==250,]$Utilization_proportion, all_result[all_result$Resource==500,]$Utilization_proportion) ##p-value < 2.2e-16(0.8990928 0.9617039)
t.test(all_result[all_result$Resource==500,]$Utilization_proportion, all_result[all_result$Resource==750,]$Utilization_proportion) ##p-value = 5.36e-05(0.9617039 0.9731320)
t.test(all_result[all_result$Resource==750,]$Utilization_proportion, all_result[all_result$Resource==1000,]$Utilization_proportion) ##p-value = 0.07016(0.9731320 0.9778099)
t.test(all_result[all_result$Resource==1000,]$Utilization_proportion, all_result[all_result$Resource==1500,]$Utilization_proportion) ##p-value = 0.3428(0.9778099 0.9802540)
t.test(all_result[all_result$Resource==1500,]$Utilization_proportion, all_result[all_result$Resource==2000,]$Utilization_proportion) ##p-value = 0.1715(0.9802540 0.9835639)
