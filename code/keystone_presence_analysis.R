library(dplyr)
library(reshape2)
group <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group.csv',header = T, row.names = 1)
env <- read.csv('../data/environment_all-norm-influx.csv', row.names = 1,header=TRUE)
env <- env[group$sample,]
env_group <- cbind(env, group)
t_test <- read.csv('./keystone/pearson/Keystone_samples_group-diffenv-t_test-sig.csv',header = T)
t_test_select <- t_test[t_test$Sig !="n.s",]
t_test_select$illu <- "keystone_big"
t_test_select[t_test_select$Keystone_mean <= t_test_select$No_keystone_mean,]$illu <- "no_keystone_big"
t_test_select_order <- arrange(t_test_select, illu)
t_test_select_order <- t_test_select_order[!t_test_select_order$para %in%c("EffBOD","EffNH4N"),]
env_group_sig <- env_group[,c(t_test_select_order$para,"group")]
env_group_sig_long <- melt(env_group_sig,id.vars = c('group'), variable.name='para', value.name='values')
env_group_sig_long$group <- factor(env_group_sig_long$group, levels = c("Keystone_sample","No_keystone_sample"))
env_group_sig_long$para <- factor(env_group_sig_long$para, levels =c("SMP","AtInfBOD","InfTN","Lat","HRT","AtHRT","SRT","Nitri","MLSS","Pop","MAT","AMMaxT","AMMinT","SMT","MIT",
                                                                     "DC","InfR","IndPer","InfNH4N","AtInfNH4N","InfBOD","ReInfBOD","InfCOD","ReInfCOD","AtInfTP","NH4N","Con","F.M","SVI") )

par(mar=c(3,5,2,2))
p<- ggplot(env_group_sig_long, aes(x=para, y=values, fill=group)) + 
  geom_boxplot(aes(fill=group), outlier.colour="white", width=0.5,size=0.2) +##异常点去除
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,100)
p + labs(y = "Normed parameters values") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"),
        legend.position="none")
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_para_diff_sig_compare-boxplot.pdf", width = 4, height = 3.2)


##基于贝叶斯网络分析的环境因子与是否包含keystone类群的因果推断
library(bnlearn)
group <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group.csv',header = T, row.names = 1)
env <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env-influx.csv', row.names = 1,header=TRUE)
env$sample <- rownames(env)
env_group <- merge(env, group, by="sample")
write.csv(env_group, "./keystone/pearson/Keystone_samples_group_meta_predata_influx-all.csv",quote = FALSE, row.names = FALSE)
env_group <- read.csv('./keystone/pearson/Keystone_samples_group_meta_predata_influx-all.csv',header = T, row.names = 1)
##分组转化为数值型变量
env_group[env_group$group=="Keystone_sample",]$group <- 1
env_group[env_group$group=="No_keystone_sample",]$group <- 0
env_group$group <- as.numeric(env_group$group)
##将字符型数据转化为数值型数据
t <- data.frame(matrix(as.numeric(unlist(env_group)),ncol = length(env_group[1,])))
rownames(t) <- rownames(env_group)
colnames(t) <- colnames(env_group)
env_group <- t
bayes_corr <- h2pc(env_group)
test <- bayes_corr$arcs
write.csv(bayes_corr$arcs,"./keystone/pearson/Keystone_samples_group_meta_predata_influx-h2pc-arcs.csv", quote=FALSE, row.names = FALSE)
bayes_corr1 <- mmhc(env_group)
test1 <- bayes_corr1$arcs
write.csv(bayes_corr1$arcs,"./keystone/pearson/Keystone_samples_group_meta_predata_influx-mmhc-arcs.csv", quote=FALSE, row.names = FALSE)
bayes_corr2 <- rsmax2(env_group,restrict = "hpc", maximize = "tabu")
test2 <- bayes_corr2$arcs
write.csv(bayes_corr2$arcs,"./keystone/pearson/Keystone_samples_group_meta_predata_influx-rsmax2-arcs.csv", quote=FALSE, row.names = FALSE)

##标准化环境因子
group <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group.csv',header = T, row.names = 1)
env <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/environment_all-norm-influx.csv', row.names = 1,header=TRUE)
env$sample <- rownames(env)
env_group <- merge(env, group, by="sample")
write.csv(env_group, "./keystone/pearson/Keystone_samples_group_meta_normed_influx-all.csv",quote = FALSE, row.names = FALSE)
env_group <- read.csv('./keystone/pearson/Keystone_samples_group_meta_normed_influx-all.csv',header = T, row.names = 1)
##分组转化为数值型变量
env_group[env_group$group=="Keystone_sample",]$group <- 1
env_group[env_group$group=="No_keystone_sample",]$group <- 0
env_group$group <- as.numeric(env_group$group)
##将字符型数据转化为数值型数据
t <- data.frame(matrix(as.numeric(unlist(env_group)),ncol = length(env_group[1,])))
rownames(t) <- rownames(env_group)
colnames(t) <- colnames(env_group)
env_group <- t
bayes_corr <- h2pc(env_group)
test <- bayes_corr$arcs
write.csv(bayes_corr$arcs,"./keystone/pearson/Keystone_samples_group_meta_normed_influx-h2pc-arcs.csv", quote=FALSE, row.names = FALSE)

##根据校正的h2pc结果，分析核心类群出现的可能原因
library(ggplot2)
group <- read.csv('./keystone/pearson/RMT-Net_0.00001_0.01_keystone_samples_group.csv',header = T, row.names = 1)
env <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env-influx.csv', row.names = 1,header=TRUE)
env <- env[group$sample,]
env_group <- cbind(env, group)
##得到SRT->F/M->contain keystone?的关系
##SRT
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=SRT, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,100)
p + labs(y = "SRT") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_SRT_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$SRT, env_group[env_group$group=="No_keystone_sample",]$SRT) ##p-value = 0.0008067(16.69897  14.69626)
##F/M
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=F.M, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "F/M") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_F.M_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$F.M, env_group[env_group$group=="No_keystone_sample",]$F.M) ##p-value = 5.221e-12(0.1467996 0.2255089)

##得到SRT->ReInfBOD->contain keystone?和InfBOD->ReInfBOD->contain keystone?的关系
##ReInfBOD
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=ReInfBOD, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,300)
p + labs(y = "ReInfBOD") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_ReInfBOD_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$ReInfBOD, env_group[env_group$group=="No_keystone_sample",]$ReInfBOD) ##p-value = 2.543e-11(97.49454 119.92551)
##InfBOD
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=InfBOD, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1000)
p + labs(y = "InfBOD") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_InfBOD_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$InfBOD, env_group[env_group$group=="No_keystone_sample",]$InfBOD) ##p-value = 6.108e-08(175.0888  211.8563)


##得到SMP->SVI->contain keystone? / SMT->SVI->contain keystone? / MIT->SVI->contain keystone? / Denitri->SVI->contain keystone? 和IndConInf->SVI->contain keystone?的关系
##SVI
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=SVI, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,400)
p + labs(y = "SVI") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_SVI_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$SVI, env_group[env_group$group=="No_keystone_sample",]$SVI) ##p-value = 5.557e-06(198.0386  219.1767)
##SMP
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=SMP, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,500)
p + labs(y = "SMP") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_SMP_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$SMP, env_group[env_group$group=="No_keystone_sample",]$SMP) ##p-value = 0.007395(117.7137  106.7961)
##SMT
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=SMT, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,30)
p + labs(y = "SMT") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_SMT_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$SMT, env_group[env_group$group=="No_keystone_sample",]$SMT) ##p-value = 1.052e-05(16.66412  18.81828)
##MIT
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=MIT, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,40)
p + labs(y = "MIT") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_MIT_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$MIT, env_group[env_group$group=="No_keystone_sample",]$MIT) ##p-value = 3.257e-05(23.33731  24.44403)
##Denitri
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=Denitri, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1.0)
p + labs(y = "Denitri") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_Denitri_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$Denitri, env_group[env_group$group=="No_keystone_sample",]$Denitri) ##p-value = 0.9845(0.6221374 0.6216012)
##IndConInf
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=group, y=IndConInf, fill=group)) + 
  geom_boxplot(alpha=0.8, show.legend=FALSE) +
  scale_fill_manual(values = c('#556B2F','#800000')) +
  ylim(0,1)
p + labs(y = "IndConInf") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/Keystone_samples_group_IndConInf_compare-boxplot.pdf", width = 4, height = 3)
t.test(env_group[env_group$group=="Keystone_sample",]$IndConInf, env_group[env_group$group=="No_keystone_sample",]$IndConInf) ##p-value = 0.08535(0.6660305 0.7122356)


##直接分析F/M与AS性能的相关性
library(ggplot2)
env <- read.csv('E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_env.csv', row.names = 1,header=TRUE)
performance <- read.csv("E:/Desktop/数据库数据下载/SRA数据库/WWTP/WWTP_16S/PRJNA509305/WWTP_meta_predata_all_performance.csv", row.names = 1,header=TRUE)
env_performance <- cbind(env, performance)
##F.M~BOD_remove_rate
par(mar=c(3,5,2,2))
p<- ggplot(env_performance, aes(x=F.M, y=BOD_remove_rate)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(BOD_remove_rate ~ F.M, data = env_performance)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 6.0, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(BOD)/(kg(MLSS)*d))", y = "BOD removal rate(g/(g(biomass)*d))") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_BOD_performance-linear.pdf", width = 4, height = 3)

##F.M~COD_remove_rate
par(mar=c(3,5,2,2))
p<- ggplot(env_performance, aes(x=F.M, y=COD_remove_rate)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(COD_remove_rate ~ F.M, data = env_performance)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 16, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(BOD)/(kg(MLSS)*d))", y = "COD removal rate(g/(g(biomass)*d))") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_COD_performance-linear.pdf", width = 4, height = 3)

##F.M~NH4N_remove_rate
par(mar=c(3,5,2,2))
p<- ggplot(env_performance, aes(x=F.M, y=NH4N_remove_rate)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(NH4N_remove_rate ~ F.M, data = env_performance)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 0.9, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(BOD)/(kg(MLSS)*d))", y = "NH4N removal rate(g/(g(biomass)*d))") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_NH4N_performance-linear.pdf", width = 4, height = 3)

##F.M~TN_remove_rate
par(mar=c(3,5,2,2))
p<- ggplot(env_performance, aes(x=F.M, y=TN_remove_rate)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(TN_remove_rate ~ F.M, data = env_performance)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 1.2, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(BOD)/(kg(MLSS)*d))", y = "TN removal rate(g/(g(biomass)*d))") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_TN_performance-linear.pdf", width = 4, height = 3)

##F.M~TP_remove_rate
par(mar=c(3,5,2,2))
p<- ggplot(env_performance, aes(x=F.M, y=TP_remove_rate)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(TP_remove_rate ~ F.M, data = env_performance)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 0.3, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(TP)/(kg(MLSS)*d))", y = "TP removal rate(g/(g(biomass)*d))") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_TP_performance-linear.pdf", width = 4, height = 3)


##F.M~ReInfBOD
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=F.M, y=ReInfBOD)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(ReInfBOD ~ F.M, data = env_group)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 300, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(TP)/(kg(MLSS)*d))", y = "ReInfBOD(mg/L)") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_ReInfBOD-linear.pdf", width = 4, height = 3)

##F.M~SVI
par(mar=c(3,5,2,2))
p<- ggplot(env_group, aes(x=F.M, y=SVI)) + 
  geom_point(size=1)  +
  geom_smooth(method="lm" , formula = y~x, color="#009E73", fill="#56B4E9", se=TRUE,linetype = "dashed",size=1)
model.lm <- lm(SVI ~ F.M, data = env_group)
l <- list(a = format(coef(model.lm)[1], digits = 4),
          b = format(abs(coef(model.lm)[2]), digits = 4),
          r2 = format(summary(model.lm)$r.squared, digits = 4),
          p = format(summary(model.lm)$coefficients[2,4], digits = 4))
eq <- substitute(italic(R)^2~"="~r2~","~italic(P)~"="~p, l)
p + xlim(0,1) +
  geom_text(aes(x = 0.25, y = 400, label = as.character(as.expression(eq))), parse = TRUE, size = 3) +
  labs(x="F/M(kg(TP)/(kg(MLSS)*d))", y = "SVI(ml/g)") +
  theme(axis.text = element_text(face = "bold", size = 10, color = "black"),
        axis.title = element_text(face = "bold", size = 10, color = "black"),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(fill = "transparent", color = "transparent"),
        axis.line = element_line(color = "black"))
dev.off()
ggsave("./keystone/pearson/AS_F.M_SVI-linear.pdf", width = 4, height = 3)
