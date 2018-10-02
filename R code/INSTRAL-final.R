
require(ggplot2)
require(reshape2)

all<-read.csv('all.hierarchical.RF.stat',sep=" ",header=T)
all <- all[which(all$Generations!= 500000 | !all$Rep %in% c(8,15,49) | all$Speciation_rate !=1e-06),]
head(all)
allm =melt(all, measure.vars=c("AST.RF","INST.RF")) 
names(allm)[9]="Method"
names(allm)[10]="RF"
head(allm)
# ggplot(data=all,aes(x=as.factor(Pruned),y=AST.RF-INST.RF,fill=as.factor(Genes)))+
#   geom_boxplot(outlier.size = 1,outlier.alpha = 0.6)+theme_bw()+
#   facet_grid(.~Generations)+
#   scale_x_discrete(labels=c("3/4","1/2","1/4"))

summary(aov(formula=RF~Method*(Genes+Generations+as.factor(Speciation_rate)+Pruned),data = allm[allm$Pruned==150,]))

b$diff=b$AST.br-b$INST.br
y=all[all$Pruned==150 & all$Genes==50,]
y=all
t.test( y$AST.RF, y$INST.RF,alternative = "less",paired = TRUE)

#----------------------Creating tables for single placement------------------------
library(Hmisc)

tableToLatex = function(d){
  is.num <- sapply(d, is.numeric)
  d[is.num] <- lapply(d[is.num], round, 4)
  options(scipen=999)
  rownames(d) <- NULL
  
  latex(d, file="",rowlabel.just="l" ,rowlabel = F ,math.row.names=T)
  
}
#tableToLatex(brd)
#tableToLatex(sum)
#-----------------------------Single placement results----------------------

require(reshape2)
instq<-read.csv('instral.qs.single.stat',sep=" ",header=T)
astq<-read.csv('astral.qs.single.stat',sep=" ",header=T)
mq=merge(instq,astq)
ising<-read.csv('instral.single.stat.all',sep=" ",header=T)
asing<-read.csv('astral.single.stat',sep=" ",header=T)
head(ising)
head(asing)
msing=merge(ising,asing)
#msing$d = msing$AST.RF-msing$INST.RF
mall=merge(mq,msing)
mall$dq=mall$AST.QS - mall$INST.QS
mall$d=mall$AST.RF - mall$INST.RF
mall<- mall[which(mall$Generations!= 500000 | !mall$Rep %in% c(8,15,49) | mall$Speciation_rate !=1e-06),]
nrow(msing)

mall[mall$Generations==10000000 & mall$Genes==200 & mall$d !=0,]

nrow(mall)
mall[mall$dq == 0 & mall$d != 0,]$d=0
mall[mall$dq == 0 & mall$d != 0,]
mall[mall$dq != 0 & mall$d == 0,]
drf=dcast(Generations~Genes,data=mall,fun.aggregate = function(x) paste(sum(x!=0),sum(x>0),sep="/"),value.var = "d")
ddq=dcast(Generations~Genes,data=mall,fun.aggregate = function(x) sum(x!=0),value.var = "dq")
tableToLatex(ddq)
tableToLatex(drf)

head(levels(as.factor(ising$Generations)))
dcast(Generations~Genes,data=msing,fun.aggregate = function(x) sum(x!=0))
dcast(Generations~Genes,data=msing,fun.aggregate = function(x) sum(x<0))
dcast(Generations~Genes,data=msing,fun.aggregate = function(x) mean(x,na.rm = T))

# dcast(Generations~Genes,data=sing,fun.aggregate = function(x) sum(x!=0))
# dcast(Generations~Genes,data=sing,fun.aggregate = function(x) max(x))
# 
# dcast(Generations~Genes,data=sing,fun.aggregate = function(x) mean(x,na.rm = T))
# brd = dcast(Generations~Genes,data=sing,fun.aggregate = function(x) mean(x*198,na.rm = T))
# sum = dcast(Generations~Genes,data=sing,fun.aggregate = function(x) sum(x!=0))
ggplot(data=sing,aes(INST.RF))+geom_histogram(binwidth = 0.01)+theme_bw()+facet_grid(Genes~Generations)+scale_fill_brewer(palette = "YlGnBu")

#--------------------------------this change should be after running t-tests-------------------------

all$Generations = as.factor(all$Generations)
levels(all$Generations)
levels(all$Generations) <- c("Very High", "High", "Moderate")


ggplot(data=all,aes(x=as.factor(Pruned),y=AST.RF-INST.RF,fill=as.factor(Genes)))+
  geom_boxplot(outlier.size = 1,outlier.alpha = 0.4)+theme_bw()+
  facet_grid(.~Generations)+scale_fill_brewer(palette = "YlGnBu")+
  scale_x_discrete(labels=c("3/4","1/2","1/4"))+ xlab("Portion of leaves in the backbone tree")+
  guides(fill=guide_legend(title="# Genes"))+ ylab(expression(Delta~RF~(ASTRAL-INSTRAL)))
ggsave("RF-final.pdf",height = 3.5 , width = 8)


#----------------------QS-RF Hierarchical--------------------
a<-read.csv('all.hierarchical.stat.final',sep=" ",header=T)
head(a)
a <- a[which(a$Generations!= 500000 | !a$Rep %in% c(8,15,49) | a$Speciation_rate !=1e-06),]
nrow(a)
require(scales)
a$Generations = as.factor(a$Generations)
levels(a$Generations)
levels(a$Generations) <- c("Very High", "High", "Moderate")
a$Pruned = as.factor(a$Pruned)
levels(a$Pruned)<- c("3/4", "1/2", "1/4")
ggplot(data=a[a$Pruned=="1/4",], aes(y=AST.QS-INST.QS, x =AST.RF-INST.RF ))+geom_point(alpha=0.33)+ geom_smooth(method="lm",se=F)+ 
  facet_grid(Generations~Genes)+ theme_bw()+ylab(expression(Delta~QS~(ASTRAL-INSTRAL)))+ xlab(expression(Delta~RF~(ASTRAL-INSTRAL)))+
  geom_rug(col=rgb(.9,0.2,0,alpha=.1),sides="bl",size=1.5)+coord_cartesian(ylim = c(-0.001,0.001))+
  scale_x_continuous(labels=percent)+
  geom_hline(color="black",linetype="dashed",yintercept = 0)+geom_vline(color="black",linetype="dashed",xintercept = 0)
ggsave("RF-QS-150-final.pdf",width = 7.1,height = 4)

ggplot(data=a, aes(y=AST.QS-INST.QS, x =AST.RF-INST.RF ))+geom_point(alpha=0.33)+ geom_smooth(method="lm",se=F)+ 
  facet_grid(Generations~Pruned)+ theme_bw()+ylab(expression(Delta~QS~(ASTRAL-INSTRAL)))+ xlab(expression(Delta~RF~(ASTRAL-INSTRAL)))+
  geom_rug(col=rgb(.9,0.2,0,alpha=.1),sides="bl",size=1.5)+coord_cartesian(ylim = c(-0.001,0.001))+
  scale_x_continuous(labels=percent)+
  geom_hline(color="black",linetype="dashed",yintercept = 0)+geom_vline(color="black",linetype="dashed",xintercept = 0)
ggsave("RF-QS-all-final.pdf",width = 7.1,height = 4)






#-----------------------------Running time of instral ----------------------------


# rt<-read.csv('runningtime.hier.stat',sep=" ",header=T)
# head(rt)
# rt$Pruned = as.factor(rt$Pruned)
# levels(rt$Pruned)
# ggplot(data=rt, aes(y=RT, x =Gene,color=as.factor(Pruned) ))+geom_point(alpha=0.33)
# ggplot(data=rt, aes(y=RT, x =factor(Pruned, levels = rev(levels(Pruned))),color=as.factor(Gene) ))+geom_jitter(alpha=0.53)+ 
#   stat_summary(fun.y = "mean", colour = "red", size = 1, geom = "point")+geom_smooth(method = "lm")+theme_bw()+xlab("N")
# 
# rta<-read.csv('runningtime.hier.stat.detailed.all',sep=" ",header=T)
# head(rta)
# nrow(rta)
# rta$N = as.factor(rta$N)
# levels(rta$N)
# head(levels(rta$RT))
# rta=droplevels(rta[rta$RT != "-",])
# rta$RT=as.numeric(as.character(rta$RT))
# ggplot(data=rta, aes(y=RT, x =N,color=as.factor(Gene) ))+theme_bw()+facet_grid(Pruned~Generation)+stat_summary(size=0.1,alpha=0.5)
#   scale_y_log10()+scale_x_log10()
#   stat_summary(fun.y = "mean", colour = "red", size = 1, geom = "point")+geom_smooth(method = "lm")+xlab("N")
#  -------------------------------------------
rt2<-read.csv('RunningTimes/instral.runningtime',sep=" ",header=T)
rt2$RT=as.numeric(as.character(rt2$RT))
rt2$N = as.factor(rt2$N)
levels(rt2$N)
rt2<-rt2[which( !rt2$Rep %in% c(8,15,49) | rt2$Speciation_rate !=1e-06),]
ggplot(data=rt2[rt2$Speciation_rate==0.000001,], aes(y=(RT), x =as.numeric(N)+51 ,color=as.factor(Gene) ))+theme_bw()+geom_smooth()+theme_bw()+xlab("Size of backbone tree")+
  ylab("Running Time (secs)")+scale_color_brewer(palette = "YlGnBu", direction=-1)+ guides(color=guide_legend(title="#Genes")) + ggsave("RT.pdf",height = 3.5,width = 5.5)
+scale_y_log10()+scale_x_log10() 
stat_summary(fun.y = "mean", colour = "red", size = 1, geom = "point")+geom_smooth(method = "lm")+xlab("N")
#  

