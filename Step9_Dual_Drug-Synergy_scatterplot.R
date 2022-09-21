argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/Dual.plot/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/Synergy.score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/RDA/")
# argv=c("~/TCGA-LINCS/LINCS.Score/",
#        "~/TCGA-LINCS/TCGA.Score/",
#        "~/TCGA-LINCS/script/",
#        "~/TCGA-LINCS/Comb.plot/",
#        "~/TCGA-LINCS/",
#        "~/TCGA-LINCS/Synergy.score/",
#        "~/TCGA-LINCS/script/RDA/")

lincs_score = read.table(paste0(argv[1],"_total_lincs_score.tsv"),
                         sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
tcga_score = read.table(paste0(argv[2],"_total_tcga_score.tsv"),
                        sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)

tf.list=read.table(paste0(argv[3],"TF_list.txt"),
                   sep="\t",stringsAsFactors = FALSE,header=TRUE)
rowName = lincs_score$Symbol
lincs_score = lincs_score[,3:ncol(lincs_score), drop = FALSE]
tcga_score = tcga_score[,3:ncol(tcga_score), drop = FALSE]
drugname_set = colnames(lincs_score)
cancername_set = colnames(tcga_score)
tot.drugset=read.table(paste0(argv[3],"_total_drug-score_results.tsv"),
                       sep='\t',header=T,stringsAsFactors=F,quote='',
                       check.names=F,comment.char="")
rownames(tot.drugset)=tot.drugset[,1]
tot.drugset=tot.drugset[,-1]
# lincs.score.1=drug 1 lincs score
# lincs.score.2=drug 2 lincs score
# tcga.score=BRCA tcga score

print("#####################################################################################")
print("#######################                                       #######################")
print("#######################    Calculating  Single Comb-score     #######################")
print("#######################                                       #######################")
print("#####################################################################################")

#rm(list=ls())
args=commandArgs(trailingOnly = TRUE)
#args[1]=Compared Drug 1 or all
#args[2]=Compared Drug 2 or all
#args[3]= cancer type TCGA
args[1]="taselisib"
args[2]="acivicin"
args[3]="BRCA"
lincs.score.1=lincs_score[,args[1]]
lincs.score.2=lincs_score[,args[2]]
tcga.score=tcga_score[,args[3]]



mixed.score=c()

for ( i in 1:nrow(lincs_score)){
  if(abs(lincs.score.1[i]) >= abs(lincs.score.2[i])){
    mixed.score[i]=lincs.score.1[i]
  }else{
    mixed.score[i]=lincs.score.2[i]
  }
}
dual.score=(-1)*(sum(mixed.score * tcga.score)/nrow(tf.list))
syn.score=dual.score/mean(tot.drugset["BRCA",args[1]],tot.drugset["BRCA",args[2]])
dual.score
syn.score


print("#####################################################################################")
print("#######################                                       #######################")
print("#######################          Plotting scatter plot        #######################")
print("#######################          with comb-scores by drug     #######################")
print("#######################                                       #######################")
print("#####################################################################################")

if(!require(ggplot2)){BiocManager::install("ggplot2")};library(ggplot2)
if(!require(reshape2)){BiocManager::install("reshape2")};library(reshape2)

#rm(list=ls())

c.type=args[3]
result_matrix = cbind.data.frame(tcga_score[,c.type,drop=FALSE],mixed.score,lincs_score)
melted_result = melt(result_matrix, id.vars=c.type)

print(args[1])
print(args[2])

png(filename = paste0(argv[4],c.type,'_',args[1],"-",args[2],'.png'),width = 750, height = 750, type = "cairo")

#"lm", "glm", "gam", "loess"
plt=ggplot(melted_result[melted_result$variable != "mixed.score",])+
  geom_point(aes(x=value,y=eval(parse(text = c.type)),col=variable),col = "#999999",size=0.1)+
  geom_point(data = subset(melted_result,variable == args[1]),aes(x=value,y=eval(parse(text = c.type))),
             col="#ff2400",size=0.9,pch=15)+
  geom_point(data = subset(melted_result,variable == args[2]),aes(x=value,y=eval(parse(text = c.type))),
             col="#ffa07a",size=0.9,pch=15)+
  geom_smooth(data=subset(melted_result,variable == args[1]),aes(x=value,y=eval(parse(text = c.type))),
               col="#ff2400",formula = y ~ x,se=FALSE)+
  geom_smooth(data=subset(melted_result,variable == args[2]),aes(x=value,y=eval(parse(text = c.type))),
               col="#ffa07a",formula = y ~ x,se=FALSE)+
  geom_smooth(data=subset(melted_result,variable == "mixed.score"),aes(x=value,y=eval(parse(text = c.type))),
              col="#8B0000",formula = y ~ x, method = "loess",fill="#FA8072")+
  theme_bw()+
  xlab('LINCS Score')+
  ylab('TCGA Score')+
  theme(legend.position = "None",axis.title = element_text(size=15),axis.text = element_text(size=10))+
  geom_vline(xintercept = 0,linetype='dashed',col='black')+
  geom_hline(yintercept = 0,linetype='dashed',col='black')+
  ggtitle(paste0(c.type,":",args[1],"-",args[2]))

print (plt)

invisible(dev.off())

print(paste0(c.type," and ",args[1],"-",args[2]," scatter plot is generated"))

print("#####################################################################################")
print("#######################                                       #######################")
print("#######################          Plotting scatter plot        #######################")
print("#######################          with Synergy-scores by drug  #######################")
print("#######################                                       #######################")
print("#####################################################################################")

if(!require(ggplot2)){BiocManager::install("ggplot2")};library(ggplot2)
if(!require(reshape2)){BiocManager::install("reshape2")};library(reshape2)

#rm(list=ls())

c.type=args[3]
result_matrix = cbind.data.frame(tcga_score[,c.type,drop=FALSE],mixed.score,lincs_score)
melted_result = melt(result_matrix, id.vars=c.type)

print(args[1])
print(args[2])

png(filename = paste0(argv[4],c.type,'_',args[1],"-",args[2],'.png'),width = 750, height = 750, type = "cairo")


#"lm", "glm", "gam", "loess"
plt=ggplot(melted_result[melted_result$variable != "mixed.score",])+
  geom_point(aes(x=value,y=eval(parse(text = c.type)),col=variable),col = "#999999",size=0.1)+
  geom_point(data = subset(melted_result,variable == args[1]),aes(x=value,y=eval(parse(text = c.type))),
             col="#ff2400",size=0.9,pch=15)+
  geom_point(data = subset(melted_result,variable == args[2]),aes(x=value,y=eval(parse(text = c.type))),
             col="#ffa07a",size=0.9,pch=15)+
  geom_smooth(data=subset(melted_result,variable == args[1]),aes(x=value,y=eval(parse(text = c.type))),
              col="#ff2400",formula = y ~ x,se=FALSE)+
  geom_smooth(data=subset(melted_result,variable == args[2]),aes(x=value,y=eval(parse(text = c.type))),
              col="#ffa07a",formula = y ~ x,se=FALSE)+
  geom_smooth(data=subset(melted_result,variable == "mixed.score"),aes(x=value,y=eval(parse(text = c.type))),
              col="#8B0000",formula = y ~ x, method = "loess",fill="#FA8072")+
  theme_bw()+
  xlab('L-TF Score')+
  ylab('T-TF Score')+
  geom_vline(xintercept = 0,linetype='dashed',col='black')+
  geom_hline(yintercept = 0,linetype='dashed',col='black')+
  labs(title=paste0(args[1],"-",args[2]),
       subtitle =paste0("Dual Drug Score = ",dual.score,"\nSynergy Score = ",syn.score))+
  scale_fill_manual(values = c("#ff2400", "#ffa07a","#8B0000"),
                    labels = c("Tiger", "Lion","Tiger-Lion"))


print (plt)

invisible(dev.off())

print(paste0(c.type," and ",args[1],"-",args[2]," scatter plot is generated"))