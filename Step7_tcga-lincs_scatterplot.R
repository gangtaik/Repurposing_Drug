print("#####################################################################################")
print("#######################                                       #######################")
print("#######################          Plotting scatter plot        #######################")
print("#######################          with enrichment files        #######################")
print("#######################                                       #######################")
print("#####################################################################################")

# argv=c("/data/MRC1_data9/gangtl95/LINCS/LINCS.Score/",
#        "/data/MRC1_data9/gangtl95/LINCS/TCGA.Score/",
#        "/data/MRC1_data9/gangtl95/LINCS/script/",
#        "/data/MRC1_data9/gangtl95/LINCS/Drug.plot/")

argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/Sinlge.plot/")

if(!require(ggplot2)){BiocManager::install("ggplot2")};library(ggplot2)
if(!require(reshape2)){BiocManager::install("reshape2")};library(reshape2)

#rm(list=ls())

lincs_score = read.table(paste0(argv[1],"_total_lincs_score.tsv"),
                         sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)
tcga_score = read.table(paste0(argv[2],"_total_tcga_score.tsv"),
                        sep='\t',header=T,stringsAsFactors=F,quote='',check.names=F)

tf.list=read.table(paste0(argv[3],"TF_list.txt"),
                   sep="\t",stringsAsFactors = FALSE,header=TRUE)
rowName = lincs_score$Symbol
colnames(lincs_score)
lincs_score = lincs_score[,3:ncol(lincs_score), drop = FALSE]
tcga_score = tcga_score[,3:ncol(tcga_score), drop = FALSE]
drugname_set = colnames(lincs_score)
cancername_set = colnames(tcga_score)

c.type="BRCA"
result_matrix = cbind.data.frame(tcga_score[,c.type,drop=FALSE],lincs_score)
melted_result = melt(result_matrix, id.vars=c.type)
print(c.type)

for(cmpd in drugname_set){
  
  #cmpd = "ZSTK-474"
  #cmpd="BGT-226"
  print(cmpd)
  
  png(filename = paste0(argv[4],c.type,'_',cmpd,'.png'),width = 750, height = 750, type = "cairo")
  
  #"lm", "glm", "gam", "loess"
  plt = ggplot(melted_result)+
    geom_point(aes(x=value,y=eval(parse(text = c.type)),col=variable),col = "#999999",size=0.1)+
    geom_point(data = subset(melted_result,variable == cmpd),aes(x=value,y=eval(parse(text = c.type))),
               col="#ff2400",size=0.9,pch=15)+
    geom_smooth(data=subset(melted_result,variable == cmpd),aes(x=value,y=eval(parse(text = c.type))),
                col="#ff2400",formula = y ~ x, method = "loess",fill="#FFA07A")+
    theme_bw()+
    xlab('LINCS Score')+
    ylab('TCGA Score')+
    theme(legend.position = "None",axis.title = element_text(size=15),axis.text = element_text(size=10))+
    geom_vline(xintercept = 0,linetype='dashed',col='black')+
    geom_hline(yintercept = 0,linetype='dashed',col='black')+
    ggtitle(paste0(c.type,":",cmpd))
  
  print (plt)
  
  invisible(dev.off())
  
  print(paste0(c.type," and ",cmpd," scatter plot is generated"))
  
  
}



print("#####################################################################################")
print("#######################                                       #######################")
print("#######################          Calculating Drug-score       #######################")
print("#######################                                       #######################")
print("#####################################################################################")

sum_total = NULL
final_result = NULL


for(i in 1:ncol(tcga_score)){
  
  for(j in 1:ncol(lincs_score)){ 
    
    print(colnames(tcga_score)[i])
    print(colnames(lincs_score)[j])
  
    #TF.score=lincs_score[,colnames(tcga_score)[i]] * tcga_score[,colnames(lincs_score)[j]]
    TF.score = lincs_score[,j] * tcga_score[,i]
    Drug.score = (-1) * {sum(TF.score, na.rm = FALSE) / nrow(tf.list)}
    #print(TF.Score)
    sum_total = append(sum_total,Drug.score)
    
  }
  
  final_result = rbind.data.frame(final_result,sum_total)
  sum_total <- NULL
}

dim(final_result)
final_result = cbind(colnames(tcga_score),final_result)
final_result = rbind(c("#",colnames(lincs_score)),final_result)

write.table(final_result,paste0(argv[3],"_total_drug-score_results.tsv"),sep='\t',col.names = F, row.names=F,quote=F)

