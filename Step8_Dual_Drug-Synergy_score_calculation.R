# argv=c("/data/MRC1_data9/gangtl95/LINCS.Score/",
#        "/data/MRC1_data9/gangtl95/TCGA.Score/",
#        "/data/MRC1_data9/gangtl95/script/",
#        "/data/MRC1_data9/gangtl95/Comb.plot/",
#        "/data/MRC1_data9/gangtl95/TCGA-LINCS/",
#        "/data/MRC1_data9/gangtl95/Synergy.score/",
#        "/data/MRC1_data9/gangtl95/script/RDA/")

argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/Dual.plot/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/Synergy.score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/RDA/")


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

print("#####################################################################################")
print("#######################                                        #######################")
print("#######################   Calculating total Dual Drug Score    #######################")
print("#######################                                        #######################")
print("#####################################################################################")
#rm(list=ls())
tot.combset=c()
tcga.score=tcga_score[,"BRCA"]
for ( i in colnames(lincs_score)){
  drug.combset=c()
  for ( j in colnames(lincs_score)){
    lincs.score.1=lincs_score[,i]
    lincs.score.2=lincs_score[,j]
    mixed.score=c()
    for ( k in 1:nrow(lincs_score)){
      if(abs(lincs.score.1[k]) >= abs(lincs.score.2[k])){
        mixed.score[k]=lincs.score.1[k]
      }else{
        mixed.score[k]=lincs.score.2[k]
      }
    }
    dd.score=(-1)*(sum(mixed.score * tcga.score)/nrow(tf.list))
    print (paste0(i," /// ",j," : ", dd.score))
    drug.combset=append(drug.combset,dd.score)
  }
  tot.combset=cbind(tot.combset, as.matrix(drug.combset));dim(tot.combset)
}
tot.combset=data.frame(tot.combset)
tot.combset=cbind(colnames(lincs_score),tot.combset)
colnames(tot.combset)=c("#",colnames(lincs_score))
#save(tot.combset,file="~/Desktop/OneDrive/Research/TCGA-LINCS/script/RDA/tot.combset.rda")
write.table(tot.combset, paste0(argv[5],"_total_Dual_Drug_Score.tsv"),sep="\t",quote=FALSE,col.names = TRUE, row.names = FALSE)

print("#####################################################################################")
print("#######################                                       #######################")
print("#######################     Calculating total Synergy-score   #######################")
print("#######################                                       #######################")
print("#####################################################################################")

tot.drugset=read.table(paste0(argv[3],"Drug-score_results.tsv"),
                       sep='\t',header=T,stringsAsFactors=F,quote='',
                       check.names=F,comment.char="")
rownames(tot.drugset)=tot.drugset[,1]
tot.drugset=tot.drugset[,-1]


load(paste0(argv[7],"tot.combset.rda"));dim(tot.combset)
rownames(tot.combset)=tot.combset[,1]
tot.combset=tot.combset[,-1]

tot.syn=c()
for ( i in colnames(tot.combset)){
  #i="10-DEBC"
  #j="1271738-62-5"
  single.score=c()
  for ( j in colnames(tot.combset)){
    print (paste0(i," /// ",j))
    single.score = rbind(single.score,
                         c(tot.drugset["BRCA",i],
                           tot.drugset["BRCA",j],
                           tot.combset[i,j],
                           tot.combset[i,j]/mean(tot.drugset["BRCA",i],tot.drugset["BRCA",j])))
  }
  
  single.save=single.score
  single.save=cbind(colnames(tot.combset),single.save)
  colnames(single.save)=c("#drug.B",paste0(i,".drug.score"), "drug.B.drug.score", "combination.score","synergy.score")
  write.table(single.save,file=paste0(argv[6],i,"_synergey.score.tsv"),sep="\t",quote=FALSE,col.names = TRUE, row.names = FALSE)

  mrg.col=data.frame(apply(single.score, 1, paste ,collapse = " , " ))
  colnames(mrg.col)=i
  tot.syn=cbind(tot.syn,as.matrix(mrg.col))
}
tot.syn=cbind(colnames(tot.combset),tot.syn)
tot.syn=rbind.data.frame(c("[Colum.Drug.score, Row.Drug.score, Combination.score, Synergy.score]",
                           colnames(tot.combset)),tot.syn)
save(tot.syn,file="~/Desktop/OneDrive/Research/TCGA-LINCS/script/RDA/tot.syn.rda")
#load(paste0(argv[7],"tot.syn.rda"))

write.table(tot.syn, paste0(argv[5],"_total_Synergy_Score.tsv"),sep="\t",quote=FALSE,col.names = FALSE, row.names = FALSE)