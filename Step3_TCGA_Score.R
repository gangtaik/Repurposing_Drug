#argv[1]=working dir
#argv[2]=result dir
#argv[3]=TF.idex.file
#argv[5]=Enrichment Score results which number of c.type normal  is more than 10 or none
#argv[6]=TF.list.file.loc
#rm(list=ls())

# argv=c("/data/MRC1_data9/gangtl95/TCGA/",
#        "/data/MRC1_data9/gangtl95/TCGA.GSVA/",
#        "/data/MRC1_data9/gangtl95/script/tft_benchmark_symbol_10TG.tsv",
#        "/data/MRC1_data9/gangtl95/TCGA.Score/",
#        "/data/MRC1_data9/gangtl95/TCGA.Score/Less10/",
#        "/data/MRC1_data9/gangtl95/script/")
argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA.GSVA/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/tft_benchmark_symbol_10TG.tsv",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/TCGA.Score/Less10/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/")

type.list=dir(path=argv[1])
type.list=type.list[type.list!="Ref"]


for (c.type in type.list){
  #c.type="COAD"
  if(!require(stringr)){install.packages("stringr")};library(stringr)
  if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
  
  print (c.type)
  
  print ("TCGA Score is being calculated")
  #TF score = expression matrix.t.value + gsva.t.value
  setwd(paste0(argv[1],c.type))
  ##column annotation (T_01, N_10, N_11, M_06)
  ##Tumor : 01, 03, 05, 09
  ##Recurrent : 02, 04, 40
  ##Metastatic : 06, 07
  ##Normal : 10, 11, 12 ,14
  gsva.scale=read.table(paste0(argv[2],"TCGA_",c.type,"_GSVA_z-score.tsv"),
                        sep="\t", check.names = FALSE, comment.char = "",stringsAsFactors = FALSE,header=TRUE)
  
  gsva.scale.T=data.frame(gsva.scale[,c(1:2)],gsva.scale[,grep("0[1-9]$",colnames(gsva.scale))])
  gsva.scale.N=data.frame(gsva.scale[,c(1:2)],gsva.scale[,grep("1[1-9]$",colnames(gsva.scale))])
  
  tf.list=read.table(paste0(argv[6],"TF_list.txt"),sep="\t",stringsAsFactors = FALSE,
                     header=TRUE)
  
  gsva.scale.T=gsva.scale.T[gsva.scale.T$Symbol %in% tf.list$Symbol,]
  gsva.scale.N=gsva.scale.N[gsva.scale.N$Symbol %in% tf.list$Symbol,]
  
  #expression files
  df.T=read.table(paste0(argv[1],c.type,"/TCGA_RNA-seq_",c.type,"_T.tsv"),
                  sep="\t", check.names = FALSE, comment.char = "",stringsAsFactors = FALSE,header=TRUE)
  tmp.T=df.T[-1,-c(1:2)]
  colnames(tmp.T)=df.T[1,-c(1:2)]
  rownames(tmp.T)=df.T[-1,2]
  
  #No value of FOXO6 #deletion
  exp.T=tmp.T[rownames(tmp.T) %in% gsva.scale.T$Symbol,]
  exp.T=exp.T[order(rownames(exp.T)),]
  exp.T=mutate_all(exp.T, function(x) as.numeric(as.character(x)))
  
  exp.list=dir(path=paste0(argv[1],c.type),pattern = ".tsv")
  if(paste0("TCGA_RNA-seq_",c.type,"_N.tsv") %in% exp.list){
    
    df.N=read.table(paste0(argv[1],c.type,"/TCGA_RNA-seq_",c.type,"_N.tsv"),
                    sep="\t", check.names = FALSE, comment.char = "",stringsAsFactors = FALSE,header=TRUE)
    tmp.N=data.frame(df.N[-1,-c(1:2)])
    colnames(tmp.N)=df.N[1,-c(1:2)]
    tmp.N=cbind(Symbol=df.N[-1,2],tmp.N)
    tmp.N=tmp.N[order(tmp.N$Symbol),]
    
    tmp2.N=data.frame(tmp.N[tmp.N$Symbol %in% gsva.scale.T$Symbol,])
    rownames(tmp2.N)=tmp2.N$Symbol
    exp.N=data.frame(tmp2.N[,-1])
    colnames(exp.N)=colnames(tmp2.N)[-1]
    rownames(exp.N)=rownames(tmp2.N)
    exp.N=mutate_all(exp.N, function(x) as.numeric(as.character(x)))
    
  }else{
    df.N=tf.list
    exp.N=data.frame(value=rep(NA,nrow(tf.list)))
    rownames(exp.N)=tf.list$Symbol
  }
  
  if(ncol(gsva.scale.N)<=3){
    print("Normal samples are not existed or only one normal existed")
    
    tcga.score=data.frame(tf.list,
                          exp.t.val=rep(NA,nrow(tf.list)),
                          gsva.t.val=rep(NA,nrow(tf.list)),
                          exp.p.val=rep(NA,nrow(tf.list)),
                          gsva.p.val=rep(NA,nrow(tf.list)),
                          TF.Score=rep(NA,nrow(tf.list)))
    write.table(tcga.score,file=paste0(argv[5],"TCGA_",c.type,"_tcga_score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
                row.names = FALSE)
    
  }else if(ncol(gsva.scale.N)>3 && ncol(gsva.scale.N)<=10){
    print(c.type)
    print("Normal samples are not enough but will include")
    tcga.score=data.frame(tf.list)
    
    for ( i in 1:nrow(tcga.score)){
      
      #i=227
      tcga.score$exp.t.val[i]=t.test(exp.T[i,],exp.N[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$statistic
      tcga.score$gsva.t.val[i]=t.test(gsva.scale.T[i,-c(1:2)],gsva.scale.N[i,-c(1:2)],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$statistic
      tcga.score$exp.p.val[i]=t.test(exp.T[i,],exp.N[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$p.value
      tcga.score$gsva.p.val[i]=t.test(gsva.scale.T[i,-c(1:2)],gsva.scale.N[i,-c(1:2)],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$p.value
      
      if(t.test(exp.T[i,],exp.N[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$stderr == 0){
        tcga.score$exp.t.val[i]=0
      }
      if(t.test(gsva.scale.T[i,-c(1:2)],gsva.scale.N[i,-c(1:2)],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$stderr == 0){
        tcga.score$gsva.t.val[i]=0
      }
    }
    tcga.score$TF.Score=tcga.score$exp.t.val+tcga.score$gsva.t.val
    
    write.table(tcga.score,file=paste0(argv[4],"TCGA_",c.type,"_tcga_score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
                row.names = FALSE)
    write.table(tcga.score,file=paste0(argv[5],"TCGA_",c.type,"_tcga_score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
                row.names = FALSE)
    
  }else{
    print("Normal samples are enough")
    tcga.score=data.frame(tf.list)
    
    for ( i in 1:nrow(tcga.score)){
      tcga.score$exp.t.val[i]=t.test(exp.T[i,],exp.N[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$statistic
      tcga.score$gsva.t.val[i]=t.test(gsva.scale.T[i,-c(1:2)],gsva.scale.N[i,-c(1:2)],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$statistic
      tcga.score$exp.p.val[i]=t.test(exp.T[i,],exp.N[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$p.value
      tcga.score$gsva.p.val[i]=t.test(gsva.scale.T[i,-c(1:2)],gsva.scale.N[i,-c(1:2)],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$p.value
      
      if(t.test(exp.T[i,],exp.N[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$stderr == 0){
        tcga.score$exp.t.val[i]=0
      }
      if(t.test(gsva.scale.T[i,-c(1:2)],gsva.scale.N[i,-c(1:2)],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$stderr == 0){
        tcga.score$gsva.t.val[i]=0
      }
    }
    tcga.score$TF.Score=tcga.score$exp.t.val+tcga.score$gsva.t.val
    write.table(tcga.score,file=paste0(argv[4],"TCGA_",c.type,"_tcga_score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
                row.names = FALSE)
  }
  list=ls()
  list=list[list!="type.list"]
  list=list[list!="argv"]
  rm(list)
}
list=ls()
list=list[list!="type.list"]
list=list[list!="argv"]
rm(list)

#Total ES score for c.type which have enough normal sample list
es.list=dir(path=argv[4],pattern="_tcga_score.tsv")
es.list=es.list[es.list!="_total_tcga_score.tsv"]
tot.es=c()
setwd(argv[4])
for ( tsv in es.list){
  df=read.table(tsv,header=TRUE,stringsAsFactors = FALSE, check.names = FALSE, sep="\t",comment.char = "")
  tot.es=cbind(tot.es, df$TF.Score)
}
es.colnm=gsub("TCGA_",replacement = "",x=es.list)
es.colnm=gsub("_tcga_score.tsv",replacement = "",x=es.colnm)
colnames(tot.es)=es.colnm
fin.es=data.frame(df[,c(1:2)],tot.es)
write.table(fin.es,file=paste0(argv[4],"/","_total_tcga_score.tsv"),sep="\t",quote=FALSE,col.names = TRUE, row.names = FALSE)





