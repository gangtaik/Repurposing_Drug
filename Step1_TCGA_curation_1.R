getwd()
file.loc="/Users/hotcoldbrew/Desktop/OneDrive/Research/TCGA-LINCS/TCGA/"
type.list=dir(path=file.loc);type.list

for ( c.type in type.list){
  print (c.type)
  setwd(paste0(file.loc,c.type))
  #install.packages("readxl") 
  #install.packages("dplyr")
  #c.type="DLBC"
  library(readxl)
  library(dplyr)
  id.set=read_excel("../../script/TCGA_mapping.xlsx",sheet = 1,col_names = T)
  #check.names =FALSE => prevent to change header name automatically
  o.df=read.table("./HiSeqV2.gz",sep="\t",header =T,stringsAsFactors = FALSE, check.names = FALSE)
  o.df=dplyr::rename(o.df,"Original"="sample")
  id.set=id.set[,c(2:4)]
  
  #Step.1 (HiSeqV2_curated.tsv)
  ##Add Ensemble, Symbol Id to HiseqV2.gz 
  print("##Add Ensemble, Symbol Id to HiseqV2.gz")
  library(tibble)
  a.df=left_join(id.set,o.df,by="Original")
  t.a.df=data.frame(t(a.df))
  t.a.df=add_column(t.a.df,X0=rownames(t.a.df),.before=1)
  ##column annotation (T_01, N_10, N_11, M_06)
  ##Tumor : 01, 03, 05, 09
  ##Recurrent : 02, 04, 40
  ##Metastatic : 06, 07
  ##Normal : 10, 11, 12 ,14
  print("##column annotation")
  smn.code=colnames(a.df[,4:ncol(a.df)])
  state=data.frame(do.call('rbind',
                           strsplit(as.character(smn.code),split='-',fixed = TRUE)))[,4]
  ann.col=data.frame(cbind(smn.code,state));ann.col
  ann.col$annoate='#'
  
  T.lis=c("01","03","05","09")
  R.lis=c("02","04","40") #259
  M.lis=c("06","07") #161
  N.lis=c("10","11","12","14")
  
  for ( i in 1:nrow(ann.col)){
    if(ann.col$state[i] %in% T.lis){
      ann.col$annoate[i]=paste0("T_",ann.col$state[i])
    }else if(ann.col$state[i] %in% R.lis){
      ann.col$annoate[i]=paste0("R_",ann.col$state[i])
    }else if(ann.col$state[i] %in% M.lis){
      ann.col$annoate[i]=paste0("M_",ann.col$state[i])
    }else{
      ann.col$annoate[i]=paste0("N_",ann.col$state[i])
    }
  }
  emty=data.frame(matrix(rep('#',9),nrow = 3))
  colnames(emty)=colnames(ann.col)
  t.a.df=add_column(t.a.df,annotate=rbind(emty,ann.col)[,3],.before = 1)
  rownames(t.a.df)=1:nrow(t.a.df)
  
  ##column sorting by 1)column annotation and, 2) column name alphabetically 
  ##Generally '#' was used to ignore header
  print("##column sorting")
  t.a.df.1=t.a.df[1:3,]
  t.a.df.2=t.a.df[4:nrow(t.a.df),]
  t.a.df.2$annotate = factor(t.a.df.2$annotate, 
                             levels=c("T_01","T_03","T_05","T_09",
                                      "N_10","N_11","N_12","N_14",
                                      "R_02","R_04","R_40",
                                      "M_06","M_07"))
  
  t.a.df.2=arrange(t.a.df.2,annotate, X0)
  cur.df.1=data.frame(t(rbind(t.a.df.1,t.a.df.2)))
  rm (cur.df)

  file.n=paste0("HiSeqV2_",c.type,"_curated_1.tsv")
  print (file.n)
  write.table(cur.df.1,file=file.n,sep="\t",quote=FALSE,col.names = FALSE,
              row.names = FALSE)
}
print ("Finished")
