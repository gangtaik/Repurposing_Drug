getwd()
file.loc="/Users/hotcoldbrew/Desktop/OneDrive/Research/TCGA-LINCS/TCGA/"
type.list=dir(path=file.loc)
type.list=type.list[type.list!="Ref"]

for ( c.type in type.list){
  #c.type="DLBC"
  file.n=paste0("HiSeqV2_",c.type,"_curated_1.tsv")
  print (c.type)
  setwd(paste0(file.loc,c.type))
  library(dplyr)
  cur.df.1=read.table(file.n,stringsAsFactors = FALSE,header=FALSE,sep="\t",
                      check.names = FALSE, comment.char="")
  colnames(cur.df.1)=paste0("V",seq(1,ncol(cur.df.1)))
  rownames(cur.df.1)=seq(1,nrow(cur.df.1))

  #Step.C (HiSeqV2_curated_2.tsv)
  ##(1) Withdrawn gene row deletion
  print("Withdrawn gene row deletion")
  df.3.1=cur.df.1[cur.df.1$V3!="Withdrawn",] #[row, col]
  
  #(2)중복된 entrez ID가 없도록 collapse한 파일(row의 average 발현값이 높은 row를 남기고, 나머지 삭제)
  # 20435 row가 나오면 맞게 계산된 것 (V1=Original, V2=Entrez ID ,V3= Gene Symbol)
  print("Duplicated entrez id row collapse")

  coln=df.3.1[c(1:2),]
  dup.df=df.3.1[which(duplicated(df.3.1$V2)|duplicated(df.3.1$V2,fromLast = TRUE)),]
  dup.df=dup.df[order(dup.df$V2),]
  remain.df=df.3.1[!df.3.1$V2 %in% unique(dup.df$V2),][-c(1:2),]
  uniq.df=c()
  for (i in unique(dup.df$V2)){
    #unique(dup.df$V2)
    #i="100134869"
    #print (i)
    dup.entz=dup.df[dup.df$V2==i , ]
    dup.entz.id=dup.entz[,1:3]
    dup.entz.id$idex=paste0("id_",seq(1,nrow(dup.entz.id)))
    dup.entz.val=mutate_all(dup.entz[,-c(1:3)], function(x) as.numeric(as.character(x)))
    rownames(dup.entz.val)=dup.entz.id$idex
    mx.entz.val=dup.entz.val[which.max(rowMeans(dup.entz.val)),]
    mx.row=cbind(dup.entz.id[dup.entz.id$idex==rownames(mx.entz.val),][,c(1:3)],mx.entz.val)
    uniq.df=rbind(uniq.df,mx.row)
  }
  uniq.df=rbind(remain.df,uniq.df)
  uniq.df=uniq.df[order(uniq.df$V2),]
  uniq.df=rbind(coln,uniq.df)
  ##Original column deletion
  print("Original column deletion")
  cur.df.2=uniq.df[,-1]
  
  file2.n=paste0("HiSeqV2_",c.type,"_curated_2.tsv")
  print (file2.n)
  write.table(cur.df.2,file=file2.n,sep="\t",quote=FALSE,col.names = FALSE,
              row.names = FALSE)
  
}
print ("Step1-C Finished")


