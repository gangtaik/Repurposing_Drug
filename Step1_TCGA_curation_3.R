getwd()
file.loc="/Users/hotcoldbrew/Desktop/OneDrive/Research/TCGA-LINCS/TCGA/"
type.list=dir(path=file.loc)
type.list=type.list[type.list!="Ref"]
output.tbl=c()

for (c.type in type.list){
  #c.type="COAD"
  library(stringr)
  library(dplyr)
  file2.n=paste0("HiSeqV2_",c.type,"_curated_2.tsv")
  print (c.type)
  setwd(paste0(file.loc,c.type))
  
  cur.df.2=read.table(file2.n,stringsAsFactors = FALSE,header=FALSE,sep="\t",
                      check.names = FALSE, comment.char="")
  colnames(cur.df.2)=paste0("V",seq(1,ncol(cur.df.2)))
  rownames(cur.df.2)=seq(1,nrow(cur.df.2))
  
  gene.id=cur.df.2[,c(1:2)]
  T.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^T")])
  N.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^N")])
  R.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^R")])
  M.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^M")])
  
  file3.n=paste0("TCGA_RNA-seq_",c.type)
  
  if (ncol(T.df) >=3){
    write.table(T.df,file=paste0(file3.n,"_T.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
    T.df.dim=paste0(nrow(T.df)-2," x ",ncol(T.df)-2," [",2,",",2,"]")
  }else{
    cat (c.type,"HiSeqV2.gz do not include Tumor sample\n",sep=" ")
    T.df.dim=""
  }
  if (ncol(N.df) >=3){
    write.table(N.df,file=paste0(file3.n,"_N.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
    N.df.dim=paste0(nrow(N.df)-2," x ",ncol(N.df)-2," [",2,",",2,"]")
  }else{
    cat (c.type,"HiSeqV2.gz do not include Normal sample\n",sep=" ")
    N.df.dim=""
  }
  if (ncol(R.df) >=3){
    write.table(R.df,file=paste0(file3.n,"_R.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
    R.df.dim=paste0(nrow(R.df)-2," x ",ncol(R.df)-2," [",2,",",2,"]")
  }else{
    cat (c.type,"HiSeqV2.gz do not include Recurrent sample\n",sep=" ")
    R.df.dim=""
  }
  if (ncol(M.df) >=3){
    write.table(M.df,file=paste0(file3.n,"_M.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
    M.df.dim=paste0(nrow(M.df)-2," x ",ncol(M.df)-2," [",2,",",2,"]")
  }else{
    cat (c.type,"HiSeqV2.gz do not include Metastated sample\n",sep=" ")
    M.df.dim=""
  }
  ##Chekc output results
  file.n=paste0("HiSeqV2_",c.type,"_curated_1.tsv")
  o.df=read.table("./HiSeqV2.gz",sep="\t",header =FALSE,stringsAsFactors = FALSE, check.names = FALSE)
  cur.df.1=read.table(file.n,stringsAsFactors = FALSE,header=FALSE,sep="\t",
                      check.names = FALSE, comment.char="")
  colnames(cur.df.1)=paste0("V",seq(1,ncol(cur.df.1)))
  rownames(cur.df.1)=seq(1,nrow(cur.df.1))
  
  
  
  #paste0(nrow(o.df)-1," x ",ncol(o.df)-1," [",1,",",1,"]")
  tbl=data.frame(paste0(nrow(o.df)-1," x ",ncol(o.df)-1," [",1,",",1,"]"),
                 paste0(nrow(cur.df.1)-2," x ",ncol(cur.df.1)-3," [",3,",",2,"]"),
                 paste0(nrow(cur.df.2)-2," x ",ncol(cur.df.2)-2," [",2,",",2,"]"),
                 T.df.dim,N.df.dim,R.df.dim,M.df.dim)
  colnames(tbl)=c("HiSeqV2.gz","HiSeqV2_curated_1.tsv","HiSeqV2_curated_2.tsv","T.tsv","N.tsv",
                  "R.tsv","M.tsv")
  output.tbl=rbind(output.tbl,tbl)
  output.tbl
  
  
}
output.tbl=cbind(type.list,output.tbl)
colnames(output.tbl)=c("","HiSeqV2.gz","HiSeqV2_curated_1.tsv","HiSeqV2_curated_2.tsv","T.tsv","N.tsv",
                       "R.tsv","M.tsv")
write.table(output.tbl,file="../../script/ouput_check.tsv",sep="\t",quote=FALSE,col.names = TRUE,
            row.names = FALSE)
print ("Step1-D is Finished")

