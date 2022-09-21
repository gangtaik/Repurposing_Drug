if(!require(stringr)){install.packages("stringr");library(stringr)}
if(!require(GSVA)){BiocManager::install("GSVA");library(GSVA)}
if(!require(dplyr)){install.packages("dplyr");library(dplyr)}

#argv[1]=working dir
#argv[2]=result dir
#argv[3]=TF.idex.file
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

#making Transcription factor interaction index list
print("making Transcription factor interaction index list")
TF.index=read.table(argv[3],header = TRUE, sep='\t',na.strings = "", fill = TRUE)
TF.index.list <- list()
for(i in 1:nrow(TF.index)){
  TF.interact=TF.index[i,c(4:ncol(TF.index))]
  TF.interact=TF.interact[!is.na(TF.interact)]
  TF.index.list[[TF.index[i,1]]]=TF.interact
}

for (c.type in type.list){
  print (c.type)
  #c.type="KIRC"
  setwd(paste0(argv[1],c.type))
  
  file.list=list.files(path="./")
  file.N=file.list[str_detect(file.list,pattern="^TCGA") & str_detect(file.list,pattern="N.tsv$")]
  file.T=file.list[str_detect(file.list,pattern="^TCGA") & str_detect(file.list,pattern="T.tsv$")]
  
  if(length(file.N)!=0){
    expression.file.N=read.table(paste0(argv[1],c.type,"/",file.N),
                                 stringsAsFactors = FALSE,header=FALSE,sep="\t",
                                 check.names = FALSE, comment.char="")
    full.N=expression.file.N[,-c(1:2)]
  }else{
    full.N=c()
  }
  if(length(file.T)!=0){
    expression.file.T=read.table(paste0(argv[1],c.type,"/",file.T),
                                 stringsAsFactors = FALSE,header=FALSE,sep="\t",
                                 check.names = FALSE, comment.char="")
    full.T=expression.file.T[,-c(1:2)]
  }else{
    full.T=c()
  }
  
  if(is.null(full.N)){
    expression.full=full.T
  }else{
    expression.full=cbind(full.N,full.T)
  }
  
  gene.id.col=expression.file.T[,c(1:2)]
  rownames(expression.full)=expression.file.T[,2]
  colnames(expression.full)=expression.full[2,]
  
  gsva.input=apply(expression.full[-c(1:2),],2,as.numeric)
  rownames(gsva.input)=expression.file.T[-c(1:2),2]
  
  print(c.type)
  print("ssgsea is calculating")
  suppressWarnings({ 
    gsva.output <- gsva(gsva.input, gset.idx.list = TF.index.list, method='ssgsea',
                        kcdf="Gaussian", abs.ranking=FALSE, min.sz=2,
                        max.sz=Inf, parallel.sz=20, mx.diff=TRUE,
                        ssgsea.norm=FALSE,
                        verbose=F)
    
    # kcdf=c("Gaussian", "Poisson", "none") 
    #"Gaussian" for which expression values are continuous
    #"Poisson" for which expression values are integer counts
    # parallel.sz => number of threads to execution to use when doing calculations
  })
  
  #gsva.output
  gsva.unscale=cbind(Symbol=rownames(gsva.output),as.data.frame(gsva.output))
  TF.symbol=data.frame(rownames(gsva.output))
  colnames(TF.symbol)="Symbol"
  gene.id=gene.id.col[3:nrow(gene.id.col),]
  colnames(gene.id)=gene.id.col[2,]
  fr.col=left_join(TF.symbol,gene.id,by="Symbol")
  gsva.unscale=cbind(Entrez=fr.col$Entrez,gsva.unscale)
  write.table(gsva.unscale,file=paste0(argv[2],"/","TCGA_",c.type,"_GSVA.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
              row.names = FALSE)
  
  #scaling by z-score, mode as "col"
  print(c.type)
  print("scaling by z-score")
  zscr <- function(x, mode = "col") {
    x <- as.matrix(x)
    z <- x
    if (mode == "row") {
      for (i in 1:nrow(x)) {
        z[i, ] <- (x[i, ] - mean(x[i, ]))/sd(x[i, ])
      }
      return(z)
    }
    else if (mode == "col") {
      for (i in 1:ncol(x)) {
        z[, i] <- (x[, i] - mean(x[, i]))/sd(x[, i])
      }
      return(z)
    }
    else if (mode == "all") {
      z <- scale(x, center = TRUE, scale = TRUE)
      return(z)
    }
    else errorCondition("\"mode\" should be one of (1) \"row\", (2) \"col\", (3) \"all\"")
  }
  gsva.output.scaling=apply(gsva.output,2,zscr)
  gsva.scale=cbind(Symbol=rownames(gsva.output),as.data.frame(gsva.output.scaling))
  gsva.scale=cbind(Entrez=fr.col$Entrez,gsva.scale)
  write.table(gsva.scale,file=paste0(argv[2],"/","TCGA_",c.type,"_GSVA_z-score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
              row.names = FALSE)
}


print ("Step2: GSVA calculation is Finished")
