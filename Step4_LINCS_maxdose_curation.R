getwd()
setwd("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer/")

#rm(list=ls())
if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
if(!require(stringr)){install.packages("stringr")};library(stringr)

# argv=c("/data/MRC1_data9/gangtl95/LINCS/LINCS_chemical_cancer/",
#        "/data/MRC1_data9/gangtl95/LINCS/LINCS_chemical_cancer_maxdose/",
#        "/data/MRC1_data9/gangtl95/LINCS/")

argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer_maxdose/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/")

chem_can.list=dir(path = argv[1], pattern=".tsv")
nomaxdose=c()
mixmaxose=c()

n=0
rm.celL=c()
for (cp_can in chem_can.list){
  #cp_can="AZD-5438_cancer.tsv"
  #cp_can=chem_can.list[3]
  #maxdose column separation
  df=read.table(paste0(argv[1],cp_can),sep="\t",check.names = FALSE, na.strings = "", 
                header = TRUE,stringsAsFactors = FALSE);print(cp_can)
  info=colnames(df)[-c(1:2)]
  info.table=data.frame(t(do.call("cbind",str_split(info,pattern = "_"))),info)
  
  if(cp_can=="_DMSO_cancer.tsv"){
    
    colnames(info.table)=c("treat","drug","time.hr","cell.line","id")
    info.table$time.hr=gsub(pattern="h",replacement = "", info.table$time.hr)
    info.table$time.hr=as.numeric(info.table$time.hr)
    info.table=info.table %>% group_by(time.hr)
    max.time=max(info.table$time.hr)
    max.coln=info.table[info.table$time.hr==max.time,]
    max.coln=max.coln[order(max.coln$cell.line),]
    maxdose.tmp=df[,colnames(df) %in% max.coln$id]
    maxdose=cbind(df[,1:2],maxdose.tmp);dim(maxdose)
    
    cp=sub("_cancer.tsv", "_max_cancer.tsv", cp_can)
    write.table(maxdose, paste0(argv[2],cp),sep="\t",quote=FALSE,col.names = TRUE, row.names = FALSE)
    
  }else{
    #View(info.table)
    colnames(info.table)=c("treat","drug","time.hr","concentration.uM","cell.line","id")
    info.table$time.hr=gsub(pattern="h",replacement = "", info.table$time.hr)
    info.table$time.hr=as.numeric(info.table$time.hr)
    info.table$concentration.uM=gsub(pattern="uM",replacement = "", info.table$concentration.uM)
    info.table$concentration.uM=as.numeric(info.table$concentration.uM)
    info.table=info.table %>% group_by(time.hr,concentration.uM,cell.line)
    celL=sort(unique(info.table$cell.line))
    
    max.coln=c()
    for ( cL in celL){
      assign(cL, info.table[info.table$cell.line==cL,])
      max.cct=max(get(cL)$concentration.uM)
      max.time=max(get(cL)$time.hr)
      max.coln.tmp=get(cL)[get(cL)$concentration.uM==max.cct & get(cL)$time.hr==max.time,]
      max.coln=rbind(max.coln,max.coln.tmp)
      
    }
    #cell line 중 maxtime-maxdose data가 혼합되어 나타난 cell line은 제거한 후 저장
    #일부 cell line들이 제거된 drug expression 들 목록
    if(nrow(max.coln)!=length(celL)){
      rm.celL=rbind(rm.celL,cp_can)
      n=n+1
    }

    if(nrow(max.coln)==0){
      print (paste0(gsub(pattern="_cancer.tsv",replacement = "",cp_can)," has no maxdose data"))
      nomaxdose=rbind(nomaxdose,gsub(pattern="_cancer.tsv",replacement = "",cp_can))

    }else if(length(unique(max.coln$time.hr))!=1 & length(unique(max.coln$concentration.uM))!=1 ){
      print (paste0(gsub(pattern="_cancer.tsv",replacement = "",cp_can)," has different maxdose data among cell lines"))
      mixmaxose=rbind(mixmaxose,gsub(pattern="_cancer.tsv",replacement = "",cp_can))

    }else{
      print (paste0(gsub(pattern="_cancer.tsv",replacement = "",cp_can)," has same maxdose data among cell lines"))
      max.coln=max.coln[order(max.coln$cell.line),]
      length(max.coln$id)

      maxdose.tmp=df[,colnames(df) %in% max.coln$id]
      maxdose=cbind(df[,1:2],maxdose.tmp)

      cp=sub("_cancer.tsv", "_max_cancer.tsv", cp_can)
      write.table(maxdose, paste0(argv[2],cp),sep="\t",quote=FALSE,col.names = TRUE, row.names = FALSE)
    }
    
  }
}

print("Maxdose detection is finished")
#cell line 중 maxtime-maxdose data가 혼합되어 나타난 일부 cell line은 제거한 후 저장
#일부 cell line들이 제거된 drug expression 들 목록
write.table(mixmaxose, paste0(argv[3],"rm_whole_list.tsv"),sep="\t",quote=FALSE,col.names = FALSE, row.names = FALSE)

#cell line 중 maxtime-maxdose data가 모든 cell line data가 혼합되었기에 cell line은 제거했기에 _max_cancer.tsv 없음
#모든 cell line들이 제거된 drug expression 들 목록
print(n)
write.table(rm.celL, paste0(argv[3],"rm_partly_list.tsv"),sep="\t",quote=FALSE,col.names = FALSE, row.names = FALSE)

#cell line 중 maxtime-maxdose data가 아예 존재를 안하는 drug 목록
write.table(nomaxdose, paste0(argv[3],"No_maxdose.tsv"),sep="\t",quote=FALSE,col.names = FALSE, row.names = FALSE)

