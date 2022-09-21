#rm(list=ls())
#argv[1]=input.file.loc
#argv[2]=ouput.fil.loc
#argv[3]=maxdose_cancer.fil.loc

# argv=c("/data/MRC1_data9/gangtl95/LINCS/LINCS.GSVA/",
#        "/data/MRC1_data9/gangtl95/LINCS/LINCS.Score/",
#        "/data/MRC1_data9/gangtl95/LINCS/LINCS_chemical_cancer_maxdose/",
#        "/data/MRC1_data9/gangtl95/LINCS/script/")


argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS.GSVA/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS.Score/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer_maxdose/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/script/")

#Calculating TF Scores by compounds
print ("Calculating TF scores by compounds")
outlier = NULL
cp.list = dir(path=argv[3],pattern="_max_cancer.tsv")
cp.list = cp.list[cp.list!="_DMSO_max_cancer.tsv"];length(cp.list)
gene.id=read.table(paste0(argv[4],"TF_list.txt"),sep="\t",stringsAsFactors = FALSE,
                   header=TRUE)
for ( cp in cp.list ){
  if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
  if(!require(stringr)){install.packages("stringr")};library(stringr)

  #GSVA.score table
  input.df = read.table(paste0(argv[1],"LINCS_maxdose_GSVA_z-score.tsv"),sep='\t',header=T,
                        stringsAsFactors=F,quote='',check.names=FALSE);dim(input.df)
  
  #setdiff(input.df$Symbol,gene.id$Symbol)
  input.tbl=input.df[input.df$Symbol %in% gene.id$Symbol,]
  input.tbl = input.tbl[,3:ncol(input.tbl)]
  dim(input.tbl)
  # Extract DMSO(Ctl) samples
  col.dbso=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = "^ctl_DMSO")];length(col.dbso)
  col.drug=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = "^trt_")];length(col.drug)
  
  ctrl.dmso = input.tbl[,colnames(input.tbl) %in% col.dbso];dim(ctrl.dmso)
  trt.drug=input.tbl[,colnames(input.tbl) %in% col.drug];dim(trt.drug)
  
  
  #cp=cp.list[length(cp.list)]
  #cp="16,16-dimethylprostaglandin-e2_max_cancer.tsv"
  #cp="2-,5--dideoxyadenosine_max_cancer.tsv"
  drug=gsub(pattern = "_max_cancer.tsv",replacement = "",cp)
  
  #Maxdose_cancer table
  trt.exp.tmp = read.table(paste0(argv[3],drug,"_max_cancer.tsv"),sep='\t',header=T,
                           stringsAsFactors=F,quote='',check.names=FALSE);dim(trt.exp.tmp)
  trt.exp.tmp = trt.exp.tmp[order(trt.exp.tmp$Symbol),]
  trt.exp2.tmp = trt.exp.tmp[trt.exp.tmp$Symbol %in% gene.id$Symbol, ];nrow(trt.exp2.tmp)
  
  trt.exp = data.frame(trt.exp2.tmp[,-c(1:2)])
  rownames(trt.exp) = trt.exp2.tmp$Symbol
  trt.exp=mutate_all(trt.exp, function(x) as.numeric(as.character(x)))
  

  
  ctrl.exp.tmp =  read.table(paste0(argv[3],"_DMSO_max_cancer.tsv"),sep='\t',header=T,
                             stringsAsFactors=F,quote='',check.names=FALSE);dim(ctrl.exp.tmp)
  ctrl.exp.tmp = ctrl.exp.tmp[order(ctrl.exp.tmp$Symbol),]
  ctrl.exp2.tmp = ctrl.exp.tmp[ctrl.exp.tmp$Symbol %in% gene.id$Symbol, ];nrow(ctrl.exp2.tmp)
  
  ctrl.exp = data.frame(ctrl.exp2.tmp[,-c(1:2)])
  rownames(ctrl.exp) = ctrl.exp2.tmp$Symbol
  ctrl.exp=mutate_all(ctrl.exp, function(x) as.numeric(as.character(x)))
  
  print("Extract treatment group sample for each drug")
  if (drug == "atenolol--+---"){
    print(drug)
    drug.col=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = paste0("^trt_atenolol"))];length(drug.col)
  }else if(drug =="butorphanol--+--tartrate"){
    print(drug)
    drug.col=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = paste0("^trt_butorphanol"))];length(drug.col)
  }else if(drug =="RHO-kinase-inhibitor-III[rockout]"){
    print(drug)
    drug.col=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = paste0("^trt_RHO-kinase"))];length(drug.col)
  }else if(drug =="SK-F-10047--+-"){
    print(drug)
    drug.col=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = paste0("^trt_RHO-kinase"))];length(drug.col)
  }else{
    print (drug)
    drug.col=colnames(input.tbl)[str_detect(string = colnames(input.tbl),pattern = paste0("^trt_",drug))];length(drug.col)
  }
  
  # Extract treatment group sample for each drug
  
  drug.tbl=input.tbl[,colnames(input.tbl) %in% drug.col];dim(drug.tbl)
  
  
  if(ncol(drug.tbl) < 2 ){    
    print(paste(drug,"data has small set.",sep=" "))
    
    outlier = append(outlier,drug)
    
    next
  }
  
  
  # Do the student t.test 
  lincs.score=data.frame(gene.id)
  #trt.exp = treatment expression
  #ctrl.exp = dmso expression
  #drug.tbl = treatment GSVA calculating expression
  #ctrl.dmso = dmso GSVA calculating expression
  
  for(i in 1:nrow(lincs.score)){
    
    lincs.score$exp.t.val[i]=t.test(trt.exp[i,],ctrl.exp[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$statistic
    lincs.score$gsva.t.val[i]=t.test(drug.tbl[i,],ctrl.dmso[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$statistic
    lincs.score$exp.p.val[i]=t.test(trt.exp[i,],ctrl.exp[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$p.value
    lincs.score$gsva.p.val[i]=t.test(drug.tbl[i,],ctrl.dmso[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$p.value
    
    if(t.test(trt.exp[i,],ctrl.exp[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$stderr == 0){
      lincs.score$exp.t.val[i]=0
    }
    if(t.test(drug.tbl[i,],ctrl.dmso[i,],alternative = "two.sided",paired = FALSE,var.equal = FALSE)$stderr == 0){
      lincs.score$gsva.t.val[i]=0
    }
  }
  lincs.score$TF.Score=lincs.score$exp.t.val+lincs.score$gsva.t.val
  
  write.table(lincs.score,paste0(argv[2],"LINCS_",drug,"_lincs_score.tsv"),sep='\t',quote=FALSE,col.names = TRUE, row.names = FALSE)
  list=ls()
  list=list[list!="argv"]
  list=list[list!="outlier"]
  list=list[list!="cp.list"]
  list=list[list!="gene.id"]
  rm(list)
}
write.table(outlier,paste0(argv[2],"lincs_score_outlier.txt"),sep='\t',quote = F)


#Total ES score for compounds which have enough treated sample number ( > 2)
lincs.es.list=dir(path=argv[2],pattern = "_lincs_score.tsv");length(lincs.es.list)
tot.es=c()
lincs.es.list=lincs.es.list[lincs.es.list!="_total_lincs_score.tsv"]
setwd(argv[2])

for ( tsv in lincs.es.list){
  #tsv="LINCS_ZSTK-474_lincs_score.tsv"
  df=read.table(tsv, header=TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep="\t",comment.char = "")
  tot.es=cbind(tot.es, df$TF.Score)
}

es.colnm=gsub("LINCS_",replacement = "", x=lincs.es.list)
es.colnm=gsub("_lincs_score.tsv",replacement = "",x=es.colnm)
colnames(tot.es)=es.colnm

fin.es=cbind(gene.id,tot.es);nrow(fin.es)

write.table(fin.es,file=paste0(argv[2],"_total_lincs_score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
            row.names = FALSE)





