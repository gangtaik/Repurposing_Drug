getwd()
setwd("~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer/")
if(!require(dplyr)){install.packages("dplyr")};library(dplyr)
if(!require(stringr)){install.packages("stringr")};library(stringr)
if(!require(GSVA)){BiocManager::install("GSVA");library(GSVA)}
if(!require(pheatmap)){install.packages("pheatmap")};library(pheatmap)

#rm(list=ls())
#argv[1]=tft benchmark ssymbol TG file loc
#argv[2]=LINCS chemical cancer data loc
#argv[3]=LINCS chemical cancer maxdose data loc
#argv[4]=LINCS GSVA output results loc

# argv=c("/data/MRC1_data9/gangtl95/script/",
#        "/data/MRC1_data9/gangtl95/LINCS/LINCS_chemical_cancer/",
#        "/data/MRC1_data9/gangtl95/LINCS/LINCS_chemical_cancer_maxdose/",
#        "/data/MRC1_data9/gangtl95/LINCS.GSVA/")

argv=c("~/Desktop/OneDrive/Research/TCGA-LINCS/script/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS/LINCS_chemical_cancer_maxdose/",
       "~/Desktop/OneDrive/Research/TCGA-LINCS/LINCS.GSVA/")


# Transcription Factor index
print("making Transcription factor interaction index list")
TF.index=read.table(paste0(argv[1],"tft_benchmark_symbol_10TG.tsv"),
                    header = TRUE, sep='\t',na.strings = "", fill = TRUE)
TF.index.list <- list()
for(i in 1:nrow(TF.index)){
  TF.interact=TF.index[i,c(4:ncol(TF.index))]
  TF.interact=TF.interact[!is.na(TF.interact)]
  TF.index.list[[TF.index[i,1]]]=TF.interact
}

# Load the gene set and restructure into list

cp.list = dir(path=argv[2],pattern="tsv")
cp.list=cp.list[cp.list!="DMSO_cancer.tsv"];length(cp.list)

# Load DMSO;the control group
df = read.table(paste0(argv[2],"/DMSO_cancer.tsv"),sep='\t',header=T,stringsAsFactors = F,check.names=FALSE);dim(df)
ctrl = df[,3:ncol(df)]
rownames(ctrl) = df[,2];dim(ctrl)


# Remove duplicated columns and leave values with 24 hour result

index = read.table(paste0(argv[1],"index_sorted_descending.txt"),
                   sep='\t',header=T,stringsAsFactors = F,na.strings = "", fill = TRUE)
selected_cell = c('A375', 'A549', 'BT20', 'HCC515', 'HELA', 'HEPG2', 'HS578T', 'HT29', 'JURKAT', 'LNCAP', 'MCF7', 'MDAMB231',  'PC3', 'SKBR3', 'YAPC')
index_sorted = index[!duplicated(index$new_drug_name),]
index_sorted = index_sorted[index_sorted$cell_id %in% selected_cell & index_sorted$time == 24,]
index_sorted$drugname <- gsub('DMSO', '_DMSO', index_sorted$drugname)

ctrl_col <- index_sorted[index_sorted$drugname == "_DMSO" ,"new_drug_name"]
ctrl <- ctrl[,colnames(ctrl) %in% ctrl_col];dim(ctrl)
no_result <- NULL
final_matrix <- data.frame(matrix(NA, nrow = 12326, ncol = 0))

# Make matrix for GSVA function
for(drug_f in cp.list){
  #drug_f="AR-A014418_cancer.tsv"
  drug <- strsplit(drug_f,'_cancer.tsv')[[1]][1]
  print(drug)
  if(drug!='DMSO'){
    drug_mat <- read.table(paste0(argv[2],drug_f),sep='\t',header = T,stringsAsFactors = F,check.names=FALSE)
    drug_mat[drug_mat$Symbol=="ADNP",]
    
    if(ncol(drug_mat)>2){
      drug_mat1 <-as.data.frame(drug_mat[,3:ncol(drug_mat)])
      rownames(drug_mat1) <- rownames(ctrl)
      new_col <- index_sorted[index_sorted$drugname == drug,"new_drug_name"]
      drug_mat1<-drug_mat1[,colnames(drug_mat1) %in% new_col]
      
      if( ncol(drug_mat1) > 1 ){
        final_matrix <- cbind.data.frame(final_matrix,drug_mat1)
      }else{
        no_result <- append(no_result,drug)
      }
    }
  }else{
    # For DMSO 
    final_matrix <- cbind.data.frame(ctrl,final_matrix)
    
  }
}


rownames(final_matrix)=rownames(ctrl)
gsva.input = as.matrix(final_matrix)
input.save = cbind(df[,1:2],final_matrix)
write.table(input.save,paste0(argv[4],"LINCS_GSVA_input.tsv"),sep='\t',quote = F)



# Apply gsva function
gsva.output = gsva(gsva.input, gset.idx.list = TF.index.list, method='ssgsea',
                   kcdf="Gaussian", abs.ranking=FALSE, min.sz=2,
                   max.sz=Inf, parallel.sz=20, mx.diff=TRUE,
                   ssgsea.norm=FALSE,
                   verbose=F) 

gsva.unscale=cbind(Symbol=rownames(gsva.output),as.data.frame(gsva.output))
TF.symbol=data.frame(rownames(gsva.output))
colnames(TF.symbol)="Symbol"
gene.id=df[,1:2]
fr.col=left_join(TF.symbol,gene.id,by="Symbol")
gsva.unscale=cbind(Entrez=fr.col$Entrez,gsva.unscale)
write.table(gsva.unscale,file=paste0(argv[4],"LINCS_GSVA.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
            row.names = FALSE)

# scaling by z-score
print("column scaling by z-score")
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
write.table(gsva.scale,file=paste0(argv[4],"LINCS_GSVA_z-score.tsv"),sep="\t",quote=FALSE,col.names = TRUE,
            row.names = FALSE)

print ("GSVA calculation is Finished")
print (no_result)
write.table(no_result,paste0(argv[4],"LINCS_maxdose_small_number.tsv"),sep='\t',quote = F, col.names=F)
