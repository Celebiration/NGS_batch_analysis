setwd("D:\\projects\\CRISPResso\\2024_10_12_xsj")
library(Biostrings)
library(stringr)
mapping<-read.csv("mapping.csv",header = T)
if (sum(duplicated(mapping$sample_id))>0) {
  stop("id重复！")
}
#去除adapter中的空格
mapping$adapter_5<-unlist(lapply(mapping$adapter_5,function(x){return(gsub(" ","",x))}))
mapping$adapter_3<-unlist(lapply(mapping$adapter_3,function(x){return(gsub(" ","",x))}))

#去除amplicon中的adapters(形如ADAPTERsequence，去除amplicon中的sequence部分)
for (i in 1:nrow(mapping)) {
  amplicon<-str_trim(toupper(mapping$amplicon_sequence[i]),side = "both")
  ad5<-toupper(str_trim(gsub("^[A-Z]+", "", mapping$adapter_5[i]),side = "both"))
  ad3<-as.character(reverseComplement(DNAString(str_trim(gsub("^[A-Z]+", "", mapping$adapter_3[i]),side = "both"))))
  ind<-0
  if (regexpr(paste0("^",ad5),amplicon)[[1]]!=-1) {
    amplicon<-gsub(pattern = paste0("^",ad5),replacement = "",x = amplicon)
    print(paste0(mapping$sample_id[i]," amplicon adapter5 trimmed."))
    ind<-1
  }else if(regexpr(paste0("^",ad5),as.character(reverseComplement(DNAString(amplicon))))[[1]]!=-1){
    amplicon<-gsub(pattern = paste0("^",ad5),replacement = "",x = as.character(reverseComplement(DNAString(amplicon))))
    print(paste0(mapping$sample_id[i]," amplicon adapter5 trimmed."))
    ind<-1
  }else{
    print(paste0("Warning: ",mapping$sample_id[i]," amplicon adapter5 untrimmed."))
  }
  if(ind==1){
    if (regexpr(paste0(ad3,"$"),amplicon)[[1]]!=-1) {
      amplicon<-gsub(pattern = paste0(ad3,"$"),replacement = "",x = amplicon)
      print(paste0(mapping$sample_id[i]," amplicon adapter3 trimmed."))
    }else{
      print(paste0("Warning: ",mapping$sample_id[i]," amplicon adapter3 untrimmed."))
    }
  }else{
    if (regexpr(paste0(ad3,"$"),amplicon)[[1]]!=-1) {
      amplicon<-gsub(pattern = paste0(ad3,"$"),replacement = "",x = amplicon)
      print(paste0(mapping$sample_id[i]," amplicon adapter3 trimmed."))
    }else if(regexpr(paste0(ad3,"$"),as.character(reverseComplement(DNAString(amplicon))))[[1]]!=-1){
      amplicon<-gsub(pattern = paste0(ad3,"$"),replacement = "",x = as.character(reverseComplement(DNAString(amplicon))))
      print(paste0(mapping$sample_id[i]," amplicon adapter3 trimmed."))
    }else{
      print(paste0("Warning: ",mapping$sample_id[i]," amplicon adapter3 untrimmed."))
    }
  }
  mapping$amplicon_sequence[i]<-amplicon
}


for (i in 1:nrow(mapping)) {
  sgRNA<-str_trim(toupper(mapping$sgRNA_sequence[i]),side = "both")
  amplicon<-str_trim(toupper(mapping$amplicon_sequence[i]),side = "both")
  if (regexpr(sgRNA,amplicon)[[1]]==-1) {
    if (regexpr(sgRNA,reverseComplement(DNAString(amplicon)))[[1]]==-1) {
      print(paste0("Error: ",mapping$sample_id[i]," sgRNA isn't in amplicon!"))
    }else{
      mapping$amplicon_sequence[i]<-as.character(reverseComplement(DNAString(amplicon)))
      print(paste0(mapping$sample_id[i]," amplicon strand fixed."))
    }
  }
  mapping$sgRNA_sequence[i]<-sgRNA
}

#更新！
write.csv(mapping,"mapping.csv",row.names = F,na="")
