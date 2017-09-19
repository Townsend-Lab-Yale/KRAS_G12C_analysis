
# Different MAFs use different conventions when determining the "tumor" allele. This attempts to bring consistency. 



tumor.allele.adder <- function(MAF){
  
  
  #Make sure all the characters in these columns are uppercase for consistency. 
  MAF$Tumor_Seq_Allele2 <- toupper(MAF$Tumor_Seq_Allele2)
  MAF$Tumor_Seq_Allele1 <- toupper(MAF$Tumor_Seq_Allele1)
  MAF$Reference_Allele <- toupper(MAF$Reference_Allele)
  
  #Delete rows with no information in either Tumor Seq columns. 
  if(length(which(MAF$Tumor_Seq_Allele1=="" & MAF$Tumor_Seq_Allele2=="")>0)){
    MAF <- MAF[-which(MAF$Tumor_Seq_Allele1=="" & MAF$Tumor_Seq_Allele2==""),]
  }
  
  #Need to determine which column contains the tumor allele. 
  tumor.allele.function <- function(MAF_input){
    ref <- MAF_input[,"Reference_Allele"]
    alt1 <- MAF_input[,"Tumor_Seq_Allele1"]
    alt2 <- MAF_input[,"Tumor_Seq_Allele2"]
    if(ref==alt1 & ref!=alt2){
      return(alt2)
    }else{
      if(ref==alt2 & ref!=alt1){
        return(alt1)
      }else{
        if(length(which(c(alt1,alt2)==""))==1){
          if(alt1 == "" & alt2 != ref){
            return(alt2)
          }else{
            if(alt2 == "" & alt1 != ref){
              return(alt1)
            }
          }
        }
      }
    }
  }
  
  MAF_allele_list <- split(x = MAF[,c("Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2")],f = rep(1:nrow(MAF),each=1))
  MAF$Tumor_allele <- unlist(lapply(X = MAF_allele_list,FUN = tumor.allele.function))
  
  return(MAF)
}