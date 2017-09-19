
# Adds column to the MAF that contains the unique tumor ID from the $Tumor_Sample_Barcode

unique.tumor.addition.function <- function(MAF.file,non.TCGA.characters.to.keep,sum.stats=T,figures=F){
  this.maf <- MAF.file
  TCGA.tumors <- grep(pattern = "TCGA",x = this.maf$Tumor_Sample_Barcode)
  other.tumors <- setdiff(1:nrow(this.maf),TCGA.tumors)
  
  first.12 <- function(string_to_12){
    return(paste(unlist(strsplit(string_to_12,split = ""))[1:12],collapse = ""))
  }
  
  if(class(non.TCGA.characters.to.keep)=="numeric"){
    first.other <- function(string_to_other){
      return(paste(unlist(strsplit(string_to_other,split = ""))[1:non.TCGA.characters.to.keep],collapse = ""))
    }
    
    this.maf$Unique_patient_identifier <- NA
    
    this.maf$Unique_patient_identifier[TCGA.tumors] <- unlist(lapply(this.maf$Tumor_Sample_Barcode[TCGA.tumors],first.12))
    this.maf$Unique_patient_identifier[other.tumors] <- unlist(lapply(this.maf$Tumor_Sample_Barcode[other.tumors],first.other))
    
    
  }else{
    this.maf$Unique_patient_identifier <- NA
    
    this.maf$Unique_patient_identifier[TCGA.tumors] <- unlist(lapply(this.maf$Tumor_Sample_Barcode[TCGA.tumors],first.12))
    this.maf$Unique_patient_identifier[other.tumors] <- this.maf$Tumor_Sample_Barcode[other.tumors]
    
  }
  
  if(sum.stats){
    unique.patients <- unique(this.maf$Unique_patient_identifier)
    
    tumor.mutation.number <- NULL
    for(i in 1:length(unique.patients)){
      tumor.mutation.number[i] <- nrow(this.maf[which(this.maf$Unique_patient_identifier==unique.patients[i]),])
      
    }
    
    if(figures){
      pdf(file = "Figures/mutations_per_tumor.pdf")  
      hist(tumor.mutation.number,breaks=100,xlab="Number of mutations per tumor",main="Histogram of the number of mutations in each tumor")
      dev.off()
    }
    
    message("Summary statistics of the number of mutations per unique tumor:")
    print(summary(tumor.mutation.number))
    
  }
  
  return(this.maf)   
}