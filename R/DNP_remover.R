
# ##Function to remove possible DNP from

DNP.remover <- function(MAF){
  message("Removing possible DNP")
  remove_it <- MAF
  final_MAF <- NULL
  tumor.list <- unique(remove_it$Unique_patient_identifier)
  counter <- 0
  possible.DNP <- NULL
  for(i in 1:length(tumor.list)){
    to.delete <- NULL
    this.tumor_full <- remove_it[which(remove_it$Unique_patient_identifier==tumor.list[i] & remove_it$Variant_Type!="SNP"),]
    this.tumor <- remove_it[which(remove_it$Unique_patient_identifier==tumor.list[i] & remove_it$Variant_Type=="SNP"),]
    for(j in 1:nrow(this.tumor)){
      if(length(which(this.tumor$Chromosome==this.tumor$Chromosome[j] &
                      (this.tumor$Start_Position==this.tumor$Start_Position[j]+1 |
                       this.tumor$Start_Position==this.tumor$Start_Position[j]-1 )))>0){
        to.delete <- c(to.delete,j)
        counter <- counter+1
      }
    }
    if(length(to.delete)>0){
      this.tumor_full <- rbind(this.tumor[-to.delete,],this.tumor_full)
      final_MAF <- rbind(final_MAF,this.tumor_full)
    }else{
      this.tumor_full <- rbind(this.tumor,this.tumor_full)
      final_MAF <- rbind(final_MAF,this.tumor_full)
    }
  }
  message(paste("Total count of potential DNP removed: ", counter))
  message("DNP removal complete")
  return(final_MAF)
}