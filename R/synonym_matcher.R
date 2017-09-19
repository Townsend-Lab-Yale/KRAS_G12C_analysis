
# Creates consistent Hugo_Symbol names when applicable. Tries to match covariate expression files used with MutSigCV.  

synonym.matcher.function <- function(input_to_be_changed,synonym.table,covariates.file,isoforms_or_MAF){
  if(length(which(colnames(covariates.file)=="hugo"))>0){
    colnames(covariates.file)[which(colnames(covariates.file)=="hugo")] <- "gene"
  }
  if(isoforms_or_MAF=="MAF"){
    this_MAF <- input_to_be_changed
    synonyms_matrix <- synonym.table
    genes.to.investigate <- unique(this_MAF$Hugo_Symbol)
    genes.in.covariate <- covariates.file$gene
    counter <- 0
    for(i in 1:length(genes.to.investigate)){
      if(!(genes.to.investigate[i] %in% genes.in.covariate)){
        if(length(which(synonyms_matrix$Symbol == genes.to.investigate[i]))>0){
          # print("Match!");print(i)
          if(any(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate)){
            if(length(which(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate))==1){
              if(all(synonyms_matrix$Chromosome[which(synonyms_matrix$Symbol==genes.to.investigate[i])[which(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate)]] == this_MAF$Chromosome[which(this_MAF$Hugo_Symbol==genes.to.investigate[i])])){ #if all the chromosomes in the MAF equal the chromosome for this synonym... change it to the synonym
                this_MAF$Hugo_Symbol[which(this_MAF$Hugo_Symbol == genes.to.investigate[i])] <- synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])[which(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate)]]
                counter <- counter + 1
              }
            }else{
              print(paste("More than 1 match in the synonym table for ",genes.to.investigate[i]," in the covariates file! Could not assign.",sep=""))
            }
          }
        }
      }
    }
    message(paste("Number of changed gene names: ",counter))
    return(this_MAF)
  }
  if(isoforms_or_MAF=="isoforms"){
    
    chr.trimmer <- function(string.to.trim){
      paste(unlist(strsplit(string.to.trim,split = ""))[4:nchar(string.to.trim)],sep="",collapse = "")
    }
    
    this_MAF <- input_to_be_changed
    synonyms_matrix <- synonym.table
    genes.to.investigate <- unique(this_MAF$geneName)
    genes.in.covariate <- covariates.file$gene
    counter <- 0
    for(i in 1:length(genes.to.investigate)){
      if(!(genes.to.investigate[i] %in% genes.in.covariate)){
        if(length(which(synonyms_matrix$Symbol == genes.to.investigate[i]))>0){
          # print("Match!");print(i)
          if(any(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate)){
            if(length(which(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate))==1){
              if(all(synonyms_matrix$Chromosome[which(synonyms_matrix$Symbol==genes.to.investigate[i])[which(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate)]] == unlist(lapply(this_MAF$chrom[which(this_MAF$geneName==genes.to.investigate[i])],chr.trimmer)))){ #if all the chromosomes in the MAF equal the chromosome for this synonym... change it to the synonym
                this_MAF$geneName[which(this_MAF$geneName == genes.to.investigate[i])] <- synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])[which(synonyms_matrix$Synonyms[which(synonyms_matrix$Symbol==genes.to.investigate[i])] %in% genes.in.covariate)]]
                counter <- counter + 1
              }
            }else{
              print(paste("More than 1 match in the synonym table for ",genes.to.investigate[i]," in the covariates file! Could not assign.",sep=""))
            }
          }
        }
      }
    }
    message(paste("Number of changed gene names: ",counter))
    return(this_MAF)
  }
}