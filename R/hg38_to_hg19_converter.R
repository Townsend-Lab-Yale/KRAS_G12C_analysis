
# Converts hg38 to hg19

# if the necessary chain isn't downloaded, download now. 
if(length(which(list.files("input_data")=="hg38Tohg19.chain"))==0){
  message("Downloading necessary chain file...")
  download.file(url = "http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/hg38ToHg19.over.chain.gz",destfile = "input_data/hg38Tohg19.chain.gz")
  install.packages("R.utils")
  require("R.utils")
  R.utils::gunzip("input_data/hg38Tohg19.chain.gz",remove=F)
}


hg38.to.hg19.converter <- function(chain,hg38_maf){
  
  ## ch = import.chain("~/Downloads/hg38ToHg19.over.chain") #chain needs to be downloaded from UCSC http://hgdownload.cse.ucsc.edu/gbdb/hg38/liftOver/
  ch = import.chain(chain)
  
  message("Loading in specified MAF...")
  hg38_maf <- hg38_maf
  
  hg38.Grange <- makeGRangesFromDataFrame(df = hg38_maf,keep.extra.columns = T,
                                          ignore.strand = F,
                                          seqinfo = NULL,
                                          seqnames.field = "Chromosome",
                                          start.field = "Start_Position",
                                          end.field = "End_Position",
                                          strand.field = "Strand",
                                          starts.in.df.are.0based = F)
  genome(hg38.Grange) <- "GRCh38"
  
  seqlevelsStyle(hg38.Grange) = "UCSC" 
  
  hg19.Grange <- liftOver(hg38.Grange,ch)
  
  hg19.Grange <- unlist(hg19.Grange)
  genome(hg19.Grange) <- "hg19"
  
  message(paste("Number of rows in the MAF that failed to convert: ", length(hg38.Grange) - length(hg19.Grange)))
  
  hg19.df <- as.data.frame(hg19.Grange)
  
  
  
  colnames(hg19.df)[which(colnames(hg19.df)=="seqnames")] <- "Chromosome"
  colnames(hg19.df)[which(colnames(hg19.df)=="start")] <- "Start_Position"
  colnames(hg19.df)[which(colnames(hg19.df)=="end")] <- "End_Position"
  hg19.df <- unique.tumor.addition.function(MAF.file = hg19.df,non.TCGA.characters.to.keep = 8,sum.stats = F,figures = T)
  
  
  hg19.df$NCBI_Build <- "Converted_from_GRCh38_to_hg19"
  
  #reduce chromosome column to just numbers/letters
  remove.3 <- function(string_to_3){
    return(paste(unlist(strsplit(as.character(string_to_3),split = ""))[4:length(unlist(strsplit(as.character(string_to_3),split = "")))],collapse = ""))
  }
  
  hg19.df$Chromosome <- unlist(lapply(X = hg19.df$Chromosome,FUN = remove.3))
  
  #To keep things consistent with always being on the "+" strand, we need to flip the genes that were just flipped. 
  to.flip <- which(hg19.df$strand=="-" & hg19.df$Variant_Type=="SNP")
  hg19.df$Reference_Allele[to.flip] <- unlist(lapply(X = hg19.df$Reference_Allele[to.flip],FUN = flip.function))
  hg19.df$Tumor_Seq_Allele1[to.flip] <- unlist(lapply(X = hg19.df$Tumor_Seq_Allele1[to.flip],FUN = flip.function))
  hg19.df$Tumor_Seq_Allele2[to.flip] <- unlist(lapply(X = hg19.df$Tumor_Seq_Allele2[to.flip],FUN = flip.function))
  
  return(hg19.df)
}