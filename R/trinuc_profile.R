

###function to output trinucleotide context using deconstructSigs

trinuc.profile.function <- function(input.MAF,
                                    signature.choice="signatures.cosmic",
                                    minimum.mutations.per.tumor=50,
                                    save.figs=F){
  
  ip <- as.data.frame(installed.packages()[,c(1,3:4)])
  if (!require("deconstructSigs")) {
    install.packages("deconstructSigs", dependencies = TRUE)
    library(deconstructSigs)
  }
  library(deconstructSigs)
  
  if(length(which(colnames(input.MAF)=="Unique_patient_identifier"))==0){
    warning("You require a column named \'Unique_patient_identifier\'with unique IDs for all tumors.\nReturning NULL")
    return(NULL)
  }
  
  input.MAF$Chromosome <- paste("chr",trimws(input.MAF$Chromosome),sep="")
  input.MAF <- input.MAF[which(input.MAF$Variant_Type=="SNP"),]
  
  #We do not want the signals for selection here, so we get rid of all recurrently mutated SNP (same start position and same chromosome)
  message("Removing all recurrent mutations... ")
  duplicated.vec.first <- duplicated(input.MAF$Start_Position)
  duplicated.vec.last <- duplicated(input.MAF$Start_Position,fromLast=T)
  
  duplicated.vec.pos <- which(duplicated.vec.first | duplicated.vec.last)
  duplicated.vec.start <- input.MAF$Start_Position[duplicated.vec.pos]
  
  
  
  # duplicated.positions <- MAF_for_analysis$Start_Position[duplicated(MAF_for_analysis$Start_Position)]
  
  remove.from.eee <- NULL
  for(eee in 1:length(duplicated.vec.pos)){
    
    positions <- which(input.MAF$Start_Position==duplicated.vec.start[eee])
    positions <- positions[-which(positions == duplicated.vec.pos[eee])]
    if(!(input.MAF$Chromosome[duplicated.vec.pos[eee]] %in% input.MAF$Chromosome[positions])){
      remove.from.eee <- c(remove.from.eee,eee)
    }
    
  }
  if(length(remove.from.eee)>0){
    duplicated.vec.pos <- duplicated.vec.pos[-remove.from.eee]
  }
  
  input.MAF <- input.MAF[-duplicated.vec.pos,]
  # no.recur.MAF <- MAF_for_analysis[include.vec,]
  
  #Finding unique tumors
  
  message("Finding the number of mutations per tumor")
  unique.patients <- unique(input.MAF$Unique_patient_identifier)
  tumor.mutation.number <- NULL
  for(i in 1:length(unique.patients)){
    tumor.mutation.number[i] <- nrow(input.MAF[which(input.MAF$Unique_patient_identifier==unique.patients[i]),])
    
  }
  
  
  patients.over.mut.minimum <- unique.patients[which(tumor.mutation.number>minimum.mutations.per.tumor)]
  
  message(paste("Number of tumors over specified minimum mutation number of ",minimum.mutations.per.tumor,": ",length(patients.over.mut.minimum),sep=""))
  
  message("Cleaning input to only contain tumors above the minimum...")
  
  #Only keeping patients above the minimum. There has to be a much faster way for this, but this works for now... 
  keep <- NULL
  for(i in 1:nrow(input.MAF)){
    if(input.MAF$Unique_patient_identifier[i] %in% patients.over.mut.minimum){keep <- c(keep,i)}
  }
  input.MAF <- input.MAF[keep,]
  
  
  message("Calculating trinucleotide mutation counts...")
  sigs.input <- mut.to.sigs.input(mut.ref = input.MAF, 
                                  sample.id = "Unique_patient_identifier", 
                                  chr = "Chromosome", 
                                  pos = "Start_Position", 
                                  ref = "Reference_Allele", 
                                  alt = "Tumor_Seq_Allele2")
  
  
  message("Calculating individual tumor mutational signatures...")
  
  if(signature.choice=="signatures.cosmic"){
    signatures.output <- list()
    for(i in 1:nrow(sigs.input)){
      signatures.output[[i]] <- whichSignatures(tumor.ref = sigs.input, 
                                                signatures.ref = signatures.cosmic, 
                                                sample.id = rownames(sigs.input)[i], 
                                                contexts.needed = TRUE,
                                                tri.counts.method = 'default')
    }
  }
  if(signature.choice=="signatures.nature2013"){
    signatures.output <- list()
    for(i in 1:nrow(sigs.input)){
      signatures.output[[i]] <- whichSignatures(tumor.ref = sigs.input, 
                                                signatures.ref = signatures.nature2013, 
                                                sample.id = rownames(sigs.input)[i], 
                                                contexts.needed = TRUE,
                                                tri.counts.method = 'default')
    }
    
  }
  
  
  # melt(rbind(signatures.output[[1]]$weights,signatures.output[[2]]$weights))
  
  if(length(signatures.output)>1){
    dir.create(path = "Figures", showWarnings = FALSE)
    weights.df <- signatures.output[[1]]$weights
    for(i in 2:length(signatures.output)){
      weights.df <- rbind(weights.df,signatures.output[[i]]$weights)
    }
    
    weights.df.melt <- melt(weights.df)
    
    if(save.figs){
      viol <- ggplot(data=weights.df.melt,aes(x=variable,y=value)) + geom_violin() + theme_bw() +  
        theme(axis.text.x = element_text(angle=90,hjust = 1)) +
        geom_jitter(width=0.25,aes(x=variable,y=value),alpha=0.2) + 
        labs(x="Cosmic signature") + labs(y="Signature weight") + ggtitle(paste("Signatures in ",tumor.name,sep="")) +
        theme(plot.title = element_text(hjust = 0.5))
      
      # viol
      ggsave(plot = viol,filename = paste("figures/Signature_Violin_plot_",tumor.name,".pdf",sep=""),units = "in",width = 11,height = 8)
    }
    
    
    
  }else{
    message("Only one signature calculated!")
    dir.create(path = "Figures", showWarnings = FALSE)
    weights.df <- signatures.output[[1]]$weights
    weights.df.melt <- melt(weights.df)
  }
  
  
  
  
  unknown.list <- rep(NA,length(signatures.output))
  for(i in 1:length(signatures.output)){
    unknown.list[i] <- signatures.output[[i]]$unknown
  }
  message("Statistical summary of the proportion of the mutational signature in each tumor sample that is \'unknown\'")
  print(summary(unknown.list))
  
  averaged.product <- signatures.output[[1]]$product
  for(i in 2:length(signatures.output)){
    averaged.product <- averaged.product + signatures.output[[i]]$product
  }
  averaged.product <- averaged.product/length(signatures.output)
  averaged.product <- averaged.product/sum(averaged.product)
  # sum(averaged.product)
  
  ###
  #Create trinucleotide data frame
  ###
  trinuc.mutation_data <- expand.grid(mutation=c("CtoA","CtoG","CtoT","TtoA","TtoC","TtoG"),Upstream=c("A","C","G","T"),Downstream=c("A","C","G","T"))
  trinuc.mutation_data$mutated_from <- rep(NA,nrow(trinuc.mutation_data)) 
  trinuc.mutation_data$mutated_to <- rep(NA,nrow(trinuc.mutation_data))   
  for(i in 1:nrow(trinuc.mutation_data)){
    trinuc.mutation_data$mutated_from[i] <- strsplit(as.character(trinuc.mutation_data$mutation[i]),split = "")[[1]][1]
  }
  for(i in 1:nrow(trinuc.mutation_data)){
    trinuc.mutation_data$mutated_to[i] <- strsplit(as.character(trinuc.mutation_data$mutation[i]),split = "")[[1]][4]
  }
  trinuc.mutation_data$total_count <- rep(0,nrow(trinuc.mutation_data)) ## 
  trinuc.mutation_data$proportion <- rep(0,nrow(trinuc.mutation_data)) ## 
  trinuc.mutation_data$section_labels <- factor(trinuc.mutation_data$mutation, labels = c("C%->%A", "C%->%G", "C%->%T","T%->%A", "T%->%C", "T%->%G"))
  
  
  
  for(i in 1:length(averaged.product)){
    splitname <- strsplit(colnames(averaged.product)[i],split="")[[1]]
    
    trinuc.mutation_data$proportion[which(trinuc.mutation_data$Upstream==splitname[1] &
                                            trinuc.mutation_data$Downstream==splitname[7] &
                                            trinuc.mutation_data$mutated_from==splitname[3] &
                                            trinuc.mutation_data$mutated_to==splitname[5])] <- averaged.product[i]
  }
  if(save.figs){
    p <- ggplot(data=trinuc.mutation_data, aes(Downstream, Upstream)) +
      geom_tile(aes(fill = proportion*100), colour = "white") + scale_fill_gradient(low = "white", high = "steelblue", name="Percent")
    p <- p + facet_grid(.~section_labels, labeller = label_parsed) 
    p <- p +  geom_text(aes(label = round(proportion, 4)*100),size=3)
    # p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    # panel.background = element_blank(), axis.line = element_line(colour = "black"))
    p <- p + theme_bw() + theme(panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.ticks = element_blank(),
                                strip.text=element_text(size=15),
                                axis.title.x = element_text(size=15),
                                axis.title.y = element_text(size=15),
                                axis.text.x = element_text(size=12),
                                axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5)) + 
      ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
    # p
    ggsave(paste("figures/",tumor.name,"_trinuc_heatmap.pdf",sep=""),height = 2.5,width = 10)
  }
  return(trinuc.mutation_data)
}