
# Calculates selection intensities


selection.intensity.calculation.function <- function(genes_for_analysis,
                                                     MAF_for_analysis,
                                                     this.substitution,
                                                     trinuc.mutation_data,
                                                     LabReference,
                                                     translations,
                                                     mut_rates,
                                                     low.mut,
                                                     mutsig_siggenes,
                                                     isoform.file=NULL,
                                                     tumor.number=NULL,
                                                     trinuc_count_matrix=NULL,
                                                     output_from_mainMAF=NULL){
  
  if(is.null(tumor.number)){
    warning("You did not provide the number of tumors in this dataset!\nUsing the unique tumors in the \'Unique_patient_identifier\' column of the MAF_for_analysis input.\nIf you broke the MAF into separate files, you will be missing the true number of tumors in this dataset!") 
    tumor.number <- length(unique(MAF_for_analysis$Unique_patient_identifier))
  }else{
    tumor.number <- tumor.number
  }
  
  # File containing isoform choices from a previous run to keep consistent variant calling. If not present ("NULL"), then make one to fill in. 
  if(is.null(isoform.file)){
    isoform.list <- as.data.frame(matrix(data=NA,ncol=3,nrow=length(unique(MAF_for_analysis$Hugo_Symbol))))
    colnames(isoform.list) <- c("Gene","Isoform","LabReference_location")
    isoform.list$Gene <- unique(MAF_for_analysis$Hugo_Symbol)
  }else{
    isoform.list <- isoform.file
  }
  
  # Loading data into the environment now so it is only called once per function. 
  LabReference <- LabReference 
  MAF_for_analysis <- MAF_for_analysis
  
  if(length(which(colnames(MAF_for_analysis)=="Start_position"))>0){
    colnames(MAF_for_analysis)[which(colnames(MAF_for_analysis)=="Start_position")] <- "Start_Position"
  }
  translations <- as.matrix(translations)
  rownames(translations) <- translations[,"Nucs"]
  
  # Creating two class-specific matrices of the trinucleotide context because it saves a lot of time looping through later
  trinuc.mutation_data.char <- as.matrix(trinuc.mutation_data[,c("mutation","Upstream","Downstream","mutated_from","mutated_to","section_labels")])
  trinuc.mutation_data.num <- as.matrix(trinuc.mutation_data[,c("total_count","proportion")])
  
  namesvec <- NULL
  for(i in 1:nrow(trinuc.mutation_data)){
    namesvec[i] <- paste0(trinuc.mutation_data$mutated_from[i],trinuc.mutation_data$Upstream[i],trinuc.mutation_data$mutated_to[i],trinuc.mutation_data$Downstream[i],collapse = "")
  }
  # namesvec
  rownames(trinuc.mutation_data.num) <- namesvec
  
  #Make sure all the characters in these columns are uppercase for consistency.
  MAF_for_analysis$Tumor_Seq_Allele2 <- toupper(MAF_for_analysis$Tumor_Seq_Allele2)
  MAF_for_analysis$Tumor_Seq_Allele1 <- toupper(MAF_for_analysis$Tumor_Seq_Allele1)
  MAF_for_analysis$Reference_Allele <- toupper(MAF_for_analysis$Reference_Allele)
  
  
  all.muts <- as.data.frame(matrix(data=NA,nrow=0,ncol=16),header=T) #Data frame to hold mutation data
  colnames(all.muts) <- c("Gene",
                          "AA_Pos",
                          "Nucleotide_position",
                          "Nuc_Ref",
                          "Nuc_Change",
                          "AA_Ref",
                          "AA_Change",
                          "gamma",
                          "gamma_epistasis",
                          "freq",
                          "mu",
                          "gene_AA_size",
                          "MutSigCV_p",
                          "MutSigCV_q",
                          "Prop_tumors_with_specific_mut",
                          "Prop_of_tumors_with_this_gene_mutated") #Column names for this data frame
  
  
  #Create a data frame to hold all the mutational data. 
  mutation.data <- as.data.frame(matrix(data=NA,nrow=0,ncol=30),header=T)
  colnames(mutation.data) <- c("Gene",
                               "Gene_isoform",
                               "Gene_size",
                               "Percent_A",
                               "Percent_T",
                               "Percent_G",
                               "Percent_C",
                               "Nucleotide_Gene_Pos",
                               "Nucleotide_chromosome_position",
                               "Chromosome",
                               "Reference_Nucleotide",
                               "Alternative_Nucleotide",
                               "Reference_Count",
                               "Alternative_Count",
                               "Tumor_origin",
                               "Nucleotide_change_tally",
                               "Nucleotide_mutation_rate",
                               "Amino_acid_position",
                               "Amino_acid_codon",
                               "Codon_position",
                               "Amino_acid_reference",
                               "Amino_acid_alternative",
                               "Amino_acid_mutation_rate",
                               "Amino_acid_change_tally",
                               "Gamma",
                               "MAF_location",
                               "Nucleotide_trinuc_context",
                               "synonymous.mu",
                               "strand",
                               "trinucs")
  
  #matrix to hold trinucleotide count data
  if(is.null(trinuc_count_matrix)){
    trinuc.count.matrix <- matrix(nrow=96,ncol=length(unique(MAF_for_analysis$Hugo_Symbol)),data=NA) 
    colnames(trinuc.count.matrix) <- unique(MAF_for_analysis$Hugo_Symbol)
    rownames(trinuc.count.matrix) <- namesvec
  }
  
  #Going to cut all substitutions out of the analysis that do not match the reference genome. Will keep track. 
  if(is.null(output_from_mainMAF)){
    mutations.to.cut <- as.data.frame(matrix(data=NA,nrow=1,ncol=3))
    colnames(mutations.to.cut) <- c("Hugo_Symbol","Chromosome","Start_Position")
  }
  
  ####Loop to iterate through every gene----
  for(zzz in 1:length(genes_for_analysis)){
    
    
    this.gene <- genes_for_analysis[zzz] #The specific gene for this run is... <- 
    
    #If the name of this gene is "Unknown", tell us! 
    if(this.gene=="Unknown"){
      print("This gene is named UNKNOWN!")
      next
    }
    
    
    ### finding isoform----
    
    
    #If there is no reference for this gene name in our LabReference file, tell us! 
    if(length(which(LabReference$geneName==this.gene))==0){
      print(paste("There is no reference for isoforms of gene ",this.gene," !"))
      next
    }
    
    
    
    
    this.gene.ref.loc <- which(LabReference$geneName==this.gene) # get locations in LabReference for this gene name
    this.gene.isoforms <- LabReference$name[this.gene.ref.loc] # get isoform names for this gene 
    
    this.gene_MAF <- MAF_for_analysis[which(MAF_for_analysis$Hugo_Symbol==this.gene & MAF_for_analysis$Variant_Type=="SNP"),] #get a gene-specific MAF file with just the SNPs we want to analyze 
    
    
    #Make sure this gene is in the MAF file, and then... 
    if(nrow(this.gene_MAF)>0){
      
      
      if(is.null(output_from_mainMAF)){
        if(is.na(isoform.list[which(isoform.list[,"Gene"]==this.gene),"Isoform"])){ #if the isoform choice does not exist, pick it
          
          length.vec <- rep(NA,length=length(this.gene.ref.loc)) #A vector that holds information on how long each isoform is
          mut.count.vec <- rep(0,length=length(this.gene.ref.loc)) #A vector that holds information on how many mutations fall within each isoform. 
          
          
          ###Looping through isoforms to find best choice ----
          ##Loop to iterate through each isoform to find the one with 
          ##the most mutations and longest length
          for(i in 1:length(this.gene.ref.loc)){
            
            #For each isoform we need to construct the sequence, and place the mutations
            
            this.gene.ref.loc.choice <- this.gene.ref.loc[i] #get the location within the reference file for this specific isoform choice 
            
            this.strand <- LabReference$strand[this.gene.ref.loc.choice] #which strand is this isoform on? 
            
            #Pull cds information
            cds.start <- LabReference$cdsStart[this.gene.ref.loc.choice]
            cds.end <- LabReference$cdsEnd[this.gene.ref.loc.choice]
            
            
            #Sometimes the cds.start == the cds.end?? If so, we do not have any information on this isoform, skip it. 
            if(cds.start==cds.end){
              print(paste("The cds.start equals the cds.end for gene: ",this.gene))
              next
            }
            
            #Pull information about all the exon start and end locations for this isoform
            exon.starts <- as.numeric(unlist(strsplit(LabReference$exonStarts[this.gene.ref.loc.choice],split=",")))
            exon.ends <- as.numeric(unlist(strsplit(LabReference$exonEnds[this.gene.ref.loc.choice],split=",")))
            
            
            start.pos <- which(exon.ends>cds.start)[1] #exon with cds start in it
            end.pos <- which(exon.starts<cds.end)[length(which(exon.starts<cds.end))] #exon with cds end in it. 
            seq.exon.start.vec <- rep(NA,length=length(start.pos:end.pos)) #vector to hold all the information about where exons start
            seq.exon.end.vec <- rep(NA,length=length(start.pos:end.pos)) #vector to hold all the information about where exons end
            
            seq.exon.start.vec[1] <- cds.start #first exon starts at the cds.start location
            seq.exon.end.vec[length(seq.exon.end.vec)] <- cds.end #and last exon ends there 
            
            #fill in the rest of the information if there is more than one exon
            if(length(seq.exon.start.vec)>1){
              seq.exon.start.vec[2:length(seq.exon.start.vec)] <- exon.starts[(start.pos+1):end.pos] 
              seq.exon.end.vec[1:(length(seq.exon.end.vec)-1)] <- exon.ends[start.pos:(end.pos-1)]
            }
            
            # print(paste("This isoform:", LabReference$name[this.gene.ref.loc.choice], "This isoform length:",sum(seq.exon.end.vec-seq.exon.start.vec),sep=" "))
            
            length.vec[i] <- sum(seq.exon.end.vec-seq.exon.start.vec) #total length of this isoform
            
            #Loop through all mutations
            this.count <- 0 
            for(ii in 1:nrow(this.gene_MAF)){
              #for each mutation, this position is... 
              this.mut.pos <- this.gene_MAF$Start_Position[ii]
              #Is this location within an exon? then add a count! 
              for(iii in 1:length(seq.exon.start.vec)){
                if(this.mut.pos>seq.exon.start.vec[iii] & this.mut.pos<seq.exon.end.vec[iii]){this.count <- this.count+1}
              }
            }
            
            mut.count.vec[i] <- this.count #Record total number of mutations within the isoform
            # print(paste("Mutations in this isoform:",this.count,sep=" "))
          }
          if(cds.start==cds.end){
            print(paste("The cds.start equals the cds.end for gene: ",this.gene))
            # next
          }
          if(mean(mut.count.vec)==0){
            print("Mutations given in the MAF do not fall into any of the isoforms we reference!")
            next
          }
          
          #If there is one isoform that has the most mutations within it, we have a winner. 
          if(length(which(mut.count.vec==max(mut.count.vec)))==1){
            this.isoform <- LabReference$name[this.gene.ref.loc[which(mut.count.vec==max(mut.count.vec))]] #this is the isoform name
            pick <- which(mut.count.vec==max(mut.count.vec)) #and the location within our LabReference locations we pulled
          }else{
            #if there is a tie for most mutations, pick largest isoform. If there is a tie for this as well, pick first
            if(length(which(mut.count.vec==max(mut.count.vec)))>1){
              max.mut.vec <- which(mut.count.vec==max(mut.count.vec))
              pick <- max.mut.vec[which(length.vec[max.mut.vec]==max(length.vec[max.mut.vec]))[1]]
              this.isoform <- LabReference$name[this.gene.ref.loc[[pick]]]
              
            }
          }
          this.gene.ref.loc.choice <- this.gene.ref.loc[pick] #get the location within the reference file for this specific gene choice
          isoform.list[which(isoform.list[,"Gene"]==this.gene),"Isoform"] <- this.isoform
          isoform.list[which(isoform.list[,"Gene"]==this.gene),"LabReference_location"] <- this.gene.ref.loc.choice
        }else{
          #If the isoform choice exists already, pick that one
          this.isoform <- isoform.list[which(isoform.list[,"Gene"]==this.gene),"Isoform"]
          this.gene.ref.loc.choice <- isoform.list[which(isoform.list[,"Gene"]==this.gene),"LabReference_location"]
        }
        
        
        print(paste("Isoform choice: ",this.isoform," Gene: ",this.gene,sep="")) #Tell us what we picked
        
        this.strand <- LabReference$strand[this.gene.ref.loc.choice] #get which strand this isoform information is on
        
        #Isoform position data 
        cds.start <- LabReference$cdsStart[this.gene.ref.loc.choice]
        cds.end <- LabReference$cdsEnd[this.gene.ref.loc.choice]
        
        exon.starts <- as.numeric(unlist(strsplit(LabReference$exonStarts[this.gene.ref.loc.choice],split=",")))
        exon.ends <- as.numeric(unlist(strsplit(LabReference$exonEnds[this.gene.ref.loc.choice],split=",")))
        
        start.pos <- which(exon.ends>cds.start)[1] #exon with cds start in it
        end.pos <- which(exon.starts<cds.end)[length(which(exon.starts<cds.end))] #exon with cds end in it. 
        seq.exon.start.vec <- rep(NA,length=length(start.pos:end.pos))
        seq.exon.end.vec <- rep(NA,length=length(start.pos:end.pos))
        
        seq.exon.start.vec[1] <- cds.start
        seq.exon.end.vec[length(seq.exon.end.vec)] <- cds.end
        if(length(seq.exon.start.vec)>1){
          seq.exon.start.vec[2:length(seq.exon.start.vec)] <- exon.starts[(start.pos+1):end.pos] 
          seq.exon.end.vec[1:(length(seq.exon.end.vec)-1)] <- exon.ends[start.pos:(end.pos-1)]
        }
        
        
        ###Getting isoform sequence data ----
        if(this.strand=="+"){
          myseq <- paste(as.character(getSeq(Hsapiens,LabReference$chrom[this.gene.ref.loc.choice],start=seq.exon.start.vec+1,end=seq.exon.end.vec,strand="+")),collapse = "")
          
          
          #If this isoform strand is "-", we need to flip the positions and strand to make it "+" 
        }else{
          
          myseq <- paste(rev(as.character(getSeq(Hsapiens,LabReference$chrom[this.gene.ref.loc.choice],start=seq.exon.start.vec+1,end=seq.exon.end.vec,strand="-"))),collapse = "")
          
          
        }
        
        
        ###replace it with the mutation----
        
        if(this.substitution[1]==this.gene){
          tmp.string <- strsplit(myseq,"")[[1]]
          tmp.string[as.numeric(this.substitution[2])] <- this.substitution[3]
          myseq <- paste(tmp.string,collapse = "")
        }
        
        
        
        
        
        tmp.string <- strsplit(myseq,"")[[1]] #temporary string with our sequence split up
        
        
        myseq.codons <- paste0(tmp.string[c(T,F,F)],tmp.string[c(F,T,F)],tmp.string[c(F,F,T)]) #every 3 as a group for the codon designations 
        
        
        AA.code <- rep(NA,length=length(myseq.codons)) #The amino acid code will go here
        
        #Loop to give the amino acid designations for each codon. Can make this a function to speed it up. TODO
        for(i in 1:length(myseq.codons)){
          AA.code[i] <- translations[which(translations[,1]==myseq.codons[i]),4]
        }
        
        #If there is no mutation rate data for this gene, tell us.
        if(length(mut_rates$r_x_X[which(mut_rates$gene==this.gene)])==0){
          print(paste("There is no mutation rate data for gene ",this.gene))
          next
        }
        
        #Get the synonymous mutation rate from the MutSigCV output. If it does not exist, use the lowest mutation rate calculated for the bunch. 
        if(mut_rates$r_x_X[which(mut_rates$gene==this.gene)]!=0){
          synonymous.mu <- mut_rates$r_x_X[which(mut_rates$gene==this.gene)]
        }else{
          synonymous.mu <- as.numeric(low.mut[which(low.mut[,1]=="Min non-zero x/X:"),2])
        }
        
        
        #What proportion of all the mutation rates is each specific mutation
        # mut.proportion <- mut.rates/sum(mut.rates)
        
        
        ##turn the sequence into a vector  
        myseq.split <- unlist(strsplit(myseq,split=""))
        
        
        #Chromosome for this gene
        chromosome <- LabReference$chrom[this.gene.ref.loc.choice]
        
        #Sometimes the same gene name is listed for multiple chromosomes. This throws off the sequencing, reduce the MAF for this gene to the gene chosen by isoform search
        this.gene_MAF <- this.gene_MAF[which(this.gene_MAF$Chromosome==paste(unlist(strsplit(chromosome,split=""))[4:length(unlist(strsplit(chromosome,split="")))],collapse ="")),]
        
        # cutting out substitutions that do not match the reference genome
        matcher.vec <- as.vector(this.gene_MAF$Reference_Allele == as.character(getSeq(Hsapiens,chromosome,start=this.gene_MAF$Start_Position,end=this.gene_MAF$End_Position,strand="+")))
        if(length(which(matcher.vec==FALSE))>0){
          not.matching <- which(matcher.vec==F)
          
          mutations.to.cut <- rbind(mutations.to.cut,this.gene_MAF[not.matching,c("Hugo_Symbol","Chromosome","Start_Position")])
          this.gene_MAF <- this.gene_MAF[-not.matching,]
          if(nrow(this.gene_MAF)==0){
            print(paste("After removing inconsistencies between the reference genome and the MAF file, there are no substitutions left to analyze."))
            next
          }
        }
        
        
        #vectors to hold the position data at the gene-level and chromosome-level
        gene.pos <- rep(NA,length=((seq.exon.end.vec[length(seq.exon.end.vec)]-seq.exon.start.vec[1])+2))
        chrom.gene.pos <- rep(NA,length=((seq.exon.end.vec[length(seq.exon.end.vec)]-seq.exon.start.vec[1])+2))
        
        #Vector of the total sequence from start to end, including introns 
        if(this.strand=="-"){
          pos.sequence  <- rev(strsplit(as.character(getSeq(Hsapiens,chromosome,start=seq.exon.start.vec[1],end=seq.exon.end.vec[length(seq.exon.end.vec)]+1,strand="-")),split="")[[1]])
        }else{      
          pos.sequence  <- strsplit(as.character(getSeq(Hsapiens,chromosome,start=seq.exon.start.vec[1],end=seq.exon.end.vec[length(seq.exon.end.vec)]+1,strand="+")),split="")[[1]]
        }
        
        #All the chromosome positions in the gene, including introns
        chrom.pos <- seq.exon.start.vec[1]:(seq.exon.end.vec[length(seq.exon.end.vec)]+1)
        
        #All the positions in the chromosome containing exons from this gene
        nuc.pos.vec <- c()
        for(i in 1:length(seq.exon.start.vec)){
          nuc.pos.vec <- c(nuc.pos.vec,(seq.exon.start.vec[i]+1):seq.exon.end.vec[i])
        }
        
        
        
        # positions within the chromosome vector containing the gene vector
        these.positions <- which(chrom.pos %in% nuc.pos.vec)
        
        #fill in the gene-level position information, based on its position in the chromosome 
        for(i in 1:length(nuc.pos.vec)){  #This could be made a function to speed it up #TODO
          # this.pos <- which(chrom.pos==nuc.pos.vec[i])
          this.pos <- these.positions[i]
          chrom.gene.pos[this.pos] <- nuc.pos.vec[i]
          if(this.strand=="+"){
            gene.pos[this.pos] <- i
          }else{
            gene.pos[this.pos] <- (length(nuc.pos.vec)+1)-i
          }
        }
        
        
        ##If the substitution is here, replace to determine new trinucleotide context. 
        if(this.substitution[1]==this.gene){
          pos.sequence[which(gene.pos==this.substitution[2])] <- this.substitution[3]
        }
        
        
        #Generate the trinucleotide context at every position. The NA at the end are placeholders and are not at postions used in calculations 
        minus.vec <- c(pos.sequence[-1],NA)
        plus.vec <- c(NA,pos.sequence[-length(pos.sequence)])
        if(this.strand=="-"){
          trinuc.pos <- paste(minus.vec,pos.sequence,plus.vec,sep="")
        }else{
          trinuc.pos <- paste(plus.vec,pos.sequence,minus.vec,sep="")
        }
        
        
        
        ##Nucleotide matrix----
        #Constructing a matrix to hold all the mutation frequency data
        mut.matrix <- matrix(nrow=4,
                             ncol=length(myseq.split),
                             data=NA)
        rownames(mut.matrix) <- c("A","T","G","C") #Rows will be the possible nucleotides
        colnames(mut.matrix) <- paste("Pos. ",seq(1,length(myseq.split),by=1)," Ref: ",myseq.split,sep="") #columns will be positions in the gene
        
        mut.matrix.trinuc <- matrix(nrow=nrow(mut.matrix),
                                    ncol=ncol(mut.matrix),
                                    data=NA)
        rownames(mut.matrix.trinuc) <- c("A","T","G","C") #Rows will be the possible nucleotides
        colnames(mut.matrix.trinuc) <- paste("Pos. ",seq(1,length(myseq.split),by=1)," Ref: ",myseq.split,sep="") #columns will be positions in the gene
        
        
        
        #Matrix that holds the mutation data for same-nucleotide mutations. This has zero observed rate
        zeros.mat <- as.data.frame(matrix(nrow=1,ncol=4,data=0)) 
        colnames(zeros.mat) <- c("AtoA","TtoT","GtoG","CtoC")
        # mut.proportion.combined <-  merge(mut.proportion.combined,zeros.mat) #add this to the previous mutation rate data frame
        
        #Fill in the matrix with all the possible mutation proportions
        position.vec <- which(!(is.na(gene.pos)))
        if(this.strand=="-"){position.vec <- position.vec[order(-position.vec)]}
        
        trinuc.mutation_data.num.this.gene <- trinuc.mutation_data.num
        # trinuc.data.this.gene <- trinuc.mutation_data
        # trinuc.data.this.gene$total_count <- 0
        trinuc.mutation_data.num.this.gene[,"total_count"] <- 0
        
        #Loop to count the trinucleotide contexts of this gene, and assign the site-specific mutation rate
        for(i in 1:length(myseq.split)){   # I feel like this could be made into a function and sped up with the lapply function. TODO
          # this.pos <- which(gene.pos==i)
          
          
          this.trinuc.split <- unlist(strsplit(trinuc.pos[position.vec[i]],split = ""))
          
          if(this.trinuc.split[2]=="C" | this.trinuc.split[2]=="T"){
            
            
            mut.matrix[this.trinuc.split[2],i] <- 0
            if(this.trinuc.split[2]=="C"){
              trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"T",this.trinuc.split[3],collapse = "")),
                                                 "total_count"] <-           trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"T",this.trinuc.split[3],collapse = "")),
                                                                                                                "total_count"]+1
              mut.matrix[c("A","G","T"),i] <- c(trinuc.mutation_data.num[c(paste0("C",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                           paste0("C",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                                           paste0("C",this.trinuc.split[1],"T",this.trinuc.split[3],collapse = "")),
                                                                         "proportion"])
              mut.matrix.trinuc[c("A","G","T"),i] <- c(paste0("C",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                       paste0("C",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""),
                                                       paste0("C",this.trinuc.split[1],"T",this.trinuc.split[3],collapse = ""))
              
            }else{
              trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                   paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = "")),
                                                 "total_count"] <-           trinuc.mutation_data.num.this.gene[c(paste0(this.trinuc.split[2],this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                                                                                  paste0(this.trinuc.split[2],this.trinuc.split[1],"G",this.trinuc.split[3],collapse = "")),"total_count"]+1
              
              # mut.matrix[c("A","C","G"),i] <- tri.nuc.flip.total.prop[positions.in.barplot] 
              mut.matrix[c("A","C","G"),i] <- c(trinuc.mutation_data.num[c(paste0("T",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                                           paste0("T",this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                                           paste0("T",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = "")),
                                                                         "proportion"])
              
              mut.matrix.trinuc[c("A","C","G"),i] <- c(paste0("T",this.trinuc.split[1],"A",this.trinuc.split[3],collapse = ""),
                                                       paste0("T",this.trinuc.split[1],"C",this.trinuc.split[3],collapse = ""),
                                                       paste0("T",this.trinuc.split[1],"G",this.trinuc.split[3],collapse = ""))
              
              
            }
          }else{
            
            mut.matrix[this.trinuc.split[2],i] <- 0
            
            if(this.trinuc.split[2]=="G"){
              trinuc.mutation_data.num.this.gene[c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"] <- trinuc.mutation_data.num.this.gene[c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"]+1
              
              mut.matrix[c("T","C","A"),i] <- c(trinuc.mutation_data.num[c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = "")),"proportion"])
              
              mut.matrix.trinuc[c("T","C","A"),i] <- c(paste0("C",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("C",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("C",flip.function(this.trinuc.split[3]),"T",flip.function(this.trinuc.split[1]),collapse = ""))
              
            }else{
              
              trinuc.mutation_data.num.this.gene[c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                   paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"] <- trinuc.mutation_data.num.this.gene[c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                                                                                                                                                                 paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = "")),"total_count"]+1 
              
              mut.matrix[c("T","G","C"),i] <- c(trinuc.mutation_data.num[c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                                           paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = "")),
                                                                         "proportion"])
              
              mut.matrix.trinuc[c("T","G","C"),i] <- c(paste0("T",flip.function(this.trinuc.split[3]),"A",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("T",flip.function(this.trinuc.split[3]),"C",flip.function(this.trinuc.split[1]),collapse = ""),
                                                       paste0("T",flip.function(this.trinuc.split[3]),"G",flip.function(this.trinuc.split[1]),collapse = ""))
            }
          }
        }
        
        trinuc.count.matrix[,this.gene] <- trinuc.mutation_data.num.this.gene[,"total_count"]
        
        # Do not need to generate the gene in the future if we calculate the mean rate in the gene once. 
        # This will speed up the bootstrap analysis later. 
        to.mean <-NULL
        for(w in 1:32){
          
          to.mean[w] <- trinuc.mutation_data.num.this.gene[w*3,"total_count"] * sum(trinuc.mutation_data.num.this.gene[c(w*3,(w*3)-1,(w*3)-2),"proportion"]) 
        }
        mean.nuc.rate.for.this.gene <- (sum(to.mean)/(sum(trinuc.mutation_data.num.this.gene[,"total_count"])/3)) 
        
        #Now that we have the mean, we can calculate the mutation rate of nucleotides without generating the gene. 
        
        #Find the sum of all rates at every position
        nuc.mutation.rate.vec <- rep(NA,length=length(myseq.split))
        for(i in 1:length(myseq.split)){
          nuc.mutation.rate.vec[i] <- sum(mut.matrix[,i])
        }
        
        #matrix that has all the synonymous mutation rates
        normalized.mut.matrix <-  (mut.matrix/mean(nuc.mutation.rate.vec))*synonymous.mu
        
        ##Amino acid matrix----
        #Matrix that stores rates for each amino acid mutation
        AA.mut.matrix <- matrix(nrow=21,
                                ncol=length(AA.code),
                                data=0)
        rownames(AA.mut.matrix) <- unique(translations[,3]) #Rows are all the possible nucleotides 
        
        #Names of the columns of this matrix will be the amino acid at that location
        aa.code.short <- rep(NA,length=length(AA.code))
        for(i in 1:length(aa.code.short)){
          aa.code.short[i] <- translations[which(translations[,4]==AA.code[i])[1],3]
        }
        colnames(AA.mut.matrix) <- paste("Pos. ",seq(1,length(AA.code),by=1)," Ref: ",aa.code.short,sep="")
        
        AA.matrix.trinuc <- matrix(nrow=nrow(AA.mut.matrix),
                                   ncol=ncol(AA.mut.matrix),
                                   data=NA)
        rownames(AA.matrix.trinuc) <- rownames(AA.mut.matrix)
        
        #########
        #Probabilities of each amino acid mutation----
        #########
        
        #Number/positional designation for each codon in terms of nucleotide position
        codon.numbers <- rep(1:length(myseq.codons),each=3)
        nuc.list <- c("A","T","G","C")
        
        #Loop to calculate mutation rates at each codon to other codons. 
        for(i in 1:length(myseq.split)){
          this.codon <- myseq.codons[codon.numbers[i]] #What is this codon?
          #see which position in the codon we are looking at 
          if(i%%3==1){pos <- 1};if(i%%3==2){pos <- 2};if(i%%3==0){pos <- 3} #gives back position within the codon
          this.codon.expanded <- unlist(strsplit(this.codon,split="")) #breaks this codon into individual nucleotides 
          # if(this.codon.expanded[2]=="G" | this.codon.expanded[2]=="A"){
          #   
          # }
          mutations <- nuc.list[which(nuc.list!=this.codon.expanded[pos])] #possible mutations at this location within the codon
          
          #Loop to go through all the possible mutations at this position in the codon and add up rates  
          for(j in 1:3){
            this.mutation <- mutations[j] #This specific mutation is... 
            this.codon.expanded[pos] <- this.mutation #Making this codon... 
            this.codon.collapsed <- paste0(this.codon.expanded,collapse="") #which we collapse... 
            new.aa <- translations[this.codon.collapsed,3] #resulting in this new amino acid... 
            AA.mut.matrix[new.aa,codon.numbers[i]] <- AA.mut.matrix[new.aa,codon.numbers[i]] + normalized.mut.matrix[this.mutation,i] #adding the rate this happens to total rates matrix
            AA.matrix.trinuc[new.aa,codon.numbers[i]] <- ifelse(is.na(AA.matrix.trinuc[new.aa,codon.numbers[i]]),mut.matrix.trinuc[this.mutation,i],paste(AA.matrix.trinuc[new.aa,codon.numbers[i]],mut.matrix.trinuc[this.mutation,i],sep=":"))
          }
        }
        
        #matrix that stores info on the complimentary strand 
        strand.switch <- as.data.frame(matrix(data=c("A","T","G","C"),ncol=1,nrow=4),stringsAsFactors = FALSE)
        rownames(strand.switch) <- c("T","A","C","G")
        
        
        ##Amino acid matrix that stores gamma information
        AA.gamma.matrix <- matrix(nrow=21,
                                  ncol=length(AA.code),
                                  data=NA)
        
        rownames(AA.gamma.matrix) <- rownames(AA.mut.matrix)
        colnames(AA.gamma.matrix) <- colnames(AA.mut.matrix)
        
        ##Amino acid that stores mutation incidence data
        AA.tally.matrix <- matrix(nrow=21,
                                  ncol=length(AA.code),
                                  data=0)
        
        rownames(AA.tally.matrix) <- rownames(AA.mut.matrix)
        colnames(AA.tally.matrix) <- colnames(AA.mut.matrix)
        
        
        #Matrix to store the incidence of mutations at the nucleotide level 
        Nuc.tally.matrix <- matrix(nrow=4,
                                   ncol=length(myseq.split),
                                   data=0)
        # rownames(Nuc.tally.matrix) <- c("A","T","G","C")
        # colnames(Nuc.tally.matrix) <- paste("Pos. ",seq(1,length(myseq.split),by=1)," Ref: ",myseq.split,sep="")
        rownames(Nuc.tally.matrix) <- rownames(mut.matrix)
        colnames(Nuc.tally.matrix) <- colnames(mut.matrix)
        
        
        #Matrix to tally up the mutations outside of exons
        outside.nucs <- as.data.frame(matrix(nrow=nrow(this.gene_MAF),ncol=5,data=NA))
        colnames(outside.nucs) <- c("nuc_pos","ref","change","tally","mutation_rate")
        outside.nucs$tally <- 0
        
        ##Setting up the data frame to store (and eventually compile) all the analysis
        mutation.data.to.add <- as.data.frame(matrix(data=NA,nrow=nrow(this.gene_MAF),ncol=31),header=T)
        colnames(mutation.data.to.add) <- c("Gene",
                                            "Gene_isoform",
                                            "Gene_size",
                                            "Percent_A",
                                            "Percent_T",
                                            "Percent_G",
                                            "Percent_C",
                                            "Nucleotide_Gene_Pos",
                                            "Nucleotide_chromosome_position",
                                            "Chromosome",
                                            "Reference_Nucleotide",
                                            "Alternative_Nucleotide",
                                            "Reference_Count",
                                            "Alternative_Count",
                                            "Tumor_origin",
                                            "Unique_patient_identifier",
                                            "Nucleotide_change_tally",
                                            "Nucleotide_mutation_rate",
                                            "Amino_acid_position",
                                            "Amino_acid_codon",
                                            "Codon_position",
                                            "Amino_acid_reference",
                                            "Amino_acid_alternative",
                                            "Amino_acid_mutation_rate",
                                            "Amino_acid_change_tally",
                                            "Gamma",
                                            "MAF_location",
                                            "Nucleotide_trinuc_context",
                                            "synonymous.mu",
                                            "strand",
                                            "trinucs")
        
        mutation.data.to.add$Gene <- this.gene
        mutation.data.to.add$Gene_isoform <- this.isoform
        mutation.data.to.add$Gene_size <- length(myseq.split)
        mutation.data.to.add$Percent_A <- length(which(myseq.split=="A"))/length(myseq.split)
        mutation.data.to.add$Percent_T <- length(which(myseq.split=="T"))/length(myseq.split)
        mutation.data.to.add$Percent_G <- length(which(myseq.split=="G"))/length(myseq.split)
        mutation.data.to.add$Percent_C <- length(which(myseq.split=="C"))/length(myseq.split)
        mutation.data.to.add$Chromosome <- this.gene_MAF$Chromosome[1]
        mutation.data.to.add$synonymous.mu <- synonymous.mu
        mutation.data.to.add$strand <- this.strand
        
        mismatched.with.ref.genome <- NULL
        
        ####Tallying up mutations ----
        ##Loop that tallys up mutations in this gene, and their positions / affect on amino acid sequence. 
        for(i in 1:nrow(this.gene_MAF)){
          
          #get the position within the gene and the amino acid code 
          
          this.pos <- gene.pos[which(chrom.gene.pos==this.gene_MAF$Start_Position[i])]
          
          this.AA.pos <- codon.numbers[this.pos]

          mutation.data.to.add$MAF_location[i] <- rownames(this.gene_MAF)[i]
          mutation.data.to.add$Nucleotide_chromosome_position[i] <- this.gene_MAF$Start_Position[i]
          
          #what was the reference and the mutation
          if(this.strand=="-"){
            this.ref <- strand.switch[this.gene_MAF$Reference_Allele[i],]
            this.mut <- strand.switch[this.gene_MAF$Tumor_allele[i],]
          }else{
            this.ref <- this.gene_MAF$Reference_Allele[i]
            this.mut <- this.gene_MAF$Tumor_allele[i]
          }
          
          mutation.data.to.add$Reference_Nucleotide[i] <- this.ref
          mutation.data.to.add$Alternative_Nucleotide[i] <- this.mut
          mutation.data.to.add$Tumor_origin[i] <- this.gene_MAF$Tumor_Sample_Barcode[i]
          mutation.data.to.add$Unique_patient_identifier[i] <- this.gene_MAF$Unique_patient_identifier[i]
          if(length(which(colnames(this.gene_MAF)=="t_ref_count"))>0){mutation.data.to.add$Reference_Count[i] <- this.gene_MAF$t_ref_count[i]}
          if(length(which(colnames(this.gene_MAF)=="t_alt_count"))>0){mutation.data.to.add$Alternative_Count[i] <- this.gene_MAF$t_alt_count[i]}
          
          #if this position exists (else, it wasn't in the isoform and tell us)...
          if(length(this.pos)>0){
            mutation.data.to.add$Nucleotide_Gene_Pos[i] <- this.pos
            mutation.data.to.add$Amino_acid_position[i] <- this.AA.pos
            mutation.data.to.add$Nucleotide_mutation_rate[i] <- normalized.mut.matrix[this.mut,this.pos]
            
            
            mutation.data.to.add$Nucleotide_trinuc_context[i] <- trinuc.pos[which(chrom.gene.pos==this.gene_MAF$Start_Position[i])]
            
            #get which position within the codon the mutation occurred in 
            if(this.pos%%3==1){pos <- 1};if(this.pos%%3==2){pos <- 2};if(this.pos%%3==0){pos <- 3}
            mutation.data.to.add$Codon_position[i] <- pos
            this.codon <- myseq.codons[this.AA.pos]
            mutation.data.to.add$Amino_acid_codon[i] <- this.codon
            this.codon.mutated <- unlist(strsplit(this.codon,split=""))
            this.codon.mutated[pos] <- this.mut
            this.codon.mutated <- paste(this.codon.mutated,collapse = "")
            
            
            this.codon.expanded <- unlist(strsplit(this.codon,split=""))
            
            this.AA.change <- translations[which(translations[,"Nucs"]==this.codon.mutated),"AA_short"]
            
            
            ##tally up the mutation in the tally matrices
            Nuc.tally.matrix[this.mut,this.pos] <- Nuc.tally.matrix[this.mut,this.pos]+1 
            AA.tally.matrix[this.AA.change,this.AA.pos] <- AA.tally.matrix[this.AA.change,this.AA.pos]+1
            
            mutation.data.to.add$Amino_acid_reference[i] <- translations[which(translations[,"Nucs"]==this.codon)[1],"AA_letter"]
            mutation.data.to.add$Amino_acid_alternative[i] <- translations[which(translations[,"AA_short"]==this.AA.change)[1],"AA_letter"]
            
          }else{
            # print(paste("A mutation is outside the isoform of choice! The gene-specific MAF file position is: ",i,sep=""))
            
            if(length(which(outside.nucs$nuc_pos==this.gene_MAF$Start_Position[i] & outside.nucs$change==this.mut))==0){
              open.spot <- which(is.na(outside.nucs$nuc_pos))[1]
              outside.nucs$nuc_pos[open.spot] <- this.gene_MAF$Start_Position[i]
              outside.nucs$ref[open.spot] <- this.ref
              outside.nucs$change[open.spot] <- this.mut
              outside.nucs$tally[open.spot] <- 1 
              
              #calculate the mutation rate at this sepecific site... 
              if(this.strand=="-"){
                outside.sequence  <- rev(strsplit(as.character(getSeq(Hsapiens,chromosome,start=this.gene_MAF$Start_Position[i]-1,end=this.gene_MAF$Start_Position[i]+1,strand="-")),split="")[[1]])
              }else{      
                outside.sequence  <- strsplit(as.character(getSeq(Hsapiens,chromosome,start=this.gene_MAF$Start_Position[i]-1,end=this.gene_MAF$Start_Position[i]+1,strand="+")),split="")[[1]]
              }
              
              mutation.data.to.add$Nucleotide_trinuc_context[i] <- paste(outside.sequence,collapse = "") 
              
              #calculating the mutation rate for this gene 
              if(outside.nucs$ref[open.spot]== "C" | outside.nucs$ref[open.spot]== "T"){
                if(outside.sequence[2]==outside.nucs$ref[open.spot]){
                  outside.nucs$mutation_rate[open.spot] <- (trinuc.mutation_data.num.this.gene[paste0(outside.sequence[2],outside.sequence[1],outside.nucs$change[open.spot],outside.sequence[3],collapse = ""),"proportion"]/mean.nuc.rate.for.this.gene)*synonymous.mu
                  mutation.data.to.add$trinucs[i] <- paste0(outside.sequence[2],outside.sequence[1],outside.nucs$change[open.spot],outside.sequence[3],collapse = "")
                }else{
                  #This means that the reference nucleotide in the MAF does not match the reference genome position. 
                  message(paste("The reference allele given in the MAF file does not match the reference genome!\nThis gene: ",this.gene," \nThis position in the gene-specific MAF file: ",i)) #Let us know
                  mismatched.with.ref.genome <- c(mismatched.with.ref.genome,i)
                }
              }else{
                if(outside.sequence[2]==outside.nucs$ref[open.spot]){
                  outside.nucs$mutation_rate[open.spot] <- (trinuc.mutation_data.num.this.gene[paste0(flip.function(outside.sequence[2]),flip.function(outside.sequence[3]),flip.function(outside.nucs$change[open.spot]),flip.function(outside.sequence[1]),collapse = ""),"proportion"]/mean.nuc.rate.for.this.gene)*synonymous.mu
                  mutation.data.to.add$trinucs[i] <- paste0(flip.function(outside.sequence[2]),flip.function(outside.sequence[3]),flip.function(outside.nucs$change[open.spot]),flip.function(outside.sequence[1]),collapse = "")
                }else{
                  #This means that the reference nucleotide in the MAF does not match the reference genome position. 
                  message(paste("The reference allele given in the MAF file does not match the reference genome!\nThis gene: ",this.gene," \nThis position in the gene-specific MAF file: ",i)) #Let us know
                  mismatched.with.ref.genome <- c(mismatched.with.ref.genome,i)
                }
              }
              mutation.data.to.add$Nucleotide_mutation_rate[i] <- outside.nucs$mutation_rate[open.spot]
            }else{
              outside.nucs$tally[which(outside.nucs$nuc_pos==this.gene_MAF$Start_Position[i] & outside.nucs$change==this.mut)] <- outside.nucs$tally[which(outside.nucs$nuc_pos==this.gene_MAF$Start_Position[i] & outside.nucs$change==this.mut)]+1
              
            }
            
          }
        }
        
        outside.nucs <- outside.nucs[-which(is.na(outside.nucs$nuc_pos)),] #reducing data frame size
        if(length(mismatched.with.ref.genome)>0){
          mutation.data.to.add <- mutation.data.to.add[-mismatched.with.ref.genome,] 
        }
        
        #codon.numbers
        
        mut.spots.nuc <- which(Nuc.tally.matrix!=0,arr.ind = T)
        if(nrow(mut.spots.nuc)==0){
          print("No mutations to tally!")
          next
        }
        for(i in 1:nrow(mut.spots.nuc)){
          these.pos <- which(mutation.data.to.add$Nucleotide_Gene_Pos==mut.spots.nuc[i,2] & mutation.data.to.add$Alternative_Nucleotide==rownames(mut.spots.nuc)[i])
          mutation.data.to.add$Nucleotide_change_tally[these.pos] <- length(these.pos)
        }
        
        
        
        #TODO, do not need to generate entire AA.gamma.matrix, may be wasting memory and computation time. 
        mut.spots <- which(AA.tally.matrix!=0,arr.ind=T)
        if(nrow(mut.spots)==0){
          print("No mutations to tally!")
          next
        }
        # as.numeric(mut.spots[1,][2])
        
        #tallying up the changes! 
        for(i in 1:nrow(mut.spots)){
          these.pos <- which(mutation.data.to.add$Amino_acid_position==mut.spots[i,2] & mutation.data.to.add$Amino_acid_alternative==translations[which(translations[,"AA_short"]==rownames(mut.spots)[i])[1],"AA_letter"])
          
          mutation.data.to.add$Amino_acid_change_tally[these.pos] <- length(these.pos)
          
          mutation.data.to.add$Amino_acid_mutation_rate[these.pos] <- AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]]
          mutation.data.to.add$trinucs[these.pos] <- AA.matrix.trinuc[mut.spots[i,1],mut.spots[i,2]]
          
          
          if(AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]]==0){message(paste("The mutation rate for this mutation is listed as ZERO"));print(mut.spots[i,])
          }else{
            AA.gamma.matrix[mut.spots[i,1],mut.spots[i,2]] <- lambda.calc(n.one = AA.tally.matrix[mut.spots[i,1],mut.spots[i,2]],n.zero = (tumor.number-AA.tally.matrix[mut.spots[i,1],mut.spots[i,2]]))/AA.mut.matrix[mut.spots[i,1],mut.spots[i,2]]
          }
          
          mutation.data.to.add$Gamma[these.pos] <-  AA.gamma.matrix[mut.spots[i,1],mut.spots[i,2]] 
        }
        
        if(nrow(outside.nucs)>0){
          outside.nucs$selection_intensity <- NA
          outside.nucs$selection_no_epistasis <- NA
          for(i in 1:nrow(outside.nucs)){
            
            outside.nucs$selection_intensity[i] <- lambda.calc(n.one = outside.nucs$tally[i],n.zero = (tumor.number - length(unique(this.gene_MAF$Unique_patient_identifier))))/outside.nucs$mutation_rate[i] 
            outside.nucs$selection_no_epistasis[i] <- lambda.calc(n.one = outside.nucs$tally[i],n.zero = (tumor.number - outside.nucs$tally[i]))/outside.nucs$mutation_rate[i]  
            
            mutation.data.to.add$Nucleotide_change_tally[which(mutation.data.to.add$Nucleotide_chromosome_position==outside.nucs$nuc_pos[i])] <- outside.nucs$tally[i]
          }
          
        }
        
        #Add the gamma assuming complete epistasis calculation to the complete_mutation_data dataframe
        mutation.data.to.add$Gamma_epistasis <- NA
        
        for(i in 1:nrow(mutation.data.to.add)){
          if(is.na(mutation.data.to.add$Amino_acid_position[i])){
            mutation.data.to.add$Gamma_epistasis[i] <- lambda.calc(n.one = mutation.data.to.add$Nucleotide_change_tally[i],n.zero = (tumor.number - length(unique(this.gene_MAF$Unique_patient_identifier))))/mutation.data.to.add$Nucleotide_mutation_rate[i]
          }else{
            mutation.data.to.add$Gamma_epistasis[i] <- lambda.calc(n.one = mutation.data.to.add$Amino_acid_change_tally[i],n.zero = (tumor.number - length(unique(this.gene_MAF$Unique_patient_identifier))))/mutation.data.to.add$Amino_acid_mutation_rate[i]
          }
        }
        
        
        this.gene.gamma <- which(!is.na(AA.gamma.matrix),arr.ind=T)
        added.rows <- as.data.frame(matrix(data=NA,nrow=(nrow(this.gene.gamma)+nrow(outside.nucs)),ncol=ncol(all.muts)),header=T)
        colnames(added.rows) <- colnames(all.muts)
        added.rows$Gene <- this.gene
        added.rows$gene_AA_size <- length(myseq.codons)
        added.rows$MutSigCV_p <- mutsig_siggenes$p[which(mutsig_siggenes$gene==this.gene)]
        added.rows$MutSigCV_q <- mutsig_siggenes$q[which(mutsig_siggenes$gene==this.gene)]
        added.rows$Prop_of_tumors_with_this_gene_mutated <- length(unique(this.gene_MAF$Unique_patient_identifier))/tumor.number
        
        counter <- 1
        for(e in 1:nrow(this.gene.gamma)){
          added.rows$AA_Pos[e] <- this.gene.gamma[e,2]
          added.rows$AA_Ref[e] <- AA.code[added.rows$AA_Pos[e]]
          added.rows$AA_Change[e] <- translations[which(translations[,"AA_short"]==rownames(AA.gamma.matrix)[this.gene.gamma[e,1]])[1],"AA_letter"]
          added.rows$gamma[e] <- AA.gamma.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]
          added.rows$mu[e] <- AA.mut.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]
          added.rows$freq[e] <- AA.tally.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]
          added.rows$Prop_tumors_with_specific_mut[e] <- AA.tally.matrix[this.gene.gamma[e,1],this.gene.gamma[e,2]]/tumor.number
          added.rows$gamma_epistasis[e] <- lambda.calc(n.one = added.rows$freq[e],n.zero = (tumor.number - length(unique(this.gene_MAF$Unique_patient_identifier))))/added.rows$mu[e]
          counter <- counter + 1
        }
        
        if(nrow(outside.nucs)>0){
          for(i in 1:nrow(outside.nucs)){
            added.rows[counter,c("Nucleotide_position","Nuc_Ref","Nuc_Change","freq","mu","gamma_epistasis","gamma")] <- 
              outside.nucs[i,c("nuc_pos","ref","change","tally","mutation_rate","selection_intensity","selection_no_epistasis")]
            added.rows[counter,"Prop_tumors_with_specific_mut"] <- added.rows[counter,"freq"]/tumor.number
            counter <- counter + 1
          }
        }

        # print(added.rows)
        all.muts <- rbind(all.muts,added.rows)
        mutation.data <- rbind(mutation.data,mutation.data.to.add)
      }else{
        #This is what happens if you already have the mutation data (isoform, trinucleotide distribution per gene, trinucleotide context per substitution)
        #output_from_mainMAF should have all possible mutations in the bootstrap MAFs
        
        if(length(which(output_from_mainMAF$complete_mutation_data$Gene==this.gene))==0){
          print(paste("There is no previous selection data for ",this.gene))
          next
        }
        
        
        main_MAF.this.gene <- output_from_mainMAF$complete_mutation_data[which(output_from_mainMAF$complete_mutation_data$Gene==this.gene),]
        this.strand <- main_MAF.this.gene$strand[1]
        
        trinuc.mutation_data.num.this.gene <- trinuc.mutation_data.num
        trinuc.mutation_data.num.this.gene[,"total_count"] <- output_from_mainMAF$trinuc_counts[,this.gene]
        
        to.mean <-NULL
        for(w in 1:32){
          to.mean[w] <- trinuc.mutation_data.num.this.gene[w*3,"total_count"] * sum(trinuc.mutation_data.num.this.gene[c(w*3,(w*3)-1,(w*3)-2),"proportion"]) 
        }
        mean.nuc.rate.for.this.gene <- (sum(to.mean)/(sum(trinuc.mutation_data.num.this.gene[,"total_count"])/3)) 
        
        #sometimes the same gene name is on different chromosomes. This filters to the isoform of choice in the previous analysis
        this.gene_MAF <- this.gene_MAF[which(this.gene_MAF$Chromosome==main_MAF.this.gene$Chromosome[1]),] 
        
        # Need to cut out substitutions that are inconsistent between the MAF file and the reference genome. These were recorded in the previous run. 
        if(length(which(output_from_mainMAF$cut.mutations$Hugo_Symbol==this.gene))>0){
          this.gene_MAF <- this.gene_MAF[-which(this.gene_MAF$Hugo_Symbol%in%output_from_mainMAF$cut.mutations$Hugo_Symbol & this.gene_MAF$Start_Position%in%output_from_mainMAF$cut.mutations$Start_Position & this.gene_MAF$Chromosome%in%output_from_mainMAF$cut.mutations$Chromosome),]
        }
        
        
        if(nrow(this.gene_MAF)==0){
          print(paste("All the data for this gene is from another chromosome, after cleaning based on the isoform picked for the original analysis"))
          next
        }
        
        #If it is on the - strand, we need to flip the variants to that strand so it is consistent with previous output. MAFs always contain "+" strand data 
        if(this.strand=="-"){
          this.gene_MAF$Reference_Allele <- unlist(lapply(X = this.gene_MAF$Reference_Allele,FUN = flip.function))
          this.gene_MAF$Tumor_allele <- unlist(lapply(X = this.gene_MAF$Tumor_allele,FUN = flip.function))
        }
        
        if(mut_rates$r_x_X[which(mut_rates$gene==this.gene)]!=0){
          synonymous.mu <- mut_rates$r_x_X[which(mut_rates$gene==this.gene)]
        }else{
          synonymous.mu <- as.numeric(low.mut[which(low.mut[,1]=="Min non-zero x/X:"),2])
        }
        
        
        
        added.rows <- as.data.frame(matrix(data=NA,nrow=nrow(this.gene_MAF),ncol=ncol(all.muts)),header=T)
        colnames(added.rows) <- colnames(all.muts)
        added.rows$Gene <- this.gene
        added.rows$MutSigCV_p <- mutsig_siggenes$p[which(mutsig_siggenes$gene==this.gene)]
        added.rows$MutSigCV_q <- mutsig_siggenes$q[which(mutsig_siggenes$gene==this.gene)]
        added.rows$Prop_of_tumors_with_this_gene_mutated <- length(unique(this.gene_MAF$Unique_patient_identifier))/tumor.number
        
        for(mafRow in 1:nrow(this.gene_MAF)){
          
          
          if(!is.na(main_MAF.this.gene$Amino_acid_position[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow])[1]])){ #within the exon
            if(length(which(added.rows$AA_Pos == main_MAF.this.gene$Amino_acid_position[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow])[1]] & added.rows$AA_Change == main_MAF.this.gene$Amino_acid_alternative[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow] & main_MAF.this.gene$Alternative_Nucleotide==this.gene_MAF$Tumor_allele[mafRow])[1]]))==0){ #If this mutation isn't already stored in the output
              this.row <- which(is.na(added.rows$mu))[1]
              
              added.rows$AA_Pos[this.row] <- main_MAF.this.gene$Amino_acid_position[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow])[1]]
              added.rows$AA_Ref[this.row] <- main_MAF.this.gene$Amino_acid_reference[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow])[1]]
              this.AA.change <- main_MAF.this.gene$Amino_acid_alternative[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF$Tumor_allele[mafRow])[1]]
              #Need the rows in the main_MAF.this.gene that are also this AA pos and change, and their unique nucleotide positions and nucleotide changes, and then 
              #if the this.gene_MAF file has these as well, to tally up! 
              added.rows$AA_Change[this.row] <- this.AA.change
              
              nuc.changes <- unlist(strsplit(main_MAF.this.gene$trinucs[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF$Tumor_allele[mafRow])[1]],split=":"))
              
              
              added.rows$mu[this.row] <- sum(trinuc.mutation_data.num.this.gene[nuc.changes,"proportion"])/mean.nuc.rate.for.this.gene*synonymous.mu
              
              
              #now, need to find all spots in the this.geneMAF that correspond to this AA change.
              #Need to find all mutation types in the main_MAF file that correspond to this AA change, and count them in the this.geneMAF file
              this.main_MAF.pos <- which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF$Tumor_allele[mafRow])[1]
              all.AA.mut.positions <- which(main_MAF.this.gene$Amino_acid_alternative==main_MAF.this.gene$Amino_acid_alternative[this.main_MAF.pos] & main_MAF.this.gene$Amino_acid_position==main_MAF.this.gene$Amino_acid_position[this.main_MAF.pos])
              
              #all nuc positions and changes that resulted in this AA pos and AA change
              these.nuc.pos <- main_MAF.this.gene$Nucleotide_chromosome_position[all.AA.mut.positions]
              these.nuc.change <- main_MAF.this.gene$Alternative_Nucleotide[all.AA.mut.positions]
              
              this.gene_MAF.nucs <- which(this.gene_MAF$Start_Position %in% these.nuc.pos)
              this.gene_MAF.change <- which(this.gene_MAF$Tumor_allele %in% these.nuc.change)
              
              added.rows$freq[this.row] <- length(intersect(this.gene_MAF.nucs,this.gene_MAF.change))
              
              added.rows$gamma_epistasis[this.row] <- lambda.calc(n.one = added.rows$freq[this.row],n.zero = (tumor.number-length(unique(this.gene_MAF$Unique_patient_identifier))))/added.rows$mu[this.row]
              added.rows$gamma[this.row] <- lambda.calc(n.one = added.rows$freq[this.row],n.zero = (tumor.number-added.rows$freq[this.row]))/added.rows$mu[this.row]
              
            } #if it is, just skip it. Should have stored all the information on previous loops
            
            
            
          }else{ #outside the exon 
            
            #only do the calculation if it is not already present in the added.rows dataframe
            if(length(which(added.rows$Nucleotide_position==this.gene_MAF$Start_Position[mafRow] & added.rows$Nuc_Change==this.gene_MAF$Tumor_allele[mafRow]))==0){
              
              this.row <- which(is.na(added.rows$mu))[1]
              added.rows$Nucleotide_position[this.row] <- this.gene_MAF$Start_Position[mafRow]
              added.rows$Nuc_Ref[this.row] <- this.gene_MAF$Reference_Allele[mafRow]
              added.rows$Nuc_Change[this.row] <- this.gene_MAF$Tumor_allele[mafRow]
              
              nuc.changes <- unlist(strsplit(main_MAF.this.gene$trinucs[which(main_MAF.this.gene$Nucleotide_chromosome_position==this.gene_MAF$Start_Position[mafRow] & main_MAF.this.gene$Alternative_Nucleotide == this.gene_MAF$Tumor_allele[mafRow])[1]],split=":"))
              added.rows$mu[this.row] <- sum(trinuc.mutation_data.num.this.gene[nuc.changes,"proportion"])/mean.nuc.rate.for.this.gene*synonymous.mu
              
              added.rows$freq[this.row] <- length(which(this.gene_MAF$Start_Position==added.rows$Nucleotide_position[this.row] & this.gene_MAF$Tumor_allele==added.rows$Nuc_Change[this.row]))
              
              added.rows$gamma_epistasis[this.row] <- lambda.calc(n.one = added.rows$freq[this.row],n.zero = (tumor.number-length(unique(this.gene_MAF$Unique_patient_identifier))))/added.rows$mu[this.row]
              
              added.rows$gamma[this.row] <- lambda.calc(n.one = added.rows$freq[this.row],n.zero = (tumor.number-added.rows$freq[this.row]))/added.rows$mu[this.row]
              
            }
          }
          
          
 
          
        }
        added.rows <- added.rows[which(!is.na(added.rows$mu)),]
        added.rows$Prop_tumors_with_specific_mut <- added.rows$freq/tumor.number
        all.muts <- rbind(all.muts,added.rows)
      }
    }else{
      print(paste("The gene ",this.gene," is not in the MAF file!",sep=""))
    }
    
    if(zzz%%1000==0){message(paste("Gene number: ",zzz," out of ",length(genes_for_analysis)))} #If zzz is a multiple of 1000 print the progress
    
  }
  if(exists('nuc.mutation.rate.vec')){ #this catch means the function generated the gene sequence. If not, just store what we have 
    output.list <- list(all_mutations=all.muts,
                        complete_mutation_data=mutation.data,
                        nucleotide_mutation_rates=nuc.mutation.rate.vec,
                        nucleotide_tally=Nuc.tally.matrix,
                        amino_acid_tally=AA.tally.matrix,
                        amino_acid_mutation_rates=AA.mut.matrix,
                        norm_mut=normalized.mut.matrix,
                        myseqsplit=myseq.split,
                        isoform.list=isoform.list,
                        trinuc_counts=trinuc.count.matrix,
                        cut.mutations=mutations.to.cut[-is.na(mutations.to.cut$Hugo_Symbol),])
  }else{
    
    all.muts$identifier <- NA
    for(i in 1:nrow(all.muts)){
      if(!(is.na(all.muts$AA_Pos[i]))){
        all.muts$identifier[i] <- paste0(all.muts$Gene[i]," ",all.muts$AA_Ref[i],all.muts$AA_Pos[i],all.muts$AA_Change[i])
      }else{
        all.muts$identifier[i] <- paste0(all.muts$Gene[i]," ",all.muts$Nuc_Ref[i],all.muts$Nucleotide_position[i],all.muts$Nuc_Change[i])
      }
    }
    
    output.list <- list(all_mutations=all.muts,
                        complete_mutation_data=mutation.data)
    
  }
  
  return(output.list)
  
}