####
# Analysis performed in the manuscript Cannataro et al., The likelihood of heterogeneity or additional mutation in KRAS or associated oncogenes to compromise targeting of oncogenic KRAS G12C
###

dir.create("output_data")


# Import NCI data ----
# File Properties
# Name	TCGA.LUAD.mutect.81ccaef3-4550-494d-882c-895fb5a3de3b.DR-7.0.somatic.maf.gz
# Access	open
# UUID	81ccaef3-4550-494d-882c-895fb5a3de3b

NCI_MAF_38 <- read.csv('input_data/NCI/gdc_download_20170919_153409/81ccaef3-4550-494d-882c-895fb5a3de3b/TCGA.LUAD.mutect.81ccaef3-4550-494d-882c-895fb5a3de3b.DR-7.0.somatic.maf',stringsAsFactors=F,skip=5,header=T,sep='\t')


# Convert NCI data to hg19, process, and save ----

require(rtracklayer)
source("R/hg38_to_hg19_converter.R")
source("R/unique_tumor_addition.R")
source("R/flip_function.R")
source("R/DNP_remover.R")
source("R/tumor_allele_add.R")

NCI_MAF_19 <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=NCI_MAF_38)
NCI_MAF_19 <- DNP.remover(MAF = NCI_MAF_19)
NCI_MAF_19$Tumor_Seq_Allele2 <- toupper(NCI_MAF_19$Tumor_Seq_Allele2)
NCI_MAF_19$Reference_Allele <- toupper(NCI_MAF_19$Reference_Allele)
NCI_MAF_19 <- tumor.allele.adder(MAF = NCI_MAF_19)

save(NCI_MAF_19,file='output_data/MAF_NCI_19.RData')
write.table(NCI_MAF_19,file='output_data/MAF_NCI_19.txt',quote = F,sep='\t',row.names = F)


# Import Yale-Gilead data, processing, combining, saving ----

source("R/merging_TCGA_and_local_MAF.R")
source("R/unique_tumor_addition.R")
source("R/flip_function.R")
source("R/DNP_remover.R")
source("R/tumor_allele_add.R")

Yale_data <- read.csv('input_data/YG/adc_inc_counts.txt',stringsAsFactors=F,header=T,sep='\t')
if(length(which(colnames(Yale_data)=="keep"))>0){
  Yale_data <- Yale_data[-which(Yale_data$keep==0),]
}
load("output_data/MAF_NCI_19.RData")

MAF_for_analysis <- merging.TCGA.and.Yale.MAF.data.function(NCI_data = NCI_MAF_19,
                                                            Yale_data = Yale_data
)

MAF_for_analysis <- unique.tumor.addition.function(MAF.file = MAF_for_analysis,non.TCGA.characters.to.keep = 'all',figures=T)
MAF_for_analysis <- DNP.remover(MAF = MAF_for_analysis)
MAF_for_analysis$Tumor_Seq_Allele2 <- toupper(MAF_for_analysis$Tumor_Seq_Allele2)
MAF_for_analysis$Reference_Allele <- toupper(MAF_for_analysis$Reference_Allele)
MAF_for_analysis <- tumor.allele.adder(MAF = MAF_for_analysis)
save(MAF_for_analysis,file='output_data/MAF_LUAD.RData')
write.table(MAF_for_analysis,file='output_data/MAF_LUAD.txt',quote = F,sep='\t',row.names = F)


# Matching gene name synonyms ----

source("R/synonym_matcher.R")

load("output_data/MAF_LUAD.RData")

MAF_for_analysis <- synonym.matcher.function(input_to_be_changed = MAF_for_analysis,synonym.table = read.csv(file='input_data/matching_table.txt',stringsAsFactors = F,sep='\t'),covariates.file = read.table(file='input_data/LUNG_expression.txt',sep='\t',header = T,stringsAsFactors = F),isoforms_or_MAF = 'MAF')
MAF_for_analysis <- tumor.allele.adder(MAF = MAF_for_analysis)

isoforms <- synonym.matcher.function(input_to_be_changed = read.csv(file='input_data/UCSC_refseq_refFlat_hg19_short.txt',header = T,sep='\t',stringsAsFactors = F),synonym.table = read.csv(file='input_data/matching_table.txt',stringsAsFactors = F,sep='\t'),covariates.file = read.table(file='input_data/LUNG_expression.txt',sep='\t',header = T,stringsAsFactors = F),isoforms_or_MAF = 'isoforms')

save(MAF_for_analysis,file='output_data/MAF_LUAD.RData')
write.table(MAF_for_analysis,file='output_data/MAF_LUAD.txt',quote = F,sep='\t',row.names = F)
save(isoforms,file='output_data/Isoforms_LUAD.RData')


# Determine trinucleotide profile ---- 

source("R/trinuc_profile.R")
require(ggplot2)
require(reshape2)
tumor.name <- "LUAD"

load("output_data/MAF_LUAD.RData")

trinuc.mutation_data <- trinuc.profile.function(input.MAF = MAF_for_analysis,save.figs=T)

save(trinuc.mutation_data,file='output_data/trinuc_data_LUAD.RData')


load("output_data/trinuc_data_LUAD.RData")

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
                            axis.text.y=element_text(size=12),plot.title = element_text(hjust = 0.5)) 
  # ggtitle(paste("Trinucleotide profile for ",tumor.name,sep=""))
p
ggsave(paste("figures/",tumor.name,"_trinuc_heatmap_noname.pdf",sep=""),height = 2.5,width = 10)





# Calculating mutation rates with MutSigCV ----

dir.create("MutSigCV") #Download and unzip MutSigCV 1.41 here 

###
# Download MutSigCV 1.41 from http://archive.broadinstitute.org/cancer/cga/mutsig_download
# Download MutSigCV 1.41 reference files from http://archive.broadinstitute.org/cancer/cga/mutsig_run#reference_files 
# Add the following lines within MutSigCV.m and run
# Outputs *.gene_rates.txt and *.overall_rates.txt, which we use to calculate mutation rate within genes
###

# %ADDED AT LINE 953 (MutSigCV.m):
#   G_rates = rmfield(G,{'expr','reptime','hic'});
# G_rates.r_x_X = G_rates.x ./ G_rates.X;
# G_rates.r_n_N = (G_rates.n_silent + G_rates.n_noncoding + G_rates.n_nonsilent) ./ (G_rates.N_silent + G_rates.N_noncoding + G_rates.N_nonsilent);
# G_rates.r_nns_Nns = (G_rates.n_nonsilent) ./ (G_rates.N_nonsilent);
# G_rates.r_ns_Ns = (G_rates.n_silent) ./ (G_rates.N_silent);
# G_rates.r_nnc_Nnc = (G_rates.n_noncoding) ./ (G_rates.N_noncoding);
# save_struct(G_rates, 'gene_rates.txt');
# %save(rate_file,'G_rates')
# 
# r_all = sum(G.n_silent + G.n_noncoding + G.n_nonsilent)/sum(G.N_silent + G.N_noncoding + G.N_nonsilent);
# r_silent = sum(G.n_silent)/sum(G.N_silent);
# r_nonsilent = sum(G.n_nonsilent)/sum(G.N_nonsilent);
# r_silent_noncoding = sum(G.n_silent + G.n_noncoding)/sum(G.N_silent + G.N_noncoding);
# r_x_X_min = min(G_rates.r_x_X(G_rates.r_x_X>0));
# r_x_X_max = max(G_rates.r_x_X);
# 
# rateId = fopen('overall_rates.txt', 'w');
# fprintf(rateId,'Overall mutation rate:\t%e\n',r_all);
# fprintf(rateId,'Silent mutation rate:\t%e\n',r_silent);
# fprintf(rateId,'Nonsilent mutation rate:\t%e\n',r_nonsilent);
# fprintf(rateId,'Silent+noncoding rate:\t%e\n',r_silent_noncoding);
# fprintf(rateId,'Min non-zero x/X:\t%e\n',r_x_X_min);
# fprintf(rateId,'Max x/X:\t%e\n',r_x_X_max);
# fclose(rateId);



# Calculating mutation rates and selection intensity for KRAS G12C  ---- 

require(BiocInstaller)
require(ggplot2)
require(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg19)


source("R/selection_intensity_calculation.R")
source("R/flip_function.R")
source("R/lamba_calculation.R")

load("output_data/trinuc_data_LUAD.RData")
load("output_data/Isoforms_LUAD.RData")
load("output_data/MAF_LUAD.RData")
LUAD.maf <- MAF_for_analysis

KRAS.mutation.output <- selection.intensity.calculation.function(genes_for_analysis = c("KRAS"),
                                                                 MAF_for_analysis = LUAD.maf,
                                                                 this.substitution = c("KRAS",34,"T"),
                                                                 trinuc.mutation_data = trinuc.mutation_data,
                                                                 LabReference =  isoforms,
                                                                 translations =  read.csv(file = "input_data/translations.csv",header = T,stringsAsFactors = F),
                                                                 mut_rates = read.csv(file="MutSigCV/MutSigCV_1.41/gene_rates.txt",header = T,stringsAsFactors = F,sep="\t"),
                                                                 low.mut = read.csv(file="MutSigCV/MutSigCV_1.41/overall_rates.txt",header = F,stringsAsFactors = F,sep="\t"),tumor.number = length(unique(LUAD.maf$Unique_patient_identifier)),mutsig_siggenes = read.csv(file="MutSigCV/MutSigCV_1.41/mutsig_output/.sig_genes.txt",header = T,stringsAsFactors = F,sep="\t"))


save(KRAS.mutation.output,file = "output_data/KRAS_mutation_output.RData")


# Figure: nucleotide level mutation rates ---- 


mutation.rates <- read.csv(file="MutSigCV/MutSigCV_1.41/gene_rates.txt",header = T,stringsAsFactors = F,sep="\t")
message("Mutation rate of KRAS:");print(mutation.rates$r_x_X[which(mutation.rates$gene=="KRAS")])



# run previous section, or
load("output_data/KRAS_mutation_output.RData")

require(reshape2)

source("R/fancy_scientific_code.R") 



normalized.mut.matrix.df <- as.data.frame(matrix(data=0,nrow=ncol(KRAS.mutation.output$norm_mut),ncol=5))
colnames(normalized.mut.matrix.df) <- c("Position","A","C","G","T")
normalized.mut.matrix.df$Position <- colnames(KRAS.mutation.output$norm_mut)
normalized.mut.matrix.df$A <- as.numeric(KRAS.mutation.output$norm_mut["A",])
normalized.mut.matrix.df$C <- as.numeric(KRAS.mutation.output$norm_mut["C",])
normalized.mut.matrix.df$G <- as.numeric(KRAS.mutation.output$norm_mut["G",])
normalized.mut.matrix.df$T <- as.numeric(KRAS.mutation.output$norm_mut["T",])
normalized.mut.matrix.df$Pos <- 1:nrow(normalized.mut.matrix.df)
normalized.mut.matrix.df$Ref <- KRAS.mutation.output$myseqsplit
head(normalized.mut.matrix.df)

normalized.mut.matrix.df.m <- melt(normalized.mut.matrix.df,id.vars = c("Position","Ref","Pos"))

head(normalized.mut.matrix.df.m)
normalized.mut.matrix.df.m <- normalized.mut.matrix.df.m[order(normalized.mut.matrix.df.m$Pos),]

positions.to.plot <- seq(from=34,to=36,by=1)
p <- ggplot(normalized.mut.matrix.df.m[normalized.mut.matrix.df.m$Pos %in% positions.to.plot,],aes(Pos, value, fill=variable)) 
p <- p + geom_bar(stat="identity",position="dodge",color="black") 
# p <- p + theme_bw() 
p <- p + theme(panel.background = element_blank())
p <- p + theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey"))
p <- p + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
p <- p + scale_y_continuous(labels=fancy_scientific)
p <- p + geom_text(aes(label=round(value*1e6,2)), vjust=-0.2, position=position_dodge(width=0.9),size=5)
p <- p + geom_text(aes(label=variable,y=-1e-7),position=position_dodge(width = 0.9),size=5)
p <- p + xlab("Nucleotide position")
p <- p + ylab("Mutation rate")
p <- p + scale_fill_discrete(name="Mutation")
p <- p + theme(axis.text.x= element_text(size=15))
p <- p + theme(axis.text.y= element_text(size=15))
p <- p + theme(axis.title.y = element_text(size = 18, angle = 90))
p <- p + theme(axis.title.x = element_text(size = 18, angle = 00))
p <- p + theme(legend.position = c(0.92, .85))
# p <- p + scale_x_discrete(labels=paste(Pos))
p

ggsave(p, file="figures/nuc_mutation_rates_KRASG12C.pdf",units = "in",height=7,width = 10)



# Mutation rates output and Figure: Amino acid mutation rates ----- 

# run previous section, or, 
load("output_data/KRAS_mutation_output.RData")

###
#Preserving the information from KRAS with the G12C mutation
###

just.12 <- KRAS.mutation.output$amino_acid_mutation_rates[,12]
ordered.12 <- sort(just.12,decreasing = T)[1:8] #Order it high --> low
ordered.12 <- ordered.12[c("Phe","Ser","Tyr","Arg","Trp","Gly","Cys","STOP")]#c(ordered.12[1:3],ordered.12[6],ordered.12[8],ordered.12[4],ordered.12[5],ordered.12[8]) #Change the ordering for the figure

message("Total mutation rate for oncogenic KRAS G12C mutations: ");print(sum(ordered.12[c("Phe","Ser","Tyr","Arg","Trp")])) #Total mutation rate 
message("Total mutation rate for all KRAS G12C mutations: ");print(sum(ordered.12)) #Total mutation rate 


KRAS.mutation.output$all_mutations[which(KRAS.mutation.output$all_mutations$AA_Pos==61 & KRAS.mutation.output$all_mutations$freq>1),]
message("Total rate of resistant mutations within KRAS") 
sum(ordered.12[c("Phe","Ser","Tyr","Arg","Trp")],KRAS.mutation.output$all_mutations$mu[which(KRAS.mutation.output$all_mutations$AA_Pos==61 & KRAS.mutation.output$all_mutations$freq>1)])

message("Sum of position 61:")
sum(KRAS.mutation.output$all_mutations$mu[which(KRAS.mutation.output$all_mutations$AA_Pos==61 &KRAS.mutation.output$all_mutations$freq>1)])

# ordered.12

ordered.12.df <- as.data.frame(matrix(NA,nrow=length(ordered.12),ncol=2))
colnames(ordered.12.df) <- c("Mutation","Rate")
ordered.12.df$Mutation <- names(ordered.12)
ordered.12.df$Rate <- as.numeric(ordered.12)
ordered.12.df
ordered.12.df$Mutation <- factor(ordered.12.df$Mutation, levels = c("Phe","Ser","Tyr","Arg","Trp","STOP","Cys","Gly"))#ordered.12.df$Mutation) #this makes Mutation an ordered factor, which ggplot uses. (http://stackoverflow.com/questions/20041136/avoid-ggplot-sorting-the-x-axis-while-plotting-geom-bar) 
ordered.12.df$alpha_vals <- c(rep(1,5),rep(1,3))
ordered.12.df$color_vals <- c(rep("black",5),rep("white",3))
# scale_fill_manual(values=cbPalette)


AA_mut_fig <- ggplot(ordered.12.df,aes(x=Mutation, y=Rate, fill=Mutation))
AA_mut_fig <- AA_mut_fig + geom_bar(stat="identity",position="dodge",color="black")#,aes(alpha=factor(alpha_vals)))
AA_mut_fig <- AA_mut_fig + theme(panel.background = element_blank())
AA_mut_fig <- AA_mut_fig + scale_fill_manual(values=ordered.12.df$color_vals)
AA_mut_fig <- AA_mut_fig + theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey"))
AA_mut_fig <- AA_mut_fig + theme(panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())
AA_mut_fig <- AA_mut_fig + scale_y_continuous(labels=fancy_scientific)
AA_mut_fig <- AA_mut_fig + geom_text(aes(label=round(Rate*1e6,2)), vjust=-0.2, position=position_dodge(width=0.9),size=5)
AA_mut_fig <- AA_mut_fig + xlab("Amino acid mutation")
AA_mut_fig <- AA_mut_fig + ylab("Mutation rate")
AA_mut_fig <- AA_mut_fig + theme(axis.text.x= element_text(size=15))
AA_mut_fig <- AA_mut_fig + theme(axis.text.y= element_text(size=15))
AA_mut_fig <- AA_mut_fig + theme(axis.title.y = element_text(size = 18, angle = 90))
AA_mut_fig <- AA_mut_fig + theme(axis.title.x = element_text(size = 18, angle = 00))
AA_mut_fig <- AA_mut_fig + theme(legend.position = 'none')
AA_mut_fig


ggsave(AA_mut_fig, file="figures/AA_mutation_rates.pdf",units = "in",height=7,width = 10)



# Calculate mutation rates and selection intensity of all substitutions of interest ----  


require(BiocInstaller)
require(ggplot2)
require(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg19)


source("R/selection_intensity_calculation.R")
source("R/flip_function.R")
source("R/lamba_calculation.R")

load("output_data/trinuc_data_LUAD.RData")
load("output_data/Isoforms_LUAD.RData")
load("output_data/MAF_LUAD.RData")
LUAD.maf <- MAF_for_analysis

Downstream.mutation.output <- selection.intensity.calculation.function(genes_for_analysis = c("HRAS",
                                                                                              "NRAS",
                                                                                              "PIK3CA",
                                                                                              "BRAF",
                                                                                              "MAP2K1",
                                                                                              "MAP2K2",
                                                                                              "MAPK1",
                                                                                              "MAPK3",
                                                                                              "AKT1",
                                                                                              "MTOR"),
                                                                       MAF_for_analysis = LUAD.maf,
                                                                       this.substitution = c("KRAS",34,"T"),
                                                                       trinuc.mutation_data = trinuc.mutation_data,
                                                                       LabReference =  isoforms,
                                                                       translations =  read.csv(file = "input_data/translations.csv",header = T,stringsAsFactors = F),
                                                                       mut_rates = read.csv(file="MutSigCV/MutSigCV_1.41/gene_rates.txt",header = T,stringsAsFactors = F,sep="\t"),
                                                                       low.mut = read.csv(file="MutSigCV/MutSigCV_1.41/overall_rates.txt",header = F,stringsAsFactors = F,sep="\t"),tumor.number = length(unique(LUAD.maf$Unique_patient_identifier)),mutsig_siggenes = read.csv(file="MutSigCV/MutSigCV_1.41/mutsig_output/.sig_genes.txt",header = T,stringsAsFactors = F,sep="\t"))

save(Downstream.mutation.output,file = "output_data/downstream_mut_data.RData")


plot(y=Downstream.mutation.output$all_mutations$freq[which(Downstream.mutation.output$all_mutations$freq>1)],x=Downstream.mutation.output$all_mutations$mu[which(Downstream.mutation.output$all_mutations$freq>1)],ylab="Frequency",xlab="Mutation rate")
abline(lm(Downstream.mutation.output$all_mutations$freq[which(Downstream.mutation.output$all_mutations$freq>1)] ~ Downstream.mutation.output$all_mutations$mu[which(Downstream.mutation.output$all_mutations$freq>1)]),lwd=2,col="red")


down.lm <- lm(Downstream.mutation.output$all_mutations$freq[which(Downstream.mutation.output$all_mutations$freq>1)] ~ Downstream.mutation.output$all_mutations$mu[which(Downstream.mutation.output$all_mutations$freq>1)])
summary(down.lm)$r.squared

message("Correlation test for Frequency vs. Mutation rate")
cor.test(Downstream.mutation.output$all_mutations$mu[which(Downstream.mutation.output$all_mutations$freq>1)],Downstream.mutation.output$all_mutations$freq[which(Downstream.mutation.output$all_mutations$freq>1)])


sum(Downstream.mutation.output$all_mutations$mu[which(Downstream.mutation.output$all_mutations$freq>1)])

# Figure: mutation rate within KRAS and downstream ---- 

source("R/fancy_scientific_code.R")
require(ggplot2)


###
#For poster, a figure showing mutation rates within and outside of KRAS 
###

# ordered.12.df <- ordered.12.df[c(1:4,8,6,7,5),] #if the order was wrong
colnames(ordered.12.df)[2] <- "mu"
rates.data.frame <- as.data.frame(matrix(data = NA, nrow=8, ncol=6))
colnames(rates.data.frame) <- c("Gene","AA_Ref","AA_Pos","AA_Change","mu","freq")

rates.data.frame[,"mu"] <- ordered.12.df[,"mu"]
rates.data.frame[,"Gene"] <- "KRAS"
rates.data.frame[,"AA_Pos"] <- 12
rates.data.frame[,"AA_Ref"] <- "C"
rates.data.frame[,"AA_Change"] <- as.vector(ordered.12.df[,"Mutation"])

translations =  read.csv(file = "input_data/translations.csv",header = T,stringsAsFactors = F)
# levels(rates.data.frame) <- droplevels(rates.data.frame)#levels(droplevels(rates.data.frame$Change))
for(i in 1:nrow(rates.data.frame)){
  rates.data.frame[i,"AA_Change"] <- translations[which(translations$AA_short==rates.data.frame[i,"AA_Change"])[1],"AA_letter"]
}

rates.data.frame <- rbind(rates.data.frame,KRAS.mutation.output$all_mutations[which(KRAS.mutation.output$all_mutations$AA_Pos==61 & KRAS.mutation.output$all_mutations$freq>1),c("Gene","AA_Ref","AA_Pos","AA_Change","mu","freq")])

rates.data.frame <- rbind(rates.data.frame,Downstream.mutation.output$all_mutations[which(Downstream.mutation.output$all_mutations$freq>1),c("Gene","AA_Ref","AA_Pos","AA_Change","mu","freq")])

rates.data.frame$Name <- paste(rates.data.frame$Gene," ",rates.data.frame$AA_Ref,rates.data.frame$AA_Pos,rates.data.frame$AA_Change,sep="")

rates.data.frame$Name <- factor(rates.data.frame$Name, levels = unique(rates.data.frame$Name))
rates.data.frame$Name <- factor(rates.data.frame$Name, levels = rev(unique(rates.data.frame$Name)))
rates.data.frame$Gene <- factor(rates.data.frame$Gene, levels = unique(rates.data.frame$Gene))

rates.data.frame$alpha_vals <- c(rep(1,5),rep(.8,3),rep(1,15))

rates.data.frame$Gene <- factor(rates.data.frame$Gene, levels = c("BRAF","PIK3CA","KRAS","NRAS","MAP2K1","AKT1"))

mutation.plot <- ggplot(data = rates.data.frame, aes(x= Name, y=mu)) + geom_bar(aes(fill=Gene,color=Gene,alpha=alpha_vals),stat = "identity")#,aes(alpha=alpha_vals))
mutation.plot <- mutation.plot + geom_text(aes(label=freq, hjust=0),show.legend =FALSE,size=13)
mutation.plot <- mutation.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + scale_y_continuous(labels=fancy_scientific)
mutation.plot <- mutation.plot + theme(legend.position = c(0.8, 0.8)) + labs(y="Mutation rate",x="Mutation") + theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(panel.grid.major = element_line(colour = "lightgray"), panel.grid.minor = element_line(colour = "lightgray"))

mutation.plot <- mutation.plot + theme(axis.text.y=element_text(size=18),axis.text.x = element_text(size=25,angle = 0,hjust = 0.5), axis.title = element_text(size=25, face="bold")) + theme(legend.text=element_text(size=22)) 
mutation.plot <- mutation.plot + coord_flip() + theme(legend.position = c(0.5,0.15)) + scale_alpha(guide = 'none')
mutation.plot

ggsave(filename = "figures/mutation_rates_within_and_downstream.pdf",plot = mutation.plot)


# Figure: sum of mutation rates -----

plot.sum.df <- as.data.frame(matrix(data=NA,nrow=2,ncol=2))
colnames(plot.sum.df) <- c("Name","mu")

plot.sum.df[1,"Name"] <- "Within KRAS"
plot.sum.df[2,"Name"] <- "Downstream of KRAS"
plot.sum.df$Name <- factor(plot.sum.df$Name, levels=c("Downstream of KRAS","Within KRAS"))
plot.sum.df[1,"mu"] <- sum(rates.data.frame[1:5,"mu"],rates.data.frame[9:10,"mu"])
plot.sum.df[2,"mu"] <- sum(rates.data.frame[11:23,"mu"])

plot.sum <- ggplot(data = plot.sum.df,aes(x=Name,y=mu)) + geom_bar(stat = "identity") + scale_y_continuous(labels=fancy_scientific) + labs(y="Mutation rate from tumorigenesis to resection",x="Mutations") + theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black")) + theme(panel.grid.major = element_line(colour = "lightgray"), panel.grid.minor = element_line(colour = "lightgray")) + theme(axis.text.y=element_text(size=20,angle=90,hjust = 0.5), axis.title = element_text(size=20), axis.text.x=element_text(size=30))
plot.sum <- plot.sum + coord_flip()
plot.sum

ggsave(filename = "figures/mutation_rates_sum.pdf",plot = plot.sum)



# Selection intensity analysis ---- 

load("output_data/downstream_mut_data.RData")
load("output_data/KRAS_mutation_output.RData")

downstream_mutations <- Downstream.mutation.output$all_mutations


require(BiocInstaller)
require(ggplot2)
require(BSgenome)
require(BSgenome.Hsapiens.UCSC.hg19)


source("R/selection_intensity_calculation.R")
source("R/flip_function.R")
source("R/lamba_calculation.R")

load("output_data/trinuc_data_LUAD.RData")
load("output_data/Isoforms_LUAD.RData")
load("output_data/MAF_LUAD.RData")
LUAD.maf <- MAF_for_analysis

KRAS.no.mut.output <- selection.intensity.calculation.function(genes_for_analysis = c("KRAS"),
                                                               MAF_for_analysis = LUAD.maf,
                                                               this.substitution = c("not_a_gene",34,"T"),
                                                               trinuc.mutation_data = trinuc.mutation_data,
                                                               LabReference =  isoforms,
                                                               translations =  read.csv(file = "input_data/translations.csv",header = T,stringsAsFactors = F),
                                                               mut_rates = read.csv(file="MutSigCV/MutSigCV_1.41/gene_rates.txt",header = T,stringsAsFactors = F,sep="\t"),
                                                               low.mut = read.csv(file="MutSigCV/MutSigCV_1.41/overall_rates.txt",header = F,stringsAsFactors = F,sep="\t"),tumor.number = length(unique(LUAD.maf$Unique_patient_identifier)),mutsig_siggenes = read.csv(file="MutSigCV/MutSigCV_1.41/mutsig_output/.sig_genes.txt",header = T,stringsAsFactors = F,sep="\t"))


save(KRAS.no.mut.output,file = "output_data/KRAS_NO_mutation_output.RData")


load("output_data/KRAS_NO_mutation_output.RData")
targeted_mut <- c("KRAS","G",12,"C")

mutation.data <- KRAS.no.mut.output$complete_mutation_data

KRAS.muts_from_analysis <- mutation.data[which(mutation.data$Gene==targeted_mut[1] & 
                                                 mutation.data$Amino_acid_reference==targeted_mut[2] &
                                                 mutation.data$Amino_acid_position==as.numeric(targeted_mut[3]) &
                                                 mutation.data$Amino_acid_alternative==targeted_mut[4])[1],]







#Without G12C
# selection.subset <- rbind(downstream_mutations[downstream_mutations$freq>1,],KRAS.no.mut.output$all_mutations[c(2,11,12),])

# Pos 61 KRAS included, from the KRAS mutation since assume G12C is already present
selection.subset <- rbind(downstream_mutations[downstream_mutations$freq>1,],KRAS.mutation.output$all_mutations[which(KRAS.mutation.output$all_mutations$AA_Pos==61),])

# KRAS.no.mut.output$all_mutations

selection.subset

selection.subset.ordered <- selection.subset[order(selection.subset$gamma_epistasis),]
selection.subset.ordered$Name <- paste(selection.subset.ordered$Gene," ",selection.subset.ordered$AA_Ref,selection.subset.ordered$AA_Pos,selection.subset.ordered$AA_Change,sep="")
selection.subset.ordered$Name <- factor(selection.subset.ordered$Name, levels=selection.subset.ordered$Name)


selection.subset.ordered$Gene <- factor(selection.subset.ordered$Gene, levels = c("BRAF","PIK3CA","KRAS","NRAS","MAP2K1","AKT1"))
#from http://stackoverflow.com/questions/18265941/two-horizontal-bar-charts-with-shared-axis-in-ggplot2-similar-to-population-pyr
library('grid')
library('gridExtra')
source("R/fancy_scientific_code.R")

g.mid <- ggplot(selection.subset.ordered,aes(x=1,y=Name)) +
  geom_text(aes(label=Name),size=5, fontface = "bold") +
  # geom_segment(aes(x=0.94,xend=0.96,yend=Name)) +
  # geom_segment(aes(x=1.04,xend=1.065,yend=Name)) +
  labs(title=" ",subtitle=" ") +# ggtitle("") +
  ylab(NULL) + 
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065)) +
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA,size=12),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm")) + theme(plot.title = element_text(hjust = 0, face="bold",size=30),
                                                        plot.subtitle = element_text(size=25,hjust=0.5))
g.mid

g1 <- ggplot(data=selection.subset.ordered,aes(x=Name,y=mu,fill=Gene)) +
  geom_bar(stat="identity") + labs(title="A",subtitle="Mutation rate") +#ggtitle("Mutation rate from tumorigenesis to resection") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        plot.margin = unit(c(1,-1,1,10), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  geom_text(aes(label=round(mu*1e6,2)), vjust=-0.5, hjust=0.5, position=position_dodge(width=0.9),size=4,angle=90) +
  geom_text(aes(label=freq,y=-0.5e-7), position=position_dodge(width=0),size=5,angle=0,color="black") +
  theme(plot.title = element_text(hjust = 0, face="bold",size=30),
        plot.subtitle = element_text(size=25,hjust=0.5)) +
  theme(legend.position = 'none',axis.text.x = element_text(size=12)) +
  scale_y_reverse(labels=fancy_scientific) + coord_flip() 
g1

# levels(selection.subset.ordered$Gene)

g2 <- ggplot(data=selection.subset.ordered, aes(x=Name,y=gamma_epistasis,fill=Gene)) +
  geom_bar(stat="identity") + labs(title="B",subtitle="Selection intensity") +# ggtitle("Selection intensity") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = unit(c(1,0,1,-1), "mm")) +
  theme(panel.background = element_blank()) +
  theme(panel.grid.major =element_line(color="lightgrey"),panel.grid.minor =element_line(color="lightgrey")) +
  theme(panel.grid.major.y = element_blank(),panel.grid.minor.y = element_blank()) +
  theme(plot.title = element_text(hjust = 0, face="bold",size=30),
        plot.subtitle = element_text(size=25,hjust=0.5),
        axis.text.x = element_text(size=12)) +
  coord_flip() +
  theme(legend.position = c(0.75, .5),legend.text = element_text(size=12))
g2

gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4.25/10,1/10,4.25/10))


gg.combined <- arrangeGrob(gg1,gg.mid,gg2,ncol=3,widths=c(4/10,1.5/10,4/10))
ggsave(gg.combined, filename = "figures/combined_mu_selection_plot_titles.pdf",units = "in",height=7,width = 10)





# Checking if there are multiple KRAS mutations or downstream mutations within seqeuencing data----

source("R/merging_TCGA_and_local_MAF.R")
source("R/hg38_to_hg19_converter.R")
require(rtracklayer)
source("R/unique_tumor_addition.R")
source("R/DNP_remover.R")
source("R/flip_function.R")
source("R/tumor_allele_add.R")
# source("R/selection_intensity_calculation.R")

# Download PAAD TCGA data: 
# File Properties
# Name	TCGA.PAAD.mutect.a2c08388-2be9-4124-83ba-ed8d686dd277.DR-7.0.somatic.maf.gz
# Access	open
# UUID	a2c08388-2be9-4124-83ba-ed8d686dd277

# Obtain Yale-Gilead PAAD data

# Download COAD TCGA data:
# File Properties
# Name	TCGA.COAD.mutect.853e5584-b8c3-4836-9bda-6e7e84a64d97.DR-7.0.somatic.maf.gz
# Access	open
# UUID	853e5584-b8c3-4836-9bda-6e7e84a64d97

# Download READ TCGA data: 
# File Properties
# Name	TCGA.READ.mutect.c999f6ca-0b24-4131-bc53-1665948f8e3f.DR-7.0.somatic.maf.gz
# Access	open
# UUID	c999f6ca-0b24-4131-bc53-1665948f8e3f

# Load in and process all data

PAAD.TCGA <- read.csv(file = "input_data/PAAD/gdc_download_20170925_182341/a2c08388-2be9-4124-83ba-ed8d686dd277/TCGA.PAAD.mutect.a2c08388-2be9-4124-83ba-ed8d686dd277.DR-7.0.somatic.maf", skip=5, sep="\t",stringsAsFactors = F,header = T)
PAAD.TCGA <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=PAAD.TCGA)
PAAD.YG <- read.csv(file = "input_data/PAAD/mutationsTN_26_Pancreatic_Cancer_YG_data.maf",stringsAsFactors = F,header = T,sep = "\t")

PAAD.TCGA.YG <- merging.TCGA.and.Yale.MAF.data.function(NCI_data = PAAD.TCGA,Yale_data = PAAD.YG)
PAAD.TCGA.YG <- unique.tumor.addition.function(MAF.file = PAAD.TCGA.YG,non.TCGA.characters.to.keep = 'all')
PAAD.TCGA.YG <- tumor.allele.adder(MAF = PAAD.TCGA.YG)




COAD.TCGA <- read.csv(file="input_data/COAD/gdc_download_20170925_182443/853e5584-b8c3-4836-9bda-6e7e84a64d97/TCGA.COAD.mutect.853e5584-b8c3-4836-9bda-6e7e84a64d97.DR-7.0.somatic.maf",header = T,sep = "\t",stringsAsFactors = F,skip=5)
COAD.TCGA <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=COAD.TCGA)
COAD.TCGA <- unique.tumor.addition.function(MAF.file = COAD.TCGA,non.TCGA.characters.to.keep = 'all')
COAD.TCGA <- tumor.allele.adder(MAF = COAD.TCGA)



READ.TCGA <- read.csv(file="input_data/READ/gdc_download_20170925_182528/c999f6ca-0b24-4131-bc53-1665948f8e3f/TCGA.READ.mutect.c999f6ca-0b24-4131-bc53-1665948f8e3f.DR-7.0.somatic.maf",header = T,sep = "\t",stringsAsFactors = F,skip=5)
READ.TCGA <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=READ.TCGA)
READ.TCGA <- unique.tumor.addition.function(MAF.file = READ.TCGA,non.TCGA.characters.to.keep = 'all')
READ.TCGA <- tumor.allele.adder(MAF = READ.TCGA)




LUAD.TCGA <- read.csv('input_data/NCI/gdc_download_20170919_153409/81ccaef3-4550-494d-882c-895fb5a3de3b/TCGA.LUAD.mutect.81ccaef3-4550-494d-882c-895fb5a3de3b.DR-7.0.somatic.maf',stringsAsFactors=F,skip=5,header=T,sep='\t')
LUAD.TCGA <- hg38.to.hg19.converter(chain='input_data/hg38Tohg19.chain',hg38_maf=LUAD.TCGA)
LUAD.YG <- read.csv('input_data/YG/adc_inc_counts.txt',stringsAsFactors=F,header=T,sep='\t')
if(length(which(colnames(LUAD.YG)=="keep"))>0){
  LUAD.YG <- LUAD.YG[-which(LUAD.YG$keep==0),]
}
LUAD.TCGA.YG <- merging.TCGA.and.Yale.MAF.data.function(NCI_data = LUAD.TCGA,
                                                        Yale_data = LUAD.YG)
LUAD.TCGA.YG <- unique.tumor.addition.function(MAF.file = LUAD.TCGA.YG,non.TCGA.characters.to.keep = 'all')
LUAD.TCGA.YG <- tumor.allele.adder(MAF = LUAD.TCGA.YG)


save(PAAD.TCGA.YG,
     COAD.TCGA,
     READ.TCGA,
     LUAD.TCGA.YG,
     file="output_data/PAAD_COAD_READ_LUAD_sequencing_data.RData")

# 21


# Importing data in Table 1. 

other.genes.Table1 <- c("KRAS",
                        "ABL1",
                        "AKT2",
                        "BRAF",
                        "EGFR",
                        "FGFR1",
                        "FGFR3",
                        "HRAS",
                        "JAK2",
                        "MET",
                        "NRAS",
                        "PDGFRA",
                        "PIK3CA")
table.1.matrix <- matrix(data=NA,nrow=18,ncol=2)

table.1.matrix[,1] <- c(249,
                        16,
                        598,
                        757,
                        123,
                        363,
                        59,
                        9,
                        607,
                        137,
                        55,
                        4,
                        1245,
                        9,
                        55,
                        560,
                        540,
                        1038)

table.1.matrix[,2] <- c(258,
                        26,
                        615,
                        761,
                        136,
                        374,
                        76,
                        20,
                        618,
                        148,
                        65,
                        15,
                        1256,
                        20,
                        67,
                        572,
                        551,
                        1049)

rownames(table.1.matrix) <- c("ABL1",
                              "AKT2",
                              "BRAF",
                              "EGFR",
                              "FGFR1",
                              "FGFR3",
                              "HRAS",
                              "HRAS",
                              "JAK2",
                              "KRAS",
                              "KRAS",
                              "KRAS",
                              "MET",
                              "NRAS",
                              "NRAS",
                              "PDGFRA",
                              "PIK3CA",
                              "PIK3CA")

colnames(table.1.matrix) <- c("Start","End")

# Function to find if Table 1 mutations are present with KRAS G12C

additional.mutations.function <- function(this.MAF,genes.to.check){
  this.tumor.file <- this.MAF
  additional.muts <- matrix(data = NA,nrow=1,ncol=ncol(this.tumor.file))
  colnames(additional.muts) <- colnames(this.tumor.file)
  for(i in 1:length(unique(this.tumor.file$Unique_patient_identifier))){
    this.tumor <-  this.tumor.file[which(this.tumor.file$Unique_patient_identifier==unique(this.tumor.file$Unique_patient_identifier)[i]),]
    if(length(which(this.tumor$Hugo_Symbol=="KRAS" & 
                    this.tumor$Start_Position==25398285  & 
                    this.tumor$Chromosome==12 &
                    this.tumor$Tumor_allele=="A"))>0){ #This is the mutation that results in KRAS G12C
      
      #Skip if this was actually a DNP (and not KRAS G12C) 
      if(length(which(this.tumor$Hugo_Symbol=="KRAS" & 
                      this.tumor$Start_Position==25398284 & 
                      this.tumor$Chromosome==12))>0){
        next
      }
      # Need to store if any downstream or within-KRAS mutations
      # If there is an additional mutation besides KRAS G12C... 
      if(length(which(this.tumor$Hugo_Symbol %in% genes.to.check))>1){
        additional.muts <- rbind(additional.muts,this.tumor[which(this.tumor$Hugo_Symbol %in% genes.to.check),])
        
      }
    }
    
  }
  if(nrow(additional.muts)>1){
    additional.muts <- (additional.muts[which(additional.muts$Variant_Classification == "Missense_Mutation"),])
  }
  return(additional.muts[-1,])# get rid of initiating row 
}

# Check for relevant tumor types

additional.muts.LUAD <- additional.mutations.function(this.MAF = LUAD.TCGA.YG,genes.to.check = other.genes.Table1)
additional.muts.PAAD <- additional.mutations.function(this.MAF = PAAD.TCGA.YG,genes.to.check = other.genes.Table1)
additional.muts.COAD <- additional.mutations.function(this.MAF = COAD.TCGA,genes.to.check = other.genes.Table1)
additional.muts.READ <- additional.mutations.function(this.MAF = READ.TCGA,genes.to.check = other.genes.Table1)


# unique(additional.muts.LUAD$Unique_patient_identifier)

# Distil into easy to read lists

LUAD.muts <- list()
for(i in 1:length(unique(additional.muts.LUAD$Unique_patient_identifier))){
  LUAD.muts[[i]] <- additional.muts.LUAD[which(additional.muts.LUAD$Unique_patient_identifier==unique(additional.muts.LUAD$Unique_patient_identifier)[i]),c("Unique_patient_identifier","Hugo_Symbol","HGVSc","HGVSp","HGVSp_Short","Variant_Classification","Variant_Type")] 
  to.delete <- NULL
  for(j in 1:nrow(LUAD.muts[[i]])){
    this.aa.pos <- as.numeric(unlist(gsub("[^0-9]", "", unlist(LUAD.muts[[i]][j,]$HGVSp_Short)), ""))
    table.1.rows <- which(rownames(table.1.matrix)==LUAD.muts[[i]][j,"Hugo_Symbol"])
    delete.this.round <- T
    for(k in 1:length(table.1.rows)){
      if((table.1.matrix[table.1.rows[k],"Start"] <= this.aa.pos &
          table.1.matrix[table.1.rows[k],"End"] >= this.aa.pos)){
        delete.this.round <- F
      }
    }
    if(delete.this.round){to.delete <- c(to.delete,j)}
  }
  
  if(length(to.delete)>0){LUAD.muts[[i]] <- LUAD.muts[[i]][-to.delete,]}
  
}
LUAD.muts


COAD.muts <- list()
for(i in 1:length(unique(additional.muts.COAD$Unique_patient_identifier))){
  COAD.muts[[i]] <- additional.muts.COAD[which(additional.muts.COAD$Unique_patient_identifier==unique(additional.muts.COAD$Unique_patient_identifier)[i]),c("Unique_patient_identifier","Hugo_Symbol","HGVSc","HGVSp","HGVSp_Short","Variant_Classification","Variant_Type")] 
  to.delete <- NULL
  for(j in 1:nrow(COAD.muts[[i]])){
    this.aa.pos <- as.numeric(unlist(gsub("[^0-9]", "", unlist(COAD.muts[[i]][j,]$HGVSp_Short)), ""))
    table.1.rows <- which(rownames(table.1.matrix)==COAD.muts[[i]][j,"Hugo_Symbol"])
    delete.this.round <- T
    for(k in 1:length(table.1.rows)){
      if((table.1.matrix[table.1.rows[k],"Start"] <= this.aa.pos &
          table.1.matrix[table.1.rows[k],"End"] >= this.aa.pos)){
        delete.this.round <- F
      }
    }
    if(delete.this.round){to.delete <- c(to.delete,j)}
  }
  
  if(length(to.delete)>0){COAD.muts[[i]] <- COAD.muts[[i]][-to.delete,]}
  
}
COAD.muts

# Number of tumors within each dataset

load("output_data/PAAD_COAD_READ_LUAD_sequencing_data.RData")

length(unique(LUAD.TCGA.YG$Unique_patient_identifier))
length(grep(pattern = "TCGA-",x = unique(LUAD.TCGA.YG$Unique_patient_identifier)))
length(unique(LUAD.TCGA.YG$Unique_patient_identifier)) - length(grep(pattern = "TCGA-",x = unique(LUAD.TCGA.YG$Unique_patient_identifier)))

length(unique(PAAD.TCGA.YG$Unique_patient_identifier))
length(grep(pattern = "TCGA-",x = unique(PAAD.TCGA.YG$Unique_patient_identifier)))
length(unique(PAAD.TCGA.YG$Unique_patient_identifier)) - length(grep(pattern = "TCGA-",x = unique(PAAD.TCGA.YG$Unique_patient_identifier)))

length(unique(COAD.TCGA$Unique_patient_identifier))

length(unique(READ.TCGA$Unique_patient_identifier))


