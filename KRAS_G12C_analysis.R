####
# Analysis performed in the manuscript Cannataro et al., The likelihood of heterogeneity or additional mutation in KRAS or associated oncogenes to compromise targeting of oncogenic KRAS G12C
###


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

Yale_data <- read.csv('input_data/YG/mutationsTN_22_Lung_adenocarcinoma.maf',stringsAsFactors=F,header=T,sep='\t')
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

MAF_for_analysis <- synonym.matcher.function(input_to_be_changed = MAF_for_analysis,synonym.table = read.csv(file='input_data/matching_table.txt',stringsAsFactors = F,sep='\t'),covariates.file = read.table(file='input_data/LUNG.txt',sep='\t',header = T,stringsAsFactors = F),isoforms_or_MAF = 'MAF')
MAF_for_analysis <- tumor.allele.adder(MAF = MAF_for_analysis)

isoforms <- synonym.matcher.function(input_to_be_changed = read.csv(file='input_data/UCSC_refseq_refFlat_hg19_short.txt',header = T,sep='\t',stringsAsFactors = F),synonym.table = read.csv(file='input_data/matching_table.txt',stringsAsFactors = F,sep='\t'),covariates.file = read.table(file='input_data/LUNG.txt',sep='\t',header = T,stringsAsFactors = F),isoforms_or_MAF = 'isoforms')

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


# Calculating mutation rates with MutSigCV ----

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
p <- p + ylab("Mutation rate from tumorigenesis to resection")
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



# Figure: 












