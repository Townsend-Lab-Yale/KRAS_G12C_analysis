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















