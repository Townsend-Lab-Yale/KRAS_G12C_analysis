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

require("rtracklayer")
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

















