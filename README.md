# KRAS_G12C_analysis
Analysis within The likelihood of heterogeneity or additional mutation in KRAS or associated oncogenes to compromise targeting of oncogenic KRAS G12C


# Guide for use: 

All analyses are contained within KRAS_G12C_analysis.R. The different sections of the script are split with "----", providing section breaks in RStudio. Detailed instructions and comments are provided within this script, however, you may need to double check the filenames to make sure they match your directory structure. Briefly, to reproduce the results within the manuscript: 

1) Download and unzip the LUAD .maf file from https://portal.gdc.cancer.gov/ and store it within the input data folder. File details: 

```
File Properties
Name:	TCGA.LUAD.mutect.81ccaef3-4550-494d-882c-895fb5a3de3b.DR-7.0.somatic.maf.gz
Access:	open
UUID:	81ccaef3-4550-494d-882c-895fb5a3de3b
```
2) Convert this to hg19, remove potential DNP listed as SNP, and process for analysis. 

3) Obtain Yale-Gilead dataset, save it within input_data, combine with NCI data to a single MAF file. 

4) Match gene names and isoform names to gene expression files used in MutSigCV. 

5) Calculate average trinucleotide context among the tumors.

6) Agree to the MutSigCV license and download MutSigCV 1.41 from http://archive.broadinstitute.org/cancer/cga/mutsig_download. Download MutSigCV 1.41 reference files from http://archive.broadinstitute.org/cancer/cga/mutsig_run#reference_files. Save files to the MutSigCV folder and unzip. Add contents of input_data/to_add_MutSigCV.txt to line 953 of MutSigCV.m. On my computer, I changed the MutSigCV.m file and function to MutSigCV_rates.m and MutSigCV(), and used this script, located in the same directory as MutSigCV_rates.m, to load in the input. We use LUNG_expression.txt as the covariate file, which is the averaged expression of all LUNG cell lines downloaded from The Cancer Cell Line Encyclopedia (Barretina, Jordi, et al. "The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity." Nature 483.7391 (2012): 603-607.)

```
addpath('.../KRAS_G12C_analysis/');
mutfile = '.../KRAS_G12C_analysis/output_data/MAF_LUAD.txt';
p_cov = '.../KRAS_G12C_analysis/MutSigCV/exome_full192.coverage.txt';
p_covar = '.../KRAS_G12C_analysis/MutSigCV/LUNG_expression.txt';
p_dict = '.../KRAS_G12C_analysis/MutSigCV/mutation_type_dictionary_file.txt';
p_chr = '.../KRAS_G12C_analysis/MutSigCV/chr_files_hg19/';

stem = strcat('mutsig_output/KRAS_G12C_');
MutSigCV_rates(mutfile, p_cov, p_covar, stem, p_dict, p_chr)
```
7) Now all the input data is processed. Run the selection sections, and the figure generating portions. Double check lines that import and export data to make sure they match your directory structure! 


