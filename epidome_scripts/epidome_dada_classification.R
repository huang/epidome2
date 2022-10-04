source("https://raw.githubusercontent.com/ssi-dk/epidome/master/scripts/epidome_functions.R")
setwd("/Volumes/data/MPV/projects/git.repositories/epidome/")

### Load amplicon table for classification ###
ST_amplicon_table = read.table("DB/epidome_ST_amplicon_frequencies.txt",sep = "\t")

### Load dada2 output for the two primers ###
epi01_table = read.table("example_data/epi01_dada_output_article.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/epi02_dada_output_article.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("/Volumes/data/MPV/ANMC/Epidome_190920/EPI_dada2_output/EPI_G216_Bayes_190920_tax_and_counts_mb35.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("/Volumes/data/MPV/ANMC/Epidome_190920/EPI_dada2_output/EPI_yycH_Bayes_190920_tax_and_counts_mb35.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/190920_run1_G216_classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/190920_run1_yycH_classified.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/190920_run2_G216_classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/190920_run2_yycH_classified.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/190920_run2_G216_classified_997.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/190920_run2_yycH_classified_997.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/article_g216_reclassified.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/article_yycH_reclassified.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/article_G216_reclassified.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/article_yycH_reclassified.csv",sep = ";",header=TRUE,row.names=1)
epi03_table = read.table("example_data/190920_run2_g216_test.csv",sep = ";",header=TRUE,row.names=1)
epi04_table = read.table("example_data/190920_run2_yycH_test.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/article_g216_test.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/article_yycH_test.csv",sep = ";",header=TRUE,row.names=1)
epi01_table = read.table("example_data/190920_run1_G216_classified_995p.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/190920_run1_yycH_classified_995p.csv",sep = ";",header=TRUE,row.names=1)


### Load metadata table
metadata_table = read.table("example_data/article_metadata.txt",header=TRUE,row.names=1)
metadata_table = read.table("../epidome_backup/sample_metadata.txt",header=TRUE,row.names=1)

### Setup an object for easy handling of epidome data
epidome_object = setup_epidome_object(epi01_table,epi02_table,metadata_table = metadata_table)
#epidome_object = setup_epidome_object(epi01_table,epi02_table)

### Check if number of sequences from each primer for each samples match up approximately ###
compare_primer_output(epidome_object)
compare_primer_output(epidome_object,"sample.type")
compare_primer_output(epidome_object,"patient.sample.site")



#### Data manipulation of various sorts ####

### Combine ASVs from dada output ###
epidome_ASV_combined = combine_ASVs_epidome(epidome_object)

### Filter lowcount samples (removes any sample that has less than X sequences from one of the two primer sets, here 500) ###
epidome_filtered_samples = filter_lowcount_samples_epidome(epidome_object,500,500)


#### Plots and figures ###

### Setup a data frame with sequence classification based on the two primer outputs. Combine ASVs first ###
epidome_ASV_combined = combine_ASVs_epidome(epidome_filtered_samples)
count_table = classify_epidome(epidome_ASV_combined,ST_amplicon_table)
### Make barplot based on classification. Set reorder=TRUE to order samples based on Bray Curtis dissimilarity and/or set normalize=FALSE to not normalize to percent ###
p = make_barplot_epidome(count_table,reorder=TRUE,normalize=TRUE)
p

p2 = p

grid.arrange(p1,p2,ncol=2)


### Data stratification based on metadata variables
epidome_object_clinical = prune_by_variable_epidome(epidome_object,"sample.type",c("Clinical"))
epidome_object_mock = prune_by_variable_epidome(epidome_object,"sample.type",c("Mock community"))


### make barplot of "clinical" samples only 
epidome_object_clinical_ASV_combined = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Clinical"))
count_table_clinical = classify_epidome_3(epidome_object_clinical_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_clinical,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples")

### make barplot of mock samples only 
epidome_object_mock_ASV_combined = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Mock community"))
count_table_mock = classify_epidome(epidome_object_mock_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_mock,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples")

count_table_mock = classify_epidome_2(epidome_object_mock_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_mock,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples")

count_table_mock = classify_epidome_3(epidome_object_mock_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_mock,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples")


count_table_mock = classify_epidome_4(epidome_object_mock_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_mock,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples")

count_table_mock = classify_epidome_4(eo2,ST_amplicon_table)
p = make_barplot_epidome(count_table_mock,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples")


epidome_object_clinical_ASV_combined = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Clinical"))
count_table_clinical = classify_epidome_4(epidome_object_clinical_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_clinical,reorder=FALSE,normalize=TRUE)
p

epidome_object_clinical_ASV_combined = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Clinical"))
count_table_clinical = classify_epidome_3(epidome_object_clinical_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_clinical,reorder=FALSE,normalize=TRUE)
p


### make PCA plots of clinical samples. Color according to variable in metadata and (optional) indicate colors to use. Set plot_ellipse=FALSE to not plot ellipse ###
epidome_object_clinical_norm = normalize_epidome_object(epidome_object_clinical) ### Normalize counts to percent
PCA_patient_colored = plot_PCA_epidome(epidome_object_clinical_norm,color_variable = "patient.ID",colors = c(),plot_ellipse = FALSE)
PCA_patient_colored + ggtitle("PCA plot of nose and skin samples colored by subject")
PCA_sample_site_colored = plot_PCA_epidome(epidome_object_clinical_norm,color_variable = "sample.site",colors = c("Red","Blue"),plot_ellipse = TRUE)
PCA_sample_site_colored + ggtitle("PCA plot of nose and skin samples colored by sampling site")



### Other stuff


### Filter lowcount ASVs (removes any ASV that does not appear with more than X% abundance in any sample) ###
epidome_filtered_ASVs = filter_lowcount_ASVs_epidome(epidome_object,percent_threshold = 1)

### Normalize data to percent ###
epidome_object_normalized = normalize_epidome_object(epidome_filtered_samples)


