source("https://raw.githubusercontent.com/ssi-dk/epidome/master/scripts/epidome_functions.R")
setwd("/Volumes/data/MPV/projects/git.repositories/epidome/")

### Load amplicon table for classification ###
ST_amplicon_table = read.table("DB/epidome_ST_amplicon_frequencies.txt",sep = "\t")

### Load dada2 output for the two primers ###
epi01_mock_table = read.table("example_data/190920_run2_G216_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_mock_table = read.table("example_data/190920_run2_yycH_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi01_clinical_table = read.table("example_data/190920_run1_G216_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_clinical_table = read.table("example_data/190920_run1_yycH_seqtab_from_dada2.csv.classified.csv",sep = ";",header=TRUE,row.names=1)

epi01_table = combine_ASV_tables(epi01_mock_table,epi01_clinical_table)
epi02_table = combine_ASV_tables(epi02_mock_table,epi02_clinical_table)


epi01_table = read.table("example_data/190920_run1_and_2_G216_seqtab_nochim.csv.classified.csv",sep = ";",header=TRUE,row.names=1)
epi02_table = read.table("example_data/190920_run1_and_2_yycH_seqtab_nochim.csv.classified.csv",sep = ";",header=TRUE,row.names=1)

### Load metadata table
metadata_table = read.table("example_data/article_metadata_with_qPCR.txt",header=TRUE,row.names=1)

### Setup an object for easy handling of epidome data
epidome_object = setup_epidome_object(epi01_table,epi02_table,metadata_table = metadata_table)
#epidome_object = setup_epidome_object(epi01_table,epi02_table)

### Check if number of sequences from each primer for each samples match up approximately ###
compare_primer_output(epidome_object)
compare_primer_output(epidome_object,"sample.type")
compare_primer_output(epidome_object,"patient.sample.site")

plot(colSums(epidome_object$p1_table),epidome_object$metadata$qPCR)
plot(colSums(epidome_object$p2_table),epidome_object$metadata$qPCR)



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
count_table_clinical = classify_epidome(epidome_object_clinical_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_clinical,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in nose and skin samples") + theme(axis.text.x = element_text(angle = 90))

### make barplot of mock samples only 
epidome_object_mock_ASV_combined = prune_by_variable_epidome(epidome_ASV_combined,"sample.type",c("Mock community"))
count_table_mock = classify_epidome(epidome_object_mock_ASV_combined,ST_amplicon_table)
p = make_barplot_epidome(count_table_mock,reorder=FALSE,normalize=TRUE)
p + ggtitle("Distribution of S. epidermidis sequence types in mock communities") + theme(axis.text.x = element_blank())



### make PCA plots of clinical samples. Color according to variable in metadata and (optional) indicate colors to use. Set plot_ellipse=FALSE to not plot ellipse ###
epidome_object_clinical_norm = normalize_epidome_object(epidome_object_clinical) ### Normalize counts to percent
PCA_patient_colored = plot_PCA_epidome(epidome_object_clinical_norm,color_variable = "patient.ID",colors = c(),plot_ellipse = FALSE)
PCA_patient_colored + ggtitle("PCA plot of nose and skin samples colored by subject")
PCA_sample_site_colored = plot_PCA_epidome(epidome_object_clinical_norm,color_variable = "sample.site",colors = c("Red","Blue"),plot_ellipse = TRUE)
PCA_sample_site_colored + ggtitle("PCA plot of nose and skin samples colored by sampling site")



norm_combined = rbind(epidome_object_clinical_norm$p1_table,epidome_object_clinical_norm$p2_table)
m = epidome_object_clinical_norm$metadata
idx = grep("_1",m$sample.ID)
norm_combined_uniq = norm_combined[,idx]
m = m[idx,]

dist = dist(t(norm_combined_uniq))
anosim(dist,m$sample.site,1000)


### Other stuff


### Filter lowcount ASVs (removes any ASV that does not appear with more than X% abundance in any sample) ###
epidome_filtered_ASVs = filter_lowcount_ASVs_epidome(epidome_object,percent_threshold = 1)

### Normalize data to percent ###
epidome_object_normalized = normalize_epidome_object(epidome_filtered_samples)


eo_clin_norm = normalize_epidome_object(epidome_object_clinical_ASV_combined)
count_table_clinical_norm = classify_epidome(eo_clin_norm,ST_amplicon_table)
plot(sort(bin_mat[which(bin_mat>0)]))
bin_mat = as.matrix(count_table_clinical_norm)
plot(sort(bin_mat[which(bin_mat>0)]))
sort(bin_mat[which(bin_mat>0)])
sort(as.vector(bin_mat[which(bin_mat>0)]))[1:100]
plot(sort(as.vector(bin_mat[which(bin_mat>0)]))[1:300])

p_theshold = 3
bin_mat[bin_mat<p_theshold] = 0
bin_mat[bin_mat>=p_theshold] = 1
sort(colSums(bin_mat))
mean(colSums(bin_mat))
