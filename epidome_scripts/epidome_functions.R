require(ggplot2)
require(pls)
require(reshape)
require(vegan)

setup_epidome_object <- function(primer1_table,primer2_table,metadata_table) {
  primer1_counts = primer1_table[,3:ncol(primer1_table)]
  primer2_counts = primer2_table[,3:ncol(primer2_table)]
  primer1_all_sample_names = colnames(primer1_counts)
  primer2_all_sample_names = colnames(primer2_counts)
  samples_with_both_primers = primer1_all_sample_names[which(primer1_all_sample_names %in% primer2_all_sample_names)]
  samples_missing_primer1_data = primer2_all_sample_names[which(!primer2_all_sample_names %in% primer1_all_sample_names)]
  samples_missing_primer2_data = primer1_all_sample_names[which(!primer1_all_sample_names %in% primer2_all_sample_names)]
  primer1_seqs = as.vector(primer1_table$Seq_number)
  primer2_seqs = as.vector(primer2_table$Seq_number)
  primer1_seqs[which(is.na(primer1_seqs))] = "seqUnclassified"
  primer2_seqs[which(is.na(primer2_seqs))] = "seqUnclassified"
  primer1_counts = primer1_table[,which(colnames(primer1_table) %in% samples_with_both_primers)]
  primer2_counts = primer2_table[,which(colnames(primer2_table) %in% samples_with_both_primers)]
  if (!missing(metadata_table)) {
    metadata_names = rownames(metadata_table)
    samples_missing_metadata = samples_with_both_primers[which(!samples_with_both_primers %in% metadata_names)]
    samples_missing_primer_data = metadata_names[which(!metadata_names %in% samples_with_both_primers)]
    include_names = metadata_names[which(metadata_names %in% samples_with_both_primers)]
    metadata_include = match(include_names,metadata_names)
    metadata_include = metadata_include[which(!is.na(metadata_include))]
    metadata_table = metadata_table[metadata_include,]
    #primer1_include = match(colnames(primer1_counts),include_names)
    primer1_include = match(include_names,colnames(primer1_counts))
    primer1_include = primer1_include[which(!is.na(primer1_include))]
    #primer2_include = match(colnames(primer2_counts),include_names)
    primer2_include = match(include_names,colnames(primer2_counts))
    primer2_include = primer2_include[which(!is.na(primer2_include))]
    epi1_table = primer1_counts[,primer1_include]
    epi2_table = primer2_counts[,primer2_include]
    epidome_object = list('p1_seqs'=primer1_seqs,'p1_table'=epi1_table,'p2_seqs'=primer2_seqs,'p2_table'=epi2_table,'metadata'=metadata_table,'sample_names'=include_names,'meta_variables'=colnames(metadata_table))
    print(paste0("Metadata loaded with ",length(metadata_names)," samples and ",ncol(metadata_table)," variables"))
    print(paste0(length(samples_missing_metadata)," samples are found in both primer sets but not in metadata: ",paste0(samples_missing_metadata,collapse = " ")))
    print(paste0(length(samples_missing_primer_data)," samples are found in metadata but is missing from one or both primer sets: ",paste0(samples_missing_primer_data,collapse = " ")))
    print(paste0(length(include_names)," samples are found in metadata and both tables and are included in epidome object"))
    
  } else {
    epi1_table = primer1_counts[,match(colnames(primer1_counts),samples_with_both_primers)]
    epi2_table = primer2_counts[,match(colnames(primer2_counts),samples_with_both_primers)]
    epidome_object = list('p1_seqs'=primer1_seqs,'p1_table'=epi1_table,'p2_seqs'=primer2_seqs,'p2_table'=epi2_table,'sample_names'=samples_with_both_primers)
    print(paste0("No metadata loaded"))
    print(paste0(length(samples_missing_primer2_data)," samples are found in p1 table but not in p2 table: ",paste0(samples_missing_primer2_data,collapse = " ")))
    print(paste0(length(samples_missing_primer1_data)," samples are found in p2 table but not in p1 table: ",paste0(samples_missing_primer1_data,collapse = " ")))
    print(paste0(length(samples_with_both_primers)," samples are found in both tables and are included in epidome object"))
  }
  return(epidome_object)
}



compare_primer_output <- function(epidome_object, color_variable = "") {
  p1_counts = colSums(epidome_object$p1_table)
  p2_counts = colSums(epidome_object$p2_table)
  cor = cor.test(p1_counts,p2_counts)
  return_df = data.frame(p1_counts,p2_counts,epidome_object$metadata)
  if (color_variable=="") {
    p =ggplot(data.frame(p1_counts,p2_counts),aes(x=p1_counts,y=p2_counts)) + geom_point() + ggtitle(paste0("Pearson correlation coefficient: ",sprintf(cor$estimate, fmt = '%#.2f'),", p = ",sprintf(cor$p.value, fmt = '%#.3f'))) +
      theme_bw() + xlab("g216 read count") + ylab("yycH read count")
  } else {
    p =ggplot(data.frame(p1_counts,p2_counts),aes(x=p1_counts,y=p2_counts,color = get_metadata(epidome_object,color_variable))) + geom_point() + ggtitle(paste0("Pearson correlation coefficient: ",sprintf(cor$estimate, fmt = '%#.2f'),", p = ",sprintf(cor$p.value, fmt = '%#.3f'))) + labs(color = color_variable) +
      theme_bw() + xlab("g216 read count") + ylab("yycH read count")
  }
  
  p
  return(list('plot'=p,'cor_test'=cor,'df'=return_df))
}


combine_ASVs_epidome = function(epidome_object) {
  return_epidome_object = epidome_object
  uniq_p1_seqs = unique(as.vector(epidome_object$p1_seqs))
  uniq_p2_seqs = unique(as.vector(epidome_object$p2_seqs))
  return_p1_table = matrix(nrow = 0, ncol = length(epidome_object$sample_names))
  return_p2_table = matrix(nrow = 0, ncol = length(epidome_object$sample_names))
  uniq_p1_seqs = uniq_p1_seqs[!uniq_p1_seqs=="seqUnclassified"]
  uniq_p2_seqs = uniq_p2_seqs[!uniq_p2_seqs=="seqUnclassified"]
  for (ASV in uniq_p1_seqs) {
    sub_d = epidome_object$p1_table[which(epidome_object$p1_seqs==ASV),]
    ASV_sums = colSums(sub_d)
    return_p1_table = rbind(return_p1_table,ASV_sums)
  }
  p1_unclassified_d = epidome_object$p1_table[which(epidome_object$p1_seqs=="seqUnclassified"),]
  return_p1_table = rbind(return_p1_table,p1_unclassified_d)
  rownames(return_p1_table) = 1:nrow(return_p1_table)
  return_p1_seqs = rep("seqUnclassified",nrow(return_p1_table))
  return_p1_seqs[1:length(uniq_p1_seqs)] = uniq_p1_seqs
  
  for (ASV in uniq_p2_seqs) {
    sub_d = epidome_object$p2_table[which(epidome_object$p2_seqs==ASV),]
    ASV_sums = colSums(sub_d)
    return_p2_table = rbind(return_p2_table,ASV_sums)
  }
  p2_unclassified_d = epidome_object$p2_table[which(epidome_object$p2_seqs=="seqUnclassified"),]
  return_p2_table = rbind(return_p2_table,p2_unclassified_d)
  rownames(return_p2_table) = 1:nrow(return_p2_table)
  return_p2_seqs = rep("seqUnclassified",nrow(return_p2_table))
  return_p2_seqs[1:length(uniq_p2_seqs)] = uniq_p2_seqs
  return_epidome_object$p1_seqs = return_p1_seqs
  return_epidome_object$p2_seqs = return_p2_seqs
  return_epidome_object$p1_table = return_p1_table
  return_epidome_object$p2_table = return_p2_table
  return(return_epidome_object)
}


filter_lowcount_samples_epidome = function(epidome_object,p1_threshold,p2_threshold) {
  original_sample_names = epidome_object$sample_names
  include_index = which(colSums(epidome_object$p1_table) >= p1_threshold & colSums(epidome_object$p2_table) >= p2_threshold)
  epidome_object$p1_table = epidome_object$p1_table[,include_index]
  epidome_object$p2_table = epidome_object$p2_table[,include_index]
  epidome_object$sample_names = epidome_object$sample_names[include_index]
  epidome_object$metadata = epidome_object$metadata[include_index,]
  filtered_samples = original_sample_names[-include_index]
  print(paste0(length(filtered_samples)," low count samples removed from data: ",paste0(filtered_samples,collapse = " ")))
  return(epidome_object)
}


filter_lowcount_ASVs_epidome = function(epidome_object,percent_threshold) {
  return_epidome_object = epidome_object
  epidome_object_norm = normalize_epidome_object(epidome_object)
  p1_max = apply(epidome_object_norm$p1_table, 1, max)
  p2_max = apply(epidome_object_norm$p2_table, 1, max)
  p1_include_idx = which(p1_max>=percent_threshold)
  p2_include_idx = which(p2_max>=percent_threshold)
  print(paste0((length(epidome_object$p1_seqs)-length(p1_include_idx)), " ASVs removed from epi01 data"))
  print(paste0(length(p1_include_idx), " ASVs remaining in epi01 data"))
  print(paste0((length(epidome_object$p2_seqs)-length(p2_include_idx)), " ASVs removed from epi02 data"))
  print(paste0(length(p2_include_idx), " ASVs remaining in epi02 data"))
  return_epidome_object$p1_seqs = epidome_object$p1_seqs[p1_include_idx]
  return_epidome_object$p2_seqs = epidome_object$p2_seqs[p2_include_idx]
  return_epidome_object$p1_table = epidome_object$p1_table[p1_include_idx,]
  return_epidome_object$p2_table = epidome_object$p2_table[p2_include_idx,]
  return(return_epidome_object)
}

normalize_epidome_object = function(epidome_object) {
  epidome_object$p1_table = apply(epidome_object$p1_table,2,function(x) x/sum(x)*100)
  epidome_object$p2_table = apply(epidome_object$p2_table,2,function(x) x/sum(x)*100)
  return(epidome_object)
}

classify_epidome = function(epidome_object,ST_amplicon_table,strict_classifier=FALSE) {
  epidome_object_norm = normalize_epidome_object(epidome_object)
  p1_seqs = unlist(lapply(as.vector(epidome_object$p1_seqs),function(x) substr(x,4,nchar(x))))
  p2_seqs = unlist(lapply(as.vector(epidome_object$p2_seqs),function(x) substr(x,4,nchar(x))))
  n_samples = length(epidome_object$sample_names)
  n_p1_seqs = length(p1_seqs)
  count_table = matrix(nrow = 0, ncol = n_samples,dimnames = list('ST'=c(),'Samples'=epidome_object$sample_names))
  match_type_table = matrix(nrow = n_p1_seqs, ncol = n_samples)
  count_table_names = c()
  unclassified_count_vec = rep(0,n_samples)
  g1_unclassified_count_vec = rep(0,n_samples)
  for (i in 1:n_p1_seqs) {
    p1_seq = p1_seqs[i]
    p1_seq_split = strsplit(p1_seq,',')[[1]]
    if (length(p1_seq_split)>1) {
      possible_STs = as.vector(ST_amplicon_table$ST)[which(ST_amplicon_table$epi01_ASV %in% p1_seq_split)]
      if (length(possible_STs)==1) {
        p1_seq = p1_seq_split[1]
      } else {
        p1_seq = "Unclassified"
      }
    }
    if (p1_seq != "Unclassified") {
      p1_seq_ST_table = ST_amplicon_table[which(ST_amplicon_table$epi01_ASV==p1_seq),]
      unique_p1_ASVs = unique(as.vector(p1_seq_ST_table$ST))
      p2_ASVs_matching_p1 = p1_seq_ST_table$epi02_ASV
      count_vec = rep(0,n_samples)
      for (j in 1:n_samples) {
        p1_percent = epidome_object_norm$p1_table[i,j]
        p1_count = epidome_object$p1_table[i,j]
        if (p1_percent > 0.01) {
          p2_seqs_present_within_difference_threshold_idx = which(epidome_object_norm$p2_table[,j]>(p1_percent-10) & epidome_object_norm$p2_table[,j]>1)
          p2_seqs_present_ASVs = p2_seqs[p2_seqs_present_within_difference_threshold_idx]
          p2_seqs_present_ASVs_matching_p1 = p2_seqs_present_ASVs[which(p2_seqs_present_ASVs %in% p2_ASVs_matching_p1)]
          p2_percent = epidome_object_norm$p2_table[which(p2_seqs %in% p2_seqs_present_ASVs_matching_p1),j]
          p2_count = epidome_object$p2_table[which(p2_seqs %in% p2_seqs_present_ASVs_matching_p1),j]
          if (length(p2_seqs_present_ASVs_matching_p1) == 0) {
            match_type_table[i,j] = "epi01 match without epi02 match"
            ST_order = order(p1_seq_ST_table$freq,decreasing = T)
            classification_group = as.vector(p1_seq_ST_table$ST)[ST_order[1]]
            if (classification_group %in% count_table_names) {
              classification_idx = which(count_table_names == classification_group)
              count_table[classification_idx,j] = count_table[classification_idx,j]+p1_count
            } else {
              count_vec = rep(0,n_samples)
              count_vec[j] = p1_count
              count_table = rbind(count_table,count_vec)
              count_table_names = c(count_table_names,classification_group)
            }
          }
          else if (length(p2_seqs_present_ASVs_matching_p1) == 1) {
            p1_p2_seq_ST_table = p1_seq_ST_table[which(p1_seq_ST_table$epi02_ASV==p2_seqs_present_ASVs_matching_p1[1]),]
            ST_order = order(p1_p2_seq_ST_table$freq,decreasing = T)
            classification_group = as.vector(p1_p2_seq_ST_table$ST)[ST_order[1]]
            if (classification_group %in% count_table_names) {
              classification_idx = which(count_table_names == classification_group)
              count_table[classification_idx,j] = count_table[classification_idx,j]+p1_count
            } else {
              count_vec = rep(0,n_samples)
              count_vec[j] = p1_count
              count_table = rbind(count_table,count_vec)
              count_table_names = c(count_table_names,classification_group)
            }
            match_type_table[i,j] = "Unique epi01 epi02 combination"
          } else {
            if (strict_classifier) {
              unclassified_count_vec[j] = unclassified_count_vec[j]+p1_count
            } else {
              p1_p2_seq_ST_table = p1_seq_ST_table[which(p1_seq_ST_table$epi02_ASV %in% p2_seqs_present_ASVs_matching_p1),]
              ST_order = order(p1_p2_seq_ST_table$freq,decreasing = T)
              classification_group = as.vector(p1_p2_seq_ST_table$ST)[ST_order[1]]
              if (classification_group %in% count_table_names) {
                classification_idx = which(count_table_names == classification_group)
                count_table[classification_idx,j] = count_table[classification_idx,j]+p1_count
              } else {
                count_vec = rep(0,n_samples)
                count_vec[j] = p1_count
                count_table = rbind(count_table,count_vec)
                count_table_names = c(count_table_names,classification_group)
              }
            }
            
            #unclassified_count_vec[j] = unclassified_count_vec[j]+p1_count
            match_type_table[i,j] = c("Non unique epi01 epi02 combination")
          }
          
        } else {
          match_type_table[i,j] = "Low counts"
          unclassified_count_vec[j] = unclassified_count_vec[j]+p1_count
        }
        
      }
    } else {
      count_vec = as.numeric(as.vector(epidome_object$p1_table[i,]))
      unclassified_count_vec = unclassified_count_vec+count_vec
      g1_unclassified_count_vec = g1_unclassified_count_vec+count_vec
      match_type_vec = rep("Unclassified epi01 match",n_samples)
      match_type_table[i,] = match_type_vec
    }
    
  }
  count_table = rbind(count_table,unclassified_count_vec)
  rownames(count_table) = c(count_table_names,"Unclassified")
  count_df = as.data.frame(count_table)
  return(count_df)
}


make_barplot_epidome_1_5 = function(count_table, reorder = FALSE, normalize = TRUE) {
	  count_df_ordered = count_table[order(rowSums(count_table),decreasing = T),]
  if (normalize) {
	      dd<-apply(count_df_ordered, 2, function(x) x/sum(x)*100)
      count_df_ordered<-as.data.frame(dd)
          ylabel = "Relative abundance (% of sequences)"
        } else {
		    ylabel = "Absolute abundance (number of sequences)"
	  }
    count_novel = count_df_ordered[which(rownames(count_df_ordered)=="-"),]
    count_unclassified = count_df_ordered[which(rownames(count_df_ordered)=="Unclassified"),]
      count_df_top12 = count_df_ordered[-which(rownames(count_df_ordered) %in% c("-","Unclassified")),]
      count_df_top12 = rbind(count_df_top12[1:5,],colSums(count_df_top12[-c(1:5),])+count_novel,count_unclassified)
        rownames(count_df_top12)[c((nrow(count_df_top12)-1),nrow(count_df_top12))] = c("Other","Unclassified")
        count_df_top12$ST = rownames(count_df_top12)
	  melt_df = melt(count_df_top12)
	  colnames(melt_df) = c("ST","Sample","Count")
	    non_ST_levels = c("Other","Unclassified")
	    ST_levels = as.vector(count_df_top12$ST)[which(!count_df_top12$ST %in% non_ST_levels)]
	      ST_levels_ordered = c(sort(as.numeric(ST_levels)),non_ST_levels)
	      if (reorder) {
		          BC = vegdist(t(count_table))
	          fit = hclust(BC, method = "ward.D")
		      melt_df$Sample = factor(as.vector(melt_df$Sample), levels = fit$labels[fit$order])
		    }
	        melt_df$ST = factor(as.vector(melt_df$ST),levels = ST_levels_ordered)
	        p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = RColorBrewer::brewer.pal(7,"Paired")) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic()
		  return(p)
}

make_barplot_epidome = function(count_table, reorder = FALSE, normalize = TRUE) {
  count_df_ordered = count_table[order(rowSums(count_table),decreasing = T),]
  if (normalize) {
    dd<-apply(count_df_ordered, 2, function(x) x/sum(x)*100)
    count_df_ordered<-as.data.frame(dd)
    ylabel = "Relative abundance (% of sequences)"
  } else {
    ylabel = "Absolute abundance (number of sequences)"
  }
  count_novel = count_df_ordered[which(rownames(count_df_ordered)=="-"),]
  count_unclassified = count_df_ordered[which(rownames(count_df_ordered)=="Unclassified"),]
  #NEED MODIFICATION 10-->27
  count_df_top12 = count_df_ordered[-which(rownames(count_df_ordered) %in% c("-","Unclassified")),]
  count_df_top12 = rbind(count_df_top12[1:28,],colSums(count_df_top12[-c(1:28),])+count_novel,count_unclassified)
  rownames(count_df_top12)[c((nrow(count_df_top12)-1),nrow(count_df_top12))] = c("Novel","Unclassified")
  count_df_top12$ST = rownames(count_df_top12)
  melt_df = melt(count_df_top12)
  colnames(melt_df) = c("ST","Sample","Count")
  non_ST_levels = c("Novel","Unclassified")
  ST_levels = as.vector(count_df_top12$ST)[which(!count_df_top12$ST %in% non_ST_levels)]
  #MODIFIED
  #ST_levels_ordered = c(sort(as.numeric(ST_levels)),non_ST_levels)
  ST_levels_ordered = c(sort(ST_levels),non_ST_levels)
  if (reorder) {
    BC = vegdist(t(count_table))
    fit = hclust(BC, method = "ward.D")
    melt_df$Sample = factor(as.vector(melt_df$Sample), levels = fit$labels[fit$order])
  }
  melt_df$ST = factor(as.vector(melt_df$ST),levels = ST_levels_ordered)
  #MODIFIED
  #https://www.r-bloggers.com/2013/09/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
  #values = getPalette(colourCount)
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))  #(12, "Paired"))  #(8,"Accent")) #(9, "Set1"))
  p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = getPalette(30)) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic() + theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))
  #p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(29)) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic()
  #p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = RColorBrewer::brewer.pal(12,"Paired")) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic()
  return(p)
}

#rename Other to Novel, Unclassified are not accounted, since each sequence couldn't find in the database. The Novel are the combination of the two amplicon couldn't find in the databases.
make_barplot_epidome_1_27 = function(count_table, reorder = FALSE, normalize = TRUE) {
  count_df_ordered = count_table[order(rowSums(count_table),decreasing = T),]
  if (normalize) {
    dd<-apply(count_df_ordered, 2, function(x) x/sum(x)*100)
    count_df_ordered<-as.data.frame(dd)
    ylabel = "Relative abundance (% of sequences)"
  } else {
    ylabel = "Absolute abundance (number of sequences)"
  }
  count_novel = count_df_ordered[which(rownames(count_df_ordered)=="-"),]
  #count_unclassified = count_df_ordered[which(rownames(count_df_ordered)=="Unclassified"),]
  #NEED MODIFICATION 10-->26
  count_df_top12 = count_df_ordered[-which(rownames(count_df_ordered) %in% c("-")),]
  count_df_top12 = rbind(count_df_top12[1:27,],colSums(count_df_top12[-c(1:27),])+count_novel)
  rownames(count_df_top12)[c( nrow(count_df_top12) )] = c("Novel")
  count_df_top12$ST = rownames(count_df_top12)
  melt_df = melt(count_df_top12)
  colnames(melt_df) = c("ST","Sample","Count")
  non_ST_levels = c("Novel")
  ST_levels = as.vector(count_df_top12$ST)[which(!count_df_top12$ST %in% non_ST_levels)]
  #MODIFIED
  #ST_levels_ordered = c(sort(as.numeric(ST_levels)),non_ST_levels)
  ST_levels_ordered = c(sort(ST_levels),non_ST_levels)
  if (reorder) {
    BC = vegdist(t(count_table))
    fit = hclust(BC, method = "ward.D")
    melt_df$Sample = factor(as.vector(melt_df$Sample), levels = fit$labels[fit$order])
  }
  melt_df$ST = factor(as.vector(melt_df$ST),levels = ST_levels_ordered)
  #MODIFIED
  #https://www.r-bloggers.com/2013/09/how-to-expand-color-palette-with-ggplot-and-rcolorbrewer/
  #values = getPalette(colourCount)
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))  #(12, "Paired"))  #(8,"Accent")) #(9, "Set1"))
  p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = getPalette(28)) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic() + theme(axis.text.x = element_text(angle = 30, vjust = 0.5, hjust=1))
  #p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Accent"))(28)) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic()
  #p = ggplot() + geom_bar(aes(y = Count, x = Sample, fill = ST), data = melt_df, stat="identity") + scale_fill_manual(values = RColorBrewer::brewer.pal(12,"Paired")) + theme(axis.text.x = element_text(angle = 90,hjust = 0.95)) + ylab(ylabel) + theme_classic()
  return(p)
}

setup_colors = function(factor_levels,colors) {
  group_count = length(factor_levels)
  if (!group_count == length(colors)) {
    if (group_count<=9) {
      colors = RColorBrewer::brewer.pal(length(factor_levels),"Set1")
    } else if (group_count<=12) {
      colors = RColorBrewer::brewer.pal(12,name="Paired")[1:group_count]
    } else {
      colors  = grDevices::rainbow(group_count)
    }
  }
  return(colors)
}

summarize_reads_epidome <- function(epidome_object,color_variable,p1_threshold=0,p2_threshold=0,colors) {
  m = epidome_object$metadata
  var_factor = m[,which(epidome_object$meta_variables==color_variable)]
  p1_read_sums = colSums(epidome_object$p1_table)
  p2_read_sums = colSums(epidome_object$p2_table)
  color_vector = setup_colors(levels(var_factor),colors)
  if (!class(var_factor) == "factor") {
    var_levels = as.vector(unique(var_factor))
  } else {
    var_levels = levels(var_factor)
  }
#  plot_df = data.frame("Group"=var_factor,"Read.count"=read_sums)
#  p = ggplot(plot_df,aes(x=Group,y=Read.count,fill=Group)) + geom_boxplot() + geom_jitter(width = 0.2) + ggtitle("Read count distribution") + theme_bw() + xlab(element_blank()) + ylab("Read count") + 
#    theme(legend.position = "none", axis.text.x = element_text(size=14), axis.title.y = element_text(size=14), plot.title = element_text(size=16)) + scale_fill_manual(values = col_vec)
  p1_mat = matrix(ncol=(length(var_levels)+1),nrow=7,dimnames = list(c("Number of samples","Samples included","Samples excluded","Median","SD","Minimum","Maximum"),c("Total",var_levels)))
  p2_mat = matrix(ncol=(length(var_levels)+1),nrow=7,dimnames = list(c("Number of samples","Samples included","Samples excluded","Median","SD","Minimum","Maximum"),c("Total",var_levels)))
#  rownames(return_mat) = c("Number of samples","Samples included","Samples excluded","Median","SD","Minimum","Maximum")
#  colnames(return_mat) = c("Total",var_levels)
  vec = c(length(p1_read_sums),length(which(p1_read_sums>=p1_threshold & p2_read_sums>=p2_threshold)),length(which(p1_read_sums<p1_threshold | p2_read_sums<p2_threshold)),
          median(p1_read_sums),sd(p1_read_sums),min(p1_read_sums),max(p1_read_sums))
  p1_mat[,1] = vec
  vec = c(length(p2_read_sums),length(which(p1_read_sums>=p1_threshold & p2_read_sums>=p2_threshold)),length(which(p1_read_sums<p1_threshold | p2_read_sums<p2_threshold)),
          median(p2_read_sums),sd(p2_read_sums),min(p2_read_sums),max(p2_read_sums))
  p2_mat[,1] = vec
  p1_compare_mat = matrix(nrow=length(var_levels),ncol=length(var_levels),dimnames = list(var_levels,var_levels))
  p2_compare_mat = matrix(nrow=length(var_levels),ncol=length(var_levels),dimnames = list(var_levels,var_levels))
  for (i in 1:length(var_levels)) {
    var_value = var_levels[i]
    p1_read_sums_1 = p1_read_sums[which(var_factor==var_value)]
    p2_read_sums_1 = p2_read_sums[which(var_factor==var_value)]
    vec = c(length(p1_read_sums_1),length(which(p1_read_sums_1>=p1_threshold & p2_read_sums_2>=p2_threshold)),length(which(p1_read_sums_1<p1_threshold | p2_read_sums_2<p2_threshold)),
            median(p1_read_sums_1),sd(p1_read_sums_1),min(p1_read_sums_1),max(p1_read_sums_1))
    p1_mat[,(i+1)] = vec
    vec = c(length(p2_read_sums_1),length(which(p1_read_sums_1>=p1_threshold & p2_read_sums_1>=p2_threshold)),length(which(p1_read_sums_1<p1_threshold | p2_read_sums_1<p2_threshold)),
            median(p2_read_sums_1),sd(p2_read_sums_1),min(p2_read_sums_1),max(p2_read_sums_1))
    p2_mat[,(i+1)] = vec
    for (j in i:length(var_levels)) {
      print(paste0(i,' ',j))
      p1_read_sums_2 = p1_read_sums[which(var_factor==var_levels[j])]
      print(p1_read_sums)
      print(p1_read_sums_2)
      wtest = wilcox.test(p1_read_sums,p1_read_sums_2)
      p1_compare_mat[i,j] = wtest$p.value
      p1_compare_mat[j,i] = wtest$p.value
      p2_read_sums_2 = p2_read_sums[which(var_factor==var_levels[j])]
      wtest = wilcox.test(p2_read_sums,p2_read_sums_2)
      p2_compare_mat[i,j] = wtest$p.value
      p2_compare_mat[j,i] = wtest$p.value
    }
  }
  return_df = as.data.frame(return_mat)
  return(list("Summary"=return_df,"p.values"=compare_mat))
}


plot_PCA_epidome = function(epidome_object,color_variable,colors,plot_ellipse = TRUE) {
  m = epidome_object$metadata
  color_variable_factor = m[,which(epidome_object$meta_variables==color_variable)]
  data_combined = rbind(epidome_object$p1_table,epidome_object$p2_table)
  pca = prcomp(t(data_combined))
  plot_df = data.frame(pca$x)
  color_vector = setup_colors(levels(color_variable_factor),colors)
  labels = c(paste0("PC1 [",sprintf("%.1f",explvar(pca)[1]),"%]"),paste0("PC2, [",sprintf("%.1f",explvar(pca)[2]),"%]"))
  #MODIFIED: size=1.5 --> size=3
  if (plot_ellipse) {
    p = ggplot(as.data.frame(pca$x),aes(x=PC1,y=PC2,color = color_variable_factor)) + labs(color = color_variable) + geom_point(size=3, alpha=1) + stat_ellipse(level=0.75) + scale_colour_manual(values = color_vector) + xlab(labels[1]) + ylab(labels[2]) + theme_bw()
  } else {
    p = ggplot(as.data.frame(pca$x),aes(x=PC1,y=PC2,color = color_variable_factor)) + labs(color = color_variable) + geom_point(size=3, alpha=1) + scale_colour_manual(values = color_vector) + xlab(labels[1]) + ylab(labels[2]) + theme_bw()
  }
  return(p)
}

prune_samples_epidome = function(epidome_object,sample_names, keep_factor_levels = FALSE) {
  return_epidome_object = epidome_object
  include_idx = which(epidome_object$sample_names %in% sample_names)
  return_epidome_object$p1_table = epidome_object$p1_table[,include_idx]
  return_epidome_object$p2_table = epidome_object$p2_table[,include_idx]
  new_m = epidome_object$metadata[include_idx,]
  if (!keep_factor_levels) {
    for (i in 1:ncol(new_m)) {
      col_factor = new_m[,i]
      if (class(col_factor) == "factor") {
        lvls = levels(col_factor)
        new_lvls = lvls[which(lvls %in% as.vector(col_factor))]
        new_factor = factor(as.vector(col_factor), levels = new_lvls)
        new_m[,i] = new_factor
      }
    }
  }
  return_epidome_object$metadata = new_m
  return_epidome_object$sample_names = epidome_object$sample_names[include_idx]
  return(return_epidome_object)
}

prune_by_variable_epidome = function(epidome_object,variable_name,variable_values, keep_factor_levels = FALSE) {
  return_epidome_object = epidome_object
  m = epidome_object$metadata
  variable_factor = m[,which(epidome_object$meta_variables==variable_name)]
  include_idx = which(variable_factor %in% variable_values)
  return_epidome_object$p1_table = epidome_object$p1_table[,include_idx]
  return_epidome_object$p2_table = epidome_object$p2_table[,include_idx]
  new_m = epidome_object$metadata[include_idx,]
  if (!keep_factor_levels) {
    for (i in 1:ncol(new_m)) {
      col_factor = new_m[,i]
      if (class(col_factor) == "factor") {
        lvls = levels(col_factor)
        new_lvls = lvls[which(lvls %in% as.vector(col_factor))]
        new_factor = factor(as.vector(col_factor), levels = new_lvls)
        new_m[,i] = new_factor
      }
    }
  }
  return_epidome_object$metadata = new_m
  return_epidome_object$sample_names = epidome_object$sample_names[include_idx]
  return(return_epidome_object)
}


combine_ASV_tables = function(ASV_table_1, ASV_table_2) {
  IDs1 = colnames(ASV_table_1)[3:ncol(ASV_table_1)]
  IDs2 = colnames(ASV_table_2)[3:ncol(ASV_table_2)]
  table_names = c("ASV","Seq_number",IDs1,IDs2)
  ASV_table_1$ASV = as.vector(ASV_table_1$ASV)
  ASV_table_2$ASV = as.vector(ASV_table_2$ASV)
  ASVs_1 = ASV_table_1$ASV
  ASVs_2 = ASV_table_2$ASV
  ASV1_match_index = which(ASVs_1 %in% ASVs_2)
  ASVs_in_both = ASVs_1[ASV1_match_index]
  ASV_table_1_match = ASV_table_1[ASV1_match_index,]
  ASV_table_2_match = ASV_table_2[match(ASVs_in_both,ASVs_2),]
  ASV_table_1_nomatch = ASV_table_1[-ASV1_match_index,]
  ASV_table_2_nomatch = ASV_table_2[which(!ASVs_2 %in% ASVs_in_both),]
  ASV_table_match_combined = cbind(ASV_table_1_match,ASV_table_2_match[,3:ncol(ASV_table_2_match)])
  ASV_table_1_nomatch_2 = cbind(ASV_table_1_nomatch,matrix(0L,nrow=nrow(ASV_table_1_nomatch),ncol=(ncol(ASV_table_2_nomatch)-2)))
  ASV_table_2_nomatch_2 = cbind(ASV_table_2_nomatch[,1:2],matrix(0L,nrow=nrow(ASV_table_2_nomatch),ncol=(ncol(ASV_table_1_nomatch)-2)),ASV_table_2_nomatch[,3:ncol(ASV_table_2_nomatch)])
  colnames(ASV_table_match_combined) = table_names
  colnames(ASV_table_1_nomatch_2) = table_names
  colnames(ASV_table_2_nomatch_2) = table_names
  ASV_table_combined = rbind(ASV_table_match_combined,ASV_table_1_nomatch_2,ASV_table_2_nomatch_2)
  return(ASV_table_combined)
}

dist_comparison = function(dist_object,group_factor) {
  group_levels = levels(group_factor)
  group_vector = as.vector(group_factor)
  g1_vec = c()
  g2_vec = c()
  p_mat = matrix(nrow=length(group_levels),ncol=length(group_levels),dimnames = list(group_levels,group_levels))
  median_IQR_mat = p_mat
  for (i in 1:length(group_vector)) {
    g1_vec = c(g1_vec,rep(group_vector[i],(length(group_vector)-i)))
    g2_vec = c(g2_vec,group_vector[(i+1):length(group_vector)])
  }
  g2_vec = g2_vec[1:length(g1_vec)]
  in_out_group_vec = rep("Between group",length(g1_vec))
  in_out_group_vec[which(g1_vec==g2_vec)] = "Within group"
  #group_factor = factor(paste0(g1_vec,'__',g2_vec),levels = paste0()
  return_df = data.frame('Group.1'=g1_vec,'Group.2'=g2_vec,'Groups'=paste0(g1_vec,'__',g2_vec),'within_or_between_group'=in_out_group_vec,'dist'=as.vector(dist_object))
  for (i in 1:length(group_levels)) {
    dist_g1  = return_df$dist[which(return_df$Group.1==group_levels[i])]
    for (j in 1:length(group_levels)) {
      if (!i==j) {
        dist_g2  = return_df$dist[which(return_df$Group.2==group_levels[i])]
        p_value = NA
        p_value = try(sprintf(wilcox.test(dist_g1,dist_g2)$p.value, fmt = '%#.5f'))
        p_mat[i,j] = p_value
        p_mat[j,i] = p_value
      }
      dist_g1_g2 = return_df$dist[which(return_df$Group.1==group_levels[i] & return_df$Group.2==group_levels[j])]
      iqr = quantile(dist_g1_g2)
      median_print = paste0('Median ',sprintf(iqr[3], fmt = '%#.2f'),', (IQR: ',sprintf(iqr[2], fmt = '%#.2f'),' - ',sprintf(iqr[4], fmt = '%#.2f'),')')
      median_IQR_mat[i,j] = median_print
      median_IQR_mat[j,i] = median_print
    }
  }
  p <- plot_ly(type="box",data=return_df,x=~Groups,y=~dist, boxpoints = "all", pointpos = -1.5)
  return_list = list('plot'=p,'dist.data'=return_df,'p.values'=as.data.frame(p_mat),'IQR.values'=as.data.frame(median_IQR_mat))
  return(return_list)
}


sample_names_epidome = function(epidome_object) {
  return(epidome_object$sample_names)
}
p1_table = function(epidome_object) {
  return(epidome_object$p1_table)
}
p2_table = function(epidome_object) {
  return(epidome_object$p2_table)
}
p1_seqs = function(epidome_object) {
  return(epidome_object$p1_seqs)
}
p2_seqs = function(epidome_object) {
  return(epidome_object$p2_seqs)
}
get_metadata = function(epidome_object,variable="") {
  if (variable=="") {
    return_object = epidome_object$metadata
  } else {
    return_object = epidome_object$metadata[,which(epidome_object$meta_variables==variable)]
  }
  return(return_object)
}


