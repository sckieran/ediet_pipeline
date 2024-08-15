args = commandArgs(trailingOnly=TRUE)

#args list: 1) prefix, 2) gene, 3) rlib location, 4) project dirr, 5) name of NCBI tax file
library("tidyverse",lib=args[3])

setwd(args[4])


rrbo <- read.delim(paste0(args[4],"/",args[1],"_",args[2],"_remote_raw_blast_out"), header=FALSE)

ncbi <- read.csv(paste0(args[4],"/'<args[5]))


colnames(rrbo) <- c("seqnum","accession","identity","length","taxon","bitscore_remote","tax_id")
rrbo$tax_id <- as.integer(rrbo$tax_id)


t3_r <- left_join(rrbo,ncbi,by="tax_id")
t3_r <- t3_r[,1:14]
t4_r <- t3_r %>% group_by(seqnum) %>% summarise(identity=mean(identity_remote),"bitscore"=mean(bitscore_remote),length=mean(length),n_species=n_distinct(species), species=paste(unique(species), collapse=','),n_genus=n_distinct(genus),genus=paste(unique(genus), collapse=','),n_family=n_distinct(family),family=paste(unique(family), collapse=','), n_order=n_distinct(order), order=paste(unique(order), collapse=',') ,n_class=n_distinct(class),class=paste(unique(class), collapse=','),n_phylum=n_distinct(phylum),phylum=paste(unique(phylum), collapse=','))
t4_r$resolution <- ifelse(!rowSums(t4_r == 1), names(t4_r)[ncol(t4_r)], names(t4_r)[max.col(t4_r == 1, 'first')])
t4_r$resolution <- gsub("n_","",t4_r$resolution)

sp_df_r <- filter(t4_r, resolution=="species")
sp_df_r$best_hit <- sp_df_r$species
gen_df_r <- filter(t4_r,resolution=="genus")
gen_df_r$best_hit <- paste(gen_df_r$genus,"sp.")
fam_df_r <- filter(t4_r,resolution=="family")
fam_df_r$best_hit <- paste(fam_df_r$family,"sp.")
ord_df_r <- filter(t4_r,resolution=="order")
ord_df_r$best_hit <- paste(ord_df_r$order,"sp.")
class_df_r <- filter(t4_r,resolution=="class")
class_df_r$best_hit <- paste(class_df_r$class,"sp.")
phy_df_r <- filter(t4_r, resolution=="phylum")
phy_df_r$best_hit <- paste(phy_df_r$phylum,"sp.")
list_of_dfs_r <- mget(ls(pattern = "_df_r$"))
best_hits_all_r <- bind_rows(list_of_dfs_r)
bbh_r <- best_hits_all_r[,c(1,18,2:4,6,8,10,12,14,16,17,5)]
colnames(bbh_r) <- c("seqnum","best_hit","identity","bitscore","length","species","genus","family","order","class","phylum","best_resolution","num_species")

stable <- read.delim(paste0(args[4],"cat_file_list.txt"),header=FALSE)
colnames(stable) <- c("sample","sequence","reads")

asvs <- read.delim(paste0(args[4],"asvs.txt"), header=FALSE)
colnames(asvs) <- c("seqnum","sequence")

seq_table <- left_join(stable,asvs,by="sequence")
seq_table$seqnum <- gsub(">","",seq_table$seqnum)

remote_taxa_table <- left_join(seq_table,bbh_r,by="seqnum")
remote_taxa_table$best_hit_remote[which(is.na(remote_taxa_table$best_hit_remote))] <- "No Hit"
write_delim(local_taxa_table, paste0(args[4],"/",args[1],"_",args[2],"_full_local_taxatable.txt"),delim="\t",quote="none")

bbh_tax_r <- unique(remote_taxa_table[,c(1,5,9:15)])
bbh_tax_r$species[bbh_tax_r$best_resolution!="species"] <- "Not Resolved"
bbh_tax_r$genus[bbh_tax_r$best_resolution=="family"] <- "Not Resolved"
bbh_tax_r$genus[bbh_tax_r$best_resolution=="order"] <- "Not Resolved"
bbh_tax_r$genus[bbh_tax_r$best_resolution=="class"] <- "Not Resolved"
bbh_tax_r$genus[bbh_tax_r$best_resolution=="phylum"] <- "Not Resolved"
bbh_tax_r$family[bbh_tax_r$best_resolution=="order"] <- "Not Resolved"
bbh_tax_r$family[bbh_tax_r$best_resolution=="class"] <- "Not Resolved"
bbh_tax_r$family[bbh_tax_r$best_resolution=="phylum"] <- "Not Resolved"
bbh_tax_r$order[bbh_tax_r$best_resolution=="class"] <- "Not Resolved"
bbh_tax_r$order[bbh_tax_r$best_resolution=="phylum"] <- "Not Resolved"
bbh_tax_r$class[bbh_tax_r$best_resolution=="phylum"] <- "Not Resolved"
bbh_tax_r <- unique(bbh_tax_r)

remote_taxa_sample_summary <- remote_taxa_table %>% group_by(sample,species) %>% summarise("reads"=sum(reads),"mean_identity"=(mean(identity)))
remote_taxa_sample_summary2 <- left_join(remote_taxa_sample_summary, bbh_tax_r,by=c("sample"="sample","best_hit"="best_hit"))
write_delim(remote_taxa_sample_summary2, paste0(args[4],"/",args[1],"_",args[2],"_sample_by_taxon_taxatable.txt"),delim="\t",quote="none")

remote_taxa_species_summary <- remote_taxa_table %>% group_by(species) %>% summarise("reads"=sum(reads),mean_identity=(mean(identity)),"n_samps"=n_distinct(sample),best_resolution=unique(best_resolution))
write_delim(remote_taxa_species_summary, paste0(args[4],"/",args[1],"_",args[2],"_species_summary_taxatable.txt"),delim="\t",quote="none")
