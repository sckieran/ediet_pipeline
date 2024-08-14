args = commandArgs(trailingOnly=TRUE)

#args list: 1) prefix, 2) gene, 3) rlib location, 4) project dirr, 5) name of NCBI tax file
library("tidyverse",lib=args[3])

setwd(args[4])

rbo <- read.delim(paste0(args[4],args[1],"_",args[2],"_raw_blast_out"), header=FALSE)
rrbo <- read.delim(paste0(args[4],args[1],"_",args[2],"_remote_raw_blast_out"), header=FALSE)

ncbi <- read.csv(paste0(args[4],args[5]))

colnames(rbo) <- c("seqnum","accession","identity","length","taxon_info","bitscore")
colnames(rrbo) <- c("seqnum","accession","identity_remote","length","taxon_remote","bitscore_remote","tax_id")
rbo$tax_id <- as.integer(str_split_fixed(rbo$taxon_info,"taxid=",2)[,2])
rrbo$tax_id <- as.integer(rrbo$tax_id)

t3 <- left_join(rbo,ncbi,by="tax_id")
t3 <- t3[,1:14]
t4 <- t3 %>% group_by(seqnum) %>% filter(bitscore==max(bitscore)) %>% summarise(identity=mean(identity),"bitscore"=mean(bitscore),length=mean(length),n_species=n_distinct(species), species=paste(unique(species), collapse=','),n_genus=n_distinct(genus),genus=paste(unique(genus), collapse=','),n_family=n_distinct(family),family=paste(unique(family), collapse=','), n_order=n_distinct(order), order=paste(unique(order), collapse=',') ,n_class=n_distinct(class),class=paste(unique(class), collapse=','),n_phylum=n_distinct(phylum),phylum=paste(unique(phylum), collapse=','))
t4$resolution <- ifelse(!rowSums(t4 == 1), names(t4)[ncol(t4)], names(t4)[max.col(t4 == 1, 'first')])
t4$resolution <- gsub("n_","",t4$resolution)

sp_df <- filter(t4, resolution=="species")
sp_df$best_hit <- sp_df$species
gen_df <- filter(t4,resolution=="genus")
gen_df$best_hit <- paste(gen_df$genus,"sp.")
fam_df <- filter(t4,resolution=="family")
fam_df$best_hit <- paste(fam_df$family,"sp.")
ord_df <- filter(t4,resolution=="order")
ord_df$best_hit <- paste(ord_df$order,"sp.")
class_df <- filter(t4,resolution=="class")
class_df$best_hit <- paste(class_df$class,"sp.")
phy_df <- filter(t4, resolution=="phylum")
phy_df$best_hit <- paste(phy_df$phylum,"sp.")
list_of_dfs <- mget(ls(pattern = "_df$"))
best_hits_all <- bind_rows(list_of_dfs)
bbh <- best_hits_all[,c(1,18,2:4,6,8,10,12,14,16,17,5)]
colnames(bbh) <- c("seqnum","best_hit","identity","bitscore","length","species","genus","family","order","class","phylum","best_resolution","num_species_in_best_hit")

rm(list = ls(pattern = "_df$"))

t3_r <- left_join(rrbo,ncbi,by="tax_id")
t3_r <- t3_r[,1:14]
t4_r <- t3_r %>% group_by(seqnum) %>% filter(bitscore_remote==max(bitscore_remote)) %>% summarise(identity=mean(identity_remote),"bitscore"=mean(bitscore_remote),length=mean(length),n_species=n_distinct(species), species=paste(unique(species), collapse=','),n_genus=n_distinct(genus),genus=paste(unique(genus), collapse=','),n_family=n_distinct(family),family=paste(unique(family), collapse=','), n_order=n_distinct(order), order=paste(unique(order), collapse=',') ,n_class=n_distinct(class),class=paste(unique(class), collapse=','),n_phylum=n_distinct(phylum),phylum=paste(unique(phylum), collapse=','))
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
colnames(bbh_r) <- c("seqnum","best_hit_remote","identity_remote","bitscore_remote","length_remote","species_remote","genus_remote","family_remote","order_remote","class_remote","phylum_remote","best_resolution_remote","num_species_remote")

bbh_full <- full_join(bbh,bbh_r,by="seqnum")

stable <- read.delim(paste0(args[4],"cat_file_list.txt"),header=FALSE)
colnames(stable) <- c("sample","sequence","reads")

asvs <- read.delim(paste0(args[4],"asvs.txt"), header=FALSE)
colnames(asvs) <- c("seqnum","sequence")

seq_table <- left_join(stable,asvs,by="sequence")
seq_table$seqnum <- gsub(">","",seq_table$seqnum)

local_taxa_table <- left_join(seq_table,bbh,by="seqnum")
local_taxa_table$best_hit[which(is.na(local_taxa_table$best_hit))] <- "No Hit"
write_delim(local_taxa_table, paste0(args[4],"/",args[1],"_",args[2],"_full_local_taxatable.txt"),delim="\t",quote="none")

remote_taxa_table <- left_join(seq_table,bbh_r,by="seqnum")
remote_taxa_table$best_hit_remote[which(is.na(remote_taxa_table$best_hit_remote))] <- "No Hit"
write_delim(local_taxa_table, paste0(args[4],"/",args[1],"_",args[2],"_full_local_taxatable.txt"),delim="\t",quote="none")

combo_taxa_table <- left_join(seq_table,bbh_full,by="seqnum")
combo_taxa_table$best_hit[which(is.na(taxa_table$best_hit))] <- "No Hit"
names(combo_taxa_table)[names(combo_taxa_table) == 'best_hit'] <- 'best_hit_local'
write_delim(local_taxa_table, paste0(args[4],"/",args[1],"_",args[2],"_full_combined_localremote_taxatable.txt"),delim="\t",quote="none")

local_taxa_sample_summary <- local_taxa_table %>% group_by(sample,species) %>% summarise("reads"=sum(reads),"mean_identity"=(mean(identity)))
write_delim(local_taxa_sample_summary, paste0(args[4],"/",args[1],"_",args[2],"_sample_by_taxon_taxatable.txt"),delim="\t",quote="none")

local_taxa_species_summary <- local_taxa_table %>% group_by(species) %>% summarise("reads"=sum(reads),mean_identity=(mean(identity)),"n_samps"=n_distinct(sample))
write_delim(t3, paste0(args[4],"/",args[1],"_",args[2],"_species_summary_taxatable.txt"),delim="\t",quote="none")

remote_taxa_sample_summary <- remote_taxa_table %>% group_by(sample,species) %>% summarise("reads"=sum(reads),"mean_identity"=(mean(identity)))
write_delim(remote_taxa_sample_summary, paste0(args[4],"/",args[1],"_",args[2],"_sample_by_taxon_taxatable.txt"),delim="\t",quote="none")

remote_taxa_species_summary <- remote_taxa_table %>% group_by(species) %>% summarise("reads"=sum(reads),mean_identity=(mean(identity)),"n_samps"=n_distinct(sample))
write_delim(remote_taxa_species_summary, paste0(args[4],"/",args[1],"_",args[2],"_species_summary_taxatable.txt"),delim="\t",quote="none")
