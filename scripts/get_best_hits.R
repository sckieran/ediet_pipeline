args = commandArgs(trailingOnly=TRUE)

#args list: 1) prefix, 2) gene, 3) rlib location, 4) project dirr, 5) name of NCBI tax file
library("tidyverse",lib=args[3])

setwd(args[1])

rbo <- read.delim(paste0(args[4],args[1],"_",args[2],"_raw_blast_out"), header=FALSE)

ncbi <- read.csv(paste0(args[4],args[5]))

colnames(rbo) <- c("seqnum","accession","identity","length","taxon_info","bitscore")
rbo$tax_id <- as.integer(str_split_fixed(rbo$taxon_info,"taxid=",2)[,2])

t3 <- left_join(rbo,ncbi,by="tax_id")
t3 <- t3[,1:14]
t4 <- t3 %>% group_by(seqnum) %>% summarise(identity=mean(identity),"bitscore"=mean(bitscore),length=mean(length),n_species=n_distinct(species), species=paste(unique(species), collapse=','),n_genus=n_distinct(genus),genus=paste(unique(genus), collapse=','),n_family=n_distinct(family),family=paste(unique(family), collapse=','), n_order=n_distinct(order), order=paste(unique(order), collapse=',') ,n_class=n_distinct(class),class=paste(unique(class), collapse=','),n_phylum=n_distinct(phylum),phylum=paste(unique(phylum), collapse=','))
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

stable <- read.delim(paste0(args[4],"cat_file_list.txt"),header=FALSE)
colnames(stable) <- c("sample","sequence","reads")

asvs <- read.delim(paste0(args[4],"asvs.txt"), header=FALSE)
colnames(asvs) <- c("seqnum","sequence")

seq_table <- left_join(stable,asvs,by="sequence")
seq_table$seqnum <- gsub(">","",seq_table$seqnum)

taxa_table <- left_join(seq_table,bbh,by="seqnum")
taxa_table$best_hit[which(is.na(taxa_table$best_hit))] <- "No Hit"

write_delim(taxa_table, paste0(args[4],args[1],"_",args[2],"_full_taxatable.txt"),delim="\t",quote="none")

taxa_sample_summary <- taxa_table %>% group_by(sample,species) %>% summarise("reads"=sum(reads),"mean_identity"=(mean(identity)))
write_delim(taxa_sample_summary, paste0(args[4],args[1],"_",args[2],"_sample_by_taxon_taxatable.txt"),delim="\t",quote="none")

taxa_species_summary <- taxa_table %>% group_by(species) %>% summarise("reads"=sum(reads),mean_identity=(mean(identity)),"n_samps"=n_distinct(sample))
write_delim(taxa_species_summary, paste0(args[4],args[1],"_",args[2],"_species_summary_taxatable.txt"),delim="\t",quote="none")

