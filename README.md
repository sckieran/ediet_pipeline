# The eDiet Metabarcoding Pipeline
Shannon Kieran Blair 2024

## Overview of Pipeline

A pipeline and tutorial to analyze metabarcoding data using pears and fastx-collapser with a limited local reference database. It is intended to be an analysis guide for students working on their own projects. It should run on any SLURM-based linux cluster that supports conda environments. A non-SLURM, bash-only (slower) version is coming.

This pipeline assumes reasonably good resolution of locus data and is designed as a "first pass" at your metabarcoding data, not a finished product. The way this pipeline assesses BLAST hits is both **naive and conservative**. If more than one equally-likely BLAST hit is available (judged by bitscore), it will attempt to assign taxonomy to each hit using the NCBI taxonomy tree file, then it will up the taxonomy tree until it finds a consensus rank. If no consensus is available at the phylum level, it calls a "no-hit". Therefore **this pipeline only works with reference libraries made using this pipeline, or with remote BLAST.** However, you can use your own FASTAs as input to step 2 of this pipeline if you want to skip querying NCBI, as long as each sequence in that FASTA begins with a header line in the format >accession genus species taxid=[Valid NCBI taxid]. You can input a list of taxa into https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi to get a list of taxids.


## Glossary of terms
- _**taxid**_ the NCBI-assigned taxonomic ID number that corresponds to a particular taxon. Every taxa at every rank that has ever had sequence information submitted to NCBI has a taxid.
- _**ASV**_ Amplicon Sequence Variant. A **unique sequence** in a sample. IE, the list of all the merged reads in a sample with duplicates removed. In order to be included in the data, an ASV must be present in >1 read.
- _**RRA**_ Relative Read Abundance. The proportion of reads in a sample that any particular ASV (or taxon) must reach to be included in the downstream analysis. Usual values run between 0.001 and 0.01, ie, between 0.1% and 1% of reads. So if a sample contains 100,000 reads total, and the RRA cutoff is 0.001, any ASV (or taxon) with less than 100 reads is excluded from downstream analysis. Taxon-level filtering, done in step 8, concatenates all the reads from any ASV that BLAST to the same taxon with a percent identity above your user-defined cutoff. ASV-level filtering, done in step 5, removes any ASV that does not reach the RRA threshold BEFORE performing BLAST. The benefit of the ASV-level filter in step 5 is that it considerably reduces noise produced by this naive collapsing of unique ASVs and allows for remote BLAST, which is impossible otherwise because of the large number of sequences produced by fastx-collapser. The drawback is that it may artificially reduce the number of reads belonging to taxa of interest, if there are many sequences at low read counts that all correspond to the same taxa, common with noisy or low-quality data, also seen commonly when you have very high overall read counts (>100K reads/sample). Try both methods and compare!
- _**percent identity**_ An NCBI metric for comparing two sequences. Literally the hamming distance between two sequences, ie, how many changes are required before they are identical.
- _**bitscore**_ An NCBI metric for assessing how good a BLAST hit is. Roughly, it is the log(2) normalized raw-score. Raw-score is the requred size of a database in order to for a match to occur by chance. 

## Quickstart Guide:
1. If necessary install miniconda for linux-64 in your home directory (`cd ~`)

`mkdir -p ~/miniconda3`

`wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh`

`bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3`

`rm -rf ~/miniconda3/miniconda.sh`

`~/miniconda3/bin/conda init bash`

`~/miniconda3/bin/conda init zsh`

2. Exit the cluster and log in again. This finalizes the conda install. Next, add conda package channels and install the pipeline dependencies:

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

`conda config --set channel_priority strict`

`conda create -n pipeline bioconda::blast bioconda::fastx_toolkit conda-forge::biopython anaconda::pandas bioconda::pear`

`conda activate pipeline` #just to check that it installed correctly

3. If necessary, go to https://account.ncbi.nlm.nih.gov/ and make an account. Then go to Account Settings, scroll down, and generate an NCBI API key. 
4. make a project directory (`mkdir your_project`) and add folders that contain your trimmed sequence files for each gene or primer set (ie ITS, 12S, trnL)
5. cd into your project directory (`cd ./your_project`) and do all of steps 5-9 from that directory.
6. set up the scripts and conda environment by typing the following:
   `module load git` #or ask your administrator about how to clone a git repository on your cluster, you may need to install git
   
   `git clone https://github.com/sckieran/ediet_pipeline/`
   
   `cp -r ./ediet_pipeline/scripts/ .` #to get the scripts into your project directory

   `cp ./ediet_pipeline/pipeline_wrapper.sh .` #to get the wrapper into your project directory

    `cp ./ediet_pipeline/slurm_template.txt .` #to move the slurm template into your project directory
   
8. add a taxa list for each gene to your project directory. They should be identically named except for the gene name (taxlist_12S, taxlist_ITS, taxlist_fwh1). They can be identical, but don't have to be. Each entry in the taxa list must contain two words, genus and species, separated by a space. If you want to search a whole genus, put "genus sp." Do not include a header line.
9. Add a list of gene terms to search for each gene/primer set. Scan NCBI for common variation on gene names. Each gene should be its own column. Columns do not have to be the same length. Each column should have a short header with no special characters or spaces, ie 12S/ITS. That header will not automatically be considered a search term, so add it to the column itself if you want it to be included. See example genelist.
10. Set your parameters by editing the pipeline_wrapper.sh script using a text editor like nano, vim or emacs (see "Parameter Setting" in this tutorial)
11. Set your slurm template by editing the slurm_template.txt file to correspond to your cluster. Internal slurm will have their own err/outfiles and jobnames, so don't add those. But do add anything your cluster needs you to submit - often a partition name, sometimes minimums for node requests, etc. No need to add tons of memory to the jobs.
13. If you wish, copy the ncbi_lineages_2024-01-23.csv.gz from ./ediet_pipeline/ to you project directory. This will use this static version of the NCBI lineages from January 2024. I will periodically update this static file, but if you're working on organisms that may have been added to NCBI very recently, you may want to skip this step and download a new version using ncbitax2lin.
15. From your project directory, run the program: `sbatch pipeline_wrapper.sh`
    
## Parameter Setting
This pipeline contains a wrapper script with all the parameters you'll need to set to run it successfully. Here is a breakdown of those parameters:

- **dir=**$PWD :this is the containing directory for your project. It should contain your scripts folder, your pipeline_wrapper.sh script, and directories called gene1...geneN with your fastqs. You can use $PWD for this parameter as long as the directory structure is correct.
- **prefix=**your_project :pick a name for your project, will be used to name your outfiles. NO SPECIAL CHARACTERS OR SPACES ALLOWED IN THIS NAME. This should be the name of your project directory, although it doesn't have to be.
- **env_name**=pipeline :the name of the conda environment you installed containing the pipeline prereqs.
- **genelist=**$PWD/genelist #the absolute path to your genelist, see example file for format. If your genelist is stored in your project_directory and named genelist, you can use $PWD/genelist for this parameter.
- **taxlist=**$PWD/taxlist :the absolute path to your taxlist, see example file for format. If your taxlist is stored in your project_directory and named taxlist, you can use $PWD/taxlist for this parameter.
- **retmax=**20 :How many sequences per taxon to return from NCBI. Default is 20. Recommended values: 5-100. Bigger values are recommended if the taxlist includes genera or higher-order taxa.
- **db_dirr=**reference_database #name (not path) of the directory you want to contain your reference database. Default is 'reference_database'.
- **key=**YOUR_NCBI_KEY #mandatory: your NCBI key for querying NCBI's entrez database. See installation instructions to obtain one. It is free.
- **R1_pattern=**"_R1.fastq" #the pattern that ends your forward read fastq files. Default is "_R1.fastq".
- **R2_pattern=**"_R2.fastq" #the pattern that ends your reverse read fastq files. Default is "_R2.fastq".
- **max_jobs=**10 #the maximum number of SLURM jobs to submit. Default is 10. Max jobs for RCDS users is 500. Play nice with others, even cutting the tasks into only 10 jobs (instead of 1) will reduce processing time to <24 hours even for very large datasets.
- **extra_seqs=**"extra_seqs" #a name for any number of files that include extra sequences not pulled from NCBI that you want included in your database. Full name **must* match the name you give here with this format: **extra_seqs_${gene}_sequences.fasta**. Leave it in your project_directory and it will be automatically included in your reference.
- **filter=**FALSE #do you want to filter your data by taxa after BLAST and classification? Default is FALSE.
- **asv_rra**=0.005 #the relative read abundance cutoff for ASV filtering. If you don't want to do any ASV filtering, set this value to 0 but you MUST set remote=FALSE and remote_comp=FALSE.
- **taxa_rra=**0.005 #what should your relative-read abundance cutoff be per taxa? Default is 0.005, or 0.5% of total reads in a sample. See relevant literature in your system/gene for recommended cutoffs.
- **identity_cutoff**=97 #what do you want the percent identity cutoff to be? See relevant literature in your system/gene for recommended cutoffs.
- **minlen=**70 #ASVs shorter than this length are discarded before BLASTing. Set to 1 for no filtering, expect processing times to increase. If you need to set separate minlens for each gene, we recommend pre-trimming in cutadapt (-m $minlen) when you perform your adapter trimming.
- **remote=FALSE** : if this is set to TRUE, the pipeline will perform ONLY remote BLAST. Not compatible with asv_rra=0. Not recommended for >10,000 ASVs.
- **remote_comp**=FALSE :if this is set to TRUE, the pipeline will perform both remote and local BLAST and return both results. Not compatible with filter=TRUE or asv_rra=0.
- **return_low=**TRUE #should BLAST return all results, or only those that hit above your $percent_identity cutoff? Default is TRUE. If FALSE, sequences with no BLAST hits above threshold are returned as "No Hit".
- **user=**your_username :your RCDS username, for checking if jobs are done.
- **email=**your_ncbi_email@email.com :the email associated with your NCBI account. Not mandatory, but throws lots of warnings if you don't include.

Once the wrapper is filled out, simply comment out any steps you want to skip (ensuring the mandatory outfiles and directories from those steps exist) and run:

`sbatch pipeline_wrapper.sh`

## Inputs and Preparation

### Inputs
This pipeline requires the following inputs:
- Demultiplexed, appropriately-named **paired-end metabarcoding data in fastq (or fastq.gz) format**, with adapters trimmed. Each gene or primer set should be analyzed separately, and should be in a different folder in your project directory called gene1, gene2...geneN.
- The output of a single run of **ncbitax2lin**. This script will attempt to install and run ncbitax2lin if it can't find a taxonomy file. If you are having issues with permissions and install, you can run it yourself by copying the commands from [this readme](https://github.com/zyxue/ncbitax2lin), the script assumes a name of "ncbi_lineages_[date_of_utcnow].csv.gz" but is agnostic to the date.
- A **taxa list file** containing the scientific names of your target taxa, one per line, with no header, if you want to run step 1 and query NCBI for your fasta.
- A **"genelist" file** with lists of common terms for your target genes. One gene per column, as many permutations on the gene as you'd like, one per line (ie, "Cytochrome Oxidase I", "COI", "COX1"). See the example files for a template. The columns should have a header (which you can repeat in the body) that is a short, human-readable name of the gene that contains no spaces, slashes, quote marks or other special characters. For example, instead of heading your column "Cytochrome Oxidase I", head it "COI". 

Data is often demultiplexed by the sequencing service company at no (or minor) cost. However, if your data has not been demultiplexed, we recommend using either [fastq-multx](https://github.com/brwnj/fastq-multx) or the demux-by-name function of [BBMap](https://github.com/BioInfoTools/BBMap). Some demultiplexers (fastq-multx, for instance) do automatic trimming, others do not, which is why this pipeline only accepts pre-trimmed data. We recommend using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) to trim. 

    
**Providing an NCBI API Key** 

Currently, you must provide Step 1 with an NCBI API key. Getting an NCBI account is free: [sign up here](https://account.ncbi.nlm.nih.gov/signup/). Once logged in, click your name in the upper right hand corner of the screen, click "Account Settings" and scroll down to the button that says "generate API key". Copy this key into the wrapper script and you'll be good to go.

**A note on organizing your data**
To maximize pipeline success, here is an example directory organization for a project that has one sample (sample_1) sequenced at COI and 16S:


![pipeline_directory_org](https://github.com/sckieran/Metabarcoding_Pipeline_Waits/assets/53580356/81b8f9a6-9251-4868-9509-a6348187704d)


You should avoid spaces and special characters (^,$,%,@,#,!,*) in the names of your files/folders and check for hidden characters like carriage returns (sometimes displayed as ^M or \r) in your filenames. Taxlist is only required if you're building a ref database. genelist is required. You can name your genes/loci/primer sets anything (again, avoiding spaces and special characters), but must provide a single-line, tab-separated list of the terms in a file. If you're building the local reference database, the program will automatically pull these term from the header of your gene search terms list. Your data files must be separated by gene in folders that correspond to the terms in the "genelist" file. This is even true if you only have one marker.

If you simply provide a project directory (folder) containing the step 1 input files and also containing a directory called gene1 (through ...geneN) that contains your fastas, the pipeline will create the rest of the directories for you within that project directory.


# FAQ
**What if I want my database to include every organism available for X gene?**

This pipeline is for limited-taxa reference databases or remote querying of NCBI. It is not designed to manage the curation of a large database that includes, say, all inverts with COI sequences in NCBI. If you want a larger database, consider using another reference database building tool, for instance [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) or [bcDatabaser](https://bcdatabaser.molecular.eco/). This pipeline is being actively developed to better improve performance with other databases. Any FASTA file that contains unique sequences each with a valid NCBI taxid will work as long as they are in the format >Accession genus species taxid=[valid NCBI taxid].


**What if I have a few of my own in-house sequences I want to add to my database?**

That's fine, but it will add some tedious work upfront. First, check if the species you are including have taxids assigned in NCBI. All taxa with at least one sequence in the NCBI nucleotide or sra database have a taxid, and you can look it up here: [here](https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi). Put your homebrewed sequences in a file called extras_gene1_sequences.fasta with the header format ">[any unique non_NCBI accession] genus species taxid=unique_taxid" and put it in your project directory. If your species has a taxid assigned in NCBI, use that as the taxid. Otherwise, you can assign one yourself, I recommend a series of letters and numbers separated by underscores. Do not use a strictly numeric value, you will almost definitely accidentally assign a real taxid to your species.  Next, in your project directory, install and run ncbitax2lin. Then, unzip and modify the ncbitax2lin file. Add the unique_taxids you assigned to the end of the file and add taxonomic information for each species following the comma-separated format of the ncbitax2lin file (described in the header of the file). Then re-zip the ncbitax2lin file and you're readdy go fo. This is tedious to do by hand, so we really recommend submitting your homebrew sequences to NCBI if they're long enough!

**What if I want to use remote BLAST to query all of NCBI?**
Select "remote=TRUE" or "remote_comp=TRUE" in the parameters section of the pipeline_wrapper.sh script. But you cannot combine this with an asv_rra=0, because you will end up with many hundreds of thousands of sequences to BLAST, and remote BLAST is very slow and picky. NCBI will be _unhappy_ if you try to remote BLAST more than about 10K ASVs. There is nothing I can do to change this. They really just don't want you to ping their servers 500,000 times.

**Why doesn't this pipeline use e-value instead of bitscore or pident?**
E-values help asses likelihood of homology and are affected by database size. The conserved regions of metabarcoding primers almost guarantee homology above any standard e-value cutoff, and these databases are generally small, so e-values are high across the board.

**I want to do 16S bacterial metagenomics. Is this pipeline right for me?**
Probably not. There is an immense, incredibly well-maintained array of resources for 16S microbial amplicon sequencing, some commercial and some open source, that will be infinitely better than this pipeline.


## Walkthrough of Each Step
### Step One: Fetch FASTAs

**Scripts:**
- step_1_get_seqs_for_database.sh
- query_rentrez.py

**What This Step Does:**

This step does the following, in order:
1. Parses your arguments and makes relevant output directories
2. Executes the script query_rentrez.py, which:
3. Executess a loop using biopython's entrez capability that searches NCBI for each taxa and each gene term you give it. So the search for the first taxa in the example files would query NCBI for "Aplodontus rufia[ORGN] AND ("12S OR 12s Ribosomal RNA OR 12S RNA OR 12S Mitochondrial"). It does this for each taxa, one a at time. It stores the following:
   * Whether or not there are sequences in genbank for that set of gene terms
   * How many sequences there are
   * The FASTA header and sequence for the first 20 (default) or retmax (set by user) matching sequences in NCBI. These sequences are written to a FASTA file in your out directory with the name "taxa"_"gene"_sequences.fasta. Each sequence has the header >Accession genus species taxid=[ncbi taxid].
4. Creates a summary file with the above information.
5. If genus_search is set to TRUE, it checks for any genera where every species of that genus in the taxlist had no hits in NCBI. It then expands the search of that genus to any species in the genus and stores the hits as "genus sp".

That's it for Step 1. We wanted to make sure that users had the opportunity to add their own, potentially off target samples to their databases For instance, there shouldn't be Homo Sapiens hits to trnL or other plant-specific loci, but they may still show up in your data from contamination in the sequencing lane. Similarly, the European Carp genome (Cyprinus carpio) and lots of COVID sequences are full of Illumina adapters, so including a few can help you identify primer dimers and untrimmed sequences in your data.

### Step Two: Validate FASTAs and Make Ref Database

**Scripts:**

- step_2_make_database.sh
- step_2_p1_rmdups.py

**What This Script Does:**

This step does the following, in order:
1. Parses arguments and makes out directories
2. Moves your "extra_seqs_gene1_sequences.fasta" files into your "reference_database" directory to be validated and incorporated into your database
3.  Concatenates all taxa FASTAs for each gene into one large FASTA.
4.  Checks for (and removes) duplicate sequences with the script step_2_p1_rmdups.sh. This is because sequences with identical headers cause makeblastdb to throw an error.
7. Builds your NCBI-formatted reference database using the `ncbi-blast` function `makeblastdb` with the following parameters: -db_type nucl (options are "nucl" and "prot") -in your_project_gene1_database_sequences.fasta -out your_project_gene1_reference -parse_seqids -blastdb_version 5
8. Moves all the extra taxa-specific FASTA files into folders to keep things tidy.

End result is a set of reference database files (10 of them, created by makeblastdb) that NCBI blast can use as a reference, along with a fasta file for each gene, all in a folder called reference_database (or a name supplied by you) inside your project folder


### Step Three: Merge forward and reverse reads

**Scripts**

- step_3_pears.sh
- pear.sh

**What This Script Does**

This script does the following, in order:
1. Parses arguments and create out directories for one gene at a time ("gene1" in this tutorial)
3. Grabs all the files in your data directory that match the forward/reverse patterns you provide (default: _R1.fastq and _R2.fastq)
4. Assesses the number of jobs it should create based on your max_jobs and the number of samples in your gene1 directory
5. Makes N lists of samples. Each list contains X samples, where X=(num_samples / max_jobs) and N=max_jobs.
6. Uses the software `pear` to merge your forward and reverse reads.
7. Checks if your jobs are done, and when they are:
8. Outputs the results of `pear`, one for each sample, to sample_paired.assembled.fastq

The end result is a single paired.assembled.fastq file containing merged reads for each sample, along with the reads that were not assembled (_paired.unassembled.forward.fastq and _paired.unassembled.reverse.fastq) and the discarded merged reads (internal pears defaults discard mergers <50 bp, to disable this, or to add additional filtering to pears, edit the file `pears.sh` to include your desired parameters, potential parameters can be found in the pears manual [here](https://cme.h-its.org/exelixis/web/software/pear/doc.html)

### Step Four: Collapse unique sequences into ASVs

**Scripts**

step_4_collapse.sh
run_collapser.sh

**What These Scripts Do**

Step 4 does the following:
1. Creates jobs as in step 3, based on your max_jobs and the number of samples found that match *_paired.assembled.fastq
2. Cleans up your unpaired, unassembled and discarded fastqs
3. Runs fastx-collapser on each sample, which functionally pulls the unique sequences out from each sample as a fasta file, naively and with no error correction, counting the number of reads per ASV.
4. Removes singletons. This prevents major ASV overflow issues later and helps control for sequencing errors.
5. Checks if your jobs are done
6. Cleans up your outfiles, collapsed outfiles have the name sample_collapsed.fasta

### Step Five: Make per-sample sequence files

**Scripts**

step_5_mk_seqfiles.sh
filter_rra.R

**What This Step Does**

Step 5 does the following, in order:
1. Reformats your collapsed fastas into per-sample ASV tables with the format "sequence  reads", tab-separated
2. For each sample, runs an R script that filters your ASVs by your asv_rra cutoff. For instance, if there were 100,000 total reads in a sample, and your asv_rra=0.001, then it would remove ASVs from a sample's seqfile if it contains <100 reads (100,000*0.001=100). It also excludes sequences shorter than your specified minimum length (minlen=, default=70 bp). Note that pears automatically throws out merged reads shorter than 50 bp, so that is a de facto minimum value for minlen.
3. Cleans up your outfiles and ensures you're ready for BLAST, stored unfiltered seqfiles in their own folder so you can always go back to them later.


### Step Six: BLAST and Classify

**Scripts**

new_step_6_blast_local.sh
new_step_6_blast_comp.sh
new_step_6_blast_remote.sh
get_best_hits.R
get_best_hits_comp.R
get_best_hits_remote.R

**What This Step Does**

Step 6 does the following, in order:
1. Copies your per-sample sequence files into a directory called "gene1_out"
2. Combines the sequences from all samples into one master file, calls uniques and filters out anything <minlen.
3. Creates a BLAST-formatted fasta from your unique ASV master file. Sequences are assigned a unique sequence number formatted as seq_000001 through seq_999999. File is called your_project_gene1_combined_ASVs.fasta
4. Parses your taxonomy file, or attempts to install and run ncbitax2lin to produce a taxonomy file.
5. BLASTs your BLAST-formatted query fasta against your reference database. If "return_low" is set to TRUE, returns all the hits with the highest bitscore (no matter how many are equally good), plus the hits from the next two highest bitscores. If return_low is set to FALSE, returns only hits with a percent identity above your identity_cutoff. Runs on 4 threads. If remote_comp is set to TRUE, also runs remote BLAST and assess best hits. If remote is set to TRUE, it will ONLY do remote BLAST.
6. Parses your raw BLAST results to add convenient taxonomy information (taxid, species) from description of the query sequence, and runs the R file get_best_hits.R (or get_best_hits_comp.R, get_best_hits_remote.R)
7. Joins the NCBI taxonomy file to the raw blast output
8. Summarises each sequence's best hits to identify how many species, genera, families, etc. are in the best scoring hits. Multiple equally-best-scoring hits occur when multiple species have identical sequences at the amplicon region.
9. Returns a "best hit" entry for each sequence that includes the sequence number, the length of the alignment, the percent identity of the top hit, the assigned taxa, the taxonomic information for that hit (phylum, order, class, family, genus), the bitscore of the top hit, the number of species present in the potential best hits, and a comma-separated list of each species that was an equally-good hit.
10. Uses the sample sequence files (_seqs.txt) and the combined_ASV.fasta file to associate sequences, assigned sequence numbers, reads, and samples.
11. Joins the sample and sequence table with the best hits table to include the blast and best hit results for every sample/ASV in the dataset.
12. Produces taxa tables for local blast, remote blast, compared remote/local blast, and summary tables that show sample by species (ie, ASVs are collapsed by species assignment) and species (ie, samples are collapsed and information about the number of reads/samples for a taxon in the entire dataset are returned).


### Step eight: Filter Taxatable -not currently functional

**Scripts**

step_8_filter_data.sh
filter_id_taxa.R

**What This Step Does**

1. Feeds your unfiltered taxatable into R
2. Loops over each sample and does the following:
  3. uses the `dplyr` package in R (pre-loaded onto the RCDS R installation, but also included in `tidyverse`) to filter your data by identity_cutoff.
  4. Groups sequences by taxa of best hit. Only sequences with percent identity > identity_cutoff are included in this grouping.
  5. Makes a list of taxa with total reads > cutoff, where cutoff=$( taxa_rra * total_reads_in_sample ), the "passing_taxa" list
  6. Filters the sample's sequence to include sequences from only taxa in the "passing_taxa" list
7. Re-builds the taxatable to include all samples
8. Produces outfile your_project_filtered_taxatable.txt and cleans up outfiles.

This step runs in sequence rather than by submitting jobs. It is usually pretty quick. 200 samples and 61k unique sequences takes about 4 minutes on the RCDS cluster.

