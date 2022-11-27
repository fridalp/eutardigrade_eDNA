# Tobias Guldberg Frøslev 2022-11-24

# script for filtering the soil eDNA metabarcoding data for 18S to contain potential tardigrade sequences for further analyses
# sequence data has previously been published and is accessed from https://github.com/tobiasgf/sample_storage
# script requires VSEARCH to be installed

# First we read in some libraries/packages we need ---- 
library(here)
library(dada2)

#download the 18S data from the Biowide project
#  https://github.com/tobiasgf/sample_storage/raw/main/asv_tables_referencedata/bw_tar_seqtab.nochim_Both.rds

#checking some basics and fixing some names and order of samples
bw <- readRDS(here::here("bw_tar_seqtab.nochim_Both.rds"))
dim(bw) #45468 OTUs (and 130 samples)
sum(bw) # 8557798 reads in total!
bw_sorted <- bw[order(as.numeric(substr(rownames(bw),3,5))),] #lets sort the table according to the sample numbers (1-130)
#actually a few of the sample names are wrong, so we fix them
rownames(bw_sorted) <-  c("NV001", "NV002", "NV003", "NV004", "NV005", "NV006", "NV007", "NV008", "NV009", "NT010", "NT011", "NT012", "NT013", "NT014", "NT015", "NT016", "NT017", "NT018", "NH019", "NH020", "NH021", "NH022", "NH023", "NH024", "NH025", "NH026", "VU027", "VU028", "VU029", "VU030", "VU031", "VU032", "VU033", "VU034", "VU035", "VO036", "VO037", "VO038", "VO039", "VO040", "VO041", "VO042", "VO043", "VD044", "VD045", "VD046", "VD047", "VD048", "VD049", "VD050", "VD051", "VD052", "EM053", "EM054", "EM055", "EM056", "EM057", "EM058", "EM059", "EM060", "EM061", "ES062", "ES063", "ES064", "ES065", "ES066", "ES067", "ES068", "ES069", "ES070", "EV071", "EV072", "EV073", "EV074", "EV075", "EV076", "EV077", "EV078", "SN079", "SN080", "SN081", "SN082", "SN083", "SN084", "SN085", "SN086", "SN087", "SV088", "SV089", "SV090", "SV091", "SV092", "SV093", "SV094", "SV095", "SM096", "SM097", "SM098", "SM099", "SM100", "SM101", "SM102", "SM103", "SM104", "FF105", "FF106", "FF107", "FF108", "FF109", "FF110", "FF111", "FF112", "FL113", "FL114", "FL115", "FL116", "FL117", "FL118", "FL119", "FL120", "FM121", "FM122", "FM123", "FM124", "FM125", "FM126", "FM127", "FM128", "FM129", "FM130")

#function for filtering the data against a reference database
filter_dada2_against_refdb <- function(dada2_table, query_coverage = 0.95, match_cutoff = 0.9, ref_db = "", out_fasta = "filtered_otus.fasta"){
  require("here")
  require("dada2")
  require("tidyverse")
  asv_file_namex <- tempfile(pattern = "", fileext = ".asvs")
  hit_file_name <- tempfile(pattern = "", fileext = ".hits")
  
  uniquesToFasta(getUniques(dada2_table), fout=asv_file_namex, ids=paste0("Seq", seq(length(getUniques(dada2_table)))))
  #perform match against refence database
  commandX <- paste0('vsearch --usearch_global "',asv_file_namex,'" --threads 40 --dbmask none --qmask none --query_cov "',query_coverage,'" --rowlen 0 --notrunclabels --userfields query+id+target --maxaccepts 0 --maxrejects 0 --top_hits_only --maxhits 1 --db "', ref_db, '" --id "', match_cutoff, '" --iddef 0 --userout "', hit_file_name, '"')
  
  system(commandX)
  
  #read in hits
  hits <- read_tsv(hit_file_name, col_names = F)
  
  #restrict table to filtered OTUs
  hits <- as.data.frame(hits)
  names(hits) <- c("otu","match","best_matching_ref")
  positions <- as.numeric(gsub("Seq","",hits[,1]))
  filt_tab <- dada2_table[,positions]
  hits$otu <-   paste0("eutardigrade_otu_", sprintf('%0.3d', 1:length(getUniques(filt_tab))))
  hits$sequence <- colnames(filt_tab)
  hits$length <- nchar(hits$sequence)
  uniquesToFasta(getUniques(filt_tab), fout=here::here(out_fasta), ids=paste0("eutardigrade_otu_", sprintf('%0.3d', 1:length(getUniques(filt_tab)))))
  filt_tab2 <- filt_tab
  colnames(filt_tab2) <- hits$otu
  result <- list(dada2_table = filt_tab, named_otu_table = filt_tab2, filter_info = hits)
  return(result)
}

# filter the biowide 18S data against a custom tardigrade reference database
filtered_data <- filter_dada2_against_refdb(bw_sorted, match_cutoff = 0.87, ref_db = "/PATH_TO_FILE/Updated_reference_Eutardigrada_Mini18S.fas")

#save the results
saveRDS(filtered_data, "full_biowide_tardigrade_filtered_v2.rds") # full result file
write.table(filtered_data$dada2_table, "full_biowide_tardigrade_filtered_v2_dada2_table.txt", sep = "\t", quote = F, row.names = F) # filtered table with sequences (dada2 format)
write.table(filtered_data$named_otu_table, "full_biowide_tardigrade_filtered_v2_otu_table.txt", sep = "\t", quote = F, row.names = F) # filtered table with otu names 
write.table(filtered_data$filter_info, "full_biowide_tardigrade_filtered_v2_filter_info.txt", sep = "\t", quote = F, row.names = F) # file with sequence, best match, length, etc...

# the rest has been done more manually, outside r.

