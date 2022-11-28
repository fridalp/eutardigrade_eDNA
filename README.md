# eutardigrade_eDNA
This repository contains all data and scripts necessary to run the analyses and produce the figures from the study Environmental DNA metabarcoding of Danish soil samples reveals new insight into the hidden diversity of eutardigrades in Denmark.

**Bioinformatic tools:**

VSEARCH v.2.9 (or later) (https://github.com/torognes/vsearch)

IQ-tree v.1.6.12  (http://www.iqtree.org/)

APPLES-2 v.2.0.9 (https://github.com/navidh86/apples)


**1.1	Resources**

There is a GitHub repository with all metabarcoding sequence data from Biowide. The file used in this study is bw_tar_seqtab.nochim_Both.rds

GitHub: https://github.com/tobiasgf/sample_storage/tree/main/asv_tables_referencedata 

The associated environmental data from all 130 sampling sites are available for download here: https://ecos.au.dk/forskningraadgivning/temasider/data 


**1.2	Reproducing the analyses**

The script used for extracting potential tardigrade sequences from the full 18S dataset is called filtering_of_sequences.r 


**1.3	Sequence data**

The 221 extracted sequences are listed here in fasta-format in the file Unfiltered_Potential_eutardigrade_MOTUs_Biowide.fas
The extracted sequences filtered to 96 MOTUs are listed in fasta-format in the file Filtered_Potential_eutardigrade_MOTUs_Biowide.fasta

Information on MOTU reads across the 130 study sites are provided in Table S4, the excel file  Table S4_BioWide_MOTU_site_data.xlsx

This table also includes sequences, sequence length, filtering, and information extracted from the environmental datasheet (see 1.1 Resources). 

MOTU-IDs are based on the numbering given in the file full_biowide_tardigrade_filtered_v2_filter_info.txt


**1.4	Reference data** 

The first reference database used in the script filtering_of_sequences.r contains 368 sequences. 

The alignment containing the 313 selected reference sequences from Genbank are provided in fasta-format in the file dataref_MSA.fasta (using MAFFT with the G-INS-i algorithm (% mafft, reorder, maxiterate 2, retree 1, globalpair input)

The alignment containing both the 313 reference sequences and the 96 MOTUs are provided in fasta-format in the file dataquery_MSA.fasta (using MAFFT (% mafft -inputorder –keeplength --addfragments fragments --auto input).


**1.5	Phylogenetic analysis** 

Maximum Likelihood of dataref_MSA.fasta performed in IQ-TREE version 1.6.12 (Nguyen et al., 2015), with the GTR+F+I+G4 as substitution model and branch support calculated with 1000 SH-aLRT (Guindon et al., 2010) and 1000 Ultrafast (UF) bootstrap replicates (Hoang et al., 2018). 
The resulting backbone tree is available in both treefile (backbone.treefile) and contree-format (backbone.contree).

In the tree-based approach, problematic MOTUs (#012, #015, #033, #050, #053, #068, #112, #102) were removed from the alignment prior to phylogenetic reconstruction, because these caused severe topological disruptions (see the file treebased_w_all_MOTUs.treefile). 

backbone.contree was used as consensus-tree, and branch support was calculated with 10.000 Ultrafast bootstrap replicates. 
The resulting tree is available in treefile-format (treebased_results.treefile)

In the phylogeny-based approach, APPLES-2 use FastTree version 2.1.11 to re-estimate branch lengths on the consensus backbone tree (backbone.contree). 
Phylogenetic placement using APPLES-2 was supported with 1000 Fast bootstrapping replicates (-F -N 1000) and printed only placements with the minimum least square error (--lse).  Minimum Least Square Error (MLSE), Pendant lengths (PL) and support values were extracted from from the output file (jplace-file, see Matsen et al., 2012 for details), and provided in supplementary Table S3. The output jplace-file is also provided (phylogenybased_results.jplace).
