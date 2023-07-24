<p align="left">
<img src="https://github.com/popphylotools/Phylosweeper/blob/main/Phylosweeper_logo.jpg" width="700" />
</p>

Created By: Carlos Congrains and Scott Geib. Email: carloscongrains@gmail.com

###################################

The Phylosweeper is a program written in python which aims to generate accurate alignments of coding sequences for phylogenetic analyses. This program analizes alignments of coding sequences to remove falsely inferred orthologs, misassembled sequences, and poorly alignment regions. This script performs three steps of filtering following by a round of alignment and trimming, which are described in the following figure. 

<p align="center">
<img src="https://github.com/popphylotools/Phylosweeper/blob/main/Phylosweeper_workflow.jpg" width="700" />
</p>

REQUIREMENTS:
-----
We developed and tested this script using Python version 3. It requires numpy and Biopython modules. The following programs must also be previously installed: 

* TranslatorX v1.1 (http://translatorx.co.uk) - Abascal et al., 2010.
* trimAl 1.2rev59 (http://trimal.cgenomics.org) - Capella-Gutierrez et al.,2009.
* distmat tool from EMBOSS:6.6.0.0 (http://emboss.open-bio.org/html/use/ch02s07.html) - Rice et al., 2000.


USAGE:
-----
The script takes an alignment in fasta format as input. The program will create a log file and the outputs will be saved in a directory named as the basename of the input file. The users must set three arguments. The --length_filter is the cutoff of the percentage of length variation. For example, 30 means that sequences with a length greater than 130% and lower than 70% of the average length will be removed. The --distance_threshold is the threshold of the average of uncorrected pairwise distance. Possible values range from 0 to 100. The --proportion_Ns_gaps is the maximum fraction of allowed missing data or gaps. For example, 0.25 means that sequences with more than 25% of missing data or gaps will be excluded.
 
```sh
python Phylosweeper.py --fasta_file ${full_path_phylosweeper}/example/cluster3580_1_1to1.fas --log_file ${full_path_log_folder}/cluster3580_1_1to1.log --out_dir ${path_output_folder} --distance_threshold 25 --proportion_Ns_gaps 0.25  --length_filter 30 --cleanup yes
```

**WARNING**: Please set the full path of the input, log and output.

Tip to run the Phylosweeper in parallel:

The following bash script will create a temporary file containing the basename of the fasta files and it will submit this information to the parallel tool. The option -j sets the number of threads to be used in parallel.

```sh
for fasta_file in ${full_path_phylosweeper}/example/.fas; do path_intermediate=$(echo "${fasta_file%.}"); ClusterID=$(basename ${path_intermediate}); echo ${ClusterID} >> temp_file; done

cat temp_file | parallel -j 2 --bar --no-notice "nice python Phylosweeper.py --fasta_file ${full_path_phylosweeper}/example/{}.fas --log_file ${full_path_log_folder}/{}.log --out_dir ${path_output_folder} --distance_threshold 25 --proportion_Ns_gaps 0.25  --length_filter 30 --cleanup yes"

rm temp_file
date
```

Citation:
-----
Congrains et al 2023.

This script is in the public domain in the United States per 17 U.S.C. § 105
