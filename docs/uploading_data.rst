Uploading data
==============

INSaFLU needs: NGS data (fastq reads) (mandatory), Sample metadata (to link each sample to the respective NGS data) (mandatory) and Reference data 
(additional user-restricted reference sequences) (optional)

Uploading Sample metadata and NGS data
++++++++++++++++++++++++++++++++++++++

# Option 1 (Batch)
------------------

A. Go to Samples menu and choose Add Samples from csv / tsv file.
.................................................................

The Sample metadata should be a comma-separated value (.csv) or tab-separated value (.tsv) table containing the columns “sample name”, 
“fastq1” and “fastq2” (mandatory columns to fulfill; NOTE: fastq2 is exceptionally not fulfilled only for single-end data) as well these 
additional variables (that may not be fulfilled): “data set”,”vaccine status”,”week”,”onset date”,”collection date”,”lab reception date”,”latitude”,
”longitude”. Examples of template table files are provided here. Samples names in your account must be unique and exclude full stops and commas 
and fastq file names must be complete and match the ones from the files that you are going to upload later. Users are encouraged to include any 
other columns (maximum X) with metadata variables to be associated with samples (see advantages here) 

B. Go to Samples menu and choose Add fastq files
................................................

Here you can simultaneously upload multiple paired-end fastq.gz files (50M maximum per file), which are automatically linked to the corresponding samples.

# Option 2 (Individual)
-----------------------

Samples’ information and respective NGS data (single-end or paired-end reads in fastq.gz format obtained through Illumina or Ion Torrent technologies) 
can be uploaded to INSaFLU as a batch (option 1) or individually (option 2):

A. Go to Samples menu and choose Add Sample metadata
....................................................

Here you can upload each sample at the time (including metadata and NGS data).


Uploading Reference data
++++++++++++++++++++++++
NSaFLU needs reference sequence files to be used for reference-based mapping (mandatory) or for extra alignment/phylogeny analyses (optional). 

In References menu, INSaFLU provides a set of ready-to-use reference sequences, all publicly available at NCBI, currently including:
i. post-pandemic (2009) vaccine/reference influenza A(H1N1)pdm2009, A(H3N2) and B viruses (from both Northern and Southern hemispheres);
ii. representative virus of multiple combinations of HA/NA subtypes (i.e., H1N1, H2N2, H5N1, H7N9, etc)

The default reference files (FASTA and GenBank formats) have been prepared to fit amplicon-based schemas capturing the whole CDS of the main eight 
genes of influenza virus (PB2, PB1, PA, HA, NP, NA, M and NS), such as the wet-lab pre-NGS protocol (here) for influenza whole genome amplification 
adapted from a RT-PCR assay described by Zhou and colleagues (Zhou et al, 2009, for Influenza A; and Zhou et al, 2014, for Influenza B; 
Zhou and Wentworth, 2012)

Note: If you are using this wet-lab pre-NGS protocol and you want to compare your sequences against a reference available at INSaFLU database, 
no further actions are needed.

Still, you may UPLOAD additional reference files (“.fasta” extension; maximum 20000 bp per file) to the user-restricted reference database. 
If you use this option, you can upload::

    i. multi-FASTA files containing the set of reference sequences that constitute the influenza “whole-genome” sequence of a particular virus 
    (e.g, the combination of the traditional 8 amplicons targeting the 8 eight influenza RNA segments). Each individual sequence must have the 
    precise size of each “intra-amplicon” target sequence that you capture by each one of the RT-PCR amplicons. INSaFLU automatically annotates
    uploaded multi-FASTA sequences upon submission, but, if you prefer, you can also upload (optionally) the respective multi-GenBank file.
    
IMPORTANT NOTES TO GENERATE ADDITIONAL REFERENCE SEQUENCES
----------------------------------------------------------

1. (multi) FASTA format is widely applied to save either nucleotide sequences or peptide sequences. An easy way to handle/generate multi-FASTA 
files is by opening a text file (e.g., NOTEPAD) and paste individual sequences after each header line. The FASTA IDs (after the '>' character) 
represent the individual sequence names. For the sake of simplicity, you may designate each sequence as 1, 2, 3, 4, 5, 6 , 7 and 8 (see example), 
following the traditional influenza segments order (keeping this numerical order is advisable). At the end, you just have to save the multi-FASTA 
file as “.fasta” (please avoid symbols or blank spaces in the file names). 

2. you may generate your multi-FASTA files in order to fit your amplicon schema by simply adjusting (i.e., the whole-genome sequences available 
for download at INSaFLU or at influenza-specific sequence repositories, such as the Influenza Research Database 
(https://www.fludb.org), NCBI Influenza Virus Resource (https://www.ncbi.nlm.nih.gov/genomes/FLU/Database/nph-select.cgi?go=database) 
and EpiFLU/GISAID (https://www.gisaid.org/).

3. INSaFLU requires reference sequences exclusively composed by non-degenerate bases (i.e. A, T, C, or G). As such, please ensure that all 
degenerated bases (e.g., R, Y, M, K, S and W) are replaced by non-degenerate sequences before uploading. The choice of the base used in 
the replacement (e.g., “A” or “G” when replacing an “R”) has no impact on the analysis. It simply means that mutations falling in the 
replaced nucleotide position will be reported taking into account the reference base selected.

ii. single FASTA files containing a particular complete or partial locus sequence (e.g., the traditionally used HA1 sequence of a virus 
representative of a particular clades/group). This can be used in “Extra Alignment/Phylogeny” projects.


Explore your Sample and Reference databases
+++++++++++++++++++++++++++++++++++++++++++

Samples menu displays all information for all loaded samples (Samples’ Names in your account must be unique). Upon submission, INSaFLU automatically 
updates samples’ information with reads quality and typing data (automate bioinformatics pipeline modules “Read quality analysis and improvement” and 
Type and sub-type detection”; see Data analysis in the Documentation). Just explore the “More info” icon next to each sample.

References menu displays all information for all reference sequences available at your confidential session. Both FASTA and GenBank files can be downloaded
by clicking on the displayed links.

