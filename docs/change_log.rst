Change log
==========

This tab includes a list (chronologically ordered) of notable changes in INSaFLU.

2022
-------

January 26, 2021
..........................


**Main changes:**

**- Settings**:

- This tab is now organized by **Sequencing technology** (ONT or Illumina/IonTorrent) and **Module** (e.g., Reads Quality Analysis and Improvement, Classification, Minor variant detection, etc)
- It is now possible to **turn ON/OFF** specific modules.
 Note: Users should turn ON/OFF specific modules and select the software settings before uploading new samples. Still, changes can always be done for specific samples afterwards
 
**- Masking consensus**

- Users can now **mask (i.e., put NNs) specific regions (or sites) of the consensus sequences for all (or individual) samples within a given Project** (check all the possibilities in the updated Project Settings button). This new feature is especially useful for masking the start/end of the sequences or known error-prone nucleotide sites. For ONT data, medaka-derived mutations with frequencies below the user-defined “minfrac” (i.e. Minimum proportion for variant evidence) are now automatically masked with an “N”. 
 Note: All user-defined masked regions are reported in the new Sample_list_settings.tsv table; As before, “Ns” are automatically introduced in low coverage regions at a user-selected coverage cut-off

**Minor changes:**

- Available hyperlinks to Nextclade (https://clades.nextstrain.org/) were updated to automatically link to specific SARS-CoV-2 or influenza (A/H3N2, A/H1N1,B/Yam or B/Vic) analysis. 
- Available options to **add new Samples (metadata and reads)** were now collapsed in a single new Button **“Add Sample”** in Samples menu.
- Tables (.tsv) listing all Samples (and respective metadata and QC statistics) and Projects in the user account can be downloaded using the **new “Download” buttons** added to the respective tabs.
- The former “Sample_list.tsv” provided for each Project is now divided in two tables: **“Sample_list.tsv”** (including metadata, Classification, etc) and  **“ Sample_list_settings.tsv” (including the software settings and user-defined cut-offs applied for each sample.). These and other Project tables can now be download using the **new “Download” button**
- A few sequences of WHO recommended vaccine influenza for the 2021-22 season were made available in the Reference menu.

This upgrade is already available in both INSaFLU free online platform (https://insaflu.insa.pt) and locally instable version https://github.com/INSaFLU/docker.
To update your local docker installation, you need to XXXX



2021
-------

December 11, 2021
..........................

**Updated Classification**: INSaFLU now detects Omicron-like Spike sequences just after reads upload (the classification is provided as “SCoV2_potential_Omicron” (this update was performed on 11 Dec 2021; more details in  https://insaflu.readthedocs.io/en/latest/data_analysis.html#influenza-type-and-sub-type-identification-and-human-betacoronavirus-classification-as-of-march-2020)

July 27, 2021
..........................

- INSaFLU online now provides **direct links for consensus sequences analysis using Nextclade (https://clades.nextstrain.org/)**. For SARS-CoV-2 projects, users just need to click in the "Nextclade" icon available next to the link for downloading individual or AllConsensus (by project) sequences. This option is not yet available in the locally instalable version (docker).

- **INSaFLU now also performs influenza type and subtype/lineage identification, as well as Human Betacoronavirus (BetaCov) identification using Oxford Nanopore Technologies (ONT) read data**. Until this update, this rapid classification (which is automatically performed after reads upload) was only available for Illumina / Ion Torrent reads. 

Other minor changes:

- Sequences markers for Human BetaCoV classification were shortened to better accomodate the classification directly from ONT reads. 

Details about the rationale behind this classification and outputs can be found in https://insaflu.readthedocs.io/en/latest/data_analysis.html#influenza-type-and-sub-type-identification-and-human-betacoronavirus-classification-as-of-march-2020 (see also the list of current genetic markers used for classification).


April 27, 2021
..........................

**INSaFLU now automatically assigns SARS-CoV-2 Pango lineages (https://pangolin.cog-uk.io/)** using Pangolin (https://github.com/cov-lineages/pangolin), as described by Rambaut and colleagues (Nat Microbiol; 5:1403-1407).

This novel feature works as follows:

- Everytime a new sample is added to a Project, the latest pangolin and pangoLEARN versions are automatically run for all samples within the Project.
- Whenever a new Pangolin / Pangolearn version is released*, a button **"Update Pango lineage"** will be automatically made available at the bottom of “Projects” tab, so that users can re-assign all samples in the project using the latest software/database versions (*INSaFLU will check every day whether a novel pangolin/pangoLearn version is available);
- Results (and software versions) are provided in the “Sample_list” and are automatically available for coloring tree nodes (and/or display colored metadata blocks next to the tree) according to the Pango lineage

Other minor changes:

- Trimmomatic version was upgraded, and ILLUMINACILP was made available for user-defined configuration;
- Downsized samples will be flagged in the “Sample_list.tsv”.

NOTE:  Users might need to do CTRL+F5 to activate this new feature.

This upgrade is already available in both INSaFLU free online platform (https://insaflu.insa.pt) and locally instable version https://github.com/INSaFLU/docker. 


March 25, 2021
..........................

**MAJOR UPGRADE – INSaFLU now also handles Oxford Nanopore Technologies (ONT) data**

Available both in INSaFLU free online (https://insaflu.insa.pt) and locally installable (https://github.com/INSaFLU/docker) versions.

In this update, we added these new main features to INSaFLU: 

- **an automate pipeline for ONT data analysis**, from raw reads to quality analysis, reference-based generation/curation of consensus sequences, mutation annotation, gene/protein/genome alignments, phylogenetic tree, metadata visualization… (details about the pipeline, including software version, default settings, etc, can be found in: https://insaflu.readthedocs.io/en/latest/data_analysis.html# ) 

- For enhanced data navigation, **two new interactive and dynamic “expand-and-collapse” panels were added to the Projects: “Mutations list” (lists all validated mutations, i.e., those inserted in the consensus sequences, for all samples); “Coverage for all samples” (provides an additional interactive color-coded coverage report, summarizing the mean depth of coverage and horizontal coverage per locus for all samples within a project)**

- As for the Illumina/IonTorrent data analysis, **INSaFLU allows users to configure key parameters for ONT reads quality analysis, mapping and consensus generation/curation**. Settings can be user-defined for the whole user account (tab “Settings”), for each project (after project creation) or for individual samples within a project (novel “Magic wand” icon) (more info in: https://insaflu.readthedocs.io/en/latest/data_analysis.html#user-defined-parameters) 

- **Mutation annotation (i.e., impact at protein level) and amino acid alignments were improved** (for SARS-CoV-2 analysis, please use the reference sequences “SARS_CoV_2_Wuhan_Hu_1_MN908947” available at the default reference database). NOTE: Protein alignments only include samples with < 10% of undefined amino acids (X).

- A new “Magic wand” icon was added to the Samples menu. It allows re-running reads’s QC for samples that are not inserted in any project (and for which the original reads have not been deleted). This feature overcomes the previous need of uploading the original fastq files to re-run the quality analysis. 

An updated summary of the main INSaFLU outputs is available here:
:download:`INSaFLU_current_outputs_25_03_2021.xlsx <_static/INSaFLU_current_outputs_25_03_2021.xlsx>`

Other minor changes include:

- Samples generated from different technologies (Illumina/Ion Torrent/ONT) can be analysed within the same Project.

- The csv/tsv file with the list of samples in a project (which compiles all samples' metadata and additional INSaFLU outputs) now also **summarizes the software settings and user-defined cut-offs applied for each sample.**

- Analysis of minor variants (Illumina data only): besides the report of a “validated_minor_iSNVs.tab” table per sample/project (listing SNV displaying intra-sample variation at frequency between 1 and 50% - minor variants), INSafLU now also reports an additional minor variants table “minor_variants_inc_indels.tab” per sample, which includes minor “indels”

- The “coverage.tsv” file was also improved.


2020
----


December 19, 2020
.......................

- Corrected an issue in “AllConsensus.fasta” file creation. We detected a bug where “red” flagged samples (not fulfilling user-selected coverage thresholds) were mistakenly included in this file (other outputs, such individual consensus sequences, variants list, alignments and trees were not affected by this bug). The issue is now solved and "AllConsensus.fasta" files were corrected by excluding “red” flagged samples. 

NOTE: If you already used individual consensus sequences (downloaded for each sample) or the alignments combining all validated locus/genome consensus sequences (Alignment_nt_locus.fasta), this bug was not a problem. If you had already downloaded the combined "AllConsensus.fasta" file,  please confirm that you exclude “red” flagged samples from your downstream analyses or, instead, please re-use the novel corrected file.


November 24, 2020
.......................

This update is available in both INSaFLU free online (https://insaflu.insa.pt) and locally installable (https://github.com/INSaFLU/docker) versions.

- Add a new button to delete fastq.gz files that are not attached to any sample ("Remove not processed files") 
- Add a new button to unlock sample metadata tables ("Unlock last file").
- As for nucleotide alignments (see update 30 Oct 2020), amino acid alignments now also include samples with incomplete locus, i.e., undefined amino acids (“X”) are automatically introduced in low coverage regions at a user-selected coverage thresholds. This update will be applied to all novel Projects. Samples within old projects (before this update) will remain unchanged unless any parameter is altered. In that case, the updated samples will be included in the amino acid alignments following the new criteria.


October 30, 2020
.......................

This important update is available in both INSaFLU free online (https://insaflu.insa.pt) and locally installable (https://github.com/INSaFLU/docker) versions.

**Main changes:**

-  INSaFLU now allows users to configure key parameters for reads quality analysis, mapping and consensus generation. Settings can be user-defined for the whole user account (tab “Settings”), for each project (after project creation) or for individual samples within a project (novel “Magic wand” icon). 

- INSaFLU now generates consensus sequences for incomplete locus, i.e., undefined nucleotides (“N”) are automatically introduced in low coverage regions at a user-selected coverage thresholds. Users can select the minimum “vertical” coverage (depth) threshold per site (mincov; default = 10) and the minimum percentage of “horizontal” coverage to generate the consensus sequence (default = 70%). 

- To better accommodate these novel features, the interactive color-coded coverage report by locus was updated to:

GREEN: % of locus size covered by at least X-fold = 100%

YELLOW: % of locus size covered by at least X-fold is ≥Y% and < 100%

RED: % of locus size covered by at least X-fold is <Y%

	X is the user-defined "mincov" value (i.e., the minimum number of reads covering a site to be considered for variant calling) selected for each project or sample (within a project) (default = 10)

	Y is the user-defined "Minimum percentage of locus horizontal coverage (with depth of coverage equal or above X) to generate consensus sequence" value selected for each project or sample (within a project) (default = 70);

**IMPORTANT NOTE:** These novel criteria will be applied to all Projects and Samples. Samples within old projects (before this update) will remain unchanged, unless the users re-run them with novel user-selected parameters. All updated samples and novel samples run from now on will be flagged ("Calendar" icon).

**Minor changes:**

- Consensus sequences can now be downloaded as a batch.

- Tabular coverage reports per sample are also provided for download.



May 06, 2020
..............

- INSaFLU local installation - a Docker version of INSaFLU, which eases the manual installation process, is now available here: https://github.com/INSaFLU/docker

- Multitasking configurations were changed, considerably speeding up the analyses. 

- A new tab “Settings” was created so that the user can change some software parameters.

All updates are available at both INSaFLU docker version and original free website (https://insaflu.insa.pt/)


March 10, 2020
..............

The following updates have been performed so that INSaFLU can better accommodate genome-based analyses of the novel coronavirus (SARS-CoV-2 / hCoV-19):

- INSaFLU now performs rapid assignment of Human Betacoronavirus (BetaCoV), including the novel coronavirus (SARS-CoV-2 / hCoV-19). Details about the rationale behind this classification and outputs can be found in https://insaflu.readthedocs.io/en/latest/data_analysis.html#influenza-type-and-sub-type-identification-and-human-betacoronavirus-classification-as-of-march-2020 (see also the list of current genetic markers used for classification).

- The publicly available SARS-CoV-2 reference genome sequence (NCBI accession number MN908947 https://www.ncbi.nlm.nih.gov/nuccore/MN908947) is available in the default INSaFLU reference database (several sequence versions with differential trimming of the sequence boundaries are available, as these regions might not be captured by your wet-lab NGS strategy). As before, the users can still insert their own reference sequences.  

- Maximum size per fastq.gz file remains 300 MB, but files will be downsized to ~150 MB before analysis (and not ~50 MB, as previously). This change minimizes the risk of losing considerable depth of coverage in your analysis, specially for SARS-CoV-2 genome analysis.


January 15, 2020
................

- INSaFLU now allows you to easily color tree nodes and to display colored metadata blocks near to the phylogenetic trees

This update largely facilitates the visualization, exploration and interpretation of your phylogenetic data, while potentiating the association/integration of relevant epidemiological and/or clinical data and pathogen genomic data towards an enhanced laboratory surveillance. See how to do it here: https://insaflu.readthedocs.io/en/latest/output_visualization.html#b-navigate-through-phylogenetic-trees-and-explore-your-metadata

- INSaFLU also allows you to “Add/update Sample metadata” at any time

To take advantage of the novel metadata visualization tools, you can now add/update the samples descriptive data by simply uploading a comma-separated (.csv) or tab-separated (.tsv or .txt) table with the updated data (a template file is provided in Samples menu / Add or Update Samples from csv / tsv file). Specific documentation can be found here:
https://insaflu.readthedocs.io/en/latest/uploading_data.html#updating-sample-metadata


January 10, 2020
................

- The INSaFLU list of genetic markers "influenza_assign_segments2contigs" was upgraded (now includes 544 sequences). This update allows the rapid assignment of additional representative virus of distinct genetic clades, which, for instance, can facilitate the sub-group HA classification and potentiate the detection of (intra-subtype) reassortments.


Latest database can be downloaded here: :download:`INSaFLU_current_genetic_markers_v5_after_10_01_2020.xlsx <_static/INSaFLU_current_genetic_markers_v5_after_10_01_2020.xlsx>`

All database versions can be found here: https://insaflu.readthedocs.io/en/latest/data_analysis.html?highlight=genetic_markers#type-and-sub-type-identification 


- The default reference database of INSaFLU was also updated. All reference sequences at INSaFLU are publicly available at NCBI (or are made available under permission of the authors). 

Download the current list here: :download:`INSaFLU_current_REFERENCE_DATABASE_10_01_2020.xlsx <_static/INSaFLU_current_REFERENCE_DATABASE_10_01_2020.xlsx>`) 

Instructions to upload additional reference sequences (e.g., "vaccine-like" sequences available in GISAID) to your confidential account can be found here: https://insaflu.readthedocs.io/en/latest/uploading_data.html#uploading-reference-data


2019
----

January 02, 2019
................

- The INSaFLU list of genetic markers "influenza_assign_segments2contigs" was upgraded (now includes 464 sequences), so, from now one, INSaFLU can assign additional representative virus of distinct genetic sub-groups of seasonal A(H3N2) viruses, not only facilitating the sub-group HA classification, but also potentiating the detection of (intra-subtype) reassortments.


Latest database can be downloaded here: :download:`INSaFLU_current_genetic_markers_v4_after_02_01_2019.xlsx <_static/INSaFLU_current_genetic_markers_v4_after_02_01_2019.xlsx>`

All database versions can be found here: https://insaflu.readthedocs.io/en/latest/data_analysis.html?highlight=genetic_markers#type-and-sub-type-identification 


2018
----

October 30, 2018 
.............

- Original reads (i.e., reads uploaded) will now be deleted after 10 days of their upload. In fact, after quality analysis and improvement, the INSaFLU pipeline does not use those original reads for any other downstream analysis (quality reports and derived quality processed reads will remain available for download).


June 29, 2018 
.............

INSaFLU now published in Genome Medicine.

Borges V, Pinheiro M et al. Genome Medicine (2018) 10:46

https://doi.org/10.1186/s13073-018-0555-0


May 14, 2018 
.............

- The INSaFLU list of genetic markers "influenza_assign_segments2contigs" was upgraded (now includes 416 sequences), so, from now one, INSaFLU can assign additional close references sequences to your viruses, such as representative virus of distinct genetic sub-groups or seasonal A(H3N2) viruses or  representative A(H5N1) sequences of distinct H5 genetic clades.


All database versions can be found here: https://insaflu.readthedocs.io/en/latest/data_analysis.html?highlight=genetic_markers#type-and-sub-type-identification 


April 9, 2018 
.............

- Maximum size per fastq.gz file was upgraded from 50 MB to 300 MB. 

	* IMPORTANT NOTE: Files between 50 - 300 MB will be downsized to ~50 MB before analysis by randomly sampling reads using fastq-sample from fastq-tools package https://github.com/dcjones/fastq-tools (developed by Daniel C. Jones dcjones@cs.washington.edu) 

- The draft assembly provided by INSaFLU (FASTA format) now additionally includes potential non-influenza specific contigs (i.e., contigs not assigned to any influenza segment / reference by INSaFLU). This feature allows users to better inspect the draft assemblies and reinforces the applicability of INSaFLU for other viruses.  


March 9, 2018 
.............

- INSaFLU now provides a draft genome assembly (FASTA format) including influenza-specific NODES/contigs. These are identified by screening the SPAdes-derived draft assemblies against an in house database using ABRIcate, which allows assigning NODES/contigs to the corresponding viral segments and to a related reference influenza virus (output: table in ".tsv" format). Please check these new outputs and guide to interpret them at the INSaFLU tab "Samples" / "Extra info" / "Type and subtype/lineage identification". Please also check software settings and parameters at the "Data analysis" tab of this Documentation. 

	This new feature reinforces the application of INSaFLU to:
	
		* analyse viruses for which a close related whole-genome sequence is not available (e.g., avian influenza) at the INSaFLU or other databses (NCBI, GISAID, etc);
		* investigate reassortments
		* disclose mixed infections
	


January 25, 2018 
................

- INSaFLU 1.0.0 is released for the scientific community at https://insaflu.insa.pt 
	
	INSaFLU ("INSide the FLU") is an bioinformatics free web-based suite that deals with primary NGS data (reads) towards the automatic generation of the output data that are actually the core first-line “genetic requests” for effective and timely influenza laboratory surveillance. While INSaFLU has indeed some influenza-specific features (e.g., automatic type/subtype identification), there is no restrictions to use it for other viruses. 

	Main highlights:
    
		* open to all, free of charge, user-restricted accounts
		* applicable to NGS data collected from any amplicon-based schema
		* allows advanced, multi-step software intensive analyses in a user-friendly manner without previous training in bioinformatics
		* automatic identification of influenza type and subtype/lineage, detection of putative mixed infections and intra-host minor variants
		* allows integrating data in a cumulative manner, thus fitting the analytical dynamics underlying the continuous epidemiological surveillance during flu epidemics
		* outputs are provided in nomenclature-stable and standardized formats and can be explored in situ or through multiple compatible downstream applications for fine-tune data analysis and visualization
