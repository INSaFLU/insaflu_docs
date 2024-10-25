**Frequently Asked Questions (FAQs)**
======================================
FAQs
....

**UPLOAD (Input, Classification, Settings, etc)**
-------------------------------------------------

- **Which sequencing technologies are compatible with the platform?**
The platform supports second- and third-generation sequencing technologies, including Illumina (single-end and paired-end), Ion Torrent, and Oxford Nanopore Technologies (ONT).

- **What read file types can be used as input?**
The platform only accepts gzipped FASTQ files with the extension “.fastq.gz”.

- **How can I convert files to GZIP on Windows?**
You can use third-party tools like 7-Zip. Here's a simple guide: https://www.freecodecamp.org/news/how-to-convert-files-to-gzip-on-windows/

- **How many “.fastq.gz” files can be uploaded per sample?**
You can upload 1 file for ONT and Illumina single-end, or 2 files for Illumina paired-end.

- **ONT technology generates several files per sample. Is there an easy way to merge them on Windows or Mac?**
Yes, you can find executable files to merge “fastq” or “fastq.gz” files with a simple double-click at this link: https://github.com/INSaFLU/INSaFLU/tree/master/files_helpful. Available both for Windows and Mac.

- **What are the upload options?**
You can either:
   1. Add One Sample (upload reads and metadata simultaneously)
   2. Add Multiple Samples (upload metadata.tsv first, then DRAG&DROP reads under “Add fastq files”)

- **What fields are mandatory in the metadata input table?**
The mandatory columns are “sample name,” “fastq1” (and “fastq2”, if paired-end Illumina).

- **What is the purpose of the “data set” column?**
The “data set” column allows you to assign a unique ID to a batch of samples, facilitating joint selection and analysis. It also helps organize files by criteria such as sequencing date or run.

- **What about the other columns? Can I add more?**

Yes, users are encouraged to fill in the existing columns or include additional columns with metadata variables related to their samples. For example, including latitude/longitude or country, region, division, and location data will allow geographic placement of samples in Nextstrain Datasets.

- **I’ve uploaded the table (but not the reads) and realized I made a mistake in the metadata. How can I delete it?**
> Samples > Add Multiple Samples > Unlock last file

- **I’ve uploaded both the table and reads but need to delete everything due to a major error (e.g., incorrect “sample-fastq” correspondence). What should I do?**
> Samples > Remove ALL samples (CAUTION: This action will delete ALL samples that have not been added to INSaFLU and TELEVIR projects. )

- **I’ve uploaded a sample but want to change the settings. Do I need to delete everything and start again?**
No. Click in the button next to the “More info” to change the sample-specific settings. 

- **I’ve uploaded everything successfully, but some reads are not attached to their respective samples. What should I do?**
> Samples > Add Fastq files > Try to relink unattached files 

- **Are the draft contigs the final curated genome sequences?**
No, the draft contigs are not the final curated genome sequences. However, they can be highly useful for identifying closely related references (e.g., using BLAST) for INSaFLU mapping, especially for highly diverse viruses like influenza A/H5N1.

- **Does the Classification result provided just after upload reflects the metagenomic identification of any virus present in the sample?**
No. This output is based on a rapid screening of draft contigs (generated right after upload) to identify/classify specific viruses of interest. Currently, it identifies influenza types A and B, all known influenza A subtypes (18 hemagglutinin and 11 neuraminidase subtypes), the two influenza B lineages (Yamagata and Victoria), five human Betacoronaviruses, RSV A/B, and the four clades of MPXV (Ia, Ib, IIa, and IIb). This classification helps in selecting the appropriate references for the INSaFLU module. For actual virus metagenomic detection, run a TELEVIR project instead.

**INSaFLU module** 
------------------

- **Can I combine samples from different sequencing technologies in the same project?**
Yes.

- **Where do I click "Run" to start the analyses?**
You don’t need to. The analysis begins automatically after you add samples (you can select and add samples at any time). Each time a sample is added or deleted, the project outputs, such as alignments and trees, will automatically be recalculated.

- **Can I have more than one reference per project?**
No. In INSaFLU projects, reads are mapped against a single reference, which will guide the position of mutations reported, alignments, etc.

- **Is it possible to add external FASTA sequences (e.g., from GenBank or GISAID) to the draft phylogenetic trees in INSaFLU projects?**
No. To build more robust phylogenetic trees with external sequences, use the Nextstrain module available on the platform.

- **What are the minimum recommended vertical and horizontal coverage values?**
The users should adjust the coverage threshold according to the context of their analysis. In Portugal, for genomic surveillance of SARS-CoV-2, we routinely use:
mincov (minimum number of reads covering a site to be considered): 30x for SARS-CoV-2 (NOTE: the default for new accounts is 10x!)
minimum horizontal coverage to generate a consensus sequence: >89% for SARS-CoV-2 (NOTE: the default for new accounts is 70%).

- **Does the analysis pipeline include primer clipping?**
Yes, INSaFLU uses the iVar strategy for primer clipping. Primer clipping is advisable. To select your primers, go to Project settings and select the appropriate primer scheme.

- **I can't find my primers. How do I add them?**
Currently, it’s not yet possible to upload primers sequences directly to the platform. Please contact us via insaflu@insa.min-saude.pt if you would like to add new primer pools to the online tool.

- **Does the platform record the parameters used for each sample in INSaFLU projects?**
Yes, they are saved in the file Sample_list_settings.tsv.

- **Does INSaFLU automatically assign SARS-CoV-2 lineages?** How does it handle constant updates?**
Yes, INSaFLU assigns lineages using Pangolin (https://github.com/cov-lineages/pangolin) (usher mode). Every time a new sample is added to a project, the latest Pangolin and database versions are automatically applied to all samples within the project. When a new version is released, an "Update Pango lineage" button will appear at the bottom of the old "Projects" tab, allowing users to re-assign lineages with the latest software/database versions. INSaFLU checks for new versions daily. Results (and software versions) are provided in the "Sample_list.tsv" and can be used to color tree nodes or display metadata blocks next to the tree based on the Pango lineage.

- **Why do I see “NNN” (undefined nucleotides) in my consensus sequences?**
INSaFLU consensus generation pipelines automatically places NNN in: 1) regions (or sites) selected to be masked by the user (in Projects settings); The beginning and end of sequences are typically set to be masked according to the primer scheme used.  2) low coverage regions (i.e., regions with coverage below --mincov); 3)
(only for ONT data) sites with mutations with frequencies between 50% and the user-defined “minfrac” (i.e. Minimum proportion for variant evidence; default: 80%).

- **Does the platform detect and report minor variants?** If so, what frequency is required for a mutation to be considered a minor variant?**
For Illumina (currently), INSaFLU reports a list of minor intra-host single nucleotide variants (iSNVs), i.e., SNVs displaying intra-sample frequencies between 1-50%. (se table “validated_minor_iSNVs.tsv”). In a dynamic manner, distinct minimum iSNV frequency cut-offs are assumed depending on the depth of coverage reached at each site, i.e., the identification of iSNV sites at frequencies of 10, 2, and 1% is only allowed if the depth of coverage at a particular site exceeds 100-fold, 500-fold, and 1000-fold, respectively. For each INSaFLU project, a plot showing the proportion of iSNVs at 1-50% (minor iSNVs) and at 50-90% for each sample are also provided to guide the identification of potential mixed infections. Mutations above 50% are considered as “main variants” and  are inserted in the consensus and reported in a separate table (validated_variants.tsv).

- **What are the criteria that trigger the “putative mixed infection” flag?**
For Illumina (currently), INSaFLU flags samples as “putative mixed infections” based on intra-host SNVs if the following cumulative criteria are met:
The ratio of the number of iSNVs at 1-50% (minor iSNVs) to those at 50-90% falls within the range 0.5-2.0.
The sum of the number of iSNVs in both categories exceeds 20. Alternatively, for mixed or co-infections involving highly divergent viruses (e.g., A/H3N2 and A/H1N1), the flag is also displayed when the sum of iSNVs from both categories exceeds 100, regardless of the first criterion.
How can I inspect if my ONT sample contains a mixture of the same-species virus?**
In the “Sample_list_settings.tsv” file, INSaFLU lists nucleotide sites with mutations having frequencies between 50% and the user-defined “minfrac” (which are automatically masked for ONT data). If this list is extensive, you should suspect a potential mixture of same-species viruses in your sample. Look at read mapping profile using IGV.

**Nextstrain module** 
----------------------

- **What is a Nextstrain build in INSaFLU?**
A nextstrain build is a workflow to perform spatiotemporal phylodynamic analysis of a given viral pathogen. We adapted available workflows to run within INSaFLU. Very briefly, input sequences are aligned into a multiple sequence alignment, and a phylogenetic tree is inferred from this alignment (most builds also include the temporal dimension in this process, as well as species-specific inferences, such as clade information).

- **What builds do we provide in INSaFLU?**
We provide a generic build with and without temporal information, and species-specific builds for SARS-CoV-2, seasonal influenza (Hemagglutinin in H3N2, H2N1pdm, B Victoria and B Yamagata), avian influenza (each of the 8 segments), monkeypox (hMPXV - clade IIb) and RSV (RSV A and B).

- **Are the INSaFLU Nextstrain builds exactly the same as the official Nextstrain builds?**
No. While we aim to make the INSaFLU Nextstrain builds as similar as possible to the official ones, there are frequent updates. We modified certain QC steps to include as many user sequences as possible, thus users are responsible for ensuring the quality of their input sequences. Since we only include user-generated sequences, some build-specific inferences, such as clade determination (particularly when dependent on tree position), may be inaccurate. You can find the code used for our builds here: INSaFLU Nextstrain Builds https://github.com/INSaFLU/nextstrain_builds

- **Can you modify parameters internal to the build?**
Currently you cannot modify any parameter related to the build. You can only choose the build.

- **Do I always need to provide a reference for the generic build?**
Yes, currently you always need to provide at least one reference for the generic build. Even though the generic build with temporal information does not explicitly require a reference to be used as root for the phylogenetic tree inference, it still requires a reference for the annotation. If you do not add a reference to a generic build, it will fail. 

- **Can I use any reference for the generic build?**
No. Currently, the generic build (with and without temporal information) only works with references with a single segment. For example, you cannot use an influenza reference (with 8 segments). If you want to use the generic build with a multi-segmented virus, you will need to separate the reference into independent segments and use one single segment as reference. 

- **Can I add more than one reference to the generic build?**
You can add more than one reference, but at the moment you cannot control which one will be used as a reference to the build (only one will be used as reference in the build). So we advise that you only add one reference, and if you need more you can add them as external sequences.

- **What metadata do I need to provide for the build?**
Most builds (except the generic build without temporal information) require that sequences be associated with dates. Spatial information such as location and latitude/longitude is often also provided but is not absolutely required for the builds, and can be added afterwards. Other variables are also not usually required for the builds and can be added even during visualization.

- **What happens if there is no metadata associated with my sequences?**
When adding sequences to a nextstrain project, metadata information will be generated for those sequences for mandatory variables (date is the most common mandatory variable)  Sequences coming from reference-based projects will bring the metadata that is already associated with them, so the date that was associated to their corresponding sample during upload is used for nextstrain. Internal INSaFLU references or external sequences (from uploaded fasta files) do not have any metadata associated with them. If the sequence does not have metadata value attributed, a default value will be used. In the case of dates, the current date (date when the sequence was added to the nextstrain project) will be used. 

- **Can I update the metadata information before running the build?**
Yes, you can and you should update the metadata information associated with sequences in the nextstrain build, particularly for those cases where default values were used. You can download the current metadata table, update any values, and upload the updated metadata table.

- **What software can I use to update the metadata?**
The metadata table is a tab-separated text file. You can edit in any text editor, or in a spreadsheet software such as Excel. If you use Excel, take care of the format of some columns such as the date, which must follow the YYYY-MM-DD format (Excel may automatically change the format). 

- **How do I visualize results from nextstrain builds?**
Nextstrain builds will provide a multiple alignment, a phylogenetic tree, and an auspice file with all the results. The alignment and phylogenetic tree can be visualized directly within INSaFLU, but the auspice file currently has to be visualized using a third party software. Namely, the auspice json files can be imported for visualization in  https://auspice.us/. Note that, although auspice.us is a website, all the data is processed locally in your browser.

- **If I want to add/modify metadata to visualize, do I need to rerun the build?**
If you do not need to modify the temporal information associated with the samples (in which case you will need to rerun the build), you can add or modify the metadata after the build is run, and drag-n-drop the modified metadata to the auspice.us visualization. 

- **I created a nextstrain dataset, added sequences, updated metadata and ran the build. Nonetheless, after finishing I could get no results (no json files in auspice, and no alignment nor phylogenetic tree). How can I find out what went wrong?**
There can be several reasons for build failure. If you do not provide a reference in a generic build, the auspice will not even be generated, and you will get an error file instead. If you get an auspice file, there will be log files inside that you can open with a text editor to find out at which step did the build fail. In many cases the build will fail at the tree building step, due to issues with metadata, such as a lack of diversity in metadata values (eg. all samples with the same date).


**TELEVIR module** 
----------------------

- **What is TELEVIR used for in INSaFLU-TELEVIR platform?** 
TELEVIR is used for hypothesis-free virus detection through classic workflows or targeted screening and validation of viral hits.

- **How do I deploy workflows in TELEVIR?** 
Workflows are deployed by turning ON/OFF workflow steps and selecting software/parameters in the Settings panel, and then clicking the Run icons for either all samples or user-selected ones.

- **Can I combine samples from different sequencing technologies in the same TELEVIR project?**
No. We encourage users to create different TELEVIR projects per different metagenomics sequencing run (including negative/positive controls) for an enhanced sample comparison and output interpretation.

- **What is the difference between Classic and Validation workflows?** 
Classic workflows focus on initial virus detection, while validation workflows verify viral hits or target specific viruses through mapping options.

- **Can TELEVIR detect all viruses with a single workflow?** 
Based on extensive benchmarking, we found that there is no "one-size-fits-all" approach. TELEVIR allows users to combine different workflows (relying on diferente software and databases) to improve and potentiate virus detection.

- **What is the purpose of the TELEVIR Reference Management Panel available for each Sample?** 
It displays classification outputs and offers options to map against panels or selected references for validation workflows.

- **How are viral hits classified in TELEVIR?** 
Viral hits are classified based on reads and contigs and undergo reference-based mapping. Hits that don't meet remapping criteria are flagged as "Unmapped."

- **How are references ranked in TELEVIR Reference Management?**
References are ranked using Global Ranking (based on the TELEVIR sorting algorithm which aggregates number of read and contig hits obtained by classifiers) and Ensemble Ranking sorts (sorts references according to their average TELEVIR sorting rankings across workflows they were identified in)

- **What is a Map Panel in TELEVIR?** 
A Map Panel refers to a set of viral references created in the TELEVIR References Management page, which can be used for mapping workflows.

- **What is TELEVIR-Focus, and how does it work?**
TELEVIR-Focus allows users to map multiple samples against a single reference, enabling detailed analysis of a given virus across multiple samples.

- **Can I run multiple validation workflows on several samples simultaneously?**
Yes, TELEVIR allows multiple samples to be mapped against specific viruses simultaneously through its "Map Added" or "Map Panels" options.

- **What are automatic warnings in TELEVIR reports?** 
TELEVIR provides automatic warnings like “Likely False Positive” and “Vestigial Mapping” to help identify potential issues in viral hits.

- **How does TELEVIR handle negative control samples?** 
Negative controls (e.g., buffers) should be included to detect contamination, and viral hits in controls are flagged in the main report.

- **What are the advantages of using multiple control samples in TELEVIR?** 
Using multiple control samples helps improve accuracy by identifying contamination sources and minimizing false positives.

- **How can I interpret the TELEVIR reports?** 

TELEVIR reports provide virus detection results at different levels: per workflow, per sample, and per project, with options, metrics, and downloadable files to further investigate flagged hits. 

- **How does TELEVIR handle closely related viral hits?** 
Cross-mapping occurs between closely related viruses. The virus with the best mapping metrics is likely the true positive.

- **Why do I see multiple hits for the same virus (TAXID)?**
This can happen due to segmented viruses (where each segment has a different accession number) or because there are multiple reference genomes for the same virus, with the closest match yielding the best mapping metrics.

- **What causes multiple hits for closely related TAXIDs?**
Cross-mapping of reads across genomes with high nucleotide homology, such as viruses from the same family, can lead to multiple hits. Interactive Heatmaps for Cross-mapping inspection allow users to easily inspect these situations.

- **What if my suspected virus isn't listed in the Main report?**
Check the Reference Management panel. If the virus is flagged as "Unmapped," request confirmatory mapping or adjust filtering steps like Viral enrichment. You can also map directly against the suspected virus, either through Map Added or Map Paneç

- **Can I manually inspect viral hits flagged by TELEVIR?** 
Yes, flagged hits can be further examined using internal tools like Integrative Genomics Viewer (IGV) or external resources, like BLAST.


Guide for pre-NGS steps
........................

Suggested pre-NGS wet-lab protocol for influenza
-------------------------------------------------

The reference-based surveillance-oriented component of INSaFLU (https://insaflu.readthedocs.io/en/latest/routine_genomic_surveillance.html#reference-based-genomic-surveillance) is highly flexible and **allows handling NGS data collected from "any" amplicon-based schema**, provided that users fit the reference files to their amplicon design and data.

The default reference database of INSaFLU includes reference sequences of:

i) post-pandemic (2009) vaccine/reference influenza A(H1N1)pdm2009, A(H3N2) and B viruses (from both Northern and Southern hemispheres); 

ii) representative virus of multiple combinations of influenza HA/NA subtypes (i.e., H1N1, H2N2, H5N1, H7N9, etc)

iii) SARS-CoV-2 reference (e.g., SARS_CoV_2_Wuhan_Hu_1_MN908947.fasta) 

etc

All reference sequences at INSaFLU  are publicly available at NCBI (or are made available under permission of the authors). Download the current list here: :download:`INSaFLU_current_REFERENCE_DATABASE_11_01_2023.xlsx <_static/INSaFLU_current_REFERENCE_DATABASE_11_01_2023.xlsx>`) 

For influenza, the reference files have been prepared to fit amplicon-based schemas capturing the whole CDS of the main eight genes of influenza virus (PB2, PB1, PA, HA, NP, NA, M and NS).

During development, INSaFLU pipeline has been tested with NGS data collected after applying the wet-lab pre-NGS protocol for influenza whole genome amplification adapted from a RT-PCR assay described by Zhou and colleagues (Zhou et al, 2009, for Influenza A; and Zhou et al, 2014, for Influenza B; Zhou and Wentworth, 2012). This protocol can be applied to simultaneously amplify the eight genomic RNA segments, irrespective of influenza virus subtype or lineage.

You can download the suggested protocol here: :download:`Suggested_RT_PCR_assay_for_influenza_WGS.pdf <_static/Suggested_RT_PCR_assay_for_influenza_WGS.pdf>`

How to design a NGS run for influenza?
---------------------------------------

According to our tests during INSaFLU development, we suggest you ask your NGS service provider to perform runs in order to yield a final output of about 300000 (2 x 150000) reads per sample, if you use the influenza RT-PCR protocol indicated above.

This will account for issues arising from both the PCR reactions (e.g., fluctuations in the percentage of influenza-specific amplicons across samples and unbalanced relative proportions of the in-sample amplicons) and the NGS run (e.g., low yield and unbalanced demultiplexing of the reads across the samples).

This approach will allow you to end-up with more than 150000 (2 x 75000) reads per sample. This cut-off yielded a success (i.e., sample with 100% of the length of the 8 influenza CDS covered by ≥ 10-fold) of 92% on our pilot study using 2 x 150 paired-end reads (300 cycles). 

.. note::
   Examples of Illumina MiSeq runs that fit this suggestion are:
   
   i) run 96 samples using Illumina V2 Standard flow cells (30 M reads total; 300 cycles); 
   
   ii) run 24 samples using Illumina Micro flow cells (4 M reads total; 300 cycles).


References:

- Zhou B, Donnelly ME, Scholes DT, St George K, Hatta M, Kawaoka Y, Wentworth DE. 2009. Single-reaction genomic amplification accelerates sequencing and vaccine production for classical and Swine origin human influenza a viruses. J Virol, 83:10309-13.

- Zhou B, Lin X, Wang W, Halpin RA, Bera J, Stockwell TB, Barr IG, Wentworth DE.  2014. Universal influenza B virus genomic amplification facilitates sequencing, diagnostics, and reverse genetics. J Clin Microbiol, 52:1330-1337. 

- Zhou B, Wentworth DE. 2012. Influenza A virus molecular virology techniques. Methods Mol Biol, 865:175-92.
