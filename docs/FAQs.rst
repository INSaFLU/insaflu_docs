**Frequently Asked Questions (FAQs)**
======================================
FAQs
....

**UPLOAD (Input, Classification, Settings, etc)**
-------------------------------------------------

- **Which sequencing technologies are compatible with the platform?**
The platform supports second- and third-generation sequencing technologies, including Illumina (single-end and paired-end), Ion Torrent, and Oxford Nanopore Technologies (ONT).

- **What file types can be used as input?**
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

*Under construction*

**Nextstrain module** 
----------------------

*Under construction*

**TELEVIR module** 
----------------------

*Under construction*



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
