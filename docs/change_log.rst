Change log
==========

This tab includes a list (chronologically ordered) of notable changes in INSaFLU.


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