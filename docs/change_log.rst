Change log
==========

This tab includes a list (chronologically ordered) of notable changes in INSaFLU.


2018
----


April 9, 2018 
.............

- Maximum size per fastq.gz file was upgraded from 50 MB to 300 MB. 

	* IMPORTANT NOTE: Files between 50 - 300 MB will be downsized to ~50 MB before analysis by randomly sampling reads.

- The draft assembly provided by INSaFLU (FASTA format) now also includes potential non-influenza specific contigs (i.e., contigs not assigned to any influenza segment / reference by INSaFLU). This feature    


March 9, 2018 
.............

- INSaFLU now provides a draft genome assembly (FASTA format) of the influenza virus genome. Influenza-specific NODES/contigs are assigned both to the corresponding viral segments and to a related reference influenza virus (output is a table in .tsv format). Check these new outputs and guide to interpret them at the INSaFLU tab "Samples / Extra info / Type and subtype/lineage identification"

	This new feature reinforces the application of INSaFLU to:
	
		* analyse viruses for which a close related whole-genome sequence is not available (e.g., avian influenza) at the INSaFLU or other databses (NCBI, GISAID, etc);
		* investigate reassortments
		* disclose mixed infections
	


January 25, 2018 
................

- INSaFLU 1.0.0 is released for the scientific community at https://insaflu.insa.pt 
	
	INSaFLU ("INSide the FLU") is an bioinformatics free web-based suite that deals with primary NGS data (reads) towards the automatic generation of the output data that are actually the core first-line “genetic requests” for effective and timely influenza laboratory surveillance. While INSaFLU has indeed some influenza-specific features (e.g., automatic type/subtype identification), there is no restrictions to use it for other viruses. 

	Highlights
	++++++++++

		* open to all, free of charge, user-restricted accounts
		* applicable to NGS data collected from any amplicon-based schema
		* allows advanced, multi-step software intensive analyses in a user-friendly manner without previous training in bioinformatics
		* automatic identification of influenza type and subtype/lineage, detection of putative mixed infections and intra-host minor variants
		* allows integrating data in a cumulative manner, thus fitting the analytical dynamics underlying the continuous epidemiological surveillance during flu epidemics
		* outputs are provided in nomenclature-stable and standardized formats and can be explored in situ or through multiple compatible downstream applications for fine-tune data analysis and visualization