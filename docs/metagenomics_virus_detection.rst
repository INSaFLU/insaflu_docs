Metagenomics virus detection
=================================

The TELEVIR  bioinformatics component of INSaFLU is a modular pipeline for the identification of viral sequences in metagenomic data (both Illumina and ONT data). 

It is composed of these main steps (detailed in https://insaflu.readthedocs.io/en/latest/bioinformatics-pipeline.html#metagenomics-virus-detection):

1. Read quality analysis and improvement [optional]
2. Viral Enrichment [optional].
3. Host Depletion [optional].
4. Assembly of the reads [optional].
5. Identification of the viral sequences.
	- Using reads.
	- Using contigs (if assembled).
	- Using several reference databases.
6. Selection of viral TAXID and representative genome sequences for confirmatory re-mapping
7. Remapping of the viral sequences against selected reference genome sequences. 

The pipeline culminates in the production of a set of summary statistics and visualizations of the results.
 
Below, you can find instructions on how to create a TELEVIR project, run samples and visualize/intrepret the results.

**TELEVIR Projects** - How to create and scale-up a metagenomics virus detection project
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
Within the *TELEVIR Projects* menu:

1. Go to *Projects* menu and choose *Create project*
....................................................

For enhanced sample comparison and output interpretation, users are encouraged to create a TELEVIR project per metagenomics sequencing run. 


2. Add a *Name* and *Description*, Choose the *Sequencing technology* and *Save*
................................................................................

3. Select the *Workflow and Software* to be run
................................................

After creating a project, and before adding/running sample, you can clicking in the "Magic wand" to select the bioinformatics workflow to be applied  to every sample added to the project. 

.. note::
   - As there is no “one-size-fits-all” bioinformatics pipeline that can detect all viruses, the **TELEVIR module was designed to allow users to easily run  complex workflows simultaneously** (covering several combinations of classification algorithms, databases and parameters, etc). The default workflows are “good-performant” workflows (selected after multiple testing and benchmarking) that together can potentiate the detection of clinical relevant viruses.
   - Nonetheless, users can always **turn ON/OFF specific steps (such as, Viral enrichment or Host depletion)** and and **select/deselect alternative software for each pipeline step**. Details about the current pipeline can be found here: https://insaflu.readthedocs.io/en/latest/bioinformatics-pipeline.html#metagenomics-virus-detection):

4. *Add samples* to the TELEVIR project
.........................................

Before running the added samples, please make sure that you have previously selected the bioinformatics workflows ("runs") to be applied  to every sample added to the project using the "Magic wand"  (see previous step)

5. *Run samples* by clicking in *Deploy Pathogen identification*
................................................................

At this time, users may start monitoring the Project progress by checking the runs (i.e., combinations of workflows) "Queued" or "Running". 

By clicking in *Run panel*, users can get an overviwe of the workflows run.


Output Visualization and Download
++++++++++++++++++++++++++++++++++

The INSaFLU-TELEVIR bioinformatics pipeline for metagenomics virus diagnostic generates multiple outputs, reflecting the multiple steps of the pipeline (detailed here: https://insaflu.readthedocs.io/en/latest/bioinformatics-pipeline.html#metagenomics-virus-detection).

The outputs are organized in dynamic 'expand-and-collapse' panels:


Main report - **Pathogen identification**
-----------------------------------------
   
This tab displays an interactive table with **summary statistics and visualizations of the end-point results of the TELEVIR metagenomic virus detection pipeline**. In summary, through this pipeline, reads and contigs (if available) are classified independently, then viral hits (TAXID) detected in both sides (and/or within the top-15 from each side) are selected for reference-based mapping against viral genome sequences present in the available databases. **This main report (interactive table) only includes viral hits that were classified at reads and/or contig level AND that had mapped reads.** 

.. note::
   - Only viral hits (TAXIDs and representative accession numbers) with mapped reads are automatically shown in this Main Report table.  
   - Other viral TAXIDs that were not automatically selected for confirmatory re-mapping step (flagged as "Unmapped") can be user-selected for mapping at any time by clicking in the "eye" icon available in the **Raw Classification and Mapping Summary** panel.
   
Below, you can find a description of the main outputs and statistics.

Mapping statistics
...................

- **Cov (%)**: horizontal coverage (i.e., percentage of the reference sequence covered)
- **Depth**: mean depth of coverage throughout the whole genome
- **DepthC**: mean depth of coverage exclusively in the covered regions
- **Mapped reads**: number of mapped reads
- **start prop (%)**:  mapped reads divided by the number of input reads (after QC)
- **mapped_prop (%)**: mapped reads divided by the number of reads used for mapping (i.e., reads retained after the "Virus enrichment" and/or "host depletion steps)
- **Gaps**: number of regions with 0 coverage
- **Windows Covered**: proportion of windows with mapped reads (the reference sequence is divided by x windows. 
- **class. success**: indication whether the TAXID was selected for mapping after reads and/or contigs classification
- **mapping success**: indication whether reads/and contigs successfully mapped against the TAXID representative references sequence
- **Warning**: 
	- *"Likely False Positive"*: when most reads map in a very small region of the reference sequence, i.e., hits with high “DepthC" but low “Depth” and low "Cov (%)". Flagged for hits with DepthC / Depth > 10 and Cov (%) > 5%.
	- *"Vestigial Mapping"*: when only a vestigial amount of reads (<= 2) mapped.
	

Mapping plots and output files
..............................

By clicking in a TAXID description, user can visualize/download multiple outputs regarding:

# READS MAPPING

- **Mapping Coverage** plot (depth of coverage thorughout the reference genome)

- **Integrative Genomics Viewer (IGV)** (interactive visualization of the mapped reads)

- **Mapped reads** in FASTA and BAM

- **Reference sequence** (".fa" format) and ".fai" index

# CONTIGS MAPPING

- **Assembly to reference dotplot** (location of the mapped contigs into the reference sequence)

- **Mapped contigs** in FASTA

- **Contigs alignment** in Pairwise mApping Format (PAF)

-**Sample remap**: statistics for reads'mapping against the set of contigs classified for a each TAXID.


Tips for Results interpretation
................................

XXXXXXXXXXXXXXXXX




Intermediate outputs 
---------------------

Multiple intermediate outputs and statistics are available by clicking in the following 'expand-and-collapse' panels:

Pre-processing: **Viral Enrichment** and/or **Host depletion**
..............................................................

This tab provides an overview on the number of reads filtered during the **Viral enrichment** and/or **Host depletion** steps of the metagenomics virus detection pipeline.

- "Viral enrichment" - retains potential viral reads based on a rapid and permissive classification of the reads against a viral sequence database.
- "Host Depletion" - remove potential host reads based on reference-based mapping against host genome sequence(s) 

The reads retained are provided for download (fastq.gz format).


**Assembly**
............

This tab provides an overview on **de novo** assembly step (this steps uses the reads retained after the "Viral enrichment" and/or "Host depletion" steps).

Filtered contigs are provided for download (fasta.gz format).

**Reads and Contigs classification**
....................................

This tab provides **reads and/or contigs classification reports** (tsv format) with the list of viral hits (TAXID and representative accession numbers) detected after the intermediate screening  against viral sequence databases. The two reports are merged to select the most "suspected" viral hits to be automatically selected for confirmatory re-mapping (see next steps). These reports are also compiled in the **Raw Classification and Mapping Summary** panel (below).


**Remapping** of the viral sequences against selected reference genome sequences. 

.................................................................................

Reads (and contigs) are mapped against  representative genome sequences of the "suspected" viral hits identified in the previous step.

This tab provides an overview on the amount of viral hits (TAXIDs and representative accession numbers) yielding mapped reads/contigs. Only viral hits with mapped reads are shown in the Main Report - Pathogen Identification.  



**Raw Classification and Mapping Summary**

.................................................................................

This table lists all viral hits (TAXID and representative accession numbers) detected during the intermediate step of **Reads and Contigs Classification** (see above), indicating if they were (or not) automatically selected for confirmatory re-mapping.

TAXIDs that were not automatically selected for confirmatory re-mapping step (flagged as "Unmapped") can be user-selected for mapping at any time by clicking in the "eye" icon.


