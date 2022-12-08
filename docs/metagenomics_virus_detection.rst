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

The INSaFLU-TELEVIR bioinformatics pipeline for metagenomics virus diagnostic generates multiple outputs, reflecting the multiple steps of the pipeline (detailed here: https://insaflu.readthedocs.io/en/latest/bioinformatics-pipeline.html#metagenomics-virus-detection):

The outputs are organized in dynamic 'expand-and-collapse' panels as follows:

- Main Table Report (**Pathogen Identification**)



- **Viral Enrichment** and/or **Host depletion** (Pre-processing)
- **Assembly**
- **Reads and Contigs classification** 	
- **Remapping** of the viral sequences against selected reference genome sequences. 





Pre-processing: **Viral Enrichment** and/or **Host depletion** (Pre-processing)
-------------------------------------------------------------------------------- 
   
This tab displays all information for all loaded samples.

A. XXXXX
..............................................................................




B. XXXX
...............................................................................

