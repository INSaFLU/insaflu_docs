Change log
==========

This tab includes a list (chronologically ordered) of notable changes in INSaFLU.

November 28, 2024
-------

New update for the reference-based genome assembly projects:

- Small modification to the snippy used by INSaFLU to enable variants being called in duplicated regions, such as the monkeypox inverted terminal repeats. Users wanting to enable variant calling in duplicated regions using snippy should set to **0** the `--mapqual` parameter in the "mutation detection and consensus generation" section.

New updates regarding monkeypox virus:

- Added new monkeypox virus references, such as the clade Ia Zaire_1979-005 (DQ011155.1)

INSaFLU-TELEVIR updates:

- Downloadable mapping files in TELEVIR Reference Focus projects.

- Add Windows Covered statistic to Reference Focus Mapping; Add Reference Accid and Description to downloaded TSV.

- Improved reporting of running reference mapping processes.


October 18, 2024
.............

New updates regarding monkeypox virus:

- INSaFLU rapid sample classification can now distinguish the different Monkeypox clades (Clade Ia, Ib, IIa and IIb). 

- The assign2contigs database was updated with representative sequences of the Monkeypox Clade Ia, Ib, IIa and IIb.

- Added Nextclade link to Monkeypox Clade I;

Televir:

- fix TELEVIR reference panels load bug

- filter panel suggestions of existing accids


August 9, 2024
.............

**INSaFLU-TELEVIR update - version 2.1.0**

# TELEVIR module #
##################

The new update aims at facilitating confirmatory mapping against specific references, beyond the hypothesis-free and automatic confirmatory mapping of classic TELEVIR workflows.

For this purpose, Update 2.1.0 introduces several new features:

**Global References Management**

INSaFLU-TELEVIR "Reference" menu now has two sections: 

- INSaFLU (for reference-based genome assembly) (kept unchanged)
- TELEVIR (for Virus detection) - **new** section to support targeted reference mapping within TELEVIR projects

The **new** TELEVIR Reference path makes available three pages:

	i. **Files**, where the user can upload references to INSaFLU-TELEVIR and view existing, global databases. Uploaded references become availble for use as targets in TELEVIR mapping workflows, individually, or as part of Reference Panels. 
	
	ii. **References**, where the user can search references in the TELEVIR databases (global and personally uploaded), and select them for transfer to the INSaFLU reference database for deeper investigation. 
	
	iii. **Panel**. In this section the user is able to build mapping panels of existing (global or personal) TELEVIR references for mapping in TELEVIR projects.

**TELEVIR References Management**

A **new** column "References" was added to the sample table of TELEVIR projects, linking to the Sample Reference Management page. This page is composed of four sections:

	i. “Workflow” : a graphical display of classification workflows deployed on this sample.
	ii. “Added References”: displays references that were manually added to a sample.
	iii. “Added Panels”: section displays reference panels that were attached to this sample.
	iv. **Actions** sidebar and **References table**.

The **References table** displays Viral hits (description, taxid, accid) attached to this sample, either as a result of a classification workflow, or through manual addition and mapping. Additionally, the References table comports three summary columns:

	- “Runs”: indicates the Classification workflows that reference was identified in.
	- **“Global Ranking”**: sorts references by applying the TELEVIR sorting algorithm to the aggregate number of read and contig hits obtained by classification software across Classification Workflows.
	- **“Ensemble Ranking”**: sorts references according to their average TELEVIR sorting rankings across workflows they were identified in.
	- “Best Mapping”: indicates whether a reference was mapped against, and whether mapping was successful, providing a link to the workflow with the best mapping statistics.

The **Actions** sidebar enables to deploy mapping workflows against Added References, Added Panels, Selected References or Combined.

	-  **Map Added**: Map sample against all manually added references as well as references selected in the select column of the main reference table.
	- **Map Panels**: Map sample against all references within each panel added to the respective Reference Management Panel Section. Note: Panel mapping workflows deployed as part of Map Panels request are deployed independently.
	- **Map Combined**: maps the sample against the top reference hits as sorted across workflows using the Ensemble Ranking described above. The number of top references is defined in the Settings under the heading Global, and is the same as for individual workflows.
	
	Note: All Mapping workflows use parameters: i) from Workflows for Pre-processing; and ii) from Validation for Request Mapping and Filtering (see below).


**TELEVIR-Focus**

This update introduces the concept **“TELEVIR-Focus”**: a reference specific project to monitor and **deploy Mapping workflows against a single reference across several samples**.

Through the **Actions** button in the TELEVIR Project Page, the user is able to select a reference from among references attached to any project sample (either manually or as a result of a Classification Workflow) to create TELEVIR-Focus Projects, which are displayed in the Project page, below the samples table. Samples Selected through the Sample Select column will automatically be added to that project.

Within the Focus page, the user is able to add mapping workflows available from combinations of currently active software in the projects (or, if those are not set, global) settings. Mapping workflows are displayed separately and permit the display of mapping statistics against the Focus reference across Samples added to the project. If not yet deployed the user can choose to deploy new mappings. After completion, stacked IGV displays are availble, as well as the possibility to generate stacked a VCF and variant-specific igv-reports.

Finally, the **TELEVIR-Focus project also bears an INSaFLU connect button**: the reference in focus is automatically transfered to the INSaFLU reference database and to an INSaFLU project using that reference and including all samples in the Focus project. The new INSaFLU project and its status are displayed within the TELEVIR-Focus project.

**TELEVIR Settings and New Software**

INSaFLU-TELEVIR Update 2.1.0 introduces a new settings configuration, whereby settings are now grouped as "Workflows", "Validation" and "Global". This expanded configuration reflects the focus of this update on confirmatory mapping.

	- **"Workflow"** settings: These parameters control the combinations of workflows deployed classicaly, and require at least one classifier to be turned ON. However, one important development relative to version 2.0.0 is that software in the Remapping step can now be turned OFF. In practice, this will run the workflows as usual, but produce no detailed report of mapping statistics against top hits at the end of the run. However, raw hits (classically displayed beneath the mapping report section), are still collected, and will appear in the Sample Reference Management page as unmapped. This allows the user to delay confirmatory mapping until evidence from several workflows has been gathered, at which point one of the Validation Workflows can be deployed.
	
	- **"Validation"** Workflows: These settings control the Request Mapping and Map filtering deployed specifically as part of "Mapping Only" requests: Map Combined, Map Added and Map Panels, described above. Of Note, Mapping requests incorporate any of the three pre-processing steps in the Workflow section (Extra Filtering, Enrichment and Host Depletion), and will deploy active software in those sections. This will result in more than one workflow being deployed if more than one software is active.
	
	- **"Global"** settings: This single step section controls Final Reporting and Remap Management across all workflows (flag type, overlap threshold for report grouping, max taxid and accids for inclusion in single workflow remapping and Combined Mapping requests).

	Note: This flexibility allow running Classic workflows (including Mapping), Classification-only (with the mapping being deployed later on taking advantage of the new "Combined Deployment" (see below)

Other updates in TELEVIR Settings:

	- Added a new “Remapping filtering” option called “dustmasker - low complexity filtering”. This option  will mask low-complexity regions (e.g., homopolymeric tracts) in the references in order to reduce false positives caused by cross-mapping exclusively in these regions.
	- Remapping can be turned OFF (specially useful when using “Map Combined”, which reduces analysis time by avoiding the repetition of mappings against references that have been identified in multiple individual workflows)


New software:

- **Host depletion**
	- Added **Bowtie2 for Illumina**

- **Viral Enrichment**
	- Added **Kraken2 for Illumina**

- **Read / Contig Classification**
	- Added **Diamond** for Read CLassification (**Illumina & ONT**)
	- Added **Kraken2 for ONT Read Classification**
	- Added **Kraken2 for Contig Classification (Illumina & ONT)**

- **Remapping**
	- Added **Bowtie2 for Illumina**

- **Remap filtering**
	- Added **“Dustmasker*** for filtering  low complexity regions in References

- **Remapping - Management**
	- Default number of Accids to map lowered to 4, applied to new accounts.

**TELEVIR Combined Deployment**

Classic Classification Workflows receive a new deployment architecture in Update 2.1.0: Tree Deployment.

In INSaFLU-TELEVIR 2.0.0, the presence of multiple active software in any single parameter section resulted in the deployment of multiple workflows, corresponding to every possible combination from the available set of software / pipeline steps. This feature remains one of the important developments of INSaFLU-TELEVIR in terms of promoting cross-validation and robust identification. However, different workflow combinations were deployed independently, possibly resulting in the repeat computation of redundant pipeline steps. Update 2.1.0 introduces a deployment architecture that branches in line with the configuration steps. Pipeline steps are then deployed sequentially, by branch. The end results are a faster overall runtime and a reduction in storage requirements.


**TELEVIR - Reporting**

	- Interative workflow diagrams are now coloured according to the step
	- Simplified display of reports with collapsed reporting groups. Within group hits are sorted by "Cov (%), with the top hit always shown (group secondary hits hide and toggle - left indicator row)
	- Interactive Heatmaps for Cross-mapping inspection are provided, both across all groups (“Read Overlap Summary”) and within groups (”Reads Overlap)
	- Added 2 new columns to the Report: Private reads and Mismatch rate 
	- Sample workflow page separates "Classic workflows" (with classification) from Mapping workflows (upon request).




# Other changes #
##################

	- When uploading a sample, you can now specify its technology (Illumina or ONT). This can be done when uploading a single sample, or when uploading in batch by adding an optional column 'technology' in the input metadata file. This is reflected in the example input tsv metadata file.  
	- For single-end reads that for some reason fail the preprocessing (e.g short ONT reads that are wrongly automatically set as being Illumina), you can swap the technology (it will rerun the preprocessing for the new technology)
	- In the samples page, we added a button that allows a user to batch delete all samples not associated to projects.
	- The clades for the Monkeypox nextclade build were updated to include the C.1.1 clade



June 21, 2024
..........................

Representative sequences of the **ongoing A/H5N1 cattle outbreak** were included in the defalut Reference database:

- **A_H5N1_A_cattle_Texas_24_008749_002_2024** (downloaded from GenBank: **A/cattle/Texas/24-008749-002/2024(H5N1)** https://www.ncbi.nlm.nih.gov/nuccore/?term=A%2Fcattle%2FTexas%2F24-008749-002%2F2024(H5N1) ). It corresponds to same reference used for genome assembly by https://github.com/andersen-lab/avian-influenza.

- **A_H5N1_A_cattle_Texas_56283_2024** (downloaded from GenBank: **A/cattle/Texas/56283/2024(H5N1)** https://www.ncbi.nlm.nih.gov/nuccore/?term=A%2Fcattle%2FTexas%2F56283%2F2024(H5N1) ). This sequence was first described by Oguzie JU et al, Emerg Infect Dis. 2024 https://doi.org/10.3201/eid3007.240717

June 4, 2024
..........................

New updates regarding influenza A/H5Nx:

- Added Nextclade links for A-H5Nx-2.3.4.4, A-H5Nx-2.3.2.1 and A-H5Nx;

- New Nextstrain builds are available for A/H5N1:  PB1, PA, NP, MP and NS segments. So, the 8 segments are now covered

ABRIcate rapid identification and/or (sub)typing is now also performed on the consensus sequences obtained in the reference-based projects (to refine the rapid classification obtained from draft contigs just after reads upload). In the case of SARS-CoV-2 projects, the pangolin-based lineage is displayed in the Classification column instead.

Clade information was updated for the SARS-CoV-2 nextstrain build to include the 24A and 24B clades.

Primer cleaning for ONT samples in reference-based projects was refined to avoid excessive read filtering, particularly in samples with higher diversity relative to the chosen reference.


April 25, 2024
..........................

**The upgraded INSaFLU-TELEVIR is now published at Genome Medicine** https://doi.org/10.1186/s13073-024-01334-3. This article describes the extensions of INSaFLU since its first release in 2018, highlighting the development and implementation of a new module for metagenomic virus detection (TELEVIR), the incorporatioon of Nextstrain, the release of findONTime (https://github.com/INSaFLU/findONTime) and algn2pheno (https://github.com/insapathogenomics/algn2pheno), among other multiple features.

If you use INSaFLU-TELEVIR, please cite:
- Santos, J. D., Sobral, D., Pinheiro, M., Isidro, J., Bogaardt, C., Pinto, M., Eusébio, R., Santos, A., Mamede, R., Horton, D. L., Gomes, J. P., TELEVIR Consortium, & Borges, V. (2024). INSaFLU-TELEVIR: an open web-based bioinformatics suite for viral metagenomic detection and routine genomic surveillance. Genome medicine, 16(1), 61. https://doi.org/10.1186/s13073-024-01334-3


2023
-------

December 13, 2023
..........................

**RSV specific features:** 

- The **RSV Nextstrain builds** were updated to follow the more recent Nextstrain implementation, so that the new lineage classification [https://github.com/rsv-lineages] is automatically shown in the interactive trees.
- The existing direct links for **rapid RSV classification** of consensus sequences already offer the new genotype nomenclature implemented by **NextClade**.
- The assign2contigs database was updated with representative sequences of the RSV A [https://github.com/rsv-lineages/lineage-designation-A] and B lineages [https://github.com/rsv-lineages/lineage-designation-B] to facilitate the **identification of closely related references sequences and improve the user selection of appropriate reference sequences** for reads mapping.


October 20, 2023
..........................

- **TELEVIR Projects (virus detection):**

New host/vector sequences available for HOST DEPLETION:

**Host/vector name**  | **Common name** | **sequence** 

- aedes_albopictus	| **mosquito** |	GCF_006496715.2_Aalbo_primary.1_genomic.fna.gz
- anas_platyrhynchos	| **duck** |	GCF_015476345.1_ZJU1.0_genomic.fna.gz
- bos_taurus	| **cow** |	GCF_002263795.3_ARS-UCD2.0_genomic.fna.gz
- canis_lupus_familiaris	| **dog** |	GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.fna.gz
- culex_pipiens	| **mosquito** |	GCF_016801865.2_TS_CPP_V2_genomic.fna.gz
- cyprinus_carpio	| **carp** |	GCF_018340385.1_ASM1834038v1_genomic.fna.gz
- felis_catus	| **cat** |	GCF_018350175.1_F.catus_Fca126_mat1.0_genomic.fna.gz
- gallus_gallus	| **chicken** |	GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz
- marmota_marmota	| **marmot** |	GCF_001458135.2_marMar_genomic.fna.gz
- neogale_vison	| **mink** |	GCF_020171115.1_ASM_NN_V1_genomic.fna.gz
- oncorhynchus_mykiss	| **rainbow_trout** |	GCF_013265735.2_USDA_OmykA_1.1_genomic.fna.gz
- phlebotomus_papatasi	| **sandfly** |	GCF_000439695.1_Ppap_1.0_genomic.fna.gz
- pipistrellus_kuhlii	| **bat** |	GCF_024763615.1_Ppap_2.1_genomic.fna.gz
- salmo_salar	| **atlantic_salmon** |	GCF_905237065.1_Ssal_v3.1_genomic.fna.gz 
- sus_scrofa	| **pig** |	GCF_000003025.6_Sscrofa11.1_genomic.fna.gz



September 8, 2023
..........................

- **TELEVIR Projects (virus detection):**
	1. **Reports** are now generated per **Workflow** (as previously), per **Sample** (**NEW REPORT** combining non-redundant hits detected across workflows) and per **Project** (combining several samples, as previously), with a decreasing level of detail.
	2. **New button to “Sort sample reports”**. Viral hits (reference accession IDs) in the main reports (at both “Workflow” and “Sample” levels) can now be grouped and sorted by the degree of overlap of cross-mapped reads. This grouping intends to place together true positive hits with their corresponding cross-mapped potential false positives, allowing for the easy identification of the latter. It can be also useful to join same-segment references (for segmented virus) and to help identifying reference sequences most closely related to the virus present in the sample. The grouping parameter (--r-overlap) is modifiable in a new “Reporting” section of the TELEVIR Settings Menu for both technologies. “Sort sample report” should be deployed everytime the grouping parameter is changed for existing projects.
	3. **New step in the Workflow - “Extra filtering”**. Low complexity regions (e.g., homopolymeric tracts or repeat regions) are a common source of false-positive bioinformatics hits, as such we added an filtering layer that targets low complexity reads using the software PrinSeq++ (Cantu et al. 2019). This additional layer is optional and disabled by default.
	4. **New step in the Workflow - “Mapping stringency”**. An optional, extra layer of “mapping stringency” was added to this step to minimize false positive hits, allowing users to set a maximum sum of the mismatch qualities before marking a read unmapped and a maximum fraction of nucleotide mismatches allowed before soft clipping from ends (Using Bamutils). This additional layer is optional and disabled by default within the settings “Remapping” section.

- **Type/Subtype identification upon ONT reads upload**
	1. Screening is now performed over a draft assembly (using Raven) instead of directly from reads to increase precision. This new feature will be turned ON by default in new accounts.

June 16, 2023
..........................

- **Mutation detection and consensus generation**: We've added an extra parameter to enable primer removal using iVar (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1618-7), for both Illumina and ONT data. The procedure is an adaptation of the iVar CookBook (https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb). This additional layer is optional within the settings “Mutation detection and consensus generation” section, for both Illumina and ONT.


- **Nextstrain DATASETS:** 
	1. A new **Generic with Time Tree** build is now available. It is similar to the Generic build, but it also builds a time tree, inferring a mutation rate from the sample dates. Like in the Generic build, one reference is required to align the dataset consensus sequences. Nonetheless, unlike in the Generic build, the reference is not specifically defined as the root, but inferred from the data instead. To make use of this build, you need to accurately specify dates associated with each sample.
	2. In the specific case of the **SARS-CoV-2 build**, when importing consensus from projects, the reference of the project is no longer included automatically in the dataset. For the other builds, the project reference is still automatically included.
	

- **Algn2pheno module**: We introduced a new database of Spike amino acid mutations in epitope residues listed in Carabelli et al, 2023, 21(3), 162–177, Nat Rev Microbiol (https://doi.org/10.1038/s41579-022-00841-7), Figure 1. This is now the report that is visualized in the project page. **Important**: for older projects, the visualization of the Algn2pheno report will fail. Nonetheless, the old reports are still available in the algn2pheno.zip file that is downloadable from the project page. To update to the new database, you can either modify the settings of at least one sample within a project, add/remove samples in the project, or create a new project with the same samples.


- **Type/Subtype identification**: We relaxed the mincov parameter in abricate from 60% to 40%. 


- **Other:** We performed internal modifications to improve the stability of the website. We also performed minor aesthetical adjustments. 


May 8, 2023
..........................

- **Nextstrain DATASETS:** new builds for the **avian influenza (A/H5N1)** are now available (HA, NA and PB2 genes), allowing phylogenetic and spatiotemporal analysis using Nextstrain workflow https://github.com/INSaFLU/nextstrain_builds/tree/main/avian-flu (adapted from https://github.com/nextstrain/avian-flu). From now on, the build is selected upon creation of a New Dataset [cannot be changed afterwards].

- **References menu:**
	1. **Vaccine-like reference sequences for the 2023-2024 season publicly available at GenBank, for A/H3N2 (A/Darwin/6/2021) and A/H1N1 (A/Wisconsin/67/2022)**, are now available in INSaFLU reference default database. This update was performed with kind support of the WHOCC Reference and Research on influenza, VIDRL, Melbourne, Australia (special thanks to Dr. Ammar Aziz and Dr. Ian Barr);
	2. All seasonal influenza sequences (A/H3N2, A/H1N1, B/Victoria and B/Yamagata) available at the default INSaFLU database were re-annotated to allow mutation annotation following the **HA1 numbering** (i.e., mutations will now be annotated for each peptide: signal peptide, HA1 and HA2 peptides, instead of the full-protein). **The new annotation (HA1 numbering) will only be applied to new projects**

- **TELEVIR Projects (virus detection):**
	1. **Controls:** user can now select “control” sample(s) within a TELEVIR project. Viral TAXID detected in the Main report of the user-selected “control” sample(s) will be flagged in the reports of samples in the same project as “Taxid found in control” in a new “Control” column. **This new functionality is designed to facilitate the background subtraction of negative controls.** Multiple controls are possible.
	2. Added a **new button to start analyses of particular samples** within a TELEVIR project. 
	3. New search tab in TELEVIR projects. Relies on Project and Sample names.  
  

- **Release of findONTime** (https://github.com/INSaFLU/findONTime)
	1. **Description:** This tool **runs concurrently with MinION sequencing** and merges (at user defined time intervals) the FASTQ files that are being generated in real-time for each sample. It can also automatically upload the files to a local docker instance of the INSaFLU-TELEVIR platform and launch the metagenomics virus detection analysis using the TELEVIR module. 
	2. **Motivation and Goal:** This development will allow users **to detect a virus in a sample as early as possible during the sequencing run**, reducing the time gap between obtaining the sample and the diagnosis, and also reducing sequencing costs (as ONT runs can be stopped at any time and the flow cells can be cleaned and reused). 
	3. **Usage:** findONTime can be used as a “start-to-end” solution or for particular tasks (e.g., merging ONT output files, metadata preparation and upload to INSaFLU-TELEVIR). See examples here: https://github.com/INSaFLU/findONTime#usage 
	

- **Local DOCKER installation:**  The new docker installation version 2.0.0 (including the TELEVIR module) is now available at https://github.com/INSaFLU/docker. To avoid incompatibilities when updating the previous local installations, **we recommend that users set up a brand new installation.** 


March 7, 2023
..........................

- **Respiratory Syncytial Virus (RSV) analysis**
	- Added multiple reference sequences (dispersed accross the RSV phylogeny) to the Default Reference Database (https://insaflu.readthedocs.io/en/latest/uploading_data.html#uploading-reference-data)
	- Added multiple RSV sequences to the assign2contigs database as a mean to faciliate the selection of closely related references for mapping.

- **SARS-CoV-2 clade/ lineage classification**:  
	- Upgraded the PANGO version; *usher* mode is now the default (instead of pangolearn).
	- Update clades of the SARS-CoV-2 nextstrain build.

- **TELEVIR** Projects:
	- Renamed the "Deploy Pathogen Identification" button to "Run".
	- Updated the coverage graphics components: coverage plots now using weighed average.
	- Corrected bug in the generation of outputs after mapping by request in the "Raw Classification and Mapping Summary"
	

For more information, please consult:

- Documentation : https://insaflu.readthedocs.io/en/latest/

- Github page: https://github.com/INSaFLU



February 2, 2023
..........................

Bug fix:

**Algn2pheno module**: solve bug in mutation count for sequences with no mutations (default 0); fix final report phenotype categories to display sets of flagged mutations instead of single draw. Update algn2pheno package to 1.1.5



January 26, 2023
..........................

**Important update:**

**New features for Respiratory Syncytial Virus (RSV) analysis**:

	- INSaFLU PROJECTS (reference-based mapping): **direct links for rapid  RSV clade/genotype classification using Nextclade (https://clades.nextstrain.org/)** are now automatically provided for RSV projects. The reference sequences used in NextClade for RSV-A (hRSV/A/England/397/2017) and RSV-B (hRSV/B/Australia/VIC-RCH056/2019) were also made available in the References database, with kind permission of the sequence authors/owners (UKHSA and WHO CCRI, respectively).
	- Nextstrain DATASETS: **two new builds (RSV_A and RSV_B) are available**, allowing RSV-specific phylogenetic and spatiotemporal analysis using Nextstrain workflow https://github.com/nextstrain/rsv. 
	- Samples menu: **RSV-A / RSV-B** was included in the typing database for **rapid classification** just after reads upload. 

Other changes:

- **TELEVIR** Projects:
	- The Run table report column (TELEVIR Projects > Project > Sample) is now dynamically updated to represent the current status of an ongoing run, by module.
	- Refinements in the Reference mapping optimization to prevent memory overflow crash in large samples.
	- **Centrifuge software was added to the Illumina Read Classification Panel**. To activate this feature, the **user must visit the mains Settings page**. For existing projects with project settings these must be reset.
	
- **Nextstrain influenza**: to allow more sequences to be inserted in the tree, we've slightly alleviated the inclusion criteria allowing more NNN and divergence in the consensus sequences (25 ambiguous positions are allowed in the HA protein and clock_filter_iqd increased to 12)

For more information, please consult:

	- Documentation : https://insaflu.readthedocs.io/en/latest/

	- Github page: https://github.com/INSaFLU


2022
-------

December 21, 2022
..........................

**Major update:**

A **New module for metagenomics virus detection (called TELEVIR)** has been released.  The main features of the TELEVIR module are:

	- handles both Illumina and ONT data;

	- allows easily running complex modular workflows, covering several combinations of steps (e.g., with/without Viral enrichment/Host depletion), classification software (e.g., Kaiju, Krakenuniq, Kraken2, Centrifuge, FastViromeExplorer), databases (NCBI RefSeq viral genome, Virosaurus, etc) and parameters;

	- includes automate “confirmatory” re-mapping against reference viral genome(s) present in the available databases;

	- culminates in user- and diagnosis-oriented  reports, including (interactive) tables and  graphs (e.g., coverage plots, Integrative Genomics Viewer visualization, Assembly to reference dotplot), as well as multiple downloadable output files (e.g., list of the software/parameters, reads/contigs classification reports, mapped reads/contigs identified per each virus; reference sequences, etc)
 

For more information about this new module (features, functionality, etc), please consult:

	- Tutorial and outputs: https://insaflu.readthedocs.io/en/latest/metagenomics_virus_detection.html#metagenomics-virus-detection

	- Pipeline details: https://insaflu.readthedocs.io/en/latest/bioinformatics_pipeline.html#metagenomics-virus-detection

	- INSaFLU Github page: https://github.com/INSaFLU

 

October 27, 2022
..........................

**Important update:**

- **New module (called “Datasets”) for Nextstrain (https://nextstrain.org/) phylogenetic and geotemporal analysis.** This user-friendly functionality will allow INSaFLU users to launch virus-specific Nextstrain builds (seasonal Influenza, SARS-CoV-2 and Monkeypox) as well as a “generic” build that can be used for other viruses.

See more details in INSaFLU documentation: https://insaflu.readthedocs.io/en/latest/data_analysis.html#nextstrain-datasets and https://insaflu.readthedocs.io/en/latest/output_visualization.html#navigate-through-nextstrain-datasets  and https://github.com/INSaFLU/nextstrain_builds


- **Integration of the “algn2pheno” (https://github.com/insapathogenomics/algn2pheno) tool within the “Projects” menu**. This new functionality screens SARS-CoV-2 Spike amino acid alignments in each SARS-CoV-2 project against two default “genotype-phenotype” databases: the COG-UK Antigenic mutations (https://sars2.cvr.gla.ac.uk/cog-uk/)  and the Pokay Database (https://github.com/nodrogluap/pokay/tree/master/data). **Align2pheno reports the full repertoire of Spike amino acid change found in each sequence, flagging for the presence of mutations of interest (and their potential impact on phenotype) included in those databases.**


See more details in INSaFLU documentation:  https://insaflu.readthedocs.io/en/latest/data_analysis.html#algn2pheno and https://insaflu.readthedocs.io/en/latest/output_visualization.html#h-explore-the-algn2pheno-report-panel-and-results

*Acknowledgements*

This important update was only possible with the contribution of several people and teams. **We would like to deeply acknowledge to:**

	- All INSaFLU developing team, with special thanks to Daniel Sobral (INSA), Miguel Pinheiro (Institute of Biomedicine - iBiMED, University of Aveiro), João Dourado Santos (INSA), Miguel Pinto (INSA), Joana Isidro (INSA) and Vítor Borges (INSA).
	- Carlijn Bogaart and Daniel Horton (University of Surrey, UK), for their key contribution to build the algn2pheno (https://github.com/insapathogenomics/algn2pheno) tool.
	- Nextstrain https://nextstrain.org/ team, for their amazing work in developing open-source tools for phylogenetic and geotemporal tracking of viral pathogens.
	- COK-UK consortium (https://www.cogconsortium.uk/) (UK) and the University of Calgary (Canada) for making available updated and comprehensive SARS-CoV-2 mutations databases (https://sars2.cvr.gla.ac.uk/cog-uk/ and https://github.com/nodrogluap/pokay/tree/master/data, respectively) for algn2pheno screening.
	- The Infraestrutura Nacional de Computação Distribuída (INCD) (https://www.incd.pt/)  for providing computational resources for testing the INSaFLU platform.
	- INSaFLU work has been supported by funding from the European Union’s Horizon 2020 Research and Innovation programme under grant agreement No 773830: One Health European Joint Programme, under the TELE-Vir project (https://onehealthejp.eu/jrp-tele-vir/) 


October 10, 2022
..........................

Users can now use trimmomatic to perform trimming of primer sequences of several predefined Primer pool sets:

– SARS-CoV-2 Primal Scheme V3 (https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv)

—SARS-CoV-2 Primal Scheme V4.1 (https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4.1)

– Monkeypox Primal Scheme from Welkers, Jonges and van den Ouden (https://www.protocols.io/view/monkeypox-virus-whole-genome-sequencing-using-comb-n2bvj6155lk5/v1)

—Monkeypox Primal Scheme from Chen et al. (https://www.protocols.io/view/monkeypox-virus-multiplexed-pcr-amplicon-sequencin-5qpvob1nbl4o/v2)

Please contact us if you want to add new Primer pools to the online tool


January 26, 2022
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

To update the local docker installation, please follow the instructions in https://github.com/INSaFLU/docker

	Note: After this update (i.e., INSaFLU versions **equal or higher 1.5.0**) users will be able to update their local installation to the latest version with a single command:
	```
	$ docker exec -it insaflu-server update-insaflu
	```


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
