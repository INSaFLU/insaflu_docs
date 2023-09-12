**Bioinformatics Pipeline**
============================

INSaFLU is an open web-based bioinformatics platform for metagenomic virus detection and routine genomic surveillance.

As such, it integrates two main pipelines / components:

- a **Reference-based genomic surveillance** pipeline (from NGS reads to quality control, mutations detection, consensus generation, virus classification, alignments, “genotype-phenotype” screening, phylogenetics, integrative phylogeographical and temporal analysis etc)
	
- a **Virus detection** pipeline (from NGS reads to quality control and metagenomics virus identification)

INSaFLU relies on a multi-software bioinformatics pipelines that will be under continuous development and improvement not only to enrich it with new features, but also to continuously shape the protocol to the best methodological advances in the field. 

The current software and default settings, which were chosen upon intensive testing, are described below, together with the list of Steps and Settings that can be turned ON/OFF or configured by the user, respectively. For additional details about the bioinformatics pipeline, please visit the INSaFLU github account: https://github.com/INSaFLU/INSaFLU (more information for each software can also be found in the official repositories; links are also provided below). 


Reference-based genomic surveillance
+++++++++++++++++++++++++++++++++++++

The “surveillance-oriented” bioinformatics component of INSaFLU allows launching:

**Projects** - From reads to reference-based generation of consensus sequences and mutations annotation, followed by gene- and genome-based alignments, amino acid alignments, Pango classification, NextClade linkage, etc.

and 

**Datasets** - From consensus sequences to advanced Nextstrain phylogenetic and genomic analysis, coupled with geographic and temporal data visualization and exploration of sequence metadata.

The bioinformatics pipeline currently consists of multiple steps (see WorkFlow) yielding multiple graphical, and sequence outputs (see *Output visualization and download* menu for details)

.. image:: _static/workflow_update_20221114.png

Simplified illustration of the main steps of the modular INSaFLU-TELEVIR bioinformatics pipeline for reference-based genomics surveillance.

*Find below the description of the main bioinformatics steps (Description and Current software )*

Read quality analysis and improvement
--------------------------------------


*Description*

This step takes the input single- or paired-end reads (fastq.gz format) and produces quality processed reads, as well as quality control reports for each file, before and after quality improvement. This module is automatically run upon reads upload (i.e., no user intervention is needed). 

*Software version/settings*

.. note::

	**## ILLUMINA / Ion Torrent data ##**
	
   	**FastQC** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.5; date 15.01.2018)

		input: single- or paired-end reads (fastq.gz format) (e.g., sample_L001_R1_001.fastq.gz and sample_L001_R2_001.fastq.gz for Illumina technology reads)
		
		--nogroup option: all reports will show data for every base in the read. 
		
	**Trimmomatic** (http://www.usadellab.org/cms/index.php?page=trimmomatic) (version 0.27; date 15.01.2018)
	
		input: single- or paired-end reads (fastq.gz format) (e.g., sample_L001_R1_001.fastq.gz and sample_L001_R2_001.fastq.gz for Illumina paired-end reads)
	
		ILLUMINACLIP: To clip adapters from the input file using user-specified adapter sequences. ILLUMINACLIP:<ADAPTER_FILE>:3:30:10:6:true
		
		HEADCROP: Cut the specified number of bases from the start of the read. Range: [0:100]. If value equal to 0 this parameter is excluded.
		
		CROP: Cut the read to a specified length. Range: [0:400]. If value equal to 0 this parameter is excluded.
	
		SLIDINGWINDOW: perform a sliding window trimming, cutting once the average quality within the window falls below a threshold (default: SLIDINGWINDOW:5:20, where 5 refers to window and 20 to the minimum average quality)
	
		LEADING: cut bases off the start of a read, if below a threshold quality (default: LEADING:3). This will allow discarding bases with very quality or N bases (quality score of 2 or less).
	
		TRAILING: cut bases off the end of a read, if below a threshold quality (default: TRAILING:3). This will allow discarding bases with very quality or N bases (quality score of 2 or less).
	
		MINLEN: drop the read if it is below a specified length (MINLEN:35)
	
		TOPHRED33:  Convert quality scores to Phred-33
		
	**## Oxford Nanopore Technologies (ONT) data ##**
		
	**NanoStat** (https://github.com/wdecoster/nanostat) (version 1.4.0)
		
		input: ONT reads (fastq.gz format) 

	**NanoFilt** (https://github.com/wdecoster/nanofilt) (version 2.6.0)
	

		**-q (QUALITY)**: Filter on a minimum average read quality score. Range: [5:30] (default: 10)
		
		**-l (LENGTH)**: Filter on a minimum read length. Range: [50:1000]. (default: 50)
		
		**--headcrop**: Trim n nucleotides from start of read. Range: [1:1000]. If value equal to 0 this parameter is excluded. (default: 70)
		
		**--tailcrop**: Trim n nucleotides from end of read. Range: [1:1000]. If value equal to 0 this parameter is excluded. (default: 70)
		
		**--maxlength**: Set a maximum read length. Range: [100:50000]. If value equal to 0 this parameter is excluded. (default: 0)
		

	**RabbitQC** (https://github.com/ZekunYin/RabbitQC)  (version 0.0.1)**
		
		input: ONT reads (fastq.gz format) pre- and post- quality improvement with NanoFilt
		
		Files between 50 - 300 MB are downsized to ~50 MB before analysis by randomly sampling reads using fastq-sample from fastq-tools package https://github.com/dcjones/fastq-tools (developed by Daniel C. Jones dcjones@cs.washington.edu)


.. note::

	**## ILLUMINA data only ##**
	
		***Users can also use trimmomatic to perform trimming of primer sequences of several predefined Primer pool sets:
		
			-- SARS-CoV-2 Primal Scheme V3 (https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv)
			
			-- SARS-CoV-2 Primal Scheme V4.1 (https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4.1)
			
			-- Monkeypox Primal Scheme from Welkers, Jonges and van den Ouden (https://www.protocols.io/view/monkeypox-virus-whole-genome-sequencing-using-comb-n2bvj6155lk5/v1)
			
			-- Monkeypox Primal Scheme from Chen et al. (https://www.protocols.io/view/monkeypox-virus-multiplexed-pcr-amplicon-sequencin-5qpvob1nbl4o/v2)
			
		Please contact us if you want to add new Primer pools to the online tool

.. important::
	INSaFLU allows users to configure key parameters for reads quality analysis in the tab **“Settings”**. 
	
	**Please check your settings before uploading new samples to your account.**
	
	See details in https://insaflu.readthedocs.io/en/latest/data_analysis.html#user-defined-parameters


Influenza type and sub-type identification (and identification of other viruses: Human Betacoronavirusm, RSV and MPXV)
-------------------------------------------------------------------------------------------------------------------------------------

*Description*
 
In this module, draft assemblies derived from post-QC reads are screened (using ABRIcate) against two INSaFLU in house sequence markers databases: 

i) "influenza_typing", which drives the discrimination of the influenza types A and B, all currently defined influenza A subtypes (18 hemagglutinin subtypes and 11 neuraminidase sub-types) and the two influenza B lineages (Yamagata and Victoria).

ii) "influenza_assign_segments2contigs", which allows the automatic assignment of the assembled contigs to both the corresponding viral segments and to a related reference influenza virus. 

The generated outputs (i.e., draft assemblies, the identified type and subtype/lineage and a table linking contigs to segments/references) are automatically provided upon reads upload (i.e., no user intervention is needed). INSaFLU flags samples as "putative mixed infections" if more than one type, HA or NA subtype or lineage is detected. In addition, specific alerts are generated if an incomplete type/subtype is assigned. 

**Since March 10, 2020, these two databases have been upgraded for rapid classification and/or contigs assignment of Human Betacoronavirus (BetaCoV) and other viruses.** Details about the rationale behind this classification and outputs can be found here: :download:`INSaFLU_current_genetic_markers_v11_after_03_03_2023.xlsx <_static/INSaFLU_current_genetic_markers_v11_after_07_03_2023.xlsx>`

Similarly to influenza classification, alerts are generated if, for instance, no BetaCoV virus is assigned or an incomplete human BetaCoV classification is obtained (for instance, due to the presence of a low number of human BetaCoV reads, etc)

*Software version/settings*

.. note::

**## ILLUMINA / Ion Torrent data ##**
	
	**SPAdes** (http://cab.spbu.ru/software/spades/) (version 3.11.1; date 15.01.2018)
   
   		--pe1-1 and --pe1.2 (for paired-end) or -s (for single-end data): define the input files, i.e, quality processed reads (e.g., sample_1P.fastq.gz and sample_2P.fastq.gz)
				
		--only-assembler: runs assembly module only and does not perform reads correction
		
				(contigs with k-mer coverage below '3' are discarded for subsequent ABRIcate analyses to avoid the classification of vestigial sequencer-derived contaminating sequences)

**## Oxford Nanopore Technologies (ONT) data ##**

	**Raven** (https://github.com/lbcb-sci/raven) (version 1.8.1; date 08.09.2018)


**Illumina and ONT**

	**ABRIcate** (https://github.com/tseemann/abricate) (version 0.8-dev; date 15.01.2018)
	
		# For type and subtype/lineage identification (and Human BetaCoV classification*):
	
		--db influeza_typing: the INSaFLU "influenza_tying" database includes a set of type- and sub-type/lineage-specific gene markers that ensure the discrimination of the influenza types A and B, all currently defined influenza A subtypes (18 hemagglutinin subtypes and 11 neuraminidase sub-types) and the two influenza B lineages (Yamagata and Victoria).
	
		--minid: minimum DNA %identity (--minid 70)
		
		--mincov: minimum DNA % coverage (--mincov 40, until 15/06/2023: --mincov 60)
		

		***As of March 10th, 2020**, samples can be classified as: 

		- "BetaCoV” if the draft assemblies contain an “M gene” with ≥70% identity and ≥40% coverage (until 15/06/2023: 60%) to one of the M (partial) gene marker sequences of the five representative Human BetaCoronavirus genomes in the database)
		
		- “SARS_CoV_2”, "SCoV2_potential_Omicron", “MERS_CoV”, “SARS_CoV”, “HCoV_HKU1” or “HCoV_OC43” if the draft assemblies contain a “S gene” with ≥70% Identity and ≥40% coverage (until 15/06/2023: 60%) coverage to one of the S (partial) gene marker sequences of the five representative Human BetaCoronavirus (the classification reflects the closest match among the five human BetaCoV listed above).

				
		# For segments/references assignment: 
		
		--db influeza_assign_segments2contigs: this database includes segment sequence markers of several seasonal human influenza [including: i) post-pandemic (2009) vaccine/reference influenza A(H1N1)pdm2009, A(H3N2) and B viruses; ii) representative viruses of specific genetic groups/lineages/clades, as defined by International Health Authorities for each season)], as well as of avian influenza from several HA/NA subtypes (i.e., H1N1, H2N2, H5N1, H7N9, etc)
	
		--minid: minimum DNA %identity (--minid 70)
		
		--mincov: minimum DNA % coverage (--mincov 30)
		
		**Draft assemblies (Illumina/Ion Torrent data or ONT data) are labeled with the closest match among the five human BetaCoV (see above) if they have ≥70% Identity and ≥30% coverage to one of the five BetaCoV full-genome sequences or partial S/M genes in the database.
		
		Important note: Since the "influeza_assign_segments2contigs" database is naturally not as exhaustive as other databases (such as, NCBI, Fludb or EpiFLU/GISAID), users may need to run the draft assemblies in these databases (or associated tools, such as BLAST) for some purposes (e.g., to detect/confirm reassortments or to infer the closest reference sequence of each segment / genome).
		


Latest list of genetic markers (version 11; 07.03.2023) can be downloaded here: 

:download:`INSaFLU_current_genetic_markers_v11_after_07_03_2023.xlsx <_static/INSaFLU_current_genetic_markers_v11_after_07_03_2023.xlsx>`
				
Previous database versions can be downloaded here:

version 10 (until 07.03.2023) :download:`INSaFLU_genetic_markers_v10_before_07_03_2023.xlsx <_static/INSaFLU_genetic_markers_v10_before_07_03_2023.xlsx>`

version 9 (until 26.01.2023) :download:`INSaFLU_genetic_markers_v9_before_26_01_2022.xlsx <_static/INSaFLU_genetic_markers_v9_before_26_01_2022.xlsx>`

version 8 (until 26.10.2022) :download:`INSaFLU_genetic_markers_v8_before_26_10_2022.xlsx <_static/INSaFLU_genetic_markers_v8_before_26_10_2022.xlsx>`

version 7 (until 11.12.2021) :download:`INSaFLU_genetic_markers_v7_before_11_12_2021.xlsx <_static/INSaFLU_genetic_markers_v7_before_11_12_2021.xlsx>`

version 6 (until 27.07.2021) :download:`INSaFLU_genetic_markers_v6_before_27_07_2021.xlsx <_static/INSaFLU_genetic_markers_v6_before_27_07_2021.xlsx>`

version 5 (until 10.03.2020) :download:`INSaFLU_genetic_markers_v5_before_10_03_2020.xlsx <_static/INSaFLU_genetic_markers_v5_before_10_03_2020.xlsx>`

version 4 (until 10.01.2020) :download:`INSaFLU_genetic_markers_v4_before_10_01_2020.xlsx <_static/INSaFLU_genetic_markers_v4_before_10_01_2020.xlsx>`

version 3 (until 02.01.2019) :download:`INSaFLU_genetic_markers_v3_before_02_01_2019.xlsx <_static/INSaFLU_genetic_markers_v3_before_02_01_2019.xlsx>`

version 2 (until 05.06.2018) :download:`INSaFLU_genetic_markers_v2_before_05_06_2018.xlsx <_static/INSaFLU_genetic_markers_v2_before_05_06_2018.xlsx>`

version 1 (until 14.05.2018) :download:`INSaFLU_genetic_markers_v1_before_14_05_2018.xlsx <_static/INSaFLU_genetic_markers_v1_before_14_05_2018.xlsx>`		

Variant detection and consensus generation
------------------------------------------

*Description*

This key module takes enables reference-based mapping, followed by SNP/indel calling and annotation and generation of consensus sequences (quality processed reads obtained through Trimmomatic analysis are used as input). Quality processed reads obtained through Trimmomatic (Illumina/IonTorrent data) NanoFilt (ONT data) are used as input. A reference sequence must be selected for each project (select one from INSaFLU default reference database or upload one of your choice).  Uploaded “.fasta” files are annotated upon submission and automatically become available at the user-restricted reference database. For influenza, each project should ideally include viruses from the same type and sub-type/lineage (this typing data is automatically determined upon reads submission to INSaFLU).

*Software version/settings*

.. note::

	**##REFERENCE ANNOTATION##**
	
	**Prokka** (https://github.com/tseemann/prokka) (version 1.12; date 15.01.2018)
   
		--kingdom: defines the Annotation mode (Viruses)
	
	
	**##ILLUMINA / Ion Torrent data##**
	
	**Snippy** (https://github.com/tseemann/snippy) (version 3.2-dev - sligthly modified (details in https://github.com/INSaFLU/INSaFLU); date 15.01.2018)
	
		--R1 (and --R2): define the reads files used as input, i.e, quality processed reads (e.g., sample_1P.fastq.gz and sample_2P.fastq.gz) obtained after Trimmomatic analysis
		
		--ref: define the reference sequence selected by the users (.fasta or gbk format) 
		
		--mapqual: minimum mapping quality to accept in variant calling(default: --mapqual 20) 
		
		--mincov: minimum coverage of variant site (default: --mincov 10)
		
		--minfrac: minimum proportion for variant evidence (default: --minfrac 0.51)
		
		--primer: defines primer sequences to be removed using iVar(version 1.4.2, available since 16/06/2023) (by default no primers are removed). The primer removal procedure was based on the iVar CookBook (https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb), but where no quality filtering is performed, and reads starting outside the primer are not excluded. Primer removal is obtained after the alignment step, but before variant calling and consensus generation.
		
		
	**## Oxford Nanopore Technologies (ONT) data ##**
	
	_Mapping:
	
	**Medaka** (https://nanoporetech.github.io/medaka/ (version 1.2.1)
		
		input: ONT quality processed reads obtained after NanoFilt analysis.
		
		medaka consensus -m model (default: r941_min_high_g360) --primer (default: empty)
			Optional primer removal using iVar follows the same procedure as described before for snippy, being applied before consensus generation.
					
		medaka variant
		
	_VCF filtering:
	
		Mutations are filtered out based on the following user-defined criteria:
		
			**Minimum depth of coverage  per site** (equivalent to --mincov in Illumina pipeline) (default: 30)
			
			**Minimum proportion  for variant evidence** (equivalent to --minfrac in Illumina pipeline) (default: 0.8)
			
			
	For each mutation, two COVERAGE values are provided in final table output: the depth of unambiguous reads spanning pos +-25 (as provided by medaka variant module) and depth per site as provided by samtools (depth -aa). Values are separated by “/”. 
	
	_Consensus generation and mutation annotation (i.e., impact at protein level):
	
	Consensus sequences are generated using bcftools (consensus -s sample.filtered.vcf.gz -f reference.fasta > sample.consensus.fasta) based on the vcf file containing the validated mutations. As for the Illumina pipeline, variant annotation is performed using snpEff 4.1l available with Snippy (see above).


.. note::

**PRIMER CLIPPING:** An extra parameter to enable primer removal using iVar (https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1618-7) is available within the settings “Mutation detection and consensus generation” section, for both Illumina and ONT. The procedure is an adaptation of the iVar CookBook (https://github.com/andersen-lab/paper_2018_primalseq-ivar/blob/master/cookbook/CookBook.ipynb). 
	
		***Users can request  trimming of primer sequences of several predefined Primer pool sets:
		
			-- SARS-CoV-2 Primal Scheme V3 (https://github.com/artic-network/artic-ncov2019/blob/master/primer_schemes/nCoV-2019/V3/nCoV-2019.tsv)
			
			-- SARS-CoV-2 Primal Scheme V4.1 (https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4.1)
			
			-- Monkeypox Primal Scheme from Welkers, Jonges and van den Ouden (https://www.protocols.io/view/monkeypox-virus-whole-genome-sequencing-using-comb-n2bvj6155lk5/v1)
			
			-- Monkeypox Primal Scheme from Chen et al. (https://www.protocols.io/view/monkeypox-virus-multiplexed-pcr-amplicon-sequencin-5qpvob1nbl4o/v2)
			
		Please contact us if you want to add new Primer pools to the online tool



**Masking low coverage regions:**

	**msa_masker.py** (https://github.com/rfm-targa/BioinfUtils/blob/master/FASTA/msa_masker.py; kind contribution of Rafael Mamede).
	
	This script substitutes positions with a low depth of coverage in a Multiple Sequence Alignment (MSA) with 'N'. The depth of coverage value below which the process masks positions is user-selected (see  “User-defined parameters”). It will not mask gaps/indels contained in the aligned sequences.
	
	-i: input FASTA file that contains a MAFFT nucleotide alignment enrolling the reference sequence (first sequence of the alignment) and consensus sequence(s) to be masked.
	
	-df: the coverage files (.depth)
	
	-r: define the reference sequence selected by the users (.fasta format) 
	
	-c: Positions with a depth value equal or below the value of this argument will be substituted by N (default= “mincov” - 1).
	
	
	**MAPPING VISUALIZATION**
					
	**Integrative Genomics Viewer** (http://software.broadinstitute.org/software/igv/)
	
		inputs: reference file (.fasta); mapping file (.bam; .bai)
		

.. important::
	INSaFLU allows users to configure key parameters for variant detection and consensus generation. **Settings** can be user-defined for the whole user account (tab “Settings”), for each project (after project creation) or for individuals samples within a project. 
	When parameters are changed for a given sample within a Project, the sample is automatically re-analysed using the novel parameters and re-inserted in the Project.
	See details in https://insaflu.readthedocs.io/en/latest/data_analysis.html#user-defined-parameters



Coverage analysis
-----------------

*Description*

This module yields a deep analysis of the coverage for each per sample by providing the following data: depth of coverage per nucleotide site, mean depth of coverage per locus, % of locus size covered by at least 1-fold and % of locus size covered by at least a user-defined "mincov" threshold (this parameter is user-selected for a Project or for a given sample within a Project). The latter constitutes the guide for consensus generation, i.e., consensus sequences are exclusively provided for locus fulfilling the criteria of having Y% of their size covered by at least X-fold (X = mincov; Y = minimum horizontal coverage) (see sections “Variant detection and consensus generation” and “User-defined parameters”). Coverage data is provided both in tabular format and interactive plots.

*Software version/settings*

.. note::
   	
	**Script used to generate Coverage statistics:**
	
	**getCoverage.py** (https://github.com/monsanto-pinheiro/getCoverage, by Miguel Pinheiro) (version v1.1; date 15.01.2018)
   
  	 	-i: define the input files, i.e, the coverage files (.depth.gz)
   
  		-r: define the reference sequence selected by the users (.fasta format) 
   
  		-o: defines the output file name (tab-separated value)
		
		
	**Script used to mask low coverage regions**

	**msa_masker.py** (https://github.com/rfm-targa/BioinfUtils/blob/master/msa_masker.py; kind contribution of Rafael Mamede)
	
	This script substitutes positions with a low depth of coverage in a Multiple Sequence Alignment (MSA) with 'N'. The depth of coverage value below which the process masks positions is user-selected (see  “User-defined parameters”). It will not mask gaps/indels contained in the aligned sequences.
	
	-i: input FASTA file that contains a MAFFT nucleotide alignment enrolling the reference sequence (first sequence of the alignment) and consensus sequence(s) to be masked.
	
	-df: the coverage files (.depth) 
	
	-r: define the reference sequence selected by the users (.fasta format) 
	
	-c: Positions with a depth value equal or below the value of this argument will be substituted by N (default= “mincov” - 1).

		

Alignment/Phylogeny
-------------------

*Description*
 
This module uses filtered nucleotide consensus sequences and performs refined nucleotide/protein sequence alignments and phylogenetic inferences. These outputs are automatically re-build and updated as more samples are added to user-restricted INSaFLU projects, making continuous data integration completely flexible and scalable. 

Users can also easily color the phylogenetic tree nodes and/or display colored metadata blocks next to the tree according to any combination of metadata variables, which facilitates the integration of relevant epidemiological and/or clinical data towards an enhanced genome-based pathogen surveillance. 

*Software version/settings*

.. note::
  	**MAUVE** (http://darlinglab.org/mauve/mauve.html) (version 2.4.0; date 15.01.2018)
   
   		progressiveMAUVE module (default settings): this algorithm is applied to perform primary draft alignments, and has the particular advantage of automatically concatenating multi-fasta input sequences during whole-genome alignments construction.
		
		input file: filtered nucleotide consensus sequences for each sample, one per each amplicon target (which are , in general, influenza CDSs) and another for the whole-genome sequence (i.e., the set of sequence targeted by the amplicon-based NGS shema, which, in general, is the pool of main 8 influenza CDSs). xmfa to fasta conversion is carried out using "convertAlignment.pl" (https://github.com/lskatz/lyve-SET/blob/master/scripts/convertAlignment.pl
		
		(default settings)
		
	**MAFFT**  (https://mafft.cbrc.jp/alignment/software/) (version 7.313; date 15.01.2018)

		For nucleotide alignments:
		
			input file: progressiveMAUVE-derived draft alignments (multifasta format), one per each locus and another for the whole-genome sequence 
		
			(default settings)
		
		For amino acid alignments:
		
			--amino: assume the sequences are in amino acid.
		
	**FastTree**  (http://www.microbesonline.org/fasttree/) (version 2.1.10 Double precision; date 15.01.2018)
	
			Double-precision mode: suitable for resolving very-short branch lengths accurately (FastTreeDbl executable)
			
			-nt: defines the input nucleotide alignment, which is a MAFFT-derived refined alignments (multifasta format). Alignments to be run include one per each locus and another for the whole-genome sequence.
			
			--gtr: defines the Generalized time-reversible (GTR) model of nucleotide evolution (CAT approximation with 20 rate categories)
			
			-boot: defines the number resample (-boot 1000)
			
	**Seqret** EMBOSS tool (http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/) (version 6.6.0.0; date 15.01.2018)
	
		input file: nucleotide alignments in FASTA (.fasta) to be converted in NEXUS (.nex) format 
	
	**MSAViewer**  (http://msa.biojs.net/) (latest; date 15.01.2018)
	
		input files: consensus nucleotide alignments for each locus and for the consensus 'whole-genome' sequence (upon concatenation of all individual locus); and amino acid alignments for the encoded proteins
		
	**Phylocanvas** (http://phylocanvas.org/) (version 2.8.1; date 15.01.2018)
	
		input files: phylogenetic tree obtained from each locus-specific nucleotide alignment and from the alignment of the 'whole-genome' sequences (upon concatenation of all individual locus)

		Metadata visualization tools were built with great contribution from Luís Rita: https://github.com/warcraft12321

Intra-host minor variant detection (and uncovering of putative mixed infections)
--------------------------------------------------------------------------------

*Description*

This module uses mapping data for the set of samples from each user-restricted INSaFLU project and provides a list of minor intra-host single nucleotide variants (iSNVs), i.e., SNV displaying intra-sample frequency between 1- 50%. This output is automatically re-build and cumulatively updated as more samples are added to each INSaFLU project, making continuous data integration completely flexible and scalable. Plots of the proportion of iSNV at frequency at 1-50%  (minor iSNVs) and at frequency 50-90% detected for each sample are also provided as mean to a guide the uncovering of putative mixed infections (exemplified in the Figure). INSaFLU flags samples as “putative mixed infections” based on intra-host SNVs if the following cumulative criteria are fulfilled: the ratio of the number of iSNVs at frequency at 1-50%  (minor iSNVs) and 50-90% falls within the range 0,5-2,0 and the sum of the number of these two categories of iSNVs exceeds 20. Alternatively, to account for mixed infections involving extremely different viruses (e.g., A/H3N2 and A/H1N1), the flag is also displayed when the the sum of the two categories of iSNVs exceeds 100, regardless of the first criterion. 

.. image:: _static/graph_mixed.png

*Software version/settings*

.. note::
   **Freebayes** (https://github.com/ekg/freebayes) (version v1.1.0-54-g49413aa; date 15.01.2018)
   
   		--min-mapping-quality: excludes read alignments from analysis if they have a mapping quality less than Q (--min-mapping-quality 20)
   		
   		--min-base-quality: excludes alleles from iSNV analysis if their supporting base quality is less than Q (--min-base-quality 20)
   		
   		--min-coverage: requires at least 100-fold of coverage to process a site (--min-coverage 100)
   		
   		--min-alternate-count: require at least 10 reads supporting an alternate allele within a single individual in order to evaluate the position (--min-alternate-count 10)
   		
   		--min-alternate-fraction: defines the minimum intra-host frequency of the alternate allele to be assumed (--min-alternate-fraction 0.01). This frequency is contingent on the depth of coverage of each processed site since min-alternate-count is set to 10, i.e., the identification of iSNV sites at frequencies of 10%, 2% and 1% is only allowed for sites with depth of coverage of at least 100-fold, 500-fold and 1000-fold, respectively.

Algn2pheno
--------------------------------------------------------------------------------

*Description*

The align2pheno module in INSaFLU performs the screening of genetic features potentially linked to specific phenotypes. Aln2pheno currently screens SARS-CoV-2 Spike amino acid alignments in each SARS-CoV-2 project against three default "genotype-phenotype" databases: the Carabelli mutations, the COG-UK Antigenic mutations and the Pokay Database (detailed below). Align2pheno reports the repertoire of mutations of interest per sequence and their potential impact on phenotype.

.. note::
   **Algn2pheno** (https://github.com/insapathogenomics/algn2pheno)
   
   		INSaFLU only runs the align2pheno module over Spike amino acid sequences with less than 10% of undefined amino acids (i.e., positions below the coverage cut-off; labelled as “X” in the protein alignments/sequences).
   		
   		Software and databases versions are provided in a log file in each run.

*Databases*

**Carabelli Database**

Description: Database of Spike amino acid mutations in epitope residues listed in Carabelli et al, 2023, 21(3), 162–177, Nat Rev Microbiol (https://doi.org/10.1038/s41579-022-00841-7), Figure 1.

Source: https://github.com/insapathogenomics/algn2pheno/blob/main/tests/DB_SARS_CoV_2_Spike_EpitopeResidues_Carabelli_2023_NatRevMic_Fig1.tsv (prepared and adapted for align2pheno based on https://doi.org/10.1038/s41579-022-00841-7)

**Pokay Database**

Description: Database of Spike amino acid mutations adapted from the curated database available through the tool Pokay, which includes a comprehensive list of SARS-CoV-2 mutations, and their associated functional impact (e.g., vaccine efficacy, pharmaceutical effectiveness, etc.) collected from literature. Made available by the CSM Center for Health Genomics and Informatics, University of Calgary.

Source: https://github.com/nodrogluap/pokay/tree/master/data


**COG-UK Antigenic Mutations Database**

Description: Database of Spike amino acid mutations adapted from the COG-UK Antigenic Mutations Database that includes “Spike amino acid replacements reported to confer antigenic change relevant to antibodies, detected in the UK data. The table lists those mutations in the spike gene identified in the UK dataset that have been associated with weaker neutralisation of the virus by convalescent plasma from people who have been infected with SARS-CoV-2, and/or monoclonal antibodies (mAbs) that recognise the SARS-CoV-2 spike protein.” Made available by the COVID-19 Genomics UK (COG-UK) Consortium through the COG-UK/Mutation Explorer.

Source: https://sars2.cvr.gla.ac.uk/cog-uk/


Nextstrain Datasets
--------------------------------------

*Description*

This module allows the creation of datasets for further in-depth phylogenetic analysis using Nextstrain (https://docs.nextstrain.org/en/latest/index.html). This provides an advanced vizualization and exploration of phylogenetic and genomic data, allowing the integration of geographic and temporal data and further user-provided metadata.

Currently, INSaFLU allows the creation of Datasets using virus-specific Nextstrain builds (seasonal Influenza, SARS-CoV-2 and Monkeypox) as well as a "generic" build that can be used for any pathogen.

More details here: https://github.com/INSaFLU/nextstrain_builds

*Builds*

**Seasonal influenza**

INSaFLU allows running four Nexstrain builds for the seasonal influenza (A/H3N2, A/H1N1/, B/Victoria and B/Yamagata), which are simplified versions of the Influenza Nextstrain builds available at https://github.com/nextstrain/seasonal-flu

So far, influenza analyses are restricted to the Hemagglutinn (HA) coding gene. The reference HA sequences used for site (nucleotide  / amino acid) numbering in the output JSON files are:

- H1N1PDM: A/California/07/2009(H1N1) (https://www.ncbi.nlm.nih.gov/nuccore/CY121680.1/)
- H3N2: A/Beijing/32/1992 (https://www.ncbi.nlm.nih.gov/nuccore/U26830.1/)
- VIC: Influenza B virus (B/Hong Kong/02/1993) (https://www.ncbi.nlm.nih.gov/nuccore/CY018813.1/)
- YAM: Influenza B virus (B/Singapore/11/1994) (https://www.ncbi.nlm.nih.gov/nuccore/CY019707.1/)

**Avian influenza** (under construction)

INSaFLU allows running Nexstrain builds for the avian influenza (A/H5N1), which are a simplified version of the Nextstrain builds available at https://github.com/nextstrain/avian-flu

So far, Nextstrain avian influenza can be launched for the Hemagglutinn (HA), Neuraminidase (NA) and polymerase protein PB2 (PB2) coding genes. The reference sequences used for site (nucleotide  / amino acid) numbering in the output JSON files are:

- HA: Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) hemagglutinin (HA) (https://www.ncbi.nlm.nih.gov/nuccore/AF144305.1/)
- NA: Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) neuraminidase (NA) (https://www.ncbi.nlm.nih.gov/nuccore/AF144304.1)
- PB2: Influenza A virus (A/Goose/Guangdong/1/96(H5N1)) polymerase (PB2)(https://www.ncbi.nlm.nih.gov/nuccore/AF144300.1)


**SARS-CoV-2**

This build is a simplified version of the SARS-CoV-2 Nextstrain build available at https://github.com/nextstrain/ncov

The reference genome used for site (nucleotide  / amino acid) numbering and genome structure in the output JSON files is:

- Wuhan-Hu-1/2019 (https://www.ncbi.nlm.nih.gov/nuccore/MN908947)


**Monkeypox virus**

This build is a simplified version of the Monkeypox virus Nextstrain build available at https://github.com/nextstrain/monkeypox

The reference genome used for site (nucleotide  / amino acid) numbering and genome structure in the output JSON files is:

- MPXV-M5312_HM12_Rivers (https://www.ncbi.nlm.nih.gov/nuccore/NC_063383)


**Respiratory Syncytial Virus (RSV)**

This build is a simplified version of the RSV virus Nextstrain build available at https://github.com/nextstrain/rsv

The reference genome used for site (nucleotide  / amino acid) numbering and genome structure in the output JSON files is:

- RSV A: RSV-A/US/BID-V8469/2001 (https://www.ncbi.nlm.nih.gov/nuccore/KJ627695.1/)
- RSV B: RSVB/Homo sapiens/USA/MCRSV_208/1980 (https://www.ncbi.nlm.nih.gov/nuccore/MG642037.1/)


**Generic**

This build is a simplified version of the Nextstrain build available at https://github.com/nextstrain/zika

This generic build uses as reference sequence (as tree root and for mutation annotation) one of the reference sequences of the projects included in the Nextstrain dataset.

Currently, the generic build does not generate a Time-Resolved Tree. To do this you need to select the Generic with TimeTree option.


**Generic with TimeTree**

This build is similar to the Generic build, but it also builds a time tree, inferring a mutation rate from the sample dates. Like in the Generic build, one reference is required to align the dataset consensus sequences. Nonetheless, unlike in the Generic build, the reference is not specifically defined as the root, but the root is inferred from the data instead. To make use of this build, you need to accurately specify dates associated with each sample.



.. important::
	**To take advantage of temporal and geographical features of Nextstrain**, please make sure you provide:
	
	- **"collection date"** for all samples added to Nextstrain datasets. If no collection date is provided, INSaFLU will automatically insert the date of the analysis as the "collection date", which might (considerably) bias (or even break) the time-scale trees generated for influenza, SARS-CoV-2 and Monkeypox.
	
	- **"latitude" and "longitude"** AND/OR **"region", "country", "division" and/or "location"** columns in the metadata. These values will be screened against a vast database of "latitude and longitude" coordinates (https://github.com/INSaFLU/nextstrain_builds/blob/main/generic/config/lat_longs.tsv) to geographically place the sequences in the Nextstrain map.
	
	To add/update the Nextstrain metadata of a given Dataset, click in **"Metadata for Nextstrain"**, download the previous table, update it with new data and upload it. Then, click in the "hourglass" icon to Rebuild the Nexstrain outputs. Note: you can also add/update the metadata of sequences previously obtained with INSaFLU (i.e., consensus sequences coming from the "Projects" module), please follow these instructions: https://insaflu.readthedocs.io/en/latest/uploading_data.html#updating-sample-metadata (this option is not available for external sequences).



Metagenomics virus detection
+++++++++++++++++++++++++++++++++++++

The TELEVIR  bioinformatics component of INSaFLU is a modular pipeline for the identification of viral sequences in metagenomic data (both Illumina and ONT data). 

It is composed of these main steps (detailed below):

1. Read quality analysis and improvement [optional]
2. Extra-filtering [optional].
3. Viral Enrichment [optional].
4. Host Depletion [optional].
5. *De novo* assembly of the reads [optional].
6. Identification of the viral sequences.
	- Using reads.
	- Using contigs (if assembled).
	- Using several reference databases.
7. Selection of viral TAXID and representative genome sequences for confirmatory re-mapping
8. Remapping of the viral sequences against selected reference genome sequences.
9. Reporting


.. image:: _static/televir_workflow_update_20221219.png

Simplified illustration of the main steps of the modular INSaFLU-TELEVIR bioinformatics pipeline for metagenomics virus diagnostics.

*Find below the description of the main bioinformatics steps (Description, Current software versions and settings)*


Read quality analysis and improvement
--------------------------------------


*Description*

This step takes the input single- or paired-end reads (fastq.gz format; ONT or Illumina) and produces quality processed reads, as well as quality control reports for each file, before and after this step. This module overlaps the two components (virus detection and surveillance) of the platform. 

*Software version/settings*

.. note::

	**## ILLUMINA data ##**
	
   	**FastQC** (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) (version 0.11.5; date 15.01.2018)

		input: single- or paired-end reads (fastq.gz format) (e.g., sample_L001_R1_001.fastq.gz and sample_L001_R2_001.fastq.gz for Illumina technology reads)
		
		--nogroup option: all reports will show data for every base in the read. 
		
	**Trimmomatic** (http://www.usadellab.org/cms/index.php?page=trimmomatic) (version 0.27; date 15.01.2018)
	
		input: single- or paired-end reads (fastq.gz format) (e.g., sample_L001_R1_001.fastq.gz and sample_L001_R2_001.fastq.gz for Illumina paired-end reads)
	
		ILLUMINACLIP: To clip adapters from the input file using user-specified adapter sequences. ILLUMINACLIP:<ADAPTER_FILE>:3:30:10:6:true
		
		HEADCROP: Cut the specified number of bases from the start of the read. Range: [0:100]. If value equal to 0 this parameter is excluded.
		
		CROP: Cut the read to a specified length. Range: [0:400]. If value equal to 0 this parameter is excluded.
	
		SLIDINGWINDOW: perform a sliding window trimming, cutting once the average quality within the window falls below a threshold (default: SLIDINGWINDOW:5:20, where 5 refers to window and 20 to the minimum average quality)
	
		LEADING: cut bases off the start of a read, if below a threshold quality (default: LEADING:3). This will allow discarding bases with very quality or N bases (quality score of 2 or less).
	
		TRAILING: cut bases off the end of a read, if below a threshold quality (default: TRAILING:3). This will allow discarding bases with very quality or N bases (quality score of 2 or less).
	
		MINLEN: drop the read if it is below a specified length (MINLEN:35)
	
		TOPHRED33:  Convert quality scores to Phred-33
		
	**## Oxford Nanopore Technologies (ONT) data ##**
		
	**NanoStat** (https://github.com/wdecoster/nanostat) (version 1.4.0)
		
		input: ONT reads (fastq.gz format) 

	**NanoFilt** (https://github.com/wdecoster/nanofilt) (version 2.6.0)
	

		**-q (QUALITY)**: Filter on a minimum average read quality score. Range: [5:30] (default: 10)
		
		**-l (LENGTH)**: Filter on a minimum read length. Range: [50:1000]. (default: 50)
		
		**--headcrop**: Trim n nucleotides from start of read. Range: [1:1000]. If value equal to 0 this parameter is excluded. (default: 70)
		
		**--tailcrop**: Trim n nucleotides from end of read. Range: [1:1000]. If value equal to 0 this parameter is excluded. (default: 70)
		
		**--maxlength**: Set a maximum read length. Range: [100:50000]. If value equal to 0 this parameter is excluded. (default: 0)
		

	**RabbitQC** (https://github.com/ZekunYin/RabbitQC)  (version 0.0.1)**
		
		input: ONT reads (fastq.gz format) pre- and post- quality improvement with NanoFilt
		
		Files between 50 - 300 MB are downsized to ~50 MB before analysis by randomly sampling reads using fastq-sample from fastq-tools package https://github.com/dcjones/fastq-tools (developed by Daniel C. Jones dcjones@cs.washington.edu)


.. important::
	INSaFLU allows users to configure key parameters for reads quality analysis in the tab **“Settings”**. 
	
	**Please check your settings before uploading new samples to your account.**
	
	See details in https://insaflu.readthedocs.io/en/latest/data_analysis.html#user-defined-parameters


Extra Filtering
------------------

*Description*

This step **remove reads enriched in low complexity regions** (e.g., homopolymeric tracts or repeat regions), which are a common source of false-positive bioinformatics hits**. This step is directly performed using over raw reads (if QC was turned OFF) or quality processed reads (if QC was turned ON).


.. note::

*Software*

	**PrinSeq++** (https://github.com/Adrian-Cantu/PRINSEQ-plus-plus) 
	
	


Viral Enrichment
------------------

*Description*

This step **retains potential viral reads** based on a rapid and permissive classification of the reads against a viral sequence database. If "Extra-filtering" is OFF, this step is directly performed using over raw reads (if QC was turned OFF) or quality processed reads (if QC was turned ON).


.. note::

*Software*

	**Centrifuge** (https://github.com/centrifuge/) 
	
	**Kraken2** (https://github.com/DerrickWood/kraken2)
	
	
*Databases*

	**Virosaurus90v 2020_4.2**  (https://viralzone.expasy.org/8676)
	
	**NCBI refseq viral genomes** release 4 (https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)



Host depletion
------------------

*Description*

This step **removes potential host reads** based on reference-based mapping against host genome sequence(s). Mapped reads are treated as potential host reads and removed. If the Extra filetring and Viral enrichment steps were turned OFF, host depletion will be directly performed over raw reads (if QC was turned OFF) or quality processed reads (if QC was turned ON).

.. note::

*Software*

	**BWA**  (https://github.com/lh3/bwa)
	
	**Minimap2** (https://github.com/lh3/minimap2)

	
	
*Host reference sequences**

	**Human reference genome hg38 - NCBI accid GCA_000001405.15**


*De novo* Assembly
------------------

*Description*

This step performs *de novo* assembly using reads retained after the "Viral enrichment" and/or "Host depletion" steps. If the latter steps were turned OFF, assembly will be directly performed using raw reads (if QC was turned OFF) or quality processed reads (if QC was turned ON).


.. note::

*Software*

	**SPAdes** (https://github.com/ablab/spades)
	
	**Raven** (https://github.com/lbcb-sci/raven)
	


Identification of the viral sequences
--------------------------------------

*Description*

This step **screens reads and contigs against viral sequence databases**, generating intermediate read and/or contig classification reports of viral hits (TAXID and representative accession numbers) potentially present in the sample. The top viral hits will be selected for confirmatory re-mapping (see next steps).


.. note::

*Software*

**Reads classification**

	**Centrifuge** (https://github.com/DaehwanKimLab/centrifuge)
	
	**FastViromeExplorer** (https://github.com/saima-tithi/FastViromeExplorer)
	
	**Kraken2** (https://github.com/DerrickWood/kraken2)
	
	**Krakenuniq** (https://github.com/fbreitwieser/krakenuniq)
	
	**Kaiju** (https://github.com/bioinformatics-centre/kaiju)
	

**Contigs classification**
	
	**Blast** (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
	
	**FastViromeExplorer** (https://github.com/saima-tithi/FastViromeExplorer)


*Databases*


	**Virosaurus90v 2020_4.2**  (https://viralzone.expasy.org/8676)
	
	**NCBI refseq viral genomes** release 4 (https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
	
	**RefSeq complete viral genomes/proteins**, as modified for the kraken2 and centrifuge databases.


Selection of viral TAXID and representative genome sequences for confirmatory re-mapping 
-----------------------------------------------------------------------------------------

*Description*

In this step, the previously identified viral hits (TAXID) are selected for confirmatory mapping against reference viral genome(s) present in the available databases. Viral TAXIDs are selected, up to a maximum number of hits*, as follows:

Viral TAXIDs are selected, up to a maximum of number of hits, as follows:

- 1º - Viral hits corresponding to phages are removed from classification reports.
- 2º - TAXIDs present in both intermediate classification reports (reads and contigs) are selected;
- 3º - additional TAXIDs are selected across the read classification report (by number of hits, in decreasing order) and contigs classification report (by number of hits and total length of matching sequences, from top-down) until reaching the defined maximum number of hits to be selected
- 4º - Representative sequences (accession ID) of the selected TAXID are queried from internal collection of databases (see next step).

*currently, this number is set as 15 as default, but it is to be user-defined

*Databases*

	**NCBI Taxonomy** (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)


Remapping of the viral sequences against selected reference genome sequences. 
--------------------------------------------------------------------------------

*Description*

This step **maps reads and/or contigs against representative genome sequences of the selected viral TAXIDs** collected in the previous step. The reference sequences are collected from available databses (see below). If a given representative TAXID/sequence is present in more than one database, priority is given to NCBI refseq viral genomes and Virosaurus.

On note, reads are also mapped against any contigs that successfully map against reference sequences. TAXIDs that were not automatically selected for this confirmatory remapping step (but that were present in the intermediate reads and/or contigs classification reports) can still user-selected for mapping at any time.


.. note::

*Software*

	**Snippy** (https://github.com/tseemann/snippy)
	
	**Bowtie2** (https://github.com/BenLangmead/bowtie2)
	
	**Minimap2** (https://github.com/lh3/minimap2)
	

*Databases*

	**Virosaurus90v 2020_4.2**  (https://viralzone.expasy.org/8676)
	
	**NCBI refseq viral genomes** release 4 (https://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
	
	**RefSeq complete viral genomes/proteins**, as modified for the kraken2 and centrifuge databases.
	
	**NCBI Taxonomy** (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/)


Reporting
----------

TELEVIR reports are generated per **Workflow**, per **Sample** (combining non-redundant hits detected across workflows) and per **Project** (combining several samples), with a decreasing level of detail.

The workflows culminate in **user-oriented reports with a list of top viral hits, each accompanied by several robust and diagnostic-oriented metrics, statistics and visualizations**, provided as (interactive) *tables* (intermediate and final reports), *graphs* (e.g., coverage plots, Integrative Genomics Viewer visualization, Assembly to reference dotplot) and *multiple downloadable output files* (e.g., list of the software parameters, reads/contigs classification reports, mapped reads/contigs identified per each virus; reference sequences, etc).

See the description of the reports and outputs here: https://insaflu.readthedocs.io/en/latest/metagenomics_virus_detection.html#televir-output-visualization-and-download

*Sorting*

In order to simplify the final reports (per **Sample** and per **Workflow**) and facilitate the identification of potential false positive hits (often arising from cross-mapping with true positive hits), **sample-specific viral references are grouped together by mapping affinity, as measured by shared mapped reads.**

In summary, reference hits selected for remapping are first pooled across sample workflows, excluding references with no mapped reads. A neighbor-joining tree based on shared reads is then constructed. A maximum of 100000 reads are used after filtering out reads shared by 5% of references or less. A filter is then applied on the inner nodes of the tree that considers the distribution of shared reads as summarized by two statistics:

	- private_reads : the proportion of reads found to map only to descendents of that node.
	- pairwise_min : the minimum reciprocal proportion of reads mapped between every pair of descendents (i.e. all samples must share at least X% (user-defined) of their reads with another sample among descendents of the same node).

Thresholds for these two statistics are defined beforehand: private_reads is set by the user through the parameter “--r-overlap” (defaults to 50 %); pairwise_min is a TELEVIR constant set to  5 %. After filtering, nodes are sorted by the total number of reads mapped to their descendents. Finally, references (tree leaves) are mapped to the filtered inner nodes and sorted accordingly. Orphaned leaves and references with no mapped reads are appended last. 


*Warnings and Flags*

TELEVIR reports provide specific Warnings for  bioinformatics “artifacts” commonly yielding false-positive taxid assignments. Calculations depend on the flag-type, a user defined variable (TELEVIR Settings – Reporting – Final Report - Flagging and Sorting – --flag-type, default: viruses), and target broad characteristics of main input types. 

Two flag-types currently exist for **viruses** (oriented to shotgun viral metagenomics) and **probes** (oriented for probe-based NGS target panels)

## Flag-type **"viruses"** (default)

- *"Likely False Positive"*: when most reads map to a very small region of the reference sequence, i.e., hits with high “DepthC" but low “Depth” and low "Cov (%)". Flagged for hits with DepthC / Depth > 10 and Cov (%) < 5%.
	
- *"Vestigial Mapping"*: when only a vestigial amount of reads (<= 2) mapped.


## Flag-type **"probes"** 

- *"Likely False Positive"*: When the reference genome is not sufficiently covered as a function of the number of the proportion of Windows Covered, calculated as above. Flagged for hits with Windows Covered <= 50 % (calculated from the fraction presented)
	
- *"Vestigial Mapping"*: when only a vestigial amount of reads (<= 2) mapped.


User-defined parameters (UNDER CONSTRUCTION)
+++++++++++++++++++++++++++++++++++++++++++++++

INSaFLU allows turning ON/OFF specific modules and user-defined configuration of key parameters for reads quality analysis, INSaFLU and TELEVIR projects. Settings can be user-defined for the whole user account (tab “Settings”), for each project (just after project creation) or for individual samples within a project (click in the "Magic wand" icon).

**Please choose your settings before uploading new samples to your account.**

Example:

.. image:: _static/01_global_settings.gif


Read quality analysis and improvement control (QC)
--------------------------------------------------

**##ILLUMINA / Ion Torrent data##**

Users can change the following **Trimmomatic** settings (see http://www.usadellab.org/cms/index.php?page=trimmomatic):

**ILLUMINACLIP**: To clip the Illumina adapters from the input file using the adapter sequences. ILLUMINACLIP:<ADAPTER_FILE>:3:30:10:6:true (default: Not apply)
		
**HEADCROP**: <length> Cut the specified number of bases from the start of the read. Range: [0:100]. If value equal to 0 this parameter is excluded. (default = 0)

**CROP**:<length> Cut the read to a specified length. Range: [0:400]. If value equal to 0 this parameter is excluded. (default = 0)

**SLIDINGWINDOW**:<windowSize> specifies the number of bases to average across Range: [3:50]. (default = 5)

**SLIDINGWINDOW**:<requiredQuality> specifies the average quality required Range: [10:100]. (default = 20)

**LEADING**:<quality> Remove low quality bases from the beginning. Range: [0:100]. If value equal to 0 this parameter is excluded. (default = 3)

**TRAILING**:<quality> Remove low quality bases from the end. Range: [0:100]. If value equal to 0 this parameter is excluded. (default = 3)

**MINLEN**:<length> This module removes reads that fall below the specified minimal length. Range: [5:500]. (default = 35)

NOTE: "Trimming occurs in the order which the parameters are listed"

**## Oxford Nanopore Technologies (ONT) data ##**

Users can change the following **NanoFilt** settings (see: https://github.com/wdecoster/nanofilt)

**QUALITY**: Filter on a minimum average read quality score. Range: [5:30] (default: 10)

**LENGTH**: Filter on a minimum read length. Range: [50:1000]. (default: 50)

**HEADCROP**:  Trim n nucleotides from start of read. Range: [1:1000]. If value equal to 0 this parameter is excluded. (default: 70)

**TAILCROP**: Trim n nucleotides from end of read. Range: [1:1000]. If value equal to 0 this parameter is excluded. (default: 70)

**MAXLENGTH:** Set a maximum read length. Range: [100:50000]. If value equal to 0 this parameter is excluded. (default: 0)


Mapping, Variant Calling
-------------------------

**##ILLUMINA / Ion Torrent data##**

Users can change the following **Snippy** settings (see also https://github.com/tseemann/snippy):

**--mapqual**: minimum mapping quality to accept in variant calling (default = 20)

**--mincov**: minimum number of reads covering a site to be considered (default = 10)

**--minfrac**: minimum proportion of reads which must differ from the reference, so that the variant is assumed in the consensus sequence (default = 0.51)


**## Oxford Nanopore Technologies (ONT) data ##**

Users can change the following settings:

**Medaka model** (default: r941_min_high_g360) (see: https://nanoporetech.github.io/medaka/)

**Minimum depth of coverage per site** (equivalent to --mincov in Illumina pipeline) (default: 30) 

**Minimum proportion for variant evidence** (equivalent to --minfrac in Illumina pipeline) (default: 0.8). Note: medaka-derived mutations with frequencies below the user-defined “minfrac” will be masked with an “N”. 


Consensus generation (horizontal coverage cut-off) and Masking
--------------------------------------------------------------
Users can select the **Minimum percentage of horizontal coverage to generate consensus**. This threshold indicates the **Minimum percentage of locus horizontal coverage** with depth of coverage equal or above –mincov (see Mapping settings) to generate a consensus sequence for a given locus. Range: [50:100] (default = 70)

In Projects setting, users can also **mask (i.e., put Ns) specific regions (or sites)** of the consensus sequences for all (or individual) samples within a given Project. This feature is especially useful for masking the start/end of the sequences or known error-prone nucleotide sites. 


.. image:: _static/masking_consensus_projects.png

**Masking summary:**

Undefined nucleotides (NNN) are automatically placed in: 
i) low coverage regions (i.e., regions with coverage below --mincov); 
ii) regions (or sites) selected to be masked by the user (in Projects settings); 
iii) for ONT data, medaka-derived mutations with frequencies below the user-defined “minfrac” (i.e. Minimum proportion for variant evidence).


