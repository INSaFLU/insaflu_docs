Frequently Asked Questions (FAQs)
=================================
FAQs
....
Under construction

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
