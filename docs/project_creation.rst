Project creation and scaling-up
===============================

One of the main goals of INSaFLU is to make data integration completely flexible and scalable in order to fulfill the analytical demands underlying laboratory surveillance throughout each flu epidemics. As such, INSaFLU allows users to create several projects (each one including multiple user-selected samples) and add more samples to each one as needed. In a dynamic manner, project outputs (e.g., gene- and genome-based alignments and phylogenetic trees) are automatically re-build and cumulatively updated as more samples are added to each project. The outputs are provided to be compatible with multiple downstream applications.

Creating and scaling-up a project
+++++++++++++++++++++++++++++++++

Within the *Projects* menu:

1. Go to *Projects* menu and choose *Create project*
....................................................

You are encouraged to create “umbrella” projects, such as projects enrolling same sub-type viruses from the same season that will be compared with the vaccine reference virus for a given flu season. 

You can designate the projects so that the name easily indicates the combination “virus sub-type/season/reference” (e.g. **A_H3N2_2017_18_vaccine_ref**)

.. image:: _static/create_project_1_create.png


2. Choose a *Project Name*, select a *Reference sequence* and *Save*
......................................................................

.. important::
   You should select a reference sequence (e.g., the vaccine strain from the current influenza season) that fits both your amplicon design (i.e., a multi-fasta file containing the set of reference sequences with the precise size of each “intra-amplicon” target sequence that you capture by each one of the RT-PCR amplicons) and the set of samples that will be compared (e.g., same sub-type viruses from the same season to be compared with the vaccine reference virus).

.. image:: _static/create_project_2_name_ref.png


3. Choose the software parameters to be applied to the project.
.................................................................

The selected parameters will be applied by default to every sample added to the project, so please set the parameters before assigning the first sample to the project.

Note: After project creation, you are still allowed to change the parameters for individual samples within the Project. The sample is automatically re-analysed using the novel parameters and re-inserted in the Project. 


4. Add the **samples** to be included in the **project**
........................................................

.. image:: _static/create_project_3_add_samples.png

Samples are processed immediately upon selection, so, at this time, users may start monitoring the Project progress by checking the number of samples in the following status: Processed (P); Waiting (W) and Error (E).

.. image:: _static/monitoring_project_status.png


5. Scale-up your **project**. 
.............................

You may add more samples to your **Project** project at any time.

.. image:: _static/create_project_4_scale_up.png


6. Modify software parameters for a given sample within a Project
..................................................................

Users can change the mapping parameters for individual samples within a Project. The sample is automatically re-analysed using the novel parameters and re-inserted in the Project (outputs are automatically re-calculated to integrate the “updated” sample). For instance, if the updated sample fulfill the criteria for consensus generation with the novel settings, it will be automatically integrated in the alignments and trees.

NOTE: Users can also re-run samples (with user-selected parameters) included in projects created before the 30 Oct 2020 update (see "Change log"). The updated samples will be flagged accordingly. 



7. Remove samples from your **project**. 
........................................

You may want to remove some samples from your project (e.g., for exclusively keeping samples with success for all 8 locus) 

.. image:: _static/create_project_5_remove_samples.png
 
  

Monitoring Projects' progress
+++++++++++++++++++++++++++++

INSaFLU projects are automatically run upon creation. So, at this time, users may start monitoring the Project progress by checking the number of samples in the following status: Processed (P); Waiting (W) and Error (E).


.. image:: _static/monitoring_project_status.png


