# Incomplete resection of the icEEG seizure onset zone is not associated with post-surgical outcomes
Gascoigne et al. (brief communication, submitted January 2024).
The code in this repository can be used to automatically detect onsets using the imprint and Epileptigenicity Index (Bartolomei et al., 2008) algortihms

### Data formatting 
For each subject, icEEG data and accompanying metadata are required. See *data* folder to see example data from SWEZ dataset (link) formatted in the appropriate way for this analysis (*STILL NEED TO ADD THIS*)

### Onset localisation and channel to ROI conversion
The **onset_detec_s_mahal.m** code creates a table (*subquestions/final_output.mat*) presenting all data required for downstream analyses:
  - Subject ID
  - Seizure IDs
  - Clinically labelled onset (channel-wise, Lausanne-120, Lausanne-250)
  - Automatically labelled onsets for each seizure (channel-wise, Lausanne-120, Lausanne-250) for each onset localisation algorithm
  - One ALO created across seizure onsets (regions included in >=50% of onsets) (channel-wise, Lausanne-120, Lausanne-250) for each onset localisation algorithm

### Downstream analysis
All scripts for performing downstream analysis are included in the **subquestions** folder. Within this folder **final.m** calls all other scripts and will produce tables and figures which were used in this paper.

In the main analysis, we looked at onsets based on the Lausanne-120 parcellation scheme, this code can be easily adapted to iterate over channel-wise, Lausanne-120, and Lausanne-250 parcellation schemes by uncommenting lines 40 and 86 (**CHECK THAT THIS IS STILL THE CASE**) 

#### Setting parameters for analysis
- *parc*: The parcellation scheme of interest (can be set to "chan", "roi_120", or "roi_250"). You only need to set this parameter if you are not iterating through parcellation schemes (i.e., lines 40 and 86 remain commented out). The code is currently set to perform analysis using the Lausanne-120 parcellation scheme (Line 17: *parc = "roi_120";*)
- *n_perm*: The number of permutations to perform to create a "chance" distribution of AUCs from which we can compute the p-value. The smaller this number is set, the faster the code will run, however we suggest that this should not be set below 500 if the p-values are to be reliable. Our work used 1000 permutations to build the distribution of AUCs (Line 18: *n_perm = 1000;*)
- *out_thresh*: The maximum ILAE values that is labelled as a 'favourable' surgical outcome. In this work we have set this value to 2 such that ILAE 1-2 and ILAE 3+ are labelled as having 'favourable' and 'unfavourable' surgical outcomes, respectively (Line 20: *out_thresh = 2;*)
- *consensus_thresh*: The threshold for inclusion of regions in the consensus onset (i.e. one onset per individual which captures regions that are prevalent across seizures). Here we will include regions present in at least half of the subject's seizures (Line 22: *consensus_thresh = 0.5;*)
- *most_resec_thresh*: The threshold above which subjects are considered as having 'most' of their onset resected. Here we have set this as 50%, so any resections removing at least half of the osnet is labelled as 'most resected' (Line 25: *most_resec_thresh = 0.5;*)
- *det_meths*: This selects which methods of onset detection will be considered (can be set to "clo", "imprint", or "EI"). In our main text, we used clinically labelled onsets (CLO) and automatically labelled onsets (ALO) based on the imprint algorithm. Supplementary analyses additionally included automatically detected onsets using the Epileptigenicity Index (Bartolomei et al., 2008) (Line 28:*det_meths* = ["clo", "imprint", "EI"];*)

#### Saving figures
The code saves figures for each step of the analysis, the user can determine if they want the figures to be saved (save_fig = 1) or not (save_fig = 0) on line 37. 
The location where the figures are stored can be changed, as this code produces the visualisations used in *Figure Three*, we have set the save location accordingly (Line 34: *save_folder = '../figures/paper_figures/Figure 3/'*).
Additionally, it is possible to change the filetype used for each of the figures (Line 37: *file_type = "svg";*), see MATLAB documentation to see the filetypes available). 

#### Step-by-step downstream analysis
##### 3.1 icEEG Seizure onset regions tend to be resected, but more complete resections are not associated with more favourable surgical outcomes
**is_onset_resected** The proportion of the seizure onset zone (CLO and IOLA) that was subsequently resected was computed and compared across outcome groups.  We elected to omit the volume calculations from this work as the volume of regions may not accurately represent the true volume of the onsets and resections, the code has been left in this repository for completeness. Output figures will include comparisons using both the number of regions and the estimated volume of onsets.

##### 3.2 Larger onsets and resections are not associated with surgical outcomes
**larger_onset** 
Here we compute the size of seizure onsets (both CLO and automatically captured) using both count of regions and volume (based on controls) and compare across surgical outcome groups. We elected to omit the volume calculations from this work as the volume of regions may not accurately represent the true volume of the onsets and resections, the code has been left in this repository for completeness. Output figures will include comparisons using both the number of regions and the estimated volume of onsets.

**larger_resection** 
Here we will compute the size of resections using both count of regions and volume (if using ROIS, based on controls) and compare across surgical outcome groups. Again, analyses based on volumes were omitted from this work but the code has been provided and figures will be produced. 
 

#### Figures
*Figure One* Panels B-F display the process from icEEG to automatically detected (consensus) onset using Lausanne-120 parcellation scheme. Panel G shows the clinically labelled onset, also using the Lausanne-120 parcellation scheme. Panels B-G can be reproduced by running **figure_1.m**
*Figure Three* is based on the figures generated when running **subquestions/final.m**
