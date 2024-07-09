# Incomplete resection of the icEEG seizure onset zone is not associated with post-surgical outcomes
Gascoigne et al. (brief communication, submitted January 2024).
The code in this repository can be used to automatically detect onsets using the imprint, Epileptigenicity Index (Bartolomei et al., 2008), and Phase-Locked High-Gamma (Weiss et al., 2015) algorithms

### Data formatting 
For each subject, icEEG data and accompanying metadata are required. The table should have one row per seizure and include the following columns:
  - *segment_id*: Unique identifier for the seizure
  - *patient_id*: Subject identification number
  - *duration*: Duration of seizure in seconds
  - *segment_pre*: Number of seconds at the start of the recording that is considered pre-ictal (i.e., before the clinically labelled seizure onset time).
  - *segment_post*: Number of seconds at the end of the recording that is considered post-ictal (i.e., after the clinically labelled seizure offset time).
  - *segment_fs*: Sampling frequency of recording in Hertz (in this work, all seizures were sampled at 512Hz).
  - *segment_channel_labels*: Channel names for each channel in the recording (this is used to localise channels to regions of interest (ROIs). The example data does not include channel labels within the table, therefore the channels are labelled as their index. 
  - *segment_data*: icEEG time series (matrix with dimension #channels x #secondsrecorded*sampling frequency).

### Onset localisation and channel to ROI conversion
The **onset_detec.m** code creates a table (*tables/final_output.mat*) presenting all data required for downstream analyses:
  - Subject ID
  - Seizure IDs
  - Channel labels
  - Automatically labelled onsets for each seizure (channel-wise)
  
The following columns are not generated in the example data (as this information was not available, please see analysis repository to see an example output table with all columns)
  - Resection localisation (channel-wise, Lausanne-120, Lausanne-250)
  - Clinically labelled onset (channel-wise, Lausanne-120, Lausanne-250)
  - Automatically labelled onsets for each seizure (Lausanne-120, Lausanne-250) for each onset localisation algorithm
  - One ALO created across seizure onsets (regions included in >=50% of onsets) (channel-wise, Lausanne-120, Lausanne-250) for each onset localisation algorithm
  - ROI names (Lausanne-120, Lausanne-250) for visualisations
  - Channels to ROI conversion matrices (Lausanne-120, Lausanne-250)
  - Subject metadata (surgery outcome, surgery year, outcome year, operation type, seizure duration, seizure type, sex, age at epilepsy onset, outcome at year one (as is used as the outcome identifier in this work)). 

