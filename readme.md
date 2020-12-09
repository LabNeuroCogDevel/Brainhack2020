# LNCD PET tat2

## Data
 * `data/wide.csv` is a wide (547 column) format CSV with a row per participant visit.
The columns include individual measures for various surveys (behavioral) and brain (MR and PET scan) data.

 * `data/mtr_striatum_258_QC.csv` - MTR ROI mean and voxel coverage count ("NZmean") and QC (manual quality check binarized 1=pass, 0=fail)

 * `data/conn_adj*.csv` row per visit CSV files. A file per MR task type: rest1, task background, rest2. Columns are the widened (2d to 1d) upper triangle of each visits fisher-z transformed correlation adjacency matrix. This is the mean BOLD time series with each of 23 ROIS correlated to one another ( 23 choose 2 = 253). Fixation from 6 consecutive task aquistions are extracted and combined to make the "background" set. 
   * also see `data/conn_*.pkl` for pre-made (and untested) pyradigm data structure:  `reloaded = RegrDataset('data/conn_rest1.pkl')`

### images
per subject tat2 3d nii.gz images are on google dirve:
   https://drive.google.com/drive/folders/1WuCqimxEBnrb8MkqXawTKLGZi7frYcy1?usp=sharing

### wide.csv

  | prefix       | desc|
  | ------       | ---- |
  | `tat2*`      | time average T2* (proxy for iron?) |
  | `frogET*`    | anti saccade task performance (esp correct: `ncor`, `corlat`, `corsd` for reward and neutral)  |
  | `RT18.Score` | Risk taking measure (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3160867/)|
  | `UPPSP*`     | impulsivity (http://www.impulsivity.org/measurement/UPPS_P )|
  | `ysrasr.*`   | DSM-oriented scales (http://fcon_1000.projects.nitrc.org/indi/enhanced/assessments/asr.html) |
  

## Neuropredict

```bash
python3 -m pip install neuropredict pyradigm
np_regress -y data/conn_bkgrd_subset.pkl --impute_strategy most_frequent -e randomforestregressor
```

## Code

`pet_data.py` quick attempt at building pyradigm data structures


```python
from pyradigm import RegressionDataset as RegrDataset
from pet_data import PET
pet = PET('data/merged_data.csv')
# see pet.widedf.shape: 384 visits with 6,987 measures

##  get sesid, age, and all the uppsp measures

# pyradigm
upps = RegrDataset()
upps.description = "Urgency, Premeditation (lack of), Perseverance (lack of), Sensation Seeking, Positive Urgency, Impulsive Behavior Scale"
pet.add_subset_to(upps, '^uppsp_') 
# Urgency, Premeditation (lack of), Perseverance (lack of), Sensation Seeking, Positive Urgency, Impulsive Behavior Scale 
# 337 samplets, 6 features

# or as dataframe
upps_df = pet.col_subset('^uppsp_')
upps_df.columns # ['sesid', 'age', 'upps_pre', 'upps_pers', 'upps_ss', 'upps_pu', 'upps_tot', 'upps_negurg']

# also include other confound columns (works same for add_subest_to and col_subset)
upps_df = pet.col_subset('^sex|^uppsp_', 'uppsp_upps')
upps_df.columns # ['sesid', 'age', 'sex', 'pre', ...]

```

## Glossary
* *lunaid*  unique identifier for a participant. Repeats across visits
* *visitnum* longitudinal study. We expect 3 distinct visits. Not every participant finishes. 
* *sessionid*/*sesid*  session identifier unique for visit year. Behavioral and scan may happen on separate days. These will still have the same sesid.
* *vdate* MR visit date, *behavedate* date of behavioral battery (RT18, UPPS, YSR/ASR), *dtbzdate* may have returned if the scanner portion did not complete. this is the date MR was resumed. few have separate vdate and dtbzdate.
* *ROI*  brain region of interest. Anatomical or functionally defined region composed of many voxels. Voxel values are often averaged within this region.
* *MR* 
* *taT2* time averaged T2\* - average of all the BOLD signal time course.
* *R2'* R2 prime
* *MTR* magnetic transfer ratio
* *BOLD* 
* *ASY/YSR* adult/youth self report. 2 sets of questions, depending on age of participant.
* *PET* 
* *DTBZ* PET radio tracker administered second. vessicular monoamine transpoter/DA
* *Rac* Raclopride PET radio tracker administered first
