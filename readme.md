# LNCD "PET" 

## Data
`data/merged_data.csv` is an ultra-wide format CSV with a row per participant visit.
The thousands of columns include individual measures for various surveys (behavioral) and brain (MR and PET scan) data.

`data/conn_adj.csv` is another row per visit CSV. Columns are the widened (2d to 1d) upper triangle of each visits fisher-z transformed correlation adjacency matrix. This is the mean BOLD time series with each of 23 ROIS correlated to one another ( 23 choose 2 = 253).

## Code

`pet_data.py` quick attempt at building pyradigm data structures


```python
from pet_data import PET
main_data = PET('data/merged_data.csv')

upps = RegrDataset()                                                                                           
upps.description = "Urgency, Premeditation (lack of), Perseverance (lack of), Sensation Seeking, Positive Urgency, Impulsive Behavior Scale"
main_data.add_subset_to(upps, '^uppsp_') 

```

## Glossary
* *lunaid*  unique identifier for a participant. Repeats across visits
* *visitnum* longitudinal study. We expect 3 distinct visits. Not every participant finishes. 
* *sessionid*/*sesid*  session identifier unique for visit year. Behavioral and scan may happen on separate days. These will still have the same sesid.
* *vdate* MR visit date, *behavedate* date of behavioral battery (RT18, UPPS, YSR/ASR), *dtbzdate* may have returned if the scanner portion did not complete. this is the date MR was resumed. few have separate vdate and dtbzdate.
* *ROI*  brain region of interest. Anatomical or functionally defined region composed of many voxels. Voxel values are often averaged within this region.
* *Rac* PET radio tracker administered first
* *taT2* time averaged T2\* - average of all the BOLD signal time course.
* *MTR* magnetic transfer ratio
* *DTBZ* PET radio tracker administered second
* *PET* 
* *MR* 
* *BOLD* 
* *RT18*
* *UPPS*
* *ASY/YSR* adult/youth self report. 2 sets of questions, depending on age of participant.