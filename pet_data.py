from glob import glob
import numpy as np
from pandas import read_excel, read_csv, DataFrame
from pandas.core.common import flatten
from pyradigm import ClassificationDataset as ClfDataset
from pyradigm import RegressionDataset as RegrDataset
from pyradigm.multiple import MultiDatasetRegress
from sspipe import p, px
import re

# Qs for pyradigm
#  * multiple ids (sub, ses)? just merge into a single string?
#      maybe `MultiDatasetRegress`?
#  * must choose prediction variable?
#  * can features be named?
#      yes - feature_names=...
#  * pandas to pyradigm?
#      no(?) - using `for i,r in df.iterrows(): ds.add_samplet(...)`


class ROI_FC:
    """used exclusively by `pyradigms_from_adj`
    ROI-ROI pair info for symmetric adjacency matrices
    and file globing
    """
    def __init__(self, fileglob, roi_names):
        self.fileglob = fileglob
        self.names = roi_names
        # roi-roi info
        self.n_roi = len(roi_names) # 23
        # in a NxN matrix, where are the unique pairs (exclude diagonal)
        idx = np.triu_indices(self.n_roi, 1) 
        self.n_feats = len(idx[0]) # 253 == scipy.special.comb(23,2)
        self.feat_names = [f'{roi_names[idx[0][i]]}-{roi_names[idx[1][i]]}'
                      for i in range(self.n_feats)]
        self.tri_idx = idx

    def all_files(self):
        return glob(self.fileglob)

    def extract_feats(self, adj_fpath):
        """ 253 roi-roi connectivity features per file (20201208) """
        # m is actually 24x24 instead of 23x23. last row and last col are all NA
        # NAs mean we need to use `genfromtxt` to read in
        # extra columns will never actually be seen.
        #  we index on the uper tri for 0-22 (discard 24th row and col)
        m = np.genfromtxt(adj_fpath, skip_header=True, missing_values="NA")
        v = m[self.tri_idx[0], self.tri_idx[1]]
        assert len(v) == self.n_feats
        return v

class PET:
    """tools around a central very wide csv file
    conveniences:
      * lookup session id OR age
      * pyradigm datasets from subset of metrics
        (e.g. only 'uppsp' related measures)
    """
    def __init__(self, csv='data/wide.csv'):
        # very very wide csv file. 
        # also has id and age lookup
        self.widedf = read_csv(csv)
        # widedf.shape # previosly (384, 6986). now (384, 547)

        # add sessid
        # will eventaully be added to merged_data.csv from R code.
        # this will become dead code
        if 'sessid' not in self.widedf.columns:
            self.widedf['sesid'] = [f"{x.lunaid}_{int(np.nan_to_num(x.visitnum))}" 
                                    for i, x in self.widedf.iterrows()]
        # quick dict lookups to
        # * go from MR ID to session a
        # * get age from sesid
        self.ses_lookup = {f"{x.lunaid}_{str(int(np.nan_to_num(x.vdate)))}": x.sesid
                           for i, x in self.widedf.iterrows()}
        self.age_lookup = {x.sesid: x.age
                           for i, x in self.widedf.iterrows()}

    def col_subset(self, patt, rmprefix=None):
        """ return a subset of the columns using a pattern 
        >>> d = PET()
        >>> len(d.col_subset('uppsp_.*tot').colums) # includes sesid, age, and uppsp_upps_tot
        3
        """
        cols = self.widedf.columns
        want = [bool(re.search('^sesid$|^age$|'+patt, x)) for x in cols]
        subset = self.widedf.iloc[:, want]
        if rmprefix is None:
            # assume
            # * the best prefix to remove is from the first non-{age,id} column 
            # * prefix is any text before the first '_'
            subset_specifc_cols = [x for x in cols[want] if x not in ['sesid','age']]
            rmprefix = subset_specifc_cols[0].split("_")[0]
        if rmprefix:
            subset.columns = [re.sub(f'^{rmprefix}_','',x) for x in subset.columns]
        return subset.dropna()
    
    def add_subset_to(self, ds, pattern, rmprefix=None):
        """
        add each row of the subsetted main dataframe to a pyradigm dataset
        
        widedf subsets always have sesid and age. use those as samplete_id and target
        the rest of dataframe columns get thrown in as features

        Parameters
        ----------
        ds : pyradigm dataset
          dataset to add samplets to
        pattern : str (regex)
          pattern for columns to keep (e.g. '^upps_*', 'upps_.*tot|upps_negurg')
        rmprefix : str
          prefix to remove from column names. defaults last columns first component (e.g. 'upps_')
          this reduces some redundency in names
         
        """
        d = self.col_subset(pattern, rmprefix)
        for (i, r) in d.iterrows():
            feats = r.drop(['sesid','age'])
            v = feats.T.tolist()
            feat_names = feats.axes[0].tolist()
            ds.add_samplet(samplet_id=r.sesid, target=r.age,
                           features=v, feature_names=feat_names,
                           overwrite=True) 


    def ld8_extract(self, text):
        """ extract id from string like 12345_20191231 """
        return re.search('\d{5}_\d{8}', text).group(0)

    def sesid(self, ld8):
        """
        make luna_visit from luna_date
        """
        return self.ses_lookup.get(ld8, '')

    def file_info(self, f):
        """ session id and age from a file """
        ld8 = self.ld8_extract(f) # get luna_date
        sid = self.sesid(ld8)     # make luna_visitnum
        age = self.age_lookup.get(sid)
        return (sid, age)
    
    def all_adj_into(self, ds, roi, outcsv=None):
        """
        functional connectivity adjecency matrices (scaled(?), center diag == 15)
        into pyradigm structure

        Parameters
        ----------
        ds : pyradigm
            dataset to add samplets to
        roi : ROI_FC
            ROI_FC object pointing to adj matrices and roi labels
        outcsv : str (path)
           write out to specified csv file
        """

        # also write to csv file
        outcsv = open(outcsv,'w')
        outline = ",".join(["sesid","age", *roi.names])
        outcsv.write(outline+"\n")

        for f in roi.all_files():
            #print(f)
            (sid, age) = self.file_info(f)
            v = roi.extract_feats(f)

            ds.add_samplet(samplet_id=sid, target=age,
                        features=v, feature_names=roi.feat_names,
                        overwrite=True) 
            outcsv.write(",".join([str(x) for x in [sid,age,*v]])+"\n")

        outcsv.close()

def pyradigms_from_adj(main_data):
    """
    This is only useful when run on LNCD server (access to Hera and Zeus)

    generate pyradigm for each connectivity set: task background, rest1, and rest2.
      background - fixiation extracted from 6 MR tasks. correlation of chunked timeseries
      rest1      - first (before task) do nothing task. correlation of continious timeseries
      rest2      - second (after task) do nothing task. correlation of continious timeseries
    
    23-23 ROI pairs each have a correlation value (then fisher z transformed)

    
    (be kind, first application of pyradigm)
    
    Parameters
    ----------
    main_data : PET
      object that knows how to get session ids and age from MR visit luna_date
    
    Return
    ----------
    dictionary of pyradigm datasets
   
    """
    # names are not in files. Hard code here instead
    roi_names = [
        'L_ant_med_OFC_1', 'L_post_med_OFC_1', 'L_ant_VMPFC', 'L_post_med_OFC_2', 'L_ant_med_OFC_2',
        'L_vent_ACC', 'L_subg_cingulate', 'L_rostr_ACC', 'R_ant_med_OFC_1', 'R_post_med_OFC_1', 'R_ant_VMPFC',
        'R_post_med_OFC_2','R_ant_med_OFC_2','R_vent_ACC', 'R_subg_cingulate', 'R_rostr_ACC', 'L_caudate',
        'L_NAC', 'L_putamen', 'R_caudate', 'R_NAC', 'R_putamen', 'VTA']

    # output of ROI_TempCor.R nested inside visit directories
    prefix="vmpfc_striatal_vta_mean_gsr_pearson_adj_pearson.txt"

    # each connectivity set has their own path to the adjacency matrix
    fc_globs = {'bkgrd': f'/Volumes/Hera/Projects/mMR_PETDA/subjs/1*_2*/func/*{prefix}',
                'rest1':  f'/Volumes/Zeus/preproc/petrest_rac1/MHRest_FM_ica/1*_2*/*{prefix}',
                'rest2':  f'/Volumes/Zeus/preproc/petrest_rac2/MHRest_FM_ica/1*_2*/*{prefix}'}

    # 'Regr' b/c age is continuous.
    # allow_nan_inf might be a mistake. easy for now (20201208)
    FCds = {k: RegrDataset(allow_nan_inf=True) for k in fc_globs.keys()}
    FCds['bkgrd'].description = "Functional connectivity w/23 ROIs from task chunked background timeseries"
    FCds['rest1'].description = "Functional connectivity w/23 ROIs from first rest timeseries"
    FCds['rest2'].description = "Functional connectivity w/23 ROIs from second rest timeseries"

    # 253 roi-roi pairs for each of background, rest1, and rest1
    # use info stored in the super wide dataframe to get age and session id from each file
    # to create each samplet that composes each (bkgrd,rest1,2) of the connectivity datasets
    for k in fc_globs.keys():
        roi = ROI_FC(fc_globs[k], roi_names)
        main_data.all_adj_into(FCds[k], roi, f'data/conn_adj_{k}.csv')
    
    return FCds

if __name__ == "__main__":
    main_data = PET('data/wide.csv')
    # functional connectivity ROIs are the same across background and rest
    # but the location of the adj files are not

    FCds = pyradigms_from_adj(main_data)

    upps = RegrDataset()
    upps.description = "Urgency, Premeditation (lack of), Perseverance (lack of), Sensation Seeking, Positive Urgency, Impulsive Behavior Scale"
    main_data.add_subset_to(upps, '^uppsp_')

    
    # error (20201208)
    # # > Differing set of IDs in two datasets. Unable to add this dataset to the MultiDataset.  

    # multi_ds = MultiDatasetRegress()
    # # bkgrd, rest1, rest2
    # for k in FCds.keys():
    #     multi_ds.append(FCds[k], k)
    # multi_ds.append(upps, 'upps')


