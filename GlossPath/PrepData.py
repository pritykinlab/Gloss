## take in raw genes , raw sample hashtag , raw biotin, cell type annotations
## make it ready for downstream processing in other files
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

class PrepData():
    def __init__(self, adata, donors_profiled = False):
        ## we need to have raw counts in the anndata.X, 
        ## and we need to have 'biotin_raw' as a field
        ## and we need to have 'sample hashtag' out
        self.donors_profiled = donors_profiled
        ## what to do in the context of donors also being profiled?
        ## separate functions 
        if type(adata) == str:
            self.adata = sc.read_h5ad(adata)
        else:
            self.adata = adata

        self._normalize_adata()
        self._prep_confounders()
        self._prep_biotin()
    
    def _prep_confounders(self):

        ## assume hashtag is in a specific field
        hto_path = '/Genomics/pritykinlab/sarah/lipstic_analysis/ulipstic-analysis/gut.extra_sample_HTO_count_info.csv'
        hto_info = pd.read_csv(hto_path, index_col=0)

        self.adata.obs['avg_sample_hto'] = hto_info['sample_counts_avgd']

        # new confounder variables
        rna_logcounts = np.log10(self.adata.obs['total_counts'] + np.percentile(self.adata.obs['total_counts'], 0.1))
        hash_logcounts = np.log10(self.adata.obs['avg_sample_hto'] + np.percentile(self.adata.obs['avg_sample_hto'], 0.1))

        mean = np.mean(rna_logcounts)
        std_dev = np.std(rna_logcounts)
        # Standardize the array
        rna_logcounts = (rna_logcounts - mean) / std_dev

        mean = np.mean(hash_logcounts)
        std_dev = np.std(hash_logcounts)
        # Standardize the array
        hash_logcounts = (hash_logcounts - mean) / std_dev

        self.adata.obs['log_sample_hashtag'] = hash_logcounts
        self.adata.obs['log_RNA_libsize'] = rna_logcounts

    def _prep_biotin(self):
        new_biotin = np.log10(self.adata.obs['raw_biotin'] + np.percentile(self.adata.obs['raw_biotin'], 5))
        self.adata.obs['new_biotin'] = new_biotin - min(new_biotin)

        new_biotin_temp = self.adata.obs['raw_biotin'] / self.adata.obs['avg_sample_hto']
        new_biotin = np.log10(new_biotin_temp + np.percentile(new_biotin_temp, 5))
        self.adata.obs['other_new_biotin'] = new_biotin - min(new_biotin)

    def _normalize_adata(self):
        if not self.donors_profiled:
            self.adata.layers['raw_counts'] = self.adata.X.copy()
            sc.pp.normalize_total(self.adata, target_sum = 10000)
            sc.pp.log1p(self.adata)
            self.adata.layers["log_lib_norm"] = self.adata.X.copy()
            sc.pp.scale(self.adata)
            self.adata.layers["log_scaled"] = self.adata.X.copy()
        else:
            pass
