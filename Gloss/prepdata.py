## take in raw genes , raw sample hashtag , raw biotin, cell type annotations
## make it ready for downstream processing in other files
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

from .pathwaysmodule import PathwaysModule
from .normalizelipstic import NormalizeLipsticRegressor

class PrepData():
    def __init__(self, adata, pathway_string, 
                 sample_hashtag = 'hash_max', libsize = 'n_counts',
                 interaction = 'biotin',
                 donors_profiled = False, only_pathways=False):
        ## we need to have raw counts in the anndata.X, 
        ## and we need to have 'biotin_raw' as a field
        ## and we need to have 'sample hashtag' out
        self.norm = 'log_scaled'
        self.donors_profiled = donors_profiled
        self.only_pathways = only_pathways
        self.sample_hashtag = sample_hashtag
        self.libsize = libsize 
        self.interaction = interaction

        if type(adata) == str:
            self.adata = sc.read_h5ad(adata)
        else:
            self.adata = adata

        self._normalize_adata()
        self._prep_confounders()
        self._prep_biotin()

        self.pathway_string = pathway_string
        self.pathways_mod = self._prep_pathways()

        self.genes = list(self.adata.var.index)

        self.X = self._prep_X()
        self.y = self._prep_y()

        self._normalize_biotin()

    def subset(self, resolution, celltype):
        index = self.adata.obs[resolution].isin([celltype])
        y = self.y[index]
        X = self.X.to_numpy()[index, :]
        return X, y
    
    def _prep_pathways(self):
        return PathwaysModule(self.pathway_string, list(self.adata.var.index))

    def _prep_confounders(self):
        # new confounder variables
        self.adata.obs['log_RNA_libsize'] = self._normalize_confounder(self.libsize)
        self.adata.obs['log_sample_hashtag'] = self._normalize_confounder(self.sample_hashtag)

    def _normalize_confounder(self, confounder_str):
        confounder = np.log10(self.adata.obs[confounder_str] + np.percentile(self.adata.obs[confounder_str], 0.1))
        mean = np.mean(confounder)
        std_dev = np.std(confounder)
        # Standardize the array
        return (confounder - mean) / std_dev

    def _prep_biotin(self):
        new_biotin = np.log10(self.adata.obs[self.interaction] + np.percentile(self.adata.obs[self.interaction], 5))
        self.adata.obs['new_biotin'] = new_biotin - min(new_biotin)

        new_biotin_temp = self.adata.obs[self.interaction] / self.adata.obs[self.sample_hashtag]
        new_biotin = np.log10(new_biotin_temp + np.percentile(new_biotin_temp, 5))
        self.adata.obs['other_new_biotin'] = new_biotin - min(new_biotin)

    def _normalize_adata(self):
        self.adata.layers['raw_counts'] = self.adata.X.copy()
        sc.pp.filter_genes(self.adata, min_counts=1)
        sc.pp.normalize_total(self.adata, target_sum = 10000)
        sc.pp.log1p(self.adata)
        self.adata.layers["log_lib_norm"] = self.adata.X.copy()
        sc.pp.scale(self.adata)
        self.adata.layers["log_scaled"] = self.adata.X.copy()
        if self.donors_profiled:
            self.adata = self.adata[~self.adata.obs['donor']].copy()

    def _prep_X(self):
        self.adata.X = self.adata.layers[self.norm].copy()
        ad_df = self.adata.to_df()
        new_df = pd.DataFrame({}, index= ad_df.index)
        mydf_list = []
        for pathway in self.pathways_mod.pathways:
            mydf = ad_df.loc[:, self.pathways_mod.pathways[pathway]]
            mydf.columns = mydf.columns + '_' + pathway
            mydf_list.append(mydf)
        
        if not self.only_pathways:
            #nogroup_df = ad_df[[col for col in ad_df.columns if col not in self.pathways_mod.pathway_genes]]
            nogroup_df = ad_df.copy()
            nogroup_df.columns = nogroup_df.columns + '_no_pathway'
            mydf_list.append(nogroup_df)
        
        new_df = pd.concat(mydf_list, axis=1)
        new_df = new_df.copy()
        
        new_df['log_sample_hashtag'] = self.adata.obs['log_sample_hashtag']
        new_df['log_RNA_libsize'] = self.adata.obs['log_RNA_libsize']

        return new_df
        
    def _prep_y(self):
        return self.adata.obs['new_biotin']
    
    def _normalize_biotin(self):
        rna_logcounts = self.adata.obs['log_RNA_libsize']
        hash_logcounts = self.adata.obs['log_sample_hashtag']

        features_matrix = np.concatenate((
                                        rna_logcounts.to_numpy().reshape(-1,1),
                                        hash_logcounts.to_numpy().reshape(-1,1),
                                        ), axis=1)
        y = self.y
        
        norm_reg = NormalizeLipsticRegressor().normalize_fit(features_matrix, y)

        my_libsize = self.adata.obs['log_RNA_libsize']
        my_sample_hashtag = self.adata.obs['log_sample_hashtag']
        my_biotin = self.adata.obs['new_biotin']

        self.adata.obs['gloss_normalized_biotin'] = my_biotin - (norm_reg.coef_[-2] * my_libsize + norm_reg.coef_[-1] * my_sample_hashtag)
