# designed to be done over a single cell type, using normalized data

from skglm import GroupLasso
import pandas as pd 
import numpy as np
import scanpy as sc

import PathwaysModule

# should there be a function that automatically does the CV loop? Probably, right?
# so this class will just load and run the regression for a given data slice, over a specific subtype.

class Regressor():
    # separate pathway database 
    def __init__(self, adata, pathway_string, group_reg=0.005, l1_reg=3, only_pathways=False):
        self.group_reg = group_reg 
        self.l1_reg = l1_reg
        self.only_pathways = only_pathways
        self.norm = 'log_lib_norm'
        self.adata = adata
        self.adata_genes = set(self.adata.var.index)
        self.pathway_string = pathway_string
        self.pathways = self._prep_pathways().pathways
        self.pathway_genes = set([ x for mylist in self.pathways.values() for x in mylist ])
        self.X = self._prep_X()
        self.y = self._prep_y()
        self.p_w = self._get_partitions_and_weights()
        self.reg = GroupLasso(self.p_w[0], self.group_reg, self.p_w[1])

    def fit(self):
        self.reg.fit(self.X, self.y)
    
    def _prep_X(self):
        self.adata.X = self.adata.layers[self.norm].copy()
        sc.pp.scale(self.adata)
        ad_df = self.adata.to_df()
        new_df = pd.DataFrame({}, index= ad_df.index)
        mydf_list = []
        for pathway in enumerate(self.pathways):
            mydf = ad_df.loc[:, self.pathways[pathway]]
            mydf.columns = mydf.columns + '_' + pathway
            mydf_list.append(mydf)
        
        if not self.only_pathways:
            nogroup_df = ad_df[[col for col in ad_df.columns if col not in self.pathway_genes]]
            nogroup_df.columns = nogroup_df.columns + '_no_pathway'
            mydf_list.append(nogroup_df)
        
        new_df = pd.concat(mydf_list, axis=1)
        new_df = new_df.copy()
        
        new_df['log_sample_hashtag'] = self.adata.obs['log_sample_hashtag']
        new_df['log_RNA_libsize'] = self.adata.obs['log_RNA_libsize']

        return new_df
        
    def _prep_y(self):
        return self.adata.obs['new_biotin']
        
    def _prep_pathways(self):
        return PathwaysModule(self.pathway_string, self.adata_genes)

    def _get_partitions_and_weights(self):
        nogroup_cols = self.adata.var.index[~self.adata.var.index.isin(self.pathway_genes)]
        partitions = []
        weights = []
        i = 0

        ## taking care of pathway-included groups
        for pathway in self.pathways:
            partitions.append(list(np.arange(i, i+len(self.pathways[pathway]))))
            i += len(self.pathways[pathway])
        for part in partitions:
            weights.append(np.sqrt(len(part)))

        # accounting for the non-pathway-included groups
        if not self.only_pathways:
            k = 0
            for j in range(len(nogroup_cols)):
                partitions.append(list(np.arange(i, i+1)))
                i += 1
                k += 1
            for _ in range(k):
                weights.append(self.l1_reg)
        
        # accounting for the two confounders
        for _ in range(2):
            partitions.append(list(np.arange(i, i+1)))
            weights.append(1)
            i += 1

        weights = np.array(weights)
        
        return partitions, weights
