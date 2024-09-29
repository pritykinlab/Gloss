from skglm import GroupLasso
import pandas as pd 
import numpy as np
import scanpy as sc

from .pathways import PathwaysModule

# should there be a function that automatically does the CV loop? Probably, right?

class Regressor():
    # separate pathway database 
    def __init__(self, adata, pathway_string, group_reg, l1_reg):
        self.norm = 'log_lib_norm'

        self.adata = adata
        self.pathway_string = pathway_string
        self.pathways = self.prep_pathways()
        self.X = self._prep_X()
        self.y = self._prep_y()

    def _prep_X(self):
        self.adata.X = self.adata.layers[self.norm].copy()
        sc.pp.scale(self.adata)
        ad_df = self.adata.to_df()

        new_df = pd.DataFrame({}, index= ad_df.index)
        mydf_list = []
        for i, pathway in enumerate(self.pathways):
            mydf = ad_df.loc[:, self.pathways[pathway]]
            mydf.columns = mydf.columns + '_' + pathway
            mydf_list.append(mydf)
        new_df = pd.concat(mydf_list, axis=1)
        new_df = new_df.copy()
        
        new_df['log_sample_hashtag'] = self.adata.obs['log_sample_hashtag']
        new_df['log_RNA_libsize'] = self.adata.obs['log_RNA_libsize']

        return new_df
        
    def _prep_y(self):
        return self.adata.obs['new_biotin']
        
    def _prep_pathways(self):
        return PathwaysModule(self.pathway_string)

    def fit(self):
        reg = GroupLasso()

    def get_partitions_and_weights(self, path_dict):
        partitions = []
        i = 0
        for pathway in path_dict:
            partitions.append(list(np.arange(i, i+len(path_dict[pathway]))))
            i += len(path_dict[pathway])
        # accounting for the two confounders
        for j in range(2):
            partitions.append(list(np.arange(i, i+1)))
            i += 1
        
        weights = []
        for part in partitions:
            weights.append(np.sqrt(len(part)))

        weights = np.array(weights)
        
        return partitions, weights

