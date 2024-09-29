from skglm import GroupLasso
import pandas as pd 
import numpy as np

from .pathways import PathwaysModule

# should there be a function that automatically does the CV loop? Probably, right?

class Regressor():
    # separate pathway database 
    def __init__(self, pathway_string):
        self.adata = None
        self.pathway_string = pathway_string
        self.pathways = {}

    def prep_data(self, adata):
        self.adata = adata 
        
    def prep_pathways(self):
        self.pathways = PathwaysModule(self.pathway_string)

    def fit():
        pass 

    def _create_duplicated_feature_matrix(path_dict, gene_ad, norm='log_lib_norm'):
        gene_ad.X = gene_ad.layers[norm].copy()
        sc.pp.scale(gene_ad)
        ad_df = gene_ad.to_df()
        new_df = pd.DataFrame({}, index= ad_df.index)
        mydf_list = []
        for i, pathway in enumerate(path_dict):
            if i % 1000 == 0:
                print(i)
            mydf = ad_df.loc[:, path_dict[pathway]]
            mydf.columns = mydf.columns + '_' + pathway
            mydf_list.append(mydf)
        new_df = pd.concat(mydf_list, axis=1)
        new_df = new_df.copy()
        
        new_df['log_sample_hashtag'] = gene_ad.obs['log_sample_hashtag']
        new_df['log_RNA_libsize'] = gene_ad.obs['log_RNA_libsize']

        return new_df

    def get_partitions_and_weights(path_dict):
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

