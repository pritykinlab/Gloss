# where you can build a distribution of coefficients for your data, for each regression
# designed to be done within a single cell type

from .regressor import Regressor
from .prepdata import PrepData

import numpy as np
import pandas as pd
from sklearn.utils import resample

class RegressBootstrap():
    def __init__(self, adata, resolutions, pathway_string, best_params, n_bootstraps = 100,
                 sample_hashtag = 'hash_max', libsize = 'n_counts',
                 interaction = 'biotin', 
                 donors_profiled=False):
        self.resolutions = resolutions
        self.best_params = best_params
        self.n_bootstraps = n_bootstraps
        self.pathway_string = pathway_string
        self.prepped_data = PrepData(adata, self.pathway_string, 
                                     sample_hashtag, libsize,
                                     interaction,
                                     donors_profiled=donors_profiled)
    
        self.response, self.intercepts = self._bootstrap()
        self.response_individual_genes = self._bootstrap_res_individual_genes()
        self.pathway_frequencies = self._get_pathway_frequencies()

    def _bootstrap(self):
        bootstrap_response = {}
        bootstrap_response_intercepts = {}
        for res in self.resolutions:
            bootstrap_response[res] = {}
            bootstrap_response_intercepts[res] = {}
            for ctype in self.resolutions[res]:
                X, y = self.prepped_data.subset(res, ctype)

                # List to store bootstrap coefficients
                bootstrap_coefficients = []
                bootstrap_intercepts = []

                # Bootstrapping
                for i in range(self.n_bootstraps):
                    X_resampled, y_resampled = resample(X, y, random_state=i)
                    myreg = Regressor(gene_list = self.prepped_data.genes, 
                                      pathway_string=self.pathway_string,
                                      group_reg=self.best_params[res][ctype]['group_reg'],
                                      single_gene_reg=self.best_params[res][ctype]['single_gene_reg']
                                      )
                    myreg.fit(X_resampled, y_resampled)
                    bootstrap_coefficients.append(myreg.coef_)
                    bootstrap_intercepts.append(myreg.intercept_)

                # Convert to DataFrame for analysis
                bootstrap_coefficients = np.array(bootstrap_coefficients)
                coefficients_df = pd.DataFrame(bootstrap_coefficients)
                coefficients_df.columns = self.prepped_data.X.columns

                # Sort the columns by median value
                medians = coefficients_df.median()
                sorted_columns = medians.sort_values().index
                sorted_df = coefficients_df[sorted_columns]
                bootstrap_response[res][ctype] = sorted_df
                bootstrap_response_intercepts[res][ctype] = bootstrap_intercepts

        return bootstrap_response, bootstrap_response_intercepts

    def _bootstrap_res_individual_genes(self):
        ind_genes_response = {}
        for res in self.resolutions:
            ind_genes_response[res] = {}
            for ctype in self.resolutions[res]:
                mydf = self.response[res][ctype].copy()
                
                grouped_data = {}
                for col in mydf:
                    # Remove replicate number (anything after '_') from the column name
                    base_feature = col.split('_')[0]
                    
                    # Sum the columns with the same base feature name
                    if base_feature in grouped_data:
                        grouped_data[base_feature] = grouped_data[base_feature] + mydf[col].copy()
                    else:
                        grouped_data[base_feature] = mydf[col].copy()
                
                # Create a new DataFrame with the summed coefficients
                summed_df = pd.DataFrame(grouped_data)
                medians = summed_df.median()
                sorted_columns = medians.sort_values().index
                summed_df = summed_df[sorted_columns]

                ind_genes_response[res][ctype] = summed_df.copy()
            
        return ind_genes_response
    
    def _get_pathway_frequencies(self):
        pathway_frequencies = {}
        for res in self.resolutions:
            pathway_frequencies[res] = {}
            for ctype in self.resolutions[res]:
                included_counts = {}
                mydf = self.response[res][ctype].copy()
                mydf = mydf.loc[:, ~mydf.columns.str.contains('_no_pathway')]
                for i in range(self.n_bootstraps):
                    test_boot_coefs = mydf.iloc[i]
                    myarray = test_boot_coefs[test_boot_coefs != 0]
                    myset = set([ pathway_names.split("_", 1)[1] for pathway_names in myarray.index  ])
                    for item in myset:
                        included_counts[item] = included_counts.get(item, 0) + 1
                my_dict = included_counts.copy()
                my_dict.pop('RNA_libsize', 0)
                my_dict.pop('sample_hashtag', 0)
                sorted_dict = sorted(my_dict.items(), key=lambda item: item[1], reverse=True)
                pathway_frequencies[res][ctype] = sorted_dict

        return pathway_frequencies

