# designed to be done over a single cell type, using normalized data
from sklearn.base import BaseEstimator, RegressorMixin
from skglm import GroupLasso
import pandas as pd 
import numpy as np
import scanpy as sc

from .pathwaysmodule import PathwaysModule

# should there be a function that automatically does the CV loop? Probably, right?
# so this class will just load and run the regression for a given data slice, over a specific subtype.

# the weights are dependent on the l1_penalty fed in
# I guess then this regressor can be initialized based off of just feeding in the pathway string, right?
# not quite - I need to feed in a geneset that defines the 'scope' of the thing

class Regressor(BaseEstimator, RegressorMixin):
    # separate pathway database 
    def __init__(self, gene_list, pathway_string, group_reg=0.005, l1_reg=3, only_pathways=False):
        super().__init__()
        self.group_reg = group_reg 
        self.l1_reg = l1_reg
        self.gene_list = gene_list
        self.only_pathways = only_pathways
        self.pathway_string = pathway_string
        self.pathways = self._prep_pathways().pathways
        self.pathway_genes = set([ x for mylist in self.pathways.values() for x in mylist ])
        self.p_w = self._get_partitions_and_weights()
        self.reg = GroupLasso(self.p_w[0], self.group_reg, self.p_w[1])

    def fit(self, X, y):
        self.reg.fit(X, y)
        return self

    def predict(self, X):
        return self.reg.predict(X)
        
    def _prep_pathways(self):
        return PathwaysModule(self.pathway_string, self.gene_list)

    def _get_partitions_and_weights(self):
        nogroup_cols = [ gene for gene in self.gene_list if gene not in self.pathway_genes ]
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
