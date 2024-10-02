## where you can figure out hyperparameters for the model and an initial pass for your data
# designed to be done over a single cell type

import Regressor
import PathwaysModule
import PrepData

from sklearn.model_selection import KFold, GridSearchCV, cross_validate
from sklearn.metrics import make_scorer, r2_score, mean_squared_error
import numpy as np
from scipy.stats import pearsonr
from skglm import GroupLasso

class RegressCV():
    # in init, I need to pass in the regressions as 
    def __init__(self, adata, resolutions, pathway_string,
                 alphas = [0.004, 0.005, 0.006, 0.007], l1s = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]):
        self.adata = adata
        self.alphas = alphas 
        self.l1s = l1s
        self.resolutions = resolutions
        self.pathway_string = pathway_string

    def _cv_loop(self, X, y, partitions_weights, alphas):
        
        # Set up outer cross-validation (for model evaluation)
        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        # Set up inner cross-validation (for hyperparameter tuning)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        
        for res in self.resolutions:
            for ctype in self.resolutions[res]:

                # need to throw in some data preparation logic here using PrepData

                myreg = Regressor(gene_list = self.adata.var.index, pathway_string='hallmark')
                param_grid = {
                    'group_reg' : self.alphas,
                    'l1_reg' : self.l1s
                }
                grid_search = GridSearchCV(myreg, param_grid, cv=inner_cv, n_jobs=-1)
                response = cross_validate(grid_search, X, y, cv=outer_cv, n_jobs=-1)

                # what to return, what to return
        
        return res