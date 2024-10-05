## where you can figure out hyperparameters for the model and an initial pass for your data
# designed to be done over a single cell type

from .regressor import Regressor
from .prepdata import PrepData

from sklearn.model_selection import KFold, GridSearchCV, cross_validate
from sklearn.metrics import make_scorer, r2_score, mean_squared_error
import numpy as np

class RegressCV():
    # in init, I need to pass in the regressions as 
    def __init__(self, adata, resolutions, pathway_string,
                 sample_hashtag = 'hash_max', libsize = 'n_counts',
                 interaction = 'biotin',
                 alphas = [0.004, 0.005, 0.006, 0.007], l1s = [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5], 
                 donors_profiled=False):
        self.pathway_string = pathway_string
        self.prepped_data = PrepData(adata, self.pathway_string, 
                                     sample_hashtag, libsize,
                                     interaction,
                                     donors_profiled=donors_profiled)
        self.alphas = alphas 
        self.l1s = l1s
        self.resolutions = resolutions
        self.response = {}
        self._cv_loop()

    def _cv_loop(self):
        
        # Set up outer cross-validation (for model evaluation)
        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        # Set up inner cross-validation (for hyperparameter tuning)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        
        for res in self.resolutions:
            self.response[res] = {}
            for ctype in self.resolutions[res]:

                # need to throw in some data preparation logic here using PrepData
                X, y = self.prepped_data.subset(res, ctype)

                myreg = Regressor(gene_list = self.prepped_data.genes, pathway_string=self.pathway_string)
                param_grid = {
                    'group_reg' : self.alphas,
                    'l1_reg' : self.l1s
                }
                grid_search = GridSearchCV(myreg, param_grid, cv=inner_cv, n_jobs=-1)
                self.response[res][ctype] = cross_validate(grid_search, X, y, cv=outer_cv, n_jobs=-1, 
                                                           return_train_score=True, return_estimator=True)
        
        return self