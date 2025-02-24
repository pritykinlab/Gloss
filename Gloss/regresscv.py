## where you can figure out hyperparameters for the model and an initial pass for your data
# designed to be done over a single cell type

from .regressor import Regressor
from .prepdata import PrepData

from sklearn.model_selection import KFold, GridSearchCV, cross_validate
from sklearn.metrics import make_scorer, r2_score, mean_squared_error
import statistics
from scipy.stats import pearsonr
import numpy as np

class RegressCV():
    # in init, I need to pass in the regressions as 
    def __init__(self, adata, resolutions, pathway_string,
                 sample_hashtag = 'hash_max', libsize = 'n_counts',
                 interaction = 'biotin',
                 group_regs = [0.004, 0.005, 0.006, 0.007, 0.008, 0.02], single_gene_regs = [1., 2., 3., 4., 5.], 
                 donors_profiled=False):
        self.pathway_string = pathway_string
        self.prepped_data = PrepData(adata, self.pathway_string, 
                                     sample_hashtag, libsize,
                                     interaction,
                                     donors_profiled=donors_profiled)
        self.group_regs = group_regs 
        self.single_gene_regs = single_gene_regs
        self.resolutions = resolutions
        self.response = self._cv_loop()
        self.best_params = self._get_best_params()

    def _get_best_params(self):
        best_params = {}
        for res in self.resolutions:
            best_params[res] = {}
            for ctype in self.resolutions[res]:
                ctype_response = self.response[res][ctype]
                best_group_regs = [model.best_params_['group_reg'] for model in ctype_response['estimator']]
                best_single_gene_regs = [model.best_params_['single_gene_reg'] for model in ctype_response['estimator']]
                best_params[res][ctype] = { 'group_reg' : statistics.median(best_group_regs),
                                            'single_gene_reg' : statistics.median(best_single_gene_regs) }
        return best_params
        


    def _scoring(self):

        def pearson_corr(y, y_pred):
            return pearsonr(y, y_pred)[0]
        def pearson_sig(y, y_pred):
            return pearsonr(y, y_pred,)[1]
        
        pearson = make_scorer(pearson_corr)
        pearson_2 = make_scorer(pearson_sig)
        r2 = make_scorer(r2_score)
        mse = make_scorer(mean_squared_error)

        score_dict = {
            'pearson' : pearson,
            'pearson_sig' : pearson_2,
            'r2' : r2,
            'mse' : mse
        } 

        return score_dict 
    
    def _cv_loop(self):
        response = {}
        # Set up outer cross-validation (for model evaluation)
        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        # Set up inner cross-validation (for hyperparameter tuning)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        
        for res in self.resolutions:
            response[res] = {}
            for ctype in self.resolutions[res]:

                # need to throw in some data preparation logic here using PrepData
                X, y = self.prepped_data.subset(res, ctype)

                myreg = Regressor(gene_list = self.prepped_data.genes, pathway_string=self.pathway_string)
                param_grid = {
                    'group_reg' : self.group_regs,
                    'single_gene_reg' : self.single_gene_regs
                }
                score_dict = self._scoring()
                grid_search = GridSearchCV(myreg, param_grid, cv=inner_cv, n_jobs=-1, scoring='neg_mean_squared_error')
                response[res][ctype] = cross_validate(grid_search, X, y, cv=outer_cv, n_jobs=-1, 
                                                           return_train_score=True, return_estimator=True,
                                                           scoring=score_dict)
        
        return response