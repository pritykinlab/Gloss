# users input an anndata, and they get back biotin for that anndata that's normalized, useful for vis, etc
# regressor (that is trained) has a function that allows the computation of the normalized biotin for a data point

from skglm import GroupLasso
import pandas as pd 
import numpy as np
import scanpy as sc
import statistics

from sklearn.linear_model import Lasso, LassoCV
from sklearn.model_selection import cross_validate, KFold

class NormalizeLipsticRegressor():
    def __init__(self):
        self.r = 42
        self.inner_cv = KFold(n_splits=5, shuffle=True, random_state=self.r)
        self.outer_cv = KFold(n_splits=5, shuffle=True, random_state=self.r) 
        self.model = LassoCV(cv=self.inner_cv)

    def normalize_fit(self, features_matrix, y):
        lasso_res = cross_validate(self.model, X=features_matrix, y=y, cv=self.outer_cv,
                return_train_score=True,
                return_estimator=True,
                n_jobs=-1
                )
        
        best_alphas = [model.alpha_ for model in lasso_res['estimator']]
        optimal_lasso_param = statistics.median(best_alphas)

        return Lasso(alpha = optimal_lasso_param).fit(features_matrix, y)


            