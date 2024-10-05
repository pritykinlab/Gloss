# where you can build a distribution of coefficients for your data, for each regression
# designed to be done within a single cell type

import Regressor

class RegressBootstrap():
    def __init__(self,
                 adata, regressor,
                 n_bootstraps=100,
                 ):
        self.n_bootstraps = n_bootstraps
        self.adata = adata
        self.regressor = regressor