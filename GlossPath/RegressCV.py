## where you can figure out hyperparameters for the model and an initial pass for your data
# designed to be done over a single cell type

import Regressor
import PathwaysModule
from sklearn.model_selection import KFold
from sklearn.metrics import make_scorer, r2_score, mean_squared_error
import numpy as np
from scipy.stats import pearsonr
from skglm import GroupLasso

class RegressCV():
    def __init__(self):
        pass

    def _cv_loop(X, y, partitions_weights, alphas):
        
        # Set up outer cross-validation (for model evaluation)
        outer_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        # Set up inner cross-validation (for hyperparameter tuning)
        inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)
        
        # Loop over outer folds
        estimators = []

        test_mse_arr = []
        train_mse_arr = []
        train_pearson_arr = []
        test_pearson_arr = []
        train_r2_score_arr = []
        test_r2_score_arr = []

        for train_index, test_index in outer_cv.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            best_alpha = None
            best_mse = np.inf  # Keep track of the best score

            # Loop over possible alpha values (hyperparameter tuning)
            for alpha in alphas:
                print(alpha)
                inner_mse = []  # Store validation MSE for each inner fold

                for inner_train_index, inner_val_index in inner_cv.split(X_train):
                    X_inner_train, X_inner_val = X_train[inner_train_index], X_train[inner_val_index]
                    y_inner_train, y_inner_val = y_train[inner_train_index], y_train[inner_val_index]

                    # Initialize and fit the Group Lasso model with current alpha
                    model = GroupLasso(alpha=alpha, groups=partitions_weights[0], weights=partitions_weights[1])
                    model.fit(X_inner_train, y_inner_train)

                    # Predict on the validation set
                    y_pred_inner = model.predict(X_inner_val)
                    mse_inner = mean_squared_error(y_inner_val, y_pred_inner)
                    inner_mse.append(mse_inner)

                # Calculate the average MSE across the inner folds for this alpha
                avg_inner_mse = np.mean(inner_mse)

                # Check if this is the best alpha so far
                if avg_inner_mse < best_mse:
                    best_mse = avg_inner_mse
                    best_alpha = alpha

            # After hyperparameter tuning, train on the full training set using the best alpha
            best_model = GroupLasso(alpha=best_alpha, groups=partitions_weights[0], weights=partitions_weights[1])
            best_model.fit(X_train, y_train)
            estimators.append(best_model)

            # Evaluate on the outer test set
            y_pred_test = best_model.predict(X_test)
            test_mse_arr.append(mean_squared_error(y_test, y_pred_test))
            test_r2_score_arr.append(r2_score(y_test, y_pred_test))
            test_pearson_arr.append(pearsonr(y_test, y_pred_test)[0])

            # Evaluate on the inner train set
            y_pred_train = best_model.predict(X_train)
            train_mse_arr.append(mean_squared_error(y_train, y_pred_train))
            train_r2_score_arr.append(r2_score(y_train, y_pred_train))
            train_pearson_arr.append(pearsonr(y_train, y_pred_train)[0])   
        
        res = {
            'estimator' : estimators,
            'train_mse' : train_mse_arr,
            'train_pearson' : train_pearson_arr,
            'train_r2_score' : train_r2_score_arr,
            'test_mse' : test_mse_arr,
            'test_pearson' : test_pearson_arr,
            'test_r2_score' : test_r2_score_arr,
        }
        
        return res