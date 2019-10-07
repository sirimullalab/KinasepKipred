# script for Xgboost model development
sys.path.append("./evaluation/")
import numpy as np
import os
import xgboost as xgb
from evaluation_metrics import rmse, pearson, spearman, ci, average_AUC, AUC
from sklearn.model_selection import train_test_split, RandomizedSearchCV, GridSearchCV
from sklearn.externals import joblib

MODEL_DIR = './models'

if not os.path.isdir(MODEL_DIR): os.makedirs(MODEL_DIR)

def train_test(dataset):
    """
    This function is used to trains the model by grid search method and evaluates the test set
    
    Parameters
    ----------
    Features : np.array
        Protein and ligand features are concatenated and used as input file
    
    Output
    ------
        Saves the best model and also prints its performance in various evaluation metrics
    """
    
    data_set = np.load(dataset)
    X = data_set[:, :-1]
    y = data_set[:, -1]
    x_train, x_test, y_train, y_test = train_test_split(X, y, test_size = 0.25, random_state = 42)
    
#    x_train = np.load('x_train.npy')
#    y_train = np.load('y_train.npy')
#    x_test = np.load('x_test.npy')
#    y_test = np.load('y_test.npy')
    
    param_grid = {'n_estimators':[1000, 500, 100], 'objective':['reg:linear'],
                  'colsample_bytree':[0.3, 0.6, 1.0], 'learning_rate':[0.1, 0.001, 0.005, 1.0],\
                  'subsample':[0.8, 1.0], 'max_depth':[3, 5, 10 ], 'alpha':[0,10], 'gamma':[0, 1, 5]}
    xgbr = xgb.XGBRegressor()
#    xgbr = RandomizedSearchCV(estimator = xgbr, param_distributions = param_grid, n_iter = 10, cv = 5)
    xgbr = GridSearchCV(estimator=xgbr, param_grid=param_grid, cv= 5)
    
    print('XGBoost mdl fitting is started')
    
    xgbr=xgbr.fit(x_train, y_train)
    print(xgbr.best_params_)
    best_model = xgbr.best_estimator_
    model_name = MODEL_DIR+"/xgb.mdl"
    joblib.dump(best_model, model_name)
    y_pred = best_model.predict(x_test)
    
    print('XGBoost model is saved')
    
    PEARSON_R = pearson(y_test, y_pred)
    SPEARMAN_R = spearman(y_test, y_pred)
    RMSE = rmse(y_test, y_pred)
    Conc_Index = ci(y_test,y_pred)
    Avg_AUC = average_AUC(y_test, y_pred)
    
    print("PEARSON_R{:0.3f}: ".format(PEARSON_R))
    print("SPEARMAN_R{:0.3f}: ".format(SPEARMAN_R))
    print("RMSE{:0.3f}: ".format(RMSE))
    print("Conc_Index{:0.3f}".format(Conc_Index))
    print("Avg_AUC{:0.3f}".format(Avg_AUC))

def external_set(external_data):
    """
    This function is used to evaluate the model on external data set
    
    Parameters
    ----------
    Features : np.array
        Protein and ligand features are concatenated and used as input file
    
    Output
    ------
        Prints model performance  on external set in various evaluation metrics
    """
    ext_data = np.load(external_data)
    x_ext = ext_data[:, :-1]
    y_ext = ext_data[:,-1:].ravel()
    model_name = MODEL_DIR+"/xgb.mdl"
    model = joblib.load(model_name)
    y_pred_ext = model.predict(x_ext)
    
    print('external set is getting evaluated')
    
    PEARSON_R = pearson(y_ext, y_pred_ext)
    SPEARMAN_R = spearman(y_ext, y_pred_ext)
    RMSE = rmse(y_ext, y_pred_ext)
    Conc_Index = ci(y_ext,y_pred_ext)
    auc = AUC(y_ext, y_pred_ext)
    
    print("PEARSON_R {:0.3f}: ".format(PEARSON_R))
    print("SPEARMAN_R {:0.3f}: ".format(SPEARMAN_R))
    print("RMSE {:0.3f}: ".format(RMSE))
    print("Conc_Index {:0.3f}: ".format(Conc_Index))
    print("Avg_AUC {:0.3f}: ".format(auc))

if __name__ == "__main__":
    
    train_test_dataset = 'total_data_features.npy'
    external_data_set = 'total_metz_pki_propy.npy'
    
    train_test(train_test_dataset)
    external_set(external_data_set)
