import numpy as np
import joblib
from sklearn.model_selection import GridSearchCV
from sklearn.neural_network import MLPRegressor
import itertools
from evaluation_metrics import rmse, pearson, spearman, ci, average_AUC, AUC
#from hypopt import GridSearch

class Train:
    def __init__(self, x_train, y_train, x_test, y_test, externalData):
        self.x_train = x_train
        self.y_train = y_train
        self.x_test = x_test
        self.y_test = y_test
        self.x_ext = externalData[:,:-1]
        self.y_ext = externalData[:,-1]
    
    def mlp_regressor(self):    
        nodesList = [64, 128, 256]
        layers = []
        
        layers = list(itertools.permutations(nodesList))
        for i in range(1, len(nodesList)):
            layers.extend([p for p in itertools.product(nodesList, repeat=i+1)])
        
        regressor = MLPRegressor(solver='adam', alpha=1e-5, early_stopping=True)
        params = {'hidden_layer_sizes': layers}
        
        mdl = GridSearchCV(estimator=regressor, param_grid = params, cv = 5)
        
        print('Fitting ANN')
        mdl.fit(self.x_train, self.y_train)
        bestmodel = mdl.best_estimator_
        print(mdl.best_params_)
    
        joblib.dump(bestmodel, 'ann_mlp.mdl')
        #bestmodel = joblib.load('models/ann_mlp.mdl')     
        self.evaluateModel(bestmodel)

    def evaluateModel(self, mdl):

        y_pred = mdl.predict(self.x_test)

        r_test = pearson(self.y_test, y_pred)
        rho_test = spearman(self.y_test, y_pred)
        rmse_test = rmse(self.y_test, y_pred)
        ci_test = ci(self.y_test, y_pred)
        auc_test = average_AUC(self.y_test, y_pred)

        y_pred_ext = mdl.predict(self.x_ext)

        r_ext = pearson(self.y_ext, y_pred_ext)
        rho_ext = spearman(self.y_ext, y_pred_ext)
        rmse_ext = spearman(self.y_ext, y_pred_ext)
        ci_ext = ci(self.y_ext, y_pred_ext)
        auc_ext = AUC(self.y_ext, y_pred_ext)

        print('Test Set Results')
        print(f'r_test: {r_test:.3f}, rho_test: {rho_test:.3f}, rmse_test: {rmse_test:.3f}, \
                ci_test: {ci_test:.3f} auc_test: {auc_test:.3f}')
        
        print('Metz dataSet results')
        print(f'r_ext: {r_ext:.3f}, rho_ext: {rho_ext:.3f}, rmse_ext {rmse_ext:.3f}, \
                ci_ext: {ci_ext:.3f}, auc_ext: {auc_ext:.3f}')
        
        

if __name__ == '__main__':
    x_train = np.load('x_train.npy')
    y_train = np.load('y_train.npy')
    x_test = np.load('x_test.npy')
    y_test = np.load('y_test.npy')

    externalData = np.load('total_metz_pki_propy.npy')
    
    train = Train(x_train, y_train, x_test, y_test, externalData)
    train.mlp_regressor()
