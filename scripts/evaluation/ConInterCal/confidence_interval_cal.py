import stat_util
import numpy as np
import matplotlib.pyplot as plt
from evaluation_metrics import rmse, pearson, spearman, ci, average_AUC, AUC
from sklearn.externals import joblib
from keras.models import load_model

class Model_Evaluation:
    def __init__(self, model, data_set):
        self.model = model
        self.data_set = data_set
    def model_load(self):
        if self.model == 'rfr':
            model_rfr = joblib.load('/Users/gvin/Desktop/kinasepkipred/models/Random_forest_gridsearch_py27.mdl')
            if self.data_set == 'test':
                x_test = np.load('/Users/gvin/pki_paper_work/features/x_test_rfr.npy')
                y_test_rfr = np.load('/Users/gvin/pki_paper_work/features/x_test_rfr.npy')
                y_pred = model_rfr.predict(x_test)
                return y_pred
            if self.data_set == 'external':
                data = np.load('/Users/gvin/pki_paper_work/features/total_metz_pki_propy.npy')
                x_test = data[:,:-1]
                y_test_ext = data[:,-1]
                y_pred = model_rfr.predict(x_test)
                return y_pred
        else:
            model_xgb = joblib.load('/Users/gvin/Desktop/kinasepkipred/models/Xgboost_rs.mdl')
            if self.data_set == 'test':
                x_test = np.load('/Users/gvin/pki_paper_work/features/pki_xgboost/x_test_xgb.npy')
                y_test_xgb = np.load('/Users/gvin/pki_paper_work/features/y_test_xgb.npy')
                y_pred = model_xgb.predict(x_test)
                return y_pred
            if self.data_set == 'external':
                data = np.load('/Users/gvin/pki_paper_work/features/total_metz_pki_propy.npy')
                x_test = data[:,:-1]
                y_test = data[:,-1]
                y_pred = model_xgb.predict(x_test)
                return y_pred

    def evaluate_metrics(self, y_pred, score_fun_used):
        y_test_rfr = np.load('/Users/gvin/pki_paper_work/features/x_test_rfr.npy')
        y_test_xgb = np.load('/Users/gvin/pki_paper_work/features/pki_xgboost/x_test_xgb.npy')
        data = np.load('/Users/gvin/pki_paper_work/features/total_metz_pki_propy.npy')
        x_test = data[:,:-1]
        y_test_ext = data[:,-1]

        if self.model=='rfr' and self.data_set == 'test':
            y_true = y_test_rfr
        elif self.model=='rfr' and self.data_set == 'external':
            y_true = y_test_ext

        elif self.model == 'xgb' and self.data_set == 'test':
            y_true = y_test_xgb
        elif self.model == 'xgb' and self.data_set == 'external':
            y_true = y_test_ext
        score, ci_lower, ci_upper, scores = stat_util.score_ci(y_true, y_pred,
                                                       score_fun=score_fun_used)

        print(name+" = {:.3f}, 95% CI: {:.3f}-{:.3f}".format(score, ci_lower, ci_upper))


if __name__ == '__main__':
    model_name=['rfr', 'xgb']
    data_sets = ['test', 'external']
    for m_name in model_name:
        for data_set in data_sets:

            mdl = Model_Evaluation(m_name, data_set)
            y_pred = mdl.model_load()
            if data_set == 'test':
                metrics = [rmse, pearson, spearman, average_AUC, AUC]
                name_metrics = ['rmse', 'pearson', 'spearman', 'average_AUC']
            else:
                metrics = [rmse, pearson, spearman, AUC]
                name_metrics = ['rmse', 'pearson', 'spearman', 'AUC']

            print()
            print('RESULTS FOR ' +mdl.model+' WITH '+mdl.data_set)
            for name, metric in zip(name_metrics, metrics):
                mdl.evaluate_metrics(y_pred, metric)

