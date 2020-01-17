# Script for Artificial neural network
# *************************************************************
# Author: Govinda KC                                          #
# UTEP, Computational Science                                 #
# Last modified: 1/17/2020                                    #
# *************************************************************
from evaluation_metrics import rmse, pearson, spearman, ci, average_AUC, AUC
import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from sklearn.metrics import r2_score
from keras.models import load_model
from sklearn.model_selection import train_test_split
import os
MODEL_DIR = '../models'
if not os.path.isdir(MODEL_DIR): os.makedirs(MODEL_DIR)

def ann():    
    R2 = []
    RMSE = []
    PEARSON_R = []
    SPEARMAN_R = []
    Conc_Index = []
    Avg_AUC = []
    R2_ext = []
    RMSE_ext = []
    PEARSON_R_ext = []
    SPEARMAN_R_ext = []
    Conc_Index_ext = []
    AUC_ext = []
    
    for i in range(5):
        # Loading file for training and test
        data = np.load('total_data_features.npy')
        x = data[:,:-1]
        y = data[:,-1]
    
        x_train, x_test, y_train, y_test = train_test_split(x, y, test_size = 0.25, random_state = None)
        regressor = Sequential()
        regressor.add(Dense(32,activation = 'relu', kernel_initializer='normal', input_dim = 4741))
        regressor.add(Dense(64,activation = 'relu'))
        regressor.add(Dense(128,activation = 'relu'))
        regressor.add(Dense(128,activation = 'relu'))
        regressor.add(Dense(32, activation = 'relu'))
        regressor.add(Dense(1, activation = 'linear'))
        regressor.summary()
        
        # compiling the ANN
        regressor.compile(optimizer = 'adam', loss = 'mean_squared_error')
        regressor.fit(x_train, y_train, batch_size=25, epochs = 100, verbose=True)
        
        model_name = MODEL_DIR+'/ann.mdl'
        regressor.save(model_name)       

        model = load_model(model_name)
        # Evaluation of test dataset
        y_pred = model.predict(x_test)
        y = y_pred.ravel()
        f = y_test
        print('y', y)
        print('f', f)
        r2 = r2_score(y,f)
        R2.append(r2)
        Rmse = rmse(y,f)
        RMSE.append(Rmse)
        Pearson_r = pearson(y,f)
        PEARSON_R.append(Pearson_r)
        Spearman_r = spearman(y,f)
        SPEARMAN_R.append(Spearman_r)
        conc_index = ci(y,f)
        Conc_Index.append(conc_index)
        avg_auc = average_AUC(y,f)
        Avg_AUC.append(avg_auc)
        
        # Evaluation of external dataset
        ext_data = np.load('/Users/gvin/pki_paper_work/pki_github/total_metz_pki_propy.npy')
        ext_x = ext_data[:,:-1]
        ext_y = ext_data[:,-1:].ravel()
        y_pred = model.predict(ext_x)
        y = y_pred.ravel()
        f = ext_y
        r2_ext = r2_score(y,f)
        R2_ext.append(r2_ext)
        Rmse_ext = rmse(y,f)
        RMSE_ext.append(Rmse_ext)
        Pearson_r_ext = pearson(y,f)
        PEARSON_R_ext.append(Pearson_r_ext)
        Spearman_r_ext = spearman(y,f)
        SPEARMAN_R_ext.append(Spearman_r_ext)
        conc_index_ext = ci(y,f)
        Conc_Index_ext.append(conc_index_ext)
        auc_ext = AUC(y,f)
        AUC_ext.append(auc_ext)
        
    mean_R2 = sum(R2)/len(R2)
    mean_RMSE = sum(RMSE)/len(RMSE)
    mean_PEARSON_R = sum(PEARSON_R)/len(PEARSON_R)
    mean_SPEARMAN_R = sum(SPEARMAN_R)/len(SPEARMAN_R)
    mean_Conc_Index = sum(Conc_Index)/len(Conc_Index)
    mean_Avg_AUC = sum(Avg_AUC)/len(Avg_AUC)
        
    mean_R2_ext = sum(R2_ext)/len(R2_ext)
    mean_RMSE_ext = sum(RMSE_ext)/len(RMSE_ext)
    mean_PEARSON_R_ext= sum(PEARSON_R_ext)/len(PEARSON_R_ext)
    mean_SPEARMAN_R_ext = sum(SPEARMAN_R_ext)/len(SPEARMAN_R_ext)
    mean_Conc_Index_ext = sum(Conc_Index_ext)/len(Conc_Index_ext)
    mean_AUC_ext = sum(AUC_ext)/len(AUC_ext)
    print('Test set results')
    print(mean_R2, mean_RMSE, mean_PEARSON_R, mean_SPEARMAN_R, mean_Conc_Index, mean_Avg_AUC)
    print('External set results')
    print(mean_R2_ext, mean_RMSE_ext, mean_PEARSON_R_ext, mean_SPEARMAN_R_ext, mean_Conc_Index_ext,
           mean_AUC_ext)
    
if __name__ == "__main__":
    ann()
