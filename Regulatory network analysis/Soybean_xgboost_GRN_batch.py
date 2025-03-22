#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import mean_absolute_error,make_scorer,r2_score
from sklearn.model_selection import KFold,train_test_split
from xgboost import XGBRegressor
import xgboost as xgb
import numpy as np


data = pd.read_table("TPM.txt",header=0,index_col=0) 
data = data.T
print(data.shape) 



TFs = pd.read_table("Gm.TFs.txt",header=None,index_col=0)
TFs.head()
TFs = TFs[1].values.tolist()
print(TFs[1:10])


xgb_params={
     'seed':0,
    'eta':0.1,
    'colsample_bytree':0.5,
    'subsample':0.5,
    'objective':'reg:linear',  
    'max_depth':3,
    'min_child_weight':3
}
n_estimator=100 

target_genes = list(data.columns)
for m in range(0,len(target_genes)):
    target_gene = target_genes[m]
    print(m)
    train_y = data[target_gene]
    train_x = data[list(set(list(data)) & set(TFs))] # using all TFs as train_x fatures 
    if target_gene in train_x.columns.values:
        train_x = train_x.drop(target_gene, 1) 
        
    # make training matrix
    dtrain=xgb.DMatrix(train_x,train_y) 
    # 模型训练和预测
    xgb_model = xgb.train(xgb_params, dtrain, 100)    
    ## output
    imp = pd.DataFrame()
    imp["TF"] = xgb_model.get_fscore().keys()
    imp["Fscore"] = xgb_model.get_fscore().values()
    imp = imp.sort_values(["Fscore"],ascending=False)[0:10]
    imp["Target"] = target_gene
    imp["Num"] = m
    imp.to_csv("TF-target.boost.results.csv",header=False,mode='a',index=False)





