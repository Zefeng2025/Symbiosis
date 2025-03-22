#!/usr/bin/env python
# coding: utf-8


import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import mean_absolute_error,make_scorer,r2_score
from sklearn.model_selection import KFold,train_test_split
from xgboost import XGBRegressor
import xgboost as xgb
import numpy as np


data = pd.read_table("20220901_MtExpressV3-Dataset/data/gene/mtexpress_v3.tmm.matrix.gene.tsv",header=0,index_col=0) 
data = data.T
print(data.shape) # 1942个样品 X 50343个基因.
#print(data.head)

# 读取TF列表文件
gene_type = pd.read_table("20220901_MtExpressV3-Dataset/annotation/MtrunA17r5.0-ANR.5.1.9.annotation.tsv",header=0,index_col=0) 
#print(gene_type.shape)
TFs = gene_type.index[gene_type["regulatorType"]=="transcription factor"]
TFs = list(TFs)
print(len(TFs),TFs[0:3]) # 3232 TFs

# 参数设置
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
for m in range(43488,len(target_genes)):
    target_gene = target_genes[m]
    print(m)
    train_y = data[target_gene]
    train_x = data[list(set(list(data)) & set(TFs))] # using all TFs as train_x fatures 
    if target_gene in train_x.columns.values:
        train_x = train_x.drop(target_gene, 1) # 如果预测的基因是TF，则从train_x去除
        
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
    imp.to_csv("/home/wuzefeng/MyResearch/1Gannong/6Medicago_reg/TF-target.boost.results.csv",header=False,mode='a',index=False)
        




