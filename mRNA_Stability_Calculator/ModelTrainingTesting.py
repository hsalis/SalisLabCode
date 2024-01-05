#!/usr/bin/env python
# coding: utf-8

import math, os, random, sys, pickle, copy, json, shutil
from itertools import product
sys.path.append('../')
sys.path.append('/usr/local/lib/')
sys.path.append('/usr/local/lib/python3.10/site-packages/')

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import lightgbm as lgb
from lightgbm import plot_importance
from sklearn.metrics import r2_score

print('LGB version: ', lgb.__version__)

NUM_TOP_ISOFORMS = 5

def loadData():
    shutil.unpack_archive('TrainingDataset.csv.zip')
    df_train = pd.read_csv('TrainingDataset.csv')

    shutil.unpack_archive('TestDataset.csv.zip')
    df_test = pd.read_csv('TestDataset.csv')

    df_train = categorizeFeatures(df_train)
    df_test = categorizeFeatures(df_test)

    stats = {}
    stats['test_data'] = df_test.describe().to_json()
    stats['train_data'] = df_train.describe().to_json()

    return (df_train, df_test, stats)

def categorizeFeatures(df):

    categoricalFeatureList = ['rppH','GQUAD','iMOTIF','CsrA']

    for i in range(NUM_TOP_ISOFORMS):
        isoform = lambda feature: '{}_iso_{}'.format(feature,i)

        for feature in categoricalFeatureList:
            df[isoform(feature)] = df[isoform(feature)].astype('category')

        for j in range(80):
            seq_feature = isoform('nt_{}'.format(j))
            df[seq_feature] = df[seq_feature].astype('category')

    return df

def r2_score_without_nan(xlist, ylist):

    xnew = []
    ynew = []
    for (x, y) in zip(xlist, ylist):
        if np.isnan(x) or np.isnan(y):
            continue
        else:
            xnew.append(x)
            ynew.append(y)

    return r2_score(xnew, ynew)

def returnKeepColumns(decayModel = False):

    keepColumns = []

    for i in range(NUM_TOP_ISOFORMS):
        isoform = lambda feature: '{}_iso_{}'.format(feature,i)

        features = []
        features += ['rppH','TX_rate','TL_rate']
        features += ['initial_sum_external_ssRNA', 'initial_dsRNA_sum']
        features += ['initial_external_ssRNA_A','initial_external_ssRNA_T','initial_external_ssRNA_G','initial_external_ssRNA_C']
        features += ['GQUAD', 'iMOTIF', 'CsrA']
        features += ['nt_{}'.format(i) for i in range(80)]
        for feature in features:
            keepColumns += [isoform(feature)]

    if decayModel:
        keepColumns += ['predicted_T0','predicted_T2','predicted_T4','predicted_T8', 'predicted_T16']
    return keepColumns

def trainTimepointModels(df_train, df_test, outcome = 'T0', modelName = 'model_T0'):

    outputData = {}
    outcome_feature = 'RNA_{}'.format(outcome)

    df_outcomes = pd.read_excel('SupplementaryData1.xlsx')

    y_train = []
    idxList = df_train['Index']
    for idx in idxList:
        DNA_counts = df_outcomes.iloc[idx-1]['DNA']
        RNA_counts = df_outcomes.iloc[idx-1][outcome_feature]

        if DNA_counts > 0.0 and RNA_counts > 0.0:
            y_train.append( np.log(RNA_counts / DNA_counts) )
        else:
            y_train.append(np.nan)
    df_train[outcome_feature] = y_train

    y_test = []
    idxList = df_test['Index']
    for idx in idxList:
        DNA_counts = df_outcomes.iloc[idx-1]['DNA']
        RNA_counts = df_outcomes.iloc[idx-1][outcome_feature]

        if DNA_counts > 0.0 and RNA_counts > 0.0:
            y_test.append( np.log(RNA_counts / DNA_counts) )
        else:
            y_test.append(np.nan)
    df_test[outcome_feature] = y_test

    keepColumns = returnKeepColumns()
    x_train = df_train[df_train.columns[df_train.columns.isin(keepColumns)]]
    x_test = df_test[df_test.columns[df_test.columns.isin(keepColumns)]]

    hyper_params = {
        'task'             : 'train',
        'boosting_type'    : 'gbdt',
        "linear_tree"      : False,
        "force_row_wise"   : False,
        "force_col_wise"   : True,
        'objective'        : 'l2',
        'metric'           : 'l2',
        "num_leaves"       : 100,
        "max_depth"        : 5,
        "max_bin"          : 1000,
        "n_estimators"     : 119,
        'learning_rate'    : 0.1,
        'verbose'          : 1,
        'bagging_fraction': 0.50,
        'bagging_freq': 5,
        'feature_fraction': 0.25,
        'min_data_in_leaf' : 50,
        'importance_type' : 'gain'
    }

    gbm = lgb.LGBMRegressor(**hyper_params)
    gbm.fit(x_train, y_train, eval_set=[(x_test, y_test),(x_train, y_train)])
    saved_model = gbm.booster_.save_model('{}_saved.txt'.format(modelName))

    lgb.plot_metric(gbm)
    plt.savefig('{}_learning_curve.eps'.format(modelName))

    pred_train = gbm.predict(x_train, num_iteration=gbm.best_iteration_, pred_contrib=False)
    pred_test = gbm.predict(x_test, num_iteration=gbm.best_iteration_)
    outputData['numDatapoints'] = (len(y_train), len(y_test))
    outputData['train_R2'] = r2_score_without_nan(y_train, pred_train)
    outputData['test_R2'] = r2_score_without_nan(y_test, pred_test)

    prediction_feature = 'predicted_{}'.format(outcome)
    df_train[prediction_feature] = pred_train
    df_test[prediction_feature] = pred_test

    total_pred = np.concatenate((pred_test,pred_train))
    total_actual = np.concatenate((y_test,y_train))

    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        np.exp(total_pred),
        np.exp(total_actual),
        alpha = 0.2,
        s=1
    )
    ax.set_xlabel('Predicted mRNA Levels')
    ax.set_ylabel('Measured mRNA Levels')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(0.03, 30)
    ax.set_xlim(0.03, 30)
    ax.tick_params(labelleft=True, labelright=True)

    plt.savefig('{}_predictions_all_data.eps'.format(modelName),dpi=300)
    plt.close()

    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        np.exp(pred_test),
        np.exp(y_test),
        alpha = 0.2,
        s=1)
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(0.03, 30)
    ax.set_xlim(0.03, 30)
    ax.tick_params(labelleft=True, labelright=True)
    ax.set_xlabel('Predicted mRNA Decay')
    ax.set_ylabel('Measured mRNA Levels')
    plt.savefig('{}_predictions_test_data.eps'.format(modelName),dpi=300)
    plt.close()

    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        np.exp(pred_train),
        np.exp(y_train),
        alpha = 0.2,
        s=1)
    ax.set_xlabel('Predicted mRNA Levels')
    ax.set_ylabel('Measured mRNA Levels')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_ylim(0.03, 30)
    ax.set_xlim(0.03, 30)
    ax.tick_params(labelleft=True, labelright=True)
    plt.savefig('{}_predictions_train_data.eps'.format(modelName),dpi=300)
    plt.close()

    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    plt.hist(np.log2(2**np.array(pred_test)/2**np.array(y_test)),bins=100)
    ax.set_xlim(left= -5)
    ax.tick_params(labelleft=True, labelright=True)
    plt.savefig('{}_error_distribution.eps'.format(modelName))
    plt.close()

    plt.rcParams.update({'font.size': 8})
    plot_importance(gbm, figsize=(12,60))
    plt.savefig('{}_feature_importances.eps'.format(modelName))

    return (df_train, df_test, outputData)

def trainDecayModel(df_train, df_test, modelName = 'model_decay_rate'):

    outputData = {}
    df_outcomes = pd.read_excel('SupplementaryData1.xlsx')

    y_train = []
    idxList = df_train['Index']
    for idx in idxList:
        y_train.append(df_outcomes.iloc[idx-1]['k_fit'])
    df_train['k_fit'] = y_train

    y_test = []
    idxList = df_test['Index']
    for idx in idxList:
        y_test.append(df_outcomes.iloc[idx-1]['k_fit'])
    df_test['k_fit'] = y_test

    keepColumns = returnKeepColumns(decayModel = True)
    x_train = df_train[df_train.columns[df_train.columns.isin(keepColumns)]]
    x_test = df_test[df_test.columns[df_test.columns.isin(keepColumns)]]

    hyper_params = {
        'task'             : 'train',
        'boosting_type'    : 'gbdt',
        "linear_tree"      : False,
        "force_row_wise"   : False,
        "force_col_wise"   : True,
        'objective'        : 'l2',
        'metric'           : 'l2',
        "num_leaves"       : 100,
        "max_depth"        : 5,
        "max_bin"          : 1000,
        "n_estimators"     : 119,
        'learning_rate'    : 0.1,
        'verbose'          : 1,
        'bagging_fraction': 0.50,
        'bagging_freq': 5,
        'feature_fraction': 0.25,
        'min_data_in_leaf' : 50,
        'importance_type' : 'gain'
        #'feature_pre_filter' : True,
        #'categorical_features': ",".join(categoricalFeatureList),
    }

    gbm = lgb.LGBMRegressor(**hyper_params)

    gbm.fit(x_train, y_train, eval_set=[(x_test, y_test),(x_train, y_train)])
    lgb.plot_metric(gbm)
    plt.savefig('{}_learning_curve.eps'.format(modelName))

    pred_train = gbm.predict(x_train, num_iteration=gbm.best_iteration_, pred_contrib=False)
    pred_test = gbm.predict(x_test, num_iteration=gbm.best_iteration_)
    outputData['numDatapoints'] = (len(y_train), len(y_test))
    outputData['train_R2'] = r2_score_without_nan(y_train, pred_train)
    outputData['test_R2'] = r2_score_without_nan(y_test, pred_test)

    prediction_feature = 'predicted_k_fit'

    # Save actual and predicted for export
    df_train[prediction_feature] = pred_train
    df_test[prediction_feature] = pred_test

    saved_model = gbm.booster_.save_model('{}_saved.txt'.format(modelName))

    total_pred = np.concatenate((pred_test,pred_train))
    total_actual = np.concatenate((y_test,y_train))

    #Convert to half-life [minutes]
    tau_pred_train = [np.log(2) / x for x in pred_train]
    tau_pred_test = [np.log(2) / x for x in pred_test]

    tau_actual_train = [np.log(2) / x for x in y_train]
    tau_actual_test = [np.log(2) / x for x in y_test]

    tau_total_pred = [np.log(2) / x for x in total_pred]
    tau_total_actual = [np.log(2) / x for x in total_actual]

    ## All Data
    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        total_pred,
        total_actual,
        alpha = 0.2,
        s=1
    )
    ax.set_xlabel('Predicted mRNA Decay')
    ax.set_ylabel('Measured mRNA Decay')
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim(0.0, 1.5)
    ax.set_ylim(0.0, 1.5)
    ax.tick_params(labelleft=True, labelright=True)
    plt.savefig('{}_predictions_all_data.eps'.format(modelName),dpi=300)
    plt.close()

    ## Test Data
    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        pred_test,
        y_test,
        alpha = 0.2,
        s=1)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim(0.0, 1.5)
    ax.set_ylim(0.0, 1.5)
    ax.tick_params(labelleft=True, labelright=True)
    ax.set_xlabel('Predicted mRNA Decay')
    ax.set_ylabel('Measured mRNA Decay')
    plt.savefig('{}_predictions_test_data.eps'.format(modelName),dpi=300)
    plt.close()

    ## Train Data

    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    ax.scatter(
        pred_train,
        y_train,
        alpha = 0.2,
        s=1)
    #ax.set_yscale('log')
    #ax.set_xscale('log')
    ax.set_xlim(0.0, 1.5)
    ax.set_ylim(0.0, 1.5)
    ax.tick_params(labelleft=True, labelright=True)
    ax.set_xlabel('Predicted mRNA Decay')
    ax.set_ylabel('Measured mRNA Decay')
    plt.savefig('{}_predictions_train_data.eps'.format(modelName),dpi=300)
    plt.close()

    plt.rcParams.update({'font.size': 14})
    fig,ax = plt.subplots(figsize=(8, 8))
    plt.hist(pred_test-y_test,bins=100)
    ax.set_xlim(left= -5)
    ax.tick_params(labelleft=True, labelright=True)
    plt.savefig('{}_error_distribution.eps'.format(modelName))
    plt.close()

    plt.rcParams.update({'font.size': 8})
    plot_importance(gbm, figsize=(12,60))
    plt.savefig('{}_feature_importances.eps'.format(modelName))

    return (df_train, df_test, outputData)

if __name__ == "__main__":

    random.seed(1)

    (df_train, df_test, stats) = loadData()
    print('# Datapoints: ', len(df_train), len(df_test) )

    outcomeList = ['T0','T2','T4','T8','T16']
    for outcome in outcomeList:
        modelName = 'model_{}'.format(outcome)
        (df_train, df_test, outputData) = trainTimepointModels(df_train, df_test, outcome = outcome, modelName = modelName)
        stats[modelName] = outputData

    (df_train, df_test, outputData) = trainDecayModel(df_train, df_test, modelName = 'decay_rate')
    stats['decay_model'] = outputData

    df_all = pd.concat([df_train, df_test])
    df_all.to_csv('SupplementaryData2.csv')

    print(stats)
    with open('stats.txt','w') as fc:
        json.dump(stats, fc)