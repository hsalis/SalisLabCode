"""
Machine learning tools for model improvement
*In development*

Copyright 2017 Alexander C. Reis, Howard M. Salis, all rights reserved.

"""

import numpy as np
import pandas as pd
import sklearn
import scipy
import scipy.stats
from itertools import combinations
from synbiomts import stats

# Feature selection
# http://scikit-learn.org/stable/modules/feature_selection.html#feature-selection
from sklearn.feature_selection import SelectKBest # select k best features
from sklearn.feature_selection import chi2 # use chi-squared to find n best
from sklearn.feature_selection import RFE # recursive feature eliminiation
from sklearn.feature_selection import f_regression

# Remove features with low variance (i.e. low functional diversity)
# http://scikit-learn.org/stable/modules/feature_selection.html#removing-features-with-low-variance
from sklearn.feature_selection import VarianceThreshold

# Pipeline
# http://scikit-learn.org/stable/modules/generated/sklearn.pipeline.Pipeline.html#sklearn.pipeline.Pipeline
from sklearn.pipeline import Pipeline

# Principle Component Analysis
# http://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html
from sklearn.decomposition import PCA


def ANOVA(X,y):
    '''Univariate linear regression tests
    Quick linear model for sequentially testing the effect of many regressors
    Using scikit learn's Feature selection toolbox
    Returns:
        F (array) = F-values for regressors
        pvalues (array) = p-values for F-scores'''

    (F,pvalues) = f_regression(X,y)
    return (F,pvalues)

fx = lambda m,x,b: m*x+b

def feature_reduction(X,y):
    #Where X is a pandas dataframe slice of the model features/predictors
    #y can be a numpy array of outcomes or pandas series

    assert isinstance(X,pd.DataFrame), "Make sure you provide X as a pandas dataframe"
    y = np.array(y)
    features = list(X)
    n = len(features)

    # Assess complete model
    x = X.sum(axis=1).values
    (m,b) = stats.fit_linear_model(x,y)
    y_predicted = fx(m,x,b)
    e1 = y - y_predicted

    reduced_models = {}

    # Recursive feature eliminiation-like approach

    check = True
    while check:

        n -= 1
        check_list = False

        for comb in combinations(features,n):

            temp = X.filter(list(comb))
            x = temp.sum(axis=1).values
            (m,b) = stats.fit_linear_model(x,y)
            e2 = fx(m,x,b)

            (h,F,pvalue) = stats.vartest2(e1,e2,alpha=0.05,test="F")

            if not h:
                reduced_models[comb] = (F,pvalue)
                check_list += True

        # Terminate if no reduced models with n features
        # are statistically indistinguishable from the complete model
        if not check_list:
            check = False

    return reduced_models





