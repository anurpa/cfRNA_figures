#Random Forest cross validation for LVQ top 10 genes 

#Import modules
from sklearn import datasets
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RepeatedKFold
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelBinarizer
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import RandomizedSearchCV
import random as r

def cfrna_rf_kfold(data,n,sample_type):
    lb = LabelBinarizer()
    y=lb.fit_transform(data.iloc[:,0])
    data.iloc[:,0]
    X=data.iloc[:,1:].values
    
    # Number of trees in random forest
    n_estimators = [int(x) for x in np.linspace(start = 10, stop = 100, num = 10)]

    # Number of features to consider at every split
    max_features = ['auto', 'sqrt']
    # Maximum number of levels in tree
    max_depth = [int(x) for x in np.linspace(10, 110, num = 11)]
    max_depth.append(None)
    # Minimum number of samples required to split a node

    min_samples_split = np.arange(1,11)
    # Minimum number of samples required at each leaf node
    min_samples_leaf = [1, 2, 4]
    # Method of selecting samples for training each tree
    bootstrap = [True, False]
    # Create the random grid
    random_grid = {'n_estimators': n_estimators,
                   'max_features': max_features,
                   'max_depth': max_depth,
                   'min_samples_split': min_samples_split,
                   'min_samples_leaf': min_samples_leaf,
                   'bootstrap': bootstrap}

    # Use the random grid to search for best hyperparameters
    # First create the base model to tune
    rf = RandomForestClassifier()

    # Set the number of fold to split the data into 
    #10% Testing
    n_splits=n

    # Set the number of times to repeat the KFold procedure
    n_repeats = 5

    # Create the cross-validation generator
    rkf = RepeatedKFold(n_splits = n_splits, n_repeats = n_repeats)

    # Random search of parameters, using 3 fold cross validation,
    # search across 20 different combinations, and use all available cores
    rf_random = RandomizedSearchCV(estimator = rf, param_distributions = random_grid, n_iter = 20, cv = rkf, verbose=2, random_state=42, n_jobs = -1)
    # Fit the random search model
    rf_random.fit(X, y)

    classifier = rf_random.best_estimator_
    
###############################################################################
    # Classification and ROC analysis

    # Run classifier with cross-validation and plot ROC curves
    cv = StratifiedKFold(n_splits=n)

    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    i = 0
    plt.figure(figsize=(15,15))
    for train, test in cv.split(X, y):
        probas_ = classifier.fit(X[train], y[train]).predict_proba(X[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y[test], probas_[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,linewidth=3.0,
                 label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

        i += 1
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Chance', alpha=.8,linewidth=3.0)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8,linewidth=3.0)
    plt.rcParams['font.family'] = "arial"
    plt.xlabel('xlabel', fontsize=60)
    plt.ylabel('ylabel', fontsize=60)
    plt.rcParams.update({'legend.fontsize':22})
    ax=plt.gca()
    ax.tick_params(axis = 'both', which = 'major', labelsize = 30)
    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')


    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('%s , %s-fold, ROC'%(sample_type,n),fontsize=60)
    #plt.legend(loc="lower right", prop={'size': 40})
    plt.legend(loc="lower right")
    
    plt.savefig("../figures/cv/ROC_%s_%sfold_v2.pdf" %(sample_type,n),dpi=500)
