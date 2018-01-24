"""
Thanks
------
https://kevinzakka.github.io/2016/07/13/k-nearest-neighbor/
http://scott.fortmann-roe.com/docs/BiasVariance.html
"""

import numpy as np
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score
import pdb

from params import par
for key,val in par.items():
    exec(key + '=val')

# make a list of possible N_want values that are odd only
Ks  =   np.arange(0,N/2)
Ks  =   2*Ks + 1

def find_optimal_KNN(positions):

    x_train,x_test,y_train,y_test   =   train_test_split( positions.T , np.arange(N) , test_size=0.33 , random_state=42 )

    # empty list that will hold cross validation scores
    cv_scores   =   np.zeros(len(Ks))

    # perform n-fold cross validation
    for i,k in enumerate(Ks):
        knn             =   KNeighborsClassifier(n_neighbors=k)
        # pdb.set_trace()
        scores          =   cross_val_score(knn, x_train, y_train, cv=10, scoring='accuracy')
        cv_scores[i]    =   scores.mean()

    # changing to misclassification error
    MSE =   1 - cv_scores

    # find index that returns smallest misclassification error
    Ki  =   ( np.abs( Ks - np.min(MSE) ) ).argmin()

    return Ks[Ki]
