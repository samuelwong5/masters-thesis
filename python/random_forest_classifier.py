import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV

from util import DataReader, partition_data


fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'


x, y = DataReader(fp).get_data()
partition = partition_data(x, y, [0.8, 0.2])

mli = lambda x: np.array(x).astype(float)
train_x = mli(partition[0][0])
train_y = mli(partition[0][1])
test_x = mli(partition[1][0])
test_y = mli(partition[1][1])

param_grid = {
    'n_estimators': [250, 500, 1000, 2000],
    'max_features': ['auto', 'log2']
}

rf = RandomForestClassifier(n_estimators=500)
rf_CV = GridSearchCV(estimator=rf, param_grid=param_grid, cv=5)
rf_CV.fit(mli(train_x), mli(train_y))
print(rf_CV.best_params_)

rf = RandomForestClassifier(n_estimators=rf_CV.best_params_['n_estimators'],
                            max_features=rf_CV.best_params_['max_features'])
rf.fit(train_x, train_y)
test_y_pred = rf.predict(test_x)
print(accuracy_score(test_y_pred, test_y))