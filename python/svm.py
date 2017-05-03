import numpy as np
from sklearn import svm
from sklearn.metrics import accuracy_score

from util import DataReader, partition_data


fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'

print('Reading file')
x, y = DataReader(fp).get_data()
argmax = lambda x: x.index(max(x))
y = list(map(argmax, y))
partition = partition_data(x, y, [0.8, 0.2])

print('Partitioning data')
#mli = lambda x: np.array(x).astype(float)
mli = lambda x: np.array(x)
train_x = mli(partition[0][0])
train_y = mli(partition[0][1])
test_x = mli(partition[1][0])
test_y = mli(partition[1][1])

print('Fitting classifier')
classifier = svm.SVC()
classifier.fit(train_x, train_y)

print('Predicting with classifier')
test_y_pred = classifier.predict(test_x)
print(accuracy_score(test_y_pred, test_y))