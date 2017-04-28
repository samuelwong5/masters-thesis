import numpy as np
from sklearn.lda import LDA
from sklearn.metrics import accuracy_score

from util import DataReader, partition_data


fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'


x, y = DataReader(fp).get_data()
argmax = lambda x: x[0] if x[0] == 1 else x[1]
y = list(map(argmax, y))
partition = partition_data(x, y, [0.8, 0.2])

mli = lambda x: np.array(x).astype(float)
train_x = mli(partition[0][0])
train_y = mli(partition[0][1])
test_x = mli(partition[1][0])
test_y = mli(partition[1][1])

lda = LDA(n_components=2, shrinkage='auto', solver='lsqr')
lda.fit(train_x, train_y)
test_y_pred = lda.predict(test_x)
print(accuracy_score(test_y_pred, test_y))