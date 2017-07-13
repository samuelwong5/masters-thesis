import numpy as np
from pathlib import Path
from sklearn import linear_model

from util import DataReader


dataset_folder = Path('../data')
dataset_names = [
    'E-GEOD-48350',
    'E-GEOD-84422',
    'E-GEOD-63063',
    'E-GEOD-63063-2',
]

datasets = list(map(lambda x: str(dataset_folder / x / '{}{}'.format(x, '-combined-data.csv')), dataset_names))

print('Reading file')
#dataset = datasets[3]
dataset = '../data/combine-48350-84422.csv'
x, y = DataReader(dataset).get_data()

print(len(x), len(x[0]))
print(len(y))
print(sum(y))

x = np.array(x)
for i in range(388):
    print(x[i,89], y[i])
y = np.array(y)
print('Fitting model')
clf = linear_model.Lasso(alpha=1.0)
clf.fit(x, y)
print(clf.coef_[1:10])
for i, x in enumerate(clf.coef_):
    if x > 0:
        print(i, x)

print(sum(clf.coef_))