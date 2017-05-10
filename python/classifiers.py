import numpy as np
from random import shuffle
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC
from sklearn.metrics import roc_auc_score
import sys

from util import DataReader, partition_data, CrossValidation

try:
    classifier = sys.argv[1]
    assert classifier in ['--rfc', '--svm'], \
        'Classifier must be --rfc or --svm'
except:
    print('Missing classifier argument')
    sys.exit(-1)

iterations = 1 if len(sys.argv) < 3 else int(sys.argv[2])


#fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'
fp = '../data/E-GEOD-84422/E-GEOD-84422-combined2.csv'
#fp = '../data/E-GEOD-63063/E-GEOD-63063-combined.csv'

print('Reading file')
x, y = DataReader(fp).get_data()
argmax = lambda x: x.index(max(x))
y = list(map(argmax, y)) # Convert labels to 1-d

def part(x, y):
    partition = partition_data(x, y, [0.8, 0.2])
    mli = lambda x: np.array(x)
    train_x = mli(partition[0][0])
    train_y = mli(partition[0][1])
    test_x = mli(partition[1][0])
    test_y = mli(partition[1][1])
    return train_x, train_y, test_x, test_y

mean_acc = 0.0
folds = 10
data = CrossValidation(x, y, folds)

for i in range(folds):
    print('[Iteration {0:2d}]'.format(i+1))
    print(' - Shuffling data')
    train_x, train_y = data.get_train_set()
    test_x, test_y = data.get_valid_set()
    classifiers = {
        '--rfc': RandomForestClassifier(n_estimators=2000, max_features='log2'),
        '--svm': SVC()
    }
    try:
        cfr = classifiers[classifier]
    except KeyError:
        print(' - Unknown classifier: {0}'.format(classifier))
        sys.exit(-1)

    print(' - Fitting classifier')
    cfr.fit(train_x, train_y)

    test_y_pred = cfr.predict(test_x)
    acc_score = accuracy_score(test_y_pred, test_y)
    print(' - Test-set accuracy: {:5f}'.format(acc_score))
    print(' - AUC: {:5f}'.format(roc_auc_score(test_y_pred, test_y)))
    mean_acc += acc_score
mean_acc /= folds
print('Average accuracy: {:5f}'.format(mean_acc))