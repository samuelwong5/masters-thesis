'''
Utility classes and functions
'''
import csv
from random import shuffle


class DataReader():
    '''
    Utility class to parse RNA array data.

    Input data is a CSV-formatted file, with the column
    headers:
        ''     - patient code
        gender - converted to 0 for female, 1 for male
        ad     - targets; 0 for healthy, >1 for braak scale
    '''
    def __init__(self, fp, delimiter=','):
        self.reader = csv.reader(open(fp, 'r'), delimiter=delimiter)
        self.headers = next(self.reader, None)
        self.fnc = {}
        for i, h in enumerate(self.headers[:10]):
            if h == '' or h == 'brainRegion' or h == 'ad':
                pass
            elif h == 'gender':
                self.fnc[i] = lambda x: 0.0 if x == 'female' else 1.0
            else:
                self.fnc[i] = lambda x: float(x)
        self.x = []
        self.y = []
        target_index = self.headers.index('ad')
        for row in self.reader:
            self.x.append(self.preprocess(row))
            if int(row[target_index]) > 0:
                self.y.append([0, 1])
            else:
                self.y.append([1, 0])

    def preprocess(self, entry):
        res = []
        for i, x in enumerate(entry[:10]):
            try:
                res.append(self.fnc[i](x))
            except:
                pass
        res += entry[10:]
        return res

    def get_features(self):
        return self.x

    def get_targets(self):
        return self.y

    def get_data(self):
        return self.x, self.y


class CrossValidation():
    '''
    Performs k-fold validation by shuffling datasets and
    returning appropriate batches for training and validation.
    '''
    def __init__(self, x, y, fold=10):
        assert len(x) == len(y), 'Length of input and targets are different!'
        self.x = x
        self.y = y
        self.fold = fold
        self.batch_size = int(len(x) / fold)
        self.shuffle_data()

    def shuffle_data(self):
        shuffle(self.x)
        shuffle(self.y)
        self.counter = 0

    def next_batch(self):
        if self.counter == self.fold:
            self.shuffle_data()
        batch_x = self.x[self.counter * self.batch_size:(self.counter + 1) * self.batch_size]
        batch_y = self.y[self.counter * self.batch_size:(self.counter + 1) * self.batch_size]
        self.counter += 1
        return batch_x, batch_y


def partition_data(x, y, pct):
    '''
    Partitions data into training/validation/test sets
    '''
    assert sum(pct) == 1, '[ERROR] Sum of data set partition percentages is not 1.'
    assert len(x) == len(y), '[ERROR] Length of features and targets vectors is not equal.'
    shuffle(x)
    shuffle(y)
    res = []
    size = len(x)
    curr = 0.0
    for i in pct:
        res.append((x[int(size*curr):int(size*(curr+i))],
                    y[int(size*curr):int(size*(curr+i))]))
        curr += i
    return res
