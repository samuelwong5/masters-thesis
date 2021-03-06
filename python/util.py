'''
Utility classes and functions
'''
import csv
from random import shuffle
import numpy as np

def str_to_float(s):
    try:
        return float(s)
    except ValueError: # Handle malformed data
        return 0.0

class DataReader():
    '''
    Utility class to parse RNA array data.

    Input data is a CSV-formatted file, with the last
    column being the targets.
    '''
    def __init__(self, fp, delimiter=','):
        self.reader = csv.reader(open(fp, 'r'), delimiter=delimiter)
        self.headers = next(self.reader, None)
        self.x = []
        self.y = []
        for row in self.reader:
            self.x.append([str_to_float(x) for x in row[2:-1]])
            self.y.append(int(row[-1]))

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
        self.size = len(x)
        self.counter = 0       # Current fold being held out
        self.pointer = 0       # Pointer to non-retrieved data if batch size != (k-1)/k * size
        self.fold = fold       # Number of folds
        self.fold_size = int(self.size / fold)
        self.x = np.array(x[:self.fold_size * fold])
        self.y = np.array(y[:self.fold_size * fold])
        self.size = self.fold_size * fold
        self.shuffle_data()


    def shuffle_data(self):
        s = list(zip(self.x, self.y))
        shuffle(s)
        s = list(zip(*s))
        self.x = np.array(list(s[0]))
        self.y = np.array(list(s[1]))

        self.counter = 0
        self.pointer = 0


    def next_batch(self, batch_size):
        if self.pointer >= self.size:
            self.pointer = 0
            self.counter += 1
        crr = self.pointer
        nxt = self.pointer + batch_size

        if self.counter == 0 and self.pointer == 0:
            self.pointer += self.fold_size
            crr = self.pointer
            nxt = self.pointer + batch_size
            batch_x = self.x[crr:nxt]
            batch_y = self.y[crr:nxt]
        #elif self.counter * self.fold_size <= nxt < (self.counter + 1) * self.fold_size:
        #    nxt += self.fold_size
        #    batch_x = self.x[crr:self.counter * batch_size] + self.x[(self.counter + 1) * batch_size:nxt]
        #    batch_y = self.y[crr:self.counter * batch_size] + self.y[(self.counter + 1) * batch_size:nxt]
        else:
            batch_x = self.x[crr:nxt]
            batch_y = self.y[crr:nxt]
        self.pointer = nxt
        return batch_x, batch_y


    def get_train_set(self):
        if self.counter == self.fold:
            self.shuffle_data()

        train_x = np.concatenate([self.x[:self.counter * self.fold_size], self.x[(self.counter + 1) * self.fold_size:]])
        train_y = np.concatenate([self.y[:self.counter * self.fold_size], self.y[(self.counter + 1) * self.fold_size:]])
        return train_x, train_y


    def get_valid_set(self):
        valid_x = self.x[self.counter * self.fold_size:(self.counter + 1) * self.fold_size]
        valid_y = self.y[self.counter * self.fold_size:(self.counter + 1) * self.fold_size]
        self.counter += 1
        return valid_x, valid_y


def partition_data(x, y, pct):
    '''
    Partitions data into sets weighted according to pct
        pct - a list of floats that represent the proportions of partitions
    '''
    assert sum(pct) == 1, '[ERROR] Sum of data set partition percentages is not 1.'
    assert len(x) == len(y), '[ERROR] Length of features and targets vectors is not equal.'
    s = list(zip(x, y))
    shuffle(s)
    s = list(zip(*s))
    x = s[0]
    y = s[1]
    res = []
    size = len(x)
    curr = 0.0
    for i in pct:
        res.append((x[int(size*curr):int(size*(curr+i))],
                    y[int(size*curr):int(size*(curr+i))]))
        curr += i
    return res
