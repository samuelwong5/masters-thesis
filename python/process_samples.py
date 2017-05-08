from os import listdir
from os.path import isfile, join
import csv
from functools import reduce

folder = '../data/E-GEOD-84422/samples/'
files = [join(folder,f) for f in listdir(folder) if isfile(join(folder, f))]

ifiles = []

print('Reading files...')
for f in files:
    reader = csv.reader(open(f, 'r'), delimiter='\t')
    next(reader, None)
    d = {}
    for row in reader:
        d[row[0]] = row[1]
    ifiles.append(d)


nodes = reduce(set.intersection, [f.keys() for f in ifiles])
print(len(nodes))