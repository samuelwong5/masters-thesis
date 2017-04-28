import numpy as np
import csv

fp = '../data/E-GEOD-48350/E-GEOD-48350-combined.csv'

with open(fp, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    headers = next(reader, None)
    data = [list(map(float, x[5:])) for x in reader]

cov_matrix = np.cov(data, rowvar=False)
np.fill_diagonal(cov_matrix, 0)
cov_matrix = np.absolute(cov_matrix)
print(np.amax(cov_matrix, axis=0)[:10])

