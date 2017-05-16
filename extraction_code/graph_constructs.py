"""Python GraphConstructs
Illustrate examples of graph constructs from pickles

"""

import os
import pickle
import pygsp
import numpy as np
from scipy import sparse as sp

os.chdir('..')

# One month graph construction

data = pickle.load(open('temporal_citation/1992-02', 'rb'))
auth_idx = data['active_authors'].nonzero()[0]
W = data['W'][np.ix_(auth_idx, auth_idx)]
G = pygsp.graphs.Graph(W)
G.set_coords(iterations=10)
G.plot(show_edges=True)

# Aggregation of several months

N = 366572
W = sp.csr_matrix((N, N))
# active_authors = sp.csr_matrix((N, 1), dtype='bool')

from_year = 2002
to_year = 2004
from_month = 01
to_month = 12

for i in xrange(from_year, to_year + 1):
    start_month = from_month if i == from_year else 1
    end_month = to_month if i == to_year else 12
    for j in xrange(start_month, end_month + 1):
        print("Read {}-{:2d}".format(i, j))
        cit_tmp = pickle.load(open('temporal_citation/{}-{:02d}.pickle'.format(i, j), 'rb'))
        # active_authors = cit_tmp['active_authors'].tocsr() + active_authors
        # diag = G['W'].diagonal().nonzero()[0] -- dont know why it served
        W += cit_tmp['W']

# auth_idx = active_authors.nonzero()[0]
# W_red = W[np.ix_(auth_idx, auth_idx)] -- keep only connected nodes
G = pygsp.graphs.Graph(W, gtype='Citation_{}-{}:{}-{}'.format(from_year, from_month, to_year, to_month), lap_type='normalized')
# W = (W + W.T) / 2. -- only if citation
G.set_coords(iterations=10)
G.plot(show_edges=True)

# Adding a month and removing the first one

print("Remove {}-{:2d}".format(from_year, from_month))
cit_tmp = pickle.load(open('temporal_citation/{}-{:02d}.pickle'.format(from_year, from_month), 'rb'))
W -= cit_tmp['W']

from_month += 1
to_month += 1
if from_month == 13:
    from_month = 1
    from_year += 1

if to_month == 13:
    to_month = 1
    to_year += 1

print("Add {}-{:2d}".format(to_year, to_month))
cit_tmp = pickle.load(open('temporal_citation/{}-{:02d}.pickle'.format(to_year, to_month), 'rb'))
W += cit_tmp['W']
