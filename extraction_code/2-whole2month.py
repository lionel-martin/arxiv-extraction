""" Whole2Month

Transform whole pickles containing data into two different graphs.
Split the information on a monthly basis (smallest piece of temporal information given).

This is step 2 for the transformation of the dataset.

"""

import os
import pickle
from collections import defaultdict
from scipy import sparse as sp
import scipy.io as sio


def def_dict_rec():
    return defaultdict(def_dict_rec)

os.chdir('..')

temporal_links = pickle.load(open(os.path.join('pickles', 'temporal_links.pickle'), 'rb'))
authors = pickle.load(open(os.path.join('pickles', 'authors.pickle'), 'rb'))

# Citation network : connect author to its references (directed)

t_graphs_citation = {}
make_dirs = ['temporal_citation', 'mat_citation']
for fold in make_dirs:
    if not os.path.isdir(fold):
        print("Creating directory '{}'".format(fold))
        os.mkdir(fold)

for y, tmp in temporal_links.items():
    for m, l in tmp.items():
        key = '{}{}-{}'.format(19 if y[0] == '9' else 20, y, m)
        print("Processing {}".format(key))
        active_authors = sp.lil_matrix((len(authors), 1), dtype='bool')

        data, row, col = [], [], []
        for p in l:
            for w in p[0]:
                active_authors[w] = True
                row += [w] * len(p[1])
                for r in p[1]:
                    active_authors[r] = True
                    col.append(r)
                    data.append(1)

        W = sp.csr_matrix((data, (row, col)), shape=(len(authors), len(authors)))
        t_graphs_citation[key] = {'W': W, 'active_authors': active_authors}
        pickle.dump(t_graphs_citation[key], open('temporal_citation/{}.pickle'.format(key), 'wb'))

        path = os.path.join('mat_citation', '{}.mat'.format(key))
        mat_data = {'i': row, 'j': col, 'data': data, 'active_authors': np.array(map(lambda x: x+1, active_authors.nonzero()[0]))}
        sio.savemat(path, mat_data)

# Coauthor network : connect authors that coauthored (undirected)

t_graphs_coauthor = {}
if not os.path.isdir('temporal_coauthor'):
    os.mkdir('temporal_coauthor')
if not os.path.isdir('mat_coauthor'):
    os.mkdir('mat_coauthor')

for y, tmp in temporal_links.items():
    for m, l in tmp.items():
        key = '{}{}-{}'.format(19 if y[0] == '9' else 20, y, m)
        print("Processing {}".format(key))
        active_authors = sp.lil_matrix((len(authors), 1), dtype='bool')

        data, row, col = [], [], []
        for p in l:
            for i in xrange(len(p[0])):
                w = p[0][i]
                activate = False
                for j in xrange(i+1, len(p[0])):
                    if p[0][j] != w:
                        row.append(w)
                        col.append(p[0][j])
                        data.append(1)
                        activate = True
                if activate:
                    active_authors[w] = True

        W = sp.csr_matrix((data, (row, col)), shape=(len(authors), len(authors)))
        W += W.T
        t_graphs_coauthor[key] = {'W': W, 'active_authors': active_authors}
        pickle.dump(t_graphs_coauthor[key], open('temporal_coauthor/{}.pickle'.format(key), 'wb'))

        path = os.path.join('mat_coauthor', '{}.mat'.format(key))
        mat_data = {'i': row, 'j': col, 'data': data, 'active_authors': np.array(map(lambda x: x+1, active_authors.nonzero()[0]))}
        sio.savemat(path, mat_data)
