""" Pickle2Mat
This code transforms the Python pickles into importable mat files

This code is optional since step 2 now also construct the mat files.
"""

import pickle
import os
import scipy.io as sio
import numpy as np

os.chdir('..')

for fold in ['temporal_citation', 'temporal_coauthor']:
    print("Starting directory {}. This might take time.".format(fold))
    mat_fold = 'mat_{}'.format(fold.split('_')[-1])
    if not os.path.isdir(mat_fold):
        os.mkdir(mat_fold)

    for fil in os.listdir(fold):
        print("Working on {}".format(fil))
        pik = pickle.load(open(os.path.join(fold, fil), 'rb'))
        i, j, data = [], [], pik['W'].data.tolist()
        for row in xrange(len(pik['W'].indptr)-1):
            for idx in xrange(pik['W'].indptr[row], pik['W'].indptr[row+1]):
                i.append(row+1)
                j.append(pik['W'].indices[idx]+1)

        path = os.path.join(mat_fold, '{}.mat'.format(fil))
        mat_data = {'i': i, 'j': j, 'data': data, 'active_authors': np.array(map(lambda x: x+1, pik['active_authors'].nonzero()[0]))}
        sio.savemat(path, mat_data)
